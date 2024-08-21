#! /usr/bin/env wolframscript
(*
Jonathan Urrutia

Maxima and minima algothims:

 - CentralFiniteDifferences
 - GoldenRatioAlgorithm
 - GoldenRatioMaxima
 - GoldenRatioMinima
 - FindExtremaInterval
 - FindExtrema


*)

(************************************
 CentralFiniteDifferences
*************************************)
CentralFiniteDifferences::usage = "CentralFiniteDifferences[{{x1,y1},{x2,y2},...{xn,yn}}, i] calculates the ith numerical derivative with an O(h^4) error.";
CentralFiniteDifferences::DerOrd = "The method for error order `1` has not been implemented yet.";
CentralFiniteDifferences[list_List, order_Integer] /; (
              If[order <= 2,
                True,
                Message[CentralFiniteDifferences::DerOrd, order]; False]
                ) := Module[{x, y, step, coeff},
                      If[order == 0, Return[list]];
                      {x, y} = Transpose[list];
                      x = x[[3 ;; Length[x] - 1]];
                      step = Differences[x]^order;
                      coeff = Which[order == 1, 1./{12, -3/2, 3/2, -12},
                                    order == 2, 1./{-12, 3/4, -2/5, 3/4, -12}];
                      y = If[Mod[order, 2] == 1,
                              coeff*Drop[#, {3}] & /@ Partition[y, 5, 1],
                              coeff*# & /@ Partition[y, 5, 1]];
                      Transpose[{Drop[x, -1], Apply[Plus, y, 2]/step}]
                    ];


(************************************
 GoldenRatioAlgorithm
*************************************)
Options[GoldenRatioAlgorithm] = {Method -> "Maxima"};
GoldenRatioAlgorithm[func_, interval_, tol_, OptionsPattern[]] :=
      Module[{gr =  .5*(Sqrt[5.] + 1.) , a, b, c, d, test, extrema },
  	         test =  Which[ OptionValue[Method] == "Maxima", Greater,
                            OptionValue[Method] == "Minima", Less] ;

            {a, b} = Sort[interval];
  	         c = b - (b - a)/gr;
  	         d = a + (b - a)/gr;
  	         While[Abs[c - d] > tol,
   		             If[test[func[c], func[d]],
                      b = d,
                      a = c];
   		             c = b - (b - a)/gr;
   		             d = a + (b - a)/gr
                   ];
  	         {#, func[#]} &[0.5 *(a + b)]
      ];

(************************************
 GoldenRatioMaxima
*************************************)
GoldenRatioMaxima::usage =  "GoldenRatioMaxima[func, intervar] calculates the maxima within the interval of the one variable function func empoloying the Golden Ratio Sectioning algorithm. It threads over a list of intervals and returns {x, f(x)}.";
GoldenRatioMaxima::interval =  "The interval should be a pair of numbers or a list of pairs.";
Options[GoldenRatioMaxima] = {Tolerance ->  1.*^-3};
GoldenRatioMaxima[func_, interval_List, OptionsPattern[]] := If[
      Length[Dimensions[interval]] == 1,
            If[Length[interval] == 2,
                GoldenRatioAlgorithm[func, interval, OptionValue[Tolerance] , Method -> "Maxima"] ,
                Message[GoldenRatioMaxima::interval]
            ],
  			    If[(Length /@ interval) == ConstantArray[2, Length[interval]],
               Map[GoldenRatioAlgorithm[func, #, OptionValue[Tolerance] ,  Method -> "Maxima"] &, interval],
              Message[GoldenRatioMaxima::interval]
            ]
        ];

(************************************
 GoldenRatioMinima
*************************************)
GoldenRatioMinima::usage =  "GoldenRatioMaxima[func, intervar] calculates the minima within the interval of the one variable function func empoloying the Golden \
                              Ratio Sectioning algorithm. It threads over a list of intervals and returns {x, f(x)}";
GoldenRatioMinima::interval =  "The interval should be a pair of numbers or a list of pairs.";
Options[GoldenRatioMinima] = {Tolerance ->  1.*^-3};
GoldenRatioMinima[func_, interval_List, OptionsPattern[]] :=  If[
      Length[Dimensions[interval]] == 1,
            If[Length[interval] == 2,
                GoldenRatioAlgorithm[func, interval, OptionValue[Tolerance] , Method -> "Minima"] ,
                Message[GoldenRatioMaxima::interval]
            ],
  			    If[(Length /@ interval) == ConstantArray[2, Length[interval]],
               Map[GoldenRatioAlgorithm[func, #, OptionValue[Tolerance] ,  Method -> "Minima"] &, interval],
              Message[GoldenRatioMaxima::interval]
            ]
        ];

(************************************
FindExtremaInterval
  interval = {a,b}
  h is a number
{maxs, mins}= intExtrema[f,{a,b}, h]
************************************)

FindExtremaInterval::usage =  "Find the interval where the maxima or minima of a f : R -> R^n for \
                                each component";
FindExtremaInterval[func_, interval_List, h_] :=
      Module[{test, var, f, df, ddf, min, max, int, len},
            var = Range @@ Flatten[{interval, h}];
            If[(len = Length[func[var[[1]]] ]) == 0,
              f = {func /@ var},
              f = Transpose[func /@  var]
            ];

            df = CentralFiniteDifferences[Transpose[{var, #}], 1] & /@  f;
            ddf =  CentralFiniteDifferences[Transpose[{var, #}], 2] & /@  f;

            min = max = ConstantArray[{}, Length[f]];
            Table[
              int = Partition[#, 2, 1] &@ df[[i, ;; , 1]];
   	          test = Partition[#, 2, 1] &@ df[[i, ;; , 2]];
   	          test = Apply[Times, test, 1];

              If[ddf[[i, #, 2]] < 0.,
                AppendTo[max[[i]], Flatten[int[[#]]]],
      					AppendTo[min[[i]], Flatten[int[[#]]]]
              ] & /@  Flatten[Position[test, n_ /; n < 0.]],
            {i, Length[f]}];

            If[len == 0,
              {min[[1]], max[[1]]},
              {min, max}
            ]
          ];

(************************************
FindExtrema
  interval = {a,b}
  h is a number
{maxs, mins}= intExtrema[f,{a,b}, h]
************************************)
FindExtrema::usage =  "Find the specific extrema values of a function f : R -> R^n for \
                      each component";
FindExtrema::extType = "Must specify a valid type of extrema point.";

Options[FindExtrema] = {Tolerance -> 1*^-3};
FindExtrema[func_, interval_List, h_, case_String, OptionsPattern[]] :=
    Module[{grFunc, min, max, len},
    len = Length[func[interval[[1]]]];
    {min, max} =  FindExtremaInterval[func, interval, h];

    Which[case == "Minima",
              grFunc = GoldenRatioMinima[#1, #2, Tolerance -> OptionValue[Tolerance]] &;
   					  If[len == 0,
                Map[ grFunc[func, # ] &, min],
    					  Table[Map[ grFunc[func[#][[i]] &, #] & , min[[i]]], {i, len}]
              ],
   	      case == "Maxima",
              grFunc = GoldenRatioMaxima[#1, #2, Tolerance -> OptionValue[Tolerance]] &;
   					  If[len == 0,
                Map[ grFunc[func, # ] &, max],
    						Table[Map[ grFunc[func[#][[i]] &, #] & , max[[i]]], {i, len}]
              ],
   	      case == "Both",
            {grFunc = GoldenRatioMinima[#1, #2, Tolerance -> OptionValue[Tolerance]] &;
    					     If[len == 0,
                      Map[ grFunc[func, # ] &, min],
       					      Table[Map[ grFunc[func[#][[i]] &, #] & , min[[i]]], {i, len}]
                    ],
    				grFunc = GoldenRatioMaxima[#1, #2, Tolerance -> OptionValue[Tolerance]] &;
    					   If[len == 0,
                      Map[ grFunc[func, # ] &, max],
                      Table[Map[ grFunc[func[#][[i]] &, #] & , max[[i]]], {i, len}]
                  ]
    				},
   		    True,
            Message[FindExtrema::extType]
        ]
    ];
