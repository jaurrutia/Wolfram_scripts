#! /usr/bin/env wolframscript
(*
Jonathan Urrutia

*)

toMap::usage  = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;

(************************************
 FresnelCoefficient
*************************************)
FresnelCoefficient::angDomain = "The incident angle `1` must be between 0 and Pi/2.";
FresnelCoefficient::refIndeces = "The length of `1` must equal the length of `2`.";

(************************************
 FresnelCoefficientS
 FresnelCoefficientP
*************************************)
FresnelReflectionS::usage = "It returns the Fresnel amplitude coefficient of reflection for S polarization. The results is an array of dimensions Length[ang] x Length[ni]";
FresnelReflectionS[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], costhi, costht, rs},
  	    If[ang < 0 || Pi/2. < ang,
            Return[Message[FresnelCoefficient::angDomain, ang]]
        ];
  	    If[Length[ni] != Length[nt],
            Return[Message[FresnelCoefficient::refIndeces, ni, nt]]
        ];

        costhi = Cos[thetai];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		rs = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {costhi, costht};
  		Which[Length[ni] == Length[ang] == 0, rs[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[rs],
   			  True, rs]
        ]

FresnelReflectionP::usage = "It returns the Fresnel amplitude coefficient of reflection for P polarization. The results is an array of dimensions Length[ang] x Length[ni]";
FresnelReflectionP[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], Ncosthi, costht, rp},
  	    If[ang < 0 || Pi/2. < ang,
            Return[Message[FresnelCoefficient::angDomain, ang]]
        ];
  	    If[Length[ni] != Length[nt],
            Return[Message[FresnelCoefficient::refIndeces, ni, nt]]
        ];

  		Ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		rp = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {Ncosthi, costht};
  		Which[Length[ni] == Length[ang] == 0, rp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[rp],
   			  True, rp]
        ]

FresnelTransmissionS::usage = "It returns the Fresnel amplitude coefficient of transmission for S polarization. The results is an array of dimensions Length[ang] x Length[ni]";
FresnelTransmissionS[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], costhi, costht, ts},
  	    If[ang < 0 || Pi/2. < ang,
            Return[Message[FresnelCoefficient::angDomain, ang]]
        ];
  	    If[Length[ni] != Length[nt],
            Return[Message[FresnelCoefficient::refIndeces, ni, nt]]
        ];

        costhi = Cos[thetai]; (*[[ang]]*)
  		costht = Sqrt[m^2 - #^2] & /@ Sin[thetai]; (*[[ang, ni]] *)
  		ts = Divide[2*#1, Plus[#1, #2]] & @@ {costhi, costht};
  		Which[Length[ni] == Length[ang] == 0, ts[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[ts],
   			  True, ts]
        ]

FresnelTransmissionP::usage = "It returns the Fresnel amplitude coefficient of transmission for P polarization. The results is an array of dimensions Length[ang] x Length[ni]";
FresnelTransmissionP[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], Ncosthi, costht, tp},
  	    If[ang < 0 || Pi/2. < ang,
            Return[Message[FresnelCoefficient::angDomain, ang]]
        ];
  	    If[Length[ni] != Length[nt],
            Return[Message[FresnelCoefficient::refIndeces, ni, nt]]
        ];

        Ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		tp = Divide[ 2*#1  , Plus[#1, #2]] & @@ {Ncosthi, costht};
  		tp = Map[#/m & , tp];
  		Which[Length[ni] == Length[ang] == 0, tp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[tp],
   			  True, tp]
        ]

FresnelTransmissionP::usage = "It returns the the superficial impedances between two diferent planar interfaces. The results is an array of dimensions Length[ang] x Length[ni]";
FresnelImpedance[ang_, {ni_, nt_}] :=
    Module[{ thetai = toMap[ang], m = toMap[nt/ni], imp},
  	    If[ang < 0 || Pi/2. < ang,
            Return[Message[FresnelCoefficient::angDomain, ang]]
        ];
  	    If[Length[ni] != Length[nt],
            Return[Message[FresnelCoefficient::refIndeces, ni, nt]]
        ];
  		imp = Re@Outer[(Sqrt[#2^2 - Sin[#1]^2]/Cos[#1]) & ,  thetai, m];
  		Which[Length[ni] == Length[ang] == 0, imp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[imp],
   			  True, imp]
        ]
