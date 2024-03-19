#!/usr/bin/env wolframscript
(* ::Package:: *)

toMap = If[Length[#]==0,{#},#]&;


BarzilaiBorwein[df_,{xi_,xf_}]:= With[{diffGrad = (df @@ xf) - (df @@ xi)},
                              Dot[(xf-xi), diffGrad]/ Dot[diffGrad,diffGrad]];
                              
BarzilaiBorwein[df_,var_,{xi_,xf_}]:= With[{diffGrad = (df /. Thread[var -> xf]) - (df /. Thread[var -> xi])},
                              Dot[(xf-xi), diffGrad]/ Dot[diffGrad,diffGrad]]


(************************************
 GradientDescent
  
  Finding a local minimum of a function whose gradient is either given explicitly or
  can be calculated symbolically by Mathematica
 
*************************************)
GradientDescent::usage = "GradientDescent[f[x1,x2,..xn], {{x1...xn},{x01,x02,...,x0n}}, tol] CentralFiniteDifferences[f[x1,x2,..xn], {{x1...xn},{x01,x02,...,x0n}}, tol]  = f[x1,x2,..xn] whose gradient is either given explicitly or
  can be calculated symbolically by Mathematica, considering an initial porint {x01,x02,...,x0n} and with a tolerance tol.
GradientDescent[f, df, {x01,x02,...,x0n}, tol] here f and df are natural functions given by the user."
GradientDescent::Lenghtx0 = "The length of `1` must equal the length of `2`."
GradientDescent::LearningRate = "The LearningRate \.b41\.b4 must be a number or set as 'Barzilai\[Dash]Borwein'. ";
Options[GradientDescent] = {LearningRate -> 1.*^-2, Tolerance -> 1.*^-5, MaxIterations -> 1000};

GradientDescent[f_, {var_,point_}, OptionsPattern[]]:=Module[{x = toMap[var], x0, x1, x2 ,grad,df, alpha,i = 1},
If[! Apply[Equal, Length /@ {var,point}],
            Return[Message[GradientDescent::Lenghtx0, var, point]]];
            
  x0 = toMap[point];
  x1 = x0 + 1.*^-2;
  grad = D[f,{x}];
  alpha = Which[NumberQ[OptionValue[LearningRate]], OptionValue[LearningRate]& ,
                OptionValue[LearningRate] == "Barzilai-Borwein", BarzilaiBorwein,
                True, Return[Message[GradientDescent::LearningRate,OptionValue[LearningRate]]]];
  df = (grad /. Thread[x -> x1]);    
  x2 = x1 - alpha[grad,x,{x0, x1}] * df;
  While[ (Norm[x2-x1]>OptionValue[Tolerance]) && (i < OptionValue[MaxIterations]),
    x0 = x1;
    x1 = x2;
    df = (grad /. Thread[x -> x1]);
    x2 = x1 - alpha[grad,x,{x0,x1}] * df;
    i++];
  Return[{x2, f /. Thread[x -> x2]}]
]

GradientDescent[f_, derivativef_, point_, OptionsPattern[]]:=Module[{x0,x1,x2,df, alpha, i = 1},
 If[! Apply[Equal, Length /@ {derivativef @@ point,point}],
            Return[Message[GradientDescent::Lenghtx0, derivativef , point]]];
  x0 = toMap[point];
  x1 = x0 + 1.*^-2;
  df = toMap[derivativef];
  alpha = Which[NumberQ[OptionValue[LearningRate]], OptionValue[LearningRate]& ,
                OptionValue[LearningRate] == "Barzilai-Borwein", BarzilaiBorwein,
                True, Return[Message[GradientDescent::LearningRate,OptionValue[LearningRate]]]];
   x2 = x1 - alpha[df,{x0, x1}] * (df @@ x1);
  While[(Norm[x2-x1]>OptionValue[Tolerance]) && (i < OptionValue[MaxIterations]),
    x0 = x1;
    x1 = x2;
    x2 = x1 - alpha[df,{x0,x1}] * (df @@ x1);
    i++];
  Return[{x2, f @@ x2}]
]


ForwardGrad[f_,x0_,step_]:=Module[{dummy,h = toMap[step],coeff, xi = toMap[x0],xx},
coeff = {-25./12, 4., -3., 4./3, -1./4};
xx = Array[(dummy = Transpose[xi &/@ Range[Length[coeff]] ]; (*[[xi,5]]*)
            dummy[[#1]] += Range[Length[coeff]]*h[[#1]];
            dummy = Transpose[dummy])&, (*[[5,xi]]*)
              Dimensions[xi]]; (*[[xi,5,xi]]*)
xx = Transpose[xx]; (*[[5,xi,xi]]*)
Plus @@ ((#/h &/@ coeff) * Apply[f,xx,{2}])]; (* [[xi]] <-Plus <-[[xi,5]]<- Transspose <-[[5,xi]] * [[5,xi]]*)

NBarzilaiBorwein[f_,{xi_,xf_}, step_]:= With[{diffGrad = ForwardGrad[f,xf,step] - ForwardGrad[f,xi,step]},
                              Dot[(xf-xi), diffGrad]/ Dot[diffGrad,diffGrad]];
                              
NGradientDescent::step = "The 'Step' parameter `1` must be either a real number or a list of the same length as point `2`.";
Options[NGradientDescent] = {LearningRate -> 1.*^-1, Tolerance -> 1.*^-3, MaxIterations -> 10000, "Step" -> 1.*^-4};
NGradientDescent[f_, point_,OptionsPattern[]]:=Module[{x0,x1,x2,df, alpha, step, i = 1},
  x0 = toMap[point];
  x1 = x0 + 1.*^-2;
  alpha = Which[NumberQ[OptionValue[LearningRate]], OptionValue[LearningRate]& ,
                OptionValue[LearningRate] == "Barzilai-Borwein", NBarzilaiBorwein,
                True, Return[Message[GradientDescent::LearningRate,OptionValue[LearningRate]]]];
  step = Which[Length[OptionValue["Step"]] == 0, Map[ OptionValue["Step"]&,x0],
                Length[OptionValue["Step"]] == Length[x0], OptionValue["Step"],
                True, Return[Message[NGradientDescent::step,OptionValue["Step"],point]]];
  x2 = x1 - alpha[f,{x0,x1},step] * ForwardGrad[f,x1, step];
   While[(Norm[x2-x1]>OptionValue[Tolerance]) && (i < OptionValue[MaxIterations]),
    x0 = x1;
    x1 = x2;
    x2 = x1 - alpha[f,{x0,x1},step] * ForwardGrad[f,x1, step];
    i++];
  Return[{x2, f @@ x2}]
]


ForwardGrad[f_,var_, x0_,step_]:=Module[{dummy,x = toMap[var],h = toMap[step],coeff, xi = toMap[x0],xx},
coeff = {-25./12,4.,-3.,4./3,-1./4};
xx = Array[(dummy = Transpose[xi &/@ Range[Length[coeff]] ]; (*[[xi,5]]*)
            dummy[[#1]] += Range[Length[coeff]]*h[[#1]];
            dummy = Transpose[dummy])&, (*[[5,xi]]*)
              Dimensions[xi]]; (*[[xi,5,xi]]*)
xx = Transpose[xx]; (*[[5,xi,xi]]*)
Plus @@ ((#/h &/@ coeff) * Map[f/.Thread[x->#] &, xx,{2}])]; 

NBarzilaiBorwein[f_,var_, {xi_,xf_}, step_]:= With[{diffGrad = ForwardGrad[f,var, xf,step] - ForwardGrad[f,var, xi,step]},
                              Dot[(xf-xi), diffGrad]/ Dot[diffGrad,diffGrad]];
                              
NGradientDescent[f_, var_, point_,OptionsPattern[]]:=Module[{x = toMap[var],x0,x1,x2,df, alpha, step, i = 1},
  x0 = toMap[point];
  x1 = x0 + 1.*^-2;
  alpha = Which[NumberQ[OptionValue[LearningRate]], OptionValue[LearningRate]& ,
                OptionValue[LearningRate] == "Barzilai-Borwein", NBarzilaiBorwein,
                True, Return[Message[GradientDescent::LearningRate,OptionValue[LearningRate]]]];
  step = Which[Length[OptionValue["Step"]] == 0, Map[ OptionValue["Step"]&,x0],
                Length[OptionValue["Step"]] == Length[x0], OptionValue["Step"],
                True, Return[Message[NGradientDescent::step,OptionValue["Step"],point]]];
  x2 = x1 - alpha[f,var,{x0,x1},step] * ForwardGrad[f,x, x1, step];
  While[(Norm[x2-x1]>OptionValue[Tolerance]) && (i < OptionValue[MaxIterations]),
    x0 = x1;
    x1 = x2;
    x2 = x1 - alpha[f,var,{x0,x1},step] * ForwardGrad[f,x,x1, step];
    i++];
  Return[{x2, f /. Thread[x->x2]}]
]
