#! /usr/bin/env wolframscript
(* ::Package:: *)

(* ---------------------------------------------------------------------------------------------
Jonathan Urrutia
Plane Interface Optics Package  (Under developement)
Standard Fresnel's Transmission and Reflection Amplitude Coefficients are define
as function of either the angle of incidence or the refractive indices of the media
formin the interface.

All functions thread over both the angle of incidence and  the refractive index.

Las Modification: 2022 / 06 / 29
------------------------------------------------------------------------------------------------*)
toMap::usage = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;

norm2 = Chop[#*Conjugate[#]]& ;


(* ---------------------------------------------------------------------------------------------
 Fresnel's Coefficients for Transmission and Refelction considering both polarizations
------------------------------------------------------------------------------------------------*)
FresnelCoefficient::IncAngleRange = "The angle of incidence `1` must be between 0 and \[Pi]/2 radins.";
FresnelCoefficient::RefIndicesUnmatch = "The length of `1` must equal the length of `2`.";
FresnelReflectionS::usage = "It returns the Fresnel's amplitude coefficient of reflection for S polarization. The results is an array of dimensions {Length[ang], Length[ni]}."
FresnelReflectionP::usage = "It returns the Fresnel's amplitude coefficient of reflection for P polarization. The results is an array of dimensions {Length[ang], Length[ni]}.";
FresnelTransmissionS::usage = "It returns the Fresnel's amplitude coefficient of transmission for S polarization. The results is an array of dimensions {Length[ang], Length[ni]}."
FresnelTransmissionP::usage = "It returns the Fresnel's amplitude coefficient of transmission for P polarization. The results is an array of dimensions {Length[ang], Length[ni]}."
FresnelImpedance::usage = "It returns surface impedance if an interface between two medi. The results is an array of dimensions {Length[ang], Length[ni]}."
(* --------------------------------------------------------------------------------------------- *)
FresnelReflectionS[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], costhi, costht, rs},
  	    If[ang < 0. || Pi/2. < ang, Return[Message[FresnelCoefficient::IncAngleRange, ang]]];
  	    If[Length[ni] != Length[nt], Return[Message[FresnelCoefficient::RefIndicesUnmatch, ni, nt]]];

      costhi = Cos[thetai]; (*[[ang]]*)
  	costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m]; (*[[ang, m]]*) (*Cosine of the transmitted angle*)
  	rs = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {costhi, costht}; (*[[ang]]*)
  	
  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[ni] == Length[ang] == 0, rs[[1, 1]], 
  		  Length[ni] == 0 || Length[ang] == 0, Flatten[rs],
   		 True, rs]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
FresnelReflectionP[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], ncosthi, costht, rp},
  	    If[ang < 0. || Pi/2. < ang, Return[Message[FresnelCoefficient::IncAngleRange, ang]]];
  	    If[Length[ni] != Length[nt], Return[Message[FresnelCoefficient::RefIndicesUnmatch, ni, nt]]];

  		ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		rp = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {ncosthi, costht};
  		
  		(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  		Which[Length[ni] == Length[ang] == 0, rp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[rp],
   			  True, rp]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
FresnelTransmissionS[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], costhi, costht, ts},
  	    If[ang < 0. || Pi/2. < ang, Return[Message[FresnelCoefficient::IncAngleRange, ang]]];
  	    If[Length[ni] != Length[nt], Return[Message[FresnelCoefficient::RefIndicesUnmatch, ni, nt]]];
  	    
  	    costhi = Cos[thetai]; (*[[ang]]*)
  		costht = Sqrt[m^2 - #^2] & /@ Sin[thetai]; (*[[ang, ni]] *)
  		ts = Divide[2*#1, Plus[#1, #2]] & @@ {costhi, costht}; (*[[ang, ni]] *)
  		Which[Length[ni] == Length[ang] == 0, ts[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[ts],
   			  True, ts]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
FresnelTransmissionP[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], ncosthi, costht, tp},
  	    If[ang < 0. || Pi/2. < ang, Return[Message[FresnelCoefficient::IncAngleRange, ang]]];
  	    If[Length[ni] != Length[nt], Return[Message[FresnelCoefficient::RefIndicesUnmatch, ni, nt]]];

          ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		tp = Divide[ 2*#1  , Plus[#1, #2]] & @@ {ncosthi, costht};
  		tp = Map[#/m & , tp];
  		Which[Length[ni] == Length[ang] == 0, tp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[tp],
   			  True, tp]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
FresnelImpedance[ang_, {ni_, nt_}] :=
    Module[{ thetai = toMap[ang], m = toMap[nt/ni], mcosthi, costht,imp},
  	    If[ang < 0. || Pi/2. < ang, Return[Message[FresnelCoefficient::IncAngleRange, ang]]];
  	    If[Length[ni] != Length[nt], Return[Message[FresnelCoefficient::RefIndicesUnmatch, ni, nt]]];
  	    
  	    imp = Re@Outer[(Sqrt[#2^2-Sin[#1]^2]/Cos[#1])& ,  thetai,m];
  		Which[Length[ni] == Length[ang] == 0, imp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[imp],
   			  True, imp]]
(* --------------------------------------------------------------------------------------------- *) 			  


(* ---------------------------------------------------------------------------------------------
Landau Three Media formula for reflection amplitude coefficient: Both polarization
------------------------------------------------------------------------------------------------*)

Landau3Media::thick = "The thickness `1` must be a scalar.";
Landau3Media::lda = "The Length of `1`, `2`, `3` and `4` must be the same.";
Landau3MediaReflectionS::usage = "It returns the reflection amplitud coeficient of a three media system, where an anisotropy refractive index can be considered.";
Landau3MediaReflectionP::usage = "It returns the reflection amplitud coeficient of a three media system, where an anisotropy refractive index can be considered.";
Landau3MediaElipsometry::usage = "It returns the elipsometric parameters {Psi,Delta} in radians given by rp/rs = Tan[Psi]*Exp[I*Delta]. An anisotropic refractive index can be taken into account.";
(* --------------------------------------------------------------------------------------------- *)


(* --------------------------------------------------------------------------------------------- *)
Landau3MediaReflectionS[ang_,{n1_,n2_,n3_},lda_,thick_]:= 
	Module[{thetai = toMap[ang], epsilon = {n1,n2,n3}^2 ,kz, dkz, r12,r23,r123,psi},
	If[Length[thick]!= 0, Return[Message[Landau3Media::thick, thick]]];
	If[ Length[lda]== 0,
		If[!Apply[Equal, Length /@ Prepend[epsilon[[{1,3}]],lda] ] && Length[epsilon[[2,1]]] == Length[lda], Return[Message[Landau3Media::lda, n1,n2,n3,lda]],
																										epsilon = Map[toMap,epsilon ,2]; epsilon[[2]] ={Flatten[epsilon[[2]]]}],
		If[!Apply[Equal, Length /@ Prepend[epsilon,lda]],Return[Message[Landau3Media::lda, n1,n2,n3,lda]]]];

	If[Length[Dimensions[epsilon[[2]]]] == 2,  (*If anisotropic Film*)
			epsilon[[2]] = epsilon[[2,;;,1]]]; (*We only consider the parallel component*)

	kz = Map[# * Sin[thetai]^2&, epsilon[[1]]]; (* (n1*Sin[angle])^2 , [[n1, angle]]*)
	kz = Sqrt[# -  kz ]&/@   epsilon  (*[[3(medium), ni, angle]]*);

	dkz = Partition[kz, 2, 1]; 
	r12 = Apply[Divide  , {Plus[#1,-#2],Plus[#1,#2]}]& @@ dkz[[1]];
	r23 =Apply[ Divide  , {Plus[#1,-#2],Plus[#1,#2]}]&@@ dkz[[2]];
	psi = thick* 2*Pi * (kz[[2]] / toMap[lda]);
	r123 = r12+r23*Exp[2*I*psi];
	r123/= 1+r12*r23*Exp[2*I*psi];

	Which[Length[n1]==Length[ang]==0, r123[[1,1]],
			Length[n1]==0 || Length[ang]== 0, Flatten[r123],
			True, Transpose[r123]]]
(* --------------------------------------------------------------------------------------------- *)			



(* --------------------------------------------------------------------------------------------- *)
Landau3MediaReflectionP[ang_,{n1_,n2_,n3_},lda_,thick_]:= 
	Module[{thetai = toMap[ang], epsilon = {n1,n2,n3}^2 ,kz, dkz, r12,r23,r123,psi, ani},
	If[Length[thick]!= 0, Return[Message[Landau3Media::thick, thick]]];
	If[ Length[lda]== 0,
		If[!Apply[Equal, Length /@ Prepend[epsilon[[{1,3}]],lda] ] && Length[epsilon[[2,1]]] == Length[lda], Return[Message[Landau3Media::lda, n1,n2,n3,lda]],
																										epsilon = Map[toMap,epsilon ,2]; epsilon[[2]] ={Flatten[epsilon[[2]]]}],
		If[!Apply[Equal, Length /@ Prepend[epsilon,lda]],Return[Message[Landau3Media::lda, n1,n2,n3,lda]]]];

	ani = ConstantArray[1.,{3,Length[epsilon[[1]]]}];
	If[Length[Dimensions[epsilon[[2]]]] == 2,  (*If anisotropic Film*)
		ani = Map[ If[Length[#]==2, Divide@@ # ,1.]& , epsilon,{2}];
		epsilon[[2]] = epsilon[[2,;;,1]]]; (*We only consider the parallel component*)
	
	kz = Map[# * Sin[thetai]^2&,epsilon[[1]]]; (* (eps1*Sin[angle])^2 , [[n1, angle]]*)
	kz = kz*# &/@ ani; (*ani*Sin[angle]*eps1  [[3,n1,angle]]*)
	kz = Sqrt[epsilon-  kz]  (*[[3(medium), ni, angle]]*);

	r12 = Apply[Divide  , {Plus[#1,-#2],Plus[#1,#2]}]& @@ (epsilon[[{2,1}]]*kz[[{1,2}]]);
	r23 =Apply[ Divide  , {Plus[#1,-#2],Plus[#1,#2]}]&@@ (epsilon[[{3,2}]]*kz[[{2,3}]]);
	psi = thick* 2*Pi * (kz[[2]] / toMap[lda]);
	r123 = r12+r23*Exp[2*I*psi];
	r123/= 1+r12*r23*Exp[2*I*psi];

	Which[Length[n1]==Length[ang]==0, r123[[1,1]],
			Length[n1]==0 || Length[ang]== 0, Flatten[r123],
			True, Transpose[r123]]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
Landau3MediaElipsometry[ang_,{n1_,n2_,n3_},lda_,thick_]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ Through[{Landau3MediaReflectionP,Landau3MediaReflectionS}[ang,{n1,n2,n3},lda,thick]] (*[[pol (2), ang, lda]]*))
(* --------------------------------------------------------------------------------------------- *) 



(* ---------------------------------------------------------------------------------------------
Transfere Matrix for the reflection and transmission amplitude coefficient: Both polarization
------------------------------------------------------------------------------------------------*)

TransferMatrix::thick = "This distances `1` must be a scalar quantity or a list of scalar qunatities.";
TransferMatrix::ind = "This length of `1` must be a the length of `2`  minus 2.";
TransferMatrix::len = "This length of `1` must equal the length of each element of `2`";
TransferMatrix::ani = "This length of `1` must equal the length of each element of `2`";

TransferMatrixS::usage ="It returns the reflection and transmission {r,t} amplitude coefficients for a stratified media; it can consider anisotropic refractive indices.";
TransferMatrixP::usage ="It returns the reflection and transmission {r,t} amplitude coefficients for a stratified media; it can consider anisotropic refractive indices.";
TransferMatrixElipsometry::usage ="IIt returns the elipsometric parameters {Psi,Delta} in radians given by rp/rs = Tan[Psi]*Exp[I*Delta] for a stratified media; it can consider anisotropic refractive indices.";
(* --------------------------------------------------------------------------------------------- *)


(* --------------------------------------------------------------------------------------------- *)

TransferMatrixS[ang_,indices_, lda_, distances_]:=
	Module[{thetai = toMap[ang], indicesPar, dist,  phase,p,Z0,mMatrix,rt={0,0} ,m},
	Z0 = 120.*Pi;(*Sqrt[mu0/eps0]=Impredance free Space (approximated value)*)
	dist = Flatten[{0.,distances,0.},1]; (*[[dist]]*)
	indicesPar = toMap /@ indices ;  (*[[dist, lda]]*)

	If[ !Apply[Equal, Length/@  dist], Return[Message[TransferMatrix::thick,distances]]];
	If[Length[indices]!= Length[dist], Return[Message[TransferMatrix::ind,distances,indices]]];
	If[Length /@ Select[indices, Length[#]== lda&]  != Length[dist], Return[0]];
	
	If[Length[lda]==0, (*We take only the parallel contribution of the refractive indices*)
			If[Flatten [Dimensions /@ indices] != {} , (*If not isotropic*)
						indicesPar = toMap /@  (toMap /@ indices)[[;;,1]];
						If[!Apply[Equal,  Dimensions /@  Append[Select[indices, Length[#]!= 0 &] ,{0,0}] ],
								 Return[Message[TransferMatrix::len,lda,indices]]]],
								If[Flatten[Length/@Dimensions/@ Prepend[indices,#]]!=ConstantArray[1,Length[#]]&[lda],(*If not isotropic *) 
									indicesPar = Map[If[Length[Dimensions[#]]==2,#[[;;,1]],#]&, indices]]];
									
	If[!Apply[Equal, Length/@ Prepend[indicesPar,toMap[lda]]], Return[Message[TransferMatrix::len,lda,indices]]];

	p = (#*Sin[thetai])^2 &/@indicesPar[[1]]; (*(n1Sin[th_i])^2   [[lda, ang]]*)
	p =  (1/Z0)*(Sqrt[#^2 - p]&/@ indicesPar) ; (*[[dist, lda, ang]] *)
	phase =((2.*Pi/toMap[lda])*#&/@  p* Z0 )* dist ; (* kz*d   [[dist, lda, ang]]*)

	m = {{Cos[#1], - I* Sin[#1]/#2},{-I*#2*Sin[#1], Cos[#1]}}&[phase, p];  (*[[ren, col, dist, lda, ang]]*)
	m = Transpose[m,{4,5,3,1,2}];  (*[[lda, ang, dist, ren, col]]*)
	m = Array[Apply[Dot,m[[#1,#2]]]&,Length/@{toMap[lda], thetai} ]; (*[[lda, ang, ren, col]]*)
	m = Transpose[m, {3,4,1,2}]; (*[[ren, col, lda, ang]]*)

	rt[[1]] =  (m[[1,1]]+m[[1,2]]*p[[-1]])*p[[1]]-(m[[2,1]]+m[[2,2]]*p[[-1]]);
	rt[[1]] /=( (m[[1,1]]+m[[1,2]]*p[[-1]])*p[[1]]+(m[[2,1]]+m[[2,2]]*p[[-1]]));
	rt[[2]] = 2 *p[[1]]/( (m[[1,1]]+m[[1,2]]*p[[-1]])*p[[1]]+(m[[2,1]]+m[[2,2]]*p[[-1]])); (*rt [[rt(2), lda, ang]]*)
	rt = Transpose[rt,{1,3,2}]; (*rt [[rt(2), ang, lda]]*)

	Which[Length[lda]==Length[ang]==0, rt[[;;, 1,1]],
		Length[lda]==0 , rt[[;;,;;,1]],
		Length[ang]==0, rt[[;;,1,;;]],
		True, rt]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
TransferMatrixImpedanceS[ang_,indices_, lda_, distances_]:= 
	Module[{thetai = toMap[ang], indicesPar, dist,  p},
	dist = Flatten[{0.,distances,0.},1]; (*[[dist]]*)
	indicesPar = toMap /@ indices ;  (*[[dist, lda]]*)

	If[ !Apply[Equal, Length/@  dist], Return[Message[TransferMatrix::thick,distances]]];
	If[Length[indices]!= Length[dist], Return[Message[TransferMatrix::ind,distances,indices]]];
	If[Length /@ Select[indices, Length[#]== lda&]  != Length[dist], Return[0]];
	
	If[Length[lda]==0, (*We take only the parallel contribution of the refractive indices*)
			If[Flatten [Dimensions /@ indices] != {} , (*If not isotropic*)
						indicesPar = toMap /@  (toMap /@ indices)[[;;,1]];
						If[!Apply[Equal,  Dimensions /@  Append[Select[indices, Length[#]!= 0 &] ,{0,0}] ],
								 Return[Message[TransferMatrix::len,lda,indices]]]],
								If[Flatten[Length/@Dimensions/@ Prepend[indices,#]]!=ConstantArray[1,Length[#]]&[lda],(*If not isotropic *) 
									indicesPar = Map[If[Length[Dimensions[#]]==2,#[[;;,1]],#]&, indices]]];
									
	If[!Apply[Equal, Length/@ Prepend[indicesPar,toMap[lda]]], Return[Message[TransferMatrix::len,lda,indices]]];

	p = (#*Sin[thetai])^2 &/@indicesPar[[1]]; (*(n1Sin[th_i])^2   [[lda, ang]]*)
	p =  (Sqrt[#^2 - p]&/@ indicesPar) ; (*[[dist, lda, ang]] *)
	p = Re [Divide @@ p[[{-1,1}]] ];

	Which[Length[lda]==Length[ang]==0, p[[1,1]],
		Length[lda]==0 , p[[1]],
		Length[ang]==0, p[[;;,1]],
		True, Transpose[p]]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
TransferMatrixP[ang_,indices_, lda_, distances_]:= 
	Module[{thetai = toMap[ang], indicesPar, dist,  phase,q,Z0,mMatrix,rt={0,0} ,m, ani},
	Z0 = 120.Pi;(*Sqrt[mu0/eps0]=Impredance free Space*)
	dist = Flatten[{0.,distances,0.},1]; (*[[dist]]*)
	indicesPar = toMap /@ indices ;  (*[[dist, lda]]*)
	ani = ConstantArray[1.,Dimensions[indicesPar]];

	If[ !Apply[Equal, Length/@  dist], Return[Message[TransferMatrix::thick,distances]]];
	If[Length[indices]!= Length[dist], Return[Message[TransferMatrix::ind,distances,indices]]];
	If[Length /@ Select[indices, Length[#]== lda&]  != Length[dist], Return[0]];
	If[Length[lda]==0, (*We take only the parallel contribution of the refractive indices*)
		If[Flatten [Dimensions /@ indices] != {} , (*If not isotropic*) 
				ani = Map[ If[Length[#]==2, Divide@@ # ,1.]& , indices];
					indicesPar = toMap /@  (toMap /@ indices)[[;;,1]];
					If[!Apply[Equal,  Dimensions /@  Append[Select[indices, Length[#]!= 0 &] ,{0,0}] ] , Return[Message[TransferMatrix::len,lda,indices]]]],
				If[Flatten[Length/@Dimensions/@ Prepend[indices,#]]!=ConstantArray[1,Length[#]]&[lda],(*If not isotropic *)  
					ani =Map[ If[Length[#]==2, Divide@@ # ,1.]& , indices, {2}];
					indicesPar = Map[If[Length[Dimensions[#]]==2,#[[;;,1]],#]&, indices]]];
					
	If[!Apply[Equal, Length/@ Prepend[indicesPar,toMap[lda]]], Return[Message[TransferMatrix::len,lda,indices]]];

	q = #*Sin[thetai] &/@indicesPar[[1]]; (*(n1Sin[th_i])   [[lda, ang]]*)
	q = #*q &/@ ani; (*a*n* Sin[thi]   [[dist, lda, ang]]*)
	q =  Z0*Sqrt[indicesPar^2 - q^2] ; (*[[dist, lda, ang]] *)
	phase =((2.*Pi/toMap[lda])*#&/@  (q/ Z0) )* dist ; (* kz*d   [[dist, lda, ang]]*)
	q /= indicesPar^2;

	m = {{Cos[#1], - I* Sin[#1]/#2},{-I*#2*Sin[#1], Cos[#1]}}&[phase, q];  (*[[ren, col, dist, lda, ang]]*)
	m  = Transpose[m,{4,5,3,1,2}];  (*[[lda, ang, dist, ren, col]]*)
	m = Array[Apply[Dot,m[[#1,#2]]]&,Length/@{toMap[lda], thetai} ]; (*[[lda, ang, ren, col]]*)
	m = Transpose[m, {3,4,1,2}]; (*[[ren, col, lda, ang]]*)

	rt[[1]] =  (m[[1,1]]+m[[1,2]]*q[[-1]])*q[[1]]-(m[[2,1]]+m[[2,2]]*q[[-1]]);
	rt[[1]] /=( (m[[1,1]]+m[[1,2]]*q[[-1]])*q[[1]]+(m[[2,1]]+m[[2,2]]*q[[-1]]));
	rt[[2]] = 2 *q[[1]]/( (m[[1,1]]+m[[1,2]]*q[[-1]])*q[[1]]+(m[[2,1]]+m[[2,2]]*q[[-1]])); (*rt [[rt(2), lda, ang]]*)
	rt = Transpose[rt,{1,3,2}]; (*rt [[rt(2), ang, lda]]*)

	Which[Length[lda]==Length[ang]==0, rt[[;;, 1,1]],
		Length[lda]==0 , rt[[;;,;;,1]],
		Length[ang]==0, rt[[;;,1,;;]],
		True, rt]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
TransferMatrixImpedanceP[ang_,indices_, lda_, distances_]:=
	Module[{thetai = toMap[ang], indicesPar, dist,  phase,q, ani},
	dist = Flatten[{0.,distances,0.},1]; (*[[dist]]*)
	indicesPar = toMap /@ indices ;  (*[[dist, lda]]*)
	ani = ConstantArray[1.,Dimensions[indicesPar]];

If[ !Apply[Equal, Length/@  dist], Return[Message[TransferMatrix::thick,distances]]];
	If[Length[indices]!= Length[dist], Return[Message[TransferMatrix::ind,distances,indices]]];
	If[Length /@ Select[indices, Length[#]== lda&]  != Length[dist], Return[0]];
	If[Length[lda]==0, (*We take only the parallel contribution of the refractive indices*)
		If[Flatten [Dimensions /@ indices] != {} , (*If not isotropic*) 
				ani = Map[ If[Length[#]==2, Divide@@ # ,1.]& , indices];
					indicesPar = toMap /@  (toMap /@ indices)[[;;,1]];
					If[!Apply[Equal,  Dimensions /@  Append[Select[indices, Length[#]!= 0 &] ,{0,0}] ] , Return[Message[TransferMatrix::len,lda,indices]]]],
				If[Flatten[Length/@Dimensions/@ Prepend[indices,#]]!=ConstantArray[1,Length[#]]&[lda],(*If not isotropic *)  
					ani =Map[ If[Length[#]==2, Divide@@ # ,1.]& , indices, {2}];
					indicesPar = Map[If[Length[Dimensions[#]]==2,#[[;;,1]],#]&, indices]]];
					
	If[!Apply[Equal, Length/@ Prepend[indicesPar,toMap[lda]]], Return[Message[TransferMatrix::len,lda,indices]]];

	q = (#*Sin[thetai])^2 &/@indicesPar[[1]]; (*(n1Sin[th_i])^2   [[lda, ang]]*)
	q = #*q &/@ ani; (*a*(n1 Sin[thi])^2   [[dist, lda, ang]]*)
	q =  Sqrt[indicesPar^2 - q] /indicesPar^2; (*[[dist, lda, ang]] *)
	q = Re[Divide @@ q[[{-1,1}]]];

	Which[Length[lda]==Length[ang]==0, q[[ 1,1]],
		Length[lda]==0 , q[[1]],
		Length[ang]==0, q[[;;,1]],
		True, Transpose[q]]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
TransferMatrixElipsometry[ang_,indices_,lda_,distances_]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ (Through[{TransferMatrixP,TransferMatrixS}[ang,indices,lda,distances]])[[;;,1]] (*[[pol (2), ang, lda]]*)) 
(* --------------------------------------------------------------------------------------------- *)
