#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ---------------------------------------------------------------------------------------------
Jonathan Urrutia
Coherent Scattering Model Package  (Under developement)
This package requires the CSM
Standard Mie Theory related functions as described by Bohren & Hufmann.

Las Modification: 2022 / 06 / 30
------------------------------------------------------------------------------------------------*)
toMap::usage = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;

norm2::usage = "The squared norm of quantity. This functions threads over all indices.";
norm2 = Chop[#*Conjugate[#]]& ;


(* --------------------------------------------------------------------------------------------- *)
MieAB[n_, x_, m_] :=
    Module[{an, bn, psiMX, dpsiMX, psiX, dpsiX, xiX, dxiX},
        {psiMX, dpsiMX} = {m*x*#, -n*# + m*x*SphericalBesselJ[n - 1, m*x]} &@  SphericalBesselJ[n, m*x];
        {psiX, dpsiX} = {x*#, -n*# + x*SphericalBesselJ[n - 1, x]} &@ SphericalBesselJ[n, x];
        {xiX, dxiX} = {psiX, dpsiX} + I*{x*#, -n*# + x*SphericalBesselY[n - 1, x]} &@ SphericalBesselY[n, x];

        an = (m*psiMX*dpsiX - psiX*dpsiMX)/(m*psiMX*dxiX - xiX*dpsiMX);
        bn = (psiMX*dpsiX - m*psiX*dpsiMX)/(psiMX*dxiX - m*xiX*dpsiMX);
  {an, bn}]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MiePI[n_, costh_] :=
    Module[{f},
        f[0] = 0;
        f[1] = 1;
        f[i_] := f[i] = ((2 i - 1)/(i - 1)) * costh * f[i - 1] - (i/(i - 1)) * f[i - 2];
    f[n]]
SetAttributes[MiePI, Listable]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieTAU[n_, costh_] :=
    Module[{f}, f[0] = 0; f[1] = 1;
        f[i_] :=
        f[i] = ((2 i - 1)/(i - 1)) * costh *f[i - 1] - (i/(i - 1)) * f[i - 2];
    n* costh *f[n] - (n + 1)*f[n - 1]]
SetAttributes[MieTAU, Listable]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieS012[x_, m_, angle_] :=
    Module[{ab, poles, pitau, coeff,  s0, s1, s2},
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]]; (*Wacombe criteria for convergence*)

        ab = MieAB[poles, x, m];
        pitau = Through[{MiePI, MieTAU}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;
		
		s0 = 0.5 * Plus @@ (((2.*#+1) &/@ poles) * Plus @@ ab);
        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  { s0, s1, s2}]
(* --------------------------------------------------------------------------------------------- *)


(* --------------------------------------------------------------------------------------------- *)
fresnelRs[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], costhi, costht, rs},
      costhi = Cos[thetai]; (*[[ang]]*)
  	costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m]; (*[[ang, m]]*) (*Cosine of the transmitted angle*)
  	rs = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {costhi, costht}; (*[[ang]]*)
  	
  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[ni] == Length[ang] == 0, rs[[1, 1]], 
  		  Length[ni] == 0 || Length[ang] == 0, Flatten[rs],
   		 True, rs]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
fresnelRp[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], ncosthi, costht, rp},
  	    ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		rp = Divide[Plus[#1, -#2], Plus[#1, + #2]] & @@ {ncosthi, costht};
  		
  		(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  		Which[Length[ni] == Length[ang] == 0, rp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[rp],
   			  True, rp]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
fresnelTs[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni]op, thetai = toMap[ang], costhi, costht, ts},
  	    costhi = Cos[thetai]; (*[[ang]]*)
  		costht = Sqrt[m^2 - #^2] & /@ Sin[thetai]; (*[[ang, ni]] *)
  		ts = Divide[2*#1, Plus[#1, #2]] & @@ {costhi, costht}; (*[[ang, ni]] *)
  		Which[Length[ni] == Length[ang] == 0, ts[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[ts],
   			  True, ts]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
fresnelTp[ang_, {ni_, nt_}] :=
    Module[{ m = toMap[nt/ni], thetai = toMap[ang], ncosthi, costht, tp},
  	    ncosthi = Outer[#2^2 * Cos[#1] & ,  thetai, m];
  		costht = Outer[Sqrt[#2^2 - Sin[#1]^2] & ,  thetai, m];
  		tp = Divide[ 2*#1  , Plus[#1, #2]] & @@ {ncosthi, costht};
  		tp = Map[#/m & , tp];
  		Which[Length[ni] == Length[ang] == 0, tp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[tp],
   			  True, tp]]
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
CSMImpedance[ang_, {ni_, nt_}] :=
    Module[{ thetai = toMap[ang], m = toMap[nt/ni], mcosthi, costht,imp},
  	    imp = Re@Outer[(Sqrt[#2^2-Sin[#1]^2]/Cos[#1])& ,  thetai,m];
  		Which[Length[ni] == Length[ang] == 0, imp[[1, 1]],
   			  Length[ni] == 0 || Length[ang] == 0, Flatten[imp],
   			  True, imp]]
(* --------------------------------------------------------------------------------------------- *) 			  


(* ---------------------------------------------------------------------------------------------
Mie Coefficients for the scattered and the internal electric field due to a spherical particle
	x -> Size Parameter
    m -> refractive indices contrast
    Dimensions[x]  == Dimensions[m] == 0
    n -> Multipole orders // threads over this variable
------------------------------------------------------------------------------------------------*)
CSM::radius = "The radius `1` must be a scalar.";
CSM::cfrac = "The coverage fraction `1` must be a scalar.";
CSM::length = "The Length of `1`, `2` and `3` must be the same.";

CSMAmplitudeReflectionTransmission::usage = "CSMAmplitudeReflectionTransmission returns the amplitude reflection and transmission coefficidnts given by the CSM. It returns an array of dimensions 2x2xLength[ang]xLenght[lda], pol x ref/tra  is the 2x2.";
(* --------------------------------------------------------------------------------------------- *)
(*CSMAmplitudeReflectionTransmission[ang, {nM,nP},lda, rad,cfrac] = {{rs,ts},{rp,tp}}*)
CSMAmplitudeReflectionTransmission[ang_, {nM_,nP_},lda_, rad_,cfrac_ ]:=
	Module[{x, alpha, m, thetai = toMap[N[ang]],S012, func,aa},
		If[Length[rad]!=0,Return[Message[CSM::radius,rad]]];
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
	
		x = toMap[2.*Pi*nM*rad/ lda];  (*[[lda]]*)
		m = toMap[nP/nM ];(* [[lda]] *)
		alpha  = 2*cfrac / (x^2 * #) &/@ Cos[thetai];  (*[[theta, lda]] *)
		S012 = Array[MieS012 @@ {x[[#1]],m[[#1]],Pi-2*thetai[[#2]]}&, Length /@ {x,thetai} ]; (*[[lda,ang, 3]]*)
		S012 = Transpose[S012,{3,2,1}]; (*[[3, ang, lda]]*)

		func = {- alpha* #2 /(1+alpha* #1 + .25 * alpha^2 * (#1^2 - #2^2)), (1.-0.25*alpha^2 *(#1^ 2-#2^2))/(1+alpha* #1 + .25 * alpha^2 * (#1^2 - #2^2))}&;
		aa = Map[ func @@ # &, S012[[{1,#}]]&/@ {2,3}];  (*[[pol (2) , rt (2), ang, lda]]*)
	Which[Length[nM]==Length[ang]==0, aa[[;;,;;,1,1]],
		Length[nM]==0 , aa[[;;,;;,;;,1]],
		Length[ang]== 0, aa[[;;,;;,1,;;]],
		True, aa]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
CSM3MediaReflectionExt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:=
	Module[{thetai = toMap[ang], wlength = toMap[lda], rtCoh, r13, beta, func,rr},
	If[Length[rad]!=0,Return[Message[CSM::radius,rad]]];
	If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
	If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
	
	rtCoh= CSMAmplitudeReflectionTransmission[thetai,toMap/@{nM,nP}, wlength, rad, cfrac]; (*[[pol, rt, ang, lda]]*)
	r13= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nM,nS}]]; (*[[pol, ang, lda]]*)
	beta = (2.*Pi/wlength)*(rad)*toMap[nM] * # &/@ Cos[thetai];  (*[[ang, lda]]*)

	func = #1[[1]] +( #2*#1[[2]]*#1[[2]]*#1[[2]]*Exp[2*I*beta])/(1.- #2 * #1[[1]] *Exp[2*I*beta])&;
	rr = MapThread[func ][ {rtCoh, r13}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
CSM3MediaReflectionInt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rCoh, r31, beta, func,rr},
		If[Length[rad]!=0,Return[Message[CSM::radius,rad]]];
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];

		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)

		rCoh = Array[ CSMAmplitudeReflectionTransmission[ thetat[[#1,#2]],Map[toMap,{nM,nP}][[;;, #2]], wlength[[#2]], rad, cfrac]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
		rCoh = Transpose[rCoh, {3,4,2,1}]; (*[[ rt, pol,  ang, lda]]*)
		rCoh = rCoh[[1]] ; (*[[pol, ang, lda]]  -  Only reflection*)
		r31= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nS,nM}]]; (*[[pol, ang, lda]]*)
		beta = (2.*Pi/wlength)*(rad)*toMap[nM]*# & /@ Cos[thetat]; (*[[ang, lda]]*)

		func = (#2 + #1* Exp[2*I*beta])/(1.+#1 * #2 * Exp[2*I*beta]) &;
		rr = MapThread[func][{rCoh, r31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
CSM3MediaTransmissionExt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:=
	Module[{thetai = toMap[ang], wlength = toMap[lda], rtCoh, rt13, beta, func,tt},
		If[Length[rad]!=0,Return[Message[CSM::radius,rad]]];
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];

		rtCoh= CSMAmplitudeReflectionTransmission[thetai,toMap/@{nM,nP}, wlength, rad, cfrac]; (*[[pol, rt, ang, lda]]*)
		rt13= Through[{fresnelRs,fresnelTs, fresnelRp,fresnelTp}[thetai, toMap/@{nM,nS}]]; (*[[4, ang, lda*)
		rt13= Partition[rt13,2];(*[[pol, rt, ang, lda*)
		beta = (2.*Pi/wlength)*(rad)*toMap[nM] * # &/@ Cos[thetai];  (*[[ang, lda]]*)

		func = ( #1[[2]]* #2[[2]] *Exp[I*beta])/(1.- #1[[1]]* #2[[1]] *Exp[2*I*beta])   & ;
		tt = MapThread[func ][ {rtCoh, rt13}] ;(*[[pol (2), ang, lda]]*)
		
	Which[Length[nM]==Length[ang]==0, tt[[;;,1,1]],
		Length[nM]==0 , tt[[;;,;;,1]],
		Length[ang]== 0, tt[[;;,1,;;]],
		True, tt]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
CSM3MediaTransmissionInt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, r13t31, beta, func,tt},
		If[Length[rad]!=0,Return[Message[CSM::radius,rad]]];
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];

		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)
		
		rtCoh = Array[ CSMAmplitudeReflectionTransmission[ thetat[[#1,#2]],Map[toMap,{nM,nP}][[;;, #2]], wlength[[#2]], rad, cfrac]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
		rtCoh = Transpose[rtCoh, {3,4,1,2}]; (*[[pol, rt, ang, lda]]*)
		r13t31 = Table[Array[rfresnel[thetat[[#1,#2]], {toMap[nM][[#2]],toMap[nS][[#2]]}]&,Dimensions[thetat]],{rfresnel,{fresnelRs, fresnelRp}}]; (*[[pol, ang, lda]]*)
		r13t31 = Partition[#,2]&@Join[r13t31,  Through[{fresnelTs, fresnelTp}[thetai, toMap/@{nS,nM}]]]; (*[[rt, pol, ang, lda]]*)
		r13t31 = Transpose[r13t31]; (*[[pol, rt, ang, lda]]*)
		beta = Transpose[(2.*Pi/wlength)*(rad)*toMap[nM]*Transpose[Cos[thetat]]]; (*[[ang, lda]]*)

		func = ( #1[[2]]* #2[[2]] *Exp[I*beta])/(1.- #1[[1]]* #2[[1]] *Exp[2*I*beta]) & ;
		tt = MapThread[func][ {rtCoh, r13t31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, tt[[;;,1,1]],
		Length[nM]==0 , tt[[;;,;;,1]],
		Length[ang]== 0, tt[[;;,1,;;]],
		True, tt]]
(* --------------------------------------------------------------------------------------------- *)		
		


(* --------------------------------------------------------------------------------------------- *)
CSM3MediaElipsometryExt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSM3MediaReflectionExt[ang,{nM,nP,nS},lda,rad,cfrac](*[[pol (2), ang, lda]]*)) 
CSM3MediaElipsometryInt[ang_,{nM_,nP_,nS_},lda_,rad_,cfrac_]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSM3MediaReflectionInt[ang,{nM,nP,nS},lda,rad,cfrac](*[[pol (2), ang, lda]]*)) 
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
(*CSMPolyAmplitudeReflectionTransmission = {{rs,ts},{rp,tp}}*)
Options[CSMPolyAmplitudeReflectionTransmission] = {"IntegrationInterval"-> "LogNormal","PartitionSize"-> 40.};
CSMPolyAmplitudeReflectionTransmission::int = "The values of `1` must be Normal, LogNormal or a list of the form {limInf,limSup}."
CSMPolyAmplitudeReflectionTransmission::step = "The values of `1` must a single numerical value."
CSMPolyAmplitudeReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}, OptionsPattern[]]:=
	Module[{km, kmz, eta, x, alpha, m, theta = toMap[N[ang]],S012,rhoS012, func,aa, lim, h, phase, int,betaF,betaC, betaB,radii},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
		If[!NumberQ[#],Return[Message[CSMPolyAmplitudeReflectionTransmission::step,#]]]&[OptionValue["PartitionSize"]];
		lim = Switch[OptionValue["IntegrationInterval"],
				"Normal", mu + {-8,8}*sigma,
				"LogNormal", mu * Power[sigma, {-1,+1}*Sqrt[18.]],
				_, If[Length[#]==2,#,Return[Message[CSMPolyAmplitudeReflectionTransmission::int,#]]]&[OptionValue["IntegrationInterval"]] ];
		
		h = Abs[Differences[lim] / OptionValue["PartitionSize"]][[1]];
		radii = Table[a,{a,lim[[1]], lim[[-1]], h}];
		
		km = toMap[(2.*Pi /lda)*nM]; (*[[lda]]*)
		kmz = km*Cos[#] &/@ theta; (*[[theta, lda]]*)
		eta = 2.*Pi /(km*km * Cos[#]) &/@  theta; (*[[theta, lda]]  Omar Version*)
		(*eta = 2*Pi*Cos[#]/km &/@  theta ; (*[[theta, lda]] Alma Version*) *)
		m = toMap[nP/nM ];(* [[lda]] *)
		
		rhoS012 = Table[x = toMap[2.Pi*nM*a/lda];
				S012 = Array[MieS012[x[[#1]],m[[#1]],Pi-2*theta[[#2]]]&, Length/@{m,theta}]; (*[[lda,ang, 3]]*)
				sizeDist[a, mu, sigma]*Transpose[S012,{3,2,1}] (*[[3, ang, lda]]*) ,
					{a,radii}]; (*[[a, 3, ang, lda ]]*)
		
		rhoS012 = Transpose[rhoS012]; (*[[3 = {S0,S1,S2}, a, ang, lda ]]*)
		phase = Table[Exp[#*2*I*kmz*a], {a,radii}]&/@{1,-1};(*[[ 2 = {+,-}, a, ang, lda ]]*)
		
		
		int = NIntegrate[Interpolation[Transpose[{radii,#1}]][a],{a,lim[[1]]*sigma,lim[[-1]]/sigma}]&;
		betaF = {#,#}&[eta * Array[int[rhoS012[[1,;;,#1,#2]]]&, Length /@ { theta, m}]]; (*[[pol, ang, lda]]*)
		betaB = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[1,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		betaC = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[2,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		
		func = {-#1 , 1.-0.25*(#2*#2-#1*#3)}/(1.+#2+0.25*(#2*#2-#1*#3)) &;
		aa = MapThread[ func , {betaC, betaF, betaB}, 3];  (*[[pol (2), ang, lda, rt (2)]]*)
		aa = Transpose[aa,{1,3,4,2}]; (*[[ pol (2), rt (2),ang, lda]]*)
	Which[Length[nM]==Length[ang]==0, aa[[;;,;;,1,1]],
		Length[nM]==0 , aa[[;;,;;,;;,1]],
		Length[ang]== 0,  aa[[;;,;;,1,;;]],
		True, aa]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
CSMPolyReflectionExt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, r13, beta, func,rr},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
	
	rtCoh= CSMPolyAmplitudeReflectionTransmission[thetai,toMap/@{nM,nP}, lda, cfrac,{mu, sigma, sizeDist}]; (*[[pol, rt, ang, lda]]*)
	r13= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nM,nS}]]; (*[[pol, ang, lda]]*)
	
	func = #1[[1]] +( #2*#1[[2]]*#1[[2]]*#1[[2]])/(1.- #2 * #1[[1]])   & ;
	rr = MapThread[func][{rtCoh, r13}]; (*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
CSMPolyReflectionInt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rCoh, r31, beta, func,rr},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
		
		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)
(*CSMPolyAmplitudeReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}]*)
		rCoh = Array[ CSMPolyAmplitudeReflectionTransmission[ thetat[[#1,#2]], Map[toMap,{nM,nP}][[;;, #2]], lda[[#2]],cfrac,{mu, sigma, sizeDist}]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
		rCoh = Transpose[rCoh, {3,4,2,1}]; (*[[ rt, pol,  ang, lda]]*)
		rCoh = rCoh[[1]] ; (*[[pol, ang, lda]]  -  Only reflection*)
		r31= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nS,nM}]]; (*[[pol, ang, lda]]*)
		(*beta = (2.*Pi/wlength)*(30.)*toMap[nM]*# & /@ Cos[thetat]; (*[[ang, lda]]*)

		func = (#2 + #1* Exp[2*I*beta])/(1.+#1 * #2 * Exp[2*I*beta]) &;*)
		func = (#2 + #1)/(1.+#1 * #2 ) &;
		rr = MapThread[func][{rCoh, r31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)




(* --------------------------------------------------------------------------------------------- *)
CSMPolyTransmissionExt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, rt13, beta, func,tt},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
	
		rtCoh= CSMPolyAmplitudeReflectionTransmission[thetai,toMap/@{nM,nP}, lda, cfrac,{mu, sigma, sizeDist}]; (*[[pol, rt, ang, lda]]*)
		rt13= Through[{fresnelRs,fresnelTs, fresnelRp,fresnelTp}[thetai, toMap/@{nM,nS}]]; (*[[4, ang, lda*)
		rt13= Partition[rt13,2];(*[[pol, rt, ang, lda*)
		
		func = ( #1[[2]]* #2[[2]])/(1.- #1[[1]]* #2[[1]] )   & ;
		tt = MapThread[func ][ {rtCoh, rt13}] ;(*[[pol (2), ang, lda]]*)
		
	Which[Length[nM]==Length[ang]==0, tt[[;;,1,1]],
		Length[nM]==0 , tt[[;;,;;,1]],
		Length[ang]== 0, tt[[;;,1,;;]],
		True, tt]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
CSMPolyTransmissionInt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, r13t31, beta, func,tt},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nP}], Return[Message[CSM::length,nP,nM,lda]]];
		
		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)
		(*CSMPolyAmplitudeReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}]*)
		rtCoh = Array[ CSMPolyAmplitudeReflectionTransmission[ thetat[[#1,#2]], Map[toMap,{nM,nP}][[;;, #2]], lda[[#2]],cfrac,{mu, sigma, sizeDist}]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
		rtCoh = Transpose[rtCoh, {3,4,1,2}]; (*[[  pol, rt,  ang, lda]]*)
		
		r13t31 = Table[Array[rfresnel[thetat[[#1,#2]], {toMap[nM][[#2]],toMap[nS][[#2]]}]&,Dimensions[thetat]],{rfresnel,{fresnelRs, fresnelRp}}]; (*[[pol, ang, lda]]*)
		r13t31 = Partition[#,2]&@Join[r13t31,  Through[{fresnelTs, fresnelTp}[thetai, toMap/@{nS,nM}]]]; (*[[rt, pol, ang, lda]]*)
		r13t31 = Transpose[r13t31]; (*[[pol, rt, ang, lda]]*)

		func = ( #1[[2]]* #2[[2]]) / (1.- #1[[1]]*#2[[1]])& ;
		tt = MapThread[func ][ {rtCoh, r13t31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, tt[[;;,1,1]],
		Length[nM]==0 , tt[[;;,;;,1]],
		Length[ang]== 0, tt[[;;,1,;;]],
		True, tt]]
(* --------------------------------------------------------------------------------------------- *)




(* --------------------------------------------------------------------------------------------- *)
CSMPolyElipsometryExt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSMPolyReflectionExt[ang,{nM,nP,nS},lda,cfrac,{mu, sigma, sizeDist}](*[[pol (2), ang, lda]]*)) 
CSMPolyElipsometryInt[ang_,{nM_,nP_,nS_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSMPolyReflectionInt[ang,{nM,nP,nS},lda,cfrac,{mu, sigma, sizeDist}](*[[pol (2), ang, lda]]*)) 
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
(*CSMPolyAmplitudeReflectionTransmissionSizeCorrected = {{rs,ts},{rp,tp}}*)
Options[CSMPolyAmplitudeReflectionTransmissionSizeCorrected] = {"IntegrationInterval"-> "LogNormal","PartitionSize"-> 40.};
CSMPolyAmplitudeReflectionTransmission::int = "The values of `1` must be Normal, LogNormal or a list of the form {limInf,limSup}."
CSMPolyAmplitudeReflectionTransmission::step = "The values of `1` must a single numerical value."
CSMPolyAmplitudeReflectionTransmissionSizeCorrected[ang_, {nMFunc_,nPFunc_},lda_, cfrac_, {mu_, sigma_, sizeDist_}, OptionsPattern[]]:=
	Module[{km, kmz, eta, x, alpha, m, nM,nP,theta = toMap[N[ang]],S012,rhoS012, func,aa, lim, h, phase, int,betaF,betaC, betaB,radii},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!NumberQ[#],Return[Message[CSMPolyAmplitudeReflectionTransmission::step,#]]]&[OptionValue["PartitionSize"]];
		lim = Switch[OptionValue["IntegrationInterval"],
				"Normal", mu + {-8,8}*sigma,
				"LogNormal", mu * Power[sigma, {-1,+1}*Sqrt[18.]],
				_, If[Length[#]==2,#,Return[Message[CSMPolyAmplitudeReflectionTransmission::int,#]]]&[OptionValue["IntegrationInterval"]] ];
		
		h = Abs[Differences[lim] / OptionValue["PartitionSize"]][[1]];
		radii = Table[a,{a,lim[[1]], lim[[-1]], h}];
		
		nM = nMFunc[lda];
		km = toMap[(2.*Pi /lda)*nM]; (*[[lda]]*)
		kmz = km*Cos[#] &/@ theta; (*[[theta, lda]]*)
		eta = 2.*Pi /(km*km * Cos[#]) &/@  theta; (*[[theta, lda]]  Omar Version*)
		(*eta = 2*Pi*Cos[#]/km &/@  theta ; (*[[theta, lda]] Alma Version*) *)
		
	
		rhoS012 = Table[x = toMap[2.Pi*nM*a/lda];
				nP = nPFunc[lda,a];
				m = toMap[nP/nM ];(* [[lda]] *) 
				S012 = Array[MieS012[x[[#1]],m[[#1]],Pi-2*theta[[#2]]]&, Length/@{m,theta}]; (*[[lda,ang, 3]]*)
				sizeDist[a, mu, sigma]*Transpose[S012,{3,2,1}] (*[[3, ang, lda]]*) ,
					{a,radii}]; (*[[a, 3, ang, lda ]]*)
		
		rhoS012 = Transpose[rhoS012]; (*[[3 = {S0,S1,S2}, a, ang, lda ]]*)
		phase = Table[Exp[#*2*I*kmz*a], {a,radii}]&/@{1,-1};(*[[ 2 = {+,-}, a, ang, lda ]]*)
		
		
		int = NIntegrate[Interpolation[Transpose[{radii,#1}]][a],{a,lim[[1]]*sigma,lim[[-1]]/sigma}]&;
		betaF = {#,#}&[eta * Array[int[rhoS012[[1,;;,#1,#2]]]&, Length /@ { theta, m}]]; (*[[pol, ang, lda]]*)
		betaB = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[1,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		betaC = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[2,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		
		func = {-#1 , 1.-0.25*(#2*#2-#1*#3)}/(1.+#2+0.25*(#2*#2-#1*#3)) &;
		aa = MapThread[ func , {betaC, betaF, betaB}, 3];  (*[[pol (2), ang, lda, rt (2)]]*)
		aa = Transpose[aa,{1,3,4,2}]; (*[[ pol (2), rt (2),ang, lda]]*)
	Which[Length[nM]==Length[ang]==0, aa[[;;,;;,1,1]],
		Length[nM]==0 , aa[[;;,;;,;;,1]],
		Length[ang]== 0,  aa[[;;,;;,1,;;]],
		True, aa]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
CSMPolyReflectionExtSizeCorrected[ang_,{nMFunc_,nPFunc_,nSFunc_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, nM,nS,wlength = toMap[lda], rtCoh, r13, beta, func,rr},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
	
	nM = nMFunc[lda];
	nS = nSFunc[lda];
	rtCoh= CSMPolyAmplitudeReflectionTransmissionSizeCorrected[thetai,{nMFunc,nPFunc}, lda, cfrac,{mu, sigma, sizeDist}]; (*[[pol, rt, ang, lda]]*)
	r13= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nM,nS}]]; (*[[pol, ang, lda]]*)
	
	func = #1[[1]] +( #2*#1[[2]]*#1[[2]]*#1[[2]])/(1.- #2 * #1[[1]])   & ;
	rr = MapThread[func][{rtCoh, r13}]; (*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
CSMPolyReflectionIntSizeCorrected[ang_,{nMFunc_,nPFunc_,nSFunc_},lda_,cfrac_,{mu_, sigma_, sizeDist_}]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rCoh, r31, beta, func,rr, nM, nS},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		
		nM = nMFunc[lda];
		nS = nSFunc[lda];
		
		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)
(*CSMPolyAmplitudeReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}]*)
		rCoh = Array[ CSMPolyAmplitudeReflectionTransmissionSizeCorrected[ thetat[[#1,#2]], {nMFunc,nPFunc}, lda[[#2]],cfrac,{mu, sigma, sizeDist}]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
		rCoh = Transpose[rCoh, {3,4,2,1}]; (*[[ rt, pol,  ang, lda]]*)
		rCoh = rCoh[[1]] ; (*[[pol, ang, lda]]  -  Only reflection*)
		r31= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nS,nM}]]; (*[[pol, ang, lda]]*)
		(*beta = (2.*Pi/wlength)*(30.)*toMap[nM]*# & /@ Cos[thetat]; (*[[ang, lda]]*)

		func = (#2 + #1* Exp[2*I*beta])/(1.+#1 * #2 * Exp[2*I*beta]) &;*)
		func = (#2 + #1)/(1.+#1 * #2 ) &;
		rr = MapThread[func][{rCoh, r31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)


