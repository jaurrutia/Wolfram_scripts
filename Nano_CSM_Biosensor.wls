#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ---------------------------------------------------------------------------------------------
Jonathan Urrutia
Coherent Scattering Model Package  (Under developement)
This package requires the CSM
Standard Mie Theory related functions as described by Bohren & Hufmann.

Las Modification: 2022 / 09 / 12
------------------------------------------------------------------------------------------------*)
(* --------------------------------------------------------------------------------------------- *)
(*CSMBioReflectionTransmission = {{rs,ts},{rp,tp}}*)
Options[CSMBioReflectionTransmission] = {"SizeDistribution"-> "LogNormal"};
CSMBioReflectionTransmission::Dist = "The values of `1 must be either 'LogNormal' or 'Normal'."
CSMBioReflectionTransmission[ang_, nM_,lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:=
	Module[{km, kmz, eta, x, alpha, m, theta = toMap[N[ang]],S012,rhoS012, func,aa, lim, h, phase, int,betaF,betaC, betaB,radii,rhost,sizeDist,nP},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM}], Return[Message[CSM::length,nP,nM,lda]]];
	
		Switch[OptionValue["SizeDistribution"],
				"LogNormal", lim = mu * Power[sigma, {-1,+1}*Sqrt[18.]];
							 h = Abs[Differences[lim] / 40][[1]];
							 rhost =  cfrac * Exp[-2*Log[sigma]^2]/(Pi * mu^2);
							 sizeDist = rhost * PDF[LogNormalDistribution[Log[mu],Log[sigma]],#]&,
				"Normal", lim = mu + 5*sigma*{-1,+1}; lim[[1]] = If[lim[[1]]< 0,0,lim[[1]]];
							 h = Abs[Differences[lim] / 40][[1]];
							 rhost =  cfrac /(Pi * mu^2);
							 sizeDist = rhost * PDF[NormalDistribution[mu,sigma],#]&,
						_, Return[Message[CSMBioReflectionTransmission::Dist,"SizeDistribution"]]];
		radii = Table[a,{a,lim[[1]], lim[[-1]], h}];
		nP = JohnsonChristyAuRefSize[toMap[lda], radii]; (*[[lda,a]]*)

		km = toMap[(2.*Pi /lda)*nM]; (*[[lda]]*)
		kmz = km*Cos[#] &/@ theta; (*[[theta, lda]]*) 
		eta = 2.*Pi /(km*km * Cos[#]) &/@  theta; (*[[theta, lda]]  Omar Version*)
		(*eta = 2*Pi*Cos[#]/km &/@  theta ; (*[[theta, lda]] Alma Version*) *)
		m = nP/toMap[nM];(* [[lda,a]] *)
		
		rhoS012 = Table[x = toMap[2.Pi*nM*radii[[i]]/lda];
				S012 = Array[MieS012[x[[#1]],m[[#1,i]],Pi-2*theta[[#2]]]&, Length/@{m[[;;,i]],theta}]; (*[[lda,ang, 3]]*)
				sizeDist[radii[[i]]]*Transpose[S012,{3,2,1}] (*[[3, ang, lda]]*) ,
					{i,Length[radii]}]; (*[[a, 3, ang, lda ]]*)
		
		rhoS012 = Transpose[rhoS012]; (*[[3 = {S0,S1,S2}, a, ang, lda ]]*)
		phase = Table[Exp[#*2*I*a], {a,radii}]&/@{kmz,-Re[kmz]+Abs[Im[kmz]]*I};(*[[ 2 = {+,-}, a, ang, lda ]]*) (*MODIFICACI\[CapitalOAcute]N EN LA FASE. 25/08/22*)
		(*phase = Transpose[phase]; (*[[2 = {+,-}, a, ang, lda ]]*)*)

		int = Switch[OptionValue["SizeDistribution"],
				"LogNormal", NIntegrate[Interpolation[Transpose[{radii,#1}]][a],{a,lim[[1]]*sigma,lim[[-1]]/sigma}]&,
				"Normal", NIntegrate[Interpolation[Transpose[{radii,#1}]][a],{a,lim[[1]]+sigma/2,lim[[-1]]-sigma/2}]&];
		
		betaF = {#,#}&@(eta * Array[int[rhoS012[[1,;;,#1,#2]]]&, Length /@ { theta, m}]); (*[[pol, ang, lda]]*)
		betaB = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[1,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		betaC = Table[eta * Array[int[rhoS012[[i,;;,#1,#2]]*phase[[2,;;,#1,#2]]]&, Length /@ {theta, m}], {i,2,3}]; (*[[pol, ang, lda]]*) 
		(*Return[{betaF,betaB,betaC}];*)
		func = {-#1 , 1.-0.25*(#2*#2-#1*#3)}/(1.+#2+0.25*(#2*#2-#1*#3)) &;
		aa = MapThread[ func , {betaC, betaF, betaB}, 3];  (*[[pol (2), ang, lda, rt (2)]]*)
		aa = Transpose[aa,{1,3,4,2}]; (*[[ pol (2), rt (2),ang, lda]]*)
	Which[Length[nM]==Length[ang]==0, aa[[;;,;;,1,1]],
		Length[nM]==0 , aa[[;;,;;,;;,1]],
		Length[ang]== 0, Transpose[aa,3<->4]; aa[[;;,;;,;;,1]],
		True, aa]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
Options[CSMBioReflectionExt] = {"SizeDistribution"-> "LogNormal"};
CSMBioReflectionExt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, r13, beta, func,rr},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM}], Return[Message[CSM::length,nM,nM,lda]]];
	
	rtCoh= CSMBioReflectionTransmission[thetai,toMap[nM], lda, cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]]; (*[[pol, rt, ang, lda]]*)
	r13= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nM,nS}]]; (*[[pol, ang, lda]]*)
	
	func = #1[[1]] +( #2*#1[[2]]*#1[[2]]*#1[[2]])/(1.- #2 * #1[[1]])   & ;
	rr = MapThread[func][{rtCoh, r13}]; (*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
Options[CSMBioReflectionInt] = {"SizeDistribution"-> "LogNormal"};
CSMBioReflectionInt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rCoh, r31, beta, func,rr},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nS}], Return[Message[CSM::length,nS,nM,lda]]];
		
		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)

(*CSMBioReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}]*)
		rCoh = Array[ CSMBioReflectionTransmission[ thetat[[#1,#2]], nM[[#2]], lda[[#2]],cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)

		rCoh = Transpose[rCoh, {3,4,2,1}]; (*[[ rt, pol,  ang, lda]]*)
		rCoh = rCoh[[1]] ; (*[[pol, ang, lda]]  -  Only reflection*)
		r31= Through[{fresnelRs, fresnelRp}[thetai, toMap/@{nS,nM}]]; (*[[pol, ang, lda]]*)

		func = (#2 + #1)/(1.+#1 * #2 ) &;
		rr = MapThread[func][{rCoh, r31}] ;(*[[pol (2), ang, lda]]*)

	Which[Length[nM]==Length[ang]==0, rr[[;;,1,1]],
		Length[nM]==0 , rr[[;;,;;,1]],
		Length[ang]== 0, rr[[;;,1,;;]],
		True, rr]]
(* --------------------------------------------------------------------------------------------- *)




(* --------------------------------------------------------------------------------------------- *)
Options[CSMBioTransmissionExt] = {"SizeDistribution"-> "LogNormal"};
CSMBioTransmissionExt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, rt13, beta, func,tt},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM}], Return[Message[CSM::length,nM,nM,lda]]];
	
		rtCoh= CSMBioReflectionTransmission[thetai,toMap[nM], lda, cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]]; (*[[pol, rt, ang, lda]]*)
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
Options[CSMBioTransmissionInt] = {"SizeDistribution"-> "LogNormal"};
CSMBioTransmissionInt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:=
	Module[{thetai = toMap[ang],thetat, wlength = toMap[lda], rtCoh, r13t31, beta, func,tt},
		If[Length[cfrac]!=0,Return[Message[CSM::cfrac,cfrac]]];
		If[!Apply[Equal, Length /@ {lda, nM,nS}], Return[Message[CSM::length,nM,nS,lda]]];
		
		thetat = ArcSin[ #* toMap[ nS/nM] ] &/@ Sin[thetai] ; (* [[ang, lda]] *)
		(*CSMBioReflectionTransmission[ang_, {nM_,nP_},lda_, cfrac_, {mu_, sigma_, sizeDist_}]*)
		rtCoh = Array[ CSMBioReflectionTransmission[ thetat[[#1,#2]], nM[[#2]], lda[[#2]],cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]]& , Dimensions[thetat]] ; (*[[ang, lda, pol, rt]]*)
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
Options[CSMBioElipsometryExt] = {"SizeDistribution"-> "LogNormal"};
Options[CSMBioElipsometryInt] = {"SizeDistribution"-> "LogNormal"};

CSMBioElipsometryExt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:= {ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSMBioReflectionExt[ang,{nM,nS},lda,cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]](*[[pol (2), ang, lda]]*))

CSMBioElipsometryInt[ang_,{nM_, nS_},lda_, cfrac_, {mu_, sigma_}, OptionsPattern[]]:={ArcTan[Abs[#]],Arg[#]}&@ (Divide@@ CSMBioReflectionInt[ang,{nM,nS},lda,cfrac,{mu, sigma},"SizeDistribution" -> OptionValue["SizeDistribution"]](*[[pol (2), ang, lda]]*))
(* --------------------------------------------------------------------------------------------- *)

