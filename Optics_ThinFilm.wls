#! /usr/bin/env wolframscript
(* ::Package:: *)

(* ---------------------------------------------------------------------------------------------
Jonathan Urrutia
Thin Film, Theory  (Under developement)
A scattering theory in the dipolar approximation to describe the optical properties of a monolayer of spheroids places on a substrate
Las Modification: 2024 / 04 / 28
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


ThinFilmGeometricParameters[aRad_,bRad_]:=Module[{vol, exc, lSta = {0,0}, geom},
	exc = 1./Sqrt[(aRad/bRad)^2 - 1.]; (*Flatness parameter instead of excentricity*)
	
	vol = (4./3.) * Pi * (aRad*aRad*bRad);
	geom = (1 + #2) * (#1 * ArcTan[1./#1] - (#2 + 3./2.)*#2*Log[1.+1./#2] + #2 )&[exc,exc*exc];  (*Dynamic geometrical corrections*)
	
	lSta[[1]] = (1. + #*#)*(1. - # * ArcTan[1./#])&[exc]; (*static depolarization factor -- perpendicular*)
	lSta[[2]] = 0.5 * (1 -  lSta[[1]]); (*static depolarization factor -- parallel*)
	Return[{vol, lSta, geom}]]
	
ThinFilmGeometricParameters[aRad_]:=Module[{vol, exc, lSta = {0,0}, geom},	
	vol = (4./3.) * Pi * (aRad*aRad*aRad);
	geom = 1.;  (*For the Sphere; see Albano*)
	lSta = {1.,1.}/3.; (*static depolarization factors *)
	Return[{vol, lSta, geom}]]



ThinFilmReflectionIntS[ang_,{epsilonMatrix_,epsilonParticle_,epsilonSubstrate_},lda_,{aRad_,bRad_},cfrac_] :=
    Module[{thetai = toMap[N[ang]], cos = {0,0}, rs, epsS, epsP, epsM, k0, vol, geom, lDyn, alpha, gamma,size, imageA},

      {epsM,epsP,epsS} = toMap /@ {epsilonMatrix, epsilonParticle, epsilonSubstrate}; (*[[lda]]*)
      k0 = 2*Pi/toMap[lda]; (*[[lda]]*)    
      {vol, lDyn, geom} = ThinFilmGeometricParameters[aRad,bRad];
      imageA = (epsS-epsM) / (epsS+epsM);
      lDyn += ( # * geom * imageA &/@  {1. , 0.5}); (*Dynamic depolarization factors [[perp/para, lda]]*)
      
      alpha = vol * epsM * (epsP - epsM) / (epsM + # * (epsP - epsM) ) &/@ lDyn; (*Polarizability [[perp/para, lda]]*)
      size = k0 * Sqrt[epsM] * aRad; (*[[lda]]*)
     
      If[ Apply[Or, Thread[size < ConstantArray[0.2,Length[k0]]]], 
         alpha[[2]] = 1. / ( 1./alpha[[2]] - #^2 / (4.*Pi*aRad) - I * #^3 / (6.*Pi) )&[size/aRad]];
      gamma = cfrac * alpha / {epsS^2, ConstantArray[1., Length[toMap[lda]]]}; (*surface Polarizability [[perp/para, lda]] -- Cambiamos matriz por sustrato*)
      
      cos[[1]] = Transpose[Cos[thetai]& /@  epsS]; (*[[ang, m]]*) (*Cosine of the incidence angle*)
      cos[[2]] = Outer[Sqrt[1. - #2 * Sin[#1]^2] & ,  thetai, (epsS/epsM)];  (*[[ang, m]]*) (*Cosine of the transmitted angle*)
      
      rs = MapThread[Sqrt[epsS]*#1 - Sqrt[epsM]*#2 + I*k0*gamma[[2]] &, cos ];
      rs /= MapThread[Sqrt[epsS]*#1 + Sqrt[epsM]*#2 - I*k0*gamma[[2]] &, cos];	

  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[epsilonMatrix] == Length[N[ang]] == 0, rs[[1, 1]],
  		  Length[epsilonMatrix] == 0 || Length[N[ang]] == 0, Flatten[rs],
   		 True, rs]]
(* --------------------------------------------------------------------------------------------- *)
ThinFilmReflectionIntS[ang_,{epsilonMatrix_,epsilonParticle_,epsilonSubstrate_},lda_,aRad_,cfrac_] :=
    Module[{thetai = toMap[N[ang]], cos = {0,0}, rs, epsS, epsP, epsM, k0, vol, geom, lDyn, alpha, gamma,size, imageA},

      {epsM,epsP,epsS} = toMap /@ {epsilonMatrix, epsilonParticle, epsilonSubstrate}; (*[[lda]]*)
      k0 = 2*Pi/toMap[lda]; (*[[lda]]*)    
      {vol, lDyn, geom} = ThinFilmGeometricParameters[aRad];
      imageA = (epsS-epsM) / (epsS+epsM);
      lDyn += ( # * geom * imageA &/@  {1. , 0.5}); (*Dynamic depolarization factors [[perp/para, lda]]*)
      
      alpha = vol * epsM * (epsP - epsM) / (epsM + # * (epsP - epsM) ) &/@ lDyn; (*Polarizability [[perp/para, lda]]*)
      size = k0 * Sqrt[epsM] * aRad; (*[[lda]]*)
     
      If[ Apply[Or, Thread[size < ConstantArray[0.2,Length[k0]]]], 
         alpha[[2]] = 1. / ( 1./alpha[[2]] - #^2 / (4.*Pi*aRad) - I * #^3 / (6.*Pi) )&[size/aRad]];
      gamma = cfrac * alpha / {epsS^2, ConstantArray[1., Length[toMap[lda]]]}; (*surface Polarizability [[perp/para, lda]] -- Cambiamos matriz por sustrato*)
      
      cos[[1]] = Transpose[Cos[thetai]& /@  epsS]; (*[[ang, m]]*) (*Cosine of the incidence angle*)
      cos[[2]] = Outer[Sqrt[1. - #2 * Sin[#1]^2] & ,  thetai, (epsS/epsM)];  (*[[ang, m]]*) (*Cosine of the transmitted angle*)
      
      rs = MapThread[Sqrt[epsS]*#1 - Sqrt[epsM]*#2 + I*k0*gamma[[2]] &, cos ];
      rs /= MapThread[Sqrt[epsS]*#1 + Sqrt[epsM]*#2 - I*k0*gamma[[2]] &, cos];	

  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[epsilonMatrix] == Length[N[ang]] == 0, rs[[1, 1]],
  		  Length[epsilonMatrix] == 0 || Length[N[ang]] == 0, Flatten[rs],
   		 True, rs]]
(* --------------------------------------------------------------------------------------------- *)


(* --------------------------------------------------------------------------------------------- *)
ThinFilmReflectionIntP[ang_,{epsilonMatrix_,epsilonParticle_,epsilonSubstrate_},lda_,{aRad_,bRad_},cfrac_] :=
    Module[{thetai = toMap[N[ang]], cos = {0,0}, sin, rp, epsS, epsP, epsM, k0, vol, geom, lDyn, alpha, gamma, kappa = {0,0},nM,nS,size, imageA},

      {epsM,epsP,epsS} = toMap /@ {epsilonMatrix, epsilonParticle, epsilonSubstrate}; (*[[lda]]*)
      {nM,nS} = Sqrt[{epsM,epsS}];
      k0 = 2*Pi/toMap[lda]; (*[[lda]]*)    
      {vol, lDyn, geom} = ThinFilmGeometricParameters[aRad,bRad];
      imageA = (epsS-epsM) / (epsS+epsM);
      lDyn += ( # * geom * imageA &/@  {1. , 0.5}); (*Dynamic depolarization factors [[perp/para, lda]]*)
      
      alpha = vol * (epsP-epsM) * epsM / (epsM + # * (epsP-epsM)) &/@ lDyn; (*Polarizability [[perp/para, lda]]*)
      size = k0 * nM * aRad; (*[[lda]]*)
      If[ Apply[Or, Thread[size < ConstantArray[0.6,Length[k0]]]], 
         alpha[[2]] = 1/(1/alpha[[2]]-#^2/(4.*Pi*aRad)-I*#^3/(6.*Pi))&[size/aRad]];
      gamma = cfrac * alpha / {epsS^2, ConstantArray[1,Length[toMap[lda]]]}; (*surface Polarizability [[perp/para, lda]]*)
      
      cos[[1]] = Transpose[Cos[thetai]& /@  epsS]; (*[[ang, m]]*) (*Cosine of the incidence angle*)
      cos[[2]] = Outer[Sqrt[1 - #2 * Sin[#1]^2] & ,  thetai, (epsS/epsM)];  (*[[ang, m]]*) (*Cosine of the transmitted angle*)
      sin = Transpose[Sin[thetai]& /@  epsS]; (*[[ang, m]]*) (*Sin of the incidence angle*)
            
      kappa[[1]] = (nM*#&/@cos[[1]]) + (nS*#&/@cos[[2]]);
      kappa[[1]] *= (1 - 0.25 * k0*k0* epsS * gamma[[1]] * gamma[[2]] * #^2 &/@ sin);(*[[ang, m ; Es matriz o sustraro??]]*)
      kappa[[2]] = (nM*#&/@cos[[1]]) - (nS*#&/@cos[[2]]);
      kappa[[2]] *= (1 - 0.25 * k0*k0* epsS * gamma[[1]] * gamma[[2]] * #^2 &/@ sin);(*[[ang, m ; Es matriz o sustraro??]]*)
      
      rp = MapThread[#1 - I * k0*gamma[[2]]*#2 + I * k0*nM*nS*epsS*gamma[[1]] * #3 &, {kappa[[2]], cos[[1]]*cos[[2]], sin*sin}];
      rp /= MapThread[#1 - I * k0*gamma[[2]]*#2 - I * k0*nM*nS*epsS*gamma[[1]] * #3 &, {kappa[[1]], cos[[1]]*cos[[2]], sin*sin}];

  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[epsilonMatrix] == Length[N[ang]] == 0, rp[[1, 1]],
  		  Length[epsilonMatrix] == 0 || Length[N[ang]] == 0, Flatten[rp],
   		 True, rp]]
(* --------------------------------------------------------------------------------------------- *)
ThinFilmReflectionIntP[ang_,{epsilonMatrix_,epsilonParticle_,epsilonSubstrate_},lda_,aRad_,cfrac_] :=
    Module[{thetai = toMap[N[ang]], cos = {0,0}, sin, rp, epsS, epsP, epsM, k0, vol, geom, lDyn, alpha, gamma, kappa = {0,0},nM,nS,size, imageA},

      {epsM,epsP,epsS} = toMap /@ {epsilonMatrix, epsilonParticle, epsilonSubstrate}; (*[[lda]]*)
      {nM,nS} = Sqrt[{epsM,epsS}];
      k0 = 2*Pi/toMap[lda]; (*[[lda]]*)    
      {vol, lDyn, geom} = ThinFilmGeometricParameters[aRad];
      imageA = (epsS-epsM) / (epsS+epsM);
      lDyn += ( # * geom * imageA &/@  {1. , 0.5}); (*Dynamic depolarization factors [[perp/para, lda]]*)
      
      alpha = vol * (epsP-epsM) * epsM / (epsM + # * (epsP-epsM)) &/@ lDyn; (*Polarizability [[perp/para, lda]]*)
      size = k0 * nM * aRad; (*[[lda]]*)
      If[ Apply[Or, Thread[size < ConstantArray[0.6,Length[k0]]]], 
         alpha[[2]] = 1/(1/alpha[[2]]-#^2/(4.*Pi*aRad)-I*#^3/(6.*Pi))&[size/aRad]];
      gamma = cfrac * alpha / {epsS^2, ConstantArray[1,Length[toMap[lda]]]}; (*surface Polarizability [[perp/para, lda]]*)
      
      cos[[1]] = Transpose[Cos[thetai]& /@  epsS]; (*[[ang, m]]*) (*Cosine of the incidence angle*)
      cos[[2]] = Outer[Sqrt[1 - #2 * Sin[#1]^2] & ,  thetai, (epsS/epsM)];  (*[[ang, m]]*) (*Cosine of the transmitted angle*)
      sin = Transpose[Sin[thetai]& /@  epsS]; (*[[ang, m]]*) (*Sin of the incidence angle*)
            
      kappa[[1]] = (nM*#&/@cos[[1]]) + (nS*#&/@cos[[2]]);
      kappa[[1]] *= (1 - 0.25 * k0*k0* epsS * gamma[[1]] * gamma[[2]] * #^2 &/@ sin);(*[[ang, m ; Es matriz o sustraro??]]*)
      kappa[[2]] = (nM*#&/@cos[[1]]) - (nS*#&/@cos[[2]]);
      kappa[[2]] *= (1 - 0.25 * k0*k0* epsS * gamma[[1]] * gamma[[2]] * #^2 &/@ sin);(*[[ang, m ; Es matriz o sustraro??]]*)
      
      rp = MapThread[#1 - I * k0*gamma[[2]]*#2 + I * k0*nM*nS*epsS*gamma[[1]] * #3 &, {kappa[[2]], cos[[1]]*cos[[2]], sin*sin}];
      rp /= MapThread[#1 - I * k0*gamma[[2]]*#2 - I * k0*nM*nS*epsS*gamma[[1]] * #3 &, {kappa[[1]], cos[[1]]*cos[[2]], sin*sin}];

  	(*The returned value is either a scalar, a list or an array according to the dimension of the given arguments.*)
  	Which[Length[epsilonMatrix] == Length[N[ang]] == 0, rp[[1, 1]],
  		  Length[epsilonMatrix] == 0 || Length[N[ang]] == 0, Flatten[rp],
   		 True, rp]]
