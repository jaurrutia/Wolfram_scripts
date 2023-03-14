#! /usr/bin/env wolframscript
(* ::Package:: *)

(* ---------------------------------------------------------------------------------------------
Jonathan Urrutia
Mie Theory Package  (Under developement)
Standard Mie Theory related functions as described by Bohren & Hufmann.

Las Modification: 2022 / 06 / 30
------------------------------------------------------------------------------------------------*)
toMap::usage = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;

norm2::usage = "The squared norm of quantity. This functions threads over all indices."
norm2 = Chop[#*Conjugate[#]]& ;

CartesianToSpherical::usage = "Transformation Matrix if the cartesian canonical vector base to the canonical spherical vector base (physisist convention)."
CartesianToSpherical = {{Sin[#1]*Cos[#2],Sin[#1]*Sin[#2], Cos[#1]},
                        {Cos[#1]*Cos[#2],Cos[#1]*Sin[#2], -Sin[#1]},
                        {-Sin[#2],Cos[#2], 0}}& ;

SphericalToCartesian::usage = "Transformation Matrix of the spherical canonical vector base to the canonical cartesina vector base (physisist convention)."
SphericalToCartesian = Transpose[CartesianToSpherical[#1,#2]]& ;



(* ---------------------------------------------------------------------------------------------
Mie Coefficients for the scattered and the internal electric field due to a spherical particle
	x -> Size Parameter
    m -> refractive indices contrast
    Dimensions[x]  == Dimensions[m] == 0
    n -> Multipole orders // threads over this variable
------------------------------------------------------------------------------------------------*)
MieCoefficient::list = "Both `1`and `2` must be scalar quantities.";
MieCoefficient::Int = "Parameter `1` must be an integer number.";
MieCoefficient::usage = "MieCoefficient[n,x,m] returns the Mie Coefficients for a given size parameter x and a relative refractive index. This functions threads over n.";
MieCoefficientInt::usage = "MieCoefficientInt[n,x,m] returns the oefficients for a given size parameter x and a relative refractive index for the internal electric field. This functions threads over n.";
(* --------------------------------------------------------------------------------------------- *)
             
(* --------------------------------------------------------------------------------------------- *)
MieCoefficient[n_, x_, m_] :=
    Module[{an, bn, psiMX, dpsiMX, psiX, dpsiX, xiX, dxiX},
        If[! Apply[Equal, Length /@ {x, m, 0}],
            Return[Message[MieCoefficient::list, x, m]]];

        {psiMX, dpsiMX} = {m*x*#, -n*# + m*x*SphericalBesselJ[n - 1, m*x]} &@  SphericalBesselJ[n, m*x];
        {psiX, dpsiX} = {x*#, -n*# + x*SphericalBesselJ[n - 1, x]} &@ SphericalBesselJ[n, x];
        {xiX, dxiX} = {psiX, dpsiX} + I*{x*#, -n*# + x*SphericalBesselY[n - 1, x]} &@ SphericalBesselY[n, x];

        an = (m*psiMX*dpsiX - psiX*dpsiMX)/(m*psiMX*dxiX - xiX*dpsiMX);
        bn = (psiMX*dpsiX - m*psiX*dpsiMX)/(psiMX*dxiX - m*xiX*dpsiMX);
  {an, bn}]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieCoefficientInt[n_, x_, m_] :=
    Module[{cn, dn, psiMX, dpsiMX, psiX, dpsiX, xiX, dxiX},
        If[! Apply[Equal, Length /@ {x, m, 0}],
            Return[Message[MieCoefficient::list, x, m]]];

        {psiMX, dpsiMX} = {m*x*#, -n*# + m*x*SphericalBesselJ[n - 1, m*x]} &@ SphericalBesselJ[n, m*x];
        {psiX, dpsiX} = {x*#, -n*# + x*SphericalBesselJ[n - 1, x]} &@ SphericalBesselJ[n, x];
        {xiX, dxiX} = {psiX, dpsiX} + I*{x*#, -n*# + x*SphericalBesselY[n - 1, x]} &@ SphericalBesselY[n, x];

        cn = (-m*dxiX*psiX + m * xiX * dpsiX)/( m*xiX*dpsiMX - psiMX*dxiX);
        dn = (m*dxiX*psiX - m * xiX * dpsiX )/(m*psiMX*dxiX - xiX*dpsiMX);
  {cn, dn}]
(* --------------------------------------------------------------------------------------------- *)



(* ---------------------------------------------------------------------------------------------
Mie Pi and Tau Angular functions
	costh -> Cosine of the angle of incidence 
    n -> Multipole orders 
This functions threads over both variables
------------------------------------------------------------------------------------------------*)
MiePi::usage = "Angular function MiPi[Cos[th]] = LegendreP[1,1,costh]/Sqrt[1-costh^2]] developed in the Mie Theory Formalism"
MieTau::usage = "Angular function MiTau[Cos[th]] = D[LegendreP[1,1,Cos[th]],th] developed in the Mie Theory Formalism"
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
MiePi[n_, costh_] :=
    Module[{f},
        f[0] = 0;
        f[1] = 1;
        f[i_] := f[i] = ((2 i - 1)/(i - 1)) * costh * f[i - 1] - (i/(i - 1)) * f[i - 2];
    f[n]]
SetAttributes[MiePi, Listable]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieTau[n_, costh_] :=
    Module[{f}, f[0] = 0; f[1] = 1;
        f[i_] :=
        f[i] = ((2 i - 1)/(i - 1)) * costh *f[i - 1] - (i/(i - 1)) * f[i - 2];
    n* costh *f[n] - (n + 1)*f[n - 1]]
SetAttributes[MieTau, Listable]
(* --------------------------------------------------------------------------------------------- *)
    


(* ---------------------------------------------------------------------------------------------
Mie Scattering Amplitude Matrix Elements
	amgle -> Angle of incidence 
	x -> Size Parameter
    m -> refractive indices contrast
This function does not thread over any arguement.
------------------------------------------------------------------------------------------------*)
MieScatteringAmplitude12::usage = "Scattering Amplitude Matrix elements S1 and S2 for a spherical scatterer. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
(* --------------------------------------------------------------------------------------------- *)
MieScatteringAmplitude12[x_, m_, angle_] :=
    Module[{ab, poles, pitau, coeff,  s1, s2},
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]]; (*Wacombe criteria for convergence*)

        ab = MieCoefficient[poles, x, m];
        pitau = Through[{MiePi, MieTau}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;

        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  {s1, s2}]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteringAmplitude12[x_, m_, angle_, pole_] :=
    Module[{ab, poles, pitau, coeff,  s1, s2},
		If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Pole orders must be integers*)    
        If[Length[pole] == 0, poles = Range[pole], poles = pole];

        ab = MieCoefficient[poles, x, m];
        pitau = Through[{MiePi, MieTau}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;

        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  {s1, s2}]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteringAmplitude012[x_, m_, angle_] :=
    Module[{ab, poles, pitau, coeff,  s0, s1, s2},
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]]; (*Wacombe criteria for convergence*)

        ab = MieCoefficient[poles, x, m];
        pitau = Through[{MiePi, MieTau}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;
		
		s0 = 0.5 * Plus @@ (((2.*#+1) &/@ poles) * Plus @@ ab);
        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  {s0, s1, s2}]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteringAmplitude012[x_, m_, angle_, pole_] :=
    Module[{ab, poles, pitau, coeff,  s0, s1, s2},
		If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Pole orders must be integers*)    
        If[Length[pole] == 0, poles = Range[pole], poles = pole];

        ab = MieCoefficient[poles, x, m];
        pitau = Through[{MiePi, MieTau}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;

		s0 = 0.5 * Plus @@ (((2.*#+1) &/@ poles) * Plus @@ ab);
        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  {s0, s1, s2}]
(* --------------------------------------------------------------------------------------------- *)



(* ---------------------------------------------------------------------------------------------
MieScatteringQ
	indices -> {nMatrix,nParticle}
This function does not thread over any arguement.
------------------------------------------------------------------------------------------------*)
MieScatteringQ::usage = "Scattering efficiency of a spherical scatterer. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
(* --------------------------------------------------------------------------------------------- *)
MieScatteringQ[indices_, wlength_, radius_] :=
	Module[{ab, sum, x, poles},
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the normSquared, drop zero-imaginary-part and the perform the sum*)
        ab = Plus @@ (Chop[#*Conjugate[#]] &@ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteringQ[indices_, wlength_, radius_, pole_] :=
	Module[{ab, sum, x, poles},
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Pole orders must be integers*)

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0, poles = Range[pole], poles = pole ];

        (*From right to left -> Calculate Mie Coefficients, ge the normSquared, drop zero-imaginary-part and the perform the sum*)
        ab = Plus @@ (Chop[#*Conjugate[#]] &@ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)
(* --------------------------------------------------------------------------------------------- *)



(* ---------------------------------------------------------------------------------------------
MieExtinctionQ
	indices -> {nMatrix,nParticle}
This function does not thread over any arguement.
------------------------------------------------------------------------------------------------*)
MieExtinctionQ::usage = "Extinction efficiency of a spherical scatterer. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
(* --------------------------------------------------------------------------------------------- *)
MieExtinctionQ[indices_, wlength_, radius_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( Re @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)

(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieExtinctionQ[indices_, wlength_, radius_, pole_] :=
    Module[{ab, sum, x, poles},
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Pole orders must be integers*)

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0, poles = Range[pole], poles = pole ];

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( Re @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent cWhiteoefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)
(* --------------------------------------------------------------------------------------------- *)


(* ---------------------------------------------------------------------------------------------
MieAbsorptionQ
	indices -> {nMatrix,nParticle}
This function does not thread over any arguement.
------------------------------------------------------------------------------------------------*)
MieAbsorptionQ::usage = "Absorption efficiency of a spherical scatterer. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
(* --------------------------------------------------------------------------------------------- *)
MieAbsorptionQ[indices_, wlength_, radius_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]] ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( (Re[#] - Chop[# * Conjugate[#]])& @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieAbsorptionQ[indices_, wlength_, radius_, pole_] :=
    Module[{ab, sum, x, poles},
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Pole orders must be integers*)

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0, poles = Range[pole], poles = pole];

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( (Re[#] - Chop[#*Conjugate[#]])& @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)
(* --------------------------------------------------------------------------------------------- *)



(* ---------------------------------------------------------------------------------------------
MieVectorSphericalHarmonics
	super -> Only valid valued: 1 and 3. This corresponds to the superscritp dennoting the SphericalBesselJ or the SphericalHankelH1 functions.
	nn -> Multipolar order(s).
	point -> Evaluation grid in the physical space for the spherical harmonic.
	wavenumber -> Magnitud of the wave vector in the media where the spherical harmonic is evaluated.
This function does not thread over any arguement.
------------------------------------------------------------------------------------------------*)
MieVectorSphericalHarmonic::PhysicalSolution = "The only radial dependency in the Mie Theory is given by SphericalBesselJ (1) and the SphericalHankelH1 (3)."
MieVectorSphericalHarmonic::CoordinateSystem = "The coordinate system must be either the cartesian or the spherical coordinate system."
MieVectorSphericalHarmonic::VectorBase = "The vector base system must be either the cartesian or spherical canonical vector system."
MieVectorSphericalHarmonic::Int = "The parameter `1 must be an integer or a list of integers."
 
MieVectorSphericalHarmonicMo1::usage = "Vector Spherical Harmonic M (Odd, m = 1, ell = n)."
MieVectorSphericalHarmonicMe1::usage = "Vector Spherical Harmonic M (Even, m = 1, ell = n)."
MieVectorSphericalHarmonicNo1::usage = "Vector Spherical Harmonic N (Odd, m = 1, ell = n)."
MieVectorSphericalHarmonicNe1::usage = "Vector Spherical Harmonic N (Even, m = 1, ell = n)."
(* --------------------------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------------------------- *)
Options[MieVectorSphericalHarmonicMo1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicMo1[super_Integer, nn_ , point_List, wavenumber_, OptionsPattern[]]:=
	Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi},
		n = toMap[nn];
		If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]]; (*Polar order must be an integer*)
		
		{etheta, ephi} = Transpose[Drop[IdentityMatrix[3],1]& /@ n]; (*Vectors in the  cannonical vector spherical base*)
		bessel = Switch[super, 
						1, SphericalBesselJ,
						3, SphericalHankelH1,
						_, Return[Message[MieVectorSphericalHarmonic::PhysicalSolution]]];
		
		{rho, theta, phi} = Switch[OptionValue["InputCoordinateSystem"],
									"Cartesian" , CoordinateTransform[ "Cartesian" -> "Spherical", point],
									"Spherical", point,
									_, Return[Message[MieVectorSphericalHarmonic::CoordinateSystem]]];
									
		VSH  = Cos[phi] * MiePi[n, Cos[theta] ] * etheta - Sin[phi] * MieTau[n, Cos[theta] ]* ephi;
		VSH *= bessel[n, wavenumber * rho];
		VSH = Switch[OptionValue["OutputVectorBase"],
					"Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                     "Spherical", VSH,
                     _, Return[Message[MieVectorSphericalHarmonic::VectorBase]]];
    If[Length[nn] == 0, VSH[[1]], VSH ]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
Options[MieVectorSphericalHarmonicMe1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicMe1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
	Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi},
		n = toMap[nn];
		If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]]; (*Polar order must be an integer*)
		
		{etheta, ephi} = Transpose[Drop[IdentityMatrix[3],1]& /@ n]; (*Vectors in the  cannonical vector spherical base*)
		bessel = Switch[super, 
						1, SphericalBesselJ,
						3, SphericalHankelH1,
						_, Return[Message[MieVectorSphericalHarmonic::PhysicalSolution]]];
		
		{rho, theta, phi} = Switch[OptionValue["InputCoordinateSystem"],
									"Cartesian" , CoordinateTransform[ "Cartesian" -> "Spherical", point],
									"Spherical", point,
									_, Return[Message[MieVectorSphericalHarmonic::CoordinateSystem]]];
		VSH  = -Sin[phi] * MiePi[n, Cos[theta]] * etheta - Cos[phi] * MieTau[n, Cos[theta]]* ephi;
		VSH *= bessel[n, wavenumber * rho];
		VSH = Switch[OptionValue["OutputVectorBase"],
					"Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                    "Spherical", VSH,
                    _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0, VSH[[1]], VSH ]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
Options[MieVectorSphericalHarmonicNo1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicNo1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
	Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi, er},
		n = toMap[nn];
		If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]]; (*Polar order must be an integer*)
		
		{er, etheta, ephi} = Transpose[IdentityMatrix[3]& /@ n];
		bessel = Switch[super, 
						1, SphericalBesselJ,
						3, SphericalHankelH1,
						_, Return[Message[MieVectorSphericalHarmonic::PhysicalSolution]]];
		{rho, theta, phi} = Switch[OptionValue["InputCoordinateSystem"],
									"Cartesian" , CoordinateTransform[ "Cartesian" -> "Spherical", point],
									"Spherical", point,
									_, Return[Message[MieVectorSphericalHarmonic::CoordinateSystem]]];
		VSH  = Sin[phi] * MieTau[n, Cos[theta]] * etheta + Cos[phi] * MiePi[n, Cos[theta]] * ephi;
	    VSH *= (-n * bessel[n, #] + # * bessel[n-1, #])&[wavenumber * rho] ;
	    VSH += Sin[phi] * bessel[n, wavenumber * rho] *(n*(n + 1)) * MiePi[n, Cos[theta]] * Sin[theta] * er;
	    VSH /= (wavenumber * rho);
		VSH = Switch[OptionValue["OutputVectorBase"],
					"Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                    "Spherical", VSH,
                    _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0, VSH[[1]], VSH ]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
Options[MieVectorSphericalHarmonicNe1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicNe1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
	Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi, er},
		n = toMap[nn];
		If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]]; (*Polar order must be an integer*)
		
		{er, etheta, ephi} = Transpose[IdentityMatrix[3]& /@ n];
		bessel = Switch[super, 
						1, SphericalBesselJ,
						3, SphericalHankelH1,
						_, Return[Message[MieVectorSphericalHarmonic::PhysicalSolution]]];
		{rho, theta, phi} = Switch[OptionValue["InputCoordinateSystem"],
									"Cartesian" , CoordinateTransform[ "Cartesian" -> "Spherical", point],
									"Spherical", point,
									_, Return[Message[MieVectorSphericalHarmonic::CoordinateSystem]]];
		VSH  = Cos[phi] * MieTau[n, Cos[theta] ] * etheta -  Sin[phi] * MiePi[n, Cos[theta] ]* ephi;
        VSH *= (-n * bessel[n, #] + # * bessel[n-1, #])&[wavenumber * rho] ;
        VSH += Cos[phi] * bessel[n, wavenumber * rho] *(n*(n + 1)) * MiePi[n, Cos[theta]] * Sin[theta] * er;
        VSH /= (wavenumber * rho);
		VSH = Switch[OptionValue["OutputVectorBase"],
					"Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                    "Spherical", VSH,
                    _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0, VSH[[1]], VSH ]]
(* --------------------------------------------------------------------------------------------- *)



(* ---------------------------------------------------------------------------------------------
MieField
	initialMesh -> Evaluating mesh: each element is a 3D vector
	indices -> {nParticcle, nMAtrix}
	wlength
	radius
This function asks for slacar quantities except for the initial mesh
------------------------------------------------------------------------------------------------*)
MieField::Mesh = "The mesh `1` is not an array of 3D vectors"
MieField::Parallel = "The Paralellize option `1` must be either True or False"
MieField::Inside = "Field is set to zero since the evaluating point `1` is smaller than the radius `2`"
MieField::Outsie = "Field is set to zero since the evaluating point `1` is larger than the radius `2`"
 
MieScatteredElectricField::usage = "Scattered  Electric Field. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
MieInternalElectricField::usage = "Scattered  Electric Field. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."	
MieScatteredHField::usage = "Scattered  H Field. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."
MieInternalHField::usage = "Scattered  H Field. The  multipolar contributions to the result can be chosen by Default (Wacombe Criteria) or up to the 'pole' \
									contribution, and if 'pole' is a list,only those multioples are considered."																	
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
Options[MieScatteredElectricField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
MieScatteredElectricField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Nen1 = map[If[Norm[#] < radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[3,#,1.]] &, mesh,{-2}];
          Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
          Mon1 = map[ If[Norm[#] < radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[3,#,1.]]&, mesh,{-2}];
          Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)
          
          field = {I, -1} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteredElectricField[initialMesh_, indices_, wlength_, radius_, pole_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		If[! Apply[Equal, IntegerQ /@ toMap[pole]], Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Polar integer must be integers*)
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          If[Length[pole] == 0, poles = Range[pole], poles = pole];
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Nen1 = map[If[Norm[#] < radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[3,#,1.]] &, mesh,{-2}];
          Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
          Mon1 = map[ If[Norm[#] < radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[3,#,1.]]&, mesh,{-2}];
          Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)
          
          field = {I, -1} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)


(* --------------------------------------------------------------------------------------------- *)
Options[MieInternalElectricField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
MieInternalElectricField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map, m},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          m = Divide @@ Reverse[indices];
          poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Nen1 = map[If[Norm[#] > radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[1,#, m]] &, mesh,{-2}];
          Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
          Mon1 = map[ If[Norm[#] > radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[1,#, m]]&, mesh,{-2}];
          Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)
          
          field = {1, -I} * MieCoefficientInt[poles, x, Divide @@ Reverse[indices]] * {Mon1, Nen1} ; (*[[ 2 = {cM,-idN}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieInternalElectricField[initialMesh_, indices_, wlength_, radius_, pole_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map, m},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		If[! Apply[Equal, IntegerQ /@ toMap[pole]], Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Polar integer must be integers*)
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          m = Divide @@ Reverse[indices];
          If[Length[pole] == 0, poles = Range[pole], poles = pole];
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Nen1 = map[If[Norm[#] > radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[1,#,m]] &, mesh,{-2}];
          Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
          Mon1 = map[ If[Norm[#] > radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[1,#,m]]&, mesh,{-2}];
          Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)
          
          field = {1, -I} * MieCoefficientInt[poles, x, Divide @@ Reverse[indices]] * {Mon1, Nen1} ; (*[[ 2 = {cM,-idN}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)




(* --------------------------------------------------------------------------------------------- *)
Options[MieScatteredHField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
MieScatteredHField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Men1 = map[If[Norm[#] < radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[3,#,1.]] &, mesh,{-2}];
          Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
          Non1 = map[ If[Norm[#] < radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[3,#,1.]]&, mesh,{-2}];
          Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)
          
          field = {-1, -I} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Men1, Non1 } ; (*[[ 2 = {-aM,-ibN}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieScatteredHField[initialMesh_, indices_, wlength_, radius_, pole_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		If[! Apply[Equal, IntegerQ /@ toMap[pole]], Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Polar integer must be integers*)
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          If[Length[pole] == 0, poles = Range[pole], poles = pole];
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Men1 = map[If[Norm[#] < radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[3,#,1.]] &, mesh,{-2}];
          Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
          Non1 = map[ If[Norm[#] < radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[3,#,1.]]&, mesh,{-2}];
          Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)
          
          field = {-1, -I} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Men1, Non1 } ; (*[[ 2 = {-aM,-ibN}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)



(* --------------------------------------------------------------------------------------------- *)
Options[MieInternalHField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
MieInternalHField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map,m},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          m = Divide @@ Reverse[indices];
          poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Men1 = map[If[Norm[#] > radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[1,#,m]] &, mesh,{-2}];
          Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
          Non1 = map[ If[Norm[#] > radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[1,#,m]]&, mesh,{-2}];
          Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)
          
          field = {-I, -1} * MieCoefficientInt[poles, x, Divide @@ Reverse[indices]] * {Non1,Men1 } ; (*[[ 2 = {-icN,-dM}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)
(* --------------------------------------------------------------------------------------------- *)
MieInternalHField[initialMesh_, indices_, wlength_, radius_, pole_,  OptionsPattern[]]:=
	Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map,m},
		If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
		If[Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
		If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}], Return[Message[MieCoefficient::list, wlength, radius]]];
		If[! Apply[Equal, IntegerQ /@ toMap[pole]], Return[Message[MieCoefficient::Int, #]] & /@ toMap[pole]]; (*Polar integer must be integers*)
		
		mesh = Switch[OptionValue["InputCoordinateSystem"],
						"Cartesian", mesh,
                        "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}] ];
         map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]];
          
          x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
          m = Divide @@ Reverse[indices];
          If[Length[pole] == 0, poles = Range[pole], poles = pole];
          coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
          
          evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
          
          Men1 = map[If[Norm[#] > radius,
                      ConstantArray[{0,0,0}, Length[poles]],
                      evenVSH[1,#,m]] &, mesh,{-2}];
          Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
          Non1 = map[ If[Norm[#] > radius,
                       ConstantArray[{0,0,0}, Length[poles]],
                       oddVSH[1,#,m]]&, mesh,{-2}];
          Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)
          
          field = {-I, -1} * MieCoefficientInt[poles, x, Divide @@ Reverse[indices]] * {Non1,Men1 } ; (*[[ 2 = {-icN,-dM}, 7, 3, grid ]]  *)
          field = Plus @@ (coeff * Plus @@ field);  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
    Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]]
(* --------------------------------------------------------------------------------------------- *)

