#! /usr/bin/env wolframscript
(*
Jonathan Urrutia
*)

(************************
 General use functions
 ************************)
toMap::usage  = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;

CartesianToSpherical::usage = "Transformation Matrix if the cartesian canonical vector base to the canonical spherical vector base (physisist convention)."
CartesianToSpherical = {{Sin[#1]*Cos[#2],Sin[#1]*Sin[#2], Cos[#1]},
                        {Cos[#1]*Cos[#2],Cos[#1]*Sin[#2], -Sin[#1]},
                        {-Sin[#2],Cos[#2], 0}}& ;

SphericalToCartesian::usage = "Transformation Matrix of the spherical canonical vector base to the canonical cartesina vector base (physisist convention)."
SphericalToCartesian = Transpose[CartesianToSpherical[#1,#2]]& ;



(************************
 Mie functions, general messages
 ************************)
 MieCoefficient::list = "Both `1`and `2` must be scalar quantities.";
 MieCoefficient::int = "Parameter `1` must be an integer number.";

(************************************
 MieCoefficients
                    x -> Size Parameter
                    m -> refractive incdices contrast
                    Dimensions[x]  == Dimensions[m] == 0
                    n -> Multipole orders // threads over this variable
*************************************)
MieCoefficient::usage = "MieCoefficient[n,x,m] give the Mie Coefficients for a given size parameter x and a relative refractive index ; this functions threads over n.";
MieCoefficient[n_, x_, m_] :=
    Module[{an, bn, psiMX, dpsiMX, psiX, dpsiX, xiX, dxiX},
        If[! Apply[Equal, Length /@ {x, m, 0}],
            Return[Message[MieCoefficient::list, x, m]]
        ];

        {psiMX, dpsiMX} = {m*x*#, -n*# + m*x*SphericalBesselJ[n - 1, m*x]} &@  SphericalBesselJ[n, m*x];
        {psiX, dpsiX} = {x*#, -n*# + x*SphericalBesselJ[n - 1, x]} &@ SphericalBesselJ[n, x];
        {xiX, dxiX} = {psiX, dpsiX} + I*{x*#, -n*# + x*SphericalBesselY[n - 1, x]} &@ SphericalBesselY[n, x];

        an = (m*psiMX*dpsiX - psiX*dpsiMX)/(m*psiMX*dxiX - xiX*dpsiMX);
        bn = (psiMX*dpsiX - m*psiX*dpsiMX)/(psiMX*dxiX - m*xiX*dpsiMX);
  {an, bn}]


  (************************************
  MieCoefficientInt
                    x -> Size Parameter
                    m -> refractive incdices contrast
                    Dimensions[x]  == Dimensions[m] == 0
                    n -> Multipole orders // threads over this variable
  *************************************)
MieCoefficientInt::usage = "MieCoefficientInt[n,x,m] give the oefficients for a given size parameter x and a relative refractive index for the internal electric field; this functions threads over n.";
MieCoefficientInt[n_, x_, m_] :=
    Module[{cn, dn, psiMX, dpsiMX, psiX, dpsiX, xiX, dxiX},
        If[! Apply[Equal, Length /@ {x, m, 0}],
            Return[Message[MieCoefficient::list, x, m]]
        ];

        {psiMX, dpsiMX} = {m*x*#, -n*# + m*x*SphericalBesselJ[n - 1, m*x]} &@ SphericalBesselJ[n, m*x];
        {psiX, dpsiX} = {x*#, -n*# + x*SphericalBesselJ[n - 1, x]} &@ SphericalBesselJ[n, x];
        {xiX, dxiX} = {psiX, dpsiX} + I*{x*#, -n*# + x*SphericalBesselY[n - 1, x]} &@ SphericalBesselY[n, x];

        cn = (m*dxiX*psiX - m * xiX*dpsiX)/(psiMX*dxiX - m*xiX*dpsiMX);
        dn = (m*dpsiMX*dpsiX - m * dxiX * dpsiX )/(m*psiMX*dxiX - xiX*dpsiMX);
  {cn, dn}]


(************************************
MiePi & MieTau: Angular functions
*************************************)

(************************************
MiePi
            n -> Multipole order (list or integer)
            costh -> Cosine of the angles to evaluate (Float or list)
*************************************)
MiePi::usage = "Angular function MiPi[Cos[th]] = LegendreP[1,1,costh]/Sqrt[1-costh^2]] developed in the Mie Theory Formalism"
MiePi[n_, costh_] :=
    Module[{f},
        f[0] = 0;
        f[1] = 1;
        f[i_] := f[i] = ((2 i - 1)/(i - 1)) * costh * f[i - 1] - (i/(i - 1)) * f[i - 2];
    f[n]]
SetAttributes[MiePi, Listable]

(************************************
MieTau
            n -> Multipole order (list or integer)
            costh -> Cosine of the angles to evaluate (Float or list)
*************************************)
MieTau::usage = "Angular function MiTau[Cos[th]] = D[LegendreP[1,1,Cos[th]],th] developed in the Mie Theory Formalism"
MieTau[n_, costh_] :=
    Module[{f}, f[0] = 0; f[1] = 1;
        f[i_] :=
        f[i] = ((2 i - 1)/(i - 1)) * costh *f[i - 1] - (i/(i - 1)) * f[i - 2];
    n* costh *f[n] - (n + 1)*f[n - 1]]
SetAttributes[MieTau, Listable]


(************************************
MieScatteringAmplitude12
            x -> Size parameter,
            m-> contrast,
            anle -> angle
*************************************)
MieScatteringAmplitude12[x_, m_, angle_] :=
    Module[{ab, poles, pitau, coeff,  s1, s2},
        (*Wacombe criteria for convergence*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];

        ab = MieCoefficient[poles, x, m];
        pitau = Through[{MiePi, MieTau}[poles, Cos[angle]] ];
        coeff = ((2.*# + 1)/((# + 1)*#)) & /@ poles;

        s1 = Plus @@ (coeff * Plus @@ (ab*pitau) );
        s2 = Plus @@ (coeff * Plus @@ (ab*Reverse[pitau]) );
  {s1, s2}]



(************************************
MieScatteringQ
indices={nMatrix,nParticle}
*************************************)
MieScatteringQ[indices_, wlength_, radius_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the normSquared, drop zero-imaginary-part and the perform the sum*)
        ab = Plus @@ (Chop[#*Conjugate[#]] &@ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)

MieScatteringQ[indices_, wlength_, radius_, pole_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0,
            poles = Range[pole],
            poles = pole
        ];

        (*From right to left -> Calculate Mie Coefficients, ge the normSquared, drop zero-imaginary-part and the perform the sum*)
        ab = Plus @@ (Chop[#*Conjugate[#]] &@ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)


(************************************
MieExtinctionQ
indices={nMatrix,nParticle}
*************************************)
MieExtinctionQ[indices_, wlength_, radius_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( Re @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)

MieExtinctionQ[indices_, wlength_, radius_, pole_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0,
            poles = Range[pole],
            poles = pole
        ];

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( Re @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)

(************************************
MieAbsorptionQ
indices={nMatrix,nParticle}
*************************************)
MieAbsorptionQ[indices_, wlength_, radius_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( (Re[#] - Chop[# * Conjugate[#]])@ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)

MieAbsorptionQ[indices_, wlength_, radius_, pole_] :=
    Module[{ab, sum, x, poles},
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
            Return[Message[MieCoefficient::list, wlength, radius]]
        ];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
            Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0,
            poles = Range[pole],
            poles = pole
        ];

        (*From right to left -> Calculate Mie Coefficients, ge the real part, and the perform the sum*)
        ab = Plus @@ ( (Re[#] - Chop[#*Conjugate[#]]) @ MieCoefficient[poles, x, Divide @@ Reverse[indices]]);
        sum = Plus @@ ((2.*# + 1 & /@ poles)*ab);       (*Muliply the past value by a pole-order dependent coefficient*)
    sum *= 2./(x^2)]       (*Normalize by the geometric cross section of the sphere*)


(************************************
MieVectorSphericalHarmonics
indices={nMatrix,nParticle}
*************************************)

 MieVectorSphericalHarmonic::PhysicalSolution = "The only radial dependency in the Mie Theory is given by SphericalBesselJ (1) and the SphericalHankelH1 (3)."
 MieVectorSphericalHarmonic::CoordinateSystem = "The coordinate system must be either the cartesian or spherical coordinate system."
 MieVectorSphericalHarmonic::VectorBase = "The vector base system must be either the cartesian or spherical canonical vector system."
 MieVectorSphericalHarmonic::Int = "The parameter `1 must be an integer or a list of int."

 Options[MieVectorSphericalHarmonicMo1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicMo1[super_Integer, nn_ , point_List, wavenumber_, OptionsPattern[]]:=
 Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi},
    n = toMap[nn];
    If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]];
    {etheta, ephi} = Transpose[Drop[IdentityMatrix[3],1]& /@ n];
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
                            _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0,
        VSH[[1]],
        VSH ]]

 Options[MieVectorSphericalHarmonicMe1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicMe1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
 Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi},
    n = toMap[nn];
    If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]];
    {etheta, ephi} = Transpose[Drop[IdentityMatrix[3],1]& /@ n];
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
    If[Length[nn] == 0,
       VSH[[1]],
       VSH ]]

 Options[MieVectorSphericalHarmonicNo1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicNo1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
 Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi, er},
    n = toMap[nn];
    If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]];
    {er, etheta, ephi} = Transpose[IdentityMatrix[3]& /@ n];
    bessel = Switch[super,
                    1, SphericalBesselJ,
                    3, SphericalHankelH1,
                    _, Return[Message[MieVectorSphericalHarmonic::PhysicalSolution]]];

    {rho, theta, phi} = Switch[OptionValue["InputCoordinateSystem"],
                              "Cartesian" , CoordinateTransform[ "Cartesian" -> "Spherical", point],
                              "Spherical", point,
                              _, Return[Message[MieVectorSphericalHarmonic::CoordinateSystem]]];

    VSH  = Sin[phi] * MieTau[n, Cos[theta] ] * etheta +  Cos[phi] * MiePi[n, Cos[theta] ]* ephi;
    VSH *= (-n * bessel[n, #] + # * bessel[n-1, #])&[wavenumber * rho] ;
    VSH += Sin[phi] * bessel[n, wavenumber * rho] *(n*(n + 1)) * LegendreP[n, Cos[theta]] * er;
    VSH /= wavenumber * rho;
    VSH = Switch[OptionValue["OutputVectorBase"],
                             "Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                             "Spherical", VSH,
                             _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0,
         VSH[[1]],
         VSH ]]

 Options[MieVectorSphericalHarmonicNe1] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Cartesian"};
MieVectorSphericalHarmonicNe1[super_Integer, nn_, point_List, wavenumber_, OptionsPattern[]]:=
 Module[{n, bessel, rho, theta, phi, VSH, etheta, ephi, er},
    n = toMap[nn];
    If[ !Apply[And, IntegerQ /@ n], Return[MieVectorSphericalHarmonic::Int,nn]];
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
    VSH += Cos[phi] * bessel[n, wavenumber * rho] *(n*(n + 1)) * LegendreP[n, Cos[theta]] * er;
    VSH /= wavenumber * rho;
    VSH = Switch[OptionValue["OutputVectorBase"],
                             "Cartesian" , Map[Dot[SphericalToCartesian[theta, phi], #]&, VSH],
                             "Spherical", VSH,
                             _, Message[MieVectorSphericalHarmonic::VectorBase]];
    If[Length[nn] == 0,
         VSH[[1]],
         VSH ]]


(************************************
MieVectorSphericalHarmonics
indices={nMatrix,nParticle}
 *************************************)
 MieField::Mesh = "The mesh `1` is not an array of 3D vectors"
 MieField::Parallel = "The Paralellize option `1` must be either True or False"
 MieField::Inside = "Field is set to zero since the evaluating point `1` is smaller than the radius `2`"
 MieField::Outsie = "Field is set to zero since the evaluating point `1` is larger than the radius `2`"


 (*Scattered Electric Field*)
  Options[MieScatteredElectricField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
 MieScatteredElectricField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
  Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
     If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
     If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
     If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
              Return[Message[MieCoefficient::list, wlength, radius]]
     ];
     mesh = Switch[OptionValue["InputCoordinateSystem"],
                              "Cartesian", mesh,
                              "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                              ];
     map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]
                      ];

     x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
     poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
     coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];

     evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
     oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

     Nen1 = map[ If[Norm[#] < radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             evenVSH[3,#,1.]]&, mesh,{-2}];
     Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
     Mon1 = map[ If[Norm[#] < radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             oddVSH[3,#,1.]]&, mesh,{-2}];
     Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)

     field = {I, -1} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
     field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
  Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
 ]


  MieScatteredElectricField[initialMesh_, indices_, wlength_, radius_, pole_, OptionsPattern[]]:=
   Module[{mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
      If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
      If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
      If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
              Return[Message[MieCoefficient::list, wlength, radius]]
      ];
      If[! Apply[Equal, IntegerQ /@ toMap[pole]],
             Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
      ];
      mesh = Switch[OptionValue["InputCoordinateSystem"],
                              "Cartesian", mesh,
                              "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                              ];
      map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"]] ]
                      ];

      x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
      If[Length[pole] == 0,
             poles = Range[pole],
             poles = pole];
      coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];

      evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
      oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

      Nen1 = map[ If[Norm[#] < radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             evenVSH[3,#,1.]]&, mesh,{-2}];
      Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
      Mon1 = map[ If[Norm[#] < radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             oddVSH[3,#,1.]]&, mesh,{-2}];
      Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)

      field = {I, -1} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
      field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
   Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
  ]


   (*Internal Electric Field*)
   Options[MieInternalElectricField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
   MieInternalElectricField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
    Module[{m, mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
      If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
      If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
      If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
               Return[Message[MieCoefficient::list, wlength, radius]]
      ];
      mesh = Switch[OptionValue["InputCoordinateSystem"],
                               "Cartesian", mesh,
                               "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                               ];
      map = Switch[OptionValue["Parallelize"],
                       False, Map,
                       True, ParallelMap,
                       _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                       ];

      x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
      poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
      coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
      m = Divide @@ Reverse[indices];

      evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
      oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

      Nen1 = map[ If[Norm[#] > radius,
                                  ConstantArray[{0,0,0}, Length[poles]],
                              evenVSH[1,#,m]]&, mesh,{-2}];
      Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
      Mon1 = map[ If[Norm[#] > radius,
                                  ConstantArray[{0,0,0}, Length[poles]],
                              oddVSH[1,#,m]]&, mesh,{-2}];
      Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)

      field = {1, -I} * MieCoefficientInt[poles, x, m] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
      field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
   Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
   ]

   MieInternalElectricField[initialMesh_, indices_, wlength_, radius_, pole_, OptionsPattern[]]:=
    Module[{m, mesh, x, poles, coeff, Nen1, Mon1, field,oddVSH,evenVSH, map},
      If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
      If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
      If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
              Return[Message[MieCoefficient::list, wlength, radius]]
      ];
      If[! Apply[Equal, IntegerQ /@ toMap[pole]],
             Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
      ];
      mesh = Switch[OptionValue["InputCoordinateSystem"],
                              "Cartesian", mesh,
                              "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                              ];
      map = Switch[OptionValue["Parallelize"],
                      False, Map,
                      True, ParallelMap,
                      _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                      ];

      x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
      If[Length[pole] == 0,
             poles = Range[pole],
             poles = pole];
      coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
      m = Divide @@ Reverse[indices];

      evenVSH = MieVectorSphericalHarmonicNe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
      oddVSH =  MieVectorSphericalHarmonicMo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

      Nen1 = map[ If[Norm[#] > radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             evenVSH[3,#, m]]&, mesh,{-2}];
      Nen1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Nen1]; (*[[7,3,grid-rev]]*)
      Mon1 = map[ If[Norm[#] > radius,
                                 ConstantArray[{0,0,0}, Length[poles]],
                             oddVSH[3,#, m]]&, mesh,{-2}];
      Mon1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Mon1];  (*[[7,3,grid-rev]]*)

      field = {1, -I} * MieCoefficientInt[poles, x, m] * {Nen1, Mon1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
      field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
   Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
   ]


    (*Magnetic Field: H*)

    (*Scattered H Field*)
     Options[MieScatteredHField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
    MieScatteredHField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
     Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
        If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
        If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
                 Return[Message[MieCoefficient::list, wlength, radius]]
        ];
        mesh = Switch[OptionValue["InputCoordinateSystem"],
                                 "Cartesian", mesh,
                                 "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                                 ];
        map = Switch[OptionValue["Parallelize"],
                         False, Map,
                         True, ParallelMap,
                         _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                         ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
        coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];

        evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
        oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

        Men1 = map[ If[Norm[#] < radius,
                                    ConstantArray[{0,0,0}, Length[poles]],
                                evenVSH[3,#,1.]]&, mesh,{-2}];
        Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
        Non1 = map[ If[Norm[#] < radius,
                                    ConstantArray[{0,0,0}, Length[poles]],
                                oddVSH[3,#,1.]]&, mesh,{-2}];
        Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)

        field = {-1, -I} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Men1, Non1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
        field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
     Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
    ]

    MieScatteredHField[initialMesh_, indices_, wlength_, radius_, pole_, OptionsPattern[]]:=
     Module[{mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
        If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
        If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
        If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
                Return[Message[MieCoefficient::list, wlength, radius]]
        ];
        If[! Apply[Equal, IntegerQ /@ toMap[pole]],
               Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
        ];
        mesh = Switch[OptionValue["InputCoordinateSystem"],
                                "Cartesian", mesh,
                                "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                                ];
        map = Switch[OptionValue["Parallelize"],
                        False, Map,
                        True, ParallelMap,
                        _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                        ];

        x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
        If[Length[pole] == 0,
               poles = Range[pole],
               poles = pole];
        coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];

        evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
        oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

        Men1 = map[ If[Norm[#] < radius,
                                    ConstantArray[{0,0,0}, Length[poles]],
                                evenVSH[3,#,1.]]&, mesh,{-2}];
        Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
        Non1 = map[ If[Norm[#] < radius,
                                    ConstantArray[{0,0,0}, Length[poles]],
                                oddVSH[3,#,1.]]&, mesh,{-2}];
        Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)

        field = {-1, -I} * MieCoefficient[poles, x, Divide @@ Reverse[indices]] * {Men1, Non1 } ; (*[[ 2 = {iaN,-bM}, 7, 3, grid ]]  *)
        field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
     Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
    ]



    (*Internal Electric Field*)
    Options[MieInternalHField] = {"InputCoordinateSystem"->"Cartesian", "OutputVectorBase" -> "Spherical", "Parallelize"-> False};
    MieInternalHField[initialMesh_, indices_, wlength_, radius_,  OptionsPattern[]]:=
     Module[{m, mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
       If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
       If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
       If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
                Return[Message[MieCoefficient::list, wlength, radius]]
       ];
       mesh = Switch[OptionValue["InputCoordinateSystem"],
                                "Cartesian", mesh,
                                "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                                ];
       map = Switch[OptionValue["Parallelize"],
                        False, Map,
                        True, ParallelMap,
                        _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                        ];

       x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
       poles = Range[Ceiling[x + 4.*x^(1./3) + 2.]];(*Wacombe criteria for convergence*)
       coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
       m = Divide @@ Reverse[indices];

       evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
       oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

       Men1 = map[ If[Norm[#] > radius,
                                   ConstantArray[{0,0,0}, Length[poles]],
                               evenVSH[1,#,m]]&, mesh,{-2}];
       Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
       Non1 = map[ If[Norm[#] > radius,
                                   ConstantArray[{0,0,0}, Length[poles]],
                               oddVSH[1,#,m]]&, mesh,{-2}];
       Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)

       field = {-I, -1} * MieCoefficientInt[poles, x, m] * {Non1, Men1 } ; (*[[ 2 = {iaN,-nM}, 7, 3, grid ]]  *)
       field = Plus @@ Plus @@ field;  (*[[ 7}, 3, grid ]] --> , [[3, grid ]]  *)
      Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
    ]

    MieInternalHField[initialMesh_, indices_, wlength_, radius_, pole_, OptionsPattern[]]:=
     Module[{m, mesh, x, poles, coeff, Men1, Non1, field,oddVSH,evenVSH, map},
       If[Length[Dimensions[initialMesh]] == 1, mesh = {initialMesh}, mesh = initialMesh];
       If[ Dimensions[mesh][[-1]] != 3, Return[Message[MieField::Mesh,initialMesh]] ];
       If[! Apply[Equal, Length /@ {wlength, radius, 0., indices[[1]], indices[[2]]}],
               Return[Message[MieCoefficient::list, wlength, radius]]
       ];
       If[! Apply[Equal, IntegerQ /@ toMap[pole]],
              Return[Message[MieCoefficient::int, #]] & /@ toMap[pole]
       ];
       mesh = Switch[OptionValue["InputCoordinateSystem"],
                               "Cartesian", mesh,
                               "Spherical", Map[CoordinateTransform["Spherical" -> "Cartesian" , #]&, mesh, {-2}]
                               ];
       map = Switch[OptionValue["Parallelize"],
                       False, Map,
                       True, ParallelMap,
                       _, Return[Message[MieField::Parallel, OptionValue["Paralellize"] ]]
                       ];

       x = (2.*Pi*radius)*indices[[1]]/wlength;    (*Size parameter*)
       If[Length[pole] == 0,
              poles = Range[pole],
              poles = pole];
       coeff = Map[ (I^#*(2.*#+1.)/(#*(#+1.)) )&, poles];
       m = Divide @@ Reverse[indices];

       evenVSH = MieVectorSphericalHarmonicMe1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;
       oddVSH =  MieVectorSphericalHarmonicNo1[#1, poles, #2, #3 * x/radius, "InputCoordinateSystem" -> OptionValue["InputCoordinateSystem"], "OutputVectorBase" -> OptionValue["OutputVectorBase"]]&;

       Men1 = map[ If[Norm[#] > radius,
                                   ConstantArray[{0,0,0}, Length[poles]],
                               evenVSH[1,#,m]]&, mesh,{-2}];
       Men1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Men1]; (*[[7,3,grid-rev]]*)
       Non1 = map[ If[Norm[#] > radius,
                                   ConstantArray[{0,0,0}, Length[poles]],
                               oddVSH[1,#,m]]&, mesh,{-2}];
       Non1 = Transpose[Transpose[#, Reverse@Range@Length@Dimensions[#]]]&[Non1];  (*[[7,3,grid-rev]]*)

       field = {-I, -1} * MieCoefficientInt[poles, x, m] * {Non1, Men1 } ; (*[[ 2 = {iaN,-nM}, 7, 3, grid ]]  *)
       field = Plus @@ Plus @@ field;  (*[[ 7, 3, grid ]] --> , [[3, grid ]]  *)
     Transpose[#, Reverse@Range@Length@Dimensions[#]]&[field]
    ]
