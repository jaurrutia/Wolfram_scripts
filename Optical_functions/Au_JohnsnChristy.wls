#! /usr/bin/env wolframscript
(*
Jonathan Urrutia

Data from: P. B. Johnson and R. W. Christy. Optical constants of the noble metals, Phys. Rev. B 6, 4370-4379 (1972)
            https://doi.org/10.1103/PhysRevB.6.4370
Retrieve from: https://refractiveindex.info/?shelf=main&book=Au&page=Johnson
*)

hPlanck = 4.135668*^-15 (*eV s  - Planks constant*);
cSpeedLight = 3*^17 (*nm s^-1   -   Speed of light*);
hcPlanckSpeed = hPlanck*cSpeedLight;


(*Drude model functions *)

    DrudeEps::usage = "Dielectric function given bty the Drudel Model as function of the incident wavelength."
    DrudeEps[wlen_, wlenPlasma_, wlenDamp_] := 1 - (wlen / wlenPlasma) ^2 * wlenDamp / (wlenDamp + I * wlen )

    DrudeEpsSize::usage = "Dielectric function given bty the Drudel Model with  a size correction on the samping constant due to surface effects, as function of the incident wavelength."
    DrudeEpsSize[wlen_, wlenPlasma_, wlenDamp_, wlenDSize_] := 1 - (wlen / wlenPlasma) ^2 * wlenDamp * wlenDSize / (wlenDamp * wlenDSize + I * (wlenDamp + wlenDSize) * wlen )


    DrudeEpsEnergy::usage = "Dielectric function given bty the Drudel Model as function of the photon energy."
    DrudeEpsEnergy[energy_, energyPlasma_, damping_] := 1 - energyPlasma^2 / (energy * (energy + I * damping))

(*Experimental data*)

    data = Rest @ Import["Au_JohnsnChristy_n_ev.csv",  "Data"];  (*[[ ev, Re[n], Im[n]  ]]*)
    JohnsonChristyAuExp::usage = "Experimental data from Johnson & Christy for the Au refractive index as function of the photon energy."
    JohnsonChristyAuExp = data;


(* Interpolation of the values: Cubic and Spline*)

    JohnsonChristyAuRef::usage = "Interpolation of the refractive inedex for Au as function of the incident energy."
    JohnsonChristyAuRef =  Interpolation[ data /. {energy_, realN_, imagN_} -> { hcPlanckSpeed / energy, realN + I *imagN}];

    JohnsonChristyAuEps::usage = "Interpolation of the dielectric function for Au as function of the incident energy."
    JohnsonChristyAuEps =  Interpolation[ data /. {energy_, realN_, imagN_} -> { hcPlanckSpeed / energy, (realN + I *imagN)^2}];

    JohnsonChristyAuEpsReal = Interpolation[ data /. {energy_, realN_, imagN_} -> { hcPlanckSpeed / energy, realN^2 - imagN^2 }, Method -> "Spline"];
    JohnsonChristyAuEpsImag = Interpolation[ data /. {energy_, realN_, imagN_} -> { hcPlanckSpeed / energy, 2 * realN * imagN } , Method -> "Spline"];
    JohnsonChristyAuEpsSpline[wlen_] := JohnsonChristyAuEpsReal[wlen] + I * JohnsonChristyAuEpsImag[wlen]

(* Size corrected dielectric function and  refractive index as function of the incident wavelength*)


    JohnsonChristyAuEpsSize::usage = "Dielectric function of Au with size correctino for nanospheres, employing the data from Johnson & Christy (1979)."
    JohnsonChristyAuEpsSize[wlen_, radius_] := ( JohnsonChristyAuEpsSpline[wlen] - DrudeEps[wlen, wlenP, gamma] + DrudeEpsSize[wlen, wlenP, gamma, gammaSize ] ) /. { wlenP -> 142.609, gamma -> 14966.2, gammaSize -> (radius*2*Pi*cSpeedLight)/(1.4*^15) }

    JohnsonChristyAuRefSize::usage = "Refractive index of Au with size correctino for nanospheres, employing the data from Johnson & Christy (1979)."
    JohnsonChristyAuRefSize[wlen_, radius_] := Sqrt[JohnsonChristyAuEpsSize[wlen, radius]]

    ClearAll[data];
