#! /usr/bin/env wolframscript
(* ::Package:: *)

(* ::Package:: *)

BeginPackage["LSUQuantumPlasmonics`"]

(* ========================================================================= *)
(* 1. PUBLIC SYMBOL DECLARATIONS                                             *)
(* ========================================================================= *)

(* Physical Constants *)
hbar::usage = "hbar is the reduced Planck constant in eV*s or J*s.";
e0::usage = "e0 is the vacuum permittivity in F/m.";
me::usage = "me is the electron rest mass in kg.";

(* Classical vs. Non-Local / Quantum Models *)
g2Analysis::usage = "Calculate the correlation from data extracted in the LSU's Quantum Harp by choosing th two chanels fo the file.";

Begin["`Private`"]

(* ========================================================================= *)
(* 2. PHYSICAL CONSTANTS (SI & eV Units)                                     *)
(* ========================================================================= *)

hbar = 6.582119569*^-16; (* eV s *)
cSpeed = 2.99792458*^17; (* nm s^-1 *)
hc = 2 * Pi * hbar * cSpeed; (* eV nm *)

(* ========================================================================= *)
(* 3. THEORY & NUMERICAL SIMULATION                                          *)
(* ========================================================================= *)



(* ========================================================================= *)
(* 4. EXPERIMENTAL DATA ANALYSIS & SIMULATION FITTING                         *)
(* ========================================================================= *)

g2Analysis[tpuFile_List, ch1_Integer:1, ch2_Integer:1]:=Module[{g2,n1,n2},
  {n1, n2} = Transpose[tpuFile][[{1+ch1,1+ch2}]];
  g2 = N[Apply[Plus,n1*n2]];
  If[ch1==ch2, 
    g2 /= #*# &@ Apply[Plus,n1],
    (*else*)
    g2 /= Apply[Plus,n1] * Apply[Plus,n2]];
  g2 *=  Length[n1]] ;

End[] (* End Private *)

EndPackage[]
