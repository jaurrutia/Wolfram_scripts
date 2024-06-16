#! /usr/bin/env wolframscript
(*
Jonathan Urrutia

Useful functions for common tasks

*)

toMap::usage  = "Changes a non-list entry to a list so Map can be used.";
toMap = If[Length[#] == 0, {#}, #] & ;


movAverage::usage = "movAverage[list, n] gives the moving average of a one dimensional list with step n";
movAverage[list_, step_] :=  Plus @@ Partition[list, Length[list] - step, 1]/(step + 1) ;

date::usage = "date[] returns a String with the date in format yyyy/mm/dd === hh:mm:ss.";
date[] :=  StringJoin @ {
              " === " , Riffle[#, " === "], " ===\n"} &@ MapThread[
                      StringJoin @ Riffle[#1, #2] &, {
                          Partition[ToString /@ Date[], 3] , {"/", ":"}
              }];
