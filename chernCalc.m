(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: chernCalc *)
(* :Context: chernCalc` *)
(* :Author: c8888 *)
(* :Date: 2016-12-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["chernCalc`"]
(* Exported symbols added here with SymbolName::usage *)
Needs["wavefunction`"];
Needs["HIOER`"];
Needs["space`"];
Needs["protocoling`"];

findCkRetr::usage =
    "findCkRetr[k_, J_, J1_, lat_, a_, rec_, RTF_, support_,nIterations_,nRepeats_,nHIO_,gamma_, pos_, neighpos_,\[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_] automates the process of finding ckrets for different k"
findCkModel::usage =
    "findCkModel[k_, J_, J1_, lat_, a_, rec_, RTF_, support_,nIterations_,nRepeats_,nHIO_,gamma_, pos_, neighpos_,\[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_] automates the process of finding ckrets for different k"
link::usage =
link::usage =
    "link[ckBZ_, i1_, j1_, i2_, j2_] returns a linking variable between two sites in Brillouin zone {i1,j1} and {i2,j2}"
FxyT::usage =
    "FxyT[ckBZ] returns a table with Fxy defined for direct calculation of Chern number"


Begin["`Private`"]

findCkRetr[k_, J_, J1_, lat_, a_, rec_, RTF_, support_,nIterations_,nRepeats_,nHIO_,gamma_, pos_, neighpos_,\[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_] := Module[{
  wavef =
      waveFunctionHarper[lat, a, J, J1, rec, RTF,
        k, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y],
  FTwavefAbs = {},
  ckRetr = {},
  ckRetrMirror = {},
  ckModel = {},
  overlapRetr,
  overlapRetrMirror
},
  ckModel = wannierBaseRectProject[wavef, lat, rec, pos,
    neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF];

  FTwavefAbs = Abs@Fourier@wavef;
  wavef = phaseRetrieveSupport[FTwavefAbs, support,
    nIterations, nRepeats, nHIO, gamma]; (*not to waste the memory*)

  ckRetr = wannierBaseRectProject[wavef, lat, rec, pos,
    neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF];

  ckRetrMirror = wannierBaseRectProject[mirrorXY@wavef, lat, rec, pos,
    neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF];

  overlapRetr = overlapWannier[ckModel, ckRetr];
  overlapRetrMirror = overlapWannier[ckModel, ckRetrMirror];
If[overlapRetr >= overlapRetrMirror,
  Return[{ckRetr, overlapRetr}],
  Return[{ckRetrMirror, overlapRetrMirror}]
];
  (* Selecting the right reflection of the image which is NECESSARY for the algorithm to calculate the Chern Number.
In real experiment the imperfections will probably not allow for such degeneracy. Anyway it needs to be analyzed further.*)
]

findCkModel[k_, J_, J1_, lat_, a_, rec_, RTF_, support_, nIterations_,
  nRepeats_, nHIO_, gamma_, pos_,
  neighpos_, \[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_] :=
    Module[{wavef =
        waveFunctionHarper[lat, a, J, J1, rec, RTF,
          k, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y]
    },
      wannierBaseRectProject[wavef, lat, rec, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF]
    ]

link[ckBZ_, i1_, j1_, i2_, j2_] :=
    Exp[I Arg@ComplexDotProduct[ckBZ[[i2, j2]], ckBZ[[i1, j1]] ] ]

FxyT[ckBZ_] := Table[Log[
  link[ckBZ, i, j, i + 1, j]*
      link[ckBZ, i + 1, j, i + 1, j + 1]*
      link[ckBZ, i + 1, j + 1, i, j + 1]*
      link[ckBZ, i, j + 1, i, j]
], {i, 1, Dimensions[ckBZ][[1]] - 1}, {j, 1, Dimensions[ckBZ][[2]] - 1}]

End[] (* `Private` *)

EndPackage[]