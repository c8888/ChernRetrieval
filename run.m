
(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: c8888 *)
(* :Date: 2016-12-11 *)

Needs["wavefunction`"];
Needs["HIOER`"];
Needs["space`"];
Needs["protocoling`"];
Needs["chernCalc`"];

t1 = DateList[];
protocolAdd[ToString[t1] <> " Program started."];
(**************************************************************)
\[Delta]x = 0.1;
\[Delta]y = 0.1;
q = 2; (* Pi-flux *)
xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;
RTF = 1;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
nIterations = 10;
nRepeats = 1;
nHIO = 20;
gamma = 0.9;
npts = 4;(*points in the 1st Brillouin zone*)
(**************************************************************)

lat = latticeProbingPoints[xmin, xmax, ymin,
  ymax, \[Delta]x, \[Delta]y];
rec = rectLatticeSites[a, RTF, xmin, xmax, ymin, ymax];
pos = rectLatticeSitesPos[lat, a, \[Delta]x, \[Delta]y, RTF];
neighpos =
    rectLatticeSitesNeighbourhood[pos, lat, \[Delta]x, \[Delta]y,
      rangeNeighbour];
BZ = latticeProbingPointsBZ[npts, a, q];
support = Map[If[Norm[#] < RTF, 1, 0] &, lat, {2}];

(**************************************************************)
\[Beta] =
    wannierNormalisationFactor[\[Sigma]w, \[Delta]x, \[Delta]y, lat];

(**************************************************************)



(**************************************************************)
(*ckModelBZ =
    ParallelMap[
      findCkModel[#, J, J1, lat, a, rec, RTF, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}];
FxyTModel = FxyT[ ckModelBZ ];
wModel = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTModel];
protocolAdd["wModel = " <> ToString[wModel]];
*)(*PutAppend[ "wModel = " <> ToString[wModel], "out/" <> ToString[$ProcessID] <> "protocol.txt"];*)(*
*)(*Export["out/ckModelBZ.dat", ckModelBZ];*)(*

*)(**************************************************************)(*

ckRetrBZ =
    ParallelMap[
      findCkRetr[#, J, J1, lat, a, rec, RTF, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}];
FxyTRetr = FxyT[ ckRetrBZ[[All, All, 1]] ];
wRetr = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTRetr];
protocolAdd[ "wRetr = " <> ToString[wRetr] <> " with mean overlap: " <> ToString[Mean@Flatten@ckRetrBZ[[All, All, 2]] ] <> " stDev: " <> ToString[StandardDeviation@Flatten@ckRetrBZ[[All, All, 2]] ] ];
*)(*Export["out/ckRetrBZ.dat", ckRetrBZ];
Export["out/" <> ToString[$ProcessID] <> "retr_overlap.pdf", ListPlot@Flatten@ckRetrBZ[[All,All,2]]];
(**************************************************************)
*)

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
