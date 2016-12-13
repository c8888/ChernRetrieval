#!/usr/local/bin/WolframScript -script

(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: c8888 *)
(* :Date: 2016-12-11 *)

Needs["wavefunction`"];
Needs["HIOER`"];
Needs["space`"];
Needs["protocoling`"];
Needs["chernCalc`"];

PutAppend[ToString[DateList[]] <> "Program started.", "out/" <> ToString[$ProcessID] <> "protocol.txt"];
(**************************************************************)
\[Delta]x = 0.1;
\[Delta]y = 0.1;
q = 2; (* Pi-flux *)
xmin = -8;
xmax = 8;
ymin = -8;
ymax = 8;
RTF = 6;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 0}; (* I know that value from experimental setup. I also \
select the energy band here (let's say its the lowest energy band). I \
don't know the tunneling amplitude ratio but can guess it *)
J = 1;
J1 = 2;
nIterations = 300;
nRepeats = 3;
nHIO = 20;
gamma = 0.9;
npts = 8;(*points in the 1st Brillouin zone*)
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
ckModelBZ =
    ParallelMap[
      findCkModel[#, J, J1, lat, a, rec, RTF, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}];
FxyTModel = FxyT[ckModelBZ];
wModel = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTModel];
PutAppend[ "wModel = " <> ToString[wModel], "out/" <> ToString[$ProcessID] <> "protocol.txt"];
Export["out/ckModelBZ.dat", ckModelBZ];

(**************************************************************)

ckRetrBZ =
    ParallelMap[
      findCkRetr[#, J, J1, lat, a, rec, RTF, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}];
FxyTRetr = FxyT[ckRetrBZ];
wRetr = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTRetr];
PutAppend[ "wRetr = " <> ToString[wRetr], "out/" <> ToString[$ProcessID] <> "protocol.txt"];
Export["out/ckRetrBZ.dat", ckRetrBZ];


