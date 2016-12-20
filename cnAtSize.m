
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
(*xmin = -8;
xmax = 8;
ymin = -8;
ymax = 8;
RTF = 6;*)
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
nIterations = 500;
nRepeats = 3;
nHIO = 20;
gamma = 0.9;
npts = 5;(*points in the 1st Brillouin zone*)
(**************************************************************)
RTFmin = 2.01; (*min and max must not be integers! TODO verify the bug*)
RTFmax = 10.01;
deltaRTF = 0.5;
RTFRepeats = 5;

margin[RTF_] := 0.3 RTF
xmin[RTF_] := -RTF-margin[RTF]
ymin[RTF_] := -RTF-margin[RTF]
xmax[RTF_] := RTF+margin[RTF]
ymax[RTF_] := RTF+margin[RTF]

protocolAdd["Parameters: "];
protocolAdd["\[Delta]x = 0.1;
\[Delta]y = 0.1;
q = 2; (* Pi-flux *)
(*xmin = -8;
xmax = 8;
ymin = -8;
ymax = 8;
RTF = 6;*)
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
nIterations = 500;
nRepeats = 3;
nHIO = 20;
gamma = 0.9;
npts = 5;(*points in the 1st Brillouin zone*)
(**************************************************************)
RTFmin = 2.01; (*min and max must not be integers! TODO verify the bug*)
RTFmax = 10.01;
deltaRTF = 0.5;
RTFRepeats = 5;

margin[RTF_] := 0.3 RTF"];

protocolAdd["Results: "];

protocolAdd[ "RTF" <> " " <> "Chern_number_retr." <> " " <> "Mean overlap" <> " " <> " Standard deviation of overlap"];

RTFTab = Table[r,{r, RTFmin, RTFmax, deltaRTF}]; (* this table will be probed*)
RTFReport = {};

BZ = latticeProbingPointsBZ[npts, a, q];

Do[ParallelMap[Module[{},
lat = latticeProbingPoints[xmin[#], xmax[#], ymin[#],
  ymax[#], \[Delta]x, \[Delta]y];
rec = rectLatticeSites[a, #, xmin[#], xmax[#], ymin[#], ymax[#]];
pos = rectLatticeSitesPos[lat, a, \[Delta]x, \[Delta]y, #];
neighpos =
    rectLatticeSitesNeighbourhood[pos, lat, \[Delta]x, \[Delta]y,
      rangeNeighbour];

support = Function[rtf, Map[If[Norm[#] < rtf, 1, 0] &, lat, {2}]][#];

(**************************************************************)
\[Beta] =
    wannierNormalisationFactor[\[Sigma]w, \[Delta]x, \[Delta]y, lat];

(**************************************************************)


ckRetrBZ =
    Function[rtf,
      Map[
      findCkRetr[#, J, J1, lat, a, rec, rtf, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}]
    ][#];
FxyTRetr = FxyT[ ckRetrBZ[[All, All, 1]] ];
wRetr = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTRetr];
AppendTo[RTFReport, {#, Re@wRetr, Mean@Flatten@ckRetrBZ[[All, All, 2]], StandardDeviation@Flatten@ckRetrBZ[[All, All, 2]]}];
protocolAdd[ ToString[#] <> " " <> ToString[Re@wRetr] <> " " <> ToString[Mean@Flatten@ckRetrBZ[[All, All, 2]] ]
    <> " " <> ToString[StandardDeviation@Flatten@ckRetrBZ[[All, All, 2]] ] ];
#]
    &,
  RTFTab
];
,
  {RTFRepeats}]

Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtRTF.dat", RTFReport];
Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtRTFPlot.pdf",
  ListPlot[RTFReport[[All,1;;2]]]];


t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
