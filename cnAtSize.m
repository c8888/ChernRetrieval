
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
k0 = {0.1, 1.67}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 0.237;
J1 = 0.125;
(*nIterations = 100;
nRepeats = 3;
nHIO = 20;
gamma = 0.9;*)
npts = 8;(*points in the 1st Brillouin zone*)
(**************************************************************)
RTFmin = 2; (*TODO: Bug. RTF must not be an odd multiple of 0.5a*)
RTFmax = 10;
deltaRTF = 1;
RTFRepeats = 1;

margin[RTF_] := 0.3*RTF
xmin[RTF_] := -RTF-margin[RTF]
ymin[RTF_] := -RTF-margin[RTF]
xmax[RTF_] := RTF+margin[RTF]
ymax[RTF_] := RTF+margin[RTF]

protocolBar[];
protocolAdd["Parameters: "];
protocolAdd["\[Delta]x = "<> ToString[\[Delta]x] ];
protocolAdd["\[Delta]y = "<> ToString[\[Delta]y] ];
protocolAdd["q = "<> ToString[q] ];
protocolAdd["xmin = "<> ToString[xmin] ];
protocolAdd["xmax = "<> ToString[xmax] ];
protocolAdd["ymin = "<> ToString[ymin] ];
protocolAdd["ymax = "<> ToString[ymax] ];
protocolAdd["RTF = "<> ToString[RTF] ];
protocolAdd["rangeNeighbour = "<> ToString[rangeNeighbour] ];
protocolAdd["a = "<> ToString[a] ];
protocolAdd["\[Sigma]w = "<> ToString[\[Sigma]w ] ];
protocolAdd["k0 = "<> ToString[k0] ];
protocolAdd["J = "<> ToString[J] ];
protocolAdd["J1 = "<> ToString[J1] ];
protocolAdd["nIterations = "<> ToString[nIterations] ];
protocolAdd["nRepeats = "<> ToString[nRepeats] ];
protocolAdd["nHIO = "<> ToString[nHIO] ];
protocolAdd["gamma = "<> ToString[gamma] ];
protocolAdd["npts = "<> ToString[npts] ];
protocolBar[];
(**************************************************************)
protocolAdd["RTFmin = " <> ToString[RTFmin]];
protocolAdd["RTFmax = " <> ToString[RTFmax]];
protocolAdd["deltaRTF = " <> ToString[deltaRTF]];
protocolAdd["RTFRepeats = " <> ToString[RTFRepeats]];
protocolAdd["margin = " <> ToString[margin[RTF]]];

protocolBar[];

protocolAdd["Results: "];

protocolAdd[ "RTF" <> " " <> "Chern_number_Model." ];

RTFTab = Table[r,{r, RTFmin, RTFmax, deltaRTF}]; (* this table will be probed*)
SetSharedVariable[RTFReport];
RTFReport = {};

BZ = latticeProbingPointsBZ[npts, a, q];



Map[Module[{},
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


ckModelBZ =
    Function[rtf,
      ParallelMap[
        findCkModel[#, J, J1, lat, a, rec, rtf, support, nIterations,
          nRepeats, nHIO, gamma, pos,
          neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}, DistributedContexts->All]
    ][#];

FxyTModel = FxyT[ ckModelBZ ];
wModel = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTModel];
AppendTo[RTFReport, {#, Re@wModel}];
protocolAdd[ ToString[#] <> " " <> ToString[wModel]  ];
#]
    &,
  RTFTab
];


Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtRTF.dat", RTFReport];


protocolBar[];

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
