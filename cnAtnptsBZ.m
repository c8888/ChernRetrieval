
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
q = 5; (* Pi-flux *)
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
RTF = 5;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
(*nIterations = 100;
nRepeats = 3;
nHIO = 20;
gamma = 0.9;*)
(*npts = 5;(*points in the 1st Brillouin zone*)*)
n = 3;
(**************************************************************)
nptsmin = 2; (* integer values *)
nptsmax = 15;
nptsRepeat = 1;

(**************************************************************)
protocolAdd[" "];
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
protocolAdd["n = "<> ToString[n]];
protocolBar[];

(***************************************************)
protocolAdd["nptsmin = "<> ToString[nptsmin] ];
protocolAdd["nptsmax = "<> ToString[nptsmax] ];
protocolAdd["nptsRepeat = "<> ToString[nptsRepeat]];
protocolBar[];
(***************************************************)

protocolAdd["Results: "];

protocolAdd[ "npts" <> " " <> "Chern_number_retr." <> " " <> "Mean overlap" <> " " <> " Standard deviation of overlap"];

nptsTab = Table[n,{n, nptsmin, nptsmax}]; (* this table will be probed*)
SetSharedVariable[nptsReport];
nptsReport = {};

lat = latticeProbingPoints[xmin, xmax, ymin,
  ymax, \[Delta]x, \[Delta]y];
rec = rectLatticeSites[a, RTF, xmin, xmax, ymin, ymax];
pos = rectLatticeSitesPos[lat, a, \[Delta]x, \[Delta]y, RTF];
neighpos =
    rectLatticeSitesNeighbourhood[pos, lat, \[Delta]x, \[Delta]y,
      rangeNeighbour];
support = Map[If[Norm[#] < RTF, 1, 0] &, lat, {2}];

Map[Module[{
  BZ = latticeProbingPointsBZ[#, a, q]
},

  (**************************************************************)
  \[Beta] =
      wannierNormalisationFactor[\[Sigma]w, \[Delta]x, \[Delta]y, lat];

  (**************************************************************)

  ckModelBZ =
      ParallelMap[
        findCkModelQ[#, J, J1, lat, a, rec, RTF, support, nIterations,
          nRepeats, nHIO, gamma, pos,
          neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, q, n] &, BZ, {2}, DistributedContexts->All];
  FxyTModel = FxyT[ ckModelBZ ];
  wModel = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTModel];

  AppendTo[nptsReport, {#, Re@wModel}];
  protocolAdd[ ToString[#] <> " " <> ToString[wModel]  ];
  #]
    &,
  nptsTab
];


Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtnptsQ.dat", nptsReport];
(*Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtnptsPlot.pdf",
  ListPlot[nptsReport[[All,1;;2]]]];*)

protocolBar[];

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
