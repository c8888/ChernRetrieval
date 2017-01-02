
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
xmin = -8;
xmax = 8;
ymin = -8;
ymax = 8;
RTF = 6;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
nIterations = 500;
nRepeats = 4;
nHIO = 20;
gamma = 0.9;
npts = 5;(*points in the 1st Brillouin zone*)
(**************************************************************)
\[Sigma]PhNoiseMin = 0.000001; (* every iterStep iterations the overlap and chern number are computed *)
\[Sigma]PhNoiseMax = 16;
\[Sigma]PhNoiseMult = 1.1;
(***************z***********************************************)

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
(*************************************************************)
(*************************************************************)

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
LaunchKernels[16];

wavefBZ = ParallelMap[waveFunctionHarper[lat, a, J, J1, rec, RTF,
  #, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y]&, BZ, {2}, DistributedContexts->All];
\[Sigma]PhNoiseTab = Table[\[Sigma]PhNoiseMin*\[Sigma]PhNoiseMult^i, {i, 0, Round@Log[\[Sigma]PhNoiseMax/\[Sigma]PhNoiseMin]/Log[\[Sigma]PhNoiseMult]}]

CkModelNoisedBZ[\[Sigma]PhNoise_] :=
      Map[wannierBaseRectProject[Abs@# * Map[Exp[I RandomVariate[NormalDistribution[Arg[#], \[Sigma]PhNoise]]]&, #, {2} ], lat, rec, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF]&, wavefBZ, {2}]

cnAtPhaseNoiseTab = Transpose[{\[Sigma]PhNoiseTab, Re@ParallelMap[1/(2 \[Pi] I )*Chop@Total@Total[FxyT[ CkModelNoisedBZ[#] ]]&, \[Sigma]PhNoiseTab, DistributedContexts->All]}];

protocolAdd[ "\[Sigma]PhaseNoise" <> " " <> "Chern number" ];


Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "cnAtPhaseNoise.dat", cnAtPhaseNoiseTab];



t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
