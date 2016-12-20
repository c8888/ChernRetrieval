(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: c8888 *)
(* :Date: 2016-12-19 *)

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
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
RTF = 3;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 2;
nIterations = 100;
nRepeats = 2;
nHIO = 20;
gamma = 0.9;
npts = 5;(*points in the 1st Brillouin zone*)
(**************************************************************)
kxmin = 0;
kxmax = 2.Pi/a;
deltaKx = 0.3;
kymin = 0;
kymax = 2.Pi/a;
deltaKy = 0.3;
(*************************************************************)
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
(***************************************************************)
protocolAdd["kxmin = "<> ToString[kxmin] ];
protocolAdd["kxmax = "<> ToString[kxmax] ];
protocolAdd["deltaKx = "<> ToString[deltaKx] ];

protocolAdd["kymin = "<> ToString[kymin] ];
protocolAdd["kymax = "<> ToString[kymax] ];
protocolAdd["deltaKy = "<> ToString[deltaKy] ];

(***************************************************************)

protocolAdd["Results: "];


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
wavef = waveFunctionHarper[lat, a, J, J1, rec, RTF,
  k0, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y];
FTwavefAbs = Abs@Fourier@wavef;

guess2[kGuess_] := Module[{
  wavefkGuess =
      waveFunctionHarper[lat, a, J, J1, rec, RTF,
        kGuess, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y]
},
  wavefkGuess =
      phaseRetrieveGuess[FTwavefAbs, Abs@wavefkGuess, support, nIterations,
        nRepeats, nHIO, gamma]; (*use the same memory*)

  Function[error,
    protocolAdd[ToString[{kGuess[[1]], kGuess[[2]], error}]];
    Return[error]
  ] [Total@Total@Abs[Abs[Fourier[wavefkGuess]]^2 - Abs[FTwavefAbs]^2]/
      Total[Total[Abs[FTwavefAbs]^2]] (*error metrics*)]
]



guess2T =
    ParallelTable[{kxGuess, kyGuess, guess2[{kxGuess, kyGuess}]}, {kxGuess, kxmin, kxmax, deltaKx}, {kyGuess, kymin, kymax, deltaKy},
      DistributedContexts -> {"space`", "wavefunction`", "HIOER`"}];
Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "kGuess.dat", Flatten[guess2T,1]];
(*Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "kGuessPlot.pdf",
  ListPlot3D[Flatten[guess2T,1], AxesLabel -> {"kx", "ky", "error"}, ColorFunction->Hue]
];*)

protocolBar[];

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
