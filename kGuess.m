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
nRepeats = 3;
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
        nRepeats, nHIO, gamma]; (*nie ma nic wspolnego z wavefkGuess,
  jest to obiekt odzyskany w tej samej pamieci*)

  Total@Total@Abs[Abs[Fourier[wavefkGuess]]^2 - Abs[FTwavefAbs]^2]/
      Total[Total[Abs[FTwavefAbs]^2]] (*metryka bledu*)
]

kxmin = -2Pi/a;
kxmax = 2Pi/a;
deltaKx = 0.6;
kymin = -2Pi/a;
kymax = 2Pi/a;
deltaKy = 0.6;

guess2T =
    ParallelTable[{kxGuess, kyGuess, guess2[{kxGuess, kyGuess}]}, {kxGuess, kxmin, kxmax, deltakx}, {kyGuess, kymin, kymax, deltaky},
      DistributedContexts -> {"space`", "wavefunction`", "HIOER`"}];
Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "kGuess.pdf",
  ListPlot3D[Flatten[guess2T,1],
  AxesLabel -> {"kx", "ky", "1-error" } ]
]

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
