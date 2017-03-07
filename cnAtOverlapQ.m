
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
q = 3; (* Pi-flux *)
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
RTF = 3.6;
rangeNeighbour = 0.6;
a = 1;
\[Sigma]w = 0.2;
k0 = {1, 2}; (* there is need to guess it from experimental data. One can use only the support too*)
J = 1;
J1 = 3;
nIterations = 2000;
nRepeats = 6;
nHIO = 20;
gamma = 0.9;
npts = 10;(*points in the 1st Brillouin zone*)
n = 3; (* band number 1...q *)
iterStepExport = 50;
(**************************************************************)
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
protocolAdd["n = "<> ToString[n] ];
protocolBar[];

(**************************************************************)

lat = latticeProbingPoints[xmin, xmax, ymin,
  ymax, \[Delta]x, \[Delta]y];
rec = rectLatticeSitesQ[a, RTF, xmin, xmax, ymin, ymax, q];
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
SetSharedVariable[ckRetrOverlap];
ckRetrOverlap = {};

phaseRetrieveSupportDebug[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_,
  iterStepExport_, wavef_, lat_, rec_, pos_,neighpos_,\[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_, RTF_, ckModel_, k_, BZ_]:= (* returns table of a retrieved object *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={},
      ckRetr = {},
      ckRetrMirror = {}
    },
      {nCol, nRow} = Dimensions[FTXAbs];

      Do[
        xi=Table[RandomComplex[], nCol, nRow]; (* random initialization, different complex numbers at each repetition *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi];
          If[Unequal[Mod[i,nHIO],0],
          (* HIO case *) xi=Table[If[support[[q,w]] == 1,
          (* inside support*) xi[[q,w]],
          (* outside support *) xiprim[[q,w]]-gamma*xi[[q,w]] ],{q,nCol},{w,nRow}]  ,
          (* ER case *)
            xi*=support];



          If[Mod[i, iterStepExport]==0,
            If[Position[ckRetrOverlap, {i,j}]=={}, AppendTo[ckRetrOverlap, {{i, j}, BZ}]];
            ckRetr = wannierBaseRectProject[xi, lat, rec, pos,
              neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF];

            ckRetrMirror = wannierBaseRectProject[mirrorXY@xi, lat, rec, pos,
              neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF];

            overlapRetr = overlapWannier[ckModel, ckRetr];
            overlapRetrMirror = overlapWannier[ckModel, ckRetrMirror];
            If[overlapRetr >= overlapRetrMirror,
              ckRetrOverlap[[First@Position[ckRetrOverlap, n_ /; n == {i, j}][[1]]]][[2]][[First@
                  First@Position[BZ, n_ /; n == k], (First@
                  Position[BZ, n_ /; n == k])[[2]]]] = {ckRetr, overlapRetr}, (*not reflected*)
            ckRetrOverlap[[First@Position[ckRetrOverlap, n_ /; n == {i, j}][[1]]]][[2]][[First@
                First@Position[BZ, n_ /; n == k], (First@
                Position[BZ, n_ /; n == k])[[2]]]] = {ckRetrMirror, overlapRetrMirror}(*1 - reflected*)
            ]
          ];
          ,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 0, retrerror=xierror;];
        If[xierror<=retrerror, retr=xi; retrerror=xierror, Null, retr=xi]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
        ,
        {j,nRepeats}
      ];
      Return[retr*support]; (* returned value *)

    ];

findCkRetrSupportDebugQ[k_, J_, J1_, lat_, a_, rec_, RTF_, support_,nIterations_,nRepeats_,nHIO_,gamma_, pos_, neighpos_,\[Sigma]w_, \[Beta]_, \[Delta]x_, \[Delta]y_, q_, n_, iterStepExport_, BZ_] := Module[{
  wavef =
      waveFunctionHarperQ[lat, a, J, J1, rec, RTF,
        k, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, q, n],
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
  wavef = phaseRetrieveSupportDebug[FTwavefAbs, support,
    nIterations, nRepeats, nHIO, gamma, iterStepExport, wavef, lat, rec, pos,
    neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, RTF, ckModel, k, BZ]; (*not to waste the memory*)

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
];

ckRetrSupportBZ =
    Map[
      findCkRetrSupportDebugQ[#, J, J1, lat, a, rec, RTF, support, nIterations,
        nRepeats, nHIO, gamma, pos,
        neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y, q, n, iterStepExport, BZ] &, BZ, {2}(*, DistributedContexts->All*)];
FxyTRetrSupport = FxyT[ ckRetrSupportBZ[[All, All, 1]] ];
wRetrSupport = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTRetrSupport];
protocolAdd[ "wRetrSupport = " <> ToString[wRetrSupport] <> " with mean overlap: " <> ToString[Mean@Flatten@ckRetrSupportBZ[[All, All, 2]] ] <> " stDev: " <> ToString[StandardDeviation@Flatten@ckRetrSupportBZ[[All, All, 2]] ] ];
Export["out/" <> ToString[$ProcessID] <>"ckRetrOverlap.mx", ckRetrOverlap];
(*Export["out/" <> ToString[$ProcessID] <>"ckRetrSupportBZ.dat", ckRetrSupportBZ];
Export["out/" <> ToString[$ProcessID] <>"FxyTRetrSupport.dat", FxyTRetrSupport];
(*Export["out/" <> ToString[$ProcessID] <> "retr_overlap.pdf", ListPlot@Flatten@ckRetrBZ[[All,All,2]]];
(**************************************************************)
*)*)





t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
