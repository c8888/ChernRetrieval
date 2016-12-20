
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
nIterations = 1000;
nRepeats = 5;
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

iterStep = 75; (* every iterStep iterations the overlap and chern number are computed *)
protocolAdd[ "Repeat" <> " " <> "Iteration" <> " " <> "Chern number retr." <> " " <> "Mean overlap" <> " " <> " Standard deviation of overlap"];
overlapIter = {};

(*modded phase retrieve to go into the process *)
phaseRetrieveSupportOverlap[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_, iterStep_]:= (* returns table of a retrieved object *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={},
      FTxi2={}
    },
      {nCol, nRow} = Dimensions[FTXAbs];
      xi=Table[RandomComplex[], nCol, nRow]; (* random initialization, different complex numbers at each repetition *)


      For[k = 0, k < nRepeats, k++,
        For[i = 0, i < nIterations, i++,
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi2 = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi2];
          If[Unequal[Mod[i,nHIO],0],
          (* HIO case *) xi=Table[If[support[[q,w]] == 1,
          (* inside support*) xi[[q,w]],
          (* outside support *) xiprim[[q,w]]-gamma*xi[[q,w]] ],{q,nCol},{w,nRow}];  ,
          (* ER case *)
            xi*=support;]
          If[ Mod[i,iterStep] == 0,
            ckRetrBZ =
                ParallelMap[
                  findCkRetr[#, J, J1, lat, a, rec, RTF, support, nIterations,
                    nRepeats, nHIO, gamma, pos,
                    neighpos, \[Sigma]w, \[Beta], \[Delta]x, \[Delta]y] &, BZ, {2}];
            FxyTRetr = FxyT[ ckRetrBZ[[All, All, 1]] ];
            wRetr = 1/(2 \[Pi] I )*Chop@Total@Total[FxyTRetr];
            AppendTo[overlapIter, {k, i, Re@wRetr, Mean@Flatten@ckRetrBZ[[All, All, 2]], StandardDeviation@Flatten@ckRetrBZ[[All, All, 2]]}];

            protocolAdd[ ToString[k] <> " " <> ToString[i] <> " " <> ToString[Re@wRetr] <> " " <> ToString[Mean@Flatten@ckRetrBZ[[All, All, 2]] ]
                <> " " <> ToString[StandardDeviation@Flatten@ckRetrBZ[[All, All, 2]] ] ];

          ];
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 0, retrerror=xierror;];
        If[xierror<=retrerror, retr=xi; retrerror=xierror;, Null;, retr=xi]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
      ];
      Return[retr*support]; (* returned value *)

    ]

phaseRetrieveSupportOverlap[FTwavefAbs, support, nIterations, nRepeats, nHIO, gamma, iterStep];

Export["out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "chernNumberAtOverlap.dat", overlapIter];



t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];
