(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: HIOER *)
(* :Context: HIOER` *)
(* :Author: c8888 *)
(* :Date: 2016-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: This package includes the methods for Hybrid Input-Output Error Reduction Phase Retrieval Algorithm.
                Wherever possible, numerical operations are parallelized.*)

BeginPackage["HIOER`"]
(* Exported symbols added here with SymbolName::usage *)
Needs["protocoling`"];
Needs["space`"];

phaseRetrieveSupport::usage =
    "phaseRetrieveSupport[FTXAbs, support, nIterations, nHIO] returns a table of retrieved complex X from its (representation in reciprocal space) absolute value FTXAbs.
    Support is the table of the same structure as FTXAbs, it has rough values {0,1} for each {kx,ky} point. Procedure performs nIterations HIOER iterations. ER is made every nHIO
iterations not to stay in local minima."
phaseRetrieveGuess::usage =
    "phaseRetrieveGuess[FTXAbs, XAbsGuess, nIterations, nHIO] returns a table of retrieved complex X from its reciprocal space representation absolute value FTXAbs.
     Uses XAbsGuess to converge faster.ER is made every nHIO iterations not to stay in local minima."

Begin["`Private`"]


phaseRetrieveGuess[FTXAbs_, wavefAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_]:= (* returns table of a retrieved object *)
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
          (* inside support*) wavefAbs[[q,w]]*Exp[I*Arg[xi[[q,w]]]],
          (* outside support *) xiprim[[q,w]]-gamma*xi[[q,w]] ],{q,nCol},{w,nRow}];  ,
          (* ER case *)
            xi=wavefAbs*Exp[I*Arg[xi]]];
          ,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 0, retrerror=xierror];
        If[xierror<=retrerror, retr=xi; retrerror=xierror, Null, retr=xi]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
      ,
        nRepeats
      ];
      Return[retr*support]; (* returned value *)

    ]
phaseRetrieveSupport[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_]:= (* returns table of a retrieved object *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={}
    },
      {nCol, nRow} = Dimensions[FTXAbs];

      Do[
        xi=Table[RandomComplex[], nCol, nRow]; (* random initialization, different complex numbers at each repetition *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = mirror2DSpace@InverseFourier[FTxi];
          If[Unequal[Mod[i,nHIO],0],
          (* HIO case *) xi=Table[If[support[[q,w]] == 1,
          (* inside support*) xi[[q,w]],
          (* outside support *) xiprim[[q,w]]-gamma*xi[[q,w]] ],{q,nCol},{w,nRow}]  ,
          (* ER case *)
            xi*=support];
          ,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 0, retrerror=xierror;];
        If[xierror<=retrerror, retr=xi; retrerror=xierror, Null, retr=xi]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
      ,
        nRepeats
      ];
      Return[retr*support]; (* returned value *)

    ]

End[] (* `Private` *)

EndPackage[]