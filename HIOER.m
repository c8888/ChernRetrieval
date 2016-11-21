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

phaseRetrieveSupport::usage =
    "phaseRetrieveSupport[FTXAbs, RTF, nIterations, nHIO] returns a table of retrieved complex X from its (representation in reciprocal space) absolute value FTXAbs.
    Support is defined to be 1 over the circle of radius RTF and 0 elsewhere. Performs nIterations HIOER iterations."
phaseRetrieveGuess::usage =
    "phaseRetrieveGuess[FTXAbs, XAbsGuess, nHIO] returns a table of retrieved complex X from its reciprocal space representation absolute value FTXAbs.
     Uses XAbsGuess to converge faster."

Begin["`Private`"]

End[] (* `Private` *)

EndPackage[]