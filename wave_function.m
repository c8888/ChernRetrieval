(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: wave_function *)
(* :Context: wave_function` *)
(* :Author: c8888 *)
(* :Date: 2016-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["wave_function`"]
(* Exported symbols added here with SymbolName::usage *)
waveFunctionHarperK::usage =
    "waveFunctionHarperK[{kx,ky}, \[Phi]k, \[Theta]k] returns a value of wave function in momentum space at {kx,ky} in the
    Harper Pi-flux model (q=2)."
waveFunctionHarper::usage =
    "waveFunctionHarper[{x, y}, rectLatticeSites, \[Phi]k, \[Theta]k] returns a value of wave function in real space at {x,y} in the
    Harper Pi-flux model (q=2)."
hamiltonianHarperK::usage =
    "hamiltonianHarperK[{kx, ky}, J, J1, \[CapitalDelta] ] returns hamiltonian matrix in momentum space at {x,y} in the
    Harper Pi-flux model (q=2)."

Begin["`Private`"]

End[] (* `Private` *)

EndPackage[]