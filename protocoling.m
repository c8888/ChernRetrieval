(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: protocoling *)
(* :Context: protocoling` *)
(* :Author: c8888 *)
(* :Date: 2016-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: This package includes functions used in order to write a human-readable protocol of a simulation while in progress. *)

BeginPackage["protocoling`"]
(* Exported symbols added here with SymbolName::usage *)

protocolSet[folder, processID]::usage =
    "protocolSetName[folder, processID] creates a new protocol file folder/processID.dat and returns the protocol name object"
protocolAdd[stringMessage]::usage =
    "protocolAdd[stringMessage] adds stringMessage at the end of protocol file"
protocolMaxMemoryUsed[]::usage =
    "protocolMaxMemoryUsed[] adds information about max memory used (in GB) in current kernel session"

Begin["`Private`"]
PutAppend[{k,i}, "out/" <> ToString[$ProcessID] <> "iterations.dat"];
End[] (* `Private` *)

EndPackage[]