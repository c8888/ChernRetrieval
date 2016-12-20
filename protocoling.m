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

protocolAdd::usage =
    "protocolAdd[stringMessage] adds stringMessage at the end of protocol file"
protocolMaxMemoryUsed::usage =
    "protocolMaxMemoryUsed[] adds information about max memory used (in GB) in current kernel session"

Begin["`Private`"]

protocolAdd[stringMessage_]:=PutAppend[stringMessage, "out/" <> ToString[Last@$CommandLine] <> "_"  <> ToString[$ProcessID] <> "protocol.txt"]
protocolMaxMemoryUsed[]:=protocolAdd["Max memory used (GB): " <> ToString[MaxMemoryUsed[]/1024.^3]]

End[] (* `Private` *)

EndPackage[]