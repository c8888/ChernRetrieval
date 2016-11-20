(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: space *)
(* :Context: space` *)
(* :Author: c8888 *)
(* :Date: 2016-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["space`"]
(* Exported symbols added here with SymbolName::usage *)

latticeProbingPoints[xmin, xmax, ymin, ymax, \[Delta]x, \[Delta]y]::usage =
    "latticeProbingPoints[xmin, xmax, ymin, ymax, \[Delta]x, \[Delta]y] generates a 2D array of points {{x_1...x_n}_1..._m} at which one measures physical variables of BEC in trap"

rectLatticeSites[ax, ay, RTF, latticeProbingPoints, rangeNeighbour]::usage =
    "rectLatticeSites[ax, ay, RTF, latticeProbingPoints] generates two lists {l1,l2}:
    l1={{x_i,y_i},type_i} are the rectangular lattice sites, where type_i is the i-th lattice site type (0 - type A or 1 - type B)
    l2={{x_i_j,y_i_j}} is the list of 1..j..k points from latticeProbingPoints which are in the distance of rangeNeighbour from i_th lattice node corresponding to l1
    "
rectLatticeBrillouinZonePoints[q, ax, ay, \[Delta]x, \[Delta]y]::usage =
    "rectLatticeBrillouinZonePoints[q, ax, ay, \[Delta]x, \[Delta]y] generates a (flattened) 1D list of the points {x_i,y_i} in the first Brillouin zone for the rectangular lattice, with spacing \[Delta]x, \[Delta]y"

honeycombLatticeSites[a, RTF, latticeProbingPoints]::usage =
    "rectLatticeSites[a RTF, latticeProbingPoints] generates two lists {l1,l2}:
    l1={{x_i,y_i},type_i} are the honeycomb lattice sites, where type_i is the i-th lattice site type (0 - type A or 1 - type B)
    l2={{x_i_j,y_i_j}} is the list of 1..j..k points from latticeProbingPoints which are in the distance of rangeNeighbour from i_th lattice node corresponding to l1
    "
honeycombLatticeBrillouinZonePoints[a, \[Delta]x, \[Delta]y]::usage =
    "honeycombLatticeBrillouinZonePoints[a, \[Delta]x, \[Delta]y] generates a (flattened) 1D list of the points {x_i,y_i} in the first Brillouin zone for the honeycomb lattice, with spacing \[Delta]x, \[Delta]y"

mirrorSpace[latticeProbingPointsValue]::usage =
    "mirrorSpace[latticeProbingPointsValue] mirrors the table with respect to x and y axes"


Begin["`Private`"]

End[] (* `Private` *)

EndPackage[]