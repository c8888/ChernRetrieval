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

latticeProbingPoints::usage =
    "latticeProbingPoints[xmin, xmax, ymin, ymax, \[Delta]x, \[Delta]y] generates a list of points {{x_i, y_i}} at which one measures physical variables of BEC in a trap"

rectLatticeSites::usage =
    "rectLatticeSites[ax, ay, RTF, latticeProbingPoints] generates a list: {{x_i,y_i},type_i} are the rectangular lattice sites, where type_i is the i-th lattice site type (-1 - type A or 1 - type B)
    "
rectLatticeSitesPos::usage =
    "returns something"
rectLatticeSitesNeighbourhood::usage =
    "returns something"
rectLatticeBrillouinZonePoints::usage =
    "rectLatticeBrillouinZonePoints[q, ax, ay, \[Delta]x, \[Delta]y] generates a (flattened) 1D list of the points {x_i,y_i} in the first Brillouin zone for the rectangular lattice, with spacing \[Delta]x, \[Delta]y"

(*)honeycombLatticeSites[a, RTF, latticeProbingPoints]::usage =
    "rectLatticeSites[a RTF, latticeProbingPoints] generates two lists {l1,l2}:
    l1={{x_i,y_i},type_i} are the honeycomb lattice sites, where type_i is the i-th lattice site type (0 - type A or 1 - type B)
    l2={{x_i_j,y_i_j}} is the list of 1..j..k points from latticeProbingPoints which are in the distance of rangeNeighbour from i_th lattice node corresponding to l1
    "
honeycombLatticeBrillouinZonePoints[a, \[Delta]x, \[Delta]y]::usage =
    "honeycombLatticeBrillouinZonePoints[a, \[Delta]x, \[Delta]y] generates a (flattened) 1D list of the points {x_i,y_i} in the first Brillouin zone for the honeycomb lattice, with spacing \[Delta]x, \[Delta]y"

mirror2DSpace[latticeProbingPointsValue]::usage =
    "mirrorSpace[latticeProbingPointsValue] mirrors the 2Dtable with respect to x and y axes"
*)

Begin["`Private`"]



latticeProbingPoints[xmin_, xmax_, ymin_, ymax_, \[Delta]x_, \[Delta]y_]:=
    Table[{x,y}, {x,xmin,xmax,\[Delta]x}, {y,ymin,ymax,\[Delta]y}]

rectLatticeSites[a_, RTF_, xmin_, xmax_, ymin_, ymax_]:=
    Module[ {
      sitesij = {}
    },
      For[i = -Floor[RTF/a], i <= Floor[RTF/a], i++,
        For[j = -Floor[RTF/a], j <= Floor[RTF/a], j++,
          If[
            Norm@N[{i*a, j*a}] <= RTF && i*a >= xmin &&  i*a <= xmax && j*a >= ymin && j*a <=ymax,
              AppendTo[sitesij, {{N[i*a], N[j*a]}, (-1)^(i+j)}]
          ]
          ]
      ];
      Return[N@sitesij]
      ]


rectLatticeSitesPos[latticeProbingPoints_, a_, \[Delta]x_, \[Delta]y_] := Module[{
  ret = {},
  istart = 0,
  iNeigh = {},
  jNeigh = {}
},
  For[i=1, i<=Length[latticeProbingPoints], i++, (*select a row closest to one of the nodes *)
    If[Abs[Abs@latticeProbingPoints[[i,1,1]] - Abs@Round[latticeProbingPoints[[i,1,1]]/ a]] <= \[Delta]x/2,
    istart=i; Break[];
    ]
  ];
  (*find nearest points of the nodes in the row*)
  For[i=istart, i<=Length[latticeProbingPoints], i++,
    If[Abs[Abs@latticeProbingPoints[[i,1,1]] - Abs@Round[latticeProbingPoints[[i,1,1]]/ a]] <= \[Delta]x/2,
      AppendTo[iNeigh, i];
    ]
  ];
  (*find neareast points of the nodes in the column*)
  For[j=1, j<=Dimensions[latticeProbingPoints][[2]], j++,
    If[Abs[Abs@latticeProbingPoints[[istart, j, 2]] - Abs@Round[latticeProbingPoints[[istart, j, 2]]/a]] <= \[Delta]y/2,
      AppendTo[jNeigh, j];
    ]
  ];

  (*namely Cartesian product of the two lists*)
  Return[Flatten[Outer[{#1, #2} &, iNeigh, jNeigh], 1]]
]

rectLatticeSitesNeighbourhood[rectLatticeSitesPos_, latticeProbingPoints_, \[Delta]x_, \[Delta]y_, rangeNeighbour_]:= Module[{
  ret={}
},
  Map[
    AppendTo[ret,
      Round@{
        {
          Max[{#[[1]] - Round[rangeNeighbour/\[Delta]x], 1}],
          Max[{#[[2]] - Round[rangeNeighbour/\[Delta]y], 1}]
        },
        {
          Min[{#[[1]] + Round[rangeNeighbour/\[Delta]x], Dimensions[latticeProbingPoints][[1]]}],
          Min[{#[[2]] + Round[rangeNeighbour/\[Delta]y], Dimensions[latticeProbingPoints][[2]]}]
        }
      }
    ]&, rectLatticeSitesPos];
  Return[ret] (**)
]


End[] (* `Private` *)

EndPackage[]