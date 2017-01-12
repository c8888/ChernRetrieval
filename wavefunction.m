(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: wavefunction *)
(* :Context: wavefunction` *)
(* :Author: c8888 *)
(* :Date: 2016-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["wavefunction`"]
Needs["space`"]
(* Exported symbols added here with SymbolName::usage *)
waveFunctionHarperQ::usage =
    "waveFunctionHarperQ[latticeProbingPoints_, a_, J_, J1_, rectLatticeSites_, RTF_, k0_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_, q_, n_] returns a value of Harper wave function of n-th energy band in momentum space at {kx,ky} in the
    Harper Pi/q-flux model."
waveFunctionHarper::usage =
    "waveFunctionHarper[{x, y}, rectLatticeSites, RTF, \[Sigma]w, \[Phi]k, \[Theta]k] returns a value of wave function in real space at {x,y} in the
    Harper Pi-flux model (q=2) for initial state \[Phi]k, \[Theta]k. Uses knowledge of neighbor positions to spped up calculations."
hamiltonianHarperK::usage =
    "hamiltonianHarperK[{kx, ky}, J, J1] returns hamiltonian matrix in momentum space at {x,y} in the
    Harper Pi-flux model (q=2). J and J1 are the tunneling factors."
hamiltonianHarperQK::usage =
    "hamiltonianHarperQK[{kx, ky}, J, J1, q] returns hamiltonian matrix in momentum space at {x,y} in the
    Harper 2Pi/q-flux model . J and J1 are the tunneling factors in x and y directions."
wannier::usage =
    "wannier[{x,y}, {xi, yi}, \[Sigma]w, \[Beta]] returns value of wannier function localized at {xj, yj}. \[Sigma]w is the Gauss peak width and \[Beta] is the normalisation factor."
wannierNormalisationFactor::usage =
    "wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, latticeProbingPoints_] returns a normalisation factor for Wannier function localized at {0,0}"
blochSphereAnglesHarper::usage =
    "blochSphereAnglesHarper[{kx,ky}] calculates the angles {\[Theta]k, \[Phi]k} in Bloch sphere representation from a Harper Hamiltonian"
wannierBaseRectProject::usage =
"wannierBaseRectProject[waveFunction_, latticeProbingPoints_, rectLatticeSites_, rectLatticeSitesPos_, rectLatticeSitesNeighbourhood_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_]
returns a wave function in a basis of orthonormal wannier functions"
ComplexDotProduct::usage =
    "ComplexDotProduct[x_, y_] := Chop[Dot[x, Conjugate[y]]]"

wannierBaseABproject::usage =
    "projects the wave function only on the two nodes A(0,0), B(0,1) in the centre"

overlapWannier::usage =
    "calculates overlap magnitude square when given wannier basis representation of wave functions"

myES::usage =
    "myES[hamiltonian] gives a sorted list of eigenvalues and eigenvectors for a hamiltonian matrix"

Begin["`Private`"]

wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, latticeProbingPoints_] := 1/Sqrt[
        Total@Total[
          Chop@Map[wannier[#, {0, 0}, \[Sigma]w, 1] &, latticeProbingPoints, {2}]^2
        ] * \[Delta]x * \[Delta]y
    ]

wannier[r_,ri_,\[Sigma]w_, wannierNormalisationFactor_]:= wannierNormalisationFactor*Exp[-Norm[N[r]-N[ri]]^2/2/\[Sigma]w^2]

hamiltonianHarperK[k_, J_, J1_] := Module[
      {
        h = {-J - J1 Cos[2 k[[1]]], J1 Sin[2 k[[1]]], -2 J Cos[k[[2]]] }
      },
      h[[1]]*PauliMatrix[1] + h[[2]]*PauliMatrix[2] + h[[3]]*PauliMatrix[3]
    ] (*no h0*Id factor, which means that we assume on-site energies equal to zero *)

Ay[m_, q_] := 2 Pi m / q

hamiltonianHarperQK[k_, J_, J1_, q_] := SparseArray[{
  {i_, i_} -> -2 J1 Cos[k[[2]] - Ay[i, q]], (*diag*)
  {i_, j_} /;
      i - j == 1 -> -J Exp[-I k[[1]]], (*pozadiag*)
  {i_, j_} /;
      i - j == -1 -> -J Exp[I k[[1]]]},(*pozadiag*)
  {q, q}] +
    SparseArray[{
      {i_, j_} /; {i, j} == {1, q} -> -J Exp[-I k[[1]]],(*róg p-
    górny*)
      {i_, j_} /; {i, j} == {q, 1} -> -J Exp[
        I k[[1]]]  (*róg l-dolny*)
    },
      {q, q}]

myES[hamiltonian_] := Module[{
  es = Transpose@
        Sort[Transpose[Eigensystem[hamiltonian]], #1[[1]] < #2[[1]] &]},
  es[[2]] = Map[Exp[-I Arg[#[[1]]]] * # &, es[[2]]];
  Return[es ] (*es[[2, i, 1]] is always real*)
]


blochSphereAnglesHarper[k_, J_, J1_] := Module[{
  h = {-J - J1 Cos[2 k[[1]] ], J1 Sin[2 k[[1]] ], -2 J Cos[ k[[2]] ] }
},
  N@CoordinateTransformData[
    "Cartesian" -> "Spherical", "Mapping", h][[2;;3]]
]

waveFunctionHarperQ[latticeProbingPoints_, a_, J_, J1_, rectLatticeSites_, RTF_, k0_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_, q_, n_] :=
    Module[{
      ret,
      um = myES[hamiltonianHarperQK[k0, J, J1, q]]
    },
      ret = Chop[
        Map[
          Function[r,
            Total[
              Map[
                If[
                  Abs[Norm[r] - Norm[#[[1]]]] < 5 \[Sigma]w && Norm[RTF]^2 > Norm[r]^2,
                    Sqrt[Norm[RTF]^2 - Norm[r]^2] * um[[2, n, #[[2]] ]] * Exp[I k0.#[[1]] ] * wannier[r, #[[1]], \[Sigma]w, wannierNormalisationFactor]
                    ,
                  0 (*When too far from a node*)
                ] &,
                rectLatticeSites
              ]
            ]
          ][#]&,
          latticeProbingPoints,{2}
        ]
      ];
      ret/Sqrt[Total[Total[Abs[ret]^2]] * \[Delta]x * \[Delta]y]
    ]


waveFunctionHarper[latticeProbingPoints_, a_, J_, J1_, rectLatticeSites_, RTF_, k0_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_] :=
    Module[{
  \[Theta]k = blochSphereAnglesHarper[k0, J, J1][[1]],
  \[Phi]k = blochSphereAnglesHarper[k0, J, J1][[2]],
      ret
    },
      ret = Chop[
      Map[
        Function[r,
          Total[
            Map[
              If[
                Abs[Norm[r] - Norm[#[[1]]]] < 5 \[Sigma]w && Norm[RTF]^2 > Norm[r]^2,
            Which[
                  #[[2]] == 1 , Sqrt[Norm[RTF]^2 - Norm[r]^2] Sin[ \[Theta]k / 2 ] * Exp[I k0.#[[1]] ] * wannier[r, #[[1]], \[Sigma]w, wannierNormalisationFactor], (* Site A *)
                  #[[2]] == -1 , Sqrt[Norm[RTF]^2 - Norm[r]^2] Cos[ \[Theta]k / 2 ] * Exp[I \[Phi]k] * Exp[I k0.#[[1]] ] * wannier[r, #[[1]], \[Sigma]w, wannierNormalisationFactor] (* Site B *)
                  ],
                0 (*When too far from a node*)
              ] &,
              rectLatticeSites
            ]
          ]
        ][#]&,
        latticeProbingPoints,{2}
      ]
    ];
    ret/Sqrt[Total[Total[Abs[ret]^2]] * \[Delta]x * \[Delta]y]
  ]

wannierBaseRectProject[waveFunction_, latticeProbingPoints_, rectLatticeSites_, rectLatticeSitesPos_, rectLatticeSitesNeighbourhood_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_,RTF_] :=
    Module[{
      ret
    },
      ret = Table[
  Sum[
    Chop[wannier[latticeProbingPoints[[i, j]], rectLatticeSites[[q,1]], \[Sigma]w, wannierNormalisationFactor ] ]*
        waveFunction[[i, j]]/If[Norm[RTF]>Norm[latticeProbingPoints[[i, j]]], Sqrt[Norm[RTF]^2-Norm[latticeProbingPoints[[i, j]]]^2], 1],
    {i, rectLatticeSitesNeighbourhood[[q, 1, 1]],
      rectLatticeSitesNeighbourhood[[q, 2, 1]]}, {j,
    rectLatticeSitesNeighbourhood[[q, 1, 2]], rectLatticeSitesNeighbourhood[[q, 2, 2]]}
  ],
  {q, Length@rectLatticeSitesPos}] * \[Delta]x * \[Delta]y;
      ret*1/Sqrt@Total[Abs[ret]^2]
      ]

ComplexDotProduct[x_, y_] := Chop[Dot[x, Conjugate[y]]]

wannierBaseABproject[waveFunction_, latticeProbingPoints_, rectLatticeSites_, rectLatticeSitesPos_, rectLatticeSitesNeighbourhood_, \[Sigma]w_, wannierNormalisationFactor_, \[Delta]x_, \[Delta]y_,RTF_] :=
    Module[{
    qApos = First@First@Position[rectLatticeSites, {0., 0.}];
      qBpos = First@First@Position[rectLatticeSites, {0., 1.}];
    },
      Map[
  Sum[
    Chop[wannier[latticeProbingPoints[[i, j]], rectLatticeSites[[#,1]], \[Sigma]w, wannierNormalisationFactor ] ]*
        waveFunction[[i, j]]/If[Norm[RTF]>Norm[latticeProbingPoints[[i, j]]], Sqrt[Norm[RTF]^2-Norm[latticeProbingPoints[[i, j]]]^2], 1],
    {i, rectLatticeSitesNeighbourhood[[#, 1, 1]],
      rectLatticeSitesNeighbourhood[[#, 2, 1]]}, {j,
    rectLatticeSitesNeighbourhood[[#, 1, 2]], rectLatticeSitesNeighbourhood[[#, 2, 2]]}
  ]&,
  {qApos, qBpos}] * \[Delta]x * \[Delta]y
]

overlapWannier[ckModel_, ckRetr_]:= Abs[ComplexDotProduct[ckModel, ckRetr]]^2

End[] (* `Private` *)

EndPackage[]