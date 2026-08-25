# Sprint: does the four-subspace (tame->wild) cross-ratio = Paper 59's elliptic lambda? (2026-08-25)

Canonical memo. Following the open lead flagged in the Berry-curvature capture (memo
`sprint_beh2_three_body_composition_memo.md` section 9): the 4-center (D~4 four-subspace) problem
carries a continuous cross-ratio parameter, and Paper 59's elliptic frontier is the
Legendre/Gamma(2) lambda-line -- are they the SAME lambda? Driver: sympy in-session +
`tests/test_paper59_four_subspace_lambda.py`.

## Current-state check (mandatory, done first)
- Paper 59 sec:modular: lambda(tau(rho)) = 1-rho exact (Legendre/X(2)); rho = c2/c1 = (a1/a2)^2 =
  ratio of the TWO density masses of the one 3-center ERI T2=(XY|XZ). Four branch points of the
  quartic model {+-1, +-i sqrt((1-rho)/rho)}.
- Memory [[cosmic-galois-elliptic-rung1]]: "landing on the universal Legendre family is EXPECTED
  for any 4-branch-point family; don't oversell." [[geovac-axis-map]]: category error -- genus /
  cross-ratio "counts densities, not nuclei."

## What was computed (exact, sympy)
1. The quartic branch points {+-1,+-ib} give COMPLEX cross-ratios (2 real + 2 imaginary pts). The
   actual elliptic curve whose periods are Paper 59's K(rho),K(1-rho) is the 2:1 reduction
   v^2 = u(u-1)(rho u + 1-rho) (u = x^2; dx/sqrt(Q) = (1/2) du/sqrt(cubic)), branch points
   **{0, 1, oo, (rho-1)/rho}**.
2. Their six-element cross-ratio orbit = **{(rho-1)/rho, 1/rho, rho/(rho-1), rho, 1-rho, -1/(rho-1)}**
   -- contains **both** period arguments rho AND 1-rho. Paper 59's lambda = 1-rho is the member
   L/(L-1). So **Paper 59's lambda IS a four-point cross-ratio** (verified exactly).

## Verdict (half-affirm, half-deflate)
- **YES (verified):** Paper 59's elliptic lambda = 1-rho is genuinely a four-point cross-ratio =
  the tame four-subspace (D~4) modulus = the cross-ratio of the four regular singular points of a
  rank-2 Fuchsian system (isomonodromy = Painleve VI; R. Fuchs elliptic representation, grounded
  via web search: PVI = rank-2 Fuchsian, 4 reg. sing. pts, cross-ratio modulus). The math bridge
  is real, not a pun.
- **BUT (the axis-map disambiguation):** the four points are the branch data of the TWO density
  masses (rho = c2/c1), NOT four nuclei. The lead's "4-CENTER" framing MISLOCATED the objects.
  The four-subspace structure lives on the DENSITY/period axis (2 masses -> 4 branch points),
  not the nuclei/geometry axis. Consistent with [[geovac-axis-map]] "counts densities not nuclei."
- **The Berry-side four-NUCLEI cross-ratio is a different object:** a genuine 4-atom configuration
  operator's cross-ratio (if definable) is a function of nuclear positions, independent of rho --
  NOT the same value. The two lambdas share only the abstract X(2) home, which ANY four-point set
  has (the memory's "expected, don't oversell").

## Captured
- Paper 59 sec:modular new paragraph "The modulus is a four-point cross-ratio" (verified identity
  + four-subspace/Fuchsian/PVI lineage + the density-not-nuclei caveat). [SYMBOLIC].
- Backing `tests/test_paper59_four_subspace_lambda.py` (2 fast: cross-ratio orbit contains
  {rho,1-rho}; branch points + period reduction). Paper 59 compiles (exit 0, 15 pp).
- Memory [[cosmic-galois-elliptic-rung1]] + [[geovac-axis-map]] updated. CHANGELOG + matrix.
- Compounds the Paper 59 Phase-4 re-review OWED (branch-local; merge/Release PI-only).

## Honest scope
No new number cracked; this is an interpretive/structural placement of an existing exact result
(lambda=1-rho) + a disambiguation of a tempting cross-corpus conflation. The generic-Hain-Brown /
integrated-MMV frontier (Rung 3) is untouched. Whether a genuine 4-nuclei config operator carries
an elliptic modulus is a separate (geometry-axis) question, not pursued.
