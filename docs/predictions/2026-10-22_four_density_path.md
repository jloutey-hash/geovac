# Pre-registered prediction: the real three-dimensional path of four densities (FP1-FP3)

**Registered:** 2026-10-22, before any computation in this file (PI-directed: "Let's do those next steps",
item two of the 2026-10-21 overview).
**Prompted by:** Paper 59 `sec:four_density` states that the real three-dimensional integral reduces into
the family of eq:tree "only for the star ... for the path of four densities the step is assumed". FD1-FD3
(`2026-10-21_four_density_star.md`) settled the star. This record tests the path.
**Extends:** `2026-10-21_four_density_star.md`, `2026-10-16_linking_ceiling.md`,
`2026-10-17_density_loop.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Take the fingerprint fibre, densities rho~_v = Delta_v^(-3), Delta_v = sqrt(c_v k^2 + 1), and a path
a - b - c - d: four densities joined by three Coulomb kernels in three dimensions, two of them internal.
Is the minimal annihilator of its one-scale slice built only from the operators of the tree ladder,
F_2, F_3 and F_4 (Paper 59 eq:four_density_op, AESZ #58), and pieces of order one, with F_4's piece
present?

## What I know before computing (disclosed)

- The star (FD1-FD3): S = (16/pi^4) t^3 int chi^3 K0(t x) dx/x, chi = int_0^x s K1(s) ds; Mellin
  operator of order 19; minimal annihilator W8 with composition factors F_4, Legendre and two of order
  one, no K3 piece. My order counts ran low twice (GS1: 17 against <= 10; FD1: 19 against <= 11), so
  no order is predicted here.
- By hand, not computed. In position space a density of scale a is proportional to a^3 K0(a r), and the
  potential of a leaf to chi(a r)/r. Put t on the internal density b, scales (a, b, c, d) = (1, t, 1, 1).
  The angular average of 1/|r - r'| is 1/max(r, r'), and r r'/max(r, r') = min(r, r'), so
  P(t) = N t^3 int_0^oo int_0^oo min(r, r') K0(t r) chi(r) K0(r') chi(r') dr dr'
       = N t^3 int_0^oo K0(t r) chi(r) H(r) dr,  H(r) = int_0^oo min(r, r') K0(r') chi(r') dr',
  with H'' = -K0 chi and H independent of t. N is fixed by the derivation that gives 8/pi^3 for the
  chain and 16/pi^4 for the star. The same steps give the chain, t^3 int K0(t r) chi(r)^2 dr, and the
  star, as checks.
- So P is a single Mellin integral against K0(t r), with g = chi H holonomic: theta(theta - 1) H =
  -r^2 K0 chi. The naive dimension of g's solution space is 3 x 8 = 24, against 10 for the star.
- In one dimension the tree identity (eq:tree) makes the holomorphic period shape-independent; in three
  dimensions the shape enters through H, an extra integration.
- Nothing in this file has been computed.

## Model

P(t) as above. The annihilator of g in theta_r from the exact relations of chi, K0 and H (tensor and
composition, as in `debug/generic_slice_A.py` and `tests/test_paper59_generic_slice.py`). The Mellin route
is `mellin_t_ode` with alpha, beta for the weight, as for the star. Exact arithmetic over Q for operators
and maps; mod p only where stated, with multi-prime CRT and an exact check before any piece is used.
Numerical checks by quadrature of theta_t^j applied under the integral, at 32 digits or more; H by its
two one-sided integrals.

## Predictions

**FP1 (the F_4 piece).** The Mellin operator L of P divides exactly into pieces, and one piece of order
four is equivalent, by an exact map rational in t, to F_4's operator (eq:four_density_op).
- Falsifier: no piece of L equivalent to F_4's operator.

**FP2 (P needs it).** The right factor of L below the F_4 piece does not annihilate P: applied to P it is
nonzero at two values of t, by more than 1e-5 of the largest term.
- Falsifier: it annihilates P to 1e-20.

**FP3 (nothing above).** The minimal annihilator of P, found as the shortest right factor of L that
annihilates P (residual at most 1e-25 of the largest term at 32 digits at four values of t, with a
perturbation of one part in 200 in one coefficient giving above 1e-5), has composition factors each of
order one or equivalent to F_2's, F_3's or F_4's operator.
- Falsifier: the minimal annihilator needs a factor of order at least two equivalent to none of
  F_2, F_3, F_4.

## Novelty status

- FP1, FP3: NOT FOUND. Paper 59's targeted search and the FD1 search found no source placing a molecular
  integral on these operators; the nearest related work is Matuzas 2026 (a different four-electron
  integral, on AESZ-330). A hit on both is a **FINDING** that removes the "assumed for the path" cap.
- FP2: PROBABLY KNOWN in kind (eq:tree); a hit is a CHECK.

## Controls

- Positive: the same reduction with one leaf removed gives the chain's registered order-seven operator;
  the code that builds g, applied to chi^3/r in place of chi H, gives the star's order-19 operator.
- Positive: P at one value of t agrees with an independent momentum-space evaluation of the path to 1e-8
  or better, which also fixes N.
- Negative: the perturbation control of FP3.

## Caps stated in advance

One slice, (1, t, 1, 1), t on an internal density. The fingerprint weights, not Slater ones. The
annihilation is measured, not derived. Minimality is bounded, not proved, unless stated otherwise at
scoring.

## Feasibility rule

If the annihilator of g, its Mellin operator or the factorization cannot be completed within the
session's memory, the computation stops and the record says STOP for feasibility, unscored, with the
stage reached. A switch to the leaf slice (t, 1, 1, 1), where the kernel is chi(t r) and g = K0 H, is
allowed only as a new registration.

## Disposition (appended 2026-10-22, after scoring)

Computed by a sub-agent from this record; scripts `debug/four_density_path_*.py`, data and logs
`debug/data/four_density_path_*`, exact pieces as strings in `debug/data/four_density_path_pieces.json`.

**Model note.** The annihilator of g = chi H was built in the differential module spanned by c_i h_j and
c_i c_j k_l, dimension 18, not 24: the chi in g and the chi in H's source are the same function. Its relation has
order 18 in theta_x and z-degree 37, so L has order 74 (t-degree 74). The minimal annihilator was found by exact
linear algebra in the x-module, using theta_t P_h = P_{(2 - theta) h} and t^2 P_h = P_{z^-1 (theta-1)^2 h}
(search mod p, rebuilt by CRT, x-identity checked exactly over Q).

**Controls.** The module code with chi^3/x gives the star's order-19 operator exactly. One leaf removed
(H -> chi, exact since chi'' = -x K0) gives the chain's registered order-seven operator exactly, and the grid code
with H -> chi gives C(t) to 2.6e-26. Normalization N = 16/pi^4; P agrees with an independent momentum-space
evaluation to all 20 digits at t = 0.7 and t = 1.6. Grid convergence 1e-40. A step-B line first integrated chi in
place of chi^2 (a bug, 4.9 percent off); corrected by `four_density_path_B2_leaf.py`. An extra probe, not
registered: with H_0 kept formal the module operator has Mellin order 15 and the chain's order-seven operator is
not an exact right factor of it.

**FP1: HELD.** L = M62 o W12 exactly (python-flint over Q; a perturbed W12 fails the same division). W12, order 12,
t-degree 78, divides exactly as Y4 o Z2 o Z2' o L1 o R2 o (theta - 2). Y4 has order four and is equivalent to F_4's
operator by an exact odd map of order three: built from B P + q P_h4 = 0, with h4 = x^6 K0^3 = (c1 - c2)^2 k0 in
the module, then R8 P = C P_h4 with R8 = Z2 o Z2' o L1 o R2 o (theta - 2), and G = t^-3 (theta-4)^2 (theta-2)^2
theta^2 from F_4's function to P_h4. Y4 o A = 0 mod Q4 and the image annihilator equals Y4, both exact. B P = -q P_h4
holds to 1e-57 at 60 digits at four t, so no boundary term. Repeated novelty search at scoring: no source places this
integral on AESZ #58; the related Matuzas instance noted in the star record stands. Class: **FINDING**, narrowed
as for the star.

**FP2: HELD.** R8, the right factor below the F_4 piece, leaves 0.10, 0.31, 0.084 and 0.053 of the largest term on P
at t = 0.5, 1.3, 2.5, 3.5. Class: **CHECK**.

**FP3: HELD (minimality bounded).** W12 annihilates P: residual 3.0e-31, 9.8e-31, 9.4e-29, 5.7e-29 at 32 digits at
the four t (4.7e-60 to 3.6e-56 at 60 digits). A 1/200 change in the theta^2 coefficient of Y4 (the F_4 piece, chosen
as in the star record before evaluation) gives 5.5e-3 to 0.50. Composition factors: Y4 (F_4), Z2, Z2' and R2 (each
equivalent to Legendre's operator, F_2, by exact odd maps; each irreducible, since its only log-free solution at
t = 0 has an order-two annihilator), L1 and theta - 2 (order one). No F_3 (K3) factor. The right factors theta - 2,
R2 o (theta - 2), R3 o (theta - 2), Z2' o R3 o (theta - 2) and R8 all fail on P (0.05 to 1.8). Cap: minimality bounded
by the x-module search, no relation of order 11 or less up to t-degree 260, mod one prime. Class: **FINDING**. The path
has the star's F_4 piece and adds two Legendre pieces; like the star it skips the K3 rung.

**Disclosed widening.** The map from Legendre's operator into Z2 gave nullity 0 in the boxes (t^2-1)^6 (t^2-9)^6 and
(t^2-1)^4 (t^2-4)^3 (t^2-9)^3 A10^2 A14^2; it was found after widening to (t^2-1)^6 (t^2-4)^4 (t^2-9)^2 A10^5 A14^3,
numerators t^-23..t^171, with A10, A14 the apparent factors of R3 and Z2'. The minimal box is (t^2-1)^6 A14^2,
numerators t^1..t^49. The map was then checked exactly.

**PM check (2026-10-22).** Rerun with the sympy helpers of `tests/test_paper59_generic_slice.py`, a code
path independent of the agent's operator module, plus an own quadrature (log
`debug/data/four_density_path_pm_check.log`). W12 / R8, R8 / (R3 o (theta - 2)) and R3 / R2 have zero
remainder, with quotients Y4, X4 and order one. X4 / Z2' has zero remainder. Its quotient is Z2 o A14^-2,
not Z2: the stored Q-forms clear denominators on the left, so the stored Z2' carries the factor A14^2.
This is the same composition factor as Z2 under a rational gauge. The three Legendre maps and the F_4 map
are exact, and a 1/1000 change in one map coefficient breaks each. An own grid (H by direct quadrature
at each node, 245 nodes, 32 digits) gives W12 P relative residual 2.3e-31 and 8.8e-29 at t = 1.3 and
2.5. The perturbed W12 gives 0.50 and 0.12, and R8 gives 0.31 and 0.084. Scored in the record's words:
FP1 HELD (FINDING, narrowed as for the star), FP2 HELD (CHECK), FP3 HELD (FINDING, minimality bounded).
