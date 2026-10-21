# Pre-registered prediction: the chain fibre off the equal-leaf slice (GS1-GS5)

**Registered:** 2026-10-21, before any computation in this file (PI-directed: "let's try to work it
down", the punchlist of the 2026-10-21 high-level review, item "the generic chain fibre").
**Prompted by:** the P3e-2W null. Integer-relation searches at four generic scale triples found no
relation in five period bases (Paper 59 `sec:three_electron`), while on the equal-leaf slice (t;1,1) the
fibre is a K3 function integrated once (Paper 59 `sec:density_ring`, eq:chain_slater_slice), a class no
finite basis of periods is expected to capture. This file asks whether the generic fibre has the same
structure, which would explain the null.
**Extends:** `2026-10-19_p3e2_widened.md`, `2026-10-20_slice_artanh.md`, `2026-10-20_slice_parity_proof.md`,
`2026-10-21_wg_closure.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

On the slice with scales (t;1,beta), beta = 2, t on the central density, is the chain's fibre at zero
separation again an elementary function plus a K3 function integrated once against an elementary kernel?

## What I know before computing (disclosed)

- The equal-leaf slice G(t) = int_0^oo x K1(t x) J(x)^2 dx, J = 3 int_0^x K0 - 3 x K0 - x^2 K1, has a
  Mellin operator of order seven dividing as N2 R3 R1 (theta+2); R3 is equivalent to F_3's operator
  4 theta^3 - t^2 (theta+1)^3; W = R3 R1 (theta+2) annihilates G (derived, record wg_closure); no rational
  map splits the K3 piece off (mod-p search, leading-term bound). Paper 59 `sec:density_ring`.
- J obeys theta J = x^3 K2 and theta ((theta-3)^2 - 4 - x^2) J = 0, an order-three operator whose solution
  space contains the constants.
- By scaling x -> x/lambda, G(lambda a; lambda b, lambda c) = lambda^(-2) G(a;b,c) for the integral
  int x K1(a x) J(b x) J(c x) dx. So the slice (t;1,2) passes through the P3e-2W triples (1;2,4) at
  t = 1/2 and (3;1,2) at t = 3. The second is a threshold triple: the central scale is the sum of the
  others, where the fibre did not collapse (P3e-4 MISS).
- The three-scale moment int K0(a x) K0(b x) K0(c x) dx is a sum of products of complete elliptic
  integrals at the Bailey moduli (BBBG 2008 Eq. 33). In the sunrise language its singular points in the
  masses are the thresholds a = b + c and pseudo-thresholds a = |b - c|.
- Counting by hand, not computed: the solution space of J(x) J(beta x)'s annihilator is the tensor
  product of two three-dimensional spaces, each containing 1. It splits as 1 + 2 + 2 + 4, so I expect the
  Mellin operator to have order at most 10, with first-order pieces for 1 x 1, J x 1 and 1 x J, and the
  4-dimensional K0-type tensor block giving the K3 piece.
- Nothing in this file has been computed.

## Model

G_b(t) = int_0^oo x K1(t x) J(x) J(beta x) dx, beta = 2. The annihilator of J(x) J(beta x) in theta_x is
the tensor product of J's operator and its rescaling (z -> beta^2 z, z = x^2), built as in
`debug/density_loop_B_ode.py` sym_power but with two different relations. The Mellin route is
`mellin_t_ode` with alpha = 2, beta = 0, r(sigma) = sigma^2 - 1 (E = M[K1]). The K3 comparison operator is
the minimal operator in t of V(t) = int_0^oo K0(t x) K0(x) K0(beta x) dx, built the same way. Factoring,
equivalence maps and split searches follow the slice scripts (`debug/p3e2w_slice_*.py`,
`debug/slice_artanh_*.py`). Exact arithmetic over Q for operators and maps; mod-p only where stated.

## Predictions

**GS1 (shape of the operator).** The Mellin operator L of G_b has order at most 10. It has a right factor
of the form E3 = R1' R1 (theta+2) or a reordering, built from three first-order pieces with algebraic
solutions: t^(-2), one singular only at t = 1, and one singular only at t = 2 (the leaf scales).
- Falsifier: order above 10, or no such right factor.

**GS2 (the K3 piece).** Above E3, L has an order-three or order-four piece R equivalent, by an exact map
rational in t, to the minimal operator of V(t).
- Falsifier: no piece equivalent to V's operator, or a piece of order above 4 needed below the top.

**GS3 (annihilation).** R E3 annihilates G_b: residual at most 1e-25 of the largest term at 32 digits at
four values of t between 0.4 and 3.5, while a change of one part in 200 in one coefficient of R gives a
residual above 1e-5. G_b needs the K3 piece: E3 G_b is nonzero at those points.
- Falsifier: residual above 1e-20 at any point, or E3 G_b = 0.

**GS4 (no split).** No rational map splits the K3 piece off: a search modulo a prime, in a box at least
the bound set by the leading terms of the local exponents, finds no split, and the same call finds a
planted split. So G_b is an elementary function plus a K3 function integrated once.
- Falsifier: a split exists (then G_b is a rational combination of V and its derivatives plus
  elementary terms, and the P3e-2W bases should have caught it at (1;2,4)).

**GS5 (thresholds).** The singular points of R in t > 0 are among {1, 3} = {|beta - 1|, beta + 1}, apart
from apparent singularities, and t = 3 is a genuine singular point of the minimal annihilator of G_b.
- Falsifier: a genuine singular point of R outside {1, 3}, or t = 3 apparent.

## Novelty status

- GS1, GS3: PROBABLY KNOWN in form. Operator factorizations for sunrise-type integrals with unequal masses
  are standard (Mueller-Stach, Weinzierl, Zayadeh 2012, the unequal-mass sunrise; not re-checked for this
  integral). Hits are CHECKs.
- GS2, GS5: PROBABLY KNOWN. The three-scale Bessel moment's family and its thresholds are BBBG 2008 and the
  sunrise literature.
- GS4: NOT FOUND. No source places a molecular chain integral off the equal-leaf slice on these operators
  (Paper 59's targeted search). A hit is a **FINDING**: it extends the slice result to a generic slice and
  explains the P3e-2W null. The search is repeated at scoring.

## Controls

- Positive: with beta = 1 the code reproduces the registered order-seven operator QG of the slice
  exactly, and V's operator reduces to F_3's.
- Positive: the split search finds a planted split (as in the slice).
- Negative: the perturbation control of GS3.
- Informational: the ratio G_b(1/2) / G_b(3), with the scaling factor 2^(-2), is compared with the P3e-2W
  grid values at (1;2,4) and (3;1,2). A mismatch is reported with its convention, not scored.

## Caps stated in advance

One slice, beta = 2. Slater weights only. The annihilation in GS3 is measured; a derived version would
follow the wg_closure route and is not part of this registration.
