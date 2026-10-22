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

## Disposition (appended 2026-10-21, after scoring)

Scripts: `debug/generic_slice_{A,A2,B,C,C2,D,D2,E,F2,G}.py`, `debug/generic_slice_num.py`; data and logs
`debug/data/generic_slice_*`. `generic_slice_F.py` (symbolic-rho exponents) was stopped for time on W7 and
replaced by the d/dt indicial route `generic_slice_F2.py`, which reproduces its R_V and R exponents.

**Controls.** b = 1 reproduces the registered QG exactly (ratio 1), and the K0 x K0 x K0 tensor at b = 1
gives F_3's operator. The split search finds the planted split (E3 o R_V, nullity 1) at the scoring box.
Informational: the P3e-2W grid values satisfy G(3;1,2) = G_b(3) and G(1;2,4) = G_b(1/2)/4 to 20 digits,
so this slice passes through two of the null triples.

**GS1: FAILED.** L has order 17 and t-degree 18, above the registered 10. The structural part held:
theta + 2 is an exact right factor, and the two first-order pieces with solutions t/(1-t^2)^4 and
t/(4-t^2)^4 are each exact right factors of L/(theta+2). They do not peel one after the other in either
order; they enter as their least common left multiple E2 (order 2), and L = M14 o E2 o (theta+2)
exactly. Class: **MISS** (the order count was wrong).

**GS2: HELD.** V's Mellin operator R_V has order 4 (exponents 0, 0, 2, 2 at t = 0; singular at t = 1, 3).
A gauge of order three, rational in t with poles only at t = 0, +-1, +-3 and numerators of degree at
most 22, maps Sol(R_V) into Sol(M14); its annihilator R has order 4 and M14 = N10 o R with zero
remainder, in exact rational arithmetic. Found mod p in the box t^-8 (t^2-1)^-6 (t^2-4)^-6 (t^2-9)^-6,
degree 60, rebuilt from three primes by CRT and rational reconstruction. A first box (poles of order 3,
degree 30) found nothing; disclosed, and expected from the b = 1 slice's map, which has poles of order 5.
Cap: the irreducibility of R_V is shown only against factors of order one and three (Fuchs exponent
count); an order-two factor is not excluded by exponents. Class: **CHECK**.

**GS3: PARTLY.** W7 = R o E2 o (theta+2) annihilates G_b at 32 digits: relative residual 2.8e-29 at
t = 0.5 (the quadrature floor there, 2.0e-29) and at most 5.2e-33 at t = 1.3, 2.5, 3.5. G_b needs the K3
piece: E3 G_b is 0.13 to 0.49 of its largest term. The perturbation clause, with the theta^1
coefficient of R times 1 + 1/200, gave 2.9e-4 and 3.2e-4 at t = 0.5, 1.3 but 4.9e-6 and 2.5e-6 at
t = 2.5, 3.5, below the registered 1e-5 at two points, though 27 orders above the residual. No other
coefficient was tried after seeing this. Also: the order-17 L annihilates G_b with no source term, at
about 1e-37 of its largest term at t = 1.3, 2.5, 3.5 (t = 0.5 unresolved on the grid). Class: **CHECK**
for the annihilation; the perturbation clause is a control shortfall.

**GS4: HELD.** Leading-term bound from the local exponents: pole order at most 4 at t = 0 and 6 at
t = +-1, +-2, +-3, numerator excess 1 at infinity. In a box beyond it (t^-8, poles of order 10, degree
65) no gauge maps Sol(R_V) into Sol(W7) (nullity 0), while the planted split is found. So on the slice
(t;1,2) the fibre is an elementary function plus a K3 function integrated once, as on (t;1,1). Caps:
mod p; the bound is the leading-term argument, not a written proof. Novelty search repeated at scoring
(2026-10-21, two extended web searches): no source places a molecular chain integral on these operators;
the nearest are Feynman-integral papers on two-parameter K3 families (arXiv:2502.15326) and the
four-mass banana (arXiv:2507.23061), neither molecular. Class: **FINDING**. It explains the P3e-2W null at
(1;2,4) and (3;1,2): no finite basis of periods is expected to capture a K3 function integrated once.

**GS5: HELD.** The singular points of R in t > 0 are t = 1 and 3, each with a half-integer exponent
(-13/2, -15/2), plus the roots of a degree-10 polynomial, apparent since R is equivalent to R_V; t = 2
is regular. W7 has exponent -9/2 at t = 1 and at t = 3. That t = 3 is genuine for the minimal
annihilator of G_b holds if R is irreducible (the GS2 cap). Class: **CHECK**.
