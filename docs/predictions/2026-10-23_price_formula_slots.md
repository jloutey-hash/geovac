# Pre-registered prediction: the price formula off the symmetric slice (PF1-PF4)

**Registered:** 2026-10-23, before any computation in this file (PI-directed: "Let's keep knocking
down these open items").
**Prompted by:** the atlas door 27 Owed line, "a correlated integral outside this family whose
transcendental the formula predicts", and Paper 59 `sec:density_ring`'s Observation: "composing unequal
decay rates costs one elliptic period ... at imaginary momentum it is set by the angle of their
triangle". Every test of that formula so far sits on the symmetric slice: base densities of equal
scale (1, 1), two equal geminals, one function L(t, 1; i gamma). On that slice the formula cannot say
which density pairs with which geminal, because the two pairings coincide. This file leaves the slice.
**Extends:** `2026-10-22_c2b_feynman_g2.md` (CF5-CF7), `2026-10-23_gamma1_mechanism.md` (MG1-MG3),
`2026-10-21_generic_slice.md` (GS1-GS5), `2026-10-22_geminal_shapes.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

The correlated triangle Y(t) = <rho_i rho_j rho_k f(r_ik) f(r_jk) / r_ij>, exponential geminals, with
the two base densities at different scales and the two geminals at different decay rates: scales
(a_i, a_j, a_k) = (1, 2, t), t on the apex k, gamma_1 = 3/2 on the edge i-k and gamma_2 = 3 on the
edge j-k. Does the one function L(a, b; c) still carry every elliptic piece, and with which slots?

This is the same integral family generalised, not a different topology. What it tests is the content
of the formula: that each elliptic piece is a bubble of two density scales at an imaginary momentum
set by a geminal, and which two scales and which geminal. A different topology is a later item.

## What I know before computing (disclosed)

- On the slice (1, 1, t) with one gamma, Y's minimal annihilator has order nine with composition
  factors 1, 1, 1, 1, 2, 3: four algebraic pieces branched at t = 2 and t = 2 gamma; one Legendre piece
  with z = (t^2 + 1 - gamma^2)/(2t) twisted by x^(1/2), x = t^2; F_3's operator twisted by t. Genuine
  singular points in t > 0: {|1 - gamma|, 1 + gamma, 2, 2 gamma} (gamma = 1/2, 3/2, 2, 3; records CF,
  MG). At gamma = 3 the Legendre map has E = ((x - 4)(x - 16))^4 (x + 8) (`tests/data/paper59_mg_2026_10_23.json`).
- The kite (route B of record QG) is the two-loop self-energy with the bubble {rho_i, rho_j}, the bubble
  {f, f} and the Coulomb line between them; the apex density closes it: Y(t) = (1/(4 pi^2)) t^3
  int z^(1/2) (z + t^2)^(-3/2) phi(z) dz. Its Feynman form has M^2 = 1 + u(1-u) z for the density pair
  and N^2 = gamma^2 + v(1-v) z for the geminal pair (`debug/ga1_0_kernel.py`, `debug/ring_ext_QG_common.py`).
  The derivation there already carries two geminal masses: c = (1-v) g1^2 + v g2^2 + v(1-v) L^2.
- The chain of three on the slice (t; 1, 2), record GS: its K3 piece is R_V, the order-four minimal
  operator of V(t) = int K0(tx) K0(x) K0(2x) dx, singular at t = 1 and 3 (`debug/generic_slice_A.py`,
  `debug/data/generic_slice_A.json`, key `QV`); at leaves (1, 1) it is F_3's operator.
- **Derived by hand while drafting, nothing computed.** Landau cuts of the kite graph. Vertices O
  (the common centre), 1, 2, 3; lines O1 = a, O2 = b, 12 = 0 (Coulomb), 13 = gamma_1, 23 = gamma_2;
  external momentum w from O to 3, z = w^2. The cuts separating O from 3 are {O1, O2}, {13, 23},
  {O1, 12, 23} and {O2, 12, 13}. Their normal and pseudo-thresholds put the kite's singular points at
  z = -(a +- b)^2, -(gamma_1 +- gamma_2)^2, -(a +- gamma_2)^2, -(b +- gamma_1)^2, and Y's at x = t^2 equal
  to the same squares. At (1, 1; gamma, gamma) this is {-4, -4 gamma^2, -(1 +- gamma)^2}, the stored set
  (MG3), and the Legendre piece's thresholds z_Leg = +-1 at t = 1 +- gamma are the cut {O2, 12, 13}:
  **the base density on the far side of the Coulomb line from the geminal.** So off the slice the
  formula's slots are crossed: the bubble (t, b) carries the momentum i gamma_1 and the bubble (t, a)
  carries i gamma_2. With (a, b, gamma_1, gamma_2) = (1, 2, 3/2, 3):
  - z_1 = (t^2 + 4 - 9/4)/(4t) = (4 t^2 + 7)/(16 t), thresholds t = 2 +- 3/2 = 7/2, 1/2;
  - z_2 = (t^2 + 1 - 9)/(2t) = (t^2 - 8)/(2t), thresholds t = 1 +- 3 = 4, 2: the same Legendre family
    as the slice at gamma = 3 (record MG), although the other scale and the other geminal differ.
  The direct pairing would give instead z = (t^2 - 5/4)/(2t) (the slice's gamma = 3/2 family, CF2) and
  z = (t^2 - 5)/(4t). The two readings differ in their singular sets too (below).
- The numbers were chosen so that the eight predicted points are distinct: (a +- b)^2 = 9, 1;
  (gamma_1 +- gamma_2)^2 = 81/4, 9/4; (a +- gamma_2)^2 = 16, 4; (b +- gamma_1)^2 = 49/4, 1/4. Under the
  direct pairing the set would be {1/4, 1, 9/4, 25/4, 9, 81/4, 25}.
- Radius: on the principal sheet M vanishes first at z = -(a+b)^2 = -9 (the pseudo-threshold -1 is not
  reached, as -(1-gamma)^2 was not on the slice), N at -81/4, D at the two thresholds; so the Cauchy
  radius is 9, larger than the slice's 4. To be checked by a second radius.
- Expected size: more algebraic pieces (two branch points for each bubble now), two Legendre pieces,
  one order-four K3 piece: order about sixteen and x-degree well above the slice's 32. Guessing Y's
  operator directly may be out of reach; the kite's operator in z (order about twelve) and the apex
  transform c_n = binom(-3/2, n) are the cheaper route (`debug/ga3_5e_modp_guess.py` uses it).

## Model

The kite of route B with M^2 = u a^2 + (1-u) b^2 + u(1-u) z, N^2 = (1-v) gamma_1^2 + v gamma_2^2 + v(1-v) z
and the kernel re-derived from `debug/ga1_0_kernel.py` with two geminal masses (the exp shape:
d/dgamma_1 d/dgamma_2 of the Yukawa-Yukawa kite), validated against route A (the closed-form loop of
`debug/ring_ext_QG_common.py` with masses (gamma_1, gamma_2, 0), itself checked against direct 3-D
quadrature) at three w. Taylor coefficients of phi exactly modulo primes (`debug/ga3_5d_modp_series.py`
generalised) and by Cauchy's formula at 600 digits; operator guessing in x = t^2 (Q-form, theta = x d/dx
conventions of the ga3 scripts); exact composition series by the recursive builder of `debug/ga3_lib.py`;
Legendre and R_V maps by the Stage E searches (`debug/ga3_map.py`, `debug/ga3_10_k3map.py`, `debug/generic_slice_D.py`);
annihilation on Y at 32 digits (`debug/ga3_3_grid.py`, `debug/ga3_9_onY.py`). Route choice is the
agent's; the kite-first route is allowed and expected.

## Predictions

**PF1 (main: the slots are crossed).** Y's minimal annihilator has exactly two composition factors of
order two, inequivalent to each other, each equivalent by an exact order-one map over Q(x), with a twist
by a power of x, to the Legendre operator of degree -1/2 in z_1 = (4 t^2 + 7)/(16 t) and in
z_2 = (t^2 - 8)/(2t) respectively.
- Falsifier: the direct pairing (z = (t^2 - 5/4)/(2t) and (t^2 - 5)/(4t)); fewer or more than two
  elliptic factors; a factor of order two whose curve is isogenous to neither (traces of Frobenius at
  several rational x and primes, as in record GA, with the two predicted families as controls).
- Confidence: moderate to high. The Landau reading fixes the thresholds; that the piece is the one
  function with those slots is the formula.

**PF2 (the singular points are the eight cuts).** The genuine singular points of Y's minimal
annihilator in x = t^2 > 0 are exactly {1/4, 1, 9/4, 4, 9, 49/4, 16, 81/4}; the kite's operator in z has
genuine points exactly minus that set (plus 0 and infinity). Kite order not predicted.
- Falsifier: a genuine point outside the set; any of the eight regular or apparent.
- Confidence: high for the six thresholds; moderate for the two bubble pseudo-thresholds 1 and 9/4,
  which on the slice sit at 0 and were not separately visible (GS saw t = 1 genuine for the chain).

**PF3 (the K3 piece is the generic-slice chain's).** Y's annihilator has a composition factor equivalent
by an exact map to R_V at beta = 2 (order four), which Y needs (the right factor below it does not
annihilate Y), and no factor above it.
- Falsifier: no such factor; Y annihilated without it; any irreducible factor of order five or more.
- Confidence: high. The apex and the two base densities are the chain (t; 1, 2).

**PF4 (no new type).** Every composition factor is of order one with an algebraic solution (branch
points among x = 1, 9, 9/4, 81/4, half-integer exponents), one of the two Legendre families of PF1, or
the R_V piece of PF3.
- Falsifier: any other factor.

## Novelty status

- PF2: KNOWN in form. Landau's equations; for the kite, Broadhurst (LL2024 notes, indico.desy.de
  event 38959) lists pairwise and triple mass sums as thresholds; pseudo-thresholds not discussed. A hit
  is a CHECK.
- PF1: NOT FOUND. Searched 2026-10-23 (web): "two-loop kite self-energy unequal masses three dimensions
  Landau pseudo-threshold elliptic Picard-Fuchs"; "explicitly correlated Slater geminal three-electron
  integrals elliptic Picard-Fuchs". The kite literature (Adams-Bogner-Schweitzer-Weinzierl arXiv:1607.01571
  and 1705.08952; "the soaring kite" arXiv:2401.14307, which lists the Landau variety of the five-mass
  kite in D = 4 - 2 epsilon) ties the elliptic curves to the sunrise subsectors in four dimensions; no
  three-dimensional case and no molecular placement. The Slater-geminal integral literature (Ten-no
  2004; Ten-no-Manby 2003) has no operator analysis. A hit is a **FINDING**: the formula's slot content
  off the slice it was read from. The search is repeated at scoring.
- PF3: PROBABLY KNOWN in form (two-parameter K3 families of the unequal-mass sunrise, arXiv:2502.15326);
  the molecular placement is record GS's. A hit is a CHECK.
- PF4: a CHECK.

## Controls

- Positive: the generalised kite at (1, 1; 2, 2) reproduces the stored gamma = 2 exponential series
  (`debug/data/cf_g2_series_g2_M1600_d600_rho2.json`) coefficient by coefficient, and at (1, 1; 3, 3) the
  record-MG series at gamma = 3 (`debug/data/mg_*`); route A against route B at the registered masses to
  1e-14 at three w; the loop against direct 3-D quadrature at five triples to 1e-8.
- Positive: `debug/generic_slice_A.py` rebuilds R_V and matches the stored `QV`; the Legendre map search
  recovers the stored gamma = 3 map on the slice's data before it is run here; the planted-split control
  of record g32 on each irreducibility claim.
- Negative: a 1/1000 change in one coefficient of each accepted map breaks it; a 1/200 change in one
  coefficient of the annihilator leaves a residual above 1e-6 of the largest term on Y (reported in the
  GA1 normalisation and in the dominant-monomial normalisation, as in record GA3).
- Minimality: no operator of lower order to the degrees searched, stated with the degrees, mod p.

## Caps stated in advance

One configuration, one shape (exponential), one placement (t on the apex). The pairing is read from
Landau cuts of the kite, a reading checked by its consequences, not a derivation of the pieces. If no
annihilator of Y is recovered with its controls passing, that is a STOP for feasibility, reported and
not scored; in that case the kite's operator alone may still score PF2's kite clause, and this is said
so. The direct pairing, if it appears, is a MISS for PF1 and a finding about the formula, recorded as
such.
