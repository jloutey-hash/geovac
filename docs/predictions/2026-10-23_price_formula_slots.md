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

## Disposition (appended 2026-10-24; predictions above unchanged)

**PF1-PF4: STOP for feasibility. All four UNSCORED.** No annihilator of the kite phi(z) was found within the
cap, so none of Y(t) was found either. By the caps above this is reported and not scored. The kite clause of
PF2 is not scored either, since it needs the kite's operator.

- The search: 300 Taylor coefficients of phi at 280 digits, tail accurate to 1.7e-216 or better. At orders 6, 9
  and 12 the residual falls smoothly as the number of unknowns grows (r = 6: 2.7e-29 at 63 unknowns to 2.1e-191
  at 259; r = 9: to 5.1e-194 at 250; r = 12: to 1.3e-196 at 247). There is no sharp drop to the data's precision.
  The earlier hits at 1e-151 to 1e-174 sat just under a 1e-150 tolerance and are threshold solutions, not
  operators. So the kite has no operator with 260 unknowns or fewer. A successor run would need about 720
  coefficients at about 600 digits.
- Partial evidence, not scored. On the principal sheet the true kite has Taylor radius 9 with leading exponent
  -1. The direct-pairing object (the same kite with its two geminals swapped) has radius 25/4 with a square-root
  branch point there. So the true kite is analytic at z = -25/4 on the principal sheet, and so are the four
  pseudo-threshold points -1/4, -1, -9/4 and -4. This leans toward the crossed pairing of PF1. It does not settle
  PF1 or PF2, because a genuine singular point of an operator can sit on another sheet. The points -25, -49/4,
  -16 and -81/4 lie beyond the radius and were not tested.
- PM recheck: an independent three-term ratio fit on the stored series gives radius 9.00025 and exponent -1.016
  at n = 297 for the true kite, and radius 6.25008 and exponent -0.503 at n = 125 for the swapped object. The
  agent's route A, built from the definition, matches the derived kite to 1.7e-21 or better at three w and rejects
  the swapped form by 1.2e-2 or more.
- Controls that passed: the loop against direct 3-D quadrature at five triples (2.4e-15); the slice members
  (1, 1; 2, 2) and (1, 1; 3, 3) against the stored series (1.9e-47 and 5e-77); the stored gamma = 3 Legendre map
  recovered and broken by a 1/1000 change; R_V's identity map exact; the kite-first route reproduces the stored
  gamma = 3 operator L9 (2.2e-362 on 800 coefficients).
- Deviations: the registered Cauchy series was stopped (it would have run more than 50 hours) and replaced by
  direct power-series arithmetic in z, checked against the Cauchy coefficients to 1.7e-76. Off the slice the
  Taylor coefficients are not in pi Q (no relation to height 1e20), so the exact mod-p route does not exist. The
  1/200, minimality and planted-split controls were not run, because nothing was accepted.
- Reading: none. The formula's slot content off the symmetric slice stays open. Novelty search not repeated, since
  nothing scored.
- Memo `debug/sprint_pf_memo.md`; data `debug/data/pf_pieces.json` (operator and map keys empty, with the
  reason); series `debug/data/pf_series_direct_explore_N300_d280_460x460_sub.json` and
  `debug/data/pf_series_swapped_g_M256_d60_rho3.json`.

## Disposition addendum (2026-10-24, after record IB; predictions above unchanged)

Record IB (`2026-10-24_ibp_kite.md`) derived the kite's operator at (1, 2; 3/2, 3) by integration by parts:
order 11, z-degree 102 to 113, exact over Q, on the stored 300-coefficient series to 1.3e-280. By the caps
above, the kite's operator alone may score PF2's kite clause, and it is scored here.

- **PF2, kite clause: HELD.** The genuine singular points of the kite's operator are exactly
  z = -1/4, -1, -9/4, -4, -9, -49/4, -16, -81/4, with 0 and infinity. Exact local exponents over Q: -1/2 at
  -1/4, -4, -49/4 and -16 (the three-line cuts); a double exponent 0, so a logarithm, at -1, -9, -9/4 and
  -81/4 (the two-line cuts); at 0 the exponents include -4 and half-integers; at infinity half-integers, two
  of them double. At the direct-pairing points -25/4 and -25 the leading coefficient does not vanish and the
  exponents are 0 to 10: ordinary points. The rest of the leading coefficient is one irreducible factor of
  degree 92 over Z, and at its roots all 11 solutions are analytic (counted modulo a prime where it has a
  root): apparent. The same count gives 8 or 10 of 11 at the eight cuts and 11 of 11 at -25/4, -25 and a
  control point.
- **PF1, PF3, PF4 and PF2's Y clause: still UNSCORED.** They are stated for Y's annihilator. The reduction on
  the tetrahedron off the slice (34 masters) needs more than 80 sample points; a run with 170 is in progress.
- Reading, not scored: at the kite level the crossed pairing holds and the direct pairing has no
  singularity at all.

## Status note (2026-10-24, later; predictions and dispositions above unchanged; nothing scored here)

The 170-sample reduction on the tetrahedron off the slice finished sampling. The fit then failed on a bug in the
PM's fitter, not in the data: the "random" normalisation vector had all 400 entries equal, so it fixed D(1), and
every denominator here vanishes at x = 1. Fits accepted earlier were each checked on every sample and on held-out
points, so they stand. Fixed, and refitted from the saved samples: 34 masters, entries of total degree at most 35,
held-out points pass. The operator for Y/t^3 in x = t^2 has order 19 and x-degree 102 to 121, modulo one 61-bit prime.

Preview at that one prime, not a score:
- The rational roots of the leading coefficient are exactly x = 0, 1/4, 1, 9/4, 4, 9, 49/4, 16, 81/4.
- At 1/4, 4, 49/4 and 16 the exponents are -1 and 0 to 17, and 18 of 19 solutions are analytic; the 19th is not
  a Laurent series, so it carries a logarithm. At 1, 9, 9/4 and 81/4 there are three half-integer exponents and
  16 of 19 analytic solutions.
- x = 25/4 and 25, the direct-pairing points, are ordinary (19 of 19).
- The remaining degree-98 part of the leading coefficient splits over F_P into irreducible factors of degrees
  8, 13, 19, 20 and 38. Counted in each extension field, all 19 solutions are analytic at their roots: apparent.
  Control: the same code at the eight cuts reproduces the counts above.

Why this is not a score: the operator comes from one prime and has not yet been checked against Y itself (the
control above requires annihilation on Y with the 1/200 perturbation), and order 19 may not be minimal (about 16
was expected; at gamma = 3/2 the reduction gave order 10 against the minimal 9). The genuine points of the minimal
annihilator lie among these, but that each of the eight is genuine for it is not yet shown. Lifting the system
needs samples at a second and a third prime: one prime leaves about 450 of 4275 coefficients at the 30-bit
reconstruction limit and 243 unreconstructed.

## Status note 2 (2026-10-25; predictions and dispositions above unchanged; nothing scored here)

The tetrahedron off the slice was sampled at two more primes (same 34 masters, order 19 at each).

- Lift. The system matrix A was rebuilt over Q from two primes and agrees with the third. The target's own
  expression in the masters, c, needed all three primes (heights up to 78 bits, against about 91 for
  reconstruction from three), so c is not held out; its independent check is the check on Y below. The operator
  for Y/t^3 in x = t^2 is exact over Q: order 19, x-degree 102 to 121, integer coefficients up to 934 bits.
- On Y (the 32-digit grid). Residual of the largest t^(2i) group: 6.1e-24, 5.1e-28, 1.2e-24, 1.8e-16 at
  t = 0.8, 1.7, 2.6, 3.8; on the grid with twice the step, 2.0e-8, 1.1e-10, 4.5e-10, 5.5e-5, so the residual
  follows the quadrature error. A 1/200 change in the dominant-monomial coefficient leaves 1.0, 1.2, 1.6, 1.1. In
  GA1's normalisation it leaves 1.0, 1.3e-5, 1.4e-19, 1.8e-16: that coefficient's monomial is negligible at large
  t, so in that normalisation the control tests only at t = 0.8 and 1.7.
- Points, exact. The leading coefficient factors over Z as x^7 (4x - 1)(x - 4)(4x - 49)(x - 16)
  [(x - 1)(x - 9)(4x - 9)(4x - 81)]^3 f, with f irreducible of degree 98. The exponents over Q are those of the
  preview. At a fresh prime, not one the operator was built from, f splits into degrees 5, 27 and 66, and all 19
  solutions are analytic at its roots: apparent. The series counts at the eight points are unchanged.
- PF2, first arm of the falsifier (a genuine point outside the set): excluded for x > 0. Y's minimal annihilator M
  divides L on the right, so every solution of M is a solution of L, and a point where all solutions of L are
  analytic is not genuine for M.
- PF2, second arm (each of the eight genuine for M): not yet decided. Method, chosen at scoring time and not named
  in the record: every analytic continuation of Y solves M, so x0 is genuine for M if some continuation of Y
  changes on a small loop around x0. L was continued numerically at 40 digits, from Y's values at t = 0.8 on the
  grid, through the upper half-plane. Two sheets were tried: Y itself, which as a convergent integral is analytic
  for Re t > 0 and so is the control, and Y continued once around x = 0. Change of the scaled vector on the loop:

  | x0 | Y itself (control) | around 0 |
  |:--|:--|:--|
  | 1/4 | 1.6e-30 | 5.9e-22 |
  | 1 | 6.0e-16 | 1.2e-16 |
  | 9/4 | 8.9e-11 | 2.5e-12 |
  | 4 | 5.7e-13 | 3.9e-15 |
  | 9 | 1.5e-8 | 2.0 |
  | 49/4 | 1.8e-11 | 0.080 |
  | 16 | 9.6e-9 | 0.14 |
  | 81/4 | 0.19 | 2.0 |
  | 3, 25/4, 25 (ordinary) | 3.6e-41, 2.2e-41, 1.2e-41 | 3.5e-41, 2.1e-41, 2.8e-41 |

  So 9, 49/4 and 16 are genuine for M. At 1/4, 1, 9/4 and 4 both sheets are analytic to the noise, so nothing is
  decided there: other sheets have not been tried. At 81/4 the control itself is 0.19, so the test does not
  decide there. The noise is the grid's error, about 1e-24, carried by solutions of L that grow along the path.
  Monodromy matrices from exact identity data were also computed at 60 digits, but the transport to the far points
  has condition number above 1e60, so they could not be used at that precision.
- In progress: a 60-digit grid at half the step, to rerun the test with about 50-digit data.

## Disposition addendum 2 (2026-10-25; predictions above unchanged)

The operator for Y off the slice is exact over Q and checked on Y (status note 2). That is enough to score
PF2's Y clause. PF1, PF3 and PF4 need the factorisation of that operator and stay UNSCORED.

- **PF2, Y clause: HELD.** The genuine singular points of Y's minimal annihilator M in x = t^2 > 0 are exactly
  1/4, 1, 9/4, 4, 9, 49/4, 16 and 81/4.
  - First arm, no genuine point outside the set. M divides the exact order-19 operator L on the right, so every
    genuine point of M is a genuine point of L. For x > 0 those are the eight cuts. The other roots of L's leading
    coefficient are apparent (status note 2).
  - Second arm, each of the eight genuine for M. Every analytic continuation of Y solves M, so x0 is genuine for
    M if some continuation of Y changes on a small loop around x0. This method was chosen at scoring time and is
    not named in the record.
- Check on Y, on a new 60-digit grid at half the step (247 nodes). Residual of the largest t^(2i) group: 3.2e-45,
  1.6e-51, 4.9e-55, 3.3e-43 at t = 0.8, 1.7, 2.6, 3.8. With the old step: 4.2e-29, 5.1e-28, 1.2e-24, 1.8e-16. So
  the residual still follows the quadrature error. A 1/200 change in the dominant monomial leaves 1.0, 1.18, 1.55,
  1.14. In GA1's normalisation it leaves 1.0, 1.3e-5, 1.4e-19, 7.3e-40; as in status note 2, that monomial is
  negligible at large t.
- The loops. L was continued by Taylor stepping from Y's values at t = 0.8 on the 60-digit grid, through the upper
  half-plane, to the base point x = 0.12. A sheet is a word of loops at the base point: "0, 9" means once around
  0, then once around 9. The number is the change of the scaled vector on a small loop around x0, relative to its
  size. A loop around an ordinary point leaves every continuation unchanged, so on such a loop the number is the
  numerical error of that path. The controls are chosen so that their path is at least as long as the deciding
  one.

  | x0 | deciding sheet | change, 64 / 90 digits | control | control, 64 / 90 digits |
  |:--|:--|:--|:--|:--|
  | 9 | 0 | 2.0 / 2.0 | Y itself | 5.9e-32 / 7.1e-37 |
  | 49/4 | 0 | 0.07967 / 0.07967 | Y itself | 3.0e-35 / 1.5e-39 |
  | 16 | 0 | 0.1332 / 0.1332 | Y itself | 3.3e-33 / 5.6e-38 |
  | 81/4 | 0 | 2.0 / 2.0 | Y itself | 5.2e-28 / 9.4e-33 |
  | 1/4 | 0, 9 | 0.4999 / 0.4999 | 0, 10 | 9.7e-15 / 8.2e-18 |
  | 1 | 0, 9 | 2.0 / 2.0 | 0, 10 | 3.3e-13 / 2.9e-16 |
  | 9/4 | 0, 16 | 2.0 / 2.0 | 0, 17 | 4.6e-4 / 9.8e-11 |
  | 4 | 0, 16 | 0.1564 / 0.1564 | 0, 17 | 6.7e-11 / 1.5e-17 |
  | 3, 25/4, 25 (ordinary) | every sheet on the 60-digit grid | 1.4e-62 or less | | |

  Each deciding value is the same to four digits at both precisions and lies at least 1e10 above its control at
  90 digits. The single-loop values also agree on the old 32-digit grid at 40 digits: 2.0, 0.07963, 0.1379, 2.0.
- Noise, disclosed. Three readings made during the run were wrong and are not used.
  - A pilot on the 32-digit grid gave 0.0042 at 1/4 on the sheet "0, 49/4". On the 60-digit grid the same sheet
    gives 2.7e-14 there. The pilot was noise.
  - On the 32-digit grid every two-loop sheet gives about 0.0026, 2.0, 2.0 and 0.27 at 1/4, 1, 9/4 and 4, whichever
    far point the second loop goes round, and so does the control around the ordinary point 25. All of it is noise.
  - On the 60-digit grid at 64 digits the sheet "0, 81/4" and the control "0, 25" are also noise (at 1/4, 0.0358
    and 0.0357). At 90 digits the control "0, 25" still leaves 1.03 at 9/4. So no two-loop sheet that goes beyond
    17 is used.
  - A first check by the PM used controls of Y itself for the two-loop sheets. These underestimate the noise there,
    since Y does not carry the components that the loop around 0 brings in. The cross-check on the second grid
    found this, and the ordinary-point controls above replace them.
  - Starting from t = 1.3 instead of 0.8 gives the same values to every printed digit, so on the 60-digit grid the
    data error does not limit. The noise falls with the working precision, by a factor 1e3 to 1e7 from 64 to 90
    digits, which is less than the 26 digits added. The cause is not pinned.
- Reading, not scored: which loop exposes which point is a statement about the pieces of Y, and that is PF1 and
  PF4's question.
- **PF1, PF3, PF4: still UNSCORED.** They need the factorisation of L and whether order 19 is minimal.
- Script `debug/ib_tetra_pf_sheets.py`. Data `debug/data/ib_tetra_pf_sheets_*.json` and logs,
  `debug/data/pf_9_onY_tetra_pf_exact_d60.json`, grid `debug/data/pf_grid_d60_h0.03125.json`.
