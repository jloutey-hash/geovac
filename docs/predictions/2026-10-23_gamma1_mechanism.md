# Pre-registered prediction: why the correlated triangle drops at gamma = 1 (MG1-MG3)

**Registered:** 2026-10-23, before any computation in this file (PI direction: "let's start knocking
out the open items", the first being the gamma = 1 mechanism).
**Prompted by:** the CF record's untested inference "gamma = 1 is a degenerate member, where the
thresholds collide" (`2026-10-22_c2b_feynman_g2.md`, CF7 note) and the geminal-shapes record's
diagnostic that at gamma = 1 the branch point i(gamma - b) sits on the momentum contour.
**Extends:** `2026-10-17_ring_extensions.md` (QG), `2026-10-22_correlated_triangle_g32.md`,
`2026-10-22_c2b_feynman.md`, `2026-10-22_c2b_feynman_g2.md`, `2026-10-22_geminal_shapes.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

With the exponential geminal e^(-gamma r) on the two edges at one density, the correlated triangle Y
has an annihilator of order nine at gamma = 3/2 and at gamma = 2, with composition factors of orders
1, 1, 1, 1, 2, 3: four algebraic pieces, one elliptic piece on the Legendre family with
z = (t^2 + 1 - gamma^2)/(2t), and one K3 piece. At gamma = 1 the minimal annihilator has order four,
R3 o R1, with no elliptic piece and one algebraic piece. The kite phi, the two densities and the two
geminals around the apex, has order six at gamma = 3/2 and order two at gamma = 1. Why?

## What I know before computing (disclosed)

- The kite in Feynman parameters (`debug/ring_ext_QG_common.py`, `phi_gem_z`): the density pair
  enters through M(u) = sqrt(1 + u(1-u) z), the geminal pair through N(v) = sqrt(gamma^2 + v(1-v) z),
  with s = M + N and k^2 = (u-v)^2 z. At gamma = 1, N is the same function of v as M is of u. So
  gamma = 1 is the member where the two bubbles have equal masses.
- The kite's genuine singular points at gamma = 3/2 (QG memo): z = -4, -9, -25/4, -1/4 and 0, that is
  minus the squares of the four pair masses 2, 2 gamma, 1 + gamma, 1 - gamma, plus the MUM point.
  At gamma = 1 these are -4, -4, -4, 0, 0: all four thresholds fall onto the density-pair threshold
  and the MUM point.
- Y's genuine singular points in t are the four pair masses: |1 - gamma|, 2, 1 + gamma, 2 gamma,
  measured at gamma = 3/2 and 2.
- The pair masses coincide two at a time at other gamma: |1 - gamma| = 2 at gamma = 3, and
  2 gamma = |1 - gamma| at gamma = 1/3. Only at gamma = 1 do all four fall onto two points, and only
  there is N the same function as M.
- The gamma = 1 Legendre family, z = t/2, is a genuine elliptic family, not isogenous to F_3's
  Clausen curve (geminal-shapes record). So the elliptic piece does not vanish at gamma = 1 by
  becoming the K3's curve.
- The elliptic piece's pseudo-threshold t = |1 - gamma| is a singular point of the operator and not of
  Y on the physical sheet, as for the kite, whose Taylor radius at gamma = 3/2 is 4, not 1/4.
- Nothing in this file has been computed.

## The mechanism, stated so it can fail

The drop at gamma = 1 is the equal-mass degeneration of the two bubbles, N = M, under which the four
pair thresholds collapse to two. It is not caused by two thresholds merely coinciding. Two readings
are therefore separated: (A) equal masses, predicting the order nine at every gamma other than 1,
including gamma = 3 where two thresholds coincide; (B) colliding thresholds, under which gamma = 3
should also lose pieces. I predict (A).

## Model

The record-g32 and cf_g2 pipeline for the exponential shape (`debug/cf_g2_*.py`, `debug/ga1_*.py`,
`debug/ring_ext_QG_*.py`): exact Taylor coefficients of the kite modulo primes, CRT, the kite's
operator in z and Y's operator in t by guessing with an exact check on float coefficients, the
composition series by exact pseudo-division, the Legendre map by the Stage E search, the K3 piece by
local data. Numerical check of each operator on the function by the QG quadrature at 30 digits or
more at four values of t, with the 1/200 control.

## Predictions

**MG1 (gamma = 1/2 is generic).** Y at gamma = 1/2 has an annihilator of order nine whose composition
series has orders 1, 1, 1, 1, 2, 3. Its genuine singular points in t > 0 are exactly {1/2, 1, 3/2, 2}.
The order-two piece is equivalent, by an exact order-one map rational in x = t^2, to the Legendre
operator of degree -1/2 in z = (t^2 + 3/4)/(2t) twisted by x^(1/2). The order-three piece has F_3's
local data twisted by t.
- Falsifier: order below nine with no elliptic piece, or an order-two piece not equivalent to that
  Legendre operator while the gamma = 2 control map is found by the same call.
- This separates the mechanism from the branch-point reading: for gamma < 1 the branch point
  i(gamma - 1) has not crossed the contour, and nothing else distinguishes gamma = 1/2 from
  gamma = 3/2.

**MG2 (main: gamma = 3 keeps the order; the mechanism is equal masses, not colliding thresholds).**
At gamma = 3 the pair masses are |1 - gamma| = 2, 2, 1 + gamma = 4 and 2 gamma = 6, so the
pseudo-threshold of the elliptic piece falls on the density-pair threshold t = 2. Y's annihilator at
gamma = 3 still has order nine with the series 1, 1, 1, 1, 2, 3; its genuine singular points in t > 0
are exactly {2, 4, 6}; the order-two piece is equivalent to the Legendre operator in
z = (t^2 - 8)/(2t) twisted by x^(1/2); the K3 piece is unchanged.
- Falsifier: an order below nine at gamma = 3, or no elliptic piece there. That would support
  reading (B) and refute the equal-mass mechanism as stated.
- Confidence: moderate. The collision at gamma = 3 puts an I2 fibre of the elliptic piece on top of
  the K3's pseudo-reflection at t = 2; a module can lose a factor at such a point. I predict it does
  not, because nothing in the integrand degenerates at gamma = 3.

**MG3 (the kite's singular points are the pair thresholds, as a rule).** At gamma = 1/2 the kite's
operator in z has order six and genuine singular points exactly {-4, -1, -9/4, -1/4, 0, infinity}; at
gamma = 3, order six and genuine singular points exactly {-4, -36, -16, 0, infinity}, with -4 now
carrying both the density pair and the pseudo-threshold.
- Falsifier: a genuine singular point of the kite's operator that is not minus the square of a pair
  mass, or a kite order other than six at either gamma.

## Novelty status

- MG3 at one gamma is known in the QG memo; the rule at new gamma is a CHECK.
- MG1: a CHECK of the one-function reading for gamma < 1.
- MG2: the mechanism statement. NOT FOUND in the literature for this integral. A hit is a FINDING
  only in that it pins the mechanism; the equal-mass degeneration of a two-loop self-energy is
  classical in the Feynman-integral literature and is credited as such at scoring.

## Controls

- Positive: the same pipeline at gamma = 2 reproduces the record's L9g, its series and the CF5 map.
- Positive: the kite at gamma = 1 reproduces the stored order-two operator.
- Negative: the 1/200 change in one coefficient of each accepted operator leaves a residual above
  1e-6 on the function.

## Caps stated in advance

One geminal shape, the exponential. Two new values of gamma. The mechanism is read from the integrand
and tested by its consequences at gamma = 1/2 and 3; no symbolic derivation of the drop from N = M is
attempted here. A feasibility STOP at either gamma is reported as a STOP and is not scored.

## Disposition (appended 2026-10-23, after computing; the predictions above are unchanged)

Pipeline `debug/mg_*.py`, data `debug/data/mg_*`, everything exact as strings in `debug/data/mg_pieces.json`.
Series of 800 coefficients at 600 digits at each gamma (Cauchy route at radius 1/2 for gamma = 1/2 and 2 for
gamma = 3), cross-checked coefficient by coefficient against an exact mod-p series at one prime (0 mismatches).

**MG1 HELD.** At gamma = 1/2 Y has an annihilator of order nine and x-degree 25; no operator of order six,
seven or eight to x-degree 63, 49, 39 (float), and none of order eight to degree 81 mod p. The composition
series is 1, 1, 1, 1, 2, 3 by exact pseudo-division. Genuine singular points in t > 0: exactly {1/2, 1, 3/2, 2}
(a degree-19 factor of the leading coefficient is apparent: exponents 0..7, 9, log-free). The order-two piece
is mapped exactly from the Legendre operator in z = (t^2 + 3/4)/(2t) twisted by x^(1/2), with
E = ((x - 1/4)(x - 9/4))^4 (x - 3/4) and numerator degrees 19, 21; the PM re-check on an independent code
path (`debug/pm_map_check.py`) confirms the identity and that a 1/1000 perturbation breaks it. The
order-three piece is mapped exactly from F_3's operator twisted by t.

**MG2 HELD (reading A).** At gamma = 3 the order stays nine (x-degree 23), the series stays 1, 1, 1, 1, 2, 3,
and the elliptic piece is present: an exact map from the Legendre operator in z = (t^2 - 8)/(2t) twisted by
x^(1/2), E = ((x - 4)(x - 16))^4 (x + 8), numerator degrees 19, 21, PM re-checked. Genuine singular points in
t > 0: exactly {2, 4, 6}. At x = 4 the local exponents are the union predicted: {-3/2, -1/2} from the K3 piece
and {-1 with a logarithm, 0, ..., 5} from the elliptic piece's I2 fibre. The K3 piece is equivalent to F_3's
operator twisted by t (exact map t (x - 4)^(-3) with numerator degree 18, against (x - 4)^(-4) and degree 20
at gamma = 1/2 and 2); its pseudo-reflection exponent at x = 4 is -9/2 rather than -11/2. So the collision of
the pseudo-threshold with the density-pair threshold merges local data but removes no factor. Colliding
thresholds do not cause the drop at gamma = 1; equal masses do.

**MG3 HELD.** The kite's operator in z has order six at both gamma (z-degree 23 and 21; no order five to
degree 122 mod p). Genuine singular points: {-4, -1, -9/4, -1/4, 0, infinity} at gamma = 1/2 and
{-4, -36, -16, 0, infinity} at gamma = 3, exactly minus the squared pair masses; at gamma = 3 the point -4
carries both the density pair (exponents -1 with a logarithm, 0, 0, 1, 2) and the pseudo-threshold (-1/2).
The operator applied to the kite by quadrature at 50 digits leaves at most 5e-41 at four values of z.

**Controls.** The same pipeline at gamma = 2 reproduces L9g (order nine, x-degree 25, residual 7e-364 on the
coefficients), its chain and factors, and the CF5 map with the same E and numerators. The kite at gamma = 1
reproduces the stored order-two operator. The 1/200 change in one coefficient leaves residuals between
1.8e-6 and 5.7e-4 on Y (the smallest at t = 5/2, gamma = 3; above the 1e-6 bar at every point) and between
5e-5 and 2e-2 on the kite, against at most 3e-24 and 5e-41 for the unperturbed operators.

**Residuals on the function.** L applied to Y by quadrature at 32 digits at t = 0.8, 1.3, 1.7, 2.5: at most
6e-27 (gamma = 1/2) and 3e-24 (gamma = 3), relative to the largest term.

**Deviations from the model.** Coefficients came from the float Cauchy route and were cross-checked mod p at
one prime, not reconstructed by CRT; the ratio-recognition tolerance was 150 digits, with the exact recheck
as the validation. A mod-p series bug (the last coefficient of a length-M run inexact) was found and fixed in
the new `mg_5d_modp_series.py`; the older ga1 route never saw it because it drops the last rows. No STOP.

**Credit.** The equal-mass degeneration of a two-loop self-energy is classical (the equal-mass sunrise);
what is new here is only that this is the mechanism for this integral's drop at gamma = 1, tested by its
consequences at gamma = 1/2 and 3. No symbolic derivation of the drop from N = M was attempted.

**Novelty at scoring.** MG1 and MG3: CHECK. MG2: FINDING (the mechanism pinned by a case that separates the
two readings).
