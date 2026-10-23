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
