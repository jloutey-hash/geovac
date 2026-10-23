# Pre-registered prediction: the cusp-correct geminal at a second gamma (GA3)

**Registered:** 2026-10-23, before any computation in this file (PI direction: "let's start knocking
out the open items").
**Prompted by:** GA1 held at gamma = 3/2 (`2026-10-22_geminal_shapes.md`): the cusp-correct geminal
f = r e^(-gamma r) on both edges brings no new type. Paper 59's cap now reads "of the LiH engine's two
shapes the cusp-correct one is tested at gamma = 3/2 only".
**Extends:** `2026-10-22_geminal_shapes.md` (GA1), `2026-10-22_c2b_feynman_g2.md` (CF5-CF7).
**Commit status:** committed on its own, before any computation in this file.

## The question

Does the GA1 result hold at a second gamma, so that it is a rule of the shape and not a feature of
gamma = 3/2?

## What I know before computing (disclosed)

- GA1 at gamma = 3/2: order nine, series 1, 1, 1, 1, 2, 3; the order-two piece of x-degree 36 is
  equivalent to the Legendre operator in z = (t^2 + 1 - gamma^2)/(2t) twisted by x^(1/2); the
  order-three piece is F_3's operator twisted by t; genuine singular points {1/2, 2, 5/2, 3}; the
  algebraic pieces carry exponents -5/2 at t = 2 and -11/2, -13/2 at t = 3.
- The exponential shape at gamma = 2 (CF5-CF7): order nine, the same series, singular points
  {1, 2, 3, 4}, Legendre in z = (t^2 - 3)/(2t).
- The linexp kernel is the second gamma-derivative of the Yukawa kite, so its singular points are the
  exponential case's, with higher exponents at t = 2 gamma.
- Nothing in this file has been computed.

## Model

The GA1 pipeline (`debug/ga1_*.py`) at gamma = 2 with the linexp shape on both edges: kernel by the
stage-0 derivation, exact mod-p Taylor coefficients, CRT, the operator in x checked on float
coefficients and on Y_lin by quadrature at 32 digits with the 1/200 control, the series by exact
pseudo-division, the Legendre map by the Stage E search with the exponential gamma = 2 map as control.

## Predictions

**GA3 (the cusp-correct shape at gamma = 2 brings no new type).** Y_lin at gamma = 2 has an
annihilator of order nine with composition factors of orders 1, 1, 1, 1, 2, 3. Its genuine singular
points in t > 0 are exactly {1, 2, 3, 4}. The order-two piece is equivalent, by an exact order-one
map rational in x, to the Legendre operator of degree -1/2 in z = (t^2 - 3)/(2t) twisted by x^(1/2),
and it is a different operator from the exponential case's order-two piece at gamma = 2. The
order-three piece has F_3's local data twisted by t. The algebraic pieces at t = 4 carry exponents
-11/2 and -13/2, as at gamma = 3/2 they did at t = 3.
- Falsifier: a composition factor of order two or more that is neither the Legendre family nor
  F_3's operator twisted by t, exhibited exactly, with the controls met; or a genuine singular point
  outside {1, 2, 3, 4}.
- Not scored: a feasibility STOP.

## Novelty status

GA1 at a second gamma. A hit is a CHECK of a FINDING already recorded; it widens Paper 59's cap from
"gamma = 3/2 only" to two values of gamma for this shape.

## Controls

- Positive: the pipeline at the exponential shape and gamma = 2 reproduces L9g and the CF5 map.
- Negative: the 1/200 change in one coefficient of the annihilator leaves a residual above 1e-6
  on Y_lin.

## Caps stated in advance

The cusp-correct shape on both edges, at gamma = 2 only in this file. The operator is guessed from
exact series and checked on the function; minimality is bounded as in record g32, not proved.
