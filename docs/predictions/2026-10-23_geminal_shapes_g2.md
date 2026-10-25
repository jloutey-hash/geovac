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

## Disposition (appended 2026-10-23, after computing; the prediction above is unchanged)

Pipeline `debug/ga3_*.py` (the GA1 stages at gamma = 2), data `debug/data/ga3_*`, everything exact as strings
in `debug/data/ga3_pieces.json`. Coefficients: exact mod-p Taylor coefficients of the kite (520 per prime, 12
primes), CRT and rational reconstruction stable from 8 primes on; an independent 600-digit float series of 800
coefficients agrees with the mod-p series on every recognisable coefficient (0 mismatches in 50).

**GA3 HELD.** Y_lin at gamma = 2 has an annihilator of order nine and x-degree 32 with integer coefficients,
exact modulo all nine primes used and leaving 3e-364 on the 800 float coefficients. Minimality bound, exact
mod p: no order six to x-degree 72, seven to 65, eight to 57. The composition series is 1, 1, 1, 1, 2, 3 by
exact pseudo-division, with x-degrees 2, 5, 9, 13, 36, 71, the same degrees as at gamma = 3/2. Genuine
singular points in t > 0: exactly {1, 2, 3, 4} (a degree-26 factor of the leading coefficient is apparent:
exponents 0..7, 9, log-free at two primes). The order-two piece is mapped exactly from the Legendre operator
in z = (t^2 - 3)/(2t) twisted by x^(1/2), with E = ((x - 1)(x - 9))^5 (x + 3) and numerator degrees 23, 25
(48 smaller boxes empty); the PM re-check on an independent code path (`debug/pm_map_check.py`) confirms the
identity and that a 1/1000 perturbation breaks it. It is a different operator from the exponential case's
C2b at gamma = 2: not proportional, x-degree 36 against 30, the same genuine points x = 1, 9. The order-three
piece is mapped exactly from F_3's operator twisted by t (A = t (x - 4)^(-4) times numerators of degree 28,
rebuilt over 12 primes, exact), with the MUM point at t = 0 and the pseudo-reflection at x = 4. The algebraic
pieces are rational, x^0 (x - 4)^(-5/2) times a quartic over x, and x^(3/2) (x - 16)^(-11/2) and
x^(5/2) (x - 16)^(-13/2) times polynomials: exponents -5/2 at t = 2 and -11/2, -13/2 at t = 4, as predicted.

**Controls.** The same map search on the exponential C2b at gamma = 2 recovers the CF5 map (same E, degree
20, equal up to scale). The mod-p guesser on the exponential series reproduces L9g (its (9, 25) null vector
equals L9g mod p; nullity 0 at (8, 39) and (9, 24)). The generic peeler on L9g gives 1, 1, 1, 1, 2, 3 with F
as its order-three factor, and the F_3 map control reproduces CF6's box. Stage-1 checks at gamma = 2 (loop
against direct 3-D, kite route A against B, the exponential kite against QG's) agree to 1.3e-12.

**Residuals on the function.** L applied to Y_lin by quadrature at 32 digits at t = 0.8, 1.3, 1.7, 2.5:
1.4e-27, 2.0e-27, 5.9e-26, 2.1e-24 of the largest term. The 1/200 change in the largest theta^0 coefficient,
the GA1 form of the control, leaves 3.5e-2, 2.2e-3, 2.0e-4 and 4.4e-8: the last point sits below the 1e-6
bar, as GA1's 2.0e-7 did at gamma = 3/2, because that coefficient's monomial is small at t = 5/2. Perturbing
instead the coefficient whose monomial is largest at each t leaves 7.4e-2, 8.7e-3, 1.5e-2, 1.5e-2
(`debug/ga3_9b_onY_control.py`). Both are reported; the control as registered is met at three of the four
points in its GA1 form and at all four in the dominant-monomial form.

**Deviations from the model.** The split-search degree cap had to be raised from 80 to 130 for the order-five
quotient to split (its K3 left factor has degree 85, as at gamma = 3/2). 520 coefficients per prime over 12
primes against GA1's 420 over 8. The stage-0 kernel is symbolic in gamma and was not re-derived; the stage-1
checks at gamma = 2 stand in. Float searches at tolerance 330 show edge fits at orders seven and eight some 30
orders above the floor; they were disposed of exactly by the mod-p scan. No STOP.

**Novelty at scoring.** CHECK, as registered. Paper 59's cap for the cusp-correct shape widens from
gamma = 3/2 only to gamma = 3/2 and 2.
