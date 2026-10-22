# Pre-registered prediction: the hyperfine x^3 ln x constant, third test (X3-1c)

**Registered:** 2026-10-22, before any computation in this file (PI-directed: "Let's do those next steps",
item one of the 2026-10-21 overview).
**Prompted by:** X3-1b (`2026-10-21_x3_constant.md`), scored UNCLEAR: its column sets A and B lacked an
x^4 ln^2 x term that the data carry. Its disposition says a decisive test is a new registration with the
log-squared columns, on data not yet fitted. This is that registration.
**Extends:** `2026-10-21_x3_constant.md`, `2026-10-19_eft_checks.md`, `2026-10-19_zemach_door.md`;
Paper 23 `sec:finite_nucleus`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Is the profile-independent constant c of the hyperfine x^3 ln x term, defined as in X3-1, equal to -64/63,
for n = 1, 2 and 3?

## What I know before computing (disclosed)

- X3-1b's data (window 1e-6 .. 1e-3, 24 points, n = 1, 2): set C gave c within 1.5e-9 to 3.6e-9 of
  -64/63 in all six cases. An unscored PM refit of the same data with x^6 and x^6 ln^2 x added gave
  1.7e-12 to 3.6e-12. The column choice below is informed by that refit, so these data cannot be reused.
- The EFT reading behind -64/63 is the memo author's and was not re-derived by the PM.
- n = 3 has not been fitted for this term. X4-1 found n-dependence at fourth order in other
  coefficients; the x^3 ln x term is third order.
- Nothing in this file has been computed.

## Model

As X3-1b: hydrogen-like, Z = 1, atomic units, uniformly charged sphere of radius R = x; bound states and
hyperfine factors from `debug/zemach_door.py`, unmodified, at mp.dps = 120, tolerance 10^-125.
Cases: point, ball and shell magnetization, n = 1, 2, 3 (nine cases).
**Fresh data:** window x = 3e-7 .. 3e-4, 30 log-spaced points. No point coincides with the X3-1b grid.

Column sets for H(x), written as (k, j) for x^k ln^j x:
- D: (1,0) (2,1) (2,0) (3,1) (3,0) (4,2) (4,1) (4,0) (5,2) (5,1) (5,0) (6,2) (6,1) (6,0);
- E: D plus (7,1) (7,0);
- F: D plus (3,2) (6,3). The (3,2) column tests that no x^3 ln^2 x term is present.

c = (x^3 ln x coefficient) - (4/5) x (x coefficient), the X3-1 definition.

## Predictions

**X3-1c.** In all nine cases and all three column sets, c = -64/63.
- Spread: the largest difference of c among D, E and F, over the nine cases, relative to 64/63.
- Pass: spread below 1e-9 and |c + 64/63| < 1e-9 x 64/63 in every case and set.
- Falsifier: spread below 1e-9 and |c + 64/63| > 1e-8 x 64/63 in any case.
- Anything else: UNCLEAR, reported as such.
- If the falsifier fires only at n = 3, that is reported as n-dependence of c, with the n = 1, 2 verdict
  stated separately.

## Novelty status

- X3-1c: NOT FOUND, as for X3-1 and X3-1b; the PM re-checks at scoring. A hit is a **FINDING** of low
  weight: the value follows from the EFT rule, which is known.

## Controls

- Positive: the x coefficient reproduces the Zemach coefficients (-3/2, -72/35, -12/5 for point, ball,
  shell) to 1e-12 in every case.
- Positive: synthetic data built from the fitted series of case (ball, n = 1, set D) with c shifted by
  1e-6 and noise at 1e-100 recovers the shifted c to 1e-12.
- Consistency: the nine cases agree with each other to 1e-10 within each set.

## Caps stated in advance

A fit, not a derivation. Three values of n, three magnetization profiles, one charge model.
