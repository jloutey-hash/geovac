# Pre-registered prediction: the differential equation of the three-centre integral in the bond angle (AO)

**Registered:** 2026-10-28, before any Legendre coefficient of a scored case is computed and before the operator
search is run on any three-centre data. PI-directed ("My mind keeps wondering back to the three-centre closed form ...
Maybe we can work on that?", then "Yep" to writing this registration).
**Commit status:** committed on its own, before any computation it scores.

## Why

Paper 59 proves the radial core N(D) has no closed form through any factorization (`sec:obstruction`) and leaves open
whether the integrated three-centre integral has one. Every test so far took the integral as a function of distance.
The water records (WA5-WA8, WX) took it as a function of the bond angle, x = cos(gamma), and mapped its singular point
there. The integral f(x) = (XY|XZ) is D-finite in x (the integrand is holonomic in the positions and in x, and
definite integration keeps that), so it obeys a linear ODE in x with polynomial coefficients. Nobody has written it
down. If it factors, that is a route to a closed form; if it does not, it is the first irreducibility statement for
the integrated object rather than its radial core. Either outcome moves the open problem.

## Method

**Data.** f(x) = sum_L F_L P_L(x), X at the origin, 1s Slater functions, |XY| = D1, |XZ| = D2. Position-space route
`debug/ao_radial.py`: F_L = (4 pi/(2L+1))^2 R_L, the radial integral of the two L-th multipoles of the densities
phi_X phi_H, each in closed form through i_L, k_L (the Barnett-Coulson / Lowdin alpha function, prior art); radial
integrals by Chebyshev on pieces split at D1, D2 and on a geometric tail, n + 2L nodes per piece, 110 digits.
Independent of the momentum route `debug/ao_engine.py` (same channel formula as `tests/test_paper59_water_bond_angle.py`).

**Operator search.** `debug/ao_guess.py`, `debug/ao_scan.py`. Ansatz sum_{i<=R, j<=Q} a_ij x^j T^i f = 0 with
T = (1 - x^2) d/dx; x and T act locally on Legendre coefficients, so each Legendre coefficient of the left side is a
finite combination of the F_L. Null vector by SVD on the fitted rows, the last 20 rows held out. Budget: R = 1..5,
Q = 1..30, (R+1)(Q+1) <= n_fit - 10, in order of (R+1)(Q+1), then R; stop at the first FOUND.
**FOUND** = gap s2/smin >= 1e15, smin <= 10 eps, held-out residual <= 1e3 eps, with eps = 1e-50 fixed for every
case: the value the controls below were run with (`ao_scan.py` argument). Control (c) requires the data itself to be
better than 1e-60, so eps is a threshold, not the data precision. Working precision 110 digits.

**Analysis.** `debug/ao_analyze.py`: standard form, finite singular points (roots of the leading coefficient, clustered,
multiplicity from the cluster), local exponents by the indicial polynomial, the Fuchs test at infinity.

## Cases

| Case | role | D1 | D2 | zX | zH | x* | D</D> |
|:--|:--|:--|:--|:--|:--|:--|:--|
| probe | unscored, measurement only | 2 | 1.5811388300841898 | 1 | 1 | 1.0277402 | 0.7905694 |
| nb | unscored neighbour of S1 | 1.85 | 2.65 | 1.3 | 1.0 | 1.0652728 | 0.6981132 |
| S1 | scored | 1.8 | 2.6 | 1.3 | 1.0 | 1.0683761 | 0.6923077 |
| S2 | scored | 1.5 | 2.4 | 1.0 | 1.2 | 1.1125 | 0.625 |

F_L for L = 0..140, at n = 100 and n = 140.

## Controls measured before registering

All on unscored cases or synthetic data; no three-centre coefficients have been put through the operator search.

| Control | Measured | Bound set |
|:--|:--|:--|
| (a) position vs momentum route, L = 0..5 (nb; momentum level 3, `data/ao_nb_mom3.json`) | max 4.3e-10 | 1e-8 |
| (b) probe sum vs the pairwise limit 0.204941652531 | 3.5e-13 | 1e-11 |
| (c) n = 100 vs n = 140, L = 0..140 (probe; nb) | 2.7e-82; 5.4e-74 (worst at L = 140) | 1e-60 |
| (d) positive control f = e^x (1.3 - x)^(7/2) + x e^(-2x) (1.3 - x)^(9/2), L = 0..140, 110 digits, eps 1e-50 | FOUND at (R, Q) = (2, 8): smin 1.2e-68, gap 2.9e51, held 7.0e-60 | must reproduce |
| (d) decoy Gamma(x + 2) (1.3 - x)^(7/2), same budget, all 107 candidates | no FOUND; worst gap 7.9e3 (smin falls to 6e-79, so the gap, not smin, discriminates) | must reproduce |
| analysis on the positive control (2, 8) | exponents at 1.3: 3.5, 4.5; at +-1: 0, 1; two other points 0, 2; infinity irregular | as found |

Earlier, the plain recurrence ansatz (r <= 30, degree <= 4) found nothing on the control at 1e-40: its residuals fell
smoothly, the curve-fit signature. That is why the ansatz is the local-operator form above. The analysis code had a
double-root bug (it read 1.3 as a simple root, exponents 0 and 4.5), fixed before registering.
Timing: F_L to L = 140 at n = 140 about 5.5 min per case on 15 cores; a full decoy scan about 15 min.

## Predictions

- **AO0 (controls; run FIRST).**
  - (a) Engine cross-check: on S1, the position route and the momentum route (`ao_engine.py`, level 3) agree for
    L = 0..5 to 1e-8 relative (measured 4.3e-10 on nb).
  - (b) The probe sum at x = 1/(2 sqrt 2.5) reproduces the pairwise limit 0.204941652531 (CHANGELOG, pairwise
    three-centre test) to 1e-11 (measured 3.5e-13).
  - (c) Per scored case, eps = max_L |F_L(n=100)/F_L(n=140) - 1| over L = 0..140 must be <= 1e-60 (measured <= 5.4e-74). A case failing is void.
  - (d) The positive control (FOUND at (R, Q) = (2, 8) on the full-size data) and the decoy (no FOUND in the whole
    budget; worst gap below 1e15) are re-run unchanged and must reproduce.
- **AO1 (existence within budget; the gate).** For each live scored case the scan returns FOUND. Probability 0.55
  (the order could exceed the budget; the kite's minimal order in its own variable was 13).
- **AO2 (the Landau map; only if AO1).** x* is a root of the leading coefficient to 1e-20, and every other finite
  singular point, x = +-1 included, has distinct non-negative integer exponents (to 1e-10). Falsifier: any other
  finite singular point with a non-integer or negative exponent, or irregular. This is the claim that the integral's
  only singularity in the finite angle plane is the complex angle where the two hydrogens touch. Probability 0.7.
- **AO3 (only if AO1).** 7/2 is a local exponent at x*, to 1e-10. Probability 0.85.
- **AO4 (only if AO1).** Infinity is an irregular singular point (the Fuchs test fails). Probability 0.85.
- **AO5 (only if AO1 in both cases).** S1 and S2 give the same minimal (R, Q). Probability 0.6.

**Reported, not scored:** R and Q; all exponents at x*; whether the operator has a right factor of order one,
tested afterwards by a method appended here before it is run.

**Gate.** STOP if a control fails or AO1 fails in both cases. GO if AO1 holds in at least one live case and AO2-AO4
hold there. BORDERLINE otherwise.

## Novelty

The D-finiteness of multicentre integrals in a geometric parameter is general holonomic theory (Zeilberger; Chyzak).
The Legendre-coefficient multipole route is Barnett-Coulson / Lowdin. An explicit ODE in the bond angle for a
three-centre two-electron Slater integral: NOT SEARCHED, treat as OPEN.

## Caps

1s functions, two scored geometries, budget R <= 5, Q <= 30; NOT FOUND within budget says nothing about the order
beyond it.
