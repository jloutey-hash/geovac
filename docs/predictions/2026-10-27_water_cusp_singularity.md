# Pre-registered prediction: the bond-angle singularity of the three-centre integral is a cusp-cusp term, derived (WA5)

**Registered:** 2026-10-27, before any computation that scores a prediction in this file. PI-directed ("Yes", to
"pre-register a new test, WA5, and run it").
**Commit status:** committed on its own, before any computation in this file.
**Rests on:** `docs/predictions/2026-10-27_water_angle_singularities.md` (WA0-WA4; WA1 missed, WA2 missed for E2,
WA3 held) and the channel formula of `tests/test_paper59_water_bond_angle.py`.

## The question

WA1 predicted that the Legendre coefficients F_l of (XY|XZ) in x = cos(gamma) decay at the rate D</D> of the
outer-centre coalescence, and missed by 2-3 percent. Is the miss the physics or the fit? This record derives the type
and the size of the singularity, so that both the power and the amplitude are fixed in advance, and nothing about the
singular term is fitted.

## The derivation (done before registering; stated here so it can be checked)

- (XY|XZ) = int int rho1(r1) rho2(r2) / r12, with rho1 = chi_X chi_Y and rho2 = chi_X chi_Z, chi = N e^{-zeta r},
  N^2 = zeta^3 / pi. As a function of the angle it can be non-analytic only where a non-smooth part of rho1 meets a
  non-smooth part of rho2 through 1/r12. A cusp convolved with 1/r12 and a smooth function is smooth in the cusp's
  position. The cusps at X are at a fixed point. So the only angle-dependent non-analytic term comes from the cusp of
  rho1 at Y meeting the cusp of rho2 at Z, and it is a function of s = |Y - Z|.
- Near Y, rho1 = chi_X(Y) N_H (1 - zeta_H u + ...), u = |r - Y|. The non-smooth part is c_Y u with
  c_Y = -zeta_H N_H chi_X(Y), and likewise c_Z at Z. Terms u times a polynomial in (r - Y) give s^6 (s . v)-type
  pieces, which are analytic; the next non-analytic piece is s^9.
- In three dimensions FT[|r|] = -8 pi / k^4 and FT[1/r] = 4 pi / k^2, so the cusp-cusp term has transform
  256 pi^3 / k^10. FT[r^7] = 161280 pi / k^10. So the leading non-analytic term is
  c_Y c_Z (pi^2 / 630) s^7 = C s^7, with C = zeta_H^5 zeta_X^3 e^{-zeta_X (D1 + D2)} / 630.
- *Check already made (disclosed).* For two 1s densities e^{-2r}/pi, the exact two-centre Coulomb integral
  J(R) = 1/R - e^{-2R}(1/R + 11/8 + 3R/4 + R^2/6) has the series 5/8 - R^2/12 + R^4/60 - R^6/126 + (2/315) R^7 + ...:
  no R, R^3 or R^5 term, and the R^7 coefficient (2/pi)^2 pi^2/630 = 2/315 exactly.
- With s^2 = 2 D1 D2 (x* - x), x* = (D1^2 + D2^2)/(2 D1 D2), the singular term is C (2 D1 D2)^{7/2} (x* - x)^{7/2}.
  Its Legendre coefficients L_l fall like rho^l l^{-4} for D1 != D2 (rho = D</D>) and like l^{-8} for D1 = D2. The s^9
  term changes F_l / L_l by a relative O(1/l) (D1 != D2) or O(1/l^2) (D1 = D2).

## What I know before computing (disclosed)

- WA1/WA2 results and the unregistered follow-up of the WA record (a free 1/l term brings r within 0.9 and 0.6
  percent; symmetric p = 8.05).
- *Method diagnosis on a synthetic function (not water data).* On the exact Legendre coefficients of (x* - x)^{7/2},
  the registered WA fit (free p) gives r/rho = 1.048 at both asymmetric pairs, with p near 6. With p fixed at 4 and
  terms 1/l, 1/l^2 it gives r/rho = 0.99919 and 0.99951. The registered fit could not have met 1 percent on a function
  whose singularity is exactly at coalescence. Script: the session scratchpad `darboux_ctrl.py`, re-created as
  `debug/wa5_cusp.py` mode `control`.
- I have not computed C, L_l, or any fit with p fixed on water coefficients.

## Setup

- Coefficients F_l: the channel formula `F_coeffs` of `tests/test_paper59_water_bond_angle.py`, copied unchanged into a
  new driver `debug/wa5_cusp.py`. Grid level 2 for values, level 1 vs level 2 for the noise floor. l from 0 to 60.
  Usable l: l >= 15 with |F_l| at least 100 times the floor, as in the WA record.
- Cases: the six of the WA record. Geometries (1.809, 1.809), (1.809, 2.5), (1.809, 1.0). Sets E1 (zeta_X = zeta_H = 1)
  and E2 (zeta_X = 2.25, zeta_H = 1.24), with zeta_X on X as in the test.
- L_l: Legendre coefficients of C (2 D1 D2)^{7/2} (x* - x)^{7/2}, by mpmath quadrature at 40 digits.

## Predictions

**WA5-0 (controls; run FIRST; WA5a-c void if any fails).**
- (a) The channel sum reproduces REFEREE_NG16 = 0.204941652260 to 1e-9 at the test geometry.
- (b) On synthetic coefficients of C' (x* - x)^{7/2} + C'' (x* - x)^{9/2}, with C'' / C' = 3, the WA5a fit gives r within
  0.1 percent of rho at both asymmetric pairs, and the WA5b fit gives A0 within 1e-2 of 1 at all three geometries.

**WA5a (the rate, power fixed by the derivation).** At the four asymmetric cases, the fit
log|F_l| + 4 log l = a + l log r + b/l + c/l^2 over the usable l gives r within 0.3 percent of D</D>.
- Falsifier: off by 0.3 percent or more in any of the four.
- Probability 0.75.

**WA5b (the amplitude, nothing fitted in the singular term).** At all six cases, the fit
F_l / L_l = A0 + A1/l + A2/l^2 over the usable l gives |A0 - 1| < 0.03.
- Falsifier: |A0 - 1| at or above 0.03 in any case. A miss with WA5a holding would locate the singularity but not its
  size (a normalization or a missed cusp term).
- Probability 0.65.

**WA5c (symmetric power).** At D1 = D2, both sets, the fit log|F_l| = a - p log l + b/l (r fixed at 1) gives p within
0.15 of 8.
- Falsifier: |p - 8| at or above 0.15 in either set.
- Probability 0.7.

## Decision gate

- **GO:** WA5-0 passes and WA5a, WA5b, WA5c hold. Then the bond-angle singularity is derived: the H-H cusp-cusp term
  C s^7, at the coalescence of the outer centres, with its size known in closed form. WA1's miss is the fit model.
- **BORDERLINE:** WA5a and WA5c hold, WA5b fails. Location and type derived, size not.
- **STOP:** WA5-0 fails, or WA5a fails.

## Caps

- At most 1 h wall, 8 processes. No change to `geovac/` or `tests/` during the run.
- A case with fewer than 10 usable l is void and reported.

## Novelty status (before computing)

- Darboux asymptotics of Legendre coefficients: KNOWN, standard.
- Cusp-generated non-analytic terms in two-centre integrals (odd powers of R): the R^7 series term is a property of the
  known closed form; the cusp-cusp reading of it and its use for the three-centre angle: NOT SEARCHED. Treat as OPEN.

## Disposition

(Appended after computing. The predictions above are not edited.)

### Disposition (appended 2026-10-27 by the PM; computed in the main session)

Driver `debug/wa5_cusp.py` (modes `control`, `cases`); logs `debug/data/wa5_{control,cases}.log`; data
`debug/data/wa5_cusp.json`.

| Item | Result | Verdict |
|:--|:--|:--|
| WA5-0(a) | channel sum 0.204941652531, diff 2.7e-10 | PASS |
| WA5-0(b) | synthetic s^7 + 3 s^9: A0 = 0.988, 0.975, 0.973 (need 1e-2); r/rho = 1.0115, 1.0551 (need 1e-3) | **FAIL** |
| WA5a-c | void as registered | VOID |

**Gate: STOP as registered** (a control failed).

**Why the control failed (diagnosed after, not registered).** The quadrature is not the cause: the synthetic
coefficients match the closed form for (1 - x)^a to all printed digits. The cause is the control's own size. The
s^9 / s^7 coefficient ratio behaves like 2(a+1)^2/l^2, about 40/l^2 for a = 7/2, so with C''/C' = 3 the correction is
about 50 percent at l = 15, and a three-term series in 1/l cannot absorb it to 1 percent. The registered fit models are
too weak for a subleading term that large. Whether water's own subleading term is that large was not known in advance.

**Water cases, reported only (they score nothing).**

| Case | usable l | r / (D</D>) (p fixed 4) | p (D1 = D2) | A0 | F_l / L_l at the two ends |
|:--|:--|:--|:--|:--|:--|
| E1 (1.809, 1.809) | 15-60 | | 8.017 | 1.0036 | 0.850, 0.990 |
| E1 (1.809, 2.5) | 15-45 | 0.99981 | | 0.9899 | 0.771, 0.940 |
| E1 (1.809, 1.0) | 15-29 | 1.00004 | | 0.9518 | 0.774, 0.875 |
| E2 (1.809, 1.809) | 15-60 | | 8.038 | 1.0068 | 0.814, 0.987 |
| E2 (1.809, 2.5) | 15-45 | 0.99976 | | 1.0670 | 0.804, 1.004 |
| E2 (1.809, 1.0) | 15-29 | 0.99996 | | 1.0145 | 0.744, 0.888 |

- With the power fixed by the derivation, the rate lands within 0.03 percent of D</D> in all four asymmetric cases
  (WA1 had missed by 2-3 percent). The symmetric power is 8.02 and 8.04.
- F_l / L_l rises toward 1 in every case and is within 1.5 percent of it at l = 60 for symmetric water. The fitted A0
  is within 7 percent of 1, but the control shows that fit is not reliable to better than a few percent here.
- Reading (not established, since unscored): the bond-angle singularity is the H-H cusp-cusp term C s^7, at the
  outer-centre coalescence, with the size C = zeta_H^5 zeta_X^3 e^{-zeta_X(D1+D2)}/630. A scored test needs either the
  s^9 coefficient derived too, so that nothing is fitted, or a control sized to the subleading term actually present.
