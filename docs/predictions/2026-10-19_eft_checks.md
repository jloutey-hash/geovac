# Pre-registered prediction: two EFT checks on the finite-nucleus doors (X4, X3)

**Registered:** 2026-10-19, before any computation in this file (PI-directed: "I want to do all this",
after the proposal to check the multi-focal wall against effective field theory).
**Prompted by:** `debug/eft_multifocal_check_memo.md`, which found the doors' mechanism to be the NRQED
matching statement and named two cheap tests. The PM verified its sources in the primary texts:
Pachucki-Patkos-Yerokhin, arXiv:2502.09455 (PPY25) Eq. (30) and (77); Karshenboim, arXiv:physics/0102085
Eq. (26) and hep-ph/9608484 Eq. (8); Hill-Paz, arXiv:1103.4617.
**Extends:** `2026-10-18_finite_nucleus.md`, `2026-10-18_nuclear_profiles.md`, `2026-10-19_zemach_door.md`;
Paper 23 `sec:finite_nucleus`; seam rows J14, J16.
**Commit status:** committed on its own, before any computation in this file.

## The questions

1. Paper 23 says the s-wave quantum defect of a uniformly charged nucleus is the same at every n. That is
   shown through x^3 and at threshold. EFT says the n-dependence starts at x^4, from the iterated contact
   term and from the <r^4> operator. Does our own bound-state code see it, with the size PPY25 implies?
2. The hyperfine door found an x^3 ln x term whose profile-independent constant, -1.0147, was not
   identified. EFT says the log comes from iterating lower-order contacts: the charge contact times the
   Zemach contact, plus the charge's Friar contact. That gives -2 delta_3 = -64/63. Does a tighter fit
   confirm it?

## What I know before computing (disclosed)

- From the memo, by its author's algebra, checked by the PM once (the -1/n consistency check): converting
  PPY25 Eq. (30) to a quantum defect with Delta E = -mu/n^3 - (3/2) mu^2/n^4 (Z = 1) gives
  mu_4(n) = mu_4(inf) - delta_0^2 [psi(n) - ln n + 1/(2n)] - <r^4>/(15 n^2), with delta_0 = -2/5 and
  <r^4> = 3R^4/7 for the uniform sphere. The PM computed by hand:
  - f(n) = psi(n) - ln n + 1/(2n): f(1) = 1/2 - gamma_E = 0.0772157 with a minus sign, i.e.
    f(1) = -0.0772157; f(2) = 3/4 - gamma_E - ln 2 = -0.0203628; f(3) = 3/2 + 1/6 - gamma_E - ln 3
    = -0.0104... (to be evaluated exactly by the script).
  - Predicted mu_4(1) - mu_4(2) = -(4/25)[f(1) - f(2)] - (1/35)(1 - 1/4) = -0.0123321.
- From the zemach-door record: the x^3 ln x coefficient was measured as (4/5) x (Zemach coefficient) +
  c with c = -1.0147, the same for point, ball and shell magnetization to 2e-5 within one column set and
  to 1e-2 across column sets. 64/63 = 1.015873. The two differ by 1.2e-3 relative, inside the spread.
  This is a numerical match found after the fact; this file is its test.
- Nothing in this file has been computed.

## Model

Hydrogen-like, Z = 1, atomic units, uniformly charged sphere of radius R = x. Bound states from
`debug/zemach_door.py` (interior power series, exterior Whittaker function, secant on nu = n - mu), at
mp.dps = 60, tolerance 10^-65. Hyperfine factors for point, ball and shell magnetization as defined there.

## Predictions

**X4-1 (n-dependence at fourth order).** D_12(x) = mu(1) - mu(2) has no x^2, x^3 or x^4 ln x term, and
its x^4 coefficient is -(4/25)[f(1) - f(2)] - (1/35)(3/4) = -0.0123321 (exact form evaluated in the
script). Likewise D_13(x) = mu(1) - mu(3), coefficient -(4/25)[f(1) - f(3)] - (1/35)(8/9).
- Fit: x = 1e-4 .. 3e-3, 14 log-spaced points; D/x^4 fitted with columns 1, x ln x, x, x^2 ln x, x^2.
  Stability: drop the last column; report both.
- Pass: the constant column agrees with the prediction to 1e-3 relative for both D_12 and D_13, and a
  fit with an added ln x column gives that column below 1e-3 of the constant.
- Falsifier: disagreement beyond 1e-3 relative, or a ln x column at x^4 above that bound. If the
  falsifier fires, the PM's next step is to recheck the memo's conversion, not the paper.
- Novelty status: KNOWN (PPY25 Eq. 30). A hit checks our code and the conversion; it is a **CHECK**. It
  also decides the scope of Paper 23's "the same at every n".

**X3-1 (the hyperfine x^3 ln x constant).** For point, ball and shell magnetization, at n = 1 and 2,
the x^3 ln x coefficient minus (4/5) x (the x coefficient) equals -64/63.
- Fit: H(x) at x = 1e-5 .. 2e-3, 16 log-spaced points, dps 60; columns x, x^2 ln x, x^2, x^3 ln x, x^3,
  x^4 ln x, x^4, x^5 ln x. Stability: drop the last column.
- Pass: |c + 64/63| < 1e-5 x 64/63 for all six (profile, n) cases in both fits.
- Falsifier: any case outside 1e-4 relative. Between 1e-5 and 1e-4: UNCLEAR, reported as such.
- Novelty status: NOT FOUND. The x^3 ln x term is not in K97, K01, KI02, EGS01 Sec. XXXIV, FP05 or
  Volotka 2006 App. A (the memo's search; the PM re-checks K01 and K97 at scoring). A hit is a
  **FINDING**, of low weight: it follows from the EFT rule, which is known.

## Controls

- The x coefficient of H reproduces the zemach-door record's Zemach coefficients to 1e-8.
- The x^2 ln x coefficient of H reproduces 4/5 to 1e-6.
- mu(1) at x = 1e-3 reproduces the finite_nucleus record's threshold series -(2/5)x^2 + (32/63)x^3 to
  the x^4 level (difference below 1e-11).

## Result

(To be appended after scoring; the predictions above are never edited.)

## Disposition (scored 2026-10-19)

Script `debug/eft_checks.py`; data `debug/data/eft_checks_x4.json`, `debug/data/eft_checks_x3.json` and
their logs. One slip in the disclosure above: the hand value of f(3) is -0.00916, not -0.0104; the script
evaluated it exactly, as the file said it would, and nothing scored depends on the hand value.

**Controls.** All pass. mu(1) at x = 1e-3 minus the threshold series through x^3 is 6.6e-13, the x^4 level.
The x coefficients of H reproduce -3/2, -72/35, -12/5. The x^2 ln x coefficients are 0.80000005 (6.5e-8
relative).

**X4-1: HELD.** D_12: predicted x^4 coefficient -0.0123321203; fitted -0.0123321204, relative error 5.4e-9
(5 columns) and 1.1e-7 (4 columns); an added ln x column is 1.2e-8 of the constant. D_13: predicted
-0.0145081249; relative error 5.3e-9 and 1.1e-7; ln x column 1.2e-8. The s-wave defect depends on n at
x^4, as PPY25 Eq. (30) implies. Class: **CHECK**. Consequence: Paper 23's "the same at every n" is scoped
to x^3; the fourth-order term is now Eq. `eq:defect_n4` there.

**X3-1: FAILED at its registered gate.** For all six (profile, n) cases the constant is -1.01486 (8
columns), the same across profiles and n to 2e-5, against -64/63 = -1.015873: relative difference 1.0e-3,
outside the 1e-4 falsifier bound. The 7-column fit gives -1.0129 to -1.0130, so the fit's own truncation
spread (2e-3) exceeds the discrepancy: the data neither confirm nor exclude -64/63. The registered
test failed; this fit does not decide -64/63, and the constant stays unidentified. Class: **MISS**. A decisive test needs an exact derivation
or a fit whose truncation spread is well below 1e-3; either would be a new pre-registration.
