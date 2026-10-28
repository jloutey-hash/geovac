# Pre-registered prediction: the zero-contact bond-angle terms on four fresh cases (WA8)

**Registered:** 2026-10-27, before any F_l in this file is computed. PI-directed ("let's do that fresh registration").
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** `docs/predictions/2026-10-27_water_zero_contact.md` (WA7: STOP on a sign clause the registration stated
wrongly; its size and remainder clauses held, reported only).

## Why a fresh record

WA7's numbers are known, so they cannot be re-scored. WA8 asks the same question on four cases never computed, with the
sign clause stated correctly and the bounds set by a stated rule, not by WA7's results.

## The derivation (done before registering)

- Script `debug/wa8_local.py`, calling `debug/wa7_local.py` unchanged (its machinery checked to R^11 exactly in WA7-0(b)).
- Cusp-free H function: phi_a - (a/b)^(5/2) phi_b. 'both': on Y and Z. 'one': 1s (exponent zeta_H) on Y, cusp-free on Z.
- S = beta_p (1 - x)^(p/2) + beta_{p+2} (1 - x)^((p+2)/2) + ..., every lower odd power exactly 0:

| Case | kind | zeta_X | zeta_H (1s) | pair (a, b) | D | p | beta_p | beta_{p+2} |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| c1 | both | 2.25 | - | (1.24, 3) | 1.809 | 11 | 0.008900840359941 | 0.019526793502478 |
| c2 | one | 2.25 | 1.24 | (1.24, 3) | 1.809 | 9 | -0.012027885734131 | -0.022171191146138 |
| c3 | both | 1 | - | (1, 2) | 2.4 | 11 | 0.027259672414193 | 0.049175250807846 |
| c4 | one | 1 | 1 | (1, 2) | 2.4 | 9 | -0.052058402179883 | -0.084104743262045 |

- L1_l, L2_l = Legendre coefficients of the leading term alone and of the two-term S, by the closed form
  (2l+1) 2^a Gamma(a+1)^2 / (Gamma(a+l+2) Gamma(a+1-l)). **Sign, stated correctly this time:** these coefficients are
  positive at large l for a = 11/2 and negative for a = 9/2, so the prediction is that F_l has the sign of L2_l, which
  is positive in all four cases.
- **Bound rule (fixed before computing).** Expected size of the omitted order at l:
  e_l = 4 (a+1)^2 (a+2)^2 / l^4 times r^2, with a = p/2 and r = beta_{p+2}/beta_p. The bound is 3 e_35:
  c1 0.090, c2 0.035, c3 0.064, c4 0.026. (Applied to WA7's known cases this rule gives 0.05 or more, above every
  WA7 deviation found; disclosed, since WA7 is known.)

## What I know before comparing (disclosed)

- No F_l for c1-c4 has been computed. WA7's four cases are known (disposition in its record).
- F_l by `debug/wa7_score.py`'s F_combo (the unchanged channel formula, bilinear in the orbital pair), with D passed in.

## Predictions

**WA8-0 (controls; run FIRST; a case failing (b) is void, and WA8a-b are scored on the rest).**
- (a) WA7's E1 'both' case recomputed with the WA8 code gives F_30 within 1e-6 relative of the stored WA7 value.
- (b) Per case, |F_l(level 2) / F_l(level 1) - 1| < 1e-3 for every l in 20..35.

**WA8a (size and sign, nothing fitted).** In every non-void case, |F_35 / L2_35 - 1| < its bound above, and F_l has
the sign of L2_l for every l in 20..35.
- Falsifier: a bound or the sign broken in any non-void case.
- Probability 0.7.

**WA8b (the remainder's order).** In every non-void case, with dev_l = |F_l / L2_l - 1|: dev_20 / dev_35 > 3 (an l^-4
remainder gives 9.4), and |F_35 / L1_35 - 1| >= 5 dev_35.
- Falsifier: either inequality fails in any non-void case.
- Probability 0.65.

## Decision gate

- **GO:** WA8-0(a) passes, at least three cases survive (b), and WA8a and WA8b hold. Then zero contact removes the
  l^-8 term and the derived next term replaces it, with nothing fitted, on cases not seen before.
- **BORDERLINE:** WA8a holds, WA8b fails.
- **STOP:** WA8-0(a) fails, fewer than three cases survive (b), or WA8a fails.

## Caps

- At most 30 min wall. No change to `geovac/` or `tests/`.
- Symmetric geometry; two cusp-free pairs; 1s on X. Not a claim about Gaussian or p functions.

## Addendum 1 (appended 2026-10-27, before computing; the predictions above are not edited)

- **A disclosure was wrong.** The bound-rule bullet says that applied to WA7's cases the rule gives "0.05 or more".
  I had not computed it. Computed now from WA7's stored betas: 0.030, 0.048 ('both', zeta_X = 1, 2.25) and 0.013, 0.019
  ('one'), against WA7's deviations at l = 35 of 0.0060, 0.0103, 0.0026, 0.0041. So on WA7's cases the rule sits 4-5
  times above what was found, not above 0.05.
- **The printed bounds are rounded loosely.** The rule gives c1 0.0915, c2 0.0347, c3 0.0618, c4 0.0267. The scored
  bounds are the ones printed in the record (0.090, 0.035, 0.064, 0.026), unchanged.

## Disposition

(Appended after comparing. The predictions above are not edited.)

### Disposition (appended 2026-10-27 by the PM; computed in the main session)

Scorer `debug/wa8_score.py`, log `debug/data/wa8_score.log`, data `debug/data/wa8_score.json`. Derivation
`debug/wa8_local.py` and `debug/wa7_local.py` unchanged since registration.

| Item | c1 | c2 | c3 | c4 | Verdict |
|:--|:--|:--|:--|:--|:--|
| WA8-0(a) | WA7 F_30 reproduced to 12 digits | | | | PASS |
| WA8-0(b) noise, l = 20..35 | 1.1e-9 | 1.9e-8 | 7.1e-10 | 1.4e-8 | PASS, 4 live |
| WA8a dev_35 (bound) | 1.88e-2 (0.090) | 7.0e-3 (0.035) | 1.22e-2 (0.064) | 5.1e-3 (0.026) | |
| WA8a sign of L2, l = 20..35 | yes | yes | yes | yes | **HELD** |
| WA8b dev_20 / dev_35 (> 3) | 14.9 | 11.5 | 13.6 | 11.3 | |
| WA8b gain over L1 at 35 (>= 5) | 7.3 | 12.2 | 9.4 | 14.7 | **HELD** |

**Gate: GO.** On four cases not computed before (two cusp-free pairs, two oxygen exponents, two bond lengths, one or
both H functions cusp-free), zero contact removes the l^-8 term, and the derived next term, (1 - x)^(11/2) with both
H functions cusp-free and (1 - x)^(9/2) with one, reproduces the coefficients to 0.5-1.9% at l = 35 with nothing
fitted, with the sign the derivation gives. The remainder falls slightly faster than l^-4 between l = 20 and 35
(ratios 11-15 against 9.4), as in WA7. Each deviation is 4-5 times inside its bound.
