# Pre-registered prediction: no contact, no l^-8 tail — the bond-angle series with a cusp-free H function (WA7)

**Registered:** 2026-10-27, before any F_l in this file is computed. PI-directed ("Yes to both", to a "no contact, no
tail" test of the contact form).
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** `docs/predictions/2026-10-27_water_cusp_next_term.md` (WA6, GO) and its derivation `debug/wa6_local.py`.

## The question

WA6 found the leading bond-angle singular term of (XY|XZ) at symmetric water. Written out, it factors as
(pi^2/630) g_Y g_Z |Y - Z|^7, where g_Y is the coefficient of |r - Y| in the density chi_X chi_Y at Y (its "contact";
by Kato's rule g = -zeta rho(Y) for a 1s function) and pi^2/630 is universal (checked symbolically: the unit-contact
singular part is exactly pi^2 |w|^7 / 630, and it rebuilds C = zeta_H^5 zeta_X^3 e^{-2 zeta_X D}/630). This is the
Coulomb analogue of Tan's contact (credit: Hofmann, Barth and Zwerger, arXiv:1304.2891, and Kimball's earlier work).
If the reading is right, an H function with zero contact must lose the l^-8 tail, and the next term, derived the same
way, must take over with nothing fitted.

## The derivation (done before registering)

- Cusp-free H function: chi = phi_1 - 2^(-5/2) phi_2 (phi_z the normalized 1s Slater function of exponent z). Its
  first-order cusp cancels, so the density's non-smooth part at that H is |a|^3 times an analytic function.
- Script `debug/wa7_local.py` (exact sympy, the WA6 machinery with Taylor degree 4). Degree-4 check: on the two-centre
  1s integral it gives (2/315) R^7 + (16/14175) R^9 + (8/93555) R^11, equal to the closed form's odd terms exactly.
- D1 = D2 = D = 1.809, S = beta_p (1 - x)^(p/2) + beta_{p+2} (1 - x)^((p+2)/2) + ..., all lower odd powers exactly 0:
  - 'both' (cusp-free on Y and Z), p = 11: zeta_X = 1: beta11 = 0.003966041796557, beta13 = 0.004938449424033;
    zeta_X = 2.25: beta11 = 0.000490690948286, beta13 = 0.000783337174502.
  - 'one' (1s, zeta_H = 1, on Y; cusp-free on Z), p = 9: zeta_X = 1: beta9 = -0.013331310330642,
    beta11 = -0.015077280231292; zeta_X = 2.25: beta9 = -0.001649390915071, beta11 = -0.002273103875760.
- Consequences: F_l ~ l^-12 ('both') and l^-10 ('one'), against l^-8 with two 1s functions; and in 'one' the series
  changes sign (F_l < 0 at large l), where every 1s case was positive.
- L1_l = Legendre coefficients of the leading term alone; L2_l = of the two-term S. Closed form
  (2l+1)/2 int (1-x)^a P_l dx = (2l+1) 2^a Gamma(a+1)^2 / (Gamma(a+l+2) Gamma(a+1-l)).

## What I know before comparing (disclosed)

- No F_l with the cusp-free function has been computed at any l or grid level.
- F_l is computed by the unchanged channel formula of `tests/test_paper59_water_bond_angle.py`, made bilinear in the
  orbital pair (the cusp-free function is a two-term sum, so each channel amplitude is a two-term sum). The two terms
  cancel at large l, so the noise floor is unknown; control (c) sizes it.
- Rough size of the omitted order, not a derivation: relative 4 (a+1)^2 (a+2)^2 / l^4 times beta_{p+4}/beta_p, about
  1e-2 ('both') and 6e-3 ('one') at l = 30 if that ratio is of order one.

## Predictions

**WA7-0 (controls; run FIRST; a case whose control (c) fails is void, and WA7a-b are scored on the rest).**
- (a) The closed form above agrees with 40-digit quadrature to a relative 1e-12 for a = 9/2, 11/2, 13/2 at l = 15, 30.
- (b) The degree-4 two-centre check reproduces 2/315, 16/14175, 8/93555 exactly (re-run in the scorer).
- (c) Per case, |F_l(level 2) / F_l(level 1) - 1| < 1e-3 for every l in 15..35.
- (d) The shared code still gives the WA6 E1 1s result: |F_30 / L2_30 - 1| within 10% of 1.6e-3.

**WA7a (size and sign, nothing fitted).** In every non-void case, |F_l / L2_l - 1| < 5e-2 ('both') or < 2e-2 ('one') at
l = 30, and F_l has the sign of beta_p for every l in 15..35.
- Falsifier: either bound or the sign broken in any non-void case.
- Probability 0.6.

**WA7b (the remainder's order).** In every non-void case, with dev_l = |F_l / L2_l - 1|, dev_15 / dev_30 > 6 (the next
order predicts about 16), and dev_30 is at least 5 times smaller than |F_30 / L1_30 - 1|.
- Falsifier: either inequality fails in any non-void case.
- Probability 0.55.

## Decision gate

- **GO:** WA7-0 (a), (b), (d) pass, at least two cases survive (c), and WA7a and WA7b hold. Then the bond-angle tail is
  set by the contacts: zero contact removes the l^-8 term, and the derived next term replaces it with nothing fitted.
- **BORDERLINE:** WA7a holds, WA7b fails.
- **STOP:** WA7-0 (a), (b) or (d) fails, fewer than two cases survive (c), or WA7a fails.

## Caps

- At most 30 min wall. No change to `geovac/` or `tests/`.
- Symmetric geometry; one cusp-free pair (exponents 1 and 2); 1s on X. Not a claim about Gaussians or p functions.

## Novelty status (before comparing)

- Contact relations for Coulomb systems (momentum tails set by the pair density at contact): KNOWN (Hofmann-Barth-Zwerger;
  Kimball). Darboux asymptotics: KNOWN.
- The contact factorization of the bond-angle term, and the zero-contact test: NOT FOUND in a short search. Treat as OPEN.

## Disposition

(Appended after comparing. The predictions above are not edited.)

### Disposition (appended 2026-10-27 by the PM; computed in the main session)

Scorer `debug/wa7_score.py` (with `debug/wa7_twocentre.py`), log `debug/data/wa7_score.log`, data
`debug/data/wa7_score.json`. Derivation `debug/wa7_local.py` unchanged since registration.

| Item | Result | Verdict |
|:--|:--|:--|
| WA7-0(a) | closed form vs 40-digit quadrature: worst relative 5e-30 | PASS |
| WA7-0(b) | 2R^7/315 + 16R^9/14175 + 8R^11/93555, exact | PASS |
| WA7-0(c) | noise floor l = 15..35: 8.6e-10, 1.0e-9 ('both'); 1.6e-8, 1.9e-8 ('one') | PASS, 4 cases live |
| WA7-0(d) | WA6 E1 1s dev_30 = 1.647e-3 | PASS |
| WA7a size | dev_30 = 1.15e-2, 1.99e-2 ('both', zeta_X = 1, 2.25; bound 5e-2); 4.9e-3, 7.7e-3 ('one'; bound 2e-2) | held |
| WA7a sign | 'both': F_l > 0, the sign of beta11, held. 'one': F_l > 0 at every l in 15..35, beta9 < 0 | **MISSED** |
| WA7b | dev_15 / dev_30 = 33.0, 44.9, 22.7, 23.8 (> 6); gain over the leading term at l = 30: 9.5, 6.8, 14.6, 11.2 (>= 5) | **HELD** |

**Gate: STOP** (WA7a missed, on its sign clause, in the two 'one' cases).

**Cause of the miss: an error in the registration, not in the derivation.** The registration said F_l would have the
sign of beta_p. But the large-l Legendre coefficients of (1 - x)^a do not keep one sign across a: they are positive for
a = 7/2 and 11/2 and negative for a = 9/2 (the sign of Gamma(a + 1 - l)). So the derived prediction for 'one' is
beta9 times a negative coefficient: F_l > 0. That is what was found, and F_l / L2_l - 1 = +4.9e-3 and +7.7e-3 at
l = 30. The registered sign clause is therefore wrong as written, and the scored gate stands as STOP. The sentence
"in 'one' the series changes sign" in the derivation section was the same error.

**Reported only (not scored as a GO).** With zero contact on one or both H functions, the l^-8 term is absent, and the
derived next term matches the coefficients with nothing fitted: to 1-2% at l = 30 with both H functions cusp-free
(F_30 about 1e-14, falling like l^-12) and to 0.5-0.8% with one (l^-10). In every case the remainder falls faster than
the next order predicts between l = 15 and 30 (ratio 23-45, against about 16), so l = 15 is not yet asymptotic.
A fresh registration on unseen cases would be needed to score this.
