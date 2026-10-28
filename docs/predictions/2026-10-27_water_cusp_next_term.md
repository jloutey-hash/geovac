# Pre-registered prediction: the symmetric-water bond-angle series from the derived cusp-cusp terms, nothing fitted (WA6)

**Registered:** 2026-10-27, before any comparison that scores a prediction in this file. PI-directed ("Yes", to a WA6
that derives the next term so that nothing is fitted).
**Commit status:** committed on its own, before any comparison in this file.
**Rests on:** `docs/predictions/2026-10-27_water_cusp_singularity.md` (WA5, STOP on a control; its water numbers
reported only) and the channel formula of `tests/test_paper59_water_bond_angle.py`.

## The question

For symmetric water (D1 = D2 = D = 1.809), the Legendre coefficients F_l of (XY|XZ) in x = cos(gamma) are set, at large
l, by the non-analytic part of the integral at gamma = 0, where the two H atoms merge. If that part is derived through
its first two orders, its own Legendre coefficients L_l should match F_l with no fitted parameter, and the mismatch
should shrink like the next order.

## The derivation (done before registering)

- The non-analytic part comes only from the H cusp of rho1 = chi_X chi_Y at Y meeting the H cusp of rho2 = chi_X chi_Z
  at Z through 1/r12 (WA5 record). Writing the non-smooth part of rho1 as |a| G1(a), a = r - Y, with
  G1(a) = -N_H N_X e^{-zeta_X |Y + a|} sinh(zeta_H |a|)/|a|, and likewise at Z, the singular part is
  S(w) = FT^-1[F1hat(-k) F2hat(k) 4 pi/k^2], w = Y - Z, built term by term from the Taylor polynomials of G1 and G2.
- **Correction to the WA5 record.** WA5 said the gradient terms of G give analytic pieces. That is wrong: |a| a_j
  against |b| gives a term w_j |w|^7, which is not analytic. At D1 = D2 it enters at the next order only (w . Y-hat =
  D(1 - x)), but at D1 != D2 it changes the leading size, since w . Y-hat is not zero at the complex singular point.
  This is a likely cause of the scatter of WA5's fitted amplitudes in the asymmetric cases. WA6 therefore scores the
  symmetric case only.
- At D1 = D2 every invariant of w (w . w, w . Y-hat, w . Z-hat) is O(1 - x), so the (1 - x)^(7/2) and (1 - x)^(9/2)
  coefficients come only from Taylor degrees dA + dB <= 2. The next non-analytic order is (1 - x)^(11/2).
- Script `debug/wa6_local.py` (sympy, exact). Its check on the two-centre Coulomb integral of e^{-2r}/pi densities
  reproduces both odd terms of the exact closed form, (2/315) R^7 + (16/14175) R^9, exactly.
- Result, S = alpha7 (1 - x)^(7/2) + alpha9 (1 - x)^(9/2) + O((1 - x)^(11/2)):
  - E1 (zeta_X = zeta_H = 1): alpha7 = 0.0305532186374231, alpha9 = 0.0285025232663417 (ratio 0.9328812).
  - E2 (zeta_X = 2.25, zeta_H = 1.24): alpha7 = 0.0110819491640684, alpha9 = 0.0134090231083288 (ratio 1.2099878).
  - alpha7 equals C (2 D^2)^(7/2), the WA5 leading term, in both sets. No (1 - x)^(1/2), ^(3/2), ^(5/2) terms.
- L0_l = Legendre coefficients of alpha7 (1 - x)^(7/2); L2_l = those of the two-term S. Both in closed form from
  int_{-1}^{1} (1 - x)^a P_l(x) dx, checked against 40-digit quadrature.

## What I know before comparing (disclosed)

- The F_l used here are the stored WA5 values (`debug/data/wa5_cusp.json`, grid level 2, l to 60), computed by the
  unchanged test formula. Their relative noise floor at D1 = D2 is below 1.4e-5 for l <= 60.
- WA5 printed F_l / L0_l at l = 15 and l = 60 for these two cases: 0.850 and 0.990 (E1), 0.814 and 0.987 (E2). I have not
  computed F_l / L2_l at any l.
- Rough size of the next order, not a derivation: the (1 - x)^(11/2) term changes the ratio by about 2400 (alpha11/alpha7)
  / l^4 relative, about 2e-4 at l = 60 and 3e-3 at l = 30 if alpha11/alpha7 is of order one.

## Predictions

**WA6-0 (controls; run FIRST; WA6a-b void if any fails).**
- (a) The closed-form Legendre coefficients of (1 - x)^(7/2) and (1 - x)^(9/2) agree with 40-digit quadrature to a
  relative 1e-12 at l = 15, 30, 45, 60.
- (b) The two-centre check above reproduces 2/315 and 16/14175 exactly (re-run in the scoring script).
- (c) The stored noise floor at D1 = D2 (level 1 vs level 2) is below 1e-4 relative for every l in 30..60. (The
  channel sum's referee control passed in WA5-0(a) on the same code.)

**WA6a (the size, nothing fitted).** For both sets, |F_l / L2_l - 1| < 1e-3 at l = 60, and below 2e-3 for every l in
50..60.
- Falsifier: either bound broken in either set.
- Probability 0.7.

**WA6b (the order of the remainder).** For both sets, with dev_l = |F_l / L2_l - 1|, dev_30 / dev_60 > 6 (the next order
predicts about 16; a missed (1 - x)^(9/2) piece would give about 4). And dev_60 is at least 10 times smaller than
|F_60 / L0_60 - 1|.
- Falsifier: either inequality fails in either set.
- Probability 0.65.

## Decision gate

- **GO:** WA6-0 passes and WA6a and WA6b hold. Then the symmetric-water bond-angle series is the derived cusp-cusp
  singular part at gamma = 0, its first two orders in closed form, with no fitted number.
- **BORDERLINE:** WA6a holds, WA6b fails (the size matches but the remainder's order is not established).
- **STOP:** WA6-0 fails, or WA6a fails.

## Caps

- At most 30 min wall. No change to `geovac/` or `tests/`.
- The asymmetric cases are not scored here. Their leading size needs the exact w . Y-hat prefactor at the complex
  singular point, which this derivation does not give.

## Novelty status (before comparing)

- Darboux asymptotics of Legendre coefficients: KNOWN.
- Cusp-generated odd powers in two-centre Slater integrals: the R^7, R^9 terms are properties of the known closed form;
  the local cusp-cusp derivation of them and of the three-centre angle series: NOT SEARCHED. Treat as OPEN.

## Disposition

(Appended after comparing. The predictions above are not edited.)
