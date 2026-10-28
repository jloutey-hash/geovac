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

### Disposition (appended 2026-10-27 by the PM; computed in the main session)

Scorer `debug/wa6_score.py`, log `debug/data/wa6_score.log`, data `debug/data/wa6_score.json`. Derivation script
`debug/wa6_local.py` unchanged since registration.

| Item | Result | Verdict |
|:--|:--|:--|
| WA6-0(a) | closed form vs 40-digit quadrature: worst relative 1.5e-29 | PASS |
| WA6-0(b) | two-centre singular part 2R^7/315 + 16R^9/14175, exact | PASS |
| WA6-0(c) | max relative noise floor, l = 30..60: 4.3e-5 (E1), 4.4e-5 (E2) | PASS |
| WA6a | dev at l = 60: 1.02e-4 (E1), 1.74e-4 (E2); max over 50..60: 2.1e-4, 3.6e-4 | **HELD** |
| WA6b | dev_30 / dev_60 = 16.21 (E1), 16.23 (E2); gain over L0 at l = 60: 101x, 76x | **HELD** |

dev_l = |F_l / L2_l - 1| for E1 at l = 15, 20, 30, 40, 50, 60: 3.1e-2, 8.8e-3, 1.6e-3, 5.1e-4, 2.1e-4, 1.0e-4; one sign
throughout. With the leading term alone: 1.5e-1, 4.0e-2, 1.0e-2 at l = 15, 30, 60.

**Gate: GO.** The symmetric-water bond-angle series is the derived cusp-cusp singular part at gamma = 0. Its first two
orders, alpha7 (1 - x)^(7/2) + alpha9 (1 - x)^(9/2), come in closed form from the H exponent, the O exponent and the bond
length, and reproduce the coefficients to 1e-4 at l = 60 with nothing fitted. The remainder falls like l^-4 (ratio 16.2,
predicted about 16), the order of the first omitted term.

**Deviation, disclosed.** The registration said the noise floor was below 1.4e-5 for l <= 60. Over l = 30..60 its maximum
is 4.4e-5. The registered control bound (1e-4) is met, and every scored deviation is at least twice the floor.

**Caps.** Symmetric geometry and 1s orbitals with two exponent sets only. The asymmetric cases are not scored: their
leading size needs the gradient prefactor at the complex singular point (the WA5 correction above).

### Novelty check (appended 2026-10-27, after the disposition; PI-directed explorer, PM-verified)

- Closest prior art, credited in Paper 59: H. J. Silverstone, J. Chem. Phys. 46, 4377 (1967). Its abstract gives the
  small-R series of two-centre Coulomb integrals with an even family and a non-analytic family, the latter's start set
  by the radial powers, plus the Poisson route. The full text was not read, so whether its non-analytic family starts
  at R^7 for 1s densities is unchecked.
- Partial, not prior art: Bouferguene and Jones, J. Chem. Phys. 109, 5718 (1998), convergence of the Barnett-Coulson
  series for three-centre nuclear-attraction integrals; Joy and Handler, J. Chem. Phys. 43, S252 (1965), convergence
  laws for one-centre H2+ expansions. Abstracts only; neither states a cusp-set rate in the abstract. Full texts owed
  before any novelty claim.
- Analogies: Schwartz, Hill, Kutzelnigg-Morgan, Wang (electron cusps); Lehmann, Froissart-Gribov (partial waves in
  cos theta).
- No source found for the cusp-cusp origin, the constant C, the bond-angle rates, or the two-term symmetric match.

**Novelty check, full text (appended 2026-10-27).** A free copy of Bouferguene and Jones 1998 was found and read in full
(PM, untrusted download kept outside the repo). Their Eq. (4.24): for the three-centre nuclear-attraction integral over
s orbitals, the Barnett-Coulson series terms behave as K (min(a,c)/max(a,c))^l / (l+1/2)^([(N2+1)/2]+2) P_l(cos),
linearly convergent for a != c and logarithmic at a = c; they credit LaBudde and Sahni, J. Chem. Phys. 33, 1015 and 1021
(1960), for asymptotic forms of such integrals, and Lowdin for tying the slow convergence to the cusp. This is prior art
for the FORM of the rate (geometric in a distance ratio times a power of l) in the one-electron case, and is credited in
Paper 59. Not in it: the two-electron (XY|XZ) integral, the cusp-cusp term C |Y - Z|^7 and its size, the bond-angle
singularity reading, and the two-term symmetric match. Silverstone 1967 and Joy-Handler 1965: no free copy found.
A secondary source reports Joy-Handler's rule as l^-4 energy increments for one-centre H2+ (not credited in the paper
until read).
