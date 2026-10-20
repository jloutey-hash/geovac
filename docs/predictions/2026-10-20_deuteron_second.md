# Pre-registered prediction: the deuteron's size, second registration (DS1-DS5)

**Registered:** 2026-10-20, before any computation in this file (PI-directed: "Both owed please", after the
five MISSes of `2026-10-20_deuteron_radius.md`). Drafted by the PM.
**Commit status:** committed on its own, before any computation that scores it.
**Extends:** `docs/predictions/2026-10-20_deuteron_radius.md` (DR1-DR5, all MISS); Paper 23 `sec:deuteron`
and `sec:finite_nucleus`; seam rows J14, J16; atlas doors 24 and 26.

## Why a second registration, and what it may not do

DR1-DR5 missed on gate placement, not on physics (disposition of the first file). Re-scoring those numbers
with moved gates would be fitting the gates to the data. So this file registers **only quantities nobody has
computed yet**: a second potential (Minnesota triplet, the corpus's own Paper 23 interaction), and a
published parameter set that may explain the one unexplained number of the first run. The MT-III numbers
of the first run appear below only as disclosed knowledge and as a regression control. They are not scored
again.

Three lessons from the first file are applied to how the gates are set:
- A gate is placed with a stated margin over the nearest known analogue, not at an estimate's edge.
- The scale lambda is the energy-optimal one. Moving lambda off it to avoid a tautology moved the basis out
  of its converged region (DR2).
- An arithmetic claim names the basis it holds in. In the unnormalized Laguerre basis the Yukawa matrices are
  rational; the orthonormal ones are not (DR3). No arithmetic claim is made here.

## What I know before computing (disclosed)

- **From the first run (MT-III, scored, `debug/deuteron_radius_results_memo.md`):** r_m = 1.9897135 fm; B =
  2.2306877 MeV; a_t = 5.5132023 fm; A_S = 0.9043067 fm^-1/2. Laguerre r_m errors at the energy-optimal
  lambda: N_r = 30, 1.6e-6; 35, 7.5e-8; 40, 2.1e-9. A_S flatness on [8, 16] fm at lambda_40: 4.2e-6.
  (r_ERE - r_m)/r_m = 2.005%. The HO basis at N_r = 40 was UV-limited by MT-III's 1/r core: B still 1.5% low
  at hbar-omega = 60 MeV. Energy-optimal lambda was 4 to 6 kappa.
- **The unexplained number.** E. Biernat, PhD thesis, arXiv:1110.3180, p. 90 (read today, Table 6.1 and
  Eq. 6.22-6.23): "close to Malfliet-Tjon's value sqrt<r^2> = 1.97625 fm", with <r^2> = (1/4) int r^2 w_0^2,
  the same definition as our r_m. Table 6.1's MT [MT69] column:
  E_B = -2.27203 MeV, a_t = 5.4739 fm, g_sigma^2/4pi = 3.22749, g_omega^2/4pi = 7.40758,
  lambda_sigma = 0.643087 fm (306.8 MeV), lambda_omega = 0.321543 fm (613.6 MeV). The potential is
  V(r) = hbar c [ -(g_sigma^2/4pi) e^{-r/lambda_sigma} + (g_omega^2/4pi) e^{-r/lambda_omega} ] / r, with
  kinetic term hbar^2/(2 m_red), 2 m_red = 938.92 MeV (Eq. 6.21). In MeV fm this is about 1461.7 and 636.9
  with ranges 3.110 and 1.555 fm^-1, against MT-III's 1438.72, 626.885, 3.11, 1.55. It binds more
  (2.272 against 2.2307 MeV), which would shrink r_m. Nothing has been computed with it.
  The thesis does not state hbar c. The table's own (mass, lambda) pairs give 197.30 to the 4 digits of the
  masses.
- **MT-V is not the explanation.** Its two-body state binds by only about 0.35-0.41 MeV (two arXiv sources
  disagree on the digit; web search today), far too weak for r_m near 2 fm.
- **Minnesota triplet.** V(r) = 200 e^{-1.487 r^2} - 178 e^{-0.639 r^2} MeV (r in fm), hbar^2/m = 41.47 MeV fm^2,
  from `geovac/nuclear/minnesota.py` and Paper 23 eq:minnesota (Thompson, Lemere and Tang 1977). In the
  3S1 channel the exchange factor is 1. Published deuteron binding 2.202 MeV (arXiv:1706.06455, a table
  entry; the primary source was not read). No published Minnesota a_t, A_S or r_m was found (one search).
  The core is a soft Gaussian, so u = r + O(r^3) near the origin, with no r^2 term.
- Nothing in this file has been computed.

## Model

Relative 3S1, l = 0, no fitted parameter. Laguerre basis, HO basis and ODE route exactly as in the first
file (`2026-10-20_deuteron_radius.md`, Model section), with these changes:
- **Laguerre elements for the Gaussian potential** are not rational. They are evaluated in closed form
  (erfc / parabolic-cylinder) or by quadrature at dps 50, and cross-checked against each other.
- **Scale rule:** lambda_N minimizes E_N at each N_r over a rational grid k/10000, as before.
- **Potentials:** P1 = Biernat's MT69 set as above, with hbar c = 197.327 MeV fm and kinetic coefficient
  hbar^2/m_N, m_N = 938.92 MeV (so hbar^2/m_N = 41.4718 MeV fm^2). P2 = Minnesota triplet as above.
- **Computed for each:** B, a_t, r_m, A_S, kappa = sqrt(B/(hbar^2/m)), r_ERE as in DR5, and for P2 the
  charge radius r_ch as in DR4 (reported, not scored).

## Predictions

**DS1 (the 0.68% gap is the parameter set).** For P1, both routes give B = 2.27203 MeV and a_t = 5.4739 fm
within 2e-3 relative, and r_m = 1.97625 fm within 1e-3 relative.
- Gate placement: r_m scales roughly as B^{-1/2}, and near threshold B moves by an estimated 10 to 20 parts per part of a
  strength change (not computed), so the undisclosed hbar c (197.327 against 197.33) can move B by about 3e-4 and r_m by
  about 2e-4. The gates are 6x and 5x that.
- Falsifier: any of the three outside its gate. If DS1 fails, the gap stays unexplained, and the record
  says so.
- Consequence if it holds: the first run's r_m and the thesis value differ by the parameter set alone; at
  the 1e-3 gate that accounts for at least 85% of the 0.68% gap.
- Novelty: KNOWN (it reproduces a published number). A hit is a **CHECK**.

**DS2 (the exponential tail converges the radius for a soft core).** For P2, at N_r = 30 and
lambda = lambda_30, |r_m(N_r) - r_m(ODE)|/r_m(ODE) <= 1e-6.
- Gate placement: MT-III, with a 1/r core, reached 1.6e-6 at N_r = 30 and 7.5e-8 at 35. The soft core
  removes the r^2 term the Laguerre functions had to build near the origin, so I expect at least one order
  below MT-III at equal N_r. Confidence moderate: the Gaussian shape of the interior could cost what the core
  saves.
- Falsifier: error above 1e-6.
- Novelty: PROBABLY KNOWN (exponential-tail bases, Caprio et al. 2012/2014; no deuteron comparison found in
  the first file's searches). A hit is a **CHECK**.

**DS3 (the first run's failed HO control was the core).** For P2, in the HO basis at N_r = 40, the
energy-optimal hbar-omega in the scan {2, 4, 6, 8, 10, 15, 20, 30, 40, 60} MeV gives B within 1e-4 relative
of the ODE B. For MT-III the same setting left B 1.5% low.
- Gate placement: at hbar-omega near 20-30 MeV, L_2 = sqrt(2(N + 7/2)) b is above 20 fm, so the IR error
  e^{-2 kappa L_2} is below 1e-4, while the UV cutoff sqrt(2(N + 7/2))/b is several times the Gaussian core's
  momentum scale. Confidence high.
- Falsifier: no hbar-omega in the scan within 1e-4.
- Novelty: KNOWN (More et al. 2013, the UV/IR picture). A hit is a **CHECK**, and it confirms the mechanism
  named in the first disposition.

**DS4 (the tail rate at the energy-optimal scale).** For P2 at N_r = 40 and lambda = lambda_40, A_S(r) =
u(r) e^{kappa_N r} varies by <= 1e-4 relative over r in [8, 16] fm, and its mean equals the ODE A_S to
1e-4.
- Gate placement: MT-III's unscored analogue was 4.2e-6 and 5.5e-8, so the gate has a margin of about 20.
- Falsifier: either quantity above 1e-4.
- Novelty: PROBABLY KNOWN (first file, DR2). A hit is a **CHECK**.

**DS5 (the effective-range overshoot tracks the core).** For P2, with kappa and a_t computed from the
potential, 0 < (r_ERE - r_m)/r_m < 2.005%, the MT-III value.
- Reasoning: the shape-independent formula treats the wavefunction as its asymptotic form all the way in.
  A harder core depletes the interior more, so the overshoot should be larger for MT-III than for a soft
  Gaussian core. Experiment's overshoot is about 0.6%.
- Confidence low. The r^2 weighting suppresses the interior, and the shape parameter can have either sign.
- Falsifier: the ratio outside (0, 2.005%).
- Novelty: PROBABLY KNOWN (Bethe's shape parameter). A hit is a **CHECK**.

## Controls

- **Regression:** the first run's MT-III numbers (B, a_t, r_m at N_r = 40 and by ODE) reproduced to 1e-8 by
  the same code that computes P1 and P2.
- **Positive, KNOWN:** P2's B against the published 2.202 MeV, within 1e-3 relative (the published value has
  four digits; a miss is reported and voids nothing).
- **Gaussian elements:** for P2, the closed-form and quadrature Laguerre elements agree to 1e-30 at
  N_r <= 5 and two lambdas.
- **Independent route:** every scored basis number is compared with the ODE route.
- **Hydrogen:** the exact checks of the first file, rerun.

## Honest caps

- Every prediction here is CHECK-class. This file closes a loose end and tests the first file's lessons on a
  second potential. It does not look for a FINDING.
- Two potentials are a sample. The Minnesota triplet is fitted to low-energy data, like MT-III, so neither
  says anything about the true deuteron beyond the potential's own error.
- DS1 depends on a thesis table whose hbar c is not stated. The gate absorbs that ambiguity, so a hit cannot
  distinguish the two conventions.

## Result

(To be appended after scoring; the predictions above are never edited.)
