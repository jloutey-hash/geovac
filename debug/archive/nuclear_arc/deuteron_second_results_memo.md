# Results memo: DS1-DS5, the deuteron's size, second registration (scored 2026-10-20)

Scores `docs/predictions/2026-10-20_deuteron_second.md` (registered in 1e9aaa71). The prediction file was
not edited, and no tracked file was modified. Scripts: `debug/deuteron_second_{laguerre,ode,ho,controls,score}.py`
(copies of the first run's scripts, generalized to Gaussian terms and a per-potential hbar^2/m), plus two
unscored diagnostics, `debug/deuteron_second_{floor,hbarc}.py`. Data: `debug/data/deuteron_second_{ode,controls,score,floor,hbarc}.json`
with `.log` captures. Run order: ODE, controls, scoring, diagnostics.

**Headline.** DS1, DS3, DS4 and DS5 HELD, so they are four CHECKs. DS2 FAILED by a factor of 4.3, so it is a
MISS. Every control passed. The independent routes agree to about 1e-8 or better on every scored number, so
no verdict depends on numerical error. The DS2 miss goes the opposite way to the registration's reasoning:
the soft Gaussian core converges r_m more slowly than MT-III's 1/r core at N_r = 30 (4.3e-6 against 1.6e-6),
not an order of magnitude faster.

## 1. Controls (run first)

| Control | Result | Gate | Verdict |
|:--|:--|:--|:--|
| Hydrogen, N_r = 1, lambda = 1 | E = -1/2 exactly (Fraction) | exact | PASS |
| Hydrogen, N_r = 2, lambda = 1/2 | det(H_u + S_u/8) = 0 exactly; eigenvalues -11/24 and -1/8 | exact | PASS |
| Regression MT-III, Laguerre N_r = 40 | lambda_40 re-found = 13949/10000 (same); B = 2.2306877450095 MeV, rel diff 7e-26; r_m = 1.9897134909168 fm, rel diff 2e-20 | 1e-8 | PASS |
| Regression MT-III, Kohn a_t | 5.5132023198067 fm, rel diff 6e-21 | 1e-8 | PASS |
| Regression MT-III, ODE | B, r_m, a_t bit-identical to the first run (rel diff 0) | 1e-8 | PASS |
| Gaussian elements, P2, N_r <= 5, lambda = 1/2 and 1: closed form vs dps-50 quadrature of <phi_m|V|phi_n> | max abs diff 2.1e-50 at both lambdas (one ulp at dps 50) | 1e-30 | PASS |
| erfc route vs parabolic-cylinder route, M_q (q <= 82), lambda = 1/2, 1, 936/625 | max rel diff 1e-115, 1e-108, 9e-103; N_r = 40 orthonormal V identical at dps 60 (also working dps 200 vs 350) | -- | PASS |
| P2 B vs published 2.202 MeV | Laguerre N_r = 40: 2.20226919 MeV, +1.22e-4; ODE: +1.22e-4 | 1e-3 | PASS |
| HO Gaussian elements, closed form vs quadrature (4 elements, hbar-omega = 20 MeV) | <= 5e-41 | -- | PASS |
| HO regression MT-III, N_r = 40, 60 MeV | B = 2.19781299937 MeV, same as the first run | -- | PASS |

ODE floor (unscored diagnostic, `deuteron_second_floor.py`). The Laguerre basis at N_r = 50 gives
P2 r_m = 1.9518278525 fm and B = 2.20226919518 MeV. Since the basis is variational, its B is a lower bound,
and it lies 5e-11 above the ODE B, so the float64 Numerov route is low by at least that much. The ODE r_m
(1.9518278469) differs from the N_r = 50 basis value by 2.8e-9. The raw ODE <r^2> values scatter at about
1e-9, and r_max = 80 fm changes nothing. So the ODE's floor is about 3e-9 in r_m and 1e-10 in B. That is
three orders of magnitude below every gate.

## 2. Predictions

### DS1: HELD -> CHECK
P1 (Biernat MT69 set, hbar c = 197.327 MeV fm, hbar^2/m_N = 197.327^2/938.92 = 41.47176 MeV fm^2).
Laguerre N_r = 40 at lambda_40 = 13987/10000; Kohn a_t at the same point; ODE route.

| Quantity | Basis | ODE | Thesis | rel. (both routes) | Gate | Margin |
|:--|:--|:--|:--|:--|:--|:--|
| B (MeV) | 2.2735331946 | 2.2735331947 | 2.27203 | +6.62e-4 | 2e-3 | 3.0x |
| a_t (fm) | 5.4723116505 | 5.4723116616 | 5.4739 | -2.90e-4 | 2e-3 | 6.9x |
| r_m (fm) | 1.9755513730 | 1.9755513770 | 1.97625 | -3.54e-4 | 1e-3 | 2.8x |

Basis vs ODE: B 3.4e-11, r_m 2.0e-9, a_t 2.0e-9. The first run's MT-III r_m was 0.681% above the thesis
value; P1 leaves 0.035%. So the parameter set accounts for 94.8% of the gap (registered consequence:
at least 85%). Unscored sensitivity (`deuteron_second_hbarc.py`, ODE): with hbar c = 197.33, B moves by
only 5.8e-5 relative (the gate placement estimated about 3e-4). The residuals become +6.04e-4 (B),
-2.55e-4 (a_t) and -3.20e-4 (r_m). Neither convention reproduces the table to its printed digits; the remaining
few-1e-4 offset is not explained here.

### DS2: FAILED -> MISS
P2 (Minnesota triplet), N_r = 30, lambda_30 = 1817/1250 = 1.4536 fm^-1 (energy-optimal on the k/10000 grid,
not at a coarse-grid edge). r_m(30) = 1.9518195437 fm against ODE 1.9518278469 fm:
|rel. error| = **4.25e-6**, against a gate of 1e-6 (4.3x over). Against the N_r = 50 basis value the error is
4.26e-6, so the ODE floor plays no part.

| N_r | lambda_N | B (MeV) | r_m (fm) | rel. err r_m vs ODE | MT-III (first run) |
|:--|:--|:--|:--|:--|:--|
| 10 | 5989/5000 | 2.19912 | 1.92230 | -1.5e-2 | -2.2e-2 |
| 15 | 6421/5000 | 2.2020936 | 1.948246 | -1.8e-3 | +1.9e-5 |
| 20 | 6133/5000 | 2.2022592 | 1.951647 | -9.3e-5 | -1.1e-4 |
| 25 | 13569/10000 | 2.20226852 | 1.9517852 | -2.2e-5 | -5.3e-6 |
| 30 | 1817/1250 | 2.202269109 | 1.95181954 | **-4.3e-6** | -1.6e-6 |
| 35 | 13539/10000 | 2.2022691889 | 1.951827811 | -1.8e-8 | -7.5e-8 |
| 40 | 936/625 | 2.2022691949 | 1.951827818 | -1.5e-8 (-1.8e-8 vs N_r = 50) | -2.1e-9 |

Convergence is not monotone (lambda_N moves). P2 is worse than MT-III at N_r = 15, 25, 30, 40 and better at
20, 35: no order-of-magnitude gain. The registration's caveat applies (the Gaussian interior costs what the
missing r^2 term saves). The optimal lambda (1.2-1.6 fm^-1, 5-7 kappa) is no smaller than MT-III's.

### DS3: HELD -> CHECK
P2, HO basis, N_r = 40, scan {2, ..., 60} MeV. The energy-optimal hbar-omega (largest B in the scan) is
**20 MeV**: B_HO = 2.2021790535 MeV against ODE 2.2022691951, rel. **-4.09e-5**, against a gate of 1e-4
(margin 2.4x). 15 MeV is also within the gate (-8.1e-5). Every HO B lies below the exact B (variational).

| hbar-omega (MeV) | 8 | 10 | 15 | 20 | 30 | 40 | 60 |
|:--|:--|:--|:--|:--|:--|:--|:--|
| rel. err B, N_r = 40 | -8.1e-3 | -2.5e-3 | -8.1e-5 | **-4.1e-5** | -3.7e-4 | -1.4e-3 | -6.7e-3 |
| L_2 (fm) | 41.1 | 36.8 | 30.0 | 26.0 | 21.2 | 18.4 | 15.0 |

The B error falls toward low hbar-omega until the UV cutoff bites, then rises again: the UV/IR trade-off
More et al. describe. For MT-III the same setting left B 1.5% low. This confirms the first disposition's
mechanism: the HO failure was the 1/r core. The HO r_m is still poor (unscored): its best error is -2.6e-4,
at 15 MeV; at 20 MeV it is -1.8e-3, and at 30 MeV -8.6e-3. That is the IR-limited radius of More et al.

### DS4: HELD -> CHECK
P2, N_r = 40, lambda_40 = 936/625. kappa_40 = sqrt(B_40/41.47) = 0.23044547 fm^-1. A_S(r) = u e^{kappa_40 r}
on 81 points in [8, 16] fm: spread (max-min)/mean = **6.25e-6** (gate 1e-4, margin 16x). Mean = 0.87374936
fm^-1/2 against ODE A_S 0.87374920 (wall-corrected [20, 30] fm plateau): **+1.80e-7** (margin 560x). The ODE
solution itself is flat on [8, 16] to 2e-10, since the Gaussian potential has died there; the basis spread is
truncation error.

### DS5: HELD -> CHECK
P2, kappa_40 = 0.23044547 fm^-1 and a_t = 5.4265163 fm (ODE; the Kohn route gives the same to 2e-10):
rho = 1.738634 fm, r_ERE = 1.9817566 fm, r_m = 1.9518278 fm, so (r_ERE - r_m)/r_m = **+1.5334%**, inside (0, 2.005%).
The margin is 0.47 percentage points to the upper edge and 1.53 to the lower. The ordering the file predicts holds:
the soft core overshoots less than MT-III (2.005%). The unscored P1 value is 2.056%.

## 3. Final classes

| ID | Verdict | Class | Number | Gate | Margin |
|:--|:--|:--|:--|:--|:--|
| DS1 | HELD | CHECK | B +6.6e-4, a_t -2.9e-4, r_m -3.5e-4 (both routes) | 2e-3, 2e-3, 1e-3 | 3.0x, 6.9x, 2.8x |
| DS2 | FAILED | MISS | 4.25e-6 | <= 1e-6 | over by 4.3x |
| DS3 | HELD | CHECK | -4.09e-5 at 20 MeV | 1e-4 | 2.4x |
| DS4 | HELD | CHECK | spread 6.25e-6; mean +1.8e-7 | 1e-4 each | 16x; 560x |
| DS5 | HELD | CHECK | 1.5334% | (0, 2.005%) | 0.47 pp to the upper edge |

## 4. Unscored numbers

| | B (MeV) | r_m (fm) | a_t (fm) | A_S (fm^-1/2) | kappa (fm^-1) | rho (fm) | r_ERE (fm) | (r_ERE-r_m)/r_m |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| MT-III (first run) | 2.2306877 | 1.9897135 | 5.5132023 | 0.9043067 | 0.2319276 | 1.879324 | 2.0296066 | 2.005% |
| P1 (MT69, Biernat) | 2.2735332 | 1.9755514 | 5.4723117 | 0.9111024 | 0.2341415 | 1.875275 | 2.0161646 | 2.056% |
| P2 (Minnesota triplet) | 2.2022692 | 1.9518279 | 5.4265163 | 0.8737492 | 0.2304455 | 1.738634 | 1.9817566 | 1.533% |

The A_S values are ODE plateaus; the basis means agree to 6e-8 (P1) and 2e-7 (P2).

**P2 charge radius**, r_ch = sqrt(r_m^2 + r_p^2 + r_n^2 + r_DF^2), with r_p = 0.84075 fm (CODATA 2022),
r_n^2 = -0.1155 fm^2, r_DF^2 = 0.0331 fm^2 and the N_r = 40 r_m: **r_ch = 2.10573 fm**, which is -1.046% against
CODATA 2018 (2.12799) and -1.036% against CODATA 2022 (2.12778). Minnesota's deuteron is too small, while
MT-III's was 0.61% too large.

For information, P1 gives 2.12774 fm (-1.2e-4 vs CODATA 2018). Not a claim: the omitted spin-orbit and
two-body terms (about 0.02 fm^2) alone shift r_ch by about 0.2%, so the agreement is below the formula's own
error.

## 5. Deviations and interpretation choices (all fixed before the scored quantity was computed)

- **D1.** Gaussian Laguerre elements: the primary route is erfc plus the upward recurrence at 200 working digits.
  It is cross-checked against the parabolic-cylinder form (all M_q, three lambdas) and against dps-50 quadrature
  (the registered control). The monomial-to-Laguerre transform runs at 200 digits and was verified against 350.
  Yukawa, kinetic, S and <r^2> are exact rationals, as before.
- **D2.** The two lambdas for the Gaussian quadrature control (1/2 and 1 fm^-1) were unregistered and chosen
  before computing. The route comparison adds lambda_40.
- **D3.** The lambda rule is the first run's D3 unchanged (coarse 0.2..3.0, step 0.2, then ternary search
  on k/10000, with a grid-local-minimum check). No optimum sat at a coarse-grid edge.
- **D4.** "Both routes" in DS1: the basis route is Laguerre N_r = 40 at lambda_40 for B and r_m, and the first
  run's Kohn functional (D2 there) at the same point for a_t. The ODE route covers all three. Both routes must
  pass each gate.
- **D5.** DS3 "energy-optimal hbar-omega" was taken as the scan point with the lowest HO energy. N_r = 40 HO
  functions, as for the Laguerre basis.
- **D6.** DS4: kappa_N = sqrt(B_40/(hbar^2/m)); spread = (max-min)/mean over 81 points; the ODE A_S is the
  wall-corrected [20, 30] fm plateau (first run's D4).
- **D7.** DS5: kappa from the N_r = 40 basis and a_t from the ODE, as in the first run. The Kohn a_t gives the
  same ratio to 8 digits. The upper bound was used as printed (2.005%); the unrounded MT-III value 2.00497%
  gives the same verdict.
- **D8.** The Kohn functional for a Gaussian potential needs the s = 1 - e^{-lambda r} couplings with the
  Gaussian. These are I_j(c, a) = int r^j e^{-cr - ar^2} dr (erfc recurrence; c = 0 is the Gamma form).
- **D9.** The ODE is the first run's float64 Numerov plus Richardson, generalized: the series start uses the
  Taylor coefficients of r V(r), and (f u)_0 = 0 for a Gaussian. Its floor is about 3e-9 in r_m (Sec. 1). It
  is therefore an independent check at the 1e-8 level, not a reference at 1e-10. No scored margin is that small.
- **D10.** P1's thesis kinetic convention "2 m_red = 938.92 MeV" was implemented, as registered, as hbar^2/m_N
  with m_N = 938.92.
