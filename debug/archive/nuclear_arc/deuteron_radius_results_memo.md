# Results memo: DR1-DR5, the deuteron's size from its own binding (scored 2026-10-20)

Scores `docs/predictions/2026-10-20_deuteron_radius.md` (registered in 9ed6f853). The prediction file was
not edited. Scripts: `debug/deuteron_radius_{laguerre,ode,ho,controls,dr3,score}.py`. Data:
`debug/data/deuteron_radius_{ode,controls,dr3,score}.json` (plus `.log` captures). The runs were
sequential: ODE, then controls, then DR3, then scoring. The old `debug/deuteron_sturmian_v2.py` was not used.

**Headline.** All five predictions FAILED as written, so the classes are five MISSes. Every control passed
except the HO consistency check (Sec. 1), and that failure changes no grade. Three of the misses are narrow:
DR1 is off by 1.6x, DR4's r_m is 0.014% under its window, and DR5's ratio is 2.005% against a <2% gate.
DR2 is not narrow. DR3 fails on a normalization oversight in the registration. The independent routes
agree to about 2e-9, so none of the misses comes from numerical error.

## 1. Controls (run first)

| Control | Result | Gate | Verdict |
|:--|:--|:--|:--|
| Hydrogen, N_r = 1, lambda = 1 | E = -1/2 exactly (Fraction) | exact | PASS |
| Hydrogen, N_r = 2, lambda = 1/2 | det(H_u + S_u/8) = 0 exactly; eigenvalues -11/24 and -1/8 | exact | PASS |
| B, ODE route (Numerov + Richardson) | 2.230687745079 MeV, rel -1.01e-6 vs 2.23069 | 1e-5 | PASS |
| B, Laguerre route (N_r = 40, lambda_40 = 13949/10000) | 2.230687745010 MeV, rel -1.01e-6 | 1e-5 | PASS |
| a_t, ODE route | 5.513202331 fm, rel +4.2e-7 vs 5.5132 | 1e-5 | PASS |
| a_t, basis route (Kohn, see deviation D2) | 5.513202320 fm (N_r = 40), rel +4.2e-7 | 1e-5 | PASS |
| MT-I singlet, Laguerre | lowest eigenvalue > 0 at every N_r = 1..40 (lambda_40), monotone non-increasing; > 0 at all 7 kappa-multiple lambdas, N_r = 10..40 | > 0, decreasing | PASS |
| MT-I singlet, ODE | zero-energy solution has no node; a_s = -23.58 fm | unbound | PASS |
| Independent route (basis vs ODE) | r_m 2.1e-9, B 3.1e-11, A_S 5.5e-8, a_t 2.1e-9 (relative, N_r = 40) | agree | PASS |
| HO consistency (More et al. L_2 curve) | **no collapse onto one curve**: at N_r = 40 B is still 1.5% low at hbar-omega = 60, and points at equal L_2 differ by up to 0.5 MeV | on the curve | **FAILED as stated** |

Orthonormality was verified exactly before use (S_u exactly diagonal, S_nn = (n+2)!/(n! 2 lambda), at every
assembled lambda). HO code: S = I to 1e-126; kinetic diagonal/off-diagonal match closed form to 1e-124/1e-153;
Yukawa elements match quadrature to 1e-43; every HO B lies above exact. ODE Richardson last-step changes
<= 1e-10 for B, <r^2>, A_S and a_t.

**The failed HO consistency control.** The More et al. IR picture assumes UV convergence. Here the 811.8/r MeV
core makes u = r + 9.79 r^2 + ... near the origin, and HO functions carry only odd powers of r. So every
(N_r, hbar-omega) point in the scan is UV-limited, and the L_2 collapse cannot appear. The registration's own
risk 3 foresaw this mechanism. Two things depend on this control: DR1's HO half and DR2's HO negative control.
Strictly, both are voided. Neither changes a grade: DR1 fails on its Laguerre half regardless, and DR2's HO
control "does not refute DR2" by registration. They are reported below for information.

## 2. Predictions

### DR1: FAILED -> MISS
- Laguerre half: N_r = 30, lambda_30 = 2763/2000. |r_m - r_m(ODE)|/r_m(ODE) = **1.62e-6**, against a gate of
  1e-6. The gate is first met between N_r = 30 and 35 (the table below). None of the kappa-multiple scan
  lambdas at N_r = 30 meets it either: the best are 4 kappa (+2.4e-6) and 6 kappa (-1.8e-6).
- HO half (voided by the consistency control; for information): at N_r = 30 the smallest error over the ten
  hbar-omega values is **1.21e-2** (at 30 MeV), so every hbar-omega is at least 1e-4 off. This half HELD.
- The falsifier fires on the Laguerre half, and the PARTLY rule covers only an HO-half failure, so DR1 is FAILED.

| N_r | lambda_N (fm^-1) | B (MeV) | r_m (fm) | rel. err r_m |
|:--|:--|:--|:--|:--|
| 10 | 47/40 | 2.22540 | 1.94673 | -2.2e-2 |
| 20 | 12969/10000 | 2.2306833 | 1.989490 | -1.1e-4 |
| 30 | 2763/2000 | 2.230687704 | 1.98971026 | -1.6e-6 |
| 40 | 13949/10000 | 2.2306877450 | 1.989713491 | -2.1e-9 |
| ODE | -- | 2.2306877451 | 1.989713495 | -- |

The energy-optimal lambda is 4-6 kappa (the core dominates the energy); convergence is not monotone in N_r
because lambda_N moves (N_r = 15: +1.9e-5; 20: -1.1e-4; 25: -5.3e-6; 35: -7.5e-8).

Novelty re-search (two web searches): Caprio et al. (arXiv:1208.4156, 1409.0108, 1409.0877) treat 4He, 6Li,
6,8He; no deuteron Laguerre-vs-HO comparison found. Stays NOT FOUND (moot: missed).

### DR2: FAILED -> MISS
Settings: N_r = 40. Lambda was set from kappa_40 = 0.2319276 fm^-1 (from the lambda_40 run). A_S(r) = u e^{kappa_N r}
was evaluated on 81 points in [8, 16] fm, with kappa_N taken from each run. The ODE A_S is 0.9043067
fm^-1/2 (wall-corrected plateau on [20, 30] fm; the raw plateau differs by 2e-7).

| lambda | B_N (MeV) | spread (max-min)/mean | mean A_S | mean vs ODE | gate |
|:--|:--|:--|:--|:--|:--|
| 2 kappa = 4639/10000 | 2.2305354 | **5.6e-3** | 0.904058 | **-2.8e-4** | 1e-4 each |
| 3 kappa = 3479/5000 | 2.2306864 | **8.2e-4** | 0.904343 | +4.0e-5 | 1e-4 each |
| (lambda_40, not scored) | 2.2306877 | 4.2e-6 | 0.9043067 | -5.5e-8 | -- |

A_S lies inside the [0.88, 0.93] window in every run. Flatness fails at both registered lambdas, and the A_S
match fails at 2 kappa. The ODE solution itself is flat to 2.2e-6 on [8, 16] fm, so the 1e-4 gate was not
limited by the physics. The cause is that lambda = 2 or 3 kappa is far from the energy-optimal 6 kappa, so at
N_r = 40 the expansion is still unconverged: B is off by 7e-5 at 2 kappa, and the error leaks into the local tail.
Moving lambda off kappa to avoid the tautology also moved it out of the converged region.
HO control (voided; for information): at N_r = 20, the energy-optimal hbar-omega in the scan is 60 MeV, the
scan's upper edge. There -u'/u deviates from kappa_N by up to 3900% (at 16 fm); the rate is 0.397 at 8 fm and
8.85 at 16 fm. The control behaved as expected (> 10%).

### DR3: FAILED -> MISS (literal clause); the other two clauses HELD
- **Laguerre clause, FAILED.** In the registered orthonormal basis phi_n = N_n f_n,
  N_n^2 = 2 lambda n!/(n+2)!, each entry is (rational)/sqrt((m+1)(m+2)(n+1)(n+2)/(4 lambda^2)).
  Checked exactly at N_r = 40, which by nesting covers N_r = 1..40, at lambda = 13949/10000 and at 1/2.
  1556 of 1600 entries of H are irrational, and 154 of 1600 entries of <r^2> are. The first case is
  <phi_0|H|phi_1>, with S_00 S_11 = 12 at lambda = 1/2. The same irrationality is already visible in the hydrogen
  control: <phi_0|1/r|phi_1> = lambda/sqrt(3). What does hold is that in the unnormalized basis f_n, H_u, S_u and
  <r^2>_u are all exactly rational, with S_u diagonal. So the secular polynomial det(H_u - E S_u) lies in Q[E],
  every truncation's eigenvalues are algebraic, and the matrices are pi-free. The orthonormal entries lie in
  Q(sqrt of a rational). The falsifier as written ("any irrational Laguerre entry") fires.
- **HO Yukawa clause, HELD.** <0s|e^{-mu r}/r|0s> = (2/(sqrt(pi) b))(1 - sqrt(pi) y e^{y^2} erfc(y)), with
  y = mu b/2 (sympy difference 0). The erfc coefficient, -2y/b, is nonzero. Closed form vs quadrature at three
  b values (1.440, 2.036, 2.880 fm), both mu, agree to <= 2e-47.
- **Minnesota clause, HELD.** HO elements n, m <= 2 are polynomials in t = (1 + kappa b^2)^{-1/2}, with
  normalization factors sqrt(5), sqrt(6), sqrt(30). Laguerre <phi_0|e^{-kappa r^2}|phi_0> carries e^{z^2} erfc(z),
  z = lambda/sqrt(kappa) (quadrature check 7e-45).

### DR4: FAILED -> MISS
- r_m(MT-III) = **1.9897135 fm** (basis N_r = 40 and ODE agree to 2e-9). The window [1.99, 2.05] is missed by
  0.00029 fm.
- r_ch = sqrt(r_m^2 + r_p^2 + r_n^2 + r_DF^2) = **2.14089 fm**, which is **+0.606%** over CODATA 2018 2.12799
  (window +0.8% to +3.5%; reaching +0.8% would need r_m = 1.99415 fm). I re-read the inputs: r_p = 0.84075(64) fm,
  CODATA 2022 (NIST CUU "rp"); r_n^2 = -0.1155(17) fm^2, PDG 2026 (pdgLive S017MCR); r_DF^2 = 0.0331 as
  registered. The from-memory values give +0.613%. Against CODATA 2022 r_d = 2.12778(27) fm the excess is +0.616%.
- The basis-error clause HELD: 2.1e-9, 5.0e6 times below the residual |r_m - 1.969|/1.969 = 1.05e-2.
- Cause: the registration's r_m of about 2.017 fm transferred experiment's +0.6% ERE overshoot. MT-III's is 2.005% (DR5).
- Novelty re-search (three searches): Biernat, PhD thesis, arXiv:1110.3180 p. 94, quotes "Malfliet-Tjon's value
  sqrt<r^2> = 1.97625 fm" (<r^2> = (1/4) ∫ r^2 w_0^2, citing MT69). Its parameter set is unverified, probably the
  original MT69 rather than the revised MT I-III, and it is 0.68% from ours. A GitHub page gives 1.98888 fm (not
  citable). It stays PROBABLY KNOWN.

### DR5: FAILED -> MISS
kappa_40 = 0.23192757 fm^-1 and a_t = 5.5132023 fm (ODE; the Kohn route agrees to 2e-9) give rho = 1.879324 fm,
r_ERE = 2.0296066 fm, and A_S(ERE) = 0.906777 against the computed 0.904307. Then
(r_ERE - r_m)/r_m = **+2.00497%**, against the gate (0, 2%). The miss is real at this precision (both routes,
1e-8). The registration's hand value r_ERE ≈ 2.030 was correct; the gate assumed the experimental
overshoot of about 0.6% would carry over.

## 3. Final classes

| ID | Verdict | Class | Margin |
|:--|:--|:--|:--|
| DR1 | FAILED (Laguerre half; HO half held, voided) | MISS | 1.62e-6 vs 1e-6 |
| DR2 | FAILED | MISS | flatness 5.6e-3 and 8.2e-4 vs 1e-4 |
| DR3 | FAILED (literal rationality clause; 2 of 3 clauses held) | MISS | 1556/1600 entries irrational |
| DR4 | FAILED (both value windows; basis-error clause held) | MISS | r_m 1.98971 vs >= 1.99; r_ch +0.61% vs >= +0.8% |
| DR5 | FAILED | MISS | 2.005% vs < 2% |

Unscored observations: no fitted parameter enters; the radius is an output of MT-III's five numbers; (B, a_t)
via shape-independent ERE fix it to 2.0%, not 1%; the Laguerre basis converges r_m to 2e-9 at N_r = 40 while
HO at equal function count stays >= 1% off.

## 4. Deviations

- **D1.** "Assembled in exact rationals" is impossible in the orthonormal basis (DR3). Assembly is exact in the
  unnormalized f_n (S_u diagonal, verified), then mapped to the orthonormal basis at dps 50; spectra unaffected.
- **D2.** No basis method for a_t was registered. I fixed the zero-energy Kohn/Schwinger functional before
  computing: u = r + w, w in span{1 - e^{-lambda r}, phi_n}, a = (<r|V|r> + <w|V|r>)/41.47, N_r = 40, lambda_40.
  It is not rational (<s|V|s> has logs). Its asymptotic coefficient (5.5132040) agrees with a to 3e-7.
- **D3.** Lambda grid k/10000 fm^-1: coarse scan 0.2..3.0, ternary search, grid-local-minimum check;
  kappa-multiples rounded to 4 decimals.
- **D4.** DR2: lambda = m kappa_40; flatness = (max-min)/mean over 81 points; ODE A_S with the r_max = 60 wall
  factor 1 - e^{-2 kappa(60-r)} divided out (raw differs by 2e-7).
- **D5.** ODE: float64 Numerov, 10-term series start, Richardson (h^4, h^6) over h = 2e-3, 1e-3, 5e-4 fm; a_t from
  the linear solution at 50 and 60 fm; nodes counted on (0, 50] fm, since the shooting root leaves a sign flip
  within 2 fm of the wall.
- **D6.** The MT-I control lambda was unregistered: I used the triplet lambda_40 (N_r = 1..40) plus the 7 kappa-multiples.
- **D7.** HO: 160 working digits for assembly, dps 30 for diagonalization. The DR2 energy-optimal hbar-omega is
  the 60 MeV scan edge.
- **D8.** Barletta and Kievsky's "a_s" was read as a_t, as the registration does.
