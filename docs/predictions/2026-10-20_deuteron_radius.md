# Pre-registered prediction: the deuteron's size from its own binding, in an exponential-tail basis (DR1-DR5)

**Registered:** 2026-10-20, before any computation in this file (PI-directed: "Yep, do it"). Drafted by a
scoping agent, reviewed by the PM. Prompted by the PI's
direction to let the nuclear skeleton supply the radius that the finite-nucleus doors consume, after
`debug/nuclear_radius_scan_memo.md` found that the HO misses <r^4> of 4He by 5.6% because of its Gaussian
tail.
**Extends:** Paper 23 `sec:deuteron` and `sec:finite_nucleus`; seam rows J14, J16 (the radius and the Zemach
moment consumed); atlas doors 24 and 26; Paper 24 `sec:ks_frame` (J15: the Coulomb Sturmians as the pushdown
of the four-dimensional oscillator). Scoping: `debug/deuteron_sturmian_scoping_memo.md`.
**Commit status:** committed on its own, before any computation that scores it.

## The question

The finite-nucleus doors consume the nuclear radius R. For the lightest nucleus, can R be an **output**?
Concretely: solve the np ³S₁ bound state for a published potential with fixed parameters, in a basis whose
tail is exponential, and compute the binding B, the point-matter radius r_m, the tail, and the charge radius.
The atlas reading to test is that the free number **moves** from R to the potential's parameters (and,
for the charge radius, the nucleon radii) and does not disappear. A sharper version: at low energy it moves
to two numbers, B and the triplet scattering length a_t.

The radius is a ground-state expectation value, <r^2> = <psi_0|r^2|psi_0>, and needs no sum over excited
or continuum states. So the failure recorded in `docs/failed_approaches_ledger.md` ("Sturmian basis for
deuteron polarizability with Minnesota NN", sum-over-states fragile for continuum response) does not
apply to this observable.

## What I know before computing (disclosed)

- **The potential.** MT-III (revised Malfliet-Tjon I-III, triplet), seen in primary text:
  V(r) = (1438.72 e^{-3.11 r} - 626.885 e^{-1.55 r})/r, with V in MeV and r in fm; hbar^2/m = 41.47 MeV fm^2;
  l = 0 only; E_2b = -2.23069 MeV and a_s = 5.5132 fm (Barletta and Kievsky, arXiv:0807.2333, Sec. IV).
  Garrido et al., arXiv:1211.1802, Sec. III: the same form, "binding energy for the deuteron of 2.2307 MeV",
  attributed to Friar et al., PRC 42, 1838 (1990). Malfliet and Tjon, Nucl. Phys. A127, 161 (1969):
  not consulted. **No published MT-III deuteron radius or A_S was found** (three web searches).
- **Prior art on the basis.** Caprio, Maris and Vary, PRC 86, 034312 (2012) and PRC 90, 034305 (2014)
  (arXiv:1208.4156, 1409.0877). Their function is S_nl(b; r) ∝ (2r/b)^{l+1} L_n^{2l+2}(2r/b) e^{-r/b}. For
  4He the CS and HO radii agree (1.4358 vs 1.4361 fm), with no gain. For 6He and 8He the CS basis gives
  "improved hbar-Omega independence" and no definitive r_m. No deuteron. More, Ekström, Furnstahl, Hagen
  and Papenbrock, PRC 87, 044326 (2013) (arXiv:1302.3815, Fig. 5): the deuteron r^2 in an HO basis is
  IR-limited, following a single curve in L_2 = sqrt(2(N + 3/2 + 2)) b.
- **Radii.** Buchmann, Henning and Sauer, nucl-th/9610018, give Eq. (3)
  r_d^2 = (r_ch^2 - r_p^2) - r_n^2 - r_DF^2 with r_DF^2 = 0.0331 fm^2, and Eq. (5), which adds the spin-orbit
  and two-body terms (r_[2]^2 ≈ 0.02 fm^2). They report r_m = 1.969(4) fm from the empirical linear relation
  at a_t = 5.419(7) fm, A_S = 0.885(2) fm^{-1/2}, and r_m ≈ 1.950 fm from scattering analyses.
  r_ch = 2.12799(74) fm (CODATA 2018), seen via a search summary. r_p = 0.8414 fm and r_n^2 = -0.116 fm^2
  are from memory and will be re-read from CODATA 2022 and PDG at scoring.
- **Hand arithmetic (no code).** kappa = sqrt(B/41.47): MT-III 0.23193 fm^-1. Shape-independent
  rho = 2(kappa - 1/a_t)/kappa^2 = 1.880 fm; A_S = sqrt(2 kappa/(1 - kappa rho)) ≈ 0.907 fm^{-1/2};
  r_ERE^2 = 1/(8 kappa^2 (1 - kappa rho)), so r_ERE ≈ 2.030 fm. The same formula on experiment
  (kappa = 0.23161, a_t = 5.419) gives 1.982 fm against 1.969, an overshoot of +0.6%. Applying it here gives
  r_m(MT-III) ≈ 2.017 fm and r_ch ≈ 2.166 fm (+1.8%).
- **Structure.** In the Laguerre basis every matrix element reduces to ∫x^k e^{-ax} dx = k!/a^{k+1}. In the
  HO basis the Yukawa element reduces to Gaussian-times-exponential integrals, which give e^{y^2} erfc(y).
- Nothing in this file has been computed.

## Model

Relative ³S₁, l = 0, hbar^2/m = 41.47 MeV fm^2, MT-III exactly as above. No parameter is fitted.

- **Laguerre basis (primary).** phi_n(r) = N_n x L_n^{(2)}(x) e^{-x/2}, x = 2 lambda r, n = 0..N_r-1,
  orthonormal in dr. The matrices are assembled in exact rationals (`fractions.Fraction`) for rational
  lambda and diagonalized with mpmath at dps 50. Scale rule: lambda_N minimizes E_N at each N_r, searched
  over a rational grid; this is the variational rule, with no physics input. Also reported: the scan
  lambda ∈ {0.5, 1, 2, 3, 4, 6, 8} x kappa_N, rounded to 4-digit rationals.
- **HO basis (comparison).** The same potential and the same function count N_r; hbar-omega ∈ {2, 4, 6, 8,
  10, 15, 20, 30, 40, 60} MeV, b = sqrt(41.47 x 2/hbar-omega) fm. Yukawa elements in closed form at dps 30,
  cross-checked by quadrature to 1e-20.
- **Reference (independent route).** Radial ODE, starting from the series u = r + c_2 r^2 at the origin;
  Numerov with Richardson extrapolation in h; r_max = 60 fm; shooting on B. It yields B, r_m, A_S (the
  plateau of u e^{kappa r} at r ≥ 20 fm) and a_t (the zero-energy solution).
- **Computed.** B_N; r_m = (1/2) sqrt(<r^2>); kappa_N = sqrt(B_N/41.47); A_S(r) = u(r) e^{kappa_N r} on
  r ∈ [8, 16] fm; r_ch = sqrt(r_m^2 + r_p^2 + r_n^2 + r_DF^2) (the spin-orbit and two-body terms are omitted,
  since MT-III has no model for them); r_ERE from the computed kappa and computed a_t.

## Predictions

**DR1 (convergence: the exponential tail wins for the weakly bound case).** At N_r = 30 and lambda = lambda_N,
|r_m(N_r) - r_m(ODE)|/r_m(ODE) ≤ 1e-6. At the same N_r, the HO r_m is off by ≥ 1e-4 relative at **every**
hbar-omega in the scan.
- Falsifier: the Laguerre error > 1e-6, or any HO hbar-omega within 1e-4. If only the HO half fails,
  record it as PARTLY.
- Novelty: the HO half is KNOWN (More et al. 2013, Fig. 5, for a different potential). The Laguerre-vs-HO
  deuteron comparison was NOT FOUND (searched Caprio 2012/2014, More 2013 and two web searches). It is low
  weight, since the expectation is standard. To be re-searched at scoring.

**DR2 (the tail rate is set by the computed binding).** At lambda = 2 kappa_N and 3 kappa_N (deliberately
not kappa_N, where the finite expansion decays as r^N e^{-kappa r} by construction) and N_r = 40, the
quantity A_S(r) = u(r) e^{kappa_N r} with kappa_N = sqrt(B_N/41.47) varies by ≤ 1e-4 relative over
r ∈ [8, 16] fm, and its mean equals the ODE A_S to 1e-4. A_S ∈ [0.88, 0.93] fm^{-1/2}. Negative control: at
N_r = 20 and the energy-optimal hbar-omega, the HO local rate -u'/u deviates from kappa_N by > 10% somewhere in
[8, 16] fm.
- Falsifier: flatness or the A_S match outside 1e-4, or A_S outside the window. The HO control failing does
  not refute DR2, but is reported.
- Novelty: PROBABLY KNOWN. Exponential bases reproducing the asymptotic normalization is the stated
  motivation of Caprio et al. 2012; no deuteron check was found.

**DR3 (the number class is set by the basis-potential pairing).** For rational lambda, every entry of
H_{N_r} (MeV) and of the <r^2> matrix (fm^2) in the Laguerre basis is rational, at N_r = 1..40, checked
exactly. The HO Yukawa element <0s|e^{-mu r}/r|0s> has a closed form whose e^{y^2} erfc(y) term (y = mu b/2)
has a nonzero coefficient, verified numerically to 1e-30 at three b. For the 2x2 table, the Minnesota
Gaussian is algebraic in HO (powers of (1 + kappa b^2)^{-1/2}) and carries e^{z^2}erfc(z) in Laguerre.
- Falsifier: any irrational Laguerre entry; a vanishing erfc coefficient; Minnesota-in-HO not algebraic.
- Novelty: KNOWN as mathematics (elementary integrals). The reading, that matching the basis tail to the
  potential's form keeps every truncation rational and π-free, was NOT FOUND as a statement. It is low
  weight: a hit is a CHECK on the math and a small atlas item.

**DR4 (the radius, and whose error the residual is).** r_m(MT-III) ∈ [1.99, 2.05] fm. r_ch(MT-III) exceeds
CODATA 2018 2.12799 fm by +0.8% to +3.5%. The basis error at N_r = 40 (vs ODE) is ≤ 1e-5 relative, at least
100x below |r_m(MT-III) - 1.969 fm|/1.969: the residual to experiment is the potential's, not the
basis's.
- Falsifier: either value outside its window, or a basis error not 100x below the residual.
- Novelty: PROBABLY KNOWN. The MT-III deuteron radius is likely tabulated in the Friar-Payne-Gibson
  benchmark literature, which was not consulted; it was not found in three searches. A hit is a CHECK.

**DR5 (the forced/free split: the radius's free content is two low-energy numbers).** With kappa_N and
a_t computed from MT-III (no experimental input), 0 < (r_ERE - r_m)/r_m < 2%, where
r_ERE^2 = 1/(8 kappa^2 (1 - kappa rho)) and rho = 2(kappa - 1/a_t)/kappa^2.
- Falsifier: the ratio outside (0, 2%).
- Novelty: KNOWN in form (Bethe 1949 effective-range theory; the empirical a_t-A_S-r_m linear relation of
  Klarsfeld et al. 1986, as reported by Buchmann et al.; Phillips, Rupak and Savage, nucl-th/9908054: at NLO
  the Z-parametrization reproduces effective-range theory for point nucleons). A hit is a CHECK. Its use
  here is the atlas price: R is released, and (B, a_t) are consumed, encoded in the potential's five
  numbers.

## Controls

- **Positive, KNOWN:** B = 2.23069 MeV and a_t = 5.5132 fm (Barletta and Kievsky) reproduced to 1e-5 relative,
  by the basis route and by the ODE route.
- **Positive, exact:** the hydrogen limit in the same Laguerre code (V = -1/r, atomic units, lambda = 1)
  gives E = -1/2 exactly at N_r = 1 (phi_0 ∝ r e^{-r} is the 1s), and at lambda = 1/2 the 2s level
  -1/8 exactly at N_r = 2. This checks the kinetic and 1/r elements.
- **Negative, KNOWN:** MT-I singlet (attractive strength 513.968, Garrido et al.) has no bound state. The
  lowest Laguerre eigenvalue is positive at every N_r ≤ 40 and decreases toward 0.
- **Independent route:** every scored basis number is compared with the ODE route (independent-route rule).
- **Consistency:** the HO energies at large N_r on the L_2 curve, as More et al. describe.

## Result

(To be appended after scoring; the predictions above are never edited.)
