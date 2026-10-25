# Pre-registered prediction: the production Gaussian fit of the Slater geminal, against an exact correlated integral (GF1-GF3)

**Registered:** 2026-10-25, before any computation in this file. PI direction: "Is there any judgment call we can
make as compared to a gaussian production method?", then, on the proposed test, "Let's make the call to do all 3."
**Prompted by:** records GY (`2026-10-23_gaussian_geminal_triangle.md`) and PF (`2026-10-23_price_formula_slots.md`).
GY showed that a single Gaussian geminal changes the type of the correlated triangle. Production explicitly correlated
codes do not use a single Gaussian. They use the Slater geminal exp(-gamma r12) and replace it, for many integrals,
by a fitted sum of Gaussians. This record measures the error of that replacement on one correlated three-electron
integral whose exact value is now known to about 50 digits.
**Extends:** `2026-10-23_gaussian_geminal_triangle.md`, `2026-10-23_price_formula_slots.md`.
**Commit status:** committed on its own, before any computation in this file.

## The object

The correlated triangle off the symmetric slice, as in record PF:
Y(t) = <rho_i rho_j rho_k f_1(r_ik) f_2(r_jk) / r_ij>, density scales a = 1 (i), b = 2 (j) and t (apex k), at the
fingerprint fibre rho~(k) = c^3 (k^2 + c^2)^(-3/2), with Slater geminals f_1 = exp(-(3/2) r) and f_2 = exp(-3 r).
Conventions of record QG, as in `debug/pf_common.py`.

The fitted object Y_n(t) replaces each geminal by its n-term fit,
  exp(-gamma r) ~ sum_i c_i exp(-alpha_i gamma^2 r^2),
with one fit {c_i, alpha_i} at gamma = 1 rescaled to each gamma. Y_n is bilinear in the two fits:
  Y_n(t) = sum_ij c_i c_j Y_GG(alpha_i (3/2)^2, alpha_j 3^2; t),
where Y_GG(A, B; t) is the same triangle with geminals exp(-A r_ik^2) and exp(-B r_jk^2).

The measured quantity is the relative error E_n(t) = Y_n(t)/Y(t) - 1.

## What I know before computing (disclosed)

- Production practice, checked against the Molpro manual (2010 PDF edition, chapter "Explicitly correlated
  methods"): the correlation factor is "a frozen expansion f12 of Gaussian type geminals. By default the geminal is
  built from six Gaussian functions, and the exponents and coefficients are optimized to obtain the best least
  squares fit to f12 = exp(-beta r12) using a suitable weight function." The manual cites A. J. May and F. R. Manby,
  J. Chem. Phys. 121, 4479 (2004) for the frozen geminal expansions. Its options name a weight exponent GEM_OMEGA
  (default derived from GEM_BETA) and a power of r in the weight, GEM_M (default 0).
- Not verified: the exact production weight function. The May-Manby and Tew-Klopper (J. Chem. Phys. 123, 074101,
  2005) papers could not be opened from this machine. So this record does not claim to reproduce the production
  fit. It uses two stated weights, below, and reports both.
- The weight W1 is plain least squares in three dimensions. That is the STO-NG problem of Hehre, Stewart and Pople
  (J. Chem. Phys. 51, 2657, 1969), which fits the 1s Slater function exp(-r) by Gaussians in the same norm.
- Recalled, not checked: that the best n-term Gaussian approximation of exp(-r) converges like exp(-c sqrt(n)), not
  geometrically, because exp(-r) has a cusp at r = 0 and Gaussians are smooth there (Kutzelnigg,
  Int. J. Quantum Chem. 51, 447, 1994). To be checked against the source before any disposition cites it.
- A Gaussian geminal has the Gaussian transform (pi/A)^(3/2) exp(-q^2/(4A)). The geminal pair at external momentum L
  combines into one Gaussian with mu = AB/(A+B), centre shift c = A/(A+B) and a factor exp(-L^2/(4(A+B))). The
  Coulomb line against one Gaussian gives a Dawson function (record GY, equal exponents). The unequal case is to be
  derived and checked, not taken from memory.
- Nothing about Y_n or E_n has been computed.

## Model

- Fits at gamma = 1, for n = 2 to 10, by nonlinear least squares in the exponents with the coefficients solved
  linearly at each step, at high precision:
  - W1: minimise int_0^oo r^2 [exp(-r) - sum c_i exp(-alpha_i r^2)]^2 dr.
  - W2: the same with the extra weight exp(-r^2), a Gaussian weight of the kind the Molpro options describe. The
    value of its exponent is our choice, not the production one.
- Y_GG(A, B; t) by the momentum route of record GY (route A), generalised to unequal exponents and unequal density
  scales: the loop in closed form with a Dawson function, then the three remaining radial integrals by quadrature in
  float64 or better, aiming at 1e-12 relative.
- The reference Y(t) from the 60-digit grid `debug/data/pf_grid_d60_h0.03125.json`.
- t in {1/2, 1, 2, 4}.

## Predictions

**GF1 (main: the six-term fit is accurate here, but not exact).** For n = 6, under both weights, at each of the
four t: 1e-7 <= |E_6(t)| <= 1e-3.
- Falsifier: |E_6| above 1e-3 at any of the four t under either weight (the fit is worse than expected), or below
  1e-7 at all four under both weights (the fit is better than expected).
- Confidence: moderate. The STO-6G analogy suggests a fit error of order 1e-3 in the norm. A smooth integral should
  do better than that, but by how much is not known.

**GF2 (the cusp sets the rate).** Under W1, at t = 1, for n = 2 to 10, ln|E_n| is fitted better by a + b sqrt(n)
than by a + b n: the residual sum of squares of the sqrt(n) fit is at most a quarter of the linear fit's.
- Falsifier: the linear form fits as well or better, or |E_n| fails to decrease at two or more steps of n.
- Confidence: moderate. It rests on the recalled theorem, which has not been checked, and on Y inheriting the
  geminal's rate.

**GF3 (diffuse apex, larger error).** Under W1, at n = 6, |E_6(1/2)| > |E_6(4)|.
- Falsifier: the reverse.
- Confidence: low. A diffuse apex density reaches large r_ik and r_jk, where a Gaussian sum falls off faster than
  the exponential. A compact apex does not.

## Novelty status

- That six Gaussians fit the Slater geminal well enough for MP2-F12 energies: KNOWN (May-Manby 2004, Tew-Klopper
  2005, as cited by the Molpro manual; not read here).
- The sqrt(n) rate for Gaussian expansions of exp(-r): KNOWN if the recalled reference holds.
- The error of the fit on an exact correlated three-electron integral, off a symmetric configuration: no prior
  measurement found in a brief search. Not a search to the standard of a novelty claim.

## Controls

- Fit code: W1 at n = 3 reproduces the STO-3G exponents of Hehre, Stewart and Pople at zeta = 1 to five digits
  (0.109818, 0.405771, 2.22766, as tabulated by Szabo and Ostlund; to be checked against a source before the run).
- Loop: the unequal-exponent Gaussian loop against direct 3-D quadrature at five triples, to 1e-8.
- Engine: at A = B and a = b = 1, Y_GG against record GY's engine (`debug/gy_common.py`) at three t, to 1e-10.
- Quadrature: Y_GG at two quadrature sizes agree to 1e-12 relative at every pair used; the worst case is reported.
- Reference: the 60-digit grid agrees with the 32-digit grid to 1e-15 at each t used. A t that fails is dropped and
  reported.
- Sanity: E_n must fall as n grows, under each weight. A normalisation error would leave it at order one.

## Caps stated in advance

- One configuration. The densities are model shapes, the fingerprint fibre, not Gaussian orbitals. Production uses
  Gaussian orbitals as well, so this is a comparison on our ground.
- The production weight function is not reproduced; two stated weights stand in for it.
- The fit replaces the geminal in this one integral. Production uses the fit only for some integral classes and
  keeps others exact. An error measured here says nothing directly about an MP2-F12 or CCSD-F12 energy.
- If the engine controls fail, that is a STOP for feasibility, reported and not scored.
