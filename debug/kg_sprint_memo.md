# Sprint KG Memo — Klein-Gordon on S^3 x R: where does pi enter?

Sprint date: 2026-05-03.

## What the sprint tested

Hypothesis: the bare relativistic Klein-Gordon spectrum on S^3 x R lives in
the same algebraic-extension ring Q[sqrt(d_1), sqrt(d_2), ...] as Paper 28's
graph-native QED, with d_i positive square-free integers. pi is NOT in the
spectrum at this stage. pi enters only when an observer compactifies the
temporal direction (finite-time observation window, finite-temperature
partition function, Casimir energy on S^3 x S^1_beta), via the Matsubara
relation omega_k^t = 2 pi k / beta. Mass m is a parameter projection that
preserves the ring (for rational m^2); time-compactification is the
observation projection that brings 2 pi into the eigenvalues themselves.

This sprint runs three computations to support or refute that picture.

## KG-1: Algebraic ring of bare KG spectrum

For omega_n^2 = n(n+2)/R^2 + m^2 on S^3 (R=1), with the rational m^2 panel
{0, 1, 1/4, 2}, n = 1..50:

- **Every omega_n in the rational panel decomposes uniquely as
  omega_n = c_n * sqrt(d_n) with c_n in Q and d_n a positive square-free
  integer.** Verified symbolically: sp.simplify(omega_n^2 - omega_n2) = 0
  for all 50*4 = 200 cases.
- **Zero transcendentals appear.** No pi, no e, no log, no zeta-values in
  any omega_n for any rational m^2 in the panel.
- The sqrt(d) generators that appear, by m^2 case, in n in [1, 50]:

  | m^2 | First few omega_n | First few d_n |
  |-----|-------------------|---------------|
  | 0 | sqrt(3), 2 sqrt(2), sqrt(15), 2 sqrt(6), sqrt(35), 4 sqrt(3) | 3, 2, 15, 6, 35, 3 |
  | 1 | 2, 3, 4, 5, 6, 7 | 1, 1, 1, 1, 1, 1 (degenerate to perfect squares!) |
  | 1/4 | sqrt(13)/2, sqrt(33)/2, sqrt(61)/2, sqrt(97)/2 | 13, 33, 61, 97 |
  | 2 | sqrt(5), sqrt(10), sqrt(17), sqrt(26) | 5, 10, 17, 26 |

  The m^2 = 1 case is special: n(n+2) + 1 = (n+1)^2 is a perfect square, so
  the spectrum collapses to integers omega_n = n+1. This is exactly the
  conformally coupled massless scalar (m^2_eff = 1 from R/6 conformal shift)
  that drives the KG-3 Casimir computation. The other cases scatter d
  generators across many distinct square-free integers in [1, 10401].

- **Full union of square-free generators across the rational panel for
  n in [1, 50]:** 134 distinct d values, ranging from 1 to 10401. See the
  JSON for the complete list.

- **Negative controls (irrational m^2):** Both m^2 = 1/pi and m^2 = sqrt(2)
  break the ring at n = 1 (the first n probed), as expected. Setting m^2
  outside Q immediately injects the corresponding irrational into omega_n^2.
  The mechanism is trivial -- omega_n^2 is a literal sum -- but it confirms
  the boundary: the spectrum is in the algebraic-extension ring iff the
  mass parameter is in it.

**Verdict: hypothesis confirmed for KG-1.** The bare KG spectrum on S^3 x R
with R = 1 and rational m^2 sits in the ring Q[{sqrt(d) : d square-free}],
with no transcendentals.

## KG-2: Where 2 pi enters

Two cases, side by side:

(a) **Open time** (t in R, no compactification). omega_n is just the
spatial spectrum from KG-1. beta is a Boltzmann parameter applied
externally; the temporal direction has no eigenvalues to quantize.
The partition function partial sum (m=0, beta=1, N=30):
Z_partial = 2.94571367781823352341530653269 (50 dps).
The leading Boltzmann factors are exp(-beta sqrt(d)) with d = 0, 3, 2, 15,
6, 35, 3, 7, 5, 11, ... -- pure-square-free arguments. **No pi anywhere
in the spectrum or in the per-mode Boltzmann factor.**

(b) **Periodic Euclidean time** (t ~ t + beta). Promotes the temporal
direction to a circle, so it acquires its own eigenmodes
omega_k^t = 2 pi k / beta, k in Z. The full spectrum is

    omega_{n,k}^2 = n(n+2) + (2 pi k / beta)^2 + m^2.

At beta = 1, m = 0, the very first non-zero-temporal eigenvalue is
omega_{0,1}^2 = 4 pi^2 (the (n=0, k=1) mode). **pi appears at exactly this
step, and not before.** Before/after:

| Step | Eigenvalue object | Transcendental content |
|------|-------------------|------------------------|
| Spatial Laplace-Beltrami on S^3 | -Delta Y_n = n(n+2) Y_n | none (integers) |
| KG dispersion on S^3 x R, t open | omega_n^2 = n(n+2) + m^2 | none for rational m^2 (KG-1) |
| Periodic time t in S^1_beta, k = 0 sector | omega_{n,0}^2 = n(n+2) + m^2 | none |
| Periodic time, k != 0 sector | omega_{n,k}^2 = n(n+2) + (2 pi k / beta)^2 + m^2 | **2 pi appears here** |

First five (n, k) eigenvalues at beta = 1, m = 0, sorted by ascending
omega^2:
  1. (n=0, k=0): omega^2 = 0
  2. (n=1, k=0): omega^2 = 3, omega = sqrt(3) -- still pi-free
  3. (n=2, k=0): omega^2 = 8, omega = 2 sqrt(2) -- still pi-free
  4. (n=3, k=0): omega^2 = 15, omega = sqrt(15) -- still pi-free
  5. (n=0, k=1): omega^2 = 4 pi^2, omega = 2 pi -- **first pi**

The injection is mechanical: the eigenvalue equation
    -d^2/dt^2 phi(t) = lambda^t phi(t)
on a circle of circumference beta has spectrum lambda^t_k = (2 pi k / beta)^2
because the smallest non-trivial wave-number on a circle of circumference
beta is 2 pi / beta. The 2 pi is the ratio (circumference) / (radian
measure). Its appearance is structurally identical to the calibration tier
in Paper 18's exchange-constant taxonomy: it is the Weyl-asymptotic ratio
attached to the compact-direction integral measure.

**Verdict: hypothesis confirmed for KG-2.** pi enters the spectrum at the
temporal-compactification step, as the ratio 2 pi / beta. No earlier step
injects pi.

## KG-3: Casimir as one-projection match

For a conformally coupled massless scalar on unit S^3 (no temporal
compactification, T = 0, "spatial-only Casimir"):

- spatial spectrum: omega_n = n + 1, multiplicity (n+1)^2, n = 0, 1, 2, ...
  (this is m^2 = 1 in the KG-1 panel, confirming the integer collapse)
- spectral zeta: zeta_X(s) = sum_n (n+1)^2 omega_n^{-s} = sum_m m^{2-s}
                = zeta_R(s - 2)
- Casimir energy: E_Cas = (1/2) zeta_X(-1) = (1/2) zeta_R(-3)
- numerical value: zeta_R(-3) = -B_4 / 4 = -(-1/30) / 4 = 1/120
  (verified symbolically via sympy: sp.zeta(-3) = 1/120, sp.bernoulli(4)
   = -1/30)
- **E_Cas = 1/240 exactly. Pure rational. No pi.**

Per-step transcendental ledger:

| Step | Quantity | Transcendentals introduced |
|------|----------|----------------------------|
| Spatial spectrum on S^3 (conformal) | omega_n = n+1, deg = (n+1)^2 | none (integers) |
| Spectral zeta zeta_X(s) = zeta_R(s-2) | analytic continuation of ordinary Riemann zeta | none in form |
| Evaluate at s = -1 -> zeta_R(-3) | rational by Bernoulli formula | NONE (1/120 is rational) |
| Casimir E_Cas = (1/2) zeta_X(-1) | (1/2) * (1/120) = 1/240 | NONE -- pure rational |

Textbook reference: Bytsenko-Cognola-Elizalde-Moretti-Zerbini, Analytic
Aspects of Quantum Fields, World Scientific 2003, sec 4.5; also Dowker &
Critchley, Phys. Rev. D 13, 3224 (1976), and Ford, Phys. Rev. D 11, 3370
(1975). All give E_Cas = 1/240 for the conformally coupled scalar on
unit S^3. Match: exact rational equality, agrees with all decimal places
of the textbook value (60 dps verified numerically).

When the temporal direction is compactified to S^1_beta, the high-temperature
expansion of the free energy density acquires the Stefan-Boltzmann term

    F(beta) / V_3 ~ -(pi^2 / 90) / beta^4 + ...,

where V_3 = 2 pi^2 (volume of unit S^3). The pi^2 / 90 is the familiar
result of the Matsubara sum collapsing to (after the s -> 0 analytic
continuation) zeta_R(4) = pi^4 / 90. This is exactly the "extra pi" the
temporal projection injects -- it lives in the Matsubara measure
(2 pi / beta)^{2s} prefactor times the zeta_R(2s) = (rational) * pi^{2s}
spatial sum, and the net pi power survives to give the SB constant.

**Verdict: spatial-only Casimir on S^3 is rational (one-projection match);
adding S^1_beta brings pi exactly through the KG-2 mechanism. The split is
clean and matches the Paper 34 dictionary.**

## Implications for Paper 34

The current Paper 34 projection list groups all of "introduce mass / observe
at finite time / finite-temperature evaluation" under broadly worded
parameter-projection rows. This sprint shows that they are structurally
different and should be split. Concrete proposal -- two distinct rows in
the projection table:

1. **Rest-mass projection** (introduces m, dimension mass = energy at
   c = 1, transcendental signature: trivial, ring-preserving for rational
   m^2). The bare-graph algebraic-extension ring is closed under this
   projection iff m^2 in Q.

2. **Observation / temporal-window projection** (introduces beta, dimension
   time = inverse energy, transcendental signature: 2 pi * Q per Matsubara
   mode; Stefan-Boltzmann constant pi^2 / 90 in the high-T limit). This
   projection is what "an observer" does when integrating over a finite
   time window or evaluating a partition function -- it lives in Paper 18's
   calibration tier.

This is a no-edits-yet proposal -- PI's call whether to amend Paper 34.
The finding is structural: the existing dictionary is *under-resolved* on
the temporal axis. Mass and observation-window are presently lumped together
as "scalar parameter introduces a constant" but they introduce categorically
different transcendental content (none vs 2 pi).

## Falsifiable prediction

If the sprint hypothesis holds, then any GeoVac observable that does not
require integrating over time should be representable in the bare-graph
algebraic ring without pi. Confirmatory examples already in the framework:

1. **Pauli counts** (Paper 14): integers, no pi. Consistent.
2. **Nuclear magic numbers** (Paper 23, Track NB): integers
   {2, 8, 20, 28, 50, 82, 126}, no pi. Consistent.
3. **Hydrogenic atomic eigenvalues on the bare graph** (Paper 7):
   lambda_n = -(n^2 - 1), pure integers. Consistent.
4. **Bargmann-Segal HO graph** (Paper 24): bit-exactly pi-free in exact
   rational arithmetic at every N_max. Consistent.
5. **Graph-native QED (Papers 28, 33)**: all matrix elements in Q[sqrt(d)],
   no transcendentals at any finite n_max. Consistent.

Counter-examples (places where pi DOES appear in observables, all of which
should be associated with a temporal / Matsubara / continuum-projection
step):

1. Vacuum-polarization Pi = 1/(48 pi^2) (Paper 28): one-loop QED on S^3,
   evaluated as a continuous spectral integral over t in (0, infinity).
   The pi^2 traces to the Seeley-DeWitt coefficient a_2 on S^3 = sqrt(pi)/8,
   squared. The continuous t-integral is the temporal-projection step.
2. QED beta function 2 alpha^2 / (3 pi) (Paper 28): same provenance.
3. Stefan-Boltzmann pi^2 / 90 in any continuum-T thermodynamic quantity.

The prediction, sharpened: **a GeoVac observable contains pi iff its
evaluation includes a continuous integration over either a temporal or a
spectral parameter that has been promoted from the discrete graph spectrum.**
The discrete graph itself is pi-free; the projection that integrates is
where pi appears.

## Honest limits

- KG-3 took the textbook value 1/240 as the reference. The script verifies
  the reduction to (1/2) zeta_R(-3) = 1/240 symbolically (sp.zeta(-3) and
  sp.bernoulli(4)), but the original derivation that the conformally coupled
  scalar's spectrum on S^3 is omega_n = n+1 with degeneracy (n+1)^2 was
  cited rather than re-derived from Lichnerowicz / conformal coupling first
  principles. The standard derivation is in any of the cited references.

- KG-2's explicit "first appearance of pi" at (n=0, k=1) uses the literal
  Matsubara dispersion. It does NOT distinguish whether pi appears at the
  level of (a) the spectrum, (b) the partition function, or (c) the
  free-energy density. The script shows (a) cleanly; (b) and (c) inherit it
  but additional projections (analytic continuation, zeta regularization)
  can reshuffle the pi-content.

- The sprint did not address whether the GeoVac framework's existing
  "pi-free certificate" for the Bargmann-Segal HO and the bare graph
  Laplacian holds for arbitrary temporal compactifications -- it tested
  S^1_beta only.

- No claim was made about the spinor sector. The Dirac-on-S^3 spectrum
  |lambda_n| = n + 3/2 (Camporesi-Higuchi) is half-integer-rational; the
  ring is Q (no sqrt extensions), so the same pi-free statement holds for
  the bare Dirac graph by inspection. But the corresponding KG-3 analog
  (Dirac Casimir on S^3 x S^1) was not computed.

- No new tests were added to tests/. The verification is in the JSON files
  and the symbolic identity assertion in kg1_algebraic_ring.py
  (sp.simplify(omega^2 - omega_n^2) == 0 for all n). Promotion to a tests/
  file would be an appropriate Paper-35-readiness step but was not in scope.

## Sprint outputs

- `debug/kg1_algebraic_ring.py` + `debug/data/kg1_algebraic_ring.json`
  (algebraic-ring decomposition of omega_n for n in [1, 50] across rational
   m^2 panel, plus irrational-m^2 negative controls)
- `debug/kg2_temporal_compactification.py` + `debug/data/kg2_temporal_compactification.json`
  (open vs periodic-time spectra, explicit pi-first-appearance flag)
- `debug/kg3_casimir_s3_s1.py` + `debug/data/kg3_casimir_s3_s1.json`
  (S^3 conformal Casimir = 1/240 + per-step transcendental ledger
   + Matsubara pi-origin diagnostic)
