# Sprint memo — Two-center decompactification as a continuous function of R

**Date:** 2026-09-05 · **Type:** DIAGNOSTIC (measurement only) · **System:** H₂⁺ (Z=1/center)
**Driver:** `debug/decompactification_R_sweep.py` · **Data:** `debug/data/decompactification_R_sweep.json`

## Verdict: **BORDERLINE**

The compact→non-compact transition at two centers is **smooth and continuous in R** —
every measured curve varies without a discontinuity or threshold, confirming the
"Boolean label / continuous coupling" reading. **But** the per-shell decompactification
front scales as **n\*(R) ∝ R^1.0**, not the predicted **n ∝ √(Z·R)** (R^0.5). The front
tracks the orbital **exponential decay length** ℓ_n = n/Z, **not** the mean radius
r_n = n²/Z the prediction assumed. Per the decision gate this is BORDERLINE:
*smooth, but window scaling ≠ √(Z·R); actual fitted exponent reported below.*

## Question and reading under test

Paper 0 "Level 2" + Paper 8: at two centers axial **m stays exactly conserved at all R**
(the φ-equation is untouched), while **ℓ stops being an exact quantum number** — it is
replaced by a separation parameter from the η-equation that reduces to ℓ only in the
united-atom limit. The reading tested: the **label** is Boolean (m exact / ℓ lost), but the
**coupling** that ℓ used to separate turns on **continuously** with R, with a per-shell
front where genuine two-center character peaks at r_n ≈ R (shells with r_n ≫ R
united-atom-like, r_n ≪ R atomic).

## Measures (stated explicitly)

- **Per-shell character (D1):** bicentric delocalization made precise as the commutator of
  the two single-center shell projectors,
  char_n(R) = ‖[P_A^(n),P_B^(n)]‖ = S_n(R)·√(1−S_n(R)²) = ½|sin 2θ_n|,
  where S_n(R) = ⟨ns_A|ns_B⟩ is the **exact** two-center overlap of two proper hydrogenic
  ns orbitals (n−1 radial nodes, mean radius ∝ n²/Z) and θ_n = arccos S_n. char_n = 0 at
  θ=0 (S=1, united, R→0) **and** θ=90° (S=0, atomic, R→∞); it **peaks at θ=45° (S=1/√2)**.
  R\*(n) ≡ R where S_n = 1/√2 is the shell's decompactification point (unique — S_n is
  monotone decreasing). This is corpus-native: the same cross-center-overlap → principal-angle
  → commutator object as the composition wall (v4.73.0) and the bond sphere (Paper 8).
- **γ(R) (D2a):** Paper 8 closed form, p_R=1/R, p₀=Z=1; cos γ=(p₀²−p_R²)/(p₀²+p_R²).
- **‖[P_A,P_B]‖ + principal angles (D2b):** SVD of the cross-center overlap in the
  orthonormalized frame; H₂⁺ s-block {1s,2s,3s} sweep + LiH σ-block {1s,2s,2p0} anchor.
- **aggregate l-mixing (D2c):** Σ_{n=1..6} char_n(R) (integral of D1).
- **seed a(R) (D2d):** Mulliken p = (R/2)(ζ_A+ζ_B), ζ=1 ⇒ a=R; the argument of the
  exponential-integral seed e^a·E₁(a) that `geovac/neumann_vee.py` (via `scipy.special.exp1`)
  carries in the closed-form two-center engine.
- **BONUS — direct η-equation l-mixing:** ground-state η-eigenvector of the homonuclear
  spheroidal angular equation (Paper 11 eq:eta_equation) in the normalized associated-Legendre
  basis, at c²(R) from the exact prolate solver (`geovac.prolate_spheroidal_lattice`).
  Metric = participation deficit 1−max_l a_l² and Shannon entropy of {a_l²}. This is the most
  literal realization of Paper 0's "ℓ replaced by a separation parameter." My η-eigenproblem
  reproduces the solver's separation constant A(c²) bit-for-bit (cross-check in driver).

## Results

**D2 — four curves, all continuous (no discontinuity anywhere):**

| curve | shape vs R | endpoints / notes |
|:--|:--|:--|
| (a) γ(R) | monotone ↓, smooth | 157°(R=0.2) → 9.5°(R=12); →π at R→0, →0 at R→∞ ✓ Paper 8 limits |
| (b) ‖[P_A,P_B]‖ H₂⁺ | non-monotone, **capped 0.500** | norm saturates at the ½ ceiling (aha_t1); **principal angles** sweep smoothly 0.8°/2.1°/7.4° (R=0.2) → 44.6°/88.2°/88.6° (R=12), monotone, max adjacent step 4.6° (grid) |
| (c) Σ char_n | rises then falls, one peak | smooth; peak R≈7.5 is N_shell-truncation-dependent |
| (d) a(R)=R ; e^a E₁(a) | linear ↑ ; monotone ↓ | smooth, analytic in R |
| BONUS η l-mixing | ↑ then saturates ~0.5 | part-deficit 0.0000(c²=0.04) → 0.50(c²≈25); the late dip is metric saturation near maximal mixing, not un-mixing; **continuous throughout** |

**LiH σ-block anchor {1s,2s,2p0}, Z=1, R=3.015 (v5.1.0 datum reproduced exactly):**
σ = 0.9913 / 0.7106 / 0.3856 → angles **7.6° / 44.7° / 67.3°**, ‖[P_A,P_B]‖ = **0.500**.

**D1 — per-shell front** (R\* where S_n = 1/√2, peak character):

| n | R\*(n) | r_n=n²/Z | R\*/r_n | ℓ_n=n/Z | **R\*/ℓ_n** |
|--:|--:|--:|--:|--:|--:|
| 1 | 1.565 | 1 | 1.565 | 1 | 1.565 |
| 2 | 4.388 | 4 | 1.097 | 2 | **2.194** |
| 3 | 6.568 | 9 | 0.730 | 3 | **2.189** |
| 4 | 8.520 | 16 | 0.533 | 4 | **2.130** |
| 5 | 10.613 | 25 | 0.425 | 5 | **2.123** |
| 6 | 12.818 | 36 | 0.356 | 6 | **2.136** |
| 7 | 15.015 | 49 | 0.306 | 7 | 2.145 |
| 8 | 17.154 | 64 | 0.268 | 8 | 2.144 |

- Fit **R\*(n) ∝ n^p: p = 0.980 (R²=0.9997, n≥2)**; all-shell p=1.114 (n=1 small-n outlier).
- ⇒ **n\*(R) ∝ R^1.02** — predicted √(ZR) is 0.5; decay-length n/Z is 1.0.
- **R\*/ℓ_n is constant ≈ 2.14** for n≥2 (drift < 3%), while R\*/r_n drifts monotonically
  1.10 → 0.27. The front is pinned to the decay length, **not** the mean radius.

**Direction of the front is as predicted qualitatively:** small R decompactifies the core
(low n) first; the window sweeps up through shells as R grows (core decompactifies last).
Only the *rate* differs — linear in n, not quadratic.

## Independent-route cross-check (memory rule: load-bearing measured numbers)

- Two overlap engines — Mulliken/Ruedenberg closed form (`aha_t1_core.two_center_s_overlap`)
  vs grid quadrature (`fast_two_center_overlap.overlap_fast`) — agree to **< 1.3e-8** across
  n=1..4, R=1..10.
- Front exponent from the independent grid route: **p = 0.968** (vs 0.980 closed-form). The
  exponent ≈ 1 is route-independent.
- η-eigenproblem A(c²) matches the prolate solver's Brent-matched A bit-for-bit.

## Interpretation

The Boolean/continuous split holds cleanly: **m is exactly conserved at all R** (structural,
not measured here — the φ-integral gives δ_{mm'}); **every coupling amplitude is a smooth
function of R** with no switch. Two independent realizations of "ℓ turns off continuously"
agree — the corpus-native cross-center commutator/principal-angles (D2b), and the direct
η-equation l-mixing (BONUS) — both rise smoothly from the united-atom limit. So the
compact→non-compact transition is a **continuous crossover**, not a Boolean event.

The front exists and is sharp per shell, but its scaling refutes the mean-radius premise:
what governs whether two copies of shell n on opposite centers "meet" is their **tail reach**
(exponential decay length ℓ_n = n/Z), not their mean radius r_n = n²/Z. The overlap of two
ns orbitals is dominated by e^{−ζ_n R} with ζ_n = Z/n, so S_n = 1/√2 at ζ_n R ≈ const ⇒
R\* ≈ 2.14 n/Z. The r_n ≈ R picture conflated two different orbital length scales.

## Caveats / not computed

- **Basis choice governs the scaling.** The front is measured in the per-center hydrogenic
  Z=1 basis (the honest basis for bicentric overlap). A measure keyed to real-space orbital
  *size* (⟨r⟩_n = R) would by construction give p=2 — but that is the prediction restated,
  not an independent measurement. The genuine two-center object (overlap) gives p≈1.
- **‖[P_A,P_B]‖ norm carries no scaling** — it saturates at the ½ ceiling (known, aha_t1
  2026-07). The principal *angles* carry the continuous signal; report those, not the norm.
- **η participation-deficit saturates** near 0.5 at large c (bounded metric); entropy (stored)
  keeps rising. Neither changes the continuity conclusion.
- **m-exactness asserted, not re-measured** (φ-integral δ_{mm'}); it is the Boolean half of the
  reading and is structural (Paper 0).
- **Aggregate (D2c) peak location is truncation-dependent** (N_shell=6); only its smoothness
  is load-bearing.
- Single electron (H₂⁺) by design — no e–e confound. Whether the n-vs-n² front distinction
  survives with correlation is untested and out of scope.
- **Guardrail (§3.5) acknowledged:** touches two-center molecular encoding (Papers 8–9). This
  is measurement of the R-dependence that theorem describes, using genuine two-center
  machinery — **not** a proposal to encode a molecule in a single-center/shared-p₀/Sturmian
  basis to bind. No binding claim made.
