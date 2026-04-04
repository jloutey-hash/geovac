# Track BI: Level 3 Hyperspherical Validation & Ab Initio PK Parameters

## 1. Level 3 Solver Validation (l_max=2, single-channel adiabatic)

| Z | Ion | E_computed (Ha) | E_exact (Ha) | Error (%) | Notes |
|:-:|:----|:----------------|:-------------|:---------:|:------|
| 2 | He  | -2.9275 | -2.9037 | 0.82 | Reference (Paper 13) |
| 3 | Li+ | -7.3216 | -7.2799 | 0.57 | Good |
| 4 | Be2+| -13.6960 | -13.6556 | 0.30 | Good |
| 5 | B3+ | -22.0067 | -22.0310 | 0.11 | Best accuracy |
| 6 | C4+ | -32.1866 | -32.4062 | 0.68 | Moderate |
| 7 | N5+ | -44.1525 | -44.7814 | 1.40 | Elevated |
| 9 | F7+ | -73.0687 | -75.5316 | 3.26 | High |

### Error analysis

The errors are NOT numerical artifacts — they are structural, from the single-channel
adiabatic approximation at l_max=2. Grid convergence tests (R_max=10-20, N_R=2000-3000,
n_alpha=100-150) confirm errors are stable to <0.01%.

The error pattern is non-monotonic: it decreases from Z=2 to Z=5, then increases
sharply for Z>=6. At Z=5, the nuclear charge is large enough that electron-electron
correlation is a smaller fraction of total energy (favorable), but at Z>=7, the
wavefunction becomes so compact that the l_max=2 angular basis is insufficient to
capture correlation structure at the cusp.

For the PK application (core screening for composed geometries), the core density
shape matters more than the absolute energy. The core eigenvector at all Z values
is >99.8% l=0, confirming that the 1s^2 core is well-described even at moderate l_max.

## 2. Ab Initio PK Parameters

All parameters computed with l_max=2, algebraic density extraction, inv2 B-method.

| Z | A (Ha*bohr^2) | B (1/bohr^2) | r_core (bohr) | E_core/e (Ha) | E_val (Ha) |
|:-:|:-------------|:-------------|:-------------|:-------------|:-----------|
| 3 | 6.927 | 6.795 | 0.271 | -3.662 | -0.198 |
| 4 | 13.014 | 12.102 | 0.203 | -6.850 | -0.343 |
| 5 | 21.402 | 18.460 | 0.165 | -11.006 | -0.305 |
| 6 | 31.366 | 25.539 | 0.140 | -16.097 | -0.414 |
| 7 | 43.093 | 33.050 | 0.123 | -22.081 | -0.534 |
| 9 | 71.804 | 48.609 | 0.101 | -36.542 | -0.640 |

### B-method candidates (all Z values)

| Z | B_rms | B_peak | B_median | B_inv2 (used) |
|:-:|:-----:|:------:|:--------:|:-------------:|
| 3 | 1.090 | 0.750 | 1.992 | 6.795 |
| 4 | 2.072 | 1.469 | 3.727 | 12.102 |
| 5 | 3.333 | 2.880 | 5.948 | 18.460 |
| 6 | 4.843 | 3.403 | 8.542 | 25.539 |
| 7 | 6.546 | 4.986 | 11.385 | 33.050 |
| 9 | 10.303 | 8.000 | 17.686 | 48.609 |

Note: Z=4 Laguerre fit fell back to spline Z_eff (monotonicity validation failed).
All other Z values used spectral Laguerre successfully.

## 3. Z^2 Scaling Analysis (anchored to Li, Z=3)

| Z | A_ab_initio | A_Z2 | A_err(%) | B_ab_initio | B_Z2 | B_err(%) |
|:-:|:-----------|:-----|:--------:|:-----------|:-----|:--------:|
| 3 | 6.927 | 6.927 | 0.0 | 6.795 | 6.795 | 0.0 |
| 4 | 13.014 | 12.315 | 5.4 | 12.102 | 12.081 | 0.2 |
| 5 | 21.402 | 19.242 | 10.1 | 18.460 | 18.876 | 2.3 |
| 6 | 31.366 | 27.709 | 11.7 | 25.539 | 27.181 | 6.4 |
| 7 | 43.093 | 37.715 | 12.5 | 33.050 | 36.997 | 11.9 |
| 9 | 71.804 | 62.344 | 13.2 | 48.609 | 61.158 | 25.8 |

### Assessment

Z^2 scaling is a **poor approximation** for both A and B parameters:
- A: 5-13% error, systematically underestimates (actual grows faster than Z^2)
- B: 0.2-26% error, with B growing slower than Z^2 at high Z

The actual scaling is approximately:
- A ~ Z^{2.25} (faster than Z^2 due to increasing core-valence energy gap)
- B ~ Z^{1.75} (slower than Z^2 because r_core ~ 1/Z, giving B ~ Z^2 only if
  the density shape were Z-independent, which it isn't)

**Conclusion: Ab initio PK parameters are required for each Z. Z^2 scaling from
Li is not accurate enough for production use, especially for B at Z>=7.**
