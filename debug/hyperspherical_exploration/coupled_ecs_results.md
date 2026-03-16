# Coupled-Channel ECS Results

**Date:** 2026-03-15
**Status:** Complete (20/20 tests passing, important structural result)
**Tests:** `tests/test_complex_scaling.py` -- TestCoupledECS (8 new tests, 20 total)

---

## 1. The Question

Can coupled-channel ECS extract the 2s^2/2s3s width ratio (experimental: 3.29) as Im(E_1)/Im(E_2)?

---

## 2. Coupled-Channel Infrastructure

### Inputs
- Adiabatic curves mu_ch(R) from `compute_adiabatic_curve` (5 channels, l_max=2, n_alpha=100)
- Non-adiabatic coupling P_mu_nu(R) from `compute_coupling_matrices` (Hellmann-Feynman)
- All interpolated as CubicSplines on a 198-point R grid (R = 0.3--20 bohr)

### Key Coupling Strengths

| Pair | max |P_mu_nu| | Physical meaning |
|------|:-----------:|:------|
| P_01 | 0.023 | s-wave channel coupling (autoionization pathway) |
| P_02 | 0.290 | s-p coupling |
| P_03 | 0.312 | s-p/d coupling |
| P_12 | 0.019 | second s-p coupling |
| P_23 | 0.961 | strong p-p coupling (nearly degenerate channels) |

### Channel Thresholds (asymptotic V_eff)

| Channel | V_eff(R=30) | Converging to | Character |
|:-------:|:-----------:|:-------------:|:---------:|
| 0 | -1.72 | -2.0 (He+ n=1) | s-wave |
| 1 | -1.72 | -2.0 (He+ n=1) | s-wave |
| 2 | -0.54 | -0.5 (He+ n=2) | p-wave |
| 3 | -0.54 | -0.5 (He+ n=2) | p/d-wave |
| 4 | -0.51 | -0.5 (He+ n=2) | d-wave |

---

## 3. Theta-Stable Resonance Candidates

### 3-Channel Results (l_max=2, R_0=12, N_R=1500)

| Feature | E (Ha) | eV above gs | Gamma (eV) | Re spread | Dominant ch | Nature |
|---------|:------:|:-----------:|:----------:|:---------:|:-----------:|:------:|
| A | -0.806 | 57.1 | 0.033 | < 10^-5 Ha | ch2 (99%) | Quasi-bound in ch2 |
| B | -0.597 | 62.8 | 0.14 | 0.002 Ha | ch2 (100%) | Quasi-bound in ch2 |

### 5-Channel Results (l_max=2, R_0=12, N_R=1200)

| Feature | E (Ha) | eV above gs | Gamma (eV) | Re spread | Dominant ch | Nature |
|---------|:------:|:-----------:|:----------:|:---------:|:-----------:|:------:|
| A' | -0.784 | 57.7 | -0.84* | 0.002 Ha | ch2 (79%) | See note |
| B' | -0.637 | 61.7 | 0.032 | 2x10^-4 Ha | ch3 (56%) | Quasi-bound in ch3 |
| C' | -0.580 | 63.2 | 0.071 | 0.002 Ha | ch4 (72%) | Quasi-bound in ch4 |

*Positive imaginary part (Im > 0) -- not a physical resonance. Likely complex conjugate partner or numerical artifact.

### Comparison with Experiment

| Resonance | E_expt (eV) | E_calc (eV) | Gamma_expt (eV) | Gamma_calc (eV) | Notes |
|-----------|:-----------:|:-----------:|:---------------:|:---------------:|:------|
| 2s^2 1S | 57.82 | 57.1 (3ch) | 0.138 | 0.033 | 4x too narrow |
| 2s3s 1S | 62.94 | 62.8 (3ch) | 0.042 | 0.14 | 3x too wide |
| 2p^2 1D | 59.90 | 61.7 (5ch) | 0.070 | 0.032 | 2x too narrow |

### Within-Sector Width Ratio

| Method | Gamma(2s^2)/Gamma(2s3s) | Experimental |
|--------|:-----------------------:|:------------:|
| 3-channel ECS | 0.033/0.14 = 0.24 | 3.29 |
| Gaunt graph | 1.13 | 3.29 |
| Avoided crossing (LZ) | fails | 3.29 |

**The ratio is inverted.** The calculation gives Gamma_A < Gamma_B, while experiment gives Gamma(2s^2) > Gamma(2s3s).

---

## 4. Why the Widths Are Wrong

### The resonance candidates are NOT the 2s^2 autoionization resonances

The theta-stable features are **quasi-bound states of the upper channels (2-4)** that sit below their own channel thresholds (~-0.5 Ha) but above the He+(n=1) threshold (-2.0 Ha). Their nature:

1. **Feature A at -0.806 Ha (ch2):** Bound state of channel 2 potential well (V_min = -1.01 Ha, threshold = -0.54 Ha). It autoionizes into channels 0,1 via P_02 coupling (0.29), giving Gamma ~ 0.033 eV.

2. **Feature B at -0.597 Ha (ch2):** Second bound state of channel 2. Higher in the well, broader coupling to continuum.

3. **Feature B' at -0.637 Ha (ch3):** Bound state of channel 3 potential well.

### These are real physics -- but not the right states

These quasi-bound states DO autoionize (they have theta-stable complex eigenvalues with Gamma > 0). But they are dominated by channels 2-4 (p-wave/d-wave character), not channels 0-1 (s-wave). The experimental 2s^2 1S resonance is a pure s-wave (L=0) state involving TWO s-electrons.

### Why can't we resolve the s-wave resonances?

The s-wave channels 0 and 1 show only **broad continuum states** (Gamma ~ 13-37 eV) in the doubly-excited region. The narrow s-wave autoionization resonances are unresolved because:

1. **Weak coupling:** P_01 = 0.023 gives Gamma ~ P^2 ~ 5x10^-4 Ha ~ 0.014 eV. This is 100x smaller than the continuum eigenvalue spacing.

2. **Channel convergence:** Both s-wave channels converge to the same He+(n=1) threshold. The second s-wave channel does NOT correlate to He+(n=2) -- it's the second adiabatic eigenstate with He+(n=1) character.

3. **Angular resolution:** With l_max=2 and n_alpha=100, the angular basis may not resolve the 2s^2 configuration accurately. The doubly-excited state involves subtle angular correlation that requires higher angular moments.

4. **Adiabatic approximation:** The adiabatic channels change character at avoided crossings. The "2s^2" state migrates between channels as a function of R. The ECS needs to capture this migration, which requires very fine R resolution near the crossing.

---

## 5. Channel Decomposition (Eigenvector Analysis)

The coupled ECS solver now returns channel norms -- the fraction of each eigenvector in each channel block. Key findings:

| Energy region | Dominant channels | Physical interpretation |
|:-------------|:-----------------:|:----------------------|
| E < -2.0 (bound) | ch0 (100%) | He ground state |
| -1.0 < E < -0.7 | ch0, ch1 (>90%) | s-wave continuum (discretized) |
| -0.7 < E < -0.5 | ch2, ch3, ch4 | Upper-channel quasi-bound states |
| E > -0.5 | mixed | Full continuum of all channels |

The narrow theta-stable features are ALL in the upper channels. The s-wave channels contribute only broad continuum states in the resonance region.

---

## 6. Convergence Study

### N_R dependence (3 channels, theta=0.3)

| N_R | Feature A Re(E) | Feature A Gamma (eV) |
|:---:|:---------------:|:-------------------:|
| 1000 | -0.8057 | 0.033 |
| 1500 | -0.8057 | 0.033 |
| 2000 | -0.8057 | 0.033 |

Feature A is fully converged in N_R. The position and width are determined by the channel 2 potential well, not the grid resolution.

### n_channels dependence (N_R=1200, theta=0.3)

| n_ch | Feature near 57 eV | Gamma (eV) | Dom. ch |
|:----:|:------------------:|:----------:|:-------:|
| 3 | -0.806 | 0.033 | ch2 |
| 5 | -0.784 (+Im!) | -0.84 | ch2 |

Going from 3 to 5 channels changes the spectrum significantly. The most stable 3-channel feature becomes unstable with 5 channels. This indicates the calculation is NOT converged in channel number.

### R_0 dependence (5 channels, theta=0.3)

| R_0 | Feature B' Re(E) | Gamma (eV) |
|:---:|:----------------:|:----------:|
| 12.0 | -0.637 | 0.032 |
| 15.0 | -0.637 | 0.017 |

Width depends on R_0 -- a sign that the feature is not a true resonance but is influenced by where the complex scaling begins.

---

## 7. Assessment: Does Coupled-Channel ECS Work?

### What works

1. **Infrastructure validated:** Ground state real (Im = 0 to 10^-15), theta-stable, matches exact energy to 1%
2. **Continuum rotates correctly:** -2theta angle after threshold separation
3. **Upper-channel quasi-bound states resolved:** Theta-stable complex eigenvalues with small positive widths
4. **Positions reasonable:** 57-63 eV above ground state, matching the doubly-excited energy region

### What doesn't work

1. **Cannot resolve s-wave autoionization:** The 2s^2 and 2s3s resonances are s-wave states that autoionize through P_01 coupling. This coupling (0.023) is too weak for the resonances to separate from the discretized s-wave continuum.

2. **Wrong channel character:** The theta-stable features are p-wave/d-wave dominated, not s-wave. They're quasi-bound states of channels 2-4, not the doubly-excited 1S resonances.

3. **Width ratio inverted:** Gamma_A/Gamma_B = 0.24 (experimental: 3.29).

4. **Not converged in n_channels:** 3-channel and 5-channel give qualitatively different spectra.

### Root cause

The coupled-channel ECS on the hyperspherical adiabatic potential fails to resolve the s-wave autoionization resonances because:

- The P_01 coupling is extremely localized (peaked at the avoided crossing) and weak (0.023)
- The adiabatic approximation requires very high angular resolution to separate the s-wave resonance from the continuum
- The channel character migration at the avoided crossing means the "resonance" is spread across channels

### What would be needed

1. **SVD (Slow Variable Discretization):** Avoids the P coupling entirely by working in a non-adiabatic representation. The full 2D problem in (R, alpha) is discretized simultaneously.

2. **Full 2D complex scaling:** Apply ECS to the 6D (or 2D in hyperspherical) coordinate space, not just the hyperradius. This correctly handles the channel-mixing physics.

3. **Much higher angular resolution:** l_max >= 6, n_alpha >= 500. The 2s^2 state requires precise angular correlation.

4. **R-matrix method:** Combine inner-region CI with outer-region channel coupling. More natural for resonances near channel thresholds.

---

## 8. Philosophical Conclusion

**The width IS an eigenvalue** -- of the complex-scaled Hamiltonian. This is mathematically certain (Simon 1979). The question is whether our computational approximation (adiabatic channels + hyperradial ECS) is sufficient to extract it.

The answer: **not with l_max=2 and n_alpha=100.** The angular basis resolves the p-wave/d-wave quasi-bound states but not the subtle s-wave autoionization resonances. The 2s^2/2s3s ratio remains out of reach for the adiabatic+ECS approach at this resolution.

### Updated hierarchy of width-determination methods

| Level | Method | s/p hierarchy | 2s^2/2s3s ratio | Status |
|:-----:|--------|:-------------:|:---------------:|:------:|
| 1 | Gaunt selection rules | YES (100x) | 1.0 | Complete |
| 2 | Weighted graph distance | YES (~7x) | 1.13 | Complete |
| 3 | Avoided crossing (LZ) | fails | fails | Complete (negative) |
| 4 | Single-channel ECS | N/A | N/A | Complete (infra only) |
| 5 | Coupled-channel ECS | YES (upper ch) | 0.24 (inverted) | Complete |
| 6 | Full P integral | -- | -- | Not attempted |
| 7 | SVD / 2D ECS | -- | -- | **Next step** |

The within-sector ratio 3.29 remains the critical test that none of the eigenvalue-based methods have passed. The fiber bundle carries irreducible content.

---

## 9. Files Modified/Created

| File | Changes |
|------|---------|
| `geovac/hyperspherical_complex_scaling.py` | Added `channel_norms` to `solve_ecs_coupled` return; existing threshold fix |
| `tests/test_complex_scaling.py` | Added `TestCoupledECS` (8 tests); 20/20 passing |
| `debug/hyperspherical_exploration/ecs_coupled_spectrum.png` | Theta-stability plot |

---

## 10. Next Steps

1. **SVD method:** Discretize the full 2D (R, alpha) problem on a product grid. Apply ECS to the hyperradius. This avoids the adiabatic approximation entirely and should resolve the s-wave resonances directly.

2. **Higher angular resolution:** If staying with the adiabatic approach, try l_max=6, n_alpha=500, n_channels=10. Computational cost: ~100x current.

3. **Diabatic representation:** Construct diabatic potentials by diagonalizing P at a reference R. The resonance is then a bound state of the diabatic potential, coupled to the continuum by the residual coupling.
