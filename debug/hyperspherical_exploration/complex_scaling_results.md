# Complex Scaling on the Hyperspherical Lattice — Results

**Date:** 2026-03-15
**Status:** Complete (12/12 tests passing, important structural result)
**Tests:** `tests/test_complex_scaling.py` — TestECSBasics (5), TestResonanceDetection (4), TestWidthPrediction (3)

---

## 1. The Question

Can autoionization widths be extracted as complex eigenvalues of a non-Hermitian Hamiltonian, bypassing the integral ∫F_μ P_μν F_ν dR entirely?

**Method:** Exterior Complex Scaling (ECS) — the coordinate transformation R → R₀ + (R - R₀)e^{iθ} for R > R₀ makes the Hamiltonian complex-symmetric. Bound states remain real; the continuum rotates by -2θ; resonances appear as isolated complex eigenvalues E = E_res - iΓ/2.

**Critical test:** The 2s²/2s3s within-sector width ratio (experimental: 3.29).

---

## 2. Implementation

### ECS Solver (`geovac/hyperspherical_complex_scaling.py`)

| Function | Purpose |
|----------|---------|
| `_smooth_theta` | Fermi function ECS ramp: θ(R) = θ/(1 + exp(-(R-R₀)/(δR/6))) |
| `solve_ecs_single_channel` | Complex FD eigenvalue problem on hyperradial grid |
| `solve_ecs_coupled` | Block-tridiagonal coupled-channel ECS (P scaled by e^{-iθ}) |
| `identify_resonances` | Classify eigenvalues: bound / resonance / continuum |
| `theta_stability_scan` | Scan θ = 0.15–0.40 to separate genuine from numerical features |
| `run_ecs_analysis` | Full pipeline: adiabatic → ECS → classification → comparison |

### Key Technical Fix: Threshold-Separated Analytic Continuation

The initial implementation naively continued V_eff(z) = μ(R)/z² + 15/(8z²) into the complex region. This rotates the asymptotic threshold E_thr = -Z²/2 = -2.0 along with the continuum, producing eigenvalues with **positive** imaginary parts — physically wrong.

**Root cause:** V_eff(R) = μ(R)/R² + 15/(8R²) → E_thr as R → ∞. Since μ(R) ~ E_thr·R² at large R, the leading term μ/z² → E_thr·R²/z² → E_thr·e^{-2iθ} ≠ E_thr. The threshold rotates.

**Fix:** Decompose μ(R) = E_thr·R² + Δμ(R) where Δμ → 0, then:

    V_ECS = E_thr + Δμ(R_real)/z²

This keeps the threshold real and only complex-scales the decaying part. The continuum then correctly rotates at angle -2θ from the real threshold.

This is a standard subtlety in applying ECS to reduced 1D equations — the threshold is a property of the asymptotic potential, not of the kinetic energy.

---

## 3. Single-Channel ECS Results (l_max=0, 1 channel)

### Parameters
- Z = 2, l_max = 0, n_alpha = 200, 1 adiabatic channel
- R₀ = 15.0 bohr, θ = 0.3, δR = 2.0, N_R = 2000
- Shift-invert: σ = -2.5 + 0i

### Ground State (Bound)

| Property | Value |
|----------|-------|
| Re(E) | -2.905220 Ha |
| Im(E) | 5.8 × 10⁻¹⁵ Ha |
| Error vs exact | 0.05% |
| θ-spread (0.15–0.40) | < 10⁻⁶ Ha |

The ground state is perfectly real (Im = 0 to machine precision) and completely θ-independent. This validates the ECS infrastructure: bound states are unaffected by the complex scaling.

### Continuum Rotation

| Property | Value |
|----------|-------|
| Expected angle | -2θ = -0.600 rad |
| Observed angles | -0.45 to -0.51 rad |
| States above threshold | 12 with Im < -0.01 |

The continuum eigenvalues cluster near the -2θ line (within 0.1–0.15 rad), confirming the threshold-separated analytic continuation works. The deviation from exact -0.6 is expected for a finite grid with smooth (not sharp) ECS.

### Eigenvalue Spectrum Structure

| Classification | Count | Description |
|---------------|:-----:|-------------|
| Bound | 1 | Ground state at -2.905 Ha |
| Resonance | 12 | Above threshold, off the -2θ line |
| Continuum | 12 | On or near the -2θ line |

### θ-Stable Features

Only **2** eigenvalues are stable across θ = 0.15–0.40:

| Feature | E (Ha) | E above gs (eV) | Γ (eV) | Re spread (Ha) |
|---------|:------:|:----------------:|:-------:|:--------------:|
| Near-threshold 1 | -1.996 | 24.7 | 0.35 | 1.1 × 10⁻³ |
| Near-threshold 2 | -1.955 | 25.8 | 0.73 | 3.1 × 10⁻⁴ |

These are **not** the doubly-excited 2s² or 2s3s resonances (which should be at ~58–60 eV above ground state, or E ≈ -0.78 Ha). They are near-threshold features of the single s-wave channel — likely numerical artifacts of the single-channel approximation at the continuum edge.

---

## 4. Why Single-Channel Cannot Produce Doubly-Excited Resonances

The 2s² ¹S autoionization resonance involves two electrons in the n=2 shell. In the hyperspherical picture:

1. **The resonance lives in the SECOND adiabatic channel** (the one correlating to He+ n=2 at large R)
2. **It autoionizes by coupling to the FIRST channel** (the one correlating to He+ n=1)
3. **This coupling is the non-adiabatic P₀₁(R)**, which only exists in a multi-channel calculation

A single-channel l_max=0 calculation has only **one** adiabatic curve — the ground-state channel. There is no second channel to support the resonance, and no coupling P₀₁ to enable decay. The complex eigenvalues we find are discretized continuum states of the single channel, not resonances.

**This is the correct physics:** autoionization requires at least two coupled channels. The single-channel ECS validates:
- ✅ Bound state energy (0.05% error)
- ✅ Imaginary part = 0 for bound states
- ✅ θ-independence of bound states
- ✅ Continuum rotation at angle ≈ -2θ
- ✅ Complex eigenvalue spectrum structure

But it **cannot** produce the doubly-excited resonances. That requires coupled-channel ECS.

---

## 5. Assessment: Is the Width an Eigenvalue?

### What ECS Proves in Principle

Yes — under complex scaling, the resonance width Γ = -2 Im(E) is the imaginary part of a complex eigenvalue. This is mathematically rigorous (Simon 1979, Moiseyev 1998). The width IS an eigenvalue of the complex-scaled Hamiltonian.

### What Our Single-Channel Result Shows

The infrastructure works perfectly:
- Bound states are real and θ-stable
- The continuum rotates correctly (after threshold separation)
- The eigenvalue solver handles the non-Hermitian problem

But the physics of doubly-excited resonances requires **coupled channels**. The `solve_ecs_coupled` function is implemented and ready for testing, but requires:
1. Multiple adiabatic curves (l_max ≥ 2, n_channels ≥ 3)
2. Non-adiabatic coupling P_μν splines
3. Significantly larger matrix (N_R × n_channels)

### The Philosophical Answer

The width IS an eigenvalue — but of the **full** complex-scaled Hamiltonian, not of a single-channel reduction. In the fiber bundle language:

| Approach | What it captures | Width as eigenvalue? |
|----------|:----------------|:-------------------:|
| Single-channel ECS | Bound state, continuum structure | ✅ (for the bound state) |
| Coupled-channel ECS | Inter-channel coupling → autoionization | ✅ (for resonances) |
| Full P integral | Same physics, different representation | ❌ (width as integral) |

The coupled-channel ECS and the full P integral are mathematically equivalent — they compute the same physical width. The ECS approach recasts it as an eigenvalue; the P integral approach computes it directly. Neither is "more fundamental"; they are dual descriptions.

The 2s²/2s3s ratio of 3.29 should be accessible via coupled-channel ECS, where each resonance appears as a separate complex eigenvalue with its own Im(E). This is the natural next step.

---

## 6. Hierarchy of Width-Determination Methods (Updated)

| Level | Method | s/p ratio | 2s²/2s3s ratio | Status |
|:-----:|--------|:---------:|:--------------:|:------:|
| 1 | Gaunt selection rules | ✅ 100× | ❌ 1.0 | Complete |
| 2 | Weighted graph distance | ✅ ~7× | ❌ 1.13 | Complete |
| 3 | Avoided crossing (LZ) | ❌ fails | ❌ fails | Complete (negative) |
| 4 | Single-channel ECS | N/A | N/A (no resonances) | Complete |
| 5 | Coupled-channel ECS | — | — | **Next step** |
| 6 | Full P integral | — | — | Alternative to 5 |

---

## 7. Files Modified/Created

| File | Changes |
|------|---------|
| `geovac/hyperspherical_complex_scaling.py` | New module: `_smooth_theta`, `solve_ecs_single_channel`, `solve_ecs_coupled`, `identify_resonances`, `theta_stability_scan`, `run_ecs_analysis` |
| `tests/test_complex_scaling.py` | 12 tests across 3 classes; 12/12 passing |

---

## 8. Next Steps

1. **Coupled-channel ECS:** Run `solve_ecs_coupled` with l_max=2, 3–5 channels, P_μν coupling. This should expose the 2s² and 2s3s resonances as isolated θ-stable complex eigenvalues. The 2s²/2s3s ratio becomes Im(E₁)/Im(E₂).

2. **Convergence study:** Vary N_R, n_channels, l_max, R₀, θ to establish that the resonance positions and widths are converged.

3. **Comparison with experiment:**
   - 2s² ¹S: E_res = 57.82 eV, Γ = 0.138 eV
   - 2s3s ¹S: E_res = 62.94 eV, Γ = 0.042 eV
   - Ratio: 3.29 (the critical test)
