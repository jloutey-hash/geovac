# Diagnostic Arc Report: Asymmetric Bond Convergence in Composed Geometry

**GeoVac v2.0.4 — March 2026**

## 1. Executive Summary

The composed geometry framework (Paper 17) constructs molecular Hamiltonians by coupling
electron-group fiber bundles, each in its own natural coordinate system. For LiH, this
combines a Level 3 hyperspherical core (Li 1s²) with a Level 4 molecule-frame
hyperspherical valence bond (Li 2s + H 1s), connected via Z_eff screening and a
Phillips-Kleinman (PK) pseudopotential.

**Problem:** The equilibrium bond length R_eq diverges linearly with the angular basis
truncation parameter l_max. At l_max=2, R_eq = 3.28 bohr (8.9% error); at l_max=5,
R_eq = 3.97 bohr (31.6% error). The drift rate is +0.23 bohr per unit l_max (R² = 0.996),
with no sign of saturation. This divergence occurs only for asymmetric bonds (Z_A ≠ Z_B);
symmetric systems like H₂ converge monotonically.

**Root cause:** Differential angular correlation. Higher-l channels add angular correlation
energy that preferentially stabilizes diffuse (large-R) geometries. The PK pseudopotential,
which confines the valence electron to the l=0 channel near the core, is spatially fixed
and cannot adapt to the l_max-dependent shift in angular weight. The wavefunction spreading
is 100% nuclear in origin (Task 5: w(0,0) = 0.9998 with e-e interaction only), not from
electron-electron correlation.

**Best solution tested:** R-dependent PK scaling, where the PK strength increases linearly
with R up to a cap: w_PK(R) = δ_{l,0} × min(1.5, R/R_ref). This reduces the drift rate
by 62% (from +0.310 to +0.119 bohr/l_max over l_max=2–4) and achieves 2.0% R_eq error
at l_max=4. Self-consistent PK (50% drift reduction) and eigenchannel rotation (0%)
were also tested.

---

## 2. Diagnostic Arc

### Task 1: Eigenchannel Rotation

**Hypothesis:** The angular coupling matrix V_12(α) has strong off-diagonal structure;
diagonalizing it (eigenchannel basis) might concentrate the wavefunction into fewer channels.

**Result: NEGATIVE.** The eigenvalue spectrum is flat — 17 out of 25 channels are needed for
90% of the Frobenius norm at the 6:1 charge ratio. The coupling matrix is already nearly
diagonal in the Legendre basis. Channel rotation cannot compress the wavefunction.

**Key numbers:** H₂ (1:1): 7/13 channels for 90%. LiH (6:1): 17/25 channels for 90%.

### Task 2: Legendre Potential Convergence

**Hypothesis:** The nuclear Coulomb potential's Legendre expansion converges slowly for
asymmetric charges, causing basis truncation error.

**Result: NEGATIVE (as bottleneck).** Convergence is exponential with decay rates
β ≈ 0.76–1.92, depending on the charge ratio and position. Even at 6:1 asymmetry, the
potential expansion converges to machine precision by k_max ≈ 10. The Legendre representation
of the potential is not the bottleneck.

### Task 3: Angular Wavefunction Decomposition

**Hypothesis:** The wavefunction itself requires many channels, regardless of potential
convergence.

**Result: CONFIRMED.** H₂ (1:1): 77% in (0,0), 3 channels for 90%. LiH (6:1 analog):
26% in (0,0), 8 channels for 90%. The asymmetric nuclear potential drives the electron
density off-center, spreading the wavefunction across many angular channels. This is the
structural bottleneck.

### Task 4: Prolate Spheroidal Compression

**Hypothesis:** Prolate spheroidal harmonics (natural basis for two-center problems) might
compress the angular expansion better than Legendre polynomials.

**Result: NEGATIVE.** At the molecular ρ values relevant to LiH (c ≈ 0.7), S₀ is 82% P₀.
No compression gain. The spheroidal basis only differs significantly from Legendre at
c > 5, which corresponds to unphysically large ρ values.

### Task 5: Nuclear vs Electron-Electron Decomposition

**Hypothesis:** The wavefunction spreading could be from either the asymmetric nuclear
potential or from electron-electron correlation.

**Result: DECISIVE.** With e-e interaction only (diagonal nuclear), w(0,0) = 0.9998 —
the wavefunction is essentially pure s-wave. With nuclear coupling only (no V_ee),
w(0,0) = 0.243 — massive spreading. The spreading is **100% nuclear** in origin.

**Eigenvalue convergence** (Task 5B): At fixed geometry (R, R_e), the lowest eigenvalue
μ₀ converges exponentially with l_max at rate β = 0.42/l_max. The problem is not that
the eigenvalue doesn't converge — it's that the converged eigenvalue at each geometry
depends on l_max differently at different R values.

### Task 6: Projected PK at Fixed Geometry

**Three PK modes compared** at fixed (R=3.015, R_e=1.5):
- **Channel-blind:** PK applied uniformly to all channels. Weakest confinement.
- **l-dependent:** PK applied only to channels with l₁=0 or l₂=0 (δ_{l,0} weighting).
- **Projected:** PK weighted by actual l=0 content of the eigenstate (iterative).

All three converge at fixed geometry. Projected PK weights ≈ 0.45 for 6:1 asymmetry
(vs 1.0 for l-dependent). The issue is not fixed-geometry convergence but R-dependent
behavior.

### Task 7: PES Sweep with Three PK Modes

**Result:** l-dependent PK gives the best R_eq at l_max=2–3. Projected PK collapses
(no equilibrium found — PES becomes monotonically attractive). Channel-blind diverges
fastest. This establishes that l-dependent PK is the baseline to improve upon.

**Key numbers (l_max=3):** Channel-blind R_eq ≈ 3.60 (19.3% error), l-dependent
R_eq ≈ 3.47 (15.1% error), projected PK: no minimum.

### Task 8: PK Saturation in Simplified Model

**Hypothesis:** R_eq might saturate (stop growing) at high l_max if the angular channel
weight redistribution eventually stabilizes.

**Result: ARTIFACT.** In the simplified model (constant Z_eff, V_NN = Z_eff/R), R_eq
appeared to saturate at ~4.45 bohr. However, this model used the wrong V_NN (Z_eff=2.69
instead of Z_bare=3.0) and omitted V_cross_nuc (core-electron attraction to nucleus B).
The apparent saturation was a consequence of these model errors.

**Diagnostic R_eq values:** l_max 2: 3.71, 3: 4.13, 4: 4.33, 5: 4.45, 6: 4.50, 7: 4.54.

### Task 9: PK Saturation in Real Solver

**Result: NO SATURATION.** Using the full ComposedDiatomicSolver with correct V_NN = 3.0/R
and V_cross_nuc, R_eq diverges linearly:

| l_max | R_eq (bohr) | ΔR_eq | % error |
|:-----:|:-----------:|:-----:|:-------:|
| 2 | 3.282 | — | 8.9% |
| 3 | 3.471 | +0.189 | 15.1% |
| 4 | 3.729 | +0.258 | 23.7% |
| 5 | 3.970 | +0.240 | 31.6% |

Linear fit: R_eq = 2.801 + 0.232 × l_max (R² = 0.996). No hint of saturation.
The simplified model's "saturation" was entirely from using Z_eff instead of Z_bare in V_NN.

### Task 10: Self-Consistent and R-Dependent PK

**Self-consistent PK:** Iteratively updates PK weights based on the l=0 content (p₁, p₂)
of the eigenstate. Converges in ~10 iterations with 0.5 damping. Reduces drift by ~50%
but still diverges.

**R-dependent PK:** Strengthens the PK barrier at large R where the valence cloud is more
diffuse: w_PK(R) = δ_{l,0} × min(cap, R/R_ref) with R_ref = 3.015 bohr, cap = 1.5.

| Mode | l_max | R_eq | % error |
|:-----|:-----:|:----:|:-------:|
| l-dependent | 2 | 2.803 | 7.0% |
| l-dependent | 3 | 3.015 | 0.0% |
| l-dependent | 4 | 3.423 | 13.5% |
| self-consistent | 2 | 2.878 | 4.5% |
| self-consistent | 3 | 3.026 | 0.4% |
| self-consistent | 4 | 3.192 | 5.9% |
| R-dependent (cap=1.5) | 2 | 2.837 | 5.9% |
| R-dependent (cap=1.5) | 3 | 2.975 | 1.3% |
| R-dependent (cap=1.5) | 4 | 3.075 | 2.0% |

**Drift rates (bohr/l_max):**
- l-dependent: +0.310
- self-consistent: +0.157 (50% reduction)
- R-dependent (cap=1.5): +0.119 (62% reduction)

**Key observation:** The self-consistent PK weights (p₁, p₂) drop strongly with R:
at l_max=4, p₁ = 0.85 at R=2.0 but p₁ = 0.33 at R=6.5. This R-dependence of the
effective PK strength is precisely what the R-dependent scaling captures.

### Task 11: Higher l_max Test (R-Dependent PK at l_max = 5, 6, 7)

Extended both l-dependent and R-dependent PK to l_max = 5, 6, 7.

**Key finding: R-dependent PK advantage erodes at higher l_max.** The R-dependent PK
drift rate was 62% lower than l-dependent at l_max 2–4 (Task 10), but this advantage
narrows dramatically at higher l_max:

| Mode | l_max | R_eq (bohr) | ΔR_eq | % error |
|:-----|:-----:|:-----------:|:-----:|:-------:|
| l-dependent | 2 | 2.803 | — | 7.0% |
| l-dependent | 3 | 3.015 | +0.212 | 0.0% |
| l-dependent | 4 | 3.423 | +0.408 | 13.5% |
| l-dependent | 5 | 3.474 | +0.051 | 15.2% |
| l-dependent | 6 | 3.562 | +0.088 | 18.1% |
| l-dependent | 7 | 3.638 | +0.076 | 20.7% |
| R-dep (cap=1.5) | 2 | 2.837 | — | 5.9% |
| R-dep (cap=1.5) | 3 | 2.975 | +0.138 | 1.3% |
| R-dep (cap=1.5) | 4 | 3.075 | +0.100 | 2.0% |
| R-dep (cap=1.5) | 5 | 3.346 | +0.271 | 11.0% |
| R-dep (cap=1.5) | 6 | 3.462 | +0.116 | 14.8% |
| R-dep (cap=1.5) | 7 | 3.610 | +0.148 | 19.7% |

**Overall drift rates (bohr/l_max):**
- l-dependent: +0.168
- R-dependent (cap=1.5): +0.160
- self-consistent: +0.157

**Critical observation:** The drift rates are nearly identical when computed over the
full l_max=2–7 range. The R-dependent PK's advantage was concentrated at l_max 2–4;
at l_max 5–7, both modes drift at similar rates (+0.05 to +0.27 bohr/step). The
l-dependent mode shows non-monotonic drift increments (large jump at l_max=4, then
smaller), while R-dependent shows a large jump at l_max=5 (+0.271).

**Interpretation:** The R-dependent PK delays the onset of divergence but does not
fundamentally change the asymptotic drift rate. The PK scaling compensates for the
l=0 weight reduction at moderate l_max but cannot counteract the deep structural
mechanism (differential angular correlation from the nuclear potential) at high l_max.

**Computational timing:** l_max=7 (32 channels) takes ~400s per PES sweep, scaling
roughly as n_ch² as expected for dense eigenvalue problems.

---

## 3. Mechanism: Differential Angular Correlation

The l_max divergence in asymmetric bonds has a clear mechanistic explanation, supported
by the diagnostic data:

### 3.1 Why Higher l Channels Shift R_eq Outward

The molecule-frame hyperspherical wavefunction is expanded in angular channels (l₁, l₂):

ψ(α, R_e) = Σ_{l₁,l₂} c_{l₁l₂}(α) f_{l₁l₂}(R_e)

For a symmetric bond (H₂), the nuclear potential is parity-symmetric, and the
wavefunction concentrates in the (0,0) channel: 77% of the weight (Task 3). Adding
higher-l channels provides angular correlation energy that is roughly R-independent,
so R_eq is stable.

For an asymmetric bond (LiH, Z_A/Z_B = 6:1 effective), the nuclear potential breaks
parity. The electron density is pulled toward the heavier nucleus, spreading across
many angular channels: only 26% in (0,0), 8 channels needed for 90% (Task 3).

The crucial asymmetry is **R-dependent**: at large R, the valence electron cloud is
more diffuse, providing more spatial structure for higher-l harmonics to capture. Each
unit increase in l_max adds roughly equal extra stabilization at large R relative to
small R (evidenced by the linear drift: +0.23 bohr/l_max, Task 9).

### 3.2 The PK Barrier Cannot Adapt

The Phillips-Kleinman pseudopotential is a fixed spatial function:

V_PK(r) = C_core × exp(-β_core × r²)

It is localized near the Li core and acts as a repulsive barrier preventing the valence
electron from collapsing into core orbitals. In the l-dependent scheme, PK is applied
only to l=0 channels (where the core orbital lives), with unit weight.

The problem: as l_max increases, the l=0 channel becomes a smaller fraction of the
total wavefunction. The PK barrier's effective strength relative to the total Hamiltonian
decreases. At large R (where the valence cloud is diffuse), this decrease is more
pronounced because:
1. More angular structure is available for higher-l channels to capture
2. The PK barrier is localized near the core, far from the valence region

This creates a systematic bias: each l_max increment deepens the PES more at large R
than at small R, shifting the minimum outward.

### 3.3 Evidence Summary

| Evidence | Source | Implication |
|:---------|:-------|:------------|
| w(0,0) = 0.9998 with e-e only | Task 5 | Spreading is 100% nuclear |
| Eigenvalue converges exponentially at fixed geometry | Task 5B | Problem is differential, not absolute |
| R_eq drift is linear (+0.23 bohr/l_max) | Task 9 | No saturation in real solver |
| p₁ drops from 0.85 → 0.33 as R: 2.0 → 6.5 (l_max=4) | Task 10 | PK weight is strongly R-dependent |
| R-dependent PK reduces drift 62% | Task 10 | Correct mechanism identified |
| Simplified model saturates but real solver doesn't | Tasks 8–9 | V_NN and V_cross_nuc are load-bearing |

---

## 4. Solutions Tested

| Approach | Result | Drift Reduction | Best R_eq Error | Reference |
|:---------|:-------|:---------------:|:---------------:|:---------:|
| Channel-blind PK | Diverges (fastest) | baseline | 10.2% (l=2) | Task 7 |
| l-dependent PK (δ_{l,0}) | Diverges | — | 0.0% (l=3, lucky) | Task 9 |
| Eigenchannel rotation | No effect | 0% | — | Task 1 |
| Spheroidal basis | No effect | 0% | — | Task 4 |
| Projected PK (iterative) | PES collapses | — | no minimum | Task 7 |
| Self-consistent PK | Reduces drift | 50% | 0.4% (l=3) | Task 10 |
| R-dependent PK (cap=1.5) | Delays onset | 62% (l≤4), ~5% (l≤7) | 2.0% (l=4), 19.7% (l=7) | Tasks 10–11 |
| Higher l_max alone | Diverges | — | — | Task 9 |

**Why some approaches fail:**

- **Eigenchannel rotation** (Task 1): The coupling matrix eigenvalues are already
  well-spread; rotating to the eigenbasis doesn't concentrate the wavefunction because
  the spreading is driven by the diagonal nuclear potential, not off-diagonal coupling.

- **Spheroidal basis** (Task 4): At the relevant ρ values, spheroidal harmonics are
  ≈ Legendre polynomials. No compression is possible because the molecular geometry
  parameter c is too small.

- **Projected PK** (Task 7): Weighting PK by the actual l=0 content makes it
  self-defeatingly weak at large R, where the l=0 fraction is already small. This
  creates positive feedback: weak PK → more spreading → even weaker PK → no minimum.

---

## 5. Recommended Implementation

**Note (Task 11 update):** R-dependent PK is the best single-parameter fix tested,
but Task 11 shows it delays divergence rather than eliminating it. Use l_max ≤ 4 for
production results with R-dependent PK. For converged results, a fundamentally
different approach (full 4-electron Level 4, or l_max-dependent R_ref) is needed.

### 5.1 R-Dependent PK Scaling

The recommended PK modification:

```
w_PK(l₁, l₂, R) = (δ_{l₁,0} × f(R), δ_{l₂,0} × f(R))

where f(R) = min(cap, max(1.0, R / R_ref))
```

**Parameters:**
- R_ref = 3.015 bohr (experimental R_eq for LiH)
- cap = 1.5 (insensitive to exact value; cap=2.0 gives identical results)

**Code change:** One modification in `level4_multichannel.py` where PK weights are
applied to the angular Hamiltonian. The scaling factor `f(R)` multiplies the existing
l-dependent weight.

### 5.2 Toward Ab Initio R_ref

The current implementation uses R_ref = R_eq(experimental), which violates the
framework's zero-parameter philosophy. Two paths to derive R_ref from atomic properties:

1. **Core screening radius:** The PK pseudopotential has a natural length scale
   r_PK = 1/√β_core ≈ 0.38 bohr. The reference distance where PK strengthening
   should begin could be set as R_ref = 2 × r_PK × (Z_A/Z_eff)^{1/3}, giving a
   purely atomic parameter.

2. **Valence cloud extent:** The hydrogenic valence orbital of Li has a characteristic
   radius ⟨r⟩ = 3n²/(2Z_eff) = 6.0 bohr for n=2, Z_eff=1. Setting R_ref as the
   distance where the valence density peaks on the internuclear axis gives a
   geometry-independent scale.

Neither has been tested. The key constraint is that R_ref must be O(R_eq) for the
scaling to have the correct magnitude.

---

## 6. Open Questions

1. **R-dependent PK drift continues at higher l_max.** Task 11 confirmed: at l_max=7,
   R_eq = 3.61 bohr (19.7% error). The overall drift rate (+0.160 bohr/l_max) is
   nearly identical to l-dependent (+0.168). The R-dependent PK delays divergence
   onset but does not change the asymptotic rate. A fundamentally different approach
   is needed for l_max > 4.

2. **Can R_ref be derived ab initio?** The R_ref parameter is currently the only
   non-ab-initio input. Deriving it from the core screening radius or valence extent
   would restore the zero-parameter property.

3. **Would full 4-electron Level 4 eliminate PK entirely?** The composed geometry
   separates core and valence electrons, introducing the PK as an approximation. A
   full 4-electron Level 4 treatment (all electrons in one hyperspherical expansion)
   would avoid PK entirely but at O(l_max^4) channel scaling.

4. **How does R-dependent PK perform for BeH₂ and H₂O?** These polyatomic systems
   use the same composed geometry with inter-fiber coupling. The PK scaling should
   apply independently to each bond pair.

5. **Is the residual drift a PK artifact or a genuine composed-geometry limitation?**
   The +0.119 bohr/l_max rate with R-dependent PK is 62% less than standard l-dependent
   but still nonzero. This could indicate either (a) the R-dependent scaling function
   needs further refinement, or (b) the composed geometry's adiabatic separation
   introduces an inherent l_max-dependent error.

6. **What is the role of V_cross_nuc?** Task 9 showed that the simplified model
   (without V_cross_nuc) saturates while the real solver doesn't. This core-electron
   attraction to nucleus B may amplify the differential correlation at large R.

---

## 7. Files Produced

### Scripts

| File | Task | Description |
|:-----|:----:|:------------|
| `debug/eigenchannel_diagnostic.py` | 1 | Eigenchannel rotation analysis |
| `debug/legendre_convergence.py` | 2 | Legendre potential convergence |
| `debug/angular_wavefunction_diagnostic.py` | 3 | Angular wavefunction decomposition |
| `debug/spheroidal_compression.py` | 4 | Prolate spheroidal basis test |
| `debug/nuclear_vs_ee_diagnostic.py` | 5 | Nuclear vs e-e decomposition |
| `debug/projected_pk.py` | 6 | Projected PK at fixed geometry |
| `debug/pes_sweep.py` | 7 | PES sweep with three PK modes |
| `debug/pk_saturation_test.py` | 8 | Simplified model PK saturation |
| `debug/pk_saturation_real.py` | 9 | Real solver PK saturation |
| `debug/pk_saturation_real_run.py` | 9 | Runner script for Task 9 |
| `debug/selfconsistent_pk.py` | 10 | Self-consistent + R-dependent PK |
| `debug/selfconsistent_pk_postprocess.py` | 10 | Post-processing for Task 10 |
| `debug/task11_higher_lmax.py` | 11 | Higher l_max PES sweeps |

### Data Files

| File | Task | Description |
|:-----|:----:|:------------|
| `debug/data/eigenchannel_diagnostic.json` | 1 | Eigenvalue spectra, cumulative fractions |
| `debug/data/legendre_convergence.json` | 2 | Convergence rates by charge ratio |
| `debug/data/angular_wavefunction_decomposition.json` | 3 | Channel weights at each l_max |
| `debug/data/spheroidal_compression.json` | 4 | Legendre vs spheroidal comparison |
| `debug/data/nuclear_vs_ee_decomposition.json` | 5 | Nuclear-only, e-e-only, full weights |
| `debug/data/projected_pk_diagnostic.json` | 6 | Three PK modes at fixed geometry |
| `debug/data/pes_sweep_pk_modes.json` | 7 | PES curves for three PK modes |
| `debug/data/pk_saturation_test.json` | 8 | Simplified model R_eq vs l_max |
| `debug/data/pk_saturation_real.json` | 9 | Real solver R_eq vs l_max |
| `debug/data/selfconsistent_pk.json` | 10 | Self-consistent + R-dependent results |
| `debug/data/task11_higher_lmax.json` | 11 | Higher l_max results |

### Plots

| File | Task | Description |
|:-----|:----:|:------------|
| `debug/plots/eigenchannel_convergence.png` | 1 | Eigenvalue spectrum comparison |
| `debug/plots/eigenchannel_coupling_matrices.png` | 1 | Coupling matrix visualization |
| `debug/plots/eigenchannel_spectra_rho1.png` | 1 | Spectra at ρ=1.0 |
| `debug/plots/eigenchannel_decay_rate.png` | 1 | Decay rate comparison |
| `debug/plots/eigenchannel_key_diagnostic.png` | 1 | Key diagnostic summary |
| `debug/plots/legendre_coefficients.png` | 2 | Legendre expansion coefficients |
| `debug/plots/legendre_convergence_loglinear.png` | 2 | Log-linear convergence |
| `debug/plots/legendre_convergence_loglog.png` | 2 | Log-log convergence |
| `debug/plots/legendre_potential_profiles.png` | 2 | Potential profiles |
| `debug/plots/angular_Veff_00.png` | 3 | Effective potential in (0,0) channel |
| `debug/plots/angular_channel_weights.png` | 3 | Channel weight distribution |
| `debug/plots/angular_cumulative_weight.png` | 3 | Cumulative weight convergence |
| `debug/plots/angular_eigenvalue_convergence.png` | 3 | μ₀ convergence with l_max |
| `debug/plots/lmax_convergence_extrapolation.png` | 3 | l_max extrapolation |
| `debug/plots/lmax_weight_evolution.png` | 3 | Weight evolution with l_max |
| `debug/plots/spheroidal_compression_cumulative.png` | 4 | Cumulative weight comparison |
| `debug/plots/spheroidal_compression_weights.png` | 4 | Weight distribution |
| `debug/plots/spheroidal_harmonics_6to1.png` | 4 | Spheroidal harmonics at 6:1 |
| `debug/plots/nuc_vs_ee_channel_weights.png` | 5 | Nuclear vs e-e weight comparison |
| `debug/plots/projected_pk_channel_weights.png` | 6 | PK mode weight comparison |
| `debug/plots/projected_pk_mu0_vs_lmax.png` | 6 | μ₀ convergence by PK mode |
| `debug/plots/pes_sweep_pk_modes.png` | 7 | PES curves by PK mode |
| `debug/plots/pes_sweep_req_vs_lmax.png` | 7 | R_eq vs l_max by PK mode |
| `debug/plots/pk_saturation_emin_vs_lmax.png` | 8 | E_min convergence (simplified) |
| `debug/plots/pk_saturation_pes_curves.png` | 8 | PES curves (simplified) |
| `debug/plots/pk_saturation_req_vs_lmax.png` | 8 | R_eq vs l_max (simplified) |
| `debug/plots/pk_saturation_differential_corr.png` | 8 | Differential correlation |
| `debug/plots/pk_saturation_drift_increments.png` | 8 | Drift increments |
| `debug/plots/richardson_extrapolation.png` | 8 | Richardson extrapolation |
| `debug/plots/pk_saturation_real_vs_diag.png` | 9 | Real vs diagnostic comparison |
| `debug/plots/pk_saturation_real_pes.png` | 9 | PES curves (real solver) |
| `debug/plots/pk_saturation_real_drifts.png` | 9 | Drift increments (real) |
| `debug/plots/selfconsistent_pk_convergence.png` | 10 | Iteration convergence |
| `debug/plots/selfconsistent_pk_weights_vs_R.png` | 10 | PK weights vs R |
| `debug/plots/selfconsistent_pk_req_vs_lmax.png` | 10 | R_eq comparison (3 modes) |
| `debug/plots/selfconsistent_pk_drift_rates.png` | 10 | Drift rate comparison |
| `debug/plots/task11_req_vs_lmax.png` | 11 | R_eq vs l_max (all modes, extended) |
| `debug/plots/task11_drift_increments.png` | 11 | Drift increments (extended) |
| `debug/plots/task11_pes_rdependent.png` | 11 | PES curves (R-dependent, extended) |
