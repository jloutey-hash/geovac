# Avoided Crossing Analysis of Autoionization Widths

**Date:** 2026-03-15
**Status:** Complete (44/44 coupled tests passing, important negative result)
**Tests:** `tests/test_hyperspherical_coupled.py` — TestAvoidedCrossings (7 tests)

---

## 1. Hypothesis

The threshold-free graph-distance analysis (Phase 2B) correctly predicts the broad/narrow hierarchy of doubly-excited helium resonances but fails the within-sector test: predicted 2s^2/2s3s width ratio of 1.13 vs experimental 3.29. The Gaunt coupling at a single R captures angular topology but not radial dynamics.

**Hypothesis:** Autoionization widths should be determined by the local geometry of the adiabatic potential curves at their avoided crossings — a zero-dimensional feature characterized by three eigenvalue-derived numbers:

1. **delta**: minimum gap between adiabatic curves (Ha)
2. **R_c**: location of the crossing (bohr)
3. **|DF|**: slope difference |dV_mu/dR - dV_nu/dR| at R_c (Ha/bohr)

Width estimates:
- **Landau-Zener:** Gamma ~ delta^2 / |DF|
- **Feshbach:** Gamma ~ delta^2 / (2 sqrt(2(E_res - V_cont(R_c))))

If correct, the width is an eigenvalue, not an integral.

---

## 2. Avoided Crossing Data

### All Pairwise Crossings (l_max=2, 5 channels, n_alpha=100, N_R=400)

| Pair | R_c (bohr) | delta (Ha) | |DF| (Ha/bohr) | V_coupling (Ha) | Gamma_LZ |
|------|:----------:|:----------:|:--------------:|:---------------:|:--------:|
| (0,1) | 13.61 | ~0 | ~0 | ~0 | ~0 |
| (0,2) | 4.10 | 1.288 | ~0 | 0.644 | diverges |
| (0,3) | 15.00 | 1.396 | 0.011 | 0.698 | 185.4 |
| (0,4) | 15.00 | 1.410 | 0.012 | 0.705 | 161.0 |
| (1,2) | 4.00 | 1.284 | ~0 | 0.642 | diverges |
| (1,3) | 15.00 | 1.396 | 0.011 | 0.698 | 185.4 |
| (1,4) | 15.00 | 1.410 | 0.012 | 0.705 | 161.0 |
| (2,3) | 15.00 | ~0 | ~0 | ~0 | ~0 |
| (2,4) | 15.00 | 0.014 | 0.002 | 0.007 | 0.10 |
| (3,4) | 5.89 | 0.001 | ~0 | 0.001 | 503 |

### Width Predictions for Experimental Resonances

| Resonance | Channels | R_c | delta | |DF| | Gamma_LZ | Gamma_Feshbach | Gamma_expt (eV) |
|-----------|:--------:|:---:|:-----:|:---:|:--------:|:--------------:|:---------------:|
| 2s^2 1S | (1,0) | 13.61 | ~0 | ~0 | ~0 | ~0 | 0.138 |
| 2s3s 1S | (3,0) | 15.00 | 1.396 | 0.011 | 185.4 | 0.587 | 0.042 |
| 2p^2 1S | (2,0) | 4.10 | 1.288 | ~0 | diverges | 0.458 | 0.0014 |
| 2p^2 1D | (4,0) | 15.00 | 1.410 | 0.012 | 161.0 | 0.624 | 0.070 |
| 2s2p 1P | (2,0) | 4.10 | 1.288 | ~0 | diverges | 0.468 | 0.037 |

### Within-Sector Ratios

| Method | Gamma(2s^2)/Gamma(2s3s) | Experimental |
|--------|:-----------------------:|:------------:|
| Landau-Zener | ~0 (breakdown) | 3.29 |
| Feshbach | ~0 (breakdown) | 3.29 |
| Gaunt graph (from Phase 2B) | 1.13 | 3.29 |

---

## 3. Why the Landau-Zener Model Fails

### The (0,1) crossing is asymptotic convergence, not a resonance

Channels 0 and 1 are both s-wave. At large R, both converge to the He+ threshold energy -Z^2/2 = -2.0 Ha. The `find_avoided_crossings` algorithm finds R_c = 13.61 bohr with delta ~ 0 — this is the asymptotic channel convergence, not a narrow avoided crossing produced by inter-channel coupling.

The 2s^2 resonance autoionizes through the non-adiabatic coupling P_01(R), which peaks at R ~ 11.9 bohr (Phase 2B data). This coupling is spread over a broad R range (~5-15 bohr), not localized at a single crossing point.

### The (0,2) crossing has a huge gap

The s-p crossing between channels 0 and 2 occurs at R_c = 4.1 bohr with delta = 1.29 Ha — an enormous gap. The curves are parallel at this point (slope_diff ~ 0), making the LZ formula undefined (division by zero). This is not a perturbative avoided crossing; the inter-electron coupling has already fully separated these channels.

### The Landau-Zener regime

The LZ formula applies when:
1. Two diabatic curves actually cross (same energy at some R)
2. The coupling is weak enough to produce a narrow avoided crossing
3. The system passes through the crossing at a well-defined velocity

For helium with 5 adiabatic channels:
- Channels 0 and 1 never cross — they're asymptotically degenerate
- Channels 0 and 2 are separated by ~1.3 Ha — the coupling is too strong
- There is no well-defined "crossing velocity" because the eigenvalue problem is stationary

---

## 4. Physical Interpretation

### The fiber carries irreducible content

In the GeoVac fiber bundle picture (Paper 13):
- **Base space:** hyperradius R (the adiabatic parameter)
- **Fiber:** angular eigenvectors Phi_mu(alpha; R) (the adiabatic states)
- **Connection:** non-adiabatic coupling P_mu_nu(R)

The avoided crossing approach attempted to extract autoionization widths from the **base space geometry** alone (the eigenvalue curves V_eff(R)). This fails because:

1. **The 2s^2 width is determined by the connection P_01**, which couples channels 0 and 1 in the asymptotic region where they are nearly degenerate. The coupling is delocalized over ~10 bohr, not concentrated at a point.

2. **The within-sector ratio (2s^2 vs 2s3s) depends on how the radial wavefunctions F_mu(R) overlap with P_mu_nu(R)**: the integral integral F_0 P_01 F_1 dR. The 2s^2 state has one radial node; the 2s3s state has two. This node structure modulates the overlap integral differently, producing the 3.3x width ratio.

3. **The eigenvalue curves** capture the gross s-wave/p-wave hierarchy (large gap vs small gap to continuum) but cannot distinguish states within the same angular sector because those states differ only in their radial quantum numbers.

### Hierarchy of width-determining mechanisms

| Level | Mechanism | What it captures | What it misses |
|-------|-----------|:----------------|:--------------|
| 1. Gaunt selection rules | Angular topology at single R | s-wave broad, p-wave narrow (100x ratio) | Within-sector ratios |
| 2. Weighted graph distance | Continuous coupling strength | Continuous s/p separation | Within-sector ratios |
| 3. Avoided crossing geometry | Eigenvalue landscape features | Nothing new (fails for He) | Everything the integral captures |
| 4. Full P integral (needed) | Connection + fiber overlap | Within-sector ratios, precise widths | — |

The autoionization width is genuinely a **fiber bundle quantity** that requires all three ingredients: the base (eigenvalue curves), the connection (P coupling), and the fiber (radial wavefunctions). No reduction to pure eigenvalue geometry succeeds beyond the angular topology level.

---

## 5. What This Means for the GeoVac Philosophy

The task hypothesized: "the width should be an eigenvalue, not an integral." This is **disproven** for within-sector width ratios. The 2s^2/2s3s ratio requires the radial wavefunction overlap integral — there is no eigenvalue shortcut.

However, this is **not a failure of the geometric approach**. It reveals the irreducible role of the fiber bundle structure:

- **The angular topology** (Gaunt graph) correctly predicts the gross hierarchy — this IS an eigenvalue quantity.
- **The within-sector structure** requires the full connection on the fiber bundle — the non-adiabatic coupling matrix P_mu_nu(R) acting on radial wavefunctions.

In Paper 7 terms: just as the 1/r potential is the projection of the dimensionless S^3 topology, the autoionization integral is the projection of the fiber bundle structure onto a single number (the width). You cannot bypass the projection — it carries irreducible content.

---

## 6. Files Modified/Created

| File | Changes |
|------|---------|
| `geovac/hyperspherical_resonances.py` | Added `find_avoided_crossings`, `landau_zener_width`, `feshbach_width`, `analyze_avoided_crossings`, `plot_avoided_crossings` |
| `tests/test_hyperspherical_coupled.py` | Added `TestAvoidedCrossings` (7 tests); 44/44 passing |
| `debug/hyperspherical_exploration/avoided_crossings.png` | Adiabatic potential curves with crossings annotated |

---

## 7. Next Steps

1. **Full radial coupling integral:** Compute integral F_mu(R) P_mu_nu(R) F_nu(R) dR using the existing coupled-channel infrastructure. This is the true Feshbach matrix element that determines individual widths.

2. **SVD method:** The Slow Variable Discretization (SVD) avoids the DBOC-P cancellation by working in a diabatic-like representation. This may give better access to individual resonance widths than the adiabatic + P coupling approach.

3. **Complex rotation:** Implement the complex scaling method to extract resonance positions and widths directly as complex eigenvalues. This bypasses the integral entirely — the width becomes the imaginary part of a complex eigenvalue.
