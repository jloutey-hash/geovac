# Sprint G4-5d F14 — 20-point quadrature panel refinement

**Date:** 2026-05-29
**Path:** Gravity arc, sprint-scale follow-on closure of the G4-5d-refined PARTIAL verdict (CLAUDE.md v3.19.0 Track 1).
**Task ID:** #24 in single-thread sprint-scale follow-on queue.
**Verdict:** **POSITIVE-F14-PHI0-PREDICTION-WITHIN-20PCT** at strict 20% gate, in fact at **sub-5%** across all nine channels.

## 1. The PARTIAL artifact (G4-5d-refined, v3.19.0)

G4-5d-refined Track 1 (v3.19.0) verified the sector-wise Mellin moment map empirically with mean deviation 18.2% and max deviation 47.9% under the φ(0) prediction. The verdict landed PARTIAL specifically because the sharp-cutoff channel at Λ ∈ {1, 2} exceeded the strict 20% gate. The poly/Gauss channel (clean quadrature comparison) already matched within 6% at every Λ.

**Diagnosed mechanism (G4-5d-refined memo §6):** the 8-point log-trapezoidal quadrature could not resolve the analytical Θ(1 − tΛ²) edge of the sharp cutoff. At Λ = 2 the cutoff turns off at t = 1/Λ² = 0.25, leaving only t ∈ {0.1, 0.2} of the 8-point panel in the active region — two points cannot capture the integrand magnitude.

## 2. Two cures applied together

**Cure 1 — denser t-grid.** 20-point log-spaced grid from t_UV = a² = 0.0025 to t_IR = 20, vs original 8 points. Distribution: 17 log-spaced base points + 3 anchor points at exact sharp-cutoff edges t = 1/Λ² for Λ ∈ {0.5, 1, 2}. Unique distinct points: 20. The G4-5d-refined panel was a subset (8 of the 20 are within 5% of the original points).

**Cure 2 — explicit edge insertion for sharp cutoff.** For sharp f(x) = Θ(1 − x):

$$
S_{\rm tip}^{\rm sharp}(\Lambda) \;=\; \frac{1}{2}\int_{t_{\rm UV}}^{1/\Lambda^2} \frac{dt}{t}\,\mathrm{tip}(t)
$$

Implementation: build sub-grid `t_sub = t_grid[t_grid < 1/Λ²]`, linearly interpolate tip(t) in log t at the exact edge, append as boundary node, then log-trapezoidal on the sub-grid. The discontinuity is handled analytically, not numerically.

Smooth cutoffs (Gauss, poly) keep the standard log-trapezoidal pass on the full 20-point grid.

## 3. Substrate panel (unchanged from G4-5d-refined)

| Parameter | Value |
|---|---|
| R | 10.0 |
| a | 0.05 |
| N_ρ | 200 |
| N_0 | 120 |
| ε = (α_+ − α_−)/2 | 0.1 |
| t_UV = a² | 0.0025 |
| t_IR | 20 |

Sanity check at t = 1: tip(1) = 0.1537672938 (bit-exact match to G4-5d at diff = 0.0).

## 4. Numerical results

### S_tip empirical (20-point panel + edge insertion)

| Λ | gauss | sharp | poly |
|---|---:|---:|---:|
| 0.50 | +0.362955 | +0.406384 | +0.383720 |
| 1.00 | +0.258255 | +0.297272 | +0.276320 |
| 2.00 | +0.162859 | +0.194684 | +0.175927 |

### Analytical regulated φ(0)[f, Λ]

| Λ | φ(0)_g | φ(0)_sh | φ(0)_po |
|---|---:|---:|---:|
| 0.50 | 6.800020 | 7.377759 | 7.089151 |
| 1.00 | 5.416747 | 5.991465 | 5.702860 |
| 2.00 | 4.037930 | 4.605170 | 4.316612 |

### Predicted vs empirical ratios

| Λ | ratio | empirical | predicted | dev % |
|---|---|---:|---:|---:|
| 0.50 | sharp/gauss | 1.1197 | 1.0850 | **+3.20%** |
| 0.50 | poly/gauss  | 1.0572 | 1.0425 | **+1.41%** |
| 0.50 | sharp/poly  | 1.0591 | 1.0407 | **+1.76%** |
| 1.00 | sharp/gauss | 1.1511 | 1.1061 | **+4.07%** |
| 1.00 | poly/gauss  | 1.0699 | 1.0528 | **+1.63%** |
| 1.00 | sharp/poly  | 1.0758 | 1.0506 | **+2.40%** |
| 2.00 | sharp/gauss | 1.1954 | 1.1405 | **+4.82%** |
| 2.00 | poly/gauss  | 1.0802 | 1.0690 | **+1.05%** |
| 2.00 | sharp/poly  | 1.1066 | 1.0668 | **+3.73%** |

**Max deviation: 4.82%. Mean deviation: 2.67%. Ordering match: True at all Λ.**

## 5. Comparison to G4-5d-refined (v3.19.0 Track 1)

| Statistic | G4-5d-refined (8-pt) | F14 (20-pt) | Improvement |
|---|---:|---:|---:|
| Max deviation | 47.9% | **4.82%** | ×9.94 |
| Mean deviation | 18.2% | **2.67%** | ×6.82 |
| Ordering match at all Λ | YES | YES | — |
| Decision-gate verdict | PARTIAL | **POSITIVE** | gate crossed |

Both cures contributed: the denser grid improved the smooth (Gauss/poly) channels uniformly, while the explicit edge insertion was the substantive cure for the sharp channel at Λ = 1, 2 (where 8 points could only place 2-4 nodes below the edge).

## 6. Structural reading — sector-wise Mellin moment map closed POSITIVE

The G4-5d structural reframing (v3.18.0 → v3.19.0 Track 1) is now empirically confirmed at the strict 20% gate, in fact at sub-5% precision on three independent ratio channels at three independent Λ values:

| Sector | Wilson coefficient | Mellin moment | Empirical confirmation |
|---|---|---|---|
| Bulk Λ⁰ (cosmological constant) | Λ_cc | φ(2) | G7 (v3.9.0) |
| Bulk Λ¹ (Einstein-Hilbert) | G_eff⁻¹ | φ(1) | G7 (v3.9.0) |
| Topological tip (S_BH) | (1/12) S-C | **φ(0)** | **F14 (this sprint), v3.19.0+** |

The φ(0) reading was first proposed in G4-5d-refined as a refinement of the naive G8 φ(2) prediction. F14 closes the empirical confirmation at all three Λ values and all three cutoff variants (Gaussian, sharp, polynomial), within sub-5%.

## 7. Honest scope

- **Quantitatively closed at sub-5%:** the φ(0) prediction holds for the topological tip sector across the three Λ values × three cutoff functions tested. This is **numerical observation**, not theorem-grade. Theorem-grade derivation is deferred to G4-6e (per `debug/g4_6_scoping_memo.md`).
- **Substrate-specific:** results computed on the G4-5a sweet-spot panel (R=10, a=0.05, N_ρ=200, N_0=120). The sector-wise moment map should be substrate-independent in the continuum limit; that's a multi-substrate verification (G4-6a, multi-month sequential gate).
- **Edge interpolation is linear in log t:** higher-order interpolation (cubic spline) would tighten further but is not needed at sub-5%.
- **Three Λ values only:** the original G4-5d panel used these three; extending to Λ ∈ {0.25, 0.5, 1, 2, 4, 8} is a 1-day follow-on if needed.

## 8. What's next in the queue

Task #25 (α > 1 numerical validation at α ∈ {4, 5, 10}) unblocks. Sub-5% closure of F14 is the cleanest possible empirical confirmation of v3.19.0 Track 1; the queue moves to validate the v3.19.0 Track 5 Möbius headline result at new α values not in the original fit set.

## 9. Files

- `debug/g4_5d_F14_20pt_panel.py` (~370 lines, driver with edge-insertion logic)
- `debug/data/g4_5d_F14_20pt_panel.json` (full ratio + φ(0) tables)
- `debug/g4_5d_F14_20pt_panel_memo.md` (this)

## 10. Cross-references

- G4-5d-refined (v3.19.0 Track 1): `debug/g4_5d_refined_phi0_prediction_memo.md`
- G4-5d original (v3.18.0): `debug/g4_5d_cutoff_dependence_memo.md`
- v3.19.0 sprint synthesis: `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md`
- Paper 51 §12.7.7 (the v3.19.0 paper extension): `papers/group5_qed_gauge/paper_51_gravity_arc.tex`
- Sprint-scale follow-on queue: task list (`/tasks list`)
- Sector-wise moment map memory: `memory/sprint_g4_5_sector_mellin_map.md`
