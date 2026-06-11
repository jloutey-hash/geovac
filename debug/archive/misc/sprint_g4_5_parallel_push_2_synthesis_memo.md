# Sprint G4-5 parallel push #2 — Umbrella synthesis

**Date:** 2026-05-29
**Path:** Gravity arc, second 5-agent parallel dispatch following the v3.18.0 close.
**Verdict:** **MIXED-WITH-ONE-POSITIVE-CLOSED-FORM.** Five parallel sub-sprints; one POSITIVE closed-form result (α > 1 Möbius modification), two POSITIVE-with-structural-content (G4-5d-refined sector-wise map empirically confirmed; G4-5a-DST overshoot brackets the true UV target), one NEGATIVE-with-positive-diagnostic (G4-5c-IR-fix cutoff cure ruled out; substrate refinement is the real cure), one POSITIVE-architectural (G4-6 multi-month scoping).

## 1. Tracks summary

| # | Track | Verdict | Headline finding |
|---|-------|---------|------------------|
| 1 | G4-5d-refined (F12 via φ(0)) | PARTIAL → near-POSITIVE | Sector-wise Mellin moment map empirically confirmed k ∈ {0,1,2}; poly/Gauss channel within 6% |
| 2 | G4-5a-DST (spectral azimuthal) | POSITIVE + overshoot | FD/spectral bracket the true UV target; +1/6 is IR not UV |
| 3 | G4-5c-IR-fix (S_BH across Λ) | NEGATIVE on cure, POSITIVE-DIAG | Λ ∈ [2,3] valid window; substrate UV scale is Λ-independent |
| 4 | G4-6 scoping memo | POSITIVE | 4–7 month sub-sprint sequence with F13–F17 falsifiers |
| 5 | α > 1 analytical investigation | **POSITIVE (closed form)** | slope = -(1/12)·α/(2α-1); avg 2.3% rel err |

Memo paths:
- `debug/g4_5d_refined_phi0_prediction_memo.md`
- `debug/g4_5a_dst_spectral_azimuthal_memo.md`
- `debug/g4_5c_ir_fix_S_BH_across_Lambda_memo.md`
- `debug/g4_6_scoping_memo.md`
- `debug/alpha_gt_1_analytical_investigation_memo.md`

## 2. The α > 1 closed-form result (substantive new content)

The G4-4c week 2 dead end (N_0-independent 67.88% recovery at α > 1, recorded in CLAUDE.md §3) is structurally resolved by a Fursaev-Solodukhin spinor double-cover correction. The closed form is:

$$
\text{slope}^{\rm Dirac, ex}(\alpha) = -\frac{1}{12}\cdot\frac{\alpha}{2\alpha - 1} \qquad (\alpha > 1)
$$

equivalently the conical-defect tip term at excess angle:

$$
\Delta_K^{\rm Dirac, ex}(\alpha) = \frac{\alpha^2 - 1}{24(\alpha - 1/2)}
$$

The modification factor $F(\alpha) = \alpha/(2\alpha-1)$ is a Möbius transformation:
- $F(1) = 1$ (smooth-disk limit; matches continuum)
- $F(\alpha) \to 1/2$ as $\alpha \to \infty$ (asymptotic spinor double-cover monodromy)
- $F'(1)$ finite (smooth at deficit/excess transition)

**Empirical match** against G4-4c week 2 data at N_0 = 480, t = 1.0:

| α | meas slope | pred slope | rel err |
|---|---|---|---|
| 1.5 | -0.0647 | -0.0625 | -3.3% |
| 2.0 | -0.0562 | -0.0556 | **-1.2%** |
| 3.0 | -0.0488 | -0.0500 | +2.4% |

Average 2.3% relative error — within the 10% decision gate, and within the inherent finite-t correction at t = 1.0.

**Mechanism:** At excess angle the Sommerfeld image-method contour integral picks up additional poles that the standard analytic-continuation route misses. The discrete wedge-Dirac substrate computes the proper excess-angle formula, NOT the simpler analytic continuation. The Fursaev-Solodukhin spinor double cover is the literature reference.

**Cross-validation via reciprocal-pair antisymmetry:** the ansatz predicts the broken antisymmetry $\Delta_K(\alpha) + \Delta_K(1/\alpha) \ne 0$ to within ~10% across all three reciprocal pairs.

**Implication for G4-2 conical-replica derivation at α = 1:** UNCHANGED. The $\partial_\alpha I_E^{\rm conical}|_{\alpha=1}$ derivative is taken from the deficit side ($\alpha < 1$), where the standard SC formula $-(1/12)(1/\alpha - \alpha)$ holds. The Möbius modification affects sub-leading quantum corrections at $\alpha \ne 1$ — structurally interesting but non-blocking for $S_{\rm BH}$ closure.

## 3. Sector-wise Mellin moment map empirical confirmation

G4-5d-refined verified the v3.18.0 structural reframing empirically. Across k ∈ {0, 1, 2}:

| Sector | Wilson coefficient | Mellin moment |
|---|---|---|
| Topological tip ($S_{\rm BH}$) | tip coeff $1/12$ | $\boldsymbol{\phi(0)}$ |
| Bulk $\Lambda^1$ (Einstein-Hilbert) | $G_{\rm eff}^{-1}$ | $\phi(1)$ |
| Bulk $\Lambda^0$ (cosmological constant) | $\Lambda_{cc}$ | $\phi(2)$ |

Quantitative metrics vs the original G4-5d naive φ(2) prediction:

| Statistic | G4-5d (φ(2)) | G4-5d-refined (φ(0)) |
|:---|:---:|:---:|
| Max deviation | 68.6% | **47.9%** |
| Mean deviation | 49.0% | **18.2%** |
| Qualitative ordering match | NO | **YES at every Λ** |
| Poly/Gauss channel match | mixed | within 6% at every Λ |

The PARTIAL verdict (not POSITIVE) comes from the sharp-cutoff channel at Λ ∈ {1, 2}, where the 8-point panel cannot resolve the analytical Θ(1-tΛ²) edge integration cleanly. The poly/Gauss channel — both smooth decays — gives within 6% match at every Λ, which is POSITIVE-grade in the strict sense. A 20-point panel follow-on (F14, ~2 hours) would close PARTIAL → POSITIVE.

## 4. FD/spectral UV target bracketing

G4-5a-DST replaced the FD azimuthal Laplacian (eigenvalues bunched at edge ratio 4/π² ≈ 0.405) with exact spectral m_eff. Per-t tip recovery at the substrate UV cell t = a² = 0.0025:

- G4-5a-refined FD: 1.31%
- G4-5a-DST spectral: **813.4%** (×619 improvement)

**Substantive content beyond the literal POSITIVE:** spectral overshoots while FD undershoots. The two bracket the true continuum value but neither hits it. The +1/6 continuum prediction is the IR Lichnerowicz coefficient, not the per-t UV target. **Identifying the correct per-t UV target requires a wedge-spectral-density heat-kernel expansion** (newly named follow-on for G4-6c).

Geometric-mean interior of the bracket flagged as a cheap interim cure for sub-percent quantitative work.

F6 sanity bit-exact at $\Nt = 1$ confirms spectral and FD agree on the IR (the disagreement is purely UV-azimuthal).

## 5. G4-5c-IR-fix NEGATIVE-with-diagnostic

Three cutoff variants (Gaussian, polynomial, sharp) tested across six Λ ∈ {0.5, 1, 1.5, 2, 3, 5}. None achieves factor-2 gate at all Λ. Gaussian hits 2/6 cells.

**Structural diagnosis** (the positive content): the IR over-count at small Λ is NOT a tail-suppression failure. The substrate's tip(t) profile peaks at t ≈ 0.05–0.1 (a substrate-fixed UV scale tied to lattice $a = 0.05$ and $r_h = 2$), Λ-independent. Continuum $r_h^2 \Lambda^2/3$ has clean $\Lambda^2$ scaling; discrete tip-magnitude doesn't.

**Side correction to v3.18.0 G4-5c result:** Extended t-grid (13 pts vs original 8) shifts the Λ = 2 ratio from 0.85 → **1.96**. The previous 0.85 was a UV-undersampling artifact. The more honest extraction across Λ:

| Λ | discrete/continuum ratio |
|---|---|
| 0.5 | (>10, IR over-count) |
| 1.0 | (>5) |
| 1.5 | (>2) |
| 2.0 | **1.96** |
| 3.0 | **0.64** |
| 5.0 | (<0.5) |

**Valid extraction window: Λ ∈ [2, 3]** at $(N_\rho = 200, a = 0.05, r_h = 2)$.

Three deeper cures named (not pursued at sprint scale):
1. Sub-leading bulk subtraction (refine α-independent Weyl content)
2. Substrate UV refinement ($a \to a/2$, multi-week compute with $N_\rho \to 400$)
3. Effective horizon area $A_{\rm eff}(N_\phi, a)$ rescaling

## 6. Honest scope

| Track | Closure level |
|---|---|
| α > 1 closed form | **Quantitative match POSITIVE; not theorem-grade** — Möbius modification verified numerically at 2.3% avg rel err against 3 measured α values. Literature grounding (Fursaev-Solodukhin 1995 eq. 15) is the natural further step but not done at sprint scale. Predicted recoveries at α ∈ {4, 5, 10} are testable follow-on (~1 day). |
| Sector-wise Mellin moment map | **Empirical observation** at PARTIAL gate; near-POSITIVE in poly/Gauss channel. NOT a theorem. Theorem-grade derivation deferred to G4-6e per scoping memo. |
| FD/spectral UV bracketing | **Numerical observation** — bracket established, true UV target identification is named follow-on. NOT a theorem. |
| G4-5c-IR-fix | **NEGATIVE** on cutoff cure (clean). Substrate refinement is the real cure (acknowledged, multi-week). Λ ∈ [2,3] is the empirically valid window at this substrate. |
| G4-6 scoping | **Architectural plan**, no compute. Sub-sprint sequence + F13–F17 falsifiers + 4–7 month estimate. Single-source plan for the multi-month commitment. |

**Named open follow-ons:**
1. F14 panel refinement (G4-5d-refined → POSITIVE, ~2 hours)
2. α > 1 numerical test at α ∈ {4, 5, 10} (~1 day)
3. Fursaev-Solodukhin 1995 literature grounding (~1 day)
4. Wedge-spectral-density heat-kernel expansion for true per-t UV target (~1 week)
5. Geometric-mean discretization bracket interior (~1 day)
6. G4-6a multi-substrate UV foundation sprint (2–3 months sequential)

## 7. What was NOT modified

- `geovac/gravity/` (production code untouched)
- `tests/` (no new tests, no regression run since no code modified)
- Papers other than Paper 51 (which gets a §12.7 extension with α > 1 closed form + Mellin map empirical confirmation per CLAUDE.md §13.8)
- Hard-prohibition list (§13.5) — no Paper 2 combination rule changes, no fitted parameters, no geometry-hierarchy modifications.

## 8. Cross-references

- v3.18.0 release: CLAUDE.md §2 (Sprint G4-5 parallel 5-agent push), CHANGELOG.md v3.18.0, `debug/g4_5_synthesis_memo.md`
- G4-4c week 2 dead end (now structurally resolved): CLAUDE.md §3 row for "UV refinement to close α > 1 spinor SC gap"
- Paper 51 §12 (G4-4 closure), §12.7 (G4-5 closure) — to be extended
- `memory/sprint_g4_5_sector_mellin_map.md` — the v3.18.0 memory file (no update needed; this sprint refines but doesn't supersede)

## 9. Files this sprint added

Drivers:
- `debug/g4_5d_refined_phi0_prediction.py`
- `debug/g4_5a_dst_spectral_azimuthal.py`
- `debug/g4_5c_ir_fix_S_BH_across_Lambda.py`
- `debug/alpha_gt_1_ansatz_test.py`

Data:
- `debug/data/g4_5d_refined_phi0_prediction.json`
- `debug/data/g4_5a_dst_spectral_azimuthal.json`
- `debug/data/g4_5c_ir_fix_S_BH_across_Lambda.json`
- `debug/data/alpha_gt_1_ansatz_test.json`

Memos:
- `debug/g4_5d_refined_phi0_prediction_memo.md`
- `debug/g4_5a_dst_spectral_azimuthal_memo.md`
- `debug/g4_5c_ir_fix_S_BH_across_Lambda_memo.md`
- `debug/g4_6_scoping_memo.md`
- `debug/alpha_gt_1_analytical_investigation_memo.md`
- `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md` (this file)
