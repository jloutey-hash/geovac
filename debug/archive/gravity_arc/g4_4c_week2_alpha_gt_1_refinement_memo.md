# Sprint G4-4c week 2 — α > 1 branch refinement

**Date:** 2026-05-29
**Path:** G4-4c sub-sprint. Tests whether the α > 1 recovery gap (60–80% of −1/12 at G4-4c first move) closes with UV refinement.
**Verdict:** **PARTIAL — STRUCTURAL FINDING.** UV refinement does not help: recovery is **bit-identical 67.88% across N_0 = 120, 240, 480**. The α > 1 gap is structurally distinct from a finite-substrate artifact. At moderate-large t (t ≈ 10), recovery improves to 73–89% but does not reach 100%.

## 1. Headline result

At fixed substrate R = 10, a = 0.05, sweep N_0 ∈ {120, 240, 480}. Compute slope = Δ_K^Dirac / (1/α − α) at t = 1.0.

| N_0 | α < 1 mean | α > 1 mean | α < 1 recovery | α > 1 recovery |
|---|---|---|---|---|
| 120 | −0.0829 | −0.0566 | **99.49%** | **67.88%** |
| 240 | −0.0829 | −0.0566 | 99.49% | 67.88% |
| 480 | −0.0829 | −0.0566 | 99.48% | 67.88% |

**Bit-identical α > 1 recovery across UV refinement.** The 32% gap is structural, not numerical.

## 2. t-extension test

Increasing t (capturing more topological content):

| α | t = 1 | t = 2 | t = 5 | t = 10 | t = 20 (IR cutoff) |
|---|---|---|---|---|---|
| 3/2 | 0.78 | 0.82 | 0.87 | **0.89** | 0.85 |
| 2 | 0.67 | 0.73 | 0.78 | **0.82** | 0.79 |
| 3 | 0.59 | 0.63 | 0.69 | **0.73** | 0.71 |

Recovery (slope / target −1/12) grows with t up to t ≈ 10, then drops as the IR cutoff R = 10 starts to bite. **At best (t = 10), α = 3/2 reaches 89%, α = 3 reaches 73%.** Still does not reach 100%.

## 3. Structural reading

Two distinct readings:
1. **Sub-leading corrections to the SC continuum formula at large α (excess angle).** The Dowker / Cheeger-Simons spinor SC result $-\frac{1}{12}(1/α - α)$ is derived for α ≤ 1 (deficit angle). For α > 1 (excess angle / saddle cone), sub-leading corrections might be significant.
2. **Different effective tip coefficient at excess vs deficit angles.** The spinor sees the conical defect via half-integer angular momentum and chirality. For α > 1, the holonomy phase around the apex changes character, possibly giving a different effective coefficient.

Both readings predict that the α > 1 branch is **structurally distinct** from the α < 1 branch, not just a finite-substrate artifact. This is consistent with the empirical finding.

## 4. Cross-references
- G4-4c first move: `debug/g4_4c_first_move_wedge_dirac_memo.md`
- G4-3c-proper scalar (also showed α > 1 asymmetry, but at 30% scalar SC recovery): `debug/g4_3c_proper_wedge_memo.md`
- Dowker 1977 / Cheeger-Simons (continuum spinor SC)

## 5. Files
- `debug/g4_4c_week2_alpha_gt_1_refinement.py`
- `debug/data/g4_4c_week2_alpha_gt_1_refinement.json`
- `debug/g4_4c_week2_alpha_gt_1_refinement_memo.md` (this)
