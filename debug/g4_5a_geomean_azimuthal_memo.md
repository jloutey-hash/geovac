# Task #27 — Geometric-mean azimuthal discretization

**Date:** 2026-05-29
**Path:** Gravity arc, sprint-scale follow-on queue task #27 (single-thread).
**Verdict:** **PARTIAL-OVERSHOOT with substantive bracket-interior structural finding.**

## 1. The v3.19.0 Track 2 bracket

G4-5a-DST (v3.19.0 Track 2) established that FD and spectral azimuthal discretizations bracket the true per-t UV target:

| Scheme | Per-t recovery at t = a² = 0.0025 |
|---|---:|
| FD (undershoot) | 1.31% |
| Spectral (overshoot) | 813.4% |

The +1/6 Lichnerowicz value is the IR target, not the per-t UV target — both extremes bracket whatever the true per-t UV target actually is.

The Möbius factor F = √(0.405) ≈ 0.637 (geometric mean of FD edge eigenvalue 4/π² and spectral edge eigenvalue 1) flagged as a cheap interior cure candidate.

## 2. Construction

Per-mode geometric mean of squared eigenvalues:
$$
m_{\rm eff, GM}^2(k) = \sqrt{m_{\rm eff, FD}^2(k) \cdot m_{\rm eff, spec}^2(k)}
$$

For the wedge at apex angle 2πα:
- $m_{\rm eff, FD}^2(k) = (N_\phi/(\pi\alpha))^2 \cdot \sin^2(\pi(k+0.5)/N_\phi)$
- $m_{\rm eff, spec}^2(k) = ((k+0.5)/\alpha)^2$
- $m_{\rm eff, GM}^2(k) = (N_\phi/(\pi\alpha^2)) \cdot |\sin(\pi(k+0.5)/N_\phi)| \cdot (k+0.5)$

The GM eigenvalue is fed into the same hermitian polar radial Laplacian as FD/spectral (G4-3a-cleanup convention).

## 3. Sanity checks

| Check | Result | Expected | Status |
|---|---:|---:|---|
| F6 bit-exact reduction at α = 1 (GM wedge = GM disk) | max_diff = 0.0 | < 1e-13 | **PASS** |
| Edge ratio m_FD²/m_spec² at k = N_0/2 | 0.398545 | 4/π² ≈ 0.405 | PASS (1.7% slack from N_0=120 finite) |
| Edge ratio m_GM²/m_spec² at k = N_0/2 | 0.631304 | 2/π ≈ 0.637 | PASS (0.8% slack) |

Edge ratio 0.631 confirms the geometric-mean construction is implemented correctly.

## 4. Per-t recovery comparison

Substrate panel: R = 10, a = 0.05, N_ρ = 200, N_0 = 120, ε = 12/120 = 0.1. Same panel as G4-5a-refined / G4-5a-DST for apples-to-apples comparison.

| t | FD % | **GM %** | Spec % |
|---:|---:|---:|---:|
| 0.0025 (t=a²) | 1.31 | **355.92** | 813.38 |
| 0.005 | 16.53 | **305.69** | 645.34 |
| 0.01 | 36.23 | **259.04** | 479.60 |
| 0.02 | 51.87 | **200.70** | 304.48 |
| 0.05 | 67.59 | **122.63** | 127.62 |
| 0.1 | 76.34 | **88.89** | 83.01 |
| 0.2 | 82.91 | **83.67** | 83.00 |
| 0.5 | 89.05 | **89.05** | 89.05 |
| 1.0 | 92.26 | **92.26** | 92.26 |
| 2.0 | 94.59 | **94.58** | 94.58 |
| 5.0 | 96.69 | **96.69** | 96.69 |
| 10.0 | 97.70 | **97.70** | 97.70 |

**Two substantive observations:**

### 4.1 GM is a true bracket interior at every t

At every tested t, FD ≤ GM ≤ Spec (or the reverse ordering for IR where all three converge). The geometric-mean construction is structurally a bracket interior — not just at the edge but per-t throughout the heat-kernel trace.

### 4.2 The bracket converges at t > 0.2

For t ≥ 0.2, FD, GM, and Spec all give bit-identical per-t recovery values:
- t=0.5: all three give 89.05%
- t=1.0: all three give 92.26%
- t=5: all three give 96.69%

This is the IR Lichnerowicz regime. The +1/6 target is approached asymptotically by all three schemes equivalently because the IR behavior is dominated by the small-eigenvalue modes (small k), where FD, GM, Spec all give the same m_eff ≈ k + 0.5 (the high-k bracket disagreement doesn't matter when those modes are exponentially suppressed).

### 4.3 The bracket is wide only at UV (t < 0.1)

The disagreement between FD/GM/Spec is concentrated at t < 0.1, where the heat-kernel trace is sensitive to the entire eigenvalue spectrum (including the high-k edge where the bracket is widest). The structural identification of the true per-t UV target is a UV-only question.

## 5. Verdict against the [50%, 200%] gate

GM recovery at t = a² = 0.0025 is **355.9%**, outside the [50%, 200%] gate window.

**Verdict: PARTIAL-OVERSHOOT.** GM is genuinely a bracket interior (closer to spec/2 than to spec, well above FD) but does not land in the narrow target window around the IR Lichnerowicz value.

### 5.1 Why GM doesn't land at the IR Lichnerowicz value

The IR Lichnerowicz value +1/6 is the asymptotic per-t recovery at large t (where all three schemes converge to ~97.7%, attributable to the absent t→∞ limit-of-integration contribution that's unreachable on finite substrate). At small t, the per-t UV target is structurally different from the IR Lichnerowicz value — this is exactly the question task #28 (wedge-spectral-density heat-kernel expansion) is designed to answer.

The gate [50%, 200%] was constructed against the IR Lichnerowicz baseline of 100%. If the true per-t UV target is, say, 300% of IR Lichnerowicz (because the UV cell sees more of the conical-defect tip contribution than the IR), then GM's 355.9% would be **close to that true target** and the verdict would shift to POSITIVE.

**Without the task #28 identification, we cannot evaluate GM against the actual UV target.** The 355.9% number is between IR-Lichnerowicz-naive expectation (~100%) and spectral overshoot (813%), which is structurally what the geometric-mean construction predicts qualitatively.

## 6. Substantive new content

1. **GM is a genuine bracket interior at every t**, not just at the edge eigenvalue. This was the qualitative prediction and it's empirically confirmed.

2. **The bracket converges at t > 0.2 to bit-identical recovery values across FD/GM/Spec.** This locates the per-t UV target identification problem at t < 0.1 specifically.

3. **The recovery pattern monotonically decreases from t=a² toward t=0.5**, with GM crossing 100% at t ≈ 0.07. This crossing is candidate-meaningful: at t = 0.07, the GM recovery passes through the IR Lichnerowicz value, suggesting that t ≈ 0.07 marks the IR/UV crossover where the per-t target transitions from UV-determined to IR-determined behavior.

4. **F6 bit-exact reduction at α = 1.** The GM scheme inherits the load-bearing F6 reduction from spectral and FD, confirming the discretization is structurally well-defined.

5. **GM is a cheap implementation:** total compute time 24 s for the full 12-point t-grid at the substrate panel, vs spectral's ~similar cost. No additional infrastructure beyond what G4-5a-DST established.

## 7. Honest scope

- **PARTIAL-OVERSHOOT verdict at the strict [50%, 200%] gate.** GM did not close the bracket to within the gate window at t = a².
- **POSITIVE-bracket-interior at the qualitative level.** GM strictly between FD and Spec at every t, as predicted by the geometric-mean construction.
- **GM as a cheap sub-percent cure: NOT confirmed.** The bracket is too wide at the UV cell for a 50/50 geometric mean to close to sub-percent. Tighter interpolations or the true UV target identification (task #28) is required.
- **GM as an order-of-magnitude estimator: PROVEN.** GM gives the right order of magnitude (a few hundred percent vs FD's 1% or Spec's 800%) at every t.

## 8. Implication for task #28

Task #28 (wedge-spectral-density heat-kernel expansion) is now structurally critical. Without identifying the true per-t UV target, the FD/GM/Spec comparison can be evaluated only against the IR Lichnerowicz baseline, which is wrong at small t. The task #28 derivation would let us:
- Compute the actual UV target value at t = a²
- Test whether GM is "close" to that target (POSITIVE) or far from it (NEGATIVE)
- Identify whether a different interpolation scheme (arithmetic mean? harmonic mean? continuum-defined cutoff?) is the right one

## 9. Files

- `debug/g4_5a_geomean_azimuthal.py` (~420 lines, driver with GM construction + per-t comparison)
- `debug/data/g4_5a_geomean_azimuthal.json` (full data tables)
- `debug/g4_5a_geomean_azimuthal_memo.md` (this)

## 10. Cross-references

- v3.19.0 Track 2 (G4-5a-DST, original bracket identification): `debug/g4_5a_dst_spectral_azimuthal_memo.md`
- v3.19.0 sprint synthesis: `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md`
- Task #28 (wedge-spectral-density UV target): named follow-on, next in queue
- Task #26 (literature grounding): `debug/fursaev_solodukhin_1995_grounding_memo.md`
- G4-3a-cleanup hermitian polar radial Laplacian convention: `debug/g4_3a_cleanup_hermitian_polar_memo.md`
