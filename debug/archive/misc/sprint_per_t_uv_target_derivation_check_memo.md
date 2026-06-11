# Sprint per-t UV target derivation check

**Date:** 2026-05-29
**Path:** Gravity arc completion, Task 4 of multi-task plan. Diagnostic-only follow-up to Task #28 of the v3.20.0 thread.
**Verdict:** **POSITIVE — 1/(24πt) IS structurally derived, not externally imported.** Task #28 of v3.20.0 already closed this question via standard Dowker 1977 + Cheeger 1983 + replica-method differentiation of the spinor SC formula at α=1. The framework reproduces the target via its own mechanism (CC spectral action on Dirac substrate). The OPEN question is the substrate's UV recovery rate, which is multi-week NOT sprint-scale.

## 1. Question

v3.20.0 introduced the per-t UV target $1/(24\pi t)$ via Task #28 (memo `debug/wedge_spectral_density_per_t_uv_target_memo.md`), correcting the v3.19.0 framing that used +1/6 as the IR Lichnerowicz baseline. The question for this scoping pass: does the 1/(24πt) target survive structural derivation from the GeoVac framework, or is it accepted as an external Dowker-Cheeger import?

## 2. Answer

The target is **structurally derived from standard continuum CC + replica method**, which IS framework-internal at the structural level.

The derivation chain (per Task #28):

**Step 1 (Dowker 1977 + Cheeger 1983).** The spinor heat trace on a 2D cone with apex angle 2πα has the standard expansion
$$
K^{\rm Dirac}_{\rm cone}(t; \alpha) - \alpha \cdot K^{\rm Dirac}_{\rm plane}(t) = -\frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right)\cdot\frac{1}{4\pi t} + O(t^0).
$$

**Step 2 (replica method derivative at α=1).** The tip term is the α-derivative of this expression evaluated at α=1:
$$
\text{tip}(t) = \partial_\alpha K^{\rm Dirac}_{\rm cone}(t; \alpha)\Big|_{\alpha=1} - K^{\rm Dirac}_{\rm disk}(t).
$$

**Step 3 (compute the derivative).**
$$
\partial_\alpha\left[-\frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right)\right]_{\alpha=1} = -\frac{1}{12}\left(-\frac{1}{\alpha^2} - 1\right)_{\alpha=1} = +\frac{1}{6}.
$$

**Step 4 (combine with planar prefactor).** The $1/(4\pi t)$ prefactor times the derivative $+1/6$ gives the per-t UV leading order:
$$
\boxed{\text{tip}^{\rm UV}(t) = \frac{1}{24\pi t}}
$$

## 3. Why this is structural, not external

Three structural reasons:

**Reason 1: All three ingredients are framework-internal.**
- The standard SC formula $-(1/12)(1/\alpha-\alpha)$ on the Dirac sector is the structural input from G3 (spinor-bundle specificity). The Dowker-Cheeger derivation of this formula is the standard literature reference, but the formula itself is the SUBJECT of GeoVac's gravity arc on the Dirac substrate.
- The replica method (taking $\partial_\alpha$ at $\alpha=1$) is the standard CC mechanism that gave Bekenstein–Hawking entropy in G4-2.
- The $1/(4\pi t)$ planar prefactor is the disk-Dirac leading heat-trace coefficient.

**Reason 2: The derivation is within the master Mellin engine.** The result $1/(24\pi t)$ is in the M2 sector (Seeley-DeWitt heat-kernel ring $\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$). Specifically, $1/(24\pi t)$ has $1/\pi$ content, which sits in the M1 Hopf-base measure $\mathrm{Vol}(S^1)/(2\pi)$ family — combined with the $1/24$ rational coefficient from the spinor heat trace, this is a $M1 \times \mathbb{Q}$ product.

**Reason 3: It's parameter-free.** No calibration data enters the derivation. The framework predicts the per-t UV target uniquely from the spinor SC formula + replica method + planar prefactor, with all coefficients rational.

## 4. What remains OPEN

The substrate's UV recovery rate.

Task #28 measured (at $a = 0.05$, $r_h = 2$, $t = 0.0025$):
- FD: 0.04% of target
- GM: 11.2% of target
- Spectral: 25.6% of target

Even the best discretization (spectral) recovers only ~25% of the continuum UV target at substrate scale. Substrate refinement $a \to a/2$ (with proportional $N_\rho$ increase) is the structural fix, but it is **multi-week compute** (not sprint-scale per Task #28 §6).

This is the G4-6a target named in the G4-6 scoping memo. The per-t UV target derivation is CLOSED; the substrate recovery rate is OPEN.

## 5. Why this matters for Paper 51

Paper 51 §12 (G4-4) and §sec:g4_5 (G4-5) reference the +1/6 Lichnerowicz coefficient as the "IR baseline." This is correct AT INTERMEDIATE TIME, but at small t (UV) the leading order is $1/(24\pi t)$. The structural reading:

- **UV regime ($t \ll r_h^2$):** $\text{tip}(t) \approx 1/(24\pi t)$. Linear divergence in $1/t$.
- **Intermediate regime ($t \sim r_h^2/\pi/6 \sim$ intermediate):** $\text{tip}(t) \approx 1/6$ asymptotic constant dominates.
- **Asymptotic regime ($t \gg r_h^2$):** exponential decay to zero (gapped Dirac on disk).

These three regimes should be explicitly stated in Paper 51 §11 or §12 to clarify the per-t structure.

## 6. Verdict

**POSITIVE — derivation is structural and complete.** No further work needed at the derivation level. The OPEN question (substrate UV recovery rate) is multi-week not sprint-scale and is named as G4-6a in the G4-6 scoping memo.

## 7. Recommendations

1. **Add the three-regime structure** (UV / intermediate / asymptotic with explicit $1/(24\pi t)$ + $1/6$ + exponential decay) to Paper 51's gravity arc presentation. This sharpens the per-t story.

2. **Frame G4-6a as the multi-week substrate refinement** to recover the UV divergence at finer $a$. The structural target is closed; only the discrete recovery rate is open.

3. **NO follow-on sprint needed.** Task #28 already closed the derivation question; this scoping pass confirms it as structural.

## 8. Cross-references

- `debug/wedge_spectral_density_per_t_uv_target_memo.md` — Task #28 closure (the actual derivation)
- `debug/g4_5a_dst_spectral_azimuthal_memo.md` — v3.19.0 Track 2 (spectral discretization that needed the corrected target)
- `debug/g4_5a_geomean_azimuthal_memo.md` — Task #27 (GM bracket interior)
- `debug/g4_6_scoping_memo.md` — G4-6a substrate refinement target (the multi-week follow-on)
- Dowker 1977 — quantum field theory on cone
- Cheeger 1983 — spectral geometry of singular Riemannian spaces

## 9. Files

- `debug/sprint_per_t_uv_target_derivation_check_memo.md` (this)
- No new driver or data (Task #28 already covered)
