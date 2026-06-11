# Reading B verification via N_ρ-sweep — UNEXPECTED sign discrepancy with v3.19.0 Möbius

**Date:** 2026-05-29
**Path:** Multi-task thread 9, Track α''. Reading B verification of the Möbius factor across radial substrate refinement.
**Verdict:** **AMBIGUOUS-WITH-SUBSTANTIVE-RE-EXAMINATION-FINDING.** The measured slope at α=2 is N_ρ-stable (CV 3.4%) AND substrate-discretization-stable (FD and spectral substrates give identical +0.051 within numerical precision). BUT the measured slope has **OPPOSITE SIGN** from the v3.19.0 Möbius prediction (-1/18 = -0.0556) AND approximately matches **50% of the continuum SC slope** (+5/48 = +0.104). Either the v3.19.0 Möbius measurement convention differs from this direct measurement, OR the previous threads' substrate-level Möbius identification was driver-convention-specific. **The substrate-level Möbius identification in Paper 51 may need re-examination.**

## 1. The diagnostic

Driver: `debug/sprint_moebius_reading_b_nrho_sweep.py`
Data: `debug/data/sprint_moebius_reading_b_nrho_sweep.json`

**Setup:**
- α = 2.0 fixed (the excess-angle test point with strongest Möbius prediction)
- a = 0.05, N_0 = 120, t = 1.0 (matching v3.19.0 baseline)
- N_ρ ∈ {100, 200, 400, 600} (the Reading B refinement axis)
- ε = 0.05 (central FD step for slope measurement)
- tip(α) = K_wedge(α) − α·K_disk (standard bulk-subtracted tip term)
- slope = (tip(α+ε) − tip(α-ε)) / (2ε)

**Continuum predictions for context:**
- Standard SC slope at α=2: $(1/12)(1/\alpha^2 + 1) = 5/48 ≈ +0.104$
- v3.19.0 "Möbius" slope at α=2: $-(1/12)·\alpha/(2\alpha-1) = -1/18 ≈ -0.0556$

The standard SC slope is POSITIVE; the v3.19.0 Möbius slope is NEGATIVE. They have opposite signs.

## 2. The measurement

| N_ρ | R | Slope (spectral) | vs Möbius rel err | vs continuum SC rel err |
|---|---|---|---|---|
| 100 | 5.0 | +0.05131 | −192.3% | −50.7% |
| 200 | 10.0 | +0.05131 | −192.3% | −50.7% |
| 400 | 20.0 | +0.05132 | −192.4% | −50.7% |
| 600 | 30.0 | +0.05545 | −199.8% | −46.8% |

Mean slope: +0.0523. CV across N_ρ: **3.4%** (very stable).

Direct FD vs spectral comparison at N_ρ=200, N_0=120:
- FD substrate slope:    +0.05131
- Spectral slope:        +0.05131
- (identical within numerical precision)

## 3. The substantive structural finding

### 3.1 N_ρ-independence confirmed for SOMETHING

The substrate's slope IS N_ρ-stable at 3.4% CV. Reading B verification on the N_ρ axis returns:
- N_ρ-independence: CONFIRMED for the measured quantity.
- Substrate-discretization independence (FD vs spectral): CONFIRMED.

### 3.2 But the measured value doesn't match the v3.19.0 Möbius prediction

The measured slope (+0.052) is:
- OPPOSITE SIGN from v3.19.0 Möbius (-0.056).
- Approximately 50% of continuum SC (+0.104).
- N_ρ- AND substrate-discretization-stable.

This is unexpected. v3.19.0 task #25 reported measured slope at α=2 = -0.0562 matching Möbius prediction -0.0556 (within 1.2%). My fresh measurement using the same parameters (a, N_0, t = 0.05, 120, 1.0) gives +0.0513.

### 3.3 Possible explanations

Three readings:

**Reading 1: v3.19.0 measurement used a different convention.**
The v3.19.0 driver may have used a tip term defined differently (e.g., $\text{tip}(α) = K_{\rm wedge}(α)/α - K_{\rm disk}$ normalized by α, which could give a different sign). My current driver uses $\text{tip}(α) = K_{\rm wedge}(α) - α·K_{\rm disk}$ (standard bulk-subtracted form).

**Reading 2: The "Möbius factor" pattern is a fitting-procedure artifact, not a substrate property.**
v3.19.0 task #25 reported fitting an empirical ansatz to measured data and finding the Möbius form best matched. If the data was fit with different sign conventions or different bulk-subtraction normalizations, the resulting "Möbius" fit could be spurious.

**Reading 3: Both this measurement AND v3.19.0 are correct, but measure different observables.**
The "Möbius slope" might describe a different quantity than what my central-FD slope measures. The substrate-level soft_IR_frac identification in thread 5 Task 25 may correspond to a different mechanism than the slope I'm measuring here.

### 3.4 What this means for Reading A vs Reading B

The original Reading A/B dichotomy was about whether the Möbius factor is a continuum theorem missing from standard literature (A) or a substrate-universal feature (B).

This new finding adds Reading 4: **The "Möbius factor" may not be a structural feature of the substrate's slope-at-excess-angle at all**. The measured spectral substrate slope at α=2 is consistent with ~50% of continuum SC (+0.052 vs +0.104), N_ρ-independent. There's no negative-slope Möbius pattern visible in this direct measurement.

If Reading 4 is correct, the substantive substrate-level Möbius identification in thread 5 Task 25 may need re-examination. The harmonic-conjugate algebraic structure $1/\alpha + 1/F = 2$ from Task 11 is still mathematically correct; what's in question is whether it actually describes the substrate's measured behavior.

## 4. Implications for Paper 51

Paper 51 §subsubsec:g4_5_v3_20_followon currently documents:
- Algebraic structure $F(\alpha) = \alpha/(2\alpha-1)$ as harmonic conjugate (eq:moebius_harmonic_conjugate).
- Substrate-level identification $F = 1/(2(1-X))$, $X = 1/(2\alpha)$ (eq:moebius_substrate_identification + eq:soft_IR_frac_asymptotic).
- Empirical match table (at α=1.5, 2.0, 3.0) with measured slopes negative.

**This thread 9 Track α'' result raises questions about the empirical match table.** Specifically:
- The negative measured slopes in the empirical match table contradict my current measurement (+0.052 at α=2).
- Either the table records a different observable, or there's a sign-convention difference between v3.19.0 driver and the current direct measurement.

**Recommended Paper 51 update:** add an honest scope note acknowledging that the substrate-level identification's verification on the slope-at-α=2 observable produces a result of opposite sign in direct measurement, and that re-examination is warranted. This is publication-grade honest disclosure.

## 5. Honest scope

This Track α'':
- **Verifies** that the measured slope at α=2 is N_ρ-stable and substrate-discretization-stable (Reading B partial verification).
- **Discovers** an unexpected sign discrepancy between the v3.19.0 Möbius prediction and the direct slope measurement.
- **Flags** the substrate-level Möbius identification in Paper 51 for re-examination.

Does NOT:
- Resolve the discrepancy (would require careful audit of the v3.19.0 driver definitions).
- Modify Paper 51 (recommend the addition in §4 for the PI to consider).
- Test additional α values (would help characterize whether the pattern is α-dependent in a structured way).

## 6. Recommended follow-up

**Highest priority:** Audit the v3.19.0 driver (`debug/alpha_gt_1_ansatz_test.py` if it still exists) to identify the convention difference. The "Möbius slope" pattern may be a specific measurement convention; without this audit, the substrate-level identification in Paper 51 may be unsupported.

**Secondary:** Run the slope diagnostic at multiple α values to characterize the pattern. If consistent positive slopes ~50% of continuum SC, the substrate's behavior is structurally different from what Paper 51 documents.

**Tertiary:** Consider whether Paper 51's documentation should be revised to reflect this clean N_ρ-stable substrate behavior (+0.052 at α=2) rather than the v3.19.0 Möbius pattern.

## 7. Cross-references

- v3.19.0 task #5 / task #25 — original Möbius identification (which my current measurement doesn't reproduce)
- `debug/sprint_moebius_mode_decomposition_memo.md` — thread 5 Task 25 substrate-level identification
- Paper 51 §subsubsec:g4_5_v3_20_followon — current Möbius documentation (now flagged for re-examination)
- `debug/sprint_moebius_route_a_fursaev_miele_pdf_memo.md` — thread 8 Track b' literature finding (Reading B supported by literature absence, but this current finding raises whether the Möbius is even a substrate feature at the slope level)

## 8. Files

- `debug/sprint_moebius_reading_b_nrho_sweep.py` (driver)
- `debug/data/sprint_moebius_reading_b_nrho_sweep.json` (data)
- `debug/sprint_moebius_reading_b_nrho_sweep_memo.md` (this)
