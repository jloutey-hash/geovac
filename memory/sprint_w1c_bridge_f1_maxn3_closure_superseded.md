---
date: 2026-05-23
type: sprint closure
sprint: β-neg / W1c × M-Z bridge / F1 max_n=3
verdict: CLEAN NEGATIVE on three-bucket; two-bucket M-Z partition stands with kernel-energetics refinement
---

# Sprint W1c bridge + F1 max_n=3 closure

## Lead-up

Three sprints landed today (2026-05-23) in sequence:

1. **F1-P1+P2** (combined W1c × multi-zeta architecture at NaH max_n=2) returned [[PARTIAL-CLOSURE]] — multi-zeta load-bearing once W1c activates Na 3s occupation, 23% descent reduction, but PES still monotonically descending. The original Layer-3 framing from the [[α arc|sprint_modular_alpha_arc_synthesis]] revised: at max_n=2 the FCI engages Na 3s at ~0.98, NOT H-dominant as inferred.

2. **W1c × M-Z bridge sprint** applied the [[M-Z partition|sprint_modular_propinquity_mZ_bethe_log]] (bimodule cross-shift / module endomorphism, from the modular propinquity sprint) to the W1c three-layer hierarchy. Proposed a **three-bucket refinement**: cross-shift / basis-closable cross-shift / endomorphism, with sub-layer 3 in the new middle class. Three falsifiable predictions for F1 max_n=3 (P1: internal min emerges; P2: D_e within 2× of experimental; P3: mz differential persists). Predictions frozen into `debug/data/sprint_w1c_mz_partition_predictions.json` before the test ran.

3. **F1 max_n=3 prediction test** ran the three predictions with the unified architecture at Q=56. **CLEAN NEGATIVE.** P1 FAIL (R_min = 2.0 bohr smallest tested, exact falsification criterion). P2 FAIL (D_e^PES = +0.7097 Ha, 10× experimental — spurious-binding signature). P3 PASS (mz differential 0.2066 Ha at R=2.0, consistent with max_n=2's 0.21 Ha — endomorphism content preserved across basis sizes).

## Three substantive new findings

### (a) Structural success / energetic failure split

The headline finding the prediction test did not anticipate. At max_n=3 under combined W1c+mz, the dominant natural orbital IS a true bonding combination:

$$\phi_\text{dom} = -0.698 \cdot \phi_\text{Na, 3s} - 0.687 \cdot \phi_\text{H, 1s} + \text{small admixtures}$$

50/50 Na/H amplitude² split, with the 2nd NO being the antibonding combination. At max_n=2 the dominant NO was 100% H-localized at every architecture — no bonding combination at all. **The basis enlargement DID close the structural gap** predicted by the bridge sprint.

But the constructed bonding orbital has higher energy than the separated configuration at every R. The PES is monotonically descending; natural occupations are [1.0, 1.0] (open-shell singlet of bonding + antibonding singly-occupied), NOT [2.0, 0.0] (closed-shell bond). **The framework can BUILD the right orbital but cannot energetically PREFER it.**

The amplitude split is stable at ~50/50 across R ∈ [2, 5] bohr — the bonding combination is not a localized feature; the orbital structure is robust, only the energy ordering is wrong.

### (b) Sub-layer 3 reclassification

The bridge sprint's hypothesis: sub-layer 3 is **BIMODULE CROSS-SHIFT (basis-closable)** — structurally cross-shift, conditionally handled at sufficient basis size. The max_n=3 test was the falsifier.

Result: the "structural" part of the cross-shift hypothesis is CORRECT (the right linear combination IS structurally a bimodule cross-shift, and the basis enlargement DOES make it constructible). The "framework handles it" part is too strong: the framework constructs the orbital but cannot energetically prefer it.

**Sub-layer 3 reclassifies to MODULE ENDOMORPHISM (basis-irreducible at tested scales)**, joining sub-layer 2 (multi-zeta basis-shape) in the endomorphism bucket. The proposed three-bucket refinement is falsified.

### (c) Two-bucket M-Z partition stands with sharper localization

The two-bucket M-Z partition (cross-shift / endomorphism) is the correct framing. Three refinements after this arc:

1. **Cross-shift bucket** broader than originally framed — includes PK cross-center (sub-layer 1), multi-zeta acting through cross-V_ne (the multiplicative effect of multi-zeta IS via the bilinear cross-V_ne kernel), cross-register V_eN, magnetization-density, Foldy-Friar.

2. **Endomorphism bucket** has chemistry-side and QED-side flavors:
   - Chemistry-side (multi-zeta basis-shape; cross-V_ne kernel energetics): intrinsically external but framework-internally computable (FrozenCore + atomic-physics fits).
   - QED-side (LS-8a counterterms; Sprint H1 Yukawa): strictly externally specified (UV-completion physics required).

3. **The "basis-closable cross-shift" sub-class evaporates.** Basis dimensionality is not the load-bearing variable; cross-V_ne kernel energetics is. **The chemistry endomorphism class now includes cross-V_ne kernel energetics** — the framework's bilinear cross-shift machinery does not autonomously favor bonding combinations over separated configurations.

## Predictions outcome

- P1 (internal minimum at R_eq ∈ [3.0, 4.5] bohr): **FAIL** — R_min = 2.0 bohr smallest tested
- P2 (D_e ∈ [0.0375, 0.150] Ha): **FAIL** (conditional on P1) — D_e^PES = +0.7097 Ha = 10× experimental
- P3 (mz differential in [0.02, 0.30] Ha at R_eq): **PASS** — |0.2066| Ha at R = 2.0

Verdict per prespecified gate: CLEAN NEGATIVE.

## Recommended next sprint

**F2 cross-V_ne kernel-shape substitution** (Track 3's named target from May 2026, 1-2 weeks). Replaces the bare cross-V_ne kernel integration on the H partner side with physical-n hydrogenic shape — structurally distinct from this sprint's basis enlargement and from multi-zeta which acts on the Na center side. Track 3 diagnostic at max_n=2 found $-0.674$ Ha differential at R=2.5 (1.9× the W1c-residual descent), suggesting F2 would over-correct from current attraction profile, possibly producing a binding minimum.

The reading from this sprint sharpens the F2 expectation: cross-V_ne kernel energetics is now identified as the load-bearing residual. If F2 produces a binding minimum, the W1c-residual wall is mitigable by within-framework kernel substitution. If not, the wall is structurally deeper than any in-bimodule mitigation.

Parallel optional: HCl/MgH₂ scoping where the W1c-residual is empirically smaller (Z-decreasing per [[CLAUDE.md §2|claude_md]]).

DO NOT pursue: further max_n enlargement on NaH (structural gap is closed at max_n=3; remaining wall is in kernel energetics, not basis size).

## Cross-references

- [[F1-P1+P2 PARTIAL-CLOSURE|sprint_f1_p1p2_combined_test_memo]] (max_n=2 baseline)
- [[Bridge sprint with three predictions|sprint_w1c_mz_partition_analysis_memo]]
- [[F1 max_n=3 CLEAN NEGATIVE|sprint_f1_maxn3_predictions_test_memo]]
- [[M-Z partition primary source|sprint_modular_propinquity_mZ_bethe_log_memo]]
- [[β synthesis memo|sprint_beta_neg_synthesis_memo]]
- [[α arc context|sprint_modular_alpha_arc_synthesis_memo]]
- [[Track 3 diagnostic for kernel-shape target|w1c_residual_nah_track3_diag]]
- [[Chemistry arc paused W1c-residual|chemistry_arc_paused_w1c_residual]]
- [[Diagnostic before engineering rule|feedback_diagnostic_before_engineering]]
- [[Structural-skeleton-scope pattern|geovac_structural_skeleton_scope_pattern]]

## Net for the framework

The three-bucket M-Z partition hypothesis is FALSIFIED but the test is structurally informative: it surfaced the structural success/energetic failure split as the substantive content, sharpened sub-layer 3's localization, and added cross-V_ne kernel energetics to the chemistry endomorphism class. The two-bucket M-Z partition is the canonical framing.

The structural-skeleton-scope statement (CLAUDE.md §1.7) is preserved and slightly sharpened: framework determines skeleton (selection rules, transcendental classes, structural orbital constructibility); does not generate calibration data (energetic preferences, parameter values, Yukawa selection). At NaH this manifests as: the bonding orbital is structurally constructible at max_n=3, but the energetic preference for that orbital over the separated configuration requires content the framework's bimodule machinery does not provide.

Next active sprint: F2 cross-V_ne kernel-shape substitution.
