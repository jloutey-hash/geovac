# Sprint β-neg synthesis — F1 max_n=3 CLEAN NEGATIVE + W1c × M-Z bridge closure

**Date:** 2026-05-23 (post-F1-max_n=3 same-day documentation pass).
**Sprint position:** Synthesis sprint capturing the negative result from the F1 max_n=3 prediction test, the falsification of the three-bucket M-Z partition hypothesis, and the substantive reclassification of W1c sub-layer 3 from "basis-closable cross-shift" to "module endomorphism (basis-irreducible at tested scales)".
**Cross-references:** `debug/sprint_f1_p1p2_combined_test_memo.md` (F1-P1+P2 max_n=2 PARTIAL-CLOSURE baseline), `debug/sprint_w1c_mz_partition_analysis_memo.md` (bridge sprint introducing the three-bucket hypothesis and three predictions), `debug/data/sprint_w1c_mz_partition_predictions.json` (the prespecified gate with measured-value verification appended), `debug/sprint_f1_maxn3_predictions_test_memo.md` (the falsifying test).

---

## §0. Executive summary

**Headline verdict: the three-bucket M-Z partition (cross-shift / basis-closable cross-shift / endomorphism) is FALSIFIED.** The two-bucket partition (cross-shift / endomorphism) stands, sharpened by where chemistry endomorphisms localize: in the cross-V_ne kernel ENERGETICS, not in basis dimensionality.

The arc closed in three sprints across a single day (2026-05-23):

1. **F1-P1+P2** built the unified W1c × multi-zeta architecture at NaH max_n=2 and surfaced PARTIAL-CLOSURE: multi-zeta is load-bearing once W1c activates Na 3s occupation, but the PES remains monotonically descending (23% reduction in descent depth, no internal minimum). The original Layer-3 framing was revised: at max_n=2 the FCI engages Na 3s at ~0.98 occupation, but the basis lacks orbital-pair flexibility for genuine bonding/antibonding partition.

2. **W1c × M-Z bridge sprint** classified the three W1c sub-layers under the modular propinquity M-Z partition (cross-shift vs endomorphism, from `debug/sprint_modular_propinquity_mZ_bethe_log_memo.md`). The substantive bridge content was the proposed THIRD bucket: "basis-closable cross-shift" — the orbital-pair flexibility wall is structurally cross-shift (the bonding combination IS a left+right linear combination in the bimodule basis) but conditionally handled at sufficient basis size. Three falsifiable predictions for NaH max_n=3 (P1: internal min emerges; P2: D_e within 2× of experimental; P3: mz differential persists).

3. **F1 max_n=3 prediction test** ran the three predictions with the F1-P1+P2 unified architecture at Q=56. **CLEAN NEGATIVE per the prespecified gate.** P1 falsified (R_min at smallest tested R=2.0 bohr, no internal minimum), P2 falsified (D_e^PES = +0.7097 Ha, ~10× experimental NaH D_e ≈ 0.075 Ha — spurious-binding signature), P3 PASSED (mz differential 0.2066 Ha at R=2.0, consistent with max_n=2's 0.21 Ha — endomorphism content preserved across basis sizes per the partition prediction).

**The substantive new structural finding** — the substantive content the test produced that the gate did not anticipate — is the **structural success/energetic failure split**. At max_n=3 under combined W1c+mz, the dominant natural orbital IS a true bonding combination (Na 3s + H 1s with 50/50 amplitude split, mixing coefficients $-0.698, -0.687$); the 2nd NO is the antibonding combination. The basis enlargement DID close the structural gap predicted by the bridge sprint. But the constructed bonding combination has higher energy than the separated configuration at every R; the framework can BUILD the right orbital but cannot energetically PREFER it. The W1c-residual wall is now sharpened: it is NOT an orbital-pair flexibility wall (the orbital pair IS constructed at max_n=3); it IS a wall in the cross-V_ne kernel ENERGETICS — specifically, the energetic asymmetry that should make the bonding combination preferred over the separated configuration.

The next sprint is F2 (cross-V_ne kernel-shape substitution on the partner side, Track 3 named target from May 2026), which addresses a structurally distinct mechanism with a $-0.674$ Ha differential at R=2.5 (1.9× the W1c-residual descent depth) — suggesting it would over-correct from the current attraction profile and possibly produce a binding minimum.

---

## §1. The path — three sprints in one day

### 1.1 F1-P1+P2 (combined architecture, PARTIAL-CLOSURE)

The F1 sprint set the stage. It diagnosed the architectural gap that prevented W1c × multi-zeta from running combined (a defensive `NotImplementedError` in `geovac/balanced_coupled.py` line 748) and removed it cleanly: the bare-Coulomb and screening-correction sub-integrals both compose at the integrand level, so the unified screened+multi-zeta path is a 25-line addition to `geovac/cross_center_screened_vne.py` with a `multi_zeta_basis: Optional[Dict] = None` kwarg.

The 8-point PES scan at max_n=2 revealed three structural findings:

(a) Multi-zeta is bit-zero on the FCI in the bare cross-V_ne case (reproducing α-PES Step 2), but **fully load-bearing when combined with W1c** (differential +0.056 Ha at R=3.5 scaling to +0.208 Ha at R=2.0). The bit-zero finding was an artifact of multi-zeta being applied to an FCI state with no Na 3s occupation.

(b) The original Layer-3 framing was revised: at W1c, the FCI natural occupations are [1.0, 1.0] (open-shell-singlet-like with Na 3s occupation ~0.98), NOT H-dominant as originally inferred from the α arc.

(c) The combined architecture does NOT close the binding gap — descent depth reduces from W1c-alone 0.898 Ha to W1c+mz 0.690 Ha (23%), but the PES is still monotonically descending; experimental NaH D_e ≈ 0.075 Ha is still ~10× smaller than the over-attraction.

The verdict was PARTIAL-CLOSURE-AT-MAX_N=2 with three named follow-ons. The bridge sprint emerged from the question: does the M-Z bimodule cross-shift / module endomorphism partition (from the modular propinquity sprint earlier in the same day) extend to the W1c three-layer hierarchy F1 surfaced?

### 1.2 W1c × M-Z bridge sprint (three-bucket hypothesis)

The bridge sprint applied the M-Z partition explicitly to each W1c sub-layer:

| Sub-layer | Mechanism | Bridge classification |
|:---|:---|:---|
| 1 — H ↔ Na core orthogonality | Phillips-Kleinman cross-center barrier | **Bimodule cross-shift** (handled at 14.6% W1c-descent reduction) |
| 2 — Na-side wavefunction shape | Multi-zeta substitution of Na 3s/3p | **Module endomorphism** (bit-zero without cross-shift activator; load-bearing once cross-shift activates Na occupancy — 23% additional reduction) |
| 3 — Orbital-pair flexibility (no bonding/antibonding mixing at max_n=2) | FCI cannot construct [H 1s ± Na 3s]_± combinations at limited basis | **Basis-closable bimodule cross-shift** — structurally cross-shift, conditionally handled at sufficient basis size |

The substantive content of the bridge was the third bucket. The M-Z partition as originally stated has two buckets; sub-layer 3's mechanism (the bonding combination is structurally a bimodule cross-shift but the FCI cannot realize it at max_n=2 due to dimensional poverty) suggested a sub-classification within the cross-shift bucket. The natural test: does max_n=3 (Q=56, adding H 2s/2p, Na 3p/3d giving ~3× more orbitals per center) provide the dimensional richness?

Three falsifiable predictions emerged:

- **P1 (internal minimum emerges):** $R_\text{eq} \in [3.0, 4.5]$ bohr. Falsifier: $R_\text{min}$ at the smallest tested R.
- **P2 (binding within 2×):** $D_e \in [0.0375, 0.150]$ Ha. Falsifier: outside this range.
- **P3 (mz differential persists):** $|E_\text{W1c} - E_\text{W1c+mz}| \in [0.02, 0.30]$ Ha. Falsifier: <0.005 Ha (mz absorbed) or >0.30 Ha (anomalous scaling).

These predictions were prespecified and frozen into `debug/data/sprint_w1c_mz_partition_predictions.json` before the F1 max_n=3 sprint ran — the curve-fit-audit discipline applied at the inter-sprint scope.

### 1.3 F1 max_n=3 (CLEAN NEGATIVE)

The F1 max_n=3 sprint executed the three predictions with the unified architecture at Q=56. Result per the prespecified gate:

| Prediction | Statement | Actual | Verdict |
|:--|:--|:--|:--:|
| P1 | Internal PES minimum at $R_\text{eq} \in [3.0, 4.5]$ bohr | $R_\text{min} = 2.0$ bohr (smallest tested) | **FAIL** |
| P2 | $D_e \in [0.0375, 0.150]$ Ha within 2× of experimental | $D_e^\text{PES} = +0.7097$ Ha (10× experimental) | **FAIL** (conditional on P1) |
| P3 | $|E_\text{W1c} - E_\text{W1c+mz}| \in [0.02, 0.30]$ Ha at $R_\text{eq}$ | $|0.2066|$ Ha at $R = 2.0$ | **PASS** |

1 of 3 predictions passed — CLEAN NEGATIVE per the bridge sprint's gate rules. The PES at max_n=3 is monotonically descending across $R \in \{2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0\}$ bohr, exactly the falsification criterion for P1.

The descent depth reduction profile is preserved across basis sizes: W1c-alone 0.916 Ha → W1c+mz 0.710 Ha (22% reduction at max_n=3, comparable to 23% at max_n=2). Multi-zeta endomorphism content is bounded per matrix element and remains visible at larger basis sizes — exactly the P3 prediction, the only prediction that survives.

---

## §2. The substantive new structural finding

### 2.1 Structural success / energetic failure split

The headline finding the prediction test did not anticipate: at max_n=3 under combined W1c+mz, the framework DOES construct the right bonding orbital. The dominant natural orbital at R=3.5 has amplitude:

$$\phi_\text{dom} = -0.698 \cdot \phi_\text{Na, 3s} - 0.687 \cdot \phi_\text{H, 1s} + \text{small admixtures}$$

with Na/H amplitude² split 0.50/0.50. The 2nd natural orbital is the antibonding combination $\phi_\text{2nd} = -0.698 \cdot \phi_\text{Na, 3s} + 0.687 \cdot \phi_\text{H, 1s}$. Compare with max_n=2, where the dominant NO at every architecture was 100% H-localized — no bonding combination at all.

The basis enlargement DID provide the orbital-mixing flexibility predicted by the bridge sprint. From a pure orbital-structure standpoint, the bridge sprint's hypothesis is structurally CONFIRMED — the bonding orbital is constructible at max_n=3 in a way it wasn't at max_n=2.

**But the constructed bonding orbital has higher energy than the separated configuration at every tested R.** The PES is monotonically descending; no internal minimum emerges. The natural occupations remain [1.0, 1.0] (open-shell singlet of bonding + antibonding singly-occupied), NOT [2.0, 0.0] (closed-shell bond where the bonding pair is doubly occupied).

This is the structural success / energetic failure split. The framework can BUILD the right orbital. The framework cannot energetically PREFER it over the separated configuration. The two halves of "binding" — orbital constructibility and energetic preference — split cleanly at this resolution.

### 2.2 The orbital amplitudes are essentially constant across R

The Na/H amplitude² split is stable at ~50/50 across most of the tested R range:

| R [bohr] | dom NO Na/H amp² |
|:--:|:--|
| 2.0 | 0.498 / 0.502 |
| 2.5 | 0.501 / 0.499 |
| 3.0 | 0.501 / 0.499 |
| 4.0 | 0.503 / 0.497 |
| 5.0 | 0.500 / 0.500 |
| 7.0 | 0.470 / 0.530 |
| 10.0 | 0.389 / 0.611 |

The bonding combination is not a localized feature of R=3.5 — it persists across the PES. The slight tilt at R=10 is the dissociation-limit weakening as the bonding combination dissolves. The orbital structure is robust; only the energy ordering is wrong.

### 2.3 Sub-layer 3 reclassification

Under the bridge sprint's three-bucket hypothesis, sub-layer 3 was classified as **basis-closable cross-shift**: structurally cross-shift (the bonding combination is a left+right bimodule linear combination), basis-conditionally handled. The max_n=3 test was the falsifier for this classification.

The result splits the classification:

- **The "structural" part of the cross-shift hypothesis is correct.** The right linear combination IS structurally a bimodule cross-shift, and the basis enlargement DOES make it constructible.
- **The "framework handles it" part is too strong.** The framework constructs the orbital but cannot make it the energetically preferred ground-state configuration.

The net reclassification: sub-layer 3 lives in the **same external-input class as sub-layer 2**, not in a separate basis-closable cross-shift sub-class. The proposed three-bucket refinement is falsified. The original two-bucket M-Z partition (cross-shift / endomorphism) is the correct framing.

### 2.4 What the two-bucket partition retains, sharpened

The two-bucket partition stands, but with a sharper understanding of where chemistry endomorphisms live:

- **Cross-shifts (handled by the framework):** Phillips-Kleinman cross-center (sub-layer 1, 14.6% descent reduction); multi-zeta substitution acting through the cross-V_ne integration kernel (the multiplicative effect of multi-zeta is via the bilinear cross-V_ne, which IS a cross-shift); cross-register V_eN recoil; magnetization-density / Zemach kernel; Foldy-Friar charge-density.
- **Endomorphisms (external input):** LS-8a renormalization counterterms $Z_2 - 1$ and $\delta m$ (QED side, strictly externally specified); Sprint H1 Yukawa $Y$ (SM unification, strictly externally specified); the multi-zeta basis-shape data itself (chemistry side, framework-internally computable via FrozenCore + atomic-physics fits); **the cross-V_ne kernel energetic asymmetry that would prefer the bonding combination over the separated configuration** (chemistry side, named structurally for the first time by this sprint).

The new entry in the endomorphism bucket is the substantive sharpening: until this sprint, the chemistry endomorphism content was localized to basis-shape data (sub-layer 2, multi-zeta). After this sprint, it extends to the cross-V_ne kernel energetics itself — the framework's bilinear cross-shift kernel does not autonomously favor bonding combinations over separated configurations.

The QED side preserves its strict-external-input character (UV-completion physics required for LS-8a counterterms). The chemistry side has both flavors: basis-shape endomorphism is mitigable via FrozenCore (computable from autonomous atomic-physics solvers); kernel-energetics endomorphism is the genuinely new external-input class.

---

## §3. What the M-Z partition retains

### 3.1 Two buckets, sharpened localization

The M-Z partition retains its two-bucket structure (cross-shift / endomorphism), but with three substantive refinements after this arc:

1. **Cross-shift bucket is broader than originally framed.** PK cross-center, multi-zeta acting through cross-V_ne, cross-register V_eN, magnetization-density, Foldy-Friar charge-density all sit in this bucket. The framework's bilinear cross-shift machinery (Roothaan 1951 cross-register kernel + Shibuya-Wulfman multi-zeta dispatch + Latrémolière bridge construction) handles all of these.

2. **Endomorphism bucket has chemistry-side and QED-side flavors.** Chemistry-side endomorphisms (multi-zeta basis-shape data, cross-V_ne kernel energetics) are intrinsically external to the framework but framework-internally computable (FrozenCore radial Schrödinger + atomic-physics basis fits). QED-side endomorphisms (LS-8a counterterms, Sprint H1 Yukawa) are strictly externally specified (UV-completion physics required, not in any GeoVac extension).

3. **The "basis-closable cross-shift" sub-class evaporates.** Even when the basis is large enough for the FCI to construct the right cross-shift orbital combination (max_n=3 at NaH), the framework cannot make that combination energetically preferred. Basis dimensionality is not the load-bearing variable; the cross-V_ne kernel energetics is.

### 3.2 Consistency with the broader framework

The CLEAN NEGATIVE verdict is fully consistent with the structural-skeleton-scope finding (CLAUDE.md §1.7): the framework determines the skeleton (selection rules, transcendental classes, scaling laws, structural orbital constructibility) but does not generate calibration data (energetic asymmetries, parameter values, Yukawa selection).

At NaH this manifests as: the bonding orbital is structurally constructible, but the energy that would make it the preferred ground-state configuration requires inputs beyond what the framework's bimodule machinery (Phillips-Kleinman + cross-V_ne + multi-zeta basis) provides.

The five-instance multi-focal-composition wall (Sprint H1 Yukawa, LS-8a Z_2/δm, HF-3 recoil, HF-4 Zemach, HF-5 multi-loop a_e) gets a sixth instance: the chemistry-side cross-V_ne kernel energetic asymmetry for second-row hydrides. The framework couples discrete labels cleanly via Wigner symbols and selection rules; the Fock projection couples space; the framework has no native composition theorem for the energetic preference between bonding and separated configurations when both are constructible from existing modules.

---

## §4. What's next

### 4.1 F2 cross-V_ne kernel-shape substitution (primary)

Track 3's named target P2 (from the May 2026 PK cross-center synthesis memo): replace the bare cross-V_ne kernel integration on the H partner side with physical-n hydrogenic shape, distinct from the multi-zeta basis on the Na center side that this sprint tested. The Track 3 diagnostic at NaH max_n=2 (`debug/w1c_residual_nah_track3_diag.py`) found a differential of $-0.674$ Ha at R=2.5 — about 1.9× the W1c-residual descent depth.

The mechanism is structurally distinct from this sprint's basis enlargement: F2 addresses the cross-V_ne kernel shape on the OPPOSITE side from where multi-zeta lives. Multi-zeta substitutes the Na-side basis function; F2 substitutes the H-partner-side basis function used INSIDE the cross-V_ne integrand. The Track 3 diagnostic suggests F2 would over-correct from the current attraction profile, possibly producing a binding minimum — at least at max_n=2 where the wall is empirically smaller.

The reading from this sprint sharpens the F2 expectation: the cross-V_ne kernel energetics is now identified as the load-bearing residual (not orbital-pair flexibility). F2 addresses kernel shape directly. If F2 produces a binding minimum, the W1c-residual wall is mitigable by within-framework kernel substitution (the M-Z partition becomes: cross-shifts handled, chemistry endomorphisms partially mitigable through both basis-shape AND kernel-shape substitution). If F2 does not produce a binding minimum, the wall is structurally deeper than any in-bimodule mitigation can address.

Compute cost: 1-2 weeks for production wiring + tests + max_n=2 verification. The diagnostic ($-0.674$ Ha differential) is at max_n=2 where the basis enlargement is not the variable; the test should run at max_n=2 first to isolate the kernel-shape mechanism from the basis-enlargement structural success.

### 4.2 HCl scoping (parallel, optional)

The W1c-residual wall is empirically Z-decreasing: NaH 5.4-6.0×, MgH₂ 2.99×, HCl 1.79× (CLAUDE.md §2). At higher Z, the wall may already be smaller fraction of the bond energy scale, so the framework may have effective binding even without complete W1c closure. A 2-day scoping run of HCl max_n=2 with the combined W1c × multi-zeta architecture would test whether the framework already supports second-row hydride binding for less-residual systems.

The multi-zeta registry would need extension to Cl (Z=17), which is a separate multi-day sub-sprint. Lower priority than F2 because the structural finding (kernel-energetics is the residual) is the same; HCl tests whether the residual magnitude crosses the bond-energy threshold, not whether the structural mechanism is mitigable.

### 4.3 Do NOT pursue: further max_n enlargement on NaH

At max_n=3 the structural gap (bonding orbital constructibility) is closed; the wall now lives in the energetic asymmetry, which is not basis-size limited at this scale. Further enlargement (max_n=4 at Q=84) would add more orbitals but would not address the energetic gap. Same structural finding would persist with marginal numerical refinement, at substantially higher FCI cost.

### 4.4 Open question: Hylleraas-class explicit-correlation for molecular bimodules?

A speculative direction surfaced by this sprint: the Hylleraas-Eckart double-α architecture (closed on He singlet at 0.0006%, He 2¹P→1¹S oscillator strength at −2.02%; CLAUDE.md §2) explicitly captures the electron-electron correlation that single-determinant FCI cannot. The NaH W1c-residual is a different mechanism (single-particle cross-V_ne kernel energetics, not electron-electron correlation), but the architectural pattern is suggestive: when single-particle bimodule machinery cannot energetically prefer a structurally constructible orbital combination, explicit-correlation methods might bridge the gap. Out of sprint scope; flagged for completeness.

---

## §5. Honest scope: confidence by claim

| Claim | Confidence | Reason |
|:---|:---:|:---|
| At max_n=3, the dominant NO under combined W1c+mz IS a true bonding combination (Na 3s + H 1s, 50/50 split) | HIGH | Direct numerical observation across 8 R-points; amplitude split stable; structural success/energetic failure split is robust |
| Three-bucket M-Z partition (cross-shift / basis-closable / endomorphism) is FALSIFIED | HIGH | Prespecified gate with predictions frozen before test ran; P1 fails on exact falsification criterion (R_min at smallest tested R) |
| Two-bucket M-Z partition (cross-shift / endomorphism) stands | HIGH | All three F1 sub-layers still classify; chemistry endomorphism class now includes both basis-shape and kernel-energetics |
| Sub-layer 3 reclassifies to genuine endomorphism class | MEDIUM-HIGH | Reclassification is the inverse of the falsified hypothesis; alternative reading (sub-layer 3 is hybrid cross-shift/endomorphism) is possible but more cumbersome |
| Cross-V_ne kernel energetics is the residual wall | MEDIUM-HIGH | Inferred from the structural success/energetic failure split; F2 test is the direct probe |
| F2 (cross-V_ne kernel-shape substitution) is the natural next sprint | MEDIUM | Track 3 diagnostic suggests F2 would over-correct from current attraction profile; possibility of producing binding minimum is high but not certain |
| HCl/MgH₂ would show qualitatively similar behavior with smaller magnitude | MEDIUM | Z-decreasing pattern is empirical (CLAUDE.md §2), not derived from first principles |
| Hylleraas-class explicit-correlation might bridge the chemistry endomorphism gap | LOW | Speculative; the mechanism is different from electron-electron correlation; flagged for completeness |

The high-confidence claims are the substantive output of this sprint. The lower-confidence claims are recommendations for next steps, not load-bearing conclusions.

---

## §6. Session summary

### Tracks
- F1 max_n=3 prediction test: **CLOSED — CLEAN NEGATIVE** (P1 FAIL, P2 FAIL, P3 PASS)
- W1c × M-Z partition bridge sprint: **CLOSED — three-bucket FALSIFIED, two-bucket stands**
- Sub-layer 3 reclassification: **CLOSED — basis-closable cross-shift → module endomorphism**

### Results
- Aggregate predictions verdict: 1/3 PASS
- Max_n=3 PES descent depth: W1c-alone 0.916 Ha → W1c+mz 0.710 Ha (22% reduction, comparable to max_n=2's 23%)
- Multi-zeta differential at max_n=3: 0.206 Ha at R=2.0, 8.5×10⁻⁵ at R=10 (consistent with max_n=2 magnitudes)
- Dominant natural orbital at max_n=3 (combined W1c+mz): $-0.698 \cdot \text{Na 3s} - 0.687 \cdot \text{H 1s}$ + small admixtures
- Natural occupations: [1.0, 1.0] (open-shell singlet), NOT [2.0, 0.0] (closed-shell bond)
- Na/H amplitude² split: stable ~50/50 across R ∈ [2, 5] bohr

### Files Modified (paper edits — see §6 of deliverable list)
- `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` — appended F1 arc subsubsection to §sec:w1c_residual
- `papers/group2_quantum_chemistry/paper_17_composed_geometries.tex` — added cross-reference to Paper 19 W1c-residual extension
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — added §V.D.9 three-bucket M-Z partition falsification entry
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` — added §VIII Sprint M-Z addendum remark
- `CLAUDE.md` — §2 sprint bullet + §3 dead-ends row (three-bucket hypothesis falsification)

### Files Created
- `memory/sprint_w1c_bridge_f1_maxn3_closure.md` — closure memory file with [[wikilinks]]
- `debug/sprint_beta_neg_synthesis_memo.md` — this synthesis memo

### Decisions
- Two-bucket M-Z partition (cross-shift / endomorphism) is the canonical framing; three-bucket refinement is dropped
- Chemistry endomorphism class extends beyond basis-shape (sub-layer 2) to cross-V_ne kernel energetics (sub-layer 3 reclassified)
- Next sprint: F2 cross-V_ne kernel-shape substitution at NaH max_n=2 (1-2 weeks; the wall is now identified as kernel-energetics, not basis-size)
- No production code modifications in this synthesis sprint (paper edits + CLAUDE.md updates only)
- DO NOT autonomously claim the W1c wall is closed or that F2 will close it — that's the next sprint's job

---

**End of Sprint β-neg synthesis memo.** Verdict: three-bucket M-Z partition hypothesis FALSIFIED via the prespecified F1 max_n=3 gate; two-bucket M-Z partition stands with sharper localization of chemistry endomorphism class (now includes cross-V_ne kernel energetics). The structural success/energetic failure split is the substantive new content; F2 cross-V_ne kernel-shape substitution is the natural next target.
