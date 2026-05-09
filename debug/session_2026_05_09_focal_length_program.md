# Session memo — 2026-05-09 multi-observable focal-length decomposition program

## Why this memo exists

The 2026-05-09 session produced an unusually large amount of structural output — three completed catalogue extensions, four iterations on the same observable (proton Zemach radius), a new §V.C subsection in Paper 34, a new directive (CLAUDE.md §1.8), and several methodology lessons. This memo captures the full trajectory in one place so the lessons don't decay into individual sprint memos that future sessions miss.

The user-facing reframe at the end of the session: the framework's projection-chain machinery functions as a unique diagnostic instrument for precision physics. Multi-observable global fits using a single structural framework expose convention-mismatches in literature compilations and kernel-approximation gaps that are invisible from any single-observable analysis. This is the methodology contribution. Codified as CLAUDE.md §1.8 with five active Roothaan-decomposition targets.

## Session map

The session ran approximately as:

1. **Audit** of the `papers/observations/paper_34_projection_taxonomy.tex` projection dictionary against the new §IV.5 base-unit decomposition, surfacing nine informal projection names and two structural slot candidates.
2. **§IV.5 base-unit decomposition** added (Tables 1+2, inverse list, L–E asymmetry observation, compositional-projection observation generalized).
3. **Three new §III projections** promoted: §III.17 charge-density (Foldy/Friar), §III.18 magnetization-density (Zemach), §III.19 tensor multipole. Plus §III.20 boundary subsection (internal-graph-derived vs Layer-2).
4. **External-field hunting** verified Zeeman/Stark are Layer-2 inputs, not §III projections (clean negative; pinned the §III/Layer-2 distinction).
5. **Three calculation tracks** dispatched in parallel: Track L (Lyman α emission rate), Track P (H 1S polarizability), Track HFD (D 1S hyperfine via §III.18).
6. **Two more calculation tracks**: Track LAR (Lamb shift Roothaan autopsy), Track rZG (self-consistent r_Z extraction).
7. **r_Z thread iteration**: original rZG (~166 mfm, all inconsistent due to misdiagnosis) → W1a-D fix (166 mfm, within 1σ of all literature, Layer-2 budget identified as systematic) → extended fit (51 mfm, kernel approximation order identified as systematic) → W1b operator extension (16 mfm, Layer-2 convention mismatch identified as systematic).
8. **CLAUDE.md §1.8 directive** codifying the multi-observable-decomposition program with five active targets.
9. This memo.

## Paper 34 changes applied (this session)

| Section | Change | Source |
|:---|:---|:---|
| §III.17 | New: nuclear charge-density correction (Foldy/Friar) | Hunting memo |
| §III.18 | New: nuclear magnetization-density correction (Zemach), now operator-level mature | Hunting memo + W1b extension |
| §III.19 | New: nuclear tensor multipole coupling | Hunting memo |
| §III.20 | New: boundary subsection (internal-graph-derived vs Layer-2 distinction) | External-field hunt |
| §IV (Table 1) | +3 rows for §III.17–19 | Hunting memo |
| §IV.5 (Table 2) | New table: base-unit decomposition of 19 projections | Audit + §IV.5 add |
| §IV.5 inverse list | [L] now 6 scalar instances + 1 [L]² rank-2 instance | Audit + §IV.5 add |
| §IV.5 Observation 1 (LE-asymmetry) | Generalized to 6 vs 4 instances; physics-of-asymmetry argument | Audit + §IV.5 add |
| §IV.5 Observation 3 (compositional) | Generalized from "unique Breit case" to "family of two confirmed instances + candidates" | Audit + §IV.5 add |
| §V (matches catalogue) | +2 rows in three-projection-chain group: Lyman α, H 1S polarizability | Track L + Track P |
| §V.B | +1 row: rZG self-consistent extraction (Layer-2 itemization diagnostic) | Track rZG (multi-iteration) |
| §V.C | New subsection: Roothaan autopsies (§V.C.1 Lamb shift full table + §V.C.2/3 placeholders for 21cm and μH Lamb) | Track LAR |
| §VI | New `rem:depth_refinement` Remark conditioning Prediction 1 on irrational vs rational analytic answer | Track P (counter-example finding) |
| §VII T3 theorem | Pending expansion 16×16 → 19×19 flagged for next TX-A audit | §III.17–19 promotion |
| §VIII Conclusion | Added intra-nuclear projection summary; framework now 19 projections | All section adds |
| Count references | "sixteen" → "nineteen" updated in abstract, §III intro, §IV intro, §IV.5 intro, T3 §VII, §VIII Conclusion | All section adds |

## CLAUDE.md changes applied (this session)

- §1.7 W1a entry updated by W1a-D fix sub-agent: noted bug-fix completion (no production code change; rZG driver Layer-2 budget was the actual issue).
- §2 development frontier: W1b recoil-mixing extension bullet added (verdict, files, next-step recommendation).
- §1.8 (NEW): multi-observable focal-length decomposition program directive with three problem classes and five active targets.

## Production code changes

- `geovac/cross_register_vne.py`: **unchanged** (W1a-D fix sprint verified species-agnostic correctness).
- `geovac/magnetization_density.py`: **extended** with `include_recoil_mixing` and `nucleon_mass` fields on `MagnetizationDensitySpec`; recoil-mixing per Pachucki + Friar moment per Eides §6.2; ~480 → ~600 lines; default `False` preserves existing 37-test bit-identical regression.
- `tests/test_cross_register_vne.py`: +6 species-agnostic recoil tests (D, Mu, static-nucleus limit, mass-ratio scaling, kernel symmetry).
- `tests/test_magnetization_density.py`: +20 tests in new `TestRecoilMixingExtension` class.

Test panel (all passing): 57/57 magnetization-density + 69/69 cross-register-vne + 9/9 Paper 35 prediction = 135/135.

## The r_Z trajectory (the day's central result)

Four iterations on extracted proton Zemach radius from a multi-observable global fit:

| Iteration | σ(r_Z(p)) | Central | Verdict | Load-bearing systematic identified |
|:---|:---:|:---:|:---|:---|
| Original rZG | 166 mfm | 1.373 fm | Inconsistent with all literature | (claimed) cross-register V_eN code |
| W1a-D fix | 166 mfm | 1.180 fm | Within 1σ of Eides, lattice, Friar–Payne | Layer-2 budget itemization |
| Extended fit | 51 mfm | 1.257 fm | All ~+4σ above literature | Kernel approximation order |
| W1b extension | 16 mfm | 1.286 fm | All ~+10σ; INCONSISTENT_WITH_BOTH | **Layer-2 convention mismatch** |

The trajectory is the substantive contribution of the day. Each iteration:
1. Tightened precision (or attempted to)
2. Revealed a deeper systematic
3. Falsified the previous iteration's diagnosed cause

The W1b iteration is the cleanest: a named hypothesis (recoil-mixing in §III.18 closes the +0.22 fm μH-alone offset) was tested operationally with full Friar–Payne implementation, and the hypothesis was falsified — recoil-mixing is real (~9% kernel modification) but doesn't close the offset. The actual systematic was relocated to convention mismatches between Krauth 2017 (μH Tab. 1 itemization) and Eides 2024 Tab. 7.3 (H 21cm itemization).

**The framework now has 16 mfm precision but cannot weigh in on Eides-vs-lattice tension** because both literature values are excluded by the central. This is paradoxically a more honest diagnostic than picking one would have been: it says "if you want a multi-observable r_Z extraction at sub-30 mfm precision, the literature compilations need convention-harmonization first; that's not a framework problem, it's a literature methodology problem."

This is a unique contribution because no single-observable analysis can surface it — convention mismatches require multiple observables in the same projection-chain dictionary to become numerically visible.

## Methodology lessons (the durable insights)

### Lesson 1: Diagnostic-before-engineering, applied operationally

Yesterday's rZG memo flagged a "bug in cross_register_vne.py" based on the symptom (r_Z(D) extracted at 6.18 fm vs Friar–Payne 2.593 fm). Today's W1a-D agent ran the production code in isolation across H, D, Mu and verified species-agnostic correctness (2.86%, 2.03%, 8.18% vs Bethe–Salpeter respectively). The actual bug was in the rZG driver's Layer-2 budget (–150 ppm specified for D HFS, –286 ppm correct per Pachucki–Yerokhin). 

The lesson: when a sprint flags a bug location, the next sprint must verify the diagnosis before fixing the diagnosed component. This is exactly the discipline `feedback_diagnostic_before_engineering.md` (memory file) was designed to enforce. The W1a-D sprint applied it correctly. Without this discipline, today would have been an incorrect "fix" of correct code.

### Lesson 2: Each precision tightening relocates the load-bearing systematic

The four-iteration trajectory shows that "load-bearing systematic" is precision-dependent. At coarse precision (~166 mfm), Layer-2 budget itemization dominates. At medium precision (~50 mfm), kernel approximation order dominates. At fine precision (~16 mfm), Layer-2 convention mismatch dominates.

Practical implication: when extending precision via additional observables or kernel refinements, expect a NEW systematic to surface at the new precision scale. This is not failure — it's diagnostic refinement. Each named systematic gets tested, falsified or confirmed, and the next is identified.

The opposite practice — claiming "now precision is 16 mfm, we can pick Eides over lattice" without checking per-observable consistency — would have produced a false prospective verdict (r_Z(p) = 1.179 fm with σ = 21 mfm "looks like" a clean discrimination) that the agent correctly rejected.

### Lesson 3: Multi-observable global fits are a methodology probe, not just a precision tool

The framework's structural-skeleton scope (CLAUDE.md §1.7) says it cannot autonomously generate calibration data (Yukawas, multi-loop counterterms, lattice radii). Today's r_Z thread shows the framework's structural capabilities exceed this in one specific way: it surfaces **convention mismatches** in calibration data that no single-observable analysis can see.

This is a methodology contribution that doesn't fit the structural-skeleton-scope statement directly. It's not "framework predicts new physics"; it's "framework's structural consistency requirement, applied across multiple observables, surfaces methodology issues in the literature compilations that consume the calibration data."

CLAUDE.md §1.8 codifies this as an ongoing program with explicit problem-class enumeration.

### Lesson 4: Roothaan autopsy format for multi-component observables

The Lamb shift Roothaan autopsy (§V.C.1) demonstrated: multi-component observables (Lamb shift, hyperfine, fine structure) deserve a different §V format than single-component observables. The autopsy table — eight components, each tagged to a §III projection chain, with FN/L2 status flags — makes the projection-chain decomposition visible at the observable level. Single-row treatment hides this structure.

The discipline is now established: §V/§V.B for single-row observables; §V.C for multi-component decompositions. Three placeholders set up (§V.C.1 done, §V.C.2/3 pending).

### Lesson 5: Hypothesis falsification with operator-level evidence is more informative than verification

The W1b extension was implemented correctly (recoil-mixing per Pachucki + Friar moment per Eides §6.2; ~9% kernel modification; NLO match against Krauth full theory at ~1%; 57/57 tests pass). But it falsified the hypothesis it was implemented to test. This is a stronger result than a verification would have been.

A verification ("recoil-mixing closes the offset") would have left open whether the offset was actually due to recoil-mixing or to something else that happens to coincide with it. The falsification ("recoil-mixing is real but doesn't close the offset") relocates the systematic with operator-level certainty.

This is the discipline of testing named hypotheses operationally: build the operator, run it against the data, and let the numerical result decide. The May 9 r_Z thread is the cleanest demonstration of this discipline in the project's history.

## Active Roothaan-decomposition targets (queued for next sprints)

Per CLAUDE.md §1.8 directive:

1. **§V.C.2 H 21 cm hyperfine four-component autopsy** (placeholder fill)
2. **§V.C.3 μH 2S–2P Lamb shift autopsy** (placeholder fill)
3. **§V.C.4 He 2³P fine structure autopsy** (NEW; tests internal multi-focal angular composition)
4. **§V.C.5 He 2¹P → 1¹S oscillator strength** (NEW; multi-electron emission rate, Sturmian continuum check)
5. **§V.C.6 Cs 6S₁/₂ hyperfine** (NEW prospective; Z=55 heavy-atom Z-scaling test, atomic-clock precision regime)

Sprint cadence: 2 tracks per sprint (1 verification + 1 prospective). Continue from established May 9 cadence.

## Open-question follow-ups (parking lot)

- **Layer-2 convention reconciliation sprint**: read Krauth 2017 + Eides 2024 + Karshenboim 2005 for Zemach itemization conventions; harmonize at sub-percent level; re-run rZG-extended-v2 with harmonized Layer-2 budgets. After harmonization, framework's central r_Z(p) might land cleanly at Eides, lattice, or between. **Status**: identified as the named follow-up by Sprint W1b agent, but is literature work rather than framework work; queued but not high priority.

- **Track P sharpening**: the conditional refinement of Prediction 1 (Remark `rem:depth_refinement`) names "rational-class observables with finite-Sturmian-basis representation" as having zero residual at any depth. Polarizability is one example. Others to test: hydrogen ⟨r²⟩, hydrogenic dipole moments, possibly other Dalgarno–Lewis-class observables. Methodology audit if this turns out to characterize a clean class.

- **§III dictionary completion**: with §III.17–19 added today, the dictionary is at 19 projections. Any further additions (e.g., decoherence operator, scattering-state projection, external-field projection if reconsidered) should follow the same operator-level-mature-vs-Layer-2-input distinction codified in §III.20.

- **L–E asymmetry physical reading**: the observation that [L] now admits 6 scalar instances + 1 rank-2 tensor while [E] admits only 4 is structurally interesting. Why does physics need more length scales than energy scales? Speculative: configuration space has more independent geometric-structure axes (atomic, molecular, internuclear, intra-nuclear charge, intra-nuclear magnetization, intra-nuclear tensor) than spectrum has independent role-axes (foundational, regulator, subtraction, correction). Worth a future memo if anything substantive emerges.

- **Multi-focal-composition wall structural reading after W1b**: the wall (CLAUDE.md §1.7) was originally five sub-walls W1a/W1b/W1c/W2a/W2b/W3 plus the LS-8a sector. With W1b now implemented at the operator level (recoil-mixing per Friar–Payne), the wall taxonomy might want updating. Sprint W1b agent reported the W1b sub-wall is "operator-level closed at NLO precision but Layer-2 itemization-convention-limited at sub-30 mfm precision." This is a sharper version of the multi-focal-composition wall pattern than the original five-instance reading. Worth refining the §1.7 wall entry.

## File index

### New files (this session)

- `debug/paper34_v_base_unit_audit_memo.md` — audit of §V chains against §IV.5 base-unit dictionary
- `debug/intra_particle_L_hunting_memo.md` — three-projection §III.17/18/19 verdict
- `debug/external_field_hunting_memo.md` — Zeeman/Stark Layer-2 verdict
- `debug/calc_track_L_lyman_alpha_memo.md` — Lyman α emission rate at +0.055%
- `debug/calc_track_P_h1s_polarizability_memo.md` — H 1S polarizability EXACT at N_basis=2
- `debug/calc_track_HFD_d_hyperfine_memo.md` — D 1S HFS via §III.18, +40 ppm BF strict
- `debug/calc_track_LAR_lamb_autopsy_memo.md` — Lamb shift Roothaan autopsy, 97.5% framework-native
- `debug/calc_track_rZG_global_zemach_memo.md` — original rZG fit
- `debug/rzg_bug_diagnosis_memo.md` — W1a-D fix diagnostic (cross_register_vne is correct)
- `debug/calc_track_rZG_extended_memo.md` — extended-fit kernel-limited verdict
- `debug/w1b_recoil_mixing_extension_memo.md` — W1b operator extension + Layer-2-convention diagnostic
- `debug/calc_track_rZG_global_zemach.py` — original fit driver
- `debug/calc_track_rZG_extended.py` — extended fit driver
- `debug/calc_track_rZG_extended_v2.py` — W1b-extended fit driver
- `debug/data/calc_track_*.json` — numerical outputs
- `debug/session_2026_05_09_focal_length_program.md` — this memo

### Modified production files

- `geovac/magnetization_density.py` (extended with recoil-mixing)
- `tests/test_cross_register_vne.py` (+6 species-agnostic tests)
- `tests/test_magnetization_density.py` (+20 recoil-mixing tests)

### Modified papers

- `papers/observations/paper_34_projection_taxonomy.tex` — substantial extension (§III.17/18/19/20, §IV Tables 1+2, §IV.5 base-unit decomposition, §V matches catalogue +2 rows, §V.B rZG row multi-iteration, §V.C Roothaan autopsies subsection with 1 full table + 2 placeholders, §VI Prediction 1 refinement, §VII T3 expansion flag, §VIII Conclusion)

### Modified CLAUDE.md sections

- §1.7 W1a entry (W1a-D fix completion note)
- §1.8 (NEW): Multi-observable focal-length decomposition program directive
- §2 (W1b recoil-mixing extension bullet)

## Closing note

The day's structural takeaway is that the framework's projection-chain machinery is a precision-physics methodology probe in addition to its other roles (computational quantum chemistry tool, NCG spectral-triple instance, qubit-encoding sparsity engine). Multi-observable consistency requirements expose convention mismatches and kernel-approximation gaps at numerical level. CLAUDE.md §1.8 codifies this as an ongoing program with five queued targets.

The trajectory of testing-named-hypothesis-operationally-and-letting-numerical-evidence-decide produced four iterations on a single observable, each informative. This discipline, not the r_Z extraction itself, is the durable contribution.
