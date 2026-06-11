# Sprint β-update — Append F4/F5/F6 to F3-maturity baseline

**Date**: 2026-05-23 (same-day continuation of F3-maturity β-comprehensive
sprint that landed F1-F2-F3 closure paperwork).
**Sprint position**: Documentation-only β-update. Three same-day diagnostic
sprints (F4, F5, F6) followed the F3-maturity β-comprehensive paperwork
and refined the W1e mechanism characterization. This memo captures the
edits applied to paper drafts + CLAUDE.md to **APPEND F4/F5/F6 findings
to the F3-maturity baseline** without superseding it. The F3-maturity
content from yesterday's β-comprehensive sprint stands; the F4/F5/F6
findings sharpen the W1e sub-layer characterization.
**Verdict line: F4-F6-APPENDED-TO-F3-BASELINE.**
**Cross-references**: `debug/sprint_f4_bonding_pk_memo.md`,
`debug/sprint_f5_explicit_core_memo.md`,
`debug/sprint_f6_maxn4_nah_memo.md`,
`debug/sprint_w1c_full_arc_synthesis_memo.md` (F3-maturity baseline that
this β-update appends to).

---

## §0. Executive summary

The F3-maturity β-comprehensive sprint (earlier 2026-05-23) landed the
following baseline narrative in CLAUDE.md §2 and Papers 19/17/34/32:
W1d (cross-block h1 architectural absence, named by F2) CLOSED at FCI
level by F3 with naturals transformation [1.0, 1.0] → [1.9991, 0.0007];
W1e (inner-region overattraction) NEWLY NAMED as F4 target under the
working hypothesis "missing Pauli repulsion against [Ne] frozen core";
F3 §6 named three candidate W1e closure mechanisms.

Three same-day diagnostic sprints tested all three F3-named candidates:

- **F4** (bonding-orbital PK extension, priority 1 by F3 §6): RULED OUT
  at Step 1. Rank-1 PK saturation gives at most 43% closure ceiling.
  W1e is NOT a single-particle Pauli-orthogonality wall.
- **F5** (explicit-core Hartree treatment, F4 §5.2 priority 1): RULED
  OUT at Step 1. Hartree J-K gives 25.7% closure ceiling (predicted).
  W1e is NOT addressable by mean-field core-bonding J-K alone.
- **F6** (basis enlargement to max_n=4, F4 §5.2 priority 2):
  PARTIAL_IMPROVEMENT. Full F3 stack at max_n=4 gives 10.2% PES
  closure (NOT 26.1% the 2-point gate suggested — 2-pt gate
  overestimates by 2.5×). W1e is NOT addressable by basis enlargement
  alone.

**Net refinement**: W1e decomposes into three structurally INDEPENDENT
sub-sub-mechanisms each with ~10-25% closure ceiling, leaving ~65-90%
of the wall as genuinely deep multi-determinant correlation requiring
multi-week+ architectural lift (Schmidt orthogonalization, fully
correlated [Ne] cores) or acceptance of structural-skeleton-scope limits.

This β-update appends these findings as **conservative additions** to the
F3-maturity paperwork — the F3-maturity content stands, the F4-F6 content
refines W1e characterization.

---

## §1. What was applied (the deliverable status table)

### Status per deliverable (per sprint prompt §7 instructions)

| # | Deliverable | Status | Location |
|:---:|:---|:---:|:---|
| 1 | Paper 19 §sec:w1c_residual append new subsubsection "Refinement of W1e mechanism (F4-F6)" | **DONE** | New §sec:w1e_refinement_f4_f6 subsubsection (~180 lines) appended after existing F1-F2-F3 narrative, before §Discussion section. Includes F4-F5-F6 sub-sprint verdicts, refined six-sub-layer W1c hierarchy, methodological 2-pt-vs-PES caveat, source documentation. |
| 2 | Paper 17 §6.10 append paragraph on F4-F6 refinement | **DONE** | New paragraph "Refinement of W1e after F4-F6 (2026-05-23 β-update)" (~40 lines) appended to existing §6.10 narrative. References F4/F5/F6 verdicts + decomposition into three sub-sub-mechanisms + 2-pt-vs-PES caveat. Cross-references to Paper 19 §sec:w1e_refinement_f4_f6 and Paper 34 §sec:conv_w1e_f4_f6_refinement. |
| 3 | Paper 34 §V.D new entry V.D.12 | **DONE** | New §sec:conv_w1e_f4_f6_refinement subsubsection (~140 lines) appended after §sec:conv_w1c_f2_f3 (the F2+F3 entry), before the Cross-pattern observation paragraph. Full five-element template (Observable / Source / Convention / Magnitude / Sprint reference). Three F4/F5/F6 verdicts; three-sub-sub-mechanism table; methodological 2-pt-vs-PES side finding; class iv (architecture-level) classification. |
| 4 | Paper 32 §VIII Sprint M-Z addendum extend with F4-F6 paragraph | **DONE** | New "F4-F6 refinement of W1e characterization (β-update, 2026-05-23)" paragraph (~100 lines) appended to existing Sprint M-Z addendum, before §sec:g3. Refines chemistry-side dissection of W1e in M-Z partition: not cross-shift (F4 rules out single-particle PK), not mean-field-extensible (F5 rules out Hartree J-K), but genuine multi-determinant correlation requiring multi-week+ architectural lift. |
| 5 | CLAUDE.md §2 append F4-F6 to F3-maturity bullet | **DONE** | Existing F3-maturity bullet (line 365) extended in-place with "F4-F6 β-update extension (2026-05-23, three same-day diagnostic sprints appending to the F3-maturity baseline above)" section. F4/F5/F6 verdicts inline + methodological finding + three-sub-sub-mechanism decomposition + natural-pause-point framing + β-update paper edit summary. F3-maturity content unchanged. |
| 6 | CLAUDE.md §3 dead-ends append three rows | **DONE** | Three new rows appended after the existing F2 kernel-shape-substitution row: (a) "Single-particle Pauli orthogonality as W1e closure mechanism" (F4); (b) "Mean-field core-bonding J-K as W1e closure mechanism" (F5); (c) "Basis enlargement to max_n=4 alone as W1e closure mechanism" (F6). Each row includes verdict, mechanism, ceiling, and file references. |
| 7 | Memory file capture | **DONE** | New `memory/sprint_w1c_arc_f4_f6_extension.md` (~110 lines) capturing F4-F5-F6 verdicts + three-sub-sub-mechanism decomposition + methodological 2-pt-vs-PES caveat + paper edit log + forward direction. Sibling to existing `memory/sprint_w1c_full_arc_f3_closure.md` (F3-maturity baseline). |
| 8 | This synthesis memo | **DONE** | `debug/sprint_beta_update_f4f6_memo.md` (this file, ~1800 words). |

### Verification

All four modified papers three-pass clean compile:

| Paper | Pages (was → now) | Compile verdict |
|:---|:---:|:---:|
| Paper 19 | 14 → 16 | 3-pass clean (only pre-existing undefined references) |
| Paper 17 | 14 → 14 | 3-pass clean (only pre-existing + expected cross-paper undefined references) |
| Paper 34 | 115 → 116 | 3-pass clean (only pre-existing + expected cross-paper undefined references) |
| Paper 32 | 54 → 55 | 3-pass clean (only pre-existing + expected cross-paper undefined references) |

Expected cross-paper undefined references (each paper compiles individually; cross-paper refs are unresolvable in individual compile):

- Paper 17 references `sec:w1c_residual` (Paper 19) — pre-existing
- Paper 17 references `sec:w1e_refinement_f4_f6` (Paper 19) — new this β-update
- Paper 17 references `sec:conv_w1e_f4_f6_refinement` (Paper 34) — new this β-update
- Paper 34 references `sec:w1e_refinement_f4_f6` (Paper 19) — new this β-update
- Paper 32 references `sec:w1e_refinement_f4_f6` (Paper 19) — new this β-update

All other warnings are pre-existing (e.g., `sec:matches`, `sec:layer2_d`,
`tab:catalog_off` in Paper 34, `sec:tensor_product_verification` in
Paper 32, `eq:Vpk` in Paper 19) and untouched by this β-update.

---

## §2. Final synthesis verdict

**F4-F6-APPENDED-TO-F3-BASELINE.**

The F3-maturity β-comprehensive sprint paperwork (earlier 2026-05-23)
stands without modification. This β-update sprint adds F4/F5/F6
findings as **conservative refinements** of the W1e mechanism
characterization that F3 named. The F1-F2-F3 narrative is preserved
verbatim in all four papers and in CLAUDE.md §2; the F4/F5/F6 content
appends to it as a "Refinement of W1e mechanism" addendum.

The substantive new content from F4-F6 (beyond the three negatives):

1. **W1e three-sub-sub-mechanism decomposition** (the structurally
   richer characterization replacing F3's working-hypothesis
   "missing Pauli repulsion"): three independent sub-sub-mechanisms
   (basis-truncation ~10% PES, Hartree-pressure ~25% 2-pt,
   deep-correlation-residual ~65-90% of wall), each with measured
   closure ceiling well below the 96% required for 2× experimental
   NaH binding.

2. **Methodological 2-point-vs-PES caveat** (a sprint-cycle-level
   refinement of the diagnostic-before-engineering rule): the 2-point
   gate at fixed reference geometry (R=3.5 vs R=10) used throughout
   F1/F2/F3 overestimates PES-derived closure by ~2.5× on
   monotonically-descending PES in the W1c-residual regime. Future
   sprint-cycle diagnostics should use PES well depth at R_min as the
   load-bearing closure metric, not 2-point differentials. The
   F3-reported D_e^F3 = +4.37 Ha is itself a 2-pt value; the
   corresponding F3 PES well depth at R_min=2.0 is likely larger
   (untested but expected by analogy to F6).

3. **Natural pause point on chemistry arc** (the strategic implication):
   the chemistry-side architectural-extension ladder
   (W1c-cross-screening → W1c-multi-zeta-basis → W1d-cross-block-h1)
   has plausibly hit a natural pause at W1e. The three closure
   candidates F3 §6 named are all empirically tested and all
   structurally insufficient; remaining candidates (Schmidt
   orthogonalization, fully-correlated [Ne] cores) are multi-week+
   architectural investments rather than sprint-scale
   diagnostic-then-implementation cycles. Chemistry arc's forward
   direction shifts from "find the right closure mechanism" to
   "characterize the framework's structural-skeleton scope limits
   empirically for second-row hydride binding," consistent with
   CLAUDE.md §1.7 multi-focal-composition wall taxonomy.

### Three sub-sub-mechanisms remaining open

Per the sprint prompt's "no autonomous claims of W1c closure" constraint,
this β-update explicitly does NOT claim W1c closure. The three W1e
sub-sub-mechanisms remain open:

- **W1e-basis-truncation**: ~10% PES closure at max_n=4; untested at
  max_n=5+; expected diminishing returns based on 6× basis-size
  increase giving only 10% closure.
- **W1e-Hartree-pressure**: ~25% 2-pt closure predicted; untested at
  PES level (F5 stopped at Step 1 prediction); PES-level closure may
  be lower per F6 2-pt-vs-PES caveat.
- **W1e-deep-correlation-residual**: ~65-90% of wall; untested by
  sprint-scale compute; requires Schmidt orthogonalization (~3-4
  weeks, mathematically definitive), fully-correlated [Ne] cores
  (multi-month, requires occupancy-restriction FCI infrastructure
  not present in framework), or acceptance of structural-skeleton-scope
  limits.

### What this β-update does NOT do

- Does not modify §1.5, §1.6, §1.7, §4, §5, §13, §14 of CLAUDE.md (per
  sprint constraints).
- Does not claim W1c closure (three sub-sub-mechanisms remain open).
- Does not introduce fitted/empirical parameters.
- Does not modify the F3-maturity narrative in Papers 19/17/34/32 —
  only appends to it.
- Does not modify production `geovac/` code or tests (no new
  implementation; this β-update is documentation-only).
- Does not run new compute (relies on the three sub-sprint memos
  from earlier 2026-05-23).

### Compute used

This β-update was documentation-only: ~0.5h paper edits + ~0.3h
LaTeX compilation verification + ~0.3h memo writing. Well under
the ~3h compute cap.

---

## §3. Forward direction

Per F6 §5 recommendation (which the PI may revise):

**Priority 1 (default)**: pivot to math.OA arc (~2-4 weeks).
The chemistry arc has produced clean intermediate results suitable
for paper updates (F1-F5 closed sub-walls + F4-F6 characterized W1e
as multi-mechanism deep correlation). The math.OA arc has more open
structural questions to pursue with sprint-scale compute. Papers
38/39/40/42/43/44/45/46/47 are in flight.

**Priority 2 (if PI continues chemistry-arc compute)**: Schmidt
orthogonalization sprint (~3-4 weeks, mathematically definitive
for the basis-orthogonalization sub-mechanism, may close
W1e-deep-correlation-residual partially).

**Priority 3 (cross-system diagnostic)**: pivot away from NaH and
test the chemistry arc on a less-extreme system (e.g., LiH at
max_n=3 where W1c is not in play; cross-system test of the W1e
characterization).

---

## §4. Files modified by this β-update

### Papers (4)
- `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` (16 pages)
- `papers/group2_quantum_chemistry/paper_17_composed_geometries.tex` (14 pages)
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` (116 pages)
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (55 pages)

### CLAUDE.md
- §2 F3-maturity bullet extended with F4-F6 β-update content
- §3 dead-ends three new rows appended after F2 kernel-shape row

### Memory file (1 new)
- `memory/sprint_w1c_arc_f4_f6_extension.md` (sibling to F3-maturity
  memory file)

### This memo (1 new)
- `debug/sprint_beta_update_f4f6_memo.md` (this file)

### NOT modified
- Production `geovac/` modules (no new implementation)
- Tests (no new implementation; baseline regression preserved)
- F3-maturity paper content (only appended, not modified)
- F3-maturity memory file `memory/sprint_w1c_full_arc_f3_closure.md`
  (kept as the F3-maturity baseline; F4-F6 memory file is a sibling
  extension)
- CLAUDE.md §1.5, §1.6, §1.7, §4, §5, §13, §14 (per sprint constraints)

---

**End of β-update memo. Verdict: F4-F6-APPENDED-TO-F3-BASELINE.**
