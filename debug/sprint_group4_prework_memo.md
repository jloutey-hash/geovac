# Sprint: group4 QA pre-work + backfill/inconsistencies — canonical memo

**Date:** 2026-06-28 · **Version:** v4.52.0 · **PI direction:** "let's look at group4 now…
do some pre-work to get ready for a qa pass" → CF-1 "quantify the full library first" + all
four pre-work items → "address the backfill/inconsistencies first" → library count "37 (ship
as-is)". Sub-records: `debug/qa/group4_cf1_library_sweep_memo.md` (CF-1 detail),
`docs/qa/group4.carryforward.md` (CF-1 carryforward), `docs/claim_test_matrix.md` (group4 rows).

## 1. What this sprint did

Prepared group4 (Paper 14 qubit encoding, 16 periodicity, 20 resource benchmarks, 23 nuclear
shell) for a `/qa group4` cert run: built the missing scaffolding (DoD, synthesis, claim-test
matrix), quantified the gating CF-1 issue, then — on PI direction — closed the backfill
coverage gaps and resolved the cross-corpus inconsistencies. **No cert run was performed**
(that is the PI-invoked next step); this is the readiness work that precedes it.

## 2. Pre-work deliverables (4/4 + deterministic pre-validation)

1. **CF-1 library re-pricing sweep** (`debug/qa/group4_cf1_library_sweep.py`, build-only, 35
   composed molecules under both ERI rules). The pair-diagonal→global-M_L re-pricing is a
   **constant factor within valence class**: **2.51× main-group, 3.25× d-block**. Two material
   consequences: (a) the O(Q^2.5)/linear scaling is ROBUST (constant prefactor ⇒ unchanged
   slope) — only the *magnitude* re-prices; the LiH-vs-STO-3G market test goes to parity; (b)
   the "d-orbitals are sparser (9.23<11.10)" claim REVERSES under the physical rule (d-block
   30.0·Q > main-group 27.9·Q) — the low pair-diagonal coefficient is an artifact of dropping
   more m-swaps, not sparser physics.
2. **DoD** `docs/qa/group4.done.md` (drafted for PI freeze; branch-defining criterion =
   QC-resource-claim honesty + CF-1 disposition).
3. **Synthesis** `papers/synthesis/group4_quantum_computing_synthesis.tex` (NEW, 4pp, three-pass
   clean; CF-1-honest spine 14→20→16→23).
4. **claim_test_matrix** group4 rows populated (was zero) by a dispatched agent; all cited
   tests RUN GREEN (254+105+78).
5. **Deterministic layer pre-validated GREEN** (C11/C13/C14/C16) incl. the new synthesis;
   C11 caught + fixed two wrong bibitem titles in the draft.

## 3. Backfill — 7 of 14 NO-TEST gaps closed

The matrix surfaced that **every keystone scaling exponent was unpinned** (the only prior test
fed `fit_pauli_scaling` synthetic data). Closed:

- **`tests/test_paper14_scaling.py`** (NEW, 4 slow tests) — exponents now fit from live He
  (n_max=2–4) / composed-LiH (max_n=1–3) data: atomic N_Pauli **Q^3.10** (paper 3.15), 1-norm
  **Q^1.67** sub-quadratic (paper 1.69), QWC **Q^3.356** (paper 3.36), composed **~Q^2.5**.
  Resolved a potential-material puzzle: composed O(Q^2.5) (within-molecule, vary max_n) and the
  exact 11.10·Q linearity (across-molecule, fixed max_n) are two different sweeps — both true.
- **`tests/test_paper23_resource_counts.py`** (NEW, 4 tests) — deuteron (592 / 80 Z-only / 512
  XY / 1-norm 342 MeV) and He-4 (712 / 1-norm 467 no-Coul / 462 with-Coul) pinned exactly.

7 gaps remain OPEN, all non-headline: ERI 1/M² (decays slower — BACKED-WEAK), double-
factorization rank, 13× QWC factor, balanced-LiH 878/0.20%, per-pair-4, Z=137, Fock rigidity,
He-4 1.20×/12.25× ratios.

## 4. Inconsistencies — all four resolved

1. **Library count** (28/30/37/38/40) — **DECIDED 37 (ship as-is, PI), applied corpus-wide.**
   `_SYSTEM_REGISTRY` = 37 systems (35 composed molecules + He + H2). The 3 organics
   CH₂O/C₂H₂/C₂H₆ Paper 20 listed are **not buildable** (bond-length entries only) → removed
   from Paper 20's abstract + multi-center table (8→5). Synced 37 across Paper 14, Paper 20,
   the synthesis, ecosystem docstring (28→37), test comment (40→37), claims_register, and
   CLAUDE.md §2 + §1.1 + §1.5 (the §1.1/§1.5 count edits applied under explicit PI "37"
   direction — normally PM-restricted; flagged).
2. **Paper 16 structure types** — was internally inconsistent (abstract/§IV "4", Table/
   conclusion/code "5"). FIXED → **5 types (A/B/C/D/E)** consistently + new Type E subsection;
   the code's E/F (p/d-block) split noted as an implementation refinement.
3. **Paper 23 1-norms** — both stale (code stable since v2.7.0, so drafting-era staleness, NOT
   a regression; term counts always correct): deuteron 227→**342 MeV**, He-4 557/552→
   **467/462 MeV**. FIXED + pinned by the new test.
4. **LiH 333 vs 334** — RESOLVED, not a defect: **334 is correct** (the shipping `hamiltonian()`
   API value = paper + CLAUDE.md); 333 is the raw `build_composed_hamiltonian` count (excludes
   the identity term).

## 5. Verification

- Deterministic gates (C11/C13/C14/C16) on group4: **PASS** (re-run after every edit batch).
- Lean regression: **173 passed**, 2 skipped (incl. the 18 topo S³ proofs + new tests).
- Full library `test_ecosystem_export --slow` + `test_paper14_scaling --slow`: **111 passed**
  (3m48s) — confirms all 37 systems build + the scaling backfills hold.
- All 5 group4 papers + the synthesis compile three-pass clean. Only code edit was the
  `ecosystem_export` docstring (no logic), covered by the passing consumer test.

## 6. Honest scope

- **Theorem grade:** none claimed this sprint (no new theorems).
- **Pinned-from-data (MEASURED, now backed):** the four Paper 14 scaling exponents (3.10/1.67/
  3.356/~2.5, fit over the tractable n_max range; paper's full-range values are 3.15/1.69/3.36/
  2.5) and the deuteron/He-4 resource counts + 1-norms. These are regression-pinned, not proofs.
- **Numerical observation:** the CF-1 re-pricing is a constant factor per valence class
  (2.51×/3.25×); the "d-block sparser" reversal is a build-only measurement.
- **Corrections (evidence-based, not new physics):** P16 5-type, P23 1-norms, library=37,
  organics dropped. All trace to live code / the registry.
- **Named open follow-ons:** (1) **CF-1 disposition A (disclose) vs B (switch)** — the
  load-bearing cert FREEZE decision, still the PI's; the synthesis is drafted to A. (2) The 7
  remaining non-headline NO-TEST gaps. (3) **Pre-existing (NOT this sprint):** Paper 20 has 4
  broken section-`\ref`s (sec:spinor_composed/composed, subsec:spinor_scope) + a missing
  `Childs2021` bibitem — surfaced incidentally during compile; for the cert run.
- **Not done:** the `/qa group4` cert run itself (PI-invoked); the DoD is drafted, not frozen.
