# Sprint: group4 A/B ERI-framing rule + corpus compliance — canonical memo

**Date:** 2026-06-28 · **Version:** v4.53.0 · **PI direction:** during the `/qa group4` step-1
freeze checkpoint, the PI chose **CF-1 option A (disclose, keep the pair-diagonal A-state)**,
directed "make very clear we left B on the table," asked me to **study the repo** to confirm
the A/B read, and then to **codify a QA rule around the framing** ("so our review process will
find framing zombies") + **run a sprint cycle** to bring the corpus into compliance. Predecessor:
`debug/sprint_group4_prework_memo.md` (the v4.52.0 pre-work that surfaced CF-1).

## 1. The decision + the repo study

**CF-1 = A (disclose).** The codebase intentionally carries two ERI angular selection rules:
**A — pair-diagonal** (`composed_qubit._ck_coefficient` & `lattice_index._ck_coefficient`,
`q=mc−ma` ⇒ nonzero only when `m_a=m_c ∧ m_b=m_d`) = the **sparsity/QC-product** rule;
**B — global-M_L** (`casimir_ci._gaunt_ck`, `q=ma−mc`, the exact `m_a+m_b=m_c+m_d`) = the
**physics-accuracy** rule. A is for sparsity, B is for accuracy — and that split is correct.

**Repo study (the PI's "does the product apply A, B, or both?"):** the QC product is **uniformly
A**. Empirically confirmed: the atomic `lattice_index` ERI table has **0 of 65** m-transfer
entries (the He n_max=2 m-swap ERI ⟨p₊₁p₋₁|p₋₁p₊₁⟩ is dropped); the composed `_ck_coefficient`
drops the same. Global-M_L (`_gaunt_ck`) is imported only by `casimir_ci`/`cusp_correction`/
`dirac_ci`/`internal_multifocal` — the precision-physics paths, **never** the shipped QC
Hamiltonians (`ecosystem_export.hamiltonian()` → composed + atomic, both A). The composed
builder imports `casimir_ci.get_rk4` for the exact **radial** R^k but uses its own pair-diagonal
**angular** rule. → confirms the "A state," uniformly.

## 2. The QA rule (the PI's primary ask)

- **`docs/qa/criteria.md` → new shared section "Dual-rule ERI framing (A = sparsity, B =
  accuracy)"** + change-log entry. Names the two rules and their purposes; defines a **framing
  zombie** (1: a pair-diagonal sparsity number presented as the *exact*/full Gaunt-selection-rule
  value without disclosure; 2: a physics-accuracy claim silently resting on A; 3: the wrong
  density/count for the context). Enforcement = `claims-reviewer` **enumerates every** sparsity/
  density/Pauli claim and checks its rule + the C16 backstop. Sharpens C3/C8/§1.5; shared
  (group4 primary, group3 P22/31 related).
- **`debug/qa/check_retracted_terms.py` → new C16 entry `pair-diagonal-as-exact-sparsity`**
  (scope group4, fail-severity). Patterns: the 2.7×-vs-STO-3G market test, the 9.23/d-block-
  sparser claim. Exempt only with a `pair-diagonal`/`approximation`/`global-M_L`/`parity` flag
  within ±5 lines. **Verified it FIRED on the 2 live zombies (P14:60, P20:997) before the fix,
  and is clean after.**

## 3. Corpus brought into compliance (option A applied)

- **Paper 14** (26pp, 3-pass clean): abstract d-block claim qualified (pair-diagonal; reverses
  under B); new **`\label{sec:eri_rule}`** subsection "The pair-diagonal ERI approximation:
  sparsity vs. accuracy" (the two rules, the constant 2.51×/3.25× re-pricing, STO-3G→parity, the
  d-block reversal, **B left on the table**, pinned by the new test); the stale "additional
  m-selection rules (m₁+m₂=m₃+m₄)" density paragraph corrected (production is pair-diagonal).
- **Paper 20** (11pp, **zero-undefined** via bibtex): new matching `sec:eri_rule` subsection; the
  2.7×-market-test line + the TM-table caption (9.23) qualified pair-diagonal; **all 5
  pre-existing undefined refs fixed** (4 cross-doc `\ref`s to Paper 14 → "Paper 14"; `Childs2021`
  resolves on clean bibtex).
- **`tests/test_paper14_eri_rule.py`** (NEW, 3 tests, pass): product realizes A (m-swap dropped,
  atomic + composed); global rule keeps it; B re-prices 2.51× (main) / 3.25× (d-block), d-block
  harder (the reversal). **Locks the A-state** — drift to B fails the test.
- **`geovac/nuclear/potential_sparsity.py::angular_zero_count` docstring** corrected (was
  mislabeled global-M_L; computes pair-diagonal D_pd; the `m_a+m_b=m_c+m_d` guard is redundant).
- **Synthesis** already CF-1-honest (drafted to A in v4.52.0; C16-clean).

## 4. Second-look sweep

- **group3 (P22/31): already compliant** — Paper 22 discloses "a stricter pair-diagonal D_pd
  1.44%" vs universal "6.06%" (done in group3 cert). No edit.
- **group2: one cross-branch flag (NOT auto-fixed — certified branch).** Paper 13 + fci_atoms
  say "exact rational Slater integrals" for the atomic CI; the atomic path is pair-diagonal
  (radial R^k exact, angular = A), so "exact integrals" doesn't disclose the pair-diagonal
  angular subset. Energy-negligible (~mHa). For a future group2 re-touch — PI's call.

## 5. Verification

- group4 deterministic gates **PASS** (C11/C13[new test ref resolves]/C14/C15/C16).
- `test_paper14_eri_rule.py` 3/3; broad regression **80 passed** (incl. topo S³ proofs +
  density/scaling consumers, 9m46s) + the prior lean 173.
- Papers compile 3-pass clean (P14 26pp, P20 11pp zero-undefined, synthesis 4pp).
- Only production code edit = the `angular_zero_count` **docstring** (no logic).

## 6. Honest scope

- **Theorem grade:** none (no new theorems).
- **Empirically established (this sprint):** the QC product is uniformly pair-diagonal (A) —
  atomic ERI table has 0 m-transfer entries; composed drops the m-swap; global-M_L is
  precision-physics-only. This is a *characterization* of the shipped code, now regression-pinned.
- **Policy/process:** the A/B dual-rule framing is now a codified QA rule (criteria.md +
  C16) — a guideline, not a physics result.
- **Disclosure (not new physics):** Papers 14/20 now disclose option A + that B was left on
  the table; the d-block-sparser claim is corrected to pair-diagonal-specific.
- **Named open follow-ons:** (1) the `/qa group4` cert RUN itself (seeded panel) — deferred to
  a fresh invocation against this committed corpus; DoD updated (CF-1=A applied), ready to
  freeze. (2) the **group2 cross-branch "exact integrals" flag** — disclose the pair-diagonal
  angular subset at a future group2 re-touch (PI's call; certified branch). (3) the 7
  remaining non-headline NO-TEST gaps from v4.52.0 (ERI 1/M², 13× QWC, balanced-878, etc.).
- **Not done:** no `/qa group4` panel ran; the group2 flag is flagged-not-fixed.
