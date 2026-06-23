# Sprint memo — `/qa group1` Batch 3 CODE-DIMENSION validation pass (2026-06-22, v4.43.5)

**Target:** Group 1 Batch 3 = Papers 29, 39, 40, 50, 52 (Ihara/Ramanujan-Hopf; tensor-product
propinquity; unified 4/π rate; CFT₃ wedge-KMS/F-theorem; Category III correspondence). **CODE
DIMENSION ONLY** (C1/C2). Per the all-dimensions rule this is a **validation pass, NOT a cert**
(1 of 5 gating dimensions ⇒ INCONCLUSIVE-for-cert). PI goal: surface NEW flavors of code-defect
to fix corpus-wide (mirroring the Batch-2 false-closure + closed-form sweeps).

## 1. Method
- Base HEAD 4650573 (v4.43.4). Worktree `../geovac-qa-seed-g1b3-code`, removed + leak-scanned
  (real corpus clean on all 4 seed sites). Seed key `debug/qa/group1_batch3_code_seed_key.json`.
- 5 code-reviewers (29/39/40/50/52), opus, blind, path-pinned, RUN the tests.
- 4 fresh vacuous-tolerance seeds (one per backing test for 29/39/40/50); Paper 52 has no test
  (code-52 = coverage-gap assessment).

## 2. Calibration — CALIBRATED (sensitivity 4/4)
| Seed | Site | Catcher | Grade |
|---|---|---|---|
| S-code-29 | `test_paper29_bound_crossing.py:131` (`> -1e9`) | code-29 | NIT (backstopped by l.132 + 3 sign-flip siblings) |
| S-code-39 | `test_gh_convergence_tensor.py:619` (`< 1e9`) | code-39 | NIT (backstopped by l.615/617) |
| S-code-40 | `test_paper40_universal_rate.py:75` (Haar `< 1e9`) | code-40 | MATERIAL (no backstop) |
| S-code-50 | `test_paper50_wedge_kms.py:138` (`< 1e9`) | code-50 | MATERIAL-NIT (backstopped by l.139) |
Specificity clean — M1–M4 controls respected (the v4.33/4.34 fixes confirmed holding). code-52
correctly assessed Paper 52 as purely conceptual (rides on Paper 38/42/50; no dedicated test needed).

## 3. NEW FLAVORS surfaced
### Flavor A — code-side withdrawn-claim zombies (SWEPT corpus-wide) — the headline
A withdrawn/false claim living in a **production-code function's returned STRING or DOCSTRING** (not
in a paper, not in a test assertion), invisible to the test suite — the code-side analog of the
Batch-2 paper zombies. code-39 Finding B/C found `gh_tensor_theorem_statement()` returning the
WITHDRAWN Pythagorean ("C₃ < 1 → 1", "1 + o(1)") **and** the Paper-38-result-as-"Latrémolière
propinquity" mislabel (should be van Suijlekom state-space GH). Swept across the contained
`gh_convergence{,_tensor}.py` pair (5 sites):
- `gh_convergence_tensor.py`: `gh_tensor_theorem_statement` (propinquity→state-space GH + Pythagorean→triangle); `tensor_L3_lipschitz` + `tensor_L5_assembly` docstrings (withdrawn "< 1"→triangle ≥1→√2).
- `gh_convergence.py`: `gh_theorem_statement` + `LimitIdentification.statement` + its docstring (propinquity→state-space GH; single-factor C₃=1 is correct, only the mislabel fixed).
- **Sub-flavor — zombie-locking tests:** `test_gh_convergence_tensor.py:487` and `test_gh_convergence.py:456` *asserted the presence* of "Latremoliere"/"propinquity" in the statement strings — pinning the zombie. Updated both to assert the corrected "state-space"/"van Suijlekom" keywords + a guard that "Latremoliere"/"propinquity" do NOT appear (so the zombie cannot re-enter). Strings/docstrings only — no logic change; the numeric functions (`c3_full_pythagorean_bound`, etc.) were already correctly annotated as withdrawn (v4.33.0) and return the triangle values.

### Flavor B — load-bearing test backing on a prunable `debug/` path with skip-on-absence (FLAGGED)
code-40: the rank-2 universality backing (the SU(3)/Sp(2)/G₂ panels + the 5641-cell C₃=1 verifier)
lives under `debug/qa/_resurrected/` (`sp2_g2_rate_constant.py`, `dirac_triangle_extended_verify.py`,
`l2_universal_rate_verify.py`) and is imported into `test_paper40_universal_rate.py` under a
`skipif(not _HAVE_RES)` guard. If `debug/` is pruned (§9 transient policy), the skip silently
disables the rank-2 tests and the suite stays green with only the SU(2) anchor. **Compounds with the
genuine coverage gap:** the general-G C₃(G)=1 claim (PRV summand / Brauer–Klimyk, "5641/5641 cells")
is verified only by `dirac_triangle_extended_verify.py::run_panel`, which `test_paper40` imports for
its constructors but **never calls** — so general-G C₃=1 has NO `tests/` assertion. (test_r25_l3_lipschitz_bound.py
backs C₃ only for the S³/SU(2) single factor, not general G.)

### Flavor C — precision overstatement vs test gate (DISSOLVED on reconcile)
code-52: p52 l.408 "matches KPS at $61+$ digits" while the shipped test gates at 1e-40. Reconcile:
code-50 *measured* the agreement at ~1e-82, so "61+ digits" is a TRUE conservative claim; the test
gate (1e-40) is merely loose, not wrong. **Not a defect** — a true measured claim is not understated
to match a conservative regression gate. No fix.

## 4. Remediation applied (Flavor A)
`geovac/gh_convergence_tensor.py` (3 sites), `geovac/gh_convergence.py` (3 sites),
`tests/test_gh_convergence_tensor.py` (zombie-lock assertion → corrected + guard),
`tests/test_gh_convergence.py` (zombie-lock keyword list → corrected + guard). No production logic
changed (strings/docstrings + test assertions only).

## 5. Flagged follow-ons (named, not done this pass)
- **Flavor B + p40 general-G C₃=1 coverage gap (MATERIAL):** add a `tests/` assertion that calls the
  panel verifier and pins C₃(G)=1 for SU(3)/Sp(2)/G₂ (a reduced panel, not the full 5641 cells), AND
  migrate the resurrected rank-2 drivers to a permanent home (or convert the `skipif` to a hard error
  when absent) so the headline rank-2 universality cannot silently de-test. → dedicated backfill sprint.
- Minor NITs: code-29 `test_galois_ihara` self-reference (Galois of hardcoded polys, not live-zeta
  factors; cross-file link via test_ihara_zeta) + Obs 2 missing inline tier; code-50 Thm 5.3
  (S⁵ log-3 cancellation) only indirectly tested. Low priority.

## 6. Honest scope
- **Validation pass, NOT a cert** (code dimension only). Code dimension is calibrated (4/4) and, beyond
  the seeds, **MATERIAL-clean except the p40 general-G C₃=1 gap** (flagged) — the new Flavor A was
  swept corpus-wide. A full Batch-3 cert (all 5 dimensions) remains the certification step.
- Flavor A is the transferable lesson: **withdrawn/descoped claims survive in production-code
  strings/docstrings AND in tests that assert their presence** — both layers must be swept, and the
  corpus-wide grep (`geovac/` for the withdrawn framing + tests/ for keyword-presence asserts) is the
  convergent fix, exactly as the Batch-2 paper sweeps were.
