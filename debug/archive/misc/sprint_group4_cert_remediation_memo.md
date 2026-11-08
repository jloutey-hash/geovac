# Sprint: /qa group4 first cert (FAIL) + full remediation — canonical memo

**Date:** 2026-06-28 · **Version:** v4.54.0 · **PI direction:** `/qa group4` (the cert gate);
then "start remediation"; then "send out an explorer agent" for the placeholder cites.
Sub-records: `docs/qa/group4.done.md` change-log (per-finding), `docs/claim_test_matrix.md`
group4 remediation note, `debug/qa/group4_seed_key.json` (answer key, gitignored).
Predecessors: v4.52.0 (pre-work), v4.53.0 (A/B framing rule + disclosure).

## 1. The cert run

First cert of group4 (Papers 14/16/20/23 + synthesis), whole-group, against the frozen DoD
(`docs/qa/group4.done.md` @ commit 923009a). Built a detached seeded worktree, planted **8
calibration seeds** (4 code S2/S3, 1 framing-zombie, 1 status-overstatement, 1 wrong-arXiv-ID,
1 synthesis-zombie) covering every gating dimension, dispatched a **fresh path-pinned panel**
(4 code-reviewers 1/paper, 2 claims chunks {14,20}/{16,23}, 1 citation, 1 synthesis; opus).

**Verdict: FAIL — fully calibrated (trustworthy).** Sensitivity **8/8** (every plant caught
by its expected reviewer), specificity **5/5** (zero false positives; the disclosed
`sec:eri_rule` correctly NOT flagged → the new framing dimension is calibrated). Worktree
removed, leak-scan clean (the McArdle wrong-ID seed was worktree-only; the real `.bib` keeps
the correct `1808.10402`).

## 2. The 9 genuine MATERIAL findings (seeds excluded) — all remediated

| # | Finding (dimension) | Fix |
|---|---|---|
| 1 | **Framing zombies** (~10 undisclosed d-block-"cheaper/sparser/economical" + market-test loci in P14/P20 that DODGED the C16 grep) — claims | Disclosed each (pair-diagonal qualifier + the reversal); **broadened the C16 `pair-diagonal-as-exact-sparsity` pattern** to catch "cheaper/economical/structurally-sparser-than-s/p/2.7×-Pauli" |
| 2 | **μ_free code bug** — `atomic_classifier.py` computed `2ν²`, not the headline SO(3N) Casimir `ν(ν+3N−2)/2 = 2(N−2)(N−1)` (N=10: code 128 vs headline 144) — code/C2 | Fixed code (2 loci); re-pinned `test_mu_free` (+row2 +3 dataclass loci) + 3 P16 secondary loci (Ne 128→144, Dirac 36450→36720, per-pair eq). **Diagnostic-only — μ_free is NOT used in any Hamiltonian/solver** (the solvers already use the correct formula); no result changed |
| 3 | **`test_balanced_row2::TestResourceTable` RED** (12 fails under `--slow`; the pre-work "all green" never ran it) — code/C1 | spec factories `composed_qubit`→`molecular_spec` + strip the R_SH/R_PH/R_SiH aliases the factory rejects. 12 fails → **39 pass** |
| 4 | **trenev2025 misattribution** — a *vibrational-spectra* paper credited with the electronic LiH/H₂O Gaussian Pauli counts/exponents (self-contradicted by its own bibitem note) — citation/C4 | 15 P14 loci + P20 caption → "GeoVac's OpenFermion recompute; Trenev = methodology + the Q^3.9–4.3 range" |
| 5 | **"~1/M²" ERI-density overstatement** — the paper's own table shows ~M^{-0.85} (≈1/M) — claims/C3 | 9 P14 loci → "~1/M (≈M^{-0.85})" |
| 6 | **P23 "fine-structure constant conjecture"** — Paper 2's K=π(B+F−Δ) is an **Observation** — C5 hard-prohibition | → "observation" |
| 7 | **Paper-38 "Latrémolière propinquity" zombie** — Paper 38 = van Suijlekom **state-space GH** — C7/C16 | Fixed P14/P20/P23/synthesis/bibitem + the P20 bibitem title |
| 8 | **P20 conclusion "38 molecules"** — stale headline (37) | → 37 |
| 9 | **Placeholder/unverifiable cites** (PI-directed explorer search) — citation/C4 | `Sunaga2025` **misattributed** → re-keyed **Swain et al. arXiv:2211.06907** (RaH 18q/47,099 exact); `caesura2025`→**PRX Quantum 6, 030337 / arXiv:2501.06165** (+"real-space BLISS"→"BLISS-THC"); `ChildsBerry`→**arXiv:1501.01715**; `MartinezYRomero2004`→**physics/0402061**; `BJL` unverifiable → **withdrawn** |

Plus a pre-existing P23 broken `\ref{ssec:cross-omega-magn}`→`ssec:omega-magn` (found incidentally).

## 3. Verification

- Deterministic gates C11/C13/C14/C15/C16 on group4: **PASS** (re-run after every batch).
- Affected tests green: `test_atomic_classifier` **197**, `test_atomic_classifier_row2` (incl.
  the corrected μ_free), `test_balanced_row2 --slow` **39**, `test_paper14_eri_rule` **3**;
  broad regression (atomic_classifier consumers + builders + topo S³ proofs) **264 passed**.
- All 5 papers compile three-pass clean (P14 26pp, P16 7pp, P20 11pp, P23 11pp, synthesis 4pp;
  **zero undefined** refs/cites). P20 via bibtex.

## 4. Honest scope

- **Theorem grade:** none (no new theorems).
- **Process verdict (trustworthy):** the cert was fully calibrated (8/8, 5/5), so the FAIL is
  trustworthy and the remediation is against verified defects.
- **Genuine bugs fixed (not cosmetic):** the μ_free code formula (diagnostic-only, no
  Hamiltonian affected — verified) and the `test_balanced_row2` RED (a real test-rot the
  pre-work's `--slow`-skipped "all green" missed).
- **Corrections (evidence-based):** all framing/citation/status fixes trace to primary text,
  the paper's own data/bibitems, or a verified external publication (the explorer search).
- **Process lesson:** even after heavy pre-work (v4.52/4.53), the FIRST cert of a fresh branch
  surfaced 9 genuine material defects — the qa.md "batch the first cert" rationale holds. The
  C16 grep was too narrow (the reviewers caught ~10 C16-dodging framing phrasings); broadened
  per the maintenance rule.
- **Named open follow-ons:** re-run `/qa group4` (the certified-PASS confirmation run) against
  the committed remediated corpus — all 9 findings are fixed, so it should PASS. The
  ChildsBerry/MartinezYRomero cite identifications are medium-confidence (explorer best-fit);
  a re-cert citation-reviewer will re-verify.
- **Not done:** the certified PASS itself (the re-cert is the next `/qa group4` invocation).
