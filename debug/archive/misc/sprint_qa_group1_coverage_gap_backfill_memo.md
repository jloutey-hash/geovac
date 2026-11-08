# Sprint memo — group1 coverage-gap test-backfill

**Date:** 2026-06-21
**Version:** v4.34.0
**Trigger:** PI — "fix the coverage gaps we didn't remediate this last round" (the four verified-true, weak/no-test gaps raised in the v4.33.0 Batch-3 re-run).
**Verdict:** all 4 gaps closed with genuine recompute-from-framework tests.

Discipline: [[feedback_validate_before_reducing]] (build genuine validation, no replacement-tautology theater) + [[feedback_resurrect_pruned_artifacts]] (recover pruned backing from git before concluding untested). Both pruned memos (`hopf_u1_block_test_memo.md`, `bh_phase0_entanglement_entropy_memo.md`) were resurrected from `56d29ddb^` (v4.0.0 hygiene refactor) for context before building.

## Gap 1 — p29 Observation 2 (m_ℓ→−m_ℓ block decomposition)

**Finding:** NOT actually untested. `tests/test_hopf_u1_block.py` exists and passes 12/12 (6 fast + 6 slow), genuinely proving — by restricting the LIVE node-level Ihara–Bass matrix M(s)=I−sA+s²Q to the ℤ₂ reflection blocks — that det(M_sym)=(s−1)(s+1)·P₂₂, det(M_antisym)=P₁₂, and det(M_sym)·det(M_antisym)=det(M) symbolically (+ the Dirac Rule-B P₂₂·P₂₄ analog). The rr3 code-reviewer (searching `test_paper29_*`) missed it because (a) it's named `test_hopf_u1_block`, (b) it imported its provider from `debug/archive/misc/` (transient).

**Fix:** `git mv debug/archive/misc/hopf_u1_block_test.py → tests/_hopf_u1_block_data.py` (underscore → not collected as a test; pytest prepend-mode makes it a bare-name sibling import); import in `test_hopf_u1_block.py` rewired `from debug.archive.misc.hopf_u1_block_test` → `from _hopf_u1_block_data`. Cited inline in Obs 2 (replacing the `debug/` memo cite — also neutralizes a dangling debug ref). 6 fast + (previously-confirmed) 12/12 pass.

## Gap 2 — p29 Alon–Boppana bound-crossing

**Finding:** backed only by `test_alon_boppana_sweep.py`, a schema/arithmetic smoke test on a frozen JSON that never recomputes spectra (a regression in any graph constructor leaves it green).

**Fix:** new `tests/test_paper29_bound_crossing.py` — recomputes the sub-Ramanujan deviation max|μ_nontrivial|−√q_max from the Hashimoto spectra of the live constructors (`GeometricLattice`, `build_bargmann_graph`, `build_dirac_s3_graph`). Verified to bit-match the JSON (S5 N=3 −0.3568, Dirac-B n=3 −0.1188, Dirac-B n=4 +1.5301). Asserts the finite-size **sign flip** for three families (S3-Coulomb, S5-Bargmann, Dirac-B): Ramanujan at small size → crosses (dev>0) at V≈30–60. 3 fast + 4 slow = 7/7.

## Gaps 3+4 — p50 wedge-KMS entropy (Sec 3) + CHM identification (Thm 4.2)

Both ride one object: the BW modular Hamiltonian K_α = J_polar on the truncated Camporesi–Higuchi triple, built by `for_bisognano_wichmann(n_max)` → `restrict_K_alpha_to_wedge()`.

**New `tests/test_paper50_wedge_kms.py` (24/24):**
- **Gap 4 (CHM):** wedge K_α distinct spectrum = {1,3,…,2n_max−1} (odd integers 2m_j) with multiplicity (n_max−k)(n_max−k+1); 2π modular-flow closure residual <1e-10 (integer spectrum → exact). Verified from the live machinery at n=2..5.
- **Gap 3 (entropy):** S(ρ_W) from the live wedge spectrum reproduces the BH-Phase0 memo values and equals the thermodynamic closed form log Z + ⟨K_α⟩ to <1e-13. The analytical shell-degeneracy entropy is validated bit-exactly against the machinery (n=2..5), then used for the large-n asymptotic (the full-Dirac build is O(dim²): n=50 → 116 GiB, so the machinery can't reach it): S − 2 log n_max → coth1−log(2sinh1)=0.45844874… (paper rounds to 0.458448), windowed slope → 2, area-law fit S∼n_max² rejected (R²<0.9) vs log fit R²>0.999. s-deformation limits (s→0 → log dim_W; s→∞ → log n_eq) confirmed.

## Verification

- C13 group1 PASS — all four inline test citations (`test_hopf_u1_block`, `test_paper29_bound_crossing`, `test_paper50_wedge_kms`) resolve to live tests.
- p29, p50 compile errors=0/undef=0.
- Matrix rows updated NO-TEST/BACKED-WEAK → BACKED-SOUND.
- Full combined verification (all 3 gap-test files --slow) — see `debug/gap_backfill_verify.log`.

## Honest scope

- Gap 1: theorem-grade (12/12, symbolic bit-exact block-det identity on live graphs).
- Gap 2: numerical/structural — recomputed deviations bit-match; the sign-flip is the finite-size content (asymptotic failure of the dense families is the claim, not a continuum theorem).
- Gap 4: theorem-grade (exact odd-integer spectrum + multiplicities + 2π closure from the live BW construction).
- Gap 3: the entropy values + thermo identity are exact from the machinery; the asymptotic constant/slope use the degeneracy formula (validated vs machinery at small n) because the full-Dirac build cannot reach large n. This is genuine (the formula IS the framework's shell structure), not a hardcode.

No paper *claims* changed — these are backing-only additions; the prose tiers were already honest.
