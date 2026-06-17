# `/qa group3` bite 2 — FAIL: findings + remediation checklist

**Date:** 2026-06-16. **Run:** `/qa group3` bite 2 (first run).
**Scope:** Papers 18, 54, 55, 56, 57 + full group3 synthesis (+ Paper 24 re-touch).
**Verdict: FAIL.** Calibrated code/claims/citation/deterministic dimensions surfaced numerous verified material defects across every bite-2 paper. The first bite (22/24/31, certified v4.18.0) was already hardened; this periods/Tannakian cluster was not, and the gate caught a lot.

> Seed key: `debug/qa/group3bite2_seed_key.json`. Worktree removed; no seed leaked to the live corpus (verified).

## Calibration scorecard
- **Sensitivity 4/5** — caught S6 (claims), S1 (citation), S2 (code), S7 (deterministic); **MISSED S8 (synthesis)**.
- **Specificity 5/5** — no known-good control flagged material.
- **Synthesis dimension NOT calibrated** — the reviewer saw S8, then re-read the *real-corpus path instead of the worktree*, got unseeded text, and reversed its own correct catch. The dimension must be **re-run** (pin reviewers to the worktree path); it cannot rescue the FAIL because four other dimensions already establish it.

## Verification status legend
- **[V]** = PM-verified against primary text this run.
- **[R]** = reviewer-reported (quoted lines / ran tests); **verify during remediation** before editing.

---

## Remediation pass 1 (2026-06-16) — status

**DONE (verified vs primary text, all edited `.tex` compile exit 0; group3 C11/C12/C13 + K-label PASS):**
- QUICK WINS — all: Paper 18 §13.5 (two `conjecture` envs → `observation` + label/refs + 2 prose lines); Paper 38 propinquity→state-space-GH sweep (Papers 54/55/57, incl. 2 bibitems retitled to the exact `\title`); citations (Paper 55 `hoffman` title + `glanois` venue→160/334–384/2016; Paper 57 `boyle` title→match-the-ID); Paper 57 C6 (discrete-vs-continuum verb); synthesis "theorem-bound boundary" softened; companion files (`forced_free_seam.md` + `principle_hunt_audit.py`) K "conjectural"→Observation.
- SUBSTANTIVE (PM-fixable): Paper 18 κ §VII.F reframe ("Derivation of κ"→"two routes", dropped "more fundamental", coincidence/Observation framing); Paper 18 zeta_d2 cite (dropped false "80-digit PSLQ", repointed to `test_qed_two_loop.py`); Paper 57 "128"→matter-sector caveat + full-axiom-260 (matches Paper 32's v4.16.1 correction).

**RE-CATEGORIZED to determinations** (verification showed these are NOT mechanical fixes; held off per the PI principle "reason out the validation before reducing the claim"):
- **Paper 40 title** — its abstract *explicitly* claims to UPGRADE state-space GH to Latrémolière propinquity via the same five-lemma chain Paper 38 was descoped on. Determine whether the propinquity claim genuinely holds (→ title stays, Paper 38 inconsistent) or descope to state-space GH like Paper 38 (group1, 78 occurrences).
- **Paper 56 §hgv dim** (9 vs 15) — §hgv uses `0≤l<n` (hydrogenic, physically correct → 9 gen); code/apparatus use `0≤l≤n` (→ 15). Cannot sync the paper to a possibly-buggy code; needs intent determination against `pro_system.py`.

**REMAINING:**
- Paper 57 60-vs-62 catalogue count (PM-fixable; needs counting the catalogue to fix the wrong number + 23-vs-24 calibration).
- "Build better validation" bucket (PI principle): Paper 56 injection-theorem genuine backing; Paper 57 P5 genuine validation; + the two determinations above.
- Coverage gaps (matrix rows; missing tests for 18/54/56/57); remaining citations (deligne-milne "Thm 2.3", EMN "G∉MT(ℤ)"); tooling (harden `check_k_label` for `\begin{conjecture}` envs + `docs/` scan; pin `/qa` reviewers to the worktree path).

## Remediation pass 2 (2026-06-16) — validation bucket (PI-approved corrections)

All three "build better validation" items were investigated (3 dispatched agents) and the PI-approved corrections applied + compile-clean; group3 deterministic gate GREEN throughout:
- **Paper 56 injection theorem** — the genuine period-map test (`tests/test_paper56_injection_g4_periodmap.py`) REFUTED C4 (M3 collapses to rank-1 like M1/M2; the `eye(n)` was unfaithful). Per PI option B, `thm:injection_g4` rewritten to the **abelianized / rank-1 image (Reading A)** across theorem, C4, proof, panel, honest-scope, abstract, and the scattered "closed immersion" mentions. The genuine test ships as a strict-xfail falsifier.
- **Paper 40** — verdict STATE-SPACE-GH-ONLY (the proof body itself concedes the dual-reach `reach_P` is a named gap — the same defect Paper 38 was descoped for). **Descoped** title + abstract + `thm:main_intro` + `thm:main` + corollaries from Latrémolière propinquity → van Suijlekom state-space GH, + general-$G$ conditionality (mirrors the completed Paper 38 fix). The Paper-40 internal cites in Papers 55/57 retitled to match (C11 group3 PASS).
- **Paper 57 P5** — verdict INHERENTLY CIRCULAR (packing_reachable ≡ F/C by construction; no F/C-blind check discriminates). **Reframed** the 98.3% as an internal-consistency check (not a discriminator) across abstract, §1, §5.5, Table caption, §5.6, §6.1; the genuine discovery (the two-family MF/period split) now carries the narrative.
- **Synthesis cascade** — group3 synthesis brought into line: Paper 56 "theorem-grade closed immersion" → abelianized homomorphism (not a closed immersion); Paper 57 98.3% "discriminator" → internal-consistency check.

**Still remaining (lower-severity):** coverage gaps (`claim_test_matrix` rows + collected tests for the 18/54/56/57 `debug/`-only claims); two citation items (Paper 56 `deligne_milne1982` "Thm 2.3" likely-misnumbered; EMN "G∉MT(ℤ)" verify); tooling (harden `check_k_label` for `\begin{conjecture}` + `docs/` scan; pin `/qa` reviewers to the worktree path); the Paper-40-retitle cascade into group5 (`paper_41`) + the group1 synthesis (out of group3 scope — for those branches' sweeps); and retiring the superseded tautological C4 (`eye(n)`) in `tests/test_paper56_injection_g4.py`.

## Remediation pass 3 (2026-06-17, v4.20.5) — lower-severity tail cleared

- **Citations (verify-then-fix; citation-reviewer verified all 4):** Paper 56 `\cite[Thm.\ 2.3]{deligne_milne1982}` → `\cite[\S5]{deligne1990}` (the exterior-tensor-product theorem is Deligne 1990 §5; the cited "Thm 2.3" does not exist in Deligne–Milne 1982); Paper 18 `deligne2010` wrong arXiv `math/0302267` removed (it is the Deligne–Goncharov 2003 paper); Paper 18 `glanois2015` venue → *J. Number Theory* **160** (2016) 334–384. (EMN 2025 + Charlton–Hoffman "Thm 2.21" verified GROUNDED — no fix.)
- **Coverage gaps:** added a Group-3-bite-2 section to `docs/claim_test_matrix.md` (14 rows) recording the load-bearing claims + the **NO-TEST** gaps honestly (Paper 18 α²-Ihara; Paper 54 selection-rules / DF; Paper 56 PS-4; Paper 57 P5 / count); the keystone Paper 56 injection row is FALSE-POSITIVE → REFUTED, with the genuine falsifier named.
- **Superseded test:** the tautological C4 tests in `test_paper56_injection_g4.py` (`gram = eye(n)` / `assert True`) skipped with a pointer to the genuine `test_paper56_injection_g4_periodmap.py` — they asserted the now-refuted injectivity (37 passed / 10 skipped).
- **Tooling — `check_k_label.py` hardened:** (1) catches the K-rule inside a tier-asserting environment (the conjecture-env blind spot) with the same NEG/NONSEL compliance escapes; (2) scans `docs/*.md` data files (`forced_free_seam.md`); + a UTF-8 stdout fix. Trunk + group3 PASS; regression `test_k_label_check.py` passes; discrimination preserved.
- **Tooling — `/qa` reviewer path-pin:** `qa.md` step 4 now mandates absolute worktree paths + forbids reading the real corpus (the S8 synthesis-reviewer miss).
- **Subsumed:** the "Paper-38-propinquity deterministic check" candidate is now covered by the C11 enforcement (v4.20.4) — propinquity titles are enforced corpus-wide, so a Paper-38/45 cite drifting from the state-space-GH title FAILS C11.
- **Found (pre-existing, out of bite-2 scope):** `paper_18` body has an undefined `\cite{paper7}` (l.176, missing bibitem) — older body debt, flagged for a future pass.

**Bite-2 remediation is now fully closed.**

## Remediation checklist

### QUICK WINS (wording / titles / labels — PM-fixable on sight)
- [ ] **Paper 38 status drift (C7), corpus pattern.** Cited as "Latrémolière propinquity" instead of "van Suijlekom state-space GH": Paper 54 (intro l.67), Paper 55 (~6 places incl. bibitem), Paper 57 (bib + l.143 prose). [R] Sweep group3 bodies; the v4.14.1 retitle never propagated here.
- [ ] **Paper 55 `hoffman2019odd` fabricated title** → "An odd variant of multiple zeta values" (ID/venue/pages correct). [R]
- [ ] **Paper 57 `boyle_farnsworth2018` wrong-ID.** Title belongs to arXiv:1910.11888 but cited ID 1604.00847 — and the DGA/Higgs-unit-vector content invoked IS the 1604 paper → **fix the TITLE to match the ID, do NOT swap the ID**. [R]
- [ ] **Paper 57 C6 (l.113)** "graph Laplacian spectrum λ_n = n²−1" → continuum/converges-to phrasing (graph L=D−A is positive-semidefinite). [R]
- [ ] **Paper 18 conjecture-environment K-labels (C5/§13.5).** Two `\begin{conjecture}` envs at l.2429 and l.2602 `[Combination rule, restated]` + "conjectural"/"derivable" prose (l.2483, 2675, 3234). [V] Demote to `observation`/plain prose; K stays an Observation.
- [ ] **`docs/forced_free_seam.md` + `debug/principle_hunt_audit.py`** still label K "conjectural per §13.5" (companion data files, gate-invisible). [R] Sweep to Observation.
- [ ] **Synthesis (C9)** "theorem-bound boundary" (Paper 57 reconvergence) → Paper 57 is an Observations paper; soften to "empirically-bounded, partially theorem-bounded by 31/55/56." [V via primary text]

### SUBSTANTIVE (author/PI-level — not wording fixes)
- [ ] **Paper 56 — headline injection theorem unbacked [LARGE, C2].** `thm:injection_g4` ("theorem-grade closed immersion") is backed by **entirely tautological tests**: C1 `simplify(x*x − x*x)==0` (no period map), C4 `gram = eye(n)` hardcoded, `assert True`. [V] No period map / Hurwitz evaluation / Gram-of-period-vectors exists in `geovac/tannakian.py`. **The keystone result of Paper 56 has no genuine backing.** Either build real backing or downgrade the theorem-grade claim.
- [ ] **Paper 56 — contradictory primitive-space dimension [LARGE, C8/C3].** §hgv defines N_sec = n(n+1)/2 (→9 at n=2); the whole injection/equality apparatus uses N = n(n+3)/2 (→15). Code (`pro_system.py`) uses 0≤l≤n (→15). Reconcile the paper definition to the code/panel convention.
- [ ] **Paper 57 — P5 98.3% tautological + unbacked [MAT, C1/C2].** `packing_reachable` is a hand-tagged column ≡ F/C status; the 98.3% measures tag-consistency, not a discovered predictor. Only artifact is `debug/principle_hunt_audit.py` (not pytest-collected, no matrix row). [R] Needs an independent witness-derivation OR an explicit tier caveat + a collected test.
- [ ] **Paper 57 — catalogue count inconsistency [MAT, C8].** 98.3% = 59/60, but inventory stated 62 (38F/1A/23C); executable catalogue = 60 (35F/1A/24C); 23-vs-24 calibration. Reconcile against `docs/forced_free_seam.md`. [R]
- [ ] **Paper 57 — "dim M(D_F)=128 per generation" [MAT, C3].** Restated at a tier its own test (`test_trunk_qa_forced_count_moduli.py`) contradicts: full-axiom moduli = 260; 128 is matter-sector only (already DOWNGRADED for Paper 32 in the matrix). Add the matter-sector caveat or correct to 260. [R]
- [ ] **Paper 18 — κ "we now derive … more fundamental" [LARGE, C3/C7/§1.5].** §VII.F titled "Derivation of κ"; asserts conformal reading "more fundamental" (ontological priority). κ is an Observation (trunk v4.13.0). [R] Retitle to "two converging routes"; drop "more fundamental."
- [ ] **Paper 18 — keystone Eq.(zeta_d2) backing [MAT, C2].** Proof-sketch cites a non-existent "80-digit PSLQ" test; symbolic half tautological (real backing is a 2-pt float spectral sum). [R] Repoint to the genuine backing; drop the PSLQ claim.

### COVERAGE GAPS (C1 — no collected test)
- [ ] **No `claim_test_matrix.md` rows for any bite-2 paper** (18/54/55/56/57). Add them.
- [ ] **Paper 18** α²-Ihara degree-44 / −35/648 headline — no test ("14/14 tests" cite covers infrastructure only). [R]
- [ ] **Paper 54** §df_thc_connection bit-exact DF claims — `debug/` drivers only (pytest-ignored). [R]
- [ ] **Paper 56** PS-4 endo-rigidity (872 residuals, ~15% of panel) — `debug/` driver only. [R]

### CITATION (C4 — verify then fix)
- [ ] Paper 56 `deligne_milne1982` "Thm 2.3" likely misnumbered (exterior tensor product is in Deligne 1990, not a numbered DM-1982 theorem) — backs the combined-equality corollary. [R]
- [ ] Paper 56 `eskandari_murty_nemoto2025` "G ∉ MT(Z)" half unverifiable from abstract (forces the level-4 refinement) — read the source. [R]
- [ ] Paper 55 Charlton–Hoffman "Thm 2.21" theorem-number unverifiable. [R]
- [ ] Paper 18 OWN bib: `glanois2015` wrong venue/vol/year; `deligne2010` appended arXiv:math/0302267 = a different (Deligne–Goncharov) paper. [R]

### TOOLING (deterministic-gate hardening — the LLM backstop caught these, harden the screen)
- [ ] **`check_k_label.py` blind spot 1:** misses K inside `\begin{conjecture}` environments (paper_18 passed `--gate group3` while wrapping K in two conjecture envs). [V]
- [ ] **`check_k_label.py` blind spot 2:** scans only `papers/**/*.tex`; misses `docs/` + `debug/` data files carrying prohibited K-labels (`forced_free_seam.md`). [R]
- [ ] **New deterministic check candidate:** Paper-38 cite says "propinquity" near a Paper-38 reference → flag (same class as internal-titles; would have caught the C7 drift mechanically).
- [ ] **`/qa` reviewer prompt fix:** dispatched reviewers must be pinned to the worktree path and forbidden from cross-reading the real corpus (the synthesis reviewer's S8 miss). Bake into the reviewer dispatch template.

---

## Honest scope
- **PM-verified [V]:** worktree S8 persistence (synthesis reviewer miss confirmed); Paper 18 conjecture envs (real corpus l.2429/2602); Paper 56 tautological injection tests (real corpus `assert True`/`gram=eye`); synthesis "theorem-bound boundary" vs Paper 57's Observations-paper status.
- **Reviewer-reported [R]:** everything else — high-confidence (quoted lines, tests run) but verify against primary text before each edit (per the §9 reconcile rule).
- **Remediation completed (v4.20.0).** Passes 1–2 above edited the bite-2 papers + synthesis; all compile-clean; group3 deterministic gate (C10–C13) GREEN. The original FAIL findings record is preserved above for provenance.

### Honest scope — post-remediation (v4.20.0)
- **Theorem-grade / settled:** Paper 56 C4 closed-immersion injectivity is **refuted** at theorem grade by a genuine computation (M3 period-vector Gram rank 1, det 0, bit-exact at n_max=2,3; falsifier `tests/test_paper56_injection_g4_periodmap.py`, 14 passed + 2 strict-xfail). The corrected theorem (abelianized / rank-1 image, Reading A) is what the construction supports. Paper 57 P5 circularity is **settled** (exact cross-tab: `packing_reachable` ≡ F/C by construction; an F/C-blind witness check scores the 60% base rate).
- **Structural determination (not a new theorem):** Paper 40 = state-space GH, not Latrémolière propinquity — established from the paper's *own* conceded `reach_P` named gap (the same defect as the Paper 38 descope), not from a fresh proof. The general-$G$ statement is *conditional* (per-group spin-window decomposition verified only for SU(2)).
- **Numerical observation:** Paper 57's 98.3% is reframed as an internal-consistency check (not a discriminator); the genuine signal is the two-family MF/period split.
- **Untouched genuine content:** Paper 40's closed-form rate γ_Λ(G), the universal 4/π theorem, C₃(G) ≤ 1; Paper 56's C1/C3 (which do hold); Paper 55's mixed-Tate classification.
- **Open follow-ons (lower-severity, not in this close):** coverage-gap matrix rows + collected tests; two Paper-56 citation verifications (deligne-milne "Thm 2.3", EMN); `check_k_label` hardening (conjecture-envs + `docs/` scan) + reviewer-path-pin; the Paper-40-retitle cascade into group5 `paper_41` + the group1 synthesis; retire the superseded tautological C4 in `tests/test_paper56_injection_g4.py`.
- **Re-cert:** bite 2 is *remediated*, not *re-certified* — a fresh PI-invoked `/qa group3` (synthesis dimension worktree-pinned) is what would certify it.
