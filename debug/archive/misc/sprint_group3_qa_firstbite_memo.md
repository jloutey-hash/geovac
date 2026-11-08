# Sprint memo — group3 QA first bite (Papers 22/24/31 + synthesis) + Paper 8 cross-branch closure

**Dates:** 2026-06-16 (runs #1–4). **Version:** v4.17.0.
**Verdict:** group3 first-bite (Papers 22/24/31 + synthesis) driven through four `/qa group3`
re-cert cycles; FAIL character converged prose → labels → test-coverage → peripheral, all
dispositioned. One real production bug fixed (Wigner-3j ×2). One stale-ref investigation
spilled into group2 and closed a dead-end (Paper 8 phase-lock probe).

> Companion records (do not duplicate): per-run detail in `memory/branch_qa_sweep_phase.md`;
> per-topic diagnostics in `debug/group3_density_diagnostic_memo.md`,
> `debug/group3_followup_3j_bug_memo.md`, `debug/group3_followup_density_layers_memo.md`;
> per-run seed keys `debug/qa/group3_seed_key{,_r2,_r3,_r4}.json`; DoD `docs/qa/group3.done.md`.

## 1. Scope & method

First bite of the branch-by-branch QA sweep beyond the trunk. Subset = the three non-trunk
group3 keystones (Paper 22 angular sparsity, Paper 24 Bargmann–Segal, Paper 31 universal/Coulomb
partition) + the group3 synthesis. Pre-registered DoD `docs/qa/group3.done.md` (frozen). Each
run: deterministic screens (C5 K-label, C6, C10 compile, C11 titles, C13 paper↔test refs) +
blind seeded worktree + fresh per-dimension reviewers (code/claims/citation/synthesis) +
calibration scorecard + PI disposition + re-run.

## 2. Run-by-run (FAIL converged, one layer per run)

- **Run #1** — panel 5/5 sensitive, 5/5 specific. Mechanical defects PM-fixed: stale K-labels in
  P24/P31 (incl. ones the deterministic screen MISSES — `strip_latex` blanks `\alpha` so the
  `α-conjecture` anchor never fires, and a bare "…constant conjecture" isn't anchored ⇒ the
  `claims-reviewer` C5 enumeration is essential, not just a backstop); C6 discrete-vs-continuum in
  24/31; citations (Dunlap/Dyall/Grant in P22, nieuviarts title in P31); overclaim rescopes (P22
  Cor-1 O(Q^2.5), rank-2 "theorem"→numerical); layer count synced to Paper 24's authoritative **six**.
  **Real bug found + fixed:** Paper 24 Wigner-3j off by ×2 (`bargmann_graph.py` denominators
  `2l+2→l+1`, `2l→l`; uniform global rescale ⇒ no published number moved; test `1/15→2/15`, 17/17 green).
- **Run #2** — fresh panel caught run-#1 residuals (incomplete dispositions): 2 stale "conjecture"
  K-labels still in P24 (incl. a §subsection *heading*, invisible to the screen); layer-count
  contradictions (P24 footnote vs list; P31 abstract "three" vs a `\section{The Five-Layer…}`
  title); P22 density table reporting D_pd as if it matched the headline D; one κ "set by Jacobian"
  overclaim in P31. All PM-fixed. **Process lesson:** code-reviewer OK'd a *gutted* planted test
  because the claim had redundant sibling backing → "score each cited test individually."
- **Run #3** — panel fully calibrated (6/6 incl. code dimension, redundant-backing blind spot
  fixed). FAIL shifted to **NO-TEST coverage**: P22 headline density, P24 HO zero-entropy, P31
  two-body numbers all lacked tests. Closed by 3 coder agents: `test_paper22_density.py`,
  `test_paper24_ho_entropy.py`, `test_paper31_two_body.py`. Test-writing **confirmed the
  `angular_zero_count` mislabel** (computes D_pd while docstring/footnote say global-M_L) —
  regression-pinned, NOT fixed (→ CF-1).
- **Run #4 (this session)** — panel calibrated: sensitivity 6/6 (all seeds, all dimensions, incl.
  code S3 = the `assert S < 1e10` planted as the *sole* entropy assertion); 1 reconciled false
  positive (P31 code-reviewer flagged the correct, tested global-D=6.06% as "not reproducible" —
  it's backed by `test_paper22_density.py`, which that P31-scoped reviewer never opened →
  **cross-paper-test blind spot**, not a real defect). Three genuine residuals, all PM-cleared:
  1. **P22 spinor Table II — NO-TEST [MATERIAL].** Spinor (jj-coupled) density column had no live
     backing; the cited driver `debug/tier2_t0_spinor_density.py` was archived. Wrote
     `tests/test_paper22_spinor_density.py` — exact `sympy` jj enumeration pinning both spinor
     columns (full-Gaunt + pair-diagonal) for l_max 0–5 + the three structural checks (Q=2(l+1)²,
     pair-diag⊆full-Gaunt, d_sp≤d_sc). Paper §"Implementation and verification" rewritten to cite
     the live test and repoint the archived driver. Matrix row added.
  2. **P24 `higgs_pickrell2025` — misattributed [SMALL].** Verified vs arXiv:2503.23549 directly:
     it's a *genuinely spherical* oscillator paper calling the **2-sphere critical**, motivated by
     2d chiral QFT — never centers d=5/SO(6), states no "d≥3 open question = our d=5." Rewrote to an
     honest distinct-construction pointer; dropped the fabricated open-question sentence.
  3. **P24 `jauchhill1940` — converse attribution [SMALL].** Confirmed Jauch–Hill 1940 proves the
     *forward* direction (oscillator HAS su(N) symmetry via N²−1 constants). Narrowed all 3 mentions
     to forward; attributed the converse/biconditional uniqueness to Wybourne/Iachello (already cited).

## 3. Density finding (closed as a D-vs-D_pd relabel)

Two distinct densities, now both pinned and named everywhere (P22/31 + synthesis + matrix):
- **D** = global-M_L Coulomb density (physical selection rule), **6.06% @ l_max=3** (14.84/8.52/6.06/4.83/3.99% l_max 1–5).
- **D_pd** = pair-diagonal *realized* density (stricter axial sub-case), **1.44% @ l_max=3** (7.81/2.76/1.44/0.90/0.62%).
- Production `composed_qubit._ck_coefficient`/`angular_zero_count` realizes **D_pd** (drops m-swap ERIs);
  `casimir_ci._gaunt_ck` realizes global **D**. The `angular_zero_count` docstring/Paper-22 footnote say
  global-M_L but the routine returns D_pd → **CF-1** (group4 carry-forward, regression-pinned, not fixed).

## 4. Paper 8 cross-branch closure (out-of-scope item that kept surfacing)

C13 flagged `Paper_8_Bond_Sphere_Sturmian.tex:1176 → test_harmonic_phase_lock.py NOT FOUND`
(group2, advisory under the group3 gate). Investigation:
- The cited path was **never in git**; the file has only ever existed at
  `debug/archive/test_harmonic_phase_lock.py`, committed once at v0.9.37 (2026-03-11 codebase lock).
- The script is a **failed speculative probe**: tests whether SO(4) Wigner-D critical points predict
  LiH R_eq≈3.015 bohr via inverse stereographic γ→R. Its own output: **0 strong/near hits** across
  5 elements × 5 p₀. Correctly archived as a dead-end (predates the §3 discipline).
- The two identities Paper 8 cites (Thm 1 σ-bond rules D²_(1,0),(1,0)≡1, D²_(0,0),(1,0)≡0) are an
  *incidental byproduct* and **genuinely true** (re-verified 4.4e-16 / 5.6e-17 on the paper's grid).
- **Fix (B):** wrote `tests/test_paper8_sigma_bond_selection.py` (5 tests: grid <1e-14 + in-paper
  closed forms cos²+sin²=1 / antisym cancellation + distinctness guard), repointed the citation.
  Paper 8 compiles fully clean (ERRORS=0, **UNDEF=0**). C13 `--gate group2` now PASS.
- **Convention (C):** added `docs/authoring_conventions.md` rule 10 — `test_*.py` citations are
  verification-backing (must be live `tests/`, gated by C13); `debug/` is provenance-only (may be
  archived). Paper 8 is the worked example.
- **§3 dead-end row** added for the phase-lock probe.

## 5. Bonus (flagged): stale hard-prohibition labels in the conventions doc

While in `docs/authoring_conventions.md`, found rule 4 saying "K… stays CONJECTURAL" and the group5
bullet "α combination rule conjectural" — both contradict the §13.5 hard prohibition (Observation
since the 2026-06-14 downgrade). Worse than stale paper prose: *guidance instructing a violation*
(the `paper_57:446` smoking-gun pattern). Corrected both to "Observation." Surfaced to PI.

## 6. Honest scope

- **Theorem grade / SYMBOLIC PROOF (closed):** P22 scalar D & D_pd densities (exact sympy + production
  `_gaunt_ck`); P22 **spinor** densities both conventions (exact jj 3j, this session); P24 Wigner-3j
  value (bug-fixed, vs sympy); P8 Thm 1 σ-bond identities (grid + closed form, this session). These
  are exact angular-combinatorial / identity facts, fully backed.
- **Panel-verified / measured:** P22 Thm 2 potential-independence (5 potentials, bit-identical mask);
  P24 HO zero-entropy (S=0, ‖[H,V]‖<1e-15); P31 two-body connected-fraction + Pearson (honest negatives).
- **Numerical observation (unchanged tier):** P31 two-body radial Pearson 0.41–0.58 decreasing — a §3
  dead-end, correctly labeled.
- **Closed dead-end:** SO(4) D-matrix critical points → LiH R_eq (0 hits); now §3-recorded.
- **Process calibration:** the panel is demonstrably calibrated (run #4: 6/6 sensitivity, 1 reconciled
  cross-paper-test false positive). The FAIL verdict is *trustworthy*, not a reviewer artifact.
- **Named open follow-ons (NOT closed this sprint):**
  1. **CF-1** (group4): `angular_zero_count` D_pd-vs-global-M_L mislabel + the composed Pauli
     advantage's dependence on the pair-diagonal approximation (LiH 333→837 Pauli = 2.5×; "2.7× fewer
     than STO-3G" → ~parity; energy ~1 mHa). PI/group4 decision pending.
  2. **CLAUDE.md §1.6 "1.44%"** tag (Phase-4 section, PM-locked) — PI-only edit, deferred.
  3. **C13 debug-vs-tests convention** now documented; other group2 papers may carry the same
     `debug/` citation pattern — to be swept when group2 is QA'd.
  4. group3 first bite was Papers 22/24/31 only; Papers 18/54/55/56/57 remain for the full group3 cert.
- **NOT a `/qa` PASS:** run #4 ended FAIL (the 3 residuals), now remediated. Per the qa.md hard rule,
  re-cert timing is the PI's; the PM does not self-trigger `/qa`.

---

## Runs #5–7 + C11 hardening + corpus title cleanup (2026-06-16, v4.18.0) — **first-bite CERTIFIED (PASS)**

The first-bite subset (Papers 22/24/31 + synthesis) reached a calibrated `/qa group3` **PASS** at run #7, after the C6 class was closed thoroughly and the deterministic title gate was hardened.

**Run #5 (FAIL → calibrated).** 6/6 sensitivity, 6/6 specificity (the v4.17.0 fixes held). 3 genuine residual **C6** discrete-vs-continuum slips found (P31 l.78/l.450, P24 l.64) — the class runs #1–2 only partly swept. PM-fixed.

**Run #6 (FAIL → calibrated).** 6/6 sensitivity (incl. S7 in the previously-blind format), 6/6 specificity. **Two MORE C6 residuals** (P31 l.816, P22 l.95–96 — they span line breaks, so the line-based sweep missed them) + a genuine **deterministic-gate blind spot**: `check_internal_titles.py` only validated bibitems ending "GeoVac Paper N (YEAR)"; P22's internal bibitems end "GeoVac Technical Report (YEAR)" and were **silently uncertified** (S7 in that format slipped past). Lesson: patch-and-rerun *leaks* on a class — sweep the whole class at once + harden the detector.

**Disposition (the fixes that earned the PASS):**
- **Thorough C6 sweep** (multiline, all 4 files): fixed P22 l.95, P31 l.816, SYN l.382 + l.668 (the 2 sites the line-grep missed + the synthesis "(D−A) produces … n²−1" / "graph whose Laplacian eigenvalues *are* n²−1"). All now attribute −(n²−1) to the continuum the graph *converges to / is conformally identified with*. Verified by the multiline scan judging every candidate (P24 was already clean; "discrete spectrum" of a Sturm–Liouville operator and the HO Bargmann diagonal N+3/2 are correct usages, not C6).
- **C11 hardened:** `KEYED` pattern resolves the paper number from the `\bibitem{GeoVac_PaperN}` key (covers the "Technical Report" format); `main_part()` now strips `:`-style subtitles (`: ` / `:\\ ` / `:\ `); `--gate <branch>` added (mirrors check_k_label / check_paper_test_refs). Self-test confirms it now catches the previously-blind P22 format.
- **borelweil** (P31) → credits Serre as expositor; **propinquity** (P31 ×2) → WH1-PROVEN substrate renamed "van Suijlekom state-space GH" (was "Latrémolière propinquity").

**Run #7 = PASS.** 6/6 sensitivity (S7 now caught in the Technical-Report format), zero false positives; the **C-G control held** — all four claims reviewers tagged every C6 sweep fix CORRECT and did not re-flag them (the sweep actually closed the class; no new C6 surfaced). All five dimensions exercised/calibrated/clean. Only genuine finding a synthesis NIT (bare "D" vs D_pd in one bullet, universality true either way).

**Corpus-wide internal-title cleanup (the "jump in on other branches" pass).** The hardened check surfaced 6 pre-existing stale internal-title cites the old check was blind to — all verified GENUINE (diagnostic prints normalized cited-vs-real), all fixed to the current `\title`:
- trunk/group1 **paper_32**: P14 ("…from Natural Geometry"→"…from Spectral Graph Theory"), P23 ("Nuclear Shell Model **Qubit** Hamiltonians…"→"…on the Hyperspherical Lattice") — *the trunk cert had the same blind spot*;
- group5 **paper_2_alpha**: P2 self-cite (v4.14.2 retitle straggler);
- group6 **paper_34**: P15, P19, P42 (P42 a real subtitle drift "…Theorem…"→"…Literal Identification…").
`check_internal_titles` now **PASSES bare, corpus-wide**; the regression test restored to the bare (corpus-wide) invariant.

**Dyall §9.3 / Grant §7.5 (P22).** The Oxford TOC confirms Dyall Ch. 9 = "Operators… under Time-Reversal Symmetry" → `§9.3` is wrong for the two-electron Coulomb reduction. Correct section unconfirmable from the web (Grant TOC didn't surface); dropped the section pointers, cite the books at work-level "(Dyall–Faegri, Grant)". Section numbers wanting the PI's physical-book check to restore.

### Honest scope (runs #5–7)
- **CERTIFIED (calibrated PASS):** group3 first-bite subset Papers 22/24/31 + the synthesis's coverage of them — survived the seed-catalog defect classes across all five dimensions, panel demonstrably discriminating (6/6 sensitivity each run).
- **Tooling hardened (not a physics result):** C11 deterministic title gate (KEYED + subtitle + --gate); the class can no longer silently recur.
- **Mechanical corpus hygiene:** 6 stale internal-title cites fixed corpus-wide; C6 wording precision.
- **Open follow-ons:** full group3 cert (Papers 18/54/55/56/57 + full synthesis); a **trunk re-touch** is warranted (the hardened C11 now covers paper_32's keyed bibitems, just cleaned, but the trunk should be re-run against the hardened gate); Dyall/Grant section numbers (PI books); CF-1 (group4) still open; the group2 sweep should expect `debug/`-citation strays (C13) + now-covered keyed bibitems (C11).
