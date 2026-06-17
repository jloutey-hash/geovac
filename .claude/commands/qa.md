---
description: The /qa gate — adversarial, calibration-tested QA review of a target (trunk or a branch) returning a trustworthy PASS / FAIL / INCONCLUSIVE. PI-invoked ONLY; the PM must never self-trigger it. Usage: /qa <target> (e.g. /qa trunk, /qa group1).
---

# /qa — the QA certification gate

This is a deliberate, PI-timed certification gate, **not** a routine pass.

> **PI-INVOKED ONLY.** The PM must never self-trigger this, run it proactively, or nudge toward it each sprint. It fires only when the PI types `/qa`. Certifying a branch "done" is a timing judgment that belongs to the PI; auto-firing would defeat that. (It is a genuine-trigger command like `/release` and `/sprint-close`, with **no** corresponding memory rule — CLAUDE.md §13.9a.)

**Why this exists.** A reviewer tasked with "find problems" always finds problems, so a "clean" can never be *earned* and a "dirty" can never be *trusted* — it is an activity, not a test. This gate converts the fault-hunt into a test with discriminating power, using four mechanisms (pre-registered criteria, materiality grading, per-run calibration controls, convergence) and returns a **three-way verdict** that distinguishes *"target not done"* from *"reviewer not trustworthy."*

## Inputs

- **Target** = `$ARGUMENTS` (`trunk`, `group1`, …). If absent, ask the PI which target.
- **Definition of done** = `docs/qa/<target>.done.md` — the pre-registered, binary, material criteria. **If it does not exist, STOP** and co-write it with the PI first. Pre-registration must precede the review; the reviewer is never allowed to define "done."
- **Seed catalog** = `docs/qa/seed_defects.md` — the known defect classes to plant.

## Protocol

1. **Load & freeze the criteria.** Read `docs/qa/<target>.done.md`, quote it back, and confirm with the PI that nothing has been added or relaxed for this run. The goalposts are fixed before any review begins. **Run the DETERMINISTIC checks the criteria name** — they certify their classes more reliably than an LLM reviewer (the run-#1 / run-#4 lessons), leaving the dispatched reviewers to *judgment* calls:
   - **C11 internal-title consistency** — `debug/qa/check_internal_titles.py` (exact string match → clean PASS/FAIL).
   - **C5 K-label backstop** — `debug/qa/check_k_label.py` (FAILs on any SUSPECT K-tier assertion in the gated **trunk** scope; the **AUDIT** list it prints is advisory non-trunk drift, surfaced for a PI sweep decision, and does **not** fail the gate). This is a *backstop*, not a replacement: the K-prohibition is semantic, so the `claims-reviewer` still enumerates K-sentences (the screen can't anchor a bare "K" without the rule name in-clause).
   - **C13 paper↔test reference integrity** — `debug/qa/check_paper_test_refs.py` (FAILs if any test cited inline in a gated paper is non-existent/renamed/archived; reports matrix-coverage gaps + non-trunk stale refs as advisory). The deterministic guard that the inline test references stay honest against the suite + `docs/claim_test_matrix.md`. `--gate <branch>` per sweep.
   - **C14 paper↔file reference integrity** — `debug/qa/check_file_refs.py` (FAILs if any `geovac/`, `benchmarks/`, or `demo/` — permanent code/artifact — file path cited inline in a gated paper does not resolve; `debug/` refs are advisory, since the transient clean-room dir is pruned by design). Closes the C13 gap (C13 covers `tests/` refs only): the group3 run-4 `geovac/jlo_chi.py` defect was a claim citing a nonexistent code module. `--gate <branch>` per sweep.

2. **Build an isolated seeded copy.** `git worktree add ../geovac-qa-seed-<target> -b qa-seed-<target>` (or reuse a temp branch). **Seeds NEVER touch the real corpus** — all planting happens in the worktree, which is deleted at the end.

3. **Plant calibration controls (blind).** From `docs/qa/seed_defects.md`, plant K (≈4–6) *realistic* defects from our actual §3 error history (fabricated arXiv ID, circular/tautological test, κ-style overclaim, stale-title drift, hard-prohibition touch, …) into the worktree copy of the target. **Cover every gating dimension:** plant at least one seed catchable by *each* dimension's reviewer type (a code-defect for `code-reviewer`, a prose-defect for `claims-reviewer`, a citation-defect for `citation-reviewer`, a synthesis-defect for the synthesis `claims-reviewer`), so no dimension that gates the verdict goes uncalibrated this run. Also designate M (≈4–6) **known-good** claims that are verified correct and must NOT be flagged material. Record the answer key — {file, location, defect-class, expected-catcher} — in `debug/qa/<target>_seed_key.json`. The reviewers must be blind to which items are seeded.

4. **Dispatch the reviewers across EVERY dimension — fresh, blind, pointed at the worktree.** The target's criteria partition into independent **review dimensions**, and a single `/qa` run MUST exercise *all* of them in one invocation. **Passing one dimension never excuses skipping another, and the run does not terminate early when a dimension comes back clean** — the dimensions are independent axes, not sequential gates, so you want the full defect picture from all of them every run. For the trunk the dimensions, their reviewer types, and the criteria they gate are:
   - **Code / test-backing** (C1–C2) — **`code-reviewer`**, one per paper that has backing tests: claim ↔ test ↔ code; *RUN the tests*; does each genuinely *prove* its claim (not tautological / false-positive / weaker than prose)?
   - **Paper claims / prose** (C3, C5, C6, C8) — **`claims-reviewer`**, one per paper: does any claim assert more than its backing? is the tier inline? are the hard prohibitions intact?
   - **External citations** (C4) — **`citation-reviewer`**, one per paper: do external cites resolve and say what we attribute?
   - **Synthesis faithfulness** (C9) — **`claims-reviewer`**, one on the branch synthesis: zombie/descoped claims? status overstatement? faithful to the papers it summarizes? (This is a *separate* dispatch from the per-paper claims review — the synthesis is its own document.)
   - **Deterministic** (C10 compiles, C11 internal titles) — the scripts in step 1, not an LLM reviewer.
   Each LLM reviewer reviews **against the pre-registered criteria**, classifies findings MATERIAL vs NIT, and attaches a counterfactual to each MATERIAL finding ("if fixed, would a *target result* change?"). They are fresh agents that never saw each other or the discovery-mode rationale (§9 QA principle 2). **If a gating dimension cannot be exercised this run, the target verdict is INCONCLUSIVE — never PASS (step 7).**

   **Path-pin (mandatory — the group3 bite-2 S8 miss).** Every reviewer prompt MUST give the ABSOLUTE **worktree** path (`../geovac-qa-seed-<target>/…`) for the files it reads, and MUST explicitly forbid reading the real corpus (`papers/…` under the project root). A reviewer that silently re-reads the real (unseeded) corpus misses every planted defect in its dimension and **de-calibrates the run with no error at all** — a synthesis `claims-reviewer` did exactly this:\ it first caught the planted synthesis seed, then *reversed* the catch after re-reading the live (unseeded) file. Before trusting a reviewer's calibration score, spot-check that its quoted line numbers/text came from the worktree, not the real corpus.

5. **Score calibration — certify the reviewers (against the answer key).**
   - **Sensitivity:** fraction of planted defects the panel collectively caught. A planted class no reviewer caught ⇒ the panel is *blind to that class*.
   - **Specificity:** any known-good control flagged MATERIAL is a *false positive*.
   - A panel that misses planted defects or manufactures faults on the controls is **not calibrated for this run**.

6. **Materiality + convergence on the genuine findings.** Set aside findings at seed/control locations. For the rest, the PM **verifies every MATERIAL finding against primary text/code before accepting it** (§9 reconcile rule). Converge: take the union of *verified* material findings across the independent reviewers.

7. **Emit the verdict — per dimension, then rolled up (AND across dimensions).** Score *each* dimension separately (was it exercised? calibrated? clean?), then combine. The target's verdict is the AND over all gating dimensions:
   - **PASS** — *every* gating dimension was **exercised**, **calibrated** (caught its plants, zero false positives on its controls), and **clean** (zero verified material defects), and independent reviewers converged. *Trustworthy clean across all dimensions.*
   - **FAIL** — every gating dimension calibrated **and** ≥1 verified material defect in *any* dimension. *Target not done* — list each defect with its **dimension** and failing criterion.
   - **INCONCLUSIVE** — any gating dimension was **not exercised**, *or* was exercised but **not calibrated** (missed its plants, or flagged a known-good control). *The verdict is untrustworthy this run.* Name the dimension and what failed. **An unexercised gating dimension is INCONCLUSIVE — never silently downgraded to a "PASS on the exercised dimensions" with a ceiling note.** Fix and re-run.

   Report the verdict as a **per-dimension scorecard** (dimension → exercised / calibrated / clean) plus the roll-up, so the PI sees which axes actually carried the verdict and which (if any) forced INCONCLUSIVE.

8. **Clean up.** `git worktree remove` the seeded copy. Confirm no seed reached the real corpus.

## Output to the PI

- **Verdict:** PASS / FAIL / INCONCLUSIVE.
- **Calibration scorecard:** sensitivity (plants caught / planted, per defect-class) and specificity (false positives / controls) — so the PI can see the panel's *demonstrated* discrimination this run.
- **Material findings** (if FAIL): each verified against primary text, with its failing criterion and counterfactual.
- **Nit log:** non-material findings (logged, do not block).
- **Honest ceiling:** "PASS" means *survived the calibrated detectors for the defect classes in the seed catalog and the pre-registered criteria* — not provably perfect. Name any criterion or defect-class that could not be exercised this run.

## Hard rules
- **PI-invoked only.** PM never self-triggers, runs proactively, or suggests it each sprint.
- Pre-registered criteria are **frozen before** the review — no goalpost-moving in either direction.
- Seeds live only in the worktree; never commit/leak them; always remove the worktree.
- A reviewer's verdict is trusted only **after** it passes calibration (caught the plants, passed the controls). An uncalibrated panel ⇒ INCONCLUSIVE, never PASS.
- **All dimensions every run.** One `/qa` run exercises *every* review dimension the criteria reference — code/test-backing, paper-claims, citations, **and** synthesis — in a single invocation. Dimensions are independent: a clean one never short-circuits another, the run never stops early on a pass, and an **unexercised gating dimension forces INCONCLUSIVE, not PASS** (the run-#2 lesson: a "PASS on the exercised dimensions" hid real defects in the two dimensions that were skipped). The verdict is the AND across all dimensions.
