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

1. **Load & freeze the criteria.** Read `docs/qa/<target>.done.md`, quote it back, and confirm with the PI that nothing has been added or relaxed for this run. The goalposts are fixed before any review begins. **Run any DETERMINISTIC checks the criteria name** (e.g. internal-title consistency via `debug/qa/check_internal_titles.py`) — string-comparison criteria are certified by a script, more reliably than an LLM reviewer (the run-#1 lesson); the dispatched reviewers cover only the *judgment* criteria.

2. **Build an isolated seeded copy.** `git worktree add ../geovac-qa-seed-<target> -b qa-seed-<target>` (or reuse a temp branch). **Seeds NEVER touch the real corpus** — all planting happens in the worktree, which is deleted at the end.

3. **Plant calibration controls (blind).** From `docs/qa/seed_defects.md`, plant K (≈4–6) *realistic* defects from our actual §3 error history (fabricated arXiv ID, circular/tautological test, κ-style overclaim, stale-title drift, hard-prohibition touch, …) into the worktree copy of the target. Also designate M (≈4–6) **known-good** claims that are verified correct and must NOT be flagged material. Record the answer key — {file, location, defect-class, expected-catcher} — in `debug/qa/<target>_seed_key.json`. The reviewers must be blind to which items are seeded.

4. **Dispatch N independent reviewers (fresh, blind), pointed at the worktree.** One of each relevant type per paper/module:
   - **`code-reviewer`** — claim ↔ test ↔ code; does each test genuinely *prove* its claim (not tautological / false-positive / weaker than prose)?
   - **`claims-reviewer`** — papers + synthesis prose: does any claim assert more than its backing? does a synthesis faithfully reflect its papers?
   - **`citation-reviewer`** — external citations resolve and say what we attribute.
   Each reviews **against the pre-registered criteria**, classifies findings MATERIAL vs NIT, and attaches a counterfactual to each MATERIAL finding ("if fixed, would a *target result* change?"). They are fresh agents that never saw each other or the discovery-mode rationale (§9 QA principle 2).

5. **Score calibration — certify the reviewers (against the answer key).**
   - **Sensitivity:** fraction of planted defects the panel collectively caught. A planted class no reviewer caught ⇒ the panel is *blind to that class*.
   - **Specificity:** any known-good control flagged MATERIAL is a *false positive*.
   - A panel that misses planted defects or manufactures faults on the controls is **not calibrated for this run**.

6. **Materiality + convergence on the genuine findings.** Set aside findings at seed/control locations. For the rest, the PM **verifies every MATERIAL finding against primary text/code before accepting it** (§9 reconcile rule). Converge: take the union of *verified* material findings across the independent reviewers.

7. **Emit the verdict.**
   - **PASS** — panel calibrated (caught the plants, zero false positives on controls) **and** zero verified material defects on the genuine content **and** independent reviewers converged. *Trustworthy clean.*
   - **FAIL** — panel calibrated **and** ≥1 verified material defect. *Target not done* — list them with the failing criterion.
   - **INCONCLUSIVE** — panel **not** calibrated (missed planted defects, or flagged known-good as material). *The verdict is untrustworthy this run.* Report which calibration check failed, fix the reviewer prompt/agent, and re-run. Do **not** report PASS or FAIL.

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
