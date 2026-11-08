# Sprint: the `/qa` certification gate — build + inaugural run (2026-06-14)

## Problem (PI-posed)
A reviewer tasked with "find problems" always finds problems → a "clean" can't be *earned*, a "dirty" can't be *trusted*. A fault-hunt is an activity, not a test (no null hypothesis). Completion-bias (the doer stretches to report done) and fault-bias (the finder always finds) are mirror failures; neither discriminates "done" from "not done." The ask: a procedure that can genuinely PASS or FAIL **and discriminate between the two**.

## Solution — four mechanisms + a three-way verdict
1. **Pre-registered criteria** — binary, material, frozen *before* the review (`docs/qa/<target>.done.md`). Turns an unbounded fault-hunt into a finite checklist with two outcomes per item; the reviewer can't redefine "done."
2. **Materiality grading** — clean = zero MATERIAL, not zero findings. Counterfactual per finding: "would a target result change?" Makes "clean" achievable despite infinite nits.
3. **Per-run calibration controls** — seed known defects *blind* into an isolated git worktree (`docs/qa/seed_defects.md`); score whether THIS run's reviewers catch them. **This is what makes a verdict trustworthy** — a "clean" is only meaningful from a detector you've just shown has detection power.
4. **Convergence** across independent fresh reviewers.

**Three-way verdict:** PASS (calibrated panel + zero material) / FAIL (calibrated + ≥1 verified material) / **INCONCLUSIVE** (panel missed plants or false-positived controls → untrustworthy, fix reviewer + re-run). The third outcome distinguishes *"target not done"* from *"reviewer not trustworthy."*

## Artifacts
- `.claude/commands/qa.md` — protocol skill. **PI-invoked only** (hard rule: PM never self-triggers). `/qa <target>`.
- `.claude/agents/claims-reviewer.md` — prose/synthesis overclaim reviewer (companion to `code-reviewer` + `citation-reviewer`).
- `docs/qa/trunk.done.md` — trunk's pre-registered 10-criterion definition-of-done.
- `docs/qa/seed_defects.md` — 9 seed-defect classes from §3 history.
- CLAUDE.md §13.9b — `/qa` registered (genuine trigger, **no memory rule** ⇒ PI controls timing).

## Inaugural run — `/qa trunk`, run #1 → **INCONCLUSIVE**
Representative panel (prose + external-citation dimensions): `claims-reviewer` ×3 (P1/P7/P32) + `citation-reviewer` ×1 (P38). 5 blind seeds + 4 known-good controls planted in a worktree (discarded after).

- **Sensitivity 4/5.** Caught: S1 fabricated arXiv ID (citation-reviewer even resolved the wrong ID to an unrelated bee-health ML paper); S4 κ "coincides"→"rigorously derived"; S5 K=π(B+F−Δ) "Observation"→"derived as a theorem" (hard-prohibition); S6 "discrete graph **directly produces** −(n²−1)". **Missed:** S7 stale internal title (Paper 0 cited by old title).
- **Specificity 4/4.** Zero false positives — every known-good control (c²/40, [T₊,L₊]=0, axiom audit, the 21 real P38 citations) correctly passed.
- **Verdict INCONCLUSIVE** — the missed planted class means the panel is not fully calibrated, so the gate returned INCONCLUSIVE *instead of a false PASS*. **The gate caught its own blind spot.** This is the intended discrimination; per-run seeding is what surfaced it (meta-validation alone never would have).
- **Real (trustworthy, given 4/4 specificity) trunk finding:** Paper 32's H1 Yukawa theorem states "128" as the *full* forced count with no matter-sector relabel — the parallel Forced-Count theorem got the v4.13.0 correction (full-axiom = 260), the H1 Yukawa twin did not. Candidate FAIL. Plus small items (P38 Perez-Sanchez stance-misattribution + `latremoliere2018` inline vol/yr slip; P1 figure captions still print stale numbers).

## Lessons (run #1 was a test of the gate itself)
1. **Criteria file wasn't in the worktree** (uncommitted; worktree = last commit). Reviewers flagged it, fell back to paraphrased criteria + the matrix/register (which ARE in the worktree). FIX: commit the qa files (or copy them in) before a run.
2. **claims-reviewer is unreliable on internal-title accuracy** (spot-checked P2's title, missed the seeded P0 title). INSIGHT: some criteria (internal-title consistency) want a **deterministic check** (grep), not an LLM reviewer — narrow LLM reviewers to judgment calls.
3. **code-reviewer / C1–C2 (test-backing) not exercised** this run (representative panel) — honest-ceiling.

## Before re-run
Commit `/qa` tooling → move internal-title check to deterministic (or add to code-reviewer) → fix H1 Yukawa "128" + small P38 items → re-run `/qa trunk`.
