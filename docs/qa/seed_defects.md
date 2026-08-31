# Seed-defect catalog (`/qa` per-run calibration controls)

The defect classes `/qa` plants — blind, on a throwaway worktree, never on the real corpus — to **measure whether this run's reviewer panel can actually detect dirt** before its verdict is trusted. Each class is drawn from GeoVac's *real* §3 error history, so catching it is genuinely informative about catching the real thing. Every class is **MATERIAL** by construction (these are exactly the things that must be caught).

**Per run:** plant a *random subset* of K≈4–6 classes (vary which, so reviewers can't pattern-learn the set), at realistic locations, blind. **Cover every gating dimension:** the subset must include at least one class catchable by each reviewer type in play this run — `code-reviewer` (S2/S3), `claims-reviewer` paper-prose (S4/S5/S6/S9), `claims-reviewer` synthesis (S8/S9), `citation-reviewer` (S1), deterministic (S7) — so no dimension that gates the verdict is left uncalibrated (the run-#3 lesson: the code and synthesis dimensions had never been seeded until they were run together). Also designate M≈4–6 **known-good controls** — verified-correct claims that must NOT be flagged material. Record the answer key in `debug/qa/<target>_seed_key.json`.

**Scoring:** sensitivity = planted classes caught / planted (a class no reviewer caught ⇒ panel blind to it ⇒ INCONCLUSIVE). specificity = 1 − (known-good controls flagged material / controls).

## Classes

| # | Defect class (historical source) | How to plant it | Expected catcher |
|---|---|---|---|
| S1 | **Fabricated / wrong external citation** (Fursaev–Solodukhin: arXiv ID resolving to a different paper) | swap a real bibitem's arXiv ID/venue for a wrong one, or invent a plausible-but-nonexistent theorem number | `citation-reviewer` |
| S2 | **Circular / tautological test** (the 4/π hardcoded-constant test) | add a test that asserts a hardcoded module constant against itself, dressed as a verification | `code-reviewer` |
| S3 | **False-positive test — passes for the wrong reason** (TC qubit-space diagonalization) | back a claim with a test that computes in the wrong space / under a cancellation that masks the real quantity | `code-reviewer` |
| S4 | **κ-style overclaim — "derived" where it only coincides** | change a paper line so a *matched* value reads "derived from …" / "SYMBOLIC PROOF" | `claims-reviewer` |
| S5 | **Hard-prohibition touch (§13.5)** | relabel K = π(B+F−Δ) as "conjecture"/"derived", or drop its Observation tag | `claims-reviewer` (tripwire) |
| S6 | **Discrete-vs-continuum imprecision** | insert "the discrete graph *produces* the −(n²−1) spectrum" | `claims-reviewer` |
| S7 | **Stale internal-title drift** | cite a GeoVac paper by an out-of-date title (e.g. Paper 0 as "Quantum State Space as a Packing Problem") | **deterministic** — `debug/qa/check_internal_titles.py` (NOT an LLM reviewer; the run-#1 lesson) |
| S8 | **Zombie / descoped claim in a synthesis** | have the synthesis assert a withdrawn result (e.g. Paper 45's retracted K⁺ propinquity theorem) as standing | `claims-reviewer` |
| S9 | **Status overstatement of a keystone** | upgrade a CONDITIONAL/PANEL-VERIFIED result to "PROVEN/unconditional" in the prose | `claims-reviewer` |

## Notes
- Seeds must be **realistic** (indistinguishable from a genuine slip in context) — not obviously-planted markers a reviewer can dismiss as a test.
- **Innocuous seed-commit message (2026-07-03 lesson, group5 delta-2):** seeds are committed onto the worktree branch (run-8 hardening), so a code reviewer running `git log` SEES the commit subject. A message like "delta calibration controls" reveals the top commit is the seed layer; use a plausible work message instead (e.g. "Tier-2 polish follow-ups"). The delta-2 code agents quoted the revealing message yet still analyzed each seed as a genuine defect, so calibration held — but don't rely on that.
- Seeds and the answer key live **only** in the worktree + `debug/qa/`. They are never committed to the real corpus; the worktree is removed at the end of every run.
- This catalog **grows**: when a genuinely new defect class slips past `/qa` and is caught some other way, add it here (and the relevant criterion to the target's `.done.md`).

## Seed-placement rule (hard, added 2026-08-28 after a two-seed failure)

**Never plant a seed on a LaTeX comment line, and assert it after planting.**

In the group6 FULL run, 2 of 9 seeds (S3 synthesis tier-promotion, S4 author-name
corruption) landed on `%`-prefixed lines — S3 on line 2 of the synthesis file
header. Neither is visible in the rendered document, so a reviewer that does not
flag them is behaving *correctly*. Consequences that run:

- the **synthesis dimension ran uncalibrated** (its only seed was void) and could
  not be certified regardless of what it found;
- the **citation dimension** fell to 1 valid seed where its Sonnet tier requires 2.

Root cause: the seeding helper selected "the first line containing X" with no
comment-line filter, and nothing checked the result.

**Required of any seeding helper:**

1. Skip lines whose stripped form starts with `%` when choosing a target line.
2. After planting, **assert** the seeded line is not a comment and, for prose
   seeds, that the changed text appears in the compiled PDF's text layer (or at
   minimum outside any `%`-comment and outside `\begin{comment}` blocks).
3. Report a per-seed placement check in the answer key, so a void seed is caught
   at planting time rather than discovered while scoring calibration.

A void seed is worse than no seed: it silently converts a calibrated dimension
into an uncalibrated one while the run still *looks* fully seeded.
