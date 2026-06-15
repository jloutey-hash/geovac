# Seed-defect catalog (`/qa` per-run calibration controls)

The defect classes `/qa` plants — blind, on a throwaway worktree, never on the real corpus — to **measure whether this run's reviewer panel can actually detect dirt** before its verdict is trusted. Each class is drawn from GeoVac's *real* §3 error history, so catching it is genuinely informative about catching the real thing. Every class is **MATERIAL** by construction (these are exactly the things that must be caught).

**Per run:** plant a *random subset* of K≈4–6 classes (vary which, so reviewers can't pattern-learn the set), at realistic locations, blind. Also designate M≈4–6 **known-good controls** — verified-correct claims that must NOT be flagged material. Record the answer key in `debug/qa/<target>_seed_key.json`.

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
| S7 | **Stale internal-title drift** | cite a GeoVac paper by an out-of-date title (e.g. Paper 0 as "Quantum State Space as a Packing Problem") | `claims-reviewer` / `code-reviewer` |
| S8 | **Zombie / descoped claim in a synthesis** | have the synthesis assert a withdrawn result (e.g. Paper 45's retracted K⁺ propinquity theorem) as standing | `claims-reviewer` |
| S9 | **Status overstatement of a keystone** | upgrade a CONDITIONAL/PANEL-VERIFIED result to "PROVEN/unconditional" in the prose | `claims-reviewer` |

## Notes
- Seeds must be **realistic** (indistinguishable from a genuine slip in context) — not obviously-planted markers a reviewer can dismiss as a test.
- Seeds and the answer key live **only** in the worktree + `debug/qa/`. They are never committed to the real corpus; the worktree is removed at the end of every run.
- This catalog **grows**: when a genuinely new defect class slips past `/qa` and is caught some other way, add it here (and the relevant criterion to the target's `.done.md`).
