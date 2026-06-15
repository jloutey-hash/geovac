# Trunk — Definition of Done (pre-registered `/qa` criteria)

**Target:** the trunk = Papers **0, 1, 7** (group3 foundations) + **32, 38** (group1 operator-algebras) + the group3 foundations synthesis.

**This file is the pre-registration.** It is frozen *before* a `/qa trunk` run; the reviewers check against it and may not redefine it. Changing it is a deliberate PI act, dated and logged below — not something done mid-review to move the goalposts.

**Verdict rule.** The trunk **PASSES** only when a *calibrated* reviewer panel (see `/qa` step 5) returns **zero MATERIAL defects** against every criterion below, and independent reviewers converge. Any verified MATERIAL defect ⇒ **FAIL**. A panel that fails calibration ⇒ **INCONCLUSIVE** (not PASS).

## Material vs nit

A finding is **MATERIAL** iff fixing it would change a result's truth value, a claim's tier, a test's validity, a citation's referent, or a reader's takeaway — i.e. it passes the counterfactual *"would a trunk result change?"*. Everything else (wording, formatting, clarity, a missing comma, a stylistic preference) is a **NIT**: logged, never blocking.

## Criteria (each binary — holds / does not)

- **C1 — Backing exists & passes.** Every load-bearing trunk claim in `docs/claim_test_matrix.md` has a backing test (or an explicit "proof-by-argument" entry), and every such test currently PASSES when run.
- **C2 — Backing proves the claim.** Each backing test genuinely *proves* its claim — not tautological, not false-positive (right answer for the wrong reason / wrong evaluation space), not weaker than the prose. (`code-reviewer` certifies per claim.)
- **C3 — Prose ≤ tier, tier inline.** Every load-bearing claim's prose stays within its provenance tier (no overclaim), and the tier is stated *inline in the paper*. "Derived / SYMBOLIC PROOF" requires a derivation test, never a matching/convergence test. (`claims-reviewer`.)
- **C4 — Citations grounded.** Every external citation attached to a load-bearing trunk claim resolves to a real work that says what we attribute (no fabricated IDs, wrong venues, or nonexistent theorem numbers). (`citation-reviewer`.)
- **C5 — Hard prohibitions intact (§13.5).** K = π(B+F−Δ) is labeled an **Observation** everywhere it appears (never derived/conjecture/conjectural/theorem); no fitted or empirical parameter has been introduced; no negative result from §3/CHANGELOG is suppressed. Backed by the deterministic **K-label screen** `debug/qa/check_k_label.py` (run in step 1; FAILs on any SUSPECT K-tier assertion in the trunk scope) **and** the `claims-reviewer` enumerate-every-K-sentence pass (the screen is a backstop, not a replacement — the run-#4 lesson, where an LLM reviewer sampled and missed a planted "K is now derived as a theorem").
- **C6 — Discrete-vs-continuum precision.** No trunk paper claims the *discrete graph* "produces / has" the −(n²−1) spectrum. The discrete graph Laplacian is positive-semidefinite; −(n²−1) is the *continuum* (S³ Laplace–Beltrami) property the graph *converges* to.
- **C7 — WH1 stated accurately.** WH1 / Paper 38 is presented as PROVEN scoped to the **van Suijlekom state-space GH distance** (translation-seminorm metrization); no residual "Latrémolière propinquity" overclaim in the body where the proved object is the state-space GH distance.
- **C8 — Headline constants honest.** κ = −1/16 is an **Observation** (coincidence, no bridge), not a derivation; 4/π is **derived (numerics-pinned)**, not asserted as full symbolic proof; the Forced-Count moduli chain is stated at its true full-axiom count, not the matter-sector subcount.
- **C9 — Synthesis faithful.** The group3 foundations synthesis traces every claim to a paper that supports it, carries no descoped/withdrawn (zombie) result, and does not overstate convergence between papers or misstate a paper's status.
- **C10 — Compiles.** Each trunk paper compiles with ERRORS=0 (pre-existing non-blocking undefined-citation warnings may be noted but are not MATERIAL unless they break a load-bearing reference).
- **C11 — Internal-citation titles (deterministic).** Every internal GeoVac citation names the cited paper by its current `\title` — certified by the deterministic check `debug/qa/check_internal_titles.py` (exit 0; also `tests/test_internal_title_consistency.py`), **not** by an LLM reviewer. This is the `/qa` run-#1 lesson: title drift is a string-comparison problem, so a script catches it more reliably than a prose reviewer. The descope-pending propinquity cluster (Papers 39/40/45–49) is *flagged*, not failed.
- **C12 — K-label cleanliness (deterministic, the run-#4 finding).** No occurrence of the Paper 2 combination rule K = π(B+F−Δ) is labeled **conjectural / conjecture / derived / a theorem / proven** anywhere in the certified scope — i.e. the **2026-06-14 conjecture→Observation downgrade is fully applied** to the papers, not just to CLAUDE.md §13.5. Certified by `debug/qa/check_k_label.py` (exit 0; also `tests/test_k_label_check.py`), scoped to the target (default = trunk; `--gate <group>` per branch). **This is a backstop, not a replacement** for the `claims-reviewer` C5 pass (the screen can't anchor a bare "K" without the rule name in-clause). Root cause of the drift: the papers were written under the pre-downgrade §13.5 and lag the (correctly updated) rule — the branch sweep clears them (see `paper_57:446` "preserved as conjectural per §13.5", the self-undermining smoking gun). **Each branch's `<branch>.done.md` inherits C12 scoped to that branch** (`check_k_label.py --gate <branch-folder>`).

## Review dimensions (all mandatory in a single run)

The criteria above partition into independent **review dimensions**. A `/qa trunk` run must exercise *every* one in one invocation — passing one dimension does **not** excuse skipping another, and an unexercised gating dimension forces **INCONCLUSIVE**, never PASS (the run-#2 lesson: a "PASS on the exercised dimensions" hid the real κ/K defects that live only in the synthesis dimension).

| Dimension | Criteria gated | Reviewer | Notes |
|---|---|---|---|
| Code / test-backing | C1, C2 | `code-reviewer` (one per paper with tests) | RUN the tests; does each *prove* its claim? |
| Paper claims / prose | C3, C5, C6, C8 | `claims-reviewer` (one per paper) | prose ≤ tier; hard prohibitions intact |
| External citations | C4 | `citation-reviewer` (one per paper) | cites resolve and say what we attribute |
| Synthesis faithfulness | C9 | `claims-reviewer` (one on the group3 synthesis) | **separate dispatch** from the per-paper claims review |
| Deterministic | C10, C11, C12 | scripts (`check_internal_titles.py`, `check_k_label.py`, compile) | not an LLM reviewer |

The target verdict is the **AND across these dimensions**: PASS only if all are exercised, calibrated, and clean; FAIL if any calibrated dimension has a verified material defect; INCONCLUSIVE if any gating dimension is unexercised or uncalibrated.

## Coverage honesty

"PASS" means the trunk **survived the calibrated detectors for the defect classes in `docs/qa/seed_defects.md` and the criteria above, across all review dimensions** — not that it is provably perfect. A `/qa` run names any criterion or defect-class it could not exercise; an unexercised gating dimension is INCONCLUSIVE, not a footnoted PASS. When a new defect class is discovered (the way §3 dead-ends grow), it is added here and to the seed catalog, and the bar quietly rises.

## Change log
- 2026-06-14 — created (co-authored PM + PI) as the first pre-registered `/qa` target.
- 2026-06-14 — added the **Review dimensions** map + all-dimensions-mandatory rule after `/qa` run #3 found that runs #1–#2 had exercised only the claims + citation dimensions; the code (C1–C2) and synthesis (C9) dimensions, run for the first time in #3, surfaced 2 real synthesis defects (κ "derivable", K "conjectural"). The gate must attempt every dimension every run.
