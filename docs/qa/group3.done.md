# Group 3 (Foundations) — Definition of Done (pre-registered `/qa` criteria)

> **FROZEN 2026-06-15** (PI-confirmed). Pre-registration is complete; criteria are locked for the run and may not be redefined mid-review.
>
> **First-pass subset.** This branch is certified in bites. The 2026-06-15 run covers Papers **22, 24, 31** (the three non-trunk group3 KEYSTONEs) only; **C9 this pass is scoped to the synthesis's coverage of those three papers.** A full group3 branch PASS additionally requires Papers 18, 54, 55, 56, 57 and the full-synthesis C9 in later passes — the criteria below are the branch standard regardless and are not relaxed for the subset.

**Target:** the **group3 foundations branch**, scoped to the papers **not already certified in the trunk**:

- **In scope:** Papers **18, 22, 24, 31, 54, 55, 56, 57** (group3) + a **re-confirmation of the group3 foundations synthesis** (`papers/synthesis/group3_foundations_synthesis.tex`) against this fuller paper set.
- **Taken as already-certified (NOT re-litigated here):** the trunk papers **0, 1, 7** and the group3 synthesis's *trunk-relevant* claims, certified at `/qa trunk` PASS (v4.16.2). A group3 paper's claim is still in scope where it **depends on or restates** a trunk result (covered by C7).

**This file is the pre-registration.** Frozen *before* the run; reviewers check against it and may not redefine it. Changing it is a deliberate PI act, dated in the change log — never done mid-review to move goalposts.

**Verdict rule.** Group3 **PASSES** only when a *calibrated* reviewer panel returns **zero MATERIAL defects** against every criterion, across every gating dimension, and independent reviewers converge. Any verified MATERIAL defect ⇒ **FAIL**. Any gating dimension unexercised or uncalibrated ⇒ **INCONCLUSIVE** (not PASS).

## Material vs nit

A finding is **MATERIAL** iff fixing it would change a result's truth value, a claim's tier, a test's validity, a citation's referent, or a reader's takeaway — it passes the counterfactual *"would a group3 result change?"*. Everything else (wording, formatting, clarity) is a **NIT**: logged, never blocking.

## Criteria (each binary — holds / does not)

- **C1 — Backing exists & passes.** Every load-bearing group3 claim in `docs/claim_test_matrix.md` has a backing test (or an explicit "proof-by-argument" entry), and every such test PASSES when run.
- **C2 — Backing proves the claim.** Each backing test genuinely *proves* its claim — not tautological, not false-positive (right answer / wrong reason / wrong evaluation space), not weaker than the prose. (`code-reviewer` per paper with tests.)
- **C3 — Prose ≤ tier, tier inline.** Every load-bearing claim's prose stays within its provenance tier (no overclaim), and the tier is stated *inline in the paper*. "Derived / SYMBOLIC PROOF" requires a derivation test, never a matching/convergence test. (`claims-reviewer`.)
- **C4 — Citations grounded.** Every external citation on a load-bearing group3 claim resolves to a real work that says what we attribute (no fabricated IDs, wrong venues, nonexistent theorem numbers). (`citation-reviewer`.) Foundations-specific watch: the periods/Tannakian apparatus of Papers 55/56 (Deligne–Milne, Brown, Fathizadeh–Marcolli, mixed-Tate / cosmic-Galois literature) — high fabrication-risk surface.
- **C5 — Hard prohibitions intact (§13.5).** K = π(B+F−Δ) is labeled an **Observation** everywhere it appears in group3 scope (never derived/conjecture/conjectural/theorem) — group3 papers and the synthesis *discuss* K even though Paper 2 lives in group5; no fitted/empirical parameter introduced; no §3/CHANGELOG negative suppressed. Backed by the deterministic **K-label screen** `debug/qa/check_k_label.py --gate group3` **and** the `claims-reviewer` enumerate-every-K-sentence pass (screen is a backstop, not a replacement).
- **C6 — Discrete-vs-continuum precision.** No group3 paper claims the *discrete graph* "produces / has" the −(n²−1) spectrum. The discrete graph Laplacian is positive-semidefinite; −(n²−1) is the *continuum* (S³ Laplace–Beltrami) property the graph *converges* to. (Watch Papers 22, 24 spectrum statements.)
- **C7 — Trunk-dependent status accurate.** Where a group3 paper cites or restates a trunk result, it states that result's *current* status — Paper 38/WH1 as **PROVEN, scoped to the van Suijlekom state-space GH distance** (no residual "Latrémolière propinquity" overclaim); κ = −1/16 as an **Observation**; no group3 paper overstates a trunk keystone it depends on.
- **C8 — Group3 headline constants honest.** Each branch headline asserts no more than its tier:
  - Paper 22 angular sparsity = **potential-independent theorem** stated at its verified scope (e.g. 1.44% ERI density at l_max=3), not extrapolated beyond what is verified;
  - Paper 24 Bargmann–Segal lattice = **bit-exactly π-free in exact rational arithmetic** at finite cutoff (a structural/arithmetic fact, not a physics derivation); HO rigidity + Coulomb/HO six-layer asymmetry stated as structural results;
  - Paper 31 A/D universal/Coulomb partition stated at its actual scope (spectral-triple language);
  - Paper 55 every GeoVac period = **cyclotomic mixed-Tate at level ≤ 4** as a *classification* (tier-appropriate), not a derivation of physics;
  - Paper 57 forced/free P5 packing-reachability discriminator at its **true 98.3%** with the one ambiguous misclassification (I3) noted, not rounded to "100%".
- **C9 — Synthesis faithful (the synthesis-layer readiness gate).** The group3 foundations synthesis traces every claim to a paper that supports it, carries **no descoped/withdrawn (zombie) result**, does not overstate convergence between papers or misstate a paper's status, and faithfully reflects Papers 18/22/24/31/54–57. (`claims-reviewer`, **separate dispatch** from the per-paper claims review.)
- **C10 — Compiles.** Each in-scope group3 paper compiles with ERRORS=0 (pre-existing non-blocking undefined-citation warnings noted, not MATERIAL unless they break a load-bearing reference).
- **C11 — Internal-citation titles (deterministic).** Every internal GeoVac citation in group3 scope names the cited paper by its current `\title` — certified by `debug/qa/check_internal_titles.py` (exit 0; also `tests/test_internal_title_consistency.py`), not an LLM reviewer. The descope-pending propinquity cluster (Papers 39/40/45–49) is *flagged*, not failed.
- **C12 — K-label cleanliness (deterministic).** No occurrence of K = π(B+F−Δ) is labeled conjectural / conjecture / derived / theorem / proven anywhere in group3 scope — the 2026-06-14 conjecture→Observation downgrade is **fully applied** to group3 papers, not just CLAUDE.md §13.5. Certified by `debug/qa/check_k_label.py --gate group3` (exit 0; also `tests/test_k_label_check.py`). **Known target:** `paper_57:446` ("preserved as conjectural per §13.5") is the smoking-gun stale label this branch sweep must clear.
- **C13 — Paper↔test reference integrity (deterministic).** Every test cited *inline* in an in-scope group3 paper resolves to a **real, live** test under `tests/` — certified by `debug/qa/check_paper_test_refs.py --gate group3` (exit 0; also `tests/test_paper_test_refs.py`). Existence is the gate; matrix↔paper coverage and non-scope stale refs are advisory.

## Review dimensions (all mandatory in a single run)

A `/qa group3` run must exercise *every* dimension in one invocation; an unexercised gating dimension forces **INCONCLUSIVE**, never PASS. Verdict = **AND across dimensions**.

| Dimension | Criteria gated | Reviewer | Notes |
|---|---|---|---|
| Code / test-backing | C1, C2 | `code-reviewer` (one per in-scope paper with tests) | RUN the tests; does each *prove* its claim? |
| Paper claims / prose | C3, C5, C6, C7, C8 | `claims-reviewer` (one per in-scope paper) | prose ≤ tier; hard prohibitions intact |
| External citations | C4 | `citation-reviewer` (one per in-scope paper) | cites resolve and say what we attribute |
| Synthesis faithfulness | C9 | `claims-reviewer` (one on the group3 synthesis) | **separate dispatch** from per-paper claims |
| Deterministic | C10, C11, C12, C13 | scripts (`check_internal_titles.py`, `check_k_label.py --gate group3`, `check_paper_test_refs.py --gate group3`, compile) | not an LLM reviewer |

## Coverage honesty

"PASS" means group3 **survived the calibrated detectors for the defect classes in `docs/qa/seed_defects.md` and the criteria above, across all review dimensions** — not that it is provably perfect. A run names any criterion or defect-class it could not exercise; an unexercised gating dimension is INCONCLUSIVE, not a footnoted PASS.

## Change log
- 2026-06-15 — **DRAFTED** by PM for PI review (second pre-registered `/qa` target, after `trunk`). Scope excludes trunk-certified Papers 0/1/7; includes group3 synthesis re-check as the synthesis-layer readiness gate.
- 2026-06-15 — **FROZEN** (PI-confirmed). First pass certifies Papers **22, 24, 31** only ("start with 3, first pass"); full-branch criteria unchanged, remaining papers (18, 54–57) + full synthesis deferred to later bites. C9 this pass scoped to the synthesis's coverage of 22/24/31.
