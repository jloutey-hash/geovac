# Branch-QA Definition of Done — shared criteria (the standard)

> **Single source of truth for the `/qa` certification criteria.** Each branch's
> `docs/qa/<target>.done.md` is now a thin **profile** that *inherits* these
> criteria and supplies only what is genuinely branch-specific: scope (papers
> in / excluded), per-criterion watch-notes, any branch-specific extra criteria
> (C14+), the deterministic-check `--gate <branch>` arg, the first-bite subset,
> and the freeze status + change log. **Change a criterion HERE once and every
> branch inherits it** — the criteria are no longer copy-pasted per branch.
> This mirrors what the deterministic checks already do (one script +
> `--gate <branch>`) and the seed catalog (one `seed_defects.md`).

> **Pre-registration / freeze is preserved.** A branch is *frozen* when its
> profile is PI-confirmed and marked FROZEN. Git pins what *this* file said at
> that commit, so "criteria frozen before the run" holds by the commit history.
> Reviewers check against {this file} + {the branch profile} and may not
> redefine either mid-review. Editing this file is a deliberate PI act, dated in
> the change log — never done mid-review to move goalposts.

## Material vs nit

A finding is **MATERIAL** iff fixing it would change a result's truth value, a
claim's tier, a test's validity, a citation's referent, or a reader's takeaway
— it passes the counterfactual *"would a result in this branch change?"*.
Everything else (wording, formatting, clarity) is a **NIT**: logged, never
blocking.

## Verdict rule

A branch **PASSES** only when a *calibrated* reviewer panel returns **zero
MATERIAL defects** against every criterion, across every gating dimension, and
independent reviewers converge. Any verified MATERIAL defect ⇒ **FAIL**. Any
gating dimension unexercised or uncalibrated ⇒ **INCONCLUSIVE** (not PASS). The
verdict is the **AND across all review dimensions** (below).

## Criteria (C1–C14; each binary — holds / does not)

- **C1 — Backing exists & passes.** Every load-bearing claim in scope listed in
  `docs/claim_test_matrix.md` has a backing test (or an explicit
  "proof-by-argument" entry), and every such test PASSES when run.
- **C2 — Backing proves the claim.** Each backing test genuinely *proves* its
  claim — not tautological, not false-positive (right answer / wrong reason /
  wrong evaluation space), not weaker than the prose. (`code-reviewer` per paper
  with tests; RUN the tests.)
- **C3 — Prose ≤ tier, tier inline.** Every load-bearing claim's prose stays
  within its provenance tier (no overclaim), and the tier is stated *inline in
  the paper*. "Derived / SYMBOLIC PROOF" requires a derivation test, never a
  matching/convergence test. (`claims-reviewer`.)
- **C4 — Citations grounded.** Every external citation on a load-bearing claim
  resolves to a real work that says what we attribute (no fabricated IDs, wrong
  venues, nonexistent theorem/def numbers). (`citation-reviewer`.) *A branch
  profile may flag a high-fabrication surface to prioritize.*
- **C5 — Hard prohibitions intact (§13.5).** K = π(B+F−Δ) is labeled an
  **Observation** everywhere it appears in scope (never derived / conjecture /
  conjectural / theorem); no fitted/empirical parameter introduced; no
  §3/CHANGELOG negative suppressed. Backed by the deterministic K-label screen
  `debug/qa/check_k_label.py --gate <branch>` **and** the `claims-reviewer`
  enumerate-every-K-sentence pass (the screen is a backstop, not a replacement —
  it can't anchor a bare "K" without the rule name in-clause).
- **C6 — Discrete-vs-continuum precision.** No paper claims the *discrete graph*
  "produces / has" the −(n²−1) spectrum (or any continuum spectrum as a bare
  graph property). The discrete object is positive-semidefinite / a truncation;
  the continuum spectrum is what it *converges* to.
- **C7 — Trunk-dependent status accurate.** Where an in-scope paper cites or
  restates a trunk result, it states that result's *current* status. *The branch
  profile names which trunk results it depends on and their current tier* (e.g.
  Paper 38/WH1 = PROVEN scoped to the van Suijlekom state-space GH distance, no
  residual "Latrémolière propinquity" overclaim; κ = −1/16 = Observation).
- **C8 — Branch headline results honest.** Each branch headline asserts no more
  than its tier. *The branch profile enumerates the per-paper headlines and their
  tier-appropriate statements.*
- **C9 — Synthesis faithful (the synthesis-layer readiness gate).** The branch's
  group synthesis traces every claim to a paper that supports it, carries **no
  descoped/withdrawn (zombie) result**, does not overstate convergence between
  papers or misstate a paper's status, and faithfully reflects the in-scope
  papers. (`claims-reviewer`, **separate dispatch** from the per-paper claims
  review.)
- **C10 — Compiles.** Each in-scope paper compiles with ERRORS=0 (pre-existing
  non-blocking undefined-citation warnings noted, not MATERIAL unless they break
  a load-bearing reference).
- **C11 — Internal-citation titles (deterministic).** Every internal GeoVac
  citation in scope names the cited paper by its current `\title` — certified by
  `debug/qa/check_internal_titles.py` (exit 0; also
  `tests/test_internal_title_consistency.py`), not an LLM reviewer. The
  descope-pending propinquity cluster (Papers 39/40/45–49) is *flagged*, not
  failed.
- **C12 — K-label cleanliness (deterministic).** No occurrence of K = π(B+F−Δ)
  is labeled conjectural / conjecture / derived / theorem / proven anywhere in
  scope. Certified by `debug/qa/check_k_label.py --gate <branch>` (exit 0; also
  `tests/test_k_label_check.py`). Backstop to C5, not a replacement.
- **C13 — Paper↔test reference integrity (deterministic).** Every test cited
  *inline* in an in-scope paper resolves to a **real, live** test under `tests/`
  — certified by `debug/qa/check_paper_test_refs.py --gate <branch>` (exit 0;
  also `tests/test_paper_test_refs.py`). Existence is the gate; matrix↔paper
  coverage and non-scope stale refs are advisory.
- **C14 — Paper↔file reference integrity (deterministic).** Every `geovac/`,
  `benchmarks/`, or `demo/` (permanent code/artifact) file path cited *inline*
  in an in-scope paper resolves to a **real** file — certified by
  `debug/qa/check_file_refs.py --gate <branch>` (exit 0; also
  `tests/test_file_ref_check.py`). Closes the C13 gap (C13 covers `tests/` refs
  only): the group3 run-4 `geovac/jlo_chi.py` defect was a claim citing a
  nonexistent code module. `debug/` refs are **advisory** — the transient
  clean-room dir (CLAUDE.md §9) is pruned by design, so dangling `debug/`
  pointers are a hygiene smell, not a cert blocker. Per the §9 policy (papers
  cite the permanent record — CHANGELOG / the paper / tests — never transient
  `debug/`), the existing ~443 dangling `debug/` refs (papers 34/32/28 …) are a
  **deferred corpus-wide sweep**, adopted-as-policy 2026-06-17, ref-removal TBD.

## Branch-specific criteria (C14+)

A branch may add criteria numbered **C14+** *in its profile* when it carries a
risk class the shared set doesn't cover (e.g. group1 **C14 — descope/partial
status accuracy**, for a branch with DESCOPED/PARTIAL papers). If a second
branch needs the same criterion, **promote it here** (the single-source rule:
don't copy a criterion into two profiles).

## Review dimensions (all mandatory in a single run)

A `/qa <branch>` run must exercise *every* dimension in one invocation; an
unexercised gating dimension forces **INCONCLUSIVE**, never PASS. Verdict = AND
across dimensions. A clean dimension never short-circuits another; the run does
not stop early on a pass.

| Dimension | Criteria gated | Reviewer | Notes |
|---|---|---|---|
| Code / test-backing | C1, C2 | `code-reviewer` (one per in-scope paper with tests) | RUN the tests; does each *prove* its claim? |
| Paper claims / prose | C3, C5, C6, C7, C8 (+ branch C14+) | `claims-reviewer` (one per in-scope paper) | prose ≤ tier; hard prohibitions intact |
| External citations | C4 | `citation-reviewer` (one per in-scope paper) | cites resolve and say what we attribute |
| Synthesis faithfulness | C9 | `claims-reviewer` (one on the branch synthesis) | **separate dispatch** from per-paper claims |
| Deterministic | C10, C11, C12, C13, C14 | scripts (`check_internal_titles.py`, `check_k_label.py --gate <branch>`, `check_paper_test_refs.py --gate <branch>`, `check_file_refs.py --gate <branch>`, compile) | not an LLM reviewer |

## Coverage honesty

"PASS" means the branch **survived the calibrated detectors for the defect
classes in `docs/qa/seed_defects.md` and the criteria above, across all review
dimensions** — not that it is provably perfect. A run names any criterion or
defect-class it could not exercise; an unexercised gating dimension is
INCONCLUSIVE, not a footnoted PASS. When a new defect class is discovered (the
way §3 dead-ends grow), it is added here and to the seed catalog, and the bar
quietly rises.

## Hard rules

- **PI-invoked only.** The PM never self-triggers `/qa`, runs it proactively, or
  suggests it each sprint.
- Pre-registered criteria (this file + the branch profile) are **frozen before**
  the review — no goalpost-moving in either direction.
- Seeds live only in the throwaway worktree; never commit/leak them; always
  remove the worktree.
- A reviewer's verdict is trusted only **after** it passes calibration (caught
  the plants, passed the controls). An uncalibrated panel ⇒ INCONCLUSIVE, never
  PASS.
- **All dimensions every run.** The verdict is the AND across all dimensions; an
  unexercised gating dimension forces INCONCLUSIVE, not PASS.

## Change log
- 2026-06-16 — **Extracted** from `trunk.done.md` / `group3.done.md` (which were
  ~90% identical) into this single shared standard, per PI direction
  ("consolidate to a single definition of done"). The three branch files become
  thin profiles. No criterion was changed in the extraction — C1–C13 are
  verbatim the trunk/group3 standard, generalized (branch-specific watch-notes
  moved into the profiles).
