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

**Code-docstring stale-prose is a NIT class (2026-06-24 PI direction, set at group1
certification).** Stale/withdrawn prose living in a `geovac/*.py` or `tests/*.py`
**docstring or comment** — e.g. a retracted convergence statement left as live prose,
or a result mislabeled "Latrémolière propinquity" where it is van Suijlekom state-space
GH — is a **fix-on-sight NIT, not a MATERIAL cert blocker**, *provided* the authoritative
source (the paper), the claim's tier, and the test's *logic/tolerances* are all correct.
Rationale: the paper — not the code docstring — is the authoritative source (§1); such
prose changes no result, claim, tier, or reader-of-the-paper takeaway; it is a long tail
(every backing-module docstring is a potential site) that fresh-adversary passes surface
one-scope-over indefinitely; and it is cheaply caught on sight by `code-reviewer` + the
**C16 docstring gate** (the propinquity / withdrawn-c3op entries now scan the backing
code modules, advisory severity). A code-docstring zombie that has drifted into a *paper*,
or that corrupts a test's *logic* (wrong tolerance, tautology, false-positive), remains
MATERIAL. This refines (does not relax) the §"Material vs nit" counterfactual: a code
docstring is not "a result, claim, tier, citation, or paper-reader takeaway."

**Secondary stale-number / provenance-attribution is a NIT class (2026-06-28 PI direction,
set at group2 certification).** A **stale or imprecise computed number on a NON-headline
secondary result** (one not listed in the branch's `<branch>.done.md` §C8 headline set —
e.g. an illustrative HeH²⁺ energy in a paper whose headline is H₂⁺), **or a provenance /
causal-attribution sentence on an otherwise-correct headline** (the headline number is
right; only the stated *reason* it changed is wrong), is a **fix-on-sight NIT, not a MATERIAL
cert blocker**, *provided* every DoD §C8 headline, its tier, and its backing test are correct.
Rationale: it passes the §"Material vs nit" counterfactual *negatively* — no DoD-listed result,
tier, or paper-reader takeaway-of-a-headline changes; and it is the long-tail class a
fresh-adversary pass surfaces ~1–2-per-run indefinitely on a large corpus (the group2 run-#4/#5
asymptote: perfect 11/11+6/6 calibration both runs, yet 2 thin DIFFERENT secondary/provenance
items each, never a headline). **A stale number or attribution that touches a DoD §C8 headline
— or that makes a headline number itself wrong — remains MATERIAL.** This refines (does not
relax) the counterfactual: a secondary number and a provenance note are not a headline result,
its tier, or its reader takeaway. (Fix them on sight regardless; they just do not gate the cert.)

**Withdrawal-chronicle placement is a NIT class (2026-07-04 PI direction, set at the
synthesis retraction-language pass).** The **chronicle** of a claim withdrawal — dates,
process narrative ("an adversarial verification pass retracted…", "this synthesis has been
brought into line…"), event-framed references ("descoped with the Paper 45 retraction") —
belongs to the **object of record** (the paper's Status note / History remark,
`docs/claims_register.md`, CHANGELOG.md), NOT to the synthesis layer. A synthesis (field
guide or group synthesis) carries **current state only**: what stands, what is
withdrawn/descoped and the *mathematical reason* (e.g. "descoped by Paper 45's degeneracy
theorem"), plus at most a one-clause flag with pointer for a headline-grade withdrawal
(stale-DOI reader protection — earlier Zenodo releases circulate the pre-withdrawal claim).
Vocabulary: **claims are "withdrawn" or "descoped"; reserve "retraction" for whole-paper
retraction, which has never occurred in the corpus.** Chronicle language found in a
synthesis is a **fix-on-sight NIT, not a MATERIAL cert blocker**, *provided* the
current-state assertion is accurate. **A synthesis that asserts a withdrawn claim as live
remains MATERIAL (the zombie rule), unchanged.** A negative *result* (e.g. a degeneracy
theorem) is content, not withdrawal language — it stays. Convention codified in
`docs/authoring_conventions.md` project-wide rule 11.

## Dual-rule ERI framing (A = sparsity, B = accuracy) — the framing-zombie rule

**The principle (2026-06-28 PI direction).** The codebase *intentionally* carries
**two electron-repulsion-integral (ERI) angular selection rules**, for two different
purposes, and a claim is honest only when it uses — and names — the rule that matches its
stated purpose:

- **Rule A — pair-diagonal** (`composed_qubit._ck_coefficient` and
  `lattice_index._ck_coefficient`, both `q = mc − ma` ⇒ nonzero only when `m_a=m_c ∧
  m_b=m_d`). This is the **sparsity / QC-product rule**: it drops genuinely-nonzero
  *m-swap* ERIs to yield a sparser qubit Hamiltonian. It is the rule the **shipped QC
  product** uses everywhere (`ecosystem_export.hamiltonian()` → composed *and* atomic). It
  is a **legitimate quality QC approximation** (energy effect ~mHa, sub-chemical) — *when
  the point is sparsity*.
- **Rule B — global-M_L** (`casimir_ci._gaunt_ck`, `q = ma − mc` ⇒ the exact Coulomb
  rule `m_a+m_b=m_c+m_d`). This is the **physics-accuracy rule**, used only in the
  precision-physics paths (`cusp_correction`, `dirac_ci`, `internal_multifocal`) — *when
  the point is accuracy*. It is **not** in the QC product.

This split is correct by design: **A where sparsity is the point, B where accuracy is the
point.** (The full study is `debug/sprint_group4_prework_memo.md` + the CF-1 carryforward.)

**A framing zombie** is a claim that breaks the A/B match:
1. a **pair-diagonal (A) resource/sparsity number presented as the *exact* / full
   Gaunt-selection-rule value** without disclosing the pair-diagonal approximation — e.g.
   the "2.7× fewer Pauli than STO-3G" market test (under B it is **parity**, 838 vs 907),
   the "d-orbital blocks are *sparser* (9.23<11.10)" claim (under B the d-block is **denser**,
   30.0>27.9), or the realized ERI density (D_pd≈1.44%) presented as *the* Gaunt
   selection-rule density (the universal D≈6.06%);
2. a **physics-accuracy claim that silently rests on the pair-diagonal (A) approximation**
   (e.g. an "exact integrals" energy framing where the angular selection is actually A);
3. the **wrong density/count for the stated context** (D=6.06% used where the realized
   product is described, or D_pd=1.44% used where the universal selection-rule density is meant).

**Enforcement.** The `claims-reviewer` (C3/C8) MUST **enumerate every sparsity / ERI-density
/ Pauli-count / "Gaunt selection rule" claim** in scope and verify each names its rule (A
pair-diagonal *disclosed as an approximation* for sparsity; B for accuracy). The deterministic
**C16** entry `pair-diagonal-as-exact-sparsity` backstops the *known* retired framings (the
market-test line, the d-block-sparser claim) — they may not re-surface without a
pair-diagonal/approximation disclosure within ±WINDOW lines. Disclosing A also means stating
that **B was considered and deliberately left on the table** (sparsity is the chosen purpose
for this corner of the corpus). This is a sharpening of C3 (prose ≤ tier) + C8 (headline
honesty) + §1.5 (benchmarking), shared across branches that make sparsity claims (group4
primary; group3 Papers 22/31 already carry the D-vs-D_pd disclosure).

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
- **C15 — Inline arXiv-ID consistency (deterministic).** No arXiv ID written
  *inline in the body prose* is a near-match transposition (same length, Hamming
  ≤ 2) of a bibitem ID without equalling it — certified by
  `debug/qa/check_inline_arxiv.py --gate <branch>` (exit 0; also
  `tests/test_inline_arxiv_check.py`). Closes the /qa group1 re-cert-1
  calibration gap:\ an inline `arXiv:2504.10830` seed (a transposition of the
  bibitem's `2504.10380`) slipped past the bibitem-focused LLM citation
  reviewers. C15 is a **complement** to the C4 LLM citation review, not a
  replacement:\ it catches transposed/typo'd IDs deterministically and cheaply;
  the reviewer still owns wrong-title / wrong-venue / wrong-attribution (right
  ID, wrong metadata — which C15 cannot see). Added 2026-06-19.
- **C16 — Retracted-claims / zombie-drift (deterministic).** No claim a prior
  /qa run retired re-surfaces *as live* — i.e. without an inline withdrawal flag
  — anywhere in the gated papers or their contained code/test modules. Certified
  by `debug/qa/check_retracted_terms.py --gate <branch>` (exit 0), whose
  `REGISTRY` holds each retired phrase + its exemption marker (the withdrawn
  Pythagorean C₃<1 form, the S⁷ "structural negative" vs its erratum, the
  Batch-2 false-closure framing, propinquity-as-achieved-metric, …).
  **Added 2026-06-23 (the lesson):** every recurring miss across the
  Batch-2/Batch-3 cert arcs was a *mechanical* drift in a low-salience region
  (status-table cell, docstring, footnote, one abstract clause) that judgment-mode
  LLM review walked past — the withdrawn Pythagorean form survived *three*
  incomplete sweeps. A grep is exhaustive and zero-variance; this moves the
  recurring zombie classes out of expensive judgment review and into a cheap
  gate. **Maintenance rule:** when a run retires/withdraws a claim, ADD its
  phrase to the registry (the [[feedback_deferral_is_churn]] single-source
  doctrine). C16 is a **complement** to the C3/C7/C9 LLM review (it catches
  *known* retracted phrases; the reviewers still own *new* overclaims and
  judgment calls). `--gate <branch>` scans the gated papers + their code modules.
- **C17 — Headline-number registry (deterministic).** Every occurrence of a
  DoD §C8 headline-number *family* carries its canonical value — certified by
  `debug/qa/check_headline_numbers.py --gate <branch>` (exit 0; also
  `tests/test_headline_numbers_check.py`, which self-tests each family against
  synthetic wrong values). **Added 2026-07-02 (PI direction; the 7th-group4-cert
  meta-lesson):** the genuine-MATERIAL classes that kept surviving calibrated
  judgment review had become *mechanical* — (a) second-locus propagation (a
  decided headline corrected at one locus, stale at another: Z=1–36 vs 56 at
  five loci; a memo-listed KH fix never applied) and (b) number-vs-source
  drift (a "190×" floor against a cited table whose own floor is 51×; a stale
  33.3 1-norm vs the live 32.6) — classes an LLM panel finds one locus at a
  time at full price, and a registry finds exhaustively for ~0. Families hold
  either known-wrong variant patterns (C16-style) or a capture+canonical pair,
  with require-context windows (so different quantities sharing a numeral,
  e.g. periodic-row "Z=1–10" prose or the ℓ-parity 0.97 Pauli-ratio cells,
  never false-positive) and exemption windows for historical mentions.
  **Maintenance rule (mirrors C16): when a run corrects or demotes a headline
  number, ADD/UPDATE its family here.**

- **C18 — Project-course duration language (deterministic).** No wall-clock
  duration (years / months / weeks / days) is attached to the project's own
  course or events — certified by `debug/qa/check_duration_language.py`
  (exit 0; `--selftest` built in). Historical-duration classes ("N years ago",
  "units of work", "across N days", "took N months") **FAIL**; forward
  effort-estimate vocabulary ("multi-month", "multi-year", "week-scale") is
  **ADVISORY standing debt** (~190 pre-existing instances in certified group
  papers) until the corpus-wide sweep retires it, then promoted to FAIL.
  Sequence language, real dates, version anchors, and unit-free effort
  vocabulary ("sprint-scale", "deep wall") are fine; external-world history is
  exempt. **Added 2026-07-04 (PI direction; the incident: the field guide
  claimed "Three years ago the project…" when the whole project is far
  younger — LLM drafting is unreliable about elapsed project time, so this is
  a grep-class check, not a judgment call).** Convention:
  `docs/authoring_conventions.md` rule 12.

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

Reviewer count per dimension follows the **granularity & scaling rule** in `qa.md` step 4 (set agents-per-dimension by depth-need, not uniformly) — the "Reviewer" column below gives the *type* and its native granularity, not a fixed per-paper count:

| Dimension | Criteria gated | Reviewer (type · granularity) | Notes |
|---|---|---|---|
| Code / test-backing | C1, C2 | `code-reviewer` · **1 paper/agent** (does not consolidate) | RUN the tests; does each *prove* its claim? |
| Paper claims / prose | C3, C5, C6, C7, C8 (+ branch C14+) | `claims-reviewer` · **chunk ~3–4 papers/agent** | prose ≤ tier; hard prohibitions intact |
| External citations | C4 | `citation-reviewer` · **chunk ~5–6 papers/agent** | cites resolve and say what we attribute |
| Synthesis faithfulness | C9 | `claims-reviewer` (one on the branch synthesis) | **separate dispatch** from per-paper claims |
| Deterministic | C10, C11, C12, C13, C14, C15, C16, C17 | scripts (`check_internal_titles.py`, `check_k_label.py`, `check_paper_test_refs.py`, `check_file_refs.py`, `check_inline_arxiv.py`, `check_retracted_terms.py`, `check_headline_numbers.py` — each `--gate <branch>`; + compile) · **whole group, first** | not an LLM reviewer |

**Enumeration mandate (the 2026-06-23 lesson).** Every LLM-reviewer prompt
demands *exhaustive enumeration*, not sampling, of the low-salience structured
regions where the recurring drifts hide — every status/catalogue-table row,
every bibitem, every occurrence of the branch's tier-terms (state-space GH vs
propinquity, derived vs match, withdrawn-vs-live), every docstring of the backing
module — one verdict per item. Prefer one enumeration-forced reviewer per
dimension over two vaguely-prompted ones (duplicate agents share blind spots).
After the panel, a one-agent **completeness-critic** names un-verified regions;
re-dispatch a focused reviewer for any gap. A full find→fix→re-scan
loop-until-dry is reserved for the *final* pre-certification convergence run.

## Run shapes: full runs vs delta-verification runs (2026-07-02, PI direction)

**The cost meta-lesson (seven group4 cycles):** the gate was priced as full
re-verification (~9–11 reviewers re-reading every paper, ~2–3M tokens/run),
but between runs the corpus changes by a ~10–20-locus remediation diff, and
from run ~4 onward the marginal yield was dominated by mechanical classes now
owned by the deterministic gates (C16/C17). The arc paid full price on "this
is probably the certifying run" four times. Run shapes are therefore:

- **Full run** — the complete protocol (all dimensions, whole-paper
  enumeration, completeness-critic). Fired at exactly two moments: the
  **first cert** of a fresh target, and the **final certifying pass**. Only a
  full run can produce **PASS** (an unexercised dimension still forces
  INCONCLUSIVE — unchanged).
- **Delta-verification run** — after any FAIL→remediation. Scope = the git
  diff since the last calibrated run: one reviewer per **affected** dimension;
  the changed loci **pasted into the prompt** with surrounding context
  (paste-don't-point — the reviewer does not re-read 3,000-line papers to see
  20 lines); **≥1 seed planted in the diff** per dispatched agent (calibration
  still per-agent); deterministic gates still run **whole-target** (they are
  free and guard the unchanged surface). Verdict: **CLEAN-DELTA / DEFECTS** —
  a delta run can never produce PASS; **a clean delta is the precondition for
  firing the full certifying run.** Skip the completeness-critic on deltas.
- **Rationale:** text that survived a calibrated panel does not lose its
  verification by sitting unchanged; the deterministic layer (now including
  the C16 zombie registry and the C17 headline-number registry) guards the
  unchanged surface between runs. Expected cost: ~0.5M tokens per delta cycle
  vs ~2.5M per full run.

**Model tiering with calibration as the safety net (2026-07-02).** The
citation and code dimensions default to **Sonnet-class** reviewers; claims
and synthesis stay **Opus-class** (they carry the judgment calls where even
Opus shows severity variance). The seeded-calibration design converts model
adequacy from a prior into a per-run *measurement*: a tiered agent that
misses its seed is de-calibrated and its dimension is re-dispatched on the
Opus tier (the run-6 P16 recovery pattern). When a dimension runs below the
Opus tier, plant **two seeds** for that agent (calibration resolution).

**Terser reporting contract (2026-07-02).** The enumeration mandate governs
the *reading*; the *report* is: defects + two-way upgrades + a compact
coverage checklist (which regions were enumerated). Reviewers do NOT
reproduce per-item verdict tables for sound items — soundness is evidenced
by the coverage checklist, the seeds, and PM spot-checks, not by report prose.

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
- 2026-07-04 — **Added C18 (project-course duration-language gate)** +
  `docs/authoring_conventions.md` rule 12 (PI direction, set at the synthesis
  time-language pass). No wall-clock units on the project's own course:
  historical durations FAIL (the "Three years ago the project…" incident);
  forward effort-estimates ("multi-month/-year") are advisory standing debt
  (~190 instances) pending the corpus sweep. Check self-tested;
  `debug/qa/check_duration_language.py`.
- 2026-07-04 — **Added the "Withdrawal-chronicle placement" NIT class** to §"Material vs
  nit" + `docs/authoring_conventions.md` project-wide rule 11 (PI direction, set at the
  synthesis retraction-language pass). Synthesis layer = current state + one-clause
  headline-grade flags with pointer; the chronicle (dates, process narrative) lives with
  the object of record (paper Status notes, claims register, CHANGELOG). Claims are
  "withdrawn/descoped," never "retracted" (no whole-paper retraction has occurred).
  Zombie rule unchanged: a withdrawn claim asserted as live is MATERIAL.
- 2026-06-16 — **Extracted** from `trunk.done.md` / `group3.done.md` (which were
  ~90% identical) into this single shared standard, per PI direction
  ("consolidate to a single definition of done"). The three branch files become
  thin profiles. No criterion was changed in the extraction — C1–C13 are
  verbatim the trunk/group3 standard, generalized (branch-specific watch-notes
  moved into the profiles).
- 2026-06-28 — **Added the "Dual-rule ERI framing (A = sparsity, B = accuracy)" section**
  (PI direction, during group4 pre-work). Codifies the framing-zombie rule: the codebase
  carries two ERI selection rules on purpose — pair-diagonal (A) for the sparsity/QC product
  (must be disclosed as an approximation wherever a sparsity claim uses it), global-M_L (B)
  for physics accuracy. A pair-diagonal number presented as the exact selection-rule value is
  a framing zombie. Enforced by the claims-reviewer (enumerate every sparsity/density/Pauli
  claim) + the new C16 entry `pair-diagonal-as-exact-sparsity`. Sharpens C3/C8/§1.5; shared
  (group4 primary, group3 P22/31 related). Repo study: `debug/sprint_group4_prework_memo.md`.
- 2026-07-02 — **Added C17 (headline-number registry gate), the "Run shapes" section
  (full vs delta-verification runs), model tiering with calibration as the safety net
  (+two-seeds-per-tiered-agent), and the terser reporting contract** (PI direction,
  post-7th-group4-cert cost review: one full run ≈ a 5-hour token session; seven
  full-price cycles on group4 with marginal yield converging to mechanical classes).
  A delta run can never produce PASS; a clean delta is the precondition for the full
  certifying run. C17 self-tests in `tests/test_headline_numbers_check.py` (6/6).
- 2026-06-28 — **Added the "Secondary stale-number / provenance-attribution is a NIT class"
  carve-out** to §"Material vs nit" (PI direction, set at group2 certification). The analog of
  the 2026-06-24 code-docstring carve-out: a stale/imprecise number on a NON-§C8-headline
  secondary result, or a wrong causal/provenance note on an otherwise-correct headline, is a
  fix-on-sight NIT, not a cert blocker — provided every DoD §C8 headline + tier + backing test
  is correct; anything touching a headline stays MATERIAL. Motivated by the group2 run-#4/#5
  thin-residual asymptote (perfect 11/11+6/6 calibration both runs, ~2 thin secondary/provenance
  items each, never a headline). Refines (does not relax) the counterfactual.
