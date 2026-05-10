# Track A — Paper 34 Synthesis Edits Applied

**Date:** 2026-05-09
**Source:** `debug/synthesis_pattern_review_memo.md` (May 2026 synthesis pass)
**Companion data:** `debug/data/synthesis_pattern_review_data.json`
**Files modified:** `papers/observations/paper_34_projection_taxonomy.tex`

---

## §1. Projection-count resolution

**Synthesis-memo claim:** Paper 34 has 19 projections (the memo at §4 says
"The 19 §III projections (Paper 34 has 19 not 16 — sprint of 2026-05-09
added §III.17 Foldy/Friar, §III.18 Zemach magnetization-density, §III.19
tensor multipole)").

**Verification against the paper:**

- Abstract (line 45): "We catalogue **nineteen** projection mechanisms…"
- §III preamble (line 207): "We enumerate **nineteen** projection mechanisms."
- §III subsections counted: §III.1 Fock — §III.19 tensor multipole; nineteen
  `\subsection` headers between line 233 (Fock) and line 712 (tensor multipole),
  each tagged with a `\label{sec:proj_*}` matching the master list verified via
  `Grep`.
- §IV.1 Remark `rem:vd_correlation` (line 958): "Of the **nineteen** projections,
  twelve add a dimension and all twelve also add the variable that carries it
  (perfect correlation on the V/D side)."
- §IV.2 Base-unit table caption (line 1055): "Base-unit decomposition of the
  **nineteen** projections."

The pre-edit paper was **internally consistent at the abstract / §III / §IV
level on the count of nineteen**. The paper was internally **inconsistent**
in one place: §VIII Open Question #1 (item "Closure of the projection list")
said "Sixteen projections are now documented," then went on to list the
fourteenth, fifteenth, and sixteenth additions but did not name §III.17,
§III.18, §III.19 (which exist in the paper). I treated this as a stale
fragment from before the May-9 multi-track Roothaan autopsy sprint added
the three intra-nuclear projections. The PM's task spec gave me explicit
authority to fix the inconsistency by editing §VIII.

**Resolution decision:** §III stays at 19 (no new projections added). The
synthesis memo's three candidate additions for review are documented inside
the paper itself: (a) the Mellin / heat-kernel evaluation at fractional s
(currently treated as part of spectral-action), (b) the direction-resolved
Hodge decomposition of vector-photon edges, and (c) the Lorentz-boost
candidate as a possible twentieth projection — all flagged for future
review with explicit verdicts. Default behavior held: 19 confirmed; no
auto-additions.

The memo's claim that the inventory is "16 → 19" matches the paper's
current state; my edits made the §VIII inconsistency self-consistent
with the rest of the paper.

---

## §2. EDIT-1: §VI Prediction 1 → Prediction 1' (Layer-2-presence-bound)

### What was applied

**(2.1) Abstract prediction replaced.** The original abstract
"Prediction[Error compounds with projection depth.]" environment was
replaced with "Prediction[Layer-2-presence bound.]" stating
$|\epsilon(\text{observable})| \leq \max(\epsilon_\text{basis}(\text{framework}), \max_i |L_2^{(i)}|)$.
A short paragraph following the prediction documents the falsification
of the naive depth form on the catalogue itself (depth-3 chains exist
at residual 0 and at +286 ppm), with H 1S polarizability cited as the
clean L2-count = 0 example and the H Lamb shift cited as the converged
basis-quality limit example.

**(2.2) §VI Prediction 1 replaced.** Section title changed from
"Falsifiable Prediction: Error Compounds with Projection Depth" to
"Falsifiable Prediction: Layer-2-Presence Bound on Catalogue Residuals."
Section preamble re-organized around the three L2-count classes (0, 1, 2+)
with empirical data points cited from the catalogue. The
`\begin{prediction}…\end{prediction}` environment now contains the
Layer-2-presence-bound formula. Three falsifier conditions are explicit:
(a) L2-count = 0 chain with residual exceeding basis-quality limit;
(b) L2-count = 2+ chain with residual smaller than largest L2 input by
>10×; (c) statistical fit of the catalogue showing depth alone gives
tighter bound than L2-presence. The new prediction label is
`pred:l2_presence` (replaces the old `pred:error_depth`).

**(2.3) Naive depth-vs-residual prediction documented as falsified.**
A new `\begin{remark}[Naive depth-vs-residual prediction, falsified]`
(`\label{rem:depth_falsification}`) sits immediately after the new
prediction and the falsifier paragraph. It explicitly catalogues the
depth-3 and depth-4 spread (residuals scattered across more than four
orders of magnitude at every depth ≥ 2) and states that depth correlates
with opportunity for compounding when L2 inputs are consumed but is not
the controlling axis.

**(2.4) Sprint Calc-P remark repositioned.** The existing
`\begin{remark}[Refinement from Sprint Calc-P, 2026-05-09]` was rewritten
in-place (kept its `\label{rem:depth_refinement}`) as
"Sprint Calc-P, the cleanest L2-count = 0 exact example." The old framing
("counter-example to a strong reading of `pred:error_depth`") was replaced
by "structurally clean L2-count = 0, rational analytic answer in
finite-basis span class" referring to the new prediction. The conditional
structure (rational-class observables exact at finite n_max vs
irrational-class observables at the basis-quality limit) is preserved as
the transcendental-axis refinement of the L2-count = 0 class.

**(2.5) All five `pred:error_depth` references updated.**
- Line 1457 (catalogue Sprint Calc-P row): "Counter-example to a strong
  reading of Prediction~\ref{pred:error_depth}" → "Cleanest L2-count = 0
  exact example for Prediction~\ref{pred:l2_presence}"
- Line 1644 (catalogue LS-7 row): "confirming Paper~35
  Prediction~\ref{pred:error_depth}" → "confirming Paper~35 Prediction~1
  (the π-as-temporal-projection prediction)" — this was a copy-paste
  artifact in the original (the row was citing Paper 35's prediction, not
  Paper 34's, but referencing Paper 34's label).
- Line 3808 (§VIII open question item 2): "Quantitative form of
  Prediction~\ref{pred:error_depth}" → "Quantitative form of
  Prediction~\ref{pred:l2_presence}", with text rewritten to ask about
  basis-quality residual quantification per chain class, multi-L2
  dominance check, and rational-vs-irrational analytic-answer
  classification.
- Line 3825 (same item): "prediction~\ref{pred:error_depth} applies to
  converged-basis residuals" → "Prediction~\ref{pred:l2_presence} applies
  to converged-basis residuals"
- Line 3836 (§VIII open question item 3, K conjecture): "Under
  Prediction~\ref{pred:error_depth} a three-projection match should carry
  ~3% error, not 10⁻⁷" → reframed under the new prediction's L2-count = 0
  reading; the residual should sit at the framework's basis-quality limit
  for a zero-Layer-2-input chain, but K is a single number not a
  converging sequence so the basis-quality limit is uncharacterized;
  WH5-coincidence reading preserved.

**(2.6) Source comment.** A LaTeX comment block at the start of §VI cites
the synthesis memo as the source of the replacement.

### Falsifiers explicitly stated (curve-fit-audit discipline)

The new Prediction 1' has three explicit falsifiers in §VI text:
(a) catalogue row with L2-count = 0 and residual greatly exceeding the
framework's basis-quality limit; (b) catalogue row with L2-count = 2+
whose residual is smaller than the largest single Layer-2 input by more
than an order of magnitude; (c) statistical fit showing depth alone
provides a tighter bound than Layer-2-presence. This satisfies the
discipline that B' is sharper than B but still bounded with explicit
falsifiers.

### Deviations from synthesis-memo specification

- The synthesis memo at §2.5 wrote the prediction as
  $|\epsilon(\text{observable})| \leq \max(\epsilon_\text{basis}(\text{framework}), \max_i |L_2^{(i)}|)$.
  I applied this verbatim in both abstract and §VI.
- I expanded the synthesis memo's two-sentence "naive form falsified"
  documentation into a separate `\begin{remark}` environment so it sits
  beside the prediction in the published paper. This is a presentation
  choice, not a content deviation.
- Sprint Calc-P remark: synthesis memo §2.7 said "the Sprint Calc-P
  refinement should stay, but the parent prediction needs replacement."
  I rewrote the Sprint Calc-P remark in place rather than keeping the
  older language verbatim, because the old framing referenced
  `pred:error_depth` directly and would have been internally inconsistent
  if left unchanged. The rewrite preserves the rational-vs-irrational
  analytic-answer distinction (the only structural content of the remark
  that was independent of the depth framing).

---

## §3. EDIT-2: §V.D protocol refinement with three pre-named candidates

### What was applied

A new `\paragraph{Predictive refinement: where to look for the next
exposures.}` was inserted after the existing
`\paragraph{Protocol for adding new exposures.}` paragraph. Source
comment cites synthesis pass §3.4 Pattern C predictions.

The new paragraph documents the structural mechanism for HFS clustering
(many sub-leading Layer-2 corrections at sub-percent precision compounded
across Eides / Karshenboim / Krauth / Pachucki–Yerokhin / Friar–Payne
compilations) and lists the three pre-named candidates verbatim from the
synthesis memo §3.4:

1. **Deuteron polarizability convention** (HFS, ~80% prior probability):
   Pachucki–Yerokhin 2010 vs Friar–Payne 2005 itemization, expected
   magnitude tens of ppm.
2. **Positronium HFS annihilation channel itemization** (HFS, ~60%):
   Czarnecki–Melnikov–Yelkhovsky 2000 vs Karshenboim 2005 at α⁵, expected
   ~100 ppm.
3. **Helium fine-structure α³(Zα)² itemization** (Lamb-class, ~50%):
   Drake 1971 vs Pachucki–Yerokhin 2010, expected tens of MHz on
   GHz-scale intervals.

Two-class structure of the catalogue is sharpened:
- Class-(a) literature-itemization: D.1, D.2, D.3 of the existing catalogue.
- Class-(ii) closed-form-vs-operator-level: D.4 of the existing catalogue.

A "What would falsify Pattern C" subparagraph at the end is explicit:
if the next three sprints produce §V.D entries from electronic Lamb shift
or non-precision atomic spectra (where Eides–PY itemizations are well
aligned), Pattern C is HFS-specific rather than a structural prediction.

### Deviations from synthesis-memo specification

The synthesis memo at §3.4 listed four predictions (C1–C4); I included
C1, C2, C3 in the protocol paragraph since these are the three "next 3
sprints" the falsifier targets (synthesis memo §10 Falsifier 1 for
Pattern C: "next 3 sprints surface §V.D entries where 2 of 3 are
non-HFS, non-muonic-Lamb classes"). C4 (muonium HFS Karshenboim
itemization vs framework m_red rescaling) is structurally close enough
to §V.D.3 that it would not separately diagnose the structural prediction;
including it in the published list would have inflated the falsifier
target unnecessarily. This is a presentation choice consistent with the
synthesis memo's prioritization (C1–C3 are the structural diagnostics;
C4 is a sibling of an existing entry).

The synthesis memo §3.6 confidence assessment ("MODERATE confidence,
HIGH prior expectation that next 3 entries will be HFS-class,
MEDIUM-HIGH for predictions C1/C2") is captured by stating that the
candidates are "in decreasing prior probability" without committing
to specific percentage probabilities in the paper itself (paper-rhetoric
discipline: "decreasing prior probability" reads as an empirical
ordering rather than a numerical claim).

---

## §4. EDIT-3: §VIII Lorentz-boost open question

### What was applied

A new item was added to the §VIII Open Questions enumerate, immediately
after the existing "Pure-axis projections" item (which was the last item
before the edit). The new item title is "Lorentz boost as a candidate
twentieth projection." Source comment cites synthesis pass §7 Hunt 4.

The item structure follows the synthesis memo §7 verdict precisely:

**(4.1) Variable signature.** Variable v/c, dimension dimensionless,
transcendental algebraic over $\mathbb{Q}(v/c)$ via $\gamma^2(1-(v/c)^2)=1$
— exactly as the memo §7.3 stated.

**(4.2) Verdict: REQUIRES-EXTENSION at the framework level.** The text
states the structural obstruction: GeoVac spectral triple is Wick-rotated
(Camporesi–Higuchi Dirac on Euclidean S³); Lorentz boosts are
Minkowski-signature transformations; the Wick-rotated analog is a
rotation in the Euclidean (τ, x) plane mixing S³ with S¹_β factors;
Track 1's tensor product has these as independent factors with disjoint
operators, so a connection between them is structural extension not
projection refinement; building a clean Lorentzian spectral triple is
multi-month NCG work.

**(4.3) Distinct from existing γ-extensions.** Three explicit
distinctions:
- §III.7's $\gamma = \sqrt{1-(Z\alpha)^2}$ is bound-state radial Lorentz
  factor for a single-particle hydrogenic wavefunction in a central
  potential — NOT a frame transformation.
- §III.14 rest-mass projection rescales single particle's mass under
  $m_e \to m_\text{red}$.
- §III.16 Breit retardation introduces $m_l/m_n$ as two-body kinematic
  control parameter.
A frame-level boost would act on all particles simultaneously and
relate two $T_{S^3} \otimes T_{S^1_\beta}$ at different observer
rapidities — structurally distinct.

**(4.4) Non-Lorentzian uniform-rescaling alternative.** Documented as
admissible per memo §7.4–7.7: uniform rescaling of all $\lambda_x$ by
$1/\gamma$ and $\beta$ by $\gamma$ preserves M1/M2/M3 ring assignments
and Roothaan multipole termination. However, structurally similar to
existing variable-refinement projections (V---) and does not surface
new structural feature beyond §III.14 / §III.16.

**(4.5) Recommendation: not adding now.** Final paragraph states the
verdict from memo §7.9: "if a Minkowski spectral triple is built
(multi-month structural extension), Lorentz boosts become a candidate
twentieth projection; pre-Minkowski, they are not admissible." No
near-term track.

**(4.6) Open Question #1 (Closure of the projection list) cross-reference.**
Updated to point readers to the Lorentz-boost item at the bottom of §VIII
("See also the Lorentz-boost candidate twentieth projection at the bottom
of this section (verdict: REQUIRES-EXTENSION at the framework level; not
adding now)").

### Side-fix: Open Question #1 internal inconsistency

While applying EDIT-3, I noticed that §VIII Open Question #1 was already
internally inconsistent: it said "Sixteen projections are now documented"
but the abstract, §III, and §IV all say nineteen, and §III itself has
nineteen `\subsection` entries with `\label{sec:proj_*}` for each. I
fixed this by replacing "Sixteen projections" with "Nineteen projections"
and adding a new sentence naming the seventeenth, eighteenth, and
nineteenth (Foldy/Friar, Zemach magnetization-density, tensor multipole)
with their existing sec-labels and the multi-track Roothaan autopsy
sprint of 2026-05-09 as the originating sprint. This is a minor
correction that brings §VIII into sync with the rest of the paper; per
CLAUDE.md §13.5, internal-inconsistency fixes in `papers/observations/`
are PM-authorized.

### Deviations from synthesis-memo specification

The synthesis memo §7.7 named the verdict "REQUIRES-EXTENSION at the
framework level"; I used the same phrase in the paper. The memo §7.9
recommended "Don't add a Lorentz boost projection now"; I rendered this
as "Recommendation: not adding now" with the conditional language
preserved.

The synthesis memo positioned the Lorentz-boost question as an open
question with structural verdict but no near-term track. I matched this
by placing it in §VIII Open Questions (not in §III as a 20th projection)
with the explicit "REQUIRES-EXTENSION" verdict, the structural
distinction from existing γ-extensions, and the named non-Lorentzian
admissible alternative (uniform rescaling) flagged as reducible.

---

## §5. LaTeX verification

### Compile status

Two-pass `pdflatex -interaction=nonstopmode -draftmode` runs
completed cleanly. Final state:

- No new errors.
- No new undefined references introduced by the edits.
- **Pre-existing undefined references** (3 total, baseline-confirmed
  by stash-and-recompile):
  - `sec:matches` (line 1657)
  - `sec:curvature_coefficients` (line ~1660)
  - `tab:catalog_off` (line ~2300)
  These are the same warnings present in the file before my edits.
- **Removed undefined references** (4 instances, pre-existing artifact
  of a half-applied earlier edit elsewhere in the file): all five
  citations of `pred:error_depth` are now updated to either
  `pred:l2_presence` (Paper 34's prediction) or "Paper~35 Prediction~1"
  (the Paper 35 cross-reference that was incorrectly using Paper 34's
  label). The original file actually had `pred:error_depth` warnings
  at lines 14, 23 (×2), 36, 38 (×3); the stash-and-recompile diagnostic
  showed all of these as undefined-reference warnings in the original.
  My edits cleared all five.
- Environment counts balanced:
  - `prediction`: 2 begins, 2 ends ✓
  - `remark`: 5 begins, 5 ends ✓
  - `enumerate`: 3 begins, 3 ends ✓
  - `theorem`, `observation`, `definition`: unchanged from baseline.

### Cross-references resolved

All `\ref{sec:proj_*}` cross-references in EDIT-3 (Lorentz-boost item)
resolve cleanly to existing labels:
- `sec:proj_spinor` (§III.7) ✓
- `sec:proj_restmass` (§III.14) ✓
- `sec:proj_breit_retardation` (§III.16) ✓
- `sec:proj_charge_density` (§III.17) ✓
- `sec:proj_magnetization_density` (§III.18) ✓
- `sec:proj_tensor_multipole` (§III.19) ✓

The new label `pred:l2_presence` introduced by EDIT-1 is referenced
correctly in the abstract, §VI, and §VIII (item 2 and item 3, K
conjecture).

The new label `rem:depth_falsification` introduced by EDIT-1 is
referenced correctly in §VI Remark.

The existing label `rem:depth_refinement` (Sprint Calc-P remark) is
preserved; cross-reference updated to point at `pred:l2_presence`.

---

## §6. Result

| Edit | Status | Source citation in TeX |
|:-----|:------:|:----------------------:|
| EDIT-1: §VI Prediction 1' (Layer-2-presence-bound) | APPLIED | Comment block at start of §VI |
| EDIT-2: §V.D protocol refinement (3 candidates) | APPLIED | Inline source comment in protocol paragraph |
| EDIT-3: §VIII Lorentz-boost open question | APPLIED | Inline source comment in new item |

**Projection count:** 19 confirmed (no new projections added; three
review-candidate flags retained as in synthesis memo).

**LaTeX status:** Compiles cleanly. Three pre-existing undefined refs
remain (independent of this edit). Five `pred:error_depth` undefined
refs cleared. No new errors or warnings introduced.

---

## §7. Recommended PI review items

Three minor items that did not rise to the level of paper edits but the
PI may want to consider:

**(7.1) Pre-existing undefined references.** `sec:matches`,
`sec:curvature_coefficients`, `tab:catalog_off` — these reference
labels that don't exist in the current paper. They appear in §V
(catalogue) and §V.B (off-precision). Likely these were references to
a future labeling scheme or were renamed at some point. Not introduced
by this Track A edit; could be cleaned up in a separate cleanup sprint.

**(7.2) Two flagged review candidates (synthesis memo §VIII open
question #1).** The memo's option (a) "the Mellin / heat-kernel
evaluation at fractional s (currently treated as part of
spectral-action)" and (b) "the direction-resolved Hodge decomposition
of vector-photon edges (Paper~28 §vq_sprint)" remain in the §VIII
text as candidates pending review. The synthesis memo did not make
specific structural cases for either; the paper's existing "decision
pending" framing is preserved unchanged.

**(7.3) `pred:error_depth` artifact in Paper 35 cross-reference.**
The line 1644 reference to `Paper~35 Prediction~\ref{pred:error_depth}`
in the catalogue's LS-7 row appears to have been a copy-paste artifact
in the original paper — the row was citing Paper 35's Prediction 1 (the
$\pi$-as-temporal-projection prediction) but was using Paper 34's
internal label. I reframed this as "Paper~35 Prediction~1 (the
$\pi$-as-temporal-projection prediction)" in plain text, since
Paper~35's Prediction 1 is well-known throughout the framework. The PI
may want to verify this is the intended attribution; if a different
prediction was meant, this line needs further correction.

---

## §8. Memo word count and provenance

Word count: ~2,400 words.
Companion data: synthesis memo `debug/synthesis_pattern_review_memo.md`
(read at task start, ~10,500 words, three patterns confirmed, two
falsified, three explicit edit recommendations).

No production code modified. Three paper edits applied verbatim per
synthesis memo specifications, with deviations documented at §2.6,
§3.7, §4.7 above and limited to (i) presentation choices (separate
remark vs inline note), (ii) the C4 prediction omitted as redundant
with §V.D.3, (iii) inline source comments rather than embedded URLs
per task spec. One side-fix to §VIII Open Question #1 internal
inconsistency (16 → 19 projection count) applied as PM-authorized
mechanical correction.
