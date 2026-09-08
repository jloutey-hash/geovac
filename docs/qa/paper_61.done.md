# Paper 61 (Bessel-Moment Period Algebra) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-61 scope + deltas + watch-notes.

> **STATUS: NEVER CERTIFIED. FULL run 2026-09-07 = FAIL (3 LARGE), remediated;
> same-day DELTA = DEFECTS (1 LARGE + 9 MATERIAL), also remediated. NOT certified.**
> Rewritten 2026-09-07 *after* the FULL run, from what the panel found; the
> DELTA record is folded in below.
>
> **What the DELTA measured, and why it is the more useful half.** All three
> LARGE fixes landed and were independently re-verified sound (the Sp₄ note's
> four steps each check; the abstract inversion is gone; the honest-scope
> sentences are all intact and all in the honest direction — 21 enumerated, 0
> inverted). The delta nonetheless found:
>
> - **1 LARGE — a goalpost ratifying the retired claim.** `docs/qa/paper_59.done.md`
>   did not merely repeat "Galois in Sp4(Z)"
>   `[retracted 2026-09-07: p61-galois-in-sp4-z]`, it *certified* it: "[SYMBOLIC] prose
>   is now correct, not an inflation". A `.done.md` is what the next certifying run
>   measures against, so the defect had a mechanism for returning. This is the
>   `.done.md`-as-vector class, first instance.
> - **1 LARGE — two remediated claims had no `check_retracted_terms.py` entry**
>   (§9 hard rule). Its first consequence was the survivor below.
> - **The retired false universal alive in an ACTIVE paper.**
>   `paper_56_tannakian_substrate.tex:1848` still said any full-modulus sweep
>   "passes through **every** CM fibre"
>   `[retracted 2026-09-07: p61-every-cm-fibre-universal]` — while citing Paper 61
>   above. Also alive in `memory/cosmic_galois_elliptic_rung1.md`, which is loaded
>   into **every session** and therefore reseeds itself, and in two live drivers.
> - **4 defects in prose the FULL-run remediation had just written**, including a
>   **new false universal** ("Every transcendental below is placed…") — the exact
>   class the same run had flagged twice. Also: it mis-tiered the companion
>   (Paper 59's own title is "Genus One at the Third Center", so this paper is at
>   the *same* rung, not one above), contradicted the paper's own π accounting by
>   routing π through temporal compactification, and named `K(½)` the lemniscate
>   constant (they differ by exactly √2; the *symbol* is the corpus convention,
>   the *name* was wrong).
> - **My two new C16 entries missed 8 of 8 live survivors.** They were "proven to
>   discriminate" against the `.tex` at git HEAD — which is not the same as
>   reaching the loci that actually survived. Causes: ASCII patterns against a
>   Unicode corpus; exemption vocabulary drawn from the surrounding *correct*
>   text, which the file's own 2026-09-03 rationale forbids in as many words; and
>   `files` lists built from memory. Fixed as a **class** (Unicode folding in
>   `_strip_markup`, en-dash → `--` so the three strict-`--` entries keep working)
>   and fire-tested 12/12, including all three legitimate near-neighbours
>   (Broadhurst–Mellit *determinant*, *monodromy* in Sp₄(ℤ), "infinitely many").
>
> **Two-way upgrades, both PM-verified before acceptance.** `W₀` went
> [OPEN] → **[SYMBOLIC]** (`eq:W0`, `W₀ = π²/ρ²`, sympy step-by-step + mpmath at
> three ρ to 40 digits), with a backing test fire-tested against 5 named wrong
> answers — including `W₀ = 1`, which is what the *existing* Wronskian test
> normalises to and the reason the false claim was invisible. And the λ identity
> is precision-limited, not `10⁻⁴¹`: 2.7e-51 at 50 digits, 1.2e-91 at 90.
>
> **Provenance warning — read before trusting this file.** Both versions were
> written by the PM. The first was drafted hours before the review, and the panel
> found it defective in five ways: three wrong backing-test names; `C5` listed as
> a deterministic criterion (the deterministic K screen is `C12`); `C22` and the
> FULL-run arithmetic audit omitted; no C7 trunk-dependency list; and a phantom
> "Katz" watch-note — Katz appears nowhere in Paper 61, direct evidence it was
> drafted from Paper 59. It also omitted the transcendental-tagging axis entirely
> for a paper made of transcendentals, and that omission hid a real defect.
> **A reviewer-authored definition of done is a known weakness of this record.**
> The criteria below are derived from the panel's findings rather than the PM's
> priors — a mitigation, not a fix.
>
> **DoD RATIFIED BY PI, 2026-09-08.** The criteria in this file are now the
> agreed goalpost for Paper 61, which retires the reviewer-authored weakness
> above as a blocker — the provenance record stays because it is why the
> criteria read as they do, not because it is still an open objection.
>
> **Ratifying the DoD is NOT certifying the paper.** Paper 61 remains
> **NOT CERTIFIED**. What is settled is *what would count* as done; what is
> not settled is whether the paper meets it. That still requires a **clean
> DELTA followed by a FULL run**, and the 2026-09-08 DELTA was not clean
> (1 LARGE + 9 MATERIAL, all remediated). A ratified goalpost measured
> against nothing is exactly the `.done.md`-as-vector failure this same
> DELTA found in `paper_59.done.md` — so the two must not be conflated.

---

## Scope

**Paper:** `papers/group3_foundations/paper_61_bessel_moment_periods.tex` —
group3, the periods / Tannakian arc, siblings 55–57. Also in the `group3` scope,
so `/qa group3` exercises it.

**The gate scope must include the seam.** Paper 59 lives in
`papers/group2_quantum_chemistry/`, Paper 61 in `papers/group3_foundations/`.
They were one paper until 2026-09-06 and describe one object. A scope of `{61}`
alone contains neither Paper 59 nor the group3 synthesis while C9 is GATING —
that gap is why the run found nine test-suite loci still crediting Paper 59 for
Paper-61-owned labels.

**C9 — GATING, and measured non-empty.** Live surface in
`group3_foundations_synthesis.tex` (L1198–1250, L1361–1378),
`group2_quantum_chemistry_synthesis.tex` (L666–682), `papers/INDEX.md:50`,
`docs/topic_to_paper_lookup.md`. Do **not** record C9 as N/A — that premise was
false for Papers 59 and 60 in the same week, and Paper 60's version concealed the
second locus of a LARGE.

---

## Branch-defining criterion: **inherited backing across a paper split**

Every Paper-61 claim is backed by tests named `test_paper59_*` or
`test_routeC_momentum.py`; **no `test_paper61_*` file exists**. Shared backing is
not wrong — duplicating would be worse — but the claim→artifact mapping is then
the only thing holding the dependency, and in this corpus that mapping has
repeatedly been wrong.

**A run PASSes this criterion only if all four hold:**

1. **Every cited test proves the *Paper-61* claim**, not merely the Paper-59
   claim it was written for. Answer per claim, explicitly.
2. **The row set is complete and correctly counted.** Measured 2026-09-07:
   **eleven** rows — lines 87, 496, 497, 498, 499, 503, 504, **505**, 506, 507,
   508. Two traps: **line 87 sits ~400 lines from the others**, in the group3
   region, and backs eq:kw + the certified 66-digit T2 + the PSLQ negative — the
   paper's most quotable claim; and **line 505 is the shared `| 59/61 |` row**,
   invisible to a `^| 61 |` grep, which is exactly the seam this criterion
   polices.
3. **`rests on:` edges reconcile.** Nine declared; row 508 deliberately has none;
   **row 505 is unaccounted for.**
4. **The `--slow` layer actually ran.** Ten `@pytest.mark.slow` tests carry the
   deep-precision witnesses — the 25-digit −π/0/2π pairing the paper cites by
   name, the eq:kw certification witness, and the M₀ monodromy the whole Sp₄
   argument rests on. CLAUDE.md §14 skips these by default. **State whether
   `--slow` was passed and report passed-vs-skipped per test.**

---

## Dimensions to exercise (FULL run)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1, **Opus**. **Run with
  `--slow`.** The nine `test_paper59_*` files the paper cites, plus
  `test_routeC_momentum.py`, plus **`test_paper59_corner_sigma2.py`** and
  **`test_paper59_four_subspace_lambda.py`** — both open "Backing test for
  Paper 61" and neither was in the first version of this list. Mandates:
  independent re-derivation; fire-test every guard; **restricted evaluation**
  (W1).
- **Paper claims / prose (C3, C5, C6, C8)** — `claims-reviewer` ×1, Opus,
  enumeration-forced over all **41** tier tags (recount 2026-09-07: 16 MEASURED,
  9 OBSERVATION, 5 SYMBOLIC, 5 OPEN, 6 MEASURED-with-inline-backing; the “35”
  predated the W6/W3 edits, which added two more).
- **External citations (C4)** — `citation-reviewer` ×1, **Opus**. 19 bibitems but
  a high-risk surface (Chowla–Selberg, Broadhurst–Roberts, Fresán–Sabbah–Yu,
  modular / cosmic-Galois). Three-way verdicts; UNVERIFIABLE is not a pass.
- **Synthesis (C9)** — `claims-reviewer` ×1. GATING; see Scope.
- **Arithmetic audit** — the shared criteria fire this on FULL runs; the first
  version omitted it. Assign it.
- **Deterministic (C10–C22)** — `--gate paper_61`. **State grounded counts, and
  record zero-surface criteria as UNMEASURED, not PASS.**
- **Completeness-critic** ×1. **Brief it from the reviewers' own reports, not
  from a PM summary** — a PM summary omitted two files in the 2026-09-07 run and
  produced two phantom coverage gaps.

### Deterministic criteria with **zero Paper-61 surface** (measured)

A criterion the target does not exercise is *unmeasured*, not passed:

| Criterion | Surface | Note |
|---|---|---|
| **C16** | **zero registry entries** name a Paper-61 locus or backing module | the gate owning the zombie class has nothing pointed here |
| **C20** | **zero `Author~Year` forms** | the bibitem-free layer is invisible **by construction** — and "Chowla–Selberg" carries a headline `[MEASURED]` result with **no bibitem at all** |
| C12 | no B/F/Δ content | vacuous |
| C6 | no graph Laplacian / S³ / −(n²−1) | vacuous |
| C15 | every arXiv ID is inside a bibitem | vacuous |
| C18 | no duration language | vacuous |
| C22-B | no `test_paper61_*` file exists | vacuous |
| C14 | one item (`benchmarks/certified_reference/`) | one-item gate |

C21 checked 2 annotations, C17 grounded 1 family — both **zero** before this run.

### C7 — cross-paper dependencies (omitted from the first version)

Two are live, load-bearing, and **unpriced**:
- **Paper 35 / WH7** — L469–471, "consistent with the compactification reading of
  Paper 35". WH7's status was rewritten **2026-09-06, the day before the run**
  (prior art found; an open hole — a compact *interval* is not a compact *circle*,
  and imaginary vs real time are identified silently).
- **Paper 56** — L234, L306–308, "answers the Paper 56 seam test T-2 in the
  negative". Paper 56's live status is Reading A: an abelianized-level-4
  *homomorphism*, not a closed immersion.

---

## C8 headline claims (frozen goalposts, backing corrected)

| # | Claim | Tier | Backing **(corrected 2026-09-07)** |
|---|---|---|---|
| 1 | `λ(τ(ρ)) = 1 − ρ` | `[MEASURED, 10⁻⁴¹; classical]` | `test_routeC_momentum.py::test_cosmic_galois_family_is_gamma2`. Re-derived to **1e-60** independently — prose understates. |
| 2 | CM fibres τ=i, τ=i√2 | `[MEASURED, ~10⁻⁵¹]` | `::test_cosmic_galois_cm_periods_are_gamma_values` — **both** legs; exact at dps=80. |
| 3 | `T2`, 66 digits | `[MEASURED, certified]` | `test_paper59_t2_value.py` pins **~21 digits**; **digits 35–66 rest on `debug/beta2_t2_*.py`, pruned by policy**. Registry stores a float64 → C21 checks 17 of 66. |
| 4 | `W₀D⁻²`; `B = πΩ`; Sp₄ | `[SYMBOLIC + MEASURED]` | `test_paper59_bessel_moment_algebra.py`. **See W2/W3 — the Galois statement is false and `W₀` is never computed.** |
| 5 | conductor-4 via `θ₃²` | `[MEASURED]` | `test_paper59_theta_chi4.py` (4 legs) |
| 6 | q-expansion `16, −256, 2112` | `[MEASURED, ~10⁻³²]` | **`test_paper59_bd_jacobian.py`** (v1 named `resurgent_skeleton` — wrong); pins only the leading 16. |
| 7 | fold `Φ(ρ)=Φ(1/ρ)/ρ²` | `[MEASURED, 2.2×10⁻⁹]` | **`test_paper59_coarea_reduction.py`** (v1 named `resurgence_corner` — wrong), and it asserts **5e-3**, ~20 orders coarser. |
| 8 | searches decisively negative, wt ≤3, h ≤10, 64 digits | `[MEASURED]` | **NO regression test.** v1 named `bessel_moment_algebra`, whose line 26 says "**No PSLQ**". Real backing: `benchmarks/certified_reference/entries_t2.py`. |
| 9 | Eisenstein/CM period, not a cusp-form L-value | `[OBSERVATION]` | `::test_paper59_gamma2_first_cusp_form_is_weight_6` |
| 10 | finite reduction remains open | `[OPEN]` | — |

---

## Watch-notes, ranked by what the run found

**W1 — restricted evaluation. Top note, and not hypothetical.** The companion
Paper 60 review found *its* headline exponent was an artifact of an undeclared
60-bohr radial box. For every Paper-61 result that is *exactly* a named constant
or *exactly* invariant, enumerate what the evaluation object excludes and exhibit
the claim on the unrestricted object.

**W2 — `Gal ⊆ Sp₄(ℤ)` is FALSE** (abstract L52–53, intro L86–87, body L619–621,
plus the test docstring and matrix row 499). Sp₄(ℤ) is discrete; a Galois group
inside it is finite, hence all solutions algebraic — contradicting the paper's own
irregularity results and the exponential torus `(ℂ*)²` Paper 59 establishes. The
test proves a **monodromy** statement. Paper 59:585 has it right as `Sp₄(ℂ)`.

**W3 — "the transcendence cancels in their determinant" is FALSE**
`[retracted 2026-09-07: p61-w0-transcendence-cancels]`, and `W₀` was
never computed anywhere. Abel fixes only the D-dependence; the test normalizes
`W₀` to 1 by construction — so no artifact in the corpus could have seen the
cancellation fail.

> **UPGRADED 2026-09-07 (same-day delta, two-way).** `W₀` is no longer
> “derived, [OPEN] pending confirmation” — it is **[SYMBOLIC]**, now in the
> paper as `eq:W0`, `W₀ = π²/ρ²`, with the branch-point derivation written
> out. PM-verified two ways before promotion:\ sympy step-by-step
> (Σx_c = 0; ∏Q′(x_c) = −16ω²; Vandermonde = 4iω/ρ²) and mpmath at
> ρ = 0.3/0.5/0.71 to 40 digits, ratio 1.0. So **π² does not cancel — it is
> squared** — and `W₀` is ρ-dependent. The withdrawal is a **net gain**: an
> exact value replaced a false cancellation. Backing test owed (§13.4a).

**W4 — the abstract inverts the body's own denial**, crediting the
integer-relation *searches* with the cosmic-Galois placement; L296–302 exists to
deny exactly that. Check against L304.

**W5 — driver-only tolerances presented as test-backed** (four tags). The claim
matrix states the split honestly; the paper does not.

**W6 — transcendental tagging.** Paper 61 cites **neither Paper 18 nor Paper 34**.
Untagged: Γ(1/4)², Γ(1/8)Γ(3/8), β(2)/Catalan, ζ, θ₂/θ₃/θ₄, e, ϖ, √π. Against
CLAUDE.md §4 and the standing memory rule. Paper 59 cites both.

**W7 — the searched-box negative is largely a PASS.** 5 of 7 loci carry an
explicit box and a reviewer independently reproduced the result with a positive
control. Record as *passed*, not as silent absence. Residual risk: the abstract
(W4) and the diagonal-A ring.

**W8 — the Paper 59 seam.** *(RESTATED 2026-09-07 — the previous wording was wrong
and would have caused a false sweep. It claimed `paper_59:1046` CONTRADICTS the row-8
finding. It does not:\ `tests/test_paper59_corner_sigma2.py:209–266` really does hold
`test_paper59_pslq_negative_disc4_guarded` — the disc-4 weight≤2 ring, height≤40,
matched decoy, and a planted `2π+3ϖ` positive control. What has NO test is Paper 61's
**headline** weight≤3 corrected-ring negative at h≤10, 64 digits. The defect is a
**scope ambiguity** in `paper_59:1047`, not a false backing — a pass acting on the old
wording would have “corrected” a true sentence.)* `paper_59:1046` asserts `corner_sigma2.py` holds "the
guarded PSLQ-negative backing the companion Paper 61", contradicting the row-8
finding; `paper_59:718–724` restates row 4 **without** Paper 61's tier-scope
caveat.

**W9 — false universals**: "sweeps the **entire** real locus" and "passes through
**every** CM fibre" `[retracted 2026-09-07: p61-every-cm-fibre-universal]`, both
contradicted by the paper's own L373–375.

> **Swept claim-wide 2026-09-07 (same-day delta).** The first remediation stopped
> inside Paper 61. The universal survived at **Paper 56:1848** — an ACTIVE paper,
> citing Paper 61 fourteen lines above — and in `memory/cosmic_galois_elliptic_rung1.md`,
> which is loaded into EVERY session and so reseeds itself. Also in two live
> drivers. All corrected to “infinitely many, not every one”; the exclusion is
> λ = 1−ρ, so only the REAL locus λ < 1 is swept (disc −3 and λ=2 are off it).
> Dated `debug/sprint_*` memos and CHANGELOG entries keep the old wording as
> historical record (§13.11) and are deliberately NOT in the gate's `files`.

**W10 — absences nobody owns.** No Status/History note recording the 2026-09-06
split (Paper 59 carries a live Zenodo DOI containing this material); `\date` is
byte-identical to Paper 59's and predates the document; no conclusion section, so
the `[OBSERVATION]`/`[OPEN]` tiers are restated compactly only in the abstract.

---

## Standing remediation owed (2026-09-07)

1. The three LARGE: W2, W4, and the row-8 backing map.
2. `W₀ = π²/ρ²`, replacing the false "transcendence cancels" (W3)
   `[retracted 2026-09-07: p61-w0-transcendence-cancels]`.
3. Broadhurst-**Mellit** → Broadhurst-**Roberts** at four loci including the
   abstract; swap in arXiv:2012.03523; drop `broadhurst2016` from that bundle.
4. Row 505's `rests on:` edge; the eleven-row count in every record.
5. A C16 entry with a Paper-61 locus; a C20-visible form (or a bibitem) for
   Chowla–Selberg.
6. `docs/claims_register.md` has **no row for Papers 54–61** (stops at 45), so
   C3's register leg is unmeasurable for this target.
7. §3 ledger row 132 names its own upgrade condition ("β(2) measured in T2 at
   ≥32 digits"); Paper 61 reports it at **64** and negative — the dead end is
   *strengthened*, but the ledger still reads pending.

## Change log

- **2026-09-07 (v2).** Rewritten from the FULL run's findings after v1 was found
  defective in five ways. Backing corrected on rows 6/7/8; `--slow` made an
  explicit gate condition; the eleven-row count and the `59/61` seam row
  recorded; zero-surface criteria enumerated as UNMEASURED; C7 added; the
  arithmetic-audit dimension restored; the phantom Katz note deleted; W1 promoted
  to restricted-evaluation on the strength of the Paper 60 box artifact. v1
  preserved in the run's scratchpad.
- **2026-09-07 (v1).** Created. Paper 61 had no `.done.md` and no scope entry; it
  was an orphan until that day.
