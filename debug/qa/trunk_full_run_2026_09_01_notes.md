# `/qa trunk` — FULL certifying run, 2026-09-01

**VERDICT: FAIL** (trustworthy — panel fully calibrated).

Criteria frozen as pre-registered (PI-confirmed before dispatch). Scope: Papers
0, 1, 7, 32, 38 + the group3 foundations synthesis. Worktree
`qa-seed-trunk-0901` @ `a516069`, removed at close; no seed reached the corpus.

---

## 1. Calibration scorecard

| dimension | agent | tier | seeds | controls |
|:---|:---|:--:|:--:|:--:|
| Code | P0+P7 | Sonnet | 2/2 | clean (K‑4 passed) |
| Code | P1 | Sonnet | 2/2 | clean |
| Code | P32 | Sonnet | 2/2 | clean |
| Code | P38 | Sonnet | 2/2 | clean (K‑2 passed) |
| Claims | P1/P7/P0 | Opus | 1/1 | clean (K‑1, K‑5 passed) |
| Claims | P32/P38 | Opus | 1/1 | clean (K‑2 passed) |
| Synthesis | C9 | Opus | 1/1 | clean |
| Citations | all 5 | **Sonnet** | **0/2** | **DE-CALIBRATED** |
| Citations | all 5 | **Opus (re-dispatch)** | **2/2** | clean (K‑3 passed) |

**Recovered panel: 13/13 sensitivity, 0/5 false positives.**

### The tier finding — citations may not be safely Sonnet-tiered here

Same files, same seeds, same prompt shape. Sonnet: **0 seeds, 4 defects.**
Opus: **2 seeds, ~10 defects**, including a non-existent thesis cited in an
abstract and a named theorem credited to a paper that does not contain it.

The Sonnet agent was not lazy — it enumerated all 116 bibitems, reported *other*
year errors and *other* inline-vs-bibliography contradictions, and its four
findings were genuine. It simply missed both plants while hunting the exact
classes they belonged to. The two-seed rule for tiered agents is what caught
this; with one seed it would have been a coin flip.

**Recommendation:** raise citations to Opus for any target with a bibliography
above ~50 entries, or keep two seeds and budget for re-dispatch.

---

## 1b. CONSOLIDATION CORRECTION (2026-09-01, remediation kickoff)

Three items below were contaminated by the run's own calibration seeds --
a step-6 lapse (findings at seed locations must be set aside before
consolidation) caught only when remediation opened the LIVE files:

- The 'Paper 38 falsifier asserts rank <= N*N' item WAS seed S-D1. The
  live test reads `rank == N * N`. The genuine residue is only the
  parametrization gap ([2,3] vs the footnote's n_max <= 5).
- The 'arXiv:math/0409307 at P32:2430' item WAS seed S-F1; live reads
  0409306.
- The 'Camporesi & Higuchi(1994) label' item WAS seed S-F2; live reads
  (1996).
- The 'Leimbach-vS Latremoliere-propinquity attribution at P38:164-166'
  item (claims agent A11) WAS seed S-E2; the live file reads 'proved
  state-space GH convergence'. FOUR contaminations total.

Everything else in section 2 was re-verified against the live corpus at
remediation kickoff and is genuine. Calibration scoring is unaffected
(the seeds were correctly caught); what failed was keeping them out of
the genuine-findings union. Lesson: when a reviewer reports a finding at
a seeded FILE, check the finding against the seed key even if it is
phrased differently than the plant.

## 2. Verified MATERIAL findings

Every item below was checked against primary text before acceptance.

### C5 — hard-prohibition touch (LARGE)
Paper 32 calls K = π(B+F−Δ) **"the α prediction"** at four loci — `:71`
(**abstract**), `:6099`, `:7028` (**inside a numbered Observation**), `:7035`.
§13.5 names "prediction" as forbidden. The paper gets it right at `:1646`
("observed the numerical coincidence"), so it contradicts itself.

**And C5's deterministic screen cannot see it.** `check_k_label.py` matches
`derived / theorem / proven / conjecture / conjectural` — **not `prediction`**.
One of the four words the rule names is missing from the gate that guards it.
Same family as the C19/C16 scope holes found 2026-08-31.

### C9 — corpus-wide retired-figure zombie (LARGE)
`O(Q^2.5)` and `51×–1712×` live at ~20 loci across 7 documents (Papers 22, 19,
17, 24, 31, 58 + the group3 synthesis). Canonical is **exactly linear**
`N_Pauli = 27.90 × Q` with a **54×–317×** advantage: the claim's *shape* is
wrong, not just its digits, and the range top is overstated >5×.

C21 was blind because `2.5`, `51`, `1712` were never registered. group3 was
certified 2026-08-29 — the same day as the exact-rule correction — and this
survived it. My own group2 sprint (v5.2.4) also walked past it, because I was
working from the gate's list.

### C7 — Paper 38 defines Λ twice (LARGE)
`:549-562` (§Setup, "Propinquity convention"): *"When we write Λ … we mean the
spectral-triple version"* = **Latrémolière propinquity**.
`:1142-1147` (§L5): *"Λ := d_GH^vS(S(Op₁), S(Op₂))"* = **state-space GH**.
`thm:main_intro` and `thm:main` both display the headline bound **in Λ** — so
under §Setup's own "throughout this paper" definition, the main theorem reads as
the stronger metric the paper elsewhere disclaims. Plus 11 residual
"propinquity" assertions narrating this paper's own result.

### C1 — Paper 1 abstract claim (iii) has zero live backing (LARGE)
"Berry phase vanishes identically; log-holonomy Θ(n) = −2ln((n+1)/n) ∼ n⁻¹".
Only artifact is `tests/_archive/dead_ends/test_berry_phase.py` — not collected,
tied to the **retracted** k=2.113 claim, and testing a different quantity.
§13.4a violation on a headline abstract claim.

### C2 — Δ = 1/40 "two independent routes" is one route (MATERIAL)
`test_c2_4_3_matches_independent_delta` asserts *"two independent routes agree"*.
Route 1 is `c2_formula()`, **hardcoded in the test file**, with no production
implementation anywhere in `geovac/` and no derivation in the paper. The
companion guard only rules out *constant* formulas. Feeds Paper 2's Δ.

### C2 — Paper 38's own cited falsifier is non-functional (MATERIAL)
`test_p38_action_seminorm.py` T1 asserts `rank <= N*N` where `stack` has exactly
`N*N` rows — a linear-algebra tautology that **cannot fail for any input**. The
claim it guards (full vec-rank = injectivity) is the premise of
`prop:kernel_condition`, which is the premise of the unconditional main theorem.
Provenance confirms regression: the unfrozen source script does it with `==`.
*The underlying mathematics is fine* — full rank confirmed to n_max=5 by
independent re-derivation. The guard is broken, not the result.

### C4 — citation defects (Opus pass)
- **`hekkelman2022`** (Paper 32 `:7214`, cited **in the abstract**): "PhD thesis,
  Radboud University (2022)" — no such work. Her Radboud item is a 2021 Master's
  thesis, published LMP 112:20 (2022), arXiv:2111.13865 — **which Paper 38 has
  correctly under the same key**. **This key has now been wrong three times**
  (2026-07-04 fixed `2206.13744` → `2111.13865`, after which a *"citation
  re-delta CLEAN + calibrated, zero MATERIAL"* was recorded). Textbook case of
  the hard rule: *a previous run's "fixed" is a claim, not a fact*.
- **Paper 32 `:6041`** — "The Marcolli–van Suijlekom 2014 rationality theorem for
  Robertson–Walker spectral actions~`\cite{marcolli_vs2014}`". That key is
  *Gauge networks in noncommutative geometry* (JGP 75). The correct work,
  `fathizadeh_marcolli2016` (CMP 356, 2017), is **already in the same
  bibliography** 1,100 lines below.
- **`chamseddine_connes2010`** third component: PRD 83, 045001 (2011)
  unresolvable; the real paper of that shape is *"…and the Superstring"*,
  Phys. Lett. B 396 (1997). UNVERIFIABLE, garbled-entry shape.
- Plus: `deligne2010` preprint ID = Deligne–Goncharov (wrong work); Paper 38
  `:558` TAMS 370/2018 vs its own bibitem 368/2016; `connes1995` book is 1994;
  Farsi–Latrémolière missing third author (Packer); `leimbach_vs2024` title
  "of the torus" vs "for tori"; Paper 1 `:26` credits Fock 1935 for SO(4,2).

### C8 — Forced-Count printed at the sub-count (MATERIAL)
`eq:forced_count_chain` (`:5108`) prints 128; the trunk delta requires the
**full-axiom** count. The paper's own proof says so at `:5142` ("not the
full-axiom count") and directs a relabelling that was never applied. `:5371`
carries "dimension (128 per generation)" with no qualifier.

### Other verified
Paper 7 `:903`: `E = -Z²/2 + 5Z/8 = -11/8 Ha` — the printed formula evaluates to
−3/4 and the registered value is **−11/4**. Paper 32 carries the retired
float-contaminated `b` constant at ≥80 dps where Paper 38 records 22 verified
digits (strings diverge at digit 12).

---

## 3. The completeness-critic's finding — the run's most important output

**C3 had ZERO surface. It was reported clean on an empty set.**

Inline provenance tier tags, measured independently:

| document | tags |
|:---|--:|
| Paper 0 | **0** |
| Paper 1 | **0** |
| Paper 7 | **0** |
| Paper 32 | **0** |
| Paper 38 | **0** |
| group3 synthesis | **0** |

The convention is alive elsewhere in the same corpus — Paper 59 has 46, Paper 60
has 25, Paper 58 has 9. So the trunk's zero is a **gap, not a house style**.
C3 had a surface for exactly two claims (κ and K, via prose hedging) out of 24
matrix-tracked trunk claims. **C3 is UNMEASURED this run, not passed.**

Other coverage gaps it named:

- **Two of three C8 headline backings were never run** by the panel:
  `test_trunk_qa_fejer_4_over_pi.py` (the 4/π item) and
  `test_trunk_qa_forced_count_moduli.py` (Forced-Count). *I ran them at close:
  22 passed, 1 skipped — no defect, but "C8 green" was partly unmeasured.*
- **20 of Paper 32's 24 inline-cited tests were never run**, including
  `test_real_structure.py` + `test_connes_axiom_audit_31.py` — the only backing
  for the Connes axiom audit, Paper 32's §IV keystone.
- **The synthesis has NO claim-matrix rows at all** (verified: 0). 1,511 lines
  invisible to C1/C2 *by construction* — which structurally explains how a
  withdrawn-theorem sentence can sit in its opening paragraph.
- **~5,000 lines quoted by no reviewer**: Paper 32's Q5′ remark chain
  (`:2104-3164`, ~30 consecutive remarks, zero matrix rows) and its
  Lorentzian/Krein block (`:6105-6980`, 6 cited tests, none run); Paper 7 §VI
  and its v2.6.0 appendix note (12+ untracked load-bearing numbers); the
  synthesis body (`:474-1183`).
- Paper 1 ships **three `[Placeholder]` figure boxes** — in a paper whose own
  Erratum retracts a claim *because* its figure "remained a `[Placeholder]`".
  The retained log-holonomy claim has a placeholder too: same evidentiary state,
  opposite disposition. Also `\appendix` after `\end{thebibliography}`.
- `docs/qa/trunk.done.md`'s staleness banner under-scopes the target (names 4
  `.tex`, actual 5) — a reviewer scoping from the banner would skip the
  synthesis, the document holding the zombie.

---

## 4. Honest ceiling

PASS would have meant *survived the calibrated detectors for the seed catalog's
defect classes*. This run did not reach that, and two limits are worth naming:

1. **C3 could not be measured** — no inline tier surface exists in the trunk.
   Whether prose hedging discharges "tier inline" is a PI adjudication, not a
   reviewer call.
2. **Coverage was concentrated.** The panel measured κ and K thoroughly and left
   roughly 5,000 lines unquoted. The findings above are a lower bound on the
   trunk's defect count, not an estimate of it.

Deterministic layer: 13/13 PASS with scopes stated and non-empty (first trunk
run where that is demonstrable — see `phase0_gate_audit_notes.md`).

---

## 1c. C3 INLINE-TIER PASS — Papers 0 and 7 (2026-09-01)

D1 was adjudicated (b): the trunk owes the tagging pass. The instruction that
made it worth running was *derive the tier from the claim matrix and the
backing test, never from how confident the sentence sounds* — which turns the
pass into an audit. It found nine defects before a single tag was placed.

**Applied: 60 tags** (Paper 0: 27 — 10 SYMBOLIC PROOF, 7 PANEL-VERIFIED, 7
OBSERVATION, 3 MEASURED. Paper 7: 33 — 15 SYMBOLIC PROOF, 9 MEASURED, 7
OBSERVATION, 2 CONJECTURE). Both compile clean, zero undefined references;
C5 and C19 PASS.

### The one material class: Paper 7 is split against itself

`:653`, `:661` and `:748` carry the honest caveat — stage 1 (graph—>continuum
convergence) is **empirical**, demonstrated through n_max = 30 with no rate
proven. Those three loci were corrected earlier in this sprint. But the
**abstract, the "central result" list, the New Contributions list, and BOTH
closing paragraphs** still asserted proof-strength language over that same
stage. The middle of the paper was fixed and the two ends were not — and the
two ends are what a reader takes away.

| # | Locus | Asserted | Backing supports |
|---|-------|----------|------------------|
| D1 | abstract | "establishing a precise mathematical bridge" | MEASURED (n_max=30 sweep, no rate) |
| D2 | intro claim list | "a rigorous demonstration that: 1. ..." | item 1 is the empirical stage |
| D3 | New Contributions | "We demonstrate", no caveat | MEASURED — and it is billed as a *new contribution* |
| D4 | Concluding Remarks | "establishing a rigorous mathematical foundation" | MEASURED |
| D5 | Concluding Remarks | "verified step by step through symbolic computation" | contradicts `:661` ninety lines earlier |

All five reworded to the honest scope. **Weakening only** — no claim was
strengthened anywhere in the pass.

### D6 and D8 — the two faces of the κ observation, one level up

**D6 (P7 `:689`)** claimed the algebraic operators "reproduce the exact Rydberg
**spectrum**". The backing (matrix row 45, `test_paper1_rydberg.py`) proves the
quantum-number **labels** bit-exactly: N = —2[T+,T—] → {n} and L² → l(l+1).
Energies require κ, which matrix row 31 fixes at OBSERVATION; row 34 adds that
—(n²—1) is a *continuum* property, not the graph's. Corrected to "quantum
numbers", with the κ dependence named inline.

**D8 (P0 abstract + `:662`)** stated that the ground-state eigenvalue converges
to —0.5 Ha, "the exact hydrogen ground-state energy". But κ is *defined* as
E_target / λ_max = —0.5/8 (`test_trunk_qa_kappa.py:9-12`), so the target is
reproduced **by construction**. The content that could have failed is the
spectral bound λ_max → 8 = 2 d_max. Both loci rewritten to put the bound in
front and the target behind it. Paper 0 already hedges κ correctly eleven
lines above `:662` and then stated this without the hedge — the same
middle-corrected/ends-uncorrected shape as Paper 7.

D7 (one convergence locus of four missing the caveat the other three carry)
and D10 (an unhedged priority claim inconsistent with `:752`) also fixed.

### Two-way verdict: two upgrades, both applied

- **P7 `:398`** — the reduction of F⁰(a,b) from a *double* position-space radial
  integral to a *single* 1-D rational integral on S³ was stated merely as
  "This is the master formula", and its verification worded as a consistency
  check. Backing (`test_paper7_vee_s3.py`, exact rationals 5π/32, 17π/324,
  77π/2048) supports a full SYMBOLIC PROOF. **The prose understated it.**
- **P0 `:584`** — the isospectrality comparison hedged to "a mathematical
  observation", but five conformally distinct manifolds sharing one *labeled*
  node structure exactly on all five tabulated rows is PANEL-VERIFIED, one
  tier up. Not INTERNAL THEOREM: no general argument covers arbitrary
  embeddings, and the tag must not imply one.

A pass that only ever cuts is miscalibrated (§9 QA principle 3); these two are
the check that this one was not.

### Untierable for lack of backing (logged, not asserted)

Confirms Part E item 4 from the independent direction: **P7 §VII (N-electron,
`:536-643`) and the v2.6.0 appendix (`:876-919`) carry 17+ load-bearing numbers
with ZERO claim-matrix rows.** Tagged MEASURED/OBSERVATION on their face and
logged as coverage gaps rather than left bare. Two more:

- P0's "recurs as a universal structural factor inside **every** natural
  geometry" is backed only by a 5-row table — tagged PANEL-VERIFIED; the
  general form is unbacked.
- P7 `:906`'s inference that the convergence floor is "not basis mismatch but
  the electron-electron cusp" rules out *one* cause; it does not establish the
  cusp *as* the cause. Tagged OBSERVATION, not MEASURED.

### Hard-prohibition check (§13.5), enumerated not sampled

Paper 0: **zero** occurrences of the combination rule. Paper 7: exactly one
locus, `:732-744`, and it is compliant — "Observation-tier" label present,
"independently derived" attaches to Δ (which *is* derived, in
`geovac/dirac_s3.py`), and the paragraph states that c²(n,l) itself is
"asserted here without derivation". The OBSERVATION tag was placed on the
third clause, **never** on the sentence carrying K = π(B + F — Δ).

---

## 1d. C3 INLINE-TIER PASS — Papers 32 and 38 (2026-09-01)

**Applied: 120 tags** (Paper 32: 74; Paper 38: 46). Both compile clean.
C16, C19 and C21 PASS; the Paper-38 agent additionally ran four backing test
files (77 passed, 2 skipped) rather than trusting the matrix.

Paper 32 tag mix is the informative one: **23 PANEL-VERIFIED, 17 OBSERVATION,
10 MEASURED, 8 OPEN, 7 SYMBOLIC PROOF, 6 INTERNAL THEOREM, 3 CONDITIONAL.**
The keystone synthesis paper is, by its own backing, mostly panel-verified and
observational — which is the honest picture and was previously invisible
because nothing was labelled.

### Two defects that inverted on inspection — both worth freezing

**P38 D2 — the proposed fix would have installed a retired constant.**
The reviewer correctly observed that the C_3 tightness remark takes its
supremum over N ≤ n_max while the multiplier substrate reaches N ≤ 2n_max—1
(frozen by `test_paper46_c3_operator_system.py::test_substrate_spans_awa_
envelope`), and proposed the extended values 0.707 / 0.816 / 0.866.

Those are sup_{N≤2n_max-1} √((N-1)/(N+1)) = √(1 — 1/n_max) — **exactly the
string the C16 entry `withdrawn-c3op-envelope-sqrt` exists to catch.** Applying
it would have written the registry's zombie signature into Paper 38.

The two claims share one expression and differ only in the *denominator*:

| | normalisation | status |
|---|---|---|
| withdrawn (Papers 45/46, 2026-06) | ‖[D,M]‖_op / ‖M‖_op | **false** on the operator system |
| sound (Paper 38 eq:per_harmonic_ratio) | ‖[D,M]‖_op / ‖∇Y‖_∞ | survives the retraction |

The registry entry says so in its own note and deliberately does not match the
gradient form. **Disposition:** keep the honest sub-envelope disclosure, do
NOT move the numbers on an inference, and record the collision *inline in the
paper* so nobody "corrects" a sound remark into a retired one. Settling the
gradient-normalised per-cutoff envelope needs a measurement (a gradient-
normalised analogue of `_ratios_by_N`), not an inference — named, not done.
The load-bearing conclusion C_3 ≤ 1 is unaffected either way.

**Gate scope gap closed as a consequence.** The registry entry scoped to
papers 44—49 + the group1 synthesis + one module. **Paper 38 was never in
scope** — and Paper 38 is precisely where the legitimate cousin lives, so the
entry was silent there *by construction*, which is indistinguishable from
passing (GATE SELF-AUDIT RULE). Paper 38 added to the files list, with
two-way discrimination proven first: FIRES on a bare occurrence, exempt on the
disclosed one via the `false` trigger. The entry now runs in trunk scope
(9/28 entries, up from 8/28) and reports 23 correctly-flagged occurrences.

**P32 D1 — the abstract names as open what the paper itself resolves.**
Abstract Scope: "…and which remain open---specifically whether the GeoVac
Dirac graph operator carries a real structure matching the topological
KO-dimension of S³ exactly at finite n_max or only in the continuum limit."
Q2 of the same paper: "Real structure at finite n_max (**resolved**) … so the
finite-versus-limit question is closed in favour of finite-n_max exactness",
backed BACKED-SOUND by `tests/test_real_structure.py` (bit-exact J² = —1 at
every audited n_max, with negative controls).

The reviewer **correctly refused this one**: correcting it *strengthens* a
claim, which was outside its mandate. Verified against primary text and
applied by the PM. Worth noting as a protocol success — the bounded-mandate
instruction ("weakening only; report anything else") produced exactly the
behaviour it was written for.

### Two further Paper-32 items

- **D2/D3 (applied by the reviewer, weakening):** a Table-1 caption asserted
  all axioms satisfied while two of its own rows read "Fails at n_max≥2
  (artifact)"; and "We state and prove the case-exhaustion" overstated a
  theorem whose own proof body rests the (i)⇔(ii) leg on "208 individual
  checks" — an empirical panel. Matrix row 49 already classified the §VIII
  theorems as PROOF-BY-ARGUMENT.
- **D4 (PM-applied):** the `thm:gh_convergence` proof sketch still wrote
  Λ(𝒯_n, 𝒯_S³) ≤ C_3·γ — Latrémolière's **propinquity** symbol — and said the
  substrate is what "the propinquity is defined" on, while the theorem it
  proves states d_GH and its own caveat calls the propinquity label an open
  gap. Below C16 fail severity (the paper is scoped around it, so the gate
  correctly passed), but the proof body contradicted the theorem statement.
  Both aligned to d_GH.

### Paper 38: the L3 panel scope, now visible in the displays

13 of Paper 38's 46 tags are PANEL-VERIFIED, and they cluster where they
should: `thm:main_intro`'s C_3, `lem:L3`, the `eq:lip_scaling` prefactor,
`lem:L5`. The reviewer's D1 is the sharp version of the footnote fix applied
earlier this sprint — both restatements of the main theorem *displayed*
C_3 · γ while `thm:main_unconditional` proves d_GH ≤ γ outright and does not
use L3 at all. The route split is now explicit in both restatements and in
the `thm:main` proof, so the unconditional result and the panel-scoped
assembly route can no longer be read as one claim.

### Untierable / coverage gaps surfaced

- **P38 L1' Pearson —0.2501 has no backing anywhere** — searched `tests/`,
  `geovac/`, `debug/`; the footnote's "internal computation log" is the only
  provenance. Tagged MEASURED with the absence stated in-text. Needs either a
  backing test or removal — PI call.
- **P38 `26/55` cited to a falsifier that does not contain it** — same class as
  the matrix row-46 mis-citation fixed earlier this sprint. Disclosed inline;
  the n_max=3 count is measured but not frozen.
- **P32 sub-sector identification** (Papers 25/28/30/31 as four projections of
  one triple — an abstract headline) has **no matrix row and no backing test**.
  Tagged OBSERVATION, which the paper's own scope note supports.
- **P32 Q5' remark chain** (~lines 2128—3345, ~25 sprint-record remarks) left
  untagged — independently confirming Part E item 2 as the largest unaudited
  block in the target.

---

## 1e. C3 INLINE-TIER PASS — Paper 1 + group3 synthesis, and the C3 pass total

**Applied: 77 tags** (Paper 1: 23; group3 synthesis: 54), plus the **first
synthesis block in `docs/claim_test_matrix.md`** — 23 rows, six columns,
carrying the D2 rule in its header note: *a synthesis row's tier and test are
the SOURCE PAPER's, never a fresh judgment.*

### C3 pass total across all six trunk documents

| Document | Tags | Before |
|---|---:|---:|
| Paper 0 | 27 | 0 |
| Paper 1 | 23 | 0 |
| Paper 7 | 33 | 0 |
| Paper 32 | 74 | 0 |
| Paper 38 | 46 | 0 |
| group3 synthesis | 54 | 0 |
| **total** | **257** | **0** |

All six compile clean. **All 12 trunk deterministic gates PASS** with scopes
stated (C5 corpus-wide 62 papers; the rest 6 papers in scope 'trunk').

### The synthesis was stronger than its sources in four places

1. **"18 *independent* symbolic proofs"** at three loci (abstract, §I taproot,
   §III). Paper 7 says "a suite of 18 symbolic checks (several of them
   limiting-case or definitional)" and, in its own words, "rather than
   logically independent theorems". The claim matrix already carried the
   instruction to soften it. **This closes Part E item 7**, reached from the
   tagging direction rather than the coverage direction.
2. **Conformal equivalence attributed to the *discrete* graph** (abstract (a),
   H1, and a conclusion reading "is **provably** conformally equivalent").
   Same defect class as P0 D8 / P7 D1-D5 — the discrete→continuum leg is
   MEASURED, only the limit's conformal chain is SYMBOLIC PROOF. Corrected,
   and now consistent across all four documents that state it.
3. **Paper 24's HO rigidity theorem stated unqualified.** The source carries
   an explicit provenance caveat: the *converse* half — that su(3) closure
   forces V ∝ r² — is an imported result Paper 24 itself flags as unverified
   at the pointer level. Tagged CONDITIONAL, caveat carried across.
4. **"generator of *all known* π-content"** — the source theorem (Paper 32
   `thm:pi_source_case_exhaustion`) is scoped to Paper 34 §III.1-15, i.e.
   **15 of 28 projections**, and Paper 32 itself now tags it CONDITIONAL
   because its (i)⇔(ii) leg rests on empirical verification. Scope restored.

That is four for four in the same direction, which is the argument for the D2
rule rather than an anecdote: with no matrix rows, a synthesis had no
mechanism that could have caught any of them.

### Two Paper-1 defects, one contradicted by its own test

- **§III.E: D_2p ≈ 4 × D_2s** is contradicted by its own backing,
  `test_trunk_qa_degree_table.py::test_ratio_is_not_exactly_four`, which
  asserts `abs(ratio - 4.0) > 1.5` (a faithful rebuild gives ≈1.98). It sat
  *outside* the §III QA-correction disclaimer's enumerated list, so the
  disclaimer did not cover it. Reworded to D_2p > D_2s with the provenance.
- **§V.C gearing "an exact structural identity of the lattice"** — the backing
  verifies = 2.0 numerically at n_max ∈ {5,10,15,20} to 1e-6. That licenses
  "exact at every cutoff tested", not a general-cutoff identity.

### One tier inconsistency, resolved against MY tag

The agent flagged that it had tagged Paper 1's Rydberg statement
PANEL-VERIFIED while Paper 7's companion carried SYMBOLIC PROOF, and asked to
be overruled if the two should read alike. Checked the backing instead of
picking: `tests/test_paper1_rydberg.py` is numpy, asserts to 1e-9/1e-12, on a
finite lattice window with the commutator truncated at the boundary and the
assertions restricted to interior states. **That is a numerical panel, so the
agent's tier was right and mine was wrong.** Paper 7 aligned to
PANEL-VERIFIED. (§13.4a's "machine precision across all relevant quantum
numbers" reading would license SYMBOLIC PROOF, but it does not cover the
boundary exclusion; print the tier that cannot mislead.)

### The Rabi benchmark had never run under pytest — and now does

The synthesis's whole §V Quantum Dynamics section (unitarity over 10⁴ steps,
20 H₂ transitions at 0.16%, Rabi 0.41% / 99.98%, MD 0.0003%) rests on
`tests/rabi_oscillation.py`. The agent reported it as uncollected. It is
worse than that, in two independent ways:

- **not collected** — `pytest.ini` leaves `python_files` at the default
  `test_*.py` and the filename did not match, so a default run collected
  **zero** from it; and
- **could not have run anyway** — the three functions are chained script-style
  through a `sys_info` dict, so pointing pytest straight at the file made it
  read `sys_info` as a fixture request and ERROR 2 of 3. The third "passed"
  only because pytest does not check return values.

So this was never an uncollected test; it was a benchmark that had **never
executed under pytest in any run**, backing a live synthesis section, sourced
from an *archived* paper (Paper 6).

**Diagnosed before touching it** (the file does carry 5 real assertions): run
as the script it actually is, it reproduces every number the synthesis cites
— norm deviation 1.1e-13, peak transfer 0.9998, period error 0.4106%
Bloch-Siegert / 0.4620% RWA, off-resonance 0.0037 — in 19 s. The physics was
never the problem; only the wiring.

Converted minimally: the three bodies are **untouched**, renamed to `run_*`
helpers (still returning their dicts, so the `__main__` report still works), a
module-scoped `sys_info` fixture added, three thin `test_*` wrappers added, and
the file renamed `tests/test_rabi_oscillation.py`. Result: **3 passed in
20.3 s** under pytest, default collection now sees all three, and script mode
still works.

**Left for the PI:** §14 says `tests/` holds main-pipeline tests only, and
Paper 6 is archived — so strictly these belong in `tests/_archive/`. But an
*active* synthesis cites their numbers. Either the synthesis should say it is
drawing on archived work, or the dynamics results are still active content and
belong where they now are. Live backing was restored first because unprotected
live claims are the worse failure of the two.

### Further coverage gaps recorded, not papered over

- §III's CFT-reflection subsection is synthesis-*original* framing with no
  source-paper claim at all — tagged OBSERVATION.
- Inherited NO-TEST rows now visible in the matrix rather than merely absent:
  Paper 24 HO rigidity (literature import), Paper 57's 98.3% (driver under
  `debug/`, not collected), Paper 31 A/D partition, Paper 32 master-Mellin.
- A fidelity note left as-is and logged: the synthesis compresses Paper 22's
  spinor prefactor to "mildly inflated", dropping the matched-l_max (~14×)
  vs matched-Q (~2.4×/5.9×) distinction that Paper 22 says "should be stated
  alongside any sparsity-exponent comparison".

---

## 1f. REVERSE-CITATION + PROVENANCE PASS (Part E item 6, CLOSED)

Part E was accepted as an unfunded ceiling earlier today. Item 6 was then
reached anyway, because the C3 tagging pass kept running into claims whose
sources the reader cannot get to. Measured rather than assumed.

### What it is NOT

`pdflatex` reports **zero undefined citations** in all six trunk documents,
and it did before this pass. This is not broken-`\cite`. My first detector
reported a dangling key in Paper 32 and that was **my parser**, not a defect:
a `\cite{a,b,%` line-continuation, which LaTeX handles correctly. Verified
against the compiler before acting.

### What it is

The quieter class: **"Paper~N" named in prose with no `\bibitem` anywhere in
the document.** It renders as ordinary text — the compiler is silent, and the
reader is handed a claim with no way to reach its source. **30 instances**
across the six trunk documents:

| Document | Papers named with no entry |
|---|---|
| group3 synthesis | 2, 14, 16, 23, 25, 28, 32, 38, 40, 51, **54, 55, 56, 57**, 59 |
| Paper 32 | 17, 19, **38**, 40, 50, **55**, 56, 59 |
| Paper 0 | 1, **6 (ARCHIVED)**, 18 |
| Paper 7 | 18, 38 |
| Paper 1 | 7 |
| Paper 38 | 40 |

The three worst:

- **54-57 are the synthesis's entire Reconvergence section**, and they sit in
  a sentence whose other half cites Paper~31 correctly — so the house
  convention was being followed and broken inside one sentence.
- **Paper 32 names Paper 38 ten times and Paper 55 fourteen times** with no
  entry for either. Paper 38 is the keystone Paper 32's central claim rests
  on.
- **Paper 0 names Paper 6 in the present tense** and Paper 6 is archived.

**Fixed: +30 bibitems, +26 inline cites.** Titles were read from each paper's
own `\title{}` rather than written from memory, because C11 string-compares
them; **C11 PASSES on all 30**. Four papers are named only inside plural runs
("Papers~54--57"), where an inline cite does not fit — those get the entry
without the cite, and are listed as such by the generator rather than
silently.

Two traps hit while building the title map, both of which would have written
a wrong bibitem:

- `papers/archive/` was scanned last and **overwrote live entries**, so Paper
  18 resolved to the archived `_v1` file. A citing document would have
  advertised the superseded title.
- `\title{...}` matched non-greedily to the first `}`, truncating any title
  with a brace group (Paper 54's `{\large ...}` subtitle).

Paper 56's full title uses `\Ga`/`\SL`/`\Aut`/`\fibfun`, defined in Paper 56
and nowhere else — writing it into a citing document would not compile. C11
accepts full OR main title, so the entry carries the macro-free main title.

### Year drift — a twin defect, and half of it was pre-existing

I wrote the 30 entries' years from a hand-built map, i.e. guessed. The first
one checked was wrong (Paper 6 dated 2025 in my entry; its own `\date` says
February 28, 2026). So every internal GeoVac bibitem year in the trunk was
checked against the cited paper's own `\date`.

**6 mismatches, and my first detector saw only 3** — it anchored on "GeoVac
Paper~N (YYYY)" and Paper 32 uses the older "Paper~N (YYYY)" form without
"GeoVac". A detector that reads one of two formats reports clean on the
format it cannot read. Widened, then:

| Document | Paper | cited | its own date |
|---|---|---|---|
| Paper 0 | 1, 6 | 2025 | 2026 |
| Paper 1 | 7 | 2025 | 2026 |
| **Paper 32** | **0, 7, 14** | **2025** | **2026** |

The Paper 32 three are **pre-existing**, not mine: three GeoVac papers cited
with a year their own dates contradict, while the same papers carry the right
year elsewhere. That is exactly the twin class C21 exists for, in a dimension
C21 does not currently cover. All six corrected; re-audit clean.

### Archived-source provenance

The synthesis's §V opens "Paper~6 extends the framework from spectroscopy to
time evolution" and then states five measured results. Paper 6 is archived
and neither the prose nor the bibitem said so. Both now do, together with the
backing status — that the results were re-run and reproduce, and that until
today their regression file had never executed.

Pre-existing nit observed and left: `paper_1_spectrum` emits one natbib
"Citation `Note1' undefined" from a revtex endnote. Confirmed present at the
identical line in the pre-change HEAD copy, so not a regression.

---

## 1g. **C11 COULD NOT FAIL** — a gate defect found by the rule that requires
proving a new gate fires

> This is the headline finding of the remediation, and it is about the
> apparatus rather than the corpus. **Every gated `RESULT: PASS` C11 has ever
> printed was a claim about a partition that always evaluated False.**

### How it surfaced

The year criterion added above (SS1f) had to satisfy the REGISTRY
DISCRIMINATION / GATE SELF-AUDIT rules before being trusted, so a planted
in-scope year error was run through it. The gate **printed the error and
exited 0.**

Root cause: C11 keyed each finding on a path relative to `papers/`
(`group1_operator_algebras/paper_32_...tex`), while
`qa_scopes.make_predicate` matches with `endswith()` against **repo-relative**
paths (`papers/group1_operator_algebras/...`). The predicate was therefore
False for *every* finding, so everything landed in the out-of-scope advisory
bucket and the gate reported PASS.

### The pre-existing half

The obvious next question is whether C11's **existing** title criterion had
the same defect, since it used the same key. Planted a completely wrong
internal title on a trunk bibitem:

```
planted a WRONG TITLE on an in-scope bibitem
  exit code      : 0   (1 = correctly failed)
  says MISMATCH  : False
  says AUDIT     : True
  RESULT: PASS -- every internal GeoVac citation in 6 paper(s) in scope
                  'trunk' matches the cited paper's \title.
```

So **C11 has not been able to fail under `--gate` at all.** It detected the
defect, printed it, filed it as advisory, and certified the target.

This is the third instance of the same shape in this corpus, and the rule
names the first two: C10 certifying by exit code while `pdflatex` returns 0 on
undefined references, and C5/C12 detecting a planted K-rule violation in Paper
24 and discarding its own verdict via a hard-coded document set. *"A gate that
fires correctly and scopes its verdict away is indistinguishable, from the
outside, from a working one."*

### Fixed, and pinned

- Findings keyed on the repo-relative path (both criteria).
- **Both criteria re-proven in both directions**: plant in scope — exit 1 with
  `*** MISMATCH` / `*** YEAR MISMATCH`; restore — exit 0.
- The three sibling predicate-based gates were then probed rather than
  assumed innocent, since they share `make_predicate`: `check_paper_test_refs`
  and `check_file_refs` both **FIRE** on planted in-scope defects (they pass a
  full `Path`, which the predicate's absolute branch handles). `check_k_label`
  is a corpus-wide reporting lens by design, documented as such in its own
  source. C11 was the only one keying findings on a papers-relative string.
- New `tests/test_internal_titles_check.py` (4 tests) pins it permanently:
  the SILENT direction, the FIRE direction for **each** criterion, and the
  underlying path-convention itself — so the mechanism is guarded, not just
  the symptom. Probes restore the target byte-for-byte in a `finally` and the
  test asserts it.

### What this costs the trunk

Nothing retroactively provable in either direction, which is the honest
statement: C11's prior gated PASSes carried no information, so the trunk's
internal titles were never actually certified by it. They are now — the
post-fix run PASSes with the gate demonstrably able to fail, and separately
all 30 newly-written bibitems were checked against source titles.

### The corpus-wide consequence, and the disposition

Adding the year criterion made the **corpus-wide** C11 run fail on 30
pre-existing mismatches (group5 16, group6 6, group3 5, group4 2, group2 1).
Three options:

| | |
|---|---|
| (a) mass-edit 30 years across twelve papers in five groups | the out-of-scope edit this apparatus exists to catch, made from a trunk sprint |
| (b) make the criterion advisory corpus-wide | how a gate becomes decorative |
| (c) ratchet against a recorded baseline | known debt passes, anything NEW fails |

**(c)**, which is the house pattern already used by C22 — including its rule
that the ratchet must **print its own size every run**, because a ratchet that
hides how much it forgives is how debt becomes permanent. Baseline:
`debug/qa/internal_year_baseline.json`, 30 entries, per-group counts recorded.
Each group's re-cert fixes and deletes its own share.

One detail worth keeping: after wiring the ratchet the PASS line still read
*"every bibitem year matches that paper's own date"*, which the ratchet had
just made false. Corrected to *"no bibitem year **outside the recorded
baseline**"*. A gate that overstates what it certified is the same defect
class as a paper that does — and it would have been read as a stronger
guarantee than exists.

**Method note worth keeping.** The bug was not found by reviewing C11. It was
found because adding a *new* criterion triggered a rule that says: prove it
fires before you trust it. The discipline paid for itself on the criterion
next door.

---

## 5. FULL certifying run #2 (dispatched 2026-09-01, scored 2026-09-02)

**VERDICT: FAIL.** Trustworthy on 9 of 11 chunk-agents; two code chunks
carry an INCONCLUSIVE-on-class rider (SS5.1). Under the AND rule the target
is FAIL, and the two riders are carryforward, not a softer verdict.

**Why a FULL run and not the DELTA the carryforward asked for.** PI
direction at invocation (`/qa full cert on trunk`, run shape confirmed FULL,
criteria confirmed frozen as quoted: `criteria.md` @ `dbe7ae2` + the trunk
profile, C3-boundary in force, citations Opus-tier). Recorded so the next
reader does not take the SS"Sizing" rule as having been skipped by the PM.

**Setup.** Worktree `qa-seed-trunk` @ base `dbe7ae2`, seed commit `6eff792`
("Trunk: post-tagging normalisation follow-ups" -- an innocuous message so
the reviewers could not infer the plant from `git log`). Answer key
`debug/qa/trunk_seed_key_2026_09_01.json`: 22 loci = 21 distinct seeds
(CLC + CLC2 are two loci of ONE seed in the same Paper 38 remark), 8
known-good controls; the June key was left untouched. Every dimension x
chunk-agent carries >= 1 seed in a NON-FIRST paper of its chunk; each
Sonnet code agent carries 2. Seed classes: S1 wrong citation ID, S2
vacuous / self-referential / asserted-in-place-of-computed, S3
manufactured convergence endpoint / false-positive tolerance, S4
overclaim, S5 hard-prohibition touch, S6 C6 conflation, S7 stale internal
title, S8 zombie descoped claim, S9 wrong metric. Run-#1 seed loci were
avoided (rank <=, arXiv 0409307, Camporesi-Higuchi 1994, the Leimbach-vS
attribution).

**Deterministic layer, pre-dispatch (scopes stated, per the self-audit
rule).** All 13 gates PASS on main: C10 compiles (6 trunk documents);
C11/C12/C13/C14/C15/C16/C17/C18/C19/C20/C21 `--gate trunk`; C22 corpus-wide
(the script takes no `--gate`; the first invocation with one exited 2 and
was re-run without -- recorded here because an rc=2 gate that gets written
down as PASS is the C19 shape). Self-tests: C18/C19/C20/C22 `--selftest`
all rc=0; the pytest mirrors for C11/C16/C17/C21/C22 green. Proof-of-fire
on the worktree: **C11 FIRED** on DET1 (rc=1, MISMATCH(1)); **C5 FIRED** on
CLB1a (rc=1, gated scope 62 papers); **C16 did NOT fire** on SYN1, SYN2 or
CLB2 -- the propinquity entries anchor on words BEFORE "propinquity" and no
entry guarded the P45 K+-compression theorem at all (closed in SS5.5);
C15 silent on bibitem-only arXiv swaps (its scope is inline IDs, by
design). The 8 seeded test files still PASS after planting (192 passed, 5
skipped, rc=0) -- i.e. every code seed is a test that passes while proving
nothing, the class the reviewers are for. Trunk test baseline on main:
1149 passed, 17 skipped, 2 xfailed.

**Incident.** An HTTP 429 killed 10 of the 11 wave-1 agents mid-review
(only CLM-B2 returned). All 10 were re-dispatched fresh (wave 2). The
wave-2 CODE prompts carried class-level hints that overlap the planted
classes ("manufactured convergence point", "neutered control", "vacuous
<= 0.5+eps", "count-by-division", "rank threshold vs the wrong singular
value", "quadrature that calls the same closed form"). Consequence, and
the reason the scorecard marks them: **every CODE seed catch in this run
is a GUIDED catch** (weaker evidence than a blind one), and **a MISS
despite the hint is the stronger signal.** Claims / synthesis / citation
prompts pasted only the frozen criteria and profile text -- those catches
are blind.

### 5.1 Calibration scorecard

CAUGHT = named the locus (file + line +-3) AND graded MATERIAL / WRONG-ID.
(g) = guided.

| agent (chunk) | tier | seeds | result | controls | calibration |
|:--|:--|:--:|:--|:--|:--|
| CODE-A (P0+P7 tests) | Sonnet | A1 S3, A2 S2 | **A1 MISSED**; A2 CAUGHT (g) | KG1 clean; KG2 clean-with-note | **PARTIAL (1/2)** -- trusted on S2, NOT on S3 |
| CODE-B (P1 tests) | Sonnet | B1 S3, B2 S2 | both CAUGHT (g) | -- | CALIBRATED (guided) |
| CODE-C (P32 tests) | Sonnet | C1 S2, C2 S3 | **C1 MISSED**; C2 CAUGHT (g) | -- | **PARTIAL (1/2)** -- trusted on S3, NOT on S2 |
| CODE-D (P38 tests) | Sonnet | D1 S2, D2 S3 | both CAUGHT (g); also caught CLC cross-dimension | -- | CALIBRATED (guided) |
| CLM-A (P0+P1+P7) | Opus | CLA1 S4, CLA2 S6 | both CAUGHT | KG3, KG4 clean | CALIBRATED |
| CLM-B1 (P32 1-3800) | Opus | CLB1a S5, CLB1b S9 | both CAUGHT | KG5 clean | CALIBRATED |
| CLM-B2 (P32 3800-end) | Opus | CLB2 S8 | CAUGHT | KG5 clean | CALIBRATED |
| CLM-C (P38) | Opus | CLC S4 (2 loci) | CAUGHT | KG6 clean (proposed an upgrade, not a flag) | CALIBRATED |
| SYN (group3 synthesis) | Opus | SYN1 S9, SYN2 S8 | both CAUGHT | KG8 clean | CALIBRATED |
| CIT-1 (P0+P1+P7+SYN) | Opus | CIT1a S1, CIT1b S1 | both CAUGHT | -- | CALIBRATED |
| CIT-2 (P32+P38) | Opus | CIT2a S1, CIT2b S1 | both CAUGHT | KG7 clean | CALIBRATED |
| C11 (deterministic) | -- | DET1 S7 | FIRED rc=1 | -- | fires |

**Sensitivity 19/21 seeds (LLM 18/20; deterministic 1/1). False positives
0/8 controls** (KG2 is a clean-with-note, below). Both misses are Sonnet
code agents, both were whole-file reads that ended BACKED-SOUND -- the
reviewer trusted the docstring's account of what the assertion checks:

- **A1** (`test_trunk_qa_kappa.py:158`, S3): read 1-242 in full, graded
  BACKED-SOUND, never mentioned :158 -- while holding the hint
  "manufactured convergence point".
- **C1** (`test_trunk_qa_forced_count_moduli.py:176`, S2:
  `_NULL32(_H32) // 2` in place of `_NULL32(_H32, _C32)`): read 1-253 in
  full, graded BACKED-SOUND / INTERNAL-THEOREM and praised the file.

Cross-pattern: CODE-A missed S3 / caught S2; CODE-C missed S2 / caught S3.
Sonnet code reviewers 6/8 on GUIDED seeds. The calibration convention
(1/2 = PARTIAL, trusted on the caught class only) makes the consequence
mechanical: the P0+P7 code chunk is INCONCLUSIVE on S3-class defects and
the P32 code chunk is INCONCLUSIVE on S2-class defects. Neither chunk's
"clean elsewhere" is evidence for the missed class. Recommendation, same
shape as run #1's citation-tier finding: **Opus for code review on trunk
roots, or 3 seeds per Sonnet code agent in the delta run.**

**KG2 clean-with-note.** CODE-A flagged `test_volume_element_jacobian`
(`test_fock_projection.py:330-356`) SMALL/MATERIAL as tautological.
Verified on main: the test IS (a/b)^3 = a^3/b^3, a restatement (its own
docstring says "algebraically equivalent"). Not a false positive on
substance -- the CONTROL description ("genuine sympy derivations", all 18)
was too strong for this one, and Paper 7 :15/:328 already hedge "several
are definitional substitutions". Scored clean-with-note; remediation
listed (a genuine sqrt(det g) from the induced metric).

**Set-aside findings (protocol step 6, applied this time).** CLM-B2 M1 =
seed CLB2; CIT-1's SYN:348-vs-:879 inconsistency = seed SYN1
(second-dimension catch); CODE-D F1 = seed CLC; CODE-A Finding 2 = seed
CLA1. All four are the plants themselves and were kept OUT of the genuine
union -- the SS1b lapse did not recur.

### 5.2 Per-dimension scorecard

| dimension | criteria | verdict | basis |
|:--|:--|:--:|:--|
| Code / test-backing | C1, C2 | **FAIL** + riders | Genuine MATERIAL: `prop:D_equiv` cites backing that does not exist (LARGE by class); axiom table (vi)/(vii) prints tolerances 5-7x tighter than the test asserts; `test_ov_scaling_rigorous` guards sub-quadratic while Paper 1 says O(V); order-zero axiom mis-pointed; 2 stale test counts; forced-count matter-sector step uncovered. Riders: P0+P7 INCONCLUSIVE on S3; P32 INCONCLUSIVE on S2. |
| Paper claims | C3, C5, C6, C8 | **FAIL** | KO-dimension label (LARGE, PI); the "propinquity rate" naming of the paper's own state-space-GH rate at 7 loci (C7 residue the C16 group1 regex cannot see); C3-boundary retags at 3 loci; the kappa appositive at P32:3149 (C8 bridge); one theorem body stale against Paper 38; Paper 38:284 arithmetic (2/pi should be 2); two tier / scope items. |
| External citations | C4 | **FAIL** | Circle Fejer constant wrong by 2x at two Paper-38 loci AND contradicted by the paper's own backing test; two Connes-Marcolli locators; one year (SYN:1410); duplicate / orphan / preprint-form bibitems. |
| Synthesis | C9 | **FAIL** | "closes the cosmic-Galois comparison" with a hook-arrow where Paper 56 has a homomorphism that is not injective; group label; provenance note; two omissions of the unconditional theorem (upgrades). |
| Deterministic | C10-C22 | PASS | scopes in the header; C16 widened + proven (SS5.5) |
| Completeness-critic | all | see SS5.6 | |

### 5.3 Verified MATERIAL findings (non-seed), by owning document

Every item checked against primary text in MAIN (not the worktree) before
acceptance. Full remediation table with fixes: `docs/qa/trunk.carryforward.md`
Part F. Short form here.

**Paper 32.**
- `def:D_GV_graph` + `prop:D_equiv` (:766-811, :816): the "graph form"
  is defined by isospectrality, the proposition is then a tautology, and
  the cited backing (`geovac/dirac_matrix_elements.py` "edge set and
  weights", `tests/test_dirac_matrix_elements.py` "108 tests") does not
  contain it -- the module has zero occurrences of "graph"; no test names
  D_equiv; no claim-matrix row. The only implemented graph-form Dirac,
  `geovac/dirac_lattice.py::DiracLattice`, is NOT isospectral to CH at
  nonzero hopping (t=1 spectrum -4.21 ... 5.93 vs diag +-1.5/+-2.5).
  LARGE by class (false backing), SMALL by consequence (no target number
  moves: every axiom verification uses the spectral form). Rescope to an
  honest Remark + a test that documents what is true. **Raised to PI.**
- KO-dimension label (:4049-4051, :4531-4534; `test_almost_commutative.py`
  docstrings): the verified sign pair (eps, eps') = (-, +) is the KO-3 pair
  in Connes' table, not KO-1; "3 + 6 = 9 == 1 hence ..." does not follow.
  Tests verify the signs at 1e-12, not the label. **LARGE, raised to PI**
  (Dabrowski-Dossena to be checked against the primary before citing).
- "propinquity rate" for the state-space GH rate at :1986, :2132,
  :2140-2141, :3198, :3569, :6874 -- C7 residue, invisible to the group1
  C16 regex (needs "Latr" adjacent). Remediate as a class.
- :3380-3382 inside `thm:gh_convergence` still says the O(log n/n) rate is
  "not rigorously proved ... deferred" -- stale against Paper 38
  (unconditional 2026-06-10).
- :3149-3150 "kappa = -1/16 (the Fock Jacobian Omega^-4)" -- the
  appositive asserts the derivation bridge C8 forbids.
- C3-boundary retags: :6477 (Frobenius-residual certificate), :3758-3762
  ("three rows verify ... to machine precision"), `thm:forced_count`
  :5113-5201 (n_max=2 enumeration + symbolic N_gen^2). Plus :5418 "128 per
  generation" vs the theorem's 128 N_gen^2 / 8 per generation.
- :1386-1388 "[MEASURED] it is an explicit Connes-style ... spectral
  triple, validated against the 21cm gap" -- the measurement backs the
  gap, not spectral-triple-hood.
- :878-882 order-zero axiom "verified in test_dirac_matrix_elements.py" --
  mis-pointed (the only order-zero test is the Lorentzian U_L one).
- :6518-6521 axiom-table (vi)/(vii) "<= 0.0675 / <= 0.101" vs
  `test_connes_axiom_audit_31.py:291-313, :328` asserting `< 0.5`.
- :4882 "provably non-overlapping" vs :5216 "seam theorem ... prove" for a
  seam the C3 pass tiered OBSERVATION; :4207-4213 untiered "forbidden ...
  vanish identically on ANY CC-compatible AC extension" (physics
  inference, not shown).
- Counts: :3715 "39/39" (49 collected: 47 + 2 skipped); :7028-7029 "67"
  (70: 66 + 4 skipped); :4053/:4167 "38 tests passing" (53 collected, 53
  pass -- found by the R2 read, 5.6). Locators: :2738-2739 "Connes-Marcolli 2008 Ch. 4"
  (the book has four chapters; the cosmic-Galois material is Ch. 1 SS1.7,
  Thm 1.100); :4032-4033 "Ch. 13" (SS1.13).
- Coverage gap: `thm:forced_count` matter sector "512 -> 256 -> 128;
  order-one / J-reality add nothing" (:5172-5176) is not computed by
  `test_trunk_qa_forced_count_moduli.py`.

**Paper 38.**
- **Circle Fejer constant** (:941-950 remark, :1702-1708 paragraph): the
  classical first absolute moment of the probability-normalised Fejer
  kernel is (2/pi) log N/N, not (4/pi). Closed form: M_n = pi/2 - (4/pi)
  sum_{k odd < N} 1/k^2 + (4/(pi N)) sum_{k odd < N} 1/k, N = n+1.
  Quadrature check: n M_n / log n = 0.9895 (n=50) -> 0.8321 (n=1600),
  doubling estimator -> 0.644 vs 2/pi = 0.6366 (4/pi = 1.273 is not in
  play). No standard normalisation gives 4/pi (Zygmund's K_n gives 2;
  integral 2 pi gives 4). **The paper's own backing test already says
  this**: `tests/test_trunk_qa_fejer_4_over_pi.py:235-268` names 2/pi
  "the circle-Fejer (unweighted) constant" and asserts the SU(2)
  constant is TWICE it. So "the constant 4/pi is the same on both sides"
  contradicts the test. Stein-Weiss 1971 SSI.1 is L^1 Fourier theory on
  R^n -- not a locus for this moment; drop that locator, keep Zygmund
  Vol. I Ch. III without "SS3.6". The SU(2) = 2 x circle relation
  (2 Vol(S^2)/Vol(SU(2)) = 4/pi mirrors 2 Vol(S^0)/Vol(S^1) = 2/pi) is an
  observation, to be stated as one. Propagates OUTSIDE trunk: Paper 40
  :966-973, :977-989 (`rem:stein_weiss_general`), :1937, :1943-1944 and
  the group1 synthesis :250, :597-599, :783, :1284 -- advisory here,
  **raised to PI** as a group1 carryforward; memory
  `l2_quantitative_rate_4_over_pi.md` to be re-read.
- :283-284 "Vol(S^2)/Vol(SU(2)) . (2/pi) = 4/pi": 4 pi / 2 pi^2 = 2/pi, so
  the factor must be 2, not 2/pi.
- :934-935 "monotonically decreasing for n >= 3 (verified at n in
  {2,...,1000})" -- no frozen test reaches 1000 (`test_central_fejer_su2.py:647-659`
  stops at 100).
- L3 (C_3 = 1, :980-1077) has no inline test citation though
  `tests/test_r25_l3_lipschitz_bound.py` exists and the matrix names it;
  Paper 40 SS3.3 may prove it at all ranks -> possible [INTERNAL THEOREM]
  upgrade, to be verified in P40 first. :553-554 whole-file cite of
  `test_p45_kplus_degeneracy.py` where one of five functions is relevant.
- Bibliography: duplicate Paper 40 bibitems (:1881 `paper40_unified` /
  :1935 `loutey_paper40`); no Paper 43 bibitem while :295 cites it inline;
  `latremoliere2018` :1840 JFA -> Trans. AMS 368 (2016) 365-411 (the body
  :591 has it right); orphans `chamseddine_connes2010`, `perez_sanchez2024`;
  preprint forms `gaudillot_vs2023` (IMRN 2025 rnaf197), `hekkelman2022`
  (LMP 112, 20); n-index clash :433 (n+1)(n+2) vs :625 n(n+1).

**Papers 0, 1, 7 and their tests.**
- `tests/test_ov_scaling_rigorous.py:169 < 1.8, :216 < 1.5` guard
  "sub-quadratic"; Paper 1 :449 says "[MEASURED] the O(V) ... scaling is
  verified by" it. Tighten the deterministic nnz exponent (:216) and state
  the measured exponent (~1.05) at P1:449.
- `test_volume_element_jacobian` restatement (KG2 note above).
- P1:26 "dynamical symmetry group SO(4,2) [barut1967, fock1935]" -- Fock
  1935 is SO(4); P1:271 Condon-Shortley phrasing; P0:749-750 Bohr vs de
  Broglie; P1 :15/:54/:59/:363 "reproduce the exact Rydberg spectrum"
  (-> quantum numbers, per P7:689); P1:393-398 K sentence inherits
  [CONJECTURE] -> own [OBSERVATION]. NIT cluster in CLM-A (P0 duplicate
  tags, "Paper 6" pointers to the archive, "~6%" vs 5.3%, "R^2 = 1.0"
  not computed).
- Orphan `loutey_paper18` in P7; inline Gaudillot-Estrada-van Suijlekom
  (arXiv:2310.14733, ID verified correct) without a bibitem in P7 :748
  and SYN :354 (add `gaudillot_vs2023` to both; re-run C20).

**group3 synthesis.**
- :164-165 "Paper 56 closes the cosmic-Galois comparison U*_GV
  hookrightarrow ..." -- Paper 56 :151-152 / :1298-1302 says the injection
  DIRECTION of the comparison, and that Phi^inj is not injective (factors
  through the abelianisation); SYN :1096-1099 itself says "a homomorphism
  (not a closed immersion)". Fix the arrow and the verb.
- :988 "Paper 25 (synthesis group)" -- it lives in group5 (:692 has it
  right). :589-594 provenance note omits the pre-rename test filename
  (`git`: renamed in 6e6ce40) and, per the R1 read (5.6), overclaims its
  scope -- the file backs unitarity + Rabi only, not applications
  (a)/(c)/(d) (matrix row 133 open); :601-603 "$10^4$ time steps" vs the
  live test's 1000 (10^4 is Paper 6's archived run). :1410 aquilanti_caligiana2003 CPL 366, 157
  -- issue is 2002 (verify vs publisher before editing). Upgrades: :176 and
  :1306-1309 omit Paper 38's unconditional state-space GH theorem (add with
  the correct metric name). Inline attributions without bibitems (:546,
  :867, :1073 Fathizadeh-Marcolli load-bearing, :1076, :1122/:1154, :1149)
  -- C20-baselined class; orphan bibitems (7).

### 5.4 What the run did NOT find (recorded so it is not mistaken for coverage)

- No hard-prohibition touch on main (C5 + CLM-B1/B2 enumeration of every
  K-sentence; KG5 clean twice).
- No live retracted-claim zombie on main under the widened C16.
- No fabricated external citation on main: the two WRONG-IDs were both
  seeds; the genuine citation defects are locators, years and hygiene.
- The kappa Observation loci (P0:652-662, P7:732-744) and the K-status
  remark (P32:2097) were each passed clean by a calibrated Opus agent.

### 5.5 Registry work done in this run (two-way proofs recorded)

**C16, two new entries** (`debug/qa/check_retracted_terms.py`), each with
the discrimination proof the 2026-08-22 rule requires, run by
`scratchpad/c16_fire_test.py` pointing `ROOT` at the worktree and then at
main:

| entry | FIRE (seeded worktree) | SILENT (main) |
|:--|:--|:--|
| `latremoliere-propinquity-named-for-gh-rate` -- anchors on `Latr...propinquity convergence` and `propinquity convergence at rate`; exempts only on an explicit denial (`no published`, `not achieved`, `strictly stronger`, ...). Files: group3, group3 synthesis, P32, P38, group1 synthesis. | SYN:348 LIVE (1/1 expected) | 0 live; P32:6112 "has no / published Latremoliere-propinquity convergence theorem" correctly exempt |
| `p45-kplus-compression-theorem-live` -- `K^+-compression theorem` (bare) or `Lorentzian quantum-metric convergence`; exempts on `degenerac / descope / withdrawn / retract / annihilat / intended to assert / open question / convergence question / open named / not a Lorentzian`. Files: P32, P38, P42-53, group1 + group3 syntheses, field guide. | SYN:693, P32 wt:6748, wt:6749 LIVE (3/3 expected) | 0 live; P46:1179, P48:452, P48:1953, P49:343 correctly exempt (all four are the descoped statement) |

Why the existing entries missed the seeds, for the record: the group3
entry exempts on `state-space` within +-5 lines, and the seeded sentence
carried "state-space" in its own next clause (the scalar half, correctly
labelled) -- so a wider pattern on that entry would have been silenced by
its own exemption; hence a separate entry with a denial-only exemption.
The first draft of the K+ entry fired on main at P49:343 ("convergence
question is an / open named research target"); the exemption was widened
and the proof re-run -- the live hit is what a registry test is for.

After the edit: `--gate trunk` PASS (11/30 entries, 38 locus patterns, 44
exempt), `--gate group1` PASS (7/30), `--gate group3` PASS (10/30);
`tests/test_retracted_terms_check.py` 5 passed; C22 check C now runs 30
retracted patterns (was 28), PASS.

**Still LLM-tier (no registry entry can see them).** CLB1b (a tier
overstatement: "established unconditionally by the Paper 35 verification"
-- semantic); and any zombie phrase that straddles a LaTeX line wrap
("Lorentzian / propinquity construction" on wt:6747-6748 was caught only
via the `K^+-compression theorem` alternative on the same line), because
`scan_entry` matches per line. Both stay with the claims reviewers.

**C16 widening still owed at remediation** (not done in-run because it
is a remediation of main, not a registry-coverage fix): a `propinquity
rate` alternative for the group1 entry once the 7 P32 loci are renamed
"state-space GH rate" -- two-way proof required then.

### 5.6 Completeness-critic pass

One fresh Opus critic (287k tokens, 45 tool uses) read the eleven
reviewer reports against the six documents on MAIN and answered "what
did nobody look at?" -- absence is not compliance. Its dispatch named
`Paper_1_Spectral_Graph_Methods.tex`, which does not exist; every
reviewer had used `paper_1_spectrum.tex`, so the slip is harmless and is
noted only so the record matches the prompt.

**Coverage gaps (regions no reviewer quoted), by document.**
- P0: Table 1/2 cells; part of the conclusions; 6/8 internal bibitems.
- P1: footnote :76; Appendix A :409-434; B.1, B.3; Table I caption.
- P7: the chordal-ansatz negative-result paragraph + table :471-493; the
  N-electron angular table :574-604; Appendix B :879-891; 5/7 internal
  bibitems.
- P32: single-reviewer coverage on :1670-3669; unquoted table cells at
  :1098-1115, :1152-1200 (whole table), :1354-1367, :1753-1776,
  :3738-3758, :4078-4092; paragraphs :176/:192/:207 and :624-713; a
  `propinquity` at :2159; the H1 block :4056-4171 + :4216-4240 (R2 below);
  ~17 unquoted `\paragraph` status strings; all 30 internal bibitems.
- P38: the Section 1 prior-art propinquity paragraph :176-203 (C7 prose,
  unquoted); footnotes :221/:603; Outline :304-319; Appendix A Steps 1-4
  :1622-1690 (the 4/pi derivation chain -- the CLC/Stein-Weiss findings
  reached it only through the body); the L3 backing test unreviewed; 7/9
  internal bibitems.
- Synthesis: abstract headlines (b)/(c) :54-60; Section I tags
  :175/:180/:182; "Status of packing-vs-physics" :271-287; Section IV tags
  :489/:569; Section V `[MEASURED]` :602/:606 (R1 below); Section VII
  `[CONDITIONAL]` :865; 24/26 internal bibitems.

**Backing tests no reviewer opened.** `tests/test_rabi_oscillation.py`
(SYN :591; the synthesis had no CODE dimension); `tests/test_paper2_corrections.py`
(matrix row 47 -- the only C5 artifact-side test for P32, audited by
nobody); `tests/test_r25_l3_lipschitz_bound.py` (P38 L3, uncited in the
paper); `debug/p38_g1g2_band_diagnostics.py`, `debug/p38_g1g2_scalar_prototype.py`;
`geovac/standard_model_triple.py` "(45 tests)" at P32:5183 (matrix row 48
already says mis-pointed); the 19 tests on synthesis matrix rows.
Partially read: `test_qed_self_energy.py`, `test_modular_hamiltonian.py:261-1074`.
Open NO-TEST rows: 133, 49, 88.

**UNMEASURED this run (per criterion x document).**

| criterion | unmeasured for |
|:--|:--|
| C7 (metric naming) | P0, P1, P7: 0 occurrences examined (the term does not occur -- vacuously fine, but no reviewer said so); synthesis sample size 1 (the seed) |
| C5 (K-label) | prose side: P0, P38 0 occurrences; P1 sample 1; P7 sample 2. Artifact side (`tests/test_paper2_corrections.py`): corpus-wide, nobody |
| C8 kappa | P1, P38, P32:3800-7527: zero occurrences examined |
| C8 4/pi | P0, P1, P7 |
| structural | footnotes, table cells, appendices, `\paragraph` strings as listed above; 83 internal bibitems had no mandated reviewer |

**Cross-reviewer contradictions (10) and how each was adjudicated.**
#1 gearing commutator (CODE-B vs CLM-A) = seed B2. #2 Berry control = seed
B1. #3 the 4/pi value (CODE-D vs CLM-C) = seed CLC. #4 CODE-C "PASS" vs
CLM-B1/B2 FAIL on the same paper -> recorded FAIL-by-content (5.2). #5
per-band injectivity = seed D2. #6 L3 "upgrade to theorem" (CLM-C) vs
"uncited backing" (CIT-2) -> both right; remediation adds the inline cite
and the tier. #7 forced_count tier: CLM-B2 is correct (no tag on the
paragraph); CODE-C's "INTERNAL-THEOREM" was the test's grade, not the
paper's. #8 SYN:348 = seed SYN1. #9 P32:1386 -- claims side examined,
backing side not; carried to remediation (the :1386-1388 item in 5.3).
#10 P7:79 "complete, machine-verified algebraic audit" vs matrix row 42
"~4 of 18 weak" -> soften P7:79 (already in 5.3 via KG2).

**The critic's one "demonstrated hit" was the seed.** It reported
P1:527 bibitem title "Deriving the Schrodinger Equation from Graph
Topology" as a live internal-title defect. That locus is seed DET1
(`trunk_seed_key_2026_09_01.json`, `deterministic_seeds`); on main
:527-529 reads "The Dimensionless Vacuum: Recovering the Schrodinger
Equation from Scale-Invariant Graph Topology," GeoVac Paper 7 (2026).
The critic had read the worktree. Recorded as a correct catch of a
planted defect by an agent that was not told about the seeds -- not as a
coverage gap with a live target.

**Recommended re-dispatches R1-R3, handled in the main session** (each was
a bounded read, not a fresh-panel job):

*R1 -- synthesis Section V `[MEASURED]` claims vs `tests/test_rabi_oscillation.py`.*
Ran the file on main: 3 passed in 16.9 s, rc=0 (`gates/rabi_main.out`).
Read all 593 lines. Findings, all genuine, none seed-adjacent:
- SMALL/MATERIAL. SYN :601-603 "Unitarity is preserved to machine
  precision over $10^4$ time steps" and :587-594 "re-run and reproduce ...
  regression backing is live in tests/test_rabi_oscillation.py". The live
  test (`run_norm_conservation`, :195) integrates **1000** steps at
  `max_n=10`, `dt=0.1`, threshold `norm_max_dev < 1e-10` (:217). The 10^4
  figure is Paper 6's original run (archive :31, :435, `<1e-14`), and
  matrix row 132 repeats 10^4. Counterfactual: the sentence's stated
  evidence base changes (10^3 in the live regression vs 10^4 archived);
  no headline changes. Fix: state both, or raise `n_steps` to 10^4 and
  re-measure.
- SMALL/MATERIAL. The provenance note :589-594 says the Section's results
  "were re-run and reproduce" with backing in that file, but the file
  backs only unitarity and application (b) Rabi (`p_tgt_peak > 0.95`,
  `period_error < 0.5`, :381-397). Applications (a) 20 H2 transitions
  0.16 % / 33 s, (c) MD 0.0003 %, (d) Langevin 300 K have NO test
  (matrix row 133, "COVERAGE GAP (still open)"). Counterfactual: three
  `[MEASURED]` results are presented as regression-backed and are not.
  Fix: scope the note to what the file backs; leave (a)/(c)/(d) tagged
  as archived Paper 6 measurements with the row-133 gap named.
- NIT. `run_off_resonance` asserts `p_tgt_max < 0.5` (:478) against a
  measured 0.0037 (135x slack); the Rabi assertions bound (0.95 / 0.5 %)
  rather than pin the quoted 99.98 % / 0.41 %; "machine precision" in the
  prose vs a 1e-10 threshold in the test. Tighten at remediation.

*R2 -- P32 H1 block :4013-4262 on main (Higgs falsifier; scope and
limitations; modular-propinquity reformulation).* Direction of the
scope-and-limitations paragraph (:4152-4163) is correct: no Higgs vev, no
M_3(C), chiralities not identified, "no claim is made that GeoVac
contains the Higgs" -- consistent with the positive-thin verdict. The
modular-propinquity paragraph (:4215-4238) uses Latremoliere's *dual
modular propinquity* as a tool on the Higgs falsifier and closes
negatively ("NOT Morita equivalence"); it is a scoped negative result,
not a C7 metric-naming zombie. Findings:
- SMALL/MATERIAL (count drift, same class as :3715 and :7028). :4053 and
  :4167 say `tests/test_almost_commutative.py` "38 tests passing"; on
  main the file collects **53** and 53 pass (3.4 s, rc=0,
  `gates/ac_collect.out`, `gates/ac_run.out`; the file grew at 23418fc,
  Sprint G3). Fix both loci to 53.
- Already in the ledger from CLM-B2: the Connes-Marcolli "Ch. 13" locator
  (:4032-4033, H4) and the KO-dimension arithmetic (:4049-4051, M4).
- NIT. The block cites CLAUDE.md by name (:4152; 15 such mentions in
  P32) and the sprint codename "Track 1 R2.5 L4" (:4156) -- audience
  register. Two `debug/` memo citations (:4235-4238) and
  `debug/data/h1_falsifier.json` (:4090) -- C14 advisory under the
  2026-06-17 policy. The reformulation paragraph and the thermal
  extension (:4171-4213) carry no inline tier (the :4207-4213 item is
  already in 5.3).

*R3 -- the 83 internal bibitems.* Title drift is deterministically
covered: C11 (`check_internal_titles.py --gate trunk`) checks every
internal bibitem's title and year against the cited paper's own
`\title`/`\date`, PASSES on main and was proven to FIRE on DET1 in this
run -- exactly the class the critic's one hit belongs to. The residual
class (a descriptor that mis-states a cited paper's *status*) was swept
by hand on main: across the six documents the only references to the
descoped Papers 45-49 are P38 :1551 ("proves this degeneracy theorem",
tagged `[OPEN]` for the Lorentzian extension) and its bibitem :1869-1873
("a degeneracy theorem", matching Paper 45's live title). Both state the
surviving negative result, not the withdrawn positive one. SYN :693 on
main carries no Paper 45 citation (the worktree's was seed SYN2). No
finding; no dispatch.

**Net effect on the verdict:** none on the FAIL (already FAIL on every
LLM dimension); R1 adds two SMALL/MATERIALs to the synthesis remediation
list and R2 adds one to P32's. The UNMEASURED table above is carried
forward verbatim so the delta run can be scoped to it.

### 5.7 Honest ceiling and what the delta run must carry

- **Two INCONCLUSIVE classes.** The delta run needs a fresh S3 seed on the
  P0/P7 backing tests and a fresh S2 seed on the P32 backing tests, blind
  (no class hint), and either Opus code reviewers or 3 seeds per Sonnet
  agent. Until then the trunk's code dimension is certified only on the
  classes its reviewers demonstrably caught.
- **Guided catches.** Six of the eight CODE seed catches were guided; the
  next run's CODE prompts must paste the frozen criteria only.
- **PI adjudications (nothing mechanical is blocked on them):** (1) the
  KO-dimension label -- fix the label to the verified KO-3 pair, or change
  the construction; (2) `prop:D_equiv` -- rescope to a Remark (recommended)
  or build the graph-form operator the definition promises; (3) the circle
  Fejer constant propagates into Paper 40 and the group1 synthesis --
  out of trunk scope, logged as group1 carryforward.
- **Remediated text is not clean text.** Run #1's delta found 4 of 11
  genuine findings were introduced by its own remediation. The next FULL
  run is gated on a clean DELTA, per the carryforward's standing rule.

## 6. Part F remediation (2026-09-02, v5.3.1)

The 30 Part F rows not covered by the v5.3.0 PI items were closed in one
sprint; each row in `docs/qa/trunk.carryforward.md` now carries its DONE
note, and the findings the remediation itself surfaced are recorded there
as F7.10-F7.16 (Paper 40's main-theorem class silently included tori; the
dim_H = 40 = g_3 = 1/Delta coincidence presented as structural at five
loci across Papers 32 and 42; c^2(3,2) vs the test's c^2(4,3); catalogue
entries vs instances in the wall count; two order-of-magnitude
restatements from the memo's own numbers; a bibitem-key collision; the
widened C16 entry finding eight more loci outside Paper 32).

Method notes for the DELTA that follows:

- **Everything measured, nothing remembered.** Every count, residual and
  exponent written into a paper this sprint was re-measured (test counts
  49/70/54/45; axiom residuals 0.0675/0.1013; nnz exponent 1.07; Rabi
  0.999756 / 0.4106 % / 0.003725; Fejer certificate 2.3 s at n = 1000).
- **Every registry change two-way proven.** The C16 widening fired on
  5/5 retired phrasings and stayed silent on 7/7 corrected ones -- and then
  fired on eight live loci the FULL run had not scoped, which is the
  argument for widening rather than hand-sweeping.
- **Gates run whole-target, scope stated:** C10 on Papers 0/1/7/18/32/38/
  40/42 + both syntheses; C11/C13/C14/C16/C18/C19/C21 on `trunk`;
  C11/C13/C14/C16/C19/C21 on `group1`; C5 corpus-wide; C15, C20, C22.
  All PASS.
- **Not done, carried:** F5 coverage debt (unchanged except App. A Step 3);
  three F1.15/F3.6 NITs the run record did not locate precisely ("H10
  sub-locators", the :1119/:1773 items, the `test_s3_eigenvalue_ground_state`
  docstring) -- named in the carryforward rows rather than guessed at.

## 7. DELTA run #1 (2026-09-02) = DEFECTS, remediated

Calibrated 9/9 (0/8 false positives) across seven blind reviewers on a
seeded worktree; verdict DEFECTS on 26 genuine items in the Part F
remediation itself (table G1 in `docs/qa/trunk.carryforward.md` Part G),
all fixed the same day. Both S3 and S2 code classes that the FULL run #2
Sonnet chunks missed were caught by Opus on the trunk-root tests. The
load-bearing non-seed catches: Paper 7's appendix still listing the
retired Jacobian identity; a Dąbrowski–Dossena equation locator wrong
(re-verified against the primary, Table 5 claim confirmed); an ε″ comment
contradicting the code; two Paper 42 twins of the dim_H = g_3 fix; five
[MEASURED] tags on bit-exact panels (register says PANEL-VERIFIED); a
19 GB test matrix in the default run. DELTA #2 with fresh seeds is the
next step; the FULL run stays gated on it.

## 8. DELTA run #2 (2026-09-02, PI-invoked) = DEFECTS, remediated

Scope = the DELTA #1 remediation (baseline reconstructed from the saved
diffs and committed on the seed branch). Calibration 8/9 with one void seed
(the CLAIMS-1 tier plant was defensible), 0/8 false positives. Fifteen
genuine items (carryforward Part H, table H1), two of them regressions of
DELTA #1's own fixes (the group1 master-theorem wording exceeding Paper 40;
the oddness pin that could not fail and whose inference the production
grading refutes). Two PI adjudications recorded (Paper 40's title; the KO
label of the finite triple, measured (-, +, +) with the finite grading).
Lesson: a delta that reviews only the remediation still finds defects in
it; the next FULL run stays gated on a clean DELTA #3.

### 8.1 PI adjudications applied (2026-09-02, v5.4.0)

The DELTA #2 H2 items were decided by the PI the same day and applied before
the FULL run: (1) Paper 40 retitled to the semisimple class, 23 loci
cascaded; (2) the finite combined triple carries **no** KO-dimension label —
the measured sign triple $(-,+,+)$ is printed (Paper 32 ×4, two modules, two
test files, matrix, C16 notes; `test_J_combined_KO3_sign` renamed
`test_J_combined_sign_pair`); (3) DELTA #3 skipped, FULL run unseeded by the
new default. Two further items the PI delegated: the Paper 7 NO-TEST row is
closed by `tests/test_paper7_graph_convergence.py` (measured 2026-09-02:
E0 = -0.41363 / -0.46375 / -0.47571 / -0.48883 / -0.49364 / -0.49713 at
n_max = 5 / 8 / 10 / 15 / 20 / 30, i.e. 17.27% -> 0.57% against -1/2,
monotone and one-sided; the six lowest eigenvalues at n_max = 30 lie within
0.06%, so the graph spectrum's bottom is dense and the claim is pinned at the
spectrum's edge), and Paper 40's §L5 / "Propinquity convention" now define
$\Lambda_{\mathrm{prop}}$ as van Suijlekom's state-space GH distance with
Latrémolière's propinquity named as not claimed (18 loci; bibitem
`vs2021_jgp` added — the compile gate caught the missing key on the first
pass). Version bumped to v5.4.0 by PI direction (gate change: seeding opt-in).

## 9. FULL run #3 (2026-09-02/03, unseeded, at fc41ec3 / v5.4.0) = FAIL

First run under the opt-in seeding default. 12 reviewers (6 code Opus, 3
claims Opus, 3 citations — P32 Opus, others Sonnet — 1 synthesis Opus) read
the committed corpus read-only; one completeness critic followed. All 12
were killed once by a session rate limit and resumed from their transcripts.
Every MATERIAL finding was verified by the PM against primary text or
recomputed by the PM's own route (the c^2 prefactor by quadrature at six
(n,l); the s/p splitting at twelve cutoffs on both lattices; the l-block
structure and top-mode weights; the Forced-count moduli by three independent
linear-algebra routes; the propagation-number null model; the Fejer moment
against the unit-S^3 geodesic distance at four cutoffs). Full ledger:
docs/qa/trunk.carryforward.md Part I (I.0 five re-pricings; I.1 P0/P1/P7;
I.2 P32; I.3 P38 + spill; I.4 gate integrity; I.5 synthesis; I.6 clean
surfaces). Scorecard: code FAIL x6, claims FAIL x3, citations CLEAN/CLEAN/FAIL,
synthesis FAIL, deterministic 13/13 PASS; roll-up FAIL with 55 verified
MATERIAL findings. Verdict rests on the standing calibration record (run #2
19/21, DELTA #1 9/9, DELTA #2 8/9; 0/8 false positives each).

Lessons (new): (1) a remediation that adds a test can CREATE the next
defect — the Paper 7 convergence test I wrote on 2026-09-02 pins kappa *
lambda_max against the target kappa was matched to, and the synthesis line I
rewrote the same day declared "no frozen test" an hour before the test
existed; (2) three reviewers converged independently on the same
false-positive from three dimensions (claims, P0 code, P7 code) — dimension
redundancy is what caught it; (3) the convention class strikes again at the
keystone: the Fejer moment's distance variable was the rotation angle, twice
the geodesic distance, and every artifact agreed because all shared the
convention — only a test against an independently defined distance can see
it; (4) a degenerate SAMPLE (18/24 zero elements) passed a genuine rank
computation for three runs — "the rank routine is sound" and "the sample
spans the algebra" are different claims; (5) C16 cannot see a line-wrapped
phrase or a module outside an entry's files list.

### 9.1 Remediation (2026-09-03, v5.4.1)

Applied in three batches (R1 Papers 0/1/7 + synthesis + coupling tests; R2
Paper 32 + gh_convergence + forced count + SM representation; R3 Paper 38 +
spill + C16). Measured facts recorded in the papers: production-lattice s/p
splitting at twelve cutoffs (37 -> 0.39 %); CG-construction non-decay
(129/68/84/139 %); lambda_max 6.618 -> 7.954; l-block structure (30
components, l = 11 block attains lambda_max, l = 0 block -> 4); c^2 prefactor
1/4 by quadrature at seven (n,l); gamma_n(module) = 2 x unit-S^3 moment at
n = 1, 2, 3, 5; Forced count 32 by three routes (matter rank 16, Majorana
16; degenerate sample -> 260); finite-algebra order-zero/one exactly 0 with
the CCM representation, combined 1.49/1.26 at n_max = 2 = GV residual
(factorisation residual 0); prop = 2 for 16/16 random *-closed subspaces;
L5 panel inequality holds at n_max 2-4 with decreasing margin. Full list:
carryforward Part I.8.

Lesson (6), added during remediation: a bash heredoc halves backslashes, so a LaTeX-bearing replacement string written through one turns `\ref` into a carriage return + `ef`; a subsequent text-mode round trip launders the CR into a newline, leaving corruption that no control-character scan can see. Both instances this session were caught only by pdflatex; C19 now checks for CR-tails at a line start. Rule reaffirmed: LaTeX-bearing edits go through script FILES, never heredocs.

## 10. DELTA #3 (2026-09-03, unseeded, on v5.4.0..v5.4.1) = DEFECTS, remediated (v5.4.2)

Five reviewers on the diff. Citations CLEAN (9/9 confirmed against primary
sources). The other four returned DEFECTS, all of one shape: the remediation
had been applied locus-by-locus, so every descoped reading survived wherever
the sweep had not looked -- the graph->S^3 convergence still asserted as
established at twelve loci across three documents, the Lorentzian "literal
identification" at four loci in Paper 42, the Forced-count 128/260 in the
field guide, Paper 57's displayed equation and Paper 32's own theorem
statement. Fixed claim-wide.

Two genuinely new findings, both PM-recomputed:
 (1) the s/p "splitting" is not a spectral gap at all -- no edge changes l,
     so the l = 0 and l = 1 blocks are disconnected components and the
     maximum-overlap mode selection is a near-tie (top-two within 1.2% at
     n_max = 30, selected eigenvalues 3.0 vs 0.011). The twelve percentages
     reproduce exactly and are not robust; they are a node-amplitude proxy.
 (2) UPGRADE: each l-block is the grid graph P_{n_max-l} x P_{2l+1}, so the
     whole spectrum is closed form and the saturation deficit is
     (42.6 + o(1))/n_max^2 -- a proven rate where the papers said none was.

Lesson (7): a locus-by-locus remediation of a claim-level defect is not a
remediation. When a run re-prices a CLAIM, sweep the claim (every assertion
of it, in every document), not the loci the reviewer happened to quote --
and check the abstract, the summary and the sibling paper first, because
that is where the survivors were every time.
Lesson (8): reproducibility is not robustness. Twelve percentages that
reproduce to five digits can still be an artifact of an arbitrary
tie-break; ask what selects the object before trusting the number.

### 10.1 L5 panel correction (2026-09-03, v5.4.3)

The fifth DELTA reviewer's last line before three server errors -- "the L5
panel inequality is violated at n_max = 7" -- was checked directly instead of
retried, and is right in substance: the panel height/Lip rises toward 1
(0.667/0.808/0.878/0.913 at n_max = 2..5, fitting 1.08 - 0.83/n_max) while
gamma falls, so the margin crosses zero near n_max = 6. The v5.4.1
"inequality holds" assertion therefore verified a small-cutoff coincidence.
Reframed as a measurement in the test, the module and both papers; the
panel-side quantity is a named open check.

Lesson (9): when replacing a tautological guard, check the replacement's
ASYMPTOTICS. Both self-inflicted defects of this arc (the Paper 7 convergence
test, the L5 panel check) were guards chosen from the same material as the
claim and true only in the regime already measured.
Lesson (10): three heredoc corruptions in one session (TAB for \times,
\nmax -> newline, halved backslashes). C19 caught two, pdflatex one. Use the
Write tool with raw strings for LaTeX-bearing edits -- no exceptions.
