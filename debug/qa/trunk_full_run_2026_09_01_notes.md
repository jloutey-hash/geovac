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
