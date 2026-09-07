# Advisory — Paper 60: The Generalized-Sturmian Secular Equation Block-Encodes Without a Metric — 2026-09-06

## The contribution, in one sentence

**How a Coulomb–Sturmian basis is *posed* sets its block-encoding cost: Avery's
isoenergetic (potential-weighted) posing makes the *atomic* secular equation a
standard, metric-free eigenproblem whose LCU 1-norm grows sublinearly — an
apparently novel, clean quantum-encoding target with both feared cost-risks of the
genre structurally absent — while for molecules the metric returns and the honest
result is a mitigated-but-not-dissolved conditioning frontier.**

That spine *exists* and is coherent (the metric is the load-bearing variable
throughout). The problem is not a missing thesis; it is a **diluted** one — see the
Verdict.

## Verdict

**DEFENSIBLE AFTER NAMED FIXES.**

The mathematics is sound and independently backed (tracked `geovac/sturmian_secular.py`,
`sturmian_sigma_law.py`; `tests/test_paper60_*`; the σ-law derived to exponent 2 with a
measured remainder; helium and H₂/H₂⁺ validations that bind and dissociate correctly).
The honesty machinery is, frankly, exemplary — the abstract's eleven tier tags cleanly
separate **[ESTABLISHED, from Avery]** (the classical isoenergetic method) from
**[MEASURED]/[OBSERVATION]** (this paper's quantum-encoding reading and resource
measurements), and the acknowledgments say so in one sentence. The novelty claim is
hedged to "appears novel" and is confirmed OPEN with no substantive collision by a
current adversarial lit scan. A committee would not send this back for a false claim.

It would send it back for two framing/scope issues that gate a *clean* defense. **(1)
Currency mismatch on the headline.** The atomic advantage is quoted as ‖M‖₁ ∼ K^0.84 in
*configuration count for a fixed atom*, then compared to Koridon's N^1.34 in *orbital
count for a growing molecule* — different axes, and a fault-tolerant-QC referee will
catch that "0.84 < 1.34" is not like-for-like. **(2) The spine is out-weighted by its
own frontier.** The clean, novel, defensible atomic win occupies roughly a third of the
paper; the molecular frontier (which has no clean win) plus a bonding-physics tangent
(decompactification/coherence fronts, lines 545–591) that carries *zero* encoding-cost
content occupy the rest. Neither issue is fatal; both are addressable by re-framing and
re-balancing, not new physics.

## What is genuinely strong

- **The core result is real, clean, and confirmed open.** The isoenergetic inversion
  makes the atomic secular equation metric-free (Eq. 5), turns the two documented
  cost-risks of a Sturmian quantum algorithm (metric conditioning; outer energy search)
  *structurally* absent (§5), and the lit scan (`debug/lit_scan/sturmian_qubit_memo.md`)
  finds no prior work combining a Sturmian basis + a quantum eigenvalue routine + the
  isoenergetic inversion. This is a genuine, narrow, defensible contribution.

- **Tier discipline is a model of the corpus's own honesty rules.** [ESTABLISHED, from
  Avery] is applied exactly to Avery's method and nothing more; the novelty is
  [OBSERVATION] ("appears novel," not "is novel"); every measured exponent is flagged as
  a numerical fit with the load-bearing claim named as the *regime/ordering*, not the
  digits (lines 281–288, 158–161). This is precisely §1.5 / MEASURED-honesty done right.

- **The "metric is one-electron, not N-electron" result is a clean second contribution.**
  The compound-matrix identity (§manyelectron, lines 519–543) — molecular-Sturmian
  configurations inherit an identity metric at every electron number — is a tidy
  [SYMBOLIC + MEASURED] structural statement with a real quantum-resource consequence
  (the metric penalty does not compound with N). It is currently under-billed.

- **The negative results are stated, not buried.** The paper answers its own hardest
  question in the negative (the molecular N-electron 1-norm is *not* sublinear, ∼n_orb^2.2,
  lines 605–637), identifies the structural reason (loss of the diagonal T⁰), and shows
  the gerade lever *fails* on water (lines 414, 656). A committee rewards a paper that
  kills its own overclaim before the examiner does.

## The hostile-examiner question

**"Your headline is the sublinear atomic 1-norm, K^0.84. But K is the configuration
count for a *fixed* atom — you are refining helium, not growing a system — so this is a
convergence/completeness axis, whereas the electronic-structure 1-norm literature you
compare against (Koridon N^1.34) scales with *system size*. In what sense is 0.84 < 1.34
a resource advantage rather than a comparison of two different axes? And since you never
report λ versus *qubits* for the isoenergetic case, and never run the algorithm, what
have you shown beyond 'a fixed atomic matrix has entries that decay with quantum
number'?"**

**Honest classification: fixable-now / name-as-scope — not a real hole.** The underlying
result (metric-free posing avoids the L² blow-up; the isoenergetic 1-norm grows
sublinearly in K where the L² one inflates) is genuine and the contrast is the load-bearing
claim, which the paper correctly identifies. What is missing is *axis discipline*: state
explicitly that K^0.84 is measured on a fixed-atom configuration/completeness axis (not a
system-size axis), that this is a different quantity from the L² §2 obstruction's λ-vs-Q
(spin-orbital count), and that the Koridon comparison is therefore heuristic, not
like-for-like. Best repair is additive: report the isoenergetic λ against *qubits* (or
log₂K + ancillas) so the payoff is stated in the same currency as the §2 obstruction and
the literature — then "the opposite of the naive inflation" becomes an apples-to-apples
sentence instead of a cross-axis juxtaposition.

## Recommendations (ranked, must-fix first)

1. **[MUST-FIX] Fix the currency mismatch on the headline (abstract; §4 lines 258–288;
   disambiguation lines 266–274).** The Babbush-2019 disambiguation *does* land its
   intended job — it cleanly separates ‖M‖₁ (fixed-atom secular matrix) from Babbush's
   plane-wave *gate count*, neutralizing the one dangerous wording collision. But the same
   paragraph then compares K^0.84 to Koridon's N^1.34 and calls it "below the best
   documented basis-choice lever," which silently swaps axes (fixed-atom config count K vs
   growing-molecule orbital count N). Either (a) drop the direct exponent comparison and
   state only that the isoenergetic posing avoids the L² superlinear inflation, or (b) keep
   it but add one clause: K is a completeness/configuration axis for a fixed atom, N is a
   system-size axis for a growing molecule, so the comparison is a heuristic regime
   contrast, not a like-for-like scaling result. This is the sentence a QC referee will
   stop on.

2. **[MUST-FIX] Rebalance for the spine — lead with the atomic win; quarantine the
   bonding-physics tangent.** Two moves. (i) The abstract opens on the *negative* setup
   (the naive Q^3.3 inflation, "the Sturmian basis is the worst choice") and the reader
   must travel four sentences before the payoff. Reorder so the clean atomic metric-free
   result is the first thing stated, with the obstruction as its foil. (ii) The
   decompactification-front (R* ≈ 2.14 n/Z, lines 545–564) and coherence-collapse (M₂,
   Eq. 16, lines 566–591) material is bonding physics with *no* block-encoding-cost
   consequence — it belongs to Paper 58's arc and rides into Paper 60 only on the shared
   compound-matrix identity. It dilutes the encoding thesis. Recommend trimming it to a
   one-paragraph pointer ("the same compound-matrix object governs the decompactification
   front; see Paper 58") or splitting it out. This is the single highest-leverage cut.

3. **[STRENGTHEN] Promote "the metric is one-electron, not N-electron" to a named
   contribution.** It is currently a mid-section [SYMBOLIC + MEASURED] paragraph
   (lines 519–543); it is one of the two genuinely clean, defensible structural results in
   the paper and deserves a billing in the abstract's lead and the conclusion on par with
   the atomic sublinearity, not below the decompactification tangent.

4. **[STRENGTHEN] Soften the abstract's Gaussian comparison to match the honest body
   caveat.** The abstract's "10³–10⁴× smaller than a standard (ratio-dependent) Gaussian
   metric" compresses a genuinely honest body discussion (lines 474–486: ratio-dependence,
   metric-vs-metric, production Gaussian absorbs the overlap classically) into a single
   hedge word. Carry one more clause of that caveat into the abstract so the number is not
   read as cherry-picked.

5. **[CONSIDER] Cite the two-projections literature at the commutator identity (line 387).**
   ‖[P_A,P_B]‖ = max_k σ_k√(1−σ_k²) saturating at ½ is classical operator theory
   (Halmos 1969; Böttcher–Spitkovsky 2010), per
   `debug/lit_scan/projector_angle_bonding_memo.md`. The principal-angle and compound-matrix
   priors are now cited (lines 369, 528) — this is the one remaining bare standard-machinery
   object. Verify whether Paper 32 (the composition-wall owner) already carries the Halmos
   cite; if not, a one-clause attribution here closes it.

6. **[CONSIDER] Tighten the "molecular payoff" language.** The gerade lever is demonstrated
   on H₂⁺ (one electron, maximal symmetry) and shown to *fail* on the first real polyatomic
   (water). The paper is honest about the equivalent-center scope (line 656), but calling it
   "the molecular payoff" in the abstract slightly outruns "the equivalent-center payoff."
   A two-word qualifier aligns the abstract with the body's own scope.

## Positioning check

- **Overclaim to watch:** the K^0.84-below-N^1.34 comparison (Rec. 1) — an axis swap, not
  a false number, but it reads as a stronger claim than the data support.
- **Overclaim, minor:** "the molecular payoff" for an H₂⁺-only, water-fails lever (Rec. 6);
  the abstract's compressed 10³–10⁴× Gaussian ratio (Rec. 4).
- **Undersell (the quiet failure):** the clean atomic result is structurally buried behind
  its own motivation in the abstract (Rec. 2i), and the one-electron-not-N-electron metric
  result is under-billed (Rec. 3). Both are real, defensible, and framed too modestly to
  land on a skim.
- **Prior art — well handled now:** the Babbush-2019 wording collision is disambiguated and
  the disambiguation lands; principal angles (Amos–Hall, King, West–Ruedenberg) and compound
  matrices (Löwdin, Prosser–Hagstrom, King, Burton) are cited where the memos said they must
  be. The one residual is Halmos/Böttcher–Spitkovsky at the commutator (Rec. 5).
- **Tier/rhetoric (§1.5):** clean. No ontological-priority drift; "apparently novel" is
  correctly hedged; exact≠accurate is respected (the residual 7 mHa is named as basis
  incompleteness, not a missing metric, lines 249–256). No tier violations found.

## What I don't know

- **Whether λ-vs-qubits is cheaply measurable for the isoenergetic atomic case.** Rec. 1's
  preferred repair (report the payoff in qubit currency) assumes the isoenergetic secular
  matrix's block-encoding admits a clean λ(Q) measurement analogous to the §2 L² sweep. The
  PM should confirm the machinery (`sturmian_secular.py`) can produce it before committing to
  that framing; if it cannot, the additive-clause repair (1b) stands alone.
- **Whether Paper 32 already carries the Halmos/two-projections citation** for the
  composition-wall commutator (Rec. 5). If it does, Paper 60 may only need a cross-reference,
  not a fresh bibitem.
- **The exact intended weight of the decompactification material (Rec. 2ii).** It was
  captured into §manyelectron deliberately (CHANGELOG v5.10.2) as the compound-matrix bridge;
  I am reading it as a distraction from the *encoding* thesis, but if the PI wants Paper 60 to
  double as the encoding-side home for that bonding result, the call is theirs — I flag the
  dilution, not the decision.
- **Certification currency:** CHANGELOG shows a "Paper 60 CERTIFIED" pass, but the standing
  staleness note (§2) marks all `*.done.md` certs stale; nothing in this brief depends on that
  cert — the backing tests are what I relied on — but a `/qa delta` after any edits is the
  right closer, not this advisory.
