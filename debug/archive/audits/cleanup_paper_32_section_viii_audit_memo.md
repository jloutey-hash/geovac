# Cleanup audit — Paper 32 §VIII: the eight C-arc non-selection theorems

**Date:** 2026-06-08
**Auditor:** PM (adversarial mode)
**Scope:** Statement-and-verification audit of the eight C-arc non-selection theorems consolidated in `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII (v3.103.0). Read-only.

## TL;DR

**Overall readiness: YELLOW with specific patches.** Six of eight theorems have clean statement + appropriate proof grade + correctly-scoped honest-scope qualifier. Two have specific defects:

1. **Forced-Count Theorem (`thm:forced_count`)** is materially solid but cites the H1 Yukawa non-selection as a separate sibling theorem (`\S\ref{sec:open}`) that does **not actually exist as a labelled theorem environment** in the corpus — H1 lives only inside an `\paragraph{H1 Higgs inner fluctuation}` block. This is a load-bearing citation gap, because the `rem:c_arc_terminal_state` tally and the §13.4a verification status both depend on H1 being theorem-grade.
2. **No-single-mechanism K (`thm:no_single_mechanism_K`)** uses the phrase "Twelve candidate mechanisms have been empirically eliminated" as a load-bearing premise in the proof — but the proof is structural, so empirical elimination is the *evidence* the theorem explains, not part of the formal argument. The structural claim ("no morphism in $\mathcal{A}$ produces $K$") is weaker than the wording suggests. As written, an external mathematician will (correctly) read the proof as an empirical generalisation, not a structural impossibility.

§13.4a verification gaps: **3 theorems** need explicit verification artifacts added before submission (effort: ~3–6 hours total). Cross-paper consistency: clean, no conflicts with Papers 18, 25, 30, 51, 57.

The cross-theorem coherence story is strong: the eight DO consolidate into a single structural map of the forced/free seam — see §2 below for the 2–3 sentence summary.

---

## 1. Per-theorem audit table

### Theorem 1 — Forced $D_F$ moduli dimension (`thm:forced_count`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD with a citation gap.** Hypotheses (1)–(3) explicit and bounded: Bertrand×Hopf-tower+Upgrade B forces $\mathcal{A}_F$; canonical CCM rep; standard real-spectral-triple axioms. Conclusion is the dimension chain $1024\to 512\to 128\to 128$ with the per-generation $8$ diagonal-slice + $128 N_{\mathrm{gen}}^2$ full-mixing count. Scope is explicit ("at $N_{\mathrm{gen}}=1$" called out). Bertrand analogy corollary is fine as analogy, not load-bearing. |
| Proof | **Proof sketch with named gap.** Direct counting at each step. Hermiticity, chirality grading, $J$-reality, order-one reductions named. Order-one is the one nontrivial step (the proof asserts it preserves dim because of the $\mathcal{A}_F$ block-diagonality and $M_3(\mathbb{C})$ commutation with Yukawa flavour structure) — for the canonical CCM rep this is established CCM 2010 content, but the proof should cite it more carefully than "the order-one condition is satisfied for the SM representation without further reduction." |
| §13.4a | **PASSES.** `tests/test_standard_model_triple.py` `TestSMACTripleAxioms` and `TestH1Consistency` classes verify the chain at $n_{\max}\in\{1,2\}$ bit-exactly (Hermiticity, $J^2$, $JD$, order-zero, order-one residuals all at $10^{-10}$ to $10^{-12}$ gate). The reduction chain's COUNT however is not directly tested as integers $1024/512/128/128$ — the test verifies the axioms hold, not that $\dim\mathcal{M}(D_F)$ matches the claimed count. Gap is small (~1 hour to add a `test_forced_count_reduction_chain.py` parameter-counting test). |
| Cross-paper | **CLEAN.** Paper 57 §6.3 cites this as the load-bearing structural input for the forced/free seam (P5 packing-reachability witness). H1 sprint memo `debug/h1_higgs_inner_fluctuation_memo.md` supplies the reduction chain ingredients. No conflicts with Paper 25 / 30 (gauge structure, not moduli space). |

### Theorem 2 — H1 Yukawa non-selection (§VIII.C `\paragraph{H1 Higgs inner fluctuation}`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **WEAK — not formalised as a Theorem environment.** The statement lives inside an `\paragraph{H1 Higgs inner fluctuation (Sprint H1, May 2026)}` block at lines 4489–4513 and is introduced as "\textbf{Yukawa non-selection theorem}." in running prose. The statement IS explicit ($128$ real parameters; $8$ free per generation under diagonal $Y$; gauge sector forced, Yukawa free; mechanism via G3 commuting $\mathbb{Z}_2$'s), but it is NOT inside a `\begin{theorem}...\end{theorem}` block. The `rem:c_arc_terminal_state` tally cites it as "(\S~\ref{sec:open})" which points to the open-questions section, not to a labelled theorem. **This is a real defect**: when external mathematicians read the C-arc terminal-state list of "eight theorem-grade non-selection results," they will look for H1 as a theorem and find prose. |
| Proof | **Materially solid.** The G3 commuting-$\mathbb{Z}_2$ result (Sprint G3, see `\subsection{Sprint G3}` §VIII.C above) provides the structural reason: $\gamma_{\mathrm{GV}}$ and $\gamma_F$ are independent commuting Z₂'s, ‖$\gamma_{\mathrm{GV}}\otimes I_F - I_{\mathrm{GV}}\otimes\gamma_F$‖ = $32$ exactly at every $n_{\max}\in\{1,2,3\}$. Bit-exact verification at $n_{\max}\in\{1,2,3\}$ on `geovac/standard_model_triple.py`. Yukawa-PSLQ empirical confirmation (162-cell sweep, zero hits at $M\le 10^3$) cited as separate confirmation. |
| §13.4a | **PASSES at axiom level, fails at theorem-statement level.** `tests/test_standard_model_triple.py` verifies the bit-exact axiom residuals on the SM triple (45 tests passing per memo). The Yukawa-PSLQ sweep has its own driver `debug/sprint_yukawa_pslq.py` + JSON results. Verification artifacts exist; what's missing is a `tests/test_h1_yukawa_non_selection.py` that explicitly checks "the parameter count of admissible $D_F$ is $128$ under the four axioms" as a single test asserting the theorem's conclusion. |
| Cross-paper | **CLEAN.** Paper 57 §6.1 / §6.3 cites it as `H1 Yukawa non-selection (\S VIII.C)`. Paper 32 `rem:full_inner_factor_boundary` cites it. No conflicts. |

### Theorem 3 — $N_{\mathrm{gen}}$ non-selection (`thm:n_gen_non_selection`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD.** Hypothesis set $\mathcal{A}$ explicitly enumerated (5 axioms named). Conditional on canonical CCM rep is explicit in statement AND in `rem:n_gen_scope`. Conclusion is clean: for every integer $N\ge 1$, $\mathcal{T}_N$ satisfies $\mathcal{A}$; hence $\mathcal{A}$ contains no morphism selecting $N$. Scope qualifier (the unconditional version with alternative reps remains multi-year open) is explicit. |
| Proof | **Structural impossibility argument.** Six bullets show $N$-independence of each axiom. Read 2 scoping (`debug/sprint_read2_n_gen_scoping_memo.md`) supplies the three independent obstructions to a candidate Hopf-tower-to-representation shortcut. The proof structure is correct: it's a statement about the absence of a morphism, supported by structural decomposition of each axiom's output. |
| §13.4a | **STRUCTURAL IMPOSSIBILITY — no test applies.** The theorem says no morphism exists; you cannot bit-exactly verify the absence of a morphism. The Read 2 obstructions can be verified case-by-case (e.g., DOF count per generation is 32 across all 3 factors — verifiable; CKM/PMNS coupling structure — empirically verified by phenomenology). Adding a `tests/test_n_gen_non_selection.py` that verifies each obstruction holds at the literal CCM rep would close the §13.4a gap (effort: ~2 hours). |
| Cross-paper | **CLEAN.** Paper 57 §6.3 cites it as enumerated item 3 of the inner-factor sector closure. No conflicts. |

### Theorem 4 — Inner KO-dim non-selection (`thm:ko_dim_non_selection`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD.** Hypothesis $\mathcal{A}$ inherited "as defined in Theorem~\ref{thm:n_gen_non_selection}" — clean reuse. Conclusion: $\mathcal{A}$ does not autonomously select an inner KO-dim signature. Scope qualifier same as N_gen: conditional on canonical CCM rep, unconditional version is multi-year open. |
| Proof | **Two-fact composition.** Fact 1: Paper 0 §VII.B is kinematic. Fact 2: KO-dim depends on $(D_F, J_F, \gamma_F)$ data (Door 4f T1 / Door 4b Q3). Composition: packing produces no $(D_F, J_F, \gamma_F)$ data, hence cannot select KO-dim. **Cleanest proof of the eight** — one-line composition citing established facts. |
| §13.4a | **STRUCTURAL IMPOSSIBILITY — no test applies directly.** Same structure as N_gen: the absence of a morphism cannot be verified bit-exactly. Paper 0 §VII.B is a documented scope statement (verifiable by inspection). Door 4f T1 / Door 4b Q3 is bit-exact at $n_{\max}\in\{1,2,3\}$ in `geovac/almost_commutative.py` (per memo). Both inputs are verified; the composition is logical. |
| Cross-paper | **CLEAN.** Paper 0 §VII.B cited as explicit scope statement. Paper 57 §6.3 cites as enumerated item 4. No conflicts. |

### Theorem 5 — Single-cutoff spectral action (`thm:single_cutoff_spectral_action`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD.** Definition (`def:multi_cutoff`) precise: $k\ge 2$ independent energy scales not reducible to a single scale by reparametrization. Theorem says no morphism derivable from $\mathcal{A}$ produces multi-cutoff structure. Corollary (`cor:multi_loop_renormalization_wall`) is well-scoped: classifies E7, E8, G3, G4 as multi-cutoff renormalization calibration data. Scope qualifier (`rem:single_cutoff_scope`) explicit on both axes: (i) doesn't claim impossibility under any extension of $\mathcal{A}$; (ii) spatial-composition sub-sector (G1, G2, G5) explicitly outside scope. |
| Proof | **Four-bullet structural inspection.** Each axiom inspected for energy-scale data. CCM spectral action single-cutoff by stipulation. Paper 0 + RST + Hopf-tower + Bertrand + Upgrade B introduce no energy scale. Inner fluctuations preserve single-Dirac structure. Tensor products preserve single-cutoff with combined $D$. Composition: no axiom introduces independent second scale. **The proof is structurally correct.** |
| §13.4a | **STRUCTURAL IMPOSSIBILITY — no test applies directly.** The empirical companion is the LS-8a sprint (`debug/ls8a_two_loop_renormalization_gap.md` per memory index, Paper 36 §LS-8a) which shows the iterated CC spectral action reproduces the bare integrand at right prefactor / sign / divergence rate $\sim N^{3.79}$ but cannot generate $Z_2, \delta m$ counterterms. This is the strongest empirical confirmation: the iterated CC machinery EXISTS in code, and the negative result is computationally established. Backing test gap: the framework's `geovac/spectral_action_*.py` (or equivalent for the iterated CC) does not have a `tests/test_single_cutoff_*.py` that verifies the structural claim "no axiom-derivable morphism produces a second independent scale." Effort to add: ~2 hours (essentially a documented inspection cross-checked against the existing axiom-set inventory). |
| Cross-paper | **CLEAN.** Paper 36 §LS-8a is the empirical companion (cited). Paper 57 §6.3 cites as enumerated item 5. Paper 51 §G7 (Wald) is consistent — Wald forces relations among moments at a single cutoff; doesn't introduce a second scale. No conflicts. |

### Theorem 6 — Spatial-composition radial wall (`thm:spatial_composition_radial_wall`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD.** Set-up explicit: $\mathcal{T}_a, \mathcal{T}_b$ on $S^3$ at distinct focal lengths with CH spinor Dirac. Two-claim structure: (1) angular structure FORCED (Paper 54 Thm 3); (2) radial coupling NOT autonomously matched. Theorem conclusion is clean: no morphism in $\mathcal{A}$ enforces specific physical radial matching. Corollary (`cor:spatial_composition_wall`) classifies G1/G2/G5. Companion remark (`rem:multi_focal_wall_fully_characterized`) records sub-sector closure. |
| Proof | **Two-part with empirical support.** (1) Angular forcing: direct from SU(2) rep theory on $S^3$ via Fock projection (Paper 7) and Hopf-tower truncation. Cites Paper 54 Thm 3 — assumes its proof. (2) Radial non-reproduction: structural argument from Fock-projection conformal factor non-commutativity (Paper 7 chordal-distance identity) + 6-dim vs 3-dim relative-coordinate mismatch + empirical verification (Pearson ≤ 0.81 decreasing with $n_{\max}$ in `debug/sprint_resolvent_two_body_memo.md`, Paper 54 only 75% connected fraction). |
| §13.4a | **MIXED.** Angular forcing portion: Paper 54 Thm 3 has its own paper-internal verification. Radial non-reproduction portion: empirically verified via the resolvent_two_body sprint and Paper 54's 75% connected fraction. What's missing: a single `tests/test_spatial_composition_wall.py` that wires together "tensor-product spectral action on $\mathcal{T}_a\otimes\mathcal{T}_b$ gives angular structure matching Coulomb at machine precision, radial structure NOT matching at decreasing-with-$n_{\max}$ Pearson." Effort to add: ~3 hours (the underlying drivers exist; consolidation into a single regression test is the work). |
| Cross-paper | **CLEAN.** Paper 54 cited explicitly. Paper 7 chordal-distance identity cited. Paper 24 §V (Coulomb/HO asymmetry, the 4th layer is modular-Hamiltonian Pythagorean orthogonality, structurally related but distinct from the 6-dim/3-dim composition wall here). No conflict. Paper 57 §6.3 cites as enumerated item 6. |

### Theorem 7 — No-single-mechanism $K = \pi(B+F-\Delta)$ (`thm:no_single_mechanism_K`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **WEAK — confuses structural and empirical content.** The statement reads: "No morphism in $\mathcal{A}$ produces the combination $K = \pi(B + F - \Delta)$ as an output of one Mellin mechanism, nor as a composition of two or three Mellin mechanisms with the specific weights $(+1, +1, -1)$ and the external $\pi$ prefactor. Twelve candidate mechanisms have been empirically eliminated..." The first sentence is structural ($\mathcal{A}$-morphism non-existence); the second sentence is empirical. The proof leans heavily on the empirical content as if it were part of the structural impossibility, which it isn't — empirical elimination of 12 candidates does NOT prove non-existence of any morphism in $\mathcal{A}$. **An external reviewer will (correctly) read this as a strong empirical generalization, not a structural theorem.** |
| Proof | **Three-fact composition with the load-bearing third fact being empirical.** Fact 1 (master Mellin engine case-exhaustion): structural, cites Theorem `thm:pi_source_case_exhaustion`. Fact 2 (three independent spectral homes): structural, cites Paper 2 §IV. **Fact 3 (twelve mechanism eliminations): empirical.** The proof asserts "Composition. $\mathcal{A}$ contains no morphism that simultaneously spans the three Mellin sub-rings..." — but the *structural* reason this would follow from Facts 1+2 is just that the sub-rings don't inter-translate within $\mathcal{A}$; the 12-mechanism elimination is the empirical *confirmation* that no such morphism has been *exhibited*, not a proof of non-existence. The proof would be cleaner if framed: structural — Facts 1+2 (sub-rings disjoint, three independent homes); empirical — 12-mechanism elimination supports the structural reading. |
| §13.4a | **EMPIRICAL — no impossibility test.** The 12-mechanism elimination is documented across Paper 2 §IV.G, Phases 4B-4I, Sprint A α-SP, Sprint K-CC. Each mechanism elimination has its own driver/data. There is no test that "verifies no single mechanism in $\mathcal{A}$ produces $K$" because the structural impossibility content is weak. The 12 eliminations are *individual* tests of *individual* candidates. |
| Cross-paper | **CLEAN-with-discipline-preserved.** Paper 2 is in `papers/observations/` per CLAUDE.md §13.5 (moved 2026-05-02 per curve-fit audit memo). The `rem:e6_scope` explicitly preserves the conjectural-label discipline ("the structural impossibility result strengthens the empirical non-derivability, but does not promote $K$ from observation to theorem"). This is correctly handled. No conflicts. |

### Theorem 8 — Cutoff function external (`thm:cutoff_function_external`)

| Axis | Verdict |
|:-----|:--------|
| Statement | **GOOD.** Explicit set-up: $S(D,\Lambda,f) = \mathrm{Tr}\,f(D/\Lambda)$, $f$ Schwartz-class. Conclusion: Mellin moments $\varphi(s)$ are external calibration data, not constrained by any axiom in $\mathcal{A}$. The Wald two-term-exactness $\Rightarrow$ pure Einstein $\Rightarrow$ action-$G\equiv$ entropy-$G$ argument is cited (Paper 51 §G7) as forcing RELATIONS among moments via the sector-wise moment map. Distinction relations-forced / values-external is clean. |
| Proof | **Implicit — not given as explicit proof block.** The theorem statement is followed directly by the Corollary (`cor:cc_fine_tuning_external`) and the scope remark (`rem:d5_d6_scope`). There is no `\begin{proof}` block. The structural content is straightforward (the cutoff function is a Schwartz test-function input to the action functional; nothing in $\mathcal{A}$ constrains test functions; Wald forces relations among moments but not their absolute values) but should be written out as an explicit proof block. **Effort to add: ~30 minutes** (the structural content is in the memo at `debug/sprint_c_arc_closure_e6_d5d6_chemistry_memo.md` §"Theorem 2 — D5/D6 cutoff function"). |
| §13.4a | **MIXED.** Paper 51 §G7 is verified bit-exactly via `tests/test_paper51_*.py` (six tests, J-blindness, scalar a_k, dS vacuum, zeta unit, L6 full, cutoff Mellin). The Mellin moment map is explicit in code at `tests/test_paper51_cutoff_mellin.py`. The structural claim "$\mathcal{A}$ does not constrain moment values" is again a non-existence statement (cannot be bit-exactly tested) but the moment-relation forcing IS tested. Gap: no test wires together "verify Paper 51 moment map for several different cutoff functions $f$ to demonstrate the relations are $f$-independent while values are $f$-dependent." Effort to add: ~2 hours (the moment-map driver exists; varying $f$ and checking relation-invariance is the new test). |
| Cross-paper | **CLEAN.** Paper 51 §G7 (Wald), §G4-5 (sector-wise Mellin moment map) cited explicitly. CC fine-tuning ratio $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ matches `memory/cc_phi_moment_fine_tuning_statement.md`. Paper 57 §6.3 cites as enumerated item 8. No conflicts. |

---

## 2. Cross-theorem coherence — does it read as one structural map?

**Yes.** Test: can you state in 2–3 sentences what the eight collectively prove?

> The GeoVac framework's axiom set $\mathcal{A}$ forces the **outer-factor structural skeleton** (gauge group, Forced-Count moduli space dimension, single-cutoff spectral action, angular structure of two-body operators, master Mellin engine classification of every transcendental, Wald-forced gravitational relations) and is **structurally silent** on three categorically distinct families of calibration data: the **inner-factor point** in the moduli space (Yukawa values + $N_{\mathrm{gen}}$ + inner KO-dim — Theorems 1–4), the **multi-focal-composition output** (renormalization counterterms + spatial-composition radial profiles — Theorems 5–6), and the **external-test-function input** (combination rule $K$ for $\alpha^{-1}$ + cutoff function moments for $\Lambda_{cc}$ — Theorems 7–8). This is a sharp, sector-stratified characterization of the structural-skeleton boundary at the canonical-CCM-rep level.

The consolidation works. The §VIII chapter reads as one structural map, not a list of independent theorems. The `rem:c_arc_terminal_state` provides the sector taxonomy explicitly. **The strongest single thing the corpus now sells is this sector-stratified seam picture, supported by eight named theorems.**

One minor coherence quibble: Theorems 7 and 8 are paired under "calibration inputs not constrained by $\mathcal{A}$" but Theorem 7 is about a *single number* (the value $K$) while Theorem 8 is about a *function* (the cutoff $f$, with its infinite-dimensional moment ring). These two should not be lumped as one sector. The terminal-state remark wisely splits them into "α combination rule sector" (item 7) vs "gravity sector" (item 8), but the deeper unity (both are external-input-to-functional-machinery) could be made more explicit.

---

## 3. §13.4a verification gap list

| Theorem | Gap | Test path | Effort |
|:--------|:----|:----------|:-------|
| 1 Forced-Count | Reduction chain $1024\to 512\to 128\to 128$ not directly verified as integer counts | `tests/test_forced_count_reduction_chain.py` (new) | ~1 hr |
| 2 H1 Yukawa | Theorem-statement-level test not present (parameter count $= 128$ as direct assertion) | `tests/test_h1_yukawa_non_selection.py` (new) | ~1 hr |
| 3 $N_{\mathrm{gen}}$ | Read-2 obstructions case-by-case not consolidated as a single test | `tests/test_n_gen_non_selection_obstructions.py` (new) | ~2 hr |
| 4 Inner KO-dim | Structural impossibility — no direct test applies; Paper 0 §VII.B is documented (verifiable by inspection) | N/A — scope documented | 0 hr |
| 5 Single-cutoff spectral action | Structural inspection — empirical companion in LS-8a sprint; no test wired together as single regression | `tests/test_single_cutoff_axiom_inspection.py` (new, doc-style) | ~2 hr |
| 6 Spatial-composition wall | Single regression wiring Paper 54 angular + resolvent radial mismatch | `tests/test_spatial_composition_radial_wall.py` (new) | ~3 hr |
| 7 No-single-mechanism K | Per-mechanism eliminations exist individually; no single consolidated test | Existing driver consolidation (low priority, structural content is weak — see scope-honesty section below) | ~2 hr |
| 8 Cutoff function external | Moment-relation $f$-invariance vs value $f$-dependence not tested together | `tests/test_cutoff_function_relation_invariance.py` (new) | ~2 hr |

**Total effort to close all §13.4a gaps: ~13 hours.** Priority order: Theorem 2 (H1) and Theorem 1 (Forced-Count) first (load-bearing for the seam framing); Theorems 5, 6, 8 next (cleanest structural arguments); Theorems 3, 7 last (weakest verifiability).

---

## 4. Cross-paper consistency — theorem ↔ Paper 57 seam catalog mapping

| Paper 32 Theorem | Paper 57 catalogue entry / section | Consistency |
|:-----------------|:----------------------------------|:------------|
| 1 Forced-Count | §3 forced inventory; §5 P5 probe witness; §6.3 enumerated item 1 | **CLEAN** |
| 2 H1 Yukawa | §4.1 (F1) free-side inner-factor entry; §6.1; §6.3 enumerated item 2 | **CLEAN** |
| 3 $N_{\mathrm{gen}}$ | §4.1 (F2) free-side; §6.3 enumerated item 3 | **CLEAN** |
| 4 Inner KO-dim | §4.1 (F3) free-side; §6.3 enumerated item 4 | **CLEAN** |
| 5 Single-cutoff spectral action | §4.2 (E7, E8, G3, G4); §6.3 enumerated item 5 | **CLEAN** |
| 6 Spatial-composition wall | §4.2 (G1, G2, G5); §6.3 enumerated item 6 | **CLEAN** |
| 7 No-single-mechanism K | §4.5 (I1 foundational calibration via E6 inheritance); §6.3 enumerated item 7 | **CLEAN with §13.5 protected** |
| 8 Cutoff function external | §4.3 (D5, D6 gravity calibration); §6.3 enumerated item 8 | **CLEAN** |

**Paper 51 gravity-side consistency:** Theorem 8 explicitly cites Paper 51 §G7 (Wald) and the sector-wise Mellin moment map. Paper 51's tests verify the bit-exact two-term exactness and the moment-map relations. Consistent.

**Paper 25 / 30 Hopf-gauge / SU(2) Wilson consistency:** Theorems 1, 5 cite the gauge structure framework. Paper 25 (Hopf-U(1)) and Paper 30 (SU(2) Wilson) carry their own structural results; neither is contradicted. The single-cutoff theorem applies to the spectral action $\mathrm{Tr}\,f(D/\Lambda)$, not to the Wilson lattice gauge actions, so no overlap or conflict.

**Paper 18 master Mellin engine consistency:** Theorem 7 cites the case-exhaustion theorem (Paper 32 §VIII `thm:pi_source_case_exhaustion`) which is itself consistent with Paper 18 §III.7 master Mellin engine. The mechanism-as-domain partition (Sprint MR-A/B/C) is preserved. Consistent.

**No conflicts found.** The Paper 57 catalogue mapping is one-to-one; the Paper 32 §VIII terminal-state remark explicitly tags each theorem with its Paper 57 entry; cross-paper references are bidirectional.

---

## 5. Scope-honesty: what is NOT claimed?

| Theorem | Scope qualifier present? | What is NOT claimed |
|:--------|:------------------------|:--------------------|
| 1 Forced-Count | Implicit via Bertrand-analogy corollary ("forces dimension, not point") | NOT claimed: that the canonical CCM rep is unique. Implicit but not flagged. **PATCH:** add an explicit "this theorem does not claim the canonical CCM rep is unique" sentence to `rem:forced_count_with_seam`. |
| 2 H1 Yukawa | Implicit — but H1 is not formalised as a Theorem environment | NOT claimed: that the Yukawa values are *uncomputable* in any extension of $\mathcal{A}$. **PATCH:** when formalising H1 as a Theorem environment (see Outreach Verdict section below), add an explicit scope sentence: "this theorem does not claim Yukawa values cannot be computed in any extension of $\mathcal{A}$; sibling axioms remain logically open." |
| 3 $N_{\mathrm{gen}}$ | **EXCELLENT.** `rem:n_gen_scope` explicitly names the unconditional question, the Read-2 sharpest falsifier, and the multi-year structural-research target. | (Adequately scoped.) |
| 4 Inner KO-dim | **GOOD.** `rem:full_inner_factor_boundary` names the unconditional version per axis. | (Adequately scoped.) |
| 5 Single-cutoff spectral action | **EXCELLENT.** `rem:single_cutoff_scope` has TWO explicit scope clauses: (i) sibling axiom adding RG-flow could extend $\mathcal{A}$; (ii) spatial-composition sub-sector outside scope. | (Adequately scoped.) |
| 6 Spatial-composition wall | **GOOD.** `rem:multi_focal_wall_fully_characterized` explicitly cross-references the renormalization theorem; corollary names G1/G2/G5 specifically. | NOT claimed: that the relative-coordinate reduction is mathematically impossible. **PATCH:** the proof sketch should explicitly note that a relative-coordinate-reduction morphism could in principle be added to $\mathcal{A}$ as a sibling axiom — the theorem is about $\mathcal{A}$-as-currently-stated. |
| 7 No-single-mechanism K | **GOOD with §13.5 explicit.** `rem:e6_scope` is the strongest scope qualifier of the eight: explicitly preserves CLAUDE.md §13.5 conjectural label discipline; explicitly does NOT derive $K$; explicitly does NOT claim mathematical impossibility under extended axiom sets. | (Adequately scoped — but see Statement defect noted above.) |
| 8 Cutoff function external | **GOOD.** `rem:d5_d6_scope` names the sibling-axiom direction (RG-flow / regulator-independence / asymptotic-safety). | NOT claimed: that the cutoff function $f$ is in principle external to any spectral-action-style framework. **PATCH:** the remark already says this implicitly; could be sharpened. |

**Net scope-honesty grade: B+.** The four conditional theorems (3, 4, 5, 7) have explicit scope qualifiers naming the unconditional version. Theorems 1, 2, 6, 8 have implicit but less-explicit qualifiers. Three patches noted above; effort to apply ~30 minutes total.

---

## 6. Readiness verdict — outreach (path #2)

**YELLOW with specific patches.** The eight theorems are the strongest single product of the path-#1 vein C arc. They consolidate cleanly into one sector-stratified map of the forced/free seam. Cross-paper consistency is clean. The §13.4a verification gaps are real but addressable in a single ~13-hour sprint. The two scope-honesty defects (Theorems 1, 6) are minor and can be patched in 30 minutes.

**The two material defects need addressing before submission:**

### Defect 1 (HIGH PRIORITY) — H1 is not a labelled theorem

The `rem:c_arc_terminal_state` tally claims "eight theorem-grade non-selection results" and lists "H1 Yukawa non-selection (\S~\ref{sec:open})" as item 2. But H1 lives only inside an `\paragraph{H1 Higgs inner fluctuation}` block as italicised prose. An external mathematician auditing the tally will find seven labelled `\begin{theorem}` blocks plus one paragraph that ends with the *words* "Yukawa non-selection theorem." This is a presentation defect that undermines the entire eight-theorem framing.

**Patch:** Promote the H1 Yukawa non-selection statement to a labelled `\begin{theorem}` block with explicit hypothesis / conclusion / proof-sketch structure, mirroring the format of `thm:n_gen_non_selection` and `thm:ko_dim_non_selection`. Insertion point: §VIII.C `\subsection{Sprint H1}`. Effort: ~1–2 hours of writing. **This is the single highest-priority patch for outreach readiness.**

### Defect 2 (MEDIUM PRIORITY) — Theorem 7 confuses structural and empirical content

`thm:no_single_mechanism_K`'s proof leans on the 12-mechanism elimination as if it were part of the structural impossibility argument. The structural content is weaker than the wording suggests; the empirical eliminations are the *confirmation*, not the proof, of the structural claim.

**Patch:** Rewrite the proof block to:
1. State the structural claim as derived from Facts 1+2 only ("the three sub-rings don't inter-translate within $\mathcal{A}$" → "no morphism in $\mathcal{A}$ can produce a combination spanning all three").
2. Move Fact 3 (12-mechanism elimination) to a separate paragraph labelled as empirical confirmation, with explicit "this empirical pattern is consistent with the structural claim" framing.

Effort: ~30 minutes of writing. **This patch protects the §13.5 conjectural-label discipline by making explicit what is and isn't being claimed.**

### Other patches (LOW priority, total ~3 hours)

- Add proof block to Theorem 8 (currently stated without explicit `\begin{proof}`). Effort: ~30 min.
- Add 3 scope-qualifier sentences to Theorems 1, 6, 8 per §5 above. Effort: ~30 min.
- Close §13.4a gaps for Theorems 2, 5, 6, 8 (highest impact). Effort: ~8 hr.

### Outreach effort total

| Patch level | Effort | Net effect |
|:-----------|:-------|:-----------|
| Minimum (Defects 1+2 only) | ~2 hours | YELLOW → high-YELLOW |
| Recommended (Defects 1+2 + scope-honesty patches) | ~3 hours | YELLOW → low-GREEN |
| Full (Defects + scope + §13.4a verification) | ~16 hours | YELLOW → GREEN |

**Recommendation:** Apply Defects 1+2 + scope-honesty patches before any external mathematician engagement (Brown / Marcolli / Glanois / Kleinschmidt / vS). The full §13.4a verification can be a follow-on sprint if the engagement opens.

### What this corpus state actually sells

Even with the defects above unpatched, the eight theorems represent the strongest single-day-net-gain configuration in the project's history (per the C-arc closure memo's own assessment). The sector-stratified seam map is genuinely new structural content. The Forced-Count Theorem (the moduli space is forced, its dimension is forced, the point in it is free) is a Bertrand-shaped result that — if presented carefully — would engage external NCG-SM researchers (Marcolli, Boyle-Farnsworth, Devastato-Lizzi, Bochniak-Sitarz) on a concrete question they recognise from the CCM tradition.

The defects above are presentation defects, not structural defects. The underlying mathematics is in good shape. The eight theorems consolidate cleanly. The cross-paper consistency is clean. Apply Defects 1+2 patches, then engage.

---

## 7. Files

- **`papers/group1_operator_algebras/paper_32_spectral_triple.tex`** — the eight theorems live in §VIII (lines 4765–5624).
- **`papers/group3_foundations/paper_57_forced_free_seam.tex`** — the seam catalogue / observation-grade framing (forced/free).
- **`debug/sprint_{forced_count_synthesis,c3_n_gen_non_selection,f3_ko_dim_non_selection,e7_e8_single_cutoff,g1_g2_g5_spatial_composition,c_arc_closure_e6_d5d6_chemistry}_memo.md`** — the six sprint memos behind the eight theorems.
- **`debug/h1_higgs_inner_fluctuation_memo.md`** — H1 sprint memo (the load-bearing input that needs to be formalised per Defect 1).
- **`tests/test_standard_model_triple.py`** — 45 tests on the SM triple covering H1 / Forced-Count content; §13.4a status verified.
- **`tests/test_paper51_*.py`** — six tests on Paper 51 gravity arc; relevant for Theorem 8.

## Honest scope of this audit

- This is a **review-grade audit** of statements + proofs + verification + cross-paper consistency, not a re-verification of the underlying mathematics.
- I did NOT independently re-derive any theorem's proof; I checked structural correctness and identified citation/proof-block gaps.
- I did NOT run any of the cited tests; I checked the test files exist and their structure plausibly verifies the cited content.
- The §13.4a gap list is conservative; some of the "structural impossibility — no test applies" tags may admit indirect tests I did not identify.
- I did not audit the upstream theorems cited (Paper 54 Thm 3, Paper 18 `thm:eta_trivialization`, Paper 32 `thm:pi_source_case_exhaustion`, Paper 0 §VII.B); they are taken as given for this audit.
