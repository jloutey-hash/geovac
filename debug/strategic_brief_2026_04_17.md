# GeoVac Strategic Brief — 2026-04-17

## Project State Summary

The project just closed a five-sprint arc (RH-1 through RH-5 + erratum) that targeted variants of the Riemann Hypothesis on GeoVac graphs. The headline: **graph-RH on GeoVac structures is finite-size-only, not asymptotic** (Ihara zeta families cross the Kotani-Sunada bound at V ~ 30-60 on three of four families), and the one clean infinite identity on the spectral side — `D_even(s) − D_odd(s) = 2^(s-1)(β(s) − β(s-2))` — ties χ_-4 content to Dirichlet β but does not globalize to a Riemann-class functional equation (13,080 templates failed, structural obstruction identified: two Hurwitz pieces at different exponents cannot share a Γ-completion). The spin-off is genuine: SU(2) Wilson lattice gauge on S^3 = SU(2) was constructed as a non-abelian sibling of Paper 25, 26/26 tests passing, module at `geovac/su2_wilson_gauge.py`. And the Sprint-5 erratum caught a real bug — S_min's 4th decimal was wrong in Paper 28 — with three independent methods now agreeing at ≥ 80 digits on the corrected value 2.47993693803... The Paper 28 irreducibility claim survives intact (it's a Q-linear-independence statement, insensitive to 4e-4 numerical shift).

The PI's instinct to step back is correct. The RH directions are closed or closing; further chasing of Ihara-zeta variants, higher-N_max GUE statistics, or functional-equation templates is now in the territory where effort-to-insight ratio has turned sharply negative. Meanwhile, the project's core assets are stronger than they have ever been: a 40-molecule library with the composed architecture (Q^2.5 Pauli scaling, isostructural invariance across rows 2-6), a three-tier QED-on-S^3 transcendental taxonomy (Paper 28), a clean Coulomb/HO structural asymmetry theorem (Paper 24), relativistic spin-ful composed Hamiltonians at Tier 2 with fine-structure verification (Sprint 5 CP closed Li/Be honest negatives), and Paper 18's exchange-constant taxonomy which has silently absorbed an enormous amount of new structure across Phases 4B-4I and the RH sprints. The last point is the most underappreciated: Paper 18 is now the single-best candidate document in the project for a "why GeoVac matters" statement, but it has not been consolidated to reflect everything that has accumulated in it.

The bottleneck is not technical — it's synthesis and outreach. The project has accumulated more genuine structural results in the last three months than it can cleanly communicate, and several core documents (Paper 18, Paper 21 synthesis, the applied-QC demo story) are behind the current state of the code and theorems. The next productive move is to consolidate before opening a new exploration frontier. That said, there are clean new openings worth tracking (the SU(2) Wilson construction is genuinely new territory, and the last continuous-numerical surfaces in Paper 18's §IV taxonomy are ripe for algebraic attack), and there is still an outstanding Phase-3 question (production-quality quantum-computing demonstrations on real hardware) that the project has capability for but has not pursued end-to-end.

---

## Recommended Directions (ranked)

### Direction 1: Paper 18 consolidation synthesis — "The Exchange-Constant Taxonomy at v2.23"

**Principle:** Transcendental Cataloging (Principle 3), with support from Algebraic Deconstruction (Principle 2).

**What:** A consolidating rewrite of Paper 18 that absorbs the accumulated taxonomic results from Phases 4B-4I (α decomposition), Sprint 4 Tier-3 Dirac-on-S^3 results (ring R_sp, odd-zeta in spinor sector, the three-axis grid), the QED-on-S^3 three-axis taxonomy (operator order × bundle × vertex topology from Paper 28), and the RH-sprint findings (χ_-4 / β(s) content as a new motivic tier, GUE signature discriminant between spectral and combinatorial sides, α²-weighted Ihara zeta as a new polynomial object over ℚ(α²)[s]). The existing Paper 18 has grown organically; the taxonomy has outgrown the current §IV structure. This sprint is a restructure + gap-filling: lay out the classification axes clearly (operator order, bundle type, vertex topology, graph/spectral side, even/odd zeta), populate the grid with what is known, and flag what is not. Include a synthesis theorem: *operator order is the primary transcendental discriminant; bundle type and vertex topology modulate the multiplicities but not the motivic weight class*. Paper 28 already contains this in condensed form for QED; Paper 18 should absorb it as the general statement.

**Why it matters:** Paper 18 is the framework's identity-critical document. Everyone who reads GeoVac closely reads it, because it's where the project turns from "we have a discretization" into "we have a classification of where transcendentals come from." It is currently behind the state of the code by a significant margin, and the three-axis taxonomy that emerged from QED and RH sprints is the cleanest statement of the project's organizing principle yet. Consolidating it is high-leverage: one sprint produces a document that every future paper can cite and that encodes the project's accumulated institutional knowledge in one place. It also directly serves the PI's interest in "understanding continuous structures of QM" — the taxonomy IS that understanding, written down.

**Positive result looks like:** A rewritten Paper 18 with a clean three-axis classification table, all structural identifications from Phases 4B-4I and Sprints 1-5 catalogued in the appropriate cells, a new §IV synthesis theorem, and pointers to every test that verifies each tier. Approximately 12-18 pages.

**Negative result looks like:** The consolidation exposes a contradiction between claimed tiers (e.g., the χ_-4 tier from Sprint 3 cannot be cleanly placed in the existing even/odd operator-order classification). This would itself be informative — it would reveal that the taxonomy needs a new axis, which is valuable.

**Risk assessment:** Very low-risk, high-expected-value. This is the single most obviously valuable thing the project can do right now. Overlap with dead ends: none — Paper 18 consolidation has never been attempted at this depth. The only risk is that the result is felt to be "just paper work" rather than new science, but the Reviewer / Decomposer agents should push back on superficial consolidation by requiring the synthesis theorem to be proved or discharged to known results.

**Estimated effort:** 1-2 sprints. This is a writeup-focused sprint but not a trivial one; the Decomposer agent should run in parallel to audit each claimed tier and each cross-reference.

**Prerequisites:** None beyond the current state. A Reviewer pass on the current Paper 18 before starting would sharpen the rewrite scope.

---

### Direction 2: Paper 30 draft — SU(2) Wilson-Hodge non-abelian gauge structure on S^3

**Principle:** Natural Geometry Search (Principle 1), supported by Algebraic Deconstruction (Principle 2).

**What:** Write up the SU(2) Wilson gauge construction as Paper 30, the non-abelian sibling of Paper 25. Sprint 4's RH-Q track produced the concrete construction: on S^3 = SU(2), the lattice gauge structure has Wilson links carrying SU(2)-valued connection, plaquette actions match the S^3 Haar-measure discrete Yang-Mills, and in the diagonal (U(1)-reduction) limit the construction reproduces Paper 25 exactly. The edge Laplacian L_1 = B^T B appears as the kinetic term at weak coupling (matching Paper 25's "edge Laplacian = gauge propagator" reading). The 26/26 tests in `geovac/su2_wilson_gauge.py` establish the structure. Paper 30 would: (i) state the construction formally (simplicial 1-complex + SU(2)-link data + plaquette action + Haar measure), (ii) prove the U(1) reduction limit = Paper 25, (iii) document the matter-free spectrum (Laplacian eigenvalues, Wilson-plaquette expectation values in the strong- and weak-coupling limits at N_max ≤ 5), (iv) compute the χ(Γ) = Euler characteristic and Betti numbers in the non-abelian Hodge decomposition, (v) flag what is not yet done (matter coupling, asymptotic-freedom probe, continuum limit).

**Why it matters:** Paper 25 reframed the α combination rule inside a lattice-gauge vocabulary but explicitly flagged SU(2) as an open direction (§VII.A). Paper 30 discharges that open direction and makes the project's lattice-gauge reading substantially stronger — the non-abelian case is what physicists think of when they hear "lattice gauge theory," and having a clean non-abelian discretization on a geometry known to discretize QM is a rare combination. It's also a natural-geometry finding: S^3 = SU(2) is structurally the simplest non-abelian compact Lie group, and the project has independently arrived at it from three directions (Fock projection, Bargmann-Segal, now Wilson gauge). The convergence is an observation worth naming.

**Positive result looks like:** Paper 30 drafted (~25-35 pages), all theorems proved symbolically in sympy or cited from Wilson-lattice standard texts, diagonal-limit reduction to Paper 25 verified to machine precision, 26 existing tests preserved + new tests for plaquette action and Haar-measure expectations. Explicit flagging of the matter-coupling gap.

**Negative result looks like:** The diagonal-limit reduction fails under some mild closer examination (e.g., a phase convention in Paper 25 that was U(1)-ambiguous but SU(2)-determined doesn't match). This would be a real physics finding and would be fed back into Paper 25.

**Risk assessment:** Low-risk. The construction exists and passes tests; writeup is discharge, not discovery. The main risk is that the paper ends up being observational rather than claiming new physics — which is fine, since Paper 25 set the template for this style.

**Estimated effort:** 1-2 sprints for the draft; potentially 3 if the continuum-limit probe (eigenvalue scaling under N_max refinement) reveals something unexpected.

**Prerequisites:** Paper 25 in mind; Hodge-1 and Paper 22 sparsity theorem available for cross-reference.

---

### Direction 3: Algebraic attack on the last continuous-numerical surfaces — Z_eff screening, rho-collapse spline, Level-3/4 radial potentials

**Principle:** Algebraic Deconstruction (Principle 2), with support from Transcendental Cataloging (Principle 3).

**What:** §12 of CLAUDE.md maintains the Algebraic Registry. The remaining "numerical-required" surfaces are few and well-identified: the Level-3 hyperradial potential V_eff(R) at l_max ≥ 1 (μ(R) is algebraic over ℚ(π,√2) per Track P1, but the integrals against the basis are non-elementary — radical obstruction at l_max=1, Cardano at l_max=2), the Level-4 angular eigenvalue sweep μ(ρ) (piecewise-smooth — Track S showed no global P(ρ,μ)=0 exists), and the rho-collapse spline cache. For composed geometries, Z_eff screening is already algebraic for Z ≤ 3 via the Laguerre spectral expansion (Track N/R) but falls back to spline at Z ≥ 4 due to polynomial cancellation. The sprint would: (a) attempt a Stieltjes-moment form for Z_eff at Z ≥ 4 by separating the polynomial-cancellation transcendental explicitly (probable outcome: a second exponential-integral seed e^{a}·E_1(a) at a different exponent), (b) characterize the radical structure of V_eff at l_max=1 (the square-root Δ(R) has known algebraic form; the integral against Laguerre basis may telescope into an elliptic integral, which is a documented new transcendental tier), (c) document μ(ρ) piecewise-smooth structure more precisely — is it smooth at the join points (the max/min boundaries), or only C^0?

**Why it matters:** Every time GeoVac has pushed harder on an "obviously numerical" surface, an algebraic structure has turned up (Gaunt integrals, Neumann expansion, split-region Legendre, Stieltjes seeds for Z_eff at Z ≤ 3). The pattern is robust enough to be a research heuristic. The payoff when it works is double: the exact arithmetic gives machine-precision accuracy and identifies exactly which transcendentals are intrinsic to the physics vs. which are quadrature artifacts. The current targets are the *last* continuous numerical surfaces in the main pipeline — closing them completes the "what GeoVac is algebraic about" story to its logical endpoint for first-row chemistry.

**Positive result looks like:** At least one of the three targets resolves cleanly — either to a new exchange-constant tier (naming a new transcendental that enters the physics for structural reasons) or to an algebraic form with identified seeds. The Algebraic Registry in §12 is updated to reflect the closures.

**Negative result looks like:** All three targets are confirmed numerical-required, with proofs. This would itself close a research direction and should be documented in §12 as "proven irreducible."

**Risk assessment:** Medium. The pattern of algebraic-structure-hidden-underneath is real, but the specific targets are the ones that have resisted multiple previous attempts — they are at the residual after easier wins have been taken. The elliptic-integral outcome for V_eff at l_max=1 is plausible but unverified; if it works, it's a named new tier. The Z ≥ 4 Z_eff attempt is the highest-probability win.

**Estimated effort:** 1-2 sprints, with Decomposer agent support. The Decomposer should draft a detailed plan for each target before code is written.

**Prerequisites:** Paper 18 for the taxonomy framing; §12 Algebraic Registry as the scoreboard; Tracks N, O, P1, P2, R, S for context.

---

### Direction 4: Phase 3 hardware VQE — end-to-end production demo on IBM Quantum

**Principle:** None of the three research principles directly; this is a benefit-to-the-community direction (mission goal d).

**What:** The project has the composed-qubit library, OpenFermion/Qiskit/PennyLane exports, and the ecosystem-export pipeline (`geovac.ecosystem_export`). The IBM Quantum demo (Track BC) converges H₂ to ~13 mHa on Aer simulator but has not been run end-to-end on real hardware to produce publishable results. Phase 3 is stated as the active frontier in CLAUDE.md §1.6, but the actual hardware arc has been paused. The sprint would: (i) pick 3-5 target molecules from the 40-molecule library (H₂, LiH composed, BeH composed, and probably CaH composed for the relativistic story), (ii) implement the hardware-efficient ansatz + error mitigation pipeline, (iii) run on IBM hardware at current noise levels (or, if hardware access is blocked, run on a high-fidelity noise-modeled Aer simulator calibrated to current IBM backend specs), (iv) report chemical accuracy vs. classical FCI reference, resource costs (Pauli terms, shot budgets, circuit depths), and head-to-head comparison with published Gaussian VQE results (e.g., Kandala et al. for BeH₂).

**Why it matters:** GeoVac's quantum-simulation pitch (51×-1,712× fewer Pauli terms than Gaussian baselines, O(Q^2.5) scaling, 190× advantage for LiH at Q ~ 30) is currently a theoretical resource-estimation argument. It needs at least one end-to-end demonstration to be taken seriously by the quantum-computing community. This is the project's single biggest outreach lever, and the composed library is ready for it. The PI's "benefit the scientific community" mission goal points directly at this work — if GeoVac is used anywhere outside this project in the next 18 months, it will be because of a hardware VQE demo paper, not a theorem paper.

**Positive result looks like:** A draft companion paper to Paper 20 ("Hardware-validated VQE with GeoVac composed Hamiltonians") reporting on 3-5 molecules, with chemical-accuracy results on real IBM hardware (or a carefully documented equivalent). The `demo/ibm_quantum_demo.py` module promoted to production with documentation.

**Negative result looks like:** Either (a) the hardware noise at current levels prevents chemical-accuracy VQE convergence for the target molecules — in which case the paper becomes a resource-requirements-for-future-hardware paper, still publishable, or (b) the composed Hamiltonians exhibit some unanticipated issue on real hardware (e.g., measurement overhead for the QWC groups is worse on hardware than in theory). Either would be informative.

**Risk assessment:** Medium-to-high. Technical risk is lower than the research risk: the infrastructure exists, so execution is engineering, not theorem-proving. The research risk is that the noise-mitigation state-of-the-art in April 2026 may not be good enough to show a clean advantage over Gaussian VQE on currently-accessible hardware. The PI has no quantum-hardware access of their own (per user context), so this may require either cloud IBM access (free tier has limits) or a collaboration. **This is the direction I'd flag most strongly as "someone should do this" rather than "we can do this alone"** — but it's important enough that identifying a collaborator is itself a valid sprint goal.

**Estimated effort:** Multi-sprint (3-5). Not a single-session task. Collaboration-sensitive.

**Prerequisites:** A path to hardware runs (IBM Quantum Experience free tier, or a collaborator with access). `geovac-hamiltonians` package polished for external use. Clear chemical-accuracy baselines computed classically for comparison.

---

### Direction 5: Phase 4 nuclear deepening — three-body systems and composed nuclear-electronic at scale

**Principle:** Natural Geometry Search (Principle 1) for the three-body structure; Algebraic Deconstruction (Principle 2) for the Moshinsky-Talmi + hyperradial machinery; Transcendental Cataloging (Principle 3) for the universal-vs-Coulomb-specific partition.

**What:** Paper 23 produced deuteron (16 qubits, 592 Pauli) and He-4 (16 qubits, 712 Pauli). Track NI produced the composed nuclear-electronic PoC for deuterium (26 qubits, 614 Pauli, 10^13 dynamic range). The next natural step is a three-body light-nucleus system — tritium ³H (2n+1p, would demonstrate the Minnesota NN interaction at n_max ≥ 2 with two-species JW at larger scale) or ⁶Li (3n+3p, closed-shell-light) — plus the composed-nuclear-electronic story for tritium itself (the T ion, T₂ molecule, or TH). This serves three complementary purposes: (a) extends Paper 23's methodology to a larger system and tests scaling, (b) tests the composed-architecture block decomposition across the nuclear/electronic scale separation at a slightly larger dynamic range, (c) potentially identifies whether the 10^13 coefficient-hierarchy issue in Track NI scales gracefully or breaks at larger systems. The Paper 24 π-free-HO certificate extends directly to tritium since it's a nuclear-HO-basis problem; verifying π-free-ness at n_max ≥ 3 for three-body is a cleaner statement than the two-body case.

**Why it matters:** Paper 23 established that GeoVac's graph methodology transfers to nuclear systems with the HO basis, and Paper 24 proved the Fock projection rigidity theorem showing this transfer is a specific structural claim, not a generic reuse. Extending to three-body is where the nuclear story gets interesting physically: the NN+NNN interaction structure in light nuclei is where modern nuclear physics has its own open problems, and GeoVac's qubit-encoding sparsity arguments might transfer. The three-body extension also tests the composed architecture in a regime (nuclear scales) very different from the chemistry scales where it was developed — confirming or breaking the "block decomposition is universal" working hypothesis is high-value.

**Positive result looks like:** A tritium (or ⁶Li) qubit Hamiltonian with Pauli-term count and 1-norm reported, π-free certificate at n_max ≥ 2 verified, and either a composed-nuclear-electronic extension (TH or similar) or a clear identification of why it's infeasible at current scale. A Paper 23 follow-up (Paper 23b or Paper 30b) could draft naturally from this.

**Negative result looks like:** The three-body NN+NNN interaction structure requires quadrature that breaks the clean Moshinsky-Talmi rationals (e.g., NNN depends on hyperradius in a way that doesn't factorize through Bargmann-Segal). This would be a real physics finding about where the framework's clean transfer to nuclear systems ends.

**Risk assessment:** Medium. Paper 23's Minnesota-NN + Moshinsky-Talmi machinery is proven for two-body; three-body is a known step up in angular-recoupling complexity (requires 9j symbols rather than 6j). The infrastructure largely exists (`geovac/nuclear/ho_two_fermion.py`). Overlap with dead ends: none — three-body nuclear work is new territory.

**Estimated effort:** Multi-sprint (2-3). The 9j-symbol machinery is the main technical risk item; cache it once and the rest follows Paper 23's template.

**Prerequisites:** Paper 23 in context; Moshinsky-Talmi module (`geovac/nuclear/`); standard 9j symbol implementations (sympy has them).

---

## Connections & Observations

**The three-axis taxonomy is more unified than Paper 18 currently shows.** The recent sprints have independently converged on the same three-axis structure:

1. **Paper 28 QED three-axis:** operator order (1st / 2nd) × bundle type (scalar / spinor) × vertex topology (no-vertex / even-parity / odd-parity).
2. **Phase 4B-4I α decomposition three-axis:** scalar-vs-spinor × finite-trace vs infinite-Dirichlet vs boundary-product × even-ζ-vs-odd-ζ content.
3. **RH sprint three-axis:** graph side (combinatorial, Poisson-distributed zeros) vs spectral side (Dirichlet-character-modulated, GUE-distributed zeros) × even / odd Dirichlet shift × abelian / non-abelian Wilson structure.

These three trichotomies are likely the same abstract structure in three instantiations. Paper 18 consolidation (Direction 1) should recognize this and state it directly. If the unifying statement holds, it would be one of the cleanest "organizing principle" results the project has produced — potentially rising to a Core paper claim. If it doesn't hold (i.e., the three triplets are genuinely different), that is informative too and should be said explicitly.

**The SU(2) Wilson construction connects Paper 25 to Paper 23 through an unexpected route.** Paper 25 frames the GeoVac Hopf graph as a lattice-gauge structure with the edge Laplacian as gauge propagator. Paper 23 uses Moshinsky-Talmi brackets + HO basis for nuclear systems. The SU(2) construction in `su2_wilson_gauge.py` is Wilson gauge on S^3 = SU(2). But S^3 is also where Paper 23's Fock-projection-rigidity theorem lives (the theorem says -Z/r is the unique potential whose one-electron Hamiltonian Fock-projects onto S^3). So the same manifold is simultaneously the spin-structure carrier (Tier-2 spinor composed), the gauge-bundle base (Paper 25 + SU(2) extension), and the Coulomb-specific projection image (Paper 23 rigidity). This four-way coincidence is not an accident — it's S^3 = SU(2) = compact simply-connected non-abelian group of rank 1 being structurally the simplest object in three different categories. Naming this explicitly (probably in Paper 30, possibly in Paper 18 consolidation) would tighten the framework's rhetoric.

**There is a hidden "no-common-generator" theme across the failed derivation attempts.** Phase 4G α-K formally documented that K = π(B + F − Δ) has B, F, Δ with categorically different origins (finite Casimir, infinite Dirichlet, boundary product — no unifying generator). Sprint 4 RH-M found that D(s)'s zeros have GUE CV while Ihara zeros have Poisson CV — categorically different. Sprint 3 RH-J / RH-K found χ_-4 content on spectral side but not on Ihara side. Sprint 4 functional-equation hunt failed because "two Hurwitz pieces at different exponents can't share a Γ-completion" — a no-common-generator statement. The project has repeatedly found that structurally-independent objects coincide numerically without unifying underneath. This pattern may itself be the project's deepest structural finding — and naming it explicitly in the Paper 18 consolidation (Direction 1) as a *meta-observation* would be a substantive contribution.

**Paper 21 synthesis update is a natural companion to Direction 1, not a separate direction.** Paper 21 is the project synthesis. Its last real update was well before the Phases 4B-4I α arc and the RH sprints. If Direction 1 (Paper 18 consolidation) produces the new axis-structure cleanly, Paper 21 §V (exchange constants) and §VI (Hopf α conjecture) both need to be refreshed to match. This is a natural downstream sprint. I did not list it as a separate Direction because doing Paper 21 before Paper 18 risks baking in a stale taxonomy. Order matters: Paper 18 first, then Paper 21.

**The "RH is closed" verdict isn't quite complete.** Sprints 1-5 closed *graph-RH* (Ihara, Kotani-Sunada) and *spectral-RH variants* (D(s) functional equation, depth-2 MZV attacks). What they did not close is whether there is a coupled *matter-gauge* RH statement on GeoVac, where the SU(2) Wilson construction provides the gauge side and Direction 5 (three-body nuclear) or Direction 4 (composed-qubit Hamiltonians as matter) provides the matter. This is one of the "open residual territory" items. I have not listed it as a separate Direction because it's too speculative for the PI's "step back to core ideas" framing, but I want it in the record. If the SU(2) Wilson construction (Direction 2) gains fermion matter coupling in a future sprint, the matter-gauge RH question would naturally re-open, and at that point it might be worth revisiting.

**The PI's "core ideas" language points at Paper 18 and Paper 25 + SU(2), not at new exploration.** Reading the PI's note carefully — "millennium prize directions are closed or mostly closed, which seems good in a sense. better to focus back on our core ideas" — the implicit preference is consolidation + proven territory over opening new speculative fronts. This is why Directions 1 and 2 are ranked above Directions 3, 4, 5. Direction 3 (algebraic attack) is core-principle-aligned and relatively low-risk but speculative in outcome. Direction 4 (hardware VQE) serves the mission but is partly outside-project-capability. Direction 5 (nuclear extension) is a natural frontier but is opening new territory rather than consolidating. The ranking reflects the PI's stated preference.

---

## What I Don't Know

**I don't know whether the three-axis taxonomy unification (discussed in Connections above) actually holds, or whether it's a pattern-matching illusion.** My recommendation for Direction 1 (Paper 18 consolidation) is contingent on this being discharged in the sprint — if it's an illusion, Direction 1 still produces value as consolidation, but the headline claim weakens.

**I don't know the PI's current time budget or availability.** Direction 1 is 1-2 sprints (doable in a few sessions). Direction 4 is multi-sprint and collaboration-sensitive. Direction 5 requires 9j-symbol machinery that I haven't audited. If the PI has one session of time, Direction 1 is the right answer; if they have four weeks, the right order is probably 1 → 2 → 3 with 4 in parallel if collaboration can be arranged.

**I don't know whether Paper 25 has been through a Reviewer pass in its current state.** Its Sec VII (open questions) is what motivated the SU(2) sprint and should be cross-checked against the actual SU(2) construction before Paper 30 is drafted. A Reviewer dispatch on Paper 25 at the top of Direction 2 would catch any wording mismatches.

**I don't know if there's an external deadline or collaboration-outreach event that changes priorities.** If there's a conference submission deadline or an outreach opportunity that favors a hardware result, Direction 4 moves up in ranking. If there's a theoretical-physics audience that wants the lattice-gauge story, Direction 2 moves up. The PI has context I don't.

**I don't know whether the `geovac/su2_wilson_gauge.py` 26/26 tests include the U(1)-reduction-to-Paper-25 match, or whether that's claimed but untested.** This is verifiable in ten minutes but I did not do so for this brief — it's a precondition check for Direction 2 and should be the first thing the PM verifies.

**I don't know if the hardware-VQE infrastructure (Direction 4) has known gaps I'm not aware of.** I noted `demo/ibm_quantum_demo.py` exists, but I haven't audited whether the noise-mitigation pipeline, the ansatz, or the measurement grouping are production-ready or are prototypes. A scoping sub-agent would be needed at the start of Direction 4 to audit readiness.

**I don't know how the PI weighs "paper count" vs "paper depth."** Direction 1 (Paper 18 consolidation) updates an existing paper; Direction 2 (Paper 30) adds a new paper. If the PI values paper count, Direction 2 moves up; if they value depth on the existing catalog, Direction 1 is the clearer answer. My ranking assumes depth-preferred, which is consistent with the "back to core ideas" framing but is an interpretation.

**I don't know whether the Sprint-5 erratum patch closed all downstream dependencies.** The erratum fixed the S_min numerical value in Paper 28 + the driver script. If any other paper or module transitively consumed the old S_min = 2.47953 value (e.g., any of Paper 28's cross-references, or a cached benchmark in `debug/data/`), those would need checking. A quick grep-pass over the papers and the benchmark registry for the old value would catch any stragglers.
