# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v4.12.1 (June 14, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (unit S3) that is mathematically equivalent to the Schrodinger equation via Fock's 1935 conformal projection. This equivalence is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Authoritative source rule:** The papers in `papers/group1_operator_algebras/`, `papers/group2_quantum_chemistry/`, `papers/group3_foundations/`, `papers/group4_quantum_computing/`, `papers/group5_qed_gauge/`, `papers/group6_precision_observations/`, and `papers/synthesis/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it. (Papers were reorganized from the previous `core/`, `methods/`, `applications/`, `synthesis/`, `standalone/`, `observations/`, `conjectures/` layout into six audience-targeted groups on 2026-05-22.)

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator (Josh Loutey — `J.~Loutey` on the papers) provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers across the six audience-targeted group folders (see §6) are written to academic standards but are not submitted to traditional journals. The project's viability case rests on its corpus of verified structural results — discretization- and encoding-structure theorems, honest negative results, and the benchmarked computational artifact (pip-installable, 40-molecule library) that demonstrates them. The tool is a research instrument, not a production-chemistry replacement (see the benchmarking rule in §1.5 and `docs/claims_register.md`). Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks. (Viability-case sentence reworded 2026-06-10 per PI direction, accessibility plan Phase 3.)

---

## 1.5. Positioning & Framing

GeoVac is a discretization framework that exploits the natural geometry of separable quantum systems to produce sparse graph Hamiltonians. It is **not** a new foundation for quantum mechanics — it works because the graph Laplacian converges to the known continuous Laplace-Beltrami operators in the continuum limit (Paper 7).

**Rhetoric rule:** The GeoVac framework is demonstrably conformally equivalent to standard continuous quantum mechanics at every level where it has been tested (Paper 7 provides 18 symbolic proofs for the S³ case; subsequent papers verify operator convergence for each natural geometry). The papers should present this conformal equivalence as the primary result. The discrete graph topology and the continuous Laplace–Beltrami operators are dual descriptions connected by proven conformal maps — the mathematics supports both readings, and the interpretation of which is more fundamental is left as a choice for the reader. Avoid language that asserts ontological priority of either description. Lead with concrete computational results (sparsity, scaling, accuracy, structural insight) rather than interpretive claims about the nature of quantum mechanics.

**Lead with concrete advantages:** O(V) sparsity, angular momentum selection rules baked into the basis, zero-parameter construction from nuclear charges and geometry alone, efficient qubit encodings (Paper 14). These structural properties are the framework's actual selling points — not philosophical claims about the nature of quantum mechanics.

**The π-free graph principle:** The GeoVac graph is π-free: all eigenvalues, degeneracies, and coupling coefficients are integers or rationals. Transcendental numbers (π, exponential integrals, spectral zeta values) enter exclusively when projecting the graph onto continuous manifolds for comparison with experiment or for computational convenience. This observation (Paper 18) motivates the design principle: stay on the graph whenever possible, and when projection is necessary, identify the minimal transcendental content (the exchange constant) required. The exchange constant taxonomy (intrinsic, calibration, embedding, flow) classifies projections by what determines them and predicts the computational cost of each departure from the graph.

**Quantum simulation positioning:** The classical solver investigation (v2.0.6-23) characterized what the natural-geometry hierarchy can and cannot achieve for ground-state PES; the primary computational value proposition is quantum simulation. The composed architecture produces qubit Hamiltonians with O(Q^2.5) Pauli scaling (vs O(Q^3.9-4.3) Gaussian), 51x-1,712x fewer Pauli terms across LiH/BeH2/H2O, structural Gaunt sparsity compatible with all downstream optimizations (tapering, grouping, tensor factorization), and a 40-molecule library across three periodic-table rows with N_Pauli = 11.11 x Q universal for composed builds. Position for VQE/NISQ (Pauli count, QWC groups) rather than fault-tolerant QPE; LiH electronic-only 1-norm matches STO-3G at 0.97x with 2.7x fewer Pauli terms. Full track-by-track chronicle: `docs/development_frontier_archive.md` + CHANGELOG.

**General composed builder foundation:** Paper 16 (Chemical Periodicity as S_N Representation Theory) provides the group-theoretic foundation for the general composed builder's atomic classification. The structure types (A/B/C/D/E), the universal angular quantum number ν = N−2, and the recursive core-valence decomposition are all derived from S_N representation theory and implemented in `atomic_classifier.py`.

**PK limitation and research path:** The composed framework's PK pseudopotential is the accuracy bottleneck (5.3-26% R_eq error); six modification attempts failed (Section 3). The balanced coupled Hamiltonian (Papers 19/20) is the PK-free alternative: cross-center V_ne via multipole expansion with exact Gaunt termination; LiH 4e FCI reaches 0.20% energy at n_max=3 with structural R_eq drift (~8.8%) — energy converges excellently, geometry drifts; optimal for single-point quantum simulation at fixed geometries. Resource anchor: composed LiH 334 Pauli vs STO-3G 907 vs cc-pVDZ 63,519. Chronicle: `docs/development_frontier_archive.md`.

**Benchmarking rule:** When comparing to other methods, always use the strongest available baseline (cc-pVTZ or better for atoms, explicitly correlated methods for molecules), not just STO-3G. If the comparison is unfavorable, say so honestly and identify what the framework offers instead (sparsity, scaling, structural insight). Position the framework as a computationally principled alternative, not as a replacement for production quantum chemistry.

---

## 1.6. Project Phase

**Phase 1 (v0.9.x-v1.x): Foundation.** Graph Laplacian equivalence proof, natural geometry hierarchy, LCAO diagnostic arc, bond sphere theory, 18 symbolic proofs (Paper 7).

**Phase 2 (v2.0.0-v2.0.23): Classical solver investigation.** Systematic exploration of all solver architectures (adiabatic, coupled-channel, 2D variational) across all levels (2-5, 4N). 30+ completed tracks, 40+ documented negative results. Key outcomes: H2 at 96.0% D_e (Level 4), LiH R_eq at 5.3% (Level 5 composed), full N-electron equilibrium without PK (Level 4N). Investigation complete: all solver x PK x basis combinations exhausted, structural accuracy ceilings characterized, exchange constant taxonomy established (Paper 18).

**Phase 3 (v2.0.24-v2.7.0): Quantum simulation.** The composed architecture's structural sparsity (O(Q^2.5) Pauli scaling, block-diagonal ERIs) is the framework's primary computational advantage. Classical solver results serve as validation benchmarks for the quantum Hamiltonians. Next steps: quantum resource estimation, hardware-aware circuit compilation, experimental collaboration.

**Phase 4 (v2.7.0+): Nuclear extension and framework delineation.** Extension of the GeoVac framework to nuclear shell model Hamiltonians using the hyperspherical (HO + spin-orbit) basis, together with a precise delineation of which parts of the framework transfer to non-Coulomb systems. Key outcomes: (1) the potential-independent angular sparsity theorem (Paper 22) — ERI density depends only on l_max, not on V(r), verified at l_max=3 as 1.44%; (2) the deuteron and He-4 qubit Hamiltonians (Paper 23, Tracks NE/NF) using the Minnesota NN potential, Moshinsky-Talmi brackets, and a two-species tensor-product JW encoding; (3) the Fock projection rigidity theorem (Track NH) — the S^3 conformal projection is unique to -Z/r and does NOT transfer to HO or Woods-Saxon; (4) the composed nuclear-electronic deuterium PoC (Track NI) demonstrating the block architecture at 26 qubits with hyperfine validation; (5) the Bargmann-Segal lattice (Paper 24, Track NK) — discrete graph encoding of the 3D HO on the holomorphic sector of S^5, bit-exactly π-free in exact rational arithmetic at every N_max (verified at N_max=5: 56 nodes, 165 edges, zero irrationals); (6) the universal vs Coulomb-specific partition as the Phase 4 conceptual result, completed by the Coulomb/HO asymmetry analysis in Paper 24 (calibration π is structurally Coulomb-specific — tied to second-order Riemannian operators with nonlinear projections, has no analog for the HO's first-order complex-analytic projection). Track NJ memo definitively shelved the nuclear → alpha connection (originally "shelve," upgraded after Sprint 2/Paper 24 to "definitively shelved" based on the structural impossibility, not just empirical absence). HO rigidity theorem (Theorem 3 of Paper 24) is the structural dual of the Fock rigidity theorem.

---

## 1.7. Working Hypotheses (Internal Register)

This section is a **bold-claim register**, distinct from the rhetoric of the papers. Papers remain cautious under §1.5 (dual-description framing, no ontological priority); this register is what sub-agents and the PM may *reason from* during synthesis work. Nothing here appears in papers unless promoted after its falsifier clears. Full status chronicles (the sprint-by-sprint evidence trail, April–June 2026) live in `docs/wh_register_history.md`; this section holds only claim, falsifier, and current status.

**Governance:**
- PM may update a WH's "Status" line based on sprint evidence — REPLACE the line and move superseded text to `docs/wh_register_history.md`; never append (§13.11 rule 9).
- Adding, retiring, or promoting a WH to paper-level claim requires explicit PI direction.
- Retired WHs move to the history doc with rationale; never silently deleted.
- No WH is a license to bypass the rhetoric rule in papers or the verification gates in §13.4.

---

**WH1 — GeoVac is an almost-commutative spectral triple.** A = functions on the Fock-projected S³ graph; H = scalar/spinor state space; D = Camporesi–Higuchi Dirac; non-abelian gauge structure enters as inner derivations of the almost-commutative extension A ⊗ M_n(ℂ) (Marcolli–van Suijlekom lineage; Papers 25/30/32).
*Falsifier:* a GeoVac observable demonstrably inconsistent with any spectral-action expansion, or a violated structural axiom (order-one, reality).
*Status:* **PROVEN — unconditional (2026-06-10).** Paper 38: the discrete truncations converge to the round-S³ spectral triple in van Suijlekom's state-space GH distance at rate (4/π + o(1))·log n_max/n_max, on the truthful CH substrate (translation-seminorm metrization; frozen falsifier `tests/test_p38_action_seminorm.py`). The Lorentzian extension (Papers 45–49) is DESCOPED (P45 annihilation theorem); repair path = Toeplitz temporal compressions (see WH7).

**WH2 — Paper 18 is the Seeley-DeWitt + ζ-invariant decomposition of this spectral triple.** The transcendental taxonomy is the structured output of spectral-action geometry, organized by operator order × bundle type.
*Falsifier:* a transcendental in a GeoVac observable that cannot be placed in the grid.
*Status:* three of four axis-quadrants filled; mixed-Tate sharpening POSITIVE (2026-06-03) — M2 on S³ sits in the pure-Tate sub-ring ⊕_k π^{2k}·ℚ (Fathizadeh–Marcolli inherited). See `debug/sprint_mixed_tate_test_memo.md`.

**WH3 — The lattice exists a priori; match to physics is evidence, not derivation.** The packing construction is independent of known physics; its persistent match (Fock S³, nuclear magic numbers, Dirac fine structure, Pauli sparsity) is evidence physics is hosted by a discrete spectral triple, not that the lattice was reverse-engineered.
*Falsifier:* the packing construction fails to force the (n, l, m, s) structure or the n²−1 spectrum; or a match turns out to depend on a hidden physics-informed parameter.
*Status:* origin-story framing permits strong ontological claim internally; papers stay under §1.5.

**WH4 (deflated 2026-05-07) — The four-way S³ unity is one Fock-projection statement plus three forced consequences.** Bertrand + SO(4) force S³; the Hopf base, the CH spinor bundle, and SU(2) Wilson all follow from S³ = SU(2) parallelizability / maximal-torus structure. Does NOT extend to inner-factor selection (Yukawas, generations remain unforced).
*Falsifier:* a construction forcing one of the four roles onto a different manifold; or a published argument decoupling a consequence from the Fock input.
*Status:* deflated to a single-input forcing statement; outer structural unity essentially closed (Sprint TS-D + Paper 38).

**WH5 — α is a projection constant, not a derivable number.** K = π(B + F − Δ) composes three structurally independent spectral objects (finite Casimir trace; Fock Dirichlet ζ(2); Dirac boundary count 1/40); the right open question is why the sum equals α⁻¹, not how to derive each piece.
*Falsifier:* a spectral-triple construction deriving K as a single coefficient of a well-defined functional.
*Status:* TWELVE mechanisms eliminated (Phases 4B–4I + Sprint A + Sprint K-CC, including the T9 algebraic obstruction: no single CC heat-kernel expansion contains B, F, Δ as terms). Standing reading: three-regime projection coincidence. Paper 2 stays in Observations; combination rule conjectural (§13.5 hard prohibition).

**WH6 — GeoVac's RH-adjacent object is the Dirac spectral zeta D(s), not classical ζ.** Internal GUE-like zero statistics (CV ≈ 0.35–0.40); the classical-RH bridge is closed by three independent walls (zeros not on one line; no spectral-triple-natural functional equation, 48 OoM; wrong Weyl class).
*Falsifier:* D(s) zeros on a single critical line at larger samples; or a natural functional equation closing the RH-O gap.
*Status:* paused (2026-04-18); if resumed, the target is D(s) and its spectral-action interpretation.

**WH7 — Time-discreteness is observer-compactification (registered 2026-06-10, PI direction).** The only temporal structure the framework can see metrically is compactified time, and compactified time is automatically discrete. Inputs: (i) Paper 35 — π enters exactly at temporal compactification (Matsubara 2πk/β); (ii) P45 annihilation theorem — ℝ-time is Lipschitz-invisible in the v1 architecture; (iii) Paper 47 — the three temporal carriers are spectrally indistinguishable. Reading: discreteness of time is supplied by the observer's compact integration window (KMS β = 2π, four-witness theorem, Connes–Rovelli thermal time) — time is the prototype free-side projection. Temporal restriction of the organizing observation below.
*Falsifier (primary):* the Toeplitz temporal-compression program — a metrically visible temporal algebra on a NON-compact carrier without compactification weakens WH7 to convention; a proven annihilation-type obstruction for the framework's whole non-compact temporal class forces it.
*Falsifier (secondary, from Paper 35):* a GeoVac observable containing π whose evaluation provably involves no temporal/spectral integration.
*Status:* REGISTERED (2026-06-10); Step-1 Toeplitz probe same day **POSITIVE-REBUILD** (`debug/sprint_wh7_toeplitz_probe_memo.md`): the translation seminorm equals the continuum Lipschitz constant exactly on Toeplitz temporal modes (err ~10⁻¹⁴); the P45 invisibility is reproduced as a momentum-diagonal-algebra artifact (input (ii) downgraded accordingly); visibility survives de-compactification at fixed bandwidth, so the primary falsifier currently leans "weakens-to-convention" on the visibility leg, pending Step-2 pointed-proper bookkeeping (follow-on B2); follow-on B1 closed same day — joint S³×S¹_T product-carrier convergence at additive rate (Paper 45 `prop:product_action_seminorm`, falsifier `tests/test_wh7_b1_joint.py`). Load-bearing content is now leg (i): only compactification makes time discrete and injects π. Sharpened reading: *time is metrically visible, but the observer's window is what makes it discrete.* Honest cap: input (iii) means discrete-vs-continuous may be empirically undecidable at every currently computable level; papers stay under §1.5.

---

**Organizing observation — discreteness is compactness (established).** The discrete, bit-exact skeleton is the closed/compact regime: compact groups have discrete spectra (Peter-Weyl), and this is the mathematical content of the S³ graph. The continuum is what remains when compactness is released. Calibration data is continuum (un-packed) data, which is why the skeleton is rational/forced and calibration is transcendental/free. Well-captured in Paper 18 §III (compactness thesis) and Paper 35 (temporal compactification injects π). The "second packing axiom" framing is retired; the 2026-05-30 confinement reframing is archived as an organizing reading (2026-05-31).

---
## 1.8. Multi-Observable Focal-Length Decomposition Program (Directive)

The May 9 r_Z extraction thread (four iterations, fully documented in §2 below) demonstrated that the framework's projection-chain machinery functions as a unique diagnostic instrument for precision physics: **a multi-observable global fit using a single structural framework exposes convention-mismatches in literature compilations and kernel approximation gaps that are invisible from any single-observable analysis.** This subsection codifies that finding as an ongoing research program with explicit targets.

**The directive.** When a sprint touches precision atomic, molecular, or nuclear observables, the PM should ask three questions in order:

1. **Is there a literature convention mismatch the framework can surface?** Different compilations (Eides, Karshenboim, Krauth, Pachucki–Yerokhin, Drake, Antognini) itemize Layer-2 corrections at sub-percent levels with different conventions for which sub-leading terms count as "Zemach" vs "polarizability" vs "recoil NLO" vs "multi-loop QED." Multi-observable global fits using GeoVac's single projection-chain dictionary expose these mismatches at numerical level (the W1b finding, 2026-05-09: ~0.01% of $\nu_F$ propagates to ~25 mfm in extracted $r_Z$, below the 30 mfm Eides-vs-lattice gap). Frameworks that consume single-observable Layer-2 inputs in their native conventions cannot see this.

2. **Is there a GeoVac kernel approximation gap the framework can identify?** Tightening precision via more observables tends to relocate the load-bearing systematic (the W1a-D and W1b iterations: cross-register V_eN was correct; kernel-leading-order at coarse precision; sub-leading recoil-mixing at 51 mfm; convention mismatch at 16 mfm). Each precision tightening uncovers a new systematic. PMs run per-observable consistency checks before accepting any global-fit result.

3. **Does the observable admit a focal-length decomposition that sharpens the §III dictionary's coverage?** Multi-component observables (Lamb shift, hyperfine, fine structure) deserve Roothaan autopsies in Paper 34 §V.C. Single-component observables get single rows in §V/§V.B. The autopsy format makes the projection-chain decomposition visible at the observable level.

**Three classes of problem the program targets.**

- **Literature convention mismatches**: differences in Layer-2 itemization between Eides 2024 / Krauth 2017 / Karshenboim 2005 / Pachucki–Yerokhin 2010 / Drake 1990 / Antognini 2013 compilations. Framework's multi-observable fits expose these structurally.
- **GeoVac kernel approximation gaps**: places where the framework's leading-order operator (Eides Zemach, Foldy–Friar contact, Friar moment) is too coarse for the precision target. Each gap is a named follow-up (W1a, W1b, etc., per CLAUDE.md §1.7 multi-focal-composition wall taxonomy).
- **General focal-length decomposition cataloguing**: §V.C Roothaan autopsies for multi-component observables. Discipline: every Layer-2 input gets a focal-length tag, every framework-native contribution gets a §III chain, the decomposition closes at sub-MHz (or sub-kHz) precision.

**Active Roothaan-decomposition targets** (queued 2026-05-09):

1. **§V.C.2 Hydrogen 21 cm hyperfine four-component autopsy** (placeholder fill). Bohr-Fermi + Schwinger $a_e$ + reduced-mass + Zemach decomposition at +18 ppm framework-residual. Tests §III.18 magnetization-density at the operator level.

2. **§V.C.3 Muonic hydrogen 2S–2P Lamb shift autopsy** (placeholder fill). Decomposition into full Uehling kernel (Antognini-style) + SE Bethe-log + Friar moment via §III.17 + Källén–Sabry two-loop VP + deuteron polarizability. Tests §III.17 + §III.18 + §III.16 in the muonic regime ($\beta = 1.475$).

3. **§V.C.4 Helium $2{}^3P$ fine structure full autopsy** (NEW). Decomposition into spin-orbit + spin-spin + spin-other-orbit (Drake combining coefficients) using bipolar harmonic $(k_1, k_2)$ decomposition. First operator-level test of Observation 3's angular compositional projection rule (internal multi-focal at $\alpha^2$). Reference: NIST + Pachucki–Yerokhin theory at sub-ppm.

4. **§V.C.5 Helium $2{}^1P \to 1{}^1S$ oscillator strength** (NEW). Multi-electron extension of Sprint Calc-L (Lyman α, +0.055% match). Tests vector-photon promotion + Wigner 3j + multi-electron correlation. Reference: Drake handbook $f \approx 0.276$. Verifies whether Sturmian's continuum-closing property (validated for hydrogen polarizability at exact at $N_\text{basis} = 2$, Sprint Calc-P) extends to multi-electron transitions.

5. **§V.C.6 Cesium $6S_{1/2}$ hyperfine (atomic clock)** (NEW, prospective). $Z = 55$ heavy-atom regime where relativistic spinor lift dominates and §III.17/§III.18 are tested at very different focal lengths than hydrogen. Atomic-clock-grade reference precision ($10^{-15}$ relative; the SI second is defined by this transition). Tests Z-scaling of the framework's projection-chain machinery in a regime where simple hydrogenic formulas break down. **Prospective angle**: heavy-atom Cs PNC (parity non-conservation) extractions have known atomic-structure uncertainties at the percent level; the framework's projection-chain decomposition might illuminate atomic-structure systematics that single-observable analyses don't expose.

These five targets together exercise §III.17, §III.18, §III.19, spinor lift (§III.7), Wigner 3j (§III.8), Wigner $D$ (§III.9), Hopf bundle (§III.2), vector-photon promotion (§III.11), spectral action (§III.6), Sturmian (§III.5), and rest-mass (§III.14) — most of the projection dictionary. Successful decomposition of all five at sub-percent framework-native precision would validate the directive's broader applicability across atomic precision physics. Where each lands also tests the directive's three problem-classes empirically: how often does multi-observable consistency reveal a literature convention issue vs a framework kernel gap vs neither?

**Sprint cadence.** 2 tracks per sprint (1 verification + 1 prospective), established in the May 9 thread. Continue this cadence for the new targets.

**Cross-references.** May 9 r_Z thread (§2 below): demonstrated all three problem classes in a single observable across four iterations. Paper 34 §V.C.1 (Lamb shift autopsy): founding instance of the cataloguing discipline. Audit memo `debug/paper34_v_base_unit_audit_memo.md`: methodological precedent for systematic analysis. Diagnostic-before-engineering memory `feedback_diagnostic_before_engineering.md`: the load-bearing discipline that made the trajectory cleanly informative.

---

## 2. Current Development Frontier

> Full sprint chronicles live in `CHANGELOG.md`. This section is a compact index. Sprint detail is in the memos linked below.

- **Honest-review: 4 roots verified (2026-06-14, v4.12.1):** Roots (7,0,32,38) content SOUND (WH1 "unconditional" verified by direct band-injectivity read; no overclaim/zombie-cites). Defects were status-drift/precision only — 11 fixes (Paper 32 WH1 label + GH-metric→vS state-space; register rows 2/17; Paper 7 κ/convergence; Paper 0 dates). Ledger `debug/honest_review_2026_06_14_ledger.md`.
- **NA-1 reconciliation + process fix (2026-06-14, v4.12.0):** Period-value irreducibility is NOT a valid Reading-A/B discriminator (irreducible≠bracket, Cartier–Milnor–Moore); JLO-Depth2 **Reading A** stands, no corpus edit. Added §9 Current-State Check + 2 standing rules. See `debug/sprint_na1_period_irreducibility_nondiscriminator_memo.md`.
- **Sprint Hodge-SL₂ (2026-06-14, v4.11.0):** GeoVac's SL₂ is NOT the Mumford–Tate group — it realizes the CM point (CM field Q(i), MT = 1-dim CM torus ⊊ SL₂); structurally explains the Hain–Brown negative; Q(i)-triangulation (Hodge=level-4=period field) spin-forced per audit. Paper 56 §sec:open_g4_hodge convention→CM-identification + KO-dim-3 correction. Tests `tests/test_paper56_hodge_sl2.py`. See `debug/sprint_hodge_sl2_memo.md`.
- **Sprint S^(4) stage-2 depth verdict (2026-06-13, v4.10.0):** realized depth ≤ 3 for S^(4) SETTLED — proven numerics-free at w=9/11/13 (stuffle spans full depth-4 space, two-prime exact rank), w=5/7 by low-weight filtration; confirms prop:depth_k at k=4. Falsifier `tests/test_s4_stage2_depth.py`. See `debug/sprint_s4_stage2_memo.md`.
- **Sprint S^(4) stage-1 (2026-06-13, v4.9.0):** k=4 rung DEFINED + rigorously bracketed S^(4) ∈ [316.443, 316.698]; exact structure mapped (15-term + o-space relations, census odd-weights 5–13); high-precision b1=2 trailing constants deferred to stage-2 symbolic (evaluator bit-exact for non-trailing only). Falsifier `tests/test_s4_stage1.py`. See `debug/sprint_s4_scoping_memo.md`.
- **Sprint S^(3) W10 identification (2026-06-12, v4.8.0):** W10 CLOSED — Charlton–Hoffman symmetry theorem (explicit products) + 4 stuffles identify all four w10 triples; the three a-doubles cancel identically in the assembly; S^(3) ∈ Q[π², ln2, ζ(odd), ζ(5,3)] explicit, sole depth-2 constant ζ(5,3) at 6π², gate 1.15e-198. Falsifier `tests/test_s3_w10_identification.py`. See `debug/sprint_s3_w10_identification_memo.md`.
- **Sprint S^(3) follow-on triple (2026-06-11, v4.7.0):** Full assembly: ALL depth-2 generators cancel (identified part ∈ depth-1 ring; ζ(11) cancels); w10 → 5-generator block via stuffle closure (saturation theorem); all w≤7 forms = verbatim Hoffman App-A rows (5 sub-agent citation errors caught by PM source verification). Falsifier extended (eq:w10_reduction). See `debug/sprint_s3_w10_symbolic_memo.md`.
- **Sprint S^(3) closure (2026-06-11, v4.6.0):** BOTH OPEN ITEMS CLOSED — S^(3) = 31.5726 in rigorous bracket [31.57063, 31.57300]; stage-1 figures 30.615/30.220 both Levin-on-log artifacts; four PSLQs ACCEPT; realized depth ≤ 2 at k=3 (≤ k−1 pattern). Falsifier `tests/test_s3_decomposition.py`. See `debug/sprint_s3_closure_memo.md`.
- **Sprint S_min identification + S^(3) decomposition (2026-06-11, v4.5.0):** S_min CLOSED — explicit Q[π²,ln2,ζ(3),ζ(5)] element; prior irreducibility = basis-coverage artifact; k=3 reaches depth-2 generators. Falsifier `tests/test_smin_decomposition.py`. See `debug/smin_dossier_round1_memo.md`.
- **Sprint B3 Phase-3 Sprint-4 band exhaustion (2026-06-10, v4.4.0):** PHASE-3 CHARTER CLOSED — 5-agent workflow: interval layer bit-exact at every window; cost layer converges in OWN power-law class (γ and thermal excluded); b-parity staircase = (−1)^{2b} grading in exhaustion dynamics; exact form ‖P_W M_{C¹}P_W‖ = √6(2j)/(2j+2). Falsifier `tests/test_wh7_band_exhaustion.py`. See `debug/sprint_wh7_band_exhaustion_memo.md`.
- **Sprint B3 Phase-3 Sprint-3b (2026-06-10, v4.3.0):** FOLD RULE CLOSED-FORM — CG lemma [RCR]=(−1)^{b+j2−j1}[C^b_{−μ′,μ}] blockwise + parity rule ε=(−1)^{b+j2−j1+μ′} (annihilation iff ε=−1, exact rational fold ratios); Sprint-3 headline facts proven j_max=1 WINDOW-EDGE ((2,1) revives 8/41 at 3/2; (2,2) stops commuting, mirror 9/25). Falsifier `tests/test_wh7_b3_fold_rule.py`. See `debug/sprint_wh7_b3_fold_rule_memo.md`.
- **Sprint B3 Phase-3 Sprint-3 (2026-06-10, v4.2.0):** ADMISSIBILITY SETTLED — folding reorganizes cone classes ((2,1) annihilated; (2,2) timelike → flow-commuting); band-limited penalties bit-exactly cutoff-INDEPENDENT n_max=2..5 (prize reduces to band exhaustion); Bures positivity REFUTED (574/2400); period-π attribution corrected (spinor grading, not folding). Falsifier `tests/test_wh7_b3_phase3_sprint3.py`. See `debug/sprint_wh7_b3_phase3_sprint3_memo.md`.
- **Sprint B3 Phase-3 Sprint-2 (2026-06-10, v4.1.0):** MECHANISM-CLOSED — evenness symmetry ([G,K]=0 ⇒ c12(ε)=c23(−ε), cost-universal) explains Sprint-1 bimodal p AND sign caveat; Umegaki chain fails generically (90/96); wedge substrate transfers bit-level, period π, operational interval ℓ. Falsifier `tests/test_wh7_b3_phase3_sprint2.py`. See `debug/sprint_wh7_b3_phase3_sprint2_memo.md`.
- **Release v4.0.0 (2026-06-10):** Major-version consolidation of the six staged units v3.110.0–v3.115.0 (repo hygiene + papers/INDEX, CLAUDE.md compaction round 2, WH7 registration + Step-1 probe, B1 product-carrier convergence, B3 Phases 1–3 Sprint 1). PI versioning convention: major versions track AI-collaborator changes — the Fable refactor of CLAUDE.md qualifies. See CHANGELOG [4.0.0].
- **Sprint B3 Phase-3 Sprint-1 state intervals (2026-06-10, v3.115.0):** FOUNDATION-LAID — wedge KMS + orbits anchored bit-level; D_max chain ("detour never costs less") 96/96 universal; bimodal kick-cost scaling by boost weight (ε² vs ε^1.2, tangent confound ruled out); caveat frozen: excess sign is reference-state dependent — D_max is the penalty layer, flow parameter is the interval. Falsifier `tests/test_wh7_b3_phase3.py`. See `debug/sprint_wh7_b3_phase3_sprint1_memo.md`.
- **Sprint B3 Phase-2 cone structure (2026-06-10, v3.114.0):** HS causal ratio = symbol classifier EXACTLY, q_F = (2m′²−b(b+1))/(b(b+1)) all-rational (CG-exact, float dev 2×10⁻¹⁶); causal form class-diagonal definite — cone is a GRADING not a signature; rate-level reverse triangle FAILS by inertia (frozen negative) → Phase 3 = state-level thermal-time intervals (Paper 49 TICI + B1 substrate). Falsifier `tests/test_wh7_b3_phase2.py`. See `debug/sprint_wh7_b3_phase2_memo.md`.
- **Sprint B3 Phase-1 boost-seminorm probe (2026-06-10, v3.113.0):** POSITIVE-STRUCTURED — boost-alone kernel = 9 (boost-invariant multipliers, the structured middle between P45 annihilation and metric condition); frame restores ℂ1; bit-exact spin-statistics grading σ_π(F)=(−1)^{2b}F; causal classifier matches symbol signs 9/9 with b=1 top weight ON the cone. Falsifier `tests/test_wh7_b3_boost.py`; Paper 45 Q1 paragraph. See `debug/sprint_wh7_b3_boost_probe_memo.md`.
- **Sprint B1 joint product-carrier convergence (2026-06-10, v3.112.0):** Paper 45 gains `prop:product_action_seminorm` — truncated S³×S¹_T system converges in vS state-space GH at additive rate γ_n+γ_K under the joint translation seminorm; 6-check panel all PASS (pure-factor exactness 0/10⁻¹⁶, Leibniz envelope, kernel condition, additive smoothing, N_t=1 reduction bit-exact); Q1 sharpened to signature-only (boost/modular-flow seminorm = named candidate); claims register row 21; P45 23pp GATE: PASS. See `debug/sprint_wh7_b1_joint_memo.md`.
- **Sprint WH7 Toeplitz probe Step 1 (2026-06-10, v3.111.0):** POSITIVE-REBUILD — translation seminorm = Lipschitz exactly on Toeplitz temporal modes (10⁻¹⁴); P45 invisibility reproduced as momentum-diagonal-algebra artifact; visibility survives de-compactification at fixed bandwidth. Compact-time Lorentzian wing rebuilt at operator level; falsifier `tests/test_wh7_toeplitz_temporal.py`. See `debug/sprint_wh7_toeplitz_probe_memo.md`.
- **Sprint CLAUDE.md compaction round 2 (2026-06-10, v3.110.1):** 320 KB → 133 KB (−59%); §1.7/§2/§6/§7/§12 + §1.5 chronicle moved verbatim to 5 docs/ files; §13.11 rule 9 (replace-don't-append, PI-authorized); repo health gate added to /release precondition 8. See CHANGELOG v3.110.1.
- **Sprint repo-hygiene + WH7 registration (2026-06-10, v3.110.0):** debug/ swept 1,824→458 top-level files into `debug/archive/<arc>/` (manifest + READMEs; test-referenced files pinned); `papers/INDEX.md` status map + README start-here block; 147 untracked LaTeX artifacts purged. WH7 (time-discreteness is observer-compactification) registered in §1.7 per PI direction. See `debug/sprint_repo_hygiene_memo.md`.
- **Sprint P38-G1/G2 closure — theorem UNCONDITIONAL (2026-06-10, v3.109.0):** Translation-seminorm reframing dissolves both named gaps: G2 kernel condition holds on the TRUTHFUL CH substrate (Schur + per-band injectivity, verified n_max=2..5); G1 dual reach closed by exact-fit spinor lifted state (CH shells ARE the V_j⊗V_{j±1/2} window blocks; defects = Fejér smoothing at γ_n). Paper 38 unconditional at rate (4/π+o(1))log n/n in vS state-space framework; cascade to P45/P32/register/README/field-guide/N1; falsifier `tests/test_p38_action_seminorm.py`; WH1 restored to PROVEN. Outreach blocker cleared. See `debug/sprint_p38_g1g2_phaseA_memo.md`.
- **Sprint repositioning — Phase 3 COMPLETE (2026-06-10, v3.108.0):** PI-directed §1 viability rewording (research instrument, not production tool); README badge + matched-qubit caveats; Paper 14 abstract caveat (P20 already compliant); corpus-wide α-mention audit CLEAN (all ~15 papers compliant, zero fixes). See CHANGELOG v3.108.0.
- **Sprint accessibility layer — Phase 2 COMPLETE (2026-06-10, v3.107.0):** N1+N2 front-door notes (docs/outreach/, send-gated), claims register (docs/claims_register.md, 20 rows), vocabulary translation (33 rows), standalone audit + rewords applied (P38/P45), README/field-guide/synthesis/.zenodo.json stale-claim fixes. Bonus: Paper 28 Thm 3 proof factor-2 error caught drafting N2, fixed, 40-digit verified. See CHANGELOG v3.107.0 + docs/corpus_accessibility_plan.md status block.
- **Sprint de-versioning + subagent budget policy (2026-06-10, v3.106.1):** PI directive: papers are a single source of truth corrected in place (git/Zenodo = version record) — all 8 papers stripped of v1/v2/erratum scaffolding, content unchanged, one History remark each in P38/P45; all GATE: PASS. New standing rule `memory/feedback_subagent_token_budget.md` (model tiering, batching, paste-don't-point, hard caps, scripted gates) + `debug/compile_3pass.sh`; sibling pass ran on sonnet at 141.7k tokens vs 366k same-scope yesterday.
- **Sprint P45-hardening + corpus descope (2026-06-09, v3.106.0):** Adversarial Phase-1 pass (accessibility plan) FALSIFIED Paper 45's main theorem — K⁺ seminorm ≡ 0 bit-exact, "Latrémolière Thm 5.5" nonexistent, L2 mass/symbol conflation. Option C executed: P45 corrected in place to the annihilation theorem, P38 to a conditional vS restatement (2 named gaps), Status notes in P40/46/47/48/49 + P32; falsifier frozen `tests/test_p45_kplus_degeneracy.py`; outreach blocked pending P38 repair. See `debug/sprint_p45_hardening_phase1_memo.md` + `docs/corpus_accessibility_plan.md`.
- **Sprint pre-outreach corpus hygiene (2026-06-08, v3.104.0):** Coordinated multi-paper cleanup pass before Brown/Kleinschmidt outreach. Five parallel audits (Papers 55, 56, 32 §VIII, LiH regression, corpus citation sweep) followed by sixteen-task execution: H1 Yukawa promoted to `\begin{theorem}` block (Paper 32 §VIII tally now 8/8); `thm:no_single_mechanism_K` proof rewritten to separate structural argument from empirical confirmation; 18 bibitem corrections across 13 papers (4 wrong arXiv IDs in P55; Perez-Sanchez titles in P38 + P56; Mondino-Sämann title in P44/45/46/47; Nieuviarts titles in P42/43/44; duplicate Paper 55 bibitem in P56; missing karamata/korevaar/tenenbaum bibitems in P55); lineage citations added to P18 + P29; field guide refreshed with C-arc closure paragraph + Papers 45/46/56 in readers' map + Deligne-Milne converse / non-commutative MS in open frontiers; README updated to v3.104.0 + new Math.OA / NCG / Periods Arc section + 17 new Paper Series entries; LiH regression memory note retired (confirmed RESOLVED v3.56.0, Paper 17 5.3% reproduces at 2.82%). 18/18 topological proofs pass; all edited papers compile three-pass clean. See `debug/cleanup_*_memo.md` (5 audit memos).
- **Sprint spin-structure moduli + Dirac-index obstruction (2026-06-08, v3.105.0):** Layer-1 flat-$\mathbb{Z}_2$ spin structures $= H^1(G;\mathbb{F}_2)$. Scalar Hopf $m\to-m$ gives 2 symmetric configs at $n_{\max}=3$; Dirac $m_j\to-m_j$ is a fixed-point-free Kramers automorphism, but chirality $\chi\to-\chi$ is not even a node bijection ($2\ell+2$ vs $2\ell$ per block, Dirac index $n_{\max}(n_{\max}+1)$) — substrate root of the relativistic-$\mathbb{Z}_2$ negatives. Paper 29 `obs:chirality_obstruction`; 11 tests. See `debug/spin_structure_moduli_memo.md`.
- **Sprint C-arc closure — E6 + D5/D6 + chemistry analog (2026-06-08, v3.103.0):** Two new theorems + one named multi-month target in Paper 32 §VIII. **E6** `thm:no_single_mechanism_K`: 3 spectral homes in 2 Mellin sub-rings; no morphism generates $K = \pi(B+F-\Delta)$; 12 mechanisms eliminated; combination rule preserved as conjectural per §13.5. **D5/D6** `thm:cutoff_function_external` + corollary: cutoff function and its φ-moments are external test-function data; Wald forces relations not values; CC fine-tuning $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ formalised as external moment selection. **Chemistry-analog** `rem:chemistry_eta_analog`: named multi-month NCG-research target (FrozenCore chirality-grading analog). **Terminal state remark** `rem:c_arc_terminal_state`: **eight theorem-grade non-selection results across four sectors**, C-arc closed. See `debug/sprint_c_arc_closure_e6_d5d6_chemistry_memo.md`.
- **Sprint G1/G2/G5 — Spatial-composition wall theorem (2026-06-08, v3.102.0):** Paper 32 §VIII gains `thm:spatial_composition_radial_wall` + `cor:spatial_composition_wall` + `rem:multi_focal_wall_fully_characterized`. Tensor-product spectral triple produces forced angular structure (Paper 54 Thm 3) but unforced radial coupling (Fock-projection conformal-factor non-commutativity). Closes **three** catalogue entries (G1, G2, G5) — the spatial-composition sub-sector of multi-focal-composition wall. **Multi-focal-composition wall now fully theorem-bound** under two structurally distinct theorems (renormalization + spatial-composition). **Six theorem-grade non-selection results now in corpus**. See `debug/sprint_g1_g2_g5_spatial_composition_memo.md`.
- **Sprint E7/E8 — Single-cutoff spectral action theorem (2026-06-08, v3.101.0):** Paper 32 §VIII gains `def:multi_cutoff` + `thm:single_cutoff_spectral_action` + `cor:multi_loop_renormalization_wall` + `rem:single_cutoff_scope`. Theorem: CC spectral action axiom is single-cutoff; no morphism in $\mathcal{A}$ produces multi-cutoff structure. Closes **four** catalogue entries (E7, E8, G3, G4) — the multi-loop QED renormalization sub-sector of the multi-focal-composition wall. **Five theorem-grade non-selection results now in corpus**: Forced-Count, H1 Yukawa, $N_{\mathrm{gen}}$, KO-dim, single-cutoff spectral action. See `debug/sprint_e7_e8_single_cutoff_memo.md`.
- **Sprint F3 — Inner KO-dim non-selection theorem (2026-06-08, v3.100.0):** Companion to C3. Paper 32 §VIII gains `thm:ko_dim_non_selection` + `rem:full_inner_factor_boundary`. Argument is one-line composition: packing is kinematic (Paper 0 §VII.B), KO-dim is real-structure data, packing produces no real-structure → $\mathcal{A}$ cannot autonomously select KO-dim. **Four theorem-grade non-selection results now characterise the full inner-factor structural-skeleton boundary at the canonical-rep level**: Forced-Count, H1 Yukawa, $N_{\mathrm{gen}}$, KO-dim. Paper 57 §3.1 + §6.3 closure updated to four. See `debug/sprint_f3_ko_dim_non_selection_memo.md`.
- **Sprint C3 — $N_{\mathrm{gen}}$ non-selection theorem (2026-06-08, v3.99.0):** Theorem-grade upgrade of Direction 2 + Read 2 NO-GO scopings. Paper 32 §VIII gains `thm:n_gen_non_selection` + `rem:n_gen_scope` (conditional on canonical CCM rep). Three theorem-grade non-selection results now in corpus: Forced-Count, H1 Yukawa, $N_{\mathrm{gen}}$. Paper 57 §3.1 + new §6.3 closure updated. See `debug/sprint_c3_n_gen_non_selection_memo.md`.
- **Sprint C2 principle hunt — P5 packing-reachability (2026-06-08, v3.98.0):** Formal test of five candidate discriminators against the v3.97.0 catalogue. P5 (packing-reachability) hits **98.3% accuracy** with a single ambiguous misclassification (I3 Higgs direction, already flagged as conditional). Two-family structure preserved as failure-mode decomposition beneath P5. Paper 57 §5.5+§5.6+§6.1+§6.2 updated. See `debug/sprint_c2_principle_hunt_memo.md`.

> Older sprint index (v2.x–v3.96.0), the long-form arc chronicles, and the RH sprint records moved verbatim to `docs/development_frontier_archive.md` (2026-06-10 compaction). CHANGELOG.md remains the canonical chronicle going forward.

**Best results by system type:**

| System | Result | Method | Paper |
|:-------|:-------|:-------|:-----:|
| He (atom) | 0.019% | 2D variational S⁵ + cusp, l_max=7 | 13 |
| He (graph-native CI) | 0.20% | Zero-parameter, exact algebraic integrals, n_max=9 | 13 |
| H⁻ | Bound, over-binds 21% | Graph-native CI, Z_c≈1.84 boundary | 13 |
| PsH | 4.1% | Level 3, sign-flipped charge | 13 |
| H₂⁺ | 0.0002% | Spectral Laguerre | 11 |
| H₂ | 96.0% D_e | Level 4, l_max=6, 61 channels | 15 |
| LiH | R_eq 5.3% | Composed, l-dependent PK, l_max=2 | 17 |
| BeH₂ | R_eq 11.7% | Composed, full 1-RDM exchange | 17 |
| H₂O | R_eq 26% | Composed, 5-block, zero parameters | 17 |
| LiH (4N) | R_eq 63.5% | Full 4e mol-frame, PK-free | 17 |
| Composed Pauli | O(Q^2.5) | 51x-1712x vs Gaussian, 28 molecules | 14 |
| Atomic Pauli | O(Q^3.15) | 1.3x-8.1x vs cc-pVDZ/cc-pVTZ | 14 |
**Key structural results (details in papers and CHANGELOG.md):**
- **κ = −1/16 derivation (v2.26.1):** Derivable from Fock projection, not fitted. Paper 18 reclassified κ to "conformal." See `debug/probe_kappa_sprint_memo.md`.
- **α structural decomposition (Phases 4B-4I, April 2026):** B=42 (Casimir), F=π²/6 (Fock Dirichlet at d_max), Δ=1/40 (Dirac degeneracy g_3). Three independent spectral homes; combination rule K=π(B+F−Δ) remains conjectural. 12 mechanisms eliminated. See Paper 2, CHANGELOG.
- **Spectral-action supertrace (v2.26.1):** SD cancellation theorem, Δ⁻¹=40 from Euler-Maclaurin, (−) sign = (-1)^F grading. Two-term exactness on S³ Dirac. See `debug/st_supertrace_sprint_memo.md`.
- **Nuclear systems (Paper 23):** Deuteron 16q/592 Pauli; He-4 16q/712 Pauli; composed nuclear-electronic deuterium 26q/614 Pauli. Fock rigidity theorem (S³ unique to −Z/r).
- **Angular sparsity theorem (Paper 22):** ERI density depends only on l_max, not V(r). Universal across potentials.
- **Bargmann-Segal lattice (Paper 24):** HO on S⁵ Hardy space, π-free, HO rigidity theorem. Coulomb/HO asymmetry = 4 layers.
- **S⁵ gauge extension (v2.26.1):** U(1) Wilson transfers; SU(3) NOT natural; CP² quotient fails. See `debug/s5_gauge_structure_memo.md`.

**Classical solver status: INVESTIGATION COMPLETE (v2.0.24).** 30+ tracks exhausted all solver × PK × basis combinations. Structural ceilings characterized. Composed at l_max=2 is the production operating point. See CHANGELOG.md for full track history.

- **Entropy/entanglement arc (Papers 26-27, v2.9.2):** Universal S_B(w̃_B/δ_B)^γ scaling with γ_∞≈1.96; HO zero-entropy rigidity exact; bond blocks break single-center curve. See Paper 27, CHANGELOG EP-2 series.
- **Cusp re-diagnosis (v2.9.2):** He 0.20% CI floor is graph-validity-boundary artifact (Z_c≈1.84), NOT cusp. TC dead at every n_max. See `debug/cusp{1,2,3}_*.py`.
- **Graph validity boundary (v2.9.0):** Variational bound violated below Z_c≈1.84; mechanism is Z-independent κ vs Z²-scaled diagonal.
- **111 Pauli count derivation (v2.9.0):** 111 = 55 direct + 56 exchange per s/p block. Universal coefficient 11.11×Q across 28 molecules.

**Quantum computing status: ACTIVE FRONTIER (Paper 14).** O(Q^2.5) composed Pauli scaling, 51x-1712x vs Gaussian, 28 molecules, ecosystem export (OpenFermion/Qiskit/PennyLane). PK classical partitioning gives 78x 1-norm reduction for H₂O. Market test: LiH 0.97× 1-norm vs STO-3G with 2.7× fewer Pauli. 40 molecules total in library. See Paper 14, Paper 20, CHANGELOG Tracks AW-CA.

**Balanced coupled (Track CD, v2.0.39+):** Cross-center V_ne via multipole; LiH 878 Pauli, 1.8% energy, 7.0% R_eq; only bound 4e config. n_max=3: 0.20% energy. BeH₂ 2,652 Pauli; H₂O 5,798 Pauli. See Paper 19, CHANGELOG.

**Precision atomic spectra (Track DI, v2.6.0):** 2D variational breaks adiabatic floor (0.022% raw, 0.004% cusp-corrected at l_max=4). Graph-native CI 0.20% at n_max=9 with exact algebraic integrals. FCI basis-invariant. See Paper 13, CHANGELOG.

**RH sprint (v2.20.0-v2.23.0, April 2026):** GeoVac Hopf graphs strictly Ramanujan (Paper 29). Six sprint-scale follow-ups: Alon-Boppana crossing at V~30-60 (finite-size statement), Hopf-U(1) block hypothesis validated, spectral χ₋₄ closed form D_even−D_odd = 2^(s−1)(β(s)−β(s−2)), GUE-like spectral-zero stats (CV≈0.35-0.40), no functional equation (48 OoM gap), SU(2) Wilson on S³ (Paper 30). Direct spectral-to-classical-RH bridge CLOSED by three independent walls (WH6). See Papers 28-30, CHANGELOG.




**Sprint 5 S_min erratum patch (v2.23.1, April 17, 2026):** the Sprint 4 RH-P side flag was independently verified by three methods (direct sum + 5-term explicit tail; mpmath.nsum Levin u-transform; Euler-Maclaurin 26-term asymptotic) all agreeing at ≥ 80 digits on the true value **S_min = 2.47993693803422255441357950082938214468792578661728845837879872655955...** Double-diagnosis: (a) Paper 28's published 2.47953699802733387 was wrong at the 4th decimal due to an erroneous tail formula in `debug/smin_identification.py` lines 45-48 which assumed T(k)~2/a² (the correct asymptotic is T(k)~2/a); (b) RH-P's 2.47993693803422255447852790477 was wrong at the 20th digit due to incorrect c_6..c_11 coefficients in `debug/compute_smin_chi_neg4.py` T_SQUARED_COEFFS. **Fixes applied**: Paper 28 §IV text updated with the corrected value at 25 digits; `debug/smin_identification.py` compute_S_min_with_tail now delegates to mpmath.nsum Levin (2-line fix, ~2 s runtime, verified to 15+ digit match); `debug/compute_smin_chi_neg4.py` T_SQUARED_COEFFS retained with documented 20-th-digit limitation (does not affect the RH-P PSLQ negative result, which was Q-linear-independence against a finite basis insensitive to 6.5e-20 numerical shifts). Irreducibility of S_min against the 47-element extended basis **unchanged** — 15 PSLQ failures stand. Key files: `debug/compute_smin_verification.py`, `debug/smin_verification_memo.md`, `debug/data/smin_verification.json`.

**QED arc (Papers 28, 33, 36, April-May 2026):**
- QED on S³: 5 theorems (T9, χ₋₄, self-energy zero, ζ(3) complementarity, product survival), S_min irreducible at 200 dps, D₅/D₆ Sommerfeld via PSLQ. See Paper 28.
- Graph-native QED (GN-1..GN-7): F₂=5√2/3 π-free; pendant-edge theorem Σ(GS)=2(n-1)/n; 1/8→4/8→7/8→8/8 selection rule recovery. See Paper 33.
- Vector-photon QED: 7/8 scalar, 8/8 Dirac via spinor phase constraint. Calibration 1/(4π) per loop. See Paper 33, CHANGELOG.
- Bound-state QED / Paper 36: Lamb shift −0.534% at one loop; LS-8a two-loop wall (counterterms not autonomous). See Paper 36.
- Sprint HF (2026-05-07): H 21cm at +18 ppm; multi-focal-composition wall crystallized (5 observables, same structural reason).
- Sprint MH (2026-05-08): μH Lamb −0.10% (full Uehling kernel); BF HFS +2 ppm; rest-mass projection verified.
- Precision catalogue (2026-05-08/09): 9 systems sub-percent across mass-hierarchy × nuclear-spin × multi-focal axes. Paper 34 §V. See CHANGELOG.

**Spectral triple / WH1 arc (Papers 32, 38-50, April-May 2026):**
- WH1 PROVEN (2026-05-06): GH-convergence theorem, five-lemma proof, Λ ≤ C₃·γ_{n_max} → 0, rate 4/π. Paper 38.
- Paper 39: Tensor-product propinquity. Paper 40: Unified GH on all compact Lie groups, 4/π universal. Master theorem subsumes 38/39/40.
- Sprint H1: AC extension POSITIVE-THIN (Higgs admitted, Yukawa not selected). Sprint G3: NEGATIVE (γ_GV ≠ γ_F).
- Sprint TS (2026-05-04): Case-exhaustion theorem (Paper 32 §VIII); master Mellin engine M1/M2/M3. See `debug/track_ts_*_memo.md`.
- TX-A/TX-B (May 2026): Paper 34 three-axis axiomatization; Paper 35 Prediction 1 graduated (208/208). See `debug/tx_{a,b}_*_memo.md`.
- ST-SU3: SU(3) Wilson on S⁵ Bargmann — gauge YES, matter NO; universal 1/(4N_c). See Paper 30 §7.7.
- WH1 R1-R3.5 (2026-05-03/04): Connes-vS operator-system alignment, prop=2, Avery-Wen-Avery, full Dirac L1'. See Paper 32 §III.

**Lorentzian arc (Papers 42-49, May 2026):**
- Sprint L0 (2026-05-16): 28-projection transfer audit (17/4/5/2 = free/Wick/Euclidean/mixed). Paper 31 §8, Paper 34 §V.E.
- Sprint L1 (2026-05-16): σ_{2π}(O)=O bit-exact at n_max=2..5; four-witness theorem at operator-system level. Paper 42.
- Sprint L2 (2026-05-16/17): Krein (3,1) + Lorentzian Dirac + axiom audit + modular Hamiltonian. Paper 43. H_local ≠ D_W signature-independent.
- Pythagorean HS-orthogonality (2026-05-17/23): ⟨H_local,D_W^L⟩=0 bit-exact, closed form, 1/π² M1 signature. Paper 43 §10.2.
- Sprint L3a-1 + Paper 44 (2026-05-17): Lorentzian operator system, prop=2/∞ envelope-dependent. Paper 44.
- Sprint L3b-2 + Paper 45 (2026-05-18): K⁺-weak Lorentzian propinquity — first in literature. Paper 45.
- Sprint L3b-2a-d + Paper 46 (2026-05-22): Strong-form closed; Λ^strong=Λ^P45 bit-exact. Paper 46.
- Sprint L3c + Paper 47 (2026-05-23): G2 norm-resolvent; three-carrier identification. Paper 47.
- Phase A + Paper 48 (2026-05-24): Krein-MS bridge; 7 newly accessible theorems. Paper 48.
- Q1' + Paper 49 (2026-05-24/25): OSLPLS strong-form bridge; twin paradox. Paper 49.
- **Sprint math.OA-arc-closure (2026-05-31):** G2-metric, G3, Q2, Q2' all CLOSED. See `debug/sprint_mathoa_arc_closure_memo.md`.

**Gauge arc (Papers 25, 30, 41, April-May 2026):**
- Paper 30 (v2.25.0): SU(2) Wilson on S³; maximal-torus → Paper 25 U(1); L₁ = kinetic term. See Paper 30.
- Paper 41 / XCWG arc (2026-05-16): Rule B Wilson U(1), seven-witness compatibility with 3D compact U(1). See Paper 41.
- Paper 18 consolidation (v2.24.0): Light-touch v1.0→v1.1; meta-pattern cross-reference.

**Dirac-on-S³ arc (v2.12.0-v2.19.4, April 2026):**
- Tier 2 (T0-T6): Spin-ful composed qubits (LiH/BeH/CaH relativistic). Spinor certificate R_sp = ℚ(α²)[γ]. Papers 14, 22 updated.
- Tier 3 (T7-T9): Darwin+MV α⁴ ladder; T9 squared Dirac ζ_{D²}(s) = π^{even} only theorem. Paper 18 4th cell filled.
- Breit SS/SOO (v2.18.0): Retarded radial, Drake combining at operator level. Sprint CP: Li +9%, Be +3%.
- One-loop QED VP (v2.18.2): a₀=a₁=√π, β(α)=2α²/(3π) from S³ spectral data.
- Compactness thesis (v2.18.1): Paper 18 Peter-Weyl rewrite. Kramers-Pasternak direct integration.

**Chemistry arc (v2.0.30-v2.3.0, Tracks BG-CX):**
- General composed builder (v2.0.30): `build_composed_hamiltonian(spec)` + atomic classifier Z=1-10.
- Second-row (v2.1.0): [Ne] frozen core, Z=11-18, Q^2.50 scaling. Third-row (v2.2.0): [Ar]/[Ar]3d¹⁰, Z=19-36.
- Multi-center (v2.3.0): 8 new molecules (LiF, CO, N₂, F₂, NaCl, CH₂O, C₂H₂, C₂H₆). N_Pauli = 11.11×Q universal.
- TM hydrides (v2.8.0): Z=21-30, d-orbital blocks, Pauli/Q=9.23. Heavy-atom (v2.12.0): SrH/BaH [Kr]/[Xe] cores.
- W1c arc (2026-05-23): F1→F2→F3 closed W1d (cross-block h1); W1e opened (inner-region overattraction). F4-F6 negatives.
- Balanced second-row FCI (v2.19.4): NaH/MgH₂ overattract — frozen-core limitation. See CHANGELOG.

**Precision catalogue / multi-focal arc (2026-05-07 through 2026-05-18):**
- Multi-focal Phase C closures: W1a cross-register V_eN, W1b Zemach operator, W1c screened V_ne, W2b tensor propinquity. See CHANGELOG.
- Roothaan autopsies (2026-05-09): 8-track sprint; 6 §V.C + 4 §V.D entries in Paper 34. See `debug/*_autopsy_*_memo.md`.
- Hylleraas r₁₂ (2026-05-09): He 1¹S at 0.0006%; Track 5 closure He 2¹P→1¹S f-value −2.02%. See Paper 34.
- Z>20 cliff diagnostic: CR67 single-zeta non-faithful; BBB93 GO with eyes open. See `debug/z_cliff_*_memo.md`.
- AdS/CFT adjacent (2026-05-25/26): Paper 50 — bit-exact F-theorem match on S³/S⁵; M2/M3 orthogonal decomposition.

**Gravity arc (Papers 51, 53, v3.4.0-v3.28.0, May 2026):**
- G1-G8: spectral action on S³ (two-term exact), thermal product, BH entropy, Newton constant, graviton diagnostic. Paper 51.
- G6-FP Fierz-Pauli: J-blindness theorem closes the graviton at spectral-action level; FP requires metric identification (propinquity). Paper 51 §FP.
- G4-3..G4-6: discrete warped substrate program (multi-month). Spinor tip −1/12 bit-exact. L6 replica-weight-harmless CLOSED.
- Paper 53: disk-with-cone propinquity (first manifold-with-boundary carrier). Bochner-Riesz plane reconstruction.
- Möbius α>1: RETIRED as finite-a substrate artifact (B4, v3.24.0). Continuum = plain Sommerfeld-Cheeger.
- Confinement reframing (v3.26.0-v3.28.0): ARCHIVED as organizing reading (2026-05-31). See §1.7.

**Scope:** See `SCOPE_BOUNDARY.md` for supported atoms/molecules (40 molecules total, Z=1-56 via frozen cores).

**Architecture locked:** The LCAO/graph-concatenation approach (v0.9.x series) is superseded. All molecular work uses natural geometry (Papers 11, 13, 15, 17).

---


## 3. Approaches That Failed

Critical institutional memory. Do not re-derive these dead ends. Full details in CHANGELOG.md. Sub-agents must check CHANGELOG.md before attempting any approach in these categories.

| Category | Count | Key Lesson |
|:---------|:-----:|:-----------|
| LCAO / single-S³ molecular encoding | 3 | Graph Laplacian kinetic energy is R-independent; need natural geometry where separation occurs |
| PK modifications (projector, spectral, self-consistent, R-dependent, l-dependent on 2D) | 6 | PK provides coordinate-space exclusion that angular projectors categorically cannot replicate; PK is the irreducible cost of composed-geometry factorization |
| Cusp treatments (alpha-only, graph absorption, θ₁₂-adapted basis) | 3 | Cusp is 2D in (α, θ₁₂); 1/r₁₂ is an embedding exchange constant; Schwartz extrapolation is the correct approach |
| Inter-group antisymmetry (node exclusion, fiber bundle, spectral) | 3 | Antisymmetry requires shared coordinate system; composed geometry factorizes into incompatible coordinates |
| Full N-electron radial solvers (adiabatic, coupled-channel, 2D variational) | 3 | Adiabatic overcounts D_e; coupled-channel numerically unstable; 2D gives unbound D_e; angular basis is the bottleneck |
| Polyatomic coupling (Z_eff partition, classical repulsion, lone pair at Z_eff>4) | 3 | Orbital-level exchange coupling needed; lone pair Slater integrals unphysical at high Z_eff |
| Geometric elevation (blow-up, Lie algebra, S³×S³) | 1 | min/max boundary is physical; per-ρ diagonalization structurally irreducible |
| Diagnostic arcs (eigenchannel rotation, spheroidal compression, enhanced Z_eff, midpoint origin) | 4 | Various; see CHANGELOG.md |
| l_max convergence via 2D solver (variational_2d in composed pipeline) | 1 | l_max divergence is structural to PK/composed architecture, NOT from adiabatic approximation. 2D solver produces identical R_eq drift (+0.15-0.22 bohr/l_max) as adiabatic. Single-point energy also diverges, violating variational bound. l_max=2 is optimal. |
| Sturmian CI (Coulomb, generalized) | 2 | Coulomb improves atomic FCI but Lowdin 1-norm inflated 2.8-4.5x; generalized loses within-config flexibility at larger basis; neither bypasses PK ceiling |
| TC Jastrow in adiabatic hyperspherical solver | 1 | TC replaces multiplicative V_ee (O(R) in angular eigenvalue) with first-derivative G_ang (O(1)); adiabatic separation requires V_ee as multiplicative operator; ~46% error; TC needs direct variational/FCI framework, not adiabatic. |
| TC in second quantization (composed qubit pipeline) | 1 | TC qubit Hamiltonians plateau at ~3.4% across n_max=2..5; standard FCI converges 5.3%→2.0%. Original benchmark used qubit-space diag (false positive). Cusp is energy-evaluation problem, not wavefunction problem. See CHANGELOG v2.9.0/v2.9.2 and `feedback_tc_correction.md`. |
| TC 2D variational cusp correction (basis doubling) | 3 | Three attempts failed (R sin(α) basis worsens, sin(α) negligible, [H,r₁₂]=0 identically). TC similarity transform succeeds non-monotonically at γ=0.11 — finite-basis correction, not true cusp removal. |
| TC angular gradient for l>0 orbitals | 1 | Adds 2.66× Pauli terms for 0.01 pp accuracy. Cost/benefit >100×. Radial-only is optimal. |
| Coupled composition (cross-block ERIs replacing PK) | 1 | Cross-block ERIs add 2.56× Pauli, 2.30× 1-norm; 29% FCI error. PK is NOT redundant in 2Q for composed basis (single-center framework lacks two-center h1). See CHANGELOG Track CB v2.0.37. |
| Single-center nested molecular LiH | 1 | R_eq 33.7% error; single Z can't represent both core and bond length scales. Track DF Sprint 4. |
| Two-center charge-center nested LiH | 1 | 48.2% energy error; truncated basis cannot represent 1s² core off-center. Track DF Sprint 4B. |
| Heterogeneous nested (per-pair Z_eff) | 1 | Löwdin orthogonalization destroys Gaunt sparsity (1711 vs 120 Pauli, 14× inflation). Track DF Sprint 5. |
| Fock energy-shell self-consistency for He | 1 | k²=−2E over-constrains 2-electron problem; SC 12.8%/5.1% vs variational 1.9%/1.6%. Single parameter insufficient for two electrons. Track DI Sprint 2. |
| SM-running origin for Δ = 1/40 (Paper 2 alpha) | 6 | Six tracks ruled out SM-running origin for Δ=1/40. Positive byproduct: Δ⁻¹ = g_3^Dirac(S³) = 40 exactly. See CHANGELOG Phase 4H. |
| Schwartz-tail hot-node patch on He Z=2 graph-native CI (CUSP-2) | 1 | Schwartz tail correction worsens accuracy; 0.20% floor is small-Z graph-validity-boundary artifact (Z_c≈1.84), not cusp. Z=10 sign flip confirms. See `cusp_rediagnosis.md`. |
| Energy graph for V_ee on S³ (Paper-12 analog search) | 1 | Pair-state graph dense (47%); cross-shell denominators don't close. Wavefunction graph diagonalizes 92-94% on its own; cusp concentrated on (1s,1s). See CHANGELOG v2.9.1. |
| Darwin+MV for He/Li/Be 2p-doublet improvement | 1 | Both 2p states share l=1; Darwin=0 for l≥1; MV cancels in splitting. Residual 66-211% errors trace to multi-electron SS/SOO. Tier 3 Track T8. |
| Balanced + frozen-core PES for second-row molecules | 2 | NaH/MgH₂ overattract monotonically at n_max=2; frozen [Ne] hides core screening from cross-center V_ne. Limited to first-row for PES. Sprint 7a v2.19.4. |
| Single-constant graph-to-continuum QED projection (C×F₂→α/(2π)) | 1 | C×F₂ grows with n_max, doesn't converge to α/(2π). Different diagrams need different projection scaling; topology-dependent. C_VP/C_SE ≈ 3/29 sub-percent stable across n_max=3,4,5. See CHANGELOG Sprint GN-QED v2.26.1. |
| σ-vertex and direction-resolved vector QED on Fock graph (VQ-1..VQ-5) | 5 | Five tracks ruled out σ-vertex approaches to vector photon selection rules. The 4 missing rules require photon to carry (L,M_L) quantum numbers, not just spin matrices at the vertex. See `vector_qed_7of8.md`. |
| Co-exact mode q-labeling for SO(4) channel count / Ward / charge conjugation | 1 | Co-exact eigenvectors have right eigenvalues but wrong support (nearest-neighbor only); triangle inequality fails at high q. 5/8 selection rules maximal without vector photon bundle. See CHANGELOG. |
| Dirac-sector lift of Paper 2 α combination rule ingredients B, F, Δ | 3 | Three obstructions: B doesn't lift (zero mode forces (m−1)), F doesn't lift (Apéry Q-linear independence), Hopf-equivariant decomp doesn't produce B/F. ζ(3) emerges natively in Dirac sector. See CHANGELOG Phase 4I Tier 1. |
| Multi-focal spatial composition (Sprint HF, 2026-05-07) | 3 | HF-3/4/5 (recoil, Zemach, multi-loop a_e) all hit one wall: framework couples discrete labels cleanly but has no native composition theorem for multiple Fock-style projections. See `multi_focal_wall_pattern.md`. |
| Phillips-Kleinman cross-center barrier for second-row chemistry binding (2026-05-08) | 1 | NaH PES still descending after W1c+PK (0.357→0.305 Ha, 14.6% only). Wall structurally deeper than orthogonality. See `chemistry_arc_paused_w1c_residual.md`. |
| Screened-Schrödinger valence basis (h1 diagonal only, Track 3, 2026-05-09) | 1 | SV diagonal correction is R-INDEPENDENT; cannot affect R-dependent descent depth. Wall lives in cross-V_ne integration shape, not eigenvalues. See `feedback_diagnostic_before_engineering.md`. |
| W3 spectral-zeta calibration-data identification (Sprint W3, 2026-05-08) | 3 | CKM Wolfenstein candidate falsified by three follow-up tracks (PMNS, lepton mass, mechanical 21,448-form basis). Selection-bias fully attributable. WH7 NOT promoted. See `w3_spectral_zeta_candidate.md`. |
| Multi-focal Path C5 saturated basis (numerical-luck f closure, 2026-05-09) | 1 | f=0.278 was cancellation artifact at cond(S)=10^10. Production result f=0.286 (+3.4%) at cond(S)=400. Diagnostic-before-engineering inside sprints. See `multifocal_phase_d_production.md`. |
| Heuristic two-zeta screening for [Xe] core (CR67 + BBB93 Kr ratios, 2026-05-09) | 1 | Cs HFS goes from −47% to −90% (wrong direction). CR67 single-zeta fits non-faithful for outer shells. Needs full BBB93/KTT or self-consistent HF. See `z_scaling_cliff_cr67_fits.md`. |
| All-positive-coefficient single-zeta hydrogenic basis lacks radial nodes (Track 2 diagnostic, 2026-05-09) | 1 | Single-zeta hydrogenic R_nl has no radial nodes; real RHF outer-shells need negative-coefficient orthogonalization. Cliff onset K Z=19. BBB93 GO with eyes open. See `z_scaling_cliff_cr67_fits.md`. |
| Multi-zeta physical Na valence basis substitution (Sprint α-Multi-zeta + α-PES, 2026-05-23) | 1 | Step 1 algebraic differential −0.135 Ha sign-consistent; Step 2 FCI shift BIT-ZERO (Na 3s sits at i=5 in eigenspectrum, unoccupied). Layer-3 FCI invisibility. Architecture retained. See `sprint_modular_alpha_arc.md`. |
| Three-bucket M-Z partition (basis-closable cross-shift) FALSIFIED at NaH max_n=3 (2026-05-23) | 1 | 1/3 predictions pass; dominant NO IS bonding combination at max_n=3 but energetically unfavored. Framework can BUILD the right orbital but cannot energetically PREFER it. See `debug/sprint_f1_maxn3_predictions_test_memo.md`. |
| Kernel-shape substitution as W1c-residual closure (Sprint F2, 2026-05-23) | 1 | Multipole expansion bit-faithful (max diff 2×10⁻⁵ Ha, 6-10 OoM below wall). KERNEL-NOT-IT. Substantive finding: cross-block h1 architecturally absent — W1d named. See `debug/sprint_f2_cross_vne_kernel_memo.md`. |
| Single-particle Pauli orthogonality (rank-1 PK on bonding orbital) as W1e closure (Sprint F4, 2026-05-23) | 1 | Rank-1 PK saturates at 43% closure ceiling regardless of barrier magnitude. W1e is multi-determinant FCI correlation, not single-particle Pauli. See `debug/sprint_f4_bonding_pk_memo.md`. |
| Mean-field core-bonding J-K (explicit-core Hartree) as W1e closure (Sprint F5, 2026-05-23) | 1 | J-K = +1.12 Ha = 25.7% of wall (right sign, insufficient magnitude). W1e is deeper than Hartree-level core-bonding interaction. See `debug/sprint_f5_explicit_core_memo.md`. |
| Basis-level Schmidt orthogonalization of H 1s against Na [Ne] core as W1e closure (Day-1 diagnostic, 2026-05-23) | 1 | Projection mass Σ|S_c|² = 1.0×10⁻³ (0.1% of H 1s in [Ne] space); diagonal differential +8.6 mHa = 0.2% of wall. H atom far outside [Ne] amplitude at R=3.5. See `debug/sprint_w1e_schmidt_diagnostic_memo.md`. |
| [Ne] core correlation as W1e closure mechanism (Day-1 literature estimate, 2026-05-23) | 1 | Binding-relevant differential ~0.001-0.003 Ha across 3 literature anchors; 14-1000× below DEFER threshold. Cannot close 3.85 Ha wall. W1e is 6th instance of multi-focal wall pattern. See `multi_focal_wall_pattern.md`. |
| Basis enlargement to max_n=4 alone as W1e closure (Sprint F6, 2026-05-23) | 1 | 10.2% PES closure ceiling at max_n=4 (2-point gate overestimates 2.5×). Diminishing returns expected at max_n=5+. See `debug/sprint_f6_maxn4_nah_memo.md`. Methodological note: PES well depth at R_min is the load-bearing metric, not 2-pt differentials. |
| UV refinement to close α > 1 spinor SC gap (G4-4c week 2, 2026-05-29) | 1 | Recovery bit-identical 67.88% across N_0 = 120, 240, 480; the 32% gap is STRUCTURAL not numerical. **Structurally RESOLVED 2026-05-29 v3.19.0 (Track 5)**: closed-form slope $-(1/12)\cdot\alpha/(2\alpha-1)$ at α > 1 (Fursaev-Solodukhin spinor double-cover correction); 2.3% rel err vs measured. See `debug/alpha_gt_1_analytical_investigation_memo.md`. |
| Cutoff-function cure for IR over-count of $S_{\rm BH}$ at small Λ (G4-5c-IR-fix, 2026-05-29) | 1 | Three cutoff variants (Gaussian, polynomial, sharp) tested across six Λ; none achieves factor-2 gate at all Λ. The IR over-count is NOT a tail-suppression issue — substrate's tip(t) peaks at substrate-fixed UV scale t ≈ 0.05–0.1 (Λ-independent), so $r_h^2 \Lambda^2/3$'s Λ² scaling cannot be inherited at fixed substrate. Real cures: sub-leading bulk subtraction, substrate UV refinement, or $A_{\rm eff}$ rescaling. See `debug/g4_5c_ir_fix_S_BH_across_Lambda_memo.md`. |
| Attributing the G7/G4-2 Newton-constant factor of 2 to a scalar-vs-Dirac cone coefficient (2026-05-30) | 1 | WRONG: 2D Dirac and scalar cones share $(1/12)(1/\alpha-\alpha)$ — c=1 both, Fursaev-Miele 1996, and our own G4-4c discrete extraction (−1/12 bit-exact). The factor of 2 is instead a forced bookkeeping artifact (Wald: two-term-exact ⇒ pure Einstein ⇒ action-G ≡ entropy-G), living in a 4D-bulk vs factorized-2D×2D normalization mismatch; close it by convention audit, NOT by tuning the cone coefficient (the one piece confirmed not to carry the 2). See `debug/wald_forces_entropy_relation_memo.md`. |
| Scalar (periodic BC) conical-defect SC extraction at sprint scale (G4-3c-proper / T1, G4-4e, 2026-05-28/29) | 1 | Scalar gets only 28-66% of SC slope at sprint scale, while spinor extracts 99.4%. Scalar zero-mode breaks α ↔ 1/α symmetry on discrete lattice; spinor (anti-periodic, half-integer m) lacks zero mode and cancels UV-asymmetry. **Anti-periodic + half-integer is structurally essential** for clean discrete-substrate SC extraction. See `debug/g4_4e_bc_sectors_memo.md`. |
| Naive K_var/K_disk → 4(l_max+1)(l_max+2) saturation prediction (G4-4b-c, 2026-05-29) | 1 | Lattice spacing $a$ fixes finite S² mass $(n+1)^2/a^2$ at apex regardless of $r_h$; naive "all-modes-massless" limit unreachable on discrete substrate. CORRECT saturation: $K_{\rm var}(r_h \to 0) = K_{\rm cone}$ (cone-Dirac heat trace) bit-exact to 6+ digits. Discrete-substrate gravity has TWO scales: $a$ (substrate UV) and $r_h$ (warp). See `debug/g4_4b_c_asymptotic_free_memo.md`. |
| Small-t polynomial fit for Seeley-DeWitt $a_0$ extraction (G4-4d, 2026-05-29) | 1 | Naive polynomial fit at $t \in [0.02, 0.5]$ gives $a_0$ = 62 vs predicted 50 (off by 24%) due to UV-overshoot bias from T2 high-m azimuthal truncation. Sweet-spot approach at $t = 0.1$ recovers $a_0 = 1.992$ at 99.6%. Continuum extrapolation or multi-window methodology needed for sub-leading coefficients. See `debug/g4_4d_seeley_dewitt_memo.md`. |
| Naive $\phi(2)$ Mellin-moment prediction for $S_{\BH}$ cutoff dependence (G4-5d, 2026-05-29) | 1 | G8 naive prediction $S_{\BH}(f) \propto \phi(2)$ rejected at 65% deviation across Gaussian/sharp/polynomial cutoffs. Structural reframing: tip term lives at $\phi(0)$ Mellin moment (log-regulated), not $\phi(2)$. Sector-wise moment map: tip$\leftrightarrow\phi(0)$ / EH$\leftrightarrow\phi(1)$ / $\Lambda_{cc}\leftrightarrow\phi(2)$. See `debug/g4_5d_cutoff_dependence_memo.md`. |
| Bulk $\Lambda^4$ cosmological-constant extraction from 2D disk-Dirac alone (G4-5b, 2026-05-29) | 1 | The disk-Dirac $K^{\rm Dirac}_{D^2}(t) \sim 2A/(4\pi t)$ has only $1/t$ leading Mellin moment; the $\Lambda^4$ cosmological term emerges only in the 4D embedding $D^2 \times S^2$. 2D substrate gives $c_4 = 0.74$ vs naive 4D prediction 2.53 (ratio 0.29). Confirmed by G4-5c (4D extraction at 0.85 ratio at UV). See `debug/g4_5b_bulk_weyl_extraction_memo.md`. |
| v3.19.0 "Fursaev-Solodukhin spinor double-cover correction" mechanism attribution for Möbius α/(2α-1) (task #26, 2026-05-29) | 1 | hep-th/9512134 is Preitschopf's "Octonions and Supersymmetry", not Fursaev-Solodukhin (fabricated arXiv ID). Correct spinor-on-cone paper is Fursaev-Miele 1996 which states spin 1/2 "resembles the scalar case" — no Möbius modification in published literature. Empirical Möbius match (task #25 sub-2% across 6 α) preserved; mechanism OPEN. See `debug/fursaev_solodukhin_1995_grounding_memo.md`. |
| Geometric-mean azimuthal as cheap UV cure for FD/Spec bracket (task #27, 2026-05-29) | 1 | GM is genuine bracket interior at every t (FD ≤ GM ≤ Spec, F6 bit-exact at α=1, edge ratio 2/π verified) but at t = a² lands at 355.9%, outside [50%, 200%] gate. Substrate refinement (a → a/2) is the real cure, not interpolation between discretization schemes. See `debug/g4_5a_geomean_azimuthal_memo.md`. |
| Single-axis N_φ-sweep for G4-6a refined A-coefficient extraction (2026-05-29) | 1 | A recovery DEGRADES with N_0 (+172.6% at N_0=60 → -18.9% at N_0=480) — opposite of expected convergence. Apparent "A signal" is finite-N_0 artifact from incomplete bulk α·K_disk subtraction. Single-axis substrate refinement does NOT close A coefficient extraction at production substrate values. See `debug/g4_6a_refined_v3_nphi_sweep_memo.md`. |
| Isotropic multi-axis refinement for A-coefficient extraction (2026-05-31) | 1 | Converges at O(a^{0.2}) at intermediate t — infeasible (10% precision requires a/14M). Per-mode Bessel correction overcorrects at t<1 (A is collective/topological, not per-mode). A = 1/(24π) supplied analytically by Theorem 1. See `debug/g4_6a_structural_closure_memo.md`. |
| soft_IR_frac → 1/(2α) as the Möbius α>1 "mechanism" (Route C, demoted by Sprint GD-2 t-audit 2026-05-29) | 1 | The Route C substrate-level "mechanism identification" (soft_IR_frac=1/(2α), sub-percent at t=1) is t≈1-tuned: (1−X)F drifts 0.31→0.59 over t∈[0.25,4], hits 1/2 only at t≈1. Demoted to a sweet-spot coincidence. The Möbius FORM (slope α/(2α−1)-shaped vs SC) is robust; the coefficient and soft-IR explanation are not. Mechanism OPEN. See `debug/sprint_gd2_moebius_t_robustness_memo.md`. |
| Spurious "sign discrepancy" Track α'' thread 9 (2026-05-29, RESOLVED same day) | 1 | The Track α'' thread 9 driver computed dΔ_K/dα = +0.052 at α=2 while v3.19.0 reported -0.0562 — flagged as a sign discrepancy. RESOLUTION: v3.19.0 convention is Δ_K(α)/(1/α-α) ratio, not derivative. Fresh spectral measurement using v3.19.0 convention reproduces -0.0562 bit-exactly. The two drivers measure different observables; substrate-level Möbius identification is supported. See Paper 51 §subsubsec:g4_5_v3_20_followon substrate-discretization-invariance paragraph. |
| Möbius α/(2α−1) as a continuum α>1 closed form (B4, 2026-05-29, RETIRES the open thread) | 1 | The Möbius FORM is a **finite-a substrate artifact**, NOT continuum physics. All prior robustness (GD-2/3/6: t-shape, azimuthal discretization, R) was at FIXED a≈0.05; the one axis never varied — radial continuum a→0 — is where it breaks. Clean spectral-azimuthal + a→0 extrapolation: at α=2,t=1 recovery climbs 0.675→0.919 (N_ρ=200→3200) AWAY from Möbius 0.667 toward Sommerfeld–Cheeger 1.0, with (1−recovery)~√a exactly (ratio √2, cone-singularity signature). Continuum α>1 = plain SC continuation −1/12. No exotic mechanism exists (FS attribution, Routes A/B/C, harmonic conjugate all found nothing because there's nothing). The α=1 replica (1/6, B1/B2) is unaffected — that's the smooth-apex O(a) case. See `debug/l6_b4_moebius_continuum_diagnostic_memo.md`. |
| Naive de-compactification of the finite-R Dirichlet disk (Paper 53 B3, 2026-05-29) | 1 | Taking R→∞ with the finite-R Dirichlet disk does NOT remove the boundary obstruction: the Markov–Cesàro ratio g=num/den is 0/0 at ρ=R (all modes vanish there), so g(R)≠f(R) regardless of f's decay — a construction artifact, not a function-content effect (collar error stays ~0.26 even where f(R−1)=5e−27). The fix is the boundaryless plane ℝ²_α directly via the 2D Bochner–Riesz/Hankel band-limit (no Dirichlet boundary). See `debug/b3_disk_backbone_memo.md`. |
| Species-II spatial fission-aperture discriminator (which-site entanglement frozen=non-binder / responsive=binder, 2026-05-30) | 1 | "Responsive vs frozen" was a LiH-only confound: MgH2 (4e, non-binder) is strongly RESPONSIVE, killing the binding reading; H2 (2e, binder) is FROZEN. The real variable is ≥2 active electron pairs (inter-pair correlation), not binding or fission topology. Spatial fission-aperture reading retired; superseded by the occupation-confinement coordinate (closed subshell = Fock-closure). See `debug/sprint_species_ii_aperture_memo.md`. |
| BW wedge entanglement entropy as Bekenstein-Hawking area law (BH-Phase0, 2026-05-31) | 1 | S ~ 2·log(n_max), NOT area law (R²=0.83 rejected). Boltzmann suppression concentrates probability on ground shell; BH entropy comes from spectral-action replica (Paper 51), not BW entanglement. See `debug/bh_phase0_entanglement_entropy_memo.md`. |
| Furnstahl IR extrapolation for E1 polarizability (2026-06-01) | 3 | Exponential, power-law, and combined models all extrapolate to 2-16 fm^3 (experiment 0.63). Multi-scale continuum response not amenable to single-scale IR correction. See `debug/sprint_cross_observable_nuclear_memo.md`. |
| Sturmian basis for deuteron polarizability with Minnesota NN (2026-06-01) | 3 | v1 grid-coarseness failure; v2 overbinds 20x (Minnesota is effective interaction for HO); v3 refit gives correct binding but alpha-unstable polarizability. Sum-over-states fragile for continuum response; LIT is the correct method. See `debug/sprint_cross_observable_nuclear_memo.md`. |
| LIT for deuteron E1 polarizability (2026-06-01) | 1 | Validated (bit-exact vs eigendecomp) but unhelpful — the +9.7% alpha_E error is from Minnesota (EWSR=152%>TRK), not the solver method. N_shells=2→3 bit-identical. See `debug/sprint_nuclear_tensor_product_memo.md`. |
| NaH Z_orb scan as chemistry-wall diagnostic (2026-06-01) | 1 | PES monotonically descending at all Z_orb ∈ {0.5–2.0}. NaH overattraction is W1e (multi-determinant FCI), not basis extent. Sturmian nuclear improvement does NOT transfer to the chemistry solver. See `debug/sprint_nuclear_tensor_product_memo.md`. |
| Resolvent (D²)⁻¹ for two-body Coulomb interaction (2026-06-01) | 4 | Four weightings tested (Laplacian, Dirac, uniform, 1/N); best Pearson 0.81 DECREASING with n_max. Mismatch is Fock projection conformal factor (Gegenbauer radial on S³ vs Slater in flat space), a function not a constant. See `debug/sprint_resolvent_two_body_memo.md`. |
| J_GV²=−1 as a forcing handle for the inner ℍ factor (Door 4c, 2026-06-01) | 1 | Over ℂ, ℍ and M₂(ℂ) are the same algebra (ℍ⊗ℂ≅M₂); the real-form distinction is an internal involution invisible to the combined J=J_GV⊗J_F (combined J²/signs/KO-dim bit-identical at n_max∈{1,2,3}). J_GV²=−1 admits but does not force ℍ — ℍ is a CCM literature import, not GeoVac-forced. See `debug/door4c_j_signtable_audit_memo.md`. |
| Gauged tensor-product spectral action (full double-sum gauge field) for two-body Coulomb radial weights (2026-06-03) | 1 | Angular selection rules recovered (Gaunt/m-cons/monopole, Paper 54 Thm 3); radial weights do NOT match Coulomb (Pearson 0.58/0.41, decreasing with n_max) — same Fock conformal-factor wall as the resolvent route. Spectral action gives metric, not Green's functions (Bochniak–Sitarz 2022 confirms). See `debug/paper54_two_body_forward_scoping_memo.md`. |
| Yukawa values in low-coefficient pure-Tate periods M1 ∪ M2 (Sprint Yukawa-PSLQ, 2026-06-03) | 1 | 162-cell PSLQ sweep (9 fermions × 3 transforms × 3 ceilings × 2 scales): zero hits at M ≤ 1000 against M1 ∪ M2 basis (basis sharpened by η-trivialization audit ruling out M3 on inner factor). Charged-lepton precision (8 digits at M_Z) honest at M=10. Empirical confirmation of Sprint H1 Yukawa non-selection theorem and Class 1 calibration-data classification. See `debug/sprint_yukawa_pslq_memo.md`. |
| Hopf-tower-to-representation extension as shortcut for forcing N_gen = 3 (Sprint Read 2 scoping, 2026-06-03) | 1 | Naive identification "same 3 (associativity wall) gives both 3 algebra factors and 3 generations" fails on representation theory: in the standard CCM SM rep every generation contains fields from ALL three algebra factors, so the two 3's are different. Multi-year deep wall stands as named (Direction 2 packing-reach NO-GO); no sprint-scale handle. See `debug/sprint_read2_n_gen_scoping_memo.md`. |
| W1e chemistry corrections as outer-factor M1/M2/M3 periods (Sprint W1e period-class, 2026-06-04) | 1 | 0/11 W1e correction terms (NaH F4/F5/F6 sprints) identify with low-coefficient M1, M2, or M3 at audit ceiling 100 or permissive 10⁶; independent random-rational null gives 0/50 on M2/M3. Wrong by structure, not precision: M3 lives at k=1 vertex-parity, W1e has zero vertex-parity content. W1e is the chemistry-side analog of the H1 Yukawa non-selection theorem — calibration-data tier (Paper 18 §IV.6 chemistry-side analog), categorically disjoint from outer-factor periods. See `debug/sprint_w1e_period_class_memo.md`. |
| Propinquity-derived Trotter bound at production parameters (Sprint Trotter propinquity, 2026-06-04) | 1 | Paper 38 GH-rate γ_{n_max} lifted via Childs-Su Duhamel + N_active-linearity gives a single-particle Lipschitz-distortion heuristic estimate of truncation error, not a rigorous bound; overshoots ε/2 budget by 3-4 OoM at production n_max=2 for LiH/BeH₂/H₂O at ε=10⁻³; LiH first becomes feasible at n_max≈5000 (computationally inaccessible). Uniformly looser than naive Suzuki-Trotter at production parameters. The 4/π M1 master-Mellin signature appears only in the (loose, non-binding) truncation budget. Tightening = multi-step research (sharper L_H from multipole structure; sharper L2 Stein-Weiss on low-harmonic subspace), not sprint-scale. See `debug/sprint_trotter_propinquity_memo.md`. |
| DMRG-on-FCIDUMP closes W1e on NaH (Sprint R3-B falsifier, 2026-06-07) | 1 | DMRG=FCI at finite sector dimension (NaH balanced n_max≤3, FCI dim ≤784): PES bit-identical to P4 baseline at every R, no interior minimum, monotone over-attraction into small R persists both n_max=2 and n_max=3 and under post-bug-fix Path B convention (W1e-Projection-Audit 2026-06-07). W1e is at the projection step from continuous Level 4 adiabatic+PK to second-quantized (h1, eri, ecore) integrals for SECOND-ROW chemistry, NOT at the determinant-expansion level — wall localized at the Hamiltonian-specification level. The bug fix raised the true overattraction depth at R=2.5 by +1.3 Ha; qualitative verdict unchanged. See `debug/sprint_{r3b_dmrg_nah_falsifier,w1e_projection_audit}_memo.md`. |
| LiH composed qubit FCI binds (R3-A diagnostic, 2026-06-07) | 1 | LiH **composed** qubit FCI gives monotone-descending PES across R ∈ [2.5, 5.0] bohr; no interior minimum. The CHANGELOG v3.56.0 2.82% R_eq comes from `ComposedDiatomicSolver.LiH_ab_initio` (CONTINUOUS Level 4 multichannel adiabatic + PK), NOT from qubit FCI on the FCIDUMP-exported integrals. Composed structurally lacks cross-block coupling (no inter-block V_ne or cross-block ERIs); its electronic energy is R-independent and only V_NN(R) varies → monotone descent. **2026-06-07 W1e-Projection-Audit correction: the original claim "LiH composed AND balanced" was wrong about balanced** — R3-A used a Path B convention that triggered a V_NN-double-count bug in the balanced corrector. Under the correct convention (or after the 2026-06-07 fix), balanced LiH binds at R_eq = 3.015 bohr with D_e = 0.158 Ha (2.4× over-binding vs continuous, but bowl-shaped). Composed row stands; balanced LiH retracted from this dead-end row. See `debug/sprint_{r3a_dmrg_lih,w1e_projection_audit}_memo.md`. |
| Existing engineering kwargs close LiH over-binding (LiH kwarg sweep, 2026-06-07) | 1 | 11 kwarg combinations tested on balanced LiH. 8 give bit-identical energy to baseline because `multi_zeta_basis`, `screened_valence_basis`, `screened_cross_center`, `pk_cross_center` are gated on `Z_nuc_center >= 11` and silently no-op for first-row; `cross_block_h1=True` makes LiH 16× over-bound (F3-style "bonding without Pauli repulsion" failure). Cheap engineering arsenal exhausted for first-row chemistry-accuracy. See `debug/lih_kwarg_sweep_log.txt` and `debug/sprint_w1e_audit_and_mvs_chemistry_open_2026_06_07_memo.md`. |
| Explicit-core HF closes NaH binding (Sprint B.1, 2026-06-07) | 1 | 12-electron RHF SCF on explicit [Ne] core NaH (30 qubits, 15 spatial orbitals at max_n=2 on Na + max_n=3 effective on bond block; DIIS + density damping + level shift) gives monotone descent into small R; 5/7 R points converge, 2 SCF-bistable in the bond-formation region. LiH HF cross-check binds at R_eq = 3.015 with D_e = 0.158 Ha — methodology validated; NaH negative is structural, not an SCF artifact. **Frozen-core projection is NOT the chemistry-accuracy wall**; un-freezing the [Ne] core does not close NaH binding. Combined with W1e-Projection-Audit + LiH kwarg sweep + F1–F6 history, the chemistry-engineering arc is empirically exhausted at second-row. See `debug/sprint_b1_nah_hf_verdict_memo.md`. |
| W1e closure via off-diagonal cross-block h1 WITHOUT orthogonalization (SO(4)-breaking diagnostic, 2026-06-07) | 1 | Per-sub-block h1 eigenvalue analysis on LiH+NaH balanced (5–6 R points each) identifies W1e as a non-orthogonality / Gaunt-sparsity tension in the per-sub-block hydrogenic basis. Diagonal cross-V_ne contribution is correct (~−Z_heavy/R on H sub-block, exact for both LiH and NaH across R sweep). Off-diagonal inter-block coupling needed for binding requires Löwdin orthogonalization, which destroys Gaunt sparsity (14× Pauli inflation per Track DF Sprint 5). **The W1e wall is structurally inseparable from the framework's qubit-encoding advantages** (angular sparsity Paper 22, gauge-network identification, Z₂ tapering schemes). Cross_block_h1 architectural extension (Sprint F3) demonstrated this: adding off-diagonal coupling on the non-orthogonal basis gives 16× over-binding on LiH. Closure within the framework's architecture requires accepting either binding-wall OR sparsity-loss; the framework chose sparsity. Cumulative dead-ends: F1–F6, Schmidt, core correlation, LiH kwarg sweep, B.1 explicit-core HF, R3-A/B DMRG. See `debug/sprint_so4_breaking_w1e_diagnostic_memo.md`. |
| Spectral action expansion has chemistry-side structural content at finite cutoff (Sprint spectral action expansion, 2026-06-07) | 1 | $S(D, \Lambda) = \mathrm{Tr}\,\exp(-D^2/\Lambda^2)$ on the M-vS-2-confirmed bit-exact chemistry Dirac (D = h₁) reduces BIT-EXACTLY to the trivial $\sum_k (-1)^k \mathrm{Tr}(D^{2k})/(k!\Lambda^{2k})$ Taylor series at large Λ. Verified across LiH (M=15), H₂ (M=10), NaH (M=10) at Λ ∈ [0.1, 100], rel. residual ≤ $10^{-15}$ for Λ ≥ 10. **No Chamseddine–Connes Seeley–DeWitt hierarchy emerges** — structurally forced by finite-dim $D$ having no UV divergence. PARTIAL: $\mathrm{Tr}(h_{\rm off}^2)$ (bond-coupling intertwiner Frobenius norm-squared) IS chemistry-meaningful (H₂ 0.37, NaH 0.31, LiH 0.10 ratio to full), but as standard linear-algebra invariant, not as emergent spectral-action content. Combined with M-vS-2 Q2 NEGATIVE (monotone $S(D)(R)$ doesn't bind, v3.87.0), **spectral action is decisively the wrong functional for chemistry observables at finite cutoff**; CC machinery requires continuum limit. See `debug/sprint_spectral_action_expansion_chemistry_diagnostic_memo.md`. |
| Camporesi-Higuchi κ-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint CH κ-parity, 2026-06-07) | 1 | $P_\kappa = \prod_{q:\kappa_q < 0} Z_q$ does NOT commute with relativistic Tier 2 chemistry $H_{\rm rel}$ on LiH_rel/BeH_rel/CaH_rel. Commutator residuals 7–9 OoM above $10^{-10}$ gate (LiH_rel $5.07\times10^{-2}$, BeH_rel $1.29\times10^{-1}$, CaH_rel $3.38\times10^{-2}$) with PK on/off and Breit on/off. **Mechanism (clean structural):** jj-coupled full-Gaunt $X_k$ angular coefficient has parity selection $(l_a + l_c + k)$ even with NO κ-dependence, so Coulomb operator freely couples $p_{3/2}$ ($\kappa = -2$) with $p_{1/2}$ ($\kappa = +1$) at same $l = 1$, flipping κ-sign. **38% of LiH_rel ERI tuples carry odd $\Delta N_{\kappa<0}$.** Dirac-Coulomb operator conserves $j$ and $m_j$ but NOT κ-branch. Relativistic-chemistry analog of H1 Yukawa non-selection theorem + W1e period-class structural disjointness. ΔQ = 0 from this tapering. Scaffold + audit module `geovac/relativistic_tapering.py` shipped for future Kramers / $m_j$-parity probes ($m_j$ IS conserved so those might fare better). See `debug/sprint_ch_kappa_parity_z2_memo.md`. |
| Direct-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity direct, 2026-06-08) | 1 | $P_{m_j} = \prod_{q:m_j(q)<0} Z_q$ in original Dirac basis also does NOT commute with $H_{\rm rel}$. Residuals 6–24 ×10⁻³ across LiH_rel/BeH_rel/CaH_rel — same order as κ-parity sprint. **Mechanism (sharper):** total $M_J$ conservation in the jj-coupled ERI is a SUM constraint ($m_j^a+m_j^b = m_j^c+m_j^d$), NOT a parity constraint. Concrete violation: $(+3/2)+(-1/2)\to(+1/2)+(+1/2)$ preserves $M_J=1$ but flips $\Delta N_{m_j<0}=+1$. Joint with κ-parity NEGATIVE: **original-basis Z-string parities in relativistic chemistry fail because the jj-coupled Coulomb freely couples sign-counted sectors under M_J/κ-branch sum-conservation rules**. The structural distinction: non-rel Hopf m_l Z₂ DOES commute because rotation to (sym, antisym) basis isolates the antisym sector cleanly; the relativistic analog (rotated-basis m_j → −m_j) is the named follow-on but requires `composed_qubit_relativistic` refactor to expose dense h1 + densify eri_sparse (~2-3 day sprint). See `debug/sprint_mj_parity_z2_memo.md`. |
| Rotated-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity rotated, 2026-06-08) | 1 | $P_{m_j}^{\rm rot}$ on the m_j ↔ −m_j (sym, antisym) rotated basis fails too (residuals 1.7–6.8×10⁻², densification bit-exact PASS rules out implementation bug); the 4-orbital product of 3-j sign factors $(-1)^{j_a+j_b+j_c+j_d}$ is not forced even for half-integer j, whereas the non-rel integer-l analog IS forced by the Gaunt parity selection rule — closes the three-sprint relativistic-Z₂ tapering thread (κ-parity v3.92.0, m_j-direct v3.94.0, m_j-rotated v3.95.0). See `debug/sprint_mj_parity_rotated_memo.md`. |
| Spectral action $S(D)(R)$ of assembled MvS Dirac binds LiH (Sprint M-vS-2 + R-sweep, 2026-06-07) | 1 | At default `lih_spec()` $n_{\max}=2$, the spectral action $S(D)(R) = \mathrm{Tr}\,\exp(-D^2/\Lambda^2)$ of the bit-exactly-assembled Marcolli-vS gauge-network Dirac is **monotone increasing** in $R$ at every $\Lambda \in \{1, 2, 4\}$ over $R \in [2.0, 8.0]$ bohr (10 panel points). $\mathrm{Tr}(D^2)$ is monotone decreasing. $E_{\mathrm{FCI}}$ is monotone decreasing (over-binding signature, W1e wall). **No interior minimum at any tested $\Lambda$ for any of the five observables.** Same shape as $E_{\mathrm{FCI}}$: the M-vS-native spectral-action observable does NOT autonomously bind chemistry where chemistry-style FCI fails to bind. **W1e is at the PROJECTION step (continuous Level-4 PK-composed → second-quantized integrals), not at the EVALUATION step (FCI vs spectral action).** Closes the reconciliation question opened mid-session and confirms the multi-focal-composition wall pattern is independent of the evaluation choice on the gauge-network data. See `debug/sprint_mvs2_lih_default_plus_rsweep_memo.md`. |
| Non-abelian M-vS gauge group reduces Pauli count beyond Z₂ Hopf-U(1) tapering (Sprint M-vS Gauge, 2026-06-07) | 2 | 5 gauge candidates tested on default LiH (M=15, Q=30): identity → 909 Pauli pre-taper, **813 post-Z₂ on 25 qubits** (the OPTIMUM). Per-sub-block h1 diag → 2069 pre/3925 post-Z₂ (2.3× denser). Per-vertex NO from FCI 1-RDM → 17173 (18.9× denser). Random block-unitary → 25317 (27.8× denser). H₂ cross-check matches pattern (514 → 1158, 2.25× denser, FCI bit-exact). **Z₂-tapered identity is optimum across all tested candidates.** Structural reason: tapering requires Hamiltonian SYMMETRY; M-vS gauge U(H_v) per vertex is a basis-change FREEDOM, not a symmetry. Hopf m_l → −m_l Z₂ is the only discrete symmetry of the GeoVac chemistry construction; the continuous non-abelian extension has no Pauli-reduction-side content. Closes one of three "M-vS upgrades chemistry" candidate paths NEGATIVE. See `debug/sprint_mvs_gauge_pauli_reduction_memo.md`. |
| K⁺ compression of the Krein Dirac as a Lorentzian quantum-metric device (P45 descope, 2026-06-09) | 1 | Krein-self-adjointness of $i\,D_{GV}\otimes I$ forces $\{J, D_{GV}\otimes I\}=0$, so $P_+ D_{GV} P_+ = 0$ exactly and the restricted Lipschitz seminorm ≡ 0 on the entire operator system (bit-exact at (2,3),(3,5); kernel 42/42, 275/275). "First Lorentzian propinquity theorem" claim withdrawn; falsifier frozen in `tests/test_p45_kplus_degeneracy.py`. See `debug/sprint_p45_hardening_phase1_memo.md`. |
| Momentum-diagonal temporal multipliers as a metric temporal algebra (P45, 2026-06-09) | 1 | $g_p(\omega_k)$ are functions of momentum, not time — they commute with the Fourier-diagonal $D_t$ exactly, so time is Lipschitz-invisible by construction (the celebrated "L3 structural identity" is the diagnosis, not a feature). Repair requires Toeplitz compressions of $e^{iqt}$ (Connes-vS S¹ pattern). Same memo. |
| Grab-bag / single-weight PSLQ bases for weight-inhomogeneous graph-QED constants (2026-06-11) | 1 | S_min's 15 "irreducibility" failures were a basis-coverage artifact (the bases held 2/8 of the weight-5 level-2 monomials); weight-grade the object first, then PSLQ against complete weight-homogeneous bases. See `debug/smin_dossier_round1_memo.md`. |
| Levin/EM series acceleration on log-modulated summands (2026-06-11) | 3 | mpmath nsum-levin silently mis-converges when the summand carries a log factor (harmonic/digamma growth), with precision-DEPENDENT error — two runs at different dps disagreeing is the tripwire. Poisoned the S^(3) trailing-1 t3 evaluator, the factorized anchor (both stage-1 figures wrong), and the published N^−1.31 "power law" ((ln N)²/N² masquerade). Fix: remove the log by Abel summation, or rigorous brackets with no acceleration; cached term arrays must honor nsum's precision contract. See `debug/sprint_s3_closure_memo.md`. |
| Nested adaptive nsum cascade for depth-≥4 multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 | A d−1-level nested `mpmath.nsum(levin)` cascade (each level's terms are themselves accelerated sums) is intractable past ONE acceleration layer: tractable at depth 3, catastrophic at depth 4 (even fast-decay t4(4,4,3,2) at 30 dps never returned; 13 CPU-h, 0/40 t4). Use the single-nsum architecture (largest var = closed Hurwitz tail, smallest = closed partial; one numeric layer). See `debug/sprint_s4_scoping_memo.md` §S1.6b. |
| `mpmath.sumem` (auto Euler-Maclaurin) for log-modulated trailing-t tails (S^(4) stage-1, 2026-06-13) | 1 | Returned a TWO-PRECISION-STABLE but WRONG value (2.4e-8 off at b1=4; ~1e+175 garbage with a shifted start) — boundary mishandling at the e_k leading zeros. **sumem two-precision-stable ≠ correct** (companion to the Levin/log row). Validate against a rigorous partial-sum ground truth, not self-consistency. See `debug/sprint_s4_scoping_memo.md` §S1.6d. |
| float64-contaminated terms inside a high-precision Levin sum (S^(4) stage-1, 2026-06-13) | 1 | `int ** (-int)` in a prefix array silently yields float64; Levin amplifies the 1e-16 jitter by its condition number (worst on the slowest series), giving a weight-graded, precision-dependent error. Every power inside an accelerated sum must be `mpf`. See `debug/sprint_s4_scoping_memo.md` §S1.6. |
| Brute high-precision summation of b1=2 high-log trailing multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 | Four methods (Abel/Levin, sumem, manual E-M, analytic-tail agent) all under-converge / are wrong / are impractically slow — this is the known-hard MZV-computation problem. Reframe: these are LOW-WEIGHT hence classically reducible, so identify them SYMBOLICALLY (stage 2) rather than sum them numerically. See `debug/sprint_s4_scoping_memo.md` §S1.6f. |
| sympy exact-Q rank for MZV-scale relation matrices (S^(4) stage-2, 2026-06-13) | 1 | sympy `Matrix.rank()` over Q is intractable at weight 11-13 (hundreds-to-thousands of admissible words; >17 min, no result). Use mod-p (large prime) streaming Gaussian elimination, cross-checked at a second prime; project to the target subspace during row construction. Standard for MZV rank computations. See `debug/sprint_s4_stage2_memo.md`. |


---

## 3.5. Guardrail Papers

Certain papers contain proven theorems or exhaustive diagnostic arcs that constrain what approaches are viable. Sub-agents MUST load the relevant guardrail paper before beginning any investigation in the constrained domain. If a proposed investigation falls within a guardrail's scope, the sub-agent must:

1. Load the guardrail paper into context
2. Identify the specific theorem or negative result that applies
3. WARN THE USER: "Paper [X] documents a proven negative result for this class of approach: [theorem statement]. The proposed investigation [specific proposal] appears to fall within scope. Do you want to proceed anyway, or modify the approach?"
4. If the user proceeds: document in the investigation prompt that the guardrail was acknowledged and state why the new attempt differs from the documented negative result

This protocol exists because Track DF (6 sprints, April 2026) re-derived at significant cost two negative results that were already proven in Papers 8-9 and FCI-M. The re-derivation produced genuine new findings (6j sparsity theorem, compact encoding, three-layer structure), but Sprints 4 and 4B could have been avoided entirely if the guardrail papers had been loaded at the investigation design stage.

### Current Guardrail Papers

| Paper | Domain | Theorem/Result | Scope |
|:------|:-------|:---------------|:------|
| Paper 8-9 | Single-center molecular encoding | Sturmian Structural Theorem: H_ij ∝ S_ij in shared-p₀ basis → eigenvalues R-independent | ANY investigation proposing to encode a heteronuclear molecule in a single-center (single-Z, single-k, single-p₀) basis, including nested hyperspherical, bond sphere, Sturmian CI, or unified orbital approaches |
| FCI-M | Graph-concatenation molecular encoding | Graph Laplacian kinetic energy is R-independent in LCAO basis → monotonically attractive PES, no equilibrium | ANY investigation proposing to concatenate atom-centered graphs into a molecular graph without natural geometry coordinates |
| Track DF record | Nested hyperspherical molecular encoding | Three molecular variants tested (single-center, charge-center, heterogeneous Löwdin) — all NEGATIVE | ANY investigation proposing to place all molecular electrons in a single S^(3N-1) Hilbert space. Extends Paper 8-9's theorem to the hyperspherical setting with additional finding: Löwdin orthogonalization of mixed-exponent bases destroys Gaunt sparsity |

### Guardrail Triggers

When a user or investigation prompt proposes work in any of the following categories, the sub-agent MUST load the indicated guardrail paper BEFORE writing any implementation plan:

- **"single center" + "molecule"** → Load Papers 8-9
- **"unified basis" + "molecule"** → Load Papers 8-9
- **"shared exponent" OR "single k" OR "single p₀" + "molecule"** → Load Papers 8-9
- **"nested" + "molecule"** → Load Papers 8-9 + Track DF record
- **"graph" + "concatenat"** → Load FCI-M
- **"LCAO" + "graph"** → Load FCI-M
- **"Sturmian" + "molecule"** → Load Papers 8-9

After loading, the sub-agent must check whether the proposed approach falls within the guardrail's scope and warn the user if it does. The warning should be constructive, not blocking — it documents what DOESN'T work, which helps refine toward what MIGHT work.

---

## 4. The Dimensionless Vacuum Principle

The graph Laplacian is dimensionless and scale-invariant, topologically equivalent to the unit three-sphere S3. Physical energies emerge only through the energy-shell constraint p0^2 = -2E, which acts as a stereographic focal length. The 1/r Coulomb potential is not an input force law -- it is the coordinate distortion from stereographic projection (chordal distance identity). The universal kinetic scale kappa = -1/16 maps graph eigenvalues to the Rydberg spectrum. Eigenvalues of the Laplace-Beltrami operator on unit S3 are pure integers: lambda_n = -(n^2 - 1).

For molecules, the natural geometry shifts from S3 (atoms) to prolate spheroidal coordinates (Paper 11), hyperspherical coordinates (Paper 13), molecule-frame hyperspherical (Paper 15), or composed fiber bundles (Paper 17). The choice of geometry is determined by where separation of variables occurs. At Level 1, the separated coordinates are all angular (on S3), making the problem fully discrete. At Levels 2+, separation produces continuous radial-like coordinates alongside discrete angular channels — a consequence of SO(4) symmetry breaking by multi-center or multi-electron potentials (Papers 8-9, negative theorem).

**Prime directive (level-aware):** The framework preserves two distinct kinds of structure, and the rules differ by level:

*Angular/symmetry structure (discrete at all levels):* The quantum number labeling (n, l, m), the per-shell degeneracy 2l+1, the selection rules from Gaunt integrals, and the channel structure are combinatorial invariants derived from the packing construction (Paper 0). These must never be modified, approximated, or bypassed at any level. They are the framework's core identity.

*Radial/parametric structure (level-dependent):* At Level 1, Fock's SO(4) symmetry converts the radial coordinate into an angular coordinate on S3, making the entire problem discrete. The graph Laplacian is the exact object and must not be modified to artificially recover continuous differential terms (like 1/r or nabla^2). At Levels 2-4, the SO(4) symmetry is broken by multi-center or multi-electron physics, and the radial-like coordinates (internuclear distance R, hyperradius R_e, prolate spheroidal xi) currently require continuous numerical methods. However, this is not necessarily a fundamental feature of the physics — it may reflect algebraic structures that haven't been found yet. The project has already demonstrated this replacement in several cases: Paper 12's Neumann expansion replaced numerical quadrature of V_ee with algebraic recurrence relations; Gaunt integrals replace angular integration at all levels; the split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. In each case, what appeared to require continuous integration was replaced by algebraic evaluation once the right structure was identified.

The places where continuous numerical methods currently survive are: Z_eff screening quadrature in composed geometries, spline caching for the rho-collapse, finite-difference radial grids in hyperspherical solvers, and the adiabatic potential curves U(R). These are computational tools for evaluating coupling matrix elements between discrete channel structures. Standard numerical methods (grid refinement, basis convergence, quadrature improvement) are appropriate for improving accuracy in these areas. Finding algebraic replacements for these numerical steps — as Paper 12 did for V_ee — is an ongoing goal of the project.

The framework's long-term aspiration is: for every natural geometry, there exists an algebraic structure that computes all coupling matrix elements without spatial integration. This has been achieved for Level 1 (S3 graph eigenvalues), for angular couplings at all levels (Gaunt integrals), and for electron-electron repulsion in prolate spheroidal coordinates (Neumann expansion). It has not yet been achieved for the hyperradial coupling in Levels 3-4 or for Z_eff screening in Level 5. What must be preserved at every level is the discrete channel structure and selection rules, not the specific numerical method used to evaluate radial amplitudes within those channels.

**Algebraic-first, observation-aware (refinement after Papers 18 / 34, May 2026):** The aspiration above is the right discipline for Layer 1 (the bare graph) and for the discrete channel structure at every level. After Paper 34's two-layer framing, the project distinguishes algebraic content (Layer 1: π-free, integer / rational / algebraic-extension matrix elements) from observation-side content (Layer 2: where projections to physical observables introduce specific transcendentals — π via the Hopf measure, π^{2k} via the spectral action, ζ(2k) via even-zeta Dirichlet series, ζ(3) via half-integer Hurwitz, Catalan G via vertex parity, 2π via temporal compactification when an observer integrates over a finite time window). Irreducible transcendentals that survive algebraic decomposition are not failures of the algebraic-first discipline — they are the *content* of a specific Paper 34 projection and should be pinned to that projection rather than chased indefinitely. The diagnostic question for any quadrature wall is therefore two-headed: (i) is the wall a missing algebraic structure (decompose, per Paper 12 / cross-block V_ne / hypergeometric Slater), or (ii) is it the irreducible signature of an observation-side projection (catalogue against Paper 34 and stop)? Both answers are valid outcomes; conflating them was the error the May-2026 curve-fit audit flagged.

*Practical test:* If a proposed modification changes which quantum numbers label the states or which transitions are allowed, it violates the prime directive. If it changes how accurately the radial amplitude is computed within a given channel, it is a legitimate numerical improvement. If it adds a transcendental, the transcendental must be tagged to a Paper 34 projection — anonymous transcendentals are not allowed in production code or papers.

**Topological integrity tests:** The 18 symbolic proofs in `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py` validate S3 conformal geometry. These must never be broken or bypassed. Run before any release.

---

## 5. Natural Geometry Hierarchy

The core organizational principle of the project. Each electron configuration has a natural coordinate system where the physics separates.

| Level | System | Natural Geometry | Best Result | Paper |
|:-----:|:-------|:-----------------|:------------|:-----:|
| 1 | H (1-center, 1e) | S3 (Fock) | < 0.1% | 7 |
| 2 | H2+ (2-center, 1e) | Prolate spheroid | 0.0002% (spectral) | 11 |
| 3 | He (1-center, 2e) | Hyperspherical | 0.019% (2D var + self-consistent cusp); 0.20% (graph-native CI n_max=9, 0 params, exact algebraic integrals) | 13 |

*Level 3 note:* The previous 0.05% result used the FD adiabatic solver (lucky error cancellation, non-variational). The adiabatic coupled-channel solver converges to a structural floor of 0.19-0.20% (v2.0.8). The 2D variational solver (Track DI, v2.6.0) breaks this floor: raw 0.022% at l_max=7 (tensor-product Laguerre × Gegenbauer basis, 8000 dim), cusp-corrected 0.004% at l_max=4. The 2D solver treats R and α simultaneously, capturing non-adiabatic R-α correlation that the adiabatic approximation misses. l_max convergence is monotonic; the per-channel angular basis (n_basis=40 Gegenbauer functions) is the convergence bottleneck, not partial-wave truncation.
| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 96.0% D_e | 15 |
| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | R_eq 63.5% (l_max=2, 2D variational; unbound D_e) | 17 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq 5.3% | 17 |
| 5 | BeH₂ (polyatomic) | Composed (Level 3 + 4) + exchange | R_eq 11.7% | 17 |
| 5 | H₂O (triatomic) | Composed (Level 3 + 4) + lone pairs | R_eq 26% | 17 |
| Nested | Be (1-center, 4e) | S¹¹ (SO(12)), H-set | Q=10, 112 Pauli, 4.9% | Track DF |

*Level 4N note:* Level 4N is the exact N-electron generalization of Level 4's mol-frame hyperspherical coordinates (SO(3N) replacing SO(6), S_N antisymmetry replacing the gerade constraint). Level 5 (composed geometry) approximates Level 4N, trading exact inter-group antisymmetry for 144x angular compression via PK. The l_max=2 result (R_eq ≈ 1.1 bohr, 63.5% error) demonstrates that equilibrium exists without PK. Track AR (v2.0.23) confirmed that D_e overcounting was an adiabatic artifact: the 2D variational solver gives E_min = -7.79 Ha (variational bound respected, above exact -8.07) with D_e = -0.19 Ha (unbound), versus the adiabatic solver's D_e = +0.49 Ha (5.3x exact). The l_max=2 S₄ [2,2] angular basis is genuinely insufficient for LiH binding. Track AS (v2.0.23) confirmed composed encoding is categorically sparser than full N-electron encoding (334 vs 3,288 Pauli terms, 20x lower 1-norm). Composed geometry's 144x angular compression and structural sparsity are essential, not approximations of convenience.

The composed geometry (Level 5) is a fiber bundle: G_total = G_nuc semi-direct G_core(R) semi-direct G_val(R, core_state). Each electron group gets its own natural coordinate system, coupled via Z_eff screening and Phillips-Kleinman pseudopotential.

**Algebraic structure:** At every level, angular matrix elements are computed from quantum number labels and Wigner 3j symbols (via Gaunt integrals), with no spatial quadrature. The split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. At Level 2, the radial solver is fully algebraic for all m: σ states (m=0) use ordinary Laguerre three-term recurrence with zero numerical integration (v2.0.9); π/δ states (m≠0) use associated Laguerre basis L_n^{|m|}(x) with partial-fraction decomposition and Stieltjes integral recurrence, reducing non-algebraic content to a single transcendental seed e^a·E₁(a) (Track J, v2.0.10). At Level 3, the angular problem is fully algebraic, and the hyperradial overlap S and kinetic K matrices are algebraic via three-term Laguerre recurrence (Track H, v2.0.10: pentadiagonal M2 for S, tridiagonal derivative expansion for K, < 1e-14 relative error, 11× build speedup). The adiabatic eigenvalues μ(R) are proven transcendental (O(R) → O(R²) regime transition, v2.0.9 Track G) — the potential V_eff(R) must stay quadrature, point-by-point diagonalization is irreducible, though spectral radial solvers achieve 95-120× speedups. At Level 4, the spectral Laguerre basis achieves 16× dimension reduction for the hyperradial coordinate (Track I), and the spectral Jacobi basis achieves 20× dimension reduction for the angular sweep (Track K, 269× speedup, 1000→50 matrix dimension). The combined spectral solver reduces the angular sweep from 99% to ~50% of total cost. The angular eigenvalues μ(ρ) remain transcendental (computed by diagonalization at each ρ-point), but the per-point cost is now a 50×50 eigensolve rather than a 1000×1000 one. n_basis_radial=20 optimal (mild conditioning at n_basis≥25). Spatial quadrature enters for: Level 3 hyperradial potential matrix elements (V_eff transcendental), Level 4 angular eigenvalue sweeps (μ(ρ) transcendental, but spectral basis reduces cost 269×), Z_eff screening in composed geometries, and rho-collapse spline caching.

---

## 6. Paper Series

> Newcomer-facing status map with one-liners: `papers/INDEX.md`. Detailed per-paper notes (the long key-result descriptions, frozen at v3.110.0): `docs/paper_notes_archive.md`. Topic → paper lookup: `docs/topic_to_paper_lookup.md`. This section keeps only what agents need at dispatch time: loading tiers, the folder map, and live status flags.

### Loading tiers

**Always load** (framework identity): Paper 0 (packing axiom, K = −1/16) · 1 (spectral graph methods) · 7 (S³ proof, 18 symbolic proofs) · 14 (qubit encoding headline) · 16 (S_N periodicity) · 22 (angular sparsity theorem) · 23 (nuclear hub, Fock rigidity) · 24 (Bargmann-Segal S⁵, Coulomb/HO asymmetry) · 27 (entropy as projection) · 31 (universal/Coulomb partition) · 32 (the spectral triple; §VIII theorems).

**Load on topic** (full list and statuses in `papers/INDEX.md`): chemistry solvers → 8–9, 11, 12, 13, 15, 17, 19, FCI-A/M; QC resources → 20; QED/gauge/gravity → 2, 25, 28, 30, 33, 36, 41, 51; math.OA arc → 29, 38, 39, 40, 42–50, 52, 53; foundations/periods → 18, 54, 55, 56, 57; precision → 26, 34, 35.

**GUARDRAIL papers** — MUST load before any investigation in their domain (trigger words and protocol in §3.5): Papers 8–9 (single-center / unified-basis / shared-exponent molecular — Sturmian structural theorem), FCI-M (graph-concatenation molecular), Track DF record (nested hyperspherical).

### Folder organization (audience groups, reorganized 2026-05-22)

| Folder | Audience | Papers (.tex) |
|:-------|:---------|:-------------:|
| `papers/group1_operator_algebras/` | math.OA / NCG | 16 |
| `papers/group2_quantum_chemistry/` | quantum chemists | 9 |
| `papers/group3_foundations/` | mathematical physicists | 11 |
| `papers/group4_quantum_computing/` | QC / NISQ / VQE | 4 |
| `papers/group5_qed_gauge/` | HEP / gauge theory | 8 |
| `papers/group6_precision_observations/` | precision AMO | 4 |
| `papers/synthesis/` | cross-group narratives + field guide | 3 |
| `papers/archive/` | historical (3, 4, 5, 6, 10, 18v1, 21) | 7 |

### Live status flags every agent must know

- **Paper 38 UNCONDITIONAL** (2026-06-10, translation-seminorm metrization; falsifier `tests/test_p38_action_seminorm.py`). The WH1 keystone.
- **Paper 45 DESCOPED + partially rebuilt** (2026-06-09 K⁺ theorem withdrawn, falsifier `tests/test_p45_kplus_degeneracy.py`; 2026-06-10 product-carrier convergence restored in the action-seminorm framework, `prop:product_action_seminorm`, falsifier `tests/test_wh7_b1_joint.py` — signature-agnostic, NOT a Lorentzian claim). Do NOT cite pre-descope claims. **Paper 46 DESCOPED; Papers 47/48/49 PARTIAL** (in-paper Status notes; the norm-resolvent arrow and the TICI/cocycle algebra survive).
- **Paper 2 is an Observation**; the combination rule K = π(B + F − Δ) stays labeled conjectural (§13.5 hard prohibition).
- **Paper 34** is the living projection catalogue (28 projections); **Paper 18 §III.7** is the master Mellin engine; tag every transcendental against both (memory rule).
- Papers are corrected **in place** (de-versioning directive 2026-06-10); git/Zenodo are the version record. No splinter files.

---

## 7. Code Architecture

> Full entry-point catalogue (~120 rows) and solver-method table: `docs/code_architecture.md` (live document — update there, not here). Most-used entry points:

| Task | Module · Entry point |
|:-----|:---------------------|
| Atomic lattice / Hamiltonian | `geovac/lattice.py` `GeometricLattice(Z, max_n)` · `geovac/hamiltonian.py` `GraphHamiltonian(lattice)` |
| Multi-electron FCI | `geovac/lattice_index.py` `LatticeIndex(Z, n_electrons, max_n)`; direct CI at N_SD ≥ 5000 |
| Molecular spec + composed builder | `geovac/molecular_spec.py` `MolecularSpec` · `geovac/composed_qubit.py` `build_composed_hamiltonian(spec)` + `*_spec()` factories |
| Balanced coupled builder | `geovac/balanced_coupled.py` `build_balanced_hamiltonian(spec, nuclei)` |
| Ecosystem export | `geovac/ecosystem_export.py` `hamiltonian(name, tapered=None/'global'/'per_block'/'extended'/'full')` → `.to_qiskit()/.to_openfermion()/.to_pennylane()` |
| Z₂ tapering | `geovac/z2_tapering.py` `apply_hopf_tapering()` · `geovac/extended_tapering.py` |
| Frozen cores / cross-center V_ne | `geovac/neon_core.py` `FrozenCore(Z)` · `geovac/shibuya_wulfman.py` `compute_cross_center_vne()` |
| Slater integrals (exact) | `geovac/hypergeometric_slater.py` `compute_rk_float()` (threshold-dispatched n ≥ 5 → exact Fraction) |
| Operator system / Connes distance / GH | `geovac/operator_system.py` `TruncatedOperatorSystem(n_max)` · `geovac/connes_distance.py` · `geovac/gh_convergence.py` |
| Physical constants | no central module; `-1/16` may be used directly (§8); `ALPHA`/`C_LIGHT` live next to their modules |

---

## 8. Coding Standards

### Sparse vs Dense: Context-Dependent

- **Hamiltonian and CI matrices (N > 100):** Always `scipy.sparse` (csr_matrix, coo_matrix). Never densify.
- **Hot-loop lookup tables (ERI, h1 in direct CI):** Use dense NumPy when array fits in memory (n_spinorb <= ~300). `scipy.sparse._validate_indices` overhead (~24us/call) is prohibitive at 100K+ lookups.

Rule of thumb: sparse for the physics matrix (N_SD x N_SD), dense for orbital-index lookup tables (n_spinorb x n_spinorb or n_spatial^4).

### Type Hints Required

All function signatures must have type hints.

```python
def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    ...
```

### Physical Constants

Import from `geovac.constants` or define at module top. No hardcoded magic numbers.
**Exception:** `-1/16` is the universal topological constant (can be used directly).

### Vectorization Over Loops

Avoid Python loops for graph operations; use NumPy masking/vectorization.

```python
mask = (n_values >= 1) & (l_values < n_values)
states_filtered = states[mask]
```

---

## 9. Workflow Protocols

### Theory Check Rule

Before implementing new physics:
1. Check `papers/` for the derivation
2. If code contradicts paper -> flag it and ask user
3. If changing physics in code -> prompt user to update papers

### Current-State Check

Before forming or reporting a verdict on any question — *especially* when resuming a thread from a `debug/` memo — verify the CURRENT state, not the snapshot. The papers (the section that *owns* the question) + CHANGELOG since the memo's date are canonical; `debug/` memos are dated snapshots the corpus moves past, and a memo's "open question / next step" may already be closed. First move when picking up a thread = read the owning paper section + post-memo CHANGELOG, THEN conclude. (Added 2026-06-14: twice in one session a verdict was formed from an 8-day-stale synthesis memo — an already-proven theorem was reported as a future sprint, and an already-settled A-vs-B question was re-answered wrong; both answers were live in Paper 56. Standing rule: `memory/feedback_verify_current_state.md`.)

### Benchmarking Rule

After any modification to production code in `geovac/`:
1. Run `/regression` (default scope `touched`) — derives the test selection from `git diff` + import graph (consumer test files of every touched module) plus the topological-integrity baseline plus a small reproducible random sample. 30s–2min typical wall.
2. If the diff spans more than 2–3 modules, prefer `/regression full` — the 10–15 min full pass. Catches cross-cutting consumers that the diff-derived selection might miss when a refactor cascades across the codebase.
3. Verify the 18 symbolic S³ proofs pass (always included in `touched` and `topo` scopes).
4. Verify H2+ < 0.1% error (topological control).
5. Verify H2 Full CI < 1.0% error (accuracy control).
6. Report any speed regression > 10%.

The previous narrow 3-file allowlist (`tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`, `tests/advanced_benchmarks.py`) silently let test rot accumulate when refactors cascaded into the dozens of test files that import from `composed_qubit`, `inter_fiber_coupling`, etc.  `/regression touched` removes the consumer-selection bottleneck by deriving it mechanically from the diff.  Use it as the standard discipline after any code edit.

### Clean Room Rule

- Generated plots -> `debug/plots/` or `papers/figures/`
- Generated data -> `debug/data/`
- Scripts -> `debug/`, `demo/`, `tests/`, or `benchmarks/` (never root)
- Documentation -> `docs/`

### Changelog Protocol

**Changelog granularity:** Each CHANGELOG entry should correspond to a releasable unit of work — a new feature, a completed diagnostic arc, a paper update, or a benchmark result. Do not create a new version entry for intermediate debugging steps, parameter sweeps, or exploratory runs. If a task requires multiple iterations to complete, document it as a single entry when the task concludes, noting the key findings and negative results. Intermediate data belongs in `debug/` with timestamped filenames, not in the CHANGELOG.

**Version numbering:** Patch versions (x.y.Z) for bug fixes and documentation; minor versions (x.Y.0) for new features, completed diagnostic arcs, or paper additions; major versions (X.0.0) for architectural changes. A diagnostic arc that tests 10 hypotheses and finds 9 negative results is one minor version, not 10 patches.

---

## 10. Validation Benchmarks

See `docs/validation_benchmarks.md` for the full benchmark table. Update that file (not this section) when adding new benchmarks.

<!-- Full table extracted to docs/validation_benchmarks.md on 2026-05-31 to reduce context-load cost. -->

---

## 11. Topic-to-Paper Lookup

See `docs/topic_to_paper_lookup.md` for the full table mapping topics to papers, sections, and loading tiers. Update that file (not this section) when adding new topic→paper mappings.

<!-- Full table extracted to docs/topic_to_paper_lookup.md on 2026-05-31 to reduce context-load cost. -->

---

## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. **Full registry** (Levels 2 / 3 / 4 / 4N / 5 + spin-ful Tier-2 tables): `docs/algebraic_registry.md` (live document — update statuses there, not here). Status vocabulary: **algebraic** (closed-form from quantum numbers) / **algebraic (implicit)** (defined by P = 0 with known coefficient ring; pointwise diagonalization is convenience, not necessity) / **algebraic-pending** (route identified, production still uses quadrature) / **numerical-required** (no known algebraic replacement). The §4 prime-directive test governs changes: anything touching quantum-number labels or selection rules is prohibited; improving radial-amplitude evaluation within a channel is legitimate.

---
## 13. Multi-Agent Protocol

The GeoVac project uses an AI-augmented agentic workflow with a formalized four-layer architecture for research direction, planning, and execution.

### 13.1 Architecture

Four layers:

1. **Research agents (strategic layer):** Four specialized agents in `agents/` that formalize the strategic planning process. The PI invokes these by loading the agent file + CLAUDE.md into a Claude session. See `agents/WORKFLOW.md` for the full protocol.

   | Agent | File | Role |
   |:------|:-----|:-----|
   | Leader | `agents/LEADER.md` | Synthesizes project state, proposes ranked research directions with reasoning. Reads CLAUDE.md + SCOPE_BOUNDARY.md + recent results. |
   | Explorer | `agents/EXPLORER.md` | Searches published literature for geometric structures that fit the natural geometry principle. Returns structured candidate reports. |
   | Decomposer | `agents/DECOMPOSER.md` | Takes mathematical expressions/code paths, catalogs transcendental content using the exchange constant taxonomy (Paper 18), proposes algebraic separations. |
   | Reviewer | `agents/REVIEWER.md` | Critiques paper drafts against GeoVac standards (rhetoric rule, benchmarking rule, transcendental cataloging) and general scientific rigor. |

   The research agents operate on the three guiding principles:
   - **Natural geometry search:** If in doubt, search for new geometries to incorporate (Explorer).
   - **Algebraic deconstruction:** When code runs long, approach algebraically. Diagonalize. Deconstruct the continuum (Decomposer).
   - **Transcendental cataloging:** Catalogue where transcendental quantities enter and their relationship to geometries (Decomposer).

   The PI remains the approval gate: agents propose, the PI decides.

2. **Plan mode (human + agents → sprint plan):** The PI reviews the Leader's Strategic Brief, selects a direction, and defines the sprint plan (track definitions, PM prompts, exit criteria). This layer produces the directive that the PM agent executes. CLAUDE.md changes originate here.

3. **PM agent (main Claude Code session, full context):** Reads CLAUDE.md, README, and all papers relevant to the current track. **Default: does the work directly in main session — coding, computation, paper edits, commits.** Sub-agents are reached for only when the PI explicitly directs dispatch OR when the work is genuinely parallelizable AND context-heavy enough that main-session would clog. Evaluates sub-agent results against verification checklists when they are used.

4. **Worker sub-agents (scoped context, opt-in):** Available for cases where parallel dispatch genuinely beats sequential main-session work (e.g. 4+ independent computational tracks, each context-heavy). Default is to NOT use them. Each sub-agent reloads CLAUDE.md and its task-relevant files, which is expensive; reserve for cases where the cost is justified by the parallelism or context-protection benefit.

**Sub-agent cost discipline (2026-05-26 policy update).** Earlier sessions defaulted to sub-agent dispatch for every sprint. The cost (CLAUDE.md re-loaded per dispatch, parallel multiplier on context budget) has been substantial, and many sprints have been done equally well or better in main session. Default is now flipped:\ main-session by default, sub-agent on PI direction. The research agents (Leader/Explorer/Decomposer/Reviewer) in `agents/` remain available for big strategic moments; they are demoted from default to opt-in.

**The research loop:**
```
Leader Brief (when invoked) → PI picks direction → PM works in main session →
PI directs sub-agent dispatch IF needed → Reviewer critiques papers (when invoked) →
Results feed back to next Leader Brief
```

### 13.2 PM Session Kickoff

Every PM session begins by reading CLAUDE.md and then executing the following:

1. Identify the current track(s) and relevant papers from the plan mode directive
2. Read those papers and any results from the previous session
3. Check the failed approaches summary (Section 3); if the current track touches PK, cusp, inter-group antisymmetry, or molecular encoding, read the full details in CHANGELOG.md before proceeding
4. Plan the session as main-session work (sequence of Read / Edit / Bash / etc.). Flag any tracks that are candidates for sub-agent dispatch (parallelizable, context-heavy) and ask the PI before launching
5. Identify which papers need updating based on the session's results (see 13.8)

**Sprint standard (2026-05-26 update):** The default workflow is **main-session work**. The PM reads, edits, computes, drafts memos, and commits directly. Sub-agent dispatch is opt-in:\ either the PI explicitly directs ("dispatch this in parallel," "run Explore for this lookup") or the PM identifies a case where dispatch genuinely beats sequential main-session work and asks first. When sub-agents are used, the standard prompt template (§13.3) applies.

Earlier convention was "one PM prompt per sprint, dispatching all tracks as parallel sub-agents." That convention is retired:\ it overspent context budget on tasks that could be done sequentially in main session at lower cost.

### 13.3 Sub-Agent Prompt Template (opt-in)

Used only when the PI directs sub-agent dispatch (see §13.2 update). Sub-agents read CLAUDE.md by default. Do not restate CLAUDE.md context inside the prompt. Keep prompts terse: target < 1500 words for diagnostics, < 2500 for implementation. The format below is the canonical shape; omit any line that's obvious from context.

```
TASK: [one sentence, the deliverable]
DECISION GATE: [what counts as GO / BORDERLINE / STOP, with thresholds]
FILES TO READ (beyond CLAUDE.md): [only files the agent wouldn't naturally find]
DO NOT MODIFY: [only if a non-obvious file is at risk]
OUTPUT:
  - [what to return in the response, brief]
  - [what files to write: debug/*.py drivers, debug/data/*.json, debug/*_memo.md]
  - [paper edits to apply, if any — apply directly per §13.8]
```

Constraints (failed approaches to avoid, structures to preserve, success criteria) are stated only when non-obvious from CLAUDE.md §3 and the paper guardrails. The agent has access to all of CLAUDE.md; trust it to find what it needs.

### 13.4 Verification Gates

Before the PM agent accepts a sub-agent result, it checks:

1. **Test gate:** Do all relevant tests pass? (Non-negotiable.)
2. **Dead-end gate:** Does the approach match any entry in the failed approaches table? If so, reject unless the sub-agent explicitly explains what is different this time.
3. **Prime directive gate:** Does the result modify any discrete structure — quantum number labeling, selection rules, channel structure, Gaunt integral coupling? If so, do NOT accept. Escalate to plan mode for human review.
4. **Consistency gate:** Does the result contradict any established result in the papers? If uncertain, flag in the session summary rather than accepting.
5. **Equation gate:** Does every equation in affected papers have a corresponding test in `tests/`? If a new equation was added to a paper in this session, was a verification test also added? If not, the session summary must flag this as an open item.

### 13.4a Equation Verification Protocol

Every equation that appears in a GeoVac paper must have a corresponding numerical verification in the codebase. This protocol bridges the research agents (which propose mathematical structures) and the PM/worker pipeline (which implements and tests them).

**The rule:** No equation goes into a paper without a test that verifies it computationally. "Verified" means one of:

| Verification type | What it checks | Example |
|:------------------|:---------------|:--------|
| **Analytical limit** | Equation reduces to a known result in a limiting case | H₂ → 2×H as R→∞ |
| **Symbolic identity** | Equation matches an equivalent form term-by-term | Gaunt integral = product of 3j symbols (symbolic, exact) |
| **Numerical cross-check** | Equation's output matches an independent computation | GeoVac He energy vs NIST reference |
| **Dimensional consistency** | All terms have consistent units | Energy terms in Ha, not mixed Ha/eV |
| **Symmetry verification** | Claimed symmetries hold numerically | Hermiticity of Hamiltonian matrix (H = H†) |

**Implementation:**

1. **When a sub-agent derives a new equation:** The PM must dispatch a verification sub-agent that implements the equation in code and checks it against at least one of the verification types above. The test goes into `tests/` alongside the production code.

2. **When the Decomposer proposes an algebraic separation:** The PM must verify that the separated form reproduces the original to machine precision (< 1e-12 relative error). Both the original and separated forms must be implemented and compared.

3. **When the Reviewer flags an unverified equation:** The PM treats this as a blocking issue — the equation is removed from the paper or a verification track is opened before the paper is finalized.

4. **Test naming convention:** Equation verification tests are named `test_paper{N}_eq{M}` or `test_paper{N}_{description}` so they can be traced back to the specific paper claim they verify.

**What counts as sufficient verification depends on the claim:**

- **Exact identities** (e.g., "the graph eigenvalue is n²-1"): must be verified symbolically or to machine precision across all relevant quantum numbers up to n_max=5.
- **Scaling laws** (e.g., "Pauli count scales as O(Q^2.5)"): must be verified across at least 4 data points with a log-log fit. Report the fit residual.
- **Accuracy claims** (e.g., "0.004% error"): must report the computed value, the reference value, the reference source, and the actual percentage to enough digits to confirm the claim.
- **Negative results** (e.g., "this approach diverges"): must show the divergence numerically (e.g., error increasing monotonically with basis size across at least 3 points).

### 13.5 Hard Prohibitions

The following changes must NEVER be made by sub-agents or the PM agent:

- Any change to the natural geometry hierarchy (new levels, changed coordinate systems)
- Introduction of any fitted or empirical parameter
- Deletion or suppression of negative results from Section 3 or CHANGELOG.md
- Removal of the "conjectural" label from the **combination rule K = π(B + F − Δ) in Paper 2**. (Paper 2 was moved Conjectures → Core on 2026-04-18 after Sprint A, then moved Core → Observations on 2026-05-02 per the curve-fit audit memo `docs/curve_fit_audit_memo.md`. The paper's surrounding structural decompositions — three independent spectral homes for B, F, Δ — are derived; the combination rule itself is a numerical observation without first-principles derivation, and the "conjectural" framing applies to the combination rule whether or not the paper's folder location changes.)

**CLAUDE.md access control:** The PM may edit CLAUDE.md for mechanical updates that keep documentation in sync with code and paper changes. The PM may NOT edit sections that define strategy, framing, or the PM's own operating rules.

| Section | PM may edit? | Examples of allowed edits |
|:--------|:-------------|:-------------------------|
| 1 (Project Identity) | Version number ONLY | Bump v2.0.25 → v2.0.26 |
| 1.5 (Positioning & Framing) | **NO** | — |
| 1.6 (Project Phase) | **NO** | — |
| 2 (Development Frontier) | Yes | Update best results, add/complete track summaries, update backlog |
| 3 (Failed Approaches) | Yes (append only) | Add new failed approach rows; never delete or modify existing entries |
| 4 (Dimensionless Vacuum) | **NO** | — |
| 5 (Natural Geometry Hierarchy) | **NO** | — |
| 6 (Paper Series) | Yes | Update file paths, loading guide descriptions, inventory tables, key results |
| 7-9 (Code/Coding/Workflow) | Yes | Add new entry points, update module paths |
| 10 (Validation Benchmarks) | Yes | Add new benchmark rows for new tests |
| 11 (Topic-to-Paper Lookup) | Yes | Add new topic → paper mappings |
| 12 (Algebraic Registry) | Yes | Update status (algebraic-pending → algebraic) when proven |
| 13 (Multi-Agent Protocol) | **NO** | — |
| 14 (Test Architecture) | **NO** | — |

### 13.6 Track Management

Active work is organized into tracks. Each track has:

- A name and one-sentence goal
- A list of relevant papers and code modules
- A current status (active / blocked / complete)
- A log of sub-agent dispatches and results (maintained by the PM in `debug/track_logs/`)

The PM agent maintains a brief track status file at `debug/track_logs/STATUS.md` that is updated at the end of each session. This file is read at the start of the next PM session to restore context.

### 13.7 General Guidance

**All changes are autonomous.** The test suite catches code errors. The PI catches framing errors in review. The PM's job is to make its best judgment, apply changes directly (including to papers), and produce a clear summary of what changed and what the results were.

The overriding principles are: **code changes are autonomous** (the test suite catches errors); **paper updates are autonomous** (adding results, tables, new-method subsections, correcting claims contradicted by new evidence, reframing based on findings). If a change is wrong, it will be caught in plan-mode review and corrected — this is faster than proposal cycles.

### 13.8 Paper Update Policy

Papers are the authoritative source for all physics (Section 1). Code that outpaces the papers creates documentation drift. PMs are expected to keep papers in sync with code results.

#### Paper edit policy

PMs may edit papers in any of the six group folders (`papers/group1_operator_algebras/`, `papers/group2_quantum_chemistry/`, `papers/group3_foundations/`, `papers/group4_quantum_computing/`, `papers/group5_qed_gauge/`, `papers/group6_precision_observations/`), the synthesis folder (`papers/synthesis/`), and the archive folder (`papers/archive/`) directly for ALL changes, including:

- Adding or updating benchmark tables and numerical results
- Adding new subsections documenting methods and results
- Correcting claims that are contradicted by new computational evidence
- Reframing results based on new findings (e.g., reclassifying transcendental → algebraic when proven)
- Updating abstracts and conclusions to reflect current best results

**PMs must still NOT:**
- Introduce fitted or empirical parameters without PI direction
- Change the natural geometry hierarchy (new levels, changed coordinates)
- Delete or suppress negative results from Section 3
- Remove the "conjectural" label from the **combination rule K = π(B + F − Δ) in Paper 2**. The prohibition is at the combination-rule level, not the paper-tier level: it applies regardless of whether Paper 2 sits in Conjectures, Core, or Observations. (History: Conjectures → Core 2026-04-18 (Sprint A), Core → Observations 2026-05-02 (curve-fit audit memo).)

**Splinter file prohibition:** PMs must edit papers in-place. Do NOT create separate .tex diff files, proposal files, or draft directories. Proposed changes go directly into the paper. If the change is wrong, `git revert` is cheaper than context-loading splinter files in plan-mode review.

### 13.9 Session Summary Format (MANDATORY)

Every PM session MUST end with a summary in the following format. This is the PI's primary review mechanism — without it, changes are invisible.

```
## Session Summary [date]

### Tracks
- Track XX: [status] — [one-line result]

### Results
[Tables of computed data, benchmark comparisons, scaling exponents,
convergence studies — whatever the track produced. Include numbers,
not just descriptions. The PI needs to see the data to evaluate
whether the result is correct.]

### Files Modified
- `path/to/file` — [one-line description of what changed]

### Files Created
- `path/to/file` — [one-line description]

### Decisions
- [Any changes to paper claims, framing, or theoretical arguments — briefly stated]
- [Any negative results or dead ends encountered]
```

The Results section goes before the file list because the data is what the PI reviews first. If a track produced no quantitative results (e.g., a pure documentation track), the Results section can be replaced with a brief description of what changed and why.

### 13.9a Two primitives for PI-desired behavior

Two distinct primitives for shaping how the PM behaves. They are NOT interchangeable, and using the wrong one is a frequent failure mode (caught and corrected 2026-05-26 in the v3.3.1 release):

**Standing rules — `memory/feedback_*.md`.** Behaviors the PM runs automatically whenever the trigger condition is met. The PI should NOT have to invoke these. Memory rules are loaded into every conversation and apply unconditionally. The body of each rule states (a) the trigger condition, (b) why the rule exists (incident or insight that produced it), and (c) how to apply it. Examples: `feedback_audit_numerical_claims.md` (any "X matches Y" claim triggers curve-fit audit); `feedback_diagnostic_before_engineering.md` (≥ 2 honest negatives on a wall triggers diagnostic-only sprint); `feedback_tag_transcendentals.md` (any transcendental appearance triggers Paper 18 + Mellin engine + Paper 34 classification).

**Slash commands — `.claude/commands/<name>.md`.** Actions the PI triggers at a specific moment chosen by the PI. Slash commands are NOT defaults — they are explicit invocations. Each command is a saved prompt that fires in current main-session context (no sub-agent dispatch, no CLAUDE.md re-load, no context multiplier).

**The distinction.** If the behavior should fire whenever a condition is met → memory rule. If the behavior fires at a specific PI-chosen moment → slash command. A behavior CAN be both:\ the memory rule is the primary mechanism (enforces the default), and the slash command is a **force-fire backup** for cases where the PM has failed to detect the trigger condition. **Putting a behavior behind a slash command alone makes it more optional, not less** — it shifts the responsibility for triggering to the PI, which is the opposite of what a discipline rule needs.

### 13.9b Project slash commands

Six slash commands defined in `.claude/commands/`:

| Command | Type | Purpose | When to fire |
|:--------|:-----|:--------|:-------------|
| `/ahha` | Trigger | Surface three unsurfaced connections from recent work; distinguish trained-flatness hedge from real epistemic caution. | PI-chosen moment after a substantive result lands. Forces the kind of reach-for-the-connection pass the PM does not do by default. |
| `/sprint-close` | Trigger | End-of-sprint protocol:\ canonical memo, CHANGELOG entry, CLAUDE.md §2 one-liner, optional §3 row and paper edits, optional MEMORY index entry, verification check, honest-scope check. | When a sprint completes and is ready for release. Standardises the close pattern. |
| `/release` | Trigger | Version bump + commit + tag + push wrapper, with precondition checks (version string bumped, §2 entry exists, papers compile clean, tests pass, hard-prohibition check). | After `/sprint-close` has staged the content; mechanical release step. |
| `/audit-claim` | **Force-fire backup** for [[feedback_audit_numerical_claims]] | Curve-fit-audit pattern from `docs/curve_fit_audit_memo.md`. The standing rule should already fire on any "X matches Y" claim; this command is a backup if it didn't. | When the PI notices the PM made a numerical-coincidence claim without running the audit. |
| `/diag` | **Force-fire backup** for [[feedback_diagnostic_before_engineering]] | Diagnostic-before-engineering pass. The standing rule should fire automatically when ≥ 2 honest negatives accumulate; this command is a backup. | When the PI notices the PM is about to launch an implementation sprint past the 2-negatives threshold without running the diagnostic. |
| `/transcendental-tag` | **Force-fire backup** for [[feedback_tag_transcendentals]] | Paper 18 + master Mellin engine + Paper 34 classification. The standing rule should fire automatically on any transcendental appearance; this command is a backup. | When the PI notices a transcendental was introduced into production code or a paper without classification. |

The first three are genuine triggers (no corresponding memory rule, fired by PI at a specific moment). The last three are force-fire backups (their primary mechanism is a memory rule that should fire automatically; the slash command exists to manually invoke the discipline when the rule failed to trigger).

**Adding more.** Memory rule:\ write `memory/feedback_<name>.md` with the standard frontmatter + body, and add an index line to `MEMORY.md` (≤ 200 chars). Slash command:\ write `.claude/commands/<name>.md` with a single `description:` frontmatter line and the prompt body. No registration, no schema, no config.

### 13.10 Research Agent Integration

**When to invoke research agents:**

| Situation | Agent | What to provide |
|:----------|:------|:----------------|
| Starting a new research phase or feeling stuck | Leader | CLAUDE.md + SCOPE_BOUNDARY.md + any recent results |
| "What geometric structure handles X?" | Explorer | The specific question + relevant papers |
| "Where do transcendentals enter this computation?" | Decomposer | The expression/code path + Paper 18 |
| Paper draft ready for critique | Reviewer | The draft + CLAUDE.md + prior papers it builds on |
| After completing a track (positive or negative) | Leader | Updated CLAUDE.md with new results |

**Handoff from research agents to PM:**

When a research agent produces a result (e.g., the Explorer identifies a candidate geometry, or the Decomposer proposes an algebraic separation), the PI translates it into a sprint plan for the PM using the standard track format:

```
Track [XX]: [Agent output → implementation goal]
  Agent source: [which agent, what it proposed]
  PM prompt: [specific implementation task]
  Exit criteria: [what constitutes success/failure]
  Papers affected: [which papers would be updated]
  Equation verification: [which new equations need tests]
```

**Research agents do not write code or modify files.** They produce analysis, proposals, and critiques that the PI evaluates and the PM/worker pipeline implements. This separation is intentional: the agents explore, the PM executes, the tests verify.

### 13.11 Content Discipline and Token Efficiency

CLAUDE.md is loaded into every PM session and every sub-agent dispatch. Its size is paid repeatedly and every edit invalidates the prompt cache. The rules below keep this cost bounded. They are hard rules, not preferences — they override the impulse to "be thorough" in CLAUDE.md, because thoroughness has homes other than CLAUDE.md.

**Where each kind of content lives:**

| Content | Home | Format |
|:--------|:-----|:-------|
| Sprint chronicle (what happened, in detail) | `CHANGELOG.md` | Full prose, no length limit |
| Sprint summary (CLAUDE.md §2 entry) | CLAUDE.md §2 | ≤ 30 words: name + date + verdict + memo path |
| Dead-end record (§3) | CLAUDE.md §3 | ≤ 2 sentences + memo link |
| Sprint memo | `debug/*.md` | ONE canonical memo per sprint, ≤ 5000 words |
| Surviving structural findings | papers/group*/*.tex | The papers, not CLAUDE.md, are the permanent record |
| Cross-session facts | `memory/*.md` | One-liner index entry in MEMORY.md, ≤ 200 chars |
| Active behavior rules | CLAUDE.md §13 | Short, dense, this section |

**Rules:**

1. **No synthesis memos.** One canonical memo per sprint. Cross-sprint synthesis lives in CHANGELOG.md or in a paper. The "comprehensive synthesis memo that supersedes earlier synthesis memos" pattern is forbidden.

2. **CLAUDE.md §2 entries are one-liners.** Format: `**Sprint NAME (YYYY-MM-DD):** Verdict in one sentence. See debug/MEMO.md.` Full sprint detail goes to CHANGELOG.md (`### Added` / `### Changed` per release entry). Existing multi-thousand-word §2 entries are technical debt to be compacted as touched.

3. **CLAUDE.md §3 dead-end rows are short.** Format: `Approach name (date) | count | One or two sentences: the lesson + path to memo.` Existing paragraph-length rows are technical debt.

4. **Memory files are for cross-session facts not derivable from CLAUDE.md or papers.** Do not auto-create memory files for sprint outcomes (sprint detail lives in CHANGELOG.md). MEMORY.md index entries strictly ≤ 200 chars — the system silently truncates beyond ~24KB and the warning has fired multiple times.

5. **Agent prompts: task + decision gate + specific files + output format.** Do not restate CLAUDE.md context inside agent prompts — sub-agents already load CLAUDE.md. Target: < 1500 words for diagnostic tasks, < 2500 words for implementation tasks. If the prompt is growing past 2500 words, the task is probably too big for one agent.

6. **Prefer Explore agent for read-only diagnostics.** Reserve general-purpose for tasks that require code modification or extensive computation. Explore has narrower context allocation.

7. **Don't dispatch sub-agents for tasks doable in the main session.** A single Read + Edit, a few-line script, or analysis under ~500 lines of context should be done by the PM directly. Sub-agents are for parallelizable work or context-heavy delegation.

8. **One canonical record per fact.** A fact appears in at most one of: CLAUDE.md, CHANGELOG.md, paper, memo, memory. Duplication is the failure mode that produced the current bloat.

9. **Status updates replace, never append (added 2026-06-10, PI-authorized).** When a WH status, a §6 status flag, or a paper-state description changes, REPLACE the existing text and move the superseded version to its history home (`docs/wh_register_history.md`, `docs/paper_notes_archive.md`, CHANGELOG.md). Chronicling-by-appending inside CLAUDE.md is the failure mode that regrew the file from 1,263 lines (2026-05-31 compaction) to 1,400 lines / 320 KB by 2026-06-10.

**Enforcement.** When the PM is tempted to write a long CLAUDE.md §2 bullet, a synthesis memo, or a verbose agent prompt — stop, move the content to its proper home, and write the short version in CLAUDE.md. Apply this to existing entries when touched, not as a one-time pass.


---

## 14. Test Architecture Policy

### Main pipeline tests (`tests/`)

The `tests/` directory contains ONLY tests that validate the main pipeline — the methods and results described in the active papers (Papers 0–18, FCI-A, FCI-M). Every test file in `tests/` should correspond to at least one paper result or a piece of infrastructure that the pipeline depends on.

**The rule:** If a paper documents a method as its primary result, that method's tests live in `tests/`. If a method is documented as failed (Section 3), superseded (architecture locked), or was a one-off diagnostic, its tests live in `tests/_archive/`.

### Archived tests (`tests/_archive/`)

Archived tests are organized into three tiers:

| Directory | Contents | When to use |
|-----------|----------|-------------|
| `tests/_archive/dead_ends/` | Tests for approaches documented as failed in Section 3 | Method confirmed dead; resolution documented |
| `tests/_archive/superseded/` | Tests for approaches replaced by a newer architecture | New method has its own tests in `tests/`; old method explicitly "architecture locked" |
| `tests/_archive/diagnostics/` | Tests from one-off investigation arcs | Diagnostic complete; resolution baked into main pipeline |

Archived tests are **not** collected by default (`--ignore=tests/_archive` in pytest config). They can be run explicitly via `pytest tests/_archive/`.

Archived tests have scientific value as institutional memory — they help prevent re-deriving failed approaches. **Do not delete archived tests.**

### Archived code (`geovac/_archive/`)

Follows the same three-tier structure (`dead_ends/`, `superseded/`, `auxiliary/`). Archived code modules are still importable from their `geovac._archive.*` paths. If a core module previously imported an archived module, a compatibility shim is left at the original location.

### When new work produces dead ends

When a track produces a negative result or a method is superseded:

1. Document the failure or supersession in CLAUDE.md Section 3 (failed approaches) or Section 2 (active frontier, marking the track COMPLETE)
2. Move the relevant tests to `tests/_archive/` in the same commit
3. Move the relevant code module to `geovac/_archive/` in the same commit
4. **Check for redirects first:** If any tests validate infrastructure that's still live (CI machinery, solver constructors, shared utilities), extract those tests and rewrite them to call the current architecture before archiving the rest
5. Update import paths in archived files

### Slow tests

Core tests that are computationally expensive are marked with `@pytest.mark.slow` and skipped by default. Run them with `pytest --slow`. A test qualifies as slow if it takes >10 seconds or performs heavy computation not needed for fast regression (e.g., cc-pVTZ integral engine, full VQE pipeline).

### Redirect-before-archive rule

Before archiving a test file, the PM must check: does this file test any infrastructure that the main pipeline still depends on? If yes, extract those tests into a core test file that calls the current API, then archive the rest. This prevents coverage gaps from forming when methods are superseded.

### Post-archive test commands

| What You Want | Command |
|---------------|---------|
| Fast regression (daily dev) | `pytest` |
| Full core + slow tests | `pytest --slow` |
| Archived tests only | `pytest tests/_archive/` |
| Everything | `pytest tests/ tests/_archive/ --slow` |
| Specific archive tier | `pytest tests/_archive/dead_ends/` |
| Topological proofs only | `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py` |
