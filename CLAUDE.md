# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v5.21.4 (October 21, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The GeoVac lattice is the finite label set (n,l,m) of the S³ hyperspherical harmonics, on which the S³ Laplace–Beltrami operator acts exactly as the diagonal SO(4) Casimir n²−1 at every cutoff; Fock's 1935 conformal projection makes that operator equivalent to the Schrödinger equation. (The lattice's edge Laplacian D−A preserves l, splits into n_max disconnected blocks, and does NOT converge to Δ_{S³} — a documented wall, 2026-09-25, `tests/test_graph_s3_nonconvergence.py`.) This structure is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Mission statement (adopted 2026-08-29, PI direction):** GeoVac charts the forced/free boundary of quantum physics. For every structure -- quantum number, selection rule, degeneracy, sparsity pattern, convergence rate, physical constant -- the program renders one of three verdicts: **FORCED** (derived from the packing construction, exactly, with a frozen falsifier), **FREE** (an exchange constant, with its projection chain named and its minimal transcendental content classified), or **WALL** (a proven obstruction with the mechanism pinned). The deliverable is the atlas of that boundary. Under this statement the corpus is one program: the chemistry/QC arc surveys how far the forced side reaches computationally; the 40+ documented negatives are the boundary itself, measured; the precision program (SS1.8) is the atlas's experimental interface; the periods/transcendence work (Papers 18/34/54-59) is the coordinate system for the free side; and the QA apparatus is what makes the atlas trustworthy. *Scope note:* this is the internal research mission (the register of SS1.7); papers remain under the SS1.5 rhetoric rule -- "forced" is atlas vocabulary, not ontology language for publication.

**Authoritative source rule:** The papers in `papers/group1_operator_algebras/`, `papers/group2_quantum_chemistry/`, `papers/group3_foundations/`, `papers/group4_quantum_computing/`, `papers/group5_qed_gauge/`, `papers/group6_precision_observations/`, and `papers/synthesis/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it. (Papers were reorganized from the previous `core/`, `methods/`, `applications/`, `synthesis/`, `standalone/`, `observations/`, `conjectures/` layout into six audience-targeted groups on 2026-05-22.)

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator (Josh Loutey — `J.~Loutey` on the papers) provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers across the six audience-targeted group folders (see §6) are written to academic standards but are not submitted to traditional journals. The project's viability case rests on its corpus of verified structural results — discretization- and encoding-structure theorems, honest negative results, and the benchmarked computational artifact (pip-installable, 37-system library) that demonstrates them. The tool is a research instrument, not a production-chemistry replacement (see the benchmarking rule in §1.5 and `docs/claims_register.md`). Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks. (Viability-case sentence reworded 2026-06-10 per PI direction, accessibility plan Phase 3.)

---

## 1.5. Positioning & Framing

GeoVac is a discretization framework that exploits the natural geometry of separable quantum systems to produce sparse graph Hamiltonians. It is **not** a new foundation for quantum mechanics — it works because its discrete labels, degeneracies and selection rules are exactly those of the known continuous operators (for S³, the Laplace–Beltrami operator is the diagonal SO(4) Casimir on the labels), with Fock's continuum projection carrying the equivalence (Paper 7). It does NOT work by a graph Laplacian converging to them: the S³ edge Laplacian D−A provably does not (2026-09-25, PI-directed re-pricing).

**Rhetoric rule:** The GeoVac framework is demonstrably conformally equivalent to standard continuous quantum mechanics through Fock's continuum projection (Paper 7: 18 symbolic checks of the continuum chain for the S³ case); the discrete side of the duality is the exact label set with its selection rules and the diagonal SO(4) Casimir, not a graph Laplacian converging to Δ_{S³} (a documented wall). The papers should present this conformal equivalence as the primary result. The discrete label lattice and the continuous Laplace–Beltrami operators are dual descriptions of the same quantum-number structure — the mathematics supports both readings, and the interpretation of which is more fundamental is left as a choice for the reader. Avoid language that asserts ontological priority of either description. Lead with concrete computational results (sparsity, scaling, accuracy, structural insight) rather than interpretive claims about the nature of quantum mechanics.

**Lead with concrete advantages:** O(V) sparsity, angular momentum selection rules baked into the basis, zero-parameter construction from nuclear charges and geometry alone, efficient qubit encodings (Paper 14). These structural properties are the framework's actual selling points — not philosophical claims about the nature of quantum mechanics.

**The π-free graph principle:** The GeoVac graph is π-free: its matrix entries, degeneracies, and coupling coefficients are integers or rationals, so its eigenvalues are algebraic numbers (roots of integer polynomials — for the S³ graph Laplacian, sums of 2−2cos(jπ/m), Paper 0 §VI), never transcendental. (The pure-integer ladder −(n²−1) is the *continuum* S³ Laplace–Beltrami spectrum, not the graph's; corrected 2026-09-25, PI-directed.) Transcendental numbers (π, exponential integrals, spectral zeta values) enter exclusively when projecting the graph onto continuous manifolds for comparison with experiment or for computational convenience. This observation (Paper 18) motivates the design principle: stay on the graph whenever possible, and when projection is necessary, identify the minimal transcendental content (the exchange constant) required. The exchange constant taxonomy (intrinsic, calibration, embedding, flow) classifies projections by what determines them and predicts the computational cost of each departure from the graph.

**Quantum simulation positioning:** The classical solver investigation (v2.0.6-23) characterized what the natural-geometry hierarchy can and cannot achieve for ground-state PES; the primary computational value proposition is quantum simulation. The composed architecture produces qubit Hamiltonians whose Pauli count is EXACTLY linear in Q across molecules at fixed basis (main-group and d-block coefficients: numeric registry `composed_coeff`, `composed_coeff_d`), against the steeper polynomial system-size scaling of Gaussian bases; a large equal-qubit Pauli advantage over Gaussian baselines for H2O at n_max = 1-3, identity-included on both sides (Paper 14); structural Gaunt sparsity compatible with all downstream optimizations (tapering, grouping, tensor factorization); and a 37-system library (35 composed molecules + He + H2) across three periodic-table rows. The within-molecule basis exponent (registry `exp_composed_3pt`) is a small-basis fit whose local slope rises with n_max — the linearity in Q, not the basis exponent, is the load-bearing advantage. Position for VQE/NISQ (Pauli count, QWC groups) rather than fault-tolerant QPE: composed LiH is near parity with minimal-basis STO-3G on raw Pauli counts, uses more terms than STO-3G's reduced count, and matches its 1-norm (registry `lih_composed_pauli`, `lih_sto3g_raw`, `lih_composed_lambda`, `lih_sto3g_lambda`); the decisive wins are the exact Q-linearity and the large reduction against cc-pVDZ (`lih_ccpvdz_pauli`). All counts are under the exact global-M_L rule (2026-08-29); the figures measured under the retired pair-diagonal sign error are blocked by C17 (debug/sprint_eri_evaluator_defects_memo.md). Full track-by-track chronicle: `docs/development_frontier_archive.md` + CHANGELOG.

**General composed builder foundation:** Paper 16 (Chemical Periodicity as S_N Representation Theory) provides the group-theoretic foundation for the general composed builder's atomic classification. The structure types (A/B/C/D/E), the minimum angular quantum number ν_min = E(n↑)+E(n↓) (the Pauli-filling count; the earlier "universal ν = N−2" holds only for the lowest spin with N ≤ 8 and was retracted 2026-09-28, PI-authorized), and the recursive core-valence decomposition are all derived from S_N representation theory and implemented in `atomic_classifier.py`.

**PK limitation and research path:** The composed framework's PK pseudopotential is the accuracy bottleneck (composed R_eq errors from several percent up to about twenty percent; Paper 17); six modification attempts failed (Section 3). The balanced coupled Hamiltonian (Papers 19/20) is the PK-free alternative: cross-center V_ne via multipole expansion with exact Gaunt termination; its LiH 4e FCI energy converges to sub-percent accuracy at n_max=3 while its R_eq drifts structurally (Paper 19) — energy converges excellently, geometry drifts; optimal for single-point quantum simulation at fixed geometries. Resource anchor: composed LiH vs STO-3G vs cc-pVDZ Pauli counts (registry `lih_composed_pauli`, `lih_sto3g_raw`, `lih_ccpvdz_pauli`; exact rule, 2026-08-29). Chronicle: `docs/development_frontier_archive.md`.

**Benchmarking rule:** When comparing to other methods, always use the strongest available baseline (cc-pVTZ or better for atoms, explicitly correlated methods for molecules), not just STO-3G. If the comparison is unfavorable, say so honestly and identify what the framework offers instead (sparsity, scaling, structural insight). Position the framework as a computationally principled alternative, not as a replacement for production quantum chemistry.

---

## 1.6. Project Phase

**Phase 1 (v0.9.x-v1.x): Foundation.** Graph Laplacian equivalence proof, natural geometry hierarchy, LCAO diagnostic arc, bond sphere theory, 18 symbolic proofs (Paper 7).

**Phase 2 (v2.0.0-v2.0.23): Classical solver investigation.** Systematic exploration of all solver architectures (adiabatic, coupled-channel, 2D variational) across all levels (2-5, 4N). 30+ completed tracks, 40+ documented negative results. Key outcomes: H2 D_e via the Level-4 hyperspherical treatment (Paper 15), LiH R_eq via Level 5 composed (Paper 17), full N-electron equilibrium without PK (Level 4N). Investigation complete: all solver x PK x basis combinations exhausted, structural accuracy ceilings characterized, exchange constant taxonomy established (Paper 18).

**Phase 3 (v2.0.24-v2.7.0): Quantum simulation.** The composed architecture's structural sparsity (exact Q-linearity of the Pauli count across molecules at fixed basis; block-diagonal ERIs) is the framework's primary computational advantage. (Corrected 2026-08-29: the retired pair-diagonal rule had named the *within-molecule* basis exponent as the advantage; under the exact rule it is the across-molecule linearity -- not the basis exponent -- that is load-bearing.) Classical solver results serve as validation benchmarks for the quantum Hamiltonians. Next steps: quantum resource estimation, hardware-aware circuit compilation, experimental collaboration.

**Phase 4 (v2.7.0+): Nuclear extension and framework delineation.** Extension of the GeoVac framework to nuclear shell model Hamiltonians using the hyperspherical (HO + spin-orbit) basis, together with a precise delineation of which parts of the framework transfer to non-Coulomb systems. Key outcomes: (1) the potential-independent angular sparsity theorem (Paper 22) — ERI density depends only on l_max, not on V(r), verified at l_max=3 under the global-M_L rule (value in Paper 22; the pair-diagonal figure is retired, 2026-08-29); (2) the deuteron and He-4 qubit Hamiltonians (Paper 23, Tracks NE/NF) using the Minnesota NN potential, Moshinsky-Talmi brackets, and a two-species tensor-product JW encoding; (3) the Fock projection rigidity theorem (Track NH) — the S^3 conformal projection is unique to -Z/r and does NOT transfer to HO or Woods-Saxon; (4) the composed nuclear-electronic deuterium PoC (Track NI) demonstrating the block architecture at 26 qubits with hyperfine validation; (5) the Bargmann-Segal lattice (Paper 24, Track NK) — discrete graph encoding of the 3D HO on the holomorphic sector of S^5, bit-exactly π-free in exact rational arithmetic at every N_max (verified at N_max=5: 56 nodes, 165 edges, zero irrationals); (6) the universal vs Coulomb-specific partition as the Phase 4 conceptual result, completed by the Coulomb/HO asymmetry analysis in Paper 24 (the asymmetry is representation-theoretic and arithmetic: π enters both sides through a measure and spectral sums, as the Kustaanheimo–Stiefel lift shows; the calibration-π layer was withdrawn 2026-10-19). Track NJ memo shelved the nuclear → alpha connection on its empirical result and the Fock rigidity narrowing (its upgrade to "definitively shelved" rested on calibration π being Coulomb-specific and was withdrawn with that layer, 2026-10-19). HO rigidity theorem (Theorem 3 of Paper 24) is the structural dual of the Fock rigidity theorem.

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
*Status:* **PROVEN for the commutative factor — unconditional (2026-06-10; constant re-priced 2026-09-03; scope narrowed 2026-09-29, PI-approved). The almost-commutative extension's convergence is OPEN.** Paper 38: the truncations of the TRUE Camporesi–Higuchi window (complete Dirac eigenspaces) converge to the round-S³ spectral triple in van Suijlekom's state-space GH distance (translation-seminorm metrization; frozen falsifiers `tests/test_p38_true_ch_window.py` (true-window injectivity) + `tests/test_p38_action_seminorm.py` (smoothing identity; its T1 is on the doubled model), re-pointed 2026-09-30); the bound d ≤ γ_n and γ_n → 0 are proved; γ_n is the Fejér-kernel moment, (4/π) log n/n (dual-Coxeter; 2/π on the unit S³, 2√2/π under Kac's form) [INTERNAL THEOREM + MEASURED], and 4/π is a kernel constant, not the convergence rate: a tapered amplitude gives d = O(1/n) [MEASURED; proof and optimal constant OPEN] (re-priced 2026-10-03, FULL #12). Convergence in the Dirac seminorm is OPEN (on the true window its kernel is the scalars, n_max ≤ 3), and so is convergence of the gauge-network / A ⊗ M_n(ℂ) structure the WH names — Paper 38 proves the commutative factor only. On the true window the J-relations hold exactly with forced KO-3 signs, the truncation is odd, and order zero/one fail and grow (a Toeplitz boundary term; `tests/test_true_window_real_structure.py`). Lemma L5 refuted 2026-09-03 (the unconditional theorem uses reach-type estimates only). Lorentzian extension DESCOPED (P45).

**WH2 — Paper 18 is the Seeley-DeWitt + ζ-invariant decomposition of this spectral triple.** The transcendental taxonomy is the structured output of spectral-action geometry, organized by operator order × bundle type.
*Falsifier:* a transcendental in a GeoVac observable that cannot be placed in the grid.
*Status:* three of four axis-quadrants filled; mixed-Tate sharpening POSITIVE (2026-06-03) — M2 on S³ sits in the pure-Tate sub-ring ⊕_k π^{2k}·ℚ (Fathizadeh–Marcolli inherited). See `debug/sprint_mixed_tate_test_memo.md`.

**WH3 — The lattice exists a priori; match to physics is evidence, not derivation.** The packing construction is independent of known physics; its persistent match (Fock S³, nuclear magic numbers, Dirac fine structure, Pauli sparsity) is evidence physics is hosted by a discrete spectral triple, not that the lattice was reverse-engineered.
*Falsifier:* the packing construction fails to force the (n, l, m, s) structure or the n²−1 spectrum; or a match turns out to depend on a hidden physics-informed parameter.
*Status:* origin-story framing permits strong ontological claim internally; papers stay under §1.5. **Graph-limit leg REFUTED (2026-09-25, PI-directed):** the edge Laplacian D−A does not converge to Δ_{S³} (n_max disconnected l-blocks, n_max zero modes, zero cross-block heat kernel); the forced (n,l,m) structure and the n²−1 spectrum survive as the exact label set plus the diagonal SO(4) Casimir L²+A² (`tests/test_graph_s3_nonconvergence.py`). The falsifier's "n²−1 spectrum" clause is therefore met by the Casimir on the labels, not by any graph limit.

**WH4 (deflated 2026-05-07) — The four-way S³ unity is one Fock-projection statement plus three forced consequences.** Bertrand + SO(4) force S³; the Hopf base, the CH spinor bundle, and SU(2) Wilson all follow from S³ = SU(2) parallelizability / maximal-torus structure. Does NOT extend to inner-factor selection (Yukawas, generations remain unforced).
*Falsifier:* a construction forcing one of the four roles onto a different manifold; or a published argument decoupling a consequence from the Fock input.
*Status:* deflated to a single-input forcing statement; outer structural unity essentially closed (Sprint TS-D + Paper 38).

**WH5 — α is a projection constant, not a derivable number.** K = π(B + F − Δ) composes three spectral objects of distinct type (finite Casimir trace; Fock Dirichlet ζ(2); Dirac boundary count 1/40); the right open question is why the sum equals α⁻¹, not how to derive each piece.
*Falsifier:* a spectral-triple construction deriving K as a single coefficient of a well-defined functional.
*Status:* TWELVE mechanisms eliminated (Phases 4B–4I + Sprint A + Sprint K-CC, including the T9 algebraic obstruction: no single CC heat-kernel expansion contains B, F, Δ as terms). Standing reading: three-regime projection coincidence. Paper 2 stays in Observations; combination rule labeled an Observation — not a conjecture or derivation (§13.5 hard prohibition; conjecture→observation downgrade 2026-06-14, PI direction — "conjecture" judged to carry unearned confidence that a derivation exists).

**WH6 — GeoVac's RH-adjacent object is the Dirac spectral zeta D(s), not classical ζ.** Internal GUE-like zero statistics (CV ≈ 0.35–0.40); the classical-RH bridge is closed by three independent walls (zeros not on one line; no spectral-triple-natural functional equation, 48 OoM; wrong Weyl class).
*Falsifier:* D(s) zeros on a single critical line at larger samples; or a natural functional equation closing the RH-O gap.
*Status:* paused (2026-04-18); if resumed, the target is D(s) and its spectral-action interpretation.

**WH7 — Time-discreteness is observer-compactification (registered 2026-06-10, PI direction).** The only temporal structure the framework can see metrically is compactified time, and compactified time is automatically discrete. Inputs: (i) Paper 35 — π enters exactly at temporal compactification (Matsubara 2πk/β); (ii) P45 annihilation theorem — ℝ-time is Lipschitz-invisible in the v1 architecture; (iii) Paper 47 — the three temporal carriers are spectrally indistinguishable. Reading: discreteness of time is supplied by the observer's compact integration window (KMS β = 2π, four-witness theorem, Connes–Rovelli thermal time) — time is the prototype free-side projection. Temporal restriction of the organizing observation below.
*Falsifier (primary):* the Toeplitz temporal-compression program — a metrically visible temporal algebra on a NON-compact carrier without compactification weakens WH7 to convention; a proven annihilation-type obstruction for the framework's whole non-compact temporal class forces it.
*Falsifier (secondary, from Paper 35):* a GeoVac observable containing π whose evaluation provably involves no temporal/spectral integration.
*Status:* REGISTERED (2026-06-10); **PRIOR ART FOUND + a technical hole opened, 2026-09-06.** Load-bearing leg is **(i)**. Prior art (verified): Dolce's Elementary Cycles (arXiv:1707.00677; Ann. Phys. 327, 1562) already postulates compact/cyclic Minkowski time as *the* quantization condition, unifying Matsubara + KK + particle-in-a-box — ontological, not observer-relative, which is now WH7's **only** distinguishing content and must be stated explicitly. Rovelli's RQM Postulate 1 ("relevant information is finite for a system with compact phase space" ⇒ "the characteristic discreteness of quantum theory", IJTP 35, 1637) anticipates the compactness→discreteness chain generally. **Technical hole (aim the falsifier here):** a compact *interval* is not a compact *circle* — a finite real-time window with no boundary condition gives a broadened but still CONTINUOUS spectrum; 2π/β spacing requires periodic identification. And the Matsubara/KMS circle is *imaginary* time while an observer's lifetime is a *real*-time interval; WH7 identifies the two silently. Unverified but flagged: bounded-region (diamond) algebras are Type III₁ with continuous modular spectrum, and the one rigorous compact-window computation (Martinetti–Rovelli, CQG 20, 4919) yields a temperature but *no* discrete spectrum — if that holds it is fatal to the strong form. Honest cap is itself published: Kempf, NJP 12, 115001 (2010), continuous↔discrete equivalent under bandlimitation — cite, do not present as an internal caveat. Prior verdicts unchanged (Toeplitz "weakens-to-convention"; Lorentzian chase structurally closed; μ₄/CM re-reading). **Claim-level restatement DONE 2026-10-16 (PI-directed, v5.20.18):** Paper 35 §XIV.B and its abstract now state the window-as-circle identification as adopted, not derived, carry both caveats (interval vs circle; Euclidean vs real time) and cite Kempf; the group6 synthesis matches. Tier unchanged; whether a real-time window forces the circle stays the open falsifier target. Scans: `debug/lit_scan/{compactness_discreteness,info_theoretic_compactness}_memo.md`. Superseded text: `docs/wh_register_history.md`.

**WH8 — The Born measure is the exchange constant of the observation projection (registered 2026-07-01, PI direction).** The compact↔continuous junction accounts for the discreteness of outcome menus (records are bound/compact ⇒ spectra discrete) and the diagonality of effective observer states (compact-window KMS, four-witness theorem): the same projection that makes time discrete and injects 2π (WH7) injects |ψ|² at observation. The selection measure itself is skeleton-external, both directions now closed: *derivation*-direction closed structurally 2026-05-26 (Paper 34 §VIII Born entry — three routes all reduce to standard GNS/Gleason inheritance; sprint `ahha_born_rule_attempt`); *generation*-direction tested negative 2026-07-01 (the graph's counting measure = Born exactly on the flat locus and only there — exact-rational TV=1/3 on the generic 4-multiplet; TV=0.86 vs the He graph-native CI ground state, 90.8% on 1s²; the flat locus is dynamically unstable under the graph's own H, TV~t², exponent 2.00). Given the projection lattice the import is unique (Gleason, dim≥3) — the residual mystery is single outcomes, not the measure's form.
*Falsifier:* a skeleton-side rule (state-independent, built from the graph's discrete invariants — degeneracies, selection rules, rational couplings) reproducing Born statistics on non-flat physical states (breaks `tests/test_wh8_born_probe.py`); or a composition of more-primitive skeleton projections that carries the measure natively (would upgrade the import to derived and falsify the exchange-constant classification).
*Status:* REGISTERED (2026-07-01); Step-1 probe same day **NEGATIVE-as-expected** (`debug/sprint_wh8_born_probe_memo.md`; falsifier `tests/test_wh8_born_probe.py`; Paper 34 §VIII Born entry updated with the quantitative converse). The external-input Class-1 placement of the Born rule is upgraded classification→tested-negative. Honest cap: this fences, not solves, the definite-outcome problem — no interpretation on the market resolves single outcomes without added structure; papers stay under §1.5 (Paper 34 §III.28 + §VIII carry the tier-honest record).

---

**Organizing observation — compactness *forces* discreteness; the converse is false (corrected 2026-09-06, PI-directed).** Compact groups have discrete spectra (Peter-Weyl), and that is the mathematical content of the S³ graph: the discrete, bit-exact skeleton is the closed/compact regime. **Do not state this as an identity.** "Discreteness *is* compactness" is false as a biconditional over spaces — Maz'ya–Shubin (Ann. of Math. 162, 919) give purely discrete spectra on non-compact ℝⁿ from a confining potential alone, and the corpus's own Paper 24 harmonic oscillator is that counterexample. The correct general statement is at the level of the **operator**: discrete spectrum with finite multiplicities ⇔ **compact resolvent** (already the definitional axiom in WH1's spectral triple). Compactness of the group or manifold is a *sufficient* mechanism and the one operative for S³; it is not necessary. **Credit:** the physics reading is not ours — Liu–Noui (CQG 34, 135008) state it causally for LQG's SU(2), Rovelli's RQM Postulate 1 (IJTP 35, 1637) derives discreteness from compact phase space, and Dolce (arXiv:1707.00677) states it for time. **What remains ours is the accounting**, not the mechanism: each released axis carries an identifiable transcendental price, so the skeleton is rational/forced and calibration (continuum, un-packed data) is transcendental/free. Captured in Paper 18 §III (with the attribution paragraph and this correction) and Paper 35 (temporal compactification injects π). The "second packing axiom" framing is retired; the 2026-05-30 confinement reframing is archived as an organizing reading (2026-05-31).

---
## 1.8. Multi-Observable Focal-Length Decomposition Program (RETIRED)

**Dormant.** Queued 2026-05-09 with five deliverables (Paper 34 §V.C.2–6);
**none was ever built** (verified 2026-09-01); successor = the v4.71.0
chemistry-error sprint. Verbatim: `docs/focal_length_program.md`. Reviving
it is a PI call.

---
## 2. Current Development Frontier

> Full sprint chronicles live in `CHANGELOG.md`. This section is a compact index. Sprint detail is in the memos linked below.

> **⚠ ORIENTATION — read before answering anything about polyatomics, accuracy, or which basis GeoVac uses.** Two facts a fresh session otherwise re-derives at the PI as a corrective lecture (PI direction 2026-08-14, after three derailed sessions): **(1) exact ≠ accurate.** The closed-form ERI arc (v4.77–79) bought *decidability*, not accuracy — the whole integral-exactness axis is worth a few mHa on H₂ and then flatlines, while basis size is worth an order of magnitude more and rising; v4.73.0 localized the chemistry defect to **100% max_n**, bit-invariant to angular and quadrature refinement. **(2) three bases are in play** — theory = Coulomb-Sturmian (shared p₀), `composed_qubit`/`two_center_eri` = **hydrogenic** (a = Z/n), `noci_engine` = STO shapes fitted by Gaussians (an *evaluator*, not a basis). Full map + the three-centre wall + "Paper 58 is not Avery and not Sturmian": `memory/polyatomic_state_of_play.md`.

> **Older entries (215 bullets, rounds 2–8; round 8 = everything from v5.15.24 back) are in `docs/development_frontier_archive.md`.** §2 is the index; CHANGELOG.md is the chronicle.

- **Ladder punchlist: second chain slice and real 3-D star (2026-10-21, v5.21.4, pre-registered):** (t;1,2) again K3 integrated once; star on F_4 rung, no K3; X3-1b UNCLEAR. CHANGELOG v5.21.4.
- **W G = 0 derived; fire_test collection errors VOID (2026-10-21, v5.21.4, pre-registered):** no Mellin source, N2 oscillates; slice and chain membership derived, contour-shift cap only; WG1-3 held. CHANGELOG v5.21.4.
- **Molecular EFT: two-centre finite-size constant; frozen core as contact series (2026-10-20/21, v5.21.3, pre-registered):** MF, CE void (controls); MFB1 held, S_mol(4) > 0, FINDING; contact series transfers worse than PK. CHANGELOG v5.21.3.
- **Slater slice: artanh term drops; P G = 0 reduced to the measured W G = 0 (2026-10-20, v5.21.2, pre-registered):** order-four factor P; parity carries W G = 0 to P G = 0, slice and chain; SA1, MP3 held; MP1-2 missed. CHANGELOG v5.21.2.
- **EFT check, novelty audit, period and deuteron tests (2026-10-19/20, v5.21.1, PI-directed):** nuclear doors = NRQED matching; 7 of 34 held predictions not in literature; X4-1, P3e-3, DS1/3-5 held; X3-1, P3e-4, DR1-5, DS2 missed. CHANGELOG v5.21.1.
- **§1.6 restated; §2 compaction round 8 (2026-10-19, v5.21.0, PI-authorized):** calibration-π reading removed, Track NJ back to "shelved"; 131 bullets moved verbatim to the frontier archive. CHANGELOG v5.21.0.
- **KS lift and hyperfine door; calibration-pi layer withdrawn (2026-10-19, v5.20.29, pre-registered, PI-delegated):** hydrogen = constrained 4D oscillator, Z = n p0, one pi/4; Zemach priced, ZM4 partly failed. CHANGELOG v5.20.29.
- **Dirac finite nucleus near Z alpha = 1 (2026-10-19, v5.20.28, pre-registered):** R^(2gamma) below; radius compulsory past it, diving radii in ratio e^(pi/sigma); Z_cr 169.7; RD3 threshold missed. CHANGELOG v5.20.28.
- **Finite-nucleus defect derived; other profiles (2026-10-18, v5.20.27, pre-registered):** one consumed moment per order (<r^2>, Friar), then universal ln(2ZR)+2gamma_E; J14 = 7th multi-focal instance. CHANGELOG v5.20.27.
- **Zero energy and a finite nucleus as doors (2026-10-18, v5.20.26, pre-registered):** E=0 contracts SO(4) to E(3); finite nucleus releases n from integrality, prices a consumed radius (J14, C). CHANGELOG v5.20.26.
- **Spin and the hydrogen continuum as released axes (2026-10-18, v5.20.25, pre-registered):** spin contraction keeps n, pi only as limits; continuum = continued bound skeleton x 1/(1-e^{-2pi eta}). J12/J13. CHANGELOG v5.20.25.
- **Seam F relative to the release; c_{6,0} drop read by its caveat (2026-10-18, v5.20.24, PI-decided):** J-row F covers surviving labels + the released class; measured parts outside. CHANGELOG v5.20.24.
- **SMALL sweep of the v5.20.21 DELTA debt (2026-10-18, v5.20.23, PI-directed):** 44 of 49 fixed, 5 left with reasons; 4 guards fire-tested 10/10; self-check 2,3,1,0. CHANGELOG v5.20.23.
- **/qa claim-search DELTA on v5.20.21 = DEFECTS, remediated (2026-10-18, v5.20.22, PI-invoked):** 2 OWNER + 10 SURFACE LARGE, math held; conductor-4 narrowed; self-check 4,5,3,2,1,0. CHANGELOG v5.20.22.
- **Rings + Ising class on period ladder (2026-10-18, v5.20.21, pre-registered):** Coulomb ring of three drops to elliptic; correlated triangle stays K3; ring of four STOP; c_{8,0} breaks parity. CHANGELOG v5.20.21.
- **Four linked densities = first rung beyond elliptic (2026-10-16, v5.20.20, PI-directed):** tree identity; F_4 operator = AESZ #58, Sp4, not elliptic-built; Hamiltonian ERIs stay elliptic. Paper 59 sec:four_density; CHANGELOG v5.20.20.
- **/qa claim-search DELTA on v5.20.16-18 = DEFECTS, remediated (2026-10-16, v5.20.19, PI-invoked):** math held; owner wording fixed; rung rule J9 F->Obs; 12 C16 + 4 C25 entries; self-check 4,4,1,0; 4 PI items. See CHANGELOG v5.20.19.
- **/walls PERIODS-TRANSCENDENCE = CRYSTALLIZED; Mellin conventions resolved; Paper 35 window-as-circle restated (2026-10-16, v5.20.18, PI-invoked):** class set by density topology, not operator/basis/nuclei; brute-T2 row MIS-SCOPED; P3e-2 filed. See CHANGELOG v5.20.18.
- **Paper 59 homes + guards for the four-centre and K3 results; Mellin engine single-sourced (2026-10-15, v5.20.17, PI-directed):** J4/J7 promoted to F; 25 restatement loci -> pointers, 6 discrepancies fixed; forcing catalogue archived; Paper 57 trimmed. See CHANGELOG v5.20.17.
- **The atlas exists as one gated register (2026-10-15, v5.20.16, PI-directed):** `docs/atlas.md`, one row per released axis with price/verdict/falsifier; seam + Paper 57 widened to 72 entries (W status, sub-domain J); field guide re-centred on it. See CHANGELOG v5.20.16.
- **Pre-registered period predictions (2026-10-13/14, v5.20.14-15):** four-centre ERI genus 1; three-electron chain = K3, period elliptic squared (Bailey); true fibre collapses at equal scales, not generically (PSLQ evidence). See docs/predictions/ + CHANGELOG.
- **Water with GeoVac-native bricks (2026-10-11, v5.20.12):** hydrogenic and Sturmian kits all short of cc-pVDZ at 24 functions (best: two-scale Sturmian, 130 mHa); §3 row added. See CHANGELOG v5.20.12.
- **Avery load-vs-generate, atoms (2026-10-11, v5.20.13):** loading angular+radial tables instead of every ERI nets 1.2x->2.5x (N=2..6) vs plain sparse; DF/THC baseline owed. See debug/sprint_avery_load_vs_generate_memo.md.
- **Avery material: PDF canonical, web page RETIRED (2026-10-11, PI):** `docs/outreach/avery_hyperspherical_map.pdf` is the one version; the private claude.ai artifact "Labels to Generalized Sturmians" is retired (never shared; do not update).
- **Water: Slater DZP vs cc-pVDZ at equal size (2026-10-11, v5.20.11):** own validated CCSD(T); Slater 41 mHa lower, identical Pauli count, 1-norm +4.7%; matched-accuracy Gaussian test owed. See CHANGELOG v5.20.11.
- **m-seam probe + water three-centre ERIs (2026-10-10, v5.20.10):** bending opens m (Renner/Walsh, prior art); shared-centre (Barnett-Coulson) channels evaluate water's three-centre ERIs to 1e-9: closed-form wall only. See CHANGELOG v5.20.10.
- **Avery-material wording made neutral before first push (2026-10-10, v5.20.9):** CHANGELOG/test/map no longer reference the shared draft or its table entries. See CHANGELOG v5.20.9.
- **Avery reading map + credit corrections (2026-10-10, v5.20.8):** docs/outreach/avery_hyperspherical_map.pdf (README-linked); P19 SW weight, P60 molecular T0/1-norm wording, Wen-Avery fixed. See CHANGELOG v5.20.8.
- **Avery-book pass: Paper 60 molecular levers are s-sector only (2026-10-09, v5.20.7):** σ basis cond exponential; no preconditioner (unitary-translate redundancy), phase-out only; P58 NaH fit sensitivity. See CHANGELOG v5.20.7.
- **LiH analytic e-n factors -8.05851 + prior-art attribution (2026-10-09, v5.20.6):** parameter-free core-orbital switches (amplitude + density) gain 6.9 mHa, VMC 0.6 sigma; one-centre reductions credited to Hylleraas-CI. See CHANGELOG v5.20.6.
- **LiH residual diagnosed; analytic CI-relaxed -8.05164 (2026-10-05, v5.20.5):** r12t2 STOP (<=3 mHa); CI relax 2.6 mHa (analytic, new best); lever = Li-core e-n factors with the geminal (VMC -8.0637). See CHANGELOG v5.20.5.
- **Owed items 1-4 (2026-10-04, v5.20.4):** P27 stale Be EP-2N numbers fixed (sign-fix fallout); P32 item 4 = J^2=-1 obstruction; heavy tests xdist-safe; LiH geminal gamma is a real parameter (golden bit-identical). See CHANGELOG v5.20.4.
- **LiH/prolate stack productionized (2026-10-04, v5.20.3):** ~11k lines debug/ -> geovac/ (prolate_allelectron/, lih_r12ci/, lih_vmc); 42-item golden oracle, all bit-identical incl. -8.049045. See CHANGELOG v5.20.3.
- **/qa trunk FULL #13 = FAIL, remediated (2026-10-03, v5.20.2):** 2 OWNER-LARGE (taper guard blind to a log rate; P32 corollary 'confirmed/closure') + 22 SURFACE (11 via self-check, 4/4/2/1/0); 4 C16 families; trunk SMALL sweep done (real_structure sign bug fixed). See CHANGELOG v5.20.2.
- **Item 5 rounds A+B + slow-test closed form (2026-10-03, v5.20.1):** sprint record relocated verbatim to docs/chronicles/ (P32 8,642 -> 6,174 lines; P18 -84); OOM Tomita tests 0.5 s via closed-form K_TT. See CHANGELOG v5.20.1.
- **/qa trunk FULL #12 = FAIL, remediated under incremental stamps + edit bar (2026-10-03, v5.20.0, PI-directed):** 4/π = Fejér-kernel constant, taper gives measured O(1/n); K=2J_z not wedge-preserving; single-cutoff → Observation. See CHANGELOG v5.20.0.
- **/qa trunk DELTA #26 = DEFECTS, remediated (2026-10-03, v5.19.2, PI-approved):** P40 4/π universality simple-G only; SD ratio-4 + production SD sign error fixed; ∇*∇ not π^even; self-check 10 rounds to zero LARGE. See CHANGELOG v5.19.2.
- **Paper 39 main-theorem coverage gap closed (2026-10-02, v5.19.1, PI-directed):** `tests/test_paper39_lifted_state.py` backs the lifted-state route (joint moment sandwich + product-state identity), fire-tested. See CHANGELOG v5.19.1.
- **/qa trunk DELTA #25 = DEFECTS, remediated (2026-10-02, v5.19.0, PI-approved):** P57 false Δ=c²(3,2); P1 ladder [SYMBOLIC PROOF]; P39/P40 C₃=1 by compression; new rule: families + self-check to zero LARGE (15 rounds); CLAUDE.md de-numbered. See CHANGELOG v5.19.0.
- **/qa trunk DELTA #24 (first claim-search DELTA) = DEFECTS, remediated (2026-10-02, v5.18.1, PI-approved):** 5 OWNER + ~20 SURFACE LARGE, no math defect; C25 hardened; 17 C16 entries swept 85 loci. See CHANGELOG v5.18.1.
- **/qa instrument: claim registry C25 + OWNER/SURFACE-LARGE stop rule + claim-search DELTA; trunk FULL #11 = FAIL, remediated (2026-09-30, v5.18.0, PI-directed):** 5 LARGE, no owner math defect; C25 found 4 front-door defects on run one. See CHANGELOG v5.18.0.
- **/qa trunk DELTA #23 = DEFECTS, remediated (2026-09-30, v5.17.17, PI-approved):** 1 LARGE + 1 borderline, both PM-fix residue ('does not force π'; P45 hybrid mislabel); 3 C16 gates; guard aliasing escape closed. See CHANGELOG v5.17.17.
- **/qa trunk DELTA #22 = DEFECTS, remediated (2026-09-30, v5.17.16, PI-approved):** 1 LARGE (10 -> 1): cb-norm label in tests/geovac; code-layer C16 gate; π-source converse 'not established', not false. See CHANGELOG v5.17.16.
- **/qa trunk DELTA #21 = DEFECTS, remediated (2026-09-30, v5.17.15, PI-approved):** ~10 LARGE (idea-level sweep reached old restatements); π-source converse (iii)⇒(i) WITHDRAWN; P39 cb-norm = 1. See CHANGELOG v5.17.15.
- **/qa trunk FULL #10 = FAIL, remediated (2026-09-30, v5.17.14, PI-approved):** 10 LARGE, no keystone math defect; P32 F=ζ_Δ(2) false, M3 carries no π, WH1 falsifier re-pointed to the true window. See CHANGELOG v5.17.14.
- **/qa trunk DELTA #20 = DEFECTS, remediated (2026-09-30, v5.17.13, PI-approved):** 4 LARGE (plateau); Paper 28 T9 false at half-integer s (fixed + test); 'M1 alone derives β' retired; next = trunk FULL run. See CHANGELOG v5.17.13.
- **/qa trunk DELTA #19 = DEFECTS, remediated (2026-09-29, v5.17.12, PI-approved):** 4 LARGE (12→8→5→4); approach change: idea-level gates + corpus-wide sweeps; 'Hopf-base measure' label retired (~77 loci). See CHANGELOG v5.17.12.
- **/qa trunk DELTA #18 = DEFECTS, remediated (2026-09-29, v5.17.11, PI-approved):** 5 LARGE (12→8→5), one-sentence unmoved-locus residue; code + citations clean; self-check caught a noise-pinning guard. See CHANGELOG v5.17.11.
- **/qa trunk DELTA #17 = DEFECTS, remediated (2026-09-29, v5.17.10, PI-approved):** 8 LARGE (12→8), all unmoved-locus residue; code + citations clean; cubic's root α_K proved transcendental (not physical α). See CHANGELOG v5.17.10.
- **Drafts 52/53 resolved (2026-09-29, v5.17.9, PI-approved):** P52 archived (KPS values are continuum; truncation never enters); P50 corrected; P53 rebuilt: flat-plane GH theorem, disk obstruction was a float artifact (negative withdrawn); P51 backbone OPEN; P50 log 2 relabelled M3-ring per P18. See CHANGELOG v5.17.9.
- **/qa trunk DELTA #16 = DEFECTS, remediated (2026-09-29, v5.17.8, PI-approved):** ~12 LARGE, framing only; family sweeps; Fock/Δ overcorrection reversed; WH1 scoped to commutative factor; K 4.8e-7 (8.8e-8 = cubic root). See CHANGELOG v5.17.8.
- **/qa trunk DELTA #15 = DEFECTS, remediated (2026-09-29, v5.17.7):** ~13 LARGE (a third PM-introduced); math all confirmed; new pre-checkpoint self-check caught 4 more LARGE. See CHANGELOG v5.17.7.
- **/qa trunk DELTA #14 = DEFECTS, remediated (2026-09-28, v5.17.6):** summaries lagging corrected owners (incl. mirror direction); Paper 43 ARCHIVED into P42 §11; P16 ν_min a theorem + prior art. See CHANGELOG v5.17.6.
- **FULL #9 open items closed (2026-09-28, v5.17.5, PI-directed):** Paper 16 ν=N−2 RETRACTED → Pauli-filling ν_min (two routes; Ne=10); 8 guards fire-tested; Pythagorean 1/π² not M1. See CHANGELOG v5.17.5.
- **/qa trunk FULL #9 = FAIL, remediated (2026-09-28, v5.17.4, PI-approved):** P32 diagonal triple spectral-not-metric; finite four-witness demoted corpus-wide (P42 retitled); P38 A2 confirmed; P7 kernel test. Paper 16 ν=N−2 flagged (PI). See CHANGELOG v5.17.4.
- **Trunk list-closing pass (2026-09-28, v5.17.3):** final re-check's ~13 LARGE closed gate-first (4 C16 entries); PI rulings applied (α approx./Observation, p-value, Δ-sign reading, P0/P7 titles); P34 III.5 real test. Next: trunk FULL run. See CHANGELOG v5.17.3.
- **/qa trunk re-check = DEFECTS, round-2 remediated (2026-09-28, v5.17.2):** ~10 LARGE residue paraphrases (<0.1% cells, Paper 8 SO(4) covariance, P31/P36/P13); C16/C17 widened; 51× not 54× (identity convention). Round 3 (whole-paper fix-mode) + final re-check: family retired; ~13 LARGE pre-existing debt → list-closing pass, then trunk FULL run. See CHANGELOG v5.17.2.
- **/qa trunk DELTA = DEFECTS, remediated (2026-09-28, v5.17.1, PI-directed):** ~19 LARGE, all one family (graph/κ credited with the hydrogen spectrum); v4.76.0 audit: refutation stands; 2 C16 classes, ~45 papers, Paper 7 retitled, AtomicSolver default → physical H. See CHANGELOG v5.17.1.
- **/qa gains RESEARCH-GRADE = no open LARGE defects (2026-09-27, v5.17.0, PI-directed):** ends the fix-SMALLs loop; LARGE/SMALL grading; SMALLs fixed on sight (no re-run); stamped + expiring; unfixable SMALLs in `docs/qa/<target>.debt.md`. DoDs pre-registered for the post-v5.16.0 DELTAs. See CHANGELOG v5.17.0.
- **LiH analytic energy −8.0490 = corpus-best (2026-09-27, v5.16.0):** π marriage + cusp-correct e^{−1.5r} geminal; gap was geminal shape (Kato cusp), not count; g_Vee triangle-weight bug fixed; VMC twin 0.16 mHa. Paper 19; CHANGELOG v5.16.0.
- **Graph→S³ convergence REFUTED; re-anchored on Fock + SO(4) Casimir (2026-09-25, v5.16.0, PI-directed):** D−A has n_max disconnected l-blocks; labels/degeneracy survive exactly. ~27 papers re-worded. `tests/test_graph_s3_nonconvergence.py`. See CHANGELOG v5.16.0.
- **He graph-native 0.216% RETIRED (2026-09-25, v5.16.0, PI-directed):** a non-physical (1/16)A one-body hop (−48 mHa); physical fixed-k CI = 1.88%. H⁻ 21%, Z_c≈1.84, "86%/10–39% correlation" were the hop. See CHANGELOG v5.16.0.
- **Sign problem = frustrated ±1 holonomy (2026-09-25, v5.16.0):** on physical He, Gaunt zeros cut cycles ~15× and Δ 43–48× vs rotated; signs no better than random; radial basis not neutral. Paper 22 §sign_holonomy.
- **π/M=16 marriage: −8.043494 at trunc 1e-3 via O(nnz) g_Vee (2026-09-25, v5.16.0):** 267 GB→8.8 GB; VMC G2 PASS to 0.14 mHa (0.56σ) after proper burn-in (an earlier 3σ gate-6 offset was a burn-in artifact). `debug/lih_pi_m16_build_plan.md`.

**Best results by system type:** values live in the papers (authoritative), `docs/claims_register.md` (one line per headline claim, with its tier) and `debug/qa/numeric_registry.py` (gated by C21); CLAUDE.md deliberately carries no measured values (PI direction 2026-10-02: a stale number here is copied into new prose with full authority). Map:

| System | Best route (qualitative) | Paper |
|:-------|:-------------------------|:-----:|
| He | 2D hyperspherical variational (raw, properly variational); cusp extrapolation marginally below exact (non-variational); fixed-k = Z hydrogenic CI with the physical Hamiltonian above the HF limit. The "graph-native" He figure is RETIRED (a non-physical (1/16)A hop, 2026-09-25) | 13 |
| H⁻ | Bound, above exact (variational), fixed-k CI with the physical H; the "over-binds" / Z_c figures were hop artifacts | 13 |
| PsH | Level 3, sign-flipped charge | 13 |
| H₂⁺ | machine precision (reference-limited) | 11 |
| H₂ | near-exact D_e via prolate explicit-r₁₂ CI with exact integrals (registry `p12_r12_de_pct`); re-based CI just below; Level 4 hyperspherical is the structural treatment | 12/15 |
| LiH | best R_eq: Route C prolate all-electron with an analytic core; best energy: analytic RI-free explicit-r₁₂ on the M=16 π CI (VMC twin agrees); composed l-dep-PK for the qubit pipeline | 19/17 |
| BeH₂, H₂O | composed (+ exchange / lone pairs), zero parameters | 17 |
| LiH (4N) | full 4e mol-frame, PK-free: equilibrium exists, D_e unbound at l_max=2 | 17 |
| Composed Pauli | exactly linear in Q across molecules (registry `composed_coeff`); large equal-qubit advantage vs Gaussian | 14 |
| Atomic Pauli | equal-qubit comparison vs cc-pVDZ/cc-pVTZ is UNFAVORABLE at small Q and favorable at larger Q (matched qubits, NOT matched accuracy: GeoVac He is tens of mHa above the Gaussian at equal Q) | 14 |
**Key structural results (details in papers and CHANGELOG.md):**
- **κ = −1/16 (v2.26.1; Observation, not derived — v4.13.0 QA):** the matching κ *coincides* with the geometric 1/16 of the Fock projection (a numerical observation; no derivation bridge — never "derived"). Not fitted. See `debug/probe_kappa_sprint_memo.md`, `memory/kappa_observation_not_derived.md`.
- **α structural decomposition (Phases 4B-4I, April 2026):** B=42 (Casimir), F=π²/6 (Fock Dirichlet at d_max), Δ=1/40 (Dirac degeneracy g_3). Three spectral homes of distinct type; combination rule K=π(B+F−Δ) is a numerical observation, not derived (conjecture→observation downgrade 2026-06-14). 12 mechanisms eliminated. See Paper 2, CHANGELOG.
- **Spectral-action supertrace (v2.26.1):** the "SD cancellation theorem" (ratio dim(S) at every order, perturbative supertrace ≡ 0) is WITHDRAWN (2026-10-02, DELTA #26: exact S³ traces give ratios 2, −1, 0), Δ⁻¹=40 from Euler-Maclaurin, (−) sign consistent with the (-1)^F grading (a reading, not a derivation; PI 2026-09-28). Two-term exactness on S³ Dirac. See `debug/st_supertrace_sprint_memo.md`.
- **Nuclear systems (Paper 23):** Deuteron and He-4 (16 qubits) and composed nuclear-electronic deuterium (26 qubits) qubit Hamiltonians (Pauli counts in Paper 23; corrected 2026-08-22, v5.0.0, when the Moshinsky N_tot guard was removed). Fock rigidity theorem (S³ unique to −Z/r).
- **Angular sparsity theorem (Paper 22):** ERI density depends only on l_max, not V(r). Universal across potentials.
- **Bargmann-Segal lattice (Paper 24):** HO on S⁵ Hardy space, π-free, HO rigidity theorem. Coulomb/HO asymmetry = 7 layers catalogued, 5 standing (Layer 1, the graph Laplacian's spectrum-computing role on S³, withdrawn 2026-09-25 with the graph→S³ re-pricing; Layer 2, calibration π, withdrawn 2026-10-19 by the Kustaanheimo–Stiefel lift; 7th: CM Hodge structure, v4.103.0; χ₋₄ periods transfer, ¼ of P28).
- **S⁵ gauge extension (v2.26.1):** U(1) Wilson transfers; SU(3) NOT natural; CP² quotient fails. See `debug/s5_gauge_structure_memo.md`.

**Classical solver status: INVESTIGATION COMPLETE (v2.0.24).** 30+ tracks exhausted all solver × PK × basis combinations. Structural ceilings characterized. Composed at l_max=2 is the production operating point. See CHANGELOG.md for full track history.

**Quantum computing status: ACTIVE FRONTIER (Paper 14).** Under the exact rule (2026-08-29): the composed Pauli count is exactly linear in Q across molecules (main-group and d-block coefficients in the registry), the within-molecule basis exponent is a small-basis fit, the equal-qubit advantage vs Gaussian is large for H₂O, and the ecosystem export covers OpenFermion/Qiskit/PennyLane. Market test: composed LiH is near parity with STO-3G on raw counts (more terms than its reduced count; matching 1-norm) at more qubits; the decisive wins are the exact Q-linearity, the large cc-pVDZ reduction and structural sparsity. Atomic path: the equal-qubit advantage INVERTS at small Q and recovers at larger Q; the 1-norm is lower but at a markedly worse He energy — disclosed honestly per the benchmarking rule. (He re-measured 2026-09-25: the library had carried a non-physical graph-hop one-body term.) (CF-1 is DISSOLVED: the pair-diagonal "convention" was a wrong-sign-q bug; every figure measured under it is retired and C17-registered.) 37-system `hamiltonian()` library (35 composed molecules + He + H2). Values: Papers 14/20 and the numeric registry. See CHANGELOG Tracks AW-CA.

**Historical arc chronicles** (entropy/EP-2; cusp; balanced coupled; precision-atomic; RH sprint + S_min erratum; and the QED / spectral-triple-WH1 / Lorentzian / gauge / Dirac-on-S³ / chemistry / precision-multi-focal / gravity arcs) moved verbatim to `docs/development_frontier_archive.md` (2026-06-16 compaction). The §6 loading tiers + `papers/INDEX.md` are the current map; CHANGELOG.md is the canonical chronicle; the §2 one-liner index above + the best-results / key-structural tables are the live current-state quick reference.

**Scope:** See `SCOPE_BOUNDARY.md` for supported atoms/molecules (the 37-system `hamiltonian()` library; Z=1-56 via frozen cores).

**Architecture locked:** The LCAO/graph-concatenation approach (v0.9.x series) is superseded. All molecular work uses natural geometry (Papers 11, 13, 15, 17).

---
## 3. Approaches That Failed

Critical institutional memory. **Do not re-derive these dead ends.** Category table below, then a name+count index of the 93 per-sprint instances; the **full account of every row is `docs/failed_approaches_ledger.md`** (moved verbatim). Refinement layer: `docs/walls/register.md`. Sub-agents must check the ledger first.

|:---------|:-----:|:-----------|
| Category | Count | Key Lesson |
| LCAO / single-S³ molecular encoding | 3 | Graph Laplacian kinetic energy is R-independent; need natural geometry where separation occurs |
| PK modifications (projector, spectral, self-consistent, R-dependent, l-dependent on 2D) | 6 | PK provides coordinate-space exclusion that angular projectors categorically cannot replicate; PK is the irreducible cost of composed-geometry factorization |
| Cusp treatments (alpha-only, graph absorption, θ₁₂-adapted basis) | 3 | Cusp is 2D in (α, θ₁₂); 1/r₁₂ is an embedding exchange constant; Schwartz extrapolation is the correct approach |
| Inter-group antisymmetry (node exclusion, fiber bundle, spectral) | 3 | Antisymmetry requires shared coordinate system; composed geometry factorizes into incompatible coordinates |
| Full N-electron radial solvers (adiabatic, coupled-channel, 2D variational) | 3 | Adiabatic overcounts D_e; coupled-channel numerically unstable; 2D gives unbound D_e; angular basis is the bottleneck |
| Polyatomic coupling (Z_eff partition, classical repulsion, lone pair at Z_eff>4) | 3 | Orbital-level exchange coupling needed; lone pair Slater integrals unphysical at high Z_eff |
| Geometric elevation (blow-up, Lie algebra, S³×S³) | 1 | min/max boundary is physical; per-ρ diagonalization structurally irreducible |
| Diagnostic arcs (eigenchannel rotation, spheroidal compression, enhanced Z_eff, midpoint origin) | 4 | Various; see CHANGELOG.md |
| l_max convergence via 2D solver (variational_2d in composed pipeline) | 1 | l_max divergence is structural to PK/composed architecture, NOT from adiabatic approximation. |
| Sturmian CI (Coulomb, generalized) | 2 | Coulomb improves atomic FCI but Lowdin 1-norm inflated 2.8-4.5x; generalized loses within-config flexibility at larger basis; neither bypasses PK ceiling |
| TC Jastrow in adiabatic hyperspherical solver | 1 | TC replaces multiplicative V_ee (O(R) in angular eigenvalue) with first-derivative G_ang (O(1)); adiabatic separation requires V_ee as multiplicative operator; ~46% error; TC needs direct variational/FCI framework, not adiabatic. |
| TC in second quantization (composed qubit pipeline) | 1 | TC qubit Hamiltonians plateau at ~3.4% across n_max=2..5; standard FCI converges 5.3%→2.0%. |
| TC 2D variational cusp correction (basis doubling) | 3 | Three attempts failed (R sin(α) basis worsens, sin(α) negligible, [H,r₁₂]=0 identically). |
| TC angular gradient for l>0 orbitals | 1 | Adds 2.66× Pauli terms for 0.01 pp accuracy. |
| Coupled composition (cross-block ERIs replacing PK) | 1 | Cross-block ERIs add 2.56× Pauli, 2.30× 1-norm; 29% FCI error. |
| Single-center nested molecular LiH | 1 | R_eq 33.7% error; single Z can't represent both core and bond length scales. |
| Two-center charge-center nested LiH | 1 | 48.2% energy error; truncated basis cannot represent 1s² core off-center. |
| Heterogeneous nested (per-pair Z_eff) | 1 | Löwdin orthogonalization destroys Gaunt sparsity (1711 vs 120 Pauli, 14× inflation). |
| Fock energy-shell self-consistency for He | 1 | k²=−2E over-constrains 2-electron problem; SC 12.8%/5.1% vs variational 1.9%/1.6%. |
| SM-running origin for Δ = 1/40 (Paper 2 alpha) | 6 | Six tracks ruled out SM-running origin for Δ=1/40. |
| Schwartz-tail hot-node patch on He Z=2 graph-native CI (CUSP-2) | 1 | Schwartz tail correction worsens accuracy; 0.20% floor is small-Z graph-validity-boundary artifact (Z_c≈1.84), not cusp. |
| Energy graph for V_ee on S³ (Paper-12 analog search) | 1 | Pair-state graph dense (47%); cross-shell denominators don't close. |
| Darwin+MV for He/Li/Be 2p-doublet improvement | 1 | Both 2p states share l=1; Darwin=0 for l≥1; MV cancels in splitting. |
| Balanced + frozen-core PES for second-row molecules | 2 | NaH/MgH₂ overattract monotonically at n_max=2; frozen [Ne] hides core screening from cross-center V_ne. |
| Single-constant graph-to-continuum QED projection (C×F₂→α/(2π)) | 1 | C×F₂ grows with n_max, doesn't converge to α/(2π). |
| σ-vertex and direction-resolved vector QED on Fock graph (VQ-1..VQ-5) | 5 | Five tracks ruled out σ-vertex approaches to vector photon selection rules. |
| Co-exact mode q-labeling for SO(4) channel count / Ward / charge conjugation | 1 | Co-exact eigenvectors have right eigenvalues but wrong support (nearest-neighbor only); triangle inequality fails at high q. |
| Dirac-sector lift of Paper 2 α combination rule ingredients B, F, Δ | 3 | Three obstructions: B doesn't lift (zero mode forces (m−1)), F doesn't lift (Apéry Q-linear independence), Hopf-equivariant decomp doesn't produce B/F. |


### Per-sprint instances (index; full rows in the ledger)

| Approach (full row in the ledger) | N |
|:---|:---:|
| TC angular-gradient ERI assembly (2026-08-30) | 1 |
| Multi-focal spatial composition (Sprint HF, 2026-05-07) | 3 |
| Phillips-Kleinman cross-center barrier for second-row chemistry binding (2026-05-08) | 1 |
| Screened-Schrödinger valence basis (h1 diagonal only, Track 3, 2026-05-09) | 1 |
| W3 spectral-zeta calibration-data identification (Sprint W3, 2026-05-08) | 3 |
| Multi-focal Path C5 saturated basis (numerical-luck f closure, 2026-05-09) | 1 |
| Heuristic two-zeta screening for [Xe] core (CR67 + BBB93 Kr ratios, 2026-05-09) | 1 |
| All-positive-coefficient single-zeta hydrogenic basis lacks radial nodes (Track 2 diagnostic, 2026-05-09) | 1 |
| Multi-zeta physical Na valence basis substitution (Sprint α-Multi-zeta + α-PES, 2026-05-23) | 1 |
| Three-bucket M-Z partition (basis-closable cross-shift) FALSIFIED at NaH max_n=3 (2026-05-23) | 1 |
| Kernel-shape substitution as W1c-residual closure (Sprint F2, 2026-05-23) | 1 |
| Single-particle Pauli orthogonality (rank-1 PK on bonding orbital) as W1e closure (Sprint F4, 2026-05-23) | 1 |
| Mean-field core-bonding J-K (explicit-core Hartree) as W1e closure (Sprint F5, 2026-05-23) | 1 |
| Basis-level Schmidt orthogonalization of H 1s against Na [Ne] core as W1e closure (Day-1 diagnostic, 2026-05-23) | 1 |
| [Ne] core correlation as W1e closure mechanism (Day-1 literature estimate, 2026-05-23) | 1 |
| Basis enlargement to max_n=4 alone as W1e closure (Sprint F6, 2026-05-23) | 1 |
| UV refinement to close α > 1 spinor SC gap (G4-4c week 2, 2026-05-29) | 1 |
| Cutoff-function cure for IR over-count of $S_{\rm BH}$ at small Λ (G4-5c-IR-fix, 2026-05-29) | 1 |
| Attributing the G7/G4-2 Newton-constant factor of 2 to a scalar-vs-Dirac cone coefficient (2026-05-30) | 1 |
| Scalar (periodic BC) conical-defect SC extraction at sprint scale (G4-3c-proper / T1, G4-4e, 2026-05-28/29) | 1 |
| Naive K_var/K_disk → 4(l_max+1)(l_max+2) saturation prediction (G4-4b-c, 2026-05-29) | 1 |
| Small-t polynomial fit for Seeley-DeWitt $a_0$ extraction (G4-4d, 2026-05-29) | 1 |
| Naive $\phi(2)$ Mellin-moment prediction for $S_{\BH}$ cutoff dependence (G4-5d, 2026-05-29) | 1 |
| Bulk $\Lambda^4$ cosmological-constant extraction from 2D disk-Dirac alone (G4-5b, 2026-05-29) | 1 |
| v3.19.0 "Fursaev-Solodukhin spinor double-cover correction" mechanism attribution for Möbius α/(2α-1) (task #26, 2026-05-29) | 1 |
| Geometric-mean azimuthal as cheap UV cure for FD/Spec bracket (task #27, 2026-05-29) | 1 |
| Single-axis N_φ-sweep for G4-6a refined A-coefficient extraction (2026-05-29) | 1 |
| Isotropic multi-axis refinement for A-coefficient extraction (2026-05-31) | 1 |
| soft_IR_frac → 1/(2α) as the Möbius α>1 "mechanism" (Route C, demoted by Sprint GD-2 t-audit 2026-05-29) | 1 |
| Spurious "sign discrepancy" Track α'' thread 9 (2026-05-29, RESOLVED same day) | 1 |
| Möbius α/(2α−1) as a continuum α>1 closed form (B4, 2026-05-29, RETIRES the open thread) | 1 |
| Naive de-compactification of the finite-R Dirichlet disk (Paper 53 B3, 2026-05-29) | 1 |
| Species-II spatial fission-aperture discriminator (which-site entanglement frozen=non-binder / responsive=binder, 2026-05-30) | 1 |
| BW wedge entanglement entropy as Bekenstein-Hawking area law (BH-Phase0, 2026-05-31) | 1 |
| Furnstahl IR extrapolation for E1 polarizability (2026-06-01) | 3 |
| Sturmian basis for deuteron polarizability with Minnesota NN (2026-06-01) | 3 |
| LIT for deuteron E1 polarizability (2026-06-01) | 1 |
| NaH Z_orb scan as chemistry-wall diagnostic (2026-06-01) | 1 |
| Resolvent (D²)⁻¹ for two-body Coulomb interaction (2026-06-01) | 4 |
| J_GV²=−1 as a forcing handle for the inner ℍ factor (Door 4c, 2026-06-01) | 1 |
| Gauged tensor-product spectral action (full double-sum gauge field) for two-body Coulomb radial weights (2026-06-03) | 1 |
| Yukawa values in low-coefficient pure-Tate periods M1 ∪ M2 (Sprint Yukawa-PSLQ, 2026-06-03) | 1 |
| Hopf-tower-to-representation extension as shortcut for forcing N_gen = 3 (Sprint Read 2 scoping, 2026-06-03) | 1 |
| W1e chemistry corrections as outer-factor M1/M2/M3 periods (Sprint W1e period-class, 2026-06-04) | 1 |
| Propinquity-derived Trotter bound at production parameters (Sprint Trotter propinquity, 2026-06-04) | 1 |
| DMRG-on-FCIDUMP closes W1e on NaH (Sprint R3-B falsifier, 2026-06-07) | 1 |
| LiH composed qubit FCI binds (R3-A diagnostic, 2026-06-07) | 1 |
| Existing engineering kwargs close LiH over-binding (LiH kwarg sweep, 2026-06-07) | 1 |
| Explicit-core HF closes NaH binding (Sprint B.1, 2026-06-07) | 1 |
| W1e closure via off-diagonal cross-block h1 WITHOUT orthogonalization (SO(4)-breaking diagnostic, 2026-06-07) | 1 |
| Spectral action expansion has chemistry-side structural content at finite cutoff (Sprint spectral action expansion, 2026-06-07) | 1 |
| Camporesi-Higuchi κ-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint CH κ-parity, 2026-06-07) | 1 |
| Direct-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity direct, 2026-06-08) | 1 |
| Rotated-basis m_j-parity Z₂ as relativistic chemistry tapering stabilizer (Sprint m_j-parity rotated, 2026-06-08) | 1 |
| Spectral action $S(D)(R)$ of assembled MvS Dirac binds LiH (Sprint M-vS-2 + R-sweep, 2026-06-07) | 1 |
| Non-abelian M-vS gauge group reduces Pauli count beyond Z₂ Hopf-U(1) tapering (Sprint M-vS Gauge, 2026-06-07) | 2 |
| K⁺ compression of the Krein Dirac as a Lorentzian quantum-metric device (P45 descope, 2026-06-09) | 1 |
| Momentum-diagonal temporal multipliers as a metric temporal algebra (P45, 2026-06-09) | 1 |
| Grab-bag / single-weight PSLQ bases for weight-inhomogeneous graph-QED constants (2026-06-11) | 1 |
| Levin/EM series acceleration on log-modulated summands (2026-06-11) | 3 |
| Nested adaptive nsum cascade for depth-≥4 multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 |
| `mpmath.sumem` (auto Euler-Maclaurin) for log-modulated trailing-t tails (S^(4) stage-1, 2026-06-13) | 1 |
| float64-contaminated terms inside a high-precision Levin sum (S^(4) stage-1, 2026-06-13) | 1 |
| Brute high-precision summation of b1=2 high-log trailing multiple-t-values (S^(4) stage-1, 2026-06-13) | 1 |
| sympy exact-Q rank for MZV-scale relation matrices (S^(4) stage-2, 2026-06-13) | 1 |
| SO(4) Wigner-D critical points predict LiH R_eq (harmonic phase-lock probe, 2026-03-11; surfaced 2026-06-16) | 1 |
| Paper 56 closed-immersion injectivity via per-sector η period map (2026-06-16) | 1 |
| Nested 2D adaptive quadrature for two-center overlap SUPPORT questions (Topos-3 attempt, 2026-07-05) | 1 |
| Graph counting/degeneracy measure as Born-rule generator (WH8 Step-1, 2026-07-01) | 1 |
| Overlap-slope law for the balanced-solver R_eq tilt (2026-07-06) | 1 |
| Löwdin retrofit of GeoVac integrals — the "sparsity-destroying option" (2026-07-07) | 1 |
| Non-orthogonal fermionic encoding to preserve Gaunt l-sparsity (2026-07-07) | 1 |
| "Exact up to one known seed" as the two-center ERI deliverable (2026-08-12) | 1 |
| Basis compactness as a quantum-resource advantage for Slater/Sturmian (QC-1, 2026-08-12) | 1 |
| Closing the *reducible* half of the three-centre ERI block first and deferring the hard half (Poly-2, 2026-08-14) | 1 |
| Dropping the three-centre ERI block to reach polyatomics with a two-centre engine (Poly-0, 2026-08-13) | 1 |
| Point-group symmetry (gerade lever) to rescue polyatomic ground-state metric conditioning (water C₂ᵥ, 2026-08-18) | 1 |
| Brute high-precision T2 by integrating the fibre over the family (2026-08-20) | 1 |
| ℚ(i) seam (P56 Dirac ↔ P59 T2) as a shared mechanism (2026-08-21) | 4 |
| Seam involution/fixed-point mechanism — τ=i pinned by a physical ℤ₂ (2026-08-21) | 1 |
| Low-rank/Woodbury flattening of the unique-center SW metric penalty (2026-08-21) | 1 |
| Analytic momentum-space factorization as a 1-norm lever ("analytic THC", 2026-08-21) | 1 |
| HO total-quanta conservation as an entanglement-rigidity mechanism (Papers 24 + 27, 2026-08-22) | 1 |
| Spectrum-aware QSVT interpolation as a scaling change for S^(−1/2) (2026-08-21) | 1 |
| TC three-body operator collapse via Gaunt/6j selection rules (2026-08-23) | 1 |
| Elliptic/CM-adapted basis for chemical accuracy (2026-08-23) | 3 |
| Symmetrising the TC operator to dodge the non-Hermitian quantum axis (2026-08-23) | 1 |
| Cheap + accurate electron-cusp on the Coulomb-Sturmian basis (2026-08-23) | 3 |
| Fixed / physics-determined GLOBAL γ for the cheap TC dressing (2026-08-24) | 1 |
| Position-dependent γ(r)∝n^{1/3} baked into the cheap non-Herm TC operator (2026-08-24) | 1 |
| Q̂₁₂ strong orthogonality as a quantum-encoding win for R12-CI (2026-08-26) | 1 |
| Per-pair γ / a second correlation length as the many-electron cusp fix (2026-08-26) | 1 |
| Explicit r₁₂ (James–Coolidge) on the PROLATE 2e basis with FD kinetic energy (2026-03-13/14; diagnosis CORRECTED 2026-09-18 — the analytical fix was already built and default, and the 94.7% it was to restore is an FD coarse-grid artifact; FD and analytical meet at ~86% at N=9) | 2 |
| Pairing explicit r₁₂ with the Neumann algebraic V_ee by calling the existing path (2026-09-18) | 1 |
| Hydrogenic per-n scaling (k_n = Z/n) as an L²-orthonormal Sturmian replacement (2026-08-26) | 1 |
| Two-block per-BLOCK radial exponent set for Paper 12's H₂ prolate CI (2026-09-18) — accuracy ~zero beyond the smallest basis (optimum is a plateau flat to 3 µHa); distinct from k_n=Z/n in costing no completeness; by-product = the "47× two-block conditioning penalty" was a RAW-cond artifact, two-block is 544× BETTER on the governing metric | 1 |
| R-adaptive shared bond exponent Z_orb(R) to heal the balanced-LiH R_eq drift (2026-09-20) | 1 | Non-variational in the parameter-free construction: h1=−Z_orb²/2n² locks Z_orb to Z_nuc, so Z_orb>1 on the Z=1 H centre drives E below exact (−9.6 vs −8.07). The drift's cure needs a variational (prolate-native) basis, not a knob on the composed recipe. See `debug/sprint_balanced_reqdrift_mechanism_memo.md`. |
| g_Vee speedup by algebra instead of sparsity: pair-space cancellation route + γ-series one-kernel route (2026-09-25) | 2 | Cancellation lever ≤9% of cost; γ-series converges exactly but ~100× dearer. The wall is sparsity (TREE4 dense field), cure = O(nnz) rewrite. See `docs/failed_approaches_ledger.md`. |
| Multi-geminal expansion / quadratic F² term for the LiH r₁₂ marriage gap (2026-09-26) | 2 | Gap was geminal SHAPE (Kato cusp), not count: one e^{−2r} geminal closes it; k≥2 saturates at +1.4 mHa; F² worth 0.013 mHa. See `docs/failed_approaches_ledger.md`. |
| r₁₂t² (e-e-n) geminal for the LiH residual (2026-10-05) | 1 | Separable and RI-free (validated vs MC) but worth ≤3 mHa: in the Li⁺ core u·t² adds 0.35 mHa; t² alone is orbital-representable radial correlation. Residual = frozen CI coefficients / orbital basis. See `docs/failed_approaches_ledger.md`. |
| GeoVac-native bricks (hydrogenic Z/n, Coulomb Sturmian) vs cc-pVDZ at equal size, water (2026-10-11) | 4 | Best native kit (two-scale Sturmian) 130 mHa above cc-pVDZ; hydrogenic fixed nodes cost 200-340 mHa on the O atom, Z_eff screening worse. Same Pauli count; 1-norm rises with accuracy. See `docs/failed_approaches_ledger.md`. |
| RE-INTERPRETATION of the CUSP-2 row (2026-09-25): its 'Z_c≈1.84 graph-validity boundary' was the non-physical (1/16)A hop; physical He has no Z_c. The negative result (Schwartz tail does not help) stands | 0 | See `docs/failed_approaches_ledger.md`. |
| RE-INTERPRETATION of the Dirichlet-disk row (2026-09-29): the 'boundary obstruction' was a float artifact at the exact grid point ρ=R (0/0 Markov ratio); off it the error decays. Negative WITHDRAWN; disk + cone OPEN | 0 | See `docs/failed_approaches_ledger.md`. |
| Generic three-electron chain fibre as K3 + elliptic + Tate periods (P3e-2, 2026-10-14; filed by `/walls` 2026-10-16) | 1 | 320-digit PSLQ at one mass triple finds no relation while the positive control recovers its exact one; the generic fibre carries more than K3 and elliptic periods. OPEN-LEANING: widen the basis (dilogs, elliptic products) before the computation. See `docs/failed_approaches_ledger.md`; walls register PERIODS-TRANSCENDENCE. |
| RE-INTERPRETATION of the Minnesota-Sturmian deuteron row (2026-10-20): the v2 '20x overbinding' came from singlet parameters labelled triplet and a non-orthonormal basis treated as orthonormal, not from Minnesota being HO-tuned. The polarizability negative stands | 0 | See `docs/failed_approaches_ledger.md`. |
| Frozen Li+ core as a Gaussian contact series matched to atomic levels only (CE, 2026-10-20; registration void on a shared control) | 1 | On LiH+ the error falls far less than 3x per order, the d channel makes it worse, and Paper 17's PK form transfers better: bound atomic levels fix only exterior phases while the molecular field probes the core interior. See `docs/failed_approaches_ledger.md`. |

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
| Paper 8-9 | Single-center molecular encoding | Sturmian Structural Theorem: in a shared-p₀ basis (H,S) share one SO(4) congruence → eigenvalues R-independent (no binding without R-dependent β_k(R)). **Scoped (v4.72.x, PI-reviewed):** the R-independence holds for the single-n D-matrix *approximation* (cross-block *imposed* = orthogonal D^(n); the backing test hard-codes it). The genuine cross-n Coulomb-Sturmian *overlap* couples n (computed: ⟨χ_1s\|χ_2s⟩=−0.47) and escapes the theorem — it binds (Paper 60's validated two-centre Sturmian H₂ CI; Avery SW closed forms; Herbst-Avery-Dreuw 2019 is atoms-only HF, cited for the cross-n overlap structure, not for molecular binding — corrected 2026-09-13) — but hits the qubit-encoding l-mixing wall (see §3 non-orthogonal-encoding row). Paper 8 §Remark + `test_paper8_overlap_cross_n.py`. | ANY investigation proposing to encode a heteronuclear molecule in a single-center (single-Z, single-k, single-p₀) basis, including nested hyperspherical, bond sphere, Sturmian CI, or unified orbital approaches |
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

The GeoVac label lattice is dimensionless and scale-invariant: its nodes are the Gel'fand–Tsetlin labels of the unit S³ harmonics, on which the Laplace–Beltrami operator is the diagonal SO(4) Casimir with pure-integer eigenvalues −(n²−1). (The lattice's edge Laplacian D−A is a different operator: l-preserving, n_max disconnected 2D-grid blocks, spectrum algebraic in [0,8]; it does not converge to Δ_{S³} — wall, 2026-09-25.) Physical energies emerge only through the energy-shell constraint p0^2 = -2E, which acts as a stereographic focal length. The 1/r Coulomb potential is not an input force law -- it is the coordinate distortion from stereographic projection (chordal distance identity). The matched scale kappa = -1/16 maps the graph's extremal eigenvalue (8 = 2 d_max, a 2D-grid band edge) to the hydrogen ground state -1/2 — an Observation; it does not map the graph spectrum onto the Rydberg ladder. Eigenvalues of the Laplace-Beltrami operator on unit S3 are pure integers: lambda_n = -(n^2 - 1).

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
| 1 | H (1-center, 1e) | S3 (Fock) | lambda_max -> 2 d_max = 8 (a spectral *bound*, not an accuracy: E_0 = kappa*lambda_max by construction; rate in Paper 0) | 7 |
| 2 | H2+ (2-center, 1e) | Prolate spheroid | machine precision (spectral; reference-limited) | 11 |
| 2 | H2 (2-center, 2e) | Prolate spheroid | near-exact D_e via explicit-r₁₂ CI with exact algebraic integrals (v5.14.10; registry `p12_r12_de_pct`); re-based CI just below (v5.13.9) | 12 |
| 3 | He (1-center, 2e) | Hyperspherical | 2D variational (raw, properly variational); cusp extrapolation (non-variational, marginally below exact); fixed-k=Z hydrogenic CI with the physical Hamiltonian (registry `he_gn_phys_n7_pct`); the graph-native figure is retired (a hop-augmented model, not He) | 13 |

*Level 3 note:* An earlier FD adiabatic result was lucky error cancellation (non-variational). The adiabatic coupled-channel solver converges to a structural floor (v2.0.8). The 2D variational solver (Track DI, v2.6.0) breaks this floor (tensor-product Laguerre × Gegenbauer basis), and a cusp extrapolation goes further but non-variationally (values: Paper 13). The 2D solver treats R and α simultaneously, capturing non-adiabatic R-α correlation that the adiabatic approximation misses. l_max convergence is monotonic; the per-channel Gegenbauer angular basis is the convergence bottleneck, not partial-wave truncation.
| 4 | H2 (2-center, 2e), hyperspherical treatment | Mol-frame hyperspherical | D_e via the hyperspherical treatment (Paper 15) — structural N-electron precursor to 4N/5; H2's *natural* geometry is Level 2 (prolate), see note | 15 |

*Level 4 note (added 2026-09-14, PI-directed):* **Level 4's placement for H2 was motivated in part by a claim that is now withdrawn, and on accuracy it is no longer the leading geometry for this system.** Paper 12 read its prolate-spheroidal residual as the electron-electron cusp being a coordinate singularity in (xi, eta) that demands a coordinate change; Papers 13 and 15 inherited that motivation. The residual was instead the sigma-only restriction — a phi-independent basis spans only m1 = m2 = 0, while a 1Sigma_g+ state constrains only the total M. Restoring the azimuthal channels in the *same* prolate-spheroidal basis recovers most of the gap, and re-conditioning that basis (v5.12.9, `geovac/prolate_recondition.py`) reaches chemical accuracy at the variational optimum (v5.13.9; registry `p12_rebased_de_pct_aopt`) — well above Level 4's Schwartz-cusp-corrected hyperspherical result; an independent Gaussian route reproduces the azimuthal-restoration step, and accurate grid-based prolate-spheroidal two-electron H2 methods exist in these same coordinates (the FE-DVR treatment of Tao-McCurdy-Rescigno, *Phys. Rev. A* **82**, 023423 (2010); a specific ground-state figure once quoted here could not be verified against that photoionization paper and was dropped 2026-09-20 — see CHANGELOG v5.14.11). What Level 4 retains, independently of the withdrawn claim, is *structural*: it is the exact N-electron generalization (SO(3N), S_N antisymmetry) that Level 4N and the composed Level 5 build on, and its cusp sits at a geometry-independent location. Those are architecture properties, not an accuracy advantage, and the hierarchy should not be read as asserting one here. **Decided (PI, 2026-09-17): H2's natural geometry is Level 2 (prolate spheroidal), alongside H2+, where its best result (explicit-r₁₂ CI with exact algebraic integrals, Paper 12 Sec. "Explicit correlation"; re-based CI just below) is obtained; H2 is now listed primarily at Level 2.** Level 4 is retained as the structural N-electron (SO(3N), S_N antisymmetry) generalization that Level 4N and Level 5 build on; its hyperspherical H2 result is documented as that structural treatment, not H2's leading geometry. See CHANGELOG v5.11.18 / v5.12.9, Paper 12 Sec. "Restoring the Azimuthal Channels" and "The monomial cap is conditioning, not a ceiling".

| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | equilibrium exists without PK, R_eq far off and D_e unbound at l_max=2 (2D variational) | 17 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq within several percent (Paper 17) | 17 |
| 5 | BeH₂ (polyatomic) | Composed (Level 3 + 4) + exchange | R_eq: Paper 17 | 17 |
| 5 | H₂O (triatomic) | Composed (Level 3 + 4) + lone pairs | R_eq: Paper 17 | 17 |
| Nested | Be (1-center, 4e) | S¹¹ (SO(12)), H-set | few-qubit encoding (Track DF record) | Track DF |

*Level 4N note:* Level 4N is the exact N-electron generalization of Level 4's mol-frame hyperspherical coordinates (SO(3N) replacing SO(6), S_N antisymmetry replacing the gerade constraint). Level 5 (composed geometry) approximates Level 4N, trading exact inter-group antisymmetry for a large angular compression via PK. The l_max=2 result demonstrates that equilibrium exists without PK (R_eq far off). Track AR (v2.0.23) confirmed that D_e overcounting was an adiabatic artifact: the 2D variational solver respects the variational bound and leaves LiH unbound, whereas the adiabatic solver overcounted D_e several-fold (values: Paper 17). The l_max=2 S₄ [2,2] angular basis is genuinely insufficient for LiH binding. Track AS (v2.0.23) confirmed composed encoding is categorically sparser than full N-electron encoding (an order of magnitude fewer Pauli terms and a much lower 1-norm; the Track AS counts predate the 2026-08-29 exact rule). Composed geometry's angular compression and structural sparsity are essential, not approximations of convenience.

The composed geometry (Level 5) is a fiber bundle: G_total = G_nuc semi-direct G_core(R) semi-direct G_val(R, core_state). Each electron group gets its own natural coordinate system, coupled via Z_eff screening and Phillips-Kleinman pseudopotential.

**Algebraic structure:** At every level, angular matrix elements are computed from quantum number labels and Wigner 3j symbols (via Gaunt integrals), with no spatial quadrature. The split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. At Level 2, the radial solver is fully algebraic for all m: σ states (m=0) use ordinary Laguerre three-term recurrence with zero numerical integration (v2.0.9); π/δ states (m≠0) use associated Laguerre basis L_n^{|m|}(x) with partial-fraction decomposition and Stieltjes integral recurrence, reducing non-algebraic content to a single transcendental seed e^a·E₁(a) (Track J, v2.0.10). At Level 3, the angular problem is fully algebraic, and the hyperradial overlap S and kinetic K matrices are algebraic via three-term Laguerre recurrence (Track H, v2.0.10: pentadiagonal M2 for S, tridiagonal derivative expansion for K, to machine precision, with an order-of-magnitude build speedup). The adiabatic eigenvalues μ(R) are proven transcendental (O(R) → O(R²) regime transition, v2.0.9 Track G) — the potential V_eff(R) must stay quadrature, point-by-point diagonalization is irreducible, though spectral radial solvers give two-order-of-magnitude speedups. At Level 4, the spectral Laguerre basis gives an order-of-magnitude dimension reduction for the hyperradial coordinate (Track I), and the spectral Jacobi basis an order-of-magnitude reduction for the angular sweep (Track K, a two-order-of-magnitude speedup). The combined spectral solver cuts the angular sweep from dominating the cost to about half of it. The angular eigenvalues μ(ρ) remain transcendental (computed by diagonalization at each ρ-point), but each point is now a small eigensolve. Values and the optimal radial basis size: the Track I/K records in CHANGELOG. Spatial quadrature enters for: Level 3 hyperradial potential matrix elements (V_eff transcendental), Level 4 angular eigenvalue sweeps (μ(ρ) transcendental, but the spectral basis reduces the cost sharply), Z_eff screening in composed geometries, and rho-collapse spline caching.

---

## 6. Paper Series

> Newcomer-facing status map with one-liners: `papers/INDEX.md`. Detailed per-paper notes (the long key-result descriptions, frozen at v3.110.0): `docs/paper_notes_archive.md`. Topic → paper lookup: `docs/topic_to_paper_lookup.md`. This section keeps only what agents need at dispatch time: loading tiers, the folder map, and live status flags.

### Loading tiers

**Always load** (framework identity): Paper 0 (packing axiom, K = −1/16) · 1 (spectral graph methods) · 7 (S³ proof, 18 symbolic proofs) · 14 (qubit encoding headline) · 16 (S_N periodicity) · 22 (angular sparsity theorem) · 23 (nuclear hub, Fock rigidity) · 24 (Bargmann-Segal S⁵, Coulomb/HO asymmetry) · 27 (entropy as projection) · 31 (universal/Coulomb partition) · 32 (the spectral triple; §VIII theorems).

**Load on topic** (full list and statuses in `papers/INDEX.md`): chemistry solvers → 8–9, 11, 12, 13, 15, 17, 19, FCI-A/M; QC resources → 20; QED/gauge/gravity → 2, 25, 28, 30, 33, 36, 41, 51; math.OA arc → 29, 38, 39, 40, 42, 44, 45, 50, 53 (43, 46–49 and 52 archived; 43's surviving content is Paper 42 §11, 52's a remark in Paper 50); foundations/periods → 18, 54, 55, 56, 57, 61; precision → 26, 34, 35.

**GUARDRAIL papers** — MUST load before any investigation in their domain (trigger words and protocol in §3.5): Papers 8–9 (single-center / unified-basis / shared-exponent molecular — Sturmian structural theorem), FCI-M (graph-concatenation molecular), Track DF record (nested hyperspherical).

### Folder organization (audience groups, reorganized 2026-05-22)

| Folder | Audience | Papers (.tex) |
|:-------|:---------|:-------------:|
| `papers/group1_operator_algebras/` | math.OA / NCG | 10 |
| `papers/group2_quantum_chemistry/` | quantum chemists | 12 |
| `papers/group3_foundations/` | mathematical physicists | 12 |
| `papers/group4_quantum_computing/` | QC / NISQ / VQE | 4 |
| `papers/group5_qed_gauge/` | HEP / gauge theory | 8 |
| `papers/group6_precision_observations/` | precision AMO | 4 |
| `papers/synthesis/` | cross-group narratives + field guide | 7 |
| `papers/archive/` | historical (3, 4, 5, 6, 10, 18v1, 21) + retired (43, 46–49, 52; `docs/retired_papers.md`) | 13 |

### Live status flags every agent must know

- **Paper 38 UNCONDITIONAL** (2026-06-10, translation-seminorm metrization; falsifiers `tests/test_p38_true_ch_window.py` + `tests/test_p38_action_seminorm.py` (smoothing only)). The WH1 keystone.
- **Paper 45 DESCOPED + partially rebuilt** (2026-06-09 K⁺ theorem withdrawn, falsifier `tests/test_p45_kplus_degeneracy.py`; 2026-06-10 product-carrier convergence restored in the action-seminorm framework, `prop:product_action_seminorm`, falsifier `tests/test_wh7_b1_joint.py` — signature-agnostic, NOT a Lorentzian claim). Do NOT cite pre-descope claims. **Paper 46 DESCOPED; Papers 47/48/49 PARTIAL** (in-paper Status notes; the norm-resolvent arrow and the TICI/cocycle algebra survive).
- **Paper 2 is an Observation**; the combination rule K = π(B + F − Δ) is labeled an Observation — never conjecture or derived (§13.5 hard prohibition; conjecture→observation downgrade 2026-06-14).
- **Paper 34** is the living projection catalogue (28 projections); **Paper 18 §III.7** is the master Mellin engine; tag every transcendental against both (memory rule).
- Papers are corrected **in place** (de-versioning directive 2026-06-10); git/Zenodo are the version record. No splinter files.

---

## 7. Code Architecture

> Full entry-point catalogue (~120 rows) and solver-method table: `docs/code_architecture.md` (live document — update there, not here). Most-used entry points:

| Task | Module · Entry point |
|:-----|:---------------------|
| Atomic lattice / Hamiltonian | `geovac/lattice.py` `GeometricLattice(max_n, nuclear_charge=Z)` · `geovac/atomic_solver.py` `AtomicSolver(max_n, Z)` (default = physical −Z²/2n² diagonal; the legacy graph model κ(D−A), whose excited states are non-physical, is `hamiltonian='graph'`) |
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
**Exception:** `-1/16` (κ, the matched energy scale — an Observation, Papers 0/7, not a derived constant) can be used directly.

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

Before forming or reporting a verdict — **or proposing a new connection or direction from corpus state** (e.g. an /aha generative pass) — on any question, *especially* when resuming a thread from a `debug/` memo, verify the CURRENT state, not the snapshot. The papers (the section that *owns* the question) + CHANGELOG since the memo's date are canonical; `debug/` memos are dated snapshots the corpus moves past, and a memo's "open question / next step" may already be closed. First move when picking up a thread = read the owning paper section + post-memo CHANGELOG, THEN conclude *or propose*. (Added 2026-06-14: twice in one session a verdict was formed from an 8-day-stale synthesis memo — an already-proven theorem was reported as a future sprint, and an already-settled A-vs-B question was re-answered wrong; both answers were live in Paper 56. Widened 2026-07-08: generation is not exempt — an /aha pass pitched a connection as novel that the Marcolli–vS chemistry arc had already established a month earlier; the `/aha` skill quarantines that verification to its Phase B (B1). Standing rule: `memory/feedback_verify_current_state.md`.)

### Branch QA Review Protocol

The corpus is QA'd branch by branch. Cycle, in order: **(1)** synthesis update; **(2)** adversarial paper review (overclaim, §1.5 rhetoric, zombie citations, status drift); **(3)** re-run 1–2 if the audit forces material synthesis changes; **(4)** adversarial *code* review, one `code-reviewer` per paper — map each claim to its backing test, **RUN the tests**, and audit whether the test actually *proves* the claim; **(5)** disposition. Run in **dependency order — trunk roots (Papers 0, 1, 7, 32, 38) before the branches**, so a finding at a root re-prices everything above it.

**Claim → artifact rule (hard).** Every load-bearing paper claim maps to a backing test, recorded in `docs/claim_test_matrix.md`. A claim with no backing test is a **coverage gap** — logged, and raised to the PI if load-bearing — never a silent omission. New equations follow the `test_paper{N}_*` convention.

**Retraction → dependents rule (hard, added 2026-09-04, PI-directed).** When a claim is retracted, corrected, or re-tiered, its `check_retracted_terms.py` entry must declare `cited_by`: the documents whose *argument* rests on it — distinct from `files`, which is only where its wording might appear. Each dependent carries a review stamp or `None`; unstamped dependents fail the gate, and a new entry declaring nothing fails too (`cited_by: {}` if genuinely none). **Stamp outcomes, not intentions** — an edit is not a review. *Why, measured:* eight of the ten recurring defect classes in the v5.4.4..v5.7.3 arc had one shape — a claim owned by one document, cited by several, corrected in the owner, left standing in the citers (“monotone from below” was corrected in Paper 0 and left false in Paper 7 *in the same commit*). Patterns cannot catch this: a citer restates the claim in its own words, so four successive pattern rebuilds each missed loci and one sweep reported clean while five loci survived a spelling difference. A dependency list is enumerated once, from the argument, and cannot be defeated by spelling. *Not covered:* the mirror direction, where an owner **strengthens** a claim and the citers keep its weaker form — a blocklist has nothing to say about a claim that is merely too weak; that needs its own registry and is owed.

**Guard-writing is a separate, separately-reviewed activity — not a step in the remediation (hard, added 2026-09-04, PI direction).** Do not write a regression guard in the same pass as the fix it protects. Complete the fix, then write the guard as its own piece of work, and review it by asking **what wrong answer would this accept?** rather than *does it pass?* Every new or rebuilt guard is fire-tested against the specific wrong answer its claim excludes (`python debug/qa/fire_test.py`), and a guard whose rejected answer cannot be named is not a guard. *Why, measured:* about ten guards in this corpus could not fail, found across DELTA #4, FULL #4, DELTA #5 and DELTA #6. Several were **replacements for guards that could not fail**, written by an author who had just been shown the failure mode; one was `test_cb_norm_is_one`, still `1 == 1` after two rebuilds; one was the `fire_test.py` selftest added to prove a cache fix worked, which could not detect that fix being reverted. The cause is structural rather than careless — a guard written in the same pass as the fix inherits the fix's reasoning, so it encodes the belief instead of testing it. Note the asymmetry in the same record: every substantive *mathematical* claim of that cycle survived independent re-derivation, while nearly every *guard* written to protect one did not.

**Disposition.** PM fixes small issues directly (status drift, cross-ref hygiene, precision, missing caveats); PM raises to the PI any load-bearing claim with no/weak/false-positive backing, a test that proves less than the prose, a suspected bug in a keystone, or anything touching a hard prohibition.

Full protocol + the three QA principles (provenance visibility / fresh adversary / two-way verdict): `docs/branch_qa_protocol.md`. The gate is `/qa` (PI-invoked only) against `docs/qa/criteria.md`.
### Summary-Surface Reading Rule (added 2026-09-12, PI direction)

**When a claim changes, reread the places that summarize it — in the same edit.**
Not a grep. A read.

1. **Every claim change:** reread the paper's **abstract**, **conclusion**, any
   paragraph labelled **Scope**, and the **Acknowledgments**. For a paper the
   size of Paper 60 that is ~6k tokens.
2. **Every session that touches claims:** read the **whole paper** once. ~30k
   tokens for Paper 60.
3. **A paper's synthesis moves with the paper.** If the claim is summarized in
   `papers/synthesis/`, that file is part of the change, not a separate errand.

*Why a reading rule and not another pattern.* The corpus has three mechanisms
for stale claims — the C16 phrase registry, `cited_by` dependents, and the
§13.8 `rests on:` edges — and **all three are document-granular**. The
2026-09-12 `/qa paper_60` FULL run found every defect **inside one file**: the
abstract claimed a 1953 theorem as derived here while the body carried
`[PRIOR ART]`; the conclusion confined a lever to the case the body had
breached; a retired energy survived as the abstract's validation endpoint under
a paragraph certifying the section converged. `cited_by` correctly reported no
dependent documents, because there were none. And no phrase registry can catch a
paraphrase — which is what a summary is, by construction.

*Why it is affordable.* Measured on Paper 60: abstract + conclusion ≈ 6.2k
tokens, whole paper ≈ 29.5k, against ≈1.5M for the review pass that found these
defects. **Reading the paper costs about 2% of reviewing it.** Cost is not a
reason to skip it.

*What it does not cover, stated so nobody over-trusts it.* Roughly half the
findings of that run were in summary surfaces; the rest were in body text, and
most of those were sentences contradicting themselves rather than stale copies.
Self-contradiction is a different failure and is caught by the reviewers'
internal-consistency mandate, not by this rule.

*The phrase registries stay.* They cost nothing to run and occasionally catch
something. They are a backstop, not the mechanism.

### Paper Retirement (added 2026-09-14, PI-directed)

**Retire claims, not papers.** The corpus's most-tested case is Paper 2: its
combination rule was demoted Conjectures -> Core -> Observations and then
relabelled conjecture -> observation, and the paper is still live and cited by
16 others. Tier demotion is the default move; retirement is the exception.

**Neither trigger is sufficient alone.**
- A withdrawn foundation is not grounds. Paper 45's main theorem was withdrawn
  and it remains load-bearing, with five healthy dependents.
- Zero dependents is not grounds. Papers 52 and 53 have none because they are
  DRAFT, which is unfinished, not supplanted, and needs finishing or dropping
  rather than archiving.

**The signal is the health of the dependents.** A paper propped up mainly by
other descoped papers is in a dying subtree. Two cheaper metrics were tried and
both failed, recorded so they are not rebuilt: plain citation count misses a
supplanted cluster entirely, because its members cite each other; and
"citations from outside the paper's own strongly connected component" collapses,
because 49 of 54 papers form one component. **C24** computes the surviving
metric and prints the profile.

**Archive, never delete.** Preservation is not the reason -- the Zenodo deposit
holds the PDF bytes, so the archival copy already exists independently of the
repo. The reason is that the repo copy is what future work greps, and pruning
has already cost this project real time
(`memory/feedback_resurrect_pruned_artifacts.md`). Deleting is not catastrophic;
archiving is simply free.

**Every retirement declares itself.** Move the `.tex` to `papers/archive/`, add
a row to `docs/retired_papers.md` (class, reason, and -- for anything with no
live successor -- trigger terms), keep the `papers/INDEX.md` row as a tombstone,
and stamp the documents that cited it, exactly as a retracted *claim* stamps its
dependents.

**Trigger terms are how the archive stays discoverable.** An archived paper's
approaches would otherwise be invisible to anyone starting new work, which is
the Sec. 3 re-derivation problem one level up. The answer is the same: record
the attempt with the phrase a future sprint would actually use, and let C24's
probe report it. No new habit, and no reminder to forget.

**The decision is the PI's.** C24 informs; it never nominates a verdict and
never blocks.

### Benchmarking Rule

After any modification to production code in `geovac/`:
1. Run `/regression` (default scope `touched`) — derives the test selection from `git diff` + import graph (consumer test files of every touched module) plus the topological-integrity baseline plus a small reproducible random sample. 30s–2min typical wall.
2. If the diff spans more than 2–3 modules, widen the `touched` selection — do **NOT** reach for `/regression full` at sprint close. **Measured 2026-08-31: the full scope is ~6.7 h serial, ~2.4 h with `-n auto`** (9,156 tests; 98 of 361 files over a 12 s budget; the cost is inherent to the Hamiltonian-building tests, not incidental). It is a **scheduled baseline**, not a close gate. Always pass `-n auto` when you do run it. See `debug/qa/test_suite_cost_memo.md`.
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

**Papers cite the permanent record, not transient `debug/` (policy, 2026-06-17 PI direction).** A paper must not cite a `debug/` sprint memo or script as backing. `debug/` is the transient clean-room dir (the bullets above) and is pruned over time, so a `\texttt{debug/...}` citation in a paper goes stale *by design* and leaves a dangling pointer (a reader following it hits a deleted file). The permanent homes are: the paper's own derivation, `CHANGELOG.md` (the sprint chronicle, §13.11), and `tests/` (the regression backing). The deterministic **C14** check (`debug/qa/check_file_refs.py`) gates the permanent code/artifact dirs (`geovac/`/`benchmarks/`/`demo/`) and reports `debug/` refs as **advisory** under this policy. *Standing debt (deferred sprint):* a corpus-wide sweep to neutralize the existing ~443 dangling `debug/` references (concentrated in papers 34/32/28; surfaced by the v4.21.0 C14 corpus sweep) is scheduled separately — the policy is adopted now, the ref-removal is its own task.

### Changelog Protocol

**Changelog granularity:** Each CHANGELOG entry should correspond to a releasable unit of work — a new feature, a completed diagnostic arc, a paper update, or a benchmark result. Do not create a new version entry for intermediate debugging steps, parameter sweeps, or exploratory runs. If a task requires multiple iterations to complete, document it as a single entry when the task concludes, noting the key findings and negative results. Intermediate data belongs in `debug/` with timestamped filenames, not in the CHANGELOG.

**Version numbering (revised 2026-08-22, PI direction).** **Default: bump the LAST number only (x.y.Z → x.y.Z+1).** The PM does not choose minor or major bumps; those are PI calls, made explicitly. The point of the change is to make version jumps *mean* something: when the second or first number moves, a reader should be able to infer that a significant corpus change happened, without having to read the CHANGELOG to find out. Under the old rule the PM bumped the minor version for every completed arc, which made minor bumps routine and therefore uninformative (v4.77 → v4.109 = 32 minor bumps).

So: patch (x.y.Z) is the standing default for everything — bug fixes, documentation, completed diagnostic arcs, paper updates, benchmark results, sprint closes. Minor (x.Y.0) and major (X.0.0) are reserved for the PI to call, and mark corpus-significant events: a retraction that moves published numbers, a change to the QA gate or the agent protocol, an architectural change, a reorganization of the paper series. If a sprint feels like it warrants more than a patch, say so in the session summary and let the PI decide — do not bump it unilaterally.

**Standing exception (added 2026-09-14, PI direction): any change to the `/qa` skill is a MINOR bump, automatically.** This one does not wait for a per-event PI call. `.claude/commands/qa.md` defines what the certification gate does, so a change there changes the meaning of every subsequent PASS — a reader seeing the second number move should be able to infer that the instrument itself moved, without reading the entry. Adding, removing or rescoping a deterministic check (C-criterion) counts; fixing a typo in the file does not. Granularity is unchanged: a diagnostic arc that tests 10 hypotheses and finds 9 negative results is ONE version entry, not 10.

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

**13.1 Architecture (moved).** Four layers: research agents (`agents/*.md` — Leader / Explorer / Decomposer / Reviewer, PI-invoked), plan mode, the PM session, opt-in worker sub-agents. Detail + cost policy: `docs/multi_agent_protocol.md`.

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

**13.4a Equation verification (moved).** **No equation enters a paper without a test that verifies it computationally.** Verification types, per-claim sufficiency bar, `test_paper{N}_*` convention: `docs/multi_agent_protocol.md`.

### 13.5 Hard Prohibitions

The following changes must NEVER be made by sub-agents or the PM agent:

- Any change to the natural geometry hierarchy (new levels, changed coordinate systems)
- Introduction of any fitted or empirical parameter
- Deletion or suppression of negative results from Section 3 or CHANGELOG.md
- Presenting the **combination rule K = π(B + F − Δ) (Paper 2)** as anything stronger than an **Observation**. The standing label is Observation — never "conjecture", "derived", "prediction", "theorem", or any stronger tier. The three ingredients B, F, Δ have derived spectral homes of distinct type; their *combination* is a numerical coincidence with no first-principles derivation (12 single-mechanism derivations eliminated). (History: Conjectures → Core 2026-04-18 Sprint A; Core → Observations 2026-05-02 per the curve-fit audit memo `docs/curve_fit_audit_memo.md`; conjecture → observation label downgrade 2026-06-14 per PI direction — "conjecture" was judged to carry unearned confidence that a derivation exists. Applies at the combination-rule level regardless of Paper 2's folder.)

**CLAUDE.md access control:** The PM may edit CLAUDE.md for mechanical updates that keep documentation in sync with code and paper changes. The PM may NOT edit sections that define strategy, framing, or the PM's own operating rules.

| Section | PM may edit? | Examples of allowed edits |
|:--------|:-------------|:-------------------------|
| 1 (Project Identity) | Version number ONLY | Bump v2.0.25 → v2.0.26 |
| 1.5 (Positioning & Framing) | **NO** | — |
| 1.6 (Project Phase) | **NO** | — |
| 2 (Development Frontier) | Yes | Update best results, add/complete track summaries, update backlog |
| 3 (Failed Approaches) | Yes (append only) | Add new failed approach rows; never delete or modify existing entries |
| 4 (Dimensionless Vacuum) | **NO** | — |
| 5 (Natural Geometry Hierarchy) | Numeric result cells ONLY (narrowed 2026-09-18, PI-directed) | Update a MEASURED result in the table's result column, or the same number where a note restates it (e.g. 99.77% → 99.81% when the optimum is re-measured). NEVER add/remove/renumber levels, change a coordinate system, or edit the framing of §5's notes — those remain prohibited by the §13.5 preamble above and §13.8 below. |
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

**Record cross-paper dependencies WHEN YOU WRITE THEM (hard, added 2026-09-04, PI direction).** When a paper claim you write or edit *rests on* a claim owned by another document, record that dependency in the same edit:

- if the owner claim already has a `check_retracted_terms.py` entry, add your document to its `cited_by` (see §9, retraction → dependents rule);
- if it does not — the common case, since most owner claims have never been retracted — note it in your claim's `docs/claim_test_matrix.md` row as `rests on: Paper N's <claim>`. Whoever later retracts that claim greps for it when building the entry.

*The operational test, because not every citation is a dependency:* **if the cited claim were withdrawn tomorrow, would this sentence have to change?** If yes, it is a dependency and gets recorded. If the citation is context, provenance, or courtesy, it is not.

*Why this sits at authoring time rather than retraction time:* at retraction time the dependency set has to be **reconstructed from memory**, and reconstruction is what failed — eight of the ten recurring defect classes in the v5.4.4..v5.7.3 arc were a claim corrected in its owner and left standing in its citers, including one corrected in Paper 0 and left false in Paper 7 *in the same commit*. When you are writing the dependent sentence you are looking straight at the thing you are relying on — that is the only moment the set is free to obtain. Patterns cannot substitute: a citing document restates the claim in its own words, which is why four successive pattern rebuilds each missed loci and one sweep reported clean while five loci survived a spelling difference.

**PMs must still NOT:**
- Introduce fitted or empirical parameters without PI direction
- Change the natural geometry hierarchy (new levels, changed coordinates)
- Delete or suppress negative results from Section 3
- Present the **combination rule K = π(B + F − Δ) in Paper 2** as anything stronger than an **Observation** (the standing label since the 2026-06-14 conjecture→observation downgrade; never "conjecture", "derived", or "theorem"). The prohibition is at the combination-rule level, not the paper-tier level: it applies regardless of whether Paper 2 sits in Conjectures, Core, or Observations. See §13.5.

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

**13.9b Slash commands (moved).** Genuine triggers: `/aha`, `/sprint-close`, `/checkpoint`, `/qa`, `/walls`. Force-fire backups for standing memory rules: `/audit-claim`, `/diag`, `/transcendental-tag`. `/qa` and `/walls` are **PI-invoked only** — never self-triggered. Each command's own file is `.claude/commands/<name>.md`; the full table is in `docs/multi_agent_protocol.md`.

**13.10 Research-agent integration (moved).** Research agents **do not write code or modify files** — they propose, the PM executes, the tests verify. Invocation + handoff format: `docs/multi_agent_protocol.md`.

### 13.11 Content Discipline and Token Efficiency

CLAUDE.md is loaded into every PM session and every sub-agent dispatch. Its size is paid repeatedly and every edit invalidates the prompt cache. The rules below keep this cost bounded. They are hard rules, not preferences — they override the impulse to "be thorough" in CLAUDE.md, because thoroughness has homes other than CLAUDE.md.

**Where each kind of content lives:**

| Content | Home | Format |
|:--------|:-----|:-------|
| Sprint chronicle (what happened, in detail) | `CHANGELOG.md` | Full prose, no length limit |
| Sprint summary (CLAUDE.md §2 entry) | CLAUDE.md §2 | ≤ 30 words: name + date + verdict + memo path |
| Dead-end record, full | `docs/failed_approaches_ledger.md` | The unabridged row. **Append here first.** |
| Dead-end record, index (§3) | CLAUDE.md §3 | Category row, or name + count in the per-sprint index |
| Superseded §2 bullets | `docs/development_frontier_archive.md` | Verbatim, by compaction round |
| Superseded WH status text | `docs/wh_register_history.md` | Verbatim (§1.7 rule 9: replace, never append) |
| Sprint memo | `debug/*.md` | ONE canonical memo per sprint, ≤ 5000 words |
| Surviving structural findings | papers/group*/*.tex | The papers, not CLAUDE.md, are the permanent record; papers cite permanent records (CHANGELOG / the paper / tests), **never transient `debug/`** (§9 policy, C14) |
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

10. **Compaction is relocation, never deletion (2026-09-01).** Content moves *verbatim* to its archive; the section keeps a pointer and enough stub that the item stays discoverable by name. The 2026-09-01 round took the file 220 KB → 99.7 KB and was verified by a conservation check (all 123 §3 rows verbatim in the ledger, all 121 dead ends still named, all 100 §2 bullets kept or archived). **Run that check on any future compaction** — one that silently drops content is worse than a large file.

11. **Apply rules 2–3 when writing the entry, not later.** At that round, 85 of 100 §2 bullets were over the 30-word budget and 76 were dated to a *single month* — never an age problem, just a rule never enforced in bulk.

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

---

## 15. Numeric Registry (the `\gvq` scheme)

**CLAUDE.md carries no measured values in its standing sections (PI direction 2026-10-02).** §1.5, §1.6, §2's results map and status paragraphs, and §5 state measured results qualitatively and point to the papers (authoritative), `docs/claims_register.md` and this registry. Reason: CLAUDE.md is loaded into every session, so a stale number here is copied into new prose with full authority (the README front page still carried the retired pair-diagonal figures at DELTA #25). The dated §2 one-liners, the §1.7 register and the §3 ledger are chronicle records and may quote the value of their date.

**Load-bearing numbers in the papers are linked objects, not free text.**

- **Artifact:** `debug/qa/numeric_registry.py` — every load-bearing quantity
  with its canonical value, the **convention** it is stated in, its
  **provenance** (measured by what route / cited from what source), and, if
  derived, the **expression** rather than a stored number.
- **Citation form:** `\gvq{registry_key}{literal}` in the `.tex`. The macro
  renders the literal only, so the PDF stays self-contained and archivable
  and the `.tex` stays hand-editable; the key is a foreign key.
- **Gate:** C21 (`debug/qa/check_numeric_consistency.py`) recomputes every
  derivation, checks every annotation against the registry (accepting
  correct display rounding), blocks retired values, and reports unregistered
  multi-document numerals.
- **Query:** `--index <key>` names every locus citing a quantity, plus the
  derivations it feeds.

### Why it exists

A census of the group4 + group6 numeric surface found **352 numerals at more
than one locus, 201 across more than one document**. Successive QA delta
passes were each finding 15–20 genuine defects — about **5% of that coupled
surface per pass** — because manual review *samples* the dependency graph
rather than traversing it. The classes that kept surviving were relational,
not value-level: **twins** (the same quantity tabulated in two papers),
**derivations** (λ/Q, N², ratios, exponents fitted from printed columns), and
**conventions** (identity-in vs identity-out; a fit's point range).

### Operating rules

1. **A value moves → edit the registry FIRST, then run C21.** The gate names
   every stale locus, including in other papers. Do not hand-sweep; that is
   the process this replaces.
2. **Then re-read the prose at each locus C21 names.** This is the step the
   scheme exists to force. A sentence around a number is a claim keyed to
   that number's *magnitude* and it does not update itself — "far below" was
   true at 5.8× and false at 0.92×; "the d-block is the sparsest" reversed
   outright. Updating a numeral without re-reading its sentence is the
   defect, not the fix.
3. **Never register a value you have not measured or cited.** `provenance`
   is not decoration; a registry of guesses launders them into authority.
4. **Annotation must not change rendered text.** Wrapping a literal is a
   no-op by construction; verify it (strip all `\gvq` and diff) rather than
   assuming it.
5. **Do not resolve values at build time.** Deliberately rejected: it would
   update numerals while leaving the surrounding claims silently false, and
   it would break the self-containedness of DOI-stamped PDFs.

Full rationale: the `numeric_registry.py` docstring and `docs/qa/criteria.md`
(C21). Self-test: `tests/test_numeric_registry.py`.
