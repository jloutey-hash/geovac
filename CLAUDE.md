# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v2.26.1 (April 18, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (unit S3) that is mathematically equivalent to the Schrodinger equation via Fock's 1935 conformal projection. This equivalence is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Authoritative source rule:** The papers in `papers/core/`, `papers/supporting/`, and `papers/observations/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it.

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers in `papers/core/`, `papers/supporting/`, and `papers/observations/` are written to academic standards but are not submitted to traditional journals. The project's viability case rests on producing a usable, benchmarked computational tool. Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks.

---

## 1.5. Positioning & Framing

GeoVac is a discretization framework that exploits the natural geometry of separable quantum systems to produce sparse graph Hamiltonians. It is **not** a new foundation for quantum mechanics — it works because the graph Laplacian converges to the known continuous Laplace-Beltrami operators in the continuum limit (Paper 7).

**Rhetoric rule:** The GeoVac framework is demonstrably conformally equivalent to standard continuous quantum mechanics at every level where it has been tested (Paper 7 provides 18 symbolic proofs for the S³ case; subsequent papers verify operator convergence for each natural geometry). The papers should present this conformal equivalence as the primary result. The discrete graph topology and the continuous Laplace–Beltrami operators are dual descriptions connected by proven conformal maps — the mathematics supports both readings, and the interpretation of which is more fundamental is left as a choice for the reader. Avoid language that asserts ontological priority of either description. Lead with concrete computational results (sparsity, scaling, accuracy, structural insight) rather than interpretive claims about the nature of quantum mechanics.

**Lead with concrete advantages:** O(V) sparsity, angular momentum selection rules baked into the basis, zero-parameter construction from nuclear charges and geometry alone, efficient qubit encodings (Paper 14). These structural properties are the framework's actual selling points — not philosophical claims about the nature of quantum mechanics.

**The π-free graph principle:** The GeoVac graph is π-free: all eigenvalues, degeneracies, and coupling coefficients are integers or rationals. Transcendental numbers (π, exponential integrals, spectral zeta values) enter exclusively when projecting the graph onto continuous manifolds for comparison with experiment or for computational convenience. This observation (Paper 18) motivates the design principle: stay on the graph whenever possible, and when projection is necessary, identify the minimal transcendental content (the exchange constant) required. The exchange constant taxonomy (intrinsic, calibration, embedding, flow) classifies projections by what determines them and predicts the computational cost of each departure from the graph.

**Quantum simulation positioning:** The classical solver investigation (v2.0.6-23) has comprehensively characterized what the natural geometry hierarchy can and cannot achieve for ground-state PES computation. The primary computational value proposition is now quantum simulation: the composed architecture produces qubit Hamiltonians with O(Q^2.5) Pauli scaling (vs O(Q^3.9-4.3) for Gaussian baselines), 51x-1,712x fewer Pauli terms across LiH/BeH2/H2O, and structural sparsity from Gaunt selection rules that is basis-intrinsic and compatible with all downstream optimizations (tapering, grouping, tensor factorization). The full N-electron quantum encoding comparison (Track AS, v2.0.23) confirmed that composed encoding is categorically sparser than direct grid encoding (334 vs 3,288 Pauli terms, 20x lower 1-norm), validating the composed architecture as the quantum computing approach. Market test (Track CA, v2.0.36): 1-norm comparison against computed Gaussian baselines shows GeoVac LiH electronic-only λ=33.3 Ha matches STO-3G λ=34.3 Ha (0.97×), with 2.7× fewer Pauli terms and 13× fewer QWC groups. He equal-qubit at Q=28: 6.8× lower λ than cc-pVTZ. Published DF/THC lambda values for Q<100 do not exist — GeoVac's competitive landscape is defined by raw JW baselines where it wins decisively. Position for VQE/NISQ (Pauli count, QWC groups) rather than fault-tolerant QPE (1-norm).

**General composed builder foundation:** Paper 16 (Chemical Periodicity as S_N Representation Theory) provides the group-theoretic foundation for the general composed builder's atomic classification. The structure types (A/B/C/D/E), the universal angular quantum number ν = N−2, and the recursive core-valence decomposition are all derived from S_N representation theory and implemented in `atomic_classifier.py`.

**PK limitation and research path:** The composed framework's PK pseudopotential is the accuracy bottleneck (5.3-26% R_eq error), and six modification attempts have failed (Section 3). Track CB (v2.0.37) confirmed that PK overcorrects: decoupled blocks (10.9% error) outperform PK (15.0%). Paper 19 documents the mathematical connection between GeoVac's S³ framework and the Coulomb Sturmian / Shibuya-Wulfman integral formalism. Track CD (v2.0.39-42) implemented the balanced coupled Hamiltonian with cross-center V_ne integrals from multipole expansion (convergence resolved: terminates exactly at L_max=2*l_max by Gaunt selection rules). 4-electron FCI: balanced coupled is the only configuration producing bound LiH (E=-7.924 Ha, 1.8% error; R_eq=3.226, 7.0% error; D_e=0.037 Ha). 878 Pauli terms (Conjecture 1 confirmed). Phase 4A with analytical V_ne integrals (incomplete gamma, machine precision): n_max=3 E(3.015)=-8.055 Ha (0.20% error). R_eq=3.280 bohr (8.8% error, structural drift +0.053 bohr/n_max, 3x smaller than PK). MIXED: energy converges excellently, R_eq drifts structurally. Optimal for single-point quantum simulation at fixed geometries. FCI accuracy characterization (Track CE, v2.0.44): balanced coupled LiH FCI PES at n_max=2 (7 R-points, Q=30) and n_max=3 (5 R-points near equilibrium, Q=84). Energy error 1.8% (n_max=2) → 0.20% (n_max=3), 9× improvement. R_eq 3.227 bohr (7.0%, n_max=2) and 3.28 bohr (8.8%, n_max=3). Gaussian resource comparison: composed 334 Pauli vs STO-3G 907 (2.7×), vs cc-pVDZ 63,519 (190×). Paper 20 written for quantum computing audience. Key file: `papers/applications/paper_20_resource_benchmarks.tex`. BeH₂ balanced coupled (Track CD, v2.0.41): 2,652 Pauli terms, 304.7 Ha 1-norm (0.86× composed — balanced is cheaper), 4-electron FCI 10.7% error vs PK's 27.4%. H₂O balanced coupled (Track CD, v2.0.42): 5,798 Pauli terms (7.45× composed), 1,509.3 Ha 1-norm (0.054× composed w/ PK, 4.18× electronic-only), 168 QWC groups. Non-collinear V_ne via Wigner D-matrix rotation. FCI infeasible (10e, Q=70). Lone pair V_ne well-behaved (all < 3 Ha/nucleus). H₂O deferred (non-collinear geometry needs general multipole) — RESOLVED v2.0.42. Second-row generalization (Tracks CH-CL, v2.1.0): frozen-core [Ne] treatment for Z=11-18 via Clementi-Raimondi Slater orbital exponents, analytical Z_eff(r) screening. 6 new molecules: NaH (Q=20), MgH₂ (Q=40), HCl (Q=50), H₂S (Q=60), PH₃ (Q=70), SiH₄ (Q=80). Pauli scaling Q^2.50 across 4 molecules (matching first-row Q^2.5). Ecosystem export for all. NaH FCI: overattraction at n_max=2 (no equilibrium), consistent with known balanced accuracy limitations. l=2 Wigner D-matrix rotation implemented (Track CM, v2.1.1) via general-l complex Wigner d-matrix converted to real SH basis — unblocks n_max=3 for all molecules. NaH n_max=3 convergence (Track CN, v2.1.1): Q=56, 5,349 balanced Pauli terms, 2-electron FCI energy improves 0.40 Ha over n_max=2, but PES overattraction persists (no equilibrium). HCl n_max=3: Q=140, 66,939 balanced Pauli. MgH2 n_max=3: Q=112, 34,473 balanced Pauli. SiH4 n_max=3: Q=224, 164,769 balanced Pauli. 1-norm cleanup (Track CO, v2.1.1): non-identity 1-norm (λ_ni) implemented in ecosystem_export — excludes identity term (nuclear repulsion + frozen-core energy, classical constants). NaH: λ_total=191 Ha → λ_ni=19 Ha (identity is 90%). Papers 14 and 20 updated with second-row data, corrected 1-norms, n_max=3 convergence. Key files: `geovac/neon_core.py` (FrozenCore), `geovac/atomic_classifier.py` (Z=11-18 entries), `geovac/shibuya_wulfman.py` (general-l rotation). Third-row extension (Tracks CQ-CT, v2.2.0): FrozenCore extended to [Ar] (18e, Z=19-20) and [Ar]3d¹⁰ (28e, Z=31-36) cores with Clementi-Raimondi exponents. Atomic classifier extended to Z=19-36 (transition metals Z=21-30 raise NotImplementedError). 6 new molecules: KH (Q=20), CaH₂ (Q=40), GeH₄ (Q=80), AsH₃ (Q=70), H₂Se (Q=60), HBr (Q=50). Isostructural invariance confirmed: frozen-core molecules with same block topology produce identical Pauli counts across rows 2-3 (NaH=KH=239, MgH₂=CaH₂=1501, etc.). Q^2.5 scaling universal across all 3 rows. Total library: 18 molecules. Ecosystem export for all 18 via `hamiltonian()` API. Multi-center molecules (Tracks CU-CX, v2.3.0): Extended composed/balanced architecture to molecules with multiple heavy-atom centers. New `OrbitalBlock.center_nucleus_idx` / `partner_nucleus_idx` fields map orbitals to specific nuclei. `MolecularSpec.nuclei` stores 3D nuclear positions. `_get_sub_block_positions()` updated with multi-center code path. 8 new molecules: LiF (Q=70), CO (Q=100), N₂ (Q=100), F₂ (Q=100), NaCl (Q=50), CH₂O (Q=120), C₂H₂ (Q=120), C₂H₆ (Q=160). Composed Pauli counts: LiF 778, CO 1111, N₂ 1111, F₂ 1111, NaCl 556, CH₂O 1333, C₂H₂ 1333, C₂H₆ 1777. Isostructural invariance extends to multi-center: CO and N₂ have identical Pauli counts (1111) despite different nuclear charges. Key finding: composed N_Pauli = 11.11 × Q across all 28 molecules (Q=10–160), coefficient stable to ±0.1. The composed Pauli count depends ONLY on Q, not on block count, structure, or atomic species. Total library: 28 molecules. Ecosystem export for all via `hamiltonian()` API. Key files: `geovac/molecular_spec.py` (extended OrbitalBlock/MolecularSpec), `geovac/composed_qubit.py` (8 new spec factories), `geovac/balanced_coupled.py` (multi-center V_ne wiring), `geovac/ecosystem_export.py` (8 new entries).

**Benchmarking rule:** When comparing to other methods, always use the strongest available baseline (cc-pVTZ or better for atoms, explicitly correlated methods for molecules), not just STO-3G. If the comparison is unfavorable, say so honestly and identify what the framework offers instead (sparsity, scaling, structural insight). Position the framework as a computationally principled alternative, not as a replacement for production quantum chemistry.

---

## 1.6. Project Phase

**Phase 1 (v0.9.x-v1.x): Foundation.** Graph Laplacian equivalence proof, natural geometry hierarchy, LCAO diagnostic arc, bond sphere theory, 18 symbolic proofs (Paper 7).

**Phase 2 (v2.0.0-v2.0.23): Classical solver investigation.** Systematic exploration of all solver architectures (adiabatic, coupled-channel, 2D variational) across all levels (2-5, 4N). 30+ completed tracks, 40+ documented negative results. Key outcomes: H2 at 96.0% D_e (Level 4), LiH R_eq at 5.3% (Level 5 composed), full N-electron equilibrium without PK (Level 4N). Investigation complete: all solver x PK x basis combinations exhausted, structural accuracy ceilings characterized, exchange constant taxonomy established (Paper 18).

**Phase 3 (v2.0.24-v2.7.0): Quantum simulation.** The composed architecture's structural sparsity (O(Q^2.5) Pauli scaling, block-diagonal ERIs) is the framework's primary computational advantage. Classical solver results serve as validation benchmarks for the quantum Hamiltonians. Next steps: quantum resource estimation, hardware-aware circuit compilation, experimental collaboration.

**Phase 4 (v2.7.0+): Nuclear extension and framework delineation.** Extension of the GeoVac framework to nuclear shell model Hamiltonians using the hyperspherical (HO + spin-orbit) basis, together with a precise delineation of which parts of the framework transfer to non-Coulomb systems. Key outcomes: (1) the potential-independent angular sparsity theorem (Paper 22) — ERI density depends only on l_max, not on V(r), verified at l_max=3 as 1.44%; (2) the deuteron and He-4 qubit Hamiltonians (Paper 23, Tracks NE/NF) using the Minnesota NN potential, Moshinsky-Talmi brackets, and a two-species tensor-product JW encoding; (3) the Fock projection rigidity theorem (Track NH) — the S^3 conformal projection is unique to -Z/r and does NOT transfer to HO or Woods-Saxon; (4) the composed nuclear-electronic deuterium PoC (Track NI) demonstrating the block architecture at 26 qubits with hyperfine validation; (5) the Bargmann-Segal lattice (Paper 24, Track NK) — discrete graph encoding of the 3D HO on the holomorphic sector of S^5, bit-exactly π-free in exact rational arithmetic at every N_max (verified at N_max=5: 56 nodes, 165 edges, zero irrationals); (6) the universal vs Coulomb-specific partition as the Phase 4 conceptual result, completed by the Coulomb/HO asymmetry analysis in Paper 24 (calibration π is structurally Coulomb-specific — tied to second-order Riemannian operators with nonlinear projections, has no analog for the HO's first-order complex-analytic projection). Track NJ memo definitively shelved the nuclear → alpha connection (originally "shelve," upgraded after Sprint 2/Paper 24 to "definitively shelved" based on the structural impossibility, not just empirical absence). HO rigidity theorem (Theorem 3 of Paper 24) is the structural dual of the Fock rigidity theorem.

---

## 1.7. Working Hypotheses (Internal Register)

This section is a **bold-claim register**, distinct from the rhetoric of the papers. Papers remain cautious under §1.5 (dual-description framing, no ontological priority, lead with computational results). This section is the register sub-agents and the PM are permitted to *reason from* when doing synthesis work. Nothing here appears in papers unless promoted after its falsifier clears.

The purpose is operational. When reasoning about what the framework *is* — what the transcendental taxonomy means, what α really is, where the RH thread actually points — the careful-paper register is under-tooled. Ramanujan worked this way: maximal-claim working notebooks, cautious proof letters sent to Hardy, both required. The working hypotheses below are the GeoVac notebook; Papers 0–30 are the letters.

**Governance:**
- PM may update a WH's "Status" field based on sprint evidence.
- Adding, retiring, or promoting a WH to paper-level claim requires explicit PI direction.
- Retired WHs are moved to an archive subsection with retirement rationale; never silently deleted.
- No WH is a license to bypass the rhetoric rule in papers or to skip the verification gates in §13.4.

---

**WH1 — GeoVac is an almost-commutative spectral triple.**
The framework is the data of a spectral triple (A, H, D): A = functions on the Fock-projected S³ graph, H = scalar/spinor state space (T1–T9), D = Camporesi-Higuchi Dirac operator (|λ_n| = n+3/2). Non-abelian gauge structure enters as inner derivations of an almost-commutative extension A × M_n(ℂ) — Paper 25 (n=1 / U(1)) and Paper 30 (n=2 / SU(2)) are the maximal-torus and full non-abelian slices of this structure. The K = π(B + F − Δ) formula has the *shape* of a finite-cutoff spectral-action coefficient: finite mode trace + regularized ζ sum + boundary mode count.
*Falsifier:* a GeoVac observable demonstrably inconsistent with any spectral-action expansion; or a structural feature of almost-commutative spectral triples (e.g., order-one condition, reality condition) that GeoVac demonstrably violates.
*Status:* **MODERATE-STRONG** (upgraded from MODERATE after supertrace sprint ST-1/ST-2, 2026-04-19). Supertrace sprint returned three paper-ready findings plus one clean negative: (F1) SD cancellation theorem — a_k^{D²}/a_k^{Δ_LB} = 4 = dim(spinor bundle) at every order k on unit S³, so the perturbative CC supertrace vanishes identically; (F2) Δ⁻¹ = 40 is the Euler-Maclaurin upper boundary term of the Dirac mode-count sum at n_CH=3; (F3) the (−) sign on Δ IS the standard (-1)^F boson-fermion grading (B,F scalar→(+), Δ spinor→(−)); (F4 CLEAN NEGATIVE) K/π is NOT the smooth-cutoff non-perturbative supertrace remainder (R_S − R_D/4 is always negative ~−1 to −3, never near K/π ≈ 43.6). Combined with Sprint A: α-PI POSITIVE PARTIAL identifies the external π in K = π(B+F−Δ) as the Hopf-bundle measure factor Vol(S²)/4 = Vol(S³)/Vol(S¹) = π, equivalently a₀² = π (squared Seeley-DeWitt zeroth coefficient); α-SP NO-MATCH-WITH-WEAK-APS-PARTIAL on direct CC-sign-rule derivation (standard CC bulk expansion gives (+,+,+) on S³, not the (+,+,−) of K; APS eta-invariant boundary subtraction is the only heat-kernel rule producing (−) next to bulk (+), and K's Δ is shape-compatible with APS but NOT literally an APS invariant since Δ⁻¹ = g_3^Dirac is a Dirac mode count); α-X PARTIAL-NEGATIVE on Paper 24 cross-check (π³ shape-match survives — every S⁵ SD coefficient factorizes as π³ × rational — but no S⁵ spectral-action coefficient produces a contribution at order α³ against 1/α, so the α-EB v2 π³α³ residual stays structural hint not derivation); α-LS MODERATE literature grounding established: Marcolli-van Suijlekom 2014 (J. Geom. Phys. 75, arXiv:1301.3480) "Gauge networks in noncommutative geometry" is a published framework of finite spectral triples on graph vertices with connection on edges whose spectral action is Wilson lattice gauge theory — structurally matches Papers 25/30 almost verbatim, with the Perez-Sanchez 2024 correction (arXiv:2401.03705, arXiv:2508.17338) clarifying that the continuum limit is Yang-Mills without Higgs. Each structural ingredient of WH1 has published precedent: finite ACG classification (Krajewski 1998, Paschke-Sitarz 2000, Ćaćić 2009), graph spectral triples (Marcolli-vS 2014, de Jong 2009), spectral truncations / operator systems (Connes-vS CMP 2021, Hekkelman 2022+2024), SM spectral action (Chamseddine-Connes 1997/2010). **However, no published framework predicts α from a discrete spectral action** — this remains GeoVac's original Paper 2 claim, not duplicated in the literature. Concrete follow-ups flagged: (a) Kluth-Litim 2020 (EPJ C 80:269, arXiv:1910.00543) gives explicit Seeley-DeWitt on spheres including odd dimensions — good independent check for α-X's S⁵ SD values; (b) Papers 25 and 30 should cite Marcolli-van Suijlekom + Perez-Sanchez correction; (c) Connes-van Suijlekom spectral truncations (CMP 2021) remains the concrete next-test target for Link 3.

---

**WH2 — Paper 18 is the Seeley-DeWitt + ζ-invariant decomposition of this spectral triple.**
The transcendental taxonomy is the structured output of spectral-action geometry, organized by operator order (2nd-order Laplacian → even-ζ via Jacobi-θ inversion; 1st-order Dirac → odd-ζ via half-integer Hurwitz; vertex-coupled → Dirichlet-L via parity characters) and bundle type (scalar, spinor, gauge). Every transcendental we have found so far sits in this three-axis grid.
*Falsifier:* a transcendental appearing in a GeoVac observable that cannot be placed in the grid.
*Status:* three of four axis-quadrants filled (Phases 4B-4I, Tiers 2-3, Sprint 4 QG/RH-J, Sprint 2 RH-M). Paper 18 v2.0 restructure should promote this from observation to theorem-with-evidence — that restructure is the operational consequence of this WH.

---

**WH3 — The lattice exists a priori; match to physics is evidence, not derivation.**
The framework originates in a geometric packing construction (pack Planck's constant → node counts give l, m in 2D → 3D lift via periodic-table row-length inductive proof gives n, s). This construction is independent of known physics. Its persistent match to physics — Fock S³ equivalence (Paper 7), nuclear magic numbers (Paper 23), Dirac fine structure (T8), Pauli sparsity scaling (Paper 14) — is therefore evidence that physics is hosted by a discrete spectral triple, not that the lattice was reverse-engineered from physics.
*Falsifier:* showing the packing construction, followed rigorously, does not force the (n, l, m, s) structure or does not produce a graph Laplacian whose spectrum is n² − 1; or, a match that turns out to depend on a hidden physics-informed parameter.
*Status:* origin-story framing permits strong ontological claim internally. Papers continue under §1.5 dual-description framing. This WH is what makes the Ramanujan-register appropriate.

---

**WH4 — The four-way S³ coincidence is one structure expressing itself four times.**
S³ appears as (a) the Fock projection image of the hydrogenic spectrum (Paper 7), (b) the base of the Hopf bundle carrying Paper 2's α construction, (c) the spin-carrier for the Dirac sector (T1, Camporesi-Higuchi), (d) the manifold of SU(2) gauge (Paper 30). Under WH1 these are not four coincidences — they are the single spectral-triple structure viewed in four projections. S³ = SU(2) is the unique rank-1 non-abelian compact Lie group where all four roles coincide.
*Falsifier:* a GeoVac construction that forces one of the four S³ roles to live on a different manifold without the framework breaking.
*Status:* strongest coincidence in the whole project; currently unnamed in papers. Motivates the spectral-triple framing more than any single computation.

---

**WH5 — α is a projection constant, not a derivable number.**
K = π(B + F − Δ) composes three structurally distinct spectral objects: finite Casimir trace (B), infinite Fock-degeneracy Dirichlet at the packing exponent (F = ζ(2)), Dirac boundary degeneracy (Δ⁻¹ = g_3^Dirac = 40). Nine mechanisms eliminated (Phases 4B–4I) confirm they do not share a common generator. α is therefore best read as the *conversion factor between three projection regimes that do not share a generator*, not as a quantity derivable from a single principle. "Why does the sum equal α⁻¹" is the correct open question, not "how do we derive each piece."
*Falsifier:* a spectral-triple construction (Connes-Chamseddine-style or otherwise) that derives K as a single coefficient of a well-defined functional.
*Status:* α-derivation program paused at this framing. Promotes Paper 2 from "conjectural derivation attempt" to "structural coincidence with three identified spectral homes." Sprint A (2026-04-18) strengthens WH5 along two axes: (i) α-MI confirmed K(m) is a single-point coincidence at m=3 — grows as Θ(m⁵), no series/asymptote, no convergent partial sum interpretation — AND found that m=3 is a TRIPLE coincidence, simultaneously satisfying (a) the B(m)/N(m) = dim(S³) = 3 selection principle, (b) the unique (+,+,−) sign pattern hit (3.5 orders of magnitude over next-best), and (c) agreement of two independent canonical Δ forms (boundary product |λ_m|·N(m−1) vs Dirac degeneracy g_m^Dirac/2) which agree only at m=3; (ii) α-SP decisively eliminated direct Connes-Chamseddine spectral-action derivation as a common generator. Both findings consistent with the projection-constant reading: no single principle generates K, but three independent selection mechanisms converge at one finite cutoff. (iii) α-EB v2 verified at 80 dps that the post-cubic residual R_predict = K − 1/α − α² = 1.2079×10⁻⁵ matches π³α³ = 1.2049×10⁻⁵ to within 0.25% (with π³ uniquely picked out among {π², π³, π⁴}); the structural reading is π³ = Vol(S⁵), which would connect Paper 2's residual to Paper 24's Bargmann-Segal S⁵ lattice as a new cross-paper hint. Next-order coefficient C ≈ 0.344 (closest to π/9 at 1.4%, then to 1/3 = 1/dim(S³) at 3.2%). Single CODATA data point cannot uniquely fix C. (iv) α-X cross-check against Paper 24 (2026-04-18) verdict: PARTIAL (negative-leaning). π³ shape-match holds at the Seeley-DeWitt coefficient level — every SD coefficient on round S⁵ factorizes as π³ × rational (a₀_scalar = π³, a₂_scalar = 10π³/3, a₄_scalar = 16π³/3; Dirac variants 4π³, −20π³/3, 14π³/3). Paper 24's π-free certificate is reconciled: it applies to the discrete graph and exact rational spectrum, while Vol(S⁵) = π³ lives in the continuum integration measure — no conflict. However, no S⁵ spectral-action coefficient produces a contribution at order α³ against 1/α structurally (CC expansion on d=5 is in odd powers of Λ with no log term, and standard CC uses a₀/a₄ ratios to fix α once, not α³). α-EB v2 finding therefore stays "structural hint, not derivation"; it should NOT be integrated into Paper 2 §IV.G as a substantive section. Recommended treatment: narrowly scoped footnote in Paper 2 Open Questions. Consistent with WH5's core thesis (α is projection constant, not derivable from a single principle). Paper 2 stays conjectural; Paper 24 unchanged.

---

**WH6 — GeoVac's RH-adjacent object is the Dirac spectral zeta D(s), not classical Riemann ζ.**
The empirical GUE-like zero statistics of D(s), D_even, D_odd (Sprint 3 RH-M, CV ≈ 0.35–0.40) are the framework's internal RMT phenomenon. The classical-RH bridge is closed by three independent walls: (i) zeros not on a single critical line (RH-M), (ii) no spectral-triple-natural functional equation (Sprint 4 RH-O, 48 orders of magnitude off), (iii) wrong Weyl law (Sprint 4 RH-N, γ_n ~ √n HO-class vs Riemann's log-density). The Weil dictionary is broken in GeoVac by design: Ihara side is Poisson (RH-G), spectral side is GUE. The right research object is D(s) itself, not a classical-RH identification.
*Falsifier:* D(s) zeros shown to lie on a single critical line at larger sample sizes; or a spectral-triple-natural functional equation discovered that closes the RH-O gap.
*Status:* three independent negatives against classical-RH, one independent positive for internal GUE. RH frontier paused (2026-04-18); if resumed, the target is D(s) and its spectral-action interpretation, not ζ_R.

---

## 2. Current Development Frontier

**Best results by system type:**
- Atoms: He at 0.019% error (2D variational on S⁵ + self-consistent cusp correction, l_max=7, n_R=35, Track DI); raw 0.022% at l_max=7; 0.004% with exact coalescence density; TC 2D cusp correction 0.001% at l_max=4, gamma=0.11 (non-Hermitian BCH, non-monotonic — sweet spot, not systematic l⁻⁸; gamma_opt ~ 0.51/(l_max+1.7), TC doubles convergence exponent ~l⁻¹ to ~l⁻², v2.9.0); graph-native CI 0.19% at n_max=7, 0.21% at n_max=8 (2,262 configs), 0.20% at n_max=9 (3,927 configs) with exact algebraic float integrals from `hypergeometric_slater.py` (zero grid error, zero free parameters). **CUSP-2 re-diagnosis (v2.9.2): the ~0.20% He graph-native CI floor at Z=2 is the small-Z graph-validity-boundary artifact, NOT the cusp — fit gives irreducible 6 mHa offset; the Schwartz cusp tail at l_max=6 is < 1 mHa.** CUSP-1 (v2.9.2): 2nd-order-RS w̃/δ screening compresses He n_max=7 CI 11× at full-CI floor (1218→109 configs) with zero accuracy loss.
- Negative ions: H⁻ (Z=1) found bound by graph-native CI at n_max≥2, but over-binds by 21% (E=-0.640 vs exact -0.528 Ha). Standard (non-graph) FCI is properly variational. Graph-native CI variational boundary at Z_c ≈ 1.84 (v2.9.0).
- Exotic atoms: PsH (positronium hydride, e⁻e⁺ in proton field) bound at 4.1% error (E=-0.756 vs exact -0.789 Ha, l_max=3) using Level 3 hyperspherical with sign-flipped charge function. Gaunt selection rules preserved. Alpha parity mixing essential (distinguishable particles). (v2.9.0)
- One-electron molecules: H2+ at 0.0002% (spectral Laguerre, Paper 11)
- Two-electron molecules: H2 at 96.0% D_e (Level 4 mol-frame hyperspherical + Schwartz cusp correction, l_max=6 σ+π 61 channels, Paper 15)
- Core-valence molecules: LiH R_eq 5.3% error with l-dependent PK at l_max=2 (composed geometry, Paper 17)
- Polyatomic molecules: BeH₂ R_eq 11.7% error (full 1-RDM exchange, composed geometry, Paper 17)
- Triatomic molecules: H₂O R_eq 26% error, uncoupled with charge-center origin (composed geometry, 5-block Level 3+4, zero parameters, Paper 17)
- Qubit encoding: O(Q^3.15) Pauli terms, O(Q^3.36) QWC groups, O(Q^1.69) 1-norm (Paper 14)
- VQE benchmark: 1.3x fewer Pauli terms at Q=10, 8.1x at Q=28 vs validated Gaussian cc-pVDZ/cc-pVTZ (Paper 14, v1.9.0)
- H₂ bond-pair qubit encoding: 112 Pauli terms at Q=10 (n_max=2), Q^3.13 scaling, R-independent sparsity, 1-norm=8.17 Ha (Paper 14, v2.0.27)
- Composed qubit encoding: consistent ~Q^2.2 within-molecule Pauli scaling across 6 composed molecules (LiH/BeH₂/H₂O/HF/NH₃/CH₄), 51x-1,712x advantage over published Gaussian baselines (Paper 14, v2.0.0-v2.0.30)
- Full N-electron molecules: LiH R_eq ~1.1 bohr at l_max=2 (63.5% error, 2D variational solver; unbound D_e confirms adiabatic overcounting was artifact) (full 4-electron mol-frame hyperspherical, PK-free, Level 4N)
- κ = -1/16 derivation (v2.26.1, April 2026): κ is derivable from the Fock projection, not merely fitted. The squared coupling between adjacent n-shells in the Gegenbauer eigenbasis is c²(n,l) = (1/16)[1 - l(l+1)/(n(n+1))]; for l=0 (s-wave), c²(n,0) = 1/16 universally for all n. Three equivalent readings: (1) squared Chebyshev transition amplitude (1/4)², (2) inverse Fock Jacobian 1/Ω⁴(0) where Ω(0) = 2 is the stereographic conformal factor at p=0, (3) l=0 base rate of the Casimir decomposition. Paper 18 reclassified κ from "calibration" to "conformal." Notable: c²(4,3) = 1/40 = Δ (Paper 2 boundary term). Four negative probes confirmed κ doesn't live in graph topology (Ollivier curvature identically zero, Ramanujan κ-independent, spectral action no extremum, S⁵ has no κ — Coulomb-specific via Fock rigidity). Key files: `debug/probe_k{1..5}_*.py`, `debug/data/probe_k{1..5}_*.json`, `debug/probe_kappa_sprint_memo.md`.
- Fine structure constant: alpha from Hopf bundle at 8.8x10^-8, zero free parameters, p-value 5.2x10^-9, universal algebraic identity B_formal/N = d, Hopf generalization negative result, circulant Hermiticity, second selection principle (Paper 2, conjectural; Phase 4 sharpening via Fock rigidity theorem in the $S^3$ specificity section); Phase 4B-4G structural decomposition (April 2026): B = 42 = finite Casimir truncation (κ↔B Fock-weight link, α-C), F = π²/6 = D_{n²}(d_max) = ζ_R(2) infinite Fock-degeneracy Dirichlet series at the packing exponent (α-J), Δ = 1/40 = |λ_{n_max}|·N(n_max−1) = c²(4,3) Fock coupling at the angular-momentum edge of the cutoff (α-K + κ sprint), three-tier composition without common generator — combination rule K = π(B + F − Δ) is now structurally decomposed but remains conjectural at the level of why the sum hits α⁻¹
- Spectral-action supertrace sprint (v2.26.1, April 2026): three paper-ready findings on K = π(B + F − Δ) in the Connes-Chamseddine framework. (F1) SD cancellation theorem: a_k^{D²}/a_k^{Δ_LB} = 4 = dim(spinor bundle) at every SD order k on unit S³ — the perturbative CC supertrace Str[f(D²/Λ²)] vanishes identically. (F2) Δ⁻¹ = 40 is the Euler-Maclaurin upper boundary term of the Dirac mode-count sum Σ g_D(n) at n_CH = 3: the EM formula gives upper boundary g_D(N)/2 = 2(N+1)(N+2) = 40 at N=3. (F3) The (−) sign on Δ is the standard (-1)^F boson-fermion grading: B and F are scalar (bosonic, +), Δ is spinor (fermionic, −). (F4 CLEAN NEGATIVE) K/π is NOT the smooth-cutoff non-perturbative supertrace remainder — R_S − R_D/4 is always negative (~−1 to −3) at all natural Λ², never near K/π ≈ 43.6. Per-shell Casimir formula verified: c(n) = n²(n²−1)/2. Key files: `debug/st_supertrace_probe.py`, `debug/st_nonperturbative_probe.py`, `debug/data/st_supertrace_probe.json`, `debug/data/st_nonperturbative_probe.json`, `debug/st_supertrace_sprint_memo.md`.
- Nuclear systems (Phase 4, Paper 23): deuteron 16 qubits / 592 Pauli terms / 227 MeV 1-norm (Track NE, Minnesota NN potential, two-species JW); He-4 16 qubits / 712 Pauli terms (Track NF, Pauli count grows only 1.20x for 12.25x larger Hilbert space); Mayer-Jensen magic numbers 2, 8, 20, 28, 50, 82, 126 from HO + spin-orbit + Nilsson l(l+1) at v_ls/hw = 0.171, d_ll/hw = 0.021 (Track NB)
- Composed nuclear-electronic deuterium (Track NI, 26 qubits): 614 Pauli terms (592 nuclear + 10 electronic + 12 hyperfine cross-register), coefficient ratio ~2e13 across nuclear/electronic/hyperfine scales, hyperfine singlet-triplet gap validated at 3*A_hf/4 = 1.62e-7 Ha (21cm line). Practical single-pass quantum simulation requires block-partitioned solving due to the 10^13 dynamic range.
- Angular sparsity theorem (Paper 22, promoted to core): ERI density depends only on l_max, not on V(r). Verified values: l_max=0: 100%, l_max=1: 7.81%, l_max=2: 2.76%, l_max=3: 1.44%, l_max=4: 0.90%, l_max=5: 0.62%. Universal across Coulomb, harmonic oscillator, Woods-Saxon, square well, Yukawa.
- Fock projection rigidity theorem (Paper 23, Track NH): the S^3 conformal projection p_0 = sqrt(-2E_n) maps a one-electron central-field Hamiltonian onto the free Laplacian on S^3 if and only if the spectrum is l-independent within each n-shell. Unique to -Z/r by SO(4) symmetry. Defines the universal/Coulomb-specific partition: angular sparsity is universal, S^3 machinery and Hopf bundle are Coulomb-only.
- Bargmann-Segal lattice (Paper 24, Track NK): discrete graph encoding of the 3D HO on the holomorphic Hardy-space sector of S^5. Built from SU(3) (N,0) symmetric reps via the Bargmann transform; nodes are (N,l,m_l), edges are SU(3) dipole transitions (ΔN=±1, Δl=±1), edge weights are exact-rational squared matrix elements. Bit-exactly π-free in exact rational arithmetic at every finite N_max (verified at N_max=5: 56 nodes, 165 edges, no irrationals). Hamiltonian is diagonal: ℏω(N+3/2). The graph adjacency encodes dipole transitions only, NOT the spectrum (in contrast to the Coulomb S^3 case where (D-A) computes the spectrum). HO rigidity theorem (dual of Fock rigidity): the 3D isotropic HO is the unique central potential whose spectrum arises from the Euler operator on the Hardy space H^2(S^5) restricted to (N,0) SU(3) irreps. Coulomb/HO asymmetry is structural: first-order complex operators give linear spectra and linear projections (no transcendentals); second-order Riemannian operators give quadratic spectra and nonlinear projections (introduce calibration π).
- S^5 gauge-structure extension (Sprint 5 Track S5, April 2026, MIXED verdict): answers Paper 25 §VII.1's open question. (a) **Abelian U(1) Wilson–Hodge structure transfers verbatim** to the Bargmann-Segal graph: at N_max=5, V=56, E=165, c=1, β_1 = E−V+c = 110; signed incidence B, node Laplacian L_0 = BB^T, edge Laplacian L_1 = B^T B, Hodge decomposition, and SVD identity between nonzero spectra of L_0/L_1 all hold; node-local ψ_v → e^(iχ_v)ψ_v acts covariantly on dipole edge amplitudes. (b) **Non-abelian SU(3) analog is NOT natural**: transitions between (N,0) and (N+1,0) irreps are Clebsch–Gordan intertwiners, not SU(3) group elements, so Wilson lattice gauge theory on the (N,0) tower fails the fixed-group-on-every-link requirement; the natural gauge group of the Bargmann graph is U(1), not SU(3). (c) **m_l-fiber-collapse quotient is NOT a CP^2 discretization**: the 12-sector quotient at N_max=5 has Laplacian eigenvalues {0, 2.22, 4.87, 4.94, 11.65, 12.11, 18.74, 28.09, 33.58, 50.61, 63.77, 99.43} that do NOT match the CP^2 Fubini–Study scalar Laplacian sequence λ_k = 4k(k+2) (best single-scale rescaling has 25% max residual, ratios vary 0.08–0.19 with no uniform fit). Sharpens Paper 24's Coulomb/HO asymmetry from two layers to THREE: (i) spectrum-computing role of L_0 (yes S^3, no S^5); (ii) calibration π (yes S^3, no S^5); (iii) Wilson physical content (full S^3, reduced S^5 — the combinatorial vocabulary is universal but there is no matter-propagator role for L_0 on S^5). Paper 2's K = π(B+F−Δ) has NO natural S^5 analog: candidate F-analog = (1/2)[ζ_R(4)+3ζ_R(5)+2ζ_R(6)] (mixed even/odd zeta, no clean π^2/6), and B/Δ analogs are Coulomb-projection-tied (Paper 23 Fock rigidity restricts calibration π to S^3). Paper 24 Corollary (§V.D) and Paper 25 §VII.1 updated. Data: `debug/s5_bargmann_segal_graph.py`, `debug/s5_edge_laplacian_analysis.py`, `debug/data/s5_graph_spectrum.json`, `debug/s5_gauge_structure_memo.md`.

**Classical solver status: INVESTIGATION COMPLETE (v2.0.24).** The systematic exploration of all classical solver architectures across the natural geometry hierarchy is complete. Over 30 investigation tracks (v2.0.6-23) exhausted all solver (adiabatic, coupled-channel, 2D variational) x PK (channel-blind, l-dependent, PK-free) x basis (FD, spectral, algebraic) combinations. Structural accuracy ceilings are characterized at every level: H2 at 96.0% D_e with CBS ~97% (Level 4), LiH R_eq at 5.3% with structural drift +0.15 bohr/l_max (Level 5 composed), full N-electron equilibrium without PK at 63.5% R_eq error (Level 4N, angular basis limited). The exchange constant taxonomy (Paper 18) classifies where transcendental content enters at each level. Remaining classical accuracy improvements require either higher l_max (diminishing returns, increasing cost) or fundamentally new coordinate systems (none identified). The composed architecture at l_max=2 represents the production operating point for classical PES.

**Entropy vs V_ee off-diagonal mass power law (v2.9.2, EP-2c/d/e/f):** Paper 27 §VII.B Prediction 2 verified, reframed, and bounded. (EP-2c) The correct predictor is **dimensionless** w̃_B = w_B / ‖H_1‖_F (raw w_B gives sign inversion from Z²/Z scaling mismatch). He-like pilot at n_max=3 across Z∈{2..10}: S_B = 8.16·w̃_B^2.383, R²=0.998. Z-scaling α_Z=-2.613 matches Paper 26's -2.56. (EP-2c multi-block) 16 single-center 2e blocks across LiH/HF/H₂O/NH₃ cores + lone pairs: combined fit S_B=7.79·w̃_B^2.374, R²=0.9983; max 7.8% deviation from He-like reference. (EP-2d) The apparent 40% core-vs-lone-pair A offset is a **Z-range fit artifact, not a block-type invariant** — identical He-like family reproduces the offset exactly when restricted to the same Z subsets. (EP-2f) Extended to Z=100: local log-log slope drifts monotonically from 2.64 (Z=2) to 1.92 (Z=100), continuing to decrease. Asymptote is close to but NOT identical to second-order RS α=2; residual ~-0.08 gap persists. The S_B(w̃_B) relation is a smooth curve, not a pure power law. Band fits: Z≤4.5 α=2.59, 5≤Z≤10 α=2.18, Z>10 α=1.95 (all R²>0.999). (EP-2e) **Bond blocks break the single-center curve**: LiH bond (10 orbitals, two-center H_1 via balanced builder) gives S=0.079 at w̃_B=0.042, **20× above** the single-center prediction 3.9e-3. Mechanism: two-center H_1 has near-degenerate Li-1s/H-1s pair → V_ee mixes at O(1) amplitude (quasi-degenerate perturbation theory), not V/ΔE perturbative. Kinetic entropy S_kin=0 preserved. (EP-2g) **Universal two-variable collapse**: adding dimensionless H_1 gap δ_B = ΔE_1/‖H_1‖_F gives S_B = 0.157·(w̃_B/δ_B)^2.228, R²=0.990 across 22 points (He-like Z∈[2,100] + LiH bond R∈[1.5,8] bohr). Single-center γ=2.11, bond γ=2.37 — same form within Δγ≈0.26. Saturation at log(2) for dissociating bonds (R=6,8) drives the 76% max deviation (physical two-site limit, not model failure). Sharpens Paper 27 §II: entropy is generated by V_ee off-diagonal between H_1 eigenspaces with amplitude controlled by w̃_B/δ_B. (EP-2h) Extension checks: (a) BeH₂ bond block at equilibrium already saturated at log(2), H₂O bond consistent with single-center curve within 2× — universal curve holds only for non-saturated bonds. (b) n_max scan Z∈{2,3,4,6,10}: γ(n_max=2)=2.63, γ(3)=2.38, γ(4)=2.32 — monotone decrease but does NOT reach 2 (residue γ-2≈0.32 persists). (c) Bounded form log(2)·tanh²(c·(w/δ)^γ) fails globally (R²=0.37 vs pure power law 0.90) — log(2) saturation is system-specific, not a generic fit feature. Combined 21-point universal fit S=0.166·(w/δ)^2.198, R²=0.897. (EP-2j) **Paper 27 §II non-degeneracy clause is essential**: Li (3e ²S, single SD GS) gives S_kin≈0.64 from open-shell antisymmetrization (not V_ee correlation; ΔS=S_full−S_kin=0.01 is the actual correlation contribution). Be (4e ¹S) gives S_kin=1.23 from genuine 3-fold H_1 degeneracy (2s/2p shell unsplit by κ·adj). Strict spatial-1-RDM floor S_kin=0 holds only for closed-shell 2e systems exhausting a complete (n,l)-degenerate sub-shell — He qualifies, Be doesn't. Paper 27 §II has new "Scope of the proposition" subsection clarifying. (EP-2k) **HF/H₂O bond R-sweeps are R-independent** (0.4%/0.2% variation across R∈[1.3,3] bohr); both sit at constant 1.20×/1.25× above EP-2g universal curve at every R. LiH is the *symmetric-Z* outlier; chemically realistic heavy-atom bonds are dominated by the heavy atom's intrinsic structure. Universal curve holds within 25% R-independent band for heavy-atom bonds. (D4) Paper 26's S~Z^-2.56 reproduces from Paper 27 machinery at n_max=4 (measured -2.546, R²=0.995). Paper 26 Sec III is operator-theoretically derivable from Paper 27 Sec VI.B; cross-reference added to Paper 26 + bibitem. (D5) Paper 24 extended with §V.C (HO entanglement rigidity corollary) tying EP-2b to the HO rigidity theorem. **Paper 27 consolidated**: 1111→900 lines, sprint chronology replaced by 2 clean subsections (HO theorem + universal scaling), abstract rewritten to 5 results. (EP-2i) Deep diagnostics resolve all three EP-2h surprises: (a) BeH₂ saturation is NOT intrinsic — composed (within-block-only) BeH₂ bond gives S=0.028 at w/δ=0.284 (within 3× of curve); balanced builder's cross-block ERIs compress δ from 0.204 to 0.027 (8×) triggering log(2) saturation. LiH composed bond is exactly degenerate (δ=0) because Li-Z_eff=H-Z=1 — only balanced cross-center V_ne breaks this. "The bond block" is not a single well-defined object. (b) H₂O factor-2 offset is INTRINSIC: S_balanced/S_composed=1.002, unaffected by cross-block ERIs. (c) γ(Z,n_max) 2D surface: local slope Z=15→30 approaches 2 from above with mild n_max downshift: 1.992 (n_max=2), 1.976 (3), 1.966 (4). The EP-2f observation of 1.92 at Z=100/n_max=3 is this surface's n_max-ridge. **True asymptote γ_∞ ≈ 1.96** (Richardson extrapolation on n_max=2,3,4,5), slightly **below** 2nd-order Rayleigh-Schrödinger γ=2. The persistent gap is attributed to multi-shell aggregation in the basis-extension limit. (EP-2L) n_max=5 local slope at Z=15→30 is 1.959, continuing the monotone decrease 1.992→1.976→1.966→1.959. (EP-2M) Two-scale bounded forms (rational saturation, stretched exponential, Hill function) all give R²≈0.985-0.989 in log-space, statistically indistinguishable from pure power law (0.990) — the log(2) saturation is *discrete* (only 1-2 LiH dissociation points) not a smooth crossover; bounded ansatze can't earn their parameter cost. (EP-2N) Be analytical degenerate-PT: 3×3 V_ee diagonalization in the H_1 GS subspace gives GS = 0.98|2s²⟩ - 0.14|2p_0²⟩ - 0.14|(2p_-1 2p_+1)_S=0⟩, occupations (2.0, 1.92, 0.02, 0.04, 0.02), S_full=0.794 nats — matches FCI EP-2j S_full=0.902 qualitatively. V_ee acts as degenerate-PT lifter selecting 2s² over 2p², leaving ~2% p-orbital leakage. Key files: `debug/ep2{c,d,e,f,g,h,i}_*.py`, `debug/data/ep2{c,d,e,f,g,h,i}_*.json`, `tests/test_paper27_entropy.py` (22/22 passing). Key files: `debug/ep2c_entropy_vs_vee_offdiag.py`, `debug/data/ep2c_entropy_vs_vee_offdiag.json`.

**HO two-fermion entanglement (v2.9.2, EP-2b):** Paper 27 §VII.A Prediction 1 verified, result stronger than predicted. On the Bargmann-Segal HO lattice with the Minnesota NN singlet kernel at ℏω=10 MeV: S_HO = 0 **exactly** at the ground state (not merely ≪ S_Coulomb), because [H_HO, V_central] = 0 to machine precision (relative commutator 4–8×10⁻¹⁶ at N_max=2,3). Mechanism: Moshinsky-Talmi transformation makes total HO quanta N_tot = N_rel + N_CM a conserved quantum number of any central two-body V, so the ground state is confined to the lowest-N_tot sector. At M_L=0, M_S=0, two-fermion closed-shell, that sector is 1-dimensional: GS = |(0s)²⟩ single Slater determinant, occupations (2,0,0,0) exactly. E_full = 3ℏω + V_00,00 ≈ 14.898 MeV (basis-size-independent at N_max=2,3). Sharpens Paper 27 §II: entanglement is generated only by the V-off-diagonal component *between H_1 eigenspaces containing the GS*, not by total V Frobenius mass (HO's V has ~93% diagonality in the H_HO eigenbasis but 0% cross-eigenspace coupling). Key files: `geovac/nuclear/ho_two_fermion.py` (new module), `debug/ep2b_ho_two_fermion.py`, `debug/data/ep2b_ho_two_fermion.json`, `tests/test_paper27_entropy.py`.

**Cusp attack sprints (v2.9.2, CUSP-1/2/3):** Three concrete attacks on the He cusp ceiling, all completed. Net: **the cusp ceiling for He graph-native CI at Z=2 is not what we thought it was.** (CUSP-1) w̃/δ-based 2nd-order-RS screening of He n_max=7 graph-native CI (1218 configs → ~100 configs via perturbation-theory scoring in the H_kin eigenbasis) achieves **11× config compression at full-CI floor, 22× at floor+0.01pp** with zero accuracy loss. The cusp lives in ~10% of configs (top-scored states are (1s,ns), (2s,2s), (2p_-1, 2p_+1)_singlet, (2p_0²) — exactly the cusp-adjacent set). Useful for VQE/qubit compactness; doesn't break the floor because the floor is structural. (CUSP-2) **Major re-diagnosis**: fitting He graph-native CI err_abs(l_max) = -A/(l+2)⁴ + c gives **irreducible floor c = 5.99 mHa**, with the Schwartz cusp tail at l_max=6 being < 1 mHa. The 0.20% ceiling long attributed to "cusp" is actually the **small-Z graph-validity-boundary artifact** (Z_c ≈ 1.84, He is just above). Confirmed at Z=10 where the error sign flips (under-binding +85 mHa, consistent with real basis truncation). **Implication: for He at Z=2 graph-native CI, cusp corrections cannot help because cusp is not the limiting error.** Real cusp work should be done at Z≥4 where the boundary artifact is suppressed. Schwartz extrapolation remains valid for Level 4 (Paper 15) molecular solvers where there's no small-Z graph issue. (CUSP-3) **TC structurally dead confirmed at high n_max**: particle-number-projected FCI benchmark at He composed n_max=2→5 shows TC error plateaus 3.40→3.47% (ratios 0.987→0.994→0.997, asymptoting to 1) while Standard continues to improve 2.48→2.02. Pauli ratio TC/Std grows 1.68→2.73. TC is worse at every accessible n_max and getting relatively more expensive. **TC is now decisively dead in second quantization on the composed basis**, not just at small n_max. Section 3 TC row updated with extended-n_max evidence. Key files: `debug/cusp{1,2,3}_*.py`, `debug/data/cusp{1,2,3}_*.json`.

**Cusp characterization (v2.9.0):** The electron-electron cusp has been fully characterized as a matrix-algebraic object: (1) h1 converges by n_max=2, diagonal V_ee (mean-field) by n_max=3 — ALL remaining convergence is off-diagonal V_ee correlation; (2) V_ee is FULL RANK in the graph eigenbasis at every n_max, with no low-rank shortcut — the cusp distributes broadly across all angular channels; (3) the TC similarity transformation doubles the convergence exponent (~l⁻¹ to ~l⁻²) with optimal gamma decreasing as ~0.51/(l_max+1.7), confirming the TC operates as a finite-basis correction (not true cusp removal). Key data: `debug/data/track_di_tc_2d.json`, `debug/data/tc_gamma_scan.json`. Energy-graph follow-up (v2.9.1) confirmed this from pure graph topology: in the wavefunction-graph eigenbasis, the unsquared Frobenius ratio ‖diag_{H₁}V_ee‖_F / ‖V_ee‖_F = 0.920 (n_max=3), 0.892 (n_max=4) — equivalently 85% / 80% of V_ee's Frobenius mass-squared is already diagonal — with ‖[H₁,V_ee]‖/‖H₁‖ = 6.1% (n_max=3) dropping only to 5.3% (n_max=4), saturating not vanishing, and the remaining off-diagonal correlation is concentrated on the (1s1s) pair-state at coalescence rather than distributed across channels. Key data: `debug/energy_graph_exploration.md`, `debug/data/energy_graph_nmax{3,4}.json`.

**Graph validity boundary (v2.9.0):** The graph-native CI with kappa=-1/16 violates the variational bound below Z_c ≈ 1.84 (CBS extrapolated ~1.87-1.89). The mechanism: the graph off-diagonal coupling (Z-independent, ~kappa) competes with the diagonal (Z²-scaled). Relative importance scales as 1/(8Z²): 12.5% at Z=1, 3.1% at Z=2, 0.13% at Z=10. Above Z_c, the graph approximation is perturbative and the CI under-binds (variational). Below Z_c, the graph topology overestimates inter-shell coupling and the CI over-binds non-variationally. Standard (non-graph) Casimir FCI is always variational. Variational-k optimization cannot fix the graph over-binding — the problem is the topology, not the exponent. Key data: `debug/data/track_di_z_sweep_variational.json`.

**111 Pauli count derivation (v2.9.0):** Each s/p block at max_n=2 (M=5, Q=10) produces exactly 111 = 55 + 56 non-identity Pauli terms. 55 = Q(Q+1)/2 direct terms (universal, all k=0 Coulomb integrals nonzero). 56 exchange terms from 3 Gaunt channels: s-s (16 Pauli from 12 ERIs), s-p cross-Coulomb (24 from 24 ERIs), k=1 dipole (16 from 16 ERIs). Total 65 ERIs per block from 7²+4²=65 allowed Gaunt quartets. Pure l-shells have zero exchange → Pauli/Q = (Q+1)/2. Universal coefficient 111/10 = 11.1.

**Quantum computing status: ACTIVE FRONTIER.** The composed architecture produces structurally sparse qubit Hamiltonians: O(Q^2.5) Pauli scaling across LiH/BeH2/H2O (exponent spread 0.02), 51x-1,712x advantage over published Gaussian baselines. Full N-electron encoding comparison (Track AS) confirmed composed is categorically sparser (334 vs 3,288 Pauli terms, 20x lower 1-norm). Equal-qubit He comparison validated against cc-pVDZ/cc-pVTZ computed integrals. Commutator-based Trotter bounds give r ~ Q^{1.47}, 7x fewer steps at Q=60. Head-to-head Gaussian comparison (Track AX, v2.0.26): CLEAR WIN on structural sparsity — 190x fewer Pauli terms for LiH at Q~30 (334 vs ~63,500 cc-pVDZ), 746x for H2O at Q=70; accuracy caveat (GeoVac 5.3%-26% R_eq error vs Gaussian <0.1%). H₂ bond-pair qubit encoding (Track AZ, v2.0.27): single-block composed encoding at Z_eff=1, Q^3.13 scaling (consistent with He atomic), 112 Pauli terms at Q=10, R-independent sparsity confirmed. Ecosystem export pipeline (Track AW): OpenFermion, Qiskit, and PennyLane export via `geovac.ecosystem_export`. `geovac-hamiltonians` PyPI package built (Track BB, v2.0.27): standalone 92 KB wheel, 6 bundled modules, all 5 systems (H2/He/LiH/BeH2/H2O). VQE validation (Track AY): H2 converges to 0.031 mHa on statevector simulator; LiH (30 qubits) infeasible for statevector (17.2 GB RAM). H₂O composed 1-norm (Track BD, v2.0.28): total 28,053 Ha at Q=70, but 98.7% from Z²-scaled PK barrier; electronic-only 1-norm 361 Ha (comparable to BeH₂ 355 Ha). PK diagonal on Z_eff=6 1s orbital is 2,387 Ha — Z² scaling is the bottleneck. PK classical partitioning (Track BF, v2.0.29): PK is a one-body operator whose energy E_PK = Tr(h1_pk · γ) can be computed classically from the VQE 1-RDM with zero additional circuits. Operator decomposition H_full = H_elec + H_pk verified to machine precision (<1e-12) for all three composed systems. Algebraic exactness: E_full = E_elec(ψ) + E_PK(ψ) with residual <1e-13 Ha. Partitioned 1-norms: LiH 33.26 Ha (PK 10.9%), BeH₂ ~355 Ha, H₂O 361 Ha (PK was 98.7% of 28,053 Ha total — 78x reduction). `composed_qubit.py` extended with `pk_in_hamiltonian` kwarg; `pk_partitioning.py` new module; `ecosystem_export.py` updated with partitioned API. 19 tests pass. IBM Quantum VQE demo (Track BC, v2.0.28): `demo/ibm_quantum_demo.py` with Aer simulator and IBM hardware modes; H₂ statevector VQE converges ~13 mHa (10 qubits, 80 params). Positioning documents (Track BE, v2.0.28): `docs/geovac_positioning.md` and `docs/geovac_onepager.md` for outreach. General composed builder and first-row generalization (Tracks BG-BJ, v2.0.30): refactored three hardcoded builders (LiH/BeH₂/H₂O) into single general `build_composed_hamiltonian(spec)` driven by `MolecularSpec` dataclass. Atomic classifier (`geovac/atomic_classifier.py`) maps Z→block decomposition for Z=1-10. Ab initio PK parameters computed for Z=3-9 via Level 3 solver (Track BI): Z² scaling has 5-26% errors, computed values required. Three new molecules: HF (Q=60, 667 Pauli), NH₃ (Q=80, 889 Pauli), CH₄ (Q=90, 1000 Pauli). Within-molecule Pauli scaling exponent ~2.2 (2-point fit max_n=1,2), consistent across all 6 composed molecules. Paper 16 promoted to Core (atomic classifier specification). Scope boundary documented (Track BK): PK s-p splitting is structurally wrong-sign (negative result); first row fully supported, second row feasible with frozen-core tabulation, transition metals out of scope. l_max convergence sprint (Tracks BP-BR, v2.0.32): 2D variational solver wired into composed pipeline (level4_method='variational_2d'), faster than adiabatic at l_max=2 (6.7s vs 9.3s/pt). l_max divergence confirmed structural to PK — NOT from adiabatic approximation. R_eq drifts +0.15-0.22 bohr/l_max with BOTH solvers; single-point E(R=3.015) violates variational bound at l_max≥2. l_max=2 is optimal operating point. Partitioned 1-norms confirmed: LiH electronic 33.26 Ha (PK 10.9%), H₂O electronic 361 Ha (PK 98.7% → 78x reduction via classical partitioning). LiH at max_n=3 (Q=84): 7,879 Pauli terms vs Gaussian cc-pVDZ 63,519 at Q=36 — sparsity advantage grows with basis. Accuracy target analysis (docs/accuracy_target_analysis.md): R_eq error is wrong metric for quantum simulation; resource estimation papers operate at fixed geometries and measure Pauli terms, 1-norm, qubit count. Sturmian CI investigation (Tracks BU-1/BU-2, v2.0.33): Coulomb Sturmian CI improves He FCI by 0.92 pp at max_n=3 (2.06% vs 2.98%). Generalized Sturmian shows crossover (better at max_n=2, worse at max_n=3 -- within-config flexibility loss). Gaunt sparsity exactly preserved in both variants. Qubit encoding: 7% fewer Pauli terms (112 vs 120 at Q=10) but 2.8-4.5x higher 1-norm from Lowdin orthogonalization. Negative result for molecular PK bypass. TC scoping (Track BX-1, v2.0.33): Transcorrelated method compatible with GeoVac angular framework. J = -(1/2)r12 preserves Gaunt selection rules. BCH terminates at 2nd order for 2 electrons. Non-Hermitian eigenproblem required. Feasibility: REQUIRES NEW INTEGRALS. TC in second quantization (Track BX-3, v2.0.35; **CORRECTED v2.9.0**): TC-modified qubit Hamiltonians implemented via `geovac/tc_integrals.py`. **Original BX-3 benchmark was incorrect**: qubit-space diagonalization (full 2^Q matrix including wrong-particle-number sectors) produced non-physical energies below the variational bound, giving a false appearance of standard FCI diverging (5.3-8.2%) and TC converging (3.3-3.6%). **Corrected particle-number-projected FCI** (v2.9.0, Track TC-V) shows standard FCI converges monotonically (5.29% → 2.48% → 2.17% at n_max=1-3) while TC plateaus at 3.3-3.4% — standard is MORE accurate than TC at n_max≥2. TC only wins at the trivial n_max=1 basis. Composed molecules: Pauli ratio exactly 1.68× (constant factor from non-Hermiticity), O(Q^2.5) scaling preserved. VarQITE (imaginary-time evolution) converges rapidly to TC ground state (16 steps from HF). BCH constant shift is -1/4 per electron pair. Key files: `geovac/tc_integrals.py`, `debug/tc_verification.py`, `debug/data/tc_verification_results.json`. Quantum resource market test (Track CA, v2.0.36): Head-to-head 1-norm comparison against Gaussian raw JW baselines. LiH STO-3G computed from OpenFermion cached integrals: Q=12, 907 Pauli terms, λ=34.3 Ha, 273 QWC groups. GeoVac LiH composed (electronic-only): Q=30, 334 Pauli, λ=33.3 Ha, 21 QWC groups. Result: 2.7× fewer Pauli terms, 13× fewer QWC groups, 0.97× 1-norm (essentially identical). He equal-qubit at Q=28: GeoVac λ=78.4 vs Gaussian cc-pVTZ λ=530.5 (6.8× lower). Literature survey: published DF/THC/SCDF lambda values for molecules at Q<100 DO NOT EXIST — QPE papers benchmark at FeMoco scale (152+ qubits). GeoVac's positioning is strongest for VQE/NISQ regime where Pauli count and QWC groups dominate runtime. Key files: `docs/market_test_results.md`, `benchmarks/gaussian_baseline_comparison.py`, `debug/data/market_test_data.json`. Balanced coupled composition (Track CD, v2.0.39): cross-center V_ne via multipole expansion of 1/|r-R_B|, terminates exactly at L_max=2*l_max by Gaunt selection rules. LiH at n_max=2: 878 Pauli terms (2.63× composed, Conjecture 1 confirmed), 74.1 Ha 1-norm (1.98×), 87 QWC groups. 4-electron FCI: 1.8% energy error, 7.0% R_eq error, only bound configuration (decoupled 10.9%, PK 15.0%, CB 29.0% — all unbound). D_e=0.037 Ha. Key files: `geovac/shibuya_wulfman.py` (compute_cross_center_vne), `geovac/balanced_coupled.py` (build_balanced_hamiltonian). Phase 4A with analytical V_ne integrals (incomplete gamma, machine precision): n_max=3 E(3.015)=-8.055 Ha (0.20% error). R_eq=3.280 bohr (8.8% error, structural drift +0.053 bohr/n_max, 3x smaller than PK). MIXED: energy converges excellently, R_eq drifts structurally. Optimal for single-point quantum simulation at fixed geometries. Nested hyperspherical investigation (Track DF, v2.5.0, 6 sprints): H-set coupled angular basis on S^(3N-1) produces 18-26% lower ERI density than uncoupled basis via 6j recoupling zeros (novel sparsity mechanism). Compact encoding: Q=10 for 4-electron systems (3× fewer qubits than composed Q=12), comparable Pauli count (112 vs 115), 6.4× lower 1-norm for atoms (PK elimination). Molecular extension: NEGATIVE for all three approaches (single-center R_eq 33.7%, charge-center 48.2% energy error, heterogeneous Löwdin destroys sparsity 1,711 vs 120 Pauli terms). Composed architecture confirmed as structurally necessary for molecular PES. Key files: `geovac/nested_hyperspherical.py`, `tests/test_nested_hyperspherical.py`. Precision atomic spectra (Track DI, v2.6.0, Sprint 1): 2D variational solver on S⁵ eliminates adiabatic approximation. Tensor-product basis: spectral Laguerre (R) × spectral Gegenbauer (α), treating hyperradius and hyperangle simultaneously. Breaks the adiabatic 0.19-0.20% floor: raw 0.022% at l_max=7 (n_basis_R=25, n_basis_alpha=40). With Schwartz cusp correction at l_max=4: 0.004% error (E=-2.90383 Ha vs exact -2.90372 Ha). Phase 1 target (<0.01%) achieved. Adiabatic floor was structural to the Born-Oppenheimer-like separation of R and α; the 2D solver captures the full R-α correlation. l_max convergence is monotonic with ~l^{-2} rate (dominated by per-channel angular basis quality, not partial-wave truncation per se). Key files: `geovac/level3_variational.py` (solve_he_variational_2d, solve_he_precision), `tests/test_level3_variational.py` (17 tests), `debug/data/track_di_he_variational_2d.json`. Algebraic Casimir CI (Track DI, v2.6.0, Sprint 2): Fully algebraic FCI matrix H(k) = Bk + Ck² with exact rational Slater integrals from Paper 7's S³ formula. 21 F^k + 14 G^k integrals verified symbolically. n_max=1: k*=27/16, E*=-729/256 (exact rationals, 1.93%). n_max=3: E*=-2.860 Ha (1.50%, 31 configurations). Self-consistency k²=-2E: NEGATIVE (over-constrains 2-electron problem, 5-13% error; variational k is always better). Three-layer structure: rational (Slater integrals) → algebraic (optimal exponents) → transcendental (e-e cusp). Off-diagonal one-body elements (k-Z)⟨a|1/r|b⟩ required for variational bound. Key files: `geovac/casimir_ci.py` (build_fci_matrix, solve_self_consistent, solve_variational), `tests/test_casimir_ci.py` (32 tests), `debug/data/track_di_casimir_ci.json`. Graph-native CI (Track DI, v2.6.0, Sprint 3C): Graph Laplacian h₁ + rational Slater V_ee, zero free parameters. He ground state 0.19% at n_max=7 (1,218 configs). Beats FCI-A paper (0.22% vs 0.35% at n_max=5) due to exact rational integrals replacing 2000-pt grid numerics. Graph off-diagonal h₁ accounts for 86% of CI correlation energy — the topology is the dominant ingredient. Three-sequence comparison: graph-native (0.19%) >> variational-k (1.41%) >> fixed-k=Z (1.85%). Key files: `geovac/casimir_ci.py` (build_graph_native_fci, _build_graph_h1), `debug/data/track_di_graph_native.json`. Basis invariance verification (Track DI, v2.6.0, Sprint 3D): FCI is invariant under orbital rotation — transforming V_ee to the graph eigenbasis gives identical energies (< 3×10⁻¹⁵ Ha) at n_max=1-4. The ~0.14% convergence floor is NOT from basis mismatch between graph h₁ and hydrogenic V_ee. Floor is from cusp + finite-n_max truncation (embedding exchange constant content, Paper 18). Three-layer structure refined: rational (graph h₁ + Slater V_ee, 98.6% of exact) → topological (T± inter-shell couplings via κ=-1/16, dominant one-body correlation) → transcendental (cusp, ~0.14% floor). Key files: `geovac/casimir_ci.py` (build_graph_consistent_fci), `tests/test_casimir_ci.py` (38 tests), `debug/data/track_di_graph_consistent.json`.

**RH sprint: GeoVac Hopf graphs are strictly Ramanujan (v2.20.0, April 17, 2026):** Ihara zeta computed on the S³ Coulomb graph (Paper 7) and S⁵ Bargmann-Segal graph (Paper 24) via the Ihara-Bass determinantal formula and the 2E × 2E Hashimoto non-backtracking edge operator. All four graphs tested (S³ at max_n=2,3; S⁵ at N_max=2,3) satisfy the Kotani-Sunada Ramanujan bound **strictly**: largest non-trivial Hashimoto eigenvalue sits strictly inside √q_max. Deviations: S³ max_n=3: −0.198; S⁵ N_max=2: −0.208; S⁵ N_max=3: −0.357 (strongly Ramanujan). Closed forms are integer-coefficient polynomial factorizations: S³ at max_n=3 decomposes per ℓ-shell (Fock adjacency preserves ℓ); S⁵ at N_max=3 factors as (s±1)²³·P_12(s²)·P_22(s²), a 12+22 dichotomy conjecturally matching Paper 25's Hopf-U(1) block decomposition (untested). No transcendentals (π, ζ(2), ζ(3)) in any zero — Paper 24's π-free certificate extends from adjacency to Ihara zeta. Standard Stark-Terras single-critical-circle functional equation does NOT hold (graphs are irregular); what does hold is s→−s reflection and Galois closure. Parallel Wick-rotation investigation (Track RH-B, memo at `debug/fock_continuation_memo.md`) is a clean negative: Fock's S³→H³ continuation (Bander-Itzykson 1966) is clean but provides no discrete Γ ⊂ SO(3,1) from framework invariants — Weyl of SO(4) is not a lattice, Hopf S¹ irrelevant for H³, Bianchi PSL(2, O_K) selected by no GeoVac invariant (Rydberg, p₀=Z/n, κ=−1/16, B=42, F=π²/6, Δ=1/40, α). Universal/Coulomb partition sharpened three ways: S³ discrete Hopf calibration-π, S⁵ holomorphic π-free, H³ continuous no-Γ. Graph-zeta is the single RH-adjacent direction intrinsic to GeoVac; Selberg-on-hydrogen is structurally closed. Paper 29 drafted (`papers/observations/paper_29_ramanujan_hopf.tex`). Literature survey (`debug/rh_literature_survey.md`, 60 entries) confirmed RH-A is not duplicated in the literature to Jan 2026; spectrum-as-zeros is ruled out for our integer Dirac spectrum (Berry-Tabor + Odlyzko) — any RH bridge must go through graph zeta or derived operators. Key files: `geovac/ihara_zeta.py`, `tests/test_ihara_zeta.py` (15/15 passing), `debug/ihara_zeta_memo.md`, `debug/data/ihara_zeta_geovac_hopf.json`. Open questions: N_max=5 extension (Paper 25 headline, 330-dim Hashimoto eigensolve), 12+22 ↔ Hopf-U(1) hypothesis test, Alon-Boppana asymptotics of the sub-Ramanujan deviation, combinatorial limit to classical Riemann zeta.

**RH-C follow-up (same sprint, v2.20.0):** Spin-ful Ihara zeta extension onto the Dirac-S³ graph (nodes labeled by DiracLabel (n_fock, κ, m_j) from `geovac/dirac_matrix_elements.py`). Two adjacency rules computed: Rule A (κ-preserving, direct lift of scalar Fock ladders) and Rule B (E1 dipole, parity-flip Δl=±1, mixes κ). Both rules Ramanujan at n_max=1,2,3. Rule A decomposes per-κ analogously to the scalar per-ℓ-shell decomposition, with tighter deviation (−0.061 at n_max=3) than the scalar S³ Coulomb at the same scale (−0.198). **Rule B at n_max=2 saturates the Ramanujan bound exactly (deviation = 0.000)**, with a (4s²+1)² factor placing zeros precisely on the critical circle |s|=1/√q_max; at n_max=3 Rule B moves strictly inside (deviation −0.119) and its Ihara zeta factorizes as (s±1)⁷⁹·(9s²+1)⁴·P_22(s²)·P_24(s²), a 22+24 dichotomy paralleling the scalar S⁵ 12+22 dichotomy. All zeros algebraic, no transcendentals, π-free certificate preserved. Key files: `geovac/ihara_zeta_dirac.py` (reuses RH-A machinery, zero-edge-graph and charpoly-variable bugs patched in-session), `tests/test_ihara_zeta_dirac.py` (35/35 passing, 2 skipped by design), `debug/ihara_zeta_dirac_memo.md`, `debug/data/ihara_zeta_dirac_s3.json`, `debug/compute_ihara_zeta_dirac.py` (reproducible driver).

**Sprint 2 RH follow-up (v2.21.0, April 17, 2026):** Four parallel tracks on the Paper 29 open questions. **(RH-D Alon-Boppana sweep):** extended Ihara zeta computation to max_n=4,5 (S³ Coulomb) and N_max=4,5 (S⁵ Bargmann-Segal, including the Paper 25 headline V=56 E=165 β₁=110 case) plus Dirac-S³ Rule A+B at n_max=4. **Three of four graph families CROSS the Kotani-Sunada Ramanujan bound at V≈30-60**, with signed deviation growing linearly in V: S³ Coulomb at max_n=5 (V=55, deviation +0.20, R²=0.994), S⁵ Bargmann-Segal at N_max=5 (V=56, deviation +0.65, R²=0.766), Dirac-S³ Rule B at n_max=4 (V=60, deviation +1.53, R²=0.828). Only Dirac-S³ Rule A remains sub-Ramanujan (|deviation| ~ V^−3.05, R²=0.994). Conjectured mechanism: mixed-density topology — q_max set by single densest vertex while ρ(T) driven by graph volume; dense near-regular sub-blocks appear with Perron > √q_max as V grows. **Graph-RH for GeoVac is a finite-size statement (secured for N_max ≤ 3), not asymptotic**. **(RH-E Hopf-U(1) block hypothesis):** **Paper 29 Hypothesis 1 VALIDATED**. The 12+22 (scalar S⁵ N_max=3) and 22+24 (Dirac-S³ Rule B n_max=3) factorizations are exactly the Z₂ m-reflection (P: m→−m) block decompositions of the node-level Ihara-Bass matrix M(s) = I − sA + s²Q. The Z₂ is the only sub-action of the Paper 25 Hopf U(1) that commutes with a real integer adjacency; the continuous U(1) phase rotation reduces to this Z₂ on real data. Verified symbolically: det(M_sym) = (s−1)(s+1)·P_22 on the 13-dim symmetric block (scalar S⁵); det(M_antisym) = P_12 on the 7-dim antisymmetric block; product equals full det(M). Same mechanism for Dirac-S³ Rule B with equal-dim 14+14 blocks under m_j→−m_j. Hypothesis → Observation in Paper 29 §5.3. **(RH-G combinatorial limit to Riemann):** STRUCTURAL REDIRECTION. The GeoVac Hopf graphs are bipartite, so tr(T^L)=0 for odd L — Ihara walk-length spectrum supported only on even integers. Paper 28's χ_−4 Dirichlet character has support only on odd integers, so χ_−4-twisted Ihara zeta is identically trivial. **Paper 28's χ_−4 structure lives on the SPECTRAL zeta D(s) = Σ g_n·|λ_n|^{−s} (Dirac Dirichlet series at quarter-integer Hurwitz shifts), NOT on the Ihara side**. Any bridge to classical Riemann ζ(s) must go through the spectral zeta, not the Ihara zeta. Confirmed by Hashimoto spectrum Poisson statistics (CV ≈ 1, Berry-Tabor, not GUE=0.42) — direct Hilbert-Pólya identification categorically excluded. **(RH-I literature re-survey with verified web):** Web access worked. 16 verified post-2024 entries appended to `debug/rh_literature_survey.md` §6. **RH-A is confirmed NOT duplicated** in the literature to April 2026. No Paper 29 bibliography errors; three citations added (Matsuura-Ohta 2025 PTEP — Kazakov-Migdal/Ihara on graphs; Yakaboylu 2024 JPhysA — new HP candidate on half-line, not compact; Huang-McKenzie-Yau 2024 — Ramanujan context). Paper 29 synced to v1.1 with all four results integrated (abstract, §5.3 Hypothesis → Observation, §6.1 Alon-Boppana resolved, §6.3 Hopf-block validated, §6.4 continuum-limit redirected, Conclusion updated). Sprint 3 recommendation: **test the spectral-zeta χ_−4 angle** — run higher-s PSLQ on D_even(s), D_odd(s) from `geovac/qed_vertex.py` for a functional-equation identification of L(s, χ_−4) as a spectral-zeta transform. This is the actual Prize-bound direction and does not require new graph code. Key Sprint 2 files: `debug/data/alon_boppana_sweep.json`, `debug/alon_boppana_memo.md`, `tests/test_alon_boppana_sweep.py`; `debug/hopf_u1_block_test.py`, `tests/test_hopf_u1_block.py` (12/12 passing), `debug/hopf_u1_block_test_memo.md`, `debug/data/hopf_u1_block_test.json`; `debug/compute_riemann_limit.py`, `debug/compute_riemann_twist_extended.py`, `debug/riemann_limit_memo.md` (3330 words), `debug/data/riemann_limit_data.json`; `debug/rh_literature_survey_sprint2_findings_memo.md`.

**Sprint 3 RH follow-up (v2.22.0, April 17, 2026):** Four parallel tracks on the Sprint 2 redirection to the spectral-zeta side. **(RH-J Spectral χ_−4 identification, POSITIVE):** exact closed-form identity **D_even(s) − D_odd(s) = 2^(s−1)·(β(s) − β(s−2))** for all integer s ≥ 2, where β(s) = L(s, χ_−4). Two Hurwitz identities cause the ζ(s, 1/4) pieces to cancel. Verified symbolically (sympy residual = 0) at s ∈ {2..10, 15} and by PSLQ at 100-digit precision. Quantitatively upgrades Paper 28's qualitative "vertex parity acts as χ_−4" to a concrete spectral-zeta-to-Dirichlet-L bridge. Paper 28 §Vertex-parity extended with new subsection "General-s closed form for D_even − D_odd" (eq:D_diff_closed). Does NOT by itself constrain Re(ρ) = 1/2: β zeros at ρ propagate to D_diff zeros at ρ and ρ+2, consistent with but not enforcing a critical-line structure. **(RH-M Spectral zero statistics, HEADLINE POSITIVE):** computed numerical zeros of D(s), D_even(s), D_odd(s) in strip Im(s) ∈ [0, 70] at 35-digit precision (n = 13, 21, 15 zeros respectively). Unfolded imaginary-part spacing CVs: 0.348, 0.343, 0.396 — all **GUE-like** (synthetic GUE CV ≈ 0.46; Poisson CV ≈ 1.05). KS rejects Poisson for D_even (p=0.004) and D_odd (p=0.030); cannot reject GUE. **First GUE signature anywhere in GeoVac.** The two sides of the Weil dictionary now give opposite verdicts: Ihara side Poisson (Sprint 2 RH-G), spectral side GUE (RH-M). Honest caveat: zeros are NOT on a single critical line — scattered in Re(s) ∈ [2.02, 3.42]. GUE-spacing yes, RH-style line no. Classical ζ has both RMT spacing AND critical-line constraint; GeoVac has only the first. Sample size n=13-21 gives CV SE ~0.15; definitive identification requires Im(s) ≤ 300 (~6 hours compute). **(RH-F Arithmetic/Galois structure):** factored Paper 29 integer polynomials over ℚ, ℚ(i), ℚ(√2), ℚ(√5), ℚ(ω), cyclotomics. Three tiers: (a) cyclotomic Φ_3/Φ_4/Φ_6 factors have ℤ/2 Galois, fingerprint specific graph cycles; (b) Ramanujan-boundary factors (4s²+1)², (9s²+1)⁴ live in ℚ(i) — consistent with the Hopf-ℤ_2 structure; (c) bulk factors have non-abelian Galois groups S_3 (cubics), S_4 (Dirac-A quartic), **S_6 (non-solvable) for P_12(s²) in scalar S⁵ N_max=3** — Abel-Ruffini: its 12 zeros have no radical expression over ℚ. By Kronecker-Weber, non-abelian Galois ⇒ not inside any cyclotomic. Paper 25's Hopf-U(1) controls the ℤ_2 block-split but does NOT abelianize the Galois group inside each block. Paper 29 Corollary 1 (π-freeness) strengthened to Corollary 2 (integer-algebraicity) plus Observation 3 (Galois structure). Paper 18 taxonomy gains new tier: "algebraic but non-radical" (P_12 zeros are the first explicit example). **(RH-K α²-weighted Ihara zeta on Dirac-S³, NEW OBJECT):** first non-0/1 edge weighting in the Ihara framework. Weight convention w(a,b) = 1 + α²·|f_SO(a) + f_SO(b)|/2 using the dimensionless Breit-Pauli SO diagonal f_SO(n,κ) = −(κ+1)/[4n³l(l+½)(l+1)]. Module `geovac/ihara_zeta_weighted.py` (14/14 tests passing). Lives in ring ℚ(α²)[s] — realizes Paper 18's spinor-intrinsic tier without γ. Four findings: (i) α² enters as a graph-invariant, not just edge-local (Rule A n_max=3: α-polynomial up to degree 44, α² coefficient has clean prefactor −35/648); (ii) three of four Ramanujan verdicts unchanged under weighting, only the boundary case Rule B n_max=2 (deviation=0 unweighted) crosses for any α > 0 — generic perturbation-of-boundary, not α²-specific; (iii) weighting pushes Hashimoto CV AWAY from GUE (Rule B n_max=3: 2.6 → 5.5), closing the Ihara-side Hilbert-Pólya route further; (iv) RH-K therefore confirms that the GUE-like behavior in RH-M is SPECTRAL, not latent in a weighted Ihara. Paper 29 v1.2 synced with both §5.4 (Corollaries + Observation) and §6 (RH-M GUE paragraph + RH-K α²-weighted paragraph). Sprint 4 scope below (Hilbert-Pólya operator reconstruction from RH-M zero list, Bianchi-like non-abelian Wilson gauge on S³=SU(2), depth-2 analog of RH-J.1 for S_min). Key Sprint 3 files: `debug/compute_spectral_chi_neg4.py`, `debug/spectral_chi_neg4_memo.md`, `debug/data/spectral_chi_neg4.json`, `tests/test_spectral_chi_neg4.py` (25 tests); `debug/compute_spectral_zero_stats.py`, `debug/spectral_zero_stats_memo.md`, `debug/data/spectral_zero_stats.json`, `tests/test_spectral_zero_stats.py` (13 tests); `debug/compute_galois_ihara.py`, `debug/galois_ihara_memo.md`, `debug/data/galois_ihara.json`, `tests/test_galois_ihara.py` (27 tests); `geovac/ihara_zeta_weighted.py`, `debug/ihara_zeta_weighted_memo.md`, `debug/data/ihara_zeta_weighted_dirac_s3.json`, `tests/test_ihara_zeta_weighted.py` (14 tests).

**Sprint 4 RH follow-up (v2.23.0, April 17, 2026):** Four parallel tracks on the Hilbert-Polya reconstruction, functional equation hunt, depth-2 closed form, and non-abelian gauge extension. **(RH-N HP operator reconstruction, structural):** four Hermitian constructions (diagonal, tridiagonal Jacobi, Toeplitz, companion) built to match the Sprint 3 RH-M zero data. Cross-construction Frobenius distances 0.40–0.81 quantify "eigenvalues don't determine H*". Headline: γ_n ≈ 18.5·√n Weyl law across D(s), D_even, D_odd (RMS residual 0.8–0.9, beats linear 2–4× and n·log n 3–5×). This is 1D-HO-class Weyl density, categorically different from Dirac-on-S³'s linear n+3/2 and from classical Riemann's log-density. **H\* is not a deformation of any native GeoVac operator.** New module `geovac/hp_operator.py` (13 tests). **(RH-O Functional equation hunt, DECISIVE CLEAN NEGATIVE):** 13,080-template grid search (6,480 single-Gamma + 6,600 Gamma-pair) across 15 candidate reflection axes at 10 complex test points, 50 dps. Classical ξ_ζ and ξ_β sanity anchors close to 10⁻⁵¹; RH-J identity closes to 10⁻⁵⁰; best GeoVac residual is 0.594 — **48 orders of magnitude worse**. Structural obstruction: each D(s), D_even(s), D_odd(s) is a linear combination of TWO Hurwitz zetas at the same shift but different exponents (s and s−2) with mismatched weights; no single Γ-product can simultaneously complete both. For D_diff, the two β-pieces demand different reflection axes (c=1 for β(s), c=5 for β(s−2)). **The missing critical-line-forcing ingredient is structurally absent, not undiscovered.** Combined with RH-M (zeros not on single line) and RH-N (wrong Weyl class), Sprint 4 **closes the direct spectral-to-classical-RH bridge** for GeoVac. **(RH-P Depth-2 χ_−4 analog, negative + side flag):** S_min split by parity of k gives S_min = S_min^even + S_min^odd; 35 PSLQ attempts across 7 basis strategies at 100 dps against {β(s), β(s−2), G, ζ(2k), D_basis, ultra-wide, depth-2 analog} — zero identifications. Structural obstruction: T(k)² is a quadratic object in Hurwitz values; the Bernoulli-asymptotic coefficients scramble the uniform prefactor structure that made RH-J.1's telescoping work at depth 1. S_min^even, S_min^odd live in weight-2 MZV sector not spanned by Dirichlet-L products. **IMPORTANT SIDE FINDING flagged for independent verification** (Sprint 5 diagnostic, NOT yet applied to Paper 28): RH-P's 100-dps Bernoulli-asymptotic computation gives S_min = 2.47993693803422255..., while Paper 28 publishes S_min = 2.47953699802733387... (discrepancy ~4×10⁻⁴ at 4th decimal). The RH-P agent claims `debug/smin_identification.py` has an incorrect tail correction. Independent verification agent dispatched (Euler-Maclaurin + direct summation to N=50000). The irreducibility claim (15 PSLQ failures) is independent of the numerical value and stands regardless. PAPER 28 NOT MODIFIED pending verification. **(RH-Q SU(2) Wilson lattice gauge on S³=SU(2), CLEAN POSITIVE):** genuine non-abelian Wilson lattice gauge theory built on the Fock-projected S³ Coulomb Hopf graph, extending Paper 25's U(1) structure. SU(2)-valued link variables U_e on directed edges; plaquettes from primitive NB walks; S_W = β·Σ(1 − ½·Re Tr U_plaq); character expansion via SU(2) Haar. Three structural results: (i) Paper 25's U(1) action is EXACTLY the maximal-torus limit of this SU(2) theory (diagonal config matches U(1) to machine precision), upgrading Paper 25 to a Cartan-subalgebra projection of a non-abelian theory; (ii) weak-coupling linearization gives exactly Paper 25's discrete Hodge-1 Laplacian L_1 = B^T B as the kinetic term — **Paper 25's L_1 is the SU(2) kinetic sector**; (iii) plaquette counts at n_max=2,3,4: (0,0,0), (2,1,1), (8,7,18) at length L=4,6,8; Monte Carlo Wilson loops on n_max=3 at β∈{0.5,1,2,5} give ⟨W⟩=0.126/0.158/0.332/0.681 (monotonic, no area/perimeter crossover at this size — only one plaquette class). HONEST FRAMING: this is a gauge-invariant Haar-normalizable non-abelian lattice gauge theory on a compact finite graph; it is NOT a continuum Yang-Mills mass-gap proof on ℝ⁴ (the Millennium Problem). New module `geovac/su2_wilson_gauge.py` (26 tests). **Potential Paper 30 material** (Wilson-Hodge SU(2) non-abelian sibling of Paper 25). Paper 29 v1.3 synced: §6 extended with four Sprint-4 paragraphs (RH-N, RH-O, RH-P, RH-Q). Paper 28 NOT MODIFIED pending S_min verification. Sprint 5 scope: (a) S_min verification diagnostic resolves Paper 28 flag; (b) pursue potential Paper 30 draft (SU(2) Wilson gauge extension); (c) χ_−4-character-twisted depth-2 sums (RH-P agent alternative suggestion); (d) larger-scale RH-M sample (n ≈ 60-100, Im s ≤ 300) for GUE identification at CV SE ~0.05. Key Sprint 4 files: `geovac/hp_operator.py`, `tests/test_hp_operator.py` (13 tests), `debug/hp_operator_memo.md`, `debug/data/hp_operator.json`; `debug/compute_functional_equation_hunt.py`, `debug/functional_equation_hunt_memo.md`, `debug/data/functional_equation_hunt.json`; `debug/compute_smin_chi_neg4.py`, `debug/smin_chi_neg4_memo.md`, `debug/data/smin_chi_neg4.json`, `tests/test_smin_chi_neg4.py` (10 tests); `geovac/su2_wilson_gauge.py`, `debug/su2_wilson_gauge_analysis.py`, `debug/su2_wilson_gauge_memo.md`, `debug/data/su2_wilson_gauge_results.json`, `tests/test_su2_wilson_gauge.py` (26 tests).

**Sprint 5 S_min erratum patch (v2.23.1, April 17, 2026):** the Sprint 4 RH-P side flag was independently verified by three methods (direct sum + 5-term explicit tail; mpmath.nsum Levin u-transform; Euler-Maclaurin 26-term asymptotic) all agreeing at ≥ 80 digits on the true value **S_min = 2.47993693803422255441357950082938214468792578661728845837879872655955...** Double-diagnosis: (a) Paper 28's published 2.47953699802733387 was wrong at the 4th decimal due to an erroneous tail formula in `debug/smin_identification.py` lines 45-48 which assumed T(k)~2/a² (the correct asymptotic is T(k)~2/a); (b) RH-P's 2.47993693803422255447852790477 was wrong at the 20th digit due to incorrect c_6..c_11 coefficients in `debug/compute_smin_chi_neg4.py` T_SQUARED_COEFFS. **Fixes applied**: Paper 28 §IV text updated with the corrected value at 25 digits; `debug/smin_identification.py` compute_S_min_with_tail now delegates to mpmath.nsum Levin (2-line fix, ~2 s runtime, verified to 15+ digit match); `debug/compute_smin_chi_neg4.py` T_SQUARED_COEFFS retained with documented 20-th-digit limitation (does not affect the RH-P PSLQ negative result, which was Q-linear-independence against a finite basis insensitive to 6.5e-20 numerical shifts). Irreducibility of S_min against the 47-element extended basis **unchanged** — 15 PSLQ failures stand. Key files: `debug/compute_smin_verification.py`, `debug/smin_verification_memo.md`, `debug/data/smin_verification.json`.

**Sprint 7 Paper 30 draft (Wilson-Hodge SU(2), v2.25.0, April 17, 2026):** Per Leader's Strategic Brief Direction 2. `papers/observations/paper_30_su2_wilson.tex` drafted (1014 lines, ~4000 words, 9 sections), formalizing Sprint 4 RH-Q's SU(2) Wilson lattice gauge construction as the non-abelian sibling of Paper 25. Two propositions: **(1) Maximal-torus reduction** — SU(2) links restricted to the Cartan torus yield Paper 25's U(1) Wilson action *exactly* (verified machine-precision at `tests/test_su2_wilson_gauge.py` lines 142-160); Paper 25 is re-read as the Cartan-subalgebra projection of a non-abelian theory. **(2) L_1 = B^T B is the weak-coupling kinetic term** — second-order expansion of S_W around the trivial vacuum gives (β/8)⟨A, L_1 A⟩ for each su(2) component; non-abelian self-couplings enter only at quartic order. Paper 25's "gauge propagator" identification is the Gaussian weak-coupling limit. **Data**: plaquette counts (0, 2, 8) at length L=(4,6,8) for n_max=(2,3,4); Monte Carlo Wilson loops on n_max=3 at β ∈ {0.5, 1, 2, 5} give ⟨W⟩ = 0.126, 0.158, 0.332, 0.681 (monotonic, matching leading-order character expansion at ~20%). **New §6: framework-wide least-action synthesis** — five action principles (Paper 0 entropy-max on 0-cochains; Rayleigh-Ritz on graph Laplacian Hamiltonian; Wilson action on gauge sector; Ihara zeta path integral from Paper 29; spectral action from Paper 28) each organize a different sector. Paper 0's max-entropy on nodes and Paper 30's min-action on gauge sector are structural duals on 0-cochains vs 1-cochains; joint variational principle flagged as future work. **Scope-honest**: Paper 30 IS a gauge-invariant, Haar-normalizable, non-abelian lattice gauge theory on a compact finite graph — the non-abelian extension of Paper 25 with three structural results. It is NOT a continuum Yang-Mills mass-gap result on ℝ⁴; NOT Lorentzian; NOT a confinement demonstration at this plaquette-class size. Conclusion names the four-way S³ coincidence (Fock projection image / Hopf graph base / Dirac spin carrier / SU(2) gauge manifold) with explicit CLAUDE.md §1.5 disclaimer against ontological claims. 26/26 tests passing. Agent-flagged review items: §6 is most interpretive (framework-wide synthesis); §5.2 L_1 scalar constant not nailed down symbolically (positive multiple via plaquette-length distribution); §8 SU(2) K-formula analog of Paper 2 framed pessimistically per rhetoric rule. Key file: `papers/observations/paper_30_su2_wilson.tex`.

**Sprint 6 Paper 18 consolidation (light-touch, v2.24.0, April 17, 2026):** Per Leader's Strategic Brief (Direction 1), Paper 18 v1.0 was behind the state of the project after Sprints 1-5 (Paper 29 Ramanujan, RH-J β(s)/β(s−2) identity, RH-M GUE signature, RH-F S_6 non-solvable Galois, RH-K α²-weighted Ihara, Phase 4G no-common-generator). Light-touch consolidation scope (v1.0 → v1.1): abstract extended with the four major findings + meta-pattern cross-reference; new §II.E "Level 5: RH-sprint exchange constants" with four sub-paragraph realizations (RH-J closed form, GUE-vs-Poisson tier, algebraic-but-non-radical tier from S_6 non-solvability, α²-weighted Ihara in ℚ(α²)[s]) plus taxonomic placement; new §VII "Structural incommensurability within GeoVac observables: a meta-pattern" naming the three independent instances (Paper 28 QED three classes, Phase 4G K=π(B+F−Δ) decomposition, RH-sprint Weil dictionary) as a meta-observation. Bibliography extended with loutey_paper24/28/29. Three-taxonomy decomposition audit verdict: **PARTIAL unification** (T1 operator-order-axis ↔ T2 sector-axis clean; T3's combinatorial-vs-spectral axis is distinct). Meta-pattern is the correct headline, not a unification theorem. Reviewer memo recommends deferring deep restructure (§IV three-axis grid rewrite, §V Phase 4G theorem treatment, new §VIII Meta-observations section) to Sprint 7-8. **Sub-agent cap:** Reviewer and Decomposer agents hit plan usage cap before producing their memos; PM executed both diagnostics from context (`debug/paper_18_review_memo.md`, `debug/three_taxonomy_decomposition_memo.md`). Key files: `papers/core/paper_18_exchange_constants.tex` (v1.0 → v1.1), `debug/strategic_brief_2026_04_17.md`, `debug/paper_18_review_memo.md`, `debug/three_taxonomy_decomposition_memo.md`.

**Scope boundary:** See `SCOPE_BOUNDARY.md` at project root for which atoms and molecules are supported, which are feasible, and which are out of scope. First row (Z=1-10) fully supported; second-row atoms (Z=11-18) supported via [Ne] frozen cores; third-row s-block (Z=19-20) via [Ar] frozen cores; third-row p-block (Z=31-36) via [Ar]3d¹⁰ frozen cores; fifth-row s-block (Z=37-38) via [Kr] frozen cores and sixth-row s-block (Z=55-56) via [Xe] frozen cores (Sprint 3 HA-A+B, v2.12.0). Full first transition series (Z=21-30) implemented as hydrides (v2.8.0): all 10 TM hydrides (ScH through ZnH) built with d-orbital blocks (l_min=2), Q=30, 277 Pauli terms, Pauli/Q=9.23 (below main-group 11.11). Atomic classifier extended to Z=1-30 with structure type F for d-block. Cr (3d⁵4s¹) and Cu (3d¹⁰4s¹) anomalous configurations handled. General `build_composed_hamiltonian(spec)` implemented in composed_qubit.py. Multi-center molecules (Track CU, v2.3.0): 8 multi-center systems. Heavy-atom alkaline-earth monohydrides SrH (Sprint 3 HA-C, v2.12.0, [Kr] core) and BaH ([Xe] core) added — scalar Q=20, 222 Pauli (isostructural with CaH/KH); relativistic Q=20, 534 Pauli (isostructural with CaH_rel, bit-identical λ_ni=13.87 Ha and QWC=52 across CaH_rel/SrH_rel/BaH_rel). RaH deferred pending [Rn] frozen core. Total library: 40 molecules (20 single-center + 8 multi-center + 10 transition metal hydrides + 2 heavy-atom monohydrides).

**Classical solver details — active frontier items (resolved or characterized):**
- l-dependent PK pseudopotential for LiH composed geometry: reduces R_eq error from 6.9% to 5.3% at l_max=2 (higher l_max with channel-blind PK diverges — see Section 3). **l_max divergence root cause confirmed (Track BQ, v2.0.32): NOT from adiabatic approximation.** Variational 2D solver (level4_method='variational_2d') produces identical drift: R_eq = 3.18, 3.33, 3.63 bohr at l_max=2,3,4 (+0.15-0.22 bohr/l_max). Single-point E(R=3.015) also diverges: -8.184, -8.204, -8.236 Ha (exact: -8.071), violating variational bound — PK overcounting worsens with angular channels. CBS extrapolation: R_eq→3.73 bohr (24% error). Root cause is structural to PK/composed architecture: higher-l channels interact with the l-dependent PK barrier in ways that systematically overcount correlation. **l_max=2 is the optimal operating point** — higher l_max degrades both R_eq and single-point energy. Diagnostic arc (v2.0.5) confirmed: divergence is linear (+0.23 bohr/l_max), angular spreading is 100% nuclear. R-dependent PK w_PK(R) = δ_{l,0} × min(cap, R/R_ref) achieves 2.0% at l_max=4 but requires empirical R_ref. Ab initio R_ref derivation attempted (Track AD, v2.0.18): all core-derived candidates (core radius, PK width, screening length: 0.27-0.82 bohr) too small to affect PES — R_ref needs to be molecular-scale (~3 bohr) but no atomic property bridges this gap without circular reasoning. Negative result.
- Polyatomic PES: BeH2 full 1-RDM exchange achieves 11.7% R_eq error (down from 20% diagonal S·F⁰, matches 12% fitted model with zero free parameters). Kinetic orthogonalization tested: +0.03 Ha uniform, no R_eq effect (negative result). Remaining 11.7% attributed to basis truncation (l_max=2) and PK structural overcounting (NOT adiabatic approximation — see Track BQ).
- H₂O composed solver: 5-block architecture (O 1s² core + 2 O–H bond pairs + 2 lone pairs), R_eq 26% error with charge-center origin. Coupling architecture validated (bond-bond ~0.5 Ha, consistent with BeH₂), but lone pair coupling unphysical at Z_eff=6 (S·F⁰ produces −28 Ha bond-lone, −15 Ha lone-lone — exceeds total electronic energy). Bottleneck is Level 4 angular basis at 6:1 charge asymmetry (R_eq error scales: 0.1% at 1:1, 18% at 2:1, 17% at 6:1 with charge-center).
- Algebraic audit confirmed: Level 4 nuclear coupling (Paper 15) uses split-region Legendre expansion with Gaunt integrals — not quadrature. Core screening (Paper 17) has algebraic density from channel coefficients. Paper 15 Section V.A corrected to match production code.
- Chemical periodicity as representation theory (Paper 16) -- computational exploitation of hierarchical structure

**Completed investigation tracks (v2.0.6-23) — see CHANGELOG.md for full details:**

The classical solver investigation comprised 30+ tracks across four phases. Key structural findings by phase:

*Spectral solvers and coupled-channel infrastructure (v2.0.6-11, Tracks A-K):* Replaced FD solvers with spectral bases at Levels 2-4, achieving 100-269× dimension reductions. Algebraic Laguerre matrix elements at Levels 2-3 (zero quadrature for m=0; single transcendental seed e^a·E₁(a) for m≠0). Level 3 coupled-channel convergence ceiling 0.19-0.20% (adiabatic truncation). Perturbation series definitive negative (branch point obstruction). 2D solver integrated into composed pipeline.

*Algebraicization (v2.0.12-13, Tracks N-S):* Level 3 eigenvalues μ(R) proven algebraic over Q(π,√2), not transcendental. Level 4 angular Hamiltonian structurally non-algebraic (piecewise-smooth). Algebraic Z_eff and Slater integrals wired as production defaults.

*Cusp attack (v2.0.14-16, Tracks U-Y):* Three negative results (alpha-only cusp factor, graph absorption of 1/r₁₂, θ₁₂-adapted basis) — cusp is 2D in (α, θ₁₂), not separable. One positive: Schwartz partial-wave extrapolation breaks through 0.19-0.20% floor (He: 0.10% at l_max=2).

*Channel convergence and diagnostics (v2.0.17-18, Tracks Z-AD):* Level 4 convergence through l_max=6: 96.0% D_e, CBS ~97%. Even-odd staircase is selection rule physics (gerade constraint). Geometric elevation structurally irreducible (three avenues negative). PK essential for equilibrium; R_ref derivation failed (scale mismatch).

*Full N-electron architecture (v2.0.19-23, Tracks AF-AS):* Built full 4-electron solver (Level 4N) as exact validation. Equilibrium exists without PK at l_max≥2. All radial solvers exhausted. Angular basis is bottleneck: l_max=2 S₄ [2,2] has 12 channels vs composed geometry's 144× compression at 5.3% R_eq. Composed quantum encoding categorically sparser (334 vs 3,288 Pauli terms, 20× lower 1-norm).

**Backlog (triaged):**
- Rebuild composed-geometry Hamiltonians with spectral Level 2 solver — DEFERRED TO QUANTUM PHASE. Low priority: spectral Level 2 would improve H2+ accuracy but composed pipeline already uses Level 3+4.
- l_max=4-5 for full 4-electron solver — SHELVED. Requires 6,250+ spectral dim (Track AK); all radial solvers exhausted; angular basis is the bottleneck, not radial treatment.
- Dense spectral PES scan (200+ R-points) for precision H₂⁺ spectroscopic constants — DEFERRED TO QUANTUM PHASE. Classical accuracy is characterized; value is in quantum Hamiltonian validation.
- **α combination rule derivation (Paper 2 → Paper 18)** — PAUSED (April 2026): seven-sprint structural decomposition complete (Phases 4B-4G), all three K ingredients (B, F, Δ) structurally identified, combination rule remains conjectural with the question reframed from 'derive each piece' to 'explain why the sum equals α⁻¹.' Detailed sprint history below. — Phase 4B sprint complete (April 2026, Tracks α-A/α-B/α-C/α-D), 3 of 4 negative or partial. K = π(B + F − Δ) remains structurally underived. (α-A) Hopf-twist spectral comparison S³ vs S¹ × S²: clean negative — neither ζ functions, regularized determinants, truncated Casimir traces, nor heat-kernel coefficients of the two manifolds reproduce K, K/π, B+F, or B+F−Δ at better than 10⁻³; cleanest near-miss ζ_{S¹×S²}(4)/ζ_{S³}(4) = 43.943 vs B+F = 43.645 (rel err 6.84×10⁻³). One more candidate mechanism eliminated. (α-B) Packing-π hypothesis: π is class-matched but NOT rigorously forced. Paper 0's σ₀ = πd₀²/2 does pin a pi^1 of ball-volume type, but the Weyl exchange constant for S² is 1/(4π), not π — Paper 18's wording "outer π = S² Weyl exchange constant" is numerically off by 4π² and should be reframed as "outer π = ω₂, the 2D ball volume". F = ζ(2) is a genuinely independent transcendental injection; any K-derivation must account for at least two independent π-type quantities. Minor observation Δ = 1/(B−2) = 1/40 (single-point coincidence, low confidence). **Paper 18 wording fix flagged for plan-mode review** — not auto-applied. (α-C) Avery-Aquilanti / Sturmian: POSITIVE PARTIAL. κ = −1/16 is structurally present in the Fock momentum-space weight ⟨n,l|(p²+p₀²)⁻²|n,l⟩ at p₀=1 (every denominator divides 32 = 2|κ|⁻¹). New exact rational identity: Σ_{(n,l), n≤3} (2l+1)·l(l+1)·⟨n,l|w|n,l⟩ = 6·B·|κ| = 252/16 = 63/4 (averaged over the 6 (n,l) cells: B·|κ| = 21/8). First exact algebraic statement linking the calibration κ and the Hopf base invariant B = 42 in the Fock formalism. F = π²/6 and Δ = 1/40 not derived. (α-D) Hopf graph morphism: clean negative. Built explicit S³ graph at n_max=3 (14 nodes, 13 edges, golden-ratio Fibonacci spectrum) and S² quotient (6 sectors, eigenvalues {0,0,0,1,3,6}); no Laplacian-spectral invariant — trace, log det, ζ, Cheeger, von Neumann entropy — hits any K target without reinserting (2l+1)l(l+1) by hand. B=42 only appears via the Casimir sum on sector labels, which is Paper 2's existing construction reframed in graph-quotient language. **Net status:** four mechanisms eliminated, one structural κ↔B link established (α-C), Paper 18 wording fix flagged for plan-mode review, Paper 2 unchanged. Data: `debug/data/track_alpha_phase4b/`. Phase 4C sprint complete (April 2026, Tracks α-E/α-F/α-G), both substantive tracks negative. (α-E) S¹ fiber Fock weight: NEGATIVE — F = ζ(2) is NOT in the discrete Hopf fiber at n_max=3. Standard Fock weight ⟨n,l,m|w|n,l,m'⟩ is m-trivial (rotation-invariant kernel collapses to δ_{m,m'}). Path-graph fiber spectral zeta sums are pure rationals under all four weighting schemes (uniform 20/3, degeneracy 28, Casimir 88/3, Hopf 136); none contain π. The continuum limit of the rescaled path-graph zeta IS π²/6 (standard identity), but the finite-size correction at n=5 (Paper 2 cutoff l_max=2) is −π²/150 ≈ −0.0658, NOT a rational multiple of Δ = 1/40 (ratio is −4π²/15, transcendental). n_max invariance fails: fiber sums grow under cutoff extension. Recommendation accepted: F is an embedding exchange constant from continuum S¹ regularization, not a graph invariant — the K formula mixes three distinct algebraic tiers (B rational base Casimir, F transcendental continuum ζ, Δ rational finite-size). (α-F) Δ identity test: NEGATIVE — Δ = 1/(B−N_init) = 1/(42−2) = 1/40 is a SINGLE-POINT COINCIDENCE at m=3, NOT a polynomial identity. LHS(m) = (m²−1)·(m−1)m(2m−1)/6 grows as m⁵/3; RHS(m) = B(m)−2 grows as m⁵/10; leading coefficients disagree 7:3. Verified table for m=1..8: agreement only at m=3. Δ remains an independent ingredient; the Phase 4B observation is logged as a false lead. Search over alternative decompositions of 1/40 in project quantities found only Paper 2's canonical form 1/(|λ_3|·N(2)) = 1/(8·5) and trivial refactorings (e.g., 1/(d_max·N_init·N(2))); no structural alternative. (α-G) Paper 18 wording correction APPLIED: three edits to Section V (lines 764-766, 769-771, 808) replacing "Weyl exchange constant for S²" with "ω₂ = π, the 2D ball volume that serves as the numerator of the S² Weyl density and as the packing-plane prefactor σ₀ = πd₀²/2 in Paper 0". Numerically correct (S² Weyl constant is 1/(4π), not π); preserves the structural Paper-0 ↔ S² connection; conjecture status of Paper 2 unchanged. **Net Phase 4C status:** F is confirmed to live outside the discrete graph (must be derived from continuum S¹ regularization, not from any graph or Fock-weight construction at n_max=3); Δ is confirmed independent; K combination rule remains conjectural with no further reduction in unknowns. Data: `debug/data/track_alpha_phase4c/`. Phase 4D sprint complete (April 2026, single Track α-H), CLEAN NEGATIVE that closes the Hopf-fiber attack on F. **Continuum-fiber base⊗fiber tensor trace:** five variants computed (uniform fiber, scaled fiber, base × continuum S¹ zeta, Mellin–Barnes/heat-kernel, resolvent), none hit any K target within 1%. (a) Uniform: T_a = (63/4)·(π²/3) = 21π²/4 ≈ 51.81 vs B+F = 43.65 (18.7% gap). (b) Scaled L=2l+1: 409/16 ≈ 25.56 (39% gap from B). (c) Base × ζ_{S¹}, no Hopf weighting, averaged over 6 cells: 47π²/288 ≈ 1.611 vs F = π²/6 ≈ 1.6449 (2.1% — by far the cleanest, but does NOT decompose as a multiple of F). (d) THE CRITICAL VARIANT — Mellin/heat-kernel: derived an EXACT closed form I(1) = Σ_{l≥0} (2l+1)·T(1+l(l+1)) where T(b) = π²·csch²(π√b)/(2b) + π·coth(π√b)/(2b^{3/2}); numerical I(1) = 3.93747 (50 dps). The asymptotic expansion via the first 5 Weyl coefficients on S² gives I_asym = π·(1 + 1/6 + 1/20 + 1/42 + 1/72) = 2119π/1680 ≈ 3.9627. **STRUCTURAL INSIGHT:** every asymptotic term is LINEAR in π, not π² — because the Jacobi inversion θ_{S¹}(t) ~ √(π/t) introduces √π, which when paired with the Schwinger weight integral yields √π·√π·ℚ = π·ℚ at every order. The expected ζ_R(2) = π²/6 contribution collapses to π·(rational) at each order. The residual I(1) − I_asym = −0.02506 does NOT PSLQ-identify with {1, π, π², ζ(3), log 2, coth(π)} at 10⁻²⁵ tolerance; it's the genuinely transcendental tail of csch²/coth at √(1+l(l+1)) for l ≥ 0. (e) Resolvent: (63/4)·π·coth(π) ≈ 49.67 (13.8% gap from B+F). **The Phase 4D structural conclusion is the strongest result of the alpha sprint series so far:** F = π²/6 CANNOT enter K through any Hopf-fiber trace, discrete OR continuous. The Jacobi inversion mechanism systematically lowers the π power by ½ at each order, which explains both Phase 4C's all-rational fiber zetas (no π at all in finite path graphs) and Phase 4D's all-linear-π asymptotics (π but never π²). Combined with α-C (B is a pure (n,l)-base property, derivable in the Fock formalism), this exhausts the Hopf bundle as a derivation route for both B and F simultaneously: B can be derived, F cannot, period. **NOTE:** corrected a typo propagated from the Phase 4B summary above and into the Phase 4C/4D sprint plans — the α-C identity sum is 252/16 = 63/4, not 21/8 (21/8 is the per-cell average). The α-C source analysis (debug/data/track_alpha_phase4b/track_c_analysis.md) had the correct values throughout. Recommendation accepted: SHELVE Hopf-fiber trace attempts for F; future work should target either (i) S⁵ spectral zeta at s=2 (a pure base computation, not a tensor trace), (ii) Eisenstein-series / L-function attached to the (n,l) lattice, or (iii) accept F as a calibration exchange constant in Paper 18's taxonomy with no microscopic discrete origin. Data: `debug/data/track_alpha_phase4d/`. Phase 4E sprint complete (April 2026, single Track α-I), CLEAN NEGATIVE that closes the entire round-sphere Laplace–Beltrami avenue for F. **S⁵ spectral geometry test:** S⁵ does NOT contain F = π²/6 in any spectral invariant. (1) Spectral zeta ζ_{S⁵}(s) converges only for s > 5/2; values at s = 3,4,5 are 0.0822, 0.0110, 0.0020 — no hit on F = 1.6449 or 2F = 3.290. ζ_{S³}(2) = 0.8850 also no hit. (2) Truncated trace B_{S⁵}(ν_max) jumps 30 → 270 → 1320 → ... — never hits B = 42 (Paper 2's value is structurally a S³ quantity, not transferable to S⁵). The S⁵ analog of Paper 2's selection principle B/N = dim fails: B_{S⁵}/N_{S⁵} jumps 30/7 ≈ 4.29 (ν_max=1) past 5 (= dim S⁵) directly to 10 (ν_max=2), never satisfying B/N = 5 at any finite ν_max. (3) Slater V_ee on S³: confirmed all 145 analytical Slater integrals in `geovac/casimir_ci.py` are pure Python Fractions (max denominator 48,828,125 = 5¹⁰), F⁰(1s,1s) = 5/8 verified, and **the entire He FCI matrix is rational + sqrt-algebraic from ⟨n,l|1/r|n',l⟩ — π² does NOT appear in V_ee by construction**. Paper 18 classification: V_ee Slater integrals are intrinsic exchange constants, not calibration. (4) Heat kernel: vol(S⁵)/vol(S³) = π³/(2π²) = π/2 (single power of π, no π²). Seeley-DeWitt coefficients (round sphere, exact sympy): a₀=1, a₁(S³)=1, a₁(S⁵)=10/3, a₂(S³)=1/2, a₂(S⁵)=16/3 — all pure rationals. Spectral determinants on odd-dimensional spheres reduce to ζ_R(3), ζ_R(5), log 2, NOT to π^even. (5) S⁵/S³ fiber contribution: ζ differences {−0.092, −0.041, −0.015} at s=3,4,5; truncated B-differences {18, 186, 996, 3756} pure integers — no target match anywhere. Cleanest near-miss is a₁(S⁵) = 10/3 ≈ 3.3333 vs 2F = π²/3 ≈ 3.2899 (rel err 1.32%, but structural coincidence between a rational and a transcendental, not an equality). **Structural insight (the headline result):** S⁵ produces (i) integer powers of π in volumes, (ii) pure rationals in heat kernel/Seeley-DeWitt/truncated traces, (iii) ζ_R(odd) + log 2 in spectral determinants. **π² = 2·ζ_R(2) never appears additively next to a rational.** This is fully consistent with Paper 24's HO rigidity theorem: the S⁵ Bargmann-Segal lattice is bit-exactly π-free in exact rational arithmetic, so if π² were natively present in any S⁵ QM construction, Paper 24 would have found it. Combined with Phases 4C (discrete Hopf S¹) and 4D (continuum Hopf S¹) negatives, this **closes the entire S³/S⁵ sphere-spectral avenue for F**. **Net Phase 4E status:** SIX mechanisms eliminated across Phases 4B-4E (Hopf-twist spectral comparison, higher Casimir traces on S³, Hopf graph quotient, discrete Hopf fiber, continuous Hopf fiber, S⁵ spectral geometry). F = ζ(2) is structurally NOT in any round-sphere Laplace-Beltrami construction — it must enter K through arithmetic/flat-lattice objects (Epstein zeta, Eisenstein E₂*, Dirichlet L at s=2) OR be accepted as a calibration exchange constant in Paper 18's taxonomy with no microscopic geometric origin. Recommendation accepted: shelve sphere-spectral attempts for F permanently. Data: `debug/data/track_alpha_phase4e/`. Phase 4F sprint complete (April 2026, single Track α-J), **POSITIVE PARTIAL — first exact identification of F = π²/6 as a graph-intrinsic Dirichlet series**. **Headline result:** D_{n²}(s) = Σ_{n=1}^∞ n²·n^{−s} = ζ_R(s−2) evaluated at s = 4 gives D_{n²}(4) = ζ(2) = π²/6 = F **EXACTLY** as a sympy symbolic equality. The weight g_n = n² is the Fock degeneracy of the n-th S³ shell (Paper 7 Sec VI). The exponent s = 4 is natural in three independent ways: (i) s = d_max from Paper 0's packing axiom, (ii) s = dim(ℝ⁴) = the ambient dimension of S³, (iii) s = 2·N_init = twice the packing initial count. This is the first positive identification of F from any graph-intrinsic quantity across the entire alpha sprint series — it softens Phase 4E's "calibration-only" reading of F and reclassifies F from "external transcendental injection" to "Dirichlet series of the Fock degeneracy lattice at s = d_max". (Subtask 1) Finite shell-lattice Epstein sanity baseline: only one trivial rational hit, Z_{n²}^{2l+1}(s=1) = 3 (the classical Fock trace = dim S³ = Paper 2's selection principle ratio); no F or B targets. (Subtask 2) D_B(s) = (1/2)[ζ(s−4) − ζ(s−2)] tabulation at s=4..12: D_B(6) = π²(15−π²)/180 ≈ 0.281 with 2·D_B(6) = F − π⁴/90 — leading term IS F but the π⁴/90 ≈ 1.08 correction is comparable; near-miss only, NOT a clean F identification. (Subtask 3 — KEY) variant Dirichlet series scan: D_{n²}, D_N, D_λ, D_lmax tested at s = 3..12; **D_{n²}(4) = ζ(2) = F is the ONLY exact symbolic equality** found across all variants and all integer/half-integer s in range. (Subtask 4) Selection principle in Dirichlet language: seeking ζ(s−4)/ζ(s−2) = 7 (the value that would give D_B(s)/D_{n²}(s) = 3) — NEGATIVE. The function is monotonic, bounded above by −4.68, asymptotes to −6 as s→∞; no sign change, no root. **Paper 2's B(3)/N(3) = 3 is a finite-truncation coincidence at m=3, NOT a Dirichlet-analytic statement.** (Subtask 5) Truncation correction: D_{n²}(4) truncated to n=1..3 is 49/36 ≈ 1.361; tail = π²/6 − 49/36 ≈ 0.2838. Tail/Δ ≈ 11.35, not a clean rational multiple. **Δ is NOT the n_max=3 truncation correction of the F-producing Dirichlet series.** **Net Phase 4F status (positive partial):** The three components of K = π(B + F − Δ) now have three separate identified origins: (1) B = 42 = finite truncated Casimir at m=3 (Paper 2 Eq. 17, with Phase 4B α-C providing the κ↔B link via Fock weight), (2) F = π²/6 = D_{n²}(d_max) = the infinite Dirichlet series of the Fock degeneracy at the packing exponent d_max = 4 (NEW), (3) Δ = 1/40 = still unknown — neither a Dirichlet truncation correction nor a function of B and N_init. The "calibration-only F" reading of Phase 4E is now obsolete: F has a graph-intrinsic origin tied to the packing axiom. The full Paper-2 closure bar (same construction giving all three of B, F, Δ) is NOT met — B comes from a finite truncation at m=3, F from an infinite series at s=4, and these are distinct mechanisms. **POSITIVE PARTIAL flagged for plan-mode review.** Paper 2 NOT auto-updated; PI to decide whether to (i) add the F = D_{n²}(d_max) identity as a remark in Paper 2 with conjecture status preserved, (ii) write it up in Paper 18's exchange-constant taxonomy as a new "arithmetic exchange constant" tier between intrinsic and calibration, or (iii) wait for Δ to be derived before any paper update. Recommendation from α-J: focus the next alpha sprint on **Δ = 1/40 alone** (try Eisenstein/Dedekind η at SL(2,ℤ) cusps, or the 1/40 = 1/(2·4·5) factorization in (N_init, d_max, l_max+2) coordinates) — Δ is the last unknown in K. Data: `debug/data/track_alpha_phase4f/`. Phase 4G sprint complete (April 2026, single Track α-K), CLEAN NEGATIVE on the arithmetic-unification hypothesis for Δ. (Subtask 1) Six finite/infinite B-F overlap candidates tested as regularization artifacts (F-partial-sum 49/36, B-summand at s=4 = 59/72, selection ratio 3, last-B-term-at-s=4 = 4/9, Dirichlet pair, F-tail π²/6 − 49/36) — all at least 10× away from 1/40. No clean overlap term reproduces K/π. (Subtask 2) Packing factorization: Paper 2's canonical form 1/Δ(m) = (m²−1)·m(m−1)(2m−1)/6 verified symbolically as equivalent to m(m−1)²(m+1)(2m−1)/6 and to |λ_m|·N(m−1). Computed B(m)·Δ(m) = 3(2m+1)(m+2)/[10(m−1)(2m−1)] symbolically — clean rational function but NOT constant in m (B·Δ = 2 at m=2, 21/20 at m=3, 27/35 at m=4, 77/120 at m=5). No cleaner factorization exists. (Subtask 3 — KEY) Brute-force ζ-combination scan, 3,228 candidates at 50 dps over ζ_R(n), 1/ζ_R(n), differences, ratios, products, (ζ_R(a)−1) tail combinations, π^k. **ZERO strict hits** within 10⁻²⁵ of 1/40. Cleanest non-tautological near-misses (all ≥ 1% relative error, no project meaning): (1/8)(ζ(3)−1) ≈ 0.02526 (1.03%), (3/10)(ζ(4)−1) ≈ 0.02470 (1.21%), 1/(4π²) ≈ 0.02533 (1.32%), (2/3)(ζ(5)−1) ≈ 0.02462 (1.53%) — accidental rational approximants to 0.025 with no structural interpretation. (Subtask 4) Laurent expansion of ζ(s−2) at s=4 and at the s=3 pole: ζ(2), ζ'(2)≈−0.9375, ζ''(2)≈1.989, ζ'''(2)≈−6.000, Stieltjes constants γ_0..5 — all combinations checked, closest 1−ζ'(2)² ≈ 0.121 (3.8× off). (Subtask 5) Hurwitz ζ(s,a) for s∈{2..6}, a∈{1..6}: closest ζ(3,5) ≈ 0.02439 at 2.4% off (no structural origin for a=5). Bernoulli check: B_6 = 1/42 is the nearest (4.76% off) — a numerological coincidence with B=42, NOT 1/40. (Subtask 6) Finite-N interpretation verified: 1/Δ(m) = |λ_m|·N(m−1), with |λ_3|=8 = the S³ Laplace-Beltrami gap to the next unused shell, N(2)=5 = cumulative state count in shells n=1,2. Both factors are intrinsic (n,l)-lattice data with no further arithmetic content. **STRUCTURAL CONCLUSION (the headline result of Phases 4B-4G):** The three components of K = π(B + F − Δ) have CATEGORICALLY DIFFERENT origins and there is NO unifying generator: (1) **B = 42** = finite Casimir trace at m = 3 — RATIONAL, COMBINATORIAL, finite truncation of the Fock degeneracy lattice (Phase 4B α-C established the κ↔B Fock-weight link); (2) **F = π²/6** = D_{n²}(s = d_max = 4) = ζ_R(2) — TRANSCENDENTAL, ARITHMETIC, infinite Dirichlet series of the Fock degeneracy at the packing exponent (Phase 4F α-J); (3) **Δ = 1/40** = |λ_{n_max}|·N(n_max−1) = (gap above cutoff) × (states below cutoff) — RATIONAL, COMBINATORIAL, finite boundary-mass invariant of the truncation (Phase 4G α-K established as irreducible). The K combination rule is therefore a sum of three categorically different objects with no common arithmetic generator: B is a finite sum, F is an infinite sum, Δ is a boundary product. K = π(B + F − Δ) is a genuine **STRUCTURAL MYSTERY** rather than a common-generator identity — the open question has shifted from "how is each piece derived" (each piece IS derived structurally now) to "why does the additive combination of these three structurally distinct objects equal α⁻¹ at 8.8×10⁻⁸." **Net status of the K combination rule after Phase 4G:** Seven mechanisms eliminated across Phases 4B-4G (six sphere-spectral plus one arithmetic-Dirichlet-for-Δ). Two positive structural identifications (κ↔B Phase 4B, F = D_{n²}(d_max) Phase 4F). One paper correction applied (Paper 18 wording, Phase 4C). All three K ingredients now have structural origins. **Phase 4G recommendation flagged for plan-mode review:** the α-K agent recommends updating Paper 2 to state explicitly that B, F, Δ have categorically different origins and the K combination rule is a structural mystery rather than a common-generator identity. This is a framing change to a conjectural paper; not auto-applied. PI to decide whether to add this remark to Paper 2 or to extend Paper 18's exchange-constant taxonomy with the three-tier (combinatorial-rational / arithmetic-transcendental / boundary-rational) classification. Data: `debug/data/track_alpha_phase4g/`. Phase 4H sprint complete (April 2026, six-track SM-origin sprint Tracks SM-A/B/C/D/E/F), **FIVE NEGATIVES + ONE STRUCTURAL POSITIVE PARTIAL**. Hypothesis tested: Δ = 1/40 is the one-loop QED vacuum-polarization shift between a "bare topological coupling" and the physical α(m_e), encoding SM particle content. PI noted the striking coincidence Σ_f N_c Q_f² = 8 = |λ_3| over three SM generations. (SM-A) Standard one-loop QED running scan: NEGATIVE — at every geometrically natural UV scale (Bohr momentum, 2 Ry, 1/r_e, m_μ, M_Z, M_Planck), Δ(1/α) misses both targets 1/40 and π/40 by ≥45×; the inverse-solve μ values that hit the targets (574.9 MeV for 1/40, 739.9 MeV for π/40) sit in the empty desert between m_e and m_π with no recognizable physics correspondence. SM running cannot reproduce 1/40 from any natural scale. (SM-B) Σ_f N_c Q_f² = 8 = |λ_3| structural map: NEGATIVE — the per-generation contribution 8/3 is constant across all three SM generations (charge-universal) while every per-shell S³ invariant (|λ_n|, g_n, |λ_n|/n) varies with n. No shell→generation map is consistent with charge-universal per-gen value. The 8 = 8 equality is a numerical coincidence with no representation-theoretic content. Other natural SM charge traces give different integers (Σ N_c Y² = 10 Weyl, 5 Dirac; Σ N_c Q⁴ = 130/27); only Σ N_c Q² hits 8, consistent with selection bias. (SM-C) U(1)_Hopf → SU(5) embedding for 1/40 as topological invariant: NEGATIVE — 1/40 does not appear naturally in SU(5)/E_8/Spin(10) as Chern number, η-invariant, or Chern-Simons level. Closest hit is the dual-Coxeter rewrite 1/40 = 1/(h^∨(SU(5))·h^∨(SO(10))) = (3/8)/15 = sin²θ_W^GUT / dim(5̄⊕10), which is a tautological renaming using the same integers 5 and 8 already present in Paper 2. The naive (n=1,n=2) shell ↔ 5̄ matching fails at the SU(3)×SU(2) branching level (1+1+3 vs 3+2). (SM-D) Vacuum polarization on S³ at n_max=3: **POSITIVE PARTIAL — the headline structural result of the sprint**. The genuine new finding: $g_3^{\text{Dirac}}(S^3) = 2 \cdot 4 \cdot 5 = 40 = \Delta^{-1}$ EXACTLY. The single-chirality degeneracy of the third Dirac eigenmode on the unit S³ (Camporesi-Higuchi spectrum $|\lambda_n|=(2n+3)/2$, $g_n=2(n+1)(n+2)$) is exactly 40 only at $n=3$ — matching Paper 2's selection-principle cutoff. This is structurally cleaner than Paper 2's current factorization $\Delta^{-1} = |\lambda_3| \cdot N(2) = 8 \cdot 5$: the Dirac decomposition is a single polynomial $g_n = 2(n+1)(n+2)$ rather than a product of two distinct objects. The perturbative one-loop self-energy itself does NOT hit 1/40 (overshoots by 45-530× at all tested normalizations), but the *combinatorial* mode-count factorization $\Delta = 1/g_3^{\text{Dirac}}$ is exact and suggests Δ is a state-counting boundary invariant (reciprocal of the charged-spinor mode count at the truncation edge), not a perturbative shift. Cross-validated against SM-A's inverse solve to 6 digits ($\mu/m_e = 1.125$). Sign check: Δ enters K with minus, corresponding to compactification *increasing* 1/α from reduced low-q photon screening — qualitatively correct for QED. (SM-E) HO/nuclear consistency check: POSITIVE — Paper 24's HO Bargmann-Segal graph is bit-exactly π-free in rational arithmetic, Paper 23's deuteron/He-4 Hamiltonians have zero π and zero calibration constants (only ℏω + Minnesota Gaussians + Moshinsky-Talmi rationals). The universal/Coulomb-specific partition holds: Δ is absent from the HO/nuclear sector, consistent with electromagnetic origin. (SM-F) Higher Hopf S⁴/S⁸ spectral invariants: NEGATIVE — the selection principle B/N = dim does not transfer (2L(L+4)/3 and 4L(L+8)/5 never hit integer dimensions for natural cutoffs). Computed B/N = 30 at L=5 on S⁴ (vs 1/α_W ≈ 29.5, 1.7% off — Diophantine integer hit, not a spectral identity); B/N = 8 at L=2 on S⁴ (vs 1/α_s(M_Z) ≈ 8.467, 5.5% off). F-analogs are mixed sums of multiple Riemann zetas with no closed-form simplification. The complex Hopf S¹→S³→S² is structurally special: only S³ has Fock projection from 1e Hamiltonian, only there does the n² degeneracy collapse F to ζ_R(2) = π²/6, only there does the B/N quadratic hit an integer. **STRUCTURAL CONCLUSION (Phase 4H headline):** The SM-running hypothesis for Δ is FALSE. But Δ has a cleaner spectral home than previously known: $\Delta^{-1} = g_3^{\text{Dirac}}(S^3) = 40$ — the third single-chirality Dirac mode degeneracy on unit S³. This is a structural identification at the same level as Phase 4F's F = D_{n²}(d_max), i.e. it gives Δ a single canonical formula but does NOT derive the K combination rule. **Eight mechanisms now eliminated across Phases 4B-4H** (six sphere-spectral, one arithmetic-Dirichlet, one SM-running). **Three positive structural identifications** (κ↔B Phase 4B, F = D_{n²}(d_max) Phase 4F, Δ = 1/g_3^Dirac Phase 4H). All three K ingredients now have clean spectral origins; the K combination rule remains a structural mystery — the open question is unchanged: "why does $\pi(B + F - 1/g_3^{\text{Dirac}})$ equal α⁻¹ at $8.8 \times 10^{-8}$." **Phase 4H recommendation flagged for plan-mode review:** the SM-D agent recommends reframing Δ in Paper 2 from $|\lambda_3| \cdot N(2) = 8 \cdot 5$ to the cleaner Dirac-degeneracy form $g_3^{\text{Dirac}} = 2(n+1)(n+2)|_{n=3}$. This is a framing change to a conjectural paper; not auto-applied. PI to decide. The SM-running hypothesis is now a documented dead end — see Section 3. Data: `debug/data/track_alpha_sm/`. Phase 4I sprint complete (April 2026, five-track Dirac-on-S³ Tier 1 sprint D1/D2/D3/D4/D5), ALL NEGATIVE with three named closed-form structural obstructions. Tests whether the Dirac sector on S³ lifts B, F, or Δ to a common-generator spectral invariant. (D1) Infrastructure module `geovac/dirac_s3.py` with Camporesi-Higuchi spectrum (|λ_n| = n+3/2, g_n^Dirac = 2(n+1)(n+2)), spinor label dataclass, and bit-exact π-free certificate (51/51 tests in exact sympy rational arithmetic, analog of Paper 24 Bargmann-Segal certificate). (D2) Dirac analog of B = 42: NEGATIVE. The scalar Laplacian has a zero mode at (n=1, l=0) whose Casimir weight l(l+1) vanishes, giving B(m) = m(m-1)(m+1)(m+2)(2m+1)/20 with factor (m-1). The Dirac spectrum is gapped (|λ_n| ≥ 3/2), every level contributes to any truncated trace, and no cumulative Dirac/Weyl sum carries the (m-1) factor. Single-point coincidence at m=3: |λ_{m-1}|·g_{m-1}^Weyl = 42 via (m-1)(m+2)=10, uniquely at m=3 and at no other integer — elementary arithmetic coincidence, not spectral identity. (D3) Dirac analog of F = π²/6: NEGATIVE. At the packing exponent s=d_max=4, D_{g_m^Dirac}(4) = Σ 2m(m+1)/m^4 = 2ζ(2) + 2ζ(3) = π²/3 + 2ζ(3). Apéry's theorem gives ζ(2), ζ(3) Q-linearly independent; isolating F requires hand-subtracting 2ζ(3), equivalent to projecting onto the m² sub-weight (i.e. returning to the scalar case). Hurwitz-regularized form D_dirac^CH(4) = π² - π⁴/12 (closed form, inseparable). F lives on the scalar Fock-degeneracy Dirichlet at d_max, not on the Dirac. **POSITIVE BYPRODUCT:** first appearance of ζ(3) as a structural transcendental in the framework, via the weight-m subchannel of the Dirac degeneracy. (D4) Orduz S¹-equivariant Hopf decomposition: NEGATIVE. The 40 states at n_CH=3 decompose under the Hopf U(1) action as a clean 20⊕20 charge-parity split (4 half-integer charges × 5 mult + 5 integer charges × 4 mult), structurally cleaner than Paper 2's scalar factorization Δ⁻¹ = |λ_3|·N(2) = 8·5 — but does NOT produce B or F. Numerical near-miss Σ |λ_n|²·g_n^Dirac = 378 = 9·B at n_CH≤2 drops at n_CH≤3 (1188, not a multiple of 42), so the factor 9 is a coincidence. Categorically distinct from Phase 4B α-D (scalar Hopf quotient): the Dirac decomposition is a first-order operator decomposition with nontrivial charge-q twists, unavailable in the bosonic sector. Rules out the Hopf-equivariant Dirac as common generator. (D5) Three-tier coincidence formally documented. **STRUCTURAL CONCLUSION:** B, F, Δ have canonical spectral homes (B: scalar Casimir trace at m=3; F: scalar Fock Dirichlet at d_max; Δ: single-level Dirac degeneracy at n=3). The three homes span TWO sectors of S³ (scalar and spinor) and THREE object types (finite trace, infinite Dirichlet, single-level degeneracy). K = π(B + F - Δ) is confirmed as a genuine cross-sector sum with no common spectral origin. NINE mechanisms eliminated across Phases 4B-4I (six sphere-spectral, one arithmetic-Dirichlet for Δ, one SM-running, one Dirac-on-S³). FOUR positive structural identifications (κ↔B Fock link α-C; F = D_{n²}(d_max) α-J; Δ⁻¹ = g_3^Dirac SM-D; ζ(3) natively in Dirac sector D3 — new framework transcendental, NOT part of K). **Recommendation per sprint plan decision table:** α structural-decomposition program closed; Paper 2 §IV rewrite applied (Tier 1 closure); NO Tier 1b proof sprint opened; move to Tier 2 (spin-ful composed qubit encoding, heavy-element relativistic chemistry) which builds on D1 infrastructure but is independent of the α program. Paper 18 taxonomy entry added: distinction between even-zeta content (second-order operators, Jacobi-θ inversion) and odd-zeta content (first-order operators, half-integer Hurwitz) — ζ(3) is the first odd-zeta framework transcendental. Data: `debug/data/dirac_d{2,3,4}_*.json`; `docs/dirac_s3_verdict.md`; `docs/paper2_section4_rewrite.tex`.
- Dirac-on-S³ Tier 2 sprint complete (April 2026, seven-track sprint T0/T1/T2/T3/T4/T5/T6), ALL POSITIVE. Engineering upgrade: spin-ful composed qubit Hamiltonians for LiH, BeH, and CaH (first instance of a relativistic molecule in the composed pipeline), built end-to-end algebraically in the (κ, m_j) native Dirac-label basis. (T0) d_spinor(l_max) table in both pair-diagonal-m and full-Gaunt conventions, l_max=0..5: ratio d_spinor/d_scalar = 1/4 (l_max=0, pure spin-dilution) → 11/20 → 553/724 → 101/118 → 2533/2820 → 9611/10396 ≈ 0.92 (l_max=5). Paper 22 theorem extends verbatim to spinor basis with sparsity-exponent preserved and prefactor inflated. (T1) `geovac/dirac_matrix_elements.py`: DiracLabel dataclass, κ ↔ (l, σ) bridge, all Szmytkowski angular matrix elements in closed form, all diagonal radial ⟨r^k⟩ and ⟨1/r^k⟩ as Bethe-Salpeter rationals; off-diagonal radial via direct sympy integration fallback; Kramers-Pasternak direct integration for ⟨r⁻²⟩ (all states) and ⟨r⁻³⟩ (|κ|≥2) added v2.18.1. 108 tests pass; 51/51 D1 regression tests preserved. (T2) `geovac/spin_orbit.py`: Breit-Pauli H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)], closed form in (n, κ) with exact Kramers cancellation at l=0 (κ=−1 forces numerator zero before the formally divergent ⟨1/r³⟩ is evaluated). Z⁴ scaling verified symbolically across Z ∈ {1, 3, 4, 38}; 2p doublet splitting α²/32 reproduces. 22 tests pass. (T3) `geovac/composed_qubit_relativistic.py`: dispatches from `build_composed_hamiltonian(spec)` when spec.relativistic=True. Three specs `lih_spec_relativistic`, `beh_spec_relativistic`, `cah_spec_relativistic` (BeH and CaH are one-bond reductions of BeH₂/CaH₂). Pauli ratio rel/scalar isostructural 1.00×/2.42×/5.89× at n_max=1/2/3 across all three molecules. λ_ni flat or decreasing vs scalar (QPE-favorable). 13 T3 regression tests + 164 pre-existing tests pass; scalar LiH/BeH₂/H₂O counts 334/556/778 bit-exact. (T4) Sunaga 2025 (PRA 111, 022817) head-to-head: only the RaH-18q cell (47,099 Pauli) is published in the main paper; per-molecule BeH/MgH/CaH/SrH/BaH at 18q are in SI Tables S1–S3, flagged DEFERRED. GeoVac native-Q ratios vs Sunaga RaH-18q: LiH 0.017×, BeH 0.017×, CaH 0.011× (at GeoVac Q ∈ {20, 30}). Projected matched-Q=18 advantage 150×–250× via Paper 14 §IV.B O(Q^2.5) scaling × rel/scalar 2.4×. Fine-structure sanity: He 2³P span −66%, Li 2²P doublet +211%, Be 2s2p ³P span −78% relative error — sign and OoM correct across all three atoms, 20–50% accuracy target NOT met (known limitation of leading-order Zα² + Slater Z_eff; Tier 3+ to lift via Darwin + mass-velocity + multi-electron SS/SOO). (T5) `geovac/spinor_certificate.py`: `verify_spinor_pi_free` walks the sympy expression tree of each T3 coefficient and enforces membership in the Paper 18 spinor-intrinsic ring R_sp := ℚ(α²)[γ]/(γ²+(Zα)²−1). T3's symbolic H_SO block passes cleanly for n_max ≤ 4 and across all κ branches; six negative-control contamination tests (bare π, π², ζ(3), log(Z), E₁, unregistered symbol) reject cleanly. 25 tests pass. Paper 18 §IV subtier drop-in applied (this commit). (T6) Paper 14 §V spinor-composed section, Paper 22 spinor-block angular sparsity section, Paper 20 Tier 2 relativistic resource table, and this CLAUDE.md bullet — all applied in this commit. **STRUCTURAL CONCLUSION:** The Paper 18 operator-order × bundle grid now has three of four cells populated (2nd-order/scalar = calibration π, 1st-order/scalar = Tier-1 odd-zeta, 1st-order/spinor = Tier-2 spinor-intrinsic α²+γ); 2nd-order/spinor filled by T9 theorem (degenerate with scalar calibration π^{even}) and confirmed by one-loop QED on S³ (v2.18.2). Tier 3 roadmap: [Kr] frozen-core for SrH/RaH, Martínez-y-Romero γ radial corrections, Darwin + mass-velocity, multi-electron SS/SOO. Data: `debug/dirac_t0_memo.md`, `docs/dirac_matrix_elements_design_memo.md`, `docs/spin_orbit_design_memo.md`, `docs/spin_ful_composed_design_memo.md`, `docs/tier2_market_test.md`, `docs/tier2_t5_verdict.md`, `docs/tier2_verdict.md`.
- Dirac-on-S³ Tier 3 sprint complete (April 2026, three-track sprint T7/T8/T9), two positive engineering results + one structural theorem. (T7) gamma = sqrt(kappa²-(Z*alpha)²) radial corrections exact for n_r=0 states via Pochhammer ratios; <1/r> exact for all states via Hellmann-Feynman; NR limits verified; R_sp ring confirmed. Full n_r>=1 for s=-2,-3 deferred (Kramers-Pasternak recursion) — PARTIALLY RESOLVED v2.18.1: direct integration gives ⟨r⁻²⟩ for all states, ⟨r⁻³⟩ for |κ|≥2; p₁/₂ limitation remains. 31 new tests pass. (T8) Darwin + mass-velocity completing the alpha^4 fine-structure ladder: E_SO + E_D + E_MV = -(Z*alpha)^4/(2n^4)*[n/(j+1/2)-3/4] verified as exact sympy symbolic equality for all 16 states through n=4. All 6 Dirac accidental degeneracies confirmed. HONEST NEGATIVE: Darwin+MV do not improve He/Li/Be 2p-doublet splittings (both states share l=1: Darwin=0 for l>=1, MV identical for same l). Residual 66-211% errors trace to multi-electron SS/SOO (Direction 3, deferred). 43 new tests pass. (T9) HEADLINE: squared Dirac D² spectral zeta theorem. zeta_{D²}(s) = 2^{2s-1}*[lambda(2s-2) - lambda(2s)] where lambda(2k) = (1-2^{-2k})*zeta_R(2k) = rational*pi^{2k}. At every integer s, zeta_{D²}(s) is a two-term polynomial in pi² with rational coefficients. No odd-zeta content at any s -- theorem, not numerical observation. Mechanism: squaring maps odd^{-s} to odd^{-2s}, always pi^{even} by Bernoulli. Paper 18's 4th cell (2nd-order x spinor-bundle) is FILLED and DEGENERATE with scalar calibration (pi^{even}). Operator order is the primary transcendental discriminant; bundle type modulates coefficients but not the transcendental class. The taxonomic grid collapses to 3 effective tiers. Data: `debug/dirac_t7_memo.md`, `debug/dirac_t8_memo.md`, `debug/dirac_t9_memo.md`, `docs/tier3_verdict.md`.
- Sprint 7 balanced second-row FCI (v2.19.4, April 2026): NaH (Q=20, 2e) and MgH₂ (Q=40, 4e) at n_max=2 balanced coupled FCI — HONEST NEGATIVE. Both show monotonic overattraction with no equilibrium. Root cause: frozen [Ne] core means valence orbitals see full Z=11/12 nuclear charge from cross-center V_ne without adequate core screening. LiH succeeds because it has an explicit core block. Structural limitation of balanced + frozen-core architecture for PES at second-row Z. Key files: `debug/sprint7_balanced_second_row.py`, `debug/data/sprint7_balanced_second_row.json`.
- Sprint 7b screened SO from FrozenCore Z_eff(r) (v2.19.4, April 2026): Fixed frozen-core SO underestimation. Three-part fix: (1) OrbitalBlock gains Z_nuc_center (nuclear charge) and n_val_offset (physical n mapping); (2) spin_orbit.py Z_wfn parameter separates operator Z from wavefunction Z; (3) relativistic builder uses screened_xi_so() which solves the radial Schrodinger equation with FrozenCore Z_eff(r) screening and computes ⟨(1/r)dV/dr⟩ from the actual screened potential. Core penetration enhancement: Ca 12×, Sr 61×, Ba 144× over hydrogenic Z_eff=2. Results: CaH 37 cm⁻¹ (-69% vs ~120 physical), SrH 178 cm⁻¹ (-70% vs ~600), BaH 433 cm⁻¹ (-61% vs ~1100). Remaining error from leading-order Breit-Pauli + Clementi-Raimondi screening approximation. Key files: `geovac/neon_core.py` (screened_r3_inverse, screened_xi_so, screened_so_splitting), `geovac/composed_qubit_relativistic.py` (screened SO path), `debug/sprint7_fine_structure.py`.
- Sprint 5 Track CP (April 2026): **closed Sprint 3 BF-E Li/Be honest negatives** to <20% via convention fix + Slater's rules. Li 2²P splitting: +8.89% err via standard single-particle SO formula zeta = alpha^2/2 * Z_val * <1/r^3>_2p with Z_val=1 (asymptotic Li 2p sees [He] core) and Z_eff=1 (full-shield hydrogenic). Be 2s2p ^3P span (E(P_0) - E(P_2)): +2.76% err via std conv, Z_val=1 (Be 2p sees [He]+2s at large r), Slater Z_eff=1.95 (0.85 per 1s + 0.35 for same-shell 2s); all three individual splittings <20% (+18.9%, +6.3%, +2.8%). Root cause of Sprint 3 negatives: Z_nuc/Z_val factor-of-Z overcount in BR_C convention (coincidentally correct for He where Z_nuc=2=2*Z_val). Core polarization (Migdalek-Bylicki PRA 57, 3456, 1998, alpha_d=0.192, r_c=0.55 a.u.) WORSENS Li from 8.9% to 25% -- attractive potential increases zeta_2p in wrong direction. Multi-config 2s2p ↔ 2p^2 mixing FORBIDDEN BY PARITY (odd vs even product parity; 1/r_12 is parity-even). 5 new tests in `tests/test_breit_integrals.py`. Paper 14 §V updated with closure narrative. Files: `debug/cp_li_core_polarization.py`, `debug/cp_be_multiconfig.py`, `debug/cp_fine_structure_memo.md`, `debug/data/cp_{li,be}_results.json`.
- Breit SS/SOO composed pipeline (v2.18.0): two-body Breit-Pauli spin-spin and spin-other-orbit interactions wired into the relativistic composed builder via `include_breit=True` kwarg in `build_composed_hamiltonian()`. Breit retarded radial integrals r_<^k / r_>^{k+3} computed by `geovac/breit_integrals.py` using exact Fraction/sympy arithmetic. Pauli counts unchanged (same Gaunt angular selection rules as Coulomb ERI). 1-norm shifts O(10⁻⁴) Ha (α²-suppressed). He 2³P J-pattern correct: Drake combining coefficients A_SS, A_SOO from retarded integrals reproduce physical splittings at O(α²Z⁴). Block diagonality preserved. Hermiticity verified. Breit vanishes identically at α=0. Scalar path ignores `include_breit` entirely (backward compatible). 11 new tests in `tests/test_breit_composed.py`. Key files: `geovac/composed_qubit_relativistic.py` (Breit ERI loop), `geovac/breit_integrals.py` (retarded radial), `geovac/composed_qubit.py` (kwarg passthrough).
- Paper 18 compactness thesis reframe (v2.18.1): abstract, introduction, and conclusion rewritten around "compactness is fundamental" — the S³ graph arises because compact Lie groups have discrete spectra, and the exchange constant taxonomy classifies departures from compactness. New Peter-Weyl subsection in Section III deriving the graph Laplacian from harmonic analysis on compact groups. Paper 0 U(1) construction-level argument added. String compactification structural parallel noted. ~115 lines added. Paper 18 §IV motivic ζ(3) subsection added: operator-order → motivic-weight correspondence (2nd-order operators produce ζ(2k) = periods of ℚ(−1)^k; 1st-order Dirac produces ζ(3) = period of ℚ(−3)). Five geometric sources of ζ(3) as periods. Zagier dimension conjecture as organizing principle. QED loop-order prediction: odd-zeta absent at one loop (T9 theorem), expected at two loops.
- Kramers-Pasternak direct integration (v2.18.1): `dirac_radial_expectation_direct()` in `geovac/dirac_matrix_elements.py`. Confluent hypergeometric polynomial expansion of Dirac-Coulomb radial wavefunctions, direct sympy integration of r^s weighted products. ⟨r⁻²⟩ works for ALL states (all n_r, all κ). ⟨r⁻³⟩ for |κ|≥2 states. p₁/₂ (|κ|=1, n_r≥1) limitation documented (sympy integral does not converge for near-singular integrand). Resolves T7 deferred item for most states. 12 new tests, 108 total in `tests/test_dirac_matrix_elements.py`.
- One-loop QED vacuum polarization on S³ (v2.18.2): new module `geovac/qed_vacuum_polarization.py`. Seeley-DeWitt heat kernel coefficients a₀=a₁=√π, a₂=√π/8 on unit S³ (exact sympy). Vacuum polarization coefficient 1/(48π²). β(α) = 2α²/(3π) — reproduces standard QED beta function from S³ spectral data alone (Camporesi-Higuchi spectrum + heat kernel expansion). No odd-zeta at one loop confirmed (T9 theorem). Transcendental classification: all one-loop quantities are calibration π (even-zeta only). Two-loop computation in progress (ζ(3) prediction from Track D3 weight-m subchannel). 33 new tests in `tests/test_qed_vacuum_polarization.py`. Key files: `geovac/qed_vacuum_polarization.py`, `debug/qed_s3_memo.md`.
- Two-loop QED vertex parity and Catalan's constant (v2.19.1): new module `geovac/qed_vertex.py`. The two-loop sunset diagram on S³ with vertex parity selection rule (n₁+n₂+n_γ odd, from γ^μ coupling) splits the Dirac Dirichlet series into even/odd sub-sums exposing Catalan's constant G = β(2) and Dirichlet β(4). Exact Hurwitz decomposition: D_even(4) = π²/2 − π⁴/24 − 4G + 4β(4); D_odd(4) = π²/2 − π⁴/24 + 4G − 4β(4); D(4) = π² − π⁴/12 (G and β(4) cancel in full sum). Mechanism: vertex parity acts as Dirichlet character χ₋₄ on the mode sum; quarter-integer Hurwitz shifts ζ(s,3/4) and ζ(s,5/4) produce L(s,χ₋₄) in the difference. This reveals a three-axis taxonomy: (1) operator order (1st vs 2nd) determines motivic weight; (2) bundle type (scalar vs spinor) modulates coefficients; (3) diagram topology (vertex parity) introduces Dirichlet L-function characters. The Dirichlet β-values are a genuinely new transcendental class from the INTERACTION vertex, not from free propagators. Verified by PSLQ at 80-digit precision. 35 new tests in `tests/test_qed_vertex.py`. Key files: `geovac/qed_vertex.py`, `debug/qed_vertex_memo.md`. Paper 18 §IV updated with new subsection "Dirichlet L-values from vertex topology."
- S_min irreducible constant confirmed (v2.19.4, extended v2.26.1): the CG-weighted two-loop sunset sum S_min = Σ_{k≥1} T(k)² = 2.47993693803422255441... (200+ verified digits at 250 dps via mpmath.nsum Levin u-transform) is confirmed irreducible — 15 PSLQ attempts across 100+ basis elements (including polylogarithms, Euler sums, Tornheim-Witten zeta, colored MZVs, double Hurwitz at quarter-integer shifts, polygamma products, and cross-products of all known constants) ALL FAILED. S_min is a depth-2 multiple Hurwitz zeta at half-integer shift a=3/2, intrinsic to the Dirac spectrum on S³. Lives at the intersection of all three Paper 18 taxonomy axes (operator order + vertex topology + CG weighting). Key files: `debug/smin_extended_pslq.py`, `tests/test_smin_extended.py` (8 fast + 2 slow tests).
- Three-loop QED on S³ (v2.19.4, factorized v2.26.1): new module `geovac/qed_three_loop.py`. Iterated sunset topology: three electron lines connected by two photon lines in chain, with SO(4) vertex selection at each vertex. **O(N⁵)→O(N³) factorization (Track Q-2):** chain topology decomposes as Total = Σ_{n2} L(n2)·g(n2)/λ(n2)^a·R(n2) where L(n2) and R(n2) are independent one-loop partial sums over (n1,q1) and (n3,q2) respectively. Factorized matches direct to 1e-25 at n_max=15; n_max=50 in <60s; n_max=100 computed. CG-weighted convergence: power-law α≈−1.31 (too slow for PSLQ — 0 reliable digits at n_max=100 even with Euler-Maclaurin tail correction). Unrestricted factorizes as D(4)³ (purely π^{even}). CG-weighted sum EXCEEDS unrestricted at n_max≥15 (W=2 amplification). Depth-3 MZV prediction structurally motivated but computationally unverified — needs analytical decomposition or dramatically faster convergence. 28 tests (25 fast, 3 slow) in `tests/test_qed_three_loop.py`. Key files: `geovac/qed_three_loop.py`.
- QED self-energy and vertex correction on S³ (v2.26.1, Track Q-1): new module `geovac/qed_self_energy.py`. One-loop electron self-energy Σ(n_ext) from spectral mode sums with SO(4) vertex selection rules. **Self-energy structural zero theorem:** Σ(n_ext=0) = 0 exactly — vertex parity forces n1+n2+q odd with n2=n_ext=0, requiring 2n1 odd, which is impossible. This is a selection-rule protection of the ground state, not a cancellation. Self-energy table: Σ(0)=0, Σ(1)=−0.389, Σ(2)=−0.215, Σ(3)=−0.144. Vertex correction infrastructure for anomalous magnetic moment extraction. Schwinger convergence toward α/(2π) at ~33% level (n_max=20). **Anomalous magnetic moment curvature expansion (v2.26.1):** F₂/Schwinger = 1 + c₁/λ² + c₂/λ⁴ + ... with c₁ = R/12 = 1/2 (Parker-Toms) and c₂ = (2 - BΔ - FΔ - F/B)/5 = 19/100 - 41π²/25200 ≈ 0.17394231, matching high-precision n_ext=1 spectral sum (n_int=0..50 + tail correction, 13-digit delta) to 8 digits (residual 1.9e-8). All three Paper 2 invariants (B=42, F=π²/6, Δ=1/40) appear in c₂ — first concrete bridge between Paper 2 (α construction) and Paper 28 (QED on S³). **c₃ = -5.946(3)×10⁻⁷ from n_int=0..50 (200.9 sigma, NONZERO)**: the expansion does not terminate at second order; |c₃x³|/|c₂x²| ≈ 5.4×10⁻⁶ indicates rapid but non-terminating convergence. PSLQ finds no identification in the T9 basis (integer coefficients ≤10⁴). T9 consistent: c₂ = rational + rational·π² only; c₃ constrained to π^{even} ring. Multi-n_ext extraction obstructed by j-dependent form factors at n_ext≥2 (j=1/2 minority component shrinks as 2/n²). Key files: `debug/g2_c2_verification.py`, `debug/g2_c2_diagnostic.py`, `debug/data/g2_c2_verification.json`, `debug/g2_c3_fast.py`, `debug/data/g2_c3_investigation.json`. 35 tests in `tests/test_qed_self_energy.py`. Key files: `geovac/qed_self_energy.py`.
- Sommerfeld fine-structure Dirichlet sums D₂–D₆ (v2.26.1): exact closed-form coefficients c_p(n) extracted from the Dirac-Coulomb energy expansion through order (Zα)¹², verified symbolically for n=1..15 at each order. D₅ = (497/128)ζ(8) − (467/16)ζ(9) + (385/32)ζ(3)ζ(6) + (75/8)ζ(4)ζ(5), identified by PSLQ at 250 dps. D₆ = (-2289/512)ζ(10) + (1589/16)ζ(11) − (1617/64)ζ(3)ζ(8) − (1785/64)ζ(4)ζ(7) − (2065/64)ζ(5)ζ(6), assembled analytically from all 7 weight-11 Euler sum decompositions (each found by sub-basis PSLQ at 400 dps), verified to 62 digits. **Product survival rule (Observation):** at order p, ζ(2)×ζ(odd) always cancels algebraically; surviving products are max(0, floor((2p-5)/2)) terms, all of the form ζ(odd≥3)×ζ(even≥4) with total weight 2p−1. Verified through D₆ (p=6): ζ(2)ζ(9) absent, 3 surviving products as predicted. Mechanism: the n² Fock degeneracy in the Euler sum coefficients drives the cancellation — the same n² that produces F = ζ(2) = D_{n²}(4) in Paper 2. **K–Sommerfeld structural separation (PSLQ negative at 200 dps):** K/π is NOT in the ℚ-span of {D₂, D₃, D₄, D₅, D₆, ζ(2), ζ(3), 1}. K lives in the spectral-zeta π-polynomial ring (T9 theorem); D_p lives in the full Euler-Zagier MZV algebra. The ζ(2) protection is structural but does NOT by itself explain why K has no odd-zeta content — the question remains why the K combination rule avoids the full MZV ring that physical energy coefficients access. Key files: `debug/fine_structure_d5_computation.py`, `debug/d6_analytical_assembly.py`, `debug/k_sommerfeld_connection.py`, `debug/data/fine_structure_d5.json`, `debug/data/d6_analytical_assembly.json`, `debug/data/k_sommerfeld_connection.json`. Paper 28 §7 updated with D₆ subsection, product survival rule extended through p=6, and K–Sommerfeld structural separation through D₆.
- **QED STRATEGIC DIRECTIVE (2026-04-27): GRAPH-NATIVE FIRST.** New QED investigation work should build out the graph-native framework (GN-1 through GN-7) rather than extending continuum spectral mode sums. Continuum QED (`qed_self_energy.py`, `qed_vertex.py`, etc.) is validated and correct but computationally expensive for precision extraction — the c3 coefficient required n_int=50 with tail corrections and Richardson extrapolation for ~3 significant digits. The graph-native approach is exact at every finite n_max, lives in an algebraic ring ℚ[√2,√3,√6,...], and the projection exchange constant C is rational at every truncation. The goal is to build out graph-native QED to higher n_max, characterize the graph limit, and use it to validate/recover continuum QED results (Paper 28 theorems, c2 formula, Sommerfeld D_p sums). This is the QED analog of the core GeoVac principle: the graph is the fundamental object, continuum physics emerges as a projection/limit. Continuum QED infrastructure remains available for cross-validation.
- Graph-native QED on the finite Fock graph (v2.26.1, Tracks GN-1 through GN-7): seven modules implementing QED entirely on the finite GeoVac graph using exact rational/algebraic arithmetic. **All graph-native QED quantities are π-free** — live in ℚ[√2, √3, √6, ...] (algebraic integers from CG coefficients), not ℚ. Key results: (1) Fock graph Hodge decomposition verified (GN-1, `geovac/fock_graph_hodge.py`); (2) electron propagator G_e = D⁻¹ in the Dirac (n,κ,m_j) basis (GN-2, `geovac/graph_qed_propagator.py`); (3) photon propagator G_γ = L₁⁺ and VP as finite matrix trace (GN-3, `geovac/graph_qed_photon.py`); (4) vertex coupling tensor V[a,b,e] via CG projection (GN-4, `geovac/graph_qed_vertex.py`); (5) self-energy Σ with Tr=44/3 at n_max=2, positive semidefinite, 5 zero eigenvalues matching photon kernel (GN-5, `geovac/graph_qed_self_energy.py`); (6) vertex correction Λ with Tr=32/9, anomalous moment F₂=5√2/3 ≈ 2.357 (irrational but π-free); (7) **broken structural zero**: Paper 28 Theorem 4 (Σ(n_ext=0)=0 by vertex parity) does NOT hold on the finite graph — **pendant-edge theorem: Σ(GS) = 2(n_max-1)/n_max → 2** exactly. The ground state v₀=|1,0,0⟩ is a leaf node (pendant vertex) at all n_max; the unique edge e₀ connecting it to the graph decouples in L₁, giving G_γ[e₀,e₀] = (n_max-1)/n_max via path-graph Laplacian inverse. Two continuum protection mechanisms are both absent: (a) vertex parity (n₁+n₂+q odd) requires vector photon quantum numbers, (b) SO(4) channel count W(0,n_int,q)=0 requires the photon to carry SU(2)_L×SU(2)_R vector structure. **Graph computes scalar QED; vector structure is calibration** (see scalar_vs_vector_qed.md); (8) VP at n_max=3 decomposes by l-sector (l=0: 2×2, l=1: 7×7, l=2: 4×4) with zero cross-block coupling (GN-6, `tests/test_graph_qed_nmax3.py`); (9) **projection exchange constant is RATIONAL** at every finite truncation — C = 50471424/1779441125 at n_max=3 — confirming Paper 18 CALIBRATION tier (GN-7, `geovac/graph_qed_continuum_bridge.py`); (10) graph captures ~35× more VP content than continuum W-weighted formula at n_max=3. **Neumann/Born series for the electron propagator:** the Dirac graph operator D = Λ + tA decomposes into diagonal eigenvalue matrix Λ and sparse adjacency A; G_e(t) = Λ⁻¹ Σ_{k≥0} (-t·Λ⁻¹A)^k terminates exactly at order N-1 (Cayley-Hamilton), making F₂(t) a rational function of t over ℚ(√2,√3,√6); spectral radius ρ(Λ⁻¹A) ~ 0.2 at t=κ=-1/16 guarantees rapid convergence of the Born series; first-order perturbation theory already captures F₂ to within ~0.2%. Implemented as production path in `graph_qed_propagator.py`. Paper 28 updated with new §Graph-Native QED. Key files: `debug/data/gn5_self_energy_vertex.json`, `debug/data/gn6_nmax3_vp.json`, `debug/data/gn7_continuum_bridge.json`, memos in `debug/gn{5,6,7}_*_memo.md`. Tests: 45/45 (GN-6), 63/63 (GN-7) passing.
- Graph-native QED three-track sprint (v2.26.1, April 2026): F₂ rational structure, selection rule census, and self-energy gap analysis. **Track 1 — F₂(t) is an EVEN rational function** of the hopping coupling t, degree 2/2 over ℚ(√2,√3,√6) at n_max=2. All odd Taylor coefficients vanish identically due to t→−t symmetry of the vertex bilinear. F₂(κ) convergence: 2.353 (n=2), 1.873 (n=3), 1.589 (n=4); effective power law F₂ ~ n_max^(−0.57). n_max=2 exact symbolic: F₂(0) = 5√2/3, Taylor c₀ = 5√2/3, c₂ ≈ −1.163 (algebraic), c₁ = c₃ = c₅ = 0. n_max=3 exact symbolic (28×28 inversion in 6.4s): F₂(κ) ≈ 1.873. n_max=4 numeric only (sympy eigenvalue bug in L₁ at E=34). **Track 2 — Selection rule survival census:** Only 1 of 8 continuum QED selection rules survives on the finite graph (Gaunt/CG sparsity). 7 broken rules partition cleanly: 3 STRUCTURAL (vertex parity, SO(4) channel count, Ward identity — require vector photon) + 4 INTRINSIC (Δm_j, spatial parity, charge conjugation, Furry — graph kinematics). Notable invariants: Tr(Σ)/Tr(Λ) = 33/8; tadpole = 16√2/15; bubble = 32/9; triangle = 3584√2/3375. **Track 3 — C × F₂ NEGATIVE RESULT:** C×F₂ grows with n_max (0.053 at n=3, 0.075 at n=4), diverging from α/(2π) = 1.16×10⁻³. The VP projection constant C is NOT the correct projection for the vertex correction. Different diagram topologies require different spectral-density factors for the graph-to-continuum projection. C(4) = 48613183465472/1024864407042525. Self-energy gap: graph Σ uniformly ~10× larger than continuum (missing SO(4) W=0 suppression). Paper 28 updated with F₂ convergence table, selection rule census table, and topology-dependent projection subsection. Key files: `debug/gn_qed_sprint_memo.md`, `debug/gn_selection_rule_census_memo.md`, `debug/data/gn_selection_rule_census.json`.
- Graph-native QED F₂ convergence sprint (v2.26.1, April 2026): Extended F₂ computation from n_max=2..4 to n_max=2..6 using numpy (bypassing sympy eigenvalue bug at n_max≥4). **F₂(κ) table:** 2.353 (n=2), 1.873 (n=3), 1.589 (n=4), 1.396 (n=5), 1.253 (n=6). **Power-law convergence:** F₂ ~ 3.507 × n^(-0.573), R² = 0.99990; local exponents increasing monotonically (-0.562 → -0.591). Richardson extrapolation: F_inf ~ -0.22 (unphysical → confirms F₂ → 0 as n_max → ∞). **Rational function structure:** n_max=2 degree 2/2 over ℚ(√2,√3,√6); n_max=3 degree 16/16 over ℚ(√2,√3,√5,√6,√10,√15); Taylor coefficients at n_max=3 all negative for k≥2 (better convergence than n_max=2's alternating signs). **F₂(κ) vs F₂(0) sign flip** at n_max=5 (zero-crossing of the difference; both sequences track within <0.1%). **Pendant-edge theorem** Σ(GS) = 2(n_max-1)/n_max verified exactly at all 5 tested n_max (mechanism: GS is pendant vertex, G_γ[e₀,e₀] = (n_max-1)/n_max from path-graph Laplacian inverse). **Self-energy trace scaling:** Tr(Σ)/N_dirac grows sublinearly (1.47 → 3.23 across n_max=2..6). **Successive ratios** F₂(n+1)/F₂(n) increasing toward 1: 0.796, 0.849, 0.878, 0.898 (consistent with slow power-law). Key files: `debug/f2_convergence_sprint_memo.md`, `debug/f2_convergence_nmax56.py`, `debug/f2_rational_nmax3.py`, `debug/data/f2_convergence_nmax56.json`, `debug/data/f2_rational_nmax3.json`.
- Native Dirac graph QED (v2.26.1, April 2026): QED built directly on the (n,κ,m_j) Dirac graph, bypassing the CG projection used in the scalar Fock graph QED (GN-1 through GN-7). Two adjacency rules tested: Rule A (κ-preserving, spinor lift of scalar ladders) and Rule B (E1 dipole, Δl=±1 parity flip). **HEADLINE: Rule B Dirac graph recovers 4/8 continuum QED selection rules** (vs 1/8 on scalar Fock graph). Selection rules partition cleanly into three categories: (1) **Always survives (1):** Gaunt/CG sparsity; (2) **Spinor-recoverable (3):** Δm_j conservation (both rules), spatial parity E1 (Rule B only), Furry's theorem (both rules) — fixed by using Dirac nodes instead of scalar Fock nodes; (3) **Vector-photon-required (4):** vertex parity, SO(4) channel count, Ward identity, charge conjugation — require promoting photon from scalar 1-cochain to vector harmonic. **GS is NOT pendant on either Dirac graph** — degree 2 (Rule A) and degree 5 (Rule B), constant across n_max=2,3,4. Structural zero Σ(GS)=0 NOT recovered; Σ is strictly positive definite (no zero eigenvalues, unlike scalar graph's 5 zeros). Exact algebraic eigenvalues: Tr(Σ_A) = 17/2, Tr(Σ_B) = 103/20; GS blocks: Rule A diag(5/8,5/8), Rule B [[103/160,21/160],[21/160,103/160]]. All quantities π-free in ℚ[√2,√17,√41,√881] — algebraic content SIMPLER than scalar Fock (ℚ[√2] for Rule A vs ℚ[√2,√3,√6]). Furry's theorem recovered because identity vertex on Dirac graph is off-diagonal (odd under graph-transposition). The Dirac graph is an **intermediate construction** between scalar Fock (1/8 rules) and continuum vector QED (8/8 rules). Paper 18 classification: spinor-recoverable rules are INTRINSIC to spinor quantum numbers; vector-required rules are CALIBRATION content from the continuum embedding. Key files: `debug/dirac_graph_qed_sprint.py`, `debug/dirac_graph_qed_exact.py`, `debug/data/dirac_graph_qed_sprint.json`, `debug/data/dirac_graph_qed_exact.json`, `debug/dirac_graph_qed_memo.md`.
- GN-QED Spectral-Zeta Projection Sprint (v2.26.1, April 2026): Three-track sprint characterizing the graph-to-continuum projection for QED on S³. **Track 1 — Per-mode topological gap:** The graph-to-continuum self-energy ratio ρ(n_ext, n_int) varies from 0 to 1.66 with no clean pattern. Root cause: Fock graph has nearest-neighbor connectivity (Δn=±1 only) while continuum SO(4) channel count W obeys a triangle inequality allowing non-nearest-neighbor contributions, creating fundamental zero/nonzero asymmetries. **Track 2 — Catalan connection CLEAN NEGATIVE:** 24 PSLQ targets against {G, β(4), π², ζ(3)} all returned null. The χ₋₄ Dirichlet character mechanism that produces Catalan's G in the continuum is NOT accessible from finite-graph parity decomposition. Graph even/odd ratios alternate with n_max (driven by degeneracy counting); continuum ratios converge to ~1. Number field shift: F₂_full ∈ ℚ(√2) but F₂_even/F₂_odd ∈ ℚ(√3). **Track 3 — Two-tier calibration structure:** Three independent projection constants C_VP, C_SE, C_F2 have structurally independent power-law scaling: C_VP ~ n^0.93 (vertex-based) or n^1.48 (photon-propagator), C_SE ~ n^1.48, C_F2 ~ n^0.57. **VP/SE constant ratio: C_VP/C_SE = 0.1035 ± 0.0009 (CV=0.83%, closest rational 3/29 = 0.10345)** — VP and SE share the same power-law exponent, differing only by a constant ~1/10 prefactor. C_F2 exponent (+0.57) matches F₂ convergence exponent (−0.57) with opposite sign, so C_F2 × F₂ → α/(2π) by construction. Two-tier structure: (a) VP/SE shared tier with exponent ~1.48 and constant ratio 3/29; (b) vertex tier with exponent ~0.57. Paper 28 updated with three-diagram projection table, VP/SE ratio finding, and new open question. Key files: `debug/spectral_projection_mode_resolved.py` (Track 1), `debug/spectral_projection_parity_split.py` (Track 2), `debug/spectral_projection_track3.py` (Track 3), `debug/spectral_projection_sprint_memo.md`, `debug/data/spectral_projection_*.json`.
- Vector QED sprint on the Fock graph (VQ-1 through VQ-5, v2.26.1, May 2026): systematic five-track investigation testing whether the 4 vector-photon-required selection rules can be recovered from the graph's intrinsic geometry without promoting the photon to a vector harmonic. **All five tracks NEGATIVE for selection rule recovery.** (VQ-1) σ·L decomposition → zero on all graph edges (intra-shell only). (VQ-2) pure σ_μ vertex → disconnected graph (2/8 rules, worse than Dirac Rule B baseline). (VQ-3) σ_μ ⊗ V_scalar vertex → trivial 3× rescaling by Pauli trace identity Σ_μ σ_μ² = 3I; Tr(Σ_vec) = 3 × Tr(Σ_scalar) = 44; F₂_vec = 3 × F₂_scalar = 5√2; number field ℚ(√2) unchanged; breaks Δm_j and Furry compared to scalar; selection rules 1/8. (VQ-4) Direction-resolved Hodge decomposition: T-channel (radial, Δn=±1) and L-channel (angular, Δm=±1) edges classified with per-channel sub-Laplacians. **POSITIVE STRUCTURAL**: T/L spectra are anisotropic (T eigenvalues {2.0} vs L {1.0, 3.0} at n_max=2); propagator block-diagonal at n_max=2, cross-channel coupled at n_max=3 (||G_TL||/||G_full|| = 5.8%). Gives scalar photon effective "polarization" from pure topology. (VQ-5) Full vector QED combining direction channels with σ vertex: all 15 σ-to-channel configurations tested (5 named + exhaustive 3×3 scan). Three distinct F₂ values (1.0×, 1.4×, 1.6× scalar), but ALL produce 1/8 surviving selection rules. **Structural conclusion:** the 4 vector-required rules (vertex parity, SO(4) channel count, Ward identity, charge conjugation) require the PHOTON to carry (L, M_L) quantum numbers — neither spin matrices at the vertex, nor direction-resolved edges, nor their combination suffices. Confirms the three-tier selection rule partition (1 always survives, 3 spinor-recoverable, 4 vector-required) from a complementary direction. σ-vertex approaches are now definitively closed (5 independent negatives). Paper 28 updated with new §VQ sprint subsection. Key files: `debug/vq{1,2,3,4,5}_*.py`, `debug/data/vq{3,4,5}_*.json`, `debug/vq{1,2,3,4,5}_*_memo.md`, `debug/vq_vector_qed_sprint_plan.md`.
- QED flat-space limit verification (v2.19.3): new module `geovac/qed_flat_limit.py`. Verifies that S³ spectral sums are structurally consistent with flat-space QED in the R→∞ limit. Key results: (1) Weyl's law confirmed — Dirac mode count N = 2(N+1)(N+2)(N+3)/3 matches Weyl prediction (2/3)R³p_max³ with O(1/n_max) convergence; (2) sunset sum scales exactly as R^{2s1+2s2+2p} (power=10 at default parameters), ratio S(R)/R^power is R-independent to machine precision; (3) D(s,R)/R^s is R-independent, confirming the even/odd parity discriminant (Paper 18) persists at all radii — the transcendental class is structural, not a finite-R artifact; (4) D(5) = 14ζ(3) − 31/2·ζ(5) verified to 80 dps via Hurwitz, ζ(3) coefficient R-independent. HONEST LIMITATION: exact coefficient matching for the full Petermann a₂ = −197/144 + π²/12 + 3ζ(3)/4 − π²ln2/2 requires matching specific diagram topologies and renormalization, beyond the scope of structural spectral sums. 28 tests (25 fast, 3 slow) in `tests/test_qed_flat_limit.py`. Key files: `geovac/qed_flat_limit.py`.
- Paper 3 CFT diagnostic (v2.18.2, research finding): spectral dimension d_s≈2 is a graph topology artifact (any tridiagonal graph gives d_s→2), central charge not genuine CFT (no conformal symmetry in graph Laplacian), Berry phase data retracted (v1.2.0), hyperradius is not holographic (no area-law entropy). Paper stays conjectural with caveats documented.
- Sprint 4 (April 2026, v2.15.0) completed three parallel tracks. **(DD) Drake 1971 fine-structure J-pattern derived from Racah 6j algebra**: f_SS(J) = (-2,+1,-1/5) and f_SOO(J) = (+2,+1,-1) at L=S=1 both derived symbolically from the identity (-1)^(L+S+J)·6j{L,S,J;S,L,k}·(-6) at k=2 (SS) and k=1 (SOO). Spin-tensor choice structurally fixed: ⟨S=1‖[s₁⊗s₂]^(2)‖S=1⟩=√5/2 (non-zero, rank-2 for SS) while ⟨S=1‖[s₁⊗s₂]^(1)‖S=1⟩=0 (vanishes, forcing SOO's s₁+2s₂ form per Bethe-Salpeter §38.15). HONEST NEGATIVE: direct/exchange mixing ratios (3/50, -2/5, 3/2, -1) NOT closed from 9j alone — requires full Varshalovich §5.17 / Brink-Satchler App. 5 bipolar harmonic machinery (future sprint). 5 new tests pass, Paper 14 §V updated. **(TR) Tier-2 T3 spinor FCI α=0 regression FIXED**: missing (-1)^(j_a+1/2) reduced-matrix-element phase in jj_angular_Xk per Grant 2007 Eq. 8.9.9/8.9.11 and Johnson Eq. 3.69. Diagnostic signature: X_0 diagonal gave -1 for j=1/2 and j=5/2 (must be +1 for unit-charge density). Single-line fix in `geovac/composed_qubit_relativistic.py`. Verification: spinor FCI at α=0 now matches scalar FCI at Z=4, Be²⁺, n_max=2,3,4 to 1-ULP (gaps < 1e-14 Ha). **Downstream impact: Pauli counts at n_max=2 increase 1.76× (533→942 for CaH/SrH/BaH rel, 805→1413 for LiH/BeH rel); at n_max=3 increase 1.92× (30,940→59,484 for CaH rel, 46,434→89,226 for LiH/BeH rel); rel/scalar ratio corrected 2.42×→4.24× (n=2) and 5.89×→11.33× (n=3)**. Isostructural invariance preserved. SO shift corrected from contaminated ~10⁻⁸ Ha to physically correct ~10⁻¹²-10⁻¹¹ Ha. GeoVac still wins Sunaga head-to-head ~80-120× at matched Q=18 (revised from 150-250× projection). Paper 14 §V, Paper 20 Tier-2 table, `docs/tier2_market_test.md` all updated with post-TR values. **(QG) Paper 25 drafted** (papers/synthesis/paper_25_hopf_gauge_structure.tex, 1023 lines, 8 sections, 23 refs): framework observation that GeoVac's Hopf graph at finite n_max is simultaneously a discretization of S³ (Paper 7), a triangle-free simplicial 1-complex carrying discrete Hodge decomposition (Lim 2020, Schaub 2020), and a Wilson-type lattice gauge structure with matter on nodes and U(1) connection on edges via ladder-operator phases. Under Paper 2's conjectural interpretation, B = Casimir matter trace on S² base, F = Fock-Dirichlet gauge zeta, Δ⁻¹ = Dirac-mode boundary count. CLAUDE.md §1.5 rhetoric enforced: sharp mathematical observation separated from conjectural physical interpretation at every step; no new α derivation. Data: `debug/dd_drake_derivation.md`, `debug/tr_fix_memo.md`, `debug/qg_paper25_memo.md`, `docs/sprint4_final_summary.md`.

- **Transcorrelated GeoVac (TC He at Level 3)** — NEGATIVE for adiabatic solver (v2.0.34, Track BX-2). BX-1 scoping was positive (Gaunt selection rules preserved, BCH terminates, V_ee cancels with correct Jastrow sign). BX-2 revealed structural incompatibility with adiabatic solver: V_ee is O(R) in angular eigenvalue but TC gradient G_ang is O(1). However, **TC in second quantization (Track BX-3) is POSITIVE**: TC-modified qubit Hamiltonians via the composed pipeline show (i) He accuracy improved from 5.3-8.2% to 3.3-3.6% across max_n=1-3, eliminating basis divergence; (ii) composed Pauli ratio exactly 1.68× (constant factor, O(Q^2.5) scaling preserved); (iii) electronic-only 1-norm overhead 9-16%, vanishing with PK dominance. Implementation: `geovac/tc_integrals.py` with `compute_tc_integrals_block()` and `build_tc_composed_hamiltonian()`. BCH constant shift is -1/4 per electron pair (sign error found and fixed during BX-3). **TC angular gradient (Track BX-4, v2.0.36) is a NEGATIVE RESULT for quantum efficiency**: angular gradient for l>0 orbitals was already implemented in `tc_integrals.py` but never benchmarked. He corrected FCI (apply_op_string): radial 3.621% → full 3.611% at max_n=2 (0.01 pp, 2.66× Pauli increase 188→500), radial 3.639% → full 3.625% at max_n=3 (0.014 pp, 2.31× ERI increase). Composed molecules: 2.67× Pauli increase (LiH 562→1498, BeH₂ 936→2496, H₂O 1310→3494). Cost/benefit ratio >100×. `include_angular=False` is the default for composed pipeline. Gaunt selection rules verified (|ΔL|=1, |Δm|≤1).

- **Coupled composition (Track CB, v2.0.37)** — NEGATIVE RESULT. Scoping investigation: replace PK with explicit cross-block ERIs in composed architecture. LiH at max_n=2: 854 Pauli terms (2.56× composed 334), 85.69 Ha 1-norm (2.30× composed 37.33), 88 QWC groups (4.19× composed 21). 4-electron FCI: coupled gives 29% error (WORST) vs composed 15%, no-PK-no-cross 10.9%. Root cause: cross-block ERIs add core-valence repulsion but the composed orbital basis (different Z per block) has no cross-center nuclear attraction (h1 terms where orbitals feel the other nucleus). Adding V_ee cross-terms without V_ne cross-terms is energetically unbalanced. Cross-center h1 requires two-center integrals not available in single-center hydrogenic basis. PK is NOT redundant in second quantization for the composed basis. Key files: `geovac/coupled_composition.py`, `tests/test_coupled_composition.py`, `docs/coupled_composition_scoping.md`.

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
| TC in second quantization (composed qubit pipeline) | 1 | Track BX-3 (v2.0.35) implemented TC-modified qubit Hamiltonians. Original benchmark used qubit-space diag (full 2^Q matrix) producing non-physical energies below variational bound — false positive. **Corrected** (v2.9.0, Track TC-V): particle-number-projected FCI shows standard converges (5.3%→2.5%→2.2%) while TC plateaus at 3.3-3.4%. **Extended to n_max=5 (v2.9.2, CUSP-3): TC error plateaus monotonically 3.40→3.44→3.46→3.47%, ratios approaching 1 (0.987→0.994→0.997); Standard continues to improve 2.48→2.17→2.07→2.02. Pauli ratio TC/Std grows 1.68→2.33→2.60→2.73.** TC is decisively dead at every accessible n_max in second quantization on the composed basis. Root cause: the cusp is an energy-evaluation problem, not a wavefunction problem (confirmed by entanglement analysis, Paper 26); TC smooths the wavefunction but doesn't help the energy integral in the angular momentum eigenbasis. VarQITE converges rapidly (16 steps) — the quantum algorithm works, the Hamiltonian itself is the issue. |
| TC 2D variational cusp correction (basis doubling) | 3 | Three approaches attempted: (1) R sin(α) basis doubling — Hcc/Scc=-0.84 >> E₀=-2.90, R factor adds 2 Ha kinetic energy, variational bound respected but energy worsens. (2) sin(α)-only angular augmentation — bound OK but adds negligible new content beyond existing Gegenbauer multi-channel basis (~5% worse). (3) Commutator-based Rayleigh quotient — <ψ₀\|[H,r₁₂]\|ψ₀>=0 identically for real wavefunctions (fundamental theorem). **NOTE: TC similarity transformation (solve_he_tc_2d) succeeds** — 0.001% at l_max=4, gamma=0.11, but non-monotonic convergence (sweet spot, not systematic l⁻⁸). Optimal gamma=0.11 << Kato gamma=1 indicates finite-basis correction, not true cusp removal. |
| TC angular gradient for l>0 orbitals | 1 | Angular gradient adds 2.66× Pauli terms (He max_n=2: 188→500) for 0.01 pp accuracy improvement (3.621%→3.611%). Composed: 2.67× increase (LiH 562→1498). Cost/benefit >100×. Gaunt selection rules preserved but coupling density defeats sparsity advantage. Radial-only is optimal. |
| Coupled composition (cross-block ERIs replacing PK) | 1 | Cross-block ERIs add 2.56× Pauli terms and 2.30× 1-norm for LiH. 4-electron FCI: 29% error (worst of three configurations). Root cause: composed orbital basis (different Z per block) means cross-block ERIs add V_ee repulsion without balancing V_ne cross-center attraction. PK is NOT redundant in second quantization for the composed basis. Two-center h1 terms required but unavailable in single-center framework. |
| Single-center nested molecular LiH | 1 | R_eq 33.7% error (2.0 bohr vs 3.015). Single Z for all electrons can't represent both core and bond length scales at max_n=2. D_e correct (4.7% error) but geometry wrong. Track DF Sprint 4. |
| Two-center charge-center nested LiH | 1 | 48.2% energy error (catastrophic). Truncated orbital basis cannot represent 1s² core density off-center. Charge-center origin incompatible with nested encoding. Track DF Sprint 4B. |
| Heterogeneous nested (per-pair Z_eff) | 1 | Löwdin orthogonalization destroys Gaunt sparsity: 1,711 vs 120 Pauli terms (14× inflation). Cross-exponent orthogonalization mixes radial indices, creating dense ERIs. R_eq improved to 17.1% but sparsity advantage eliminated. Track DF Sprint 5. |
| Fock energy-shell self-consistency for He | 1 | k²=-2E over-constrains 2-electron problem. SC gives 12.8% at n_max=1, 5.1% at n_max=2 vs variational 1.9%, 1.6%. Single projection parameter insufficient for two electrons — transition from calibration to embedding exchange constants. Track DI Sprint 2. |
| SM-running origin for Δ = 1/40 (Paper 2 alpha) | 6 | Phase 4H sprint, six tracks. (SM-A) One-loop QED running gives Δ(1/α)=1/40 only at μ=574.9 MeV (1.125·m_e), with no recognizable physics correspondence; all natural geometric scales (Bohr, 2 Ry, m_μ, M_Z, Planck) miss by ≥45×. (SM-B) Σ N_c Q_f² = 8 = \|λ_3\| is a numerical coincidence — per-generation 8/3 is charge-universal while every per-shell S³ invariant varies with n. (SM-C) 1/40 does not appear naturally in SU(5)/E_8/Spin(10) as Chern number, η-invariant, or CS level; dual-Coxeter rewrite is tautological. (SM-F) Higher Hopf S⁴/S⁸: selection principle B/N=dim does not transfer; closest near-misses are Diophantine (B/N=30 vs 1/α_W=29.5 at 1.7%, B/N=8 vs 1/α_s=8.467 at 5.5%), no spectral identity. POSITIVE SIDE-RESULT (SM-D): the perturbative vacuum-pol calculation fails (overshoots 45-530×) but produces the cleaner combinatorial identity Δ⁻¹ = g_3^Dirac(S³) = 2(n+1)(n+2)\|_{n=3} = 40, the third single-chirality Dirac mode degeneracy on unit S³ — exactly matches Paper 2's selection-principle cutoff. SM-E confirms Δ is absent from HO/nuclear sector (Papers 23, 24), consistent with Coulomb-specific. Phase 4H. |
| Schwartz-tail hot-node patch on He Z=2 graph-native CI (CUSP-2) | 1 | CUSP-2 (v2.9.2) tested adding Schwartz cusp tail `-A/(l_max+2)^4` to the (1s,1s) singlet pair-state diagonal of He graph-native CI. Fit across l_max=3..6 gave err_abs = -1.48/(l+2)^4 + 5.99e-3 Ha: the residual after any Schwartz subtraction is a CONSTANT 6 mHa offset that does not decrease with l. The actual Schwartz tail at l_max=6 is <1 mHa; the 0.20% ceiling is the **small-Z graph-validity-boundary artifact** (Z_c≈1.84, He at Z=2 is just above). Z=10 confirmation: error sign flips to under-binding +85 mHa, consistent with real basis truncation. Patching with a Schwartz tail that doesn't match the data actively worsens accuracy. **Implication: for He at Z=2 graph-native CI, cusp corrections cannot improve accuracy because cusp is not the limiting error.** Schwartz extrapolation remains valid for Level 4 (Paper 15) molecular solvers where there's no small-Z graph artifact. |
| Energy graph for V_ee on S³ (Paper-12 analog search) | 1 | Pair-state graph at n_max=3 (31 nodes) and n_max=4 (101 nodes) characterized. Within-parity blocks ~47% dense — orbital Gaunt sparsity does NOT project to pair-states. Diagonal V_ii are exact rationals but cross-shell denominators (2^k vs 3^j) don't close — spectrum is non-rational. No three-term recurrence in (n, l) for Slater integrals; Neumann-style separability is specific to prolate spheroidal, not S³. POSITIVE PARTIAL: wavefunction graph diagonalizes 92-94% of V_ee Frobenius mass on its own (‖[H₁,V]‖/‖V‖ = 6.1% at n_max=3, 5.3% at n_max=4, saturating); cusp is concentrated on (1s1s) pair-state, not distributed. v2.9.1. |
| Darwin+MV for He/Li/Be 2p-doublet improvement | 1 | NEGATIVE. Both 2p states share l=1, so Darwin=0 for l>=1 and MV cancels in the splitting. Residual 66-211% errors trace to multi-electron SS/SOO. Tier 3 Track T8. |
| Balanced + frozen-core PES for second-row molecules | 2 | NaH (Q=20, 2e) and MgH₂ (Q=40, 4e) both show monotonic overattraction with no equilibrium at n_max=2. Root cause: frozen [Ne] core hides core screening from cross-center V_ne; valence sees full Z=11/12 nuclear charge from other center. LiH works because explicit core block provides screening. Balanced + frozen-core is limited to first-row for PES (single-point energy at fixed R still valid). Sprint 7a, v2.19.4. |
| Single-constant graph-to-continuum QED projection (C×F₂→α/(2π)) | 1 | NEGATIVE. C×F₂ grows with n_max (0.053 at n=3, 0.075 at n=4), diverging from α/(2π)=1.16×10⁻³ by 46-65×. The VP projection constant C grows faster than F₂ decreases. Root cause: different diagram topologies (VP vs vertex correction) require different spectral-density matching factors; no single multiplicative constant converts graph-native scalar QED to physical vector QED. The graph-to-continuum projection is topology-dependent, not a universal focal-point constant. Sprint GN-QED, v2.26.1. **Extended (Spectral-Zeta Projection Sprint, v2.26.1):** three independent projection constants C_VP, C_SE, C_F2 computed at n_max=2..5 confirm the single-constant projection is definitively dead across ALL three QED diagram types, not just VP×F₂. However, VP and SE share a common projection scaling: C_VP/C_SE = 0.1035 ± 0.0009 (CV=0.83%, closest rational 3/29) is constant to sub-percent across n_max=3,4,5. The ~10× prefactor difference is a pure topology-of-diagram effect. Only the vertex correction (F₂) requires a fundamentally different projection (exponent 0.57 vs 1.48). Two-tier calibration structure: VP/SE shared tier + vertex-specific tier. Per-mode projection ratio varies 0 to 1.66 (topological, not multiplicative). Catalan/β connection from graph parity: 24 PSLQ targets all null (clean negative). |
| σ-vertex and direction-resolved vector QED on Fock graph (VQ-1 through VQ-5) | 5 | Five tracks testing whether vector photon selection rules can be recovered from graph-intrinsic geometry. (VQ-1) σ·L zero on all edges (intra-shell only). (VQ-2) pure σ_μ disconnected graph, 2/8 rules (worse than baseline). (VQ-3) σ_μ⊗V_scalar trivial 3× rescaling by Pauli trace identity, 1/8 rules, breaks Δm_j and Furry. (VQ-4) T/L direction-resolved Hodge: POSITIVE STRUCTURAL (anisotropic spectra, cross-channel coupling at n_max=3), but insufficient for selection rules. (VQ-5) full vector QED (15 configurations): three distinct F₂ values (1.0×, 1.4×, 1.6× scalar), all produce 1/8 rules. Root cause: the 4 missing rules require the PHOTON to carry (L, M_L) quantum numbers, not just directional edge labels or spin matrices at the vertex. σ-vertex approaches definitively closed. |
| Dirac-sector lift of Paper 2 α combination rule ingredients B, F, Δ | 3 | Phase 4I Tier 1 sprint (D1-D5, 2026). Three obstructions proved closed-form: (1) B does not lift because the scalar Laplacian zero mode at (n=1,l=0) forces a (m-1) factor in B(m) that no Dirac/Weyl cumulative trace carries; single-point hit |λ_{m-1}|·g_{m-1}^Weyl = 42 at m=3 is the elementary coincidence (m-1)(m+2)=10. (2) F does not lift because g_m^Dirac = 2m²+2m mixes two homogeneity classes, giving D_dirac^Fock(4) = 2ζ(2) + 2ζ(3); Apéry Q-linear independence forbids isolating F without projecting back to the m² sub-weight (scalar case). Hurwitz form: D_dirac^CH(4) = π² - π⁴/12, inseparable. (3) Hopf-equivariant Orduz decomposition reproduces Δ⁻¹ = 40 as 20⊕20 charge-parity split (cleaner than Paper 2's 8·5) but does NOT produce B or F; does NOT factorize as 8·5 (the scalar "8" = \|λ_3\| of -Δ_LB is not a Dirac quantity). Positive byproduct: ζ(3) is the first non-π transcendental in the framework, appearing natively in the Dirac sector at s=d_max via the weight-m subchannel. K = π(B + F - Δ) is now a formally-documented cross-sector coincidence (scalar + spinor) with three spectral homes. NO Tier 1b opened per sprint plan D5 decision table; move to Tier 2 spin-ful composed qubit encoding. |


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

*Practical test:* If a proposed modification changes which quantum numbers label the states or which transitions are allowed, it violates the prime directive. If it changes how accurately the radial amplitude is computed within a given channel, it is a legitimate numerical improvement.

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

### Context Loading Guide

Papers are grouped by how frequently they need to be in a sub-agent's context. **Always-load** papers define the framework's identity — nodes are quantum numbers, the graph is the physics. **Load-on-topic** papers are needed when working in their specific area. **Archive** papers are historical or scaffolding; load only on explicit request.

| Paper | Priority | Load when... |
|:------|:---------|:-------------|
| Paper 0 | **Always** | Packing construction, K = -1/16, nodes = quantum numbers — the axiom set |
| Paper 1 | **Always** | Spectral graph theory, O(N) eigenvalue methods — the computational identity |
| Paper 7 | **Always** | S³ proof, conformal equivalence, dimensionless vacuum — the theoretical foundation |
| Paper 14 | **Always** | Qubit encoding, Pauli scaling, composed sparsity — the headline result |
| Paper 22 | **Always** | Angular sparsity theorem, potential-independent ERI density, universal/Coulomb partition |
| Paper 24 | **Always** | Bargmann-Segal lattice for the 3D HO, π-free graph, HO rigidity theorem, Coulomb/HO asymmetry |
| Paper 2 | On-topic | α conjecture: K = π(B + F − Δ) three-sector spectral coincidence on the Fock-projected S³; Marcolli-vS gauge-network lineage (WH1); Sprint A (2026-04-18) structural verifications. Promoted to core 2026-04-18; combination rule itself remains conjectural per §13.5. |
| Paper 18 | On-topic | Exchange constants, transcendental taxonomy, observable classification, Weyl-Selberg |
| Paper 15 | On-topic | Level 4 solver, H₂ accuracy, channel convergence, cusp correction |
| Paper 17 | On-topic | Composed geometry, molecules, polyatomics, PK |
| Paper 13 | On-topic | Hyperspherical, He, fiber bundle, coupled-channel |
| Paper 11 | On-topic | Prolate spheroidal, H₂⁺, spectral Laguerre |
| Paper 12 | On-topic | Algebraic V_ee, Neumann expansion |
| Paper 16 | **Always** | Chemical periodicity, S_N representation theory, atomic classifier — the generalization specification |
| Paper 6 | On-topic | Quantum dynamics, spectroscopy, AIMD |
| FCI-A | On-topic | Multi-electron atom benchmarks (He, Li, Be) |
| SCOPE_BOUNDARY.md | On-topic | Atomic/molecular scope boundary, first-row limitations, second-row feasibility, transition metals |
| Papers 8-9 | On-topic | Bond sphere geometry, Sturmian structural theorem (H ∝ S, eigenvalues R-independent), D-matrix selection rules, harmonic phase locking negative result. **GUARDRAIL:** load before ANY single-center, unified-basis, or shared-exponent molecular investigation. Predicts failure of shared-p₀ approaches (proven theorem, not empirical). Track DF Sprints 4-5 re-derived this at significant cost. |
| Paper 10 | Archive | Nuclear lattice draft; rovibrational spectra |
| FCI-M | On-topic | LCAO molecular FCI, R-independent graph Laplacian kinetic energy, 29-version diagnostic arc (v0.9.8-0.9.37). **GUARDRAIL:** load before ANY graph-concatenation or multi-center graph molecular investigation. Documents why natural geometry was necessary. |
| Paper 19 | On-topic | Coupled composition, PK alternatives, two-center integrals, Shibuya-Wulfman connection |
| Paper 20 | On-topic | Applications paper: resource benchmarks, FCI accuracy, Gaussian comparison, installation |
| Paper 23 | On-topic | Nuclear shell model qubit Hamiltonians, deuteron and He-4, Fock rigidity theorem, composed nuclear-electronic |
| Paper 21 | On-topic | Synthesis: S³ equivalence, N-electron generalization, exchange constants, Hopf α, atomic spectra program |
| Paper 25 | On-topic | Synthesis observation: the GeoVac Hopf graph as a discrete lattice-gauge structure (node Laplacian = matter propagator, edge Laplacian L₁ = B^T B = discrete Hodge-1, S² quotient = Hopf base). Places Paper 2's (B, F, Δ) inside the graph-Hodge / Wilson vocabulary. Paper 2 stays conjectural — Paper 25 is framing, not a new α derivation. |
| Paper 26 | On-topic | Entanglement structure: energy-entanglement decoupling, basis-intrinsic sparsity, hub migration, core-valence decoupling, S ~ Z^{-2.56} |
| Paper 27 | On-topic | Entropy as projection artifact: one-body entanglement-inert, HO zero-entropy rigidity, universal S_B scaling |
| Paper 28 | On-topic | QED on S³ (v2.1): 5 theorems (T9, parity discriminant, χ_{-4} identity, self-energy structural zero, ζ(3) complementarity) + depth-k tower proposition + 2 observations (Euler sum cancellation, product survival rule); D₅/D₆ Sommerfeld sums; K-Sommerfeld structural separation; S_min at 200 digits irreducible; three-loop O(N³) factorization |
| Paper 29 | On-topic | Ihara zeta, Ramanujan property, graph-RH for the GeoVac Hopf graphs, per-ℓ-shell decomposition, 12+22 factorization, Wick-rotation scope boundary on Selberg-on-hydrogen |
| Paper 30 | On-topic | SU(2) Wilson lattice gauge on S^3=SU(2); maximal-torus reduction to Paper 25 U(1); L_1=B^TB as weak-coupling kinetic term; five-action-structure synthesis |

### Paper Inventory

#### Core (`papers/core/`) — Always load

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 0 | `Paper_0_Geometric_Packing.tex` | Universal constant K = -1/16 |
| Paper 1 | `paper_1_spectrum.tex` | Spectral graph theory, O(N) eigenvalue methods. Berry phase retracted v1.2.0 |
| Paper 7 | `Paper_7_Dimensionless_Vacuum.tex` | S3 proof (18/18 symbolic), Schrodinger recovery, SO(3N) generalization |
| Paper 14 | `paper_14_qubit_encoding.tex` | Structurally sparse qubits: O(Q^3.15) atoms, O(Q^2.5) composed; Trenev et al. Gaussian baselines |
| Paper 16 | `paper_16_periodicity.tex` | Chemical periodicity from S_N representation theory; atomic classifier for composed geometry |
| Paper 22 | `paper_22_angular_sparsity.tex` | Potential-independent angular sparsity theorem; ERI density 1.44% at l_max=3; universal/Coulomb-specific partition |
| Paper 24 | `paper_24_bargmann_segal.tex` | Bargmann-Segal lattice: π-free discretization of 3D HO on S^5 Hardy space; Coulomb/HO structural asymmetry; HO rigidity theorem |
| Paper 2 | `paper_2_alpha.tex` | Fine structure constant from Hopf bundle (8.8×10⁻⁸, p = 5.2×10⁻⁹, zero free parameters): K = π(B + F − Δ) as three-sector spectral coincidence on the Fock-projected S³ with three canonical spectral homes (B = finite Casimir truncation at m=3; F = D_{n²}(d_max) = ζ_R(2); Δ⁻¹ = g_3^Dirac = 40). Sprint A (2026-04-18, v2.26.1) structural verifications: triple m=3 selection, Hopf-measure π identification (Vol(S²)/4 = Vol(S³)/Vol(S¹)), APS-shape on Δ minus sign, post-cubic residual ≈ π³α³ to 0.25% (structural hint, not derivation). Marcolli-vS gauge-network lineage (WH1). Promoted from Conjectures to Core 2026-04-18; **combination rule K itself remains conjectural** per §13.5. |
| Paper 27 | `paper_27_entropy_projection.tex` | Entropy as projection artifact: one-body operators are entanglement-inert on sparse lattices; S_kin/S_full ~ 1e-14; HO zero-entropy rigidity; universal S_B(w̃_B/δ_B) scaling with γ_∞ ≈ 1.96 |

#### Supporting (`papers/core/` and `papers/methods/`) — Load on-topic

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 6 | `Paper_6_Quantum_Dynamics.tex` | O(V) dynamics: Rabi, spectroscopy, AIMD |
| Paper 11 | `paper_11_prolate_spheroidal.tex` | Prolate spheroidal lattice: H2+ 0.70% E_min |
| Paper 12 | `paper_12_algebraic_vee.tex` | Neumann V_ee: H2 92.4% D_e, cusp diagnosis (7.6% gap) |
| Paper 13 | `paper_13_hyperspherical.tex` | Hyperspherical lattice: He 0.05%, fiber bundle, ab initio spectroscopy |
| Paper 15 | `paper_15_level4_geometry.tex` | Level 4: H2 96.0% D_e, HeH+ 93.1% D_e |
| Paper 17 | `paper_17_composed_geometries.tex` | Composed geometry: LiH R_eq 6.4%, BeH+ bound, ab initio PK |
| Paper 18 | `paper_18_exchange_constants.tex` | Spectral-geometric exchange constants: Weyl-Selberg identification of κ, e^a E₁(a), μ(R); observable classification (Sec VI); α connection |
| FCI-A | `paper_fci_atoms.tex` | He 0.35%, Li 1.10%, Be 0.90% |
| Paper 19 | `paper_19_coupled_composition.tex` | Balanced coupled: 0.20% energy, 3-molecule census, PK-free, regime-dependent 1-norm advantage |

#### Applications (`papers/applications/`) — Load on-topic

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 20 | `paper_20_resource_benchmarks.tex` | Resource benchmarks for quantum computing community: FCI PES, Gaussian comparison, 51-1,712× Pauli advantage, pip install |
| Paper 23 | `paper_23_nuclear_shell.tex` | Nuclear shell model qubit Hamiltonians: deuteron (Minnesota, 16 qubits, 592 Pauli), He-4 (712 Pauli, 1.20× for 12.25× Hilbert), composed nuclear-electronic PoC, Fock projection rigidity theorem |

#### Methods (`papers/methods/`) — Load on-topic (GUARDRAIL)

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Papers 8-9 | `Paper_8_Bond_Sphere_Sturmian.tex` | Bond sphere (positive), Sturmian structural theorem (H ∝ S, eigenvalues R-independent — **GUARDRAIL** for single-center molecular), SO(4) selection rules, multi-electron Sturmian CI negative result (v2.0.33) |
| FCI-M | `paper_fci_molecules.tex` | LCAO molecular FCI, R-independent graph Laplacian kinetic energy, 29-version diagnostic arc — **GUARDRAIL** for graph-concatenation molecular |

#### Archive (`papers/archive/`) — Load only on explicit request

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 10 | `paper_10_nuclear_lattice.tex` | Rovibrational spectra from SU(2) algebraic chains |

#### Synthesis (`papers/synthesis/`) — Load on-topic

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 21 | `paper_21_geometric_vacuum_synthesis.tex` | Synthesis: S³ proof chain, N-electron S^(3N-1), exchange constant taxonomy, 6j recoupling sparsity, Hopf α conjecture, atomic spectra research program |
| Paper 25 | `paper_25_hopf_gauge_structure.tex` | Synthesis observation: GeoVac Hopf graph as a discrete lattice-gauge structure. Node Laplacian = matter propagator; edge Laplacian L₁ = B^T B = discrete Hodge-1 = gauge propagator (Paper 2 reading). Hopf quotient S³→S² yields base graph of Phase 4B α-D. Places Paper 2's (B, F, Δ) as Casimir-weighted S² matter trace / infinite Fock-Dirichlet gauge zeta / Dirac-boundary mode count. No new α claim; Paper 2 unchanged. |

#### Observations (`papers/observations/`)

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 26 | `paper_26_entanglement.tex` | Entanglement structure of the angular momentum eigenbasis: energy-entanglement decoupling (V_ee generates 100% of entropy), basis-intrinsic sparsity (step-function ERI density transition at identity), hub migration pattern, core-valence decoupling, S ~ Z^{-2.56} |
| Paper 28 | `paper_28_qed_s3.tex` | QED on S³ (v2.1): 5 theorems + 1 proposition + 2 observations; self-energy structural zero; S_min at 200 digits irreducible; three-loop O(N³) factorization; depth-k tower; c₂ = (2-BΔ-FΔ-F/B)/5 Paper 2 bridge; D₅/D₆ Sommerfeld sums and product survival rule; K-Sommerfeld structural separation (K not in Q-span of D_p); comprehensive tables |
| Paper 29 | `paper_29_ramanujan_hopf.tex` | GeoVac Hopf graphs are Ramanujan: Ihara zeta, graph-RH, scope boundary on Selberg-on-hydrogen |
| Paper 30 | `paper_30_su2_wilson.tex` | SU(2) Wilson lattice gauge on the Hopf graph — non-abelian sibling of Paper 25; maximal-torus reduction; L_1 as kinetic term; framework-wide least-action synthesis |

#### Conjectures (`papers/conjectures/`)

| Paper | File | Key Topic |
|:------|:-----|:----------|
| Paper 3 | `paper_3_holography.tex` | Holographic entropy, spectral dimension, central charge |
| Paper 4 | `Paper_4_Universality.tex` | Mass-independence, universality, muonic hydrogen |
| Paper 5 | `Paper_5_Geometric_Vacuum.tex` | Comprehensive geometric vacuum framework (synthesis) |


---

## 7. Code Architecture

### Key Entry Points

| Task | Module | Entry Point |
|:-----|:-------|:------------|
| Atomic lattice | `geovac/lattice.py` | `GeometricLattice(Z, max_n)` |
| Atomic Hamiltonian | `geovac/hamiltonian.py` | `GraphHamiltonian(lattice)` |
| Multi-electron FCI | `geovac/lattice_index.py` | `LatticeIndex(Z, n_electrons, max_n)` |
| Direct CI (large systems) | `geovac/direct_ci.py` | `DirectCISolver(lattice_index)` |
| Molecular FCI (LCAO) | `geovac/lattice_index.py` | `MolecularLatticeIndex(atoms, R)` |
| Prolate spheroidal (H2+) | `geovac/prolate_spheroidal.py` | Prolate spheroidal lattice for diatomics |
| Hyperspherical (He) | `geovac/hyperspherical.py` | Two-electron hyperspherical solver |
| Level 4 (H2) | `geovac/level4_multichannel.py` | Molecule-frame hyperspherical |
| Core screening | `geovac/core_screening.py` | `CoreScreening(Z).solve()` |
| Ab initio PK | `geovac/ab_initio_pk.py` | `AbInitioPK(core, n_core)` |
| Composed diatomic | `geovac/composed_diatomic.py` | `ComposedDiatomicSolver.LiH_ab_initio(l_max)` |
| Rho-collapse cache | `geovac/rho_collapse_cache.py` | `AngularCache`, `FastAdiabaticPES` |
| Quantum dynamics | `geovac/dynamics.py` | O(V) time evolution |
| Hopf bundle | `geovac/hopf_bundle.py` | S3 embedding, Hopf projection, fiber analysis |
| VQE benchmark | `geovac/vqe_benchmark.py` | `build_geovac_he()`, `run_vqe()`, `collect_static_metrics()` |
| Gaussian reference | `geovac/gaussian_reference.py` | `h2_sto3g()`, `he_sto3g()`, `he_cc_pvdz()`, `he_cc_pvtz()` |
| Measurement grouping | `geovac/measurement_grouping.py` | `qwc_groups()`, `analyze_measurement_cost()` |
| Trotter bounds | `geovac/trotter_bounds.py` | `pauli_1norm()`, `trotter_steps_required()` |
| Composed qubit | `geovac/composed_qubit.py` | `build_composed_lih()`, `build_composed_beh2()`, `build_composed_h2o()`, `build_h2_bond_pair()` |
| Algebraic angular solver | `geovac/algebraic_angular.py` | `AlgebraicAngularSolver(Z, l_max)` |
| Algebraic coupled-channel | `geovac/algebraic_coupled_channel.py` | `solve_hyperspherical_algebraic_coupled()` |
| Hyperradial solver | `geovac/hyperspherical_radial.py` | `solve_radial()`, `solve_radial_spectral()` |
| Spectral angular (Level 4) | `geovac/level4_spectral_angular.py` | `SpectralAngularSolver`, `angular_method='spectral'` |
| Perturbation series (Level 3) | `geovac/algebraic_angular.py` | `perturbation_series_mu()`, `pade_approximant()` |
| Algebraic curve (Level 3) | `geovac/algebraic_curve.py` | `characteristic_polynomial()`, `algebraic_veff_matrix_lmax0()` |
| Algebraic curve (Level 4) | `geovac/algebraic_curve_level4.py` | `extract_level4_matrices()`, `probe_rho_dependence()` |
| Algebraic Z_eff | `geovac/algebraic_zeff.py` | `LaguerreZeffExpansion`, `fit_density_laguerre()` |
| Algebraic Slater integrals | `geovac/algebraic_slater.py` | `slater_fk_integral_algebraic()`, `slater_f0_algebraic()` |
| Cusp factor (Level 4) | `geovac/cusp_factor.py` | `CuspFactorSolver` (negative result — alpha-only factor insufficient) |
| Cusp graph analysis | `geovac/cusp_graph.py` | `analyze_cusp_graph()` (theory: 1/r₁₂ cannot be absorbed into S⁵ graph) |
| Cusp correction | `geovac/cusp_correction.py` | `he_cusp_correction()`, `h2_cusp_correction()`, `cusp_pes_correction()` |
| Cusp angular basis | `geovac/cusp_angular_basis.py` | `verify_so6_casimir()`, `extract_channel_coefficients()` (Track Y negative result) |
| Geometric elevation | `geovac/geometric_elevation.py` | Blow-up, Lie algebra, S³×S³ analysis (Track Z negative result) |
| PK R_ref derivation | `geovac/pk_rref.py` | `compute_rref_candidates()`, `select_rref()` (Track AD negative result) |
| N-electron solver | `geovac/n_electron_solver.py` | `scan_pes_4e_lih_multichannel()` |
| N-electron scoping | `geovac/n_electron_scope.py` | `four_electron_channel_count_molecular()` |
| N-electron spectral | `geovac/n_electron_spectral.py` | Spectral compression analysis |
| N-electron 2D solver | `geovac/n_electron_2d.py` | 2D variational solver for 4-electron system (Track AR) |
| N-electron qubit encoding | `geovac/n_electron_qubit.py` | Full N-electron quantum encoding comparison (Track AS) |
| Ecosystem export | `geovac/ecosystem_export.py` | `GeoVacHamiltonian`, `hamiltonian()`, `.to_openfermion()`, `.to_qiskit()`, `.to_pennylane()` |
| PK partitioning | `geovac/pk_partitioning.py` | `pk_classical_energy()`, `reconstruct_1rdm_from_statevector()`, `validate_pk_partitioning()` |
| Atomic classifier | `geovac/atomic_classifier.py` | `classify_atom(Z)`, `AtomClassification`, `pk_params_z2_scaled()` |
| Molecular spec | `geovac/molecular_spec.py` | `MolecularSpec`, `OrbitalBlock` |
| General composed builder | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec)`, `lih_spec()`, `beh2_spec()`, `h2o_spec()`, `hf_spec()`, `nh3_spec()`, `ch4_spec()` |
| TC integrals | `geovac/tc_integrals.py` | `compute_tc_integrals_block()`, `build_tc_composed_hamiltonian()` |
| Coupled composition | `geovac/coupled_composition.py` | `build_coupled_hamiltonian()`, `coupled_fci_energy()`, `run_coupled_scoping()` (Track CB, negative result) |
| Shibuya-Wulfman V_ne | `geovac/shibuya_wulfman.py` | `compute_cross_center_vne()`, `compute_cross_center_vne_element()` |
| Balanced coupled builder | `geovac/balanced_coupled.py` | `build_balanced_hamiltonian(spec, nuclei)` |
| Frozen core (Ne-like) | `geovac/neon_core.py` | `FrozenCore(Z)` |
| NaH spec | `geovac/composed_qubit.py` | `nah_spec()` |
| MgH₂ spec | `geovac/composed_qubit.py` | `mgh2_spec()` |
| HCl spec | `geovac/composed_qubit.py` | `hcl_spec()` |
| H₂S spec | `geovac/composed_qubit.py` | `h2s_spec()` |
| PH₃ spec | `geovac/composed_qubit.py` | `ph3_spec()` |
| SiH₄ spec | `geovac/composed_qubit.py` | `sih4_spec()` |
| KH spec | `geovac/composed_qubit.py` | `kh_spec()` |
| CaH₂ spec | `geovac/composed_qubit.py` | `cah2_spec()` |
| GeH₄ spec | `geovac/composed_qubit.py` | `geh4_spec()` |
| AsH₃ spec | `geovac/composed_qubit.py` | `ash3_spec()` |
| H₂Se spec | `geovac/composed_qubit.py` | `h2se_spec()` |
| HBr spec | `geovac/composed_qubit.py` | `hbr_spec()` |
| SrH spec ([Kr] frozen core, Sprint 3) | `geovac/molecular_spec.py` | `srh_spec()`, `srh_spec_relativistic()` |
| BaH spec ([Xe] frozen core, Sprint 3) | `geovac/molecular_spec.py` | `bah_spec()`, `bah_spec_relativistic()` |
| Alkaline-earth monohydride builder | `geovac/molecular_spec.py` | `_alkaline_earth_monohydride_spec(Z, name, ...)` |
| LiF spec (multi-center) | `geovac/composed_qubit.py` | `lif_spec()` |
| CO spec (multi-center) | `geovac/composed_qubit.py` | `co_spec()` |
| N₂ spec (multi-center) | `geovac/composed_qubit.py` | `n2_spec()` |
| F₂ spec (multi-center) | `geovac/composed_qubit.py` | `f2_spec()` |
| NaCl spec (multi-center) | `geovac/composed_qubit.py` | `nacl_spec()` |
| CH₂O spec (multi-center) | `geovac/composed_qubit.py` | `ch2o_spec()` |
| C₂H₂ spec (multi-center) | `geovac/composed_qubit.py` | `c2h2_spec()` |
| C₂H₆ spec (multi-center) | `geovac/composed_qubit.py` | `c2h6_spec()` |
| Multi-center nuclear repulsion | `geovac/composed_qubit.py` | `_compute_multi_center_nuclear_repulsion()` |
| TM hydride spec (Z=21-30) | `geovac/molecular_spec.py` | `transition_metal_hydride_spec(Z)`, `sch_spec()`, ..., `znh_spec()` |
| General composed builder | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec)` |
| General-l Wigner D rotation | `geovac/shibuya_wulfman.py` | `_build_rotation_matrix_real_sh(l, theta, phi)`, `_build_block_rotation_matrix()` |
| Nested hyperspherical | `geovac/nested_hyperspherical.py` | `build_nested_hamiltonian()`, H-set coupled basis |
| Level 3 variational 2D | `geovac/level3_variational.py` | `solve_he_variational_2d()`, `solve_he_precision()` |
| Algebraic Casimir CI | `geovac/casimir_ci.py` | `build_fci_matrix()`, `solve_self_consistent()`, `solve_variational()`, `convergence_table()` |
| Graph-native CI | `geovac/casimir_ci.py` | `build_graph_native_fci()`, `build_graph_consistent_fci()`, `_build_graph_h1()` |
| Hypergeometric Slater integrals | `geovac/hypergeometric_slater.py` | `compute_rk_float()`, `compute_rk_exact()`, exact Fraction-arithmetic R^k evaluator for arbitrary n_max |
| DUCC downfolding | `geovac/downfolding.py` | Exact (2J-K) core potential, PK divergence root-cause analysis |
| Physical constants | `geovac/constants.py` | `HBAR`, `C`, `ALPHA`, etc. |
| Breit retarded integrals | `geovac/breit_integrals.py` | `compute_radial()`, exact Fraction/sympy r_<^k / r_>^{k+3} kernel |
| Breit composed pipeline | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec, include_breit=True)` passthrough to relativistic builder |
| QED vacuum polarization | `geovac/qed_vacuum_polarization.py` | `seeley_dewitt_coefficients_s3()`, `vacuum_polarization_coefficient()`, `beta_function_qed()`, `spectral_zeta_massive()`, `spectral_zeta_derivative_at_zero()`, `classify_transcendental_content()` |
| QED vertex coupling | `geovac/qed_vertex.py` | `two_loop_odd_even_split()`, `two_loop_sunset_vertex_restricted()`, `two_loop_transcendental_classification()`, `vertex_allowed_triples()`, `reduced_coupling_squared()` |
| QED flat-space limit | `geovac/qed_flat_limit.py` | `weyl_density_check()`, `two_loop_sunset_R_dependent()`, `flat_space_limit_scaling()`, `two_loop_vacuum_energy_R_dependent()`, `extract_zeta3_coefficient()`, `verify_weyl_zeta3()` |
| QED self-energy | `geovac/qed_self_energy.py` | `self_energy_spectral()`, `self_energy_table()`, `self_energy_convergence()`, `self_energy_transcendental_class()`, `vertex_correction_spectral()`, `schwinger_convergence()` |
| QED three-loop | `geovac/qed_three_loop.py` | `three_loop_sunset_s3()`, `three_loop_unrestricted()`, `three_loop_vertex_restricted()`, `three_loop_cg_weighted()`, `three_loop_factorized()`, `three_loop_factorized_convergence()`, `decompose_three_loop_mzv()`, `three_loop_euler_maclaurin_tail()` |
| Graph-native Fock graph/Hodge | `geovac/fock_graph_hodge.py` | `build_fock_graph()`, `compute_hodge_decomposition()`, `verify_hodge_identity()` |
| Graph-native electron propagator | `geovac/graph_qed_propagator.py` | `DiracGraphOperator`, `electron_propagator()`, `neumann_electron_propagator()` |
| Graph-native photon propagator | `geovac/graph_qed_photon.py` | `build_fock_graph()`, `compute_photon_propagator()`, `compute_vacuum_polarization()` |
| Graph-native vertex tensor | `geovac/graph_qed_vertex.py` | `build_projection_matrix()`, `build_vertex_tensor()`, `vertex_tensor_to_matrices()` |
| Graph-native self-energy/vertex | `geovac/graph_qed_self_energy.py` | `compute_self_energy()`, `compute_vertex_correction()`, `extract_anomalous_moment()`, `self_energy_structural_zero_check()`, `self_energy_pi_free_certificate()` |
| Graph-native continuum bridge | `geovac/graph_qed_continuum_bridge.py` | `compute_projection_exchange_constant()`, `compute_continuum_bridge()`, `bridge_self_energy()` |

### Solver Methods

| Method | Access | Use Case |
|:-------|:-------|:---------|
| Mean-field | `LatticeIndex(method='mean_field')` | Quick atomic energies |
| Full CI (matrix) | `LatticeIndex(method='full_ci')` | Small systems (N_SD < 5000) |
| Full CI (direct) | `DirectCISolver` or `fci_method='direct'` | Large systems (N_SD >= 5000) |
| Auto | `fci_method='auto'` | Switches at N_SD = 5000 |
| Frozen core | `FrozenCoreLatticeIndex` | Active-space CI for core-valence |
| Locked shell | `LockedShellMolecule` | Extreme SD reduction for heavy atoms |

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

### Benchmarking Rule

After any modification to `hamiltonian.py`, `lattice.py`, or `solver.py`:
1. Run topological integrity: `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py -v`
2. Run validation: `pytest tests/advanced_benchmarks.py`
3. Verify 18/18 symbolic proofs pass (topological foundation)
4. Verify H2+ < 0.1% error (topological control)
5. Verify H2 Full CI < 1.0% error (accuracy control)
6. Report any speed regression > 10%

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

| Test | Max Error | Purpose |
|:-----|:---------:|:--------|
| Symbolic proofs (18 tests) | 0 failures | Topological foundation |
| H (hydrogen) | < 0.1% | Basic validation |
| He+ (helium ion) | < 0.1% | Z-scaling check |
| H2+ (ionized H2, FD) | < 0.1% | Topological control |
| H2+ (prolate spheroidal, spectral) | < 0.001% | Spectral Laguerre accuracy control |
| He (hyperspherical) | < 0.1% | Multi-electron control |
| H2 Full CI | < 1.0% | Accuracy control |
| H2 Neumann V_ee | 92.4% D_e | Algebraic integral accuracy |
| H2 Level 4 (2D solver) | 96.0% D_e | Molecule-frame hyperspherical |
| HeH+ Level 4 | 93.1% D_e | Heteronuclear extension |
| LiH Composed (ab initio PK) | R_eq 6.4% | Composed geometry |
| BeH+ Composed | Bound, physical | Transferability |
| Hyperspherical (20 tests) | 0 failures | Angular + adiabatic + radial |
| Muonic H energy ratio | < 0.01% | Mass-independence |
| V_ee S3 overlap (1s-1s, 1s-2s, 2s-2s) | < 0.01% | Topological integrity |
| Direct CI vs matrix CI | < 1e-8 Ha | Algorithmic consistency |
| JW eigenvalue consistency | < 1e-4 Ha | Qubit encoding correctness |
| H2 STO-3G Pauli terms | exactly 15 | Gaussian reference validation |
| He cc-pVDZ Pauli terms | exactly 156 | Computed Gaussian baseline |
| He cc-pVTZ Pauli terms | exactly 21,607 | Computed Gaussian baseline |
| He cc-pVDZ FCI energy | < 0.001 Ha vs published | Integral engine validation |
| He cc-pVTZ FCI energy | < 0.0002 Ha vs published | Integral engine validation |
| QWC grouping correctness | 0 violations | Measurement group integrity |
| LiH composed Pauli terms (Q=30) | exactly 334 | Composed qubit validation |
| BeH2 composed Pauli terms (Q=50) | exactly 556 | Composed qubit validation |
| H2O composed Pauli terms (Q=70) | exactly 778 | Composed qubit validation |
| H2 bond-pair Pauli terms (Q=10) | exactly 112 | Bond-pair qubit validation |
| H2 bond-pair R-independence | 112 at all R | Selection rule sparsity |
| Composed cross-block ERIs | exactly 0 | Block-diagonal integrity |
| HF composed Pauli terms (Q=60) | exactly 667 | New molecule validation |
| NH₃ composed Pauli terms (Q=80) | exactly 889 | New molecule validation |
| CH₄ composed Pauli terms (Q=90) | exactly 1000 | New molecule validation |
| General builder exact reproduction | < 1e-14 | Refactored builder matches hardcoded |
| Atomic classifier Z=1-10 | 101 tests pass | First-row classification |
| TC radial-only LiH Pauli (Q=30) | exactly 562 | TC composed radial validation (BX-3) |
| TC angular Gaunt selection rules | |ΔL|=1, |Δm|≤1 | Angular gradient preserves Gaunt rules |
| TC angular He max_n=1 identity | rad = full | No l>0 orbitals → identical |
| Algebraic angular Casimir (R=0) | < 10⁻⁶ | SO(6) eigenvalue correctness |
| Algebraic a₁ perturbation | < 10⁻⁸ | First-order perturbation theory |
| Algebraic GL vs quad consistency | < 10⁻¹⁰ | Quadrature correctness |
| Algebraic l_max monotonic convergence | monotonic | Centrifugal singularity elimination |
| Spectral vs FD consistency | < 0.01 Ha | Spectral and FD agree within FD error |
| Spectral dimension reduction | > 100× | 20 basis functions vs 5000 FD grid |
| Spectral speedup | > 10× | Wall time reduction control |
| l_max=5 coupled-channel floor | 0.15%–0.25% | Adiabatic ceiling characterization |
| n_channels=5 vs 3 consistency | < 1 mHa | Channel truncation validation |
| Spectral PES (H₂⁺) | R_eq < 0.5%, E_min < 0.001% | Spectral PES scan accuracy |
| Algebraic vs quadrature matrices | < 1e-12 elementwise | Laguerre recurrence correctness |
| Algebraic PES energy | < 1e-14 Ha vs quadrature | End-to-end algebraic validation |
| Level 3 spectral-FD consistency | < 0.0001 Ha | Hyperradial spectral solver |
| Level 3 spectral coupled-channel | 0.15%–0.25% at l_max=3 | Coupled-channel ceiling preserved |
| Level 3 spectral dimension reduction | > 100× | 25 basis functions vs 3000 FD grid |
| Level 3 spectral coupled speedup | > 50× | Coupled-channel wall time |
| Perturbation a₁ vs Paper 13 | < 10⁻⁹ | RS series first-order validation |
| Perturbation series convergence | converges R < 2 bohr | Convergence radius characterization |
| 2D solver composed LiH | drift < +0.15 bohr/l_max | 2D vs adiabatic drift comparison |
| Level 3 algebraic vs quadrature S, K | < 1e-10 elementwise | Laguerre recurrence correctness (Track H) |
| Level 3 algebraic energy consistency | < 1e-10 Ha | Algebraic S+K gives identical physics |
| Level 3 algebraic ceiling unchanged | 0.15%–0.25% at l_max=3 | Coupled-channel ceiling preserved with algebraic |
| Level 4 spectral vs FD consistency | < 0.001 Ha | Level 4 spectral hyperradial solver |
| Level 4 spectral D_e% match | within 0.5% | D_e% preserved under spectral substitution |
| Level 4 spectral convergence plateau | n_basis 20-30 within 0.0001 Ha | Convergence plateau verified |
| Level 4 spectral angular vs FD U_min | < 2e-5 Ha | Angular spectral solver accuracy (Track K) |
| Level 4 spectral angular speedup | > 200× | Angular sweep speedup control (Track K) |
| Level 4 spectral angular dimension | 20× reduction | 1000 → 50 matrix dimension (Track K) |
| Cusp factor baseline reproduction | < 1e-14 eigenvalue diff | f=1 matches standard solver (Track U) |
| Cusp factor D_e degradation | D_e decreases with γ | Negative result verified (Track U) |
| S⁵ Green singularity order | 1/d³ (not 1/d¹) | Dimensionality mismatch proof (Track W) |
| He cusp correction sign | ΔE < 0 | Cusp lowers energy (Track X) |
| He cusp correction l_max=2 | error < 0.15% | Breaks through 0.19-0.20% floor (Track X) |
| H₂ cusp correction R-dependent | varies with R | Required for D_e correction (Track X) |
| Cusp correction convergence | → 0 as l_max → ∞ | Basis-dependent, not systematic error (Track X) |
| 4-electron LiH equilibrium at l_max≥2 | structural | PK-free equilibrium validation (Track AJ) |
| S₄ [2,2] channel reduction | verified | N-electron antisymmetry machinery (Track AJ) |
| N-electron spectral compression | ≥ 100× | FD-to-spectral ratio (Track AK) |
| N-electron 2D variational bound | E_min > exact | Variational principle respected (Track AR) |
| N-electron 2D D_e sign | D_e < 0 (unbound) | Adiabatic overcounting confirmed as artifact (Track AR) |
| Full vs composed Pauli terms | full > composed | Full encoding categorically denser (Track AS) |
| Balanced coupled LiH Pauli terms (Q=30) | exactly 878 | Balanced coupled qubit validation |
| Balanced coupled LiH 4e FCI bound | D_e > 0 | Only bound 4e config |
| Balanced coupled LiH n_max=3 energy | -8.055 Ha (0.20% err) | Convergence validation |
| Balanced coupled variational bound | E > exact at all R | Variational principle |
| Balanced coupled BeH₂ Pauli terms (Q=50) | exactly 2,652 | Polyatomic balanced validation |
| Balanced coupled BeH₂ 1-norm < composed | 304.7 < 354.9 Ha | 1-norm advantage verification |
| Balanced coupled BeH₂ 4e FCI | 10.7% error | Polyatomic energy validation |
| Balanced coupled H₂O Pauli terms (Q=70) | exactly 5,798 | Non-collinear balanced validation |
| Balanced coupled H₂O 1-norm | 1,509.3 Ha | Non-collinear 1-norm |
| Balanced coupled non-collinear V_ne | direction=(0,0,-1) matches nuc_parity to 1e-14 | Wigner D rotation validation |
| FrozenCore Z_eff asymptotic | Z_eff(0)≈Z, Z_eff(∞)≈Z-10 | Ne-like screening validation |
| FrozenCore density normalization | integral = 10 ± 1% | Core electron count |
| Second-row Pauli scaling | Q^2.50 | O(Q^2.5) universality |
| NaH balanced Pauli (Q=20) | exactly 239 | Second-row qubit validation |
| Atomic classifier Z=11-18 | 97 tests pass | Second-row classification |
| Wigner D l=2 orthogonality | R^T R = I to 1e-12 | l=2 rotation validation |
| Wigner D l=2 determinant | det(R) = +1 | Proper rotation check |
| Wigner D l=1 legacy match | < 1e-13 | General vs legacy consistency |
| Block rotation l=0,1,2 | R^T R = I to 1e-12 | Mixed-l block validation |
| NaH n_max=3 build | succeeds | l=2 rotation unblock |
| NaH n_max=3 FCI bound | E(3) < E(2) at all R | Variational convergence |
| LiF composed Pauli terms (Q=70) | exactly 778 | Multi-center qubit validation |
| CO composed Pauli terms (Q=100) | exactly 1111 | Multi-center qubit validation |
| N₂ composed Pauli terms (Q=100) | exactly 1111 | Multi-center isostructural invariance |
| CO = N₂ Pauli count | identical | Isostructural invariance (multi-center) |
| F₂ composed Pauli terms (Q=100) | exactly 1111 | Multi-center qubit validation |
| NaCl composed Pauli terms (Q=50) | exactly 556 | Mixed frozen-core multi-center |
| CH₂O composed Pauli terms (Q=120) | exactly 1333 | Multi-center polyatomic validation |
| C₂H₂ composed Pauli terms (Q=120) | exactly 1333 | Multi-center polyatomic validation |
| C₂H₆ composed Pauli terms (Q=160) | exactly 1777 | Multi-center polyatomic validation |
| Composed Pauli/Q ratio | 11.11 ± 0.1 | Universal linear scaling law |
| All TM hydrides composed Pauli terms (Q=30) | exactly 277 (non-identity) | Transition metal qubit validation (isostructural: all 10 identical). Note: Track CZ/DA reported 278 including the identity term; v2.8.0 standardized on excluding identity from N_pauli across all builders. |
| SrH / BaH composed Pauli terms (Q=20) | exactly 222 (non-identity) | Heavy-atom alkaline-earth monohydride validation ([Kr], [Xe] frozen cores); isostructural with KH, NaH, CaH (Sprint 3 HA-C, v2.12.0). |
| SrH_rel / BaH_rel composed Pauli terms (Q=20) | exactly 942 (non-identity) | Relativistic heavy-atom monohydride validation; isostructural with CaH_rel (post-TR, Sprint 4 v2.15.0; pre-TR was 534). |
| SrH_rel / BaH_rel / CaH_rel λ_ni | bit-identical 13.87 Ha | Cross-species relativistic 1-norm invariance (frozen core screens Z, spin-orbit sees Z_eff=2 uniformly). |
| SrH_rel / BaH_rel / CaH_rel QWC | bit-identical 52 groups | Cross-species QWC structural invariance. |
| d-only block ERI density | 4.0% | d-orbital Gaunt sparsity |
| d-block Pauli/Q (composed) | 9.23 (< main-group 11.11) | d-orbital sparsity advantage |
| Nested Be Pauli terms (Q=10) | exactly 112 | Nested encoding qubit validation |
| Nested Be 1-norm < composed | 18.95 < 121.35 Ha | PK elimination 1-norm advantage |
| H-set ERI density < uncoupled | 9.2% < 12.5% (l_max=1) | 6j recoupling sparsity |
| LiH/BeH₂/H₂O Pauli unchanged | 334/556/778 | Backward compatibility regression |
| He 2D variational bound | E > exact at all l_max | Variational principle |
| He 2D l_max monotonic | E decreasing with l_max | Convergence validation |
| He 2D breaks adiabatic floor | < 0.10% at l_max=5 | Non-adiabatic improvement |
| He 2D cusp-corrected | < 0.01% at l_max=4 | Track DI Phase 1 target |
| He 2D radial converged | stable at n_R ≥ 20 | Radial basis validation |
| Casimir CI n_max=1 k* | 9/4 (exact) | SC Hartree screening |
| Casimir CI n_max=1 E_var | -729/256 (exact) | Variational Hartree screening |
| Casimir CI polynomial structure | residual < 1e-10 | H(k) = Bk + Ck² |
| Casimir CI variational bound | E_var > exact at all n_max | Variational principle |
| Casimir CI n_max=3 error | < 2% (variational) | Algebraic CI accuracy |
| Graph-native CI n_max=7 error | < 0.20% | Graph-native CI accuracy |
| Graph-native CI n_max=8 error | 0.207% (2,262 configs) | Exact algebraic float integrals |
| Graph-native CI n_max=9 error | 0.201% (3,927 configs) | Exact algebraic float integrals |
| He 2D variational self-consistent cusp | < 0.020% | Self-consistent cusp correction (l_max=7, n_R=35) |
| Algebraic Slater R^k machine precision | < 1.5e-12 | `compute_rk_float()` vs exact Fraction |
| casimir_ci F²(2p,2p) corrected | 45/512 (was 43/512) | Typo fix validated by hypergeometric evaluator |
| Graph-native CI beats diagonal | error_graph < error_diag | Graph topology dominance |
| Graph-consistent = graph-native | < 1e-10 Ha | FCI basis invariance |
| Graph-consistent full spectrum | all evals match to 1e-8 | FCI orbital rotation invariance |
| Graph-native CI Z_c crossover | Z_c ≈ 1.84 at n_max=7 | Graph validity boundary |
| Graph-native H⁻ over-binding | E_CI < E_exact at n_max≥2 | Non-variational characterization |
| Standard FCI H⁻ variational | E_CI > E_exact at all n_max | Variational bound verification |
| He energy decomposition verified | <h1>+<V_ee>=E to 1e-12 | Decomposition correctness |
| V_ee full rank (graph eigenbasis) | rank = dim at all n_max | Cusp rank characterization |
| Diagonal V_ee converged by n_max=3 | < 0.01 mHa change | Mean-field convergence |
| PsH bound (adiabatic) | V_min < -0.5 Ha | Exotic atom binding |
| PsH energy l_max=3 | 4.1% error | Exotic atom accuracy |
| 111 Pauli per s/p block | exactly 111 | Composed Pauli derivation |
| 111 = 55 + 56 decomposition | 55 direct + 56 exchange | Pauli channel decomposition |
| TC gamma_opt(l=3) ≈ 0.10 | sweet spot | TC optimization |
| Dirac-on-S³ π-free certificate (Weyl sector) | 0 non-rationals, n_max=6 | Analog of Paper 24 §III Bargmann-Segal certification; every \|λ_n\| is exact sympy Rational (n + 3/2), every g_n^Weyl is positive int ((n+1)(n+2)) |
| Dirac-on-S³ π-free certificate (Dirac sector) | 0 non-rationals, n_max=6 | Every \|λ_n\| is Rational, every g_n^Dirac is positive int (2(n+1)(n+2)) |
| Dirac-on-S³ label generator | exactly g_n labels per level | `spinor_labels_at_n` generates exactly (n+1)(n+2) Weyl or 2(n+1)(n+2) Dirac labels; stronger invariant than bare π-free |
| Dirac Δ⁻¹ identity | g_3^Dirac = 40 exactly | Phase 4H SM-D identity (Δ = 1/(g_3^Dirac)) reproduced in D1 API as `delta_inverse_identity() == (40, Rational(1,40))` |
| Fock ↔ CH convention conversion | invertible at all n | `fock_to_ch(ch_to_fock(n)) == n` for n = 1..10; label compatibility with scalar graph (n_Fock = n_CH + 1) |
| D2 cumulative Dirac trace closed form | exact Rational | Σ g_n^Dirac = (N+1)(N+2)(2N+3)/3 symbolic identity, verified for N = 0..5 |
| D2 \|λ_{m-1}\|·g_{m-1}^Weyl = 42 at m=3 | (m-1)(m+2) = 10 uniquely | Single-point Dirac/Weyl coincidence with B = 42 occurs only at m=3, sympy exact |
| D3 Dirac Dirichlet at s=4 | 2ζ(2) + 2ζ(3) | `summation(2*m*(m+1)*m**(-4), (m,1,oo)) == pi**2/3 + 2*zeta(3)` symbolic |
| D3 Weyl Dirichlet at s=4 | ζ(2) + ζ(3) | Factor of 2 difference from Dirac, sympy exact |
| D4 Hopf charge partition sum | Σ mult = g_n^Dirac | Sympy-exact at n_CH = 0..5; 40 = 20 half-integer-charge + 20 integer-charge at n_CH = 3 |
| D4 Dirac Hurwitz spectral zeta at s=4 | π² − π⁴/12 | `summation(2*(n+1)*(n+2)/(n+Rational(3,2))**4, (n,0,oo))` symbolic closed form |
| T0 d_spinor pair-diag at l_max=0..5 | exact rational | 1/4, 11/20, 553/724, 101/118, 2533/2820, 9611/10396 from sympy Wigner 3j |
| T0 d_spinor full-Gaunt at l_max=0..5 | exact rational | 1/4, ..., 0.923 (sec'dry table, physically correct Coulomb rule) |
| T0 d_spinor ≤ d_scalar | monotonic ∀ l_max | Spinor basis sparsity bounded by scalar sparsity |
| T0 scalar density reproduces Paper 22 | bit-exact | Pair-diag convention matches Paper 22 Table §III |
| `dirac_matrix_elements` module tests | 108 tests pass | Angular (Szmytkowski) + diagonal radial (Bethe-Salpeter) + κ↔(l,σ) bridge + Kramers-Pasternak direct integration |
| σ·r̂ reduction identity | `(−κ, m_j, −1)` exact | Szmytkowski Eq. 2.7 verification |
| Diagonal ⟨1/r³⟩_{n,l} hydrogenic | Z³/[n³·l(l+½)(l+1)] | T1 closed form, verified to sympy machine precision |
| T1 ⟨1s\|r\|2s⟩ at Z=1 | −32√2/81 | Off-diagonal radial sympy integration sanity |
| `spin_orbit` module tests | 22 tests pass | H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)] closed form |
| SO Kramers cancellation at l=0 | H_SO = 0 exact | κ=−1 forces numerator zero before ⟨1/r³⟩ evaluated |
| SO Z⁴ scaling | (Z/Z_ref)⁴ symbolic | Verified at Z ∈ {1, 3, 4, 38} |
| 2p doublet splitting (Z=1) | α²/32 exact | Breit-Pauli fine-structure benchmark |
| `spin_ful_composed` module tests | 13 tests pass | Full composed-rel pipeline regression |
| LiH/BeH₂/H₂O scalar regression | 334/556/778 Pauli preserved | Bit-exact scalar path unchanged when relativistic=False |
| LiH relativistic Pauli at n_max=1 | exactly 9 | Matches scalar (no spin-orbit at l=0) |
| LiH rel/scalar Pauli ratio at n_max=2 | ∈ [3.7, 4.9] | Pinned at 4.24× in regression suite (post-TR, Sprint 4 v2.15.0) |
| Spinor FCI at α=0 matches scalar FCI | \|ΔE\| < 1e-10 Ha | TR regression test (Sprint 4): jj reduced-matrix-element phase fix |
| Spinor composed block-diagonal ERI | zero cross-block entries | Factorization preserved in relativistic path |
| α → 0 zeroes H_SO diagonal | exact zero | Non-relativistic limit verification |
| `spinor_certificate` module tests | 25 tests pass | Ring R_sp = ℚ(α²)[γ]/(γ²+(Zα)²−1) enforcement |
| Contamination rejection (π, π², ζ(3), log, E₁, unregistered) | raises SpinorTaxonomyError | Six negative controls |
| T3 H_SO block R_sp membership | passes n_max ≤ 4 | Every κ-branch coefficient in ring |
| Sunaga RaH-18q baseline | 47,099 Pauli (published) | Single calibrated cell from PRA 111, 022817 |
| GeoVac rel/Sunaga RaH-18q ratio | 0.011×–0.017× native-Q | Resource advantage verified |
| Fine-structure Li 2²P splitting | sign + OoM correct | Breit-Pauli + Z_eff sanity |
| Fine-structure He 2³P span | sign + OoM correct | 66% relative error (accepted) |
| Fine-structure Be 2s2p ³P span | sign + OoM correct | 78% relative error (accepted) |
| Sprint 5 CP: Li 2²P doublet splitting | < 20% | +8.89% err with std conv (Z_val=1, Z_eff=1) |
| Sprint 5 CP: Be 2s2p ³P span (P₀-P₂) | < 20% | +2.76% err with std conv + Slater Z_eff=1.95 |
| Sprint 5 CP: Be 2s2p ³P individual (3 splittings) | < 20% | all three pass (+18.9%, +6.3%, +2.8%) |
| Sprint 5 CP: 2s2p ↔ 2p² parity-forbidden coupling | exactly 0 | symbolic parity verification |
| Sprint 5 CP: Li core polarization worsens accuracy | err_cp > err_bare | MB α_d/r_c increases Δζ in wrong direction |
| Breit zero at α=0 | breit_eri_count = 0 | α² prefactor kills all Breit terms when alpha_num=0 |
| Breit 1-norm order of magnitude | rel_shift < 1% | Breit 1-norm shift is O(α²) ~ 10⁻⁴ relative to Coulomb |
| Breit Pauli count unchanged | N_pauli(breit) ≤ 2× N_pauli(no breit) | Same Gaunt selection rules; bounded increase |
| Breit block diagonality | cross_block_eri_count = 0 | Breit ERI remains block-diagonal |
| Breit Hermiticity | max(|imag(coeff)|) < 1e-10 | All qubit operator coefficients are real |
| Breit He 2³P splittings | O(10⁻⁵) Ha total | Drake J-pattern physical, nonzero splittings |
| Breit radial Z³ scaling | ratio = 8.0 (Z=2/Z=1) | Retarded integrals scale as Z³ |
| Kramers-Pasternak matches Pochhammer n_r=0 | < 1e-10 | Direct integration reproduces T7 Pochhammer results |
| Kramers-Pasternak matches Hellmann-Feynman ⟨1/r⟩ | < 1e-10 | Direct integration matches all-state HF formula |
| Kramers-Pasternak ⟨r⁰⟩ = 1 normalization | exact | Normalization integral equals 1 for all states |
| Seeley-DeWitt a₀ on unit S³ | √π (exact sympy) | Heat kernel coefficient from Dirac D² spectrum |
| Seeley-DeWitt a₁ on unit S³ | √π (exact sympy) | Curvature correction (R_scalar/6 = 1) |
| Seeley-DeWitt a₂ on unit S³ | √π/8 (exact sympy) | Vacuum polarization coefficient source |
| Vacuum polarization coefficient | 1/(48π²) (exact) | Standard Dirac fermion VP from S³ spectral data |
| β(α) QED one-loop | 2α²/(3π) (exact) | Reproduces standard QED beta function |
| No odd-zeta at one loop | structural theorem | T9 guarantees ζ_{D²}(s) is polynomial in π² |
| D_even(4) decomposition | π²/2 − π⁴/24 − 4G + 4β(4) (PSLQ 80 digits) | Vertex parity exposes Catalan G and Dirichlet β(4) |
| D_odd(4) opposite sign | π²/2 − π⁴/24 + 4G − 4β(4) (PSLQ 80 digits) | Opposite Dirichlet content in odd-n sub-sum |
| D_even(4) + D_odd(4) cancellation | = D(4) to 1e-60 | Dirichlet L-values cancel in full (unrestricted) sum |
| Vertex selection rule consistency | matches hodge1_s3 exactly | SO(4) triangle + parity: n₁+n₂+n_γ odd |
| Fine-structure Dirac formula (n<=4) | exact symbolic | Dirac formula verified for all (n,l,j) through n=4 |
| Dirac accidental degeneracy | 6/6 pairs confirmed | E_FS depends on (n,j) only |
| gamma radial <1/r> NR limit | Z/n^2 | Hellmann-Feynman all-state formula |
| gamma radial <1/r^2> n_r=0 NR limit | Z^2/(n^3(l+1/2)) | Pochhammer ratio |
| gamma radial <1/r^3> n_r=0 NR limit | Z^3/(n^3 l(l+1/2)(l+1)) | Pochhammer ratio |
| zeta_{D^2}(2) | pi^2 - pi^4/12 | Squared Dirac spectral zeta (T9) |
| zeta_{D^2}(s) pi^{even} only | theorem | No odd-zeta content at any integer s |
| Weyl density Dirac S³ | O(1/n_max) convergence | Spectral-to-momentum correspondence |
| Sunset R-scaling power | exact R^10 at (s1,s2,p)=(2,2,1) | Dimensional analysis verification |
| D(s,R)/R^s R-independence | ratio constant to 1e-30 | Transcendental class is R-independent |
| D(5) = 14ζ(3) − 31/2·ζ(5) | exact to 80 dps | ζ(3) structural identification |
| Even-s stays π^{even} at all R | structural theorem | One-loop parity persists |
| Screened radial solver hydrogenic limit | < 1% at Z_eff=const | Reproduces analytical ⟨1/r³⟩ |
| Screened SO enhancement Ca/Sr/Ba | 12×/61×/144× over Z_eff=2 | Core penetration quantified |
| Screened SO CaH/SrH/BaH splitting | within 70% of physical | Leading-order Breit-Pauli |
| NaH balanced FCI overattraction | no equilibrium | Frozen-core cross-V_ne limitation |
| MgH₂ balanced FCI overattraction | no equilibrium | Same mechanism as NaH |
| Self-energy Σ(n_ext=0) structural zero | exactly 0 | Vertex parity selection rule proof |
| Self-energy Σ(n_ext=1) sign | < 0 | Physical (ground state protected, excited states shifted) |
| Self-energy convergence monotonic | monotonic with n_max | Spectral sum convergence |
| Factorized matches direct three-loop | < 1e-25 at n_max=15 | O(N³) vs O(N⁵) algorithmic equivalence |
| Factorized n_max=50 speed | < 60 seconds | O(N³) performance benchmark |
| Factorized convergence monotonic | vals increasing with n_max | Three-loop sum convergence |
| S_min 200 digits verified | 3 independent methods agree | mpmath.nsum, Euler-Maclaurin, direct sum |
| S_min PSLQ irreducibility | 15 failures across 100+ basis | Extended basis including Tornheim-Witten, colored MZV |
| c₂ cross-invariant 8-digit match | |c₂_cross - c₂_apparent| < 2e-8 | Paper 2↔28 bridge verification |
| c₂ symbolic identity | 19/100 - 41π²/25200 (exact sympy) | Rational+π² decomposition |
| c₃ from n_int=0..50 | -5.946(3)×10⁻⁷ at 200.9 sigma | Nonzero; expansion does not terminate at c₂ |
| c₂ T9 consistency | rational + rational·π² only | One-loop π^{even} constraint |
| Paper 28 theorem count | ≥ 4 theorems with proofs | Q-4 exit criterion |
| D₅ c₅(n) exact match n=1..15 | rational equality | Sommerfeld closed-form verification |
| D₅ PSLQ decomposition | residual < 1e-190 | Weight-9 MZV decomposition |
| D₅ ζ(2)ζ(7) cancellation | coefficient = 0 | Product survival rule verification |
| D₆ analytical assembly | 62 matching digits | Weight-11 MZV decomposition from 7 Euler sums |
| D₆ ζ(2)ζ(9) cancellation | coefficient = 0 | Product survival rule verification at p=6 |
| D₆ surviving products | exactly 3 | z(3)z(8), z(4)z(7), z(5)z(6) as predicted |
| Product survival rule | max(0, floor((2p-5)/2)) | Verified D₂..D₆ |
| K/π not in Q-span of D₂..D₆ | PSLQ null | K-Sommerfeld structural separation |
| Graph-native QED tests (non-slow) | 72 pass | GN-5 self-energy + vertex pipeline |
| Graph-native Σ trace (t=0, n_max=2) | exactly 44/3 | Self-energy trace validation |
| Graph-native Λ trace (t=0, n_max=2) | exactly 32/9 | Vertex correction trace validation |
| Graph-native F₂ (t=0, n_max=2) | exactly 5√2/3 | Anomalous moment (irrational, π-free) |
| Graph-native Σ ground-state block ≠ 0 | [[1,1],[1,1]] | Broken structural zero (CG opens couplings) |
| Graph-native Σ eigenvalues (t=0) | {0(×5), 4/3, 2(×2), 4, 16/3} | Self-energy spectrum |
| Graph-native VP Π trace (n_max=2) | exactly 32 | Vacuum polarization trace (GN-4) |
| Graph-native VP Π trace (n_max=3) | exactly 3192/5 | VP extension (GN-6) |
| Graph-native projection C (n_max=3) | 50471424/1779441125 | Rational projection exchange constant (GN-7) |
| Graph-native all π-free | structural | All graph QED quantities algebraic, no transcendentals |
| GN-6 n_max=3 VP tests | 45 pass | Extended VP at n_max=3 |
| GN-7 continuum bridge tests | 63 pass | Projection exchange constant verification |
| F₂(t) even-function | c₁ = c₃ = c₅ = 0 exact | Graph-native anomalous moment parity |
| F₂(κ) convergence | 2.353, 1.873, 1.589 | Monotonic decrease at n_max=2,3,4 |
| F₂(κ) at n_max=5 | 1.39581063 | Extended convergence (numpy) |
| F₂(κ) at n_max=6 | 1.25321124 | Extended convergence (numpy) |
| F₂ power-law exponent | -0.573 (R²=0.99990) | F₂ ~ 3.507 × n^(-0.573) |
| Pendant-edge n_max=5 | Σ(GS) = 1.60000000 = 8/5 | Exact match 2(n-1)/n |
| Pendant-edge n_max=6 | Σ(GS) = 1.66666667 = 5/3 | Exact match 2(n-1)/n |
| F₂(t) n_max=3 rational degree | 16/16 | Even function over ℚ(√2,√3,√5,...) |
| F₂ successive ratios → 1 | 0.796, 0.849, 0.878, 0.898 | Power-law consistency |
| Selection rule census | 1/8 survives | Only Gaunt/CG sparsity survives on graph |
| C × F₂ divergence | grows with n_max | VP projection ≠ vertex projection |
| Dirac graph Rule B selection rules | 4/8 survive | Spinor-recoverable vs vector-required partition |
| Dirac graph GS NOT pendant | degree 2 (A), 5 (B) | Constant across n_max=2,3,4 |
| Dirac graph Σ(GS) ≠ 0 | nonzero both rules | Structural zero not recovered |
| Dirac graph Σ strictly PSD | 0 zero eigenvalues | Unlike scalar graph (5 zeros) |
| Dirac graph π-free | ℚ[√2,√17,√41,√881] | Algebraic, no transcendentals |
| Dirac graph Tr(Σ_A) | exactly 17/2 | Exact algebraic trace |
| Dirac graph Tr(Σ_B) | exactly 103/20 | Exact algebraic trace |
| Dirac graph Furry recovered | tadpole = 0 | Off-diagonal identity vertex parity |
| VP/SE projection ratio constant | CV < 1% across n_max=3,4,5 | Two-tier calibration structure |
| VP/SE ratio ≈ 3/29 | 0.1035 ± 0.0009 | Closest simple rational |
| C_F2 exponent matches F₂ | +0.57 vs −0.57 | Complementary scaling |
| Catalan/β from graph parity | 24 PSLQ null | Clean negative |
| Per-mode projection topological | ρ varies 0 to 1.66 | Not multiplicative |
| Speed regression | < 10% | Performance control |

---

## 11. Topic-to-Paper Lookup

| Topic | Paper | Section | Tier |
|:------|:-----:|:--------|:----:|
| Universal constant -1/16 | 0 | Sec 2 | Core |
| Graph Laplacian method | 1 | Sec 3 | Core |
| O(V) quantum dynamics | 6 | All | Core |
| Rabi oscillations | 6 | -- | Core |
| Delta-kick spectroscopy | 6 | -- | Core |
| AIMD / Langevin thermostat | 6 | -- | Core |
| Dimensionless vacuum proof | 7 | All | Core |
| Schrodinger recovery | 7 | Sec 4 | Core |
| S3 conformal geometry | 7 | Sec 3 | Core |
| V_ee on S3 (node overlap) | 7 | Sec VI | Core |
| Slater F0 master formula | 7 | Sec VI.B | Core |
| SO(3N) generalization | 7 | Sec VI | Core |
| Nuclear rovibrational spectra | 10 | All | Core |
| Prolate spheroidal lattice | 11 | All | Core |
| Neumann V_ee expansion | 12 | Sec III-V | Core |
| Prolate spheroidal CI (H2) | 12 | Sec VI | Core |
| Cusp diagnosis (7.6% gap) | 12 | Sec VII | Core |
| Hyperspherical coordinates | 13 | Sec II | Core |
| Angular eigenvalue (Gaunt) | 13 | Sec III | Core |
| Adiabatic potential curves | 13 | Sec IV | Core |
| Fiber bundle structure | 13 | Sec VII | Core |
| Natural geometry hierarchy | 13 | Sec VIII | Core |
| Ab initio spectroscopy | 13 | Sec IX | Core |
| Algebraic structure (SO(6)) | 13 | Sec XII | Core |
| Algebraic angular solver | 13 | Sec XIII | Core |
| Coupled-channel convergence | 13 | Sec XIII.C | Core |
| Hellmann-Feynman P-matrix | 13 | Sec XIII.A | Core |
| Adiabatic bottleneck (l_max) | 13 | Sec XIII.B | Core |
| Perturbation series μ(R) | 13 | Sec XII-XIII | Core |
| Transcendental boundary (μ(R)) | 13 | Sec XII.B | Core |
| Spectral Laguerre (Level 2) | 11 | Sec V | Core |
| Algebraic Laguerre moments | 11 | Sec V | Core |
| Spectral Laguerre (Level 3) | 13 | — | Core |
| Algebraic Laguerre S, K (Level 3) | 13 | — | Core |
| Spectral Laguerre (Level 4) | 15 | Sec VI.H | Core |
| Spectral Jacobi angular (Level 4) | 15 | Sec VI.I | Core |
| 2D solver in composition | 15, 17 | Sec VI.D, — | Core |
| Qubit Pauli scaling | 14 | All | Core |
| Structural sparsity | 14 | Sec III | Core |
| QWC measurement groups | 14 | Sec III.E | Core |
| Pauli 1-norm scaling | 14 | Sec III.F | Core |
| Trotter error bounds | 14 | Sec III.F | Core |
| Genuine Gaussian baselines | 14 | Sec III.G | Core |
| VQE head-to-head benchmark | 14 | Sec III.G | Core |
| Equal-qubit Gaussian comparison | 14 | Sec III.G | Core |
| Composed qubit Hamiltonians | 14 | Sec IV | Core |
| Composed Pauli scaling (Q^2.5) | 14 | Sec IV.B | Core |
| Trenev et al. Gaussian baselines | 14 | Sec IV.C | Core |
| Level 4 mol-frame hyperspherical | 15 | All | Core |
| Mol-frame charge function | 15 | Sec III | Core |
| Multichannel expansion | 15 | Sec V | Core |
| Heteronuclear extension | 15 | Sec V.D | Core |
| Variational 2D solver | 15 | Sec VI.D | Core |
| HeH+ convergence | 15 | Sec VI.E | Core |
| Double-adiabatic fiber bundle | 15 | Sec VII.C | Core |
| Chemical periodicity (S_N reps) | 16 | All | Observation |
| Structure types A/B/C/D/E | 16 | Sec III | Observation |
| Hierarchical decomposition | 16 | Sec IV | Observation |
| Dirac limit (Z~137) | 16 | Sec VI | Observation |
| Exchange constants (Weyl-Selberg) | 18 | All | Core |
| Stieltjes seed e^a E₁(a) | 18 | Sec II.B | Core |
| Transcendence hierarchy | 18 | Sec IV, VI | Core |
| α as exchange constant | 18 | Sec V | Core |
| 1/r₁₂ embedding exchange constant | 18 | Sec II (new) | Core |
| Observable classification by transcendental content | 18 | Sec VI | Core |
| Cusp dimensionality obstruction | 7, 15 | — | Core |
| Cusp correction (Schwartz) | 15 | — | Core |
| Cusp diagnosis (7.6% gap) | 12 | Sec VII | Core |
| Composed natural geometries | 17 | All | Core |
| Core-valence fiber bundle | 17 | Sec II | Core |
| Ab initio Phillips-Kleinman | 17 | Sec IV | Core |
| Rho-collapse cache | 17 | Sec V | Core |
| l_max divergence root cause | 17 | Sec VI.A | Core |
| LiH/BeH+ benchmarks | 17 | Sec VI | Core |
| Bond sphere theory | 8-9 | All | Core |
| Sturmian structural theorem | 8-9 | Sec IV | Core |
| SO(4) selection rules | 8-9 | Sec III | Core |
| Fine structure alpha (Hopf bundle) | 2 | Sec 3-5 | Conjecture |
| Hopf fibration spectral invariants | 2 | Sec 3 | Conjecture |
| Cubic self-consistency (alpha) | 2 | Sec 4 | Conjecture |
| Spectral dimension d_s | 3 | Sec 4 | Conjecture |
| Holographic entropy S | 3 | Sec 5 | Conjecture |
| Central charge c | 3 | Sec 6 | Conjecture |
| Mass-independence | 4 | Sec 3-4 | Conjecture |
| Muonic hydrogen | 4 | Sec 5 | Conjecture |
| Contact geometry | 4 | Sec 5 | Conjecture |
| Comprehensive framework | 5 | All | Conjecture |
| Level 4N full N-electron solver | 17 | Sec VIII | Core |
| S₄ antisymmetry (4-electron) | 17 | Sec VIII | Core |
| N-electron spectral compression | 17 | Sec VIII | Core |
| PK as composition exchange constant (N-electron validation) | 18 | Sec IV | Core |
| N-electron 2D variational solver | 17 | Sec VIII | Core |
| N-electron quantum encoding comparison | 14 | Sec IV | Core |
| Two-center integrals | 19 | Sec III-IV | Supporting |
| Shibuya-Wulfman integrals | 19 | Sec IV | Supporting |
| Coulomb Sturmian (molecular) | 19, 8-9 | Sec IV | Supporting |
| PK elimination | 19 | Sec I, VI | Supporting |
| Coupled composition | 19 | All | Supporting |
| Balanced coupled composition | 19 | Sec V–VI | Supporting |
| Cross-center V_ne (Shibuya-Wulfman) | 19 | Sec III-IV | Supporting |
| Multipole expansion termination | 19 | Sec III.C | Supporting |
| Non-collinear V_ne | 19 | Sec V | Supporting |
| Wigner D-matrix rotation | 19 | Sec V | Supporting |
| Resource benchmarks (quantum) | 20 | All | Applications |
| FCI PES characterization | 20 | Sec III | Applications |
| Gaussian comparison (matched accuracy) | 20 | Sec IV | Applications |
| Installation and ecosystem | 20 | Sec V | Applications |
| Frozen-core tabulation | 17, 19 | Sec II, Sec VI | Core |
| Second-row generalization | 14, 20 | — | Core |
| Clementi-Raimondi Z_eff | 17 | — | Supporting |
| Multi-center molecules | 14, 20 | — | Core |
| Multi-center isostructural invariance | 14 | — | Core |
| Nuclear scaling (Pauli vs atoms) | 14, 20 | — | Core |
| Composed linear scaling law (11.11×Q) | 14 | — | Core |
| 6j recoupling sparsity | 14 | Sec III.D | Core |
| Nested hyperspherical encoding | 14 | Sec IV.M | Core |
| Composed factorization necessity (nested validation) | 17 | Sec VIII.C | Core |
| 2D variational solver (Level 3 He) | 13 | — | Core |
| Non-adiabatic R-α correlation | 13 | — | Core |
| Schwartz cusp correction (Level 3) | 13 | — | Core |
| Synthesis: geometric vacuum framework | 21 | All | Synthesis |
| Packing axioms (synthesis) | 21 | Sec II | Synthesis |
| S³ proof chain (synthesis) | 21 | Sec III | Synthesis |
| N-electron S^(3N-1) (synthesis) | 21 | Sec IV | Synthesis |
| Exchange constant taxonomy (synthesis) | 21 | Sec V | Synthesis |
| Hopf α conjecture (synthesis) | 21 | Sec VI | Synthesis |
| Atomic spectra research program | 21 | Sec VII | Synthesis |
| Graph-native CI (He) | 7, 13, 18 | Note added v2.6.0 | Core |
| FCI basis invariance | 7, 13, 18 | Note added v2.6.0 | Core |
| Three-layer exchange structure (He) | 18, 21 | Note added v2.6.0 | Core |
| Molecular limitation (honest assessment) | 21 | Sec VIII | Synthesis |
| Hopf graph as lattice gauge structure (synthesis) | 25 | All | Synthesis |
| Graph Hodge theory (node/edge Laplacian, Betti numbers) | 25 | Sec II | Synthesis |
| Edge Laplacian L₁ = B^T B as gauge propagator (conjecture) | 25 | Sec III, IV | Synthesis |
| S² Hopf quotient graph (lattice gauge base) | 25 | Sec III.B | Synthesis |
| U(1) gauge phases from ladder operators | 25 | Sec III.C | Synthesis |
| (B, F, Δ) as gauge-theoretic invariants (reframing) | 25 | Sec IV.A | Synthesis |
| Wilson/graph-Hodge dictionary | 25 | Sec II.E | Synthesis |
| Lattice gauge theory on Fock-projected S³ | 25 | Sec III, VI | Synthesis |
| Why the lattice-gauge / graph-Hodge gap remained | 25 | Sec VI.D | Synthesis |
| SU(3) gauge on Bargmann-Segal S⁵ (open question) | 25 | Sec VII.A | Synthesis |
| Rank-2 Breit as 2-form gauge (open question) | 25 | Sec VII.B | Synthesis |
| Potential-independent angular sparsity | 22 | All | Core |
| ERI density vs l_max (theorem + data) | 22 | Sec III, IV | Core |
| Universal vs Coulomb-specific partition | 22 | Sec VI; 23 | Core |
| Nuclear magic numbers from graph | 23 | Sec II | Applications |
| Mayer-Jensen spin-orbit in HO basis | 23 | Sec II | Applications |
| Deuteron qubit Hamiltonian (Minnesota) | 23 | Sec IV | Applications |
| He-4 qubit Hamiltonian (2p+2n) | 23 | Sec V | Applications |
| Two-species tensor-product JW encoding | 23 | Sec IV | Applications |
| Composed nuclear-electronic deuterium | 23 | Sec VI | Applications |
| Fock projection rigidity theorem | 23 | Sec III | Applications |
| Hyperfine I·S cross-register coupling | 23 | Sec VI | Applications |
| Multi-scale coefficient hierarchy | 23 | Sec VI | Applications |
| Bargmann-Segal lattice (HO discretization) | 24 | All | Core |
| π-free graph certificate (HO) | 24 | Sec III | Core |
| HO rigidity theorem (dual of Fock rigidity) | 24 | Sec V | Core |
| Coulomb/HO structural asymmetry | 24 | Sec IV | Core |
| Calibration π is Coulomb-specific | 24 | Sec IV | Core |
| First-order vs second-order operator distinction | 24 | Sec IV | Core |
| S^5 Hardy space / Bargmann space / SU(3) (N,0) | 24 | Sec II | Core |
| Linear-affine HO projection | 24 | Sec II.B | Core |
| Cusp characterization (energy decomposition) | 13, 18 | — | Core |
| V_ee spectral structure (graph eigenbasis) | 13, 18 | — | Core |
| Graph validity boundary (Z_c ≈ 1.84) | 7, 13 | — | Core |
| 111 Pauli count derivation | 14 | — | Core |
| TC gamma optimal scan | 13 | — | Core |
| PsH exotic atom | 13 | — | Core |
| H⁻ negative ion (graph-native CI) | 13 | — | Core |
| Single-center molecular encoding | 8-9 (GUARDRAIL), FCI-M (GUARDRAIL), Track DF | — | Methods |
| Shared exponent molecular | 8-9 (GUARDRAIL — Sturmian theorem) | Sec IV | Methods |
| Graph concatenation molecular | FCI-M (GUARDRAIL — R-independent h1) | All | Methods |
| Nested hyperspherical molecular | Track DF record, 8-9, 17 §F.2 | — | Methods |
| Unified basis molecular | 8-9 (GUARDRAIL), Track DF Sprint 5 (Löwdin) | — | Methods |
| Sturmian basis (molecular) | 8-9 (negative theorem), Track BU (1-norm inflation) | — | Methods |
| Hypergeometric Slater integrals | 7 | Sec VI.B | Core |
| DUCC downfolding (core potential) | 17, 18 | — | Core |
| PK l_max divergence root cause | 17, 18 | — | Core |
| Energy graph for V_ee on S³ (characterization) | 18 | Note added v2.9.1 | Core |
| Dirac-on-S³ infrastructure (Camporesi-Higuchi spectrum, spinor harmonics) | 2 (§IV rewrite), 23 | — | Conjecture/Applications |
| Camporesi-Higuchi spectrum on S³ | 2, 23 | §IV | Conjecture |
| Spinor spherical harmonics on S³ | 2 | §IV | Conjecture |
| π-free certificate for Dirac-on-S³ | 2 | §IV | Conjecture |
| Orduz Hopf-equivariant Dirac decomposition | 2 | §IV | Conjecture |
| 20⊕20 charge-parity split of Δ⁻¹ | 2 | §IV | Conjecture |
| ζ(3) in the Dirac sector at s = d_max | 2 | §IV.5 (new) | Conjecture |
| Odd-zeta vs even-zeta content (spinor vs scalar Laplacian) | 2 | §IV.5 | Conjecture |
| Three-homes theorem for B, F, Δ | 2 | §IV.1 (new) | Conjecture |
| Cross-sector structural coincidence (K combination rule) | 2 | §IV.6 | Conjecture |
| Dirac Fock rigidity (deferred) | 23 | — | Applications |
| First-order vs second-order spectral operators (π content) | 24 | §IV-V | Core |
| Spinor composed encoding (Tier 2) | 14 | §V (new) | Core |
| Dirac matrix elements in (κ, m_j) basis | 14 | §V | Core |
| Szmytkowski spinor matrix elements | 14 | §V | Core |
| Martínez-y-Romero radial recursion (reserved) | 14, 18 | §V, §IV | Core |
| jj-coupled Coulomb ERI (Dyall §9 separable form) | 14, 22 | §V, §III | Core |
| Breit-Pauli spin-orbit in (κ, m_j) | 14, 18 | §V, §IV | Core |
| Kramers cancellation at l=0 | 14 | §V | Core |
| d_spinor(l_max) sparsity density | 22 | spinor section (new) | Core |
| Spinor-block angular sparsity theorem | 22 | spinor section | Core |
| Pair-diagonal vs full-Gaunt convention | 22 | spinor section | Core |
| Universal partition (sharpened to spinor) | 22 | spinor section | Core |
| Sunaga 2025 head-to-head | 14, 20 | §V, Tier-2 table (new) | Core/Applications |
| OpenFermion-Dirac baseline (deferred SI cells) | 20 | Tier-2 table | Applications |
| Relativistic resource benchmarks (LiH/BeH/CaH) | 20 | Tier-2 table | Applications |
| Fine-structure sanity (He/Li/Be) | 14, 20 | §V, Tier-2 table | Core/Applications |
| α² and γ taxonomy (spinor-intrinsic) | 18 | §IV subtier (new) | Core |
| Ring R_sp = ℚ(α²)[γ]/(γ²+(Zα)²−1) | 18 | §IV | Core |
| Operator-order × bundle grid (4 cells) | 18 | §IV | Core |
| π-free spinor certificate | 18 | §IV | Core |
| First-order operator on spinor bundle | 18 | §IV | Core |
| Lichnerowicz formula on S³ | 18 | §IV (T9 result) | Core |
| Squared Dirac spectral zeta | 18 | §IV (T9 result) | Core |
| Darwin term | 14 | §V | Core |
| Mass-velocity term | 14 | §V | Core |
| Operator-order transcendental discriminant | 18, 24 | §IV, §V | Core |
| Dirac fine-structure formula | 14 | §V | Core |
| Kramers-Pasternak direct integration | 14, 18 | §V, §IV | Core |
| One-loop QED vacuum polarization on S³ | 18 | §IV (new) | Core |
| Seeley-DeWitt heat kernel coefficients | 18 | §IV | Core |
| QED beta function from spectral data | 18 | §IV | Core |
| Motivic weight and operator order | 18 | §IV (new) | Core |
| ζ(3) as period of Q(−3) | 18 | §IV | Core |
| Zagier dimension conjecture | 18 | §IV | Core |
| Catalan constant G from vertex parity | 18 | §IV (new) | Core |
| Dirichlet beta function β(s) on S³ | 18 | §IV (new) | Core |
| Vertex parity selection rule (two-loop sunset) | 18 | §IV (new) | Core |
| Two-loop sunset diagram on S³ | 18 | §IV (new) | Core |
| Three-axis transcendental taxonomy | 18 | §IV (new) | Core |
| Dirichlet L-function L(s, χ₋₄) | 18 | §IV (new) | Core |
| Compactness thesis (Peter-Weyl) | 18 | §III (new) | Core |
| Weyl asymptotic / flat-space limit | 18 | §IV | Core |
| R-dependent spectral sums on S³ | 18 | §IV | Core |
| Petermann two-loop coefficient (structural match) | 18 | §IV | Core |
| Energy-entanglement decoupling | 26 | All | Observation |
| Basis-intrinsic ERI sparsity (step-function transition) | 26 | Sec III | Observation |
| Entanglement network topology (hub migration) | 26 | Sec IV | Observation |
| Core-valence entanglement decoupling | 26 | Sec V | Observation |
| S ~ Z^{-2.56} entanglement scaling | 26 | Sec II | Observation |
| Composed block R-independent entanglement | 26 | Sec VI | Observation |
| Entropy as projection artifact | 27 | All | Core |
| One-body entanglement-inert theorem | 27 | Sec II | Core |
| HO zero-entropy rigidity | 27 | Sec VII.A | Core |
| Universal S_B(w̃_B/δ_B) scaling | 27 | Sec VII.B | Core |
| Cusp hot-node concentration | 27 | Sec IV | Core |
| n^4 area-law as pair-counting | 27 | Sec III | Core |
| One-loop self-energy on S³ | 28 | §6 | Observation |
| Self-energy structural zero theorem | 28 | §6 (Theorem 4) | Observation |
| Vertex correction / anomalous magnetic moment | 28 | §6 | Observation |
| Curvature expansion c₂ = (2-BΔ-FΔ-F/B)/5 (Paper 2 bridge) | 28 | §6.6 | Observation |
| Three-loop O(N³) factorization | 28 | §9 | Observation |
| S_min extended PSLQ (200 digits, 100-element basis) | 28 | §8 | Observation |
| Depth-k tower proposition | 28 | §10 | Observation |
| Graph-native QED construction (finite graph formulation) | 28 | §graph_native_qed | Observation |
| Graph-native vacuum polarization (Π as finite matrix trace) | 28 | §graph_native_qed.2 | Observation |
| Graph-native self-energy and vertex correction | 28 | §graph_native_qed.3 | Observation |
| Graph-native broken structural zero (CG opens couplings) | 28 | §graph_native_qed.4 | Observation |
| Graph-native continuum bridge / projection exchange constant | 28 | §graph_native_qed.5 | Observation |
| Graph-native π-free certificate (all QED algebraic) | 28 | §graph_native_qed.6 | Observation |
| Fock graph Hodge theory (GN-1) | 25, 28 | §graph_native_qed | Observation |
| Graph QED propagators (GN-2, electron + photon) | 28 | §graph_native_qed | Observation |
| F₂(t) even rational function | 28 | §f2_even | Observation |
| F₂(κ) convergence table | 28 | §f2_convergence | Observation |
| Selection rule survival census (1/8 survives) | 28 | §selection_rule_census | Observation |
| Topology-dependent projection (C×F₂ diverges) | 28 | §topology_dependent_projection | Observation |
| Native Dirac graph QED (4/8 selection rules recovered) | 28 | §dirac_graph_qed | Observation |
| Spinor-recoverable vs vector-photon-required partition | 18, 28 | §dirac_graph_qed, §IV | Core/Observation |
| Dirac graph Hodge decomposition (Rule A/B) | 25, 28 | §dirac_graph_qed | Observation |
| GS pendant status (scalar vs Dirac graph) | 28 | §dirac_graph_qed | Observation |
| Furry's theorem recovery on Dirac graph | 28 | §dirac_graph_qed | Observation |
| F₂ power-law convergence (n_max=2..6) | 28 | §f2_convergence | Observation |
| F₂(t) even rational function structure | 28 | §f2_even, §f2_convergence | Observation |
| Pendant-edge theorem (n_max=2..6) | 28 | §pendant_edge | Observation |
| F₂(κ) vs F₂(0) sign flip at n_max=5 | 28 | §f2_convergence | Observation |
| Self-energy trace Tr(Σ)/N_dirac scaling | 28 | §f2_convergence | Observation |
| Three-diagram projection constants (C_VP, C_SE, C_F2) | 28 | §three_diagram_projection | Observation |
| VP/SE constant ratio (two-tier calibration) | 28 | §three_diagram_projection | Observation |
| Per-mode topological gap (graph vs continuum self-energy) | 28 | §three_diagram_projection | Observation |
| Catalan/β from graph parity (clean negative) | 28 | §three_diagram_projection | Observation |
| Vector QED sprint (VQ-1 through VQ-5) | 28 | §vq_sprint | Observation |
| σ-weighted vertex (Pauli trace identity 3×) | 28 | §vq_sprint | Observation |
| Direction-resolved Hodge decomposition (T vs L channels) | 28 | §vq_sprint | Observation |
| Spectral anisotropy (T/L different L₁ spectra) | 28 | §vq_sprint | Observation |
| Cross-channel coupling in photon propagator | 28 | §vq_sprint | Observation |
| Three-tier selection rule partition (graph/spinor/vector) | 28 | §vq_sprint | Observation |
| Full vector QED (direction × σ_μ, 15 configs) | 28 | §vq_sprint | Observation |

---

## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. Status: **algebraic** (closed-form from quantum numbers), **algebraic (implicit)** (defined by polynomial equation P=0 with known coefficient ring; pointwise diag is computational convenience, not mathematical necessity), **algebraic-pending** (algebraic route identified but production code still uses quadrature), **numerical-required** (no known algebraic replacement).

### Level 3 (Hyperspherical — He)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| SO(6) Casimir eigenvalues | algebraic | Exact: ν(ν+4)/2 from representation theory |
| Centrifugal barrier | algebraic | Diagonal in Gegenbauer spectral basis (confirmed Track B, v2.0.6) |
| Nuclear coupling | algebraic | Partial harmonic sums (Eqs. 31-32, Paper 13) |
| V_ee coupling | algebraic-pending | Split-region Legendre structure confirmed; GL quadrature still used in production |
| Centrifugal matrix elements | algebraic | Diagonal in Gegenbauer basis (v2.0.6) |
| Hyperradial overlap (S) | algebraic | Pentadiagonal M2 moment matrix from three-term Laguerre recurrence. S = M2/(8α³). Machine-precision agreement with quadrature (< 1e-14 relative). (Track H, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Derivative kernel B_n = -n/2 L_{n-1} + 1/2 L_n + (n+1)/2 L_{n+1} (tridiagonal). K = bbᵀ/(4α). Machine-precision agreement (< 1e-14 relative). 11× build speedup. (Track H, v2.0.10) |
| Hyperradial potential (V_eff) | algebraic at l_max=0; numerical-required at l_max≥1 | l_max=0: V_eff algebraic via Stieltjes moment decomposition (single transcendental seed e^a·E₁(a), Track P2 v2.0.13). l_max≥1: μ(R) is algebraic over Q(π,√2) (Track P1), but V_eff integrals are non-elementary (radical obstruction: √Δ(R) at l_max=1, Cardano at l_max=2). Quadrature required. |
| Hellmann-Feynman P-matrix | algebraic | Exact from R-independent dH/dR (v2.0.6) |
| Q-matrix (second derivative coupling) | algebraic | Exact Q = PP + dP/dR computed from Hellmann-Feynman quantities (v2.0.6) |
| Coupled-channel radial solve | numerical-required | Coupled ODE integration on R grid |
| Adiabatic potential curves U(R) | algebraic (implicit) | μ(R) satisfies P(R,μ)=0, polynomial in both R and μ, coefficients in Q(π,√2)[R]. Degree L+1 in both variables at l_max=L. Point-by-point diagonalization is computational convenience. (Track P1, v2.0.12) |
| TC double commutator | algebraic | [[H, R sin(α)], R sin(α)] = -kron(S_R, I + cos²α). Exact closed form, terminates at 2nd order for 2 electrons. The cos²α matrix is block-diagonal in l (Y_l orthogonality). |
| TC single commutator | algebraic | [H, R sin(α)] = (-d/dR)⊗sin(α) + (1/R)⊗[Λ²+15/8, sin(α)] + I⊗[C, sin(α)]. D_R via Laguerre derivative (x weight, anti-symmetric). Angular commutator via IBP + V_extra decomposition. |
| Angular Casimir decomposition | algebraic | V_extra = diag(casimir) - IBP_free, where IBP_free[j,k] = <χ_j'|χ_k'>. Extracts the S⁴ metric cotangent-derivative terms algebraically from the known eigenvalues. Angular derivatives from d/du C_n^λ(u) = 2λ C_{n-1}^{λ+1}(u). |
| sin/sin²/cos² angular matrices | algebraic | Block-diagonal in l by Y_l(θ₁₂) orthogonality. Cross-channel elements are zero because sin(α) preserves θ₁₂ quantum number l. (v2.8.2) |

### Level 2 (Prolate Spheroidal — H₂⁺)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Angular η-eigenfunctions | algebraic | Legendre spectral basis (n_basis=50), same as Mitnik et al. angular component |
| Azimuthal m | algebraic | Exact separation of variables, integer quantum number |
| Radial ξ-solver (m=0) | algebraic | Ordinary Laguerre basis (n_basis=20, α-adapted). 250× dimension reduction vs FD. Three-term recurrence for ALL matrix elements — zero quadrature. Machine-precision agreement with quadrature (< 1e-14). (v2.0.9) |
| Radial ξ-solver (m≠0) | algebraic | Associated Laguerre basis L_n^{|m|}(x) with weight x^|m|·e^{-x}. Partial-fraction decomposition of centrifugal 1/x singularity into lowered moment M_{-1} (algebraic, DLMF 18.9.13) + Stieltjes integral J (three-term recurrence). Single transcendental seed e^a·E₁(a); all other elements algebraic. Associated basis converges faster than ordinary for m=1 (stable by N=10). (Track J, v2.0.10) |
| Separation parameter c² | numerical-required | Iterative root-finding (Brent method) for self-consistency between angular and radial equations |
| Coupled-channel ceiling | characterized | Error floor 0.19-0.20% from adiabatic approximation; algebraic convergence ~l_max⁻²; 3 channels sufficient; n_basis and R-grid converged. Sub-0.1% requires non-adiabatic (2D variational) solver. (v2.0.8) |

### Level 4 (Mol-Frame Hyperspherical — H₂)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Hyperradial overlap (S) | algebraic | Same Laguerre three-term recurrence as Level 3. Pentadiagonal M2 moment matrix. (Track I, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Same derivative kernel as Level 3. Tridiagonal expansion. (Track I, v2.0.10) |
| Hyperradial potential (U_eff) | numerical-required | Adiabatic curve from angular sweep. μ(ρ) is piecewise-smooth in ρ (NOT algebraic — split-region Legendre expansion creates kinks at min/max boundaries). Structurally different from Level 3. (Track S, v2.0.13) |
| Angular eigenvalue sweep (FD) | numerical-required | Point-by-point diagonalization of H_ang at ~130 ρ-values. FD: 1000×1000 matrices, 99% of wall time. Structurally irreducible: no global P(ρ,μ)=0 exists (Track S). |
| Angular eigenvalue sweep (spectral) | numerical-required | Jacobi polynomial basis: 50×50 matrices, 269× speedup, ~50% of wall time. SO(6) Casimir free spectrum + precomputed V_ee coupling. `angular_method='spectral'`. Optimal approach given Track S's structural finding. (Track K, v2.0.11) |
| 2D tensor product assembly | numerical-required | H_ang evaluated at each quadrature point for (R_e, α) tensor product. Dense kronecker assembly. |
| Nuclear coupling (split-region Legendre) | algebraic | Gaunt integrals, exact via 3j triangle inequality (Paper 15). |
| Cusp factor f(α) | negative result | Alpha-only cusp factor f(α) = 1 + (R_e/2)sin(2α) does not improve D_e. Cusp is 2D in (α, θ₁₂), not separable in α alone. Slow convergence dominated by θ₁₂ Gegenbauer expansion. (Track U, v2.0.14) |
| 1/r₁₂ graph absorption | structural obstruction | Cannot absorb 1/r₁₂ into conformally weighted S⁵ Laplacian. Green's function singularity mismatch: 1/d³ (S⁵) vs 1/d¹ (Coulomb). Cusp is an embedding exchange constant. (Track W, v2.0.14) |
| Cusp energy correction | corrective (post-processing) | Schwartz partial-wave extrapolation: ΔE_cusp ~ -A/(l_max+2)⁴. A = (10/π)⟨δ³(r₁₂)⟩. He: 0.10% at l_max=2 (from 0.24%). H₂: R-dependent, ~1.7 mHa differential, ~1.0 pp D_e. Basis-dependent, not exchange constant. (Track X, v2.0.14) |

### Level 4N (Full N-Electron Mol-Frame Hyperspherical — LiH)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| SO(12) Casimir eigenvalues | algebraic | Exact from SO(3N) representation theory. Free spectrum validated. (Track AK, v2.0.20) |
| S₄ Young projector | algebraic | Character-based projection onto [2,2] singlet irrep. Optimized channel-space eigendecompose. (Track AJ, v2.0.20) |
| Cross-pair V_ee (Gaunt direction) | algebraic | Gaunt integral structure for angular coupling between electron pairs. (Track AK, v2.0.20) |
| Cross-pair V_ee (hyperangular) | numerical-required | 3D numerical hyperangular integration. ρ-independent (precomputed once). (Track AK, v2.0.20) |
| Angular eigenvalue sweep | numerical-required | Point-by-point diagonalization at each ρ. Spectral basis: 750 dim at l_max=2 (1000× compression from FD). (Track AK, v2.0.20) |
| Hyperradial solve | numerical-required | Same Laguerre spectral basis as Level 4. (Track AJ, v2.0.20) |

### Level 5 (Composed Geometries — LiH, BeH₂, H₂O)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Z_eff screening function | algebraic | Laguerre spectral expansion: n(r) projected onto L_k(2βr)·r²·e^{-2βr} basis. N_core(r) = C_inf − e^{-X}·P(X), single transcendental seed e^{-2βr}. Production default `zeff_method='spectral_laguerre'`. Auto-fallback to spline for Z≥4 (polynomial cancellation). (Track N, v2.0.12; Track R wiring, v2.0.13) |
| Slater F^k integrals | algebraic | Two implementations: (1) Laguerre moment decomposition (Track O, v2.0.12): F^k = c_A^T · R^k · c_B, single transcendental seed γ(1,x). (2) Hypergeometric R^k evaluator (v2.8.1): `geovac/hypergeometric_slater.py` with exact Fraction-arithmetic for arbitrary n_max. `compute_rk_float()` gives machine precision (1.5e-12 max error) at 25x speedup over Fraction path. Validated 144/145 table entries in `casimir_ci.py`; found+fixed F²(2p,2p) typo (43/512→45/512). Corrected systematic grid bias (0.06-0.44% per integral). 8x speedup over grid quadrature. Production default for graph-native CI. |
| Phillips-Kleinman PK | algebraic | Ab initio from core eigenvector projection (Paper 17 Sec IV). Gaussian PK parameters A, B from core screening. |
| Inter-fiber exchange coupling | algebraic | Full 1-RDM exchange via channel overlap S(R) and Slater F^0 (algebraic Laguerre). Bond-bond coupling validated (~0.5 Ha). Lone pair coupling disabled at Z_eff>4 (physics limitation). |
| Cross-center V_ne (multipole) | algebraic | Exact termination at L = 2*l_max by Gaunt selection rules. 33 terms at n_max=2, 168 at n_max=3. m-diagonal. Analytical evaluation via incomplete gamma functions (machine precision, zero quadrature). Replaces grid-based trapezoid. (Track CD, v2.0.39-40) |

### Level 5 — Spin-ful composed (Tier 2)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Szmytkowski angular σ·r̂, J², L·S, L², σ² | algebraic | Exact eigenvalues in (κ, m_j); integer / half-integer Kronecker deltas |
| Diagonal hydrogenic ⟨r^k⟩, ⟨1/r^k⟩ | algebraic | Bethe-Salpeter rationals; ⟨1/r⟩=Z/n², ⟨1/r³⟩=Z³/[n³l(l+½)(l+1)] (Z³ diverges l=0, suppressed by Kramers) |
| Off-diagonal ⟨n'l'\|r^k\|n l⟩ | algebraic | Direct sympy integration of assoc_laguerre; ~0.1–1s per call |
| Breit-Pauli H_SO in (κ, m_j) | algebraic | Closed form H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)], diagonal; exact Kramers at l=0 |
| jj-coupled angular coefficient X_k(κ_a, m_a, κ_c, m_c) | algebraic | sympy wigner_3j, cached, full-Gaunt selection |
| Spinor composed two-body ERI | algebraic | X_k·X_k·R^k factorization, radial via hypergeometric_slater (exact Fraction or machine-float) |
| Relativistic γ = √(1−(Zα)²) | algebraic | Degree-2 algebraic over ℚ(α²) via γ² + (Zα)² = 1; reserved symbol in T1/T2, not bound at Tier 2 |
| Dirac-Coulomb radial ⟨r^{-1}⟩ (all states) | algebraic | Hellmann-Feynman: Z(γn_r+κ²)/(γN_D³), exact for all n_r. NR limit Z/n². (T7, v2.12.0) |
| Dirac-Coulomb radial ⟨r^{-2}⟩, ⟨r^{-3}⟩ (n_r=0) | algebraic | Pochhammer ratios from single-term wavefunction structure. (T7, v2.12.0) |
| Dirac-Coulomb radial ⟨r^{-2}⟩ (all states) | algebraic | Kramers-Pasternak direct integration via confluent hypergeometric polynomial expansion. Works for all n_r, all κ. (v2.18.1) |
| Dirac-Coulomb radial ⟨r^{-3}⟩ (|κ|≥2) | algebraic | Kramers-Pasternak direct integration. p₁/₂ (|κ|=1, n_r≥1) limitation: sympy integral does not converge for near-singular integrand. (v2.18.1) |
| Darwin + mass-velocity α⁴ corrections | algebraic | Closed-form rationals in (Z, n, l) times α²; verified Dirac formula exact. (T8, v2.12.0) |
| Breit retarded radial R^k_BP(n₁l₁,n₂l₂;n₃l₃,n₄l₄) | algebraic | Exact Fraction/sympy arithmetic via r_<^k / r_>^{k+3} kernel; Z³ scaling verified; `geovac/breit_integrals.py`. (v2.18.0) |
| Seeley-DeWitt coefficients a₀, a₁, a₂ on S³ | algebraic | Exact sympy Rational on unit S³ (a₀=a₁=√π, a₂=√π/8). `geovac/qed_vacuum_polarization.py`. (v2.18.2) |
| QED vertex even/odd Dirac split D_even(s), D_odd(s) | algebraic | Exact Hurwitz zeta at quarter-integer shifts: D_even via ζ(s,3/4), D_odd via ζ(s,5/4). At s=4: D_even = π²/2 − π⁴/24 − 4G + 4β(4), D_odd = π²/2 − π⁴/24 + 4G − 4β(4). PSLQ-verified at 80 digits. `geovac/qed_vertex.py`. (v2.19.1) |
| Screened ⟨1/r³⟩ from FrozenCore Z_eff(r) | numerical-required | Radial Schrodinger equation with Z_eff(r)/r potential solved on uniform grid (FD, eigh_tridiagonal). ⟨1/r³⟩ = ∫\|u\|²/r³ dr from normalized wavefunction. Enhancement 12-144× over hydrogenic at Z_eff=2 for Ca/Sr/Ba 2p. `geovac/neon_core.py` (screened_r3_inverse, screened_xi_so). (v2.19.4) |

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

3. **PM agent (main Claude Code session, full context):** Reads CLAUDE.md, README, and all papers relevant to the current track. Decomposes tasks into sub-agent work units. Evaluates sub-agent results against verification checklists. Does NOT write production code itself — it plans, delegates, and synthesizes.

4. **Worker sub-agents (scoped context):** Execute specific coding, computation, or drafting tasks. Each receives only the files it needs. Reports results back to the PM agent.

**The research loop:**
```
Leader Brief → PI picks direction → PM decomposes into tracks →
Workers execute → PM verifies → Reviewer critiques papers →
Results feed back to Leader for next brief
```

### 13.2 PM Session Kickoff

Every PM session begins by reading CLAUDE.md and then executing the following:

1. Identify the current track(s) and relevant papers from the plan mode directive
2. Read those papers and any results from the previous session
3. Check the failed approaches summary (Section 3); if the current track touches PK, cusp, inter-group antisymmetry, or molecular encoding, read the full details in CHANGELOG.md before proceeding
4. Decompose the session goal into sub-agent tasks
5. Draft sub-agent prompts using the standard template (13.3)
6. Identify which papers need updating based on the session's results (see 13.8)

**Sprint standard:** The default workflow is one PM prompt per sprint, dispatching all tracks as parallel sub-agents. The PI provides the sprint plan (track definitions, PM prompts, exit criteria) in plan-mode chat. The PM executes without further approval cycles.

### 13.3 Sub-Agent Prompt Template

Every sub-agent dispatch uses this format:

```
CONTEXT FILES: [minimal set of files to read — list explicitly]
TASK: [one clear deliverable, stated in one sentence]
CONSTRAINTS:
  - Failed approaches to avoid: [list relevant entries from Section 3 / CHANGELOG.md]
  - Structures that must be preserved: [quantum numbers, selection rules, etc.]
  - Do NOT modify: [list any files that are off-limits]
SUCCESS CRITERIA:
  - Tests to pass: [specific test files or assertions]
  - Numerical targets: [specific error thresholds if applicable]
  - Consistency check: [which papers or results to verify against]
OUTPUT FORMAT: [what to report back — specific data, not just pass/fail]
PAPER UPDATES:
  - Papers affected: [list paper numbers]
  - Changes to make: [list specific edits — apply directly to papers]
```

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
- Removal of the "conjectural" label from the **combination rule K = π(B + F − Δ) in Paper 2**. (Paper 2 itself was promoted from Conjectures to Core on 2026-04-18; the paper's surrounding theorems — three homes, three obstructions, Sprint A structural verifications — are not conjectural, but the combination rule observation remains conjectural until a first-principles derivation is demonstrated. The conjectural label must stay attached to the observation itself, not the paper's tier.)

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

PMs may edit papers in `papers/core/`, `papers/supporting/`, `papers/observations/`, and `papers/conjectures/` directly for ALL changes, including:

- Adding or updating benchmark tables and numerical results
- Adding new subsections documenting methods and results
- Correcting claims that are contradicted by new computational evidence
- Reframing results based on new findings (e.g., reclassifying transcendental → algebraic when proven)
- Updating abstracts and conclusions to reflect current best results

**PMs must still NOT:**
- Introduce fitted or empirical parameters without PI direction
- Change the natural geometry hierarchy (new levels, changed coordinates)
- Delete or suppress negative results from Section 3
- Remove the "conjectural" label from the **combination rule K = π(B + F − Δ) in Paper 2** (even though Paper 2 is now Core tier; the prohibition narrowed on 2026-04-18 from paper-level to combination-rule-level)

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