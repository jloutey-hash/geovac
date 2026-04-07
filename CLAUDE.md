# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v2.4.0 (April 6, 2026)
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

**Phase 3 (v2.0.24+): Quantum simulation.** The composed architecture's structural sparsity (O(Q^2.5) Pauli scaling, block-diagonal ERIs) is the framework's primary computational advantage. Classical solver results serve as validation benchmarks for the quantum Hamiltonians. Next steps: quantum resource estimation, hardware-aware circuit compilation, experimental collaboration.

---

## 2. Current Development Frontier

**Best results by system type:**
- Atoms: He at 0.05% error (hyperspherical coordinates, Paper 13)
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
- Fine structure constant: alpha from Hopf bundle at 8.8x10^-8, zero free parameters, p-value 5.2x10^-9, universal algebraic identity B_formal/N = d, Hopf generalization negative result, circulant Hermiticity, second selection principle (Paper 2, conjectural)

**Classical solver status: INVESTIGATION COMPLETE (v2.0.24).** The systematic exploration of all classical solver architectures across the natural geometry hierarchy is complete. Over 30 investigation tracks (v2.0.6-23) exhausted all solver (adiabatic, coupled-channel, 2D variational) x PK (channel-blind, l-dependent, PK-free) x basis (FD, spectral, algebraic) combinations. Structural accuracy ceilings are characterized at every level: H2 at 96.0% D_e with CBS ~97% (Level 4), LiH R_eq at 5.3% with structural drift +0.15 bohr/l_max (Level 5 composed), full N-electron equilibrium without PK at 63.5% R_eq error (Level 4N, angular basis limited). The exchange constant taxonomy (Paper 18) classifies where transcendental content enters at each level. Remaining classical accuracy improvements require either higher l_max (diminishing returns, increasing cost) or fundamentally new coordinate systems (none identified). The composed architecture at l_max=2 represents the production operating point for classical PES.

**Quantum computing status: ACTIVE FRONTIER.** The composed architecture produces structurally sparse qubit Hamiltonians: O(Q^2.5) Pauli scaling across LiH/BeH2/H2O (exponent spread 0.02), 51x-1,712x advantage over published Gaussian baselines. Full N-electron encoding comparison (Track AS) confirmed composed is categorically sparser (334 vs 3,288 Pauli terms, 20x lower 1-norm). Equal-qubit He comparison validated against cc-pVDZ/cc-pVTZ computed integrals. Commutator-based Trotter bounds give r ~ Q^{1.47}, 7x fewer steps at Q=60. Head-to-head Gaussian comparison (Track AX, v2.0.26): CLEAR WIN on structural sparsity — 190x fewer Pauli terms for LiH at Q~30 (334 vs ~63,500 cc-pVDZ), 746x for H2O at Q=70; accuracy caveat (GeoVac 5.3%-26% R_eq error vs Gaussian <0.1%). H₂ bond-pair qubit encoding (Track AZ, v2.0.27): single-block composed encoding at Z_eff=1, Q^3.13 scaling (consistent with He atomic), 112 Pauli terms at Q=10, R-independent sparsity confirmed. Ecosystem export pipeline (Track AW): OpenFermion, Qiskit, and PennyLane export via `geovac.ecosystem_export`. `geovac-hamiltonians` PyPI package built (Track BB, v2.0.27): standalone 92 KB wheel, 6 bundled modules, all 5 systems (H2/He/LiH/BeH2/H2O). VQE validation (Track AY): H2 converges to 0.031 mHa on statevector simulator; LiH (30 qubits) infeasible for statevector (17.2 GB RAM). H₂O composed 1-norm (Track BD, v2.0.28): total 28,053 Ha at Q=70, but 98.7% from Z²-scaled PK barrier; electronic-only 1-norm 361 Ha (comparable to BeH₂ 355 Ha). PK diagonal on Z_eff=6 1s orbital is 2,387 Ha — Z² scaling is the bottleneck. PK classical partitioning (Track BF, v2.0.29): PK is a one-body operator whose energy E_PK = Tr(h1_pk · γ) can be computed classically from the VQE 1-RDM with zero additional circuits. Operator decomposition H_full = H_elec + H_pk verified to machine precision (<1e-12) for all three composed systems. Algebraic exactness: E_full = E_elec(ψ) + E_PK(ψ) with residual <1e-13 Ha. Partitioned 1-norms: LiH 33.26 Ha (PK 10.9%), BeH₂ ~355 Ha, H₂O 361 Ha (PK was 98.7% of 28,053 Ha total — 78x reduction). `composed_qubit.py` extended with `pk_in_hamiltonian` kwarg; `pk_partitioning.py` new module; `ecosystem_export.py` updated with partitioned API. 19 tests pass. IBM Quantum VQE demo (Track BC, v2.0.28): `demo/ibm_quantum_demo.py` with Aer simulator and IBM hardware modes; H₂ statevector VQE converges ~13 mHa (10 qubits, 80 params). Positioning documents (Track BE, v2.0.28): `docs/geovac_positioning.md` and `docs/geovac_onepager.md` for outreach. General composed builder and first-row generalization (Tracks BG-BJ, v2.0.30): refactored three hardcoded builders (LiH/BeH₂/H₂O) into single general `build_composed_hamiltonian(spec)` driven by `MolecularSpec` dataclass. Atomic classifier (`geovac/atomic_classifier.py`) maps Z→block decomposition for Z=1-10. Ab initio PK parameters computed for Z=3-9 via Level 3 solver (Track BI): Z² scaling has 5-26% errors, computed values required. Three new molecules: HF (Q=60, 667 Pauli), NH₃ (Q=80, 889 Pauli), CH₄ (Q=90, 1000 Pauli). Within-molecule Pauli scaling exponent ~2.2 (2-point fit max_n=1,2), consistent across all 6 composed molecules. Paper 16 promoted to Core (atomic classifier specification). Scope boundary documented (Track BK): PK s-p splitting is structurally wrong-sign (negative result); first row fully supported, second row feasible with frozen-core tabulation, transition metals out of scope. l_max convergence sprint (Tracks BP-BR, v2.0.32): 2D variational solver wired into composed pipeline (level4_method='variational_2d'), faster than adiabatic at l_max=2 (6.7s vs 9.3s/pt). l_max divergence confirmed structural to PK — NOT from adiabatic approximation. R_eq drifts +0.15-0.22 bohr/l_max with BOTH solvers; single-point E(R=3.015) violates variational bound at l_max≥2. l_max=2 is optimal operating point. Partitioned 1-norms confirmed: LiH electronic 33.26 Ha (PK 10.9%), H₂O electronic 361 Ha (PK 98.7% → 78x reduction via classical partitioning). LiH at max_n=3 (Q=84): 7,879 Pauli terms vs Gaussian cc-pVDZ 63,519 at Q=36 — sparsity advantage grows with basis. Accuracy target analysis (docs/accuracy_target_analysis.md): R_eq error is wrong metric for quantum simulation; resource estimation papers operate at fixed geometries and measure Pauli terms, 1-norm, qubit count. Sturmian CI investigation (Tracks BU-1/BU-2, v2.0.33): Coulomb Sturmian CI improves He FCI by 0.92 pp at max_n=3 (2.06% vs 2.98%). Generalized Sturmian shows crossover (better at max_n=2, worse at max_n=3 -- within-config flexibility loss). Gaunt sparsity exactly preserved in both variants. Qubit encoding: 7% fewer Pauli terms (112 vs 120 at Q=10) but 2.8-4.5x higher 1-norm from Lowdin orthogonalization. Negative result for molecular PK bypass. TC scoping (Track BX-1, v2.0.33): Transcorrelated method compatible with GeoVac angular framework. J = -(1/2)r12 preserves Gaunt selection rules. BCH terminates at 2nd order for 2 electrons. Non-Hermitian eigenproblem required. Feasibility: REQUIRES NEW INTEGRALS. TC in second quantization (Track BX-3, v2.0.35): TC-modified qubit Hamiltonians implemented via `geovac/tc_integrals.py`. He benchmark: accuracy improved from 5.3-8.2% (standard, diverging) to 3.3-3.6% (TC, converging) across max_n=1-3. Composed molecules: Pauli ratio exactly 1.68× (constant factor from non-Hermiticity), O(Q^2.5) scaling preserved. Electronic-only 1-norm overhead 9-16% (LiH/BeH₂), vanishing with PK dominance (H₂O: 1.00×). BCH constant shift is -1/4 per electron pair (sign error found and corrected). Current implementation: radial gradient only; angular gradient for p/d orbitals pending. Key files: `geovac/tc_integrals.py` (compute_tc_integrals_block, build_tc_composed_hamiltonian), `debug/data/tc_he_qubit_benchmark.json`, `debug/data/tc_composed_benchmark.json`. Quantum resource market test (Track CA, v2.0.36): Head-to-head 1-norm comparison against Gaussian raw JW baselines. LiH STO-3G computed from OpenFermion cached integrals: Q=12, 907 Pauli terms, λ=34.3 Ha, 273 QWC groups. GeoVac LiH composed (electronic-only): Q=30, 334 Pauli, λ=33.3 Ha, 21 QWC groups. Result: 2.7× fewer Pauli terms, 13× fewer QWC groups, 0.97× 1-norm (essentially identical). He equal-qubit at Q=28: GeoVac λ=78.4 vs Gaussian cc-pVTZ λ=530.5 (6.8× lower). Literature survey: published DF/THC/SCDF lambda values for molecules at Q<100 DO NOT EXIST — QPE papers benchmark at FeMoco scale (152+ qubits). GeoVac's positioning is strongest for VQE/NISQ regime where Pauli count and QWC groups dominate runtime. Key files: `docs/market_test_results.md`, `benchmarks/gaussian_baseline_comparison.py`, `debug/data/market_test_data.json`. Balanced coupled composition (Track CD, v2.0.39): cross-center V_ne via multipole expansion of 1/|r-R_B|, terminates exactly at L_max=2*l_max by Gaunt selection rules. LiH at n_max=2: 878 Pauli terms (2.63× composed, Conjecture 1 confirmed), 74.1 Ha 1-norm (1.98×), 87 QWC groups. 4-electron FCI: 1.8% energy error, 7.0% R_eq error, only bound configuration (decoupled 10.9%, PK 15.0%, CB 29.0% — all unbound). D_e=0.037 Ha. Key files: `geovac/shibuya_wulfman.py` (compute_cross_center_vne), `geovac/balanced_coupled.py` (build_balanced_hamiltonian). Phase 4A with analytical V_ne integrals (incomplete gamma, machine precision): n_max=3 E(3.015)=-8.055 Ha (0.20% error). R_eq=3.280 bohr (8.8% error, structural drift +0.053 bohr/n_max, 3x smaller than PK). MIXED: energy converges excellently, R_eq drifts structurally. Optimal for single-point quantum simulation at fixed geometries.

**Scope boundary:** See `SCOPE_BOUNDARY.md` at project root for which atoms and molecules are supported, which are feasible, and which are out of scope. First row (Z=1-10) fully supported; second-row atoms (Z=11-18) supported via [Ne] frozen cores; third-row s-block (Z=19-20) via [Ar] frozen cores; third-row p-block (Z=31-36) via [Ar]3d¹⁰ frozen cores. Transition metals Z=21-22 scoping demonstrated (Track CZ/DA, v2.4.0): ScH and TiH Hamiltonians built with d-orbital blocks (l_min=2), d-block Pauli/Q = 9.27 (lower than main-group 11.11). Z=23-30 architecturally feasible but not yet tested. Multi-center molecules (Track CU, v2.3.0): 8 multi-center systems. Total library: 30 molecules (20 single-center + 8 multi-center + 2 transition metal hydrides).

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
- **α combination rule derivation (Paper 2 → Paper 18)** — SHELVED. Paper 18's exchange constant framework suggests K = π(B + F − Δ) may be derivable. Research question: does spectral geometry of fiber bundles constrain K? Not yet investigated. Upstream of Paper 2, downstream of Paper 18.
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
| TC Jastrow in adiabatic hyperspherical solver | 1 | TC replaces multiplicative V_ee (O(R) in angular eigenvalue) with first-derivative G_ang (O(1)); adiabatic separation requires V_ee as multiplicative operator; ~46% error; TC needs direct variational/FCI framework, not adiabatic. **NOTE: TC in second quantization (Track BX-3) succeeds** — the qubit pipeline doesn't use the adiabatic solver |
| TC angular gradient for l>0 orbitals | 1 | Angular gradient adds 2.66× Pauli terms (He max_n=2: 188→500) for 0.01 pp accuracy improvement (3.621%→3.611%). Composed: 2.67× increase (LiH 562→1498). Cost/benefit >100×. Gaunt selection rules preserved but coupling density defeats sparsity advantage. Radial-only is optimal. |
| Coupled composition (cross-block ERIs replacing PK) | 1 | Cross-block ERIs add 2.56× Pauli terms and 2.30× 1-norm for LiH. 4-electron FCI: 29% error (worst of three configurations). Root cause: composed orbital basis (different Z per block) means cross-block ERIs add V_ee repulsion without balancing V_ne cross-center attraction. PK is NOT redundant in second quantization for the composed basis. Two-center h1 terms required but unavailable in single-center framework. |


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
| 3 | He (1-center, 2e) | Hyperspherical | 0.05% | 13 |

*Level 3 note:* The 0.05% result uses the original FD adiabatic solver at l_max=0. The algebraic spectral basis (v2.0.6) gives 0.16% single-channel at l_max=0 (removing lucky FD error cancellation). However, the algebraic solver enables coupled-channel integration, which shows convergent behavior: 0.37% → 0.27% at l_max 1-3, compared to divergent single-channel (0.56% → 0.65%). The coupled-channel solver with exact Q-matrix converges algebraically (~l_max⁻²) to a structural floor of 0.19-0.20% (v2.0.8, l_max=0-5 study). Sub-0.1% requires the variational 2D solver (Paper 15) to bypass the adiabatic approximation.
| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 96.0% D_e | 15 |
| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | R_eq 63.5% (l_max=2, 2D variational; unbound D_e) | 17 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq 5.3% | 17 |
| 5 | BeH₂ (polyatomic) | Composed (Level 3 + 4) + exchange | R_eq 11.7% | 17 |
| 5 | H₂O (triatomic) | Composed (Level 3 + 4) + lone pairs | R_eq 26% | 17 |

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
| Papers 8-9 | Archive | Bond sphere (positive), Sturmian structural theorem (negative), multi-electron Sturmian CI negative result (v2.0.33); historical exploration |
| Paper 10 | Archive | Nuclear lattice draft; rovibrational spectra |
| FCI-M | Archive | LCAO molecular FCI scaffold; superseded by natural geometry |
| Paper 19 | On-topic | Coupled composition, PK alternatives, two-center integrals, Shibuya-Wulfman connection |
| Paper 20 | On-topic | Applications paper: resource benchmarks, FCI accuracy, Gaussian comparison, installation |

### Paper Inventory

#### Core (`papers/core/`) — Always load

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Paper 0 | `Paper_0_Geometric_Packing.tex` | Universal constant K = -1/16 |
| Paper 1 | `paper_1_spectrum.tex` | Spectral graph theory, O(N) eigenvalue methods. Berry phase retracted v1.2.0 |
| Paper 7 | `Paper_7_Dimensionless_Vacuum.tex` | S3 proof (18/18 symbolic), Schrodinger recovery, SO(3N) generalization |
| Paper 14 | `paper_14_qubit_encoding.tex` | Structurally sparse qubits: O(Q^3.15) atoms, O(Q^2.5) composed; Trenev et al. Gaussian baselines |
| Paper 16 | `paper_16_periodicity.tex` | Chemical periodicity from S_N representation theory; atomic classifier for composed geometry |

#### Supporting (`papers/supporting/`) — Load on-topic

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

#### Archive (`papers/archive/`) — Load only on explicit request

| Paper | File | Key Result |
|:------|:-----|:-----------|
| Papers 8-9 | `Paper_8_Bond_Sphere_Sturmian.tex` | Bond sphere (positive), Sturmian structural theorem (negative), SO(4) selection rules, multi-electron Sturmian CI negative result (v2.0.33) |
| Paper 10 | `paper_10_nuclear_lattice.tex` | Rovibrational spectra from SU(2) algebraic chains |
| FCI-M | `paper_fci_molecules.tex` | LCAO LiH results and diagnostic arc |

#### Observations (`papers/observations/`)

(No active observation papers — Paper 16 promoted to Core.)

#### Conjectures (`papers/conjectures/`)

| Paper | File | Key Topic |
|:------|:-----|:----------|
| Paper 2 | `paper_2_alpha.tex` | Fine structure constant from Hopf bundle (8.8x10^-8, p = 5.2x10^-9, circulant Hermiticity, second selection, universal B_formal/N = d identity, Hopf generalization negative result, zero params) |
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
| LiF spec (multi-center) | `geovac/composed_qubit.py` | `lif_spec()` |
| CO spec (multi-center) | `geovac/composed_qubit.py` | `co_spec()` |
| N₂ spec (multi-center) | `geovac/composed_qubit.py` | `n2_spec()` |
| F₂ spec (multi-center) | `geovac/composed_qubit.py` | `f2_spec()` |
| NaCl spec (multi-center) | `geovac/composed_qubit.py` | `nacl_spec()` |
| CH₂O spec (multi-center) | `geovac/composed_qubit.py` | `ch2o_spec()` |
| C₂H₂ spec (multi-center) | `geovac/composed_qubit.py` | `c2h2_spec()` |
| C₂H₆ spec (multi-center) | `geovac/composed_qubit.py` | `c2h6_spec()` |
| Multi-center nuclear repulsion | `geovac/composed_qubit.py` | `_compute_multi_center_nuclear_repulsion()` |
| ScH spec (transition metal) | `geovac/composed_qubit.py` | `sch_spec()` |
| TiH spec (transition metal) | `geovac/composed_qubit.py` | `tih_spec()` |
| General-l Wigner D rotation | `geovac/shibuya_wulfman.py` | `_build_rotation_matrix_real_sh(l, theta, phi)`, `_build_block_rotation_matrix()` |
| Physical constants | `geovac/constants.py` | `HBAR`, `C`, `ALPHA`, etc. |

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
| ScH composed Pauli terms (Q=30) | exactly 278 | Transition metal qubit validation |
| TiH composed Pauli terms (Q=30) | exactly 278 | Transition metal qubit validation |
| d-only block ERI density | 4.0% | d-orbital Gaunt sparsity |
| d-block Pauli/Q | 5.60 (< s+p 11.20) | d-orbital sparsity advantage |
| LiH/BeH₂/H₂O Pauli unchanged | 334/556/778 | Backward compatibility regression |
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
| Slater F^k integrals | algebraic | Laguerre moment decomposition: F^k = c_A^T · R^k · c_B where R^k matrix uses incomplete gamma recurrence. Single transcendental seed γ(1,x) = 1−e^{-x}. Production default `slater_method='algebraic_laguerre'`. Agreement within 1.77% of numerical for compact Z≥2 densities. (Track O, v2.0.12; Track R wiring, v2.0.13) |
| Phillips-Kleinman PK | algebraic | Ab initio from core eigenvector projection (Paper 17 Sec IV). Gaussian PK parameters A, B from core screening. |
| Inter-fiber exchange coupling | algebraic | Full 1-RDM exchange via channel overlap S(R) and Slater F^0 (algebraic Laguerre). Bond-bond coupling validated (~0.5 Ha). Lone pair coupling disabled at Z_eff>4 (physics limitation). |
| Cross-center V_ne (multipole) | algebraic | Exact termination at L = 2*l_max by Gaunt selection rules. 33 terms at n_max=2, 168 at n_max=3. m-diagonal. Analytical evaluation via incomplete gamma functions (machine precision, zero quadrature). Replaces grid-based trapezoid. (Track CD, v2.0.39-40) |

---


## 13. Multi-Agent Protocol

The GeoVac project uses an AI-augmented agentic workflow with a formalized planner–worker architecture for Claude Code sessions.

### 13.1 Architecture

Three layers:

1. **Plan mode (human + chat):** Strategic direction, framing decisions, result synthesis. CLAUDE.md changes originate here. This layer is not automated.
2. **PM agent (main Claude Code session, full context):** Reads CLAUDE.md, README, and all papers relevant to the current track. Decomposes tasks into sub-agent work units. Evaluates sub-agent results against verification checklists. Does NOT write production code itself — it plans, delegates, and synthesizes.
3. **Worker sub-agents (scoped context):** Execute specific coding, computation, or drafting tasks. Each receives only the files it needs. Reports results back to the PM agent.

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

### 13.5 Hard Prohibitions

The following changes must NEVER be made by sub-agents or the PM agent:

- Any change to the natural geometry hierarchy (new levels, changed coordinate systems)
- Introduction of any fitted or empirical parameter
- Deletion or suppression of negative results from Section 3 or CHANGELOG.md
- Removal of the "conjectural" label from Paper 2

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
- Remove the "conjectural" label from Paper 2

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