# Changelog

All notable changes to GeoVac will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Note:** the CHANGELOG is currently behind the `CLAUDE.md` version cursor (CLAUDE.md tracks v2.10–v2.26 range; CHANGELOG below jumps from v2.9.2 to v2.26.1). Intermediate version entries for the RH sprint series (v2.20–v2.25, Papers 28/29/30) are in `git log` commit messages but have not yet been back-filled into this CHANGELOG. A consolidation sprint is flagged for future work.

## [2.26.1] - 2026-04-18

### Changed — Paper 2 promoted from Conjectures to Core

- **File move:** `papers/conjectures/paper_2_alpha.tex` → `papers/core/paper_2_alpha.tex` (via `git mv`; history preserved).
- **Rationale:** Sprint A (2026-04-18) applied the structural reframe (see v2.26.0 notes + `debug/paper_2_reframe_skeleton.md` + `debug/alpha_sprint_a/`). Three structural verifications landed with partial positive / partial negative verdicts (α-PI, α-MI, α-SP, α-EB v2, α-X, α-LS); Marcolli-van Suijlekom 2014 gauge-network framework identified as published lineage (WH1 upgraded to MODERATE); Paper 2 body extended by 7 substantive edits (+168 lines) naming the triple m=3 selection, Hopf-measure π, APS-shape on Δ minus sign, π³α³ residual footnote (as hint, not derivation), and spectral-triple setting with Marcolli-vS citation.
- **Conjectural-label prohibition narrowed** in `CLAUDE.md` §13.5 and §13.8: "Removal of the 'conjectural' label from Paper 2" → "Removal of the 'conjectural' label from the combination rule K = π(B + F − Δ) in Paper 2." The paper's surrounding theorems (three-homes, three obstructions, Sprint A verifications) are not conjectural; the combination rule observation remains conjectural until a first-principles derivation lands.
- **Papers 25 and 30 updated** (same date) with Marcolli-van Suijlekom 2014 + Perez-Sanchez 2024/2025 citations establishing published gauge-network lineage.
- **CLAUDE.md §6 Paper Inventory** updated: Paper 2 added to Core tier, removed from Conjectures tier. Context Loading Guide entry added (On-topic).
- **README.md Paper Series** table updated: Paper 2 bolded and description reflects reframe.

### Notes

- Sprint B (g−2 from QED-on-S³) plan drafted as `debug/sprint_b_g_minus_2_plan.md`; NOT auto-dispatched; awaits PI go/no-go. Sprint B will test whether the framework predicts the anomalous magnetic moment at one loop — the Layer 2 out-of-sample prediction test per the promotion criterion.

## [2.26.0] - 2026-04-18

### Added — CLAUDE.md §1.7 Working Hypotheses (Internal Register)

- **New section:** CLAUDE.md §1.7 "Working Hypotheses (Internal Register)" — bold-claim register distinct from paper rhetoric. Six WHs (almost-commutative spectral triple, Paper 18 as Seeley-DeWitt + ζ-invariant decomposition, lattice a priori from packing origin, four-way S³ coincidence as one structure, α as projection constant, D(s) not classical Riemann ζ), each with falsifier and status fields. Governance specifies PM may update Status only; add/retire/promote requires PI direction.
- **Purpose:** papers remain cautious under §1.5 rhetoric rule; internal thinking gets Ramanujan-register license. Papers are Hardy-letters; the WH register is the notebook.
- **Memory note:** `memory/wh_register_april2026.md` records the register's addition and conversational context.

## [2.9.2] - 2026-04-15

### Added — Paper 27 (entropy as projection artifact)

- **New core paper:** `papers/core/paper_27_entropy_projection.tex` ("Entropy as a Projection Artifact: One-Body Operators are Entanglement-Inert on Sparse Lattices"). Five results:
  1. Operator-theoretic floor (theorem): non-degenerate one-body GS has S_kin=0 in its natural-orbital basis. Verified S_kin/S_full ~ 1e-14 for He at n_max=2,3 (EP-1).
  2. Area-law identity: A_n = g_n² = 4n⁴ pair count corrects Paper 5's one-particle mislabeling.
  3. Cusp localization: (1s,1s) is the unique hot-node on the V_ee pair graph.
  4. **HO zero-entropy rigidity (theorem, EP-2b):** on Bargmann-Segal HO basis, [H_HO, V_central]=0 exactly via Moshinsky-Talmi total-quanta conservation. Closed-shell 2-fermion GS is a single Slater determinant for ANY central two-body V. S=0 identically. New module `geovac/nuclear/ho_two_fermion.py`.
  5. **Universal scaling (EP-2c→2N):** S_B = A(w̃_B/δ_B)^γ with γ_∞ ≈ 1.96 (Richardson on n_max=2..5). Below second-order Rayleigh-Schrödinger γ=2; residue attributed to multi-shell aggregation. Verified across He-like Z∈[2,100], LiH R-sweep, HF/H₂O R-sweeps, Be analytical degenerate-PT.
- **Paper 24 extension:** §V.C HO entanglement rigidity corollary — structural dual of Fock rigidity, completing the π-free partition.
- **Paper 26 cross-reference:** S~Z^{-2.56} (Paper 26 §III) derives from Paper 27 §VI.B's universal curve. Reproduced at n_max=4: measured α_Z=-2.546, R²=0.995.
- **CLAUDE.md Section 2** consolidated with the full EP-1 → EP-2N sprint arc.

### Added — Cusp attack sprints

- **CUSP-1** (`debug/cusp1_screened_ci.py`): w̃/δ-based 2nd-order-RS configuration screening for He graph-native CI at n_max=7. **11× config compression** (1218→109) at the full-CI floor, 22× (→56) at floor+0.01pp. Top-scored states are exactly the cusp-relevant (1s,ns), (2s,2s), (2p_-1, 2p_+1) singlet configs. Useful for VQE/qubit compactness; doesn't break the accuracy floor (which is structural).
- **CUSP-2** (`debug/cusp2_hotnode_patch.py`): Tested Schwartz-tail patch on the (1s,1s) hot-node. **Decisive negative diagnosis of a long-standing misunderstanding**: the 0.20% floor for graph-native CI on He at Z=2 is NOT the cusp. Linear fit err_abs = -A/(l+2)⁴ + c gives irreducible c ≈ 6 mHa — the small-Z graph-validity-boundary artifact (Z_c ≈ 1.84). Z=10 confirms: error sign flips (under-binding +85 mHa). **For He graph-native CI, the cusp contribution is below 1 mHa at l_max≥4; the reported ceiling is basis-internal.** Real cusp work should be conducted at Z≥4 where the small-Z artifact is suppressed.
- **CUSP-3** (`debug/cusp3_tc_high_lmax.py`): TC benchmarked with proper particle-number-projected FCI at He composed n_max=2,3,4,5. **Decisive structural negative**: TC error plateaus at 3.40→3.47%, ratios approaching 1 (0.987→0.994→0.997). Standard continues to improve 2.48→2.02 and is still decreasing. Pauli ratio TC/Std worsens 1.68→2.73. **TC is now confirmed dead at every accessible n_max in second quantization on the composed basis**, not just small n_max. Section 3 updated with this extended-n_max evidence.

### Added — Paper 27 EP-2 sprint chain (drivers, data, tests)

- Drivers: `debug/ep2{b,c,d,e,f,g,h,i,j,k,L,M,N}_*.py` (13 scripts, one per sprint).
- Data: `debug/data/ep2{b,c,c_multi_block,d,e,f,g,h,i,j,k,L,M,N}_*.json`.
- Tests (`tests/test_paper27_entropy.py`): 27 tests total (26 passing + 1 `--slow`-gated). Includes:
  - EP-1 reproduction at n_max=2,3
  - Area-law pair counting (combinatorial + leading-order + one-body counter-check)
  - Energy-graph structural invariants at n_max=3,4 (commutator ratio, diagonality, (1s,1s)=5/8)
  - EP-2b HO commutator + single-SD ground state
  - EP-2c/2d/2g/2h/2i/2L dimensionless power law + Z-range artifact + asymptote below 2
  - EP-2n Be analytical degenerate-PT
  - Paper 26 / Paper 27 Z-scaling consistency
  - Proposition non-degeneracy qualifier (Li/Be)

### Fixed

- **Paper 27 abstract + conclusion:** γ → 2 asymptote corrected to γ_∞ ≈ 1.96 based on n_max=5 Richardson extrapolation. Previous "γ → 2 from above" claim was an artifact of limited n_max coverage.
- **Paper 27 Section II.C (Scope of the proposition):** added non-degeneracy qualifier, open-shell qualifier, and Be analytical degenerate-PT subsection.
- **Paper 27 Section V (wording fix):** ‖V_ee^diag(H_1)‖_F / ‖V_ee‖_F is unsquared Frobenius ratio; previous text incorrectly stated squared formula. n_max=4 value 0.94 corrected to 0.892.
- **Paper 27 Section VII consolidated:** 1111→900 lines, sprint chronology collapsed into two clean subsections (HO rigidity + Universal scaling). Methods/robustness collected into one paragraph.

### Status

- Paper 27 core content complete and regression-locked.
- Cusp problem: fundamentally re-diagnosed at He Z=2 (not cusp, graph-validity). No breakthrough on the cusp itself at Z≥4, but three structural facts established (see Section 3 updates).

---

## [2.9.1] - 2026-04-14

### Added (documentation only — no production code changes)

- **Energy graph characterization for V_ee on S³ (exploratory sprint):** Tested the hypothesis that a "second graph" — nodes = electron pair-states |(n₁l₁m₁)(n₂l₂m₂)⟩, edges = ⟨ab|1/r₁₂|cd⟩ — could be a Paper-12-style algebraic generator for V_ee on S³ (analog of the Neumann expansion for H₂ in prolate spheroidal coordinates).
- Findings:
  1. Pair-state graph at n_max=3 has 31 nodes, n_max=4 has 101 nodes (singlet M_L=0 sector).
  2. Within-parity blocks are ~47% dense at both sizes — orbital-level Gaunt sparsity (Paper 22) does NOT survive projection to pair-states.
  3. Diagonal V_ii are exact rationals (⟨1s1s|V|1s1s⟩ = 5/8, ⟨2s2s|V|2s2s⟩ = 77/512), but cross-shell entries mix 2^k and 3^j denominators — V's spectrum is non-rational. Consistent with Paper 18's classification of 1/r₁₂ as embedding exchange constant.
  4. No three-term recurrence in (n,l) for the Slater integrals: V(1s·ns) ratios drift {2.70, 2.20}, V(ns·ns) drift {4.16, 2.26}. Neumann-style separability is specific to prolate spheroidal coordinates, not S³.
  5. **Main new finding:** In the H₁ eigenbasis, 92% of V's Frobenius mass is diagonal at n_max=3, 94% at n_max=4 (residual ‖[H₁,V]‖/‖V‖ = 6.1% → 5.3%, saturating not vanishing). The wavefunction graph nearly diagonalizes V_ee.
  6. **Cusp signature:** The (1s1s) pair-state is simultaneously the highest-diagonal node, largest weighted hub, and head of the hottest off-diagonal edge. The cusp is concentrated on one pair-state at coalescence, NOT smeared across all angular channels.
- Assessment: NEGATIVE on the Paper-12 analog hypothesis (no clean algebraic generator for V_ee on S³); POSITIVE PARTIAL on wavefunction-graph near-diagonalization of V_ee. The "energy graph as distinct object" framing is less useful than hoped — the two graphs share most eigenstructure.
- Documentation: Paper 18 updated with a new paragraph in the embedding-constants section making the graph-theoretic characterization of 1/r₁₂ explicit. CLAUDE.md Section 3 (failed approaches) records the Paper-12 analog negative result; Section 2 cusp paragraph extended with the saturation finding; Section 11 topic lookup updated.
- Pointers: `debug/energy_graph_exploration.md`, `debug/data/energy_graph_nmax3.json`, `debug/data/energy_graph_nmax4.json`.
- No production code changed; no new tests added.

---

## [2.9.0] - 2026-04-13

### Added

- **111 Pauli count derivation (Investigation 1):** Analytically derived N_Pauli = 55 + 56 per s/p block at max_n=2. 55 = Q(Q+1)/2 direct terms (universal, from k=0 monopole). 56 exchange terms from 3 Gaunt channels: s-s (16), s-p cross-Coulomb (24), k=1 dipole (16). Pure l-shells (e.g., d-only) have zero exchange → Pauli/Q = (Q+1)/2 = 5.5. Universal coefficient 11.1 = 111/10. Key files: `debug/` analysis scripts.
- **He energy decomposition (Investigation 2):** `<h1>` converges by n_max=2 (within 0.04 Ha of limit). `<V_ee>` is the slow component — ALL remaining convergence struggle is off-diagonal V_ee correlation. `<V_ee>/|E|` ratio: 0.4545 (n_max=1) → 0.3290 (n_max=7) → 0.3257 (exact). Finite basis overestimates electron coalescence probability.
- **V_ee spectral analysis (Investigation 5):** V_ee is FULL RANK in the graph eigenbasis at every n_max (2-5). No low-rank shortcut exists for the cusp. Diagonal V_ee (mean-field) converges instantly (locked by n_max=3 at -2.787 Ha). Off-diagonal V_ee correlation converges as ~n_max^(-1) toward ~0.109 Ha. Cusp distributes broadly across all angular channels.
- **Graph-native CI variational boundary (Z-sweep):** The graph Laplacian CI crosses from non-variational (over-binding) to variational at Z_c ≈ 1.84 (n_max=7, extrapolated CBS ~1.87-1.89). Below Z_c, kappa=-1/16 off-diagonal coupling dominates Z²-scaled diagonal, creating artificial stabilization. Relative importance scales as 1/(8Z²): 12.5% at Z=1, 3.1% at Z=2, 0.13% at Z=10.
- **H⁻ variational-k investigation:** Standard (non-graph) FCI is ALWAYS variational for H⁻. Graph-native FCI over-binds by 21% (E=-0.640 vs exact -0.528 Ha). Variational-k optimization makes it WORSE (k_opt=1.16, energy drops further). Over-binding is specific to graph h1 topology, not orbital exponent.
- **TC gamma optimal scan (Investigation 4):** gamma_opt decreases with l_max, approximately gamma_opt ~ 0.51/(l_max+1.7). l_max=3 sweet spot (0.001% error, 33x improvement) is a near-cancellation, not systematic. TC roughly doubles convergence exponent (~l⁻¹ to ~l⁻²) but does not achieve theoretical l⁻⁸. Large-gamma floor ~0.35% from non-Hermitian contamination.
- **PsH exotic atom solver:** Standalone positronium hydride (e⁻e⁺ in proton field) solver using Level 3 hyperspherical framework with sign-flipped charge function. PsH IS bound: V_eff minimum at R=3.62 bohr, well depth 0.042 Ha below H+e⁺ threshold (l_max=3). Energy -0.756 Ha (4.1% error vs exact -0.789 Ha). Alpha parity mixing essential (distinguishable particles). Gaunt selection rules preserved.
- **H⁻ binding via graph-native CI:** H⁻ found bound at n_max≥2 but with massive over-binding (21% error). Confirms correlation required for H⁻ (n_max=1 correctly unbound).
- Analysis scripts in `debug/`: `he_energy_decomposition.py`, `vee_spectral_analysis.py`, `z_sweep_variational.py`, `h_minus_variational_k.py`, `tc_gamma_scan.py`, `psh_solver.py`, `positronium_analysis.py`
- Plots in `debug/plots/`: `vee_graph_eigenbasis.png`, `vee_sv_decay.png`, `psh_adiabatic_curves.png`, `psh_mu_comparison.png`, `psh_charge_function.png`

### Key Scientific Findings

- **Cusp characterization:** The electron-electron cusp is a full-rank, broadly-distributed, off-diagonal V_ee object in the graph eigenbasis. h1 and diagonal V_ee converge fast; ALL convergence struggle is off-diagonal V_ee. No rank-reduction or channel-dominant shortcut exists. TC similarity transformation is the structurally correct response (modifies effective V_ee).
- **Graph validity boundary:** Z_c ≈ 1.84. The graph Laplacian with kappa=-1/16 is reliable when Z² >> kappa (nuclear potential dominates graph topology). He (Z=2) barely qualifies. H⁻ (Z=1) does not. The graph's rigid inter-shell coupling overestimates correlation for asymmetric/weakly-bound systems.
- **Exotic atoms:** PsH demonstrates the hyperspherical framework extends to matter-antimatter systems via sign flips. The composed geometry insight applies: systems with disparate particle-nucleus distance scales need composed blocks, not single-center coordinates.

---

## [2.8.2] - 2026-04-13

### Added

- **5 multi-center diatomic molecules:** LiF (Q=70), CO (Q=100), N₂ (Q=100), F₂ (Q=100), NaCl (Q=50)

### Fixed

- Paper 14 corrections: composed coefficient 11.11→11.10 (exact), TM Pauli/Q 9.27→9.23, 1-norm table corrected to electronic-only
- Paper 20 corrections: same coefficient fixes, table caption clarifies composed vs balanced per row
- Balanced coupled fix for `spec.nuclei` attribute

---

## [2.8.1] - 2026-04-12

### Added

- **Algebraic Slater integrals:** `geovac/hypergeometric_slater.py` with exact Fraction-arithmetic R^k evaluator for arbitrary n_max. Validated 144/145 table entries, found+fixed F²(2p,2p) typo (43/512→45/512). 8x speedup.
- **Float algebraic path:** `compute_rk_float()` gives machine-precision (1.5e-12) at 25x faster than Fraction. Corrected systematic grid bias (0.06-0.44% per integral).
- **DUCC downfolding:** `geovac/downfolding.py` computes exact (2J-K) core potential. Root cause of l_max divergence: PK underestimates p-orbital potential by 109x. H₂O 1-norm 9% lower with downfolding.
- **Ecosystem export:** 30→35 molecules via `hamiltonian()` API
- **He graph-native FCI convergence** to n_max=8 (0.207%, 2262 configs) and n_max=9 (0.201%, 3927 configs) with exact algebraic integrals
- **He 2D variational best:** 0.019% self-consistent cusp correction, 0.004% with exact coalescence density
- **Wigner 3j caching, property caching, double-build elimination**
- 418 files restored from OneDrive migration
- 7 broken test imports fixed, 2908 tests collect cleanly

### Fixed

- O (Z=8) PK parameter gap fixed
- `casimir_ci.py` F²(2p,2p) typo corrected (43/512→45/512)
- Pip reinstalled to correct directory

---

## [2.8.0] - 2026-04-12

### Added

- **Full first transition series (Z=21-30)** as hydrides: ScH, TiH, VH, CrH, MnH, FeH, CoH, NiH, CuH, ZnH
- **General `build_composed_hamiltonian(spec)`** in `composed_qubit.py` — MolecularSpec-driven builder consumed by balanced_coupled and coupled_composition
- **Atomic classifier extended to Z=1-30** — second row (Z=11-18), K/Ca (Z=19-20), and all first-row transition metals with structure type F
- **`l_min` field on `OrbitalBlock`** — restricts angular momentum enumeration for d-only blocks (l_min=2)
- **`_v_cross_nuc_frozen_core`** — frozen-core electrostatic potential for multi-shell cores
- **`transition_metal_hydride_spec(Z)`** — spec factory for all 10 TM hydrides with convenience aliases
- **10 TM hydrides in ecosystem export** — accessible via `hamiltonian('ScH')` etc.
- **48 new tests** in `tests/test_transition_metals.py`
- **Cr/Cu anomalous configurations** correctly handled (3d⁵4s¹ and 3d¹⁰4s¹)

### Key Results

- All 10 TM hydrides: Q=30 qubits, 277 Pauli terms, Pauli/Q = 9.23
- Isostructural invariance confirmed: identical block topology → identical Pauli count
- d-block ERI density sparser than s/p (confirming Track CZ: 4.0% vs 8.9%)
- Pauli/Q = 9.23 < 11.11 main-group coefficient — transition metals are cheaper per qubit
- Library expanded from 30 to 38 molecules

### Changed

- `SCOPE_BOUNDARY.md` updated to v2.8.0 — transition metals now "Fully Implemented"
- Paper 20 (Resource Benchmarks) updated with full 10-molecule TM hydride table
- Paper 20 future directions: TM classifier item (iv) removed (completed)

---

## [2.7.1] - 2026-04-12

### Changed

- Outreach-ready documentation correction (tag only)

---

## [2.7.0] - 2026-04-12

### Added

- Papers 22 (Angular Sparsity Theorem), 23 (Nuclear Shell Hamiltonians), 24 (Bargmann-Segal Lattice)
- Paper 21 (Geometric Vacuum Synthesis)
- Precision He: 2D variational solver (0.004%), graph-native CI (0.19%), excited states
- Nuclear shell model: deuteron (16Q/592 Pauli), He-4 (16Q/712 Pauli), composed nuclear-electronic (26Q/614 Pauli)
- `geovac/casimir_ci.py`, `geovac/level3_variational.py`, `geovac/nuclear/` package
- Alpha structural decomposition phases 4B-4H (paused — combination rule open)
- Papers 8-9 promoted from archive to methods tier

---

## [0.4.0] - 2026-02-15

### 🌟 Major Scientific Breakthrough

#### Global Metric Scaling for Isoelectronic Series
- **BREAKTHROUGH:** Conformal transformation approach for multi-electron Z-scaling
- **PHYSICS FIX:** Resolved virial mismatch from previous Jacobian scaling
- **METHOD:** Solve Helium-equivalent system, scale eigenvalues by γ = (Z/2)²
- **VALIDATION:** Li+ 10.87% error, Be2+ 15.22% error (improved from 31.6%/44.5%)

#### Theoretical Significance
- **CONFORMAL INVARIANCE:** Z-scaling is a metric transformation, not parameter change
- **VIRIAL THEOREM:** Both T and V scale uniformly by Z², preserving <T> = -<V>/2
- **UNIVERSALITY:** Lattice topology is universal, only metric (energy scale) changes with Z
- **PHYSICAL LIMIT:** Remaining 10-15% error attributed to relativistic corrections (Z⁴)

### Added

- **Global metric scaling implementation** in isoelectronic tests
- **`docs/GLOBAL_METRIC_SCALING_SUCCESS.md`** - Complete technical analysis
- **`docs/JACOBIAN_SCALING_RESULTS.md`** - Historical context (archived)
- **`debug/plots/create_isoelectronic_plot.py`** - Visualization script
- **`debug/plots/isoelectronic_scaling.png`** - Validation plot
- **`tests/test_isoelectronic.py`** - Comprehensive isoelectronic test suite
- **E/Z² ratio analysis** - Validates near-constant scaling
- **Transition state test** - Linear H3 (19.94% error)

### Changed

- **Version:** 0.3.2 → **0.4.0**
- **README:** Added v0.4.0 section with global metric scaling results
- **README:** Updated benchmarks table with isoelectronic series
- **README:** Updated roadmap (v0.4.0 current, v0.5.0 planned)
- **Scaling approach:** Jacobian (kinetic-only) → Global conformal transformation
- **Isoelectronic accuracy:** 31-45% → **10-15%** (20-30 point improvement)

### Validated

- ✅ Global metric scaling preserves virial theorem
- ✅ E/Z² ratio nearly constant (GeoVac: -0.713 to -0.723)
- ✅ Li+ (Z=3, 2e): -6.489 Ha (10.87% error)
- ✅ Be2+ (Z=4, 2e): -11.572 Ha (15.22% error)
- ✅ Linear H3 transition state: -1.321 Ha (19.94% error)
- ✅ Conformal transformation theory validated

### Deprecated

- **Jacobian scaling** (scaling only kinetic energy by Z²) - causes virial mismatch
- Use **global metric scaling** instead for isoelectronic series

### Documentation

- [RELEASE_NOTES_v0.4.0.md](RELEASE_NOTES_v0.4.0.md) - Detailed release notes
- [docs/GLOBAL_METRIC_SCALING_SUCCESS.md](docs/GLOBAL_METRIC_SCALING_SUCCESS.md) - Technical analysis

## [0.3.2] - 2026-02-14

### Added

- **AtomicSolver class** - Pure geometric formulation for single-electron atoms
- **Z²-scaling for hydrogenic ions** - Automatic scaling for H, He+, Li2+, etc.
- **`solve_atom()` convenience function** - Quick single-electron calculations
- Comprehensive benchmark suite for validation
- Complete documentation for universal kinetic scale

### Validated

- ✅ H (Z=1): -0.497 Ha (0.57% error at max_n=30)
- ✅ He+ (Z=2): -1.989 Ha (0.57% error at max_n=30)
- ✅ Li2+ (Z=3): -4.474 Ha (0.57% error at max_n=30)
- ✅ Universal kinetic scale -1/16 works for all single-electron systems
- ✅ Z²-scaling formula exact: `kinetic_scale_eff = -1/16 × Z²`

### Changed

- Version: 0.3.1 → 0.3.2
- README updated with AtomicSolver examples
- Documentation expanded for single-electron systems

## [0.3.1] - 2026-02-13

### Added

- **Multi-solver architecture** - Mean-Field, Geometric-DFT, Full CI, Dirac
- **Geometric-DFT** - Fast correlation functional (5.7% error, 79% recovery)
- **Full CI for 2-electron systems** - Exact correlation (<1% with optimization)
- **Dirac relativistic solver** - Spinor formalism with relativistic corrections
- **Geometry optimization** - PES scanning and bond length optimization

### Validated

- ✅ H₂ Mean-Field: -0.980 Ha (16.5% error)
- ✅ H₂ Geometric-DFT: -1.108 Ha (5.7% error)
- ✅ H₂ Full CI (R=1.40): -1.142 Ha (2.8% error)
- ✅ H₂ Full CI (R=1.30 optimized): -1.169 Ha (0.43% error) ⭐

## [0.2.1] - 2026-02-13

### 🔬 Major Scientific Discoveries

#### Universal Constant Discovery
- **DISCOVERED:** `kinetic_scale = -1/16` is a fundamental topological invariant, not a fitting parameter
- **VALIDATED:** Across H (Z=1), He⁺ (Z=2), and H₂⁺ with <0.1% error
- **PHYSICAL MEANING:** Dimensionless ground state eigenvalue of vacuum lattice is exactly 8

#### H₂⁺ Control Experiment
- **PROVEN:** Graph topology correctly models covalent bonding (0% error for H₂⁺)
- **CONFIRMED:** 17% H₂ discrepancy is correlation energy, not topological flaw
- **LITMUS TEST:** Single-electron H₂⁺ validates mean-field framework

#### Mean-Field Classification
- **CLASSIFIED:** GeoVac as Topological Hartree-Fock solver
- **SINGLE-ELECTRON:** Exact accuracy (0% error)
- **MULTI-ELECTRON:** Mean-field quality (~17% correlation error, expected)

#### Bridge Scaling Physics
- **MECHANISM:** Super-linear scaling (α≈1.1) from angular momentum recruitment
- **EVIDENCE:** 90% high-l states (f,g,h,i) participate at n=25
- **PHYSICAL:** Mimics d/f orbital chemistry in heavy elements

### Added

- `UNIVERSAL_KINETIC_SCALE = -1/16` constant in `geovac/__init__.py`
- `HYDROGEN_GROUND_STATE`, `H2_PLUS_USES_UNIVERSAL_SCALE`, `H2_CORRELATION_ERROR` constants
- Physics classification section in package docstring
- `validate_universal_constant.py` - Comprehensive validation tool for H/He⁺/H₂⁺
- `analyze_bridge_distribution.py` - Physical analysis of bridge scaling
- `CORE_PRODUCT_STATUS.md` - Complete status report
- H₂⁺ control experiment documentation in README
- Molecular bonding correlation test section in README
- Universal constant section in README with validation data
- Mean-field classification documentation throughout
- Paper 5 appendix: H₂⁺ experiment and bridge scaling physics

### Changed

- **DEFAULT PARAMETER:** `MoleculeHamiltonian(..., kinetic_scale)`: `-0.075551` → `-1/16`
- **PACKAGE DESCRIPTION:** From empirical to "Topological Hartree-Fock solver"
- **PERFORMANCE CLAIMS:** "~35% error for H₂" → "0% H₂⁺, ~17% H₂ (correlation)"
- **BRIDGE SCALING:** Updated from static N=16 to dynamic N≈4×max_n
- **ERROR ATTRIBUTION:** Clarified correlation vs topology
- `demo_h2.py` to use universal constant with validation references
- `geovac/__init__.py` docstring to reflect mean-field nature
- `geovac/hamiltonian.py` documentation and examples

### Fixed

- Theoretical foundation: Framework now has first-principles basis
- Error attribution: Clear separation of topology (exact) vs correlation (missing)
- Bridge scaling mechanism: Physical origin identified and validated
- Documentation: Proper classification and realistic performance claims

### Validated

- ✅ Universal constant convergence (H, He⁺, H₂⁺)
- ✅ Single-electron topology (0% error for H₂⁺)
- ✅ Multi-electron mean-field behavior (17% correlation in H₂)
- ✅ Angular momentum recruitment in bridge scaling
- ✅ All existing tests pass with new constant

### Backward Compatibility

- ✅ **MAINTAINED:** Existing code with explicit `kinetic_scale` still works
- ✅ **NEW DEFAULT:** Code without explicit parameter uses universal constant
- ✅ **API STABLE:** No breaking changes to method signatures

## [0.2.0] - 2026-02-12

### Added

- `MoleculeHamiltonian` class for molecular bonding
- `GeometricLattice.stitch_lattices()` method for bridge connections
- Spectral delocalization bonding mechanism
- `demo_h2.py` - Complete H₂ molecule demonstration
- Bridge priority ranking system
- Wavefunction delocalization analysis
- Binding energy calculations
- Performance benchmarks for molecules

### Changed

- README: Updated to "First Topological Quantum Chemistry Solver"
- Documentation: Added molecular bonding examples
- Examples: Updated with H₂ demonstrations

### Fixed

- Matrix sparsity maintenance in molecular systems
- Bridge connectivity for optimal bonding

## [0.1.0] - 2026-02-01

### Added

- Initial release
- `GeometricLattice` class for atomic systems
- `HeliumHamiltonian` class for two-electron atoms
- `DiracHamiltonian` class (experimental)
- Graph Laplacian based kinetic energy
- Sparse matrix eigenvalue solver
- Basic documentation and examples

---

## Version Naming Convention

- **Major (X.0.0):** Breaking API changes
- **Minor (0.X.0):** New features, backward compatible
- **Patch (0.0.X):** Bug fixes, documentation updates

## Links

- [v0.2.1 Release Notes](RELEASE_NOTES_v0.2.1.md) - Detailed release documentation
- [v0.2.0 Release Notes](RELEASE_NOTES_v0.2.0.md) - Previous release
- [Core Product Status](CORE_PRODUCT_STATUS.md) - Complete validation report

[0.2.1]: https://github.com/your-org/geovac/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/your-org/geovac/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/your-org/geovac/releases/tag/v0.1.0
