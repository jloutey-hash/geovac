# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v1.7.7 (March 22, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (unit S3) that is mathematically equivalent to the Schrodinger equation via Fock's 1935 conformal projection. This equivalence is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Authoritative source rule:** The core papers in `papers/core/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it.

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers in `papers/core/` are written to academic standards but are not submitted to traditional journals. The project's viability case rests on producing a usable, benchmarked computational tool. Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks.

---

## 2. Current Development Frontier

**Best results by system type:**
- Atoms: He at 0.05% error (hyperspherical coordinates, Paper 13)
- One-electron molecules: H2+ at 0.70% (prolate spheroidal, Paper 11)
- Two-electron molecules: H2 at 94.1% D_e (Level 4 mol-frame hyperspherical, Paper 15)
- Core-valence molecules: LiH R_eq = 3.21 bohr, 6.4% error (composed geometry, Paper 17)
- Fine structure constant: alpha from Hopf bundle at 8.8x10^-8, zero free parameters, p-value 5.2x10^-9, universal algebraic identity B_formal/N = d, Hopf generalization negative result, circulant Hermiticity, second selection principle (Paper 2, conjectural)

**Active frontier:**
- Extending composed geometry (Paper 17) to higher l_max with sigma+pi channels for LiH
- Polyatomic composition pattern (BeH2, H2O) via fiber bundle generalization
- Improving Level 4 D_e recovery beyond 94% (channel convergence, cusp corrections)
- Chemical periodicity as representation theory (Paper 16) -- computational exploitation of hierarchical structure

**Architecture locked:** The LCAO/graph-concatenation approach (v0.9.x series) is superseded. All molecular work uses natural geometry (Papers 11, 13, 15, 17).

---

## 3. Approaches That Failed

Critical institutional memory. Do not re-derive these dead ends.

| Approach | Why It Fails | Resolution | Reference |
|:---------|:-------------|:-----------|:----------|
| LCAO graph concatenation for molecules | Graph Laplacian kinetic energy is R-independent -> monotonically attractive PES, no equilibrium | Prolate spheroidal lattice (Paper 11) or composed geometry (Paper 17) | FCI-M, 29-version diagnostic arc |
| Sturmian basis with shared p0 | H proportional to S -> eigenvalues R-independent | Prolate spheroidal separation introduces R-dependence via beta_k(R) | Papers 8-9, Structural Theorem |
| Berry phase from lattice plaquettes | arg() = 0 identically for real SU(2)/SU(1,1) operators | Log-holonomy Theta(n) = -2ln((n+1)/n) ~ 1/n is the valid geometric quantity | Paper 1 v1.2.0 erratum |
| Numerical V_ee quadrature on prolate spheroid | Coulomb singularity at r12=0 causes slow convergence, saturates at ~80% D_e | Algebraic Neumann expansion eliminates quadrature entirely | Paper 12 |
| Single-S3 molecular encoding | One p0 cannot encode R-dependent bonding physics | Natural geometry principle: use coordinate system where separation occurs | Papers 8-9, 11 |
| Orbital exponent relaxation (zeta) | Shifts PES uniformly, not differentially; R_eq unchanged | Not a mechanism for equilibrium geometry | v0.9.36 |
| Mulliken cross-nuclear diagonal | Too strong at short R, drives R_eq inward | Bond sphere (Paper 8) or natural geometry approach | v0.9.35 diagnosis |

---

## 4. The Dimensionless Vacuum Principle

The graph Laplacian is dimensionless and scale-invariant, topologically equivalent to the unit three-sphere S3. Physical energies emerge only through the energy-shell constraint p0^2 = -2E, which acts as a stereographic focal length. The 1/r Coulomb potential is not an input force law -- it is the coordinate distortion from stereographic projection (chordal distance identity). The universal kinetic scale kappa = -1/16 maps graph eigenvalues to the Rydberg spectrum. Eigenvalues of the Laplace-Beltrami operator on unit S3 are pure integers: lambda_n = -(n^2 - 1).

For molecules, the natural geometry shifts from S3 (atoms) to prolate spheroidal coordinates (Paper 11), hyperspherical coordinates (Paper 13), molecule-frame hyperspherical (Paper 15), or composed fiber bundles (Paper 17). The choice of geometry is determined by where separation of variables occurs.

**Prime directive:** Never modify the discrete graph Laplacian to artificially recover continuous differential terms (like 1/r or nabla^2). The graph is an exact, dimensionless S3 topology. The Schrodinger equation is its flat-space projection. If eigenvalues do not match expected physics, the issue is in the projection or energy-shell constraint, not the graph.

**Topological integrity tests:** The 18 symbolic proofs in `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py` validate S3 conformal geometry. These must never be broken or bypassed. Run before any release.

---

## 5. Natural Geometry Hierarchy

The core organizational principle of the project. Each electron configuration has a natural coordinate system where the physics separates.

| Level | System | Natural Geometry | Best Result | Paper |
|:-----:|:-------|:-----------------|:------------|:-----:|
| 1 | H (1-center, 1e) | S3 (Fock) | < 0.1% | 7 |
| 2 | H2+ (2-center, 1e) | Prolate spheroid | 0.70% | 11 |
| 3 | He (1-center, 2e) | Hyperspherical | 0.05% | 13 |
| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 94.1% D_e | 15 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq 6.4% | 17 |

The composed geometry (Level 5) is a fiber bundle: G_total = G_nuc semi-direct G_core(R) semi-direct G_val(R, core_state). Each electron group gets its own natural coordinate system, coupled via Z_eff screening and Phillips-Kleinman pseudopotential.

---

## 6. Paper Series

### Reading Guide

1. **Start here:** Paper 7 (the theoretical foundation -- graph Laplacian = S3 = Schrodinger)
2. **Atoms:** Papers 0, 1 (graph construction, eigenvalue methods), then FCI-A (multi-electron)
3. **Multi-electron atoms:** Paper 13 (hyperspherical lattice, He at 0.05%, fiber bundle)
4. **Dynamics:** Paper 6 (time evolution, spectroscopy, AIMD on graph Hamiltonians)
5. **Molecules -- the problem:** Papers 8-9 (bond sphere geometry, why single-S3 fails)
6. **Molecules -- the solution:** Paper 11 (prolate spheroidal lattice, H2+ zero free params)
7. **Molecules -- two-electron:** Paper 12 (algebraic V_ee, Neumann expansion, H2 92.4% D_e)
8. **Molecules -- natural geometry:** Paper 15 (Level 4 hyperspherical, H2 94.1% D_e)
9. **Molecules -- core-valence:** Paper 17 (composed geometry, LiH 6.4%, BeH+ bound)
10. **Periodicity:** Paper 16 (periodic table from S_N representation theory on S^(3N-1))
11. **Ab initio spectroscopy:** Paper 13 Sec IX (PES -> Morse -> nuclear lattice)
12. **Quantum computing:** Paper 14 (qubit encoding, O(Q^3.15) Pauli scaling)

### Paper Inventory

#### Core (`papers/core/`)

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Paper 0 | `Paper_0_Geometric_Packing.tex` | Active | Universal constant K = -1/16 |
| Paper 1 | `paper_1_spectrum.tex` | Active | Spectral graph theory, O(N) eigenvalue methods. Berry phase retracted v1.2.0 |
| Paper 6 | `Paper_6_Quantum_Dynamics.tex` | Active | O(V) dynamics: Rabi, spectroscopy, AIMD |
| Paper 7 | `Paper_7_Dimensionless_Vacuum.tex` | Active | S3 proof (18/18 symbolic), Schrodinger recovery, SO(3N) generalization |
| Paper 10 | `paper_10_nuclear_lattice.tex` | Draft | Rovibrational spectra from SU(2) algebraic chains |
| Paper 11 | `paper_11_prolate_spheroidal.tex` | Draft | Prolate spheroidal lattice: H2+ 0.70% E_min |
| Paper 12 | `paper_12_algebraic_vee.tex` | Active | Neumann V_ee: H2 92.4% D_e, cusp diagnosis (7.6% gap) |
| Paper 13 | `paper_13_hyperspherical.tex` | Active | Hyperspherical lattice: He 0.05%, fiber bundle, ab initio spectroscopy |
| Paper 14 | `paper_14_qubit_encoding.tex` | Active | Structurally sparse qubits: O(Q^3.15) vs O(Q^4.60) Gaussian |
| Paper 15 | `paper_15_level4_geometry.tex` | Active | Level 4: H2 94.1% D_e, HeH+ 93.1% D_e |
| Paper 16 | `paper_16_periodicity.tex` | Active | Periodic table from S_N representation theory on S^(3N-1) |
| Paper 17 | `paper_17_composed_geometries.tex` | Active | Composed geometry: LiH R_eq 6.4%, BeH+ bound, ab initio PK |
| FCI-A | `paper_fci_atoms.tex` | Draft | He 0.35%, Li 1.10%, Be 0.90% |
| FCI-M | `paper_fci_molecules.tex` | Scaffold | LCAO LiH results and diagnostic arc |

#### Methods (`papers/methods/`)

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Papers 8-9 | `Paper_8_Bond_Sphere_Sturmian.tex` | Draft | Bond sphere (positive), Sturmian structural theorem (negative), SO(4) selection rules |

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

---

## 10. Validation Benchmarks

| Test | Max Error | Purpose |
|:-----|:---------:|:--------|
| Symbolic proofs (18 tests) | 0 failures | Topological foundation |
| H (hydrogen) | < 0.1% | Basic validation |
| He+ (helium ion) | < 0.1% | Z-scaling check |
| H2+ (ionized H2) | < 0.1% | Topological control |
| He (hyperspherical) | < 0.1% | Multi-electron control |
| H2 Full CI | < 1.0% | Accuracy control |
| H2 Neumann V_ee | 92.4% D_e | Algebraic integral accuracy |
| H2 Level 4 (2D solver) | 94.1% D_e | Molecule-frame hyperspherical |
| HeH+ Level 4 | 93.1% D_e | Heteronuclear extension |
| LiH Composed (ab initio PK) | R_eq 6.4% | Composed geometry |
| BeH+ Composed | Bound, physical | Transferability |
| Hyperspherical (20 tests) | 0 failures | Angular + adiabatic + radial |
| Muonic H energy ratio | < 0.01% | Mass-independence |
| V_ee S3 overlap (1s-1s, 1s-2s, 2s-2s) | < 0.01% | Topological integrity |
| Direct CI vs matrix CI | < 1e-8 Ha | Algorithmic consistency |
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
| Qubit Pauli scaling | 14 | All | Core |
| Structural sparsity | 14 | Sec III | Core |
| Level 4 mol-frame hyperspherical | 15 | All | Core |
| Mol-frame charge function | 15 | Sec III | Core |
| Multichannel expansion | 15 | Sec V | Core |
| Heteronuclear extension | 15 | Sec V.D | Core |
| Variational 2D solver | 15 | Sec VI.D | Core |
| HeH+ convergence | 15 | Sec VI.E | Core |
| Double-adiabatic fiber bundle | 15 | Sec VII.C | Core |
| Chemical periodicity (S_N reps) | 16 | All | Core |
| Structure types A/B/C/D/E | 16 | Sec III | Core |
| Hierarchical decomposition | 16 | Sec IV | Core |
| Dirac limit (Z~137) | 16 | Sec VI | Core |
| Composed natural geometries | 17 | All | Core |
| Core-valence fiber bundle | 17 | Sec II | Core |
| Ab initio Phillips-Kleinman | 17 | Sec IV | Core |
| Rho-collapse cache | 17 | Sec V | Core |
| LiH/BeH+ benchmarks | 17 | Sec VI | Core |
| Bond sphere theory | 8-9 | All | Methods |
| Sturmian structural theorem | 8-9 | Sec IV | Methods |
| SO(4) selection rules | 8-9 | Sec III | Methods |
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
