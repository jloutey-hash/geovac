# GeoVac: Computational Quantum Chemistry via Spectral Graph Theory

![Status](https://img.shields.io/badge/Status-Production-brightgreen) ![Version](https://img.shields.io/badge/Version-2.4.0-blue) ![License](https://img.shields.io/badge/License-MIT-orange)

**Version 2.4.0** - Transition metal scoping, d-orbital sparsity, 30-molecule library

GeoVac constructs **structurally sparse qubit Hamiltonians** for quantum simulation of molecular systems, achieving O(Q^2.5) Pauli term scaling -- a 51x to 1,712x advantage over published Gaussian baselines (Paper 14). The sparsity is intrinsic to the basis: angular momentum selection rules in the Gaunt integrals enforce block-diagonal electron repulsion integrals, producing Hamiltonians that are sparse in structure and concentrated in weight (O(Q^1.69) 1-norm).

The framework rests on a proven mathematical equivalence: the discrete graph Laplacian on $S^3$ is conformally equivalent to the Schrodinger equation for hydrogen (Fock 1935, Paper 7, 18 symbolic proofs). For molecules, the **natural geometry principle** extends this to prolate spheroidal, hyperspherical, and composed fiber-bundle coordinates (Papers 11-17). Classical solver results across 30+ investigation tracks serve as validation benchmarks for the quantum Hamiltonians.

---

## Development Methodology

GeoVac is developed using an AI-augmented research workflow. The principal investigator (J. Loutey) provides scientific direction, physical intuition, and quality control. Implementation, numerical exploration, and documentation drafting are performed collaboratively with large language models (Anthropic Claude). All physics results are validated against known analytical solutions, NIST reference data, and a symbolic proof suite (18 independent tests verifying the S³ conformal geometry). No result is accepted on the basis of LLM output alone.

This workflow is itself a research contribution — an experiment in whether agentic AI tools can accelerate independent scientific research. The project has no institutional affiliation. Primary dissemination is via GitHub and Zenodo (DOI-stamped releases).

---

## What's New in v2.0.42

### Balanced Coupled Composition (v2.0.39-42, Track CD — MIXED)

- **v2.0.39-42 — Track CD: Balanced coupled composition.** Cross-center nuclear attraction integrals via multipole expansion replace PK pseudopotential. Non-collinear geometries via Wigner D-matrix rotation (Sprint 5). Complete polyatomic census:

| System | Blocks | Q | Composed Pauli | Balanced Pauli | Ratio | λ_comp (Ha) | λ_bal (Ha) | λ ratio |
|--------|:------:|:-:|:--------------:|:--------------:|:-----:|:-----------:|:----------:|:-------:|
| LiH | 2 | 30 | 334 | 878 | 2.63× | 37.3 | 74.1 | 1.98× |
| BeH₂ | 3 | 50 | 556 | 2,652 | 4.77× | 354.9 | 304.7 | 0.86× |
| H₂O | 5 | 70 | 778 | 5,798 | 7.45× | 361/28,053 | 1,509.3 | 4.18×/0.054× |

- LiH is the only bound 4e FCI config (1.8% energy, 7.0% R_eq). n_max=3 analytical: 0.20% energy.
- BeH₂: balanced 1-norm < composed (0.86×) — balanced cheaper for QPE.
- H₂O: balanced eliminates PK (18.6× cheaper than composed w/ PK). FCI infeasible (10e).
- Pauli ratio grows with block count: 2.63× → 4.77× → 7.45×.

### Transcorrelated Integrals (v2.0.35-36, Tracks BX-3/BX-4)

TC-modified qubit Hamiltonians via the composed pipeline bypass the adiabatic solver entirely:

| System | Standard err% | TC err% | Pauli ratio |
|--------|:------------:|:-------:|:-----------:|
| He (max_n=1) | 5.3% | 3.3% | 1.68× |
| He (max_n=3) | 8.2% (diverging) | 3.6% (converging) | 1.68× |

O(Q^2.5) scaling preserved. TC angular gradient (Track BX-4) is a **negative result**: 2.67× Pauli increase for <0.02 pp accuracy. Radial-only is the default.

### Quantum Resource Market Test (v2.0.36, Track CA)

Head-to-head 1-norm comparison against computed Gaussian baselines:

| Metric | GeoVac LiH (Q=30) | Gaussian STO-3G (Q=12) | Advantage |
|--------|:-----------------:|:---------------------:|:---------:|
| Pauli terms | 334 | 907 | 2.7× |
| QWC groups | 21 | 273 | 13× |
| 1-norm (λ) | 33.3 Ha | 34.3 Ha | 0.97× (match) |

He equal-qubit at Q=28: GeoVac λ=78.4 vs Gaussian cc-pVTZ λ=530.5 (6.8× lower).

### General Composed Builder (v2.0.30, Tracks BG-BK)

Single `build_composed_hamiltonian(spec)` for all molecules, driven by `MolecularSpec` dataclass and atomic classifier:

```python
from geovac.ecosystem_export import hamiltonian
h = hamiltonian('LiH', R=3.015, l_max=2)  # 334 Pauli terms, 30 qubits
op_qiskit = h.to_qiskit()    # SparsePauliOp
op_of = h.to_openfermion()   # QubitOperator
op_pl = h.to_pennylane()     # qml.Hamiltonian
```

Six composed molecules: LiH (334), BeH₂ (556), H₂O (778), HF (667), NH₃ (889), CH₄ (1000 Pauli terms). Consistent ~Q^2.2 within-molecule scaling.

### Coupled Composition (v2.0.37, Track CB — NEGATIVE)

Replacing PK with cross-block ERIs *without* cross-center V_ne: 29% error (worst). Decoupled blocks (10.9%) outperform PK (15.0%), confirming PK overcorrection. Root cause: missing two-center nuclear attraction integrals.

### Paper 19 — Balanced Coupled Composition (v2.0.39-43, Track CD — Supporting)

Balanced coupled architecture: cross-center V_ne via multipole expansion (exact termination at L_max=2*l_max by Gaunt selection rules), analytical integrals via incomplete gamma functions. LiH 0.20% energy error at n_max=3, 878 Pauli terms. 3-molecule census: LiH/BeH₂/H₂O. Only bound 4-electron FCI configuration. PK-free, regime-dependent 1-norm advantage over composed framework.

### Ecosystem & Infrastructure (v2.0.27-29)

- **Ecosystem export:** OpenFermion, Qiskit, PennyLane via `geovac.ecosystem_export`
- **PK partitioning:** Classical 1-RDM evaluation removes PK from quantum circuit (78× 1-norm reduction for H₂O)
- **geovac-hamiltonians:** Standalone PyPI package (92 KB wheel, 5 systems)
- **IBM Quantum demo:** VQE on Aer simulator and IBM hardware modes
- **VQE validation:** H₂ converges to 0.031 mHa on statevector simulator

### Quantum Encoding Summary

| Metric | GeoVac (Composed) | Gaussian Baselines | Advantage |
|--------|:-----------------:|:------------------:|:---------:|
| Pauli scaling (molecules) | O(Q^2.5) | O(Q^3.9-4.3) | Structural |
| LiH (Q=30) | 334 terms | 17,065 (Trenev) | 51x |
| BeH₂ (Q=50) | 556 terms | 256,000 (Trenev) | 460x |
| H₂O (Q=70) | 778 terms | 1.33M (Trenev) | 1,712x |
| He 1-norm (Q=28) | 78.4 Ha | 530.5 Ha (cc-pVTZ) | 6.8x |
| Full vs composed (Q=10) | 334 terms | 3,288 terms (full N-e) | 10x sparser |

### Classical Validation Benchmarks (Complete)

- H₂ at 96.0% D_e (Level 4, l_max=6 with cusp correction, CBS ~97%)
- LiH R_eq at 5.3% (Level 5 composed, l-dependent PK, l_max=2)
- BeH₂ R_eq at 11.7% (full 1-RDM exchange, zero free parameters)
- H₂O R_eq at 26% (5-block composed, charge-center origin)
- Full 4-electron LiH equilibrium without PK (Level 4N structural validation)

### Prior Releases
- **v2.0.26:** Ecosystem export pipeline, head-to-head Gaussian comparison (190× fewer Pauli for LiH)
- **v2.0.21:** Track AM (l_max=3 convergence), Track AN (full documentation update)
- **v2.0.14-20:** Cusp attack (3 tracks), N-electron architecture (Tracks AF-AS), spectral solvers
- **v2.0.6-13:** Algebraic infrastructure, spectral bases, exchange constants (Paper 18)
- **v2.0.0-5:** Composed multi-center qubit Hamiltonians, diagnostics, polyatomics
- **v1.7-1.9:** VQE benchmarks, measurement cost analysis, Paper 2 rewrite, Paper 17
- **v1.0-1.6:** Foundation papers (0, 1, 6-16), solvers, Direct CI, Bond Sphere theory

---

## Precision Benchmarks

### Atoms (Graph Laplacian on $S^3$)

| System | Method | Error | Status |
|--------|--------|:-----:|:------:|
| H, He+, H$_2^+$ (1e) | Graph Laplacian | < 0.1% | Exact mean-field |
| He (2e, FCI) | slater_full + hybrid h1 | **0.35%** | Monotonic convergence |
| **He (2e, hyperspherical)** | **Adiabatic 1 ch.** | **0.05%** | **Best GeoVac He** |
| Li (3e, FCI) | slater_full + exact h1 | **1.10%** | Monotonic convergence |
| Be (4e, FCI) | Direct CI, nmax=4 | **0.90%** | Sub-1% |
| Li+ (2e, Z=3) | Conformal torsion | **0.039%** | Isoelectronic scaling |

### Molecules

| System | Method | Error | Status |
|--------|--------|:-----:|:------:|
| H$_2^+$ (prolate spheroidal) | Separated solver | **0.70%** | Zero free params |
| H$_2$ (prolate CI, relaxed) | Eckart SCF + 2×2 CI | 58% $D_e$ | Zero free params |
| **H$_2$ (Neumann $V_{ee}$)** | **Algebraic CI, 27 bf** | **92.4% $D_e$** | **Paper 12** |
| **H$_2$ (Level 4, 2D+cusp)** | **Mol.-frame hypersp.** | **96.0% $D_e$** | **Paper 15** |
| **HeH$^+$ (Level 4)** | **Charge-center hypersp.** | **93.1% $D_e$** | **Paper 15** |
| H$_2$ (Full CI) | Tensor product | 1.72% | Cross-nuclear ceiling |
| **LiH (composed, ab initio PK)** | **Composed Level 3+4** | **R_eq: 5.3%** | **Paper 17** |
| **BeH⁺ (composed, ab initio PK)** | **Composed Level 3+4** | **Bound** | **Paper 17** |
| **BeH₂ (composed, exchange)** | **Composed Level 3+4, exchange coupling** | **R_eq: 11.7%** | **Paper 17** |
| **LiH (Level 4N, full 4e)** | **Full mol-frame hypersp. (SO(12))** | **R_eq: 63.5% (structural; 2D unbound D_e)** | **Paper 17** |
| **H₂O (composed, uncoupled)** | **Composed Level 3+4, charge-center** | **R_eq: 26%** | **Paper 17** |
| LiH ($D_e$, CP-corrected) | LCAO FCI | **1.0%** | nmax=3 |

### Ab Initio Spectroscopy

| Benchmark | Ab Initio | Experiment | Error |
|-----------|-----------|------------|:-----:|
| H$_2$ $R_e$ (bohr) | 1.418 | 1.401 | +1.2% |
| H$_2$ $\omega_e$ (cm$^{-1}$) | 4435 | 4401 | +0.8% |
| H$_2$ $B_e$ (cm$^{-1}$) | 59.49 | 60.85 | -2.2% |
| H$_2$ $\nu_{01}$ (cm$^{-1}$) | 4157 | 4161 | -0.1% |
| LiH $R_{\rm eq}$ (bohr) | 3.18 | 3.015 | +5.3% |
| LiH $\omega_e$ (cm$^{-1}$) | 1471 | 1406 | +4.6% |
| LiH $B_e$ (cm$^{-1}$) | 6.63 | 7.51 | -11.7% |

### Dynamics & Spectroscopy

| Benchmark | Result | Status |
|-----------|--------|:------:|
| Rabi oscillation | 0.41% period error | PASS |
| H$_2$ delta-kick spectroscopy | 20/35 peaks, 0.16% | PASS |
| AIMD energy conservation | 0.0003% max drift | PASS |
| Langevin bond dissociation | Correct thermal behavior | PASS |

### Qubit Encoding (Paper 14)

| Metric | GeoVac Scaling | R² | Gaussian (conventional) |
|--------|:--------------:|:--:|:-----------------------:|
| Pauli terms (atoms) | O(Q^3.15) | 0.9995 | O(Q^4.60) molecules / O(Q^3.56) atoms |
| Pauli terms (composed) | O(Q^2.5) | 0.991+ | O(Q^3.92–4.25) (Trenev et al.) |
| QWC measurement groups | O(Q^3.36) | 1.0000 | — |
| Pauli 1-norm (λ) | O(Q^1.69) | 0.9972 | — |
| Trotter steps (ε=10⁻³) | O(Q^1.69) | 0.9972 | — |

#### Composed Multi-Center Systems (v2.0.0)

| Molecule | Blocks | Q (nmax=2) | Pauli Terms | vs Gaussian | Scaling α |
|----------|:------:|:----------:|:-----------:|:-----------:|:---------:|
| LiH | 3 | 30 | 334 | 51× fewer | 2.50 |
| BeH₂ | 5 | 50 | 556 | 460× fewer | 2.51 |
| H₂O | 7 | 70 | 778 | 1,712× fewer | 2.52 |

#### Equal-Qubit Comparison (He atom, validated v1.9.0)

| Q | GeoVac Encoding | Error | Pauli Terms | Gaussian Encoding | Error | Pauli Terms | Advantage |
|:-:|:----------------|:-----:|:-----------:|:------------------|:-----:|:-----------:|:---------:|
| 10 | nmax=2 (5 spatial) | 0.56% | 120 | cc-pVDZ (5 spatial) | 0.56% | 156 | 1.3× |
| 28 | nmax=3 (14 spatial) | 0.45% | 2,659 | cc-pVTZ (14 spatial) | 0.12% | 21,607 | 8.1× |

Gaussian Pauli counts are actual JW term counts from computed MO integrals (not scaling estimates).

### GeoVac vs PySCF (CI-validated)

| System | GeoVac (best) | PySCF/STO-3G | PySCF/cc-pVQZ |
|--------|:------:|:------------:|:-------------:|
| H | 0.57% | 6.68% | 0.01% |
| He | **0.05%** (hyperspherical) | 3.30% | 1.45% |
| H$_2$ | 1.73% | 3.17% | 0.06% |

GeoVac beats PySCF/STO-3G on all tested atoms. He (hyperspherical) beats PySCF/cc-pVQZ.

---

## The Universal Kinetic Constant (-1/16)

$$\text{Scale} = -\frac{1}{16} = -0.0625$$

The unique kinetic scale mapping dimensionless graph eigenvalues to the physical hydrogen spectrum, determined by requiring agreement with the Rydberg formula (Paper 7, Sec. V). Universal across all tested systems (H, He+, Li2+, He, H2, H2+). The Z-scaling formula `kinetic_scale_eff = -1/16 * Z^2` is exact.

---

## Quick Start

### Installation
```bash
git clone https://github.com/jloutey-hash/geovac.git
cd geovac
pip install -e .
```

### Single-Electron Atoms
```python
from geovac import AtomicSolver, solve_atom

# Hydrogen
E, psi = solve_atom(Z=1, max_n=30)
print(f"H ground state: {E:.6f} Ha")  # -0.497 Ha (0.57% error)

# Any hydrogenic ion (automatic Z^2 scaling)
E, psi = solve_atom(Z=3, max_n=30)  # Li2+
print(f"Li2+ ground state: {E:.6f} Ha")  # -4.474 Ha
```

### Multi-Electron Atoms (FCI)
```python
from geovac import LatticeIndex

# Helium - Full CI with Slater integrals
idx = LatticeIndex(n_electrons=2, max_n=5, nuclear_charge=2,
                   vee_method='slater_full', h1_method='hybrid')
H = idx.assemble_hamiltonian()
E, psi = idx.compute_ground_state()
print(f"He FCI: {E[0]:.6f} Ha")  # -2.894 Ha (0.35% error)
```

### Helium (Hyperspherical Adiabatic)
```python
from geovac.hyperspherical_radial import solve_helium

result = solve_helium(Z=2.0, l_max=0, n_alpha=200)
print(f"He ground state: {result['energy']:.6f} Ha")  # -2.9052 Ha (0.05% error)
```

### H2+ Molecule (Prolate Spheroidal)
```python
from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice

lattice = ProlateSpheroidalLattice(R=2.0, Z_A=1, Z_B=1, N_xi=5000)
E_total = lattice.total_energy()
print(f"H2+ at R=2.0: {E_total:.6f} Ha")  # -0.597 Ha (0.70% error)
```

### LiH Molecule (Composed Natural Geometry)
```python
from geovac.composed_diatomic import ComposedDiatomicSolver

solver = ComposedDiatomicSolver.LiH_ab_initio(l_max=2)
results = solver.run_all()
print(f"LiH R_eq: {results['spectro']['R_eq']:.3f} bohr")  # 3.21 bohr (expt: 3.015)
print(f"LiH ω_e: {results['spectro']['omega_e']:.0f} cm⁻¹")  # 1471 cm⁻¹ (expt: 1406)
```

### Time Evolution (Rabi Oscillations)
```python
from geovac import AtomicSolver, TimePropagator

solver = AtomicSolver(max_n=5, Z=1)
E, psi = solver.compute_ground_state(n_states=5)

# Build dipole operator and propagate
V_z = TimePropagator.build_dipole_z(solver.lattice)
prop = TimePropagator(solver.H, dt=0.01)
psi_t = prop.evolve(psi[:, 0].astype(complex), n_steps=1000)
```

---

## Solver Methods

| Method | Complexity | Accuracy | Use Case |
|--------|-----------|----------|----------|
| `mean_field` | O(N) | Exact for 1e, ~17% for 2e | Fast screening |
| `full_ci` | O(N_SD) | < 1% for 2-4e | Quantitative benchmarks |
| `direct` | O(N_SD * N_conn) | < 1% (same as full_ci) | Large determinant spaces |
| `dirac` | O((2N)^2) | Relativistic corrections | Heavy atoms |
| `frozen_core` | O(N_active^4) | < 2% | Core-valence separation |
| `locked_shell` | O(N_active^4) | < 2% | Tensor product of shell states |
| Prolate spheroidal | O(N_xi) | < 1% for 1e diatomics | Molecular PES |
| Composed (Level 3+4) | O(N_rho * N_ch^2 + N_R * N_Re) | ~6% R_eq | Core-valence diatomics |

---

## Project Structure

```
geovac/                 Core package
  lattice.py              GeometricLattice (quantum state graph)
  lattice_index.py        LatticeIndex / MolecularLatticeIndex (N-electron FCI)
  direct_ci.py            DirectCISolver (excitation-driven FCI)
  frozen_core.py          FrozenCoreLatticeIndex (frozen-core CI)
  locked_shell.py         LockedShellMolecule (locked-shell CI)
  hamiltonian.py          MoleculeHamiltonian (molecular bonding)
  prolate_spheroidal_lattice.py   Prolate spheroidal solver (H2+)
  prolate_scf.py            Prolate spheroidal SCF + CI (H2)
  neumann_vee.py            Algebraic V_ee via Neumann expansion
  hyperspherical_angular.py Hyperangular solver (Gaunt + Liouville FD)
  hyperspherical_adiabatic.py Adiabatic potential curves
  hyperspherical_radial.py  Hyperradial solver + solve_helium()
  hyperspherical_complex_scaling.py  ECS for resonance detection
  hyperspherical_coupling.py  Coupled-channel adiabatic solver
  hyperspherical_resonances.py  Resonance analysis tools
  level4_multichannel.py  Level 4 mol.-frame hyperspherical solver
  level4_spectral_angular.py  Spectral Jacobi angular solver (Level 4)
  core_screening.py       Core electron screening Z_eff(r)
  ab_initio_pk.py         Ab initio Phillips-Kleinman pseudopotential
  rho_collapse_cache.py   ρ-collapsed angular cache + fast PES
  pauli_projector.py      Core-valence Pauli projection
  composed_diatomic.py    General composed diatomic solver
  composed_triatomic.py   BeH₂ composed solver (exchange coupling)
  composed_water.py       H₂O composed solver (5-block, lone pairs)
  lone_pair.py            Lone pair solver (Level 3 non-bonding pairs)
  algebraic_angular.py    Algebraic angular solver (Gegenbauer spectral basis)
  algebraic_angular_sturmian.py  Sturmian angular variant (research artifact)
  algebraic_coupled_channel.py   Coupled-channel integration (algebraic P/Q matrices)
  inter_fiber_coupling.py Inter-fiber exchange coupling (monopole, exchange, 1-RDM)
  lih_composed.py         LiH composed solver (original)
  berry_phase.py            Log-holonomy on geometric lattice
  nuclear_lattice.py      Nuclear vibration/rotation lattice
  coupled_en_lattice.py   Coupled electronic-nuclear lattice
  dynamics.py             TimePropagator (Crank-Nicolson)
  aimd.py                 VelocityVerlet, LangevinThermostat
  vqe_benchmark.py        VQE pipeline (OpenFermion→Qiskit, head-to-head comparison)
  gaussian_reference.py   Gaussian baselines (H₂ STO-3G, He STO-3G/cc-pVDZ/cc-pVTZ)
  measurement_grouping.py QWC measurement group analysis
  trotter_bounds.py       Pauli 1-norm and Trotter error bounds
  composed_qubit.py       General composed qubit Hamiltonians (6 molecules via MolecularSpec)
  atomic_classifier.py    Atomic classification Z→block decomposition (Z=1-10)
  molecular_spec.py       MolecularSpec dataclass for general composed builder
  ecosystem_export.py     OpenFermion, Qiskit, PennyLane export pipeline
  pk_partitioning.py      Classical PK energy from VQE 1-RDM
  tc_integrals.py         Transcorrelated integrals for composed pipeline
  coupled_composition.py  Coupled composition (Track CB, negative result)
  cross_block_mp2.py      Cross-block ERI computation
  n_electron_2d.py        2D variational solver for full N-electron (Track AR)
  n_electron_qubit.py     Full N-electron quantum encoding comparison (Track AS)
  cusp_factor.py          Cusp factor solver (Track U, negative result)
  cusp_graph.py           Cusp graph analysis (Track W, negative result)
  cusp_correction.py      Schwartz cusp correction (Track X)
  cusp_angular_basis.py   θ₁₂-adapted angular basis (Track Y, negative result)
  geometric_elevation.py  Geometric elevation analysis (Track Z, negative result)
  algebraic_zeff.py       Algebraic Z_eff (Laguerre spectral, Track N)
  algebraic_slater.py     Algebraic Slater integrals (Track O)
  algebraic_curve.py      Level 3 algebraic curve (Track P1)
  algebraic_curve_level4.py  Level 4 algebraic curve (Track S)
  n_electron_solver.py    Full N-electron mol-frame hyperspherical solver (Level 4N)
  n_electron_scope.py     N-electron channel counting and feasibility analysis
  n_electron_spectral.py  N-electron spectral compression analysis
  hopf_bundle.py          Hopf fibration analysis (S³→S², α derivation)
  wigner_so4.py           SO(4) Wigner D-matrix elements
  atomic_solver.py        Single-electron atom solver
  dirac_hamiltonian.py    Relativistic Dirac solver

papers/
  core/                 Defensible foundations (Papers 0, 1, 6-16, 17-18, FCI)
    Paper 0:  Geometric packing & universal constant
    Paper 1:  Spectral graph theory & eigenvalue methods
    Paper 6:  Quantum dynamics & thermodynamics
    Paper 7:  The Dimensionless Vacuum (S3 proof, 18/18 proofs)
    Paper 11: Molecular Fock Projection (prolate spheroidal)
    Paper 12: Algebraic V_ee (Neumann expansion)
    Paper 13: Hyperspherical lattice (multi-electron atoms)
    Paper 14: Structurally sparse qubit Hamiltonians
    Paper 15: Level 4 mol.-frame hyperspherical (H2 96.0%, HeH+ 93.1%)
    Paper 16: Chemical periodicity from S_N representation theory (atomic classifier)
    Paper 17: Composed natural geometries (LiH 5.3%, BeH⁺)
    Paper 18: Spectral-geometric exchange constants (Weyl-Selberg, α connection)
    FCI-A:    Multi-electron atoms (He, Li, Be)
    Paper 19: Balanced coupled composition (0.20% energy, PK-free, 3-molecule census)
  conjectures/          Theoretical explorations (Papers 2-5)
    Paper 2:  Fine structure α from Hopf bundle (8.8×10⁻⁸ error)
  archive/              Historical papers, superseded versions (Papers 8-9, 10, FCI-M)
    Paper 8-9: Bond Sphere & Sturmian negative theorem, SO(4) selection rules
    Paper 10: Nuclear lattice for vibration/rotation
    FCI-M:    Heteronuclear diatomics (LiH benchmark)

tests/                  Unit tests (pytest) + validation
  test_fock_projection.py   10 proofs (stereographic geometry)
  test_fock_laplacian.py    8 proofs (conformal Laplacian)
  test_prolate_h2plus.py    11 tests (prolate spheroidal solver)
  test_hyperspherical_he.py 20 tests (hyperspherical He solver)
  test_neumann_vee.py       Neumann V_ee algebraic integrals
  test_direct_ci.py         6 tests (Direct CI consistency)
  test_frozen_core.py       14 tests (frozen-core CI)
  test_locked_shell.py      11 tests (locked-shell CI)
  test_lih_fci.py           50 tests (LiH heteronuclear FCI)
  test_hopf_bundle.py       41 tests (Hopf fibration, α formula)
  test_berry_phase.py       7 tests (log-holonomy, power-law)
  test_complex_scaling.py   ECS resonance tests
  test_hyperspherical_coupled.py  Coupled-channel tests
  test_core_screening.py    21 tests (core screening Z_eff)
  test_rho_collapse.py      9 tests (ρ-collapse cache)
  test_z_eff_injection.py   7 tests (Z_eff injection)
  test_lih_composed.py      12 tests (LiH composed solver)
  test_composed_diatomic.py 22 tests (general composed solver)
  test_ab_initio_pk.py      24 tests (ab initio PK)
  test_ab_initio_pk_v2.py   7 tests (PK v2 validation)
  test_gaussian_reference.py  23 tests (Gaussian baseline integrals, cc-pVDZ, cc-pVTZ)
  test_vqe_benchmark.py     12 tests (VQE pipeline, system builders, metrics)
  test_measurement_grouping.py 20 tests (QWC grouping correctness)
  test_trotter_bounds.py    16 tests (1-norm, Trotter steps)
  test_composed_qubit.py    25 tests (composed LiH qubit Hamiltonian)
  test_composed_beh2.py     14 tests (composed BeH₂ qubit Hamiltonian)
  test_composed_h2o.py      17 tests (composed H₂O qubit Hamiltonian)
  test_composed_water.py    21 tests (H₂O composed PES solver)
  test_bond_angle.py        29 tests (bond pair charge asymmetry diagnostic)
  test_lone_pair.py         13 tests (lone pair solver)
  test_1rdm_exchange.py     13 tests (full 1-RDM inter-fiber exchange)
  test_paper14_revision.py  5 tests (Paper 14 revision validation)
  test_algebraic_angular.py 28 tests (algebraic angular solver)
  test_algebraic_angular_sturmian.py 12 tests (Sturmian variant)
  test_algebraic_coupled_channel.py  10 tests (coupled-channel integration)
  test_cusp_angular_basis.py  Cusp angular basis (Track Y negative result)
  test_cusp_correction_2d.py  2D solver cusp validation (Track X2)
  test_cusp_correction_l4.py  Level 4 cusp correction tests
  test_composed_cusp.py       Composed geometry cusp wiring (Track AB)
  test_geometric_elevation.py Geometric elevation tests (Track Z)

demo/                   Demonstrations (H2, spectroscopy, AIMD)
debug/                  Development scratchpad + generated data/plots
benchmarks/             Performance tracking (PES, scaling, PySCF comparison)
ADSCFT/                 AdS/CFT correspondence research (retained, tested)
```

---

## Paper Series

| # | Title | Key Result |
|:-:|-------|------------|
| **2** | **Fine Structure Constant** | **α from Hopf bundle, 8.8×10⁻⁸, p = 5.2×10⁻⁹, universal B_formal/N = d, Hopf generalization negative, zero params** |
| 0 | Geometric Packing | Universal constant K = -1/16 |
| 1 | Spectral Graph Theory | Eigenvalue methods, O(N) scaling, Berry phase correction |
| 6 | Quantum Dynamics | Rabi, spectroscopy, AIMD at O(V) |
| 7 | **Dimensionless Vacuum** | **S3 proof (18/18 symbolic), SO(3N) generalization** |
| 8-9 | Bond Sphere + Sturmian | SO(4) D-matrix, negative theorem |
| 10 | Nuclear Lattice | Vibration/rotation graph structures |
| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H2+ 0.70% |
| **12** | **Algebraic V_ee** | **Neumann expansion, H2 92.4% D_e** |
| **13** | **Hyperspherical Lattice** | **He 0.05%, SO(6) Casimir, algebraic structure, Li SO(9)** |
| **14** | **Qubit Hamiltonians** | **O(Q^3.15) atoms, O(Q^2.5) composed; 51×–1,712× vs Gaussian** |
| **15** | **Level 4 Geometry** | **H₂ 96.0%, HeH⁺ 93.1%, 2D solver + cusp** |
| **16** | **Chemical Periodicity** | **μ_free = 2(N-2)², S_N irreps, periodic law** |
| **17** | **Composed Geometries** | **LiH R_eq 5.3%, BeH₂ 11.7%, H₂O 26%, ab initio PK, zero molecular fitting** |
| **18** | **Exchange Constants** | **Weyl--Selberg identification of κ, e^a E₁(a), μ(R); α connection** |
| FCI-A | Full CI (Atoms) | He 0.35%, Li 1.10%, Be 0.90% |
| FCI-M | LCAO FCI (Molecules) | LiH D_e 1.0% (CP-corrected) |
| **19** | **Coupled Composition** | **Balanced coupled: 0.20% energy, PK-free, 3-molecule census** |
| **20** | **Resource Benchmarks** | **30 molecules, FCI PES, 2.7×–190× Pauli advantage vs Gaussian, d-orbital sparsity** |

---

## Scope and Limitations

### What GeoVac Does Well
- **Structurally sparse qubit Hamiltonians:** O(Q^2.5) Pauli scaling for composed molecules (vs O(Q^3.9-4.3) Gaussian); 51x-1,712x fewer Pauli terms across LiH/BeH2/H2O; O(Q^1.69) 1-norm for fault-tolerant simulation cost
- **Block-diagonal electron repulsion integrals:** Gaunt selection rules enforce basis-intrinsic sparsity, compatible with all downstream optimizations (tapering, grouping, tensor factorization)
- **d-Orbital sparsity advantage:** d-orbital blocks have 4.0% ERI density (vs 8.9% for s/p), producing fewer Pauli terms per qubit. Transition metal hydrides (ScH, TiH) are cheaper to encode than main-group molecules at the same qubit count
- Classical validation benchmarks: H2 96.0% D_e (Level 4), LiH R_eq 5.3% (Level 5), BeH2 R_eq 11.7% (exchange coupling), He 0.05% (hyperspherical)
- Ab initio molecular spectroscopy (electron lattice -> PES -> rovibrational lines)
- O(V) scaling for all single-particle operations
- Core-valence diatomics (LiH, BeH+) via composed Level 3+4 geometry with ab initio pseudopotential
- Time-dependent dynamics via sparse Crank-Nicolson propagation

### Current Limitations
- **Core-valence diatomics:** The composed geometry achieves R_eq within 5.3% for LiH with l-dependent ab initio pseudopotential at l_max=2 (Paper 17). The l_max divergence is now understood to be an adiabatic approximation artifact (v2.0.6 diagnostic): bare HeH⁺ Level 4 drifts at +0.262 bohr/l_max with no PK or Z_eff, and PK-induced symmetry breaking adds +0.303 for composed LiH. Three mitigation attempts failed (algebraic PK projector, enhanced Z_eff, single-channel DBOC). The variational 2D solver (Paper 15) is the identified fix; integration into the composition pipeline is the next milestone.
- **H2 electron-electron cusp:** Level 4 molecule-frame hyperspherical solver recovers 96.0% D_e at l_max=6 with Schwartz cusp correction (Paper 15), surpassing prolate spheroidal CI (92.4%, Paper 12). CBS extrapolation estimates 96-97% with frozen π channels; reaching >99% requires higher π channels (m_max>=2) or δ channels.
- **Basis convergence:** At nmax=3, BSSE (0.115 Ha) exceeds LiH experimental binding energy (0.092 Ha). Convergence at larger nmax not characterized.
- **Transition metals:** ScH and TiH demonstrated with [Ar] frozen cores and d-orbital valence blocks. General transition metal classification (Z=23-30) is not yet automated. d-Orbital blocks are structurally sparser than s/p (4.0% vs 8.9% ERI density).
- **Polyatomics:** The composition pattern (core + bond pairs + lone pairs) produces PES for BeH₂ (`composed_triatomic.py`, R_eq 11.7%) and H₂O (`composed_water.py`, R_eq 26%). Qubit Hamiltonians with O(Q^2.5) Pauli scaling in `composed_qubit.py`.
- **Polyatomic accuracy (BeH₂):** Full 1-RDM exchange reduces R_eq error from 31% (block-diagonal) to 11.7%, matching the fitted model with zero free parameters. Residual error attributed to basis truncation (l_max=2) and adiabatic approximation.
- **Triatomic accuracy (H₂O):** R_eq = 1.34 bohr (26% error vs 1.809 expt). The bottleneck is the Level 4 angular basis at 6:1 charge asymmetry (O–H bonds), not the coupling framework. Bond-bond coupling is validated (~0.5 Ha, consistent with BeH₂); lone pair coupling is unphysical at Z_eff=6 (S·F⁰ breakdown). The Pauli scaling advantage (Q^2.5, Paper 14) is unaffected.

### What GeoVac Does NOT Replace
- Production quantum chemistry for general molecules
- Gaussian basis integral technology for large systems
- Density functional theory for extended systems

---

## Citation

```
@software{geovac2026,
  author = {J. Loutey},
  title = {GeoVac: Computational Quantum Chemistry via Spectral Graph Theory},
  year = {2026},
  version = {2.0.0},
  url = {https://github.com/jloutey-hash/geovac}
}
```

### Acknowledgment for Papers

If adapting the acknowledgment for academic papers in this series:

> Computational implementation and documentation drafting were performed with the assistance of large language models (Anthropic Claude) under the author's scientific direction. All results were validated against analytical benchmarks and symbolic proof suites.

---

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Acknowledgments

- **Spectral Graph Theory** foundations (Chung, 1997)
- **NIST Atomic Spectra Database** for validation data
- **SciPy/NumPy** for sparse matrix infrastructure

**Contact:** Issues and contributions welcome at [https://github.com/jloutey-hash/geovac/issues](https://github.com/jloutey-hash/geovac/issues)
