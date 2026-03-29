# GeoVac: Computational Quantum Chemistry via Spectral Graph Theory

![Status](https://img.shields.io/badge/Status-Production-brightgreen) ![Version](https://img.shields.io/badge/Version-2.0.6-blue) ![License](https://img.shields.io/badge/License-MIT-orange)

**Version 2.0.6** - Algebraic Angular Solver & Coupled-Channel Integration (Level 3)

GeoVac models quantum mechanics on **discrete, dimensionless graph topologies**. The discrete graph Laplacian -- a sparse matrix with O(V) nonzero entries -- is *mathematically equivalent* to the Schrodinger equation for hydrogen via Fock's 1935 conformal projection, as formally proven via 18 independent symbolic proofs (Paper 7). This equivalence is the computational foundation: by working directly on the graph topology, expensive continuous integration is replaced by O(N) sparse matrix eigenvalue problems that produce the same physics.

For molecules, the natural geometry shifts from $S^3$ to prolate spheroidal coordinates, where the two-center Coulomb problem separates exactly. The **natural geometry principle** -- *every separable quantum system has a natural lattice* -- unifies the atomic and molecular solvers (Paper 11).

---

## Development Methodology

GeoVac is developed using an AI-augmented research workflow. The principal investigator (J. Loutey) provides scientific direction, physical intuition, and quality control. Implementation, numerical exploration, and documentation drafting are performed collaboratively with large language models (Anthropic Claude). All physics results are validated against known analytical solutions, NIST reference data, and a symbolic proof suite (18 independent tests verifying the S³ conformal geometry). No result is accepted on the basis of LLM output alone.

This workflow is itself a research contribution — an experiment in whether agentic AI tools can accelerate independent scientific research. The project has no institutional affiliation. Primary dissemination is via GitHub and Zenodo (DOI-stamped releases).

---

## What's New in v2.0.6

Algebraic spectral basis replaces finite-difference angular solver at Level 3 (He). Coupled-channel integration resolves the l_max convergence problem.

- **Algebraic angular solver:** Gegenbauer spectral basis reduces angular matrix dimension from 200-800 (FD) to 10-30 (algebraic). All matrix elements algebraic (SO(6) Casimir, Gaunt integrals, partial harmonic sums). Exact Hellmann-Feynman P-matrix from R-independent dH/dR.
- **Coupled-channel convergence:** Single-channel error diverges with l_max (0.16% → 0.65%); coupled-channel error converges (0.37% → 0.27% at l_max 1-3). P-only and P+Q bracket exact energy from below and above.
- **DBOC analysis:** +0.035 Ha correction, 23x the adiabatic error, with 97% cancelled by off-diagonal coupling. Confirms adiabatic approximation (not angular basis) is the accuracy bottleneck.
- **Sturmian variant:** 1.6-3.6x faster angular convergence at l_max=0; counterproductive at l_max≥1 (converges faster to wrong adiabatic limit). Preserved as research artifact.
- **Paper 13:** New section on algebraic angular solver and coupled-channel convergence

### Prior Releases
- **v2.0.5:** Asymmetric Bond Diagnostic Arc — l_max divergence diagnosis (+0.23 bohr/l_max), R-dependent PK scaling (2.0% at l_max=4), negative results (eigenchannel rotation, spheroidal compression, self-consistent PK, projected PK)
- **v2.0.4:** H₂O Composed Solver — 5-block architecture, R_eq 26%, charge asymmetry diagnostic, lone pair coupling (negative result)
- **v2.0.3:** Commutator-Based Trotter Bounds — r ~ Q^{1.47}, 7× fewer steps at Q=60
- **v2.0.2:** Full 1-RDM Exchange — BeH₂ R_eq 11.7%, algebraic methods, l-dependent PK (LiH 5.3%)
- **v2.0.1:** BeH₂ Polyatomic Diagnostic — inter-fiber exchange coupling (20% R_eq), failed approaches documented
- **v2.0.0:** Composed Multi-Center Qubit Hamiltonians — universal Q^2.5 Pauli scaling, 51×–1,712× Gaussian advantage
- **v1.9.0:** VQE Benchmark Infrastructure — validated Gaussian baselines, equal-qubit comparison, 35 tests
- **v1.8.0:** Quantum Simulation Cost Analysis — O(Q^1.69) Pauli 1-norm, measurement/simulation cost modules, 48 tests
- **v1.7.7:** Epistemic calibration & rhetorical tightening (Papers 0, 7, 16, FCI-A, terminology audit)
- **v1.7.6:** Paper 2 universal algebraic identity B_formal/N = d, quaternionic/octonionic Hopf negative results, residual analysis
- **v1.7.4:** Paper 2 algebraic selection principle (B/N = 3(m+2)(m-1)/10), S³ specificity proof, spectral geometry survey (negative results)
- **v1.7.3:** Paper 2 rewrite — statistical validation (p = 5.2×10⁻⁹), circulant Z₃ structure, derivation chain assessment
- **v1.7.2:** Documentation review & epistemic tightening (Papers 0, 1, 16)
- **v1.7.1:** Paper 2 rewrite — α from Hopf bundle (8.8×10⁻⁸, zero params, `hopf_bundle.py`)
- **v1.7.0:** Composed Natural Geometries (Paper 17, LiH R_eq 6.4%, BeH⁺ bound)
- **v1.6.1:** Hierarchical Molecular Solvers (FrozenCoreLatticeIndex, LockedShellMolecule)
- **v1.6.0:** Chemical Periodicity from S_N Representation Theory (Paper 16)
- **v1.5.0:** Algebraic Structure & SO(3N) Generalization (Papers 7 & 13)
- **v1.4.0:** Papers 14-15, Level 4 Solver (H₂ 94.1%, HeH⁺ 93.1%)
- **v1.2.0:** QA & Corrections (Berry phase retraction, autoionization sections)
- **v1.1.0:** Multi-Particle Natural Geometry (Papers 12-13, He 0.05%, H$_2$ 92.4% $D_e$)
- **v1.0.x:** Paper 11, Direct CI, LiH benchmark, Bond Sphere theory

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
| **H$_2$ (Level 4, 2D)** | **Mol.-frame hypersp.** | **94.1% $D_e$** | **Paper 15** |
| **HeH$^+$ (Level 4)** | **Charge-center hypersp.** | **93.1% $D_e$** | **Paper 15** |
| H$_2$ (Full CI) | Tensor product | 1.72% | Cross-nuclear ceiling |
| **LiH (composed, ab initio PK)** | **Composed Level 3+4** | **R_eq: 5.3%** | **Paper 17** |
| **BeH⁺ (composed, ab initio PK)** | **Composed Level 3+4** | **Bound** | **Paper 17** |
| **BeH₂ (composed, exchange)** | **Composed Level 3+4, exchange coupling** | **R_eq: 11.7%** | **Paper 17** |
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
  composed_qubit.py       Composed multi-center qubit Hamiltonians (LiH, BeH₂, H₂O)
  hopf_bundle.py          Hopf fibration analysis (S³→S², α derivation)
  wigner_so4.py           SO(4) Wigner D-matrix elements
  atomic_solver.py        Single-electron atom solver
  dirac_hamiltonian.py    Relativistic Dirac solver

papers/
  core/                 Defensible foundations (Papers 0, 1, 6-15, 17, FCI)
    Paper 0:  Geometric packing & universal constant
    Paper 1:  Spectral graph theory & eigenvalue methods
    Paper 6:  Quantum dynamics & thermodynamics
    Paper 7:  The Dimensionless Vacuum (S3 proof, 18/18 proofs)
    Paper 8-9: Bond Sphere & Sturmian negative theorem, SO(4) selection rules
    Paper 10: Nuclear lattice for vibration/rotation
    Paper 11: Molecular Fock Projection (prolate spheroidal)
    Paper 12: Algebraic V_ee (Neumann expansion)
    Paper 13: Hyperspherical lattice (multi-electron atoms)
    Paper 14: Structurally sparse qubit Hamiltonians
    Paper 15: Level 4 mol.-frame hyperspherical (H2 94.1%, HeH+ 93.1%)
    Paper 17: Composed natural geometries (LiH 5.3%, BeH⁺)
    FCI-A:    Multi-electron atoms (He, Li, Be)
    FCI-M:    Heteronuclear diatomics (LiH benchmark)
  observations/         Empirical observations & pattern recognition (Paper 16)
    Paper 16: Chemical periodicity from S_N representation theory
  conjectures/          Theoretical explorations (Papers 2-5)
    Paper 2:  Fine structure α from Hopf bundle (8.8×10⁻⁸ error)
  archive/              Superseded paper versions and diagnostic reports

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
| **15** | **Level 4 Geometry** | **H₂ 94.1%, HeH⁺ 93.1%, 2D solver** |
| **16** | **Chemical Periodicity** | **μ_free = 2(N-2)², S_N irreps, periodic law** |
| **17** | **Composed Geometries** | **LiH R_eq 5.3%, BeH₂ 11.7%, H₂O 26%, ab initio PK, zero molecular fitting** |
| FCI-A | Full CI (Atoms) | He 0.35%, Li 1.10%, Be 0.90% |
| FCI-M | LCAO FCI (Molecules) | LiH D_e 1.0% (CP-corrected) |

---

## Scope and Limitations

### What GeoVac Does Well
- Single-center systems (atoms, atomic ions) via S3 graph Laplacian
- Two-electron atoms (He) via hyperspherical adiabatic solver (0.05% error)
- One-electron diatomics (H2+, HeH2+) via prolate spheroidal lattice
- Two-electron diatomics (H2) via Neumann algebraic V_ee (92.4% D_e)
- Few-electron atoms (2-4 electrons) via sparse Full CI
- Ab initio molecular spectroscopy (electron lattice -> PES -> rovibrational lines)
- O(V) scaling for all single-particle operations
- Efficient qubit encoding: O(Q^3.15) Pauli terms (atoms), O(Q^2.5) (composed molecules); 51×–1,712× fewer Pauli terms vs published Gaussian baselines
- Core-valence diatomics (LiH, BeH⁺) via composed Level 3+4 geometry with ab initio pseudopotential
- Time-dependent dynamics via sparse Crank-Nicolson propagation

### Current Limitations
- **Core-valence diatomics:** The composed geometry achieves R_eq within 5.3% for LiH with l-dependent ab initio pseudopotential at l_max=2 (Paper 17). The l_max divergence is now understood to be an adiabatic approximation artifact (v2.0.6 diagnostic): bare HeH⁺ Level 4 drifts at +0.262 bohr/l_max with no PK or Z_eff, and PK-induced symmetry breaking adds +0.303 for composed LiH. Three mitigation attempts failed (algebraic PK projector, enhanced Z_eff, single-channel DBOC). The variational 2D solver (Paper 15) is the identified fix; integration into the composition pipeline is the next milestone.
- **H2 electron-electron cusp:** Level 4 molecule-frame hyperspherical solver recovers 94.1% D_e (Paper 15), surpassing prolate spheroidal CI (92.4%, Paper 12). The remaining ~6% gap requires higher partial waves ($l_{\max} \ge 5$).
- **Basis convergence:** At nmax=3, BSSE (0.115 Ha) exceeds LiH experimental binding energy (0.092 Ha). Convergence at larger nmax not characterized.
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
