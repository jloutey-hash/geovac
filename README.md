# GeoVac: Computational Quantum Chemistry via Spectral Graph Theory

![Status](https://img.shields.io/badge/Status-Production-brightgreen) ![Version](https://img.shields.io/badge/Version-1.1.0-blue) ![License](https://img.shields.io/badge/License-MIT-orange)

**Version 1.1.0** - The Multi-Particle Release

GeoVac models quantum mechanics on **discrete, dimensionless graph topologies**. The discrete graph Laplacian -- a sparse matrix with O(V) nonzero entries -- is *mathematically equivalent* to the Schrodinger equation for hydrogen via Fock's 1935 conformal projection, as formally proven via 18 independent symbolic proofs (Paper 7). This equivalence is the computational foundation: by working directly on the graph topology, expensive continuous integration is replaced by O(N) sparse matrix eigenvalue problems that produce the same physics.

For molecules, the natural geometry shifts from $S^3$ to prolate spheroidal coordinates, where the two-center Coulomb problem separates exactly. The **natural geometry principle** -- *every separable quantum system has a natural lattice* -- unifies the atomic and molecular solvers (Paper 11).

---

## What's New in v1.1.0

### The Multi-Particle Natural Geometry (Paper 13)
The hyperspherical lattice extends GeoVac to multi-electron systems:

| System | Method | Error | Matrix dim. |
|--------|--------|:-----:|:-----------:|
| He (hyperspherical, 1 ch.) | Adiabatic | **0.05%** | 600 |
| H$_2$ (Neumann $V_{ee}$, 27 bf) | Algebraic CI | **92.4%** $D_e$ | 27 |
| H$_2$ $\nu_{01}$ (ab initio) | Full pipeline | **+1.3%** | --- |

### Paper 12: Algebraic Two-Electron Integrals
- Neumann expansion replaces numerical $V_{ee}$ quadrature with algebraic evaluation
- 12-20 percentage point improvement over numerical integration at all basis sizes
- Remaining 7.6% gap diagnosed as electron-electron cusp (basis limit)

### Paper 13: Hyperspherical Lattice for He
- Single-channel adiabatic solver: E = -2.9052 Ha (0.05% vs exact -2.9037 Ha)
- First non-trivial fiber bundle in GeoVac (angular structure varies with R)
- Structural prototype for electron-nuclear coupling in molecules
- 20 new tests, all passing

### Ab Initio Molecular Spectroscopy
- Full pipeline: electron lattice -> PES -> Morse fit -> nuclear lattice -> rovibrational spectrum
- H$_2$: $R_e$ = 1.42 bohr (+1.5%), $\omega_e$ = 4500 cm$^{-1}$ (+2.2%), $B_e$ = 59.1 cm$^{-1}$ (-2.9%)
- Zero experimental spectroscopic input at any stage

### Qubit Encoding Benchmarks
- GeoVac Pauli terms scale as O(Q^3.15) vs O(Q^4.60) for conventional Gaussian bases
- ERI density drops as ~1/M^2 (angular momentum selection rules)

### Prior Milestones (v1.0.x, v0.9.x)
- **Paper 11:** Prolate spheroidal lattice for H$_2^+$ (0.70% error, zero free params)
- **Direct CI:** Excitation-driven FCI solver, O(N_SD^1.0) scaling, Be nmax=4 in 357s
- **LiH benchmark:** D_e_CP = 0.093 Ha (1.0% error), BSSE diagnosed and corrected
- **Bond Sphere theory (Paper 8-9):** SO(4) selection rules, Sturmian negative theorem

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
| H$_2$ (Full CI) | Tensor product | 1.72% | Cross-nuclear ceiling |
| LiH ($D_e$, CP-corrected) | LCAO FCI | **1.0%** | nmax=3 |

### Ab Initio Spectroscopy

| Benchmark | Ab Initio | Experiment | Error |
|-----------|-----------|------------|:-----:|
| H$_2$ $R_e$ (bohr) | 1.422 | 1.401 | +1.5% |
| H$_2$ $\omega_e$ (cm$^{-1}$) | 4500 | 4401 | +2.2% |
| H$_2$ $B_e$ (cm$^{-1}$) | 59.09 | 60.85 | -2.9% |
| H$_2$ $\nu_{01}$ (cm$^{-1}$) | 4214 | 4161 | +1.3% |

### Dynamics & Spectroscopy

| Benchmark | Result | Status |
|-----------|--------|:------:|
| Rabi oscillation | 0.41% period error | PASS |
| H$_2$ delta-kick spectroscopy | 20/35 peaks, 0.16% | PASS |
| AIMD energy conservation | 0.0003% max drift | PASS |
| Langevin bond dissociation | Correct thermal behavior | PASS |

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
| Prolate spheroidal | O(N_xi) | < 1% for 1e diatomics | Molecular PES |

---

## Project Structure

```
geovac/                 Core package
  lattice.py              GeometricLattice (quantum state graph)
  lattice_index.py        LatticeIndex / MolecularLatticeIndex (N-electron FCI)
  direct_ci.py            DirectCISolver (excitation-driven FCI)
  hamiltonian.py          MoleculeHamiltonian (molecular bonding)
  prolate_spheroidal_lattice.py   Prolate spheroidal solver (H2+)
  prolate_scf.py            Prolate spheroidal SCF + CI (H2)
  neumann_vee.py            Algebraic V_ee via Neumann expansion
  hyperspherical_angular.py Hyperangular solver (Gaunt + Liouville FD)
  hyperspherical_adiabatic.py Adiabatic potential curves
  hyperspherical_radial.py  Hyperradial solver + solve_helium()
  nuclear_lattice.py      Nuclear vibration/rotation lattice
  coupled_en_lattice.py   Coupled electronic-nuclear lattice
  dynamics.py             TimePropagator (Crank-Nicolson)
  aimd.py                 VelocityVerlet, LangevinThermostat
  wigner_so4.py           SO(4) Wigner D-matrix elements
  atomic_solver.py        Single-electron atom solver
  dirac_hamiltonian.py    Relativistic Dirac solver

papers/
  core/                 Defensible foundations (Papers 0, 1, 6, 7, 10-13, FCI)
    Paper 0:  Geometric packing & universal constant
    Paper 1:  Spectral graph theory & eigenvalue methods
    Paper 6:  Quantum dynamics & thermodynamics
    Paper 7:  The Dimensionless Vacuum (S3 proof, 18/18 proofs)
    Paper 10: Nuclear lattice for vibration/rotation
    Paper 11: Molecular Fock Projection (prolate spheroidal)
    Paper 12: Algebraic V_ee (Neumann expansion)
    Paper 13: Hyperspherical lattice (multi-electron atoms)
    FCI-A:    Multi-electron atoms (He, Li, Be)
    FCI-M:    Heteronuclear diatomics (LiH benchmark)
  zenodo/               Publication cluster (Papers 7, 11, 12, 13)
  methods/              Paper 8-9: Bond Sphere & Sturmian negative theorem
  conjectures/          Theoretical explorations (Papers 2-5; speculative)

tests/                  Unit tests (pytest) + validation
  test_fock_projection.py   10 proofs (stereographic geometry)
  test_fock_laplacian.py    8 proofs (conformal Laplacian)
  test_prolate_h2plus.py    11 tests (prolate spheroidal solver)
  test_hyperspherical_he.py 20 tests (hyperspherical He solver)
  test_neumann_vee.py       Neumann V_ee algebraic integrals
  test_direct_ci.py         6 tests (Direct CI consistency)
  test_lih_fci.py           50 tests (LiH heteronuclear FCI)

demo/                   Demonstrations (H2, spectroscopy, AIMD)
debug/                  Development scratchpad + generated data/plots
benchmarks/             Performance tracking (PES, scaling, PySCF comparison)
ADSCFT/                 AdS/CFT correspondence research (retained, tested)
```

---

## Paper Series

| # | Title | Key Result |
|:-:|-------|------------|
| 0 | Geometric Packing | Universal constant K = -1/16 |
| 1 | Spectral Graph Theory | Eigenvalue methods, O(N) scaling |
| 6 | Quantum Dynamics | Rabi, spectroscopy, AIMD at O(V) |
| 7 | **Dimensionless Vacuum** | **S3 proof (18/18 symbolic), Schrodinger recovery** |
| 8-9 | Bond Sphere + Sturmian (`methods/`) | SO(4) D-matrix, negative theorem |
| 10 | Nuclear Lattice | Vibration/rotation graph structures |
| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H2+ 0.70% |
| **12** | **Algebraic V_ee** | **Neumann expansion, H2 92.4% D_e** |
| **13** | **Hyperspherical Lattice** | **He 0.05%, fiber bundle, ab initio spectroscopy** |
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
- Efficient qubit encoding: O(Q^3.15) Pauli terms vs O(Q^4.60) conventional
- Time-dependent dynamics via sparse Crank-Nicolson propagation

### Current Limitations
- **Multi-electron molecules:** The LCAO approach achieves correct-order binding energies but lacks R-dependent kinetic repulsion. Equilibrium geometry requires a Fock-weighted correction with adjustable lambda (see FCI-M paper, Sec. V).
- **H2 electron-electron cusp:** Prolate spheroidal CI saturates at 92.4% D_e due to missing non-analytic r_12 terms (Paper 12 diagnosis). The hyperspherical approach (Paper 13) resolves this for atoms but has not yet been combined with molecular coordinates.
- **Basis convergence:** At nmax=3, BSSE (0.115 Ha) exceeds LiH experimental binding energy (0.092 Ha). Convergence at larger nmax not characterized.
- **Polyatomics:** No natural geometry identified beyond two-center systems.
- **Level 4 geometry:** The combined two-center + two-electron natural coordinate system remains an open problem.

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
  version = {1.1.0},
  url = {https://github.com/jloutey-hash/geovac}
}
```

---

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Acknowledgments

- **Spectral Graph Theory** foundations (Chung, 1997)
- **NIST Atomic Spectra Database** for validation data
- **SciPy/NumPy** for sparse matrix infrastructure

**Contact:** Issues and contributions welcome at [https://github.com/jloutey-hash/geovac/issues](https://github.com/jloutey-hash/geovac/issues)
