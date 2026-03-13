# GeoVac: Computational Quantum Chemistry via Spectral Graph Theory

![Status](https://img.shields.io/badge/Status-Production-brightgreen) ![Version](https://img.shields.io/badge/Version-1.0.0-blue) ![License](https://img.shields.io/badge/License-MIT-orange)

**Version 1.0.0** - The Natural Geometry Release

GeoVac models quantum mechanics on **discrete, dimensionless graph topologies**. The discrete graph Laplacian -- a sparse matrix with O(V) nonzero entries -- is *mathematically equivalent* to the Schrodinger equation for hydrogen via Fock's 1935 conformal projection, as formally proven via 18 independent symbolic proofs (Paper 7). This equivalence is the computational foundation: by working directly on the graph topology, expensive continuous integration is replaced by O(N) sparse matrix eigenvalue problems that produce the same physics.

For molecules, the natural geometry shifts from $S^3$ to prolate spheroidal coordinates, where the two-center Coulomb problem separates exactly. The **natural geometry principle** -- *every separable quantum system has a natural lattice* -- unifies the atomic and molecular solvers (Paper 11).

---

## What's New in v1.0.0

### The Natural Geometry Principle (Paper 11)
The prolate spheroidal lattice extends GeoVac from atoms to molecules with zero free parameters:

| System | $R_{\mathrm{eq}}$ error | $E_{\min}$ error | Free params |
|--------|:-----------------------:|:-----------------:|:-----------:|
| H$_2^+$ (prolate spheroidal) | **0.21%** | **0.70%** | **0** |
| LiH (LCAO FCI, CP-corrected) | shifted | **1.0%** $D_e$ | 0 |

### Paper Series Complete (11 papers)
- **Paper 10:** Nuclear lattice for molecular vibration and rotation
- **Paper 11:** Molecular Fock projection -- prolate spheroidal lattice for diatomics
- Forward references added across Papers 7, 8-9, FCI, LCAO-FCI, and 10

### Prior Milestones (v0.9.x)
- **Direct CI:** Excitation-driven FCI solver, O(N_SD^1.0) scaling, Be nmax=4 in 357s
- **LiH benchmark:** D_e_CP = 0.093 Ha (1.0% error), BSSE diagnosed and corrected
- **Bond Sphere theory (Paper 8-9):** SO(4) selection rules, Sturmian negative theorem
- **Numba-accelerated DirectCI, Fock-weighted kinetic correction**

---

## Precision Benchmarks

### Atoms (Graph Laplacian on $S^3$)

| System | Method | Error | Status |
|--------|--------|:-----:|:------:|
| H, He+, H$_2^+$ (1e) | Graph Laplacian | < 0.1% | Exact mean-field |
| He (2e, FCI) | slater_full + hybrid h1 | **0.35%** | Monotonic convergence |
| Li (3e, FCI) | slater_full + exact h1 | **1.10%** | Monotonic convergence |
| Be (4e, FCI) | Direct CI, nmax=4 | **0.90%** | Sub-1% |
| Li+ (2e, Z=3) | Conformal torsion | **0.039%** | Isoelectronic scaling |

### Molecules

| System | Method | Error | Status |
|--------|--------|:-----:|:------:|
| H$_2^+$ (prolate spheroidal) | Separated solver | **0.70%** | Zero free params |
| H$_2$ (prolate CI, relaxed) | Eckart SCF + 2×2 CI | 58% $D_e$ | Zero free params |
| H$_2$ (Full CI) | Tensor product | 1.72% | Cross-nuclear ceiling |
| LiH ($D_e$, CP-corrected) | LCAO FCI | **1.0%** | nmax=3 |

### Dynamics & Spectroscopy

| Benchmark | Result | Status |
|-----------|--------|:------:|
| Rabi oscillation | 0.41% period error | PASS |
| H$_2$ delta-kick spectroscopy | 20/35 peaks, 0.16% | PASS |
| AIMD energy conservation | 0.0003% max drift | PASS |
| Langevin bond dissociation | Correct thermal behavior | PASS |

### GeoVac vs PySCF (CI-validated)

| System | GeoVac | PySCF/STO-3G | PySCF/cc-pVQZ |
|--------|:------:|:------------:|:-------------:|
| H | 0.57% | 6.68% | 0.01% |
| He | **0.35%** | 3.30% | 1.45% |
| H$_2$ | 1.73% | 3.17% | 0.06% |

GeoVac beats PySCF/STO-3G on all tested atoms with no basis-set optimization.

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
  nuclear_lattice.py      Nuclear vibration/rotation lattice
  coupled_en_lattice.py   Coupled electronic-nuclear lattice
  dynamics.py             TimePropagator (Crank-Nicolson)
  aimd.py                 VelocityVerlet, LangevinThermostat
  wigner_so4.py           SO(4) Wigner D-matrix elements
  atomic_solver.py        Single-electron atom solver
  dirac_hamiltonian.py    Relativistic Dirac solver

papers/
  core/                 Defensible foundations (Papers 0, 1, 6, 7, 10, 11, FCI)
    Paper 0:  Geometric packing & universal constant
    Paper 1:  Spectral graph theory & eigenvalue methods
    Paper 6:  Quantum dynamics & thermodynamics
    Paper 7:  The Dimensionless Vacuum (S3 proof, 18/18 proofs)
    Paper 10: Nuclear lattice for vibration/rotation
    Paper 11: Molecular Fock Projection (prolate spheroidal)
    FCI-A:    Multi-electron atoms (He, Li, Be)
    FCI-M:    Heteronuclear diatomics (LiH benchmark)
  methods/              Paper 8-9: Bond Sphere & Sturmian negative theorem
  conjectures/          Theoretical explorations (Papers 2-5; speculative)

tests/                  Unit tests (pytest) + validation
  test_fock_projection.py   10 proofs (stereographic geometry)
  test_fock_laplacian.py    8 proofs (conformal Laplacian)
  test_prolate_h2plus.py    11 tests (prolate spheroidal solver)
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
| **11** | **Molecular Fock Projection** | **Prolate spheroidal lattice, H2+ 0.21%/0.70%** |
| FCI-A | Full CI (Atoms) | He 0.35%, Li 1.10%, Be 0.90% |
| FCI-M | LCAO FCI (Molecules) | LiH D_e 1.0% (CP-corrected) |

---

## Scope and Limitations

### What GeoVac Does Well
- Single-center systems (atoms, atomic ions) via S3 graph Laplacian
- One-electron diatomics (H2+, HeH2+) via prolate spheroidal lattice
- Few-electron atoms (2-4 electrons) via sparse Full CI
- O(V) scaling for all single-particle operations
- Time-dependent dynamics via sparse Crank-Nicolson propagation

### Current Limitations
- **Multi-electron molecules:** The LCAO approach achieves correct-order binding energies but lacks R-dependent kinetic repulsion. Equilibrium geometry requires a Fock-weighted correction with adjustable lambda (see FCI-M paper, Sec. V).
- **Basis convergence:** At nmax=3, BSSE (0.115 Ha) exceeds LiH experimental binding energy (0.092 Ha). Convergence at larger nmax not characterized.
- **Two-electron diatomics:** H2 via prolate spheroidal CI achieves D_e = 0.101 Ha (58% of exact) with relaxed orbitals. HeH+ binding recovered with per-atom Z_eff optimization.
- **Polyatomics:** No natural geometry identified beyond two-center systems.

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
  version = {1.0.0},
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
