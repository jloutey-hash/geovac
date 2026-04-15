# GeoVac: Structurally Sparse Qubit Hamiltonians from Graph Theory

![Status](https://img.shields.io/badge/Status-Production-brightgreen) ![Version](https://img.shields.io/badge/Version-2.9.0-blue) ![License](https://img.shields.io/badge/License-MIT-orange)

GeoVac constructs **structurally sparse qubit Hamiltonians** for molecular quantum simulation. The angular momentum selection rules of the hyperspherical harmonic basis enforce block-diagonal electron repulsion integrals, producing Hamiltonians with **O(Q^2.5) Pauli term scaling** — a **51x to 1,712x advantage** over published Gaussian baselines across LiH, BeH₂, and H₂O (Paper 14).

```python
from geovac.ecosystem_export import hamiltonian
h = hamiltonian('LiH', R=3.015, l_max=2)  # 334 Pauli terms, 30 qubits
op_qiskit = h.to_qiskit()    # SparsePauliOp
op_of = h.to_openfermion()   # QubitOperator
op_pl = h.to_pennylane()     # qml.Hamiltonian
```

Install the standalone Hamiltonian package: `pip install geovac-hamiltonians`

---

## Headline Numbers

| Metric | Value |
|:-------|:------|
| Molecules | **38** (H₂ through C₂H₆, 3 periodic table rows, 10 transition metals); **35** via `hamiltonian()` API |
| Pauli scaling | **O(Q^2.5)** composed, universal coefficient 11.10 × Q |
| Advantage vs Gaussian | **51×–1,712×** fewer Pauli terms (LiH/BeH₂/H₂O) |
| 1-norm (LiH) | **33.3 Ha** (matches STO-3G 34.3 Ha, 13× fewer QWC groups) |
| He accuracy | **0.019%** (2D variational + self-consistent cusp, zero free parameters) |
| H₂ accuracy | **96.0% D_e** (molecule-frame hyperspherical, l_max=6) |
| Algebraic integrals | Exact hypergeometric R^k evaluator eliminates grid quadrature for Slater integrals |
| Mathematical foundation | Fock 1935 S³ conformal equivalence, **18 symbolic proofs** |

---

## Why This Matters for Quantum Computing

The Pauli term count and 1-norm are the dominant cost factors for near-term (VQE/NISQ) and fault-tolerant (QPE) quantum simulation respectively. GeoVac's basis-intrinsic sparsity — from Gaunt selection rules, not post-hoc optimization — produces qubit Hamiltonians that are structurally cheaper to simulate than Gaussian-basis alternatives at the same qubit count. The sparsity is compatible with all downstream optimizations (tapering, grouping, tensor factorization). See Paper 14 for the encoding theory and Paper 20 for resource benchmarks.

**Universal vs Coulomb-specific (Paper 22):** The angular sparsity guarantees are universal across spherical fermion systems — they hold for Coulomb, harmonic oscillator, Woods-Saxon, and any other radial potential, depending only on l_max. ERI density at l_max=3 is verified at 1.44% regardless of V(r). The S³ conformal projection and Hopf bundle structure are Coulomb-specific (by the Fock rigidity theorem, Paper 23). The sparsity extends to nuclear shell model Hamiltonians (Paper 23: deuteron 16 qubits / 592 Pauli, He-4 16 qubits / 712 Pauli); the conformal machinery does not.

**The HO has its own discretization (Paper 24):** The 3D harmonic oscillator has a discrete graph encoding on the holomorphic sector of S⁵ via the Bargmann-Segal transform, parallel to Fock's discretization of the Coulomb problem on S³. The Bargmann graph is **bit-exactly π-free** in exact rational arithmetic at every finite N_max. Calibration π is Coulomb-specific, not a generic feature of quantum discretization.

---

## Quantum Encoding Summary

| Metric | GeoVac (Composed) | Gaussian Baselines | Advantage |
|--------|:-----------------:|:------------------:|:---------:|
| Pauli scaling (molecules) | O(Q^2.5) | O(Q^3.9-4.3) | Structural |
| LiH (Q=30) | 334 terms | 17,065 (Trenev) | 51× |
| BeH₂ (Q=50) | 556 terms | 256,000 (Trenev) | 460× |
| H₂O (Q=70) | 778 terms | 1.33M (Trenev) | 1,712× |
| He 1-norm (Q=28) | 78.4 Ha | 530.5 Ha (cc-pVTZ) | 6.8× |

### Market Test (Track CA): GeoVac vs Computed Gaussian Baselines

| Metric | GeoVac LiH (Q=30) | Gaussian STO-3G (Q=12) | Advantage |
|--------|:-----------------:|:---------------------:|:---------:|
| Pauli terms | 334 | 907 | 2.7× |
| QWC groups | 21 | 273 | 13× |
| 1-norm (λ) | 33.3 Ha | 34.3 Ha | 0.97× (match) |

---

## Molecule Library (38 systems)

| Row | Molecules | Q range | Pauli range |
|:---:|:----------|:-------:|:----------:|
| 1st (Z=1-10) | H₂, He, LiH, BeH₂, H₂O, HF, NH₃, CH₄ | 10-90 | 112-1,000 |
| 2nd (Z=11-18) | NaH, MgH₂, HCl, H₂S, PH₃, SiH₄ | 20-80 | 239-7,273 |
| 3rd (Z=19-36) | KH, CaH₂, GeH₄, AsH₃, H₂Se, HBr | 20-80 | 239-7,273 |
| Multi-center | LiF, CO, N₂, F₂, NaCl, CH₂O, C₂H₂, C₂H₆ | 50-160 | 556-1,777 |
| Transition metals | ScH, TiH, VH, CrH, MnH, FeH, CoH, NiH, CuH, ZnH | 30 | 277 |

Isostructural invariance: molecules with the same block topology produce identical Pauli counts regardless of atomic species. The universal composed coefficient is N_Pauli = 11.11 × Q. Transition metal hydrides fall *below* this coefficient at Pauli/Q = 9.23, reflecting sparser d-orbital ERIs (4.0% density vs 8.9% for s/p).

### Nuclear Systems (Paper 23)

| System | Qubits | Pauli terms | 1-norm | Notes |
|:-------|:------:|:-----------:|:------:|:------|
| Deuteron (1p+1n) | 16 | 592 | 227 MeV | Minnesota potential, two-species JW |
| He-4 (2p+2n) | 16 | 712 | 557 MeV | Same Q, 1.20× Pauli for 12.25× Hilbert |
| He-4 + Coulomb | 16 | 712 | 552 MeV | pp Coulomb +0.7 MeV |
| Composed deuterium (1p+1n+1e) | 26 | 614 | 5e5 Ha | Nuclear + electronic block architecture |

---

## Classical Validation Benchmarks

| System | Method | Result | Paper |
|:-------|:-------|:-------|:-----:|
| He (2e) | 2D variational + self-consistent cusp | **0.019%** error | 13 |
| He (2e) | Graph-native CI (n_max=9) | **0.20%** (zero parameters, exact algebraic integrals) | 7, 13 |
| H₂ (2e) | Mol-frame hyperspherical | **96.0% D_e** | 15 |
| HeH⁺ (2e) | Charge-center hyperspherical | **93.1% D_e** | 15 |
| LiH (4e) | Composed Level 3+4, ab initio PK | **R_eq 5.3%** | 17 |
| BeH₂ (6e) | Composed + exchange coupling | **R_eq 11.7%** | 17 |
| H₂O (10e) | Composed 5-block | **R_eq 26%** | 17 |
| LiH (4e) | Balanced coupled, n_max=3 | **0.20% energy** | 19 |

---

## Development Methodology

GeoVac is developed using an **AI-augmented agentic research workflow**. The principal investigator (J. Loutey) provides scientific direction, physical intuition, and quality control. Implementation, numerical exploration, and documentation drafting are performed collaboratively with large language models (Anthropic Claude). All physics results are validated against known analytical solutions, NIST reference data, and a symbolic proof suite (18 independent tests verifying the S³ conformal geometry). No result is accepted on the basis of LLM output alone.

This workflow is itself a research contribution — an experiment in whether agentic AI tools can accelerate independent scientific research. The project has no institutional affiliation. Primary dissemination is via GitHub and Zenodo (DOI-stamped releases).

---

## Quick Start

### Installation
```bash
git clone https://github.com/jloutey-hash/geovac.git
cd geovac
pip install -e .
```

### Qubit Hamiltonians (Main API)
```python
from geovac.ecosystem_export import hamiltonian

# Any of the 35 API-accessible molecules (38 total library)
h = hamiltonian('LiH', R=3.015, l_max=2)
print(f"Qubits: {h.n_qubits}, Pauli terms: {h.n_terms}")
op = h.to_qiskit()  # or .to_openfermion(), .to_pennylane()

# Transition metal hydrides (new in v2.8)
h_fe = hamiltonian('FeH')
print(f"FeH: Q={h_fe.n_qubits}, Pauli={h_fe.n_terms}")  # Q=30, 277 Pauli
```

### Atomic Calculations
```python
from geovac import LatticeIndex

# Helium Full CI
idx = LatticeIndex(n_electrons=2, max_n=5, nuclear_charge=2,
                   vee_method='slater_full', h1_method='hybrid')
E, psi = idx.compute_ground_state()
print(f"He FCI: {E[0]:.6f} Ha")  # -2.894 Ha (0.35% error)
```

### Precision He (2D Variational)
```python
from geovac.level3_variational import solve_he_variational_2d

result = solve_he_variational_2d(Z=2, n_basis_R=25, n_basis_alpha=40, l_max=4)
print(f"He: {result['energies'][0]:.6f} Ha, error: {result['error_pct']:.4f}%")
```

---

## Paper Series

| # | Title | Key Result |
|:-:|-------|------------|
| 0 | Geometric Packing | Universal constant K = -1/16 |
| 1 | Spectral Graph Theory | Eigenvalue methods, O(N) scaling |
| **7** | **Dimensionless Vacuum** | **S³ proof (18/18 symbolic), SO(3N) generalization** |
| **14** | **Qubit Hamiltonians** | **O(Q^2.5) composed; 51×–1,712× vs Gaussian** |
| **16** | **Chemical Periodicity** | **S_N representation theory, atomic classifier** |
| 6 | Quantum Dynamics | Rabi, spectroscopy, AIMD at O(V) |
| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H₂⁺ 0.0002% |
| **12** | Algebraic V_ee | Neumann expansion, H₂ 92.4% D_e |
| **13** | Hyperspherical Lattice | He 0.004%, fiber bundle, algebraic structure |
| **15** | Level 4 Geometry | H₂ 96.0% D_e, HeH⁺ 93.1% |
| **17** | Composed Geometries | LiH 5.3%, BeH₂ 11.7%, ab initio PK |
| **18** | Exchange Constants | Weyl-Selberg taxonomy, π-free graph principle |
| **22** | **Angular Sparsity Theorem** | **Potential-independent ERI density (1.44% at l_max=3)** |
| **24** | **Bargmann-Segal Lattice** | **π-free HO discretization on S⁵ Hardy space** |
| FCI-A | Full CI (Atoms) | He 0.35%, Li 1.10%, Be 0.90% |
| **19** | Balanced Coupled | 0.20% energy, PK-free, 3-molecule census |
| **20** | Resource Benchmarks | 38 molecules, Gaussian comparison, `pip install` |
| **23** | **Nuclear Shell Hamiltonians** | **Deuteron/He-4 qubit encoding, Fock rigidity theorem** |
| 21 | Synthesis | S³ proof chain, exchange constants, research program |
| 2 | Fine Structure Constant | α from Hopf bundle, 8.8×10⁻⁸ (conjectural) |
| 8-9 | Bond Sphere + Sturmian | Structural theorem (guardrail), SO(4) selection rules |

---

## Scope and Limitations

### What GeoVac Does Well
- **Structurally sparse qubit Hamiltonians:** O(Q^2.5) Pauli scaling, 51×–1,712× fewer terms vs Gaussian
- **Block-diagonal ERIs:** Gaunt selection rules enforce basis-intrinsic sparsity
- **d-Orbital sparsity:** d-blocks have 4.0% ERI density (vs 8.9% s/p) — transition metals are cheaper per qubit
- **Classical benchmarks:** H₂ 96.0% D_e, LiH R_eq 5.3%, He 0.019% (self-consistent cusp)
- **Algebraic Slater integrals:** Exact hypergeometric R^k evaluator eliminates grid quadrature (machine precision, 8x speedup)
- **DUCC downfolding:** Identifies PK l_max divergence root cause (109x p-orbital underestimate); H₂O 1-norm 9% lower
- **Zero molecular fitting parameters** — all from nuclear charges and geometry

### Current Limitations
- **Classical PES accuracy:** R_eq errors of 5-26% for composed molecules (PK structural overcounting). Optimal for fixed-geometry quantum simulation, not geometry optimization.
- **Transition metals:** Full first series implemented as hydrides (v2.8.0); non-hydride TM molecules not yet built
- **Polyatomic accuracy:** H₂O R_eq 26% (charge asymmetry bottleneck, not coupling framework)

### What GeoVac Does NOT Replace
- Production quantum chemistry for general molecules
- Gaussian basis integral technology for large systems
- Density functional theory for extended systems

---

## Citation

```
@software{geovac2026,
  author = {J. Loutey},
  title = {GeoVac: Structurally Sparse Qubit Hamiltonians from Graph Theory},
  year = {2026},
  version = {2.8.2},
  url = {https://github.com/jloutey-hash/geovac}
}
```

If adapting the acknowledgment for academic papers:

> Computational implementation and documentation drafting were performed with the assistance of large language models (Anthropic Claude) under the author's scientific direction. All results were validated against analytical benchmarks and symbolic proof suites.

---

## License

MIT License - See [LICENSE](LICENSE) for details.

**Contact:** Issues and contributions welcome at [https://github.com/jloutey-hash/geovac/issues](https://github.com/jloutey-hash/geovac/issues)
