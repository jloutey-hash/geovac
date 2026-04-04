# GeoVac: Structurally Sparse Qubit Hamiltonians for Molecular Quantum Simulation

## What GeoVac Does

GeoVac is a discretization framework that exploits the natural geometry of quantum systems (hyperspherical coordinates, prolate spheroidal coordinates, composed fiber bundles) to produce qubit Hamiltonians with far fewer Pauli terms than standard Gaussian-basis Jordan-Wigner encodings. The sparsity is structural -- it comes from angular momentum selection rules baked into the basis, not from post-hoc truncation or compression. The result is a direct reduction in the measurement budget for variational quantum eigensolvers (VQE) and other near-term quantum algorithms.

## The Headline Numbers

Head-to-head Pauli term counts, GeoVac composed encoding vs Gaussian Jordan-Wigner:

| System | GeoVac qubits | GeoVac Pauli terms | Gaussian basis | Gaussian Pauli terms | Advantage |
|--------|:---:|:---:|--------|:---:|:---:|
| He | 10 | 118 | cc-pVDZ (10 qubits) | 156 | 1.3x |
| He | 28 | 2,659 | cc-pVTZ (28 qubits) | 21,607 | 8.1x |
| LiH | 30 | 334 | cc-pVDZ (~30 qubits) | ~63,500 | 190x |
| BeH2 | 50 | 556 | estimated | ~200,000 | ~360x |
| H2O | 70 | 778 | estimated | ~560,000 | ~720x |

GeoVac Pauli terms scale as Q^2.5 across all composed systems (LiH, BeH2, H2O), with an exponent spread of just 0.02. Gaussian Jordan-Wigner scaling is Q^3.9-4.3.

## Why It's Sparse

GeoVac's composed architecture assigns each electron group (core pair, bond pair, lone pair) its own natural coordinate system, then couples them via screened Coulomb integrals and Phillips-Kleinman pseudopotentials. Within each block, the angular basis is built from spherical harmonics on the natural geometry, so Gaunt integral selection rules enforce exact zeros in the electron repulsion integrals. The result is block-diagonal ERIs with no cross-block two-electron terms. This structural sparsity is intrinsic to the basis and survives all downstream transformations (Jordan-Wigner encoding, qubit tapering, measurement grouping, tensor factorization).

## Where It Fits

### Near-term / VQE

For variational algorithms on near-term hardware, the dominant cost is the measurement budget: the number of distinct Pauli terms times the number of shots per term. GeoVac's 190x-720x reduction in Pauli terms translates directly to a proportional reduction in total measurements. This advantage is complementary to measurement optimization techniques (qubit-wise commuting grouping, classical shadows, derandomization) -- GeoVac reduces the number of terms those techniques have to work with.

Note on comparisons: density fitting (DF) and tensor hypercontraction (THC) are powerful Gaussian-basis compression techniques, but they target block-encoding representations for fault-tolerant algorithms. They do not reduce the number of Pauli terms in a Jordan-Wigner Hamiltonian. GeoVac's Pauli term advantage is against raw Jordan-Wigner encoding, which is the representation used by VQE and other near-term variational methods.

1-norm data (relevant for shot budgets and Trotter bounds):
- He at Q=10: 11.29 Ha (vs Gaussian cc-pVDZ 42.95 Ha, 3.8x advantage)
- He at Q=28: 78.36 Ha (vs Gaussian cc-pVTZ 530.47 Ha, 6.8x advantage)
- He atomic 1-norm scaling: Q^1.69
- LiH at Q=30: 33.26 Ha (electronic-only; PK classical correction: +4.08 Ha)
- BeH2 at Q=50: 354.89 Ha
- H2O at Q=70: 361 Ha (electronic-only; PK classical correction resolves the 28,053 Ha total)

### Fault-tolerant / QPE

GeoVac's 1-norm scaling (Q^1.69 for atoms) and commutator-based Trotter bounds (Q^1.47, 7x fewer Trotter steps at Q=60) suggest competitive resource estimates for fault-tolerant quantum phase estimation. However, block-encoding circuits for the GeoVac Hamiltonian structure have not yet been built. This is a concrete collaboration opportunity: the composed block-diagonal structure should map naturally to LCU or QROM encodings, but the circuit construction and resource estimation remain open.

## Honest Limitations

GeoVac's classical accuracy is substantially below production quantum chemistry. The best classical results are H2 at 96.0% dissociation energy, LiH equilibrium bond length at 5.3% error, BeH2 at 11.7%, and H2O at 26% -- compared to sub-0.1% for Gaussian methods with correlating basis sets. The framework currently covers five systems (H2, He, LiH, BeH2, H2O). The value proposition is structural sparsity for quantum simulation, not classical accuracy.

Note on H2O 1-norm: The Phillips-Kleinman pseudopotential (PK) is a one-body operator whose energy contribution can be computed classically from the 1-RDM measured during VQE, with zero additional quantum circuits. This quantum-classical partitioning (v2.0.29) reduces the H2O quantum Hamiltonian 1-norm from 28,053 Ha to 361 Ha (78x reduction), with no approximation. The partitioned 1-norms across all composed systems are comparable: LiH 33 Ha, BeH2 355 Ha, H2O 361 Ha.

## Try It

```bash
pip install geovac-hamiltonians
```

```python
from geovac_hamiltonians import hamiltonian
H = hamiltonian("LiH")  # 334 Pauli terms, 30 qubits
print(f"{H.n_terms} Pauli terms, {H.n_qubits} qubits")
qiskit_op = H.to_qiskit()  # SparsePauliOp for VQE
```

The `geovac-hamiltonians` package (v0.1.0, 92 KB) is a standalone wheel with no heavy dependencies. It bundles pre-computed Hamiltonians for H2, He, LiH, BeH2, and H2O. Export to OpenFermion, Qiskit, and PennyLane is available in the full `geovac` package.

## Collaboration Opportunities

- **Circuit compilation:** Block-encoding circuits (LCU, QROM) for the composed block-diagonal Hamiltonian structure. The structural sparsity should enable efficient encodings, but this has not been built.
- **Resource estimation:** End-to-end fault-tolerant resource estimates (T-count, logical qubit count) comparing GeoVac Hamiltonians to DF/THC Gaussian Hamiltonians at matched accuracy.
- **Hardware demonstrations:** VQE or ADAPT-VQE runs on current hardware using GeoVac's compact Hamiltonians. H2 at 10 qubits and LiH at 30 qubits are natural targets.
- **Measurement optimization:** Benchmarking classical shadows, derandomized shadows, or other measurement protocols on GeoVac's structurally sparse Hamiltonians.
- **Extending the molecule set:** The composed architecture is systematic (core pairs + bond pairs + lone pairs). New molecules require computing the Level 3/4 building blocks and coupling integrals.
- **Improving classical accuracy:** The 5-26% equilibrium bond length errors are dominated by basis truncation (l_max=2) and the adiabatic approximation. Higher angular momentum or non-adiabatic solvers would improve the Hamiltonian quality.

## Papers and Citation

The GeoVac paper series is available at [github.com/jlouthan/geometric-vacuum](https://github.com/jlouthan/geometric-vacuum) with DOI-stamped releases on Zenodo.

Key papers:
- **Paper 14** (qubit encoding): Pauli scaling analysis, composed sparsity, Gaussian baselines
- **Paper 7** (foundation): S3 conformal equivalence, 18 symbolic proofs
- **Paper 17** (composed geometries): LiH, BeH2, H2O molecular solvers
- **Paper 15** (Level 4): H2 molecule-frame hyperspherical solver
- **Paper 18** (exchange constants): Spectral-geometric classification of transcendental content

If you use GeoVac Hamiltonians in your work, please cite Paper 14 and the Zenodo DOI for the version you used.
