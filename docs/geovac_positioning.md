# GeoVac: Structurally Sparse Qubit Hamiltonians for Molecular Quantum Simulation

## What GeoVac Does

GeoVac is a discretization framework that exploits the natural geometry of quantum systems (hyperspherical coordinates, prolate spheroidal coordinates, composed fiber bundles) to produce qubit Hamiltonians with far fewer Pauli terms than standard Gaussian-basis Jordan-Wigner encodings. The sparsity is structural -- it comes from angular momentum selection rules baked into the basis, not from post-hoc truncation or compression. The result is a direct reduction in the measurement budget for variational quantum eigensolvers (VQE) and other near-term quantum algorithms.

## The Headline Numbers

**Figures withdrawn 2026-10-02 (/qa trunk DELTA #26).** The table that stood
here (LiH 334, BeH₂ 556, H₂O 778 Pauli terms; "N_Pauli = 11.11 × Q"; He
counts) was measured under the retired pair-diagonal ERI rule (a wrong-sign
bug, corrected 2026-08-29) and before the He one-body hop was removed
(2026-09-25). Current values live in Paper 14 (authoritative) and the numeric
registry (`composed_coeff`, `lih_composed_pauli`, `lih_sto3g_raw`,
`lih_ccpvdz_pauli`, `lih_composed_lambda`). Under the exact rule the Pauli
count is exactly linear in Q across molecules, composed LiH is near parity
with STO-3G on raw counts, and the decisive reduction is against cc-pVDZ.

## Why It's Sparse

GeoVac's composed architecture assigns each electron group (core pair, bond pair, lone pair) its own natural coordinate system, then couples them via screened Coulomb integrals and Phillips-Kleinman pseudopotentials. Within each block, the angular basis is built from spherical harmonics on the natural geometry, so Gaunt integral selection rules enforce exact zeros in the electron repulsion integrals. The result is block-diagonal ERIs with no cross-block two-electron terms. This structural sparsity is intrinsic to the basis and survives all downstream transformations (Jordan-Wigner encoding, qubit tapering, measurement grouping, tensor factorization).

## Where It Fits

### Near-term / VQE

For variational algorithms on near-term hardware, the dominant cost is the measurement budget: the number of distinct Pauli terms times the number of shots per term. GeoVac's reduction in Pauli terms (values: Paper 14) translates directly to a proportional reduction in total measurements. This advantage is complementary to measurement optimization techniques (qubit-wise commuting grouping, classical shadows, derandomization) -- GeoVac reduces the number of terms those techniques have to work with.

Note on comparisons: density fitting (DF) and tensor hypercontraction (THC) are powerful Gaussian-basis compression techniques, but they target block-encoding representations for fault-tolerant algorithms. They do not reduce the number of Pauli terms in a Jordan-Wigner Hamiltonian. GeoVac's Pauli term advantage is against raw Jordan-Wigner encoding, which is the representation used by VQE and other near-term variational methods.

1-norm data (relevant for shot budgets and Trotter bounds): see Paper 14 and the numeric registry (the per-system 1-norms once listed here predate the 2026-08-29 exact rule and the 2026-09-25 He correction).

### Fault-tolerant / QPE

GeoVac's 1-norm scaling and commutator-based Trotter bounds (exponents in Paper 14; the figures once quoted here predate the 2026-08-29 exact rule) suggest competitive resource estimates for fault-tolerant quantum phase estimation. However, block-encoding circuits for the GeoVac Hamiltonian structure have not yet been built. This is a concrete collaboration opportunity: the composed block-diagonal structure should map naturally to LCU or QROM encodings, but the circuit construction and resource estimation remain open.

## Honest Limitations

GeoVac's classical accuracy is substantially below production quantum chemistry. The best classical results are H2 at 99.97% dissociation energy (explicit-r12 CI in prolate coordinates, Paper 12), LiH equilibrium bond length at 5.3% error, BeH2 at 11.7%, and H2O at 19.4% -- compared to sub-0.1% for Gaussian methods with correlating basis sets. For atoms, the 2D variational solver reaches 0.022% for He (variational), and a Schwartz cusp extrapolation 0.004% (non-variational, marginally below exact; zero free parameters). The framework's library covers 35 composed molecules across 3 periodic table rows. The value proposition is structural sparsity for quantum simulation, not classical accuracy.

Note on H2O 1-norm: The Phillips-Kleinman pseudopotential (PK) is a one-body operator whose energy contribution can be computed classically from the 1-RDM measured during VQE, with zero additional quantum circuits. This quantum-classical partitioning (v2.0.29) removes the PK contribution from the quantum 1-norm with no approximation. Current partitioned 1-norms: Paper 14/20 and the numeric registry (the per-system figures once quoted here predate the 2026-08-29 exact rule).

## Try It

```bash
pip install "geovac-hamiltonians @ git+https://github.com/jloutey-hash/geovac.git#subdirectory=geovac-hamiltonians"
```

```python
from geovac_hamiltonians import hamiltonian
H = hamiltonian("LiH")  # 30 qubits; Pauli count: Paper 14
print(f"{H.n_terms} Pauli terms, {H.n_qubits} qubits")
qiskit_op = H.to_qiskit()  # SparsePauliOp for VQE
```

The `geovac-hamiltonians` package is a standalone wheel with no heavy dependencies. Export to OpenFermion, Qiskit, and PennyLane is available. The full `geovac` package supports all 30 molecules with `hamiltonian('LiH', R=3.015)` API.

## Collaboration Opportunities

- **Circuit compilation:** Block-encoding circuits (LCU, QROM) for the composed block-diagonal Hamiltonian structure. The structural sparsity should enable efficient encodings, but this has not been built.
- **Resource estimation:** End-to-end fault-tolerant resource estimates (T-count, logical qubit count) comparing GeoVac Hamiltonians to DF/THC Gaussian Hamiltonians at matched accuracy.
- **Hardware demonstrations:** VQE or ADAPT-VQE runs on current hardware using GeoVac's compact Hamiltonians. H2 at 10 qubits and LiH at 30 qubits are natural targets.
- **Measurement optimization:** Benchmarking classical shadows, derandomized shadows, or other measurement protocols on GeoVac's structurally sparse Hamiltonians.
- **Extending the molecule set:** The composed architecture is systematic (core pairs + bond pairs + lone pairs) and now covers 30 molecules. New molecules require only a `MolecularSpec` definition and atomic classifier entry.
- **Improving classical accuracy:** The 5.3-19.4% equilibrium bond length errors are dominated by basis truncation (l_max=2) and the pseudopotential approximation. The balanced coupled architecture (Paper 19) achieves 0.20% energy error at n_max=3 for LiH.

## Papers and Citation

The GeoVac paper series is available at [github.com/jloutey-hash/geovac](https://github.com/jloutey-hash/geovac) with DOI-stamped releases on Zenodo.

Key papers:
- **Paper 14** (qubit encoding): Pauli scaling analysis, composed sparsity, Gaussian baselines
- **Paper 7** (foundation): S3 conformal equivalence, 18 symbolic proofs
- **Paper 17** (composed geometries): LiH, BeH2, H2O molecular solvers
- **Paper 15** (Level 4): H2 molecule-frame hyperspherical solver
- **Paper 18** (exchange constants): Spectral-geometric classification of transcendental content

If you use GeoVac Hamiltonians in your work, please cite Paper 14 and the Zenodo DOI for the version you used.
