# geovac-hamiltonians

Structurally sparse qubit Hamiltonians for quantum chemistry, built from the
GeoVac composed natural geometry framework. The discrete graph Laplacian on
the three-sphere encodes angular momentum selection rules directly into the
basis, producing Hamiltonians with 51x--1,712x fewer Pauli terms than
equivalent Gaussian-basis encodings and O(Q^2.5) Pauli scaling (vs O(Q^4) for
Gaussian bases).

## Install

```bash
pip install geovac-hamiltonians
```

Optional exports:

```bash
pip install geovac-hamiltonians[qiskit]      # Qiskit SparsePauliOp export
pip install geovac-hamiltonians[pennylane]    # PennyLane Hamiltonian export
pip install geovac-hamiltonians[all]          # both
```

## Quick start

```python
from geovac_hamiltonians import hamiltonian

# Build a 30-qubit LiH Hamiltonian (334 Pauli terms)
H = hamiltonian('LiH')
print(H.n_qubits, H.n_terms, f"{H.one_norm:.2f} Ha")
# 30  334  37.33 Ha

# Export to OpenFermion (always available)
of_op = H.to_openfermion()

# Export to Qiskit (requires pip install qiskit)
qk_op = H.to_qiskit()

# Export to PennyLane (requires pip install pennylane)
pl_op = H.to_pennylane()
```

## Benchmark

| System | Qubits | Pauli Terms | Gaussian Pauli Terms | Advantage |
|--------|--------|-------------|---------------------|-----------|
| H2     | 10     | 112         | 15 (STO-3G, Q=4)   | --        |
| He     | 10     | 120         | 156 (cc-pVDZ, Q=10) | 1.3x     |
| LiH    | 30     | 334         | 63,519 (cc-pVDZ, Q=36) | **190x** |
| BeH2   | 50     | 556         | --                  | --        |
| H2O    | 70     | 778         | 107,382 (cc-pVDZ, Q=46) | **138x** |

Scaling exponents (Pauli terms vs qubits):
- GeoVac composed: **Q^2.50** (LiH/BeH2/H2O, exponent spread 0.02)
- Gaussian (Trenev et al. 2025): Q^3.9--4.3

The sparsity advantage grows with system size. At Q=84 (LiH n_max=3), the
estimated advantage exceeds 300x.

**Accuracy caveat:** GeoVac classical accuracy is R_eq 5.3% (LiH), 11.7%
(BeH2), 26% (H2O). Gaussian cc-pVDZ achieves < 0.1%. The structural sparsity
is the framework's primary advantage for quantum simulation, not classical
accuracy.

## Supported systems

| System | Function call | Qubits | Default geometry |
|--------|--------------|--------|-----------------|
| H2     | `hamiltonian('H2')` | 10 | R = 1.4 bohr |
| He     | `hamiltonian('He')` | 10 | atomic |
| LiH    | `hamiltonian('LiH')` | 30 | R = 3.015 bohr |
| BeH2   | `hamiltonian('BeH2')` | 50 | R = 2.502 bohr |
| H2O    | `hamiltonian('H2O')` | 70 | R = 1.809 bohr |

Parameters:
- `R` (float): internuclear distance in bohr
- `l_max` (int): angular momentum cutoff for composed systems (default 2)
- `max_n` (int): principal quantum number cutoff for atomic/H2 (default 2)

## How it works

GeoVac replaces Gaussian basis functions with the natural coordinate system
for each electron group in a molecule. The discrete graph Laplacian on S^3
(Fock's 1935 stereographic projection) produces integer eigenvalues with
built-in angular momentum selection rules. For molecules, composed natural
geometries (core + valence fiber bundles) factorize the many-electron problem
into independent blocks connected by one-body screening potentials, yielding
block-diagonal two-electron integrals and O(Q^2.5) Pauli scaling.

See the [GeoVac paper series](https://github.com/jlouthan/geometric-vacuum)
for full derivations, particularly:
- **Paper 7**: S^3 conformal equivalence proof (18 symbolic proofs)
- **Paper 14**: Qubit encoding and Pauli scaling analysis
- **Paper 17**: Composed natural geometries for molecules

## Citation

```bibtex
@software{geovac2026,
  title  = {GeoVac: Structurally Sparse Qubit Hamiltonians from Geometric Quantum Chemistry},
  author = {GeoVac Development Team},
  year   = {2026},
  url    = {https://github.com/jlouthan/geometric-vacuum},
}
```

## License

MIT
