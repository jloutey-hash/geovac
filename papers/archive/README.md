# The Geometric Vacuum: Quantum Mechanics from Discrete Information Geometry

**Author:** J. Loutey, Independent Researcher, Kent, Washington

## Overview

The Geometric Vacuum (GeoVac) program constructs quantum mechanics from discrete graph topology. The central insight, proven via 18 symbolic proofs (Paper 7), is that the discrete graph Laplacian on a lattice homologous to the unit three-sphere S^3 is mathematically equivalent to the Schrodinger equation: the 1/r Coulomb potential, the hydrogen spectrum, and the continuous wavefunction all emerge from a single dimensionless, scale-invariant combinatorial topology via Fock's 1935 stereographic projection.

This equivalence is exploited computationally by identifying the *natural geometry* of each quantum system -- the coordinate system where the Schrodinger equation separates -- and discretizing it as a sparse graph. The resulting lattice Hamiltonians are O(N) sparse matrices that replace expensive continuous integration with algebraic operations on quantum number labels.

This upload contains four papers forming a coherent cluster: from the foundational proof of the dimensionless vacuum principle (Paper 7), through the molecular prolate spheroidal lattice (Paper 11) and algebraic two-electron integrals (Paper 12), to the multi-particle hyperspherical lattice (Paper 13). Together, they establish the natural geometry hierarchy: S^3 for atoms, prolate spheroid for diatomics, and hyperspherical coordinates for multi-electron systems.

## Papers

| Paper | Title | Key Result |
|:------|:------|:-----------|
| **Paper 7** | The Dimensionless Vacuum: Recovering the Schrodinger Equation from Scale-Invariant Graph Topology | 18/18 symbolic proofs: graph Laplacian = S^3 Laplace-Beltrami = Schrodinger equation |
| **Paper 11** | The Molecular Fock Projection: Diatomic Molecules as Prolate Spheroidal Lattices | H2+ to 0.70% error with zero free parameters |
| **Paper 12** | Algebraic Two-Electron Integrals on the Prolate Spheroidal Lattice | Neumann expansion eliminates V_ee integration error; H2 at 92.4% D_e |
| **Paper 13** | The Hyperspherical Lattice: Two-Electron Atoms as Coupled Channel Graphs | He ground state to 0.05% error; first non-trivial fiber bundle in GeoVac; ab initio H2 rovibrational spectrum |

## Software

All results are produced by the GeoVac open-source Python package:

- **Repository:** [github.com/jloutey-hash/geovac](https://github.com/jloutey-hash/geovac)
- **Version:** v1.1.0
- **Language:** Python 3.10+ (NumPy, SciPy)
- **License:** MIT

## Citation

See `CITATION.cff` for machine-readable citation metadata.

## License

These papers are distributed under the Creative Commons Attribution 4.0 International License (CC-BY-4.0).
