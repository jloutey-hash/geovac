# Paper 14 Figures

Figures for "Structurally Sparse Qubit Hamiltonians from Spectral Graph Theory"

## Figure 1: pauli_scaling.png
Log-log plot of Pauli term count vs number of qubits Q.
- Red squares: Gaussian bases (STO-3G, 6-31G, cc-pVDZ) with O(Q^4.60) fit
- Blue circles: GeoVac lattice (nmax=2..5) with O(Q^3.15) fit
- Dashed reference lines at Q^2, Q^4
- Generate: `python benchmarks/pauli_term_scaling.py` (reuses existing plot logic)
- Source data: `benchmarks/qubit_encoding/results.md`

## Figure 2: eri_density.png
ERI density (fraction of nonzero integrals) vs number of spatial orbitals M.
- GeoVac: 10.4% (M=5), 3.9% (M=14), 2.1% (M=30), 1.3% (M=55) -- decays as ~1/M^2
- Gaussian: 50-100% -- flat or increasing
- Overlay: dashed 1/M^2 reference curve
- Generate: add plot function to benchmarks/pauli_term_scaling.py

## Figure 3: gaunt_sparsity.png (optional)
Visualization of the Gaunt coefficient sparsity pattern.
- Heat map of |c^k(l1,l2,l3)| for k=0..4, showing which (l1,l2,l3) triples are nonzero
- Illustrates the angular momentum selection rule origin of ERI sparsity
- Generate: standalone script in debug/
