# GeoVac: 190x Fewer Pauli Terms for Molecular Quantum Simulation

I built a framework that produces qubit Hamiltonians with 190x-720x fewer Pauli terms than standard Gaussian-basis Jordan-Wigner encodings for molecular quantum simulation.

The concrete result: a lithium hydride Hamiltonian on 30 qubits with 334 Pauli terms, where a comparable Gaussian cc-pVDZ encoding produces roughly 63,500. For water on 70 qubits, it's 778 terms versus an estimated 560,000. The scaling exponent is Q^2.5 across all molecules tested, compared to Q^3.9-4.3 for Gaussian Jordan-Wigner.

The sparsity is structural, not approximate. GeoVac assigns each electron group in a molecule -- core pairs, bond pairs, lone pairs -- its own natural coordinate system (hyperspherical or prolate spheroidal), then couples them through screened Coulomb integrals. Angular momentum selection rules enforced by the basis produce exact zeros in the electron repulsion integrals, giving block-diagonal two-electron terms with no cross-block entries. This structure survives Jordan-Wigner encoding and is compatible with all downstream optimizations.

For near-term quantum computing, this matters because the measurement budget in VQE scales with the number of Pauli terms. Fewer terms means proportionally fewer measurements to estimate the energy. The 1-norm (which controls shot noise and Trotter error) also shows significant advantages: 3.8x at matched qubit count for helium, scaling as Q^1.69.

I want to be upfront about limitations. GeoVac's classical accuracy is well below production quantum chemistry -- 5% equilibrium bond length error for LiH, 12% for BeH2, 26% for water, compared to sub-0.1% for Gaussian methods. The framework currently covers five systems (H2, He, LiH, BeH2, H2O). And the Pauli advantage is against raw Jordan-Wigner -- density fitting and tensor hypercontraction use block-encoding representations that bypass the Pauli picture entirely, so the comparison is specifically relevant to VQE and other variational near-term methods. The pseudopotential barrier contribution (98.7% of the raw H2O 1-norm) is computed classically from the VQE 1-RDM with zero additional circuits, reducing the quantum 1-norm to 361 Ha -- comparable to BeH2's 355 Ha.

The Hamiltonians are available as a standalone Python package (`pip install geovac-hamiltonians`, 92 KB, no heavy dependencies) with export to OpenFermion, Qiskit, and PennyLane. The underlying theory is documented in a series of papers on GitHub with DOI-stamped Zenodo releases.

I'm looking for collaborators on three fronts: running VQE or ADAPT-VQE on current hardware with these compact Hamiltonians, building block-encoding circuits to benchmark against DF/THC for fault-tolerant algorithms, and extending the molecule set beyond the current five systems. If any of this is relevant to your work, I'd welcome the conversation.
