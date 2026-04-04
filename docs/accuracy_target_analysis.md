# Accuracy Requirements for Quantum Simulation Hamiltonians

## What the quantum computing community actually measures

The fault-tolerant quantum phase estimation (FT-QPE) literature does not use equilibrium geometry error as an accuracy metric. The landmark resource estimation papers — Reiher et al. (2017) on FeMoco, Lee et al. (2021) on double-factorized Hamiltonians, Goings et al. (2022) on cytochrome P450 — all operate at a single fixed molecular geometry and ask: how many logical qubits and T-gates are needed to resolve the ground-state energy to within chemical accuracy (1.6 mHa, or 1 kcal/mol) at that geometry? The Hamiltonian's job is to encode the physics at a given nuclear configuration. Whether that Hamiltonian's potential energy surface reproduces the experimental bond length is a separate question about basis set quality, and it is not the question these papers address.

The metrics that determine quantum simulation cost are: (i) the number of Pauli terms, which sets the measurement overhead for VQE and the gate count per Trotter step; (ii) the Pauli 1-norm, which controls QPE circuit depth via Trotter step count or qubitization query complexity; and (iii) the number of qubits, which determines hardware requirements. A Hamiltonian that is structurally sparser — fewer Pauli terms, lower 1-norm — is cheaper to simulate regardless of whether its classical PES matches experiment to 0.1% or 5%.

This distinction matters for GeoVac. The R_eq errors (5.3% for LiH, 11.7% for BeH2, 26% for H2O) measure how well the composed-geometry PES reproduces the experimental equilibrium distance. They do not measure how accurately the Hamiltonian represents the electronic structure at a given geometry. An R_eq error of 5% means the PES minimum is shifted, not that the electronic energy at any particular bond length is 5% wrong. For quantum simulation benchmarking, the relevant comparison is: at the same geometry and the same qubit count, how do GeoVac and Gaussian Hamiltonians compare on simulation cost?


## Error decomposition for GeoVac composed systems

The table below decomposes the current accuracy limitations by system and source. "Basis-fixable" means the error decreases with increasing l_max; "structural" means it persists at any l_max within the current architecture.

| System | R_eq error | Dominant source | Category | Notes |
|:-------|:-----------|:----------------|:---------|:------|
| LiH | 5.3% | Basis truncation (l_max=2) | Basis-fixable | H2 standalone converges from 83% to 96% D_e over l_max=0-6; LiH drift is +0.15 bohr/l_max (linear, slow) |
| BeH2 | 11.7% | Basis truncation + adiabatic approx. | Mostly basis-fixable | Full 1-RDM exchange validated; remaining error from l_max=2 and single-channel treatment |
| H2O | 26% | Level 4 angular basis at 6:1 charge asymmetry | Structural + basis | Lone pair coupling disabled (unphysical at Z_eff >= 6); charge-center origin helps but does not resolve the asymmetry |

Three structural limitations apply across all composed systems. First, the Phillips-Kleinman pseudopotential produces wrong-sign s-p orbital splitting: PK raises the energy of orbitals that overlap the core, but in real atoms s-orbitals are more bound than p due to core penetration. This is a structural impossibility, not a parameter issue — the sign is wrong by construction. For first-row atoms, orbital filling order is trivial and manually assigned, so this does not affect production results, but it prevents extension to systems where filling order matters. Second, lone-pair Slater integrals become unphysical at Z_eff >= 6, producing coupling magnitudes that exceed total electronic energies. Bond-bond coupling (~0.5 Ha) is validated and physical. Third, cross-block electron repulsion integrals are set to zero by the composed approximation; selection-rule counting bounds their contribution at no more than 2x the reported Pauli counts.

The critical observation is that none of these limitations affect the Pauli term count or its scaling. The O(Q^2.5) scaling and the 51x-1,712x Pauli count advantage over Gaussian baselines are structural properties of the block-diagonal architecture and Gaunt selection rules. They hold independently of the energy accuracy.


## The accuracy bar for publishable quantum resource advantage

Three possible accuracy standards exist, in decreasing order of stringency.

The first is chemical accuracy at a fixed geometry: the Hamiltonian must reproduce the FCI energy of a benchmark Gaussian calculation to within 1.6 mHa at the same geometry. This is the gold standard for quantum chemistry, but it is not what resource estimation papers require of the Hamiltonian encoding. It is a property of the basis set, not of the quantum simulation framework. GeoVac's hydrogenic basis is not designed to compete with cc-pVQZ Gaussian bases on absolute energy — it is designed to produce sparser Hamiltonians.

The second is correct qualitative physics: the Hamiltonian must produce bound states, correct symmetry, and physically reasonable energy scales. GeoVac meets this bar for all six first-row molecules tested (LiH, BeH2, H2O, HF, NH3, CH4). The composed Hamiltonians produce bound states with the correct orbital structure, and the PES shapes are qualitatively correct even where R_eq is shifted.

The third is structural resource advantage with honest caveats: the Hamiltonian encoding demonstrably reduces quantum simulation cost (Pauli terms, 1-norm, measurement groups) relative to established baselines, with the accuracy limitations stated explicitly. This is the standard that the quantum computing literature actually applies. Papers routinely compare Hamiltonian encodings — double factorization vs. tensor hypercontraction vs. plane waves vs. Gaussians — on resource cost, acknowledging that different bases have different accuracy profiles. The question is always: for a given computational budget, which encoding gets you closer to the answer?

GeoVac's value proposition fits squarely in the third category. The composed architecture produces 190x fewer Pauli terms for LiH at Q~30, 746x fewer for H2O at Q=70, with O(Q^2.5) scaling versus O(Q^3.9-4.3). The PK partitioning result (v2.0.29) further reduces the effective 1-norm by 78x for H2O by isolating the PK barrier as a classically computable one-body operator. These are large, robust advantages on the metrics that determine quantum simulation cost. The accuracy limitations are real but orthogonal: they affect how well the classical PES reproduces experiment, not how efficiently the Hamiltonian can be simulated.


## Recommendation for Paper 14

The current accuracy caveat in Paper 14 (Section IV.E) is thorough and honest, listing seven specific limitations. However, it frames accuracy and sparsity as competing concerns, which implicitly concedes a weakness that is not relevant to the paper's central claim. The framing should be restructured around three points.

First, state clearly that the comparison is a resource estimation comparison, not a basis-set quality comparison. This is already present in limitation (6) but should be elevated to the opening of the limitations subsection. The sentence "This comparison should not be interpreted as a basis-set quality comparison" should appear before the enumerated list, not buried as item (6).

Second, separate the accuracy discussion into "affects quantum simulation cost" and "does not affect quantum simulation cost." Cross-block ERIs being zero affects both accuracy and Pauli count (by at most 2x); this belongs in a simulation-cost discussion. PK Z^2-scaling, constant Z_eff, and inter-bond coupling affect energy accuracy only; these are basis-quality caveats and should be grouped as such.

Third, note that R_eq error — the most frequently cited accuracy figure — measures PES shape fidelity, which is relevant for geometry optimization but not for single-point energy calculations. For quantum simulation benchmarking at a fixed geometry, the relevant question is whether the Hamiltonian encodes the correct physics at that configuration. The composed Hamiltonians produce bound states with correct symmetry and physically reasonable energies for all tested systems.

The honest summary is: GeoVac's composed Hamiltonians are 2-3 orders of magnitude sparser than Gaussian alternatives on the metrics that determine quantum simulation cost; the basis produces larger energy errors than production Gaussian bases; and the two facts are largely independent because sparsity comes from selection-rule structure, not from accuracy of the radial wavefunctions.
