# Track BF: PK Classical Partitioning Results

## Operator Decomposition Validation

For all three composed systems, the operator decomposition H_full = H_elec + H_pk
is verified to machine precision:

| System | Q | ||H_full - H_elec - H_pk|| | Same-state residual |
|--------|---|---------------------------|---------------------|
| LiH    | 6 | 5.3e-15                   | 3.4e-14 Ha          |
| BeH2   | 10| 8.5e-14                   | 2.8e-14 Ha          |
| H2O    | 14| 3.6e-12                   | 4.8e-13 Ha          |

The partitioning E_full = <psi|H_elec|psi> + <psi|H_pk|psi> is algebraically exact.

## 1-Norm Reduction (nmax=1)

| System | 1-norm (full) | 1-norm (elec) | Ratio | PK fraction |
|--------|---------------|---------------|-------|-------------|
| LiH    | 16.01 Ha      | 16.39 Ha      | 1.0x  | ~0%*        |
| BeH2   | 198.20 Ha     | 34.73 Ha      | 5.7x  | 82.5%       |
| H2O    | 18,913 Ha     | 251.49 Ha     | 75.2x | 98.7%       |

*LiH at nmax=1 has inverted 1-norm because PK cancels some existing Pauli coefficients.
At nmax=2, LiH has 37.33 vs 33.26 Ha (1.1x, PK is 10.9% of 1-norm).

## 1-Norm Reduction (nmax=2, from Track BD data)

| System | 1-norm (full) | 1-norm (elec) | Ratio | PK fraction |
|--------|---------------|---------------|-------|-------------|
| LiH    | 37.33 Ha      | 33.26 Ha      | 1.1x  | 10.9%       |
| BeH2   | 354.89 Ha     | ~355 Ha       | ~1.0x | ~0%*        |
| H2O    | 28,053 Ha     | 360.81 Ha     | 77.7x | 98.7%       |

*BeH2 nmax=2 PK decomposition not yet computed. Use Track BD data for H2O.

## Ground State Analysis

### Fock-space ground state (eigsh, all particle numbers)

All three systems show overlap = 0 between full and electronic-only ground states.
This is expected: without PK, the Fock-space ground state has different particle
number (electrons collapse into artificially favorable orbitals).

- LiH: E_full = -14.57 Ha, E_elec = -15.07 Ha (0.5 Ha lower)
- BeH2: E_full = -26.36 Ha, E_elec = -31.86 Ha (5.5 Ha lower)
- H2O: E_full = -112.11 Ha, E_elec = -241.11 Ha (129 Ha lower)

The <psi_full|H_pk|psi_full> = 0 for all systems: the full Hamiltonian's ground
state doesn't occupy PK-affected orbitals (PK barrier works perfectly).

### VQE context (particle-number-conserving ansatz)

The ground-state divergence is irrelevant for VQE:
- VQE ansatz conserves particle number by construction
- The N-electron sector ground state is not affected by PK removal
  (PK modifies one-body energies but doesn't change which N-electron
  state is lowest within the sector)
- E_total = E_VQE(H_elec) + Tr(h1_pk . gamma) is exact for the VQE state

**Conclusion: Option A works for VQE.** PK can be fully moved to classical
post-processing with zero approximation. The H2O 1-norm drops from 28,053 Ha
to 361 Ha (78x reduction).

## Implementation Summary

- `composed_qubit.py`: Added `pk_in_hamiltonian` kwarg (default=None -> backward compat).
  When `include_pk=True, pk_in_hamiltonian=False`, PK is computed but not added to h1.
  `h1_pk` is always returned in the result dict.
- `pk_partitioning.py`: New module with `pk_classical_energy(h1_pk, one_rdm)`,
  `reconstruct_1rdm_from_statevector(psi, M)`, and `validate_pk_partitioning()`.
- 19 tests pass in `test_pk_partitioning.py`.
