"""
PK Classical Partitioning — Quantum-Classical Splitting of Phillips-Kleinman
============================================================================

The Phillips-Kleinman pseudopotential (PK) is a one-body operator that prevents
valence electron collapse into core orbitals. Under Jordan-Wigner encoding,
one-body terms map to Pauli strings that encode 1-RDM elements
γ_{pq} = ⟨ψ|a†_p a_q|ψ⟩. The PK energy contribution E_PK = Tr(h1_PK · γ)
can be computed classically from measurement data, with zero additional
quantum circuits.

This partitioning is algebraically exact:
    E_total = E_quantum(pk=False) + E_PK_classical

It is the same class of technique as frozen-core partitioning, mean-field
embedding, and the symmetry-shift method (Loaiza et al. JCTC 2024).

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np


def pk_classical_energy(
    h1_pk: np.ndarray,
    one_rdm: np.ndarray,
) -> float:
    """
    Compute the PK energy contribution from the 1-RDM.

    E_PK = Tr(h1_PK · γ)

    where h1_PK is the PK one-electron integral matrix (spatial orbitals)
    and γ is the 1-RDM in the same basis.

    Parameters
    ----------
    h1_pk : ndarray of shape (M, M)
        PK one-electron integral matrix in spatial orbital basis.
    one_rdm : ndarray of shape (M, M)
        One-particle reduced density matrix γ_{pq} = ⟨ψ|a†_p a_q|ψ⟩
        in the same spatial orbital basis.

    Returns
    -------
    float
        E_PK in Hartree.
    """
    return float(np.trace(h1_pk @ one_rdm))


def reconstruct_1rdm_from_statevector(
    statevector: np.ndarray,
    n_spatial: int,
) -> np.ndarray:
    """
    Reconstruct the spatial 1-RDM from a statevector in the JW basis.

    The 1-RDM is γ_{pq} = Σ_σ ⟨ψ|a†_{p,σ} a_{q,σ}|ψ⟩ where σ sums
    over spin-up and spin-down.

    Under Jordan-Wigner encoding with 2M qubits (M spatial orbitals),
    spin-orbital index j = 2*p + σ where σ=0 (up) or 1 (down).

    Parameters
    ----------
    statevector : ndarray of shape (2^Q,)
        State vector in computational basis (JW encoding).
    n_spatial : int
        Number of spatial orbitals M. Q = 2*M qubits.

    Returns
    -------
    ndarray of shape (M, M)
        Spatial 1-RDM: γ_{pq} = Σ_σ ⟨ψ|a†_{p,σ} a_{q,σ}|ψ⟩.
    """
    Q = 2 * n_spatial
    dim = 2 ** Q
    assert len(statevector) == dim, (
        f"Statevector length {len(statevector)} != 2^{Q} = {dim}"
    )

    rdm = np.zeros((n_spatial, n_spatial))

    for p in range(n_spatial):
        for q in range(n_spatial):
            for sigma in range(2):
                # Spin-orbital indices under JW
                j_p = 2 * p + sigma
                j_q = 2 * q + sigma
                rdm[p, q] += _compute_1rdm_element_jw(
                    statevector, j_p, j_q, Q,
                )

    return rdm


def _compute_1rdm_element_jw(
    psi: np.ndarray,
    j: int,
    k: int,
    n_qubits: int,
) -> float:
    """
    Compute ⟨ψ|a†_j a_k|ψ⟩ in the Jordan-Wigner representation.

    Under JW, a†_j a_k acts on computational basis states |b⟩ as:
    - If j == k: projects onto occupation n_j, gives ⟨n_j⟩
    - If j != k: annihilates qubit k, creates qubit j, with JW string phase

    Parameters
    ----------
    psi : ndarray
        Statevector in computational basis.
    j, k : int
        Spin-orbital indices.
    n_qubits : int
        Total number of qubits.

    Returns
    -------
    float
        Real part of ⟨ψ|a†_j a_k|ψ⟩.
    """
    dim = len(psi)

    if j == k:
        # Number operator: ⟨ψ|n_j|ψ⟩ = Σ_b |ψ_b|² × bit_j(b)
        result = 0.0
        for b in range(dim):
            if (b >> j) & 1:
                result += abs(psi[b]) ** 2
        return result

    # Off-diagonal: a†_j a_k
    # a_k annihilates qubit k (requires bit k = 1, clears it)
    # a†_j creates qubit j (requires bit j = 0, sets it)
    # JW string: (-1)^{sum of occupations between min(j,k)+1 and max(j,k)-1}
    lo, hi = min(j, k), max(j, k)

    result = 0.0 + 0.0j
    for b in range(dim):
        # a_k requires bit k occupied
        if not ((b >> k) & 1):
            continue
        # a†_j requires bit j unoccupied
        if (b >> j) & 1:
            continue

        # JW parity string: count bits set in range (lo, hi) exclusive
        parity = 0
        for bit in range(lo + 1, hi):
            if (b >> bit) & 1:
                parity += 1

        # Apply: annihilate k, create j
        b_new = (b ^ (1 << k)) | (1 << j)
        phase = (-1) ** parity

        result += np.conj(psi[b_new]) * phase * psi[b]

    return float(np.real(result))


def validate_pk_partitioning(
    builder_func: Any,
    builder_kwargs: Dict[str, Any],
    n_electrons: int,
    system_name: str = "system",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Validate PK partitioning for a composed system.

    Builds the Hamiltonian with and without PK in the quantum part,
    diagonalizes both, reconstructs the 1-RDM, and verifies:
        |E_full - (E_elec + E_PK_classical)| < threshold

    Also checks whether the electronic-only ground state is the physical
    valence state or a collapsed core state.

    Parameters
    ----------
    builder_func : callable
        One of build_composed_lih, build_composed_beh2, build_composed_h2o.
    builder_kwargs : dict
        Keyword arguments for the builder (R, l_max, etc.). Should NOT
        include include_pk or pk_in_hamiltonian.
    n_electrons : int
        Number of valence electrons in the quantum simulation.
    system_name : str
        Label for output messages.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with validation results.
    """
    from openfermion import get_sparse_operator, jordan_wigner

    if verbose:
        print(f"\n{'='*60}")
        print(f"PK PARTITIONING VALIDATION — {system_name}")
        print(f"{'='*60}")

    # Build with PK in Hamiltonian (reference)
    kw_full = {**builder_kwargs, 'include_pk': True, 'pk_in_hamiltonian': True,
               'verbose': False}
    result_full = builder_func(**kw_full)
    M = result_full['M']
    Q = result_full['Q']

    # Build without PK in Hamiltonian (partitioned)
    kw_elec = {**builder_kwargs, 'include_pk': True, 'pk_in_hamiltonian': False,
               'verbose': False}
    result_elec = builder_func(**kw_elec)

    h1_pk = result_elec['h1_pk']

    if verbose:
        print(f"  M = {M} spatial orbitals, Q = {Q} qubits")
        print(f"  N_pauli (full)  = {result_full['N_pauli']}")
        print(f"  N_pauli (elec)  = {result_elec['N_pauli']}")
        pk_norm = sum(abs(c) for c in result_full['qubit_op'].terms.values())
        elec_norm = sum(abs(c) for c in result_elec['qubit_op'].terms.values())
        print(f"  1-norm (full)   = {pk_norm:.4f} Ha")
        print(f"  1-norm (elec)   = {elec_norm:.4f} Ha")
        print(f"  1-norm ratio    = {pk_norm / max(elec_norm, 1e-15):.1f}x")

    # Exact diagonalization (only feasible for small Q)
    if Q > 20:
        if verbose:
            print(f"  Q={Q} too large for exact diag. Skipping eigenvalue check.")
        return {
            'system': system_name,
            'M': M, 'Q': Q,
            'N_pauli_full': result_full['N_pauli'],
            'N_pauli_elec': result_elec['N_pauli'],
            'one_norm_full': sum(abs(c) for c in result_full['qubit_op'].terms.values()),
            'one_norm_elec': sum(abs(c) for c in result_elec['qubit_op'].terms.values()),
            'h1_pk': h1_pk,
            'exact_diag': False,
            'note': f'Q={Q} exceeds exact diag limit',
        }

    # Sparse matrix diagonalization
    from scipy.sparse.linalg import eigsh

    if verbose:
        print(f"  Exact diagonalization (2^{Q} = {2**Q} dim)...")

    H_full_sparse = get_sparse_operator(result_full['qubit_op'])
    H_elec_sparse = get_sparse_operator(result_elec['qubit_op'])

    # Build PK-only sparse operator for direct validation
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    eri_zero = np.zeros((M, M, M, M))
    pk_fermion = build_fermion_op_from_integrals(h1_pk, eri_zero, 0.0)
    pk_qubit = jordan_wigner(pk_fermion)
    H_pk_sparse = get_sparse_operator(pk_qubit, n_qubits=Q)

    # Verify operator decomposition: H_full = H_elec + H_pk
    op_diff = abs(H_full_sparse - H_elec_sparse - H_pk_sparse).max()
    if verbose:
        print(f"  ||H_full - H_elec - H_pk|| = {op_diff:.2e}")

    # Get ground states
    E_full, psi_full = eigsh(H_full_sparse, k=1, which='SA')
    E_full = float(E_full[0])
    psi_full = psi_full[:, 0]

    E_elec, psi_elec = eigsh(H_elec_sparse, k=1, which='SA')
    E_elec = float(E_elec[0])
    psi_elec = psi_elec[:, 0]

    if verbose:
        print(f"  E_full (with PK) = {E_full:.10f} Ha")
        print(f"  E_elec (no PK)   = {E_elec:.10f} Ha")

    # Same-state validation using DIRECT sparse operator evaluation:
    # E_full = <psi_full|H_full|psi_full>
    #        = <psi_full|H_elec|psi_full> + <psi_full|H_pk|psi_full>
    # This is algebraically exact regardless of particle number sector.
    E_elec_of_full_gs = float(np.real(
        np.dot(psi_full.conj(), H_elec_sparse @ psi_full)
    ))
    E_pk_of_full_gs = float(np.real(
        np.dot(psi_full.conj(), H_pk_sparse @ psi_full)
    ))
    residual_same_state = abs(E_full - (E_elec_of_full_gs + E_pk_of_full_gs))

    # Check if ground states match (overlap)
    overlap = float(abs(np.dot(psi_full.conj(), psi_elec)))

    # PK energy of the electronic ground state (for cross-state comparison)
    E_pk_of_elec_gs = float(np.real(
        np.dot(psi_elec.conj(), H_pk_sparse @ psi_elec)
    ))

    if verbose:
        print(f"\n  --- Same-state validation (direct operator) ---")
        print(f"  <psi_full|H_elec|psi_full> = {E_elec_of_full_gs:.10f} Ha")
        print(f"  <psi_full|H_pk|psi_full>   = {E_pk_of_full_gs:.10f} Ha")
        print(f"  Sum                        = {E_elec_of_full_gs + E_pk_of_full_gs:.10f} Ha")
        print(f"  Residual                   = {residual_same_state:.2e} Ha")

        print(f"\n  --- Ground state comparison ---")
        print(f"  |<psi_full|psi_elec>|      = {overlap:.6f}")
        print(f"  <psi_elec|H_pk|psi_elec>   = {E_pk_of_elec_gs:.10f} Ha")
        print(f"  E_elec + E_pk(psi_elec)    = {E_elec + E_pk_of_elec_gs:.10f} Ha")
        print(f"  vs E_full                  = {E_full:.10f} Ha")
        diff = abs(E_full - (E_elec + E_pk_of_elec_gs))
        print(f"  Difference                 = {diff:.6f} Ha")

        if overlap > 0.99:
            print(f"\n  [OK] Ground states match (overlap {overlap:.4f})")
            print(f"    Option A works: PK can be fully moved to classical post-processing")
        else:
            print(f"\n  [!!] Ground states DIFFER (overlap {overlap:.4f})")
            print(f"    Removing PK changes the ground state -- partial PK may be needed")
            print(f"    Note: eigsh searches all particle-number sectors.")
            print(f"    In VQE, the ansatz constrains particle number, so Option A")
            print(f"    may still work if the N-electron ground state is stable.")

    results = {
        'system': system_name,
        'M': M, 'Q': Q,
        'N_pauli_full': result_full['N_pauli'],
        'N_pauli_elec': result_elec['N_pauli'],
        'one_norm_full': sum(abs(c) for c in result_full['qubit_op'].terms.values()),
        'one_norm_elec': sum(abs(c) for c in result_elec['qubit_op'].terms.values()),
        'E_full': E_full,
        'E_elec': E_elec,
        'E_pk_of_full_gs': E_pk_of_full_gs,
        'E_pk_of_elec_gs': E_pk_of_elec_gs,
        'E_elec_of_full_gs': E_elec_of_full_gs,
        'residual_same_state': residual_same_state,
        'operator_decomposition_error': float(op_diff),
        'ground_state_overlap': overlap,
        'ground_states_match': overlap > 0.99,
        'h1_pk': h1_pk,
        'exact_diag': True,
    }

    return results
