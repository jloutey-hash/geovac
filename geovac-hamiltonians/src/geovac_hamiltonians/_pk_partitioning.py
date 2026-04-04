"""
PK Classical Partitioning — Quantum-Classical Splitting of Phillips-Kleinman
============================================================================

The Phillips-Kleinman pseudopotential (PK) is a one-body operator that prevents
valence electron collapse into core orbitals. Under Jordan-Wigner encoding,
one-body terms map to Pauli strings that encode 1-RDM elements
gamma_{pq} = <psi|a+_p a_q|psi>. The PK energy contribution E_PK = Tr(h1_PK . gamma)
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

    E_PK = Tr(h1_PK . gamma)

    where h1_PK is the PK one-electron integral matrix (spatial orbitals)
    and gamma is the 1-RDM in the same basis.

    Parameters
    ----------
    h1_pk : ndarray of shape (M, M)
        PK one-electron integral matrix in spatial orbital basis.
    one_rdm : ndarray of shape (M, M)
        One-particle reduced density matrix gamma_{pq} = <psi|a+_p a_q|psi>
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

    The 1-RDM is gamma_{pq} = Sum_sigma <psi|a+_{p,sigma} a_{q,sigma}|psi>
    where sigma sums over spin-up and spin-down.

    Under Jordan-Wigner encoding with 2M qubits (M spatial orbitals),
    spin-orbital index j = 2*p + sigma where sigma=0 (up) or 1 (down).

    Parameters
    ----------
    statevector : ndarray of shape (2^Q,)
        State vector in computational basis (JW encoding).
    n_spatial : int
        Number of spatial orbitals M. Q = 2*M qubits.

    Returns
    -------
    ndarray of shape (M, M)
        Spatial 1-RDM: gamma_{pq} = Sum_sigma <psi|a+_{p,sigma} a_{q,sigma}|psi>.
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
    Compute <psi|a+_j a_k|psi> in the Jordan-Wigner representation.

    Under JW, a+_j a_k acts on computational basis states |b> as:
    - If j == k: projects onto occupation n_j, gives <n_j>
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
        Real part of <psi|a+_j a_k|psi>.
    """
    dim = len(psi)

    if j == k:
        # Number operator: <psi|n_j|psi> = Sum_b |psi_b|^2 * bit_j(b)
        result = 0.0
        for b in range(dim):
            if (b >> j) & 1:
                result += abs(psi[b]) ** 2
        return result

    # Off-diagonal: a+_j a_k
    # a_k annihilates qubit k (requires bit k = 1, clears it)
    # a+_j creates qubit j (requires bit j = 0, sets it)
    # JW string: (-1)^{sum of occupations between min(j,k)+1 and max(j,k)-1}
    lo, hi = min(j, k), max(j, k)

    result = 0.0 + 0.0j
    for b in range(dim):
        # a_k requires bit k occupied
        if not ((b >> k) & 1):
            continue
        # a+_j requires bit j unoccupied
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
