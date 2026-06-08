"""
Total-spin (S^2) sector machinery for GeoVac qubit Hamiltonians.
================================================================

This module provides the total-spin operators ``S^2`` and ``S_z`` as
JW-encoded ``QubitOperator``s on the GeoVac composed-builder spin-orbital
basis, plus helpers for:

1.  Verifying ``[S^2, H] = 0`` on a constructed GeoVac qubit Hamiltonian.
2.  Building a VQE objective ``H + lambda * S^2`` that penalizes triplet
    (and higher-spin) leakage.
3.  Counting the singlet (``S = 0``) subspace dimension of the
    ``S_z = 0`` sector, exposing the FCI dimension reduction that
    a singlet-pure VQE ansatz captures.
4.  Diagonalising a (tapered) qubit Hamiltonian and reporting
    ``<S^2>`` per low-lying eigenstate.

Honest scope
------------
``S^2`` is an **integer-valued** Casimir of an SU(2) representation,
not a Z_2 stabilizer.  Standard "qubit tapering" handles only Z_2
stabilizers and therefore cannot remove ``log2(deg(S^2))`` qubits.
What this module delivers is **Hilbert-space-dimension reduction** and
**VQE convergence improvement**: the singlet sector is roughly a
factor of 2-3 smaller than the full ``S_z = 0`` sector for closed-shell
molecules, and a singlet-projected ansatz cannot leak amplitude into
the spurious triplets that share ``S_z = 0``.

Conventions
-----------
GeoVac composed builders use the standard OpenFermion JW convention with
interleaved spin orbitals.  For ``M`` spatial orbitals there are
``Q = 2*M`` qubits; spatial orbital ``p`` carries

    alpha (spin-up)   <-> qubit  2*p
    beta  (spin-down) <-> qubit  2*p + 1

This is exactly the convention used by ``openfermion.s_squared_operator``
and ``openfermion.sz_operator``, both of which we delegate to.

Per CLAUDE.md positioning
-------------------------
``[S^2, H] = 0`` is a property of the *physical* Hamiltonian on the
full spin-orbital basis (the alpha-beta interleaved JW encoding here).
GeoVac's standard composed builder, the per-block Hopf-Z_2 tapering of
:mod:`geovac.z2_tapering`, and the global tapering all preserve total
spin symmetry, so the gate :func:`verify_s2_commutes` is expected to pass
to machine precision on those Hamiltonians.  Aggressive transformations
that mix alpha/beta orbitals (for example a Bogoliubov rotation of the
JW frame) would break this; the gate detects such breakage.

References
----------
- OpenFermion ``hamiltonians.special_operators``: definitions of
  ``s_squared_operator``, ``sz_operator``, ``s_plus_operator``,
  ``s_minus_operator``.
- Paper 14 §sec:hopf_tapering (per-block Hopf-Z_2 tapering).
- CLAUDE.md §1.5 "Quantum simulation positioning" — composed Pauli
  scaling and ``S^2`` are complementary structural optimizations.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    from openfermion import (
        QubitOperator,
        s_squared_operator,
        sz_operator,
        jordan_wigner,
    )
    from openfermion.utils import commutator
    from openfermion.linalg import get_sparse_operator
    _HAS_OPENFERMION = True
except ImportError:  # pragma: no cover
    _HAS_OPENFERMION = False


__all__ = [
    "compute_s2_operator",
    "compute_sz_operator",
    "verify_s2_commutes",
    "s2_penalty_hamiltonian",
    "singlet_sector_dim",
    "sz_zero_sector_dim",
    "compute_fci_ground_state_s2",
]


# ---------------------------------------------------------------------------
# S^2 and S_z as QubitOperators
# ---------------------------------------------------------------------------

def compute_s2_operator(M: int) -> "QubitOperator":
    r"""Return the total-spin-squared operator on ``Q = 2*M`` JW qubits.

    Constructs ``S^2 = S^- S^+ + S_z (S_z + 1)`` as a
    ``FermionOperator`` (via :func:`openfermion.s_squared_operator`),
    Jordan-Wigner transforms it to a ``QubitOperator`` Hermitian Pauli
    sum on the interleaved alpha/beta JW basis.

    Parameters
    ----------
    M : int
        Number of spatial orbitals.

    Returns
    -------
    QubitOperator
        Hermitian operator on ``2*M`` qubits, with at most ``O(M^2)``
        Pauli strings.

    Notes
    -----
    For closed-shell singlet ground states ``<S^2> = 0`` exactly.
    The eigenvalues of ``S^2`` are ``S * (S + 1)`` for non-negative
    integer or half-integer ``S``.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required for compute_s2_operator")
    if not isinstance(M, int) or M <= 0:
        raise ValueError(f"M must be a positive integer; got {M!r}")
    return jordan_wigner(s_squared_operator(M))


def compute_sz_operator(M: int) -> "QubitOperator":
    r"""Return the total-spin-z operator on ``Q = 2*M`` JW qubits.

    ``S_z = (1/2) * sum_p (n_{p alpha} - n_{p beta})`` in the
    interleaved JW convention.  Used for sanity checks and as the
    second Casimir labelling SU(2) sectors together with ``S^2``.

    Parameters
    ----------
    M : int
        Number of spatial orbitals.

    Returns
    -------
    QubitOperator
        Hermitian operator on ``2*M`` qubits, a sum of diagonal
        ``Z_{2p}`` and ``Z_{2p+1}`` terms.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required for compute_sz_operator")
    if not isinstance(M, int) or M <= 0:
        raise ValueError(f"M must be a positive integer; got {M!r}")
    return jordan_wigner(sz_operator(M))


# ---------------------------------------------------------------------------
# Commutation gate
# ---------------------------------------------------------------------------

def verify_s2_commutes(
    qubit_op: "QubitOperator",
    M: int,
    atol: float = 1e-10,
) -> Tuple[bool, float]:
    """Verify ``[S^2, H] = 0`` for a GeoVac qubit Hamiltonian.

    Computes the commutator at the ``QubitOperator`` level (exact Pauli
    arithmetic, no matrix construction) and reports its operator-1-norm
    -- the sum of absolute values of all surviving Pauli coefficients.

    The commutator is exactly zero for a spin-symmetric Hamiltonian.
    Floating-point construction can leave a residue around
    ``1e-12 .. 1e-14``; the default ``atol = 1e-10`` is comfortably
    above that.

    Parameters
    ----------
    qubit_op : QubitOperator
        The Hamiltonian to test.  Must act on ``2*M`` qubits in the
        standard interleaved JW convention.
    M : int
        Number of spatial orbitals.  Must match the encoding of
        ``qubit_op``.
    atol : float
        Absolute tolerance on the commutator 1-norm.

    Returns
    -------
    (commutes, residual) : (bool, float)
        ``commutes = True`` iff ``residual < atol``.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required for verify_s2_commutes")

    s2 = compute_s2_operator(M)
    comm = commutator(s2, qubit_op)
    residual = float(sum(abs(c) for c in comm.terms.values()))
    return residual < atol, residual


# ---------------------------------------------------------------------------
# VQE penalty objective
# ---------------------------------------------------------------------------

def s2_penalty_hamiltonian(
    qubit_op: "QubitOperator",
    M: int,
    lam: float,
) -> "QubitOperator":
    r"""Return ``H_eff = H + lambda * S^2`` for spin-projected VQE.

    For ``lam > 0`` the singlet sector (``S^2 = 0``) is unshifted,
    triplet states are lifted by ``2 * lam``, quintets by ``6 * lam``,
    etc.  Used as a VQE objective when targeting a closed-shell
    singlet ground state: any spurious triplet contamination raises
    the variational energy by ``2 * lam`` and is filtered out.

    Parameters
    ----------
    qubit_op : QubitOperator
        The physical Hamiltonian (which is taken to commute with ``S^2``;
        callers should verify with :func:`verify_s2_commutes`).
    M : int
        Number of spatial orbitals.
    lam : float
        Penalty strength in Hartree.  Typical values:
        ``lam = 1.0 .. 10.0`` Ha are generous; ``lam`` should exceed
        any splitting one cares to resolve between the singlet ground
        state and the lowest triplet.

    Returns
    -------
    QubitOperator
        ``qubit_op + lam * compute_s2_operator(M)``.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required for s2_penalty_hamiltonian")
    if lam == 0.0:
        # Return a copy to avoid aliasing the input.
        return qubit_op + QubitOperator()
    return qubit_op + lam * compute_s2_operator(M)


# ---------------------------------------------------------------------------
# Sector dimensions (closed-form)
# ---------------------------------------------------------------------------

def _binom(n: int, k: int) -> int:
    """Integer binomial coefficient ``C(n, k)`` with the usual conventions."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    # Build with integers only (avoid math.comb on edge cases).
    from math import comb
    return comb(n, k)


def sz_zero_sector_dim(M: int, n_electrons: int) -> int:
    r"""Dimension of the ``S_z = 0`` sector of the FCI Hilbert space.

    With ``n_alpha = n_beta = N/2`` electrons, the sector dimension is
    ``C(M, N/2)^2``.

    Parameters
    ----------
    M : int
        Number of spatial orbitals.
    n_electrons : int
        Total electron count.  Must be even (otherwise ``S_z = 0`` does
        not exist exactly).

    Returns
    -------
    int
        Dimension of the ``S_z = 0`` subspace.
    """
    if n_electrons % 2 != 0:
        raise ValueError(
            f"S_z = 0 requires even n_electrons; got {n_electrons}"
        )
    n_half = n_electrons // 2
    return _binom(M, n_half) ** 2


def singlet_sector_dim(M: int, n_electrons: int) -> int:
    r"""Dimension of the singlet (``S = 0``) sector for ``n_electrons`` in
    ``M`` spatial orbitals.

    Uses the Weyl branching identity for ``S = 0`` on ``N`` electrons in
    ``M`` orbitals,

    .. math::

        \dim_{S = 0}(M, N) = \binom{M}{N/2}^2 - \binom{M}{N/2 - 1}\binom{M}{N/2 + 1},

    a special case of the well-known formula

    .. math::

        \dim_S(M, N) = \frac{2S + 1}{M + 1}\binom{M + 1}{N/2 - S}
                       \binom{M + 1}{N/2 + S + 1}

    (see e.g. Pauncz, *Spin Eigenfunctions*, 1979; or the Weyl
    dimension formula in any group theory textbook).  At ``S = 0`` and
    after the standard simplification this reduces to the
    binomial-difference form above.

    Parameters
    ----------
    M : int
        Number of spatial orbitals.
    n_electrons : int
        Total electron count.  Must be even.

    Returns
    -------
    int
        Dimension of the singlet sector.  Equals
        :func:`sz_zero_sector_dim` minus the dimension of all higher-
        spin (triplet, quintet, ...) ``S_z = 0`` projections.
    """
    if n_electrons % 2 != 0:
        raise ValueError(
            f"Closed-shell singlet requires even n_electrons; got {n_electrons}"
        )
    if n_electrons < 0 or n_electrons > 2 * M:
        raise ValueError(
            f"n_electrons={n_electrons} out of range for M={M} "
            f"(must be in [0, {2 * M}])"
        )
    n_half = n_electrons // 2
    dim_sz0 = _binom(M, n_half) ** 2
    cross = _binom(M, n_half - 1) * _binom(M, n_half + 1)
    return dim_sz0 - cross


# ---------------------------------------------------------------------------
# FCI diagonaliser with <S^2> per state
# ---------------------------------------------------------------------------

def compute_fci_ground_state_s2(
    qubit_op: "QubitOperator",
    M_effective: int,
    k: int = 3,
    n_qubits: Optional[int] = None,
) -> List[Dict[str, Any]]:
    r"""Diagonalise (a tapered) qubit Hamiltonian, return lowest ``k``
    eigenvalues with their ``<S^2>`` expectation values.

    Builds the sparse matrix representations of ``qubit_op`` and
    ``S^2_{M_effective}`` on the **same** ``2 * M_effective`` JW basis,
    diagonalises to find the lowest ``k`` eigenvalues, and computes
    ``<psi | S^2 | psi>`` per eigenstate.

    .. warning::
        ``M_effective`` is the number of spatial orbitals **in the
        encoding currently used by ``qubit_op``**.  For an untapered
        composed Hamiltonian this equals ``spec.M``.  After per-block
        Hopf tapering the qubit count changes but the spatial orbital
        index of ``S^2`` does not align trivially -- you should call
        this on the *untapered* Hamiltonian (or wrap the tapered one
        only when you have constructed the tapered ``S^2`` as well).
        The honest scope check :func:`verify_s2_commutes` should
        always be run first.

    Parameters
    ----------
    qubit_op : QubitOperator
        The Hamiltonian.
    M_effective : int
        Number of spatial orbitals matching ``qubit_op``'s encoding.
    k : int
        Number of low-lying eigenvalues to return.
    n_qubits : int, optional
        Override the inferred qubit count.  Defaults to
        ``2 * M_effective``.

    Returns
    -------
    list of dict
        Each entry has keys ``'energy'``, ``'s2_expectation'``,
        ``'s_estimate'``.  ``'s_estimate'`` is the inferred total
        spin ``S`` (from ``S^2 = S(S+1)``); not necessarily an
        integer/half-integer if the eigenstate is mixed-spin (which
        only happens when ``[S^2, H] != 0``).
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for compute_fci_ground_state_s2"
        )
    if n_qubits is None:
        n_qubits = 2 * M_effective

    H_sparse = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    S2_sparse = get_sparse_operator(
        compute_s2_operator(M_effective), n_qubits=n_qubits,
    )

    dim = H_sparse.shape[0]
    if k >= dim:
        # ARPACK requires k < dim; fall back to dense for tiny problems.
        H_dense = H_sparse.toarray()
        S2_dense = S2_sparse.toarray()
        eigvals, eigvecs = np.linalg.eigh(H_dense)
        results = []
        for i in range(min(k, dim)):
            psi = eigvecs[:, i]
            s2_exp = float(np.real(np.vdot(psi, S2_dense @ psi)))
            results.append(_pack_eigen_result(eigvals[i], s2_exp))
        return results

    from scipy.sparse.linalg import eigsh
    # ``which='SA'`` -> smallest algebraic eigenvalues (most negative first).
    eigvals, eigvecs = eigsh(H_sparse, k=k, which='SA')
    # eigsh does not guarantee sorted output; sort here.
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    results = []
    for i in range(k):
        psi = eigvecs[:, i]
        s2_exp = float(np.real(np.vdot(psi, S2_sparse @ psi)))
        results.append(_pack_eigen_result(eigvals[i], s2_exp))
    return results


def _pack_eigen_result(energy: float, s2_exp: float) -> Dict[str, Any]:
    """Build one ``compute_fci_ground_state_s2`` result entry."""
    # Invert S * (S + 1) = s2_exp.  Robust against tiny negative noise.
    s2_clip = max(0.0, s2_exp)
    s_est = 0.5 * (-1.0 + (1.0 + 4.0 * s2_clip) ** 0.5)
    return {
        'energy': float(energy),
        's2_expectation': float(s2_exp),
        's_estimate': float(s_est),
    }
