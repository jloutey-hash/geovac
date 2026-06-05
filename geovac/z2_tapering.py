"""
Hopf-U(1) Z2 tapering for GeoVac composed qubit Hamiltonians.
=============================================================

Production module for the Hopf-U(1) ``m_l -> -m_l`` reflection tapering
identified in Paper 29 §5.3 (Observation 5.1) as the only sub-action of
the continuous Paper 25 Hopf U(1) that commutes with a real-integer
adjacency.  Together with the standard alpha/beta parity Z-strings this
gives a tapering qubit reduction of:

    Delta Q = 3                  (mode='global')
    Delta Q = 2 + n_sub_blocks   (mode='per_block', n_sb >= 1)

documented in Paper 14 §sec:hopf_tapering.

Mechanism
---------
In the native ``(n, l, m)`` orbital basis the reflection
``P : (n, l, m) -> (n, l, -m)`` is a permutation, not a Pauli Z-string.
An orthogonal rotation to the symmetric / antisymmetric combinations

    phi_+(n, l, |m|) = [phi_{n,l,+m} + phi_{n,l,-m}] / sqrt(2)   (P = +1)
    phi_-(n, l, |m|) = [phi_{n,l,+m} - phi_{n,l,-m}] / sqrt(2)   (P = -1)

turns ``P`` into a number operator counting electrons on antisymmetric
spatial orbitals, which under Jordan-Wigner becomes a Pauli Z-string.

In the standard composed builder the ERI tensor is block-diagonal across
sub-blocks (``composed_qubit.py:727-745``), so the ``m_l -> -m_l``
reflection acts independently within each sub-block, giving one
``Z_2`` stabiliser ``P_i`` per sub-block with at least one ``l >= 1``
shell.  The 'per_block' mode replaces the global ``P = prod_i P_i``
with the per-sub-block components, which are stronger because they
share their qubit-reduction power as ``2 + n_sb`` rather than ``3``.

This module is the production version of the same-day sprint drivers
``debug/sprint_z2_tapering.py`` (global) and
``debug/sprint_z2_per_subblock.py`` (per-sub-block).  See the sprint
memos for the verification panel (37/37 molecules, machine-precision
spectrum preservation on 6/37 Q<=20 cases).

Public API
----------
- :func:`build_pm_rotation` -- orthogonal rotation to P-eigenbasis.
- :func:`rotate_h1_eri` -- apply rotation to one- and two-body integrals.
- :func:`build_stabilizers` -- emit Z-string stabilisers for tapering.
- :func:`apply_hopf_tapering` -- end-to-end tapering of a JW QubitOperator.
- :func:`hopf_tapered_from_spec` -- full pipeline starting from a
  MolecularSpec.

Notes
-----
This module requires ``openfermion`` (for ``QubitOperator``,
``jordan_wigner``, and ``taper_off_qubits``) and ``numpy``.

When the underlying composed builder has non-zero cross-block ERIs
(``cross_block_h1=True``; or balanced-coupled multipole), the
per-sub-block ``P_i`` may not commute with the cross-block coupling.
This module's :func:`build_stabilizers` audits commutation per ``P_i``
and drops any that fail, so the construction degrades gracefully to a
smaller ``Delta Q`` rather than producing an incorrect spectrum.
"""

from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple

import numpy as np

try:
    from openfermion import QubitOperator, jordan_wigner  # noqa: F401
    from openfermion.transforms import taper_off_qubits
    from openfermion.utils import commutator
    _HAS_OPENFERMION = True
except ImportError:  # pragma: no cover
    _HAS_OPENFERMION = False


# ---------------------------------------------------------------------------
# Orbital-table -> P-eigenbasis rotation
# ---------------------------------------------------------------------------

def build_pm_rotation(
    orbital_table: Sequence[Tuple[Any, int, int, int]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Build the orthogonal rotation U mapping native ``(n, l, m)`` orbitals
    to the symmetric / antisymmetric eigenbasis of the per-sub-block
    ``m -> -m`` reflection.

    Parameters
    ----------
    orbital_table : sequence of ``(sub_block_key, n, l, m)`` tuples
        The native orbital ordering.  ``sub_block_key`` must be hashable
        and distinguish orbitals across different sub-blocks (the
        standard composed builder uses ``(blk_idx, label, side)``).

    Returns
    -------
    U : ``(M, M)`` ndarray
        Orthogonal rotation, ``U U^T = I``.  Row index ``i`` corresponds
        to native orbital ``i``; row ``i`` is the symmetric combination
        if ``parity[i] = +1``, antisymmetric if ``parity[i] = -1``.
    parity : ``(M,)`` int ndarray
        ``+1`` if the rotated orbital is in the symmetric ``P = +1``
        sector, ``-1`` if antisymmetric.  Orbitals with ``m = 0`` are
        always ``+1`` (P-fixed).
    """
    M = len(orbital_table)
    U = np.zeros((M, M))
    parity = np.zeros(M, dtype=int)

    idx_of = {orb: i for i, orb in enumerate(orbital_table)}
    handled = [False] * M

    for i, (sb, n, l, m) in enumerate(orbital_table):
        if handled[i]:
            continue
        if m == 0:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        partner_key = (sb, n, l, -m)
        if partner_key not in idx_of:
            raise RuntimeError(
                f"orbital {orbital_table[i]} has no m -> -m partner in "
                f"sub-block {sb!r}"
            )
        j = idx_of[partner_key]
        if j == i:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        if m > 0:
            i_plus, i_minus = i, j
        else:
            i_plus, i_minus = j, i
        inv_sqrt2 = 1.0 / np.sqrt(2.0)
        U[i_plus, i_plus] = inv_sqrt2
        U[i_plus, i_minus] = inv_sqrt2
        U[i_minus, i_plus] = inv_sqrt2
        U[i_minus, i_minus] = -inv_sqrt2
        parity[i_plus] = +1
        parity[i_minus] = -1
        handled[i_plus] = True
        handled[i_minus] = True

    err = float(np.max(np.abs(U @ U.T - np.eye(M))))
    if err > 1e-12:
        raise RuntimeError(f"rotation is not orthogonal: err={err:.2e}")

    return U, parity


def rotate_h1_eri(
    h1: np.ndarray, eri: np.ndarray, U: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply an orthogonal basis rotation to one- and two-body integrals.

    Parameters
    ----------
    h1 : ``(M, M)`` ndarray
        One-body integrals (real symmetric).
    eri : ``(M, M, M, M)`` ndarray
        Two-body integrals in chemist notation ``(pq|rs)``.
    U : ``(M, M)`` ndarray
        Orthogonal rotation from :func:`build_pm_rotation`.

    Returns
    -------
    h1_rot, eri_rot : ndarray, ndarray
        Integrals in the rotated (P-eigenbasis) ordering.
    """
    h1_rot = U @ h1 @ U.T
    eri_rot = np.einsum(
        'pa,qb,rc,sd,abcd->pqrs', U, U, U, U, eri, optimize='optimal',
    )
    return h1_rot, eri_rot


# ---------------------------------------------------------------------------
# Stabiliser construction
# ---------------------------------------------------------------------------

def _z_string_from_indices(spin_orbital_indices: Sequence[int]) -> "QubitOperator":
    op = QubitOperator(())
    for q in spin_orbital_indices:
        op *= QubitOperator(((int(q), 'Z'),))
    return op


def build_stabilizers(
    parity: np.ndarray,
    orbital_table: Optional[Sequence[Tuple[Any, int, int, int]]] = None,
    mode: str = 'per_block',
) -> Tuple["QubitOperator", "QubitOperator", List["QubitOperator"]]:
    """Construct Z-string stabilisers for Hopf tapering.

    Parameters
    ----------
    parity : ``(M,)`` int ndarray
        Per-orbital P-parity from :func:`build_pm_rotation` (``+1`` or
        ``-1``).
    orbital_table : sequence of ``(sub_block_key, n, l, m)``, optional
        Required when ``mode='per_block'``.  Used to group antisymmetric
        spin-orbitals by sub-block.  When ``mode='global'`` the orbital
        table is ignored.
    mode : ``'global'`` or ``'per_block'``
        ``'global'`` returns one combined ``P = prod_i P_i`` (Delta Q = 3
        total with alpha/beta parity).  ``'per_block'`` returns one
        ``P_i`` per sub-block that contains at least one antisymmetric
        orbital (Delta Q = 2 + n_sb).

    Returns
    -------
    Z_alpha, Z_beta : QubitOperator
        Standard alpha/beta parity Z-strings (always present).
    P_list : list of QubitOperator
        Either ``[P_global]`` (``mode='global'``) or ``[P_1, ..., P_{n_sb}]``
        (``mode='per_block'``).  Empty list if no antisymmetric orbital
        exists in the basis (e.g. He at ``max_n=2`` with only s-orbitals).
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for build_stabilizers. "
            "Install with: pip install openfermion"
        )
    if mode not in ('global', 'per_block'):
        raise ValueError(f"mode must be 'global' or 'per_block', got {mode!r}")
    if mode == 'per_block' and orbital_table is None:
        raise ValueError("orbital_table is required when mode='per_block'")

    M = len(parity)
    parity = np.asarray(parity)

    # alpha-parity Z_alpha = prod_p Z_{2p}, beta-parity Z_beta = prod_p Z_{2p+1}
    Z_alpha = _z_string_from_indices(2 * p for p in range(M))
    Z_beta = _z_string_from_indices(2 * p + 1 for p in range(M))

    antisym_indices = [p for p in range(M) if parity[p] == -1]
    if len(antisym_indices) == 0:
        return Z_alpha, Z_beta, []

    if mode == 'global':
        spin_qubits: List[int] = []
        for p in antisym_indices:
            spin_qubits.append(2 * p)
            spin_qubits.append(2 * p + 1)
        P_global = _z_string_from_indices(spin_qubits)
        return Z_alpha, Z_beta, [P_global]

    # mode == 'per_block': one P_i per sub-block with antisym orbitals.
    from collections import OrderedDict
    sb_to_antisym: "OrderedDict[Any, List[int]]" = OrderedDict()
    for p, (sb, n, l, m) in enumerate(orbital_table):
        sb_to_antisym.setdefault(sb, [])
        if parity[p] == -1:
            sb_to_antisym[sb].append(p)

    P_list: List["QubitOperator"] = []
    for sb, anti_orbitals in sb_to_antisym.items():
        if len(anti_orbitals) == 0:
            continue
        spin_qubits = []
        for p in anti_orbitals:
            spin_qubits.append(2 * p)
            spin_qubits.append(2 * p + 1)
        P_list.append(_z_string_from_indices(spin_qubits))

    return Z_alpha, Z_beta, P_list


def _commutes(op_a: "QubitOperator", op_b: "QubitOperator",
              atol: float = 1e-10) -> Tuple[bool, float]:
    c = commutator(op_a, op_b)
    max_coef = max((abs(v) for v in c.terms.values()), default=0.0)
    return (max_coef < atol), float(max_coef)


# ---------------------------------------------------------------------------
# End-to-end tapering of a JW QubitOperator
# ---------------------------------------------------------------------------

def apply_hopf_tapering(
    qubit_op_rotated: "QubitOperator",
    parity: np.ndarray,
    orbital_table: Optional[Sequence[Tuple[Any, int, int, int]]] = None,
    mode: str = 'per_block',
    sector_signs: Optional[Sequence[int]] = None,
    drop_noncommuting: bool = True,
    audit_tol: float = 1e-10,
) -> Tuple["QubitOperator", List[int]]:
    """Apply Hopf-U(1) tapering to a JW-encoded qubit Hamiltonian.

    The input ``qubit_op_rotated`` must already be in the P-eigenbasis
    (rotated via :func:`rotate_h1_eri` + JW).  This function builds the
    stabilisers, audits their commutation with the Hamiltonian, drops any
    that fail (defensive against cross-block extensions), and applies
    ``openfermion.transforms.taper_off_qubits``.

    Parameters
    ----------
    qubit_op_rotated : QubitOperator
        JW-encoded Hamiltonian in the P-eigenbasis.
    parity : ``(M,)`` int ndarray
        Per-orbital parity from :func:`build_pm_rotation`.
    orbital_table : sequence of ``(sub_block_key, n, l, m)``, optional
        Required when ``mode='per_block'``.
    mode : ``'global'`` or ``'per_block'``
        Tapering mode (see :func:`build_stabilizers`).
    sector_signs : sequence of ``+1`` / ``-1``, optional
        One sign per stabiliser in the order
        ``[Z_alpha, Z_beta, P_1, ..., P_k]``.  Default is all ``+1``
        which selects the fully-symmetric sector (the physical
        ground-state sector for closed-shell chemistry in the standard
        composed builder).
    drop_noncommuting : bool
        If ``True`` (default) drop any ``P_i`` that does not commute with
        the Hamiltonian at tolerance ``audit_tol``.  When the standard
        composed builder is used this should never fire; it is a
        defensive gate against cross-block extensions.
    audit_tol : float
        Tolerance for the per-stabiliser commutator audit.

    Returns
    -------
    qubit_op_tapered : QubitOperator
        Tapered Hamiltonian on ``Q_naive - n_stabilisers`` qubits.
    dropped_indices : list of int
        Indices into the un-filtered ``P_list`` that were dropped by the
        commutation audit.  Empty list when no extension breaks the
        construction.
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for apply_hopf_tapering. "
            "Install with: pip install openfermion"
        )

    Z_alpha, Z_beta, P_list = build_stabilizers(
        parity, orbital_table=orbital_table, mode=mode,
    )

    # Audit each P_i against the rotated Hamiltonian.  Z_alpha and Z_beta
    # always commute with a particle-number / S_z-conserving operator and
    # are skipped.
    kept_P: List["QubitOperator"] = []
    dropped: List[int] = []
    for i, P_i in enumerate(P_list):
        ok, _ = _commutes(qubit_op_rotated, P_i, atol=audit_tol)
        if ok or not drop_noncommuting:
            kept_P.append(P_i)
        else:
            dropped.append(i)

    stabilisers = [Z_alpha, Z_beta] + kept_P
    n_stab = len(stabilisers)

    if sector_signs is None:
        signs = [1] * n_stab
    else:
        signs = list(sector_signs)
        if len(signs) != n_stab:
            raise ValueError(
                f"sector_signs length {len(signs)} does not match number "
                f"of stabilisers {n_stab} (after dropping {len(dropped)} "
                f"non-commuting P_i)"
            )
        if any(s not in (-1, 1) for s in signs):
            raise ValueError("sector_signs entries must be +1 or -1")

    signed_stabs = [s * sgn for s, sgn in zip(stabilisers, signs)]
    qubit_op_tapered = taper_off_qubits(qubit_op_rotated, signed_stabs)
    return qubit_op_tapered, dropped


# ---------------------------------------------------------------------------
# Full pipeline from a MolecularSpec
# ---------------------------------------------------------------------------

def _enumerate_orbitals(spec: Any) -> List[Tuple[Any, int, int, int]]:
    """Mirror ``composed_qubit.build_composed_hamiltonian``'s orbital
    enumeration.  Returns the orbital_table as ``[(sb_key, n, l, m), ...]``.
    """
    from geovac.composed_qubit import _enumerate_states

    orbital_table: List[Tuple[Any, int, int, int]] = []
    for blk_idx, blk in enumerate(spec.blocks):
        l_min = getattr(blk, 'l_min', 0)
        sb_key_center = (blk_idx, blk.label, 'center')
        for (n, l, m) in _enumerate_states(blk.max_n, l_min=l_min):
            orbital_table.append((sb_key_center, n, l, m))
        if blk.has_h_partner:
            partner_max_n = (
                blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            )
            sb_key_partner = (blk_idx, blk.label, 'partner')
            for (n, l, m) in _enumerate_states(partner_max_n):
                orbital_table.append((sb_key_partner, n, l, m))
    return orbital_table


def hopf_tapered_from_spec(
    spec: Any,
    mode: str = 'per_block',
    pk_in_hamiltonian: bool = True,
    sector_signs: Optional[Sequence[int]] = None,
    verbose: bool = False,
) -> dict:
    """End-to-end Hopf-tapering pipeline starting from a MolecularSpec.

    Builds the composed Hamiltonian, rotates to the P-eigenbasis,
    JW-encodes, and applies tapering.

    Parameters
    ----------
    spec : MolecularSpec
        From ``geovac.molecular_spec`` (e.g., ``hydride_spec(8)`` for H2O).
    mode : ``'global'`` or ``'per_block'``
        Tapering mode.  See :func:`build_stabilizers`.
    pk_in_hamiltonian : bool
        Passed through to ``build_composed_hamiltonian``.
    sector_signs : sequence of ``+1`` / ``-1``, optional
        See :func:`apply_hopf_tapering`.  Default is the fully-symmetric
        sector.
    verbose : bool
        Print build progress.

    Returns
    -------
    dict with keys:
        qubit_op_naive : QubitOperator (no tapering)
        qubit_op_tapered : QubitOperator (after tapering)
        Q_naive, Q_tapered : int
        delta_Q : int (positive number of qubits removed)
        n_sub_blocks : int
        n_sub_blocks_with_antisym : int
        dropped_P_indices : list of int (empty in standard composed builder)
        orbital_table, parity, U : the intermediate objects
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for hopf_tapered_from_spec. "
            "Install with: pip install openfermion"
        )
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.qubit_encoding import build_fermion_op_from_integrals

    result = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
    )
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']

    orbital_table = _enumerate_orbitals(spec)
    M = len(orbital_table)
    if h1.shape != (M, M):
        raise RuntimeError(
            f"orbital_table length {M} does not match h1.shape {h1.shape}; "
            f"sub-block enumeration is out of sync with build_composed_hamiltonian"
        )

    # Naive JW
    fop_naive = build_fermion_op_from_integrals(h1, eri, nuc_rep)
    qop_naive = jordan_wigner(fop_naive)
    Q_naive = 2 * M

    # Rotation to P-eigenbasis
    U, parity = build_pm_rotation(orbital_table)
    h1_rot, eri_rot = rotate_h1_eri(h1, eri, U)
    fop_rot = build_fermion_op_from_integrals(h1_rot, eri_rot, nuc_rep)
    qop_rot = jordan_wigner(fop_rot)

    # Tapering
    qop_tap, dropped = apply_hopf_tapering(
        qop_rot, parity, orbital_table=orbital_table, mode=mode,
        sector_signs=sector_signs,
    )

    # Compute Q_tapered from the operator (largest qubit index + 1).
    Q_tapered = 0
    for term in qop_tap.terms:
        for q, _ in term:
            if q + 1 > Q_tapered:
                Q_tapered = q + 1

    sb_keys = sorted(set(orb[0] for orb in orbital_table), key=lambda k: str(k))
    n_sb = len(sb_keys)
    antisym_by_sb: dict = {}
    for p, (sb, n, l, m) in enumerate(orbital_table):
        antisym_by_sb.setdefault(sb, 0)
        if parity[p] == -1:
            antisym_by_sb[sb] += 1
    n_sb_anti = sum(1 for c in antisym_by_sb.values() if c > 0)

    return {
        'qubit_op_naive': qop_naive,
        'qubit_op_tapered': qop_tap,
        'Q_naive': Q_naive,
        'Q_tapered': Q_tapered,
        'delta_Q': Q_naive - Q_tapered,
        'n_sub_blocks': n_sb,
        'n_sub_blocks_with_antisym': n_sb_anti,
        'dropped_P_indices': dropped,
        'mode': mode,
        'orbital_table': orbital_table,
        'parity': parity,
        'U': U,
    }
