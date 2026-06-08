"""
Symmetry-adapted basis: per-orbital Z_2^3 irrep labelling.
==========================================================

Combines the three commuting Z_2 symmetries of the GeoVac composed-chemistry
Hamiltonian into one basis rotation U_total = U_swap @ U_hopf and tags each
orbital with its triple-Z_2 irrep label

    (p_Hopf, p_ell, p_swap) in {+1, -1}^3.

ℓ-parity is a sign on l-odd orbitals (no rotation), and inversion is built as
swap × (-1)^l (not an independent rotation when the other three are present).

Public API
----------
- :func:`build_symmetry_adapted_rotation` — returns (U_total, sector_labels)
- :func:`decompose_hamiltonian_by_sector` — rotated h1, eri grouped by sector
- :func:`find_hidden_z2_in_sector` — scan each sector for accidental
  block-diagonal sub-structure (candidate hidden Z_2) and audit against H
- :func:`extended_plus_hidden_tapered_from_spec` — v3.89.0 extended tapering
  plus any validated hidden stabilizers (additive shim on
  ``extended_tapering.extended_tapered_from_spec``)

Honest scope (sprint Symmetry-Adapted Basis, 2026-06-07)
--------------------------------------------------------
The Hopf and atom-swap rotations COMMUTE STRUCTURALLY because the Hopf
rotation only mixes orbitals within a sub-block at fixed (n, l) (it acts on
the m index), while the swap rotation pairs equivalent sub-blocks at fixed
(n, l, m). They act on orthogonal subspaces of orbital-index space, so
U_swap and U_hopf commute as matrices.  ℓ-parity is then a signed projection
onto the (-1)^l eigenspace and is naturally compatible.

The hidden-Z_2 scan is verdict-driven:
  - POSITIVE: a sector-block diagonalisation reveals two further sub-blocks
    not coupled by h1 OR by eri, AND a Z-string stabilizer with the matching
    parity commutes with the full H (h1 + eri + ecore).
  - NEGATIVE: every sector block is fully populated (no accidental splits).
  - PARTIAL: pattern suggests a new structure but the candidate stabilizer is
    not of the form ``∏ Z_q``.

This module does NOT modify ``geovac/extended_tapering.py``.
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    from openfermion import QubitOperator, jordan_wigner  # noqa: F401
    from openfermion.utils import commutator
    _HAS_OPENFERMION = True
except ImportError:  # pragma: no cover
    _HAS_OPENFERMION = False


# ---------------------------------------------------------------------------
# Helpers (shared shape with z2_tapering / extended_tapering)
# ---------------------------------------------------------------------------

def _z_string_from_indices(spin_orbital_indices: Sequence[int]) -> "QubitOperator":
    op = QubitOperator(())
    for q in spin_orbital_indices:
        op *= QubitOperator(((int(q), 'Z'),))
    return op


def _commutes(op_a: "QubitOperator", op_b: "QubitOperator",
              atol: float = 1e-10) -> Tuple[bool, float]:
    c = commutator(op_a, op_b)
    max_coef = max((abs(v) for v in c.terms.values()), default=0.0)
    return max_coef <= atol, float(max_coef)


# ---------------------------------------------------------------------------
# Combined rotation + sector labels
# ---------------------------------------------------------------------------

def build_symmetry_adapted_rotation(
    spec: Any,
    nuclei: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[Tuple[Any, int, int, int]],
           Dict[str, Any]]:
    """Build U_total = U_swap @ U_hopf and per-orbital Z_2^3 sector labels.

    Parameters
    ----------
    spec : MolecularSpec
    nuclei : list of dicts, optional
        Required for atom-swap.  If omitted, U_swap = I and p_swap = +1 for
        every orbital.

    Returns
    -------
    U_total : (M, M) ndarray
        Orthogonal rotation, U_total @ U_total.T = I to < 1e-12.
    sector_labels : (M, 3) int ndarray
        Each row is ``(p_Hopf, p_ell, p_swap)`` for orbital ``i`` in the
        ROTATED basis.
    parity_hopf, parity_swap : (M,) int ndarrays
        Raw parities returned by :func:`build_pm_rotation` and
        :func:`build_atom_swap_rotation_and_stabilizers` respectively.
    orbital_table : list
        From :func:`_enumerate_orbitals(spec)`; orbital ``i`` has unrotated
        index ``i``.  ℓ-parity uses ``orbital_table[i][2]``.
    diagnostics : dict
        Keys: ``U_hopf``, ``U_swap``, ``swap_diag``, ``orthogonality_err``.

    Notes
    -----
    U_total acts on the orbital index space (M-by-M), and the sector-labels
    are read off in the ROTATED basis.  Because ℓ-parity (-1)^l is a per-
    orbital sign that depends only on the l label (which is preserved by
    BOTH U_hopf and U_swap — Hopf mixes m at fixed n, l; swap pairs sub-
    blocks at fixed n, l, m), the same (-1)^l applies to row i of the
    rotated basis as to orbital i.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")

    from geovac.z2_tapering import _enumerate_orbitals, build_pm_rotation
    from geovac.extended_tapering import (
        build_atom_swap_rotation_and_stabilizers,
    )

    orbital_table = _enumerate_orbitals(spec)
    M = len(orbital_table)

    U_hopf, parity_hopf = build_pm_rotation(orbital_table)

    if nuclei is None or len(nuclei) == 0:
        U_swap = np.eye(M)
        parity_swap = np.ones(M, dtype=int)
        swap_diag: List[Dict[str, Any]] = []
    else:
        U_swap, parity_swap, swap_diag = (
            build_atom_swap_rotation_and_stabilizers(
                spec, orbital_table, nuclei,
            )
        )

    # Hopf acts on m-pairs within (sub_block, n, l); swap acts on sub-block
    # pairs at fixed (n, l, m).  Their supports overlap structurally but the
    # algebraic action commutes because they act on different *labels* of the
    # same index set (m-label vs. sub-block-label).
    U_total = U_swap @ U_hopf

    err = float(np.max(np.abs(U_total @ U_total.T - np.eye(M))))
    if err > 1e-12:
        raise RuntimeError(
            f"combined rotation not orthogonal: max|UU^T - I| = {err:.2e}"
        )

    # Per-orbital ℓ-parity: (-1)^l of the original orbital.
    parity_ell = np.array(
        [(-1) ** int(orbital_table[i][2]) for i in range(M)], dtype=int
    )

    sector_labels = np.stack(
        [parity_hopf, parity_ell, parity_swap], axis=1
    ).astype(int)

    diagnostics = {
        'U_hopf': U_hopf,
        'U_swap': U_swap,
        'swap_diag': swap_diag,
        'orthogonality_err': err,
        'commute_hopf_swap': float(
            np.max(np.abs(U_hopf @ U_swap - U_swap @ U_hopf))
        ),
    }

    return U_total, sector_labels, parity_hopf, parity_swap, orbital_table, diagnostics


# ---------------------------------------------------------------------------
# Sector decomposition
# ---------------------------------------------------------------------------

def decompose_hamiltonian_by_sector(
    spec: Any,
    nuclei: Optional[List[Dict[str, Any]]] = None,
    builder: str = 'composed',
    pk_in_hamiltonian: bool = True,
    R: Optional[float] = None,
    verbose: bool = False,
) -> Dict[Tuple[int, ...], Dict[str, Any]]:
    """Rotate h1, eri to the symmetry-adapted basis and group orbital
    indices by Z_2^3 sector label.

    Parameters
    ----------
    spec : MolecularSpec
    nuclei : list of dicts, optional
    builder : ``'composed'`` or ``'balanced'``

    Returns
    -------
    dict mapping ``(p_Hopf, p_ell, p_swap) -> {...}``
        Each value contains
          ``orbital_indices`` : list of int
          ``h1_block``        : (k, k) ndarray
          ``eri_block``       : (k, k, k, k) ndarray (intra-sector slice)
          ``dim``             : int (= k)
        Plus the meta key ``('__meta__',)`` carrying
          ``U_total``, ``sector_labels``, ``orbital_table``,
          ``h1_rot``, ``eri_rot``, ``nuc_rep``, ``max_offsector_h1``,
          ``max_offsector_eri``.
    """
    from geovac.z2_tapering import rotate_h1_eri

    if builder == 'composed':
        from geovac.composed_qubit import build_composed_hamiltonian
        result = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
        )
    elif builder == 'balanced':
        from geovac.balanced_coupled import build_balanced_hamiltonian
        result = build_balanced_hamiltonian(
            spec, R=R, nuclei=nuclei, verbose=verbose,
        )
    else:
        raise ValueError(f"builder must be 'composed' or 'balanced'")

    h1 = np.asarray(result['h1'])
    eri = np.asarray(result['eri'])
    nuc_rep = float(result['nuclear_repulsion'])

    U_total, sector_labels, parity_hopf, parity_swap, orbital_table, diag = (
        build_symmetry_adapted_rotation(spec, nuclei)
    )

    h1_rot, eri_rot = rotate_h1_eri(h1, eri, U_total)

    M = len(orbital_table)
    sector_of: Dict[Tuple[int, ...], List[int]] = OrderedDict()
    for i in range(M):
        key = tuple(int(x) for x in sector_labels[i])
        sector_of.setdefault(key, []).append(i)

    # Cross-sector residuals (block-diagonality gate)
    max_offsector_h1 = 0.0
    for key_a, idx_a in sector_of.items():
        for key_b, idx_b in sector_of.items():
            if key_a == key_b:
                continue
            block = h1_rot[np.ix_(idx_a, idx_b)]
            if block.size > 0:
                max_offsector_h1 = max(
                    max_offsector_h1, float(np.max(np.abs(block)))
                )

    # ERI off-sector check: every 4-index term (pq|rs) must have the product
    # of sector triples preserve symmetry pq vs rs (chemist notation).
    max_offsector_eri = 0.0
    keys_arr = np.array([tuple(int(x) for x in s) for s in sector_labels])
    # For chemist ERI (pq|rs), nonzero entries must have
    # sector_labels[p] * sector_labels[q] * sector_labels[r] * sector_labels[s] == (+1,+1,+1)
    # i.e., for each of the three Z_2 axes, the product of the four labels is +1.
    prod = keys_arr[:, None, None, None, :] * keys_arr[None, :, None, None, :] * \
           keys_arr[None, None, :, None, :] * keys_arr[None, None, None, :, :]
    # Where any axis is -1, ERI must be zero
    mask_offsector = np.any(prod == -1, axis=-1)
    if np.any(mask_offsector):
        max_offsector_eri = float(np.max(np.abs(eri_rot[mask_offsector])))

    output: Dict[Tuple[int, ...], Dict[str, Any]] = OrderedDict()
    for key, idx in sector_of.items():
        idx_arr = np.array(idx, dtype=int)
        h1_block = h1_rot[np.ix_(idx_arr, idx_arr)]
        eri_block = eri_rot[
            np.ix_(idx_arr, idx_arr, idx_arr, idx_arr)
        ]
        output[key] = {
            'orbital_indices': list(int(x) for x in idx),
            'h1_block': h1_block,
            'eri_block': eri_block,
            'dim': int(len(idx)),
        }

    output[('__meta__',)] = {
        'U_total': U_total,
        'sector_labels': sector_labels,
        'parity_hopf': parity_hopf,
        'parity_swap': parity_swap,
        'orbital_table': orbital_table,
        'h1_rot': h1_rot,
        'eri_rot': eri_rot,
        'nuc_rep': nuc_rep,
        'max_offsector_h1': max_offsector_h1,
        'max_offsector_eri': max_offsector_eri,
        'commute_hopf_swap': diag['commute_hopf_swap'],
    }
    return output


# ---------------------------------------------------------------------------
# Hidden-Z_2 detector
# ---------------------------------------------------------------------------

def _connected_components_in_sector(
    h1_block: np.ndarray,
    eri_block: np.ndarray,
    atol: float = 1e-10,
) -> List[List[int]]:
    """Identify connected components within a sector block.

    Two orbitals i, j (within the sector) are connected if
      h1_block[i, j] is nonzero  OR
      eri_block has any nonzero entry with both i and j in the
      (p, q, r, s) tuple.

    Returns a list of components; each component is a list of LOCAL indices
    (0..dim-1) into the sector's orbital_indices list.
    """
    k = h1_block.shape[0]
    if k <= 1:
        return [[i] for i in range(k)]

    adj = np.zeros((k, k), dtype=bool)
    np.fill_diagonal(adj, True)
    # h1 connections
    adj |= np.abs(h1_block) > atol
    # ERI connections: (p, q, r, s) all live in same sector by construction.
    # Project onto pair-co-occurrence.
    eri_abs = np.abs(eri_block)
    # Mark (p, q) connected if any (r, s) makes eri large
    pq_active = (eri_abs > atol).any(axis=(2, 3))
    rs_active = (eri_abs > atol).any(axis=(0, 1))
    adj |= pq_active | pq_active.T
    adj |= rs_active | rs_active.T

    # Connect pairs that co-occur in any active ERI quartet
    # If eri[p, q, r, s] is nonzero, then p–q, r–s, p–r, q–s, p–s, q–r are all
    # connected (they appear in the same term).
    nz = np.argwhere(eri_abs > atol)
    for (p, q, r, s) in nz:
        adj[p, q] = adj[q, p] = True
        adj[r, s] = adj[s, r] = True
        adj[p, r] = adj[r, p] = True
        adj[q, s] = adj[s, q] = True
        adj[p, s] = adj[s, p] = True
        adj[q, r] = adj[r, q] = True

    # BFS for components
    visited = [False] * k
    comps: List[List[int]] = []
    for start in range(k):
        if visited[start]:
            continue
        stack = [start]
        comp = []
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for w in range(k):
                if adj[v, w] and not visited[w]:
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def find_hidden_z2_in_sector(
    spec: Any,
    nuclei: Optional[List[Dict[str, Any]]] = None,
    builder: str = 'composed',
    pk_in_hamiltonian: bool = True,
    R: Optional[float] = None,
    audit_tol: float = 1e-10,
    verbose: bool = False,
) -> List[Dict[str, Any]]:
    """Scan each sector for accidental block-diagonal sub-structure
    suggesting a hidden Z_2 symmetry, and audit candidates against H.

    For each sector with dim >= 2, compute connected components of its
    intra-sector h1 + eri.  If there are >= 2 components AND a Z-string
    stabilizer that puts opposite signs on different components commutes
    with the full Hamiltonian, emit a candidate.

    Returns
    -------
    list of dicts (one per candidate, may be empty):
      'sector'                : tuple (p_Hopf, p_ell, p_swap)
      'orbital_indices'       : list of int — orbitals in this sector
      'components'            : list of list-of-int — each component is
                                global orbital indices for one sub-block
      'candidate_stabilizer'  : QubitOperator (Z-string)
      'audit_residual'        : float — max commutator entry vs full H
      'is_valid'              : bool — True iff audit_residual <= audit_tol
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")

    from geovac.qubit_encoding import build_fermion_op_from_integrals

    sectors = decompose_hamiltonian_by_sector(
        spec, nuclei=nuclei, builder=builder,
        pk_in_hamiltonian=pk_in_hamiltonian, R=R, verbose=verbose,
    )

    meta = sectors[('__meta__',)]
    h1_rot = meta['h1_rot']
    eri_rot = meta['eri_rot']
    nuc_rep = meta['nuc_rep']
    M = h1_rot.shape[0]

    # Build the rotated Hamiltonian as a QubitOperator once for auditing.
    fop_rot = build_fermion_op_from_integrals(h1_rot, eri_rot, nuc_rep)
    qop_rot = jordan_wigner(fop_rot)

    candidates: List[Dict[str, Any]] = []
    for key, info in sectors.items():
        if key == ('__meta__',):
            continue
        if info['dim'] < 2:
            continue
        comps = _connected_components_in_sector(
            info['h1_block'], info['eri_block'], atol=audit_tol,
        )
        if len(comps) < 2:
            continue

        # We have at least two components. Build a Z-string stabilizer that
        # flips the sign of one component and not the other(s). If there are
        # exactly 2 components we get one candidate; if more, pair them off
        # as (comp_0 vs rest).
        # Map local-component indices → global orbital indices in rotated
        # basis (= same index space as orbital_indices in the sector dict).
        global_idx = info['orbital_indices']

        # Try every component as the "−1" group
        for c_idx, comp in enumerate(comps):
            anti_orbitals = [global_idx[i] for i in comp]
            spin_qubits: List[int] = []
            for p in anti_orbitals:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            P = _z_string_from_indices(spin_qubits)
            ok, residual = _commutes(qop_rot, P, atol=audit_tol)
            candidates.append({
                'sector': key,
                'component_index': c_idx,
                'n_components': len(comps),
                'orbital_indices': global_idx,
                'components': [
                    [global_idx[i] for i in cc] for cc in comps
                ],
                'candidate_stabilizer': P,
                'audit_residual': residual,
                'is_valid': bool(ok),
            })
            if verbose:
                print(f"  sector {key}, comp {c_idx}: residual {residual:.2e}, "
                      f"valid={ok}")

    return candidates


# ---------------------------------------------------------------------------
# End-to-end pipeline (additive on extended_tapered_from_spec)
# ---------------------------------------------------------------------------

def extended_plus_hidden_tapered_from_spec(
    spec: Any,
    nuclei: Optional[List[Dict[str, Any]]] = None,
    use_hopf: bool = True,
    use_ell_parity: bool = True,
    use_atom_swap: bool = False,
    use_inversion: bool = False,
    pk_in_hamiltonian: bool = True,
    builder: str = 'composed',
    R: Optional[float] = None,
    audit_tol: float = 1e-10,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run the v3.89.0 extended tapering pipeline, then look for additional
    hidden-Z_2 stabilizers in the rotated sectors and append the validated
    ones to the tapering stack.

    Returns
    -------
    dict with the standard ``extended_tapered_from_spec`` keys plus
      ``hidden_z2_candidates`` : list (from :func:`find_hidden_z2_in_sector`)
      ``hidden_z2_kept``       : list of dict for validated additions
      ``Q_tapered_hidden``     : int — Q after extended + hidden
      ``delta_Q_hidden``       : int — extra qubits saved by hidden
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")

    from geovac.extended_tapering import extended_tapered_from_spec
    from openfermion.transforms import taper_off_qubits

    base = extended_tapered_from_spec(
        spec,
        use_hopf=use_hopf, use_ell_parity=use_ell_parity,
        use_atom_swap=use_atom_swap, use_inversion=use_inversion,
        pk_in_hamiltonian=pk_in_hamiltonian, builder=builder, R=R,
        nuclei=nuclei, verbose=verbose,
    )

    candidates = find_hidden_z2_in_sector(
        spec, nuclei=nuclei, builder=builder,
        pk_in_hamiltonian=pk_in_hamiltonian, R=R,
        audit_tol=audit_tol, verbose=verbose,
    )

    # Filter to validated candidates and de-duplicate by the (component_index,
    # sector) pair: when n_components > 2, multiple candidates are emitted
    # (one per component-as-anti); pairs with the same set of anti-orbitals
    # produce the same stabilizer. We just take the first valid representative
    # per (sector,c_idx) pair and rely on the audit gate for correctness.
    base_qop = base['qubit_op_tapered']

    # We append further stabilizers ON TOP OF the already-tapered operator.
    # The hidden candidate Z-strings live on the rotated (pre-taper) qubits;
    # since taper_off_qubits removes a fixed set of qubits, the surviving
    # qubits keep their indices contiguous in the operator returned by
    # extended_tapered_from_spec (openfermion's convention).  We therefore
    # cannot directly re-use the hidden Z-strings without remapping qubit
    # indices.  Simpler approach: redo tapering on the original rotated
    # operator with the union of extended + hidden stabilizers.

    # Rebuild the rotated operator and stabilizers exactly as
    # extended_tapered_from_spec did, then add validated hidden stabilizers
    # and call taper_off_qubits once at the end.
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from geovac.z2_tapering import rotate_h1_eri
    from geovac.extended_tapering import (
        apply_extended_tapering,
    )

    if builder == 'composed':
        result = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=False,
        )
    else:
        result = build_balanced_hamiltonian(
            spec, R=R, nuclei=nuclei, verbose=False,
        )
    h1 = np.asarray(result['h1'])
    eri = np.asarray(result['eri'])
    nuc_rep = float(result['nuclear_repulsion'])

    U_total = base['U_total']
    h1_rot, eri_rot = rotate_h1_eri(h1, eri, U_total)
    fop_rot = build_fermion_op_from_integrals(h1_rot, eri_rot, nuc_rep)
    qop_rot = jordan_wigner(fop_rot)

    # Run extended_tapering audit + GF(2) drop pass but with the extra
    # hidden Z-strings included.  Easiest path: call apply_extended_tapering
    # to build the standard stabilizer list, then mutate the kept list with
    # validated hidden stabilizers, then re-run linear-independence audit and
    # final taper.
    qop_tap_ext, audit_ext = apply_extended_tapering(
        qop_rot,
        parity_hopf=base['parity_hopf'],
        orbital_table=base['orbital_table'],
        parity_swap=base['parity_swap'],
        swap_diag=base['swap_diag'],
        use_hopf_per_block=use_hopf,
        use_ell_parity=use_ell_parity,
        use_atom_swap=use_atom_swap,
        use_inversion=use_inversion,
        inversion_stabilizers=None,
        verbose=verbose,
    )

    # Build the union stabilizer list manually so we can pass it to
    # taper_off_qubits.
    from geovac.extended_tapering import (
        _z_string_from_indices as _zs_ext,
        _drop_linearly_dependent_z_stabilizers,
    )
    M = len(base['orbital_table'])
    n_qubits = 2 * M
    Z_alpha = _zs_ext(2 * p for p in range(M))
    Z_beta = _zs_ext(2 * p + 1 for p in range(M))
    stabilizers: List[Tuple[str, "QubitOperator"]] = [
        ('alpha', Z_alpha), ('beta', Z_beta),
    ]

    # Replicate the same per_block stabilizer construction used internally.
    # Easiest: walk through audit_ext['kinds_kept'][2:] (skipping alpha/beta)
    # and re-emit by re-running build_stabilizers from the extended pipeline.
    # apply_extended_tapering already filtered + dependency-dropped the list,
    # but we want to take its 'kinds_kept' set as a starting point.
    # Pragmatic: we will redo the stabilizer build manually here.
    # However, the cleanest implementation reuses the internal builder by
    # calling apply_extended_tapering twice — once for the standard set,
    # then once with hidden additions stitched on top.

    # The hidden candidate Z-strings live on the ROTATED basis qubits.  We
    # add them after the standard ones and re-run the GF(2) drop + final
    # taper.  We use apply_extended_tapering's drop logic by manually
    # building the union stabilizer list and calling taper_off_qubits.

    # Reproduce the union list shape by re-emitting per-block / atom-swap /
    # inversion stabilizers in the same order as apply_extended_tapering.
    joint_sb_map: Dict[Any, Any] = {}
    if base.get('swap_diag'):
        from geovac.extended_tapering import _build_joint_sb_map
        joint_sb_map = _build_joint_sb_map(
            base['orbital_table'], base['swap_diag'],
        )

    def _joint_key(sb: Any) -> Any:
        return joint_sb_map.get(sb, sb)

    if use_hopf:
        sb_to_anti: "OrderedDict[Any, List[int]]" = OrderedDict()
        for p, (sb, n, l, m) in enumerate(base['orbital_table']):
            jsb = _joint_key(sb)
            sb_to_anti.setdefault(jsb, [])
            if base['parity_hopf'][p] == -1:
                sb_to_anti[jsb].append(p)
        for jsb, anti in sb_to_anti.items():
            if not anti:
                continue
            spin_qubits = []
            for p in anti:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(
                (f'hopf:{jsb}', _zs_ext(spin_qubits))
            )

    if use_ell_parity:
        sb_to_lodd: "OrderedDict[Any, List[int]]" = OrderedDict()
        for p, (sb, n, l, m) in enumerate(base['orbital_table']):
            jsb = _joint_key(sb)
            sb_to_lodd.setdefault(jsb, [])
            if l % 2 == 1:
                sb_to_lodd[jsb].append(p)
        for jsb, l_odd in sb_to_lodd.items():
            if not l_odd:
                continue
            spin_qubits = []
            for p in l_odd:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(
                (f'ell_parity:{jsb}', _zs_ext(spin_qubits))
            )

    if use_atom_swap and base['parity_swap'] is not None:
        anti = [i for i in range(M) if base['parity_swap'][i] == -1]
        if anti:
            spin_qubits = []
            for p in anti:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(
                ('atom_swap', _zs_ext(spin_qubits))
            )

    # Append validated hidden stabilizers (skip duplicates by orbital set).
    valid = [c for c in candidates if c['is_valid']]
    # Deduplicate by frozenset of (sector, component-as-anti-orbital set)
    seen_orbset = set()
    hidden_kept: List[Dict[str, Any]] = []
    for c in valid:
        # The c['component_index']-th component is the −1 group; map to
        # orbital indices.
        anti_orbs = c['components'][c['component_index']]
        key_set = frozenset(anti_orbs)
        if key_set in seen_orbset:
            continue
        # Also skip the trivial "all orbitals in sector" (which equals an
        # already-present stabilizer)
        all_orbs = frozenset(c['orbital_indices'])
        if key_set == all_orbs:
            continue
        seen_orbset.add(key_set)
        spin_qubits = []
        for p in anti_orbs:
            spin_qubits.append(2 * p)
            spin_qubits.append(2 * p + 1)
        P = _zs_ext(spin_qubits)
        stabilizers.append(
            (f'hidden_z2:sector={c["sector"]}:c={c["component_index"]}', P)
        )
        hidden_kept.append({
            'sector': c['sector'],
            'component_index': c['component_index'],
            'anti_orbitals': anti_orbs,
            'audit_residual': c['audit_residual'],
        })

    # Drop linearly dependent Z-string stabilizers
    kept_with_kind, dropped = _drop_linearly_dependent_z_stabilizers(
        stabilizers, n_qubits=n_qubits,
    )

    # Audit each non-alpha/beta against the rotated operator and drop
    # non-commuting ones (defensive — the hidden candidates already passed
    # auditing, but the union may have dependencies we haven't checked
    # numerically).
    final_kept: List[Tuple[str, "QubitOperator"]] = []
    audit_failures: List[Tuple[str, float]] = []
    for (kind, P_i) in kept_with_kind:
        if kind in ('alpha', 'beta'):
            final_kept.append((kind, P_i))
            continue
        ok, r = _commutes(qop_rot, P_i, atol=audit_tol)
        if ok:
            final_kept.append((kind, P_i))
        else:
            audit_failures.append((kind, r))

    kept_ops = [op for (_, op) in final_kept]
    qop_tap = taper_off_qubits(qop_rot, kept_ops)

    Q_naive = base['Q_naive']
    Q_tapered_hidden = 0
    for term in qop_tap.terms:
        for q, _ in term:
            if q + 1 > Q_tapered_hidden:
                Q_tapered_hidden = q + 1

    return {
        **base,
        'qubit_op_tapered_plus_hidden': qop_tap,
        'Q_tapered_hidden': int(Q_tapered_hidden),
        'delta_Q_hidden': int(base['Q_tapered'] - Q_tapered_hidden)
                          if Q_tapered_hidden > 0 else int(base['delta_Q']),
        'hidden_z2_candidates': candidates,
        'hidden_z2_kept': hidden_kept,
        'hidden_z2_audit_failures': audit_failures,
        'hidden_kinds_kept': [k for (k, _) in final_kept],
    }
