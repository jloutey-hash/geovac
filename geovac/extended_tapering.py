"""
Extended Z₂ tapering for GeoVac composed qubit Hamiltonians.
============================================================

Production module extending the Hopf-U(1) ``m_l -> -m_l`` Z₂ tapering
(`geovac/z2_tapering.py`, v3.52.0) with three additional Z₂ symmetries
identified by the sprint M-vS Gauge + symmetry-explorer cycle (v3.88.0):

  - **ℓ-parity Z₂** (Candidate 4): per sub-block, the operator
    ``P_ℓ = (-1)^(N_lodd)`` where ``N_lodd`` counts electrons in
    l-odd orbitals.  Commutes with H by Gaunt selection rule
    (ERI requires ``(l_a+l_b+l_c+l_d)`` even; within-block
    hydrogenic h1 is l-diagonal).  Per sub-block ΔQ = ``+n_sub_blocks``.

  - **Atom-swap Z₂** (Candidate 1): for molecules with equivalent
    atoms (BeH₂ H_1 ↔ H_2, NH₃ H_1 ↔ H_2, CH₄ pairwise, N₂/F₂ atom
    swap, etc.), the permutation of equivalent sub-blocks generates
    a Z₂ symmetry.  Klein four-group subgroup for CH₄/SiH₄/GeH₄.
    Adds ΔQ = ``+1`` (or ``+2`` for tetrahedral) per molecule.

  - **Inversion Z₂** (Candidate 2): for centrosymmetric molecules
    (BeH₂ linear, N₂, F₂, MgH₂, CaH₂, C₂H₂, C₂H₆), spatial inversion
    combines the atom-swap with a ``(-1)^l`` sign factor.  Independent
    Z₂ generator above atom-swap.  Adds ΔQ = ``+1``.

The "extended" tapering applies all of these on top of the
existing Hopf ``m_l -> -m_l`` rotation and the standard alpha/beta
parity Z-strings.  Each stabilizer is audited for commutation with
the actual Hamiltonian before tapering; non-commuting stabilizers are
gracefully dropped (the audit gate from `apply_hopf_tapering`).

Production API
--------------
- :func:`build_ell_parity_stabilizers`
- :func:`build_atom_swap_rotation_and_stabilizers`
- :func:`build_inversion_rotation_and_stabilizers`
- :func:`apply_extended_tapering`
- :func:`extended_tapered_from_spec`

Backward compatibility
----------------------
The existing :func:`geovac.z2_tapering.hopf_tapered_from_spec` is
unchanged; this module is additive.  Setting all extended kwargs
to ``False`` in :func:`extended_tapered_from_spec` reproduces the
Hopf-only behavior bit-exactly.

Honest scope (sprint M-vS Symmetries, 2026-06-07)
-------------------------------------------------
ℓ-parity is valid for the **composed** path (`build_composed_hamiltonian`
in `composed_qubit.py`).  In the **balanced** path
(`build_balanced_hamiltonian`) the cross-center V_ne adds off-diagonal
h1 elements with Δl ≠ 0 (multipole L=1 dipole term), which would
break per-block ℓ-parity.  The production `ecosystem_export.hamiltonian()`
uses the composed path, so production tapering is ℓ-parity-safe;
balanced FCI users should set ``use_ell_parity=False`` or rely on the
audit gate to drop the broken stabilizers.

Atom-swap and inversion stabilizers are valid in both composed and
balanced paths (the chemistry Hamiltonian is symmetric under
geometric equivalence of nuclei regardless of cross-center treatment).
"""

from __future__ import annotations

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    from openfermion import QubitOperator, jordan_wigner  # noqa: F401
    from openfermion.transforms import taper_off_qubits
    from openfermion.utils import commutator
    _HAS_OPENFERMION = True
except ImportError:  # pragma: no cover
    _HAS_OPENFERMION = False


# ---------------------------------------------------------------------------
# Helpers (mirror z2_tapering.py)
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
    return max_coef <= atol, max_coef


def _z_string_to_bitvec(op: "QubitOperator", n_qubits: int) -> np.ndarray:
    """Convert an all-Z Pauli string (single term) to a length-n_qubits
    binary vector with bit q = 1 if Z_q appears in the term."""
    vec = np.zeros(n_qubits, dtype=np.uint8)
    # Each term of op is a tuple of (qubit_index, Pauli) pairs;
    # for stabilizers we expect a single term.
    for term in op.terms:
        for (q, p) in term:
            if p != 'Z':
                raise ValueError(
                    f"_z_string_to_bitvec expects all-Z Pauli strings; "
                    f"got {p} on qubit {q}"
                )
            vec[q] ^= 1
    return vec


def _drop_linearly_dependent_z_stabilizers(
    stabilizers: List[Tuple[str, "QubitOperator"]],
    n_qubits: int,
) -> Tuple[List[Tuple[str, "QubitOperator"]], List[Tuple[str, int]]]:
    """Greedily keep stabilizers whose Z-bit vectors are GF(2)-independent.

    Stabilizers are kept in the order provided; the first occurrence
    in each linear-dependence orbit is retained. Returns (kept, dropped)
    where each entry is the original (kind, op) tuple.

    NOTE: alpha/beta parity stabilizers have Z on ALL alpha (even) or
    ALL beta (odd) qubits respectively; they should normally be
    independent of any sub-block Z-string. We preserve them first.
    """
    # Compute bit vectors for all
    vecs: List[Tuple[int, str, "QubitOperator", np.ndarray]] = []
    for i, (kind, op) in enumerate(stabilizers):
        try:
            v = _z_string_to_bitvec(op, n_qubits)
        except ValueError:
            v = None
        vecs.append((i, kind, op, v))

    # GF(2) row-reduction to identify independent set
    basis: List[np.ndarray] = []
    kept: List[Tuple[str, "QubitOperator"]] = []
    dropped: List[Tuple[str, int]] = []
    for (i, kind, op, v) in vecs:
        if v is None:
            # non-Z stabilizer; pass through unchanged (don't drop)
            kept.append((kind, op))
            continue
        # Reduce v against current basis
        v_red = v.copy()
        for b in basis:
            # pivot column = first 1 in b
            pivot_cols = np.where(b == 1)[0]
            if len(pivot_cols) == 0:
                continue
            pivot = pivot_cols[0]
            if v_red[pivot] == 1:
                v_red = np.bitwise_xor(v_red, b)
        if np.any(v_red):
            basis.append(v_red)
            kept.append((kind, op))
        else:
            dropped.append((kind, i))
    return kept, dropped


# ---------------------------------------------------------------------------
# Candidate 4: ℓ-parity Z₂ stabilizers
# ---------------------------------------------------------------------------

def build_ell_parity_stabilizers(
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    mode: str = 'per_block',
) -> List["QubitOperator"]:
    """Build P_ℓ Z₂ stabilizers tracking (−1)^(N_lodd).

    For each sub-block with at least one l-odd orbital (l ≥ 1, odd),
    emit ``P_ℓ^block = Π_{q: l(q) odd, sb(q) = block} Z_q`` over the
    α and β spin-orbitals of each l-odd spatial orbital (so each
    spatial l-odd orbital contributes 2 Z-factors).

    Parameters
    ----------
    orbital_table : sequence of ``(sub_block_key, n, l, m)`` tuples
        From `geovac.z2_tapering._enumerate_orbitals(spec)` or analog.
    mode : ``'global'`` or ``'per_block'``
        ``'global'`` returns one stabilizer over all l-odd orbitals
        (smaller ΔQ but more robust).  ``'per_block'`` returns one
        stabilizer per sub-block with l-odd content (largest ΔQ).

    Returns
    -------
    list of QubitOperator
        One per sub-block (per_block mode) or [P_global] (global mode).
        Empty list if no l-odd orbitals exist in the basis.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")
    if mode not in ('global', 'per_block'):
        raise ValueError(f"mode must be 'global' or 'per_block'")

    sb_to_lodd: "OrderedDict[Any, List[int]]" = OrderedDict()
    for i, (sb, n, l, m) in enumerate(orbital_table):
        sb_to_lodd.setdefault(sb, [])
        if l % 2 == 1:
            sb_to_lodd[sb].append(i)

    P_list: List["QubitOperator"] = []
    if mode == 'per_block':
        for sb, l_odd_orbitals in sb_to_lodd.items():
            if not l_odd_orbitals:
                continue
            spin_qubits: List[int] = []
            for p in l_odd_orbitals:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            P_list.append(_z_string_from_indices(spin_qubits))
    else:  # global
        all_l_odd: List[int] = []
        for orbs in sb_to_lodd.values():
            all_l_odd.extend(orbs)
        if all_l_odd:
            spin_qubits = []
            for p in all_l_odd:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            P_list.append(_z_string_from_indices(spin_qubits))

    return P_list


# ---------------------------------------------------------------------------
# Candidate 1: atom-swap Z₂
# ---------------------------------------------------------------------------

def find_equivalent_atom_pairs(
    spec: Any,
    nuclei: List[Dict[str, Any]],
    atol: float = 1e-6,
) -> List[Tuple[int, int]]:
    """Identify pairs of equivalent atoms (i, j) in the spec.

    Two atoms are 'equivalent' if:
      - Same Z (atomic number)
      - The sub-block structure attached to each is the same (same
        block_type, same Z_center, max_n, l_min, pk_A, pk_B)
      - Their geometric positions are related by a reflection or
        rotation that preserves the molecule (checked via the
        canonical-position helper below).

    Returns
    -------
    list of (i, j) tuples
        Pairs of nucleus indices (into `nuclei`) that are equivalent.
        Each pair is emitted once with i < j.  For 3-equivalent
        molecules like NH₃ we pick a single pair (any transposition);
        for 4-equivalent like CH₄ we pick the Klein V_4 generators
        ((12)(34) and (13)(24) on the H indices).
    """
    n_nuc = len(nuclei)
    # Group by Z
    z_to_indices: "OrderedDict[float, List[int]]" = OrderedDict()
    for i, nuc in enumerate(nuclei):
        z_to_indices.setdefault(float(nuc['Z']), []).append(i)

    pairs: List[Tuple[int, int]] = []
    for z, idx_list in z_to_indices.items():
        if len(idx_list) < 2:
            continue
        positions = [np.array(nuclei[i]['position'], dtype=float)
                     for i in idx_list]
        # Distance matrix to all OTHER nuclei
        # Two atoms are equivalent if the multiset of distances to
        # all other nuclei (excluding themselves) is identical.
        all_other_pos = [np.array(nuclei[k]['position'], dtype=float)
                         for k in range(n_nuc)]

        def distance_signature(idx: int) -> Tuple[float, ...]:
            p = np.array(nuclei[idx]['position'], dtype=float)
            dists = []
            for k in range(n_nuc):
                if k == idx:
                    continue
                q = np.array(nuclei[k]['position'], dtype=float)
                d = float(np.linalg.norm(p - q))
                dists.append(round(d, 8))
            return tuple(sorted(dists))

        sigs = [distance_signature(i) for i in idx_list]
        # Group indices by signature
        sig_groups: "OrderedDict[Tuple[float, ...], List[int]]" = OrderedDict()
        for ii, sig in zip(idx_list, sigs):
            sig_groups.setdefault(sig, []).append(ii)
        # For each group of ≥2 equivalents, emit pairs
        for sig, group in sig_groups.items():
            if len(group) < 2:
                continue
            if len(group) == 2:
                pairs.append((group[0], group[1]))
            elif len(group) == 3:
                # S_3: one Z₂ transposition (pick first pair)
                pairs.append((group[0], group[1]))
            elif len(group) == 4:
                # S_4: Klein V_4 generators
                # (12)(34) and (13)(24); emit as two paired-swap stabilizers
                # For simplicity, we emit two SINGLE pair swaps that
                # mutually commute: (12) is NOT compatible with (13) alone
                # under Klein V_4. Instead emit the double transpositions:
                # represent (12)(34) as a single stabilizer covering both
                # swaps simultaneously. See `build_atom_swap_rotation`.
                pairs.append((group[0], group[1]))
                pairs.append((group[2], group[3]))
                # Second Klein generator: (13)(24)
                pairs.append((group[0], group[2]))
                pairs.append((group[1], group[3]))
            else:
                # More than 4 equivalents: emit just one pair (lower bound)
                pairs.append((group[0], group[1]))

    return pairs


def _match_equivalent_blocks(
    spec: Any,
    nuc_pair: Tuple[int, int],
    block_nuc_map: Dict[int, Dict[str, int]],
) -> List[Tuple[int, int]]:
    """For an equivalent-nucleus pair (i, j), find block pairs (b_a, b_b)
    where block b_a's nuclei are swap-equivalents of block b_b's.

    Specifically, a block with (center=c, partner=p) matches another
    with (center=c', partner=p') if:
      - block parameters identical (same block_type, Z_center, max_n, ...)
      - (c, p) ↔ (c', p') under the swap (i, j) (applied as a transposition
        of nucleus indices)

    Returns list of (block_idx_a, block_idx_b) with a < b.
    """
    i, j = nuc_pair

    def apply_swap(k: int) -> int:
        if k == i:
            return j
        if k == j:
            return i
        return k

    used: set = set()
    matched: List[Tuple[int, int]] = []
    for a_idx, blk_a in enumerate(spec.blocks):
        if a_idx in used:
            continue
        c_a = block_nuc_map.get(a_idx, {}).get('center_idx', -1)
        p_a = block_nuc_map.get(a_idx, {}).get('partner_idx', -1)
        c_a_swap = apply_swap(c_a)
        p_a_swap = apply_swap(p_a)
        # Look for a block b with (center, partner) == swapped(a)
        for b_idx, blk_b in enumerate(spec.blocks):
            if b_idx <= a_idx or b_idx in used:
                continue
            c_b = block_nuc_map.get(b_idx, {}).get('center_idx', -1)
            p_b = block_nuc_map.get(b_idx, {}).get('partner_idx', -1)
            if c_b != c_a_swap or p_b != p_a_swap:
                continue
            # Verify parameter equality
            if blk_a.block_type != blk_b.block_type:
                continue
            if abs(blk_a.Z_center - blk_b.Z_center) > 1e-9:
                continue
            if blk_a.max_n != blk_b.max_n:
                continue
            if getattr(blk_a, 'l_min', 0) != getattr(blk_b, 'l_min', 0):
                continue
            if getattr(blk_a, 'has_h_partner', False) != getattr(blk_b, 'has_h_partner', False):
                continue
            # Match
            matched.append((a_idx, b_idx))
            used.add(a_idx)
            used.add(b_idx)
            break
    return matched


def build_atom_swap_rotation_and_stabilizers(
    spec: Any,
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    nuclei: List[Dict[str, Any]],
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """Build rotation to atom-swap (+/-) eigenbasis.

    For each equivalent-nucleus pair (i, j), find block pairs that map
    to each other under the swap, then swap BOTH center and partner
    sub-blocks of those blocks (this handles e.g. BeH₂ where bond_1
    and bond_2 are swap-paired, requiring swap of bond_1_center ↔
    bond_2_center AND bond_1_partner ↔ bond_2_partner).

    Returns
    -------
    U_swap : (M, M) ndarray
        Orthogonal rotation to the atom-swap eigenbasis.
    parity_swap : (M,) int ndarray
        +1 for symmetric, −1 for antisymmetric.
    diagnostics : list of dict
    """
    M = len(orbital_table)
    block_nuc_map = _get_block_nucleus_map(spec, nuclei)
    pairs = find_equivalent_atom_pairs(spec, nuclei)

    U = np.eye(M)
    parity = np.ones(M, dtype=int)
    diagnostics: List[Dict[str, Any]] = []

    used_orbital_indices: set = set()
    for (i, j) in pairs:
        block_pairs = _match_equivalent_blocks(spec, (i, j), block_nuc_map)
        if not block_pairs:
            continue

        pair_anti_orbs: List[int] = []
        sb_pair_records: List[Tuple[Any, Any]] = []
        for (a_idx, b_idx) in block_pairs:
            blk_a = spec.blocks[a_idx]
            blk_b = spec.blocks[b_idx]
            # Swap both center and partner sides
            for side in ('center', 'partner'):
                if side == 'partner' and not getattr(blk_a, 'has_h_partner', False):
                    continue
                label_a = (a_idx, blk_a.label, side)
                label_b = (b_idx, blk_b.label, side)
                sb_pair_records.append((label_a, label_b))

                idx_in_a: Dict[Tuple[int, int, int], int] = {}
                idx_in_b: Dict[Tuple[int, int, int], int] = {}
                for orb_idx, (sb, n, l, m) in enumerate(orbital_table):
                    if sb == label_a:
                        idx_in_a[(n, l, m)] = orb_idx
                    elif sb == label_b:
                        idx_in_b[(n, l, m)] = orb_idx

                for nlm, a_glob in idx_in_a.items():
                    b_glob = idx_in_b.get(nlm)
                    if b_glob is None:
                        continue
                    if a_glob in used_orbital_indices or b_glob in used_orbital_indices:
                        continue
                    inv_sqrt2 = 1.0 / np.sqrt(2.0)
                    a_low, b_high = min(a_glob, b_glob), max(a_glob, b_glob)
                    U[a_low, :] = 0.0
                    U[b_high, :] = 0.0
                    U[a_low, a_low] = inv_sqrt2
                    U[a_low, b_high] = inv_sqrt2
                    U[b_high, a_low] = inv_sqrt2
                    U[b_high, b_high] = -inv_sqrt2
                    parity[a_low] = +1
                    parity[b_high] = -1
                    used_orbital_indices.add(a_low)
                    used_orbital_indices.add(b_high)
                    pair_anti_orbs.append(b_high)

        if pair_anti_orbs:
            diagnostics.append({
                'pair': (i, j),
                'block_pairs': block_pairs,
                'sub_block_pairs': sb_pair_records,
                'n_orbital_pairs': len(pair_anti_orbs),
            })

    err = float(np.max(np.abs(U @ U.T - np.eye(M))))
    if err > 1e-12:
        raise RuntimeError(f"swap rotation not orthogonal: err={err:.2e}")

    return U, parity, diagnostics


def _auto_assign_nucleus_indices(
    spec: Any, nuclei: List[Dict[str, Any]],
) -> Dict[int, Dict[str, int]]:
    """Auto-assign block → nucleus indices using the hydride convention.

    For hydride_spec-generated molecules (LiH, BeH₂, H₂O, NH₃, CH₄, etc.),
    the spec.blocks don't have center_nucleus_idx/partner_nucleus_idx
    populated (they default to -1).  Convention:
      - Heavy atom at nucleus index 0
      - H atoms at nucleus indices 1, 2, ..., n_H
      - Core, lone_pair, and bond_center blocks: center on heavy (0)
      - bond_partner blocks: partner on H_i (incrementing per bond)

    Returns
    -------
    dict mapping block_idx → {'center_idx': int, 'partner_idx': int}
        Uses explicit values from spec.blocks where set;
        falls back to the hydride convention otherwise.
    """
    out: Dict[int, Dict[str, int]] = {}
    if not nuclei:
        return out

    # Identify heavy nucleus (first with Z > 1.5)
    heavy_idx = next(
        (i for i, n in enumerate(nuclei) if float(n['Z']) > 1.5),
        0,
    )
    # H nuclei in order
    h_nuc_idx_list = [
        i for i, n in enumerate(nuclei) if float(n['Z']) <= 1.5
    ]

    h_counter = 0
    for b_idx, blk in enumerate(spec.blocks):
        # Use explicit values if set
        c_idx = getattr(blk, 'center_nucleus_idx', -1)
        p_idx = getattr(blk, 'partner_nucleus_idx', -1)
        if c_idx < 0:
            # Apply convention based on block_type
            c_idx = heavy_idx
        if p_idx < 0 and getattr(blk, 'has_h_partner', False):
            if h_counter < len(h_nuc_idx_list):
                p_idx = h_nuc_idx_list[h_counter]
                h_counter += 1
        out[b_idx] = {'center_idx': c_idx, 'partner_idx': p_idx}
    return out


def _get_block_nucleus_map(
    spec: Any, nuclei: List[Dict[str, Any]],
) -> Dict[int, Dict[str, int]]:
    """Resolve block → {center_idx, partner_idx} mapping using explicit
    values where set, else hydride-convention auto-assignment."""
    return _auto_assign_nucleus_indices(spec, nuclei)


def _sub_blocks_on_nucleus(
    spec: Any, nuc_idx: int,
    block_nuc_map: Optional[Dict[int, Dict[str, int]]] = None,
) -> List[Tuple[int, str, str]]:
    """Return list of (block_idx, label, side) for sub-blocks attached to
    nucleus `nuc_idx`.

    Uses `block_nuc_map` if provided (from `_get_block_nucleus_map`);
    otherwise falls back to spec.blocks[i] attribute values.
    """
    out: List[Tuple[int, str, str]] = []
    for b_idx, blk in enumerate(spec.blocks):
        if block_nuc_map is not None:
            c_idx = block_nuc_map.get(b_idx, {}).get('center_idx', -1)
            p_idx = block_nuc_map.get(b_idx, {}).get('partner_idx', -1)
        else:
            c_idx = getattr(blk, 'center_nucleus_idx', -1)
            p_idx = getattr(blk, 'partner_nucleus_idx', -1)
        if c_idx == nuc_idx:
            out.append((b_idx, blk.label, 'center'))
        if p_idx == nuc_idx:
            out.append((b_idx, blk.label, 'partner'))
    return out


def _match_sub_blocks_for_swap(
    spec: Any,
    sb_on_i: List[Tuple[int, str, str]],
    sb_on_j: List[Tuple[int, str, str]],
) -> List[Tuple[Any, Any]]:
    """Match sub-blocks on nucleus i with their swap-partners on j.

    Returns list of (sb_key_i, sb_key_j) pairs where each sb_key is
    the tuple ``(block_idx, label, side)`` consumed by the orbital
    enumerator.
    """
    matched: List[Tuple[Any, Any]] = []
    used_j = [False] * len(sb_on_j)
    for (bi, label_i, side_i) in sb_on_i:
        blk_i = spec.blocks[bi]
        for k, (bj, label_j, side_j) in enumerate(sb_on_j):
            if used_j[k]:
                continue
            blk_j = spec.blocks[bj]
            if side_i != side_j:
                continue
            if blk_i.block_type != blk_j.block_type:
                continue
            if abs(blk_i.Z_center - blk_j.Z_center) > 1e-9:
                continue
            if blk_i.max_n != blk_j.max_n:
                continue
            if getattr(blk_i, 'l_min', 0) != getattr(blk_j, 'l_min', 0):
                continue
            # Match!  Construct sb_key tuples used by orbital enumerator.
            sb_key_i = (bi, label_i, side_i)
            sb_key_j = (bj, label_j, side_j)
            matched.append((sb_key_i, sb_key_j))
            used_j[k] = True
            break
    return matched


# ---------------------------------------------------------------------------
# Candidate 2: spatial inversion Z₂ (centrosymmetric molecules)
# ---------------------------------------------------------------------------

def is_centrosymmetric(nuclei: List[Dict[str, Any]], atol: float = 1e-6) -> bool:
    """Check whether the nuclear arrangement has a center of inversion.

    Centrosymmetric means: for every nucleus at position r, there is
    another nucleus of the same Z at position -r (relative to the
    center of mass / inversion center, which we take as the centroid).
    """
    if not nuclei:
        return False
    # Compute centroid weighted by Z (to find the inversion center)
    total_z = sum(float(n['Z']) for n in nuclei)
    if total_z == 0:
        return False
    centroid = np.zeros(3)
    for n in nuclei:
        centroid += float(n['Z']) * np.array(n['position'], dtype=float)
    centroid /= total_z

    used = [False] * len(nuclei)
    for i, n in enumerate(nuclei):
        if used[i]:
            continue
        p = np.array(n['position'], dtype=float) - centroid
        # Look for the inversion partner at -p relative to centroid
        target = -p + centroid
        found = False
        for j, m in enumerate(nuclei):
            if used[j]:
                continue
            if abs(float(m['Z']) - float(n['Z'])) > 1e-9:
                continue
            q = np.array(m['position'], dtype=float)
            if np.linalg.norm(q - target) < atol:
                used[i] = True
                used[j] = True
                found = True
                break
        if not found:
            return False
    return True


def build_inversion_rotation_and_stabilizers(
    spec: Any,
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    nuclei: List[Dict[str, Any]],
    U_swap: np.ndarray,
    parity_swap: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, List["QubitOperator"]]:
    """For centrosymmetric molecules, build the inversion Z₂ stabilizer
    above the atom-swap rotation.

    Inversion i: r → −r maps orbital (n, l, m) → (−1)^l × (same orbital
    on the inverted atom).  After applying U_swap (which already maps
    swap-equivalent orbitals to ±-combinations), inversion further
    multiplies by (−1)^l per orbital.

    The combined parity is ``parity_inversion[i] = parity_swap[i] *
    (−1)^l(i)``.  Only orbitals with parity_inversion = -1 (in the
    POST-swap basis) carry the inversion Z₂ stabilizer's Z-factor.

    Note: this is the independent generator IN ADDITION TO swap;
    if it gives no new content (all l-even orbitals OR no swap
    structure), the stabilizer is trivial and dropped.
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")
    if not is_centrosymmetric(nuclei):
        return U_swap, parity_swap, []

    M = len(orbital_table)
    # In the U_swap basis, each orbital is either (a_low, +) or (b_high, -)
    # in a pair, OR unpaired (parity_swap = +1).  Inversion acts as
    # additional sign (-1)^l on each orbital.
    parity_inv = np.ones(M, dtype=int)
    for i, (sb, n, l, m) in enumerate(orbital_table):
        sign = (-1) ** l
        parity_inv[i] = parity_swap[i] * sign

    # Build the inversion stabilizer Z-string on orbitals with
    # parity_inv = -1.
    anti_orbs = [i for i in range(M) if parity_inv[i] == -1]
    if not anti_orbs:
        return U_swap, parity_swap, []

    # Check independence vs swap: if parity_inv == parity_swap (or
    # == -parity_swap up to sign convention), inversion is not an
    # independent generator.  Specifically, if (−1)^l is constant
    # (= +1 for all orbitals, which happens when all orbitals are
    # l-even), inversion equals swap.
    all_l_even = all(orb[2] % 2 == 0 for orb in orbital_table)
    if all_l_even:
        return U_swap, parity_swap, []  # not independent

    spin_qubits: List[int] = []
    for p in anti_orbs:
        spin_qubits.append(2 * p)
        spin_qubits.append(2 * p + 1)
    P_inv = _z_string_from_indices(spin_qubits)

    return U_swap, parity_swap, [P_inv]


# ---------------------------------------------------------------------------
# Combined extended tapering
# ---------------------------------------------------------------------------

def _build_joint_sb_map(
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    swap_diag: List[Dict[str, Any]],
) -> Dict[Any, Any]:
    """Build a map: original_sub_block_key → joint_sub_block_key, where
    swap-paired sub-blocks map to a common joint key.  Un-swapped
    sub-blocks map to themselves.

    Used to merge per-sub-block Hopf and ℓ-parity stabilizers across
    swap pairs (otherwise the per-sub-block stabilizer doesn't commute
    with the swap-rotated Hamiltonian, since the sub-block identity
    has been mixed by the swap).
    """
    sb_map: Dict[Any, Any] = {}
    for diag in swap_diag:
        for (label_a, label_b) in diag.get('sub_block_pairs', []):
            joint = ('SWAP_JOINT', label_a, label_b)
            sb_map[label_a] = joint
            sb_map[label_b] = joint
    return sb_map


def apply_extended_tapering(
    qubit_op_rotated: "QubitOperator",
    parity_hopf: np.ndarray,
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    parity_swap: Optional[np.ndarray] = None,
    swap_diag: Optional[List[Dict[str, Any]]] = None,
    use_hopf_per_block: bool = True,
    use_ell_parity: bool = True,
    use_atom_swap: bool = True,
    use_inversion: bool = True,
    inversion_stabilizers: Optional[List["QubitOperator"]] = None,
    sector_signs: Optional[Sequence[int]] = None,
    drop_noncommuting: bool = True,
    audit_tol: float = 1e-10,
    verbose: bool = False,
) -> Tuple["QubitOperator", Dict[str, Any]]:
    """Apply combined Hopf + ℓ-parity + atom-swap + inversion tapering.

    `qubit_op_rotated` is the JW-encoded Hamiltonian AFTER all the
    rotations (Hopf + swap) have been applied to (h1, eri).  This
    function builds the stabilizers (Hopf m_l, ℓ-parity, swap,
    inversion), audits commutation with the Hamiltonian, drops any
    that fail, and applies `taper_off_qubits` to the survivors.

    Returns
    -------
    qubit_op_tapered : QubitOperator
    meta : dict with keys:
        - n_stabs_total : int
        - n_stabs_kept : int
        - dropped : list of (kind, index, residual)
        - kinds_kept : list of strings
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")

    M = len(orbital_table)
    Z_alpha = _z_string_from_indices(2 * p for p in range(M))
    Z_beta = _z_string_from_indices(2 * p + 1 for p in range(M))

    stabilizers: List[Tuple[str, "QubitOperator"]] = [
        ('alpha', Z_alpha), ('beta', Z_beta),
    ]

    # Build joint sub-block map: when swap rotation merges sub-blocks
    # (sb_A, sb_B) into a joint pair, the per-sub-block Hopf and
    # ℓ-parity stabilizers must be MERGED across the joint pair
    # (otherwise they don't commute with the swap-rotated Hamiltonian).
    joint_sb_map: Dict[Any, Any] = {}
    if swap_diag:
        joint_sb_map = _build_joint_sb_map(orbital_table, swap_diag)

    def _joint_key(sb: Any) -> Any:
        return joint_sb_map.get(sb, sb)

    # Hopf Z₂: per joint_block stabilizers from parity_hopf
    if use_hopf_per_block:
        from collections import OrderedDict
        sb_to_anti: "OrderedDict[Any, List[int]]" = OrderedDict()
        for p, (sb, n, l, m) in enumerate(orbital_table):
            jsb = _joint_key(sb)
            sb_to_anti.setdefault(jsb, [])
            if parity_hopf[p] == -1:
                sb_to_anti[jsb].append(p)
        for jsb, anti_orbitals in sb_to_anti.items():
            if not anti_orbitals:
                continue
            spin_qubits = []
            for p in anti_orbitals:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(
                (f'hopf:{jsb}', _z_string_from_indices(spin_qubits))
            )

    # ℓ-parity Z₂: per joint_block
    if use_ell_parity:
        from collections import OrderedDict
        sb_to_lodd: "OrderedDict[Any, List[int]]" = OrderedDict()
        for p, (sb, n, l, m) in enumerate(orbital_table):
            jsb = _joint_key(sb)
            sb_to_lodd.setdefault(jsb, [])
            if l % 2 == 1:
                sb_to_lodd[jsb].append(p)
        for jsb, l_odd_orbitals in sb_to_lodd.items():
            if not l_odd_orbitals:
                continue
            spin_qubits = []
            for p in l_odd_orbitals:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(
                (f'ell_parity:{jsb}', _z_string_from_indices(spin_qubits))
            )

    # Atom-swap Z₂: stabilizers come from parity_swap (parity = -1
    # on the antisymmetric orbital of each swap pair).  We don't
    # build this here; instead, the caller passes parity_swap and
    # we read off the Z-string.  This avoids re-doing the rotation
    # logic.  If parity_swap is None, no swap stabilizer.
    if use_atom_swap and parity_swap is not None:
        # Group anti-orbitals by their pair structure: each pair_anti
        # was emitted with its own stabilizer in the rotation builder.
        # Here we re-emit a single combined stabilizer per pair.
        # Since the caller wrote parity_swap[i] = -1 for the higher
        # orbital index of each pair, the Z-string is over all such.
        anti = [i for i in range(M) if parity_swap[i] == -1]
        if anti:
            spin_qubits = []
            for p in anti:
                spin_qubits.append(2 * p)
                spin_qubits.append(2 * p + 1)
            stabilizers.append(('atom_swap', _z_string_from_indices(spin_qubits)))

    # Inversion Z₂: passed in directly
    if use_inversion and inversion_stabilizers:
        for k, P in enumerate(inversion_stabilizers):
            stabilizers.append((f'inversion:{k}', P))

    # Audit each stabilizer (skip alpha/beta which are particle-number).
    # Three filtering passes:
    #   1. Commutation audit (drop if [H, P] != 0 at audit_tol)
    #   2. Linear-independence audit over GF(2) on the Z-bit vectors
    #   3. openfermion's own check at taper_off_qubits time (last resort)
    n_qubits = 2 * len(orbital_table)
    kept_with_kind: List[Tuple[str, "QubitOperator"]] = []
    dropped: List[Tuple[str, int, float]] = []
    for k_idx, (kind, P_i) in enumerate(stabilizers):
        if kind in ('alpha', 'beta'):
            kept_with_kind.append((kind, P_i))
            continue
        ok, residual = _commutes(qubit_op_rotated, P_i, atol=audit_tol)
        if ok or not drop_noncommuting:
            kept_with_kind.append((kind, P_i))
        else:
            dropped.append((kind, k_idx, residual))
            if verbose:
                print(f"  Dropping stabilizer {kind!r} (noncommuting): "
                      f"residual {residual:.2e}")

    # Linear-independence pass
    kept_with_kind, dep_dropped = _drop_linearly_dependent_z_stabilizers(
        kept_with_kind, n_qubits=n_qubits,
    )
    for (kind, k_idx) in dep_dropped:
        dropped.append((kind, k_idx, -1.0))  # -1 marker = dependence drop
        if verbose:
            print(f"  Dropping stabilizer {kind!r} (linearly dependent on prior)")

    kinds_kept = [k for (k, _) in kept_with_kind]
    kept = [op for (_, op) in kept_with_kind]

    if sector_signs is None:
        signs = [1] * len(kept)
    else:
        signs = list(sector_signs)
        if len(signs) != len(kept):
            raise ValueError(
                f"sector_signs length {len(signs)} != n_kept {len(kept)}"
            )

    signed_stabs = [s * sgn for s, sgn in zip(kept, signs)]
    qubit_op_tapered = taper_off_qubits(qubit_op_rotated, signed_stabs)

    return qubit_op_tapered, {
        'n_stabs_total': len(stabilizers),
        'n_stabs_kept': len(kept),
        'kinds_kept': kinds_kept,
        'dropped': dropped,
    }


# ---------------------------------------------------------------------------
# End-to-end pipeline from a MolecularSpec
# ---------------------------------------------------------------------------

def extended_tapered_from_spec(
    spec: Any,
    use_hopf: bool = True,
    use_ell_parity: bool = True,
    use_atom_swap: bool = False,
    use_inversion: bool = False,
    pk_in_hamiltonian: bool = True,
    builder: str = 'composed',
    R: Optional[float] = None,
    nuclei: Optional[List[Dict[str, Any]]] = None,
    sector_signs: Optional[Sequence[int]] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """End-to-end extended-tapering pipeline.

    Parameters
    ----------
    spec : MolecularSpec
    use_hopf : bool, default True
        Apply Hopf m_l → −m_l reflection Z₂ tapering (v3.52.0).
    use_ell_parity : bool, default True
        Apply ℓ-parity Z₂ tapering.  Strict improvement vs Hopf-only
        in the composed path: saves ``+n_sub_blocks`` qubits AND
        reduces Pauli count (~3–15% in the verification panel).
    use_atom_swap : bool, default False
        Apply equivalent-atom permutation Z₂ tapering.  Saves
        additional qubits on polyatomics (BeH₂ +1, H₂O +1, NH₃ +1,
        CH₄ +2) BUT typically inflates the Pauli count 2–4× because
        the rotation mixes orbitals across sub-blocks.  Off by
        default; opt in if minimizing qubit count is more important
        than minimizing Pauli/measurement count.
    use_inversion : bool, default False
        Apply spatial inversion Z₂ on top of atom-swap.  Often
        redundant with atom-swap + ℓ-parity (the inversion stabilizer
        becomes linearly dependent after these).  Off by default.
    pk_in_hamiltonian : bool
        Forwarded to the builder.
    builder : ``'composed'`` or ``'balanced'``
        Which Hamiltonian to taper.  Composed is the production
        default for `ecosystem_export.hamiltonian()`.  Balanced is
        used for FCI benchmarking; ℓ-parity is typically dropped by
        the audit gate in balanced (cross-center V_ne breaks it).
    R, nuclei : optional
        Forwarded to the builder for balanced.

    Returns
    -------
    dict with:
        - qubit_op_naive, qubit_op_tapered : QubitOperator
        - Q_naive, Q_tapered, delta_Q : int
        - n_stabs_kept, dropped : tapering audit info
        - U_total, parity_hopf, parity_swap : rotation diagnostics
        - kinds_kept : list of stabilizer-kind strings retained
    """
    if not _HAS_OPENFERMION:
        raise ImportError("openfermion is required")

    from geovac.z2_tapering import (
        _enumerate_orbitals, build_pm_rotation, rotate_h1_eri,
    )
    from geovac.qubit_encoding import build_fermion_op_from_integrals

    orbital_table = _enumerate_orbitals(spec)
    M = len(orbital_table)

    # Build the un-rotated Hamiltonian
    if builder == 'composed':
        from geovac.composed_qubit import build_composed_hamiltonian
        result = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
        )
        nuclei_list = nuclei or []
    elif builder == 'balanced':
        from geovac.balanced_coupled import build_balanced_hamiltonian
        result = build_balanced_hamiltonian(
            spec, R=R, nuclei=nuclei, verbose=verbose,
        )
        nuclei_list = nuclei or []
    else:
        raise ValueError(f"builder must be 'composed' or 'balanced'")

    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']

    # Compute Q_naive from naive JW
    fop_naive = build_fermion_op_from_integrals(h1, eri, nuc_rep)
    qop_naive = jordan_wigner(fop_naive)
    Q_naive = 2 * M

    # Build rotations
    U_hopf = np.eye(M)
    parity_hopf = np.ones(M, dtype=int)
    if use_hopf:
        U_hopf, parity_hopf = build_pm_rotation(orbital_table)

    U_swap = np.eye(M)
    parity_swap = np.ones(M, dtype=int)
    swap_diag: List[Dict[str, Any]] = []
    inv_stabs: List["QubitOperator"] = []
    if use_atom_swap and nuclei_list:
        U_swap, parity_swap, swap_diag = (
            build_atom_swap_rotation_and_stabilizers(
                spec, orbital_table, nuclei_list,
            )
        )
        if use_inversion:
            _, _, inv_stabs = build_inversion_rotation_and_stabilizers(
                spec, orbital_table, nuclei_list, U_swap, parity_swap,
            )

    # Hopf and swap rotations commute (Hopf: within-block m mixing;
    # swap: cross-block label mixing with same (n,l,m)).  Apply both.
    U_total = U_swap @ U_hopf

    # Rotate h1, eri
    h1_rot, eri_rot = rotate_h1_eri(h1, eri, U_total)

    # JW the rotated Hamiltonian
    fop_rot = build_fermion_op_from_integrals(h1_rot, eri_rot, nuc_rep)
    qop_rot = jordan_wigner(fop_rot)

    # Apply combined tapering
    qop_tap, audit = apply_extended_tapering(
        qop_rot, parity_hopf=parity_hopf,
        orbital_table=orbital_table,
        parity_swap=parity_swap,
        swap_diag=swap_diag,
        use_hopf_per_block=use_hopf,
        use_ell_parity=use_ell_parity,
        use_atom_swap=use_atom_swap,
        use_inversion=use_inversion,
        inversion_stabilizers=inv_stabs,
        sector_signs=sector_signs,
        verbose=verbose,
    )

    # Compute Q_tapered
    Q_tapered = 0
    for term in qop_tap.terms:
        for q, _ in term:
            if q + 1 > Q_tapered:
                Q_tapered = q + 1

    return {
        'qubit_op_naive': qop_naive,
        'qubit_op_tapered': qop_tap,
        'Q_naive': Q_naive,
        'Q_tapered': Q_tapered,
        'delta_Q': Q_naive - Q_tapered,
        'n_stabs_kept': audit['n_stabs_kept'],
        'kinds_kept': audit['kinds_kept'],
        'dropped': audit['dropped'],
        'orbital_table': orbital_table,
        'parity_hopf': parity_hopf,
        'parity_swap': parity_swap,
        'swap_diag': swap_diag,
        'U_total': U_total,
        'builder': builder,
    }
