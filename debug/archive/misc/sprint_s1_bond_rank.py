"""Sprint S1 — Direct MPO bond-rank measurement on GeoVac Hamiltonians.

The operator Schmidt rank chi_k of an operator H at cut k is the rank of
the matrix M[(L, R)] = coefficient_of_pauli_string_with_left_part_L_right_part_R.

This is the exact MPO bond dimension of the operator after maximally
compact decomposition (Schollwock 2011 review §4.3). For chemistry
Hamiltonians it can be computed directly from the Pauli decomposition
without running DMRG.

Theorem 3.2.A prediction (Connes-vS / MPO scoping memo):
  GeoVac MPOs should show chi_k drop dramatically at every sub-block
  boundary because cross-block ERIs vanish bit-exactly.

Sprint S1 falsifier: at the cuts corresponding to sub-block boundaries in
the qubit ordering, chi_k should be substantially smaller than at cuts
inside sub-blocks. If observed, Theorem 3.2.A is empirically supported.
If not observed, Theorem 3.2.A is falsified at this level.

Test panel: LiH (3 sub-blocks), BeH2 (5 sub-blocks), H2O (7 sub-blocks).
Plus a single-sub-block atomic He as control (no boundary, no drop).
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Tuple

import numpy as np

from openfermion import QubitOperator

from geovac.composed_qubit import build_composed_hamiltonian, _enumerate_states
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec


def operator_schmidt_rank(qop: QubitOperator, cut: int,
                          n_qubits: int,
                          rel_thr: float = 1e-12) -> Tuple[int, np.ndarray]:
    """Compute the operator Schmidt rank of `qop` at qubit cut `cut`.

    qop = sum_{P_left, P_right} M[P_left, P_right] (P_left otimes P_right)
    where P_left acts on qubits [0, cut) and P_right acts on [cut, n_qubits).
    chi_k = rank(M).

    Returns (rank, singular_values).
    """
    if cut <= 0 or cut >= n_qubits:
        return 1, np.array([1.0])

    coef = {}
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
        coef[(left, right)] = c

    left_keys = sorted(set(k[0] for k in coef))
    right_keys = sorted(set(k[1] for k in coef))
    if not left_keys or not right_keys:
        return 1, np.array([1.0])

    left_idx = {k: i for i, k in enumerate(left_keys)}
    right_idx = {k: i for i, k in enumerate(right_keys)}

    M = np.zeros((len(left_keys), len(right_keys)))
    for (l, r), v in coef.items():
        M[left_idx[l], right_idx[r]] = float(v.real if hasattr(v, 'real') else v)

    sv = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))
    return rank, sv


def get_subblock_boundaries(spec) -> List[int]:
    """Return the list of qubit boundaries between sub-blocks.

    Each sub-block has `n_orbitals` spatial orbitals -> 2 * n_orbitals qubits
    under JW with alpha/beta interleaving. The boundary AFTER sub-block i
    is at qubit position 2 * (cumulative_orbitals_up_to_and_including_i).
    """
    boundaries = []
    cumulative = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        center_states = _enumerate_states(blk.max_n, l_min=l_min)
        cumulative += len(center_states)
        boundaries.append(2 * cumulative)
        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_max_n)
            cumulative += len(partner_states)
            boundaries.append(2 * cumulative)
    # Last boundary is at end; the meaningful interior boundaries are
    # all but the final one.
    return boundaries[:-1]


def analyze_molecule(spec, label: str) -> Dict:
    print(f"\n{'='*70}\n{label}\n{'='*70}")
    res = build_composed_hamiltonian(spec, verbose=False)
    qop = res['qubit_op']
    M = res['h1'].shape[0]
    Q = 2 * M

    # Sub-block boundaries (qubit positions)
    boundaries = get_subblock_boundaries(spec)
    print(f"  M = {M}, Q = {Q}, sub-block boundaries: {boundaries}")
    print(f"  Total Pauli terms: {len(qop.terms)}")

    # Compute bond rank at every cut
    profile = []
    for cut in range(1, Q):
        rank, sv = operator_schmidt_rank(qop, cut, Q)
        is_boundary = cut in boundaries
        profile.append({
            'cut': cut,
            'bond_rank': rank,
            'is_subblock_boundary': is_boundary,
            'top_5_sv': [float(x) for x in sv[:5]],
        })

    # Sub-block boundary ranks
    boundary_ranks = [
        (cut, rank)
        for cut, rank in [(p['cut'], p['bond_rank']) for p in profile]
        if cut in boundaries
    ]
    # Interior ranks
    interior_ranks = [
        (cut, rank)
        for cut, rank in [(p['cut'], p['bond_rank']) for p in profile]
        if cut not in boundaries
    ]

    # Print profile
    print(f"\n  Bond-rank profile (chi_k at every cut):")
    print(f"    {'cut':>4s}  {'rank':>6s}  {'boundary?':>10s}")
    for p in profile:
        marker = "  <-- BLOCK BOUNDARY" if p['is_subblock_boundary'] else ""
        print(f"    {p['cut']:>4d}  {p['bond_rank']:>6d}  "
              f"{'YES' if p['is_subblock_boundary'] else '   ':>10s}  {marker}")

    boundary_max = max((r for _, r in boundary_ranks), default=0)
    boundary_min = min((r for _, r in boundary_ranks), default=0)
    interior_max = max((r for _, r in interior_ranks), default=0)
    interior_min = min((r for _, r in interior_ranks), default=0)
    interior_mean = np.mean([r for _, r in interior_ranks]) if interior_ranks else 0

    print(f"\n  SUMMARY:")
    print(f"    Bond rank at sub-block boundaries:")
    print(f"      min: {boundary_min}, max: {boundary_max}")
    print(f"    Bond rank at interior cuts:")
    print(f"      min: {interior_min}, max: {interior_max}, mean: {interior_mean:.1f}")
    if boundary_ranks and interior_ranks:
        ratio = interior_max / max(boundary_max, 1)
        print(f"    Ratio (interior_max / boundary_max): {ratio:.2f}x")

    return {
        'molecule': label,
        'M': M, 'Q': Q,
        'n_pauli': len(qop.terms),
        'subblock_boundaries': boundaries,
        'profile': profile,
        'boundary_ranks': boundary_ranks,
        'interior_ranks': interior_ranks,
        'boundary_max_rank': boundary_max,
        'boundary_min_rank': boundary_min,
        'interior_max_rank': interior_max,
        'interior_min_rank': interior_min,
        'interior_mean_rank': float(interior_mean),
    }


def main():
    results = {}
    for spec_fn, label in [
        (lih_spec, 'LiH'),
        (beh2_spec, 'BeH2'),
        (h2o_spec, 'H2O'),
    ]:
        try:
            results[label] = analyze_molecule(spec_fn(), label)
        except Exception as e:
            import traceback
            traceback.print_exc()
            results[label] = {'error': str(e)}

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s1_bond_rank.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")

    # Aggregate verdict
    print(f"\n{'='*70}\nVERDICT\n{'='*70}")
    for label, res in results.items():
        if 'error' in res:
            print(f"  {label}: ERROR")
            continue
        ratio = res['interior_max_rank'] / max(res['boundary_max_rank'], 1)
        verdict = (
            "STRONG drop at boundaries"
            if ratio >= 3.0 else
            "MODERATE drop at boundaries"
            if ratio >= 1.5 else
            "NO drop at boundaries"
        )
        print(f"  {label:6s}: interior_max = {res['interior_max_rank']}, "
              f"boundary_max = {res['boundary_max_rank']}, "
              f"ratio = {ratio:.2f}x   -- {verdict}")


if __name__ == '__main__':
    main()
