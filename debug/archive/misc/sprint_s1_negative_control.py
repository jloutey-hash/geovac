"""Sprint S1 negative-control: bond rank on Paper 19 balanced-coupled builder.

The standard composed builder has cross-block ERIs = 0 bit-exactly (F4 from
the DF memo, this morning). Sprint S1 found chi_k = 2 at every sub-block
boundary as a consequence.

The Paper 19 balanced-coupled builder has NON-ZERO cross-center V_ne and
cross-block ERIs. By the structural argument in Sprint S1 memo §S1-4, the
chi_k = 2 floor should BREAK at would-be sub-block boundaries when cross-block
coupling is present.

This is the negative control. If chi_k > 2 at would-be boundaries in balanced-
coupled, the structural reading is confirmed.
"""

from __future__ import annotations

import json
import os

import numpy as np

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sprint_s1_bond_rank import operator_schmidt_rank, get_subblock_boundaries


def main():
    print("="*70)
    print("Negative control: LiH balanced-coupled bond rank")
    print("="*70)

    spec = lih_spec()
    boundaries = get_subblock_boundaries(spec)
    print(f"  Sub-block boundaries (would-be in qubit order): {boundaries}")

    # Build balanced-coupled LiH (cross-block V_ne ON)
    print("  Building balanced-coupled Hamiltonian...")
    res = build_balanced_hamiltonian(spec, R=3.015, verbose=False)
    qop = res['qubit_op']
    M = res['M'] if 'M' in res else int(np.sqrt(res['eri'].size ** (1/2)))
    Q = res['Q'] if 'Q' in res else 2 * M
    print(f"  M = {M}, Q = {Q}")
    print(f"  Total Pauli terms: {len(qop.terms)}")

    # Compute bond rank at every cut
    profile = []
    for cut in range(1, Q):
        rank, _sv = operator_schmidt_rank(qop, cut, Q)
        is_boundary = cut in boundaries
        profile.append({
            'cut': cut,
            'bond_rank': rank,
            'is_subblock_boundary': is_boundary,
        })

    boundary_ranks = [p['bond_rank'] for p in profile if p['is_subblock_boundary']]
    interior_ranks = [p['bond_rank'] for p in profile if not p['is_subblock_boundary']]

    print(f"\n  Bond-rank profile:")
    print(f"    {'cut':>4s}  {'rank':>6s}  {'boundary?':>10s}")
    for p in profile:
        marker = "  <-- BLOCK BOUNDARY" if p['is_subblock_boundary'] else ""
        print(f"    {p['cut']:>4d}  {p['bond_rank']:>6d}  "
              f"{'YES' if p['is_subblock_boundary'] else '   ':>10s}  {marker}")

    print(f"\n  At sub-block boundaries:")
    print(f"    min: {min(boundary_ranks)}, max: {max(boundary_ranks)}")
    print(f"  At interior cuts:")
    print(f"    min: {min(interior_ranks)}, max: {max(interior_ranks)}, "
          f"mean: {np.mean(interior_ranks):.1f}")

    composed_boundary = 2  # from Sprint S1 main result
    print(f"\n  COMPARISON:")
    print(f"    Composed builder boundary rank (Sprint S1):       {composed_boundary}")
    print(f"    Balanced-coupled builder boundary rank (here):    "
          f"min={min(boundary_ranks)}, max={max(boundary_ranks)}")
    if max(boundary_ranks) > composed_boundary:
        print(f"    -> NEGATIVE CONTROL CONFIRMED: cross-block coupling lifts chi_k")
        print(f"       from {composed_boundary} (composed) to {max(boundary_ranks)} (balanced)")
    elif max(boundary_ranks) == composed_boundary:
        print(f"    -> UNEXPECTED: cross-block coupling does not lift chi_k.")
        print(f"       Either cross-block ERIs are still vanishing structurally or")
        print(f"       balanced-coupled has special structure that preserves chi=2.")
    else:
        print(f"    -> Anomalous: balanced rank < composed rank. Investigate.")

    out = {
        'molecule': 'LiH',
        'builder': 'balanced_coupled',
        'M': M, 'Q': Q,
        'n_pauli': len(qop.terms),
        'subblock_boundaries': boundaries,
        'profile': profile,
        'boundary_ranks': boundary_ranks,
        'interior_ranks': interior_ranks,
        'composed_boundary_rank_reference': composed_boundary,
    }
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s1_negative_control.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
