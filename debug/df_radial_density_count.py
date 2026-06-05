"""Follow-up: predict DF rank from distinct-radial-density count.

The observation from df_vs_multipole_rank_test.py:
  DF rank = 7 for every sub-block with orbital basis (1s, 2s, 2p_{-1,0,+1}).

Hypothesis: rank = number of distinct unordered radial-density products
   rho_alpha(r) = R_{n_p l_p}(r) * R_{n_r l_r}(r)
weighted by Gaunt-allowed L's.

For the basis {1s, 2s, 2p}:
  L=0: (1s)^2, (1s)(2s), (2s)^2, (2p)^2 = 4 densities
  L=1: (1s)(2p), (2s)(2p)             = 2 densities
  L=2: (2p)^2                          = 1 density
  Total = 7  -- match.

Verify this prediction holds for n_max=3 (basis includes 3s, 3p, 3d).
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple, Dict

import numpy as np

from geovac.composed_qubit import (
    build_composed_hamiltonian, _enumerate_states, _compute_rk_integrals_block,
    _build_eri_block,
)
from geovac.molecular_spec import MolecularSpec, OrbitalBlock


def count_distinct_radial_densities(states: List[Tuple[int, int, int]]) -> Dict:
    """
    Count distinct radial densities rho_{n_p l_p, n_r l_r}(r) = R_{n_p l_p} * R_{n_r l_r}
    weighted by Gaunt-allowed L for each (n_p l_p, n_r l_r) pair.

    Because rho_{a,b} = rho_{b,a} as a function of r, we count UNORDERED radial
    pairs. The total predicted DF rank is the sum over L of the number of
    distinct unordered (n_p l_p, n_r l_r) radial-density products that have
    a Gaunt-allowed coupling at L.
    """
    nl_set = sorted(set((n, l) for n, l, _m in states))

    # Unordered radial pairs
    unordered_radial_pairs = set()
    for i, (n_a, l_a) in enumerate(nl_set):
        for j, (n_b, l_b) in enumerate(nl_set):
            if (n_b, l_b, n_a, l_a) in unordered_radial_pairs:
                continue
            unordered_radial_pairs.add((n_a, l_a, n_b, l_b))

    densities_per_L = {}
    for (n_a, l_a, n_b, l_b) in unordered_radial_pairs:
        L_min = abs(l_a - l_b)
        L_max = l_a + l_b
        for L in range(L_min, L_max + 1):
            if (l_a + l_b + L) % 2 != 0:
                continue
            densities_per_L.setdefault(L, set()).add(
                tuple(sorted([(n_a, l_a), (n_b, l_b)]))
            )

    counts = {L: len(d) for L, d in densities_per_L.items()}
    total = sum(counts.values())

    return {
        'total_predicted_DF_rank': total,
        'densities_per_L': {L: sorted(d) for L, d in densities_per_L.items()},
        'counts_per_L': counts,
        'nl_set': nl_set,
    }


def compute_df_rank(eri_sub: np.ndarray, rel_thr: float = 1e-10) -> int:
    n = eri_sub.shape[0]
    G2 = eri_sub.reshape(n * n, n * n)
    G2 = 0.5 * (G2 + G2.T)
    sv = np.linalg.svd(G2, compute_uv=False)
    return int(np.sum(sv > rel_thr * max(sv[0], 1e-30))), [float(x) for x in sv]


def build_atomic_block_eri(Z: float, max_n: int) -> np.ndarray:
    """Build a single-block ERI for an atom at charge Z with quantum-number cutoff."""
    states = _enumerate_states(max_n, l_min=0)
    rk_cache = _compute_rk_integrals_block(Z, states)
    eri_phys = _build_eri_block(Z, states, rk_cache)

    n = len(states)
    eri = np.zeros((n, n, n, n))
    # Convert (a, b, c, d) physicist -> (a, c, b, d) chemist
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] += val
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))
    return eri, states


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'df_radial_density_count.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    results = {}

    # Test for n_max = 2 (default composed) and n_max = 3
    for max_n in [2, 3, 4]:
        for Z in [3.0, 6.0]:  # Li-like and O-like charges
            label = f"max_n={max_n}_Z={Z}"
            print(f"\n{'='*60}\n{label}\n{'='*60}")

            try:
                eri_sub, states = build_atomic_block_eri(Z, max_n)
                n = len(states)
                print(f"Orbitals: {n}, states: {states}")

                pred = count_distinct_radial_densities(states)
                print(f"\nPredicted DF rank (distinct radial densities): {pred['total_predicted_DF_rank']}")
                print(f"Counts per L: {pred['counts_per_L']}")
                for L, d in pred['densities_per_L'].items():
                    print(f"  L={L}: {[list(p) for p in d]}")

                rank, sv = compute_df_rank(eri_sub, rel_thr=1e-10)
                print(f"\nMeasured DF rank @ s/s_max > 1e-10: {rank}")
                print(f"Top SV: {[f'{x:.3e}' for x in sv[:min(20, len(sv))]]}")
                # Print also rank at looser threshold to spot near-zero modes
                rank_1e6, _ = compute_df_rank(eri_sub, rel_thr=1e-6)
                print(f"DF rank @ s/s_max > 1e-6: {rank_1e6}")

                match = (rank == pred['total_predicted_DF_rank'])
                print(f"\n*** PREDICTION MATCH: {match} ***")

                results[label] = {
                    'max_n': max_n,
                    'Z': Z,
                    'n_orbitals': n,
                    'states': [list(s) for s in states],
                    'predicted_rank': pred['total_predicted_DF_rank'],
                    'predicted_counts_per_L': pred['counts_per_L'],
                    'predicted_densities_per_L': {
                        L: [list(p) for p in d] for L, d in pred['densities_per_L'].items()
                    },
                    'measured_rank_at_1e-10': rank,
                    'measured_rank_at_1e-6': rank_1e6,
                    'top_singular_values': sv[:20],
                    'prediction_match': match,
                }
            except Exception as e:
                import traceback
                traceback.print_exc()
                results[label] = {'error': str(e)}

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
