"""Probe cross-L linear dependences in the multipole basis.

Hypothesis update: per-L rank counting (sum_L rank(R^L)) works at n_max=2
(7=7), but at n_max=3 predicts 28 vs measured 24. The R^L matrices are all
full rank, so the missing 4 ranks must come from linear dependences between
v-vectors at DIFFERENT (L, density) pairs that share orbital-pair support.

For each (L, density alpha), build the (pr)-space vector
    v^{L,alpha}(p, r) = c^L(l_p, m_p; l_r, m_r)
                       for orbital pairs (p, r) with (n_p, l_p, n_r, l_r) ~ alpha,
                       zero otherwise.

Stack all v^{L,alpha} vectors. Compute the rank of the stacked matrix.

If rank(stacked) == DF rank, F3 is fully characterized: the missing ranks
come from cross-L (L=0 vs L=2 at same density, etc.) linear dependences.

Output the null space basis as linear combinations of (L, density) channels.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple, Dict

import numpy as np

from geovac.composed_qubit import (
    _enumerate_states, _compute_rk_integrals_block, _build_eri_block,
    _ck_coefficient,
)


def build_atomic_block_eri(Z: float, max_n: int):
    states = _enumerate_states(max_n, l_min=0)
    rk_cache = _compute_rk_integrals_block(Z, states)
    eri_phys = _build_eri_block(Z, states, rk_cache)
    n = len(states)
    eri = np.zeros((n, n, n, n))
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] += val
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))
    return eri, states


def build_v_vectors(states):
    """For each (L, density alpha), build v^{L,alpha}(p, r) vector in (pr)-space.

    Density alpha = unordered (n_a, l_a, n_b, l_b).

    Returns:
        labels: list of (L, alpha) tuples in the order they were inserted
        V: matrix of shape (M^2, len(labels)) where each column is a v-vector
    """
    n = len(states)
    label_list = []
    v_dict = {}

    for p, (n_p, l_p, m_p) in enumerate(states):
        for r, (n_r, l_r, m_r) in enumerate(states):
            density = tuple(sorted([(n_p, l_p), (n_r, l_r)]))
            L_min = abs(l_p - l_r)
            L_max = l_p + l_r
            for L in range(L_min, L_max + 1):
                if (l_p + l_r + L) % 2 != 0:
                    continue
                M = m_p - m_r
                if abs(M) > L:
                    continue
                ck = _ck_coefficient(l_p, m_p, l_r, m_r, L)
                if abs(ck) < 1e-15:
                    continue
                key = (L, density)
                if key not in v_dict:
                    v_dict[key] = np.zeros(n * n)
                    label_list.append(key)
                v_dict[key][p * n + r] = ck

    V = np.column_stack([v_dict[k] for k in label_list])
    return label_list, V


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'df_cross_L_dependences.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    results = {}

    for max_n in [2, 3, 4]:
        label_outer = f"max_n={max_n}"
        print(f"\n{'='*60}\n{label_outer}\n{'='*60}")

        eri, states = build_atomic_block_eri(3.0, max_n)
        n = len(states)

        G2 = eri.reshape(n * n, n * n)
        G2 = 0.5 * (G2 + G2.T)
        sv_g = np.linalg.svd(G2, compute_uv=False)
        df_rank = int(np.sum(sv_g > 1e-10 * max(sv_g[0], 1e-30)))

        labels, V = build_v_vectors(states)
        n_channels = len(labels)

        sv_v = np.linalg.svd(V, compute_uv=False)
        # Rank of V at machine threshold (V is a 0/Gaunt matrix, no scaling issues)
        rank_V_tight = int(np.sum(sv_v > 1e-12 * max(sv_v[0], 1e-30)))
        rank_V_loose = int(np.sum(sv_v > 1e-6 * max(sv_v[0], 1e-30)))

        print(f"  n orbitals = {n}, M^2 = {n*n}")
        print(f"  V shape: {V.shape} (n_channels = {n_channels})")
        print(f"  DF rank @ 1e-10: {df_rank}")
        print(f"  rank(V) @ 1e-12: {rank_V_tight}")
        print(f"  rank(V) @ 1e-6: {rank_V_loose}")
        print(f"  Top 15 SV(V): {[f'{x:.3e}' for x in sv_v[:15]]}")
        if rank_V_tight < n_channels:
            print(f"  ** {n_channels - rank_V_tight} cross-L linear dependences in V **")

        # Print null space of V if rank-deficient
        nulls = []
        if rank_V_tight < n_channels:
            U_V, sv_V_full, Vh_V = np.linalg.svd(V, full_matrices=True)
            # Null space of V corresponds to right singular vectors with zero SV
            null_basis = Vh_V[rank_V_tight:, :]
            for k in range(null_basis.shape[0]):
                nv = null_basis[k]
                # Find significant components (relative to max)
                max_c = np.max(np.abs(nv))
                components = [(lab, c) for lab, c in zip(labels, nv) if abs(c) > 0.01 * max_c]
                components.sort(key=lambda x: -abs(x[1]))
                print(f"\n  Null vector {k}:")
                for lab, c in components[:8]:
                    print(f"    {c:+.4f} * v^{{L={lab[0]}, density={lab[1]}}}")
                nulls.append({
                    'index': k,
                    'top_components': [
                        {'L': lab[0], 'density': list(lab[1]), 'coefficient': float(c)}
                        for lab, c in components
                    ],
                })

        results[label_outer] = {
            'n_orbitals': n,
            'states': [list(s) for s in states],
            'df_rank_g2': df_rank,
            'n_multipole_channels': n_channels,
            'rank_V_tight': rank_V_tight,
            'rank_V_loose': rank_V_loose,
            'df_rank_equals_rank_V': df_rank == rank_V_tight,
            'rank_deficit': n_channels - rank_V_tight,
            'top_singular_values_V': [float(x) for x in sv_v[:30]],
            'null_vectors_of_V': nulls,
            'labels': [{'L': lab[0], 'density': list(lab[1])} for lab in labels],
        }

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
