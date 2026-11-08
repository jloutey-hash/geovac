"""Inspect the structure of the rank-7 DF factors at n_max=2.

For each of the 7 leading singular vectors of G2 = ERI reshaped to (M^2, M^2),
check whether it lives entirely within a single (L, density) channel.

Concretely: each left-singular-vector u in (p,r)-space should be supported only
on orbital pairs (p, r) where c^L(l_p, m_p; l_r, m_r) is nonzero for the
expected (L, density) of that mode.

Output: per-mode support distribution across (L, radial-density-pair) channels.
If each mode is supported on exactly one channel, that is unambiguous evidence
that DF on the GeoVac ERI tensor IS the multipole expansion (same basis, same
basis vectors).
"""

from __future__ import annotations

import os
import json
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


def build_multipole_projectors(states: List[Tuple[int, int, int]]):
    """For each (L, M, ordered radial-pair), build the (p,r)-space vector
        v^{L,M,radpair}_{(p,r)} = c^L(l_p, m_p; l_r, m_r) * indicator((n_p,l_p,n_r,l_r) = radpair)
        (with delta(m_p - m_r, M) implicit in the Gaunt coefficient)

    Sum over m within each radial-density-orbit gives the canonical (L, density)
    basis vector. These are the analytic multipole basis vectors.
    """
    n = len(states)
    projectors = {}  # key: (L, density_pair) where density_pair is unordered

    for p, (n_p, l_p, m_p) in enumerate(states):
        for r, (n_r, l_r, m_r) in enumerate(states):
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
                density = tuple(sorted([(n_p, l_p), (n_r, l_r)]))
                key = (L, density)
                if key not in projectors:
                    projectors[key] = {
                        'vectors_by_M': {},
                        'L': L, 'density': density,
                    }
                vec = projectors[key]['vectors_by_M'].setdefault(
                    M, np.zeros(n * n)
                )
                vec[p * n + r] = ck

    # Normalize per (L, density) by summing-over-M-then-Frobenius-normalizing
    # Actually, what we want is: for each (L, density), there is a single
    # "channel" but it can carry up to (2L+1) M components.  The DF mode
    # picks out specific linear combinations.  We'll just report support.
    return projectors


def analyze_mode_support(u_mode: np.ndarray, projectors: Dict, n: int) -> Dict:
    """For a single (p,r)-space singular vector u_mode, compute its overlap with
    each (L, density) projector.  Return overlap distribution.
    """
    u = u_mode / max(np.linalg.norm(u_mode), 1e-30)
    overlaps = {}
    for key, info in projectors.items():
        L = info['L']
        density = info['density']
        # Form the SUBSPACE spanned by all M-components of this (L, density)
        M_vectors = np.array(list(info['vectors_by_M'].values()))
        if M_vectors.size == 0:
            continue
        # Orthonormalize this subspace
        Q, _ = np.linalg.qr(M_vectors.T)
        # Project u onto subspace
        proj_norm = np.linalg.norm(Q.T @ u)
        overlaps[f"L={L}_density={density}"] = float(proj_norm ** 2)
    return overlaps


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'df_factor_structure.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    results = {}

    Z = 3.0
    for max_n in [2, 3, 4]:
        label = f"max_n={max_n}_Z={Z}"
        print(f"\n{'='*60}\n{label}\n{'='*60}")

        eri, states = build_atomic_block_eri(Z, max_n)
        n = len(states)
        G2 = eri.reshape(n * n, n * n)
        G2 = 0.5 * (G2 + G2.T)
        U, sv, Vh = np.linalg.svd(G2)
        rank = int(np.sum(sv > 1e-10 * max(sv[0], 1e-30)))
        print(f"n_orbitals={n}, DF rank @ 1e-10 = {rank}")

        projectors = build_multipole_projectors(states)
        print(f"Multipole (L, density) channels: {len(projectors)}")
        for k in sorted(projectors.keys()):
            print(f"  L={k[0]}, density={k[1]}, |M-vectors|={len(projectors[k]['vectors_by_M'])}")

        mode_overlaps = []
        for k in range(min(rank, n * n)):
            u_mode = U[:, k]
            overlaps = analyze_mode_support(u_mode, projectors, n)
            total = sum(overlaps.values())
            top_channels = sorted(overlaps.items(), key=lambda x: -x[1])[:3]
            mode_overlaps.append({
                'mode': k,
                'sv': float(sv[k]),
                'total_overlap': float(total),
                'top_3_channels': top_channels,
            })
            print(f"\nMode {k}: sv={sv[k]:.4e}, total overlap with multipole subspace={total:.4f}")
            for ch, ov in top_channels:
                if ov > 1e-4:
                    print(f"  {ch}: {ov:.4f}")

        results[label] = {
            'n_orbitals': n,
            'states': [list(s) for s in states],
            'df_rank': rank,
            'n_multipole_channels': len(projectors),
            'top_singular_values': [float(x) for x in sv[:min(20, len(sv))]],
            'mode_overlaps': mode_overlaps,
        }

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
