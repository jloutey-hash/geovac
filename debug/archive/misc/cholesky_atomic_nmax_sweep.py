"""Cholesky test extended to larger basis (n_max in {2, 3, 4}).

Use the atomic-block ERI builder (same path as df_radial_density_count.py)
to extend the production-basis n_max=2 test to larger basis where the
F3 graded compression appears.

Confirms F1 (multipole-subspace containment) for Cholesky vectors, not
just DF vectors.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple

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


def pivoted_cholesky(G2: np.ndarray, eps: float = 1e-12) -> Tuple[np.ndarray, list]:
    """Standard chemistry pivoted Cholesky (Beebe-Linderberg 1977)."""
    n = G2.shape[0]
    work = G2.astype(np.float64).copy()
    diag = np.diag(work).copy()
    L = np.zeros((n, 0))
    pivots = []
    rank = 0
    while rank < n:
        max_idx = int(np.argmax(diag))
        max_val = float(diag[max_idx])
        if max_val < eps:
            break
        col = work[:, max_idx] / np.sqrt(max_val)
        L = np.column_stack([L, col]) if L.size else col.reshape(-1, 1)
        pivots.append(max_idx)
        rank += 1
        work = work - np.outer(col, col)
        diag = np.maximum(np.diag(work).copy(), 0.0)
    return L, pivots


def build_multipole_subspace_basis(states):
    n = len(states)
    v_dict = {}
    label_list = []
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
    Q, _ = np.linalg.qr(V)
    return Q, label_list


def main():
    results = {}

    for max_n in [2, 3, 4]:
        for Z in [3.0]:
            label = f"max_n={max_n}_Z={Z}"
            print(f"\n{'='*60}\n{label}\n{'='*60}")

            eri, states = build_atomic_block_eri(Z, max_n)
            n = len(states)
            G2 = eri.reshape(n * n, n * n)
            G2 = 0.5 * (G2 + G2.T)

            sv = np.linalg.svd(G2, compute_uv=False)
            df_rank_strict = int(np.sum(sv > 1e-10 * max(sv[0], 1e-30)))
            df_rank_loose = int(np.sum(sv > 1e-6 * max(sv[0], 1e-30)))

            # Cholesky at two thresholds
            L_tight, _ = pivoted_cholesky(G2, eps=1e-10 * max(sv[0], 1e-30))
            L_loose, _ = pivoted_cholesky(G2, eps=1e-6 * max(sv[0], 1e-30))
            cholesky_rank_tight = L_tight.shape[1]
            cholesky_rank_loose = L_loose.shape[1]

            # Reconstruction
            recon_err_tight = (
                np.linalg.norm(G2 - L_tight @ L_tight.T) /
                max(np.linalg.norm(G2), 1e-30)
            )

            # Multipole subspace
            Q_mp, ch_labels = build_multipole_subspace_basis(states)
            n_channels = len(ch_labels)

            # Subspace overlap test
            min_overlap = 1.0
            max_residual = 0.0
            for k in range(cholesky_rank_tight):
                v = L_tight[:, k]
                v_unit = v / max(np.linalg.norm(v), 1e-30)
                proj_norm_sq = float(np.linalg.norm(Q_mp.T @ v_unit) ** 2)
                residual = float(np.linalg.norm(v_unit - Q_mp @ (Q_mp.T @ v_unit)))
                min_overlap = min(min_overlap, proj_norm_sq)
                max_residual = max(max_residual, residual)

            # Subspace match: Cholesky span vs SVD top-rank span
            Q_chol, _ = np.linalg.qr(L_tight)
            U_svd, _, _ = np.linalg.svd(G2)
            U_top = U_svd[:, :df_rank_strict]
            subspace_overlap = float(np.trace((Q_chol @ Q_chol.T) @ (U_top @ U_top.T)))

            print(f"  n orbitals: {n}")
            print(f"  DF rank @ 1e-10: {df_rank_strict}")
            print(f"  DF rank @ 1e-6: {df_rank_loose}")
            print(f"  Cholesky rank @ 1e-10: {cholesky_rank_tight}")
            print(f"  Cholesky rank @ 1e-6: {cholesky_rank_loose}")
            print(f"  Multipole channels: {n_channels}")
            print(f"  Cholesky reconstruction error: {recon_err_tight:.4e}")
            print(f"  Min Cholesky-multipole overlap squared: {min_overlap:.12f}")
            print(f"  Max Cholesky-multipole residual: {max_residual:.4e}")
            print(f"  Subspace overlap trace (Cholesky vs SVD top-{df_rank_strict}): {subspace_overlap:.6f}")
            print(f"    Expected (same subspace, ranks match): {float(min(df_rank_strict, cholesky_rank_tight)):.6f}")
            match = abs(subspace_overlap - float(min(df_rank_strict, cholesky_rank_tight))) < 1e-8
            print(f"    Match: {match}")

            results[label] = {
                'max_n': max_n,
                'Z': Z,
                'n_orbitals': n,
                'df_rank_1e-10': df_rank_strict,
                'df_rank_1e-6': df_rank_loose,
                'cholesky_rank_1e-10': cholesky_rank_tight,
                'cholesky_rank_1e-6': cholesky_rank_loose,
                'multipole_channels': n_channels,
                'reconstruction_error': float(recon_err_tight),
                'min_overlap_squared': float(min_overlap),
                'max_residual': float(max_residual),
                'subspace_overlap_trace': subspace_overlap,
                'subspace_match': bool(match),
            }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'cholesky_atomic_nmax_sweep.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
