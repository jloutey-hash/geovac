"""Directly verify the multipole decomposition reconstructs G2.

Build the (L, density)-indexed v vectors and R^L matrices, compute
G2_predicted = sum_L A_L R^L A_L^T where A_L stacks v^{L, alpha} as columns,
compare to G2 from the actual ERI tensor.

Also compute rank carefully. The standard symmetric-positive-semidefinite
result says rank(V K V^T) = rank(K) when V has full column rank. If the
predicted rank doesn't match the empirical rank, something in the
decomposition is off.

Specifically check:
1. Does G2_predicted == G2_actual to machine precision?
2. If yes, why does rank(G2) < sum_L rank(R^L)?
3. The likely answer: the v vectors at different L are not linearly
   independent in the (pr) space, i.e., V (with M-pooling) has rank < n_channels.
"""

from __future__ import annotations

import numpy as np

from geovac.composed_qubit import (
    _enumerate_states, _compute_rk_integrals_block, _build_eri_block,
    _ck_coefficient,
)
from geovac.hypergeometric_slater import compute_rk_float


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


def enumerate_densities_per_L(states):
    """Distinct unordered (n_a, l_a)-(n_b, l_b) density pairs at each Gaunt-allowed L."""
    nl_set = sorted(set((n, l) for n, l, _m in states))
    by_L = {}
    for i, (n_a, l_a) in enumerate(nl_set):
        for j, (n_b, l_b) in enumerate(nl_set):
            if j < i:
                continue
            L_min = abs(l_a - l_b)
            L_max = l_a + l_b
            for L in range(L_min, L_max + 1):
                if (l_a + l_b + L) % 2 != 0:
                    continue
                by_L.setdefault(L, []).append((n_a, l_a, n_b, l_b))
    for L in by_L:
        by_L[L] = sorted(by_L[L])
    return by_L


def build_v_matrix_per_L(states, L, densities):
    """Build matrix A_L of shape (M^2, len(densities)) where
    A_L[:, alpha] = v^{L, alpha}(p, r) = c^L(l_p, m_p; l_r, m_r) when
    density (p,r) matches alpha (with both orderings), else 0.

    No M-pooling distinction: the column accumulates Gaunt values across
    all M (different m_p, m_r) for the given density.
    """
    n = len(states)
    A = np.zeros((n * n, len(densities)))
    density_idx = {tuple(sorted([(d[0], d[1]), (d[2], d[3])])): i for i, d in enumerate(densities)}

    for p, (n_p, l_p, m_p) in enumerate(states):
        for r, (n_r, l_r, m_r) in enumerate(states):
            density_key = tuple(sorted([(n_p, l_p), (n_r, l_r)]))
            if density_key not in density_idx:
                continue
            i_d = density_idx[density_key]
            L_min_pr = abs(l_p - l_r)
            L_max_pr = l_p + l_r
            if not (L_min_pr <= L <= L_max_pr):
                continue
            if (l_p + l_r + L) % 2 != 0:
                continue
            ck = _ck_coefficient(l_p, m_p, l_r, m_r, L)
            if abs(ck) > 1e-15:
                A[p * n + r, i_d] = ck

    return A


def build_R_L_matrix(L, densities):
    """R^L(alpha, beta) with alpha = (n_a, l_a)-(n_b, l_b) on electron 1,
       beta = (n_c, l_c)-(n_d, l_d) on electron 2.

       compute_rk_float convention: first 4 args = electron 1 pair (1,3),
       next 4 args = electron 2 pair (2,4).
    """
    n_d_count = len(densities)
    R = np.zeros((n_d_count, n_d_count))
    for i, (n_a, l_a, n_b, l_b) in enumerate(densities):
        for j, (n_c, l_c, n_d_, l_d_) in enumerate(densities):
            try:
                val = compute_rk_float(n_a, l_a, n_b, l_b, n_c, l_c, n_d_, l_d_, L)
                R[i, j] = val
            except Exception:
                R[i, j] = 0.0
    return 0.5 * (R + R.T)


def main():
    for max_n in [2, 3]:
        print(f"\n{'='*70}\nmax_n = {max_n}\n{'='*70}")
        eri, states = build_atomic_block_eri(1.0, max_n)
        n = len(states)
        G2 = eri.reshape(n * n, n * n)
        G2 = 0.5 * (G2 + G2.T)

        sv_g = np.linalg.svd(G2, compute_uv=False)
        df_rank = int(np.sum(sv_g > 1e-10 * max(sv_g[0], 1e-30)))
        print(f"  M={n}, M^2={n*n}, DF rank @ 1e-10 = {df_rank}")
        print(f"  Top 15 sv of G2: {[f'{x:.3e}' for x in sv_g[:15]]}")

        densities_per_L = enumerate_densities_per_L(states)

        # Build V K V^T predicted
        G2_predicted = np.zeros_like(G2)
        all_A_columns = []
        all_K_diag = []

        for L in sorted(densities_per_L):
            densities = densities_per_L[L]
            A_L = build_v_matrix_per_L(states, L, densities)  # shape (M^2, n_d)
            R_L = build_R_L_matrix(L, densities)  # shape (n_d, n_d)

            contrib = A_L @ R_L @ A_L.T
            G2_predicted += contrib

            print(f"\n  L={L}: n_d={len(densities)}, ||A_L||_F={np.linalg.norm(A_L):.4f}, "
                  f"rank(A_L) @ 1e-10 = {int(np.sum(np.linalg.svd(A_L, compute_uv=False) > 1e-10 * np.linalg.svd(A_L, compute_uv=False)[0]))}")
            print(f"        rank(R_L) @ 1e-10 = {int(np.sum(np.linalg.svd(R_L, compute_uv=False) > 1e-10 * np.linalg.svd(R_L, compute_uv=False)[0]))}")
            sv_contrib = np.linalg.svd(contrib, compute_uv=False)
            rank_contrib = int(np.sum(sv_contrib > 1e-10 * max(sv_contrib[0], 1e-30)))
            print(f"        rank(A_L R_L A_L^T) @ 1e-10 = {rank_contrib}")

            all_A_columns.append(A_L)
            for d_i, density in enumerate(densities):
                all_K_diag.append((L, density, R_L[d_i, d_i]))

        # Check decomposition matches G2
        err = G2 - G2_predicted
        rel_err = np.linalg.norm(err) / max(np.linalg.norm(G2), 1e-30)
        print(f"\n  || G2 - V K V^T ||_F / ||G2||_F = {rel_err:.4e}")

        # Stack ALL A_L's into full V
        V_full = np.hstack(all_A_columns)
        sv_V = np.linalg.svd(V_full, compute_uv=False)
        rank_V = int(np.sum(sv_V > 1e-10 * max(sv_V[0], 1e-30)))
        rank_V_strict = int(np.sum(sv_V > 1e-6 * max(sv_V[0], 1e-30)))
        print(f"\n  Total V shape: {V_full.shape}")
        print(f"  rank(V_full) @ 1e-10 = {rank_V}")
        print(f"  rank(V_full) @ 1e-6 = {rank_V_strict}")
        print(f"  Top 15 sv(V_full): {[f'{x:.3e}' for x in sv_V[:15]]}")
        print(f"  Smallest 10 sv(V_full): {[f'{x:.3e}' for x in sv_V[-10:]]}")

        print(f"\n  EMPIRICAL DF rank = {df_rank}")
        print(f"  Number of channels (cols of V) = {V_full.shape[1]}")
        print(f"  rank(V_full) numerical = {rank_V}")


if __name__ == '__main__':
    main()
