"""Per-L rank decomposition of the GeoVac ERI tensor.

Hypothesis (F3 grounding):
  DF rank = sum_L rank(R^L_density_by_density)
where R^L is the radial-integral matrix between distinct unordered radial
densities at angular momentum L. The "extra compression" at higher n_max
comes from R^L matrices that are NOT full rank.

This script computes rank(R^L) at each L by:
1. Enumerating distinct unordered radial densities at each L.
2. Building the explicit R^L matrix using the algebraic R^k_(...) integral
   evaluator (compute_rk_float / compute_rk_exact).
3. Computing the SVD of each R^L matrix.
4. Comparing sum(rank(R^L)) to the empirical DF rank.

If the empirical DF rank matches sum(rank(R^L)) EXACTLY at every n_max,
then F3 is grounded: the rank deficiency is fully explained by R^L rank.

We then characterize the null space of R^L to identify the algebraic source.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple, Dict

import numpy as np

from geovac.composed_qubit import (
    _enumerate_states, _compute_rk_integrals_block, _build_eri_block,
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


def compute_df_rank(eri_sub: np.ndarray, rel_thr: float = 1e-10) -> Tuple[int, List[float]]:
    n = eri_sub.shape[0]
    G2 = eri_sub.reshape(n * n, n * n)
    G2 = 0.5 * (G2 + G2.T)
    sv = np.linalg.svd(G2, compute_uv=False)
    rank = int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))
    return rank, sv


def enumerate_distinct_densities_per_L(states):
    """Return dict L -> list of unordered (n_a, l_a, n_b, l_b) tuples (a <= b)
    that produce Gaunt-allowed contributions at angular momentum L.
    """
    nl_set = sorted(set((n, l) for n, l, _m in states))
    by_L = {}
    for i, (n_a, l_a) in enumerate(nl_set):
        for j, (n_b, l_b) in enumerate(nl_set):
            if j < i:
                continue  # unordered
            L_min = abs(l_a - l_b)
            L_max = l_a + l_b
            for L in range(L_min, L_max + 1):
                if (l_a + l_b + L) % 2 != 0:
                    continue
                by_L.setdefault(L, []).append((n_a, l_a, n_b, l_b))
    # Sort each list lexicographically
    for L in by_L:
        by_L[L] = sorted(by_L[L])
    return by_L


def build_R_L_matrix(L: int, densities: List[Tuple[int, int, int, int]]):
    """Build R^L matrix between distinct radial densities.

    Each density is an unordered orbital pair (n_a, l_a, n_b, l_b)
    representing rho_alpha(r) = R_{n_a l_a}(r) * R_{n_b l_b}(r).

    compute_rk_float(n1, l1, n3, l3, n2, l2, n4, l4, k) computes R^k at the
    hydrogenic Z=1 normalization; rank structure is Z-independent so we omit
    Z entirely.  Chemist Slater convention: R^k(a, b; c, d) with electrons
    1 -> rows (a, c), 2 -> cols (b, d).
    """
    n = len(densities)
    R = np.zeros((n, n))
    for i, (n_a, l_a, n_b, l_b) in enumerate(densities):
        for j, (n_c, l_c, n_d, l_d) in enumerate(densities):
            try:
                # compute_rk_float(n1, l1, n3, l3, n2, l2, n4, l4, k):
                # electron 1 pair = positions (1, 3) = (n_a, l_a)-(n_b, l_b)
                # electron 2 pair = positions (2, 4) = (n_c, l_c)-(n_d, l_d)
                val = compute_rk_float(n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, L)
                R[i, j] = val
            except Exception as e:
                # Gaunt-forbidden quartet -> contribution is zero, not an error
                R[i, j] = 0.0
    R_sym = 0.5 * (R + R.T)
    asym = np.linalg.norm(R - R_sym) / max(np.linalg.norm(R), 1e-30)
    return R_sym, asym


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'df_per_L_rank.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    results = {}

    for max_n in [2, 3, 4]:
        for Z in [3.0]:
            label = f"max_n={max_n}_Z={Z}"
            print(f"\n{'='*60}\n{label}\n{'='*60}")

            try:
                eri, states = build_atomic_block_eri(Z, max_n)
                n = len(states)
                df_rank, sv = compute_df_rank(eri, 1e-10)
                df_rank_loose, _ = compute_df_rank(eri, 1e-6)
                print(f"Orbitals: {n}")
                print(f"DF rank @ 1e-10: {df_rank}")
                print(f"DF rank @ 1e-6: {df_rank_loose}")
                print(f"Top 15 SV: {[f'{x:.3e}' for x in sv[:15]]}")

                densities_per_L = enumerate_distinct_densities_per_L(states)
                print(f"\nDistinct unordered radial densities per L:")
                for L in sorted(densities_per_L):
                    print(f"  L={L}: {len(densities_per_L[L])} densities")
                    for d in densities_per_L[L]:
                        print(f"    {d}")

                ranks_per_L = {}
                rank_sum = 0
                R_matrices = {}
                for L in sorted(densities_per_L):
                    densities = densities_per_L[L]
                    R_L, asym = build_R_L_matrix(L, densities)
                    sv_R = np.linalg.svd(R_L, compute_uv=False)
                    rank_tight = int(np.sum(sv_R > 1e-10 * max(sv_R[0], 1e-30)))
                    rank_loose = int(np.sum(sv_R > 1e-6 * max(sv_R[0], 1e-30)))

                    print(f"\n  L={L}: {len(densities)} densities, |R^L| shape ({len(densities)}, {len(densities)})")
                    print(f"    R^L symmetric residual: {asym:.2e}")
                    print(f"    rank(R^L) @ 1e-10: {rank_tight}")
                    print(f"    rank(R^L) @ 1e-6: {rank_loose}")
                    print(f"    SV of R^L: {[f'{x:.3e}' for x in sv_R]}")

                    ranks_per_L[L] = {
                        'n_densities': len(densities),
                        'rank_1e-10': rank_tight,
                        'rank_1e-6': rank_loose,
                        'singular_values': [float(x) for x in sv_R],
                    }
                    rank_sum += rank_tight
                    R_matrices[L] = R_L

                print(f"\n  Sum of rank(R^L) @ 1e-10: {rank_sum}")
                print(f"  Empirical DF rank @ 1e-10: {df_rank}")
                print(f"  Match: {rank_sum == df_rank}")

                # Save null vectors for L's with rank deficiency
                null_vectors_per_L = {}
                for L, R_L in R_matrices.items():
                    n_d = len(densities_per_L[L])
                    if ranks_per_L[L]['rank_1e-10'] < n_d:
                        U_R, sv_R, Vh_R = np.linalg.svd(R_L)
                        rk = ranks_per_L[L]['rank_1e-10']
                        null_vecs = Vh_R[rk:, :].tolist()
                        null_vectors_per_L[L] = {
                            'densities': [list(d) for d in densities_per_L[L]],
                            'null_basis': null_vecs,
                            'null_singular_values': [float(x) for x in sv_R[rk:]],
                        }
                        print(f"\n  L={L} null space ({n_d - ranks_per_L[L]['rank_1e-10']} dim):")
                        for k, nv in enumerate(null_vecs):
                            print(f"    Null vector {k} (sv={sv_R[rk+k]:.3e}):")
                            for d, c in zip(densities_per_L[L], nv):
                                if abs(c) > 1e-6:
                                    print(f"      {c:+.4f} * rho{d}")

                results[label] = {
                    'max_n': max_n,
                    'Z': Z,
                    'n_orbitals': n,
                    'df_rank_1e-10': df_rank,
                    'df_rank_1e-6': df_rank_loose,
                    'predicted_rank_from_R_L': rank_sum,
                    'match': rank_sum == df_rank,
                    'per_L': ranks_per_L,
                    'null_vectors_per_L': null_vectors_per_L,
                    'top_sv_G2': [float(x) for x in sv[:min(30, len(sv))]],
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
