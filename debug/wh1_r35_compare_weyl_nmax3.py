"""Comparison: truthful CH Weyl vs full Dirac at n_max=3.

The R3.5 truthful CH at n_max=3 surprisingly produced 96 finite distances
(out of 780 - 60 forced zeros - 624 inf = 96 finite). To check whether
the structural prediction (full == 2 * Weyl) holds, we run truthful CH
Weyl-only at n_max=3 and compare.

Predicted: Weyl alone at n_max=3 should give 96 / 4 ≈ 24 finite distances
(since chirality doubling + symmetry should multiply by ~4).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from geovac.spinor_operator_system import (
    SpinorTruncatedOperatorSystem,
    camporesi_higuchi_dirac_matrix,
    spinor_basis,
    spinor_dim,
    spinor_graph_distance_matrix,
    spinor_label_strings,
)
from geovac.connes_distance import compute_distance_matrix


print("=" * 70)
print("WH1-R3.5: truthful CH Weyl-only at n_max = 3 (comparison reference)")
print("=" * 70)

n_max = 3
O = SpinorTruncatedOperatorSystem(n_max)
print(f"  dim_H_Weyl = {O.dim_H}")
print(f"  dim(O_Weyl) = {O.dim}")

D = camporesi_higuchi_dirac_matrix(O.basis)
print(f"  D eigenvalues: {sorted(set(np.round(np.diag(D).real, 6).tolist()))}")

print("\n  Computing distance matrix...")
t0 = time.time()
dist = compute_distance_matrix(O, D=D, progress=False)
elapsed = time.time() - t0
print(f"  done in {elapsed:.1f} s")

N = O.dim_H
n_pairs = N * (N - 1) // 2
n_inf = sum(
    1 for i in range(N) for j in range(i + 1, N) if np.isinf(dist[i, j])
)
n_zero = sum(
    1 for i in range(N) for j in range(i + 1, N)
    if (not np.isinf(dist[i, j]) and dist[i, j] < 1e-7)
)
n_finite_nz = n_pairs - n_inf - n_zero
finite_values = [
    dist[i, j]
    for i in range(N) for j in range(i + 1, N)
    if not np.isinf(dist[i, j]) and dist[i, j] >= 1e-7
]
print(f"\n  pairs: total={n_pairs}, inf={n_inf}, zero={n_zero}, "
      f"finite={n_finite_nz}")
if finite_values:
    print(f"  range of finite nonzero: [{min(finite_values):.4f}, "
          f"{max(finite_values):.4f}]")

graph = spinor_graph_distance_matrix(O)
nz_pairs_d, nz_pairs_g = [], []
for i in range(N):
    for j in range(i + 1, N):
        d_ij = dist[i, j]
        if np.isinf(d_ij) or d_ij < 1e-7:
            continue
        nz_pairs_d.append(d_ij)
        nz_pairs_g.append(graph[i, j])

pearson_nz, spearman_nz = float("nan"), float("nan")
if len(nz_pairs_d) >= 3:
    pearson_nz, _ = pearsonr(nz_pairs_g, nz_pairs_d)
    spearman_nz, _ = spearmanr(nz_pairs_g, nz_pairs_d)
    print(f"  Pearson nz: {pearson_nz:.4f}, Spearman nz: {spearman_nz:.4f}")


def _serialize_matrix(M):
    rows = []
    for row in M:
        out_row = []
        for x in row:
            v = float(x)
            if np.isinf(v):
                out_row.append("inf")
            else:
                out_row.append(round(v, 8))
        rows.append(out_row)
    return rows


out = {
    "n_max": n_max,
    "construction": "weyl_only_truthful_CH_for_R3_5_comparison",
    "dim_H_Weyl": O.dim_H,
    "dim_O_Weyl": O.dim,
    "wall_time_seconds": elapsed,
    "n_pairs_total": n_pairs,
    "n_pairs_inf": n_inf,
    "n_pairs_zero_forced": n_zero,
    "n_pairs_finite_nz": n_finite_nz,
    "finite_range": (
        [float(min(finite_values)), float(max(finite_values))]
        if finite_values else None
    ),
    "pearson_nz_vs_graph": pearson_nz,
    "spearman_nz_vs_graph": spearman_nz,
    "labels": spinor_label_strings(O),
    "distance_matrix": _serialize_matrix(dist),
}
out_path = (
    project_root / "debug" / "data"
    / f"wh1_r35_weyl_truthful_nmax{n_max}.json"
)
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(out, indent=2))
print(f"\n  saved -> {out_path.relative_to(project_root)}")
