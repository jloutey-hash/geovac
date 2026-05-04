"""WH1-R3.5 — full Dirac (both chiralities) Connes distance computation.

Implements deliverable 2 of the WH1-R3.5 sprint plan: compute the
Connes distance matrix on the full-Dirac truncated operator system at
n_max=2 (and n_max=3 if time permits) under (a) the truthful CH
full Dirac (eigenvalues +/-(n_fock + 1/2) on the two chirality
sectors); and (b) the offdiag CH full Dirac (cross-chirality
couplings via standard SO(4) E1 selection rules).

The hypothesis (per debug/track_ts_a_gh_convergence_memo.md §6.2):
adding both chiralities introduces native cross-chirality off-diagonal
multiplier structure that breaks the n-degeneracy obstruction
without requiring an artificial perturbation.

Output: data dumps to debug/data/wh1_r35_full_dirac_*.json with the
distance matrices, infinities count, Pearson correlation with the
graph distance, and a one-paragraph verdict.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

# Make sure debug paths exist
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
    full_dirac_graph_distance_matrix,
    full_dirac_label_strings,
)
from geovac.connes_distance import compute_distance_matrix


def _serialize_matrix(M):
    """Convert numpy float matrix to nested list with 'inf' as JSON literal."""
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


def run_n_max_2_truthful():
    """Connes distance under truthful CH full Dirac at n_max=2."""
    print("=" * 70)
    print("WH1-R3.5: full Dirac (truthful CH) at n_max = 2")
    print("=" * 70)

    n_max = 2
    O = FullDiracTruncatedOperatorSystem(n_max)
    print(f"  dim_H_full = {O.dim_H} (expected {full_dirac_dim(n_max)})")
    print(f"  dim(O_full) = {O.dim}")
    print(f"  num multipliers = {len(O.multiplier_matrices)}")

    D = camporesi_higuchi_full_dirac_matrix(O.basis)
    eigvals = np.diag(D).real
    print(f"  D eigenvalues: {sorted(set(np.round(eigvals, 6).tolist()))}")

    print("\n  Computing distance matrix...")
    t0 = time.time()
    dist = compute_distance_matrix(O, D=D, progress=False)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f} s")

    N = O.dim_H
    n_pairs = N * (N - 1) // 2
    n_inf = 0
    n_zero = 0
    n_finite_nz = 0
    finite_values = []
    for i in range(N):
        for j in range(i + 1, N):
            d = dist[i, j]
            if np.isinf(d):
                n_inf += 1
            elif d < 1e-7:
                n_zero += 1
            else:
                n_finite_nz += 1
                finite_values.append(d)

    print(f"\n  Pair statistics:")
    print(f"    total off-diagonal pairs: {n_pairs}")
    print(f"    +infinity pairs: {n_inf}")
    print(f"    zero (forced) pairs: {n_zero}")
    print(f"    finite nonzero pairs: {n_finite_nz}")
    if finite_values:
        print(f"    range of finite nonzero: [{min(finite_values):.4f}, "
              f"{max(finite_values):.4f}]")

    # Graph distance
    graph = full_dirac_graph_distance_matrix(O)
    nz_pairs_d = []
    nz_pairs_g = []
    for i in range(N):
        for j in range(i + 1, N):
            d_ij = dist[i, j]
            if np.isinf(d_ij) or d_ij < 1e-7:
                continue
            nz_pairs_d.append(d_ij)
            nz_pairs_g.append(graph[i, j])

    pearson_nz = float("nan")
    spearman_nz = float("nan")
    if len(nz_pairs_d) >= 3:
        pearson_nz, _ = pearsonr(nz_pairs_g, nz_pairs_d)
        spearman_nz, _ = spearmanr(nz_pairs_g, nz_pairs_d)
        print(f"  Pearson nz vs graph distance: {pearson_nz:.4f}")
        print(f"  Spearman nz vs graph distance: {spearman_nz:.4f}")

    out = {
        "n_max": n_max,
        "construction": "full_dirac_truthful_CH",
        "dim_H_full": O.dim_H,
        "dim_O_full": O.dim,
        "num_multipliers": len(O.multiplier_matrices),
        "dirac_eigenvalues": sorted(set(np.round(eigvals, 6).tolist())),
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
        "labels": full_dirac_label_strings(O),
        "distance_matrix": _serialize_matrix(dist),
    }
    out_path = (
        project_root / "debug" / "data"
        / f"wh1_r35_full_dirac_truthful_nmax{n_max}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n  saved -> {out_path.relative_to(project_root)}")
    return out


def run_n_max_2_offdiag():
    """Connes distance under offdiag CH full Dirac at n_max=2.

    This Dirac has cross-chirality couplings AND within-chirality
    couplings, naturally breaking the n-degeneracy. It is the natural
    full-Dirac analog of R3.2 'CH+offdiag' mode."""
    print("=" * 70)
    print("WH1-R3.5: full Dirac (offdiag CH) at n_max = 2")
    print("=" * 70)

    n_max = 2
    O = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_offdiag_dirac_matrix(
        O.basis,
        diag_lifters=(1.0, 0.1, 0.005),
        offdiag_alpha=1.0,
        chirality_coupling=1.0,
    )
    print(f"  D shape: {D.shape}")
    print(f"  D Hermitian: ||D - D*|| = {np.linalg.norm(D - D.conj().T):.2e}")

    print("\n  Computing distance matrix...")
    t0 = time.time()
    dist = compute_distance_matrix(O, D=D, progress=False)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f} s")

    N = O.dim_H
    n_pairs = N * (N - 1) // 2
    n_inf = 0
    n_zero = 0
    n_finite_nz = 0
    finite_values = []
    for i in range(N):
        for j in range(i + 1, N):
            d = dist[i, j]
            if np.isinf(d):
                n_inf += 1
            elif d < 1e-7:
                n_zero += 1
            else:
                n_finite_nz += 1
                finite_values.append(d)

    print(f"\n  Pair statistics:")
    print(f"    total off-diagonal pairs: {n_pairs}")
    print(f"    +infinity pairs: {n_inf}")
    print(f"    zero (forced) pairs: {n_zero}")
    print(f"    finite nonzero pairs: {n_finite_nz}")
    if finite_values:
        print(f"    range of finite nonzero: [{min(finite_values):.4f}, "
              f"{max(finite_values):.4f}]")

    # Graph distance
    graph = full_dirac_graph_distance_matrix(O)
    nz_pairs_d = []
    nz_pairs_g = []
    for i in range(N):
        for j in range(i + 1, N):
            d_ij = dist[i, j]
            if np.isinf(d_ij) or d_ij < 1e-7:
                continue
            nz_pairs_d.append(d_ij)
            nz_pairs_g.append(graph[i, j])

    pearson_nz = float("nan")
    spearman_nz = float("nan")
    if len(nz_pairs_d) >= 3:
        pearson_nz, _ = pearsonr(nz_pairs_g, nz_pairs_d)
        spearman_nz, _ = spearmanr(nz_pairs_g, nz_pairs_d)
        print(f"  Pearson nz vs graph distance: {pearson_nz:.4f}")
        print(f"  Spearman nz vs graph distance: {spearman_nz:.4f}")

    out = {
        "n_max": n_max,
        "construction": "full_dirac_offdiag_CH",
        "diag_lifters": [1.0, 0.1, 0.005],
        "offdiag_alpha": 1.0,
        "chirality_coupling": 1.0,
        "dim_H_full": O.dim_H,
        "dim_O_full": O.dim,
        "num_multipliers": len(O.multiplier_matrices),
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
        "labels": full_dirac_label_strings(O),
        "distance_matrix": _serialize_matrix(dist),
    }
    out_path = (
        project_root / "debug" / "data"
        / f"wh1_r35_full_dirac_offdiag_nmax{n_max}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n  saved -> {out_path.relative_to(project_root)}")
    return out


def run_n_max_3_truthful():
    """Connes distance under truthful CH full Dirac at n_max=3."""
    print("=" * 70)
    print("WH1-R3.5: full Dirac (truthful CH) at n_max = 3")
    print("=" * 70)

    n_max = 3
    O = FullDiracTruncatedOperatorSystem(n_max)
    print(f"  dim_H_full = {O.dim_H}")

    D = camporesi_higuchi_full_dirac_matrix(O.basis)
    print("\n  Computing distance matrix... (this may be slow)")
    t0 = time.time()
    dist = compute_distance_matrix(O, D=D, progress=False)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f} s")

    N = O.dim_H
    n_pairs = N * (N - 1) // 2
    n_inf = sum(1 for i in range(N) for j in range(i + 1, N) if np.isinf(dist[i, j]))
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
    print(f"  pairs: total={n_pairs}, inf={n_inf}, zero={n_zero}, "
          f"finite={n_finite_nz}")

    graph = full_dirac_graph_distance_matrix(O)
    nz_pairs_d = []
    nz_pairs_g = []
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

    out = {
        "n_max": n_max,
        "construction": "full_dirac_truthful_CH",
        "dim_H_full": O.dim_H,
        "dim_O_full": O.dim,
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
        "labels": full_dirac_label_strings(O),
        "distance_matrix": _serialize_matrix(dist),
    }
    out_path = (
        project_root / "debug" / "data"
        / f"wh1_r35_full_dirac_truthful_nmax{n_max}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"  saved -> {out_path.relative_to(project_root)}")
    return out


def run_n_max_3_offdiag():
    """Connes distance under offdiag CH full Dirac at n_max=3."""
    print("=" * 70)
    print("WH1-R3.5: full Dirac (offdiag CH) at n_max = 3")
    print("=" * 70)

    n_max = 3
    O = FullDiracTruncatedOperatorSystem(n_max)
    print(f"  dim_H_full = {O.dim_H}")

    D = camporesi_higuchi_offdiag_dirac_matrix(
        O.basis,
        diag_lifters=(1.0, 0.1, 0.005),
        offdiag_alpha=1.0,
        chirality_coupling=1.0,
    )

    print("\n  Computing distance matrix... (this may be slow)")
    t0 = time.time()
    dist = compute_distance_matrix(O, D=D, progress=False)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f} s")

    N = O.dim_H
    n_pairs = N * (N - 1) // 2
    n_inf = sum(1 for i in range(N) for j in range(i + 1, N) if np.isinf(dist[i, j]))
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
    print(f"  pairs: total={n_pairs}, inf={n_inf}, zero={n_zero}, "
          f"finite={n_finite_nz}")

    graph = full_dirac_graph_distance_matrix(O)
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

    out = {
        "n_max": n_max,
        "construction": "full_dirac_offdiag_CH",
        "diag_lifters": [1.0, 0.1, 0.005],
        "offdiag_alpha": 1.0,
        "chirality_coupling": 1.0,
        "dim_H_full": O.dim_H,
        "dim_O_full": O.dim,
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
        "labels": full_dirac_label_strings(O),
        "distance_matrix": _serialize_matrix(dist),
    }
    out_path = (
        project_root / "debug" / "data"
        / f"wh1_r35_full_dirac_offdiag_nmax{n_max}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"  saved -> {out_path.relative_to(project_root)}")
    return out


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mode", choices=["t2", "o2", "t3", "o3", "all"], default="all",
        help="t2/o2/t3/o3 = truthful/offdiag at n_max=2/3; all = all four"
    )
    args = p.parse_args()

    results = {}
    if args.mode in ("t2", "all"):
        results["truthful_nmax2"] = run_n_max_2_truthful()
    if args.mode in ("o2", "all"):
        results["offdiag_nmax2"] = run_n_max_2_offdiag()
    if args.mode in ("t3", "all"):
        results["truthful_nmax3"] = run_n_max_3_truthful()
    if args.mode in ("o3", "all"):
        results["offdiag_nmax3"] = run_n_max_3_offdiag()

    print("\n\n" + "=" * 70)
    print("WH1-R3.5 SUMMARY")
    print("=" * 70)
    for key, r in results.items():
        if r is None:
            continue
        print(f"  {key}: dim_H={r['dim_H_full']}, "
              f"inf={r['n_pairs_inf']}/{r['n_pairs_total']}, "
              f"finite={r['n_pairs_finite_nz']}, "
              f"Pearson={r.get('pearson_nz_vs_graph', 'NaN'):.4f}")
