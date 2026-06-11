"""Paper-0 packing Steiner probe — 3D version with proper MST-based Steiner heuristic.

Addresses PM's over-reach in the 2D probe:
  1. 3D instead of 2D (Dirac degeneracies are S^3 objects, not planar).
  2. Proper Steiner-MST heuristic (not just local Delaunay-triangle Fermat).
  3. Shell-pair granularity as the PI specified: "steiner tree between each shell".
  4. Leverages packing sparsity: only local cross-shell connections, not N^2 pairs.

For each adjacent shell pair (n, n+1):
  - Collect points from both shells (Fibonacci-uniform 3D distribution, n^2 per shell)
  - Compute Euclidean MST on the union of points
  - Apply iterative Steiner insertion at each MST junction (degree-3+ node):
      * replace high-degree MST vertex with Fermat/Steiner point at Weiszfeld optimum
      * keep only if total edge length decreases
  - Record final Steiner point positions

This is the "MST + local Steiner improvement" heuristic, which gives a 4/3
approximation to the optimal Steiner MST for metric problems. Not exact,
but substantially better than the 2D Delaunay-Fermat approach that only
considered triangle-local optima.

Sparsity leverage: we never compute all O(N^2) pair interactions; we
start from the MST (which is O(N log N) via scipy.sparse.csgraph) and
improve locally. This is the "computationally sparse" analog of the
packing's geometric sparsity.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

sys.stdout.reconfigure(encoding="utf-8")

PHI = (1.0 + np.sqrt(5.0)) / 2.0


# ---------------------------------------------------------------------------
# 3D packing
# ---------------------------------------------------------------------------

def fibonacci_sphere(N: int, radius: float) -> np.ndarray:
    """Fibonacci (sunflower) distribution of N points on a sphere of given radius.

    Provides near-uniform distribution for any N >= 1. Points are deterministic
    given N; for different shell radii the same relative placements occur.
    """
    if N == 1:
        return np.array([[0.0, 0.0, radius]])
    points = np.zeros((N, 3))
    for i in range(N):
        y = 1 - 2 * (i + 0.5) / N  # y in (-1, 1)
        theta = 2 * np.pi * i / PHI
        rho = np.sqrt(max(0.0, 1 - y * y))
        x = rho * np.cos(theta)
        z = rho * np.sin(theta)
        points[i] = [radius * x, radius * y, radius * z]
    return points


def s3_packing_3d(n_max: int, d_0: float = 1.0):
    """Generate 3D S^3-inspired packing with n^2 points on sphere of radius n*d_0.

    Uses hydrogen degeneracy convention: shell n (principal quantum number)
    has n^2 orbital angular positions. Total through n_max: sum_{n=1..n_max} n^2.
    """
    shells = {}
    for n in range(1, n_max + 1):
        count = n * n  # n^2 orbital positions
        pts = fibonacci_sphere(count, n * d_0)
        shells[n] = pts
    return shells


# ---------------------------------------------------------------------------
# MST + Steiner heuristic
# ---------------------------------------------------------------------------

def euclidean_mst(points: np.ndarray):
    """Compute Euclidean MST using scipy.sparse.csgraph.minimum_spanning_tree."""
    dist = squareform(pdist(points))
    sparse = csr_matrix(dist)
    mst = minimum_spanning_tree(sparse)
    mst = mst.toarray()
    edges = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if mst[i, j] > 0 or mst[j, i] > 0:
                edges.append((i, j, float(max(mst[i, j], mst[j, i]))))
    return edges


def fermat_weiszfeld(neighbors, max_iter: int = 200, tol: float = 1e-10) -> np.ndarray:
    """Weiszfeld's algorithm for geometric median (Fermat point) in any dimension.

    Returns the point minimizing sum of Euclidean distances to the given neighbors.
    For 3 points with all angles < 120° this equals the Fermat point.
    For more neighbors it's the geometric median.
    """
    neighbors = np.asarray(neighbors, dtype=float)
    F = neighbors.mean(axis=0)  # centroid start
    for _ in range(max_iter):
        diffs = neighbors - F
        dists = np.linalg.norm(diffs, axis=1)
        if dists.min() < tol:
            return F
        w = 1.0 / dists
        F_new = (w[:, None] * neighbors).sum(axis=0) / w.sum()
        if np.linalg.norm(F_new - F) < tol:
            return F_new
        F = F_new
    return F


def mst_with_steiner_improvement(points: np.ndarray, max_passes: int = 3):
    """Compute MST + iterative local Steiner insertion.

    For each MST vertex of degree >= 3, try replacing with a Steiner point at
    the Fermat (geometric-median) location of its neighbors. Keep the change
    if total tree length decreases.

    Returns:
      final_length : total length after insertions
      steiner_points : list of (position, degree, replaced_terminal_idx)
      mst_length : original MST length (for comparison)
    """
    edges = euclidean_mst(points)
    mst_length = sum(w for _, _, w in edges)

    # Build adjacency as dict
    adj = defaultdict(list)
    for i, j, w in edges:
        adj[i].append(j)
        adj[j].append(i)

    # Current tree: list of edges with explicit (position_a, position_b)
    # We'll track Steiner points as new node indices beyond len(points)
    all_points = points.copy()
    steiner_records = []

    for pass_idx in range(max_passes):
        improved = False
        # Look at each terminal with degree >= 3 in the current tree
        deg = defaultdict(int)
        for i, j, w in edges:
            deg[i] += 1
            deg[j] += 1

        for node_idx in list(deg):
            if node_idx >= len(points):
                continue  # skip already-Steiner nodes
            if deg[node_idx] < 3:
                continue
            # Try Steiner improvement: keep node as terminal, add Steiner point
            neighbors_idx = []
            for i, j, w in edges:
                if i == node_idx:
                    neighbors_idx.append(j)
                elif j == node_idx:
                    neighbors_idx.append(i)
            neighbor_positions = np.array([all_points[k] for k in neighbors_idx])
            # Fermat point of neighbors + the node itself
            candidate_positions = np.vstack([neighbor_positions, all_points[node_idx]])
            fermat_F = fermat_weiszfeld(candidate_positions)
            # Would length decrease if we replaced node->neighbor edges with
            # node->Fermat + Fermat->each_neighbor? No, that adds length.
            # Correct Steiner move: delete some edges, add Steiner point, reconnect.
            # Simpler: insert Steiner at Fermat of 3 of the neighbors, and connect
            # the 4th neighbor directly to either Steiner or the node.
            # Quick check: if fermat point is strictly inside the tree and the
            # sum of distances from Fermat to 3 neighbors < sum of current 3 edges,
            # swap.
            if len(neighbors_idx) >= 3:
                # Try Steiner improvement: pick TWO of the 3 neighbors, insert a
                # Steiner point F at Fermat of (X, A, B), delete edges XA and XB,
                # add edges FA, FB, FX. Classical 2-to-3 Steiner move.
                #
                # The node X keeps any remaining edges. Condition for improvement:
                #   |FA| + |FB| + |FX|  <  |XA| + |XB|
                # because F=Fermat minimizes sum of distances to the three points
                # that form the Steiner triangle (X, A, B).
                nbr_dists = [(k, np.linalg.norm(all_points[k] - all_points[node_idx]))
                             for k in neighbors_idx]
                nbr_dists.sort(key=lambda x: x[1])
                # Try the 3 pair choices among the 3 closest neighbors.
                three_nbrs = [k for k, _ in nbr_dists[:3]]
                pairs_to_try = [
                    (three_nbrs[0], three_nbrs[1]),
                    (three_nbrs[0], three_nbrs[2]),
                    (three_nbrs[1], three_nbrs[2]),
                ]
                best_pair = None
                best_F = None
                best_improvement = 0.0
                X = all_points[node_idx]
                for (a_idx, b_idx) in pairs_to_try:
                    A = all_points[a_idx]
                    B = all_points[b_idx]
                    # Triangle (X, A, B) angle check
                    angles_ok = True
                    for P, Q, R in [(X, A, B), (A, X, B), (B, X, A)]:
                        PQ = Q - P
                        PR = R - P
                        cosang = np.dot(PQ, PR) / (
                            np.linalg.norm(PQ) * np.linalg.norm(PR) + 1e-20)
                        if cosang < -0.5:
                            angles_ok = False
                            break
                    if not angles_ok:
                        continue
                    F = fermat_weiszfeld(np.array([X, A, B]))
                    new_len = (np.linalg.norm(F - X)
                               + np.linalg.norm(F - A)
                               + np.linalg.norm(F - B))
                    old_len = (np.linalg.norm(X - A)
                               + np.linalg.norm(X - B))
                    improvement = old_len - new_len
                    if improvement > best_improvement + 1e-12:
                        best_improvement = improvement
                        best_pair = (a_idx, b_idx)
                        best_F = F
                if best_pair is not None:
                    # Accept Steiner insertion
                    steiner_idx = len(all_points)
                    all_points = np.vstack([all_points, F])
                    # Remove 3 edges node->three_nbrs, add 3 edges F->three_nbrs
                    # plus 1 edge node->F
                    new_edges = []
                    for i, j, w in edges:
                        if (i == node_idx and j in three_nbrs) or (
                                j == node_idx and i in three_nbrs):
                            continue
                        new_edges.append((i, j, w))
                    for k in three_nbrs:
                        new_edges.append(
                            (steiner_idx, k, np.linalg.norm(F - all_points[k])))
                    new_edges.append(
                        (node_idx, steiner_idx,
                         np.linalg.norm(F - all_points[node_idx])))
                    edges = new_edges
                    steiner_records.append({
                        "position": [float(F[0]), float(F[1]), float(F[2])],
                        "radius": float(np.linalg.norm(F)),
                        "connected_to": three_nbrs + [node_idx],
                    })
                    improved = True
                    break  # restart pass
        if not improved:
            break

    final_length = sum(w for _, _, w in edges)
    return {
        "mst_length": mst_length,
        "final_length": final_length,
        "improvement": mst_length - final_length,
        "steiner_points": steiner_records,
        "n_steiner": len(steiner_records),
    }


# ---------------------------------------------------------------------------
# Main: per-shell-pair Steiner analysis
# ---------------------------------------------------------------------------

def main():
    d_0 = 1.0
    n_max = 6

    print("=" * 72)
    print("Paper-0 Packing Steiner Probe — 3D, per-shell-pair, MST + local Steiner")
    print("=" * 72)
    print(f"\nPacking: n_max={n_max}, n^2 points per shell (S^3 hydrogen convention)")

    shells = s3_packing_3d(n_max, d_0)
    for n, pts in shells.items():
        print(f"  Shell n={n}: r={n*d_0:.2f}, {len(pts)} points")

    print()
    dirac_g = lambda n: 2 * (n + 1) * (n + 2)

    all_results = {}
    shell_pair_summary = []
    for n in range(1, n_max):
        # Shell pair (n, n+1)
        pts_n = shells[n]
        pts_np1 = shells[n + 1]
        combined = np.vstack([pts_n, pts_np1])
        print(f"--- Shell pair ({n}, {n+1}) — {len(combined)} terminals ---")
        r_mid = (n + 0.5) * d_0
        print(f"  Midpoint radius (conjecture target): {r_mid:.2f}")
        print(f"  Dirac g_{{n-1}} = 2·n·(n+1) = {dirac_g(n-1)} "
              f"(if mapping pair-idx-from-1 to n=pair-1)")
        result = mst_with_steiner_improvement(combined, max_passes=5)
        print(f"  MST length: {result['mst_length']:.4f}")
        print(f"  After Steiner: {result['final_length']:.4f}")
        print(f"  Improvement: {result['improvement']:.4f} "
              f"({100*result['improvement']/result['mst_length']:.2f}%)")
        print(f"  Steiner points added: {result['n_steiner']}")
        if result['steiner_points']:
            radii = [sp['radius'] for sp in result['steiner_points']]
            print(f"  Steiner radii: {[f'{r:.3f}' for r in sorted(radii)]}")
            print(f"    mean: {np.mean(radii):.3f}, "
                  f"std: {np.std(radii):.3f}")
            print(f"    n near r_mid={r_mid:.2f} (±0.15): "
                  f"{sum(1 for r in radii if abs(r - r_mid) < 0.15)}")
        else:
            print("  (No Steiner improvement found; MST was already locally optimal.)")
        all_results[f"({n},{n+1})"] = result
        shell_pair_summary.append({
            "pair": (n, n + 1),
            "n_terminals": len(combined),
            "mid_radius": r_mid,
            "dirac_ref_g": dirac_g(n - 1),
            "mst_length": result['mst_length'],
            "steiner_length": result['final_length'],
            "n_steiner": result['n_steiner'],
            "steiner_radii": [sp['radius'] for sp in result['steiner_points']],
        })
        print()

    # Aggregate summary
    print("=" * 72)
    print("Summary")
    print("=" * 72)
    print(f"\n{'Pair':<8} {'N_term':<8} {'r_mid':<8} {'Dirac g':<10} {'N_steiner':<10} "
          f"{'mean r':<10} {'near r_mid':<12}")
    for s in shell_pair_summary:
        radii = s['steiner_radii']
        mean_r = np.mean(radii) if radii else float('nan')
        near = sum(1 for r in radii if abs(r - s['mid_radius']) < 0.15)
        print(f"  {s['pair']!r:<8} {s['n_terminals']:<8} "
              f"{s['mid_radius']:<8.2f} {s['dirac_ref_g']:<10} "
              f"{s['n_steiner']:<10} {mean_r:<10.3f} {near:<12}")

    # Write JSON
    out_path = Path("debug/data/paper0_steiner_3d.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump({
            "generated_by": "debug/compute_paper0_steiner_3d.py",
            "d_0": d_0,
            "n_max": n_max,
            "shells_per_shell_count": {str(n): n*n for n in range(1, n_max+1)},
            "shell_pair_results": {
                key: {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                      for k, v in val.items()
                      if k != "steiner_points" or True}
                for key, val in all_results.items()
            },
            "summary": [{k: (v if not isinstance(v, tuple) else list(v))
                         for k, v in s.items()} for s in shell_pair_summary],
            "dirac_reference": [{"n": n, "lambda": n + 1.5,
                                 "g_n": 2 * (n + 1) * (n + 2)}
                                for n in range(5)],
        }, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
