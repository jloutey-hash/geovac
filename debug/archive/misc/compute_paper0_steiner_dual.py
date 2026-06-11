"""Paper-0 packing Steiner / Fermat / Voronoi dual probe.

Tests the PI's conjecture: does the Steiner-tree / Fermat-point dual of
the Paper-0 packing (integer shells r_k = k*d_0) land at half-integer
radii (k+1/2)*d_0, and does its count match the Dirac degeneracy
g_n = 2(n+1)(n+2)?

If yes: 2D geometric derivation of the Dirac spectrum's n+3/2 shift.
If no: characterize the actual pattern.

Three independent computations:
  1. Delaunay-triangle Fermat points (3-terminal Steiner minima of
     natural nearest-neighbor triples from the packing).
  2. Voronoi vertices of the packing (3+ Voronoi cells meeting).
  3. Midpoint-bisector radii (the radial-only Voronoi boundary between
     shell k and shell k+1 lives exactly at (k+1/2)*d_0 by construction;
     this is the trivial sanity check).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.spatial import Voronoi, Delaunay

sys.stdout.reconfigure(encoding="utf-8")


# ---------------------------------------------------------------------------
# Paper-0 packing generator
# ---------------------------------------------------------------------------

def paper0_packing(k_max: int, d_0: float = 1.0):
    """Generate Paper-0 packing.

    Shell k at radius r_k = k*d_0.
    Shell k has 2k-1 distinct angular positions (spin doubling collapsed;
    spin is a Z/2 label, not a spatial coordinate).
    Shell 1 initialized with N_init = 2 points (Axiom 1). For k=1 the
    angular count 2*1-1 = 1 is doubled to 2 by the Z/2 initialization
    convention, and we place the 2 points at angles 0 and pi (chord = 2*d_0,
    but the initialization is structurally distinct from later shells).

    Returns:
      points : (V, 2) array of 2D positions
      labels : list of (k, i) tuples
    """
    points: list[list[float]] = []
    labels: list[tuple[int, int]] = []
    for k in range(1, k_max + 1):
        r_k = k * d_0
        if k == 1:
            # N_init = 2: angles 0, pi.
            n_angular = 2
            for i in range(n_angular):
                theta = np.pi * i
                points.append([r_k * np.cos(theta), r_k * np.sin(theta)])
                labels.append((k, i))
        else:
            n_angular = 2 * k - 1  # 2*ell + 1 per Paper 0 Eq. 5
            for i in range(n_angular):
                theta = 2 * np.pi * i / n_angular
                points.append([r_k * np.cos(theta), r_k * np.sin(theta)])
                labels.append((k, i))
    return np.array(points), labels


# ---------------------------------------------------------------------------
# Fermat-point computation (Weiszfeld)
# ---------------------------------------------------------------------------

def fermat_point(P1, P2, P3, max_iter: int = 500, tol: float = 1e-13):
    """Fermat (1-median) point of a triangle.

    If any triangle angle >= 120 degrees, the Fermat point coincides with
    that vertex; returns (vertex, False). Otherwise returns (interior
    Fermat point, True) via Weiszfeld's iteration.
    """
    P1, P2, P3 = map(np.asarray, (P1, P2, P3))
    # Check for >= 120 degree angles
    for A, B, C in [(P2, P1, P3), (P1, P2, P3), (P1, P3, P2)]:
        # Angle at B between edges B->A and B->C
        BA = A - B
        BC = C - B
        cosang = np.dot(BA, BC) / (np.linalg.norm(BA) * np.linalg.norm(BC) + 1e-20)
        if cosang < -0.5:  # angle > 120°
            return B, False
    # Interior Fermat point via Weiszfeld
    F = (P1 + P2 + P3) / 3.0
    for _ in range(max_iter):
        d1, d2, d3 = (np.linalg.norm(F - P) for P in (P1, P2, P3))
        if min(d1, d2, d3) < tol:
            return F, True
        w = np.array([1/d1, 1/d2, 1/d3])
        F_new = (w[0]*P1 + w[1]*P2 + w[2]*P3) / w.sum()
        if np.linalg.norm(F_new - F) < tol:
            return F_new, True
        F = F_new
    return F, True


# ---------------------------------------------------------------------------
# Analysis: Delaunay-triangle Fermat dual
# ---------------------------------------------------------------------------

def delaunay_fermat_dual(points, d_0: float = 1.0):
    """For each Delaunay triangle of the packing, compute its Fermat point.

    Returns list of (r, interior_bool, triangle_indices) for all triangles.
    """
    tri = Delaunay(points)
    results = []
    for simplex in tri.simplices:
        P1, P2, P3 = points[simplex[0]], points[simplex[1]], points[simplex[2]]
        F, is_interior = fermat_point(P1, P2, P3)
        r = np.linalg.norm(F) / d_0
        results.append({
            "triangle": [int(s) for s in simplex],
            "fermat_x": float(F[0]),
            "fermat_y": float(F[1]),
            "radius_over_d0": float(r),
            "is_interior": bool(is_interior),
        })
    return results


# ---------------------------------------------------------------------------
# Analysis: Voronoi vertices
# ---------------------------------------------------------------------------

def voronoi_vertices_analysis(points, d_0: float = 1.0, max_radius: float = 10.0):
    """Compute Voronoi vertices of the packing and tabulate their radii.

    Filters out Voronoi vertices at infinity and those outside the ball
    of radius max_radius*d_0.
    """
    vor = Voronoi(points)
    vertices = vor.vertices  # (n_vertices, 2)
    results = []
    for i, v in enumerate(vertices):
        r = np.linalg.norm(v) / d_0
        if r < max_radius and np.isfinite(r):
            results.append({
                "index": int(i),
                "x": float(v[0]),
                "y": float(v[1]),
                "radius_over_d0": float(r),
            })
    return results


# ---------------------------------------------------------------------------
# Histogram helpers
# ---------------------------------------------------------------------------

def radii_histogram(radii, d_0: float = 1.0, bin_width: float = 0.1, r_max: float = 6.0):
    """Bin radii and return (bin_center, count) list."""
    bins = np.arange(0, r_max + bin_width, bin_width)
    hist, edges = np.histogram(radii, bins=bins)
    return [{"bin_low": float(edges[i]), "bin_high": float(edges[i+1]), "count": int(hist[i])}
            for i in range(len(hist)) if hist[i] > 0]


def cluster_near_half_integers(radii, d_0: float = 1.0, width: float = 0.1):
    """For each half-integer k+0.5, count radii within +/- width."""
    cluster = {}
    for k in range(1, 8):
        target = k + 0.5
        count = sum(1 for r in radii if abs(r - target) < width)
        cluster[f"r={target}"] = count
    # Also count at integer radii (baseline)
    for k in range(1, 8):
        target = float(k)
        count = sum(1 for r in radii if abs(r - target) < width)
        cluster[f"r={target} (integer)"] = count
    return cluster


# ---------------------------------------------------------------------------
# Dirac degeneracy reference
# ---------------------------------------------------------------------------

def dirac_degeneracies(n_max: int = 6):
    """Camporesi-Higuchi Dirac degeneracies on S^3: g_n = 2(n+1)(n+2)."""
    return [{"n": n, "lambda": n + 1.5, "g_n": 2 * (n + 1) * (n + 2)}
            for n in range(n_max + 1)]


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    d_0 = 1.0
    k_max = 6

    print("=" * 72)
    print("Paper-0 Packing Steiner / Fermat / Voronoi Dual Probe")
    print("=" * 72)

    # Generate packing
    points, labels = paper0_packing(k_max, d_0)
    print(f"\nPacking: k_max={k_max}, {len(points)} total points")
    per_shell = defaultdict(int)
    for k, i in labels:
        per_shell[k] += 1
    for k in sorted(per_shell):
        expected = 2 if k == 1 else (2 * k - 1)
        print(f"  Shell k={k}: r={k*d_0:.2f}, N={per_shell[k]} points "
              f"(expected {expected})")

    # Delaunay Fermat dual
    print(f"\n--- Delaunay-triangle Fermat points ---")
    delaunay_results = delaunay_fermat_dual(points, d_0)
    interior = [r["radius_over_d0"] for r in delaunay_results if r["is_interior"]]
    vertex = [r["radius_over_d0"] for r in delaunay_results if not r["is_interior"]]
    print(f"  Total Delaunay triangles: {len(delaunay_results)}")
    print(f"  Interior Fermat points (all angles < 120°): {len(interior)}")
    print(f"  Degenerate (>=120° triangles, F = vertex): {len(vertex)}")

    print("\n  Interior Fermat radii (sorted):")
    for r in sorted(set(round(x, 4) for x in interior))[:40]:
        count = sum(1 for x in interior if abs(x - r) < 1e-3)
        print(f"    r/d_0 = {r:.4f}  (count = {count})")

    print("\n  Histogram of interior Fermat radii:")
    hist_i = radii_histogram(interior, d_0, bin_width=0.1)
    for b in hist_i:
        bar = "#" * b["count"]
        print(f"    [{b['bin_low']:.2f}, {b['bin_high']:.2f}): {b['count']:>3d} {bar}")

    print("\n  Clustering near half-integers:")
    clusters = cluster_near_half_integers(interior, d_0, width=0.1)
    for key, count in clusters.items():
        print(f"    {key}: {count}")

    # Voronoi vertices
    print(f"\n--- Voronoi vertices ---")
    vor_results = voronoi_vertices_analysis(points, d_0)
    vor_radii = [v["radius_over_d0"] for v in vor_results]
    print(f"  Total Voronoi vertices (inside r < 10): {len(vor_radii)}")

    print("\n  Histogram of Voronoi vertex radii:")
    hist_v = radii_histogram(vor_radii, d_0, bin_width=0.1)
    for b in hist_v[:40]:
        bar = "#" * b["count"]
        print(f"    [{b['bin_low']:.2f}, {b['bin_high']:.2f}): {b['count']:>3d} {bar}")

    print("\n  Voronoi-vertex clustering near half-integers:")
    vor_clusters = cluster_near_half_integers(vor_radii, d_0, width=0.1)
    for key, count in vor_clusters.items():
        print(f"    {key}: {count}")

    # Compare to Dirac degeneracies
    print("\n--- Dirac degeneracy reference (Camporesi-Higuchi) ---")
    dirac = dirac_degeneracies(k_max)
    print("  |lambda_n| = n + 3/2, g_n = 2(n+1)(n+2)")
    for d in dirac:
        print(f"    n={d['n']}: |lambda|={d['lambda']:.1f}, g_n={d['g_n']}")

    # Save everything to JSON
    out = {
        "generated_by": "debug/compute_paper0_steiner_dual.py",
        "d_0": d_0,
        "k_max": k_max,
        "packing_points": [{"shell": k, "idx": i, "x": float(p[0]), "y": float(p[1])}
                           for (k, i), p in zip(labels, points)],
        "delaunay_fermat": delaunay_results,
        "voronoi_vertices": vor_results,
        "interior_fermat_histogram": hist_i,
        "voronoi_vertex_histogram": hist_v,
        "interior_fermat_half_int_clustering": clusters,
        "voronoi_vertex_half_int_clustering": vor_clusters,
        "dirac_reference": dirac,
    }
    out_path = Path("debug/data/paper0_steiner_dual.json")
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
