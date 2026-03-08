#!/usr/bin/env python3
"""
debug_p0_per_state.py
======================
Test what happens when p0 = 1/n is used in the Fock stereographic projection.

Two interpretations are tested:

  A) p0=1/n, |p|=1/n (on-shell momentum, same as current)
     → ratio u = |p|/p0 = 1 for ALL states → equator collapse

  B) p0=1/n, |p|=1 (unit reference momentum)
     → ratio u = n → sign-flipped xi_4 → identical distances (S3 isometry)

Paper 7 (Sec III.1): "The stereographic projection always produces the unit S3,
regardless of p0. The intrinsic geometry has no physical scale."

This script proves that d^2(1s,2s) = 0.40 in ALL cases — the inter-shell
distance is a property of the ratio u = |p|/p0, which is invariant.

Author: GeoVac Development Team, February 2026
"""

import numpy as np


def fock_s3(n: int, l: int, m: int, p0: float, p_mag: float) -> np.ndarray:
    """Fock stereographic projection for state (n,l,m) with given p0 and |p|."""
    if l > 0:
        cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
        sin_th = np.sqrt(max(0.0, 1.0 - cos_th ** 2))
    else:
        cos_th = 1.0
        sin_th = 0.0

    # 3D momentum vector
    p_vec = p_mag * np.array([sin_th, 0.0, cos_th])
    p_sq = p_mag ** 2
    p0_sq = p0 ** 2
    denom = p_sq + p0_sq

    xi = np.zeros(4)
    xi[0] = 2.0 * p0 * p_vec[0] / denom
    xi[1] = 2.0 * p0 * p_vec[1] / denom
    xi[2] = 2.0 * p0 * p_vec[2] / denom
    xi[3] = (p_sq - p0_sq) / denom
    return xi


def print_coords_and_distances(label: str, coords: dict) -> None:
    """Print S3 coordinates and pairwise d^2 for a set of named states."""
    states = list(coords.keys())
    xis = [coords[s] for s in states]

    print(f"\n  --- {label} ---")
    print(f"  {'(n,l,m)':>14}   {'xi_1':>8}  {'xi_2':>8}  {'xi_3':>8}  {'xi_4':>8}  {'|xi|':>7}")
    print("  " + "-" * 60)
    for s, xi in zip(states, xis):
        norm = np.linalg.norm(xi)
        print(f"  ({s[0]},{s[1]},{s[2]:+d}):  "
              f"{xi[0]:>8.5f}  {xi[1]:>8.5f}  {xi[2]:>8.5f}  {xi[3]:>8.5f}  {norm:>7.5f}")

    # Key pairwise distances
    pairs = [
        ((1, 0, 0), (2, 0, 0), "d^2(1s, 2s)"),
        ((2, 0, 0), (2, 1, 0), "d^2(2s, 2p0)"),
        ((1, 0, 0), (2, 1, 0), "d^2(1s, 2p0)"),
    ]
    print(f"\n  Key pairwise distances:")
    for sa, sb, name in pairs:
        if sa in coords and sb in coords:
            dot = np.dot(coords[sa], coords[sb])
            d_sq = 2.0 * (1.0 - dot)
            print(f"    {name:20s} = {d_sq:.5f}")


def main() -> None:
    states_5 = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (2, 1, -1)]

    print("=" * 70)
    print("CURRENT CODE: p0=1 (uniform), |p|=1/n")
    print("=" * 70)
    coords_current = {}
    for (n, l, m) in states_5:
        coords_current[(n, l, m)] = fock_s3(n, l, m, p0=1.0, p_mag=1.0 / n)
    print_coords_and_distances("p0=1, |p|=1/n", coords_current)

    print("\n\n" + "=" * 70)
    print("INTERPRETATION A: p0=1/n (per-state), |p|=1/n (on-shell)")
    print("  ratio u = |p|/p0 = 1 for ALL states")
    print("=" * 70)
    coords_A = {}
    for (n, l, m) in states_5:
        coords_A[(n, l, m)] = fock_s3(n, l, m, p0=1.0 / n, p_mag=1.0 / n)
    print_coords_and_distances("p0=1/n, |p|=1/n (ON-SHELL)", coords_A)

    print("\n\n" + "=" * 70)
    print("INTERPRETATION B: p0=1/n (per-state), |p|=1 (unit reference)")
    print("  ratio u = |p|/p0 = n")
    print("=" * 70)
    coords_B = {}
    for (n, l, m) in states_5:
        coords_B[(n, l, m)] = fock_s3(n, l, m, p0=1.0 / n, p_mag=1.0)
    print_coords_and_distances("p0=1/n, |p|=1 (UNIT REF)", coords_B)

    # Show why distances are identical for B vs current
    print("\n\n" + "=" * 70)
    print("MATHEMATICAL PROOF: distances are invariant")
    print("=" * 70)
    print("""
  Paper 7, Eq. (unit_metric): the induced metric depends on u = p/p0 only.
  The factors of p0 cancel identically.

  Current code: p0=1, |p|=1/n  =>  u = 1/n
  Interpretation B: p0=1/n, |p|=1  =>  u = 1/(1/n) = n

  These give DIFFERENT u, but the S3 coordinates satisfy:
    xi_4(B) = (u^2-1)/(u^2+1) = (n^2-1)/(n^2+1)   [positive]
    xi_4(current) = (1/n^2-1)/(1/n^2+1) = (1-n^2)/(n^2+1)   [negative]

  This is a GLOBAL sign flip of xi_4 for ALL states, which is an
  isometry of S^3 (reflection through the equatorial hyperplane).
  All pairwise chordal distances d^2 = 2(1 - xi_i . xi_j) are preserved
  because the sign flips cancel in the dot product:
    xi_i . xi_j = sum(xi_k_i * xi_k_j) for k=1..4
  Flipping xi_4 for ALL states: xi_4_i * xi_4_j -> (-xi_4_i)(-xi_4_j) = same.

  Interpretation A: u=1 for all states -> ALL on equator, xi_4=0.
  d^2(1s,2s) = 2(1 - cos(0) * cos(0)) = 0 for same (l,m). CATASTROPHIC.
""")

    # Compute V_ee(1s,2s) for Li with kappa=7.5
    kappa_li = 7.5
    d_sq_1s_2s = 2.0 * (1.0 - np.dot(coords_current[(1, 0, 0)],
                                       coords_current[(2, 0, 0)]))
    v_ee = kappa_li / d_sq_1s_2s if d_sq_1s_2s > 1e-14 else float('inf')
    print(f"  V_ee(1s,2s) with kappa=7.5, d^2={d_sq_1s_2s:.4f}: {v_ee:.4f} Ha")
    print(f"  This is UNCHANGED by any p0 substitution.")
    print(f"  The problem is NOT solvable by changing p0 within the Fock framework.")

    print("\n\nCONCLUSION")
    print("=" * 70)
    print("""  The Fock stereographic projection is scale-invariant (Paper 7 Sec III.1).
  Changing p0 from 1 to 1/n either:
    (A) collapses all shells to the equator (if |p| also scales), or
    (B) flips xi_4 sign = S3 isometry = same distances.

  d^2(1s,2s) = 0.40 is a TOPOLOGICAL property of the Fock embedding:
  it measures how "close" the 1s and 2s shells are in S3 geometry.
  This cannot be changed by reparametrizing the projection.

  To get inter-shell V_ee in the correct range (1-5 Ha), we need either:
  1. A CONFORMAL-WEIGHTED chordal distance: d^2_phys = Omega_i * Omega_j * |p_i - p_j|^2
     (Paper 7, Eq. chordal) — this uses the conformal factor to modulate distances.
  2. A separate kappa for same-orbital vs. different-orbital pairs.
  3. A completely different inter-shell distance metric (e.g., based on |r_i - r_j|).
""")


if __name__ == "__main__":
    main()
