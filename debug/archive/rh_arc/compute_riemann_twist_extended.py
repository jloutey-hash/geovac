"""
RH-G extended: the chi_{-4} twist is structurally null because all GeoVac
Hopf graphs are bipartite (tr(T^L) = 0 for all odd L). Bipartiteness is
mod-2; chi_{-4} is mod-4. The two parities disagree on support.

Question: is there a natural character that DOES align with graph walks?

Strategy:
  (a) Reparametrize L = 2m (walks only on even L). Then the "quarter-length"
      character chi_{-4}(m) becomes the odd-integer test on m.
  (b) Try chi_{-3}(n) = Legendre symbol (n|3), mod-3 non-principal.
      chi_{-3}(n) = 1 if n = 1 mod 3, -1 if n = 2 mod 3, 0 if n = 0 mod 3.
  (c) Try the sign (-1)^L (trivial once we restrict to even L) and
      (-1)^(L/2) (which gives +1,-1,+1,-1,... on L = 2,4,6,8,...).

The goal is to find a character that gives non-trivial sums.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

from geovac.ihara_zeta import (
    _count_closed_nonbacktracking_walks,
    _mobius_to_primitive,
    hashimoto_matrix,
)
from geovac.lattice import GeometricLattice
from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


def chi_minus3(n: int) -> int:
    m = n % 3
    return 1 if m == 1 else (-1 if m == 2 else 0)


def chi_sign_half(L: int) -> int:
    """(-1)^(L/2) when L even, else 0."""
    if L % 2 != 0:
        return 0
    return (-1) ** (L // 2)


def chi_mod4_on_half(L: int) -> int:
    """chi_{-4}(L/2) when L even, else 0.
    On L = 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 this is chi_{-4}(m)
    for m = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 = +1, 0, -1, 0, +1, 0, -1, 0, +1, 0
    """
    if L % 2 != 0:
        return 0
    m = L // 2
    r = m % 4
    return 1 if r == 1 else (-1 if r == 3 else 0)


def A_s3(max_n):
    lat = GeometricLattice(max_n=max_n)
    A_raw = lat.adjacency
    if hasattr(A_raw, "toarray"):
        A = np.asarray(A_raw.toarray(), dtype=int)
    else:
        A = np.asarray(A_raw, dtype=int)
    A = (np.abs(A) > 0).astype(int); np.fill_diagonal(A, 0); return A


def A_s5(N_max):
    bg = build_bargmann_graph(N_max=N_max)
    Aw = np.asarray(bg.adjacency_dense())
    A = (np.abs(Aw) > 1e-12).astype(int); np.fill_diagonal(A, 0); return A


def A_dir(n_max, rule):
    A, _, _, _ = build_dirac_s3_graph(n_max=n_max, adjacency_rule=rule)
    return np.asarray(A, dtype=int)


def investigate(A, name, L_MAX=20):
    print(f"\n### {name}  V={A.shape[0]}  E={int(A.sum())//2}")
    tr = _count_closed_nonbacktracking_walks(A, L_MAX)

    # Build trace table
    print(f"{'L':>4}{'tr(T^L)':>16}{'chi_-3':>8}{'chi_s':>8}{'chi_m4_half':>12}"
          f"{'chi_-3*tr':>18}{'chi_s*tr':>14}{'chi_m4_h*tr':>14}")
    sums = {'chi_-3': 0, 'chi_s': 0, 'chi_m4_half': 0}
    for L in range(1, L_MAX + 1):
        c3 = chi_minus3(L)
        cs = chi_sign_half(L)
        cm4 = chi_mod4_on_half(L)
        trL = int(tr[L])
        sums['chi_-3'] += c3 * trL
        sums['chi_s'] += cs * trL
        sums['chi_m4_half'] += cm4 * trL
        print(f"{L:>4d}{trL:>16d}{c3:>8d}{cs:>8d}{cm4:>12d}"
              f"{c3 * trL:>18d}{cs * trL:>14d}{cm4 * trL:>14d}")

    # Partial sums for Dirichlet L-series on graph side
    print(f"Partial sums up to L={L_MAX}: {sums}")
    return sums


def main():
    graphs = [
        ("S3 Coulomb max_n=3", A_s3(3)),
        ("S5 Bargmann-Segal N_max=2", A_s5(2)),
        ("S5 Bargmann-Segal N_max=3", A_s5(3)),
        ("Dirac-S3 rule A n_max=3", A_dir(3, 'A')),
        ("Dirac-S3 rule B n_max=3", A_dir(3, 'B')),
    ]
    all_sums = {}
    for name, A in graphs:
        all_sums[name] = investigate(A, name, L_MAX=20)

    # Save to JSON
    out = Path(__file__).parent / "data" / "riemann_twist_extended.json"
    with out.open("w") as f:
        json.dump(all_sums, f, indent=2, default=str)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
