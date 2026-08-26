"""Probe B control -- the momentum grid vs the STANDARD finite-rank factorizations.

Three exact factorizations of the SAME Loewdin-basis ERI tensor are compared on
one convention (the DF/SF spectral 1-norm, sum over leaves of the squared nuclear
norm, no 1/2, matching the v4.94.0 lambda_std convention on the other side):

  * sparse / standard    lambda_2 = sum_pqrs |(pq|rs)|
  * DF (eigen-Cholesky)  (pq|rs) = sum_m g_m U^m_pq U^m_rs  (n_orb^2 leaves)
                         lambda_2 = sum_m |g_m| (sum_i |eig U^m|)^2
  * momentum grid        lambda_2 = sum_nu w_nu (||C||_*^2 + ||S||_*^2)   [M leaves]

The point of the control: the momentum decomposition is a CONTINUUM-index single
factorization -- many leaves with infinitesimal weights -- whereas DF is a
rank-n_orb^2 one.  Squaring penalises concentration, so the refined decomposition
can have a smaller spectral 1-norm than the finite-rank one; this quantifies it.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from debug.probeB_thc_density import feynman_grid, lowdin_X, system  # noqa: E402
from debug.probeB_thc_lambda import (GRIDS, NELEC, build_grid,  # noqa: E402
                                     eri_from_nodes, grid_densities, lambdas,
                                     node_count, one_electron)

GRID = "g5"


def lambda_df(eri: np.ndarray) -> tuple:
    """Exact eigen-Cholesky (DF) spectral 1-norm of the ERI supermatrix."""
    n = eri.shape[0]
    W = eri.reshape(n * n, n * n)
    W = 0.5 * (W + W.T)
    g, U = np.linalg.eigh(W)
    lam = 0.0
    rank = 0
    for m in range(n * n):
        if abs(g[m]) < 1e-14:
            continue
        rank += 1
        Um = U[:, m].reshape(n, n)
        Um = 0.5 * (Um + Um.T)
        lam += abs(g[m]) * float(np.abs(np.linalg.eigvalsh(Um)).sum()) ** 2
    return lam, rank


if __name__ == "__main__":
    tgrid = feynman_grid(y_max=26.0, panel=0.30, n_g=12)
    print("Probe B control -- momentum grid vs finite-rank DF, at grid " + GRID)
    print("system n_orb  M nodes  DF rank   lam2_std   lam2_DF   lam2_momentum   "
          "DF/std  mom/std   mom/DF")
    for sname in ("H2", "LiH", "H2_4o", "LiH_5o"):
        orbs, nuclei = system(sname)
        n = len(orbs)
        dP = max(max(abs(o1.z - o2.z) for o1 in orbs for o2 in orbs), 1e-6)
        S, _h = one_electron(orbs, nuclei, tgrid)
        X = lowdin_X(S)
        bands = build_grid(*GRIDS[GRID][:3], dP, pad=GRIDS[GRID][3])
        om, C, Si = grid_densities(orbs, bands, tgrid, X)
        eri = eri_from_nodes(om, C, Si)
        l2_std = float(np.abs(eri).sum())
        l2_df, rank = lambda_df(eri)
        l2_mom = lambdas(om, C, Si, NELEC[sname])["spec"]
        print(f"{sname:>7} {n:>5}  {node_count(bands):>7}  {rank:>7}   "
              f"{l2_std:>8.5f}  {l2_df:>8.5f}  {l2_mom:>12.5f}   "
              f"{l2_df/l2_std:>6.3f}  {l2_mom/l2_std:>7.3f}  {l2_mom/l2_df:>7.3f}")
        sys.stdout.flush()
