"""Probe B follow-up -- the small-k (identity-shift) audit, done honestly.

As k -> 0 the momentum leaf tends to the overlap matrix, which in an orthonormal
basis is the IDENTITY, so the leaf operator becomes the number operator N and its
square is a c-number on a fixed particle-number sector.  The small-k 1-norm mass
of the factorized encoding is therefore partly spurious.  Splitting each leaf

    Chat_nu = c_nu N  +  Ctilde_nu ,   c_nu = tr(C_nu)/n ,  tr(Ctilde_nu) = 0

gives, on the N_e sector,

    sum_nu w_nu Chat_nu^2 = [sum_nu w_nu c_nu^2] N_e^2                (constant)
                          + 2 N_e sum_nu w_nu c_nu Ctilde_nu         (ONE one-body op)
                          + sum_nu w_nu Ctilde_nu^2                  (two-body remainder)

The middle term is a single one-body operator summed over nu WITH SIGNS, so it is
folded into h before taking a 1-norm (cancellation across nu is real and large).
The honest shifted 1-norm is therefore

    lambda_shift = sum_pq |h_pq + dh_pq|  +  sum_nu w_nu (||Ctilde||_*^2 + ||Stilde||_*^2)

which is what this script reports, alongside the pessimistic variant that
1-norms the cross term node-by-node (the ``lam2_shift`` column of the main run).
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

GRID = "g5"   # lambda is converged to <1e-4 relative here (see main run)


def shift_report(sysname: str, tgrid) -> str:
    orbs, nuclei = system(sysname)
    n = len(orbs)
    ne = NELEC[sysname]
    dP = max(max(abs(o1.z - o2.z) for o1 in orbs for o2 in orbs), 1e-6)
    S, h_raw = one_electron(orbs, nuclei, tgrid)
    X = lowdin_X(S)
    hm = X.T @ h_raw @ X
    bands = build_grid(*GRIDS[GRID][:3], dP, pad=GRIDS[GRID][3])
    om, C, Si = grid_densities(orbs, bands, tgrid, X)
    lam = lambdas(om, C, Si, ne)
    lam["std"] = float(np.abs(eri_from_nodes(om, C, Si)).sum())

    eye = np.eye(n)[None, :, :]
    cC = np.trace(C, axis1=1, axis2=2) / n
    cS = np.trace(Si, axis1=1, axis2=2) / n
    Ct = C - cC[:, None, None] * eye
    St = Si - cS[:, None, None] * eye
    ntC = np.abs(np.linalg.eigvalsh(Ct)).sum(axis=1)
    ntS = np.abs(np.linalg.eigvalsh(St)).sum(axis=1)

    lam2_tl = float(np.sum(om * (ntC ** 2 + ntS ** 2)))
    dh = 2 * ne * np.einsum("n,n,npq->pq", om, cC, Ct) \
        + 2 * ne * np.einsum("n,n,npq->pq", om, cS, St)
    lam1p = float(np.abs(hm + dh).sum())
    lam1 = float(np.abs(hm).sum())
    const = float(np.sum(om * (cC ** 2 + cS ** 2))) * ne ** 2

    tot_std = lam1 + lam["std"]
    tot_spec = lam1 + lam["spec"]
    tot_shift_pess = lam1 + lam["shift"]
    tot_shift = lam1p + lam2_tl
    return (f"{sysname:>7} {n:>4} {ne:>4}  {node_count(bands):>7}  "
            f"{tot_std:>10.5f} {tot_spec:>11.5f} {tot_shift_pess:>13.5f} "
            f"{tot_shift:>12.5f}   {tot_spec/tot_std:>7.3f} {tot_shift/tot_std:>9.3f}"
            f"   [lam1 {lam1:.4f} -> {lam1p:.4f}; lam2_tl {lam2_tl:.4f}; "
            f"absorbed const {const:.4f}]")


if __name__ == "__main__":
    tgrid = feynman_grid(y_max=26.0, panel=0.30, n_g=12)
    print("Probe B -- identity-shift (small-k) audit at grid " + GRID)
    print("system n_orb n_el  M nodes     lam_std   lam_spec  lam_shift_pess  "
          "lam_shift    spec/std shift/std")
    for s in ("H2", "LiH", "H2_4o", "LiH_5o"):
        print(shift_report(s, tgrid))
        sys.stdout.flush()
