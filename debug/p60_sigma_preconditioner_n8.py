"""Follow-on to debug/p60_sigma_preconditioner_diag.py (2026-10-09):
 (a) canonical-orthogonalisation accuracy cost at N = 8 (36 functions/centre), H2+ R = 2;
 (b) quadrature-floor check of lambda_min(S) at the cond ~1e7..1e15 points (160x160 vs 176x320);
 (c) involution test (C Pi)^2 -> I on the 1s,2s,2p block out to N = 9.
Output: debug/data/p60_sigma_preconditioner_n8.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geovac.sturmian_sigma_basis as sb  # noqa: E402

H2P_R2 = -1.102634214


def top(W, S, tol):
    ev, V = np.linalg.eigh(S)
    keep = ev > tol * ev[-1]
    X = V[:, keep] / np.sqrt(ev[keep])
    return float(np.linalg.eigvalsh(X.T @ (0.5 * (W + W.T)) @ X)[-1]), int((~keep).sum())


def swap_quadrature(nx, ne):
    sb.NX, sb.NE = nx, ne
    sb._TX, sb._WX = np.polynomial.laguerre.laggauss(nx)
    sb._EE, sb._WE = np.polynomial.legendre.leggauss(ne)


def main():
    out = {"canon_N8": [], "floor": [], "involution": []}
    for tol in (1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10):
        k = brentq(lambda kk: top(*sb.build(8, kk, 2.0)[1::-1], tol)[0] - kk, 0.9, 2.5, xtol=1e-12)
        S, W, lab = sb.build(8, k, 2.0)
        _, nd = top(W, S, tol)
        out["canon_N8"].append({"tol": tol, "err": -0.5 * k * k - H2P_R2, "n_dropped": nd,
                                "fns_total": 2 * len(lab)})
        print("N8", tol, -0.5 * k * k - H2P_R2, nd, flush=True)
    for N in range(2, 10):
        S, _, lab = sb.build(N, 1.0, 3.0)
        m = len(lab)
        CP = S[:m, m:] @ np.diag([(-1.0) ** l for (_, l) in lab])
        out["involution"].append({"N": N, "dev": float(np.abs((CP @ CP)[:3, :3] - np.eye(3)).max())})
    pts = [(1.0, 3), (1.0, 4), (3.0, 4), (3.0, 5), (3.0, 6), (8.0, 6), (8.0, 7)]
    base = {p: np.linalg.eigvalsh(sb.build(p[1], 1.0, p[0])[0])[:3] for p in pts}
    swap_quadrature(176, 320)
    for p in pts:
        e = np.linalg.eigvalsh(sb.build(p[1], 1.0, p[0])[0])[:3]
        out["floor"].append({"kR": p[0], "N": p[1], "lam_low3_160": base[p].tolist(),
                             "lam_low3_176x320": e.tolist()})
    swap_quadrature(160, 160)
    with open("debug/data/p60_sigma_preconditioner_n8.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
