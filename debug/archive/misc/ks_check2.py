"""KS phase-1 feasibility check 2 (draft record docs/predictions/2026-10-28_kite_speed_test.md). NOT a timing.
float64 tanh-sinh in u (both ends graded, u and 1-u kept separately) and in v on the two panels split at the ridge
v = 1-u.  Relative error of phi against the 60-digit grid; and the Y-weighted error budget per node at t = 0.8, 2.5.
No clock is read."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ks_check import kernel, A, B, G1, G2

DATA = Path(__file__).parent / "data"


def ts01(h, kmax):
    k = np.arange(-kmax, kmax + 1)
    t = k * h
    s = np.pi / 2 * np.sinh(t)
    # x = tanh(s) on (-1,1); u = (1+x)/2 = 1/(1+e^{-2s}), 1-u = 1/(1+e^{2s})
    u = 1 / (1 + np.exp(-2 * s))
    um = 1 / (1 + np.exp(2 * s))
    w = h * np.pi / 2 * np.cosh(t) / np.cosh(s) ** 2 / 2
    keep = (u > 0) & (um > 0) & (w > 1e-300)
    return u[keep], um[keep], w[keep]


def phi_ts64(wn, h, kmax):
    z = wn * wn
    u, um, wu = ts01(h, kmax)
    U = u * um
    M = np.sqrt(A * A * u + B * B * um + U * z)
    cu = wu * np.sqrt(U) / M ** 3
    r, rm, wr = ts01(h, kmax)
    tot = 0.0
    for side in (0, 1):
        if side == 0:      # v in (0, 1-u): v = (1-u) r
            L = um[:, None]
            v = L * r
            vm = 1 - v
            d = L * rm      # 1-u-v = (1-u)(1-r)
        else:              # v in (1-u, 1): v = 1 - u*(1-r) ... 1-v = u*rm
            L = u[:, None]
            vm = L * rm
            v = 1 - vm
            d = -L * r      # 1-u-v = u - u r ... = -(v-(1-u)) = -(u r)
        wv = L * wr
        V = v * vm
        N = np.sqrt(G1 * G1 * vm + G2 * G2 * v + V * z)
        cv = wv * V / N ** 5
        K = kernel(M[:, None], N, d * d * z)
        tot += np.sum(cu[:, None] * cv * K)
    return 2 * G1 * G2 * A ** 3 * B ** 3 * tot


def main():
    rec = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    ws = np.array([float(x) for x in rec["w"]])
    dws = np.array([float(x) for x in rec["dw"]])
    ph = np.array([float(x) for x in rec["phi"]])
    Y = {}
    for t in (0.8, 2.5):
        Y[t] = np.sum(ws ** 2 * ph * dws * (1 + ws ** 2 / t ** 2) ** -1.5)
    for h, kmax in ((1 / 8, 40), (1 / 16, 80), (1 / 32, 160)):
        worst_rel, worst_w = {}, {}
        line = []
        for i in range(0, len(ws), 6):
            e = abs(phi_ts64(ws[i], h, kmax) / ph[i] - 1)
            for t in Y:
                c = e * ws[i] ** 2 * ph[i] * dws[i] * (1 + ws[i] ** 2 / t ** 2) ** -1.5 / abs(Y[t])
                if c > worst_rel.get(t, 0):
                    worst_rel[t], worst_w[t] = c, ws[i]
            if i % 30 == 0:
                line.append(f"w={ws[i]:.1e}:{e:.0e}")
        print(f"h={h} kmax={kmax}", " ".join(line), flush=True)
        print("   worst node contribution to Y error (rel to Y):",
              {t: (f"{worst_rel[t]:.1e}", f"w={worst_w[t]:.1e}") for t in Y}, flush=True)


if __name__ == "__main__":
    main()
