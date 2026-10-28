"""KS phase-1 feasibility check (draft record docs/predictions/2026-10-28_kite_speed_test.md). NOT a timing.

Does a plain float64 tensor Gauss rule on the route-B double integral of pf_common reproduce the 60-digit grid's
phi at a few nodes to the 1e-9 .. 1e-12 band?  u: Gauss-Chebyshev 2nd kind (weight sqrt(u(1-u))), v: Gauss-Legendre
split at the ridge v = 1 - u (where k2 = (1-u-v)^2 z vanishes; at large z the kernel peaks there).
Prints relative errors only.  No clock is read."""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

DATA = Path(__file__).parent / "data"
A, B, G1, G2 = 1.0, 2.0, 1.5, 3.0


def kernel(M, N, k2):
    N2 = N * N
    s = M + N
    s2 = s * s
    iD = 1.0 / (s2 + k2)
    iD2 = iD * iD
    return 3 * iD + 6 * s2 * iD2 + 2 * (N2 + 3 * M * N) * (3 * s2 - k2) * iD2 * iD \
        + 24 * M * N2 * s * (s2 - k2) * iD2 * iD2


def cheb2(n):
    k = np.arange(1, n + 1)
    th = k * np.pi / (n + 1)
    x = np.cos(th)                      # on (-1, 1), weight sqrt(1-x^2)
    w = np.pi / (n + 1) * np.sin(th) ** 2
    # map to u in (0,1): sqrt(u(1-u)) du = (1/4) sqrt(1-x^2) dx / 2 ... u = (1+x)/2, du = dx/2, sqrt(u(1-u)) = sqrt(1-x^2)/2
    return (1 + x) / 2, w / 4


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return lo[..., None] + (hi - lo)[..., None] * (x + 1) / 2, (hi - lo)[..., None] * w / 2


def phi64(w, nu, nv):
    z = w * w
    u, wu = cheb2(nu)
    M = np.sqrt(A * A * u + B * B * (1 - u) + u * (1 - u) * z)
    cu = wu / M ** 3
    tot = 0.0
    for lo, hi in ((np.zeros_like(u), 1 - u), (1 - u, np.ones_like(u))):
        v, wv = gl(nv, lo, hi)
        V = v * (1 - v)
        N = np.sqrt(G1 * G1 * (1 - v) + G2 * G2 * v + V * z)
        cv = wv * V / N ** 5
        d = 1 - u[:, None] - v
        K = kernel(M[:, None], N, d * d * z)
        tot += np.sum(cu[:, None] * cv * K)
    return 2 * G1 * G2 * A ** 3 * B ** 3 * tot


def main():
    rec = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    ws = [float(x) for x in rec["w"]]
    ph = [float(x) for x in rec["phi"]]
    picks = [i for i, w in enumerate(ws) if w > 1e-3][::20] + [len(ws) - 1]
    for nu, nv in ((40, 40), (80, 80), (160, 160)):
        errs = []
        for i in picks:
            errs.append((ws[i], abs(phi64(ws[i], nu, nv) / ph[i] - 1)))
        print(nu, nv, " ".join(f"w={w:.2e}:{e:.1e}" for w, e in errs), flush=True)


if __name__ == "__main__":
    main()
