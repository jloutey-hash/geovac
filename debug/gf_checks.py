"""GF controls (record 2026-10-25_gaussian_fit_geminal.md): (1) the unequal-exponent Gaussian loop against direct 3-D
quadrature at five triples; (2) the engine at A = B = 1, a = b = 1 against record GY's grid; (3) quadrature convergence
of Y_GG at two node sets; (4) the reference grid, 60 against 32 digits, at the registered t."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import gf_engine as F

DATA = Path(__file__).parent / "data"
out = {}


def loop_closed(A, B, k, L):
    mu = A * B / (A + B); c = A / (A + B)
    pre = (np.pi ** 2 / (A * B)) ** 1.5 * np.exp(-L @ L / (4 * (A + B)))
    R = np.array([np.linalg.norm(k - c * L)])
    # without the (pi/A)^(3/2)(pi/B)^(3/2) normalisation of g~: divide it out again
    return (pre * F.G_loop(R, mu)[0]) / (np.pi ** 2 / (A * B)) ** 1.5


def loop_direct(A, B, k, L, Nang=200, Nr=400):
    """int d^3q exp(-q^2/(4A)) exp(-(q-L)^2/(4B)) / |k - q|^2, spherical coordinates centred at k (1/s^2 cancels)."""
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, 2 * Nang, endpoint=False); wph = np.pi / Nang
    st = np.sqrt(1 - uu ** 2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    sr, wr = np.polynomial.legendre.leggauss(Nr)
    scale = 4 * np.sqrt(max(A, B)) + np.linalg.norm(k) + np.linalg.norm(L)
    tot = 0.0
    for lo, hi in [(0, scale / 4), (scale / 4, scale), (scale, 3 * scale)]:
        r = (hi - lo) / 2 * sr + (hi + lo) / 2; wg = (hi - lo) / 2 * wr
        for ri, wi in zip(r, wg):
            q = k + ri * n
            f = np.exp(-np.sum(q * q, -1) / (4 * A) - np.sum((q - L) ** 2, -1) / (4 * B))
            tot += wi * np.sum(wu[:, None] * f) * wph
    return tot


rng = np.random.default_rng(20261025)
rows = []
for _ in range(5):
    A, B = 10 ** rng.uniform(-0.5, 1.5, 2)
    k, L = rng.normal(size=3) * 1.5, rng.normal(size=3) * 1.5
    lc, ld = loop_closed(A, B, k, L), loop_direct(A, B, k, L)
    rows.append({"A": A, "B": B, "closed": lc, "direct": ld, "rel": abs(lc / ld - 1)})
    print("loop", f"A={A:.4g} B={B:.4g}", f"closed {lc:.15g} direct {ld:.15g} rel {abs(lc/ld-1):.2e}", flush=True)
out["loop"] = rows

sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_3_values as V
mp.mp.dps = 30
gy = json.loads((DATA / "gy_grid_g1_d32_h0.0078125.json").read_text())
ts = [0.8, 1.7, 2.6]
yE = F.Y_GG(1.0, 1.0, ts, a=1.0, b=1.0)
rows = []
for t, ye in zip(ts, yE):
    yg = float(V.Y_theta(gy, str(t), 0)[0])
    rows.append({"t": t, "engine": ye, "gy_grid": yg, "rel": abs(ye / yg - 1)})
    print("engine vs GY", t, f"{ye:.15g} {yg:.15g} rel {abs(ye/yg-1):.2e}", flush=True)
out["engine_vs_gy"] = rows

TS = [0.5, 1.0, 2.0, 4.0]
rows = []
for A, B in [(0.05 * 2.25, 0.05 * 9), (1.0 * 2.25, 1.0 * 9), (30 * 2.25, 30 * 9), (300 * 2.25, 300 * 9), (0.05 * 2.25, 300 * 9), (2216.25 * 2.25, 2216.25 * 9), (0.127672 * 2.25, 2216.25 * 9)]:
    y1 = F.Y_GG(A, B, TS)
    y2 = F.Y_GG(A, B, TS, nw=36, ns=36, nd=36)
    rel = np.max(np.abs(y1 / y2 - 1))
    rows.append({"A": A, "B": B, "Y": list(y2), "rel_two_node_sets": rel})
    print("convergence", f"A={A:.4g} B={B:.4g}", [f"{v:.12g}" for v in y2], f"rel {rel:.2e}", flush=True)
out["convergence"] = rows

import ring_ext_QG_8_onY as E
mp.mp.dps = 70
g60 = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
g32 = json.loads((DATA / "pf_grid_d32_h0.0625.json").read_text())
rows = []
for t in ["0.5", "1", "2", "4"]:
    a, b = E.thetaY(g60, t, 0)[0], E.thetaY(g32, t, 0)[0]
    rows.append({"t": t, "Y60": mp.nstr(a, 30), "rel_60_vs_32": mp.nstr(abs(b / a - 1), 3)})
    print("reference", t, mp.nstr(a, 25), "60 vs 32:", mp.nstr(abs(b / a - 1), 3), flush=True)
out["reference"] = rows
(DATA / "gf_checks.json").write_text(json.dumps(out, indent=1, default=float))
print("saved")
