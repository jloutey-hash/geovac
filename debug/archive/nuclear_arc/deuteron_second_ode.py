"""Independent ODE route for DS1-DS5: copied from deuteron_radius_ode.py, generalized to Yukawa +
Gaussian potentials and a per-potential hbar^2/m.  Numerov + Richardson (h^4, h^6), series start
u = r + c_2 r^2 + ... (coefficients from the Taylor series of r V(r)), r_max = 60 fm, shooting on B.
Outputs B, <r^2>, r_m, A_S (wall-corrected plateau of u e^{kappa r} on [20, 30] fm, as in the first
run's D4), the [8, 16] fm ODE spread, and a_t (zero-energy solution at 50 and 60 fm).
"""
from __future__ import annotations

import json
import math
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq

RMAX = 60.0
HS = [4e-3, 2e-3, 1e-3, 5e-4]


def fpot(name: str):
    """(float potential list, hb2m) from the exact registry in deuteron_second_laguerre."""
    import deuteron_second_laguerre as L
    pot, hb = L.POTS[name]
    return [(k, float(g), float(p)) for k, g, p in pot], float(hb)


def rV_taylor(pot, order: int) -> List[float]:
    """d_m: coefficient of r^m in r V(r)."""
    d = [0.0] * (order + 2)
    for kind, g, p in pot:
        for m in range(order + 2):
            if kind == "y":
                d[m] += g * (-p) ** m / math.factorial(m)
            elif m % 2 == 1:
                j = (m - 1) // 2
                d[m] += g * (-p) ** j / math.factorial(j)
    return d


def series_coeffs(pot, hb: float, E: float, order: int = 10) -> List[float]:
    d = rV_taylor(pot, order)
    c = [0.0, 1.0]
    for n in range(0, order - 1):
        s = sum(d[n + 1 - l] * c[l] for l in range(1, n + 2) if l < len(c))
        s -= E * (c[n] if n < len(c) else 0.0)
        c.append(s / (hb * (n + 1) * (n + 2)))
    return c


def Varr(pot, r: np.ndarray) -> np.ndarray:
    V = np.zeros_like(r)
    for kind, g, p in pot:
        if kind == "y":
            V[1:] += g * np.exp(-p * r[1:]) / r[1:]
        else:
            V += g * np.exp(-p * r * r)
    return V


def integrate(pot, hb: float, E: float, h: float, rmax: float = RMAX):
    M = int(round(rmax / h))
    r = np.arange(M + 1) * h
    V = Varr(pot, r)
    f = (V - E) / hb
    d0 = rV_taylor(pot, 1)[0]  # (f u)_0 = d_0 / hb
    c = series_coeffs(pot, hb, E)
    u = np.zeros(M + 1)
    u[1] = sum(ci * h ** i for i, ci in enumerate(c))
    h2 = h * h / 12.0
    w_prev = -h2 * d0 / hb
    w_cur = (1 - h2 * f[1]) * u[1]
    ff = f.tolist()
    uu = u.tolist()
    for i in range(1, M):
        w_next = 2 * w_cur - w_prev + 12 * h2 * ff[i] * uu[i]
        uu[i + 1] = w_next / (1 - h2 * ff[i + 1])
        w_prev, w_cur = w_cur, w_next
    return r, np.array(uu)


def simpson(y: np.ndarray, h: float) -> float:
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def bound_state(pot, hb: float, h: float, Elo: float = -3.0, Ehi: float = -1.5) -> Dict:
    def end(E):
        return integrate(pot, hb, E, h)[1][-1]
    E = brentq(end, Elo, Ehi, xtol=1e-15, rtol=1e-15, maxiter=200)
    r, u = integrate(pot, hb, E, h)
    if u[int(len(u) * 0.3)] < 0:
        u = -u
    m50 = int(round(50.0 / h))
    nodes = int(np.sum(np.sign(u[1:m50]) != np.sign(u[2:m50 + 1])))
    norm = simpson(u * u, h)
    u = u / math.sqrt(norm)
    r2 = simpson(r * r * u * u, h)
    B = -E
    kap = math.sqrt(B / hb)
    AS = u * np.exp(kap * r)
    wall = 1 - np.exp(-2 * kap * (RMAX - r))
    sel = (r >= 20) & (r <= 30)
    sel816 = (r >= 8) & (r <= 16)
    return {"B": B, "r2": r2, "nodes": nodes, "AS_raw_20_30": float(AS[sel].mean()),
            "AS_corr_20_30": float((AS[sel] / wall[sel]).mean()),
            "AS_corr_spread_20_30": float(np.ptp(AS[sel] / wall[sel]) / (AS[sel] / wall[sel]).mean()),
            "AS_816_mean": float((AS[sel816] / wall[sel816]).mean()),
            "AS_816_spread": float(np.ptp(AS[sel816] / wall[sel816]) / (AS[sel816] / wall[sel816]).mean())}


def scattering_length(pot, hb: float, h: float) -> Tuple[float, int]:
    r, u = integrate(pot, hb, 0.0, h)
    i1 = int(round(50.0 / h))
    i2 = len(r) - 1
    a = (r[i1] * u[i2] - r[i2] * u[i1]) / (u[i2] - u[i1])
    nodes = int(np.sum(np.sign(u[1:-1]) != np.sign(u[2:])))
    return a, nodes


def richardson(vals: Sequence[float]) -> Tuple[float, float]:
    a1 = [(16 * vals[i + 1] - vals[i]) / 15 for i in range(len(vals) - 1)]
    if len(a1) >= 2:
        a2 = [(64 * a1[i + 1] - a1[i]) / 63 for i in range(len(a1) - 1)]
        return a2[-1], abs(a2[-1] - a1[-1])
    return a1[-1], abs(a1[-1] - vals[-1])


def run(name: str) -> Dict:
    pot, hb = fpot(name)
    out = {"potential": name, "hb2m": hb, "h": HS, "rmax": RMAX}
    runs = [bound_state(pot, hb, h) for h in HS]
    for key in ["B", "r2", "AS_raw_20_30", "AS_corr_20_30", "AS_816_mean"]:
        vals = [x[key] for x in runs]
        ext, err = richardson(vals[-3:])
        out[key] = {"raw": vals, "richardson": ext, "rich_last_change": err}
    out["nodes"] = [x["nodes"] for x in runs]
    out["AS_corr_spread_20_30"] = runs[-1]["AS_corr_spread_20_30"]
    out["AS_816_spread_ode"] = runs[-1]["AS_816_spread"]
    out["r_m"] = 0.5 * math.sqrt(out["r2"]["richardson"])
    out["kappa"] = math.sqrt(out["B"]["richardson"] / hb)
    at = [scattering_length(pot, hb, h) for h in HS]
    vals = [x[0] for x in at]
    ext, err = richardson(vals[-3:])
    out["a_t"] = {"raw": vals, "richardson": ext, "rich_last_change": err,
                  "nodes_zero_energy": [x[1] for x in at]}
    return out


if __name__ == "__main__":
    res = {name: run(name) for name in ("MT3", "P1", "P2")}
    for name, o in res.items():
        print(name, "B", o["B"]["richardson"], o["B"]["rich_last_change"], "r_m", o["r_m"],
              "dr2", o["r2"]["rich_last_change"], "a_t", o["a_t"]["richardson"], o["a_t"]["rich_last_change"],
              "A_S", o["AS_corr_20_30"]["richardson"], "nodes", o["nodes"], o["a_t"]["nodes_zero_energy"])
    with open("data/deuteron_second_ode.json", "w") as fh:
        json.dump(res, fh, indent=1)
