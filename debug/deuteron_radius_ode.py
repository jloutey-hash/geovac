"""Independent ODE route for DR1-DR5: Numerov + Richardson, series start, shooting on B.

Registered: u = r + c_2 r^2 + ... at the origin; Numerov; Richardson in h; r_max = 60 fm;
shooting on B.  Outputs B, <r^2>, r_m, A_S (plateau of u e^{kappa r} at r >= 20 fm), a_t.
float64 throughout (an accuracy target of ~1e-10 relative is checked by the Richardson table).
"""
from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq

HB2M = 41.47
MT3 = [(1438.72, 3.11), (-626.885, 1.55)]
MT1S = [(1438.72, 3.11), (-513.968, 1.55)]
RMAX = 60.0


def series_coeffs(pot, E: float, order: int = 10) -> List[float]:
    d = [sum(g * (-mu) ** m / math.factorial(m) for g, mu in pot) for m in range(order + 2)]
    c = [0.0, 1.0]
    for n in range(0, order - 1):
        s = sum(d[n + 1 - l] * c[l] for l in range(1, n + 2) if l < len(c))
        s -= E * (c[n] if n < len(c) else 0.0)
        c.append(s / (HB2M * (n + 1) * (n + 2)))
    return c


def integrate(pot, E: float, h: float, rmax: float = RMAX) -> Tuple[np.ndarray, np.ndarray]:
    M = int(round(rmax / h))
    r = np.arange(M + 1) * h
    V = np.zeros(M + 1)
    for g, mu in pot:
        V[1:] += g * np.exp(-mu * r[1:]) / r[1:]
    f = (V - E) / HB2M
    A = sum(g for g, mu in pot)
    c = series_coeffs(pot, E)
    u = np.zeros(M + 1)
    u[1] = sum(ci * h ** i for i, ci in enumerate(c))
    h2 = h * h / 12.0
    w_prev = -h2 * A / HB2M  # w_0 = u_0 - h^2 (f u)_0 / 12 with (f u)_0 = A/hb2m
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


def bound_state(pot, h: float, Elo: float = -3.0, Ehi: float = -1.5) -> Dict:
    def end(E):
        return integrate(pot, E, h)[1][-1]
    E = brentq(end, Elo, Ehi, xtol=1e-15, rtol=1e-15, maxiter=200)
    r, u = integrate(pot, E, h)
    if u[int(len(u) * 0.3)] < 0:
        u = -u
    m50 = int(round(50.0 / h))
    nodes = int(np.sum(np.sign(u[1:m50]) != np.sign(u[2:m50 + 1])))  # nodes on (0, 50] fm
    norm = simpson(u * u, h)
    u = u / math.sqrt(norm)
    r2 = simpson(r * r * u * u, h)
    B = -E
    kap = math.sqrt(B / HB2M)
    AS = u * np.exp(kap * r)
    wall = 1 - np.exp(-2 * kap * (RMAX - r))
    sel = (r >= 20) & (r <= 30)
    sel816 = (r >= 8) & (r <= 16)
    return {"B": B, "r2": r2, "nodes": nodes, "AS_raw_20_30": float(AS[sel].mean()),
            "AS_corr_20_30": float((AS[sel] / wall[sel]).mean()),
            "AS_corr_spread_20_30": float(np.ptp(AS[sel] / wall[sel]) / (AS[sel] / wall[sel]).mean()),
            "AS_816_mean": float((AS[sel816] / wall[sel816]).mean()),
            "AS_816_spread": float(np.ptp(AS[sel816] / wall[sel816]) / (AS[sel816] / wall[sel816]).mean())}


def scattering_length(pot, h: float) -> Tuple[float, int]:
    r, u = integrate(pot, 0.0, h)
    i1 = int(round(50.0 / h))
    i2 = len(r) - 1
    a = (r[i1] * u[i2] - r[i2] * u[i1]) / (u[i2] - u[i1])
    nodes = int(np.sum(np.sign(u[1:-1]) != np.sign(u[2:])))
    return a, nodes


def richardson(vals: Sequence[float]) -> Tuple[float, float]:
    """vals at h, h/2, h/4 with error series h^4, h^6. Returns (extrapolated, |last step change|)."""
    a1 = [(16 * vals[i + 1] - vals[i]) / 15 for i in range(len(vals) - 1)]
    if len(a1) >= 2:
        a2 = [(64 * a1[i + 1] - a1[i]) / 63 for i in range(len(a1) - 1)]
        return a2[-1], abs(a2[-1] - a1[-1])
    return a1[-1], abs(a1[-1] - vals[-1])


if __name__ == "__main__":
    import json
    hs = [4e-3, 2e-3, 1e-3, 5e-4]
    out = {"h": hs, "rmax": RMAX}
    runs = [bound_state(MT3, h) for h in hs]
    for key in ["B", "r2", "AS_raw_20_30", "AS_corr_20_30", "AS_816_mean"]:
        vals = [x[key] for x in runs]
        ext, err = richardson(vals[-3:])
        out[key] = {"raw": vals, "richardson": ext, "rich_last_change": err}
    out["nodes"] = [x["nodes"] for x in runs]
    out["AS_corr_spread_20_30"] = runs[-1]["AS_corr_spread_20_30"]
    out["AS_816_spread_ode"] = runs[-1]["AS_816_spread"]
    r2 = out["r2"]["richardson"]
    out["r_m"] = 0.5 * math.sqrt(r2)
    out["kappa"] = math.sqrt(out["B"]["richardson"] / HB2M)
    at = [scattering_length(MT3, h) for h in hs]
    vals = [x[0] for x in at]
    ext, err = richardson(vals[-3:])
    out["a_t"] = {"raw": vals, "richardson": ext, "rich_last_change": err,
                  "nodes_zero_energy": [x[1] for x in at]}
    asg = [scattering_length(MT1S, h) for h in hs]
    vals = [x[0] for x in asg]
    ext, err = richardson(vals[-3:])
    out["singlet_MT1"] = {"a_s": ext, "rich_last_change": err,
                          "nodes_zero_energy": [x[1] for x in asg]}
    print(json.dumps(out, indent=1))
    with open("data/deuteron_radius_ode.json", "w") as fh:
        json.dump(out, fh, indent=1)
