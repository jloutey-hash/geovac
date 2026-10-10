"""Fit quality of the shapes Paper 58's NaH ladder actually uses (debug only).

<fit|STO> overlap and the hydrogenic one-electron energy error of each shape at zeta = 1
(scales as zeta^2): STO-6G 1s and STO-3G 2s (hardcoded in noci_engine, used for NaH's
1s and 2s) vs the engine's own LSQ fits at 6 / 12 Gaussians.  Radial quadrature on a fine
log grid (shape-quality diagnostic only).  Writes debug/data/sto_sturmian_eval_fitquality.json
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

from geovac import noci_engine as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                   "sto_sturmian_eval_fitquality.json")
r = np.geomspace(1e-7, 80.0, 400001)


def radial_from(a, d, l):
    f = np.zeros_like(r)
    for ai, di in zip(a, d):
        prim = r ** l * np.exp(-ai * r * r)
        f += di * prim / math.sqrt(np.trapezoid(prim * prim * r * r, r))
    return f / math.sqrt(np.trapezoid(f * f * r * r, r))


def quality(a, d, l, nr):
    f = radial_from(a, d, l)
    t = r ** (nr - 1) * np.exp(-r)
    t /= math.sqrt(np.trapezoid(t * t * r * r, r))
    ov = np.trapezoid(f * t * r * r, r)
    # one-electron energy in the hydrogen-like field Z = nr (exact STO energy -1/2 * ... below)
    df = np.gradient(f, r)
    T = 0.5 * np.trapezoid((df * df + l * (l + 1) * f * f / r ** 2) * r * r, r)
    dt = np.gradient(t, r)
    Tt = 0.5 * np.trapezoid((dt * dt + l * (l + 1) * t * t / r ** 2) * r * r, r)
    V = -np.trapezoid(f * f * r, r)
    Vt = -np.trapezoid(t * t * r, r)
    return {"overlap": float(ov), "1-ov": float(1 - ov), "dT": float(T - Tt), "dV": float(V - Vt),
            "dE(T+V, Z=1)": float(T + V - Tt - Vt)}


def main() -> None:
    out = {"STO-6G 1s": quality(*E.STO6G_1S, 0, 1), "STO-3G 2s": quality(*E.STO3G_2S, 0, 2)}
    for kind, (l, nr) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2)), ("3s", (0, 3))):
        for ng in (6, 12):
            a, d, q = E.fit_sto_shape(l, nr, n_gauss=ng)
            out[f"LSQ{ng} {kind}"] = quality(a, d, l, nr)
            out[f"LSQ{ng} {kind}"]["engine_quality"] = q
    for k, v in out.items():
        print(k, {kk: f"{vv:.3e}" for kk, vv in v.items()})
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
