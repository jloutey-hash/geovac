"""Instrument toys for the proposed records (synthetic integrands with exactly known constants; not new shapes).

TOY1 (B1's instrument): S(eps) = int_1^oo I2((2-eps)x) x e^{-2x} dx  -> eps^(3/2) S -> A (1 + kappa eps)
TOY2 (B2's instrument, same cancellation as B2): S(eps) = int_1^oo [(I0(y) + I1(y)/y) - I2(y)] x e^{-2x} dx,
     y = (2-eps)x  -> eps^(1/2) S -> A (1 + kappa eps) + odd powers of u = sqrt(eps).
Constants derived by p63odd_tail.branch_expansion (H = x e^{-2x} exactly).  Evaluated float64 (log-grid
trapezoid, scaled Bessels, the P63L1 scheme) and mp30 (composite Gauss-Legendre, p63l1_engine.branch_mp);
fits: TOY1 {1, u^2, u^3, u^4}; TOY2 {1, u, u^2, u^3, u^4}; eps = 1e-2 .. 1e-4, 9 points.
Output: debug/data/p63odd_toy.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.special import i0e, i1e, ive

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63l1_engine as L1  # noqa: E402
import p63odd_tail as TL  # noqa: E402


def consts():
    H = TL.Ser(sp.Integer(1), 2, 1, [1] + [0] * (TL.N - 1))
    P1 = TL.powers(TL.branch_expansion([(1, TL.kernel_I(2), H)]))
    P2 = TL.powers(TL.branch_expansion([(1, TL.kernel_I(0, True), H), (-1, TL.kernel_I(2), H)], order=3))
    A1, k1 = P1[sp.Rational(-3, 2)], P1[sp.Rational(-1, 2)] / P1[sp.Rational(-3, 2)]
    A2, k2 = P2[sp.Rational(-1, 2)], P2[sp.Rational(1, 2)] / P2[sp.Rational(-1, 2)]
    return {"A1": sp.nsimplify(A1), "k1": sp.nsimplify(k1), "A2": sp.nsimplify(A2), "k2": sp.nsimplify(k2),
            "lead2_eps^-3/2": P2.get(sp.Rational(-3, 2), 0)}


def s_float(eps: float, which: int, npts: int = 400001) -> float:
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), npts)
    x = np.exp(y)
    yy = (2 - eps) * x
    k = ive(2, yy) if which == 1 else (i0e(yy) + i1e(yy) / yy - ive(2, yy))
    return float(np.trapezoid(k * x * np.exp(-eps * x) * x, y))


def s_mp(eps, which: int):
    # kernel inside branch_mp is I0 e^{-y}; supply h = x * [K(y)/I0(y)] by folding the ratio into h
    a2 = 2 - mp.mpf(eps)

    def h(x):
        yv = a2 * x
        if which == 1:
            r = mp.besseli(2, yv) / mp.besseli(0, yv)
        else:
            r = (mp.besseli(0, yv) + mp.besseli(1, yv) / yv - mp.besseli(2, yv)) / mp.besseli(0, yv)
        return x * r
    return L1.branch_mp(mp.mpf(eps), h, ngl=40)


def fit(eps_list, vals, powers):
    u = np.sqrt(np.array(eps_list))
    M = np.vstack([u**p for p in powers]).T
    c, *_ = np.linalg.lstsq(M, np.array(vals), rcond=None)
    return dict(zip(powers, c))


def main() -> None:
    C = consts()
    print(C)
    mp.mp.dps = 30
    eps_list = list(np.logspace(-2, -4, 9))
    out = {"consts": {k: str(v) for k, v in C.items()}, "rows": []}
    t0 = time.time()
    for eps in eps_list:
        row = {"eps": eps}
        for w in (1, 2):
            f = s_float(eps, w)
            m = s_mp(eps, w)
            row[f"toy{w}_float"] = f
            row[f"toy{w}_mp"] = mp.nstr(m, 20)
            row[f"toy{w}_rel"] = float(abs(f / m - 1))
        out["rows"].append(row)
        print(row, flush=True)
    out["seconds"] = time.time() - t0
    A1, k1, A2, k2 = (float(C[k]) for k in ("A1", "k1", "A2", "k2"))
    for w, pw, powers, A, kap in ((1, 1.5, [0, 2, 3, 4], A1, k1), (2, 0.5, [0, 1, 2, 3, 4], A2, k2)):
        vals = [float(mp.mpf(r[f"toy{w}_mp"])) * r["eps"] ** pw for r in out["rows"]]
        for tag, sl in (("all9", slice(0, 9)), ("smallest6", slice(3, 9))):
            c = fit(eps_list[sl], vals[sl], powers)
            out[f"toy{w}_fit_{tag}"] = {"A_ratio_minus_1": c[0] / A - 1, "u2_over_Akappa": c[2] / (A * kap)}
        logS = np.log([float(mp.mpf(r[f"toy{w}_mp"])) for r in out["rows"]])
        sl = np.polyfit(np.log(eps_list[4:]), logS[4:], 1)[0]
        out[f"toy{w}_slope_last5"] = sl
        print(w, out[f"toy{w}_fit_all9"], out[f"toy{w}_fit_smallest6"], sl)
    json.dump(out, open(HERE / "data" / "p63odd_toy.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
