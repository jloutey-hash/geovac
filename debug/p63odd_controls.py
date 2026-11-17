"""Controls for the (B) design (Paper 63, l > 0 attainment shapes B1, B2).  No branch integral of B1 or B2 is
evaluated at eps < 0.1.

C1  generic tail route (p63odd_engine.phix_f, L = 1, k = 1, source f_p) reproduces the closed form e^x E(x) of
    record P63L1 (p63l1_engine), x = 1, 3, 10, 40, 80.
C2  positive control on the scored P63L1 object: its branch S(eps) rebuilt from the generic route
    (h = (2/(3pi)) x^2 E_generic^2) vs the committed scored float values (debug/data/p63l1_run.json), 9 eps.
C3  reduction: B1 and B2 fibres at t = 0.5, 1.3 rooted at the root vs at the far leaf (different Gaunt paths).
C4  integrand instrument: e^{2x} Phi_Lx (L = 0, 2) by the Laguerre route vs a cumulative-tail grid route.
C5  local-rule tail: e^{2x} Phi_Lx vs the derived series (p63odd_tail.local_potential), x = 20..60, L = 0, 2.
C6  B1, B2 branch integrals at eps = 0.5, 0.2, 0.1 only: 10k vs 20k trapezoid nodes agree; wall time (the
    400k-node P63L1 grid costs ~80x more per node here, since each node needs the Laguerre R inside a Laguerre
    potential; the full C6 at 400k/800k was stopped after 45 min).
Output: debug/data/p63odd_controls.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import k0, k1, k1e, kv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63l1_engine as L1  # noqa: E402
import p63odd_engine as EN  # noqa: E402
import p63odd_tail as TL  # noqa: E402


def c1() -> dict:
    xs = np.array([1.0, 3.0, 10.0, 40.0, 80.0])
    w = lambda u: -k1e(u) / (2 * np.pi**2)  # noqa: E731  e^{u} f_p
    gen = EN.phix_f(xs, 1, w=w, k=1.0)
    closed = 2 / np.pi * (k1(xs) * np.exp(xs) + L1.Rscaled_f(xs) / xs**2)
    return {"x": xs.tolist(), "max_rel": float(np.max(np.abs(gen / closed - 1)))}


def c2() -> dict:
    ref = json.load(open(HERE / "data" / "p63l1_run.json"))
    rows = ref["rows"] if "rows" in ref else ref["points"]
    w = lambda u: -k1e(u) / (2 * np.pi**2)  # noqa: E731

    def h(x):
        Eg = EN.phix_f(x, 1, w=w, k=1.0)
        return 2 / (3 * np.pi) * x**2 * Eg**2
    out = []
    for r in rows:
        eps = r["eps"]
        y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), 400001)
        x = np.exp(y)
        f = np.empty_like(x)
        for i in range(0, x.size, 20000):
            xs = x[i:i + 20000]
            f[i:i + 20000] = L1.i0e((2 - eps) * xs) * h(xs) * np.exp(-eps * xs)
        s = float(np.trapezoid(f * x, y))
        sref = float(r["S_mp"])
        out.append({"eps": eps, "S_generic": s, "S_committed_mp30": sref, "rel": abs(s / sref - 1)})
    return {"rows": out, "max_rel": max(o["rel"] for o in out)}


# ---------------------------------------------------------------- C3 reduction (grid, as the tree_rung_l guard)
_u = np.linspace(np.log(1e-9), np.log(80.0), 200001)
_r = np.exp(_u)
_f = -k1(_r) / (2 * np.pi**2)


def _cum0(y):
    return cumulative_simpson(y * _r, x=_u, initial=0.0)


def _tail(y):
    zz = (y * _r)[::-1]
    return cumulative_simpson(zz, x=-_u[::-1], initial=0.0)[::-1]


def _phi(L, w):
    return 4 * np.pi / (2 * L + 1) * (_r ** (-(L + 1)) * _cum0(_r ** (L + 2) * w) + _r**L * _tail(_r ** (1 - L) * w))


def c3() -> dict:
    out = {}
    g1 = _phi(1, _f)
    for t in (0.5, 1.3):
        rr = t * _r
        root0 = t**5 / (2 * np.pi**2) * (k0(rr) - k1(rr) / rr)      # channel 0 radial of lap t^3 K0/(2pi^2)
        root2 = t**5 / (2 * np.pi**2) * kv(2, rr)                     # h_2 part: radial of P2
        # B1 rooted at the d root and at the far leaf (guard route)
        I1a = 4 * np.pi / 5 * simpson(_r**2 * root2 * _phi(2, 2 / 3 * _f * g1) * _r, x=_u)
        I1b = 4 * np.pi / 3 * simpson(_r**2 * _f * _phi(1, 2 / 5 * _f * _phi(2, root2)) * _r, x=_u)
        # B2: root = root0 - (5/2) root2 P2
        I2a = (4 * np.pi * simpson(_r**2 * root0 * _phi(0, _f * g1 / 3) * _r, x=_u)
               - 5 / 2 * 4 * np.pi / 5 * simpson(_r**2 * root2 * _phi(2, 2 / 3 * _f * g1) * _r, x=_u))
        src_P1 = _f * (_phi(0, root0) + 2 / 5 * (-5 / 2) * _phi(2, root2))   # cos * (Phi0 + Phi2 P2) -> P1 part
        I2b = 4 * np.pi / 3 * simpson(_r**2 * _f * _phi(1, src_P1) * _r, x=_u)
        out[t] = {"B1_root": I1a, "B1_leaf": I1b, "B1_rel": abs(I1a / I1b - 1),
                  "B2_root": I2a, "B2_leaf": I2b, "B2_rel": abs(I2a / I2b - 1)}
    return out


# ---------------------------------------------------------------- C4 / C5 integrand
def _grid_phix(L: int, xs: np.ndarray) -> np.ndarray:
    """cumulative-tail route on a fine linear grid in [1, 200], scaled by e^{2x}."""
    u = np.linspace(1.0, 200.0, 2_000_001)
    w = EN.what_f(u) * np.exp(-2 * (u - 1.0))          # e^{2} * source, avoids underflow
    def tail(y):
        return cumulative_simpson(y[::-1], x=-u[::-1], initial=0.0)[::-1]
    A = tail(u ** (1 - L) * w)
    B = tail(u ** (L + 2) * w)
    phi = 4 * np.pi / (2 * L + 1) * (u**L * A - u ** (-L - 1) * B) * np.exp(2 * (u - 1.0))
    return np.interp(xs, u, phi)


def c4() -> dict:
    xs = np.array([1.0, 3.0, 10.0, 40.0])
    out = {}
    for L in (0, 2):
        a = EN.phix_f(xs, L)
        b = _grid_phix(L, xs)
        out[L] = {"x": xs.tolist(), "laguerre": a.tolist(), "grid": b.tolist(), "max_rel": float(np.max(np.abs(a / b - 1)))}
    return out


def c5() -> dict:
    f_p = TL.K_ser(1).scale(-1 / (2 * sp.pi**2))
    E_full = TL.local_potential(f_p, 1)
    w0 = f_p * E_full
    out = {}
    xs = np.array([20.0, 30.0, 40.0, 50.0, 60.0])
    for L in (0, 2):
        Ph = TL.local_potential(w0, L)
        pref = float(Ph.pref)
        cs = [float(c) for c in Ph.c]
        ser = pref * xs ** float(Ph.a) * sum(c * xs ** (-j) for j, c in enumerate(cs))
        num = EN.phix_f(xs, L)
        out[L] = {"x": xs.tolist(), "rel_resid": np.abs(num / ser - 1).tolist(), "terms": len(cs)}
    return out


def c6() -> dict:
    out = {}
    for eps in (0.5, 0.2, 0.1):
        row = {}
        for name, fn in (("B1", EN.branch1_f), ("B2", EN.branch2_f)):
            t0 = time.time()
            a = fn(eps, 10001)
            ta = time.time() - t0
            b = fn(eps, 20001)
            row[name] = {"S_10k": a, "S_20k": b, "rel": abs(a / b - 1), "seconds_10k": ta}
        out[eps] = row
    return out


def main() -> None:
    res = {}
    only = sys.argv[1:] or ["C1", "C2", "C3", "C4", "C5", "C6"]
    allf = {"C1": c1, "C2": c2, "C3": c3, "C4": c4, "C5": c5, "C6": c6}
    for name, fn in ((k, allf[k]) for k in only):
        t0 = time.time()
        try:
            res[name] = fn()
        except Exception as e:  # report, don't hide
            res[name] = {"error": repr(e)}
        res[name + "_seconds"] = time.time() - t0
        print(name, json.dumps(res[name], default=str)[:600], flush=True)
    path = HERE / "data" / ("p63odd_controls.json" if not sys.argv[1:] else "p63odd_controls_" + "_".join(only) + ".json")
    json.dump(res, open(path, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
