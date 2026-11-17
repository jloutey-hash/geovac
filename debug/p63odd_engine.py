"""Paper 63 owed item (B): engine for two l > 0 attainment shapes (design pass; no scored evaluation here).

OBJECTS (conventions of debug/p63l1_engine.py and tests/test_paper63_tree_rung_l.py).
  p_z density (scale 1):  f_p(r) cos(theta),  f_p = -K1(r)/(2 pi^2).
  p_z leaf potential:     g(r) cos(theta),    g = -1/r^2 + E,  E = (2/pi)[K1 + R/r^2]  (record P63L1).
  middle p_z times leaf potential: W = f_p g cos^2 = f_p g [1/3 + (2/3) P2]; its e^{-2r} part is f_p E cos^2.
  B1 root (d):     h_2(d) t^3 K0(t r)/(2 pi^2) = t^5 K2(t r) P2/(2 pi^2).
     FIBRE_1(t) = (4pi/5) t^5/(2pi^2) int K2(t r) r^2 Phi_2[(2/3) f_p g] dr.
  B2 root (mixed): (lap - (5/2) h_2(d)) t^3 K0(t r)/(2 pi^2) = t^5/(2pi^2) [(K0 - K1/(tr)) - (5/2) K2 P2].
     FIBRE_2(t) = (4pi/3) t^5/(2pi^2) int r^2 [(K0 - K1/(tr)) Phi_0[f_p g] - K2 Phi_2[f_p g]] dr.
  Branch integrals at t = -(2 - eps) + i0 (only the e^{-2r} parts of the potentials reach -2):
     S1(eps) = int_1^oo I2(y) x^2 Phi2x[(2/3) f_p E](x) dx,
     S2(eps) = int_1^oo x^2 [(I0(y) + I1(y)/y) Phi0x[f_p E](x) - I2(y) Phi2x[f_p E](x)] dx,   y = (2 - eps) x,
  Phi_Lx[w](x) = (4pi/(2L+1)) [x^L int_x^oo u^{1-L} w - x^{-L-1} int_x^oo u^{L+2} w]
               = (4pi/(2L+1)) x int_0^oo (1+v/x)^{L+2} expm1(-(2L+1) log1p(v/x)) w(x+v) dv   (stable form).
The design memo (debug/p63odd_design_memo.md) derives the constants; this file only evaluates integrands and the
branch integral, and is to be run at eps < 0.1 only in the scored phase of a registered record.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Callable, Dict, Tuple

import mpmath as mp
import numpy as np
from scipy.special import i0e, i1e, ive, k1e, roots_laguerre

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63l1_engine as L1  # noqa: E402  (R, E closed forms, float and mp; frozen in record P63L1)

_UL, _WL = roots_laguerre(120)          # weight e^{-U}; v = U/2 for e^{-2v}


# ------------------------------------------------------------------ float64
def what_f(u: np.ndarray) -> np.ndarray:
    """e^{2u} f_p(u) E(u)."""
    Es = 2 / np.pi * (k1e(u) + L1.Rscaled_f(u) / u**2)
    return -k1e(u) * Es / (2 * np.pi**2)


def phix_f(x: np.ndarray, L: int, w: Callable = what_f, k: float = 2.0) -> np.ndarray:
    """e^{k x} Phi_Lx[w](x) for a source with e^{-k u} removed (w = e^{k u} source)."""
    V = _UL[None, :] / k
    rho = V / x[:, None]
    lp = np.log1p(rho)
    ker = np.exp((L + 2) * lp) * np.expm1(-(2 * L + 1) * lp)
    X = x[:, None] + V
    vals = w(X.ravel()).reshape(X.shape)
    return 4 * np.pi / (2 * L + 1) * x * ((ker * vals) @ _WL) / k


def h1_f(x: np.ndarray) -> np.ndarray:
    return x**2 * (2 / 3) * phix_f(x, 2)


def h2_pair_f(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    return x**2 * phix_f(x, 0), x**2 * phix_f(x, 2)


def branch1_f(eps: float, npts: int = 400001) -> float:
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), npts)
    x = np.exp(y)
    out = np.empty_like(x)
    for i in range(0, x.size, 20000):
        xs = x[i:i + 20000]
        out[i:i + 20000] = ive(2, (2 - eps) * xs) * h1_f(xs) * np.exp(-eps * xs)
    return float(np.trapezoid(out * x, y))


def branch2_f(eps: float, npts: int = 400001) -> float:
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), npts)
    x = np.exp(y)
    out = np.empty_like(x)
    for i in range(0, x.size, 20000):
        xs = x[i:i + 20000]
        yy = (2 - eps) * xs
        h0, h2 = h2_pair_f(xs)
        out[i:i + 20000] = ((i0e(yy) + i1e(yy) / yy) * h0 - ive(2, yy) * h2) * np.exp(-eps * xs)
    return float(np.trapezoid(out * x, y))


# ------------------------------------------------------------------ mpmath
_GLAG: Dict = {}


def _laguerre_mp(n: int):
    key = (n, mp.mp.dps)
    if key not in _GLAG:
        # Golub-Welsch would need eig at high precision; use mpmath's polyroots on L_n via recurrence
        xs = sorted(mp.polyroots([mp.mpf((-1) ** k) * mp.binomial(n, k) / mp.factorial(k) for k in range(n, -1, -1)],
                                 maxsteps=400, extraprec=4 * n))
        ws = []
        for xk in xs:
            xk = mp.re(xk)
            Lnp1 = mp.laguerre(n + 1, 0, xk)
            ws.append(xk / ((n + 1) ** 2 * Lnp1**2))
        _GLAG[key] = ([mp.re(x) for x in xs], ws)
    return _GLAG[key]


WDPS = 30


def what_mp(u) -> mp.mpf:
    """evaluated at the fixed working precision WDPS (mp.quad raises the context precision internally, beyond what
    the frozen R series of p63l1_engine supports at x >= 80)."""
    with mp.workdps(WDPS):
        u = mp.mpf(u)
        return -mp.besselk(1, u) * mp.exp(u) * L1.Escaled_mp(u) / (2 * mp.pi**2)


def phix_mp(x, L: int, n: int = 60) -> mp.mpf:
    """same as phix_f, adaptive mp.quad in v on [0, oo) (cached nodes would be faster; this is the referee)."""
    x = mp.mpf(x)

    def f(v):
        rho = v / x
        lp = mp.log1p(rho)
        return mp.exp((L + 2) * lp) * mp.expm1(-(2 * L + 1) * lp) * what_mp(x + v) * mp.exp(-2 * v)
    return 4 * mp.pi / (2 * L + 1) * x * mp.quad(f, [0, 1, 4, 12, 40, mp.inf], maxdegree=8)


def phix_local_series(x, L: int) -> mp.mpf:
    """not used by the engine: the local-rule series lives in p63odd_tail.py (derivation side)."""
    raise NotImplementedError


if __name__ == "__main__":
    # smoke + timing only (no branch evaluation at eps < 0.1)
    xs = np.array([1.0, 3.0, 10.0, 40.0])
    t0 = time.time()
    print("phi2 float", phix_f(xs, 2), time.time() - t0)
    mp.mp.dps = 30
    t0 = time.time()
    print("phi2 mp", [mp.nstr(phix_mp(x, 2), 15) for x in xs], time.time() - t0)
