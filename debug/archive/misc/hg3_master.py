"""Record HG (2026-10-27_h2_geminal_elliptic.md), HG3 addendum: the Lesiuk-Moszynski master integral
f(r) = r Int e^{-u3 r1A - u2 r1B - w2 r2A - w3 r2B - w1 r12} / (16 pi^2 r1A r1B r2A r2B r12) by an independent route,
momentum space (no use of LM's r-equation or w1-series).  Each Yukawa-pair density, by Feynman parameters, is
2 pi Int_0^1 dt e^{-i k.c(t)} e^{-R sqrt(Delta)}/sqrt(Delta), Delta = t(1-t)k^2 + t a^2 + (1-t) b^2 (a at A, b at B),
and the geminal kernel is 4 pi/(k^2 + w1^2), so
    f = (R / 2 pi) Int_0^inf k^2 dk/(k^2 + w1^2) Int Int dt ds j0(k R (t - s)) E1(t, k) E2(s, k),  E = e^{-R sqrt D}/sqrt D.
Double-exponential (tanh-sinh) rules in t, s and in u, k = c u/(1 - u); refinement levels show the convergence.
    python debug/hg3_master.py
"""
import math

import numpy as np


def ts(n, h):
    """Nodes t, 1 - t (each computed directly, no cancellation) and weights of the tanh-sinh rule on (0, 1)."""
    j = np.arange(-n, n + 1) * h
    y = np.pi * np.sinh(j)
    t, tc = 1 / (1 + np.exp(-y)), 1 / (1 + np.exp(y))
    wt = np.pi * np.cosh(j) * t * tc * h
    keep = (t > 0) & (tc > 0) & (wt > 1e-300)
    return t[keep], tc[keep], wt[keep]


def j0(x):
    out = np.ones_like(x)
    m = np.abs(x) > 1e-4
    out[m] = np.sin(x[m]) / x[m]
    xs = x[~m]
    out[~m] = 1 - xs ** 2 / 6 + xs ** 4 / 120
    return out


def master(u2, u3, w2, w3, w1, R, level, c=4.0):
    h = 2.0 ** -level
    t, tc, wt = ts(int(4.5 / h), h)
    u, uc, wu = ts(int(4.5 / h), h)
    k = c * u / uc
    wk = wu * c / uc ** 2
    T, S = t[:, None], t[None, :]
    total = []
    for kk, ww in zip(k, wk):
        D1 = t * tc * kk ** 2 + t * u3 ** 2 + tc * u2 ** 2      # density 1: u3 at A, u2 at B
        D2 = t * tc * kk ** 2 + t * w2 ** 2 + tc * w3 ** 2      # density 2: w2 at A, w3 at B
        E1 = wt * np.exp(-R * np.sqrt(D1)) / np.sqrt(D1)
        E2 = wt * np.exp(-R * np.sqrt(D2)) / np.sqrt(D2)
        inner = E1 @ j0(kk * R * (T - S)) @ E2
        total.append(ww * kk ** 2 / (kk ** 2 + w1 ** 2) * inner)
    return R / (2 * math.pi) * math.fsum(total)


if __name__ == '__main__':
    LM = {1: 2.39121170158262983284e-4}
    for w1 in (1.0,):
        prev = None
        for level in (3, 4, 5, 6):
            v = master(1.0, 2.0, 1.0, 2.0, w1, 2.0, level)
            print('w1', w1, 'level', level, repr(v), 'rel to LM', (v - LM[1]) / LM[1],
                  'step', None if prev is None else (v - prev) / v, flush=True)
            prev = v
