"""WA5 driver (record docs/predictions/2026-10-27_water_cusp_singularity.md).

Modes: control (WA5-0), cases (WA5a-c).  F_coeffs is copied unchanged from tests/test_paper59_water_bond_angle.py.
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss, legval
from scipy.special import spherical_jn

REFEREE_NG16 = 0.204941652260
WATER = 1.809
CASES = [(WATER, WATER), (WATER, 2.5), (WATER, 1.0)]
SETS = {'E1': (1.0, 1.0), 'E2': (2.25, 1.24)}
LMAX = 60
OUT = Path(__file__).parent / "data" / "wa5_cusp.json"


def _tanh_sinh(n, h):
    j = np.arange(-n, n + 1) * h
    x = np.tanh(0.5 * np.pi * np.sinh(j))
    w = 0.5 * np.pi * np.cosh(j) / np.cosh(0.5 * np.pi * np.sinh(j)) ** 2 * h
    t, wt = 0.5 * (x + 1), 0.5 * w
    keep = (t > 0) & (t < 1)
    return t[keep], wt[keep]


def _G(k, t, D, zA, zB):
    N = np.sqrt(zA ** 3 / np.pi) * np.sqrt(zB ** 3 / np.pi)
    Dl = np.sqrt(k ** 2 * t * (1 - t) + t * zA ** 2 + (1 - t) * zB ** 2)
    return N * 2 * np.pi * t * (1 - t) * zA * zB * np.exp(-D * Dl) * (D ** 2 * Dl ** 2 + 3 * D * Dl + 3) / Dl ** 5


def F_coeffs(D1, D2, zX, zH, lmax, level=0, k_scale=3.0):
    t, wt = _tanh_sinh(60 * 2 ** level, 0.07 / 2 ** level)
    u, wu = leggauss(400 * 2 ** level)
    u, wu = 0.5 * (u + 1), 0.5 * wu
    k = k_scale * u / (1 - u)
    wk = wu * k_scale / (1 - u) ** 2
    G1 = _G(k[:, None], t[None, :], D1, zX, zH)
    G2 = _G(k[:, None], t[None, :], D2, zX, zH)
    out = []
    for l in range(lmax + 1):
        A1 = (spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt
        A2 = (spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt
        out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))
    return np.array(out)


def leg_coeffs(f, lmax, xs):
    """Legendre coefficients of f on [-1, 1], 40 digits; xs = breakpoints."""
    mp.mp.dps = 40
    return np.array([float((2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), xs))
                     for l in range(lmax + 1)])


def singular_L(D1, D2, zX, zH, lmax, power=mp.mpf(7) / 2):
    C = mp.mpf(zH) ** 5 * mp.mpf(zX) ** 3 * mp.e ** (-mp.mpf(zX) * (D1 + D2)) / 630
    xs_ = (mp.mpf(D1) ** 2 + mp.mpf(D2) ** 2) / (2 * mp.mpf(D1) * D2)
    pref = C * (2 * mp.mpf(D1) * D2) ** power
    return leg_coeffs(lambda x: pref * (xs_ - x) ** power, lmax, mp.linspace(-1, 1, 13))


def fit_rate_fixed_p(F, ls, p):
    """WA5a: log|F_l| + p log l = a + l log r + b/l + c/l^2."""
    l = ls.astype(float)
    A = np.vstack([np.ones_like(l), l, 1 / l, 1 / l ** 2]).T
    c = np.linalg.lstsq(A, np.log(np.abs(F[ls])) + p * np.log(l), rcond=None)[0]
    return float(np.exp(c[1]))


def fit_amp(ratio, ls):
    """WA5b: F_l / L_l = A0 + A1/l + A2/l^2."""
    l = ls.astype(float)
    A = np.vstack([np.ones_like(l), 1 / l, 1 / l ** 2]).T
    return float(np.linalg.lstsq(A, ratio[ls], rcond=None)[0][0])


def fit_power(F, ls):
    """WA5c: log|F_l| = a - p log l + b/l, r fixed at 1."""
    l = ls.astype(float)
    A = np.vstack([np.ones_like(l), -np.log(l), 1 / l]).T
    return float(np.linalg.lstsq(A, np.log(np.abs(F[ls])), rcond=None)[0][1])


def control():
    Y, Z = np.array([0, 0, 2.0]), np.array([1.5, 0, 0.5])
    D1, D2 = np.linalg.norm(Y), np.linalg.norm(Z)
    F = F_coeffs(D1, D2, 1.0, 1.0, 20)
    v = legval(Y @ Z / (D1 * D2), F)
    print(f"WA5-0(a) channel sum {v:.12f}  diff {v - REFEREE_NG16:.2e}  {'PASS' if abs(v - REFEREE_NG16) < 1e-9 else 'FAIL'}")
    ls = np.arange(15, LMAX + 1)
    for (D1, D2) in CASES:
        xs_ = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        mp.mp.dps = 40
        xsm = mp.mpf(xs_)
        f = lambda x: (xsm - x) ** (mp.mpf(7) / 2) + 3 * (xsm - x) ** (mp.mpf(9) / 2)
        S = leg_coeffs(f, LMAX, mp.linspace(-1, 1, 13))
        L = leg_coeffs(lambda x: (xsm - x) ** (mp.mpf(7) / 2), LMAX, mp.linspace(-1, 1, 13))
        A0 = fit_amp(S / L, ls)
        msg = f"WA5-0(b) ({D1}, {D2}): A0 = {A0:.5f} {'PASS' if abs(A0 - 1) < 1e-2 else 'FAIL'}"
        if D1 != D2:
            rho = min(D1, D2) / max(D1, D2)
            r = fit_rate_fixed_p(S, ls, 4.0)
            msg += f"; r/rho = {r / rho:.6f} {'PASS' if abs(r / rho - 1) < 1e-3 else 'FAIL'}"
        print(msg, flush=True)


def cases():
    res = {}
    for name, (zX, zH) in SETS.items():
        for (D1, D2) in CASES:
            F1 = F_coeffs(D1, D2, zX, zH, LMAX, level=1)
            F2 = F_coeffs(D1, D2, zX, zH, LMAX, level=2)
            floor = np.abs(F2 - F1)
            ls = np.array([l for l in range(15, LMAX + 1) if abs(F2[l]) >= 100 * floor[l]])
            # keep the contiguous run from 15
            run = [l for i, l in enumerate(ls) if l == 15 + i]
            ls = np.array(run)
            L = singular_L(D1, D2, zX, zH, LMAX)
            key = f"{name} ({D1}, {D2})"
            row = dict(usable=[int(ls[0]), int(ls[-1])] if len(ls) else [], n=len(ls),
                       F=F2.tolist(), floor=floor.tolist(), L=L.tolist())
            if len(ls) >= 10:
                row["A0"] = fit_amp(F2 / L, ls)
                row["ratio_at_ends"] = [float(F2[ls[0]] / L[ls[0]]), float(F2[ls[-1]] / L[ls[-1]])]
                if D1 != D2:
                    rho = min(D1, D2) / max(D1, D2)
                    row["r_over_rho"] = fit_rate_fixed_p(F2, ls, 4.0) / rho
                else:
                    row["p"] = fit_power(F2, ls)
            res[key] = row
            print(key, {k: v for k, v in row.items() if k not in ("F", "floor", "L")}, flush=True)
    OUT.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"control": control, "cases": cases}[sys.argv[1]]()
