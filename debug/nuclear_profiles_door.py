"""Nuclear profiles: tests NP1-NP3 of docs/predictions/2026-10-18_nuclear_profiles.md.

Registered in commit 8292bf51 before this script ran. Writes debug/data/nuclear_profiles_door.json.
Both profiles use Seaton's limit: the s-wave defect is the zero-energy phase shift over pi, from
matching the interior solution to the exterior pair sqrt(2r) J_1(sqrt(8r)), -sqrt(2r) Y_1(sqrt(8r)).
"""
from __future__ import annotations

import json
import sys
from typing import Dict

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
from finite_nucleus_threshold_series import Lg, exterior, x  # noqa: E402

OUT: Dict[str, object] = {}


def shell_series(order: int = 5):
    w = sp.sqrt(2 * x)
    lam = sp.series(w * sp.cos(w) / sp.sin(w), x, 0, order + 1).removeO()
    f, g, xf, xg = exterior()
    t = (xf - lam * f) / (lam * g - xg)
    tser = sp.series(t, x, 0, order + 1).removeO()
    mu = sp.expand(sp.series(sp.atan(tser) / sp.pi, x, 0, order + 1).removeO())
    return sp.collect(mu, x)


def match_mu(lam_times_r, rZ):
    """mu = atan(t)/pi with t = (r f' - lam r f... ) using r = rZ (Z = 1 units) and lam_times_r = r u'/u."""
    fN = lambda r: mp.sqrt(2 * r) * mp.besselj(1, mp.sqrt(8 * r))
    gN = lambda r: -mp.sqrt(2 * r) * mp.bessely(1, mp.sqrt(8 * r))
    rf = rZ * mp.diff(fN, rZ)
    rg = rZ * mp.diff(gN, rZ)
    t = (rf - lam_times_r * fN(rZ)) / (lam_times_r * gN(rZ) - rg)
    return mp.atan(t) / mp.pi


def shell_exact(xv):
    w = mp.sqrt(2 * xv)
    return match_mu(w * mp.cot(w), xv)


def gauss_exact(yv, smax=12):
    """Integrate u'' = -2 y erf(s/sqrt2)/s u in s = r/b from 0 to smax; Z = 1."""
    def V(s):
        if s == 0:
            return mp.sqrt(2 / mp.pi)
        return mp.erf(s / mp.sqrt(2)) / s
    sol = mp.odefun(lambda s, Y: [Y[1], -2 * yv * V(s) * Y[0]], 0, [mp.mpf(0), mp.mpf(1)])
    u, du = sol(smax)
    return match_mu(smax * du / u, smax * yv)


def fit(xs, vals, cols):
    A = mp.matrix([[c(t) for c in cols] for t in xs])
    sol = mp.lu_solve(A.T * A, A.T * mp.matrix(vals))
    return [sol[i] for i in range(len(cols))], mp.norm(A * sol - mp.matrix(vals))


def run() -> None:
    mu = shell_series()
    OUT["shell_series"] = str(mu)
    gam = sp.EulerGamma
    c4 = sp.expand(mu.coeff(x, 4))
    pred4_universal = -sp.Rational(4, 9) * (Lg + 2 * gam)
    OUT["shell_NP1_coef_x2"] = str(mu.coeff(x, 2))
    OUT["shell_NP2_coef_x3"] = str(mu.coeff(x, 3))
    OUT["shell_NP3_x4_minus_universal"] = str(sp.simplify(c4 - pred4_universal))
    # exact check of the shell series
    mp.mp.dps = 50
    chk = {}
    fmu = sp.lambdify((x, Lg), mu, "mpmath")
    for xv in (mp.mpf(1) / 100, mp.mpf(1) / 300, mp.mpf(1) / 1000):
        ser = fmu(xv, mp.log(2 * xv))
        ex = shell_exact(xv)
        chk[mp.nstr(xv, 3)] = {"exact_minus_series": mp.nstr(ex - ser, 4),
                               "over_x6": mp.nstr((ex - ser) / xv**6, 5)}
    OUT["shell_exact_check"] = chk

    # Gaussian, numerically
    mp.mp.dps = 30
    ys = [mp.mpf(1) / 1000 * 2 ** (mp.mpf(j) / 2) for j in range(11)]  # 1e-3 .. 3.2e-2
    red = []
    for yv in ys:
        d = gauss_exact(yv)
        red.append(d / (-2 * yv**2))
    P = lambda k: (lambda t: t**k)
    PL = lambda k: (lambda t: t**k * mp.log(t))
    co0, r0 = fit(ys, red, [P(k) for k in range(6)])
    co1, r1 = fit(ys, red, [P(k) for k in range(6)] + [PL(2), PL(3)])
    a1_pred = -mp.mpf(16) / (3 * mp.sqrt(mp.pi))
    OUT["gauss"] = {"no_log": {"a0": mp.nstr(co0[0], 12), "a1": mp.nstr(co0[1], 12), "resid": mp.nstr(r0, 3)},
                    "with_logs": {"a0": mp.nstr(co1[0], 12), "a1": mp.nstr(co1[1], 12), "y2log": mp.nstr(co1[6], 8),
                                  "resid": mp.nstr(r1, 3)},
                    "a1_pred_-16/(3 sqrt pi)": mp.nstr(a1_pred, 12),
                    "a1_rel_err": mp.nstr(abs(co1[1] - a1_pred) / abs(a1_pred), 3),
                    "y2log_pred": 2, "y2log_rel_err": mp.nstr(abs(co1[6] - 2) / 2, 3)}
    print(json.dumps(OUT, indent=1))
    with open("debug/data/nuclear_profiles_door.json", "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=1)


if __name__ == "__main__":
    run()
