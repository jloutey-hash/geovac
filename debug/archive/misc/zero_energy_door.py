"""Zero-energy door: tests ZE1-ZE3 of docs/predictions/2026-10-18_zero_energy.md.

Registered in commit bd2435d1 before this script ran. Writes debug/data/zero_energy_door.json.
"""
from __future__ import annotations

import json
import random
from typing import Dict, List

import mpmath as mp
import sympy as sp

OUT: Dict[str, object] = {}
X = sp.symbols("x y z", real=True)
x, y, z = X
rr = sp.sqrt(x**2 + y**2 + z**2)
EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def eps(i, j, k):
    return EPS.get((i, j, k), 0)


def p(i, f):
    return -sp.I * sp.diff(f, X[i])


def Lop(i, f):
    return sum(eps(i, j, k) * X[j] * p(k, f) for j in range(3) for k in range(3) if eps(i, j, k))


def Aop(i, f, Z):
    pxL = sum(eps(i, j, k) * p(j, Lop(k, f)) for j in range(3) for k in range(3) if eps(i, j, k))
    Lxp = sum(eps(i, j, k) * Lop(j, p(k, f)) for j in range(3) for k in range(3) if eps(i, j, k))
    return (pxL - Lxp) / 2 - Z * X[i] / rr * f


KTERMS = 18
HARM = {(0, 0): sp.Integer(1), (1, 0): z, (1, 1): x + sp.I * y, (2, 0): 2 * z**2 - x**2 - y**2}


def ze1() -> Dict[str, object]:
    mp.mp.dps = 30
    random.seed(7)
    pts = [tuple(sp.Rational(random.randint(3, 12), 100) * s for s in (1, -1, 1)) for _ in range(3)]
    res: Dict[str, object] = {}
    for Z in (1, 2):
        for (l, m), P in HARM.items():
            # J_{2l+1}(sqrt(8 Z r)) / (sqrt(r) r^l) is (2Z)^(l+1/2) times the entire series
            # sum_k (-2Zr)^k / (k! (k+2l+1)!); truncated at K terms, exact to ~ (2Zr)^K/(K!)^2 at r ~ 0.2
            nu = 2 * l + 1
            g = sum((-2 * Z * rr) ** kk / (sp.factorial(kk) * sp.factorial(kk + nu)) for kk in range(KTERMS))
            psi = g * P
            Ax, Ay, Az = (Aop(i, psi, Z) for i in range(3))
            comm = Aop(0, Ay, Z) - Aop(1, Ax, Z)
            A2 = Aop(0, Ax, Z) + Aop(1, Ay, Z) + Aop(2, Az, Z) - Z**2 * psi
            AL = sum(Aop(i, Lop(i, psi), Z) for i in range(3))
            worst = 0
            for pt in pts:
                sub = dict(zip(X, pt))
                scale = abs(sp.N(psi.subs(sub), 30))
                for e in (comm, A2, AL):
                    v = abs(sp.N(e.subs(sub), 30)) / scale
                    worst = max(worst, float(v))
            res[f"Z{Z}_l{l}_m{m}"] = worst
    # control at E < 0: 2p (m = 1), E = -Z^2/8; [A_x, A_y] psi = -2i E L_z psi, A^2 psi = (Z^2/4) psi
    Z = 1
    psi = rr**0 * sp.exp(-Z * rr / 2) * (x + sp.I * y)
    E = -sp.Rational(Z * Z, 8)
    Ax, Ay, Az = (Aop(i, psi, Z) for i in range(3))
    comm = Aop(0, Ay, Z) - Aop(1, Ax, Z)
    A2 = Aop(0, Ax, Z) + Aop(1, Ay, Z) + Aop(2, Az, Z)
    sub = dict(zip(X, pts[0]))
    c_val = sp.N(comm.subs(sub), 30)
    c_pred = sp.N((-2 * sp.I * E * Lop(2, psi)).subs(sub), 30)
    a_val = sp.N(A2.subs(sub), 30)
    a_pred = sp.N((sp.Rational(Z * Z, 4) * psi).subs(sub), 30)
    res["control_2p_comm_nonzero"] = str(sp.N(abs(c_val), 10))
    res["control_2p_comm_matches_-2iE_Lz"] = float(abs(c_val - c_pred) / abs(c_pred))
    res["control_2p_A2_matches_Z2/4"] = float(abs(a_val - a_pred) / abs(a_pred))
    return res


# ---------------------------------------------------------------- ZE2
def R_bound(n: int, l: int, Z: int, r):
    xx = 2 * mp.mpf(Z) * r / n
    lognorm = 3 * mp.log(2 * mp.mpf(Z) / n) + mp.loggamma(n - l) - mp.log(2 * n) - mp.loggamma(n + l + 1)
    return mp.exp(lognorm / 2) * xx**l * mp.e ** (-xx / 2) * mp.laguerre(n - l - 1, 2 * l + 1, xx)


def target(l: int, Z: int, r):
    return 2 * mp.besselj(2 * l + 1, mp.sqrt(8 * Z * r)) ** 2 / r


def fit(xs: List[mp.mpf], vals: List[mp.mpf], powers: List[int]) -> List[mp.mpf]:
    A = mp.matrix([[xv ** pw for pw in powers] for xv in xs])
    sol = mp.lu_solve(A, mp.matrix(vals))
    return [sol[i] for i in range(len(powers))]


def ze2() -> Dict[str, object]:
    mp.mp.dps = 60
    Ns = [200, 250, 300, 400, 500, 650, 800, 1000]
    res: Dict[str, object] = {}
    for Z in (1, 2):
        for l in (0, 1, 2):
            for r in (mp.mpf("0.5"), mp.mpf(2), mp.mpf(10)):
                T = target(l, Z, r)
                vals = [R_bound(n, l, Z, r) ** 2 * mp.mpf(n) ** 3 / Z**2 / T - 1 for n in Ns]
                co = fit([1 / mp.mpf(n) for n in Ns], vals, list(range(len(Ns))))
                res[f"bound_Z{Z}_l{l}_r{r}"] = {"a0": mp.nstr(co[0], 3), "a1": mp.nstr(co[1], 3),
                                                "a2": mp.nstr(co[2], 8), "a3": mp.nstr(co[3], 3)}
    mp.mp.dps = 30
    ks = [mp.mpf(1) / 100, mp.mpf(1) / 200, mp.mpf(1) / 400, mp.mpf(1) / 800]
    for Z in (1, 2):
        for l in (0, 1):
            r = mp.mpf(2)
            T = target(l, Z, r)
            vals = [(2 / (mp.pi * k)) * (mp.coulombf(l, -Z / k, k * r) / r) ** 2 / T - 1 for k in ks]
            co = fit(ks, vals, [0, 1, 2, 3])
            res[f"cont_Z{Z}_l{l}_r2"] = {"a0": mp.nstr(co[0], 3), "a1(k)": mp.nstr(co[1], 3),
                                         "a2(k^2)": mp.nstr(co[2], 8)}
    return res


# ---------------------------------------------------------------- ZE3
def ze3() -> Dict[str, object]:
    xs = sp.symbols("x", positive=True)
    res: Dict[str, object] = {}
    for l in (0, 1, 2):
        nu = 2 * l + 1
        ser = sum((-1) ** k * xs ** (k + sp.Rational(nu, 2)) / (sp.factorial(k) * sp.factorial(k + nu)) for k in range(12))
        expr = sp.expand(ser**2 / xs**nu)
        coeffs = sp.Poly(sp.series(expr, xs, 0, 10).removeO(), xs).all_coeffs()
        res[f"l{l}_coeffs_rational"] = all(c.is_rational for c in coeffs)
        # cross-check the series against sympy's besselj at x = 3/7
        v1 = sp.N((sp.besselj(nu, 2 * sp.sqrt(sp.Rational(3, 7)))) ** 2 / sp.Rational(3, 7) ** nu, 30)
        v2 = sp.N(expr.subs(xs, sp.Rational(3, 7)), 30)
        res[f"l{l}_series_vs_besselj"] = float(abs(v1 - v2))
    return res


def run() -> None:
    OUT["ZE1"] = ze1()
    print("ZE1", json.dumps(OUT["ZE1"], indent=1))
    OUT["ZE2"] = ze2()
    print("ZE2", json.dumps(OUT["ZE2"], indent=1))
    OUT["ZE3"] = ze3()
    print("ZE3", OUT["ZE3"])
    with open("debug/data/zero_energy_door.json", "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=1, default=str)


if __name__ == "__main__":
    run()
