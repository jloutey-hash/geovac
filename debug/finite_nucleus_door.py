"""Finite-nucleus door: tests FN0-FN4 of docs/predictions/2026-10-18_finite_nucleus.md.

Registered in commit bd2435d1 (amendment 5ba129aa) before this script ran.
Writes debug/data/finite_nucleus_door.json.

Uniformly charged sphere of radius R: inside, the electron's radial equation is an isotropic
oscillator (omega^2 = Z/R^3, energy shifted by 3Z/(2R)), solved by a Kummer function of r^2;
outside it is pure Coulomb (Whittaker W for bound states, F and G for the continuum). Every
quantity below is one log-derivative matching at r = R.
"""
from __future__ import annotations

import json
from typing import Dict, List

import mpmath as mp
import sympy as sp

OUT: Dict[str, object] = {}


def L_in(E, l: int, Z, R):
    """u'/u at r = R for the regular interior solution."""
    w = mp.sqrt(Z / R**3)
    eps = E + 3 * Z / (2 * R)
    b = l + mp.mpf(3) / 2
    a = b / 2 - eps / (2 * w)
    xx = w * R**2
    M = mp.hyp1f1(a, b, xx)
    dM = a / b * mp.hyp1f1(a + 1, b + 1, xx)
    return (l + 1) / R - w * R + 2 * w * R * dM / M


def L_out_bound(nu, l: int, Z, R):
    """u'/u at r = R for W_{nu, l+1/2}(2 Z r / nu); x W' = (x/2 - kappa) W - W_{kappa+1}."""
    xx = 2 * Z * R / nu
    W = mp.whitw(nu, l + mp.mpf(1) / 2, xx)
    W1 = mp.whitw(nu + 1, l + mp.mpf(1) / 2, xx)
    dWdx = ((xx / 2 - nu) * W - W1) / xx
    return dWdx * (2 * Z / nu) / W


def nu_bound(n: int, l: int, Z, R):
    f = lambda nu: L_in(-Z**2 / (2 * nu**2), l, Z, R) - L_out_bound(nu, l, Z, R)
    d0 = lead(n, l, 1, Z * R)  # start the secant at the leading-order estimate (the function is steep)
    if d0 == 0:
        d0 = -mp.mpf(10) ** -30
    return mp.findroot(f, (mp.mpf(n) - d0, mp.mpf(n) - mp.mpf("1.001") * d0))


def phase(k, l: int, Z, R):
    """Short-range phase shift relative to Coulomb: u = F + tan(Delta) G outside."""
    E = k**2 / 2
    Li = L_in(E, l, Z, R)
    eta = -Z / k
    F = mp.coulombf(l, eta, k * R)
    G = mp.coulombg(l, eta, k * R)
    dF = k * mp.diff(lambda rho: mp.coulombf(l, eta, rho), k * R)
    dG = k * mp.diff(lambda rho: mp.coulombg(l, eta, rho), k * R)
    t = (dF - Li * F) / (Li * G - dG)
    return mp.atan(t)


def c_l(l: int) -> sp.Rational:
    k = sp.Rational(1, 2 * l + 2) - sp.Rational(3, 2 * (2 * l + 3)) + sp.Rational(1, 2 * (2 * l + 5))
    return 4 * 2 ** (2 * l) * k / sp.factorial(2 * l + 1) ** 2


def lead(n_or_E, l: int, Z, x, bound: bool = True):
    """Registered leading defect -c_l x^(2l+2) prod(1 - s^2/n^2), or its continuation in E."""
    c = mp.mpf(sp.Rational(c_l(l)).p) / sp.Rational(c_l(l)).q
    if bound:
        pr = mp.fprod([1 - mp.mpf(s * s) / n_or_E**2 for s in range(1, l + 1)])
    else:
        pr = mp.fprod([1 + 2 * mp.mpf(s * s) * n_or_E / Z**2 for s in range(1, l + 1)])
    return -c * x ** (2 * l + 2) * pr


def fit(xs: List, vals: List, cols) -> List:
    A = mp.matrix([[cf(xv) for cf in cols] for xv in xs])
    sol = mp.lu_solve(A.T * A, A.T * mp.matrix(vals))
    resid = mp.norm(A * sol - mp.matrix(vals))
    return [sol[i] for i in range(len(cols))], resid


def run() -> None:
    mp.mp.dps = 50
    # FN0: reduced mass, symbolic
    n_, Z_, M_ = sp.symbols("n Z M", positive=True)
    mu = M_ / (M_ + 1)
    nu_eff = Z_ / sp.sqrt(-2 * (-mu * Z_**2 / (2 * n_**2)) / mu)  # in units where the focal length absorbs mu
    OUT["FN0_reduced_mass_defect"] = str(sp.simplify(n_ - nu_eff))
    OUT["c_l"] = {l: str(c_l(l)) for l in range(3)}

    # FN2: defects vs ZR
    xs = [mp.mpf(1) / 1000 * 2**j for j in range(6)]  # ZR = 1e-3 .. 3.2e-2
    fn2: Dict[str, object] = {}
    scaling_worst = mp.mpf(0)
    for l in (0, 1, 2):
        for n in (l + 1, l + 2, 5, 10, 20, 40):
            if n <= l:
                continue
            ds = []
            for xv in xs:
                nu = nu_bound(n, l, mp.mpf(1), xv)
                ds.append(n - nu)
                if n == 5 and xv == xs[2]:  # Z-scaling check: Z = 2, R = x/2
                    nu2 = nu_bound(n, l, mp.mpf(2), xv / 2)
                    scaling_worst = max(scaling_worst, abs((n - nu2) - (n - nu)) / abs(n - nu))
            red = [d / lead(n, l, 1, xv) for d, xv in zip(ds, xs)]  # -> 1 + a1 x + ...
            cols = [lambda t: 1, lambda t: t, lambda t: t**2, lambda t: t**3]
            co, res = fit(xs, red, cols)
            colsL = [lambda t: 1, lambda t: t, lambda t: t**2, lambda t: t**2 * mp.log(t), lambda t: t**3]
            coL, resL = fit(xs, red, colsL)
            fn2[f"l{l}_n{n}"] = {"delta_at_ZR=1e-3": mp.nstr(ds[0], 12), "lead_ratio_a0": mp.nstr(co[0], 12),
                                 "a1": mp.nstr(co[1], 12), "a2": mp.nstr(co[2], 8),
                                 "fit_resid": mp.nstr(res, 3), "with_x2logx_coef": mp.nstr(coL[3], 6),
                                 "with_log_resid": mp.nstr(resL, 3)}
    fn2["Z_scaling_rel_diff(Z=2,R/2 vs Z=1,R)"] = mp.nstr(scaling_worst, 3)
    OUT["FN2"] = fn2

    # FN1: node count unchanged at small ZR; the exterior W is multiplied by the sign that makes it
    # continuous with the interior at r = R (without that, a convention sign flip counts as a node)
    def nodes(n, l, xv):
        nu = nu_bound(n, l, mp.mpf(1), xv)
        E = -1 / (2 * nu**2)
        w = mp.sqrt(1 / xv**3)
        eps = E + 3 / (2 * xv)
        b = l + mp.mpf(3) / 2
        inside = lambda r: mp.hyp1f1(b / 2 - eps / (2 * w), b, w * r * r) * r ** (l + 1) * mp.e ** (-w * r * r / 2)
        outside = lambda r: mp.whitw(nu, l + mp.mpf(1) / 2, 2 * r / nu)
        sgn = mp.sign(inside(xv)) * mp.sign(outside(xv))
        cnt, prev = 0, None
        for i in range(1, 6000):
            r = mp.mpf(i) * (4 * n * n) / 6000
            v = inside(r) if r < xv else sgn * outside(r)
            s_ = mp.sign(v)
            if prev is not None and s_ != prev and s_ != 0:
                cnt += 1
            prev = s_
        return cnt
    mp.mp.dps = 25
    OUT["FN1_nodes"] = {f"n{n}_l{l}": [nodes(n, l, mp.mpf(1) / 100), n - l - 1]
                        for n, l in ((1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2), (4, 1), (5, 0))}
    mp.mp.dps = 60

    # FN4: the s-wave defect beyond leading order, with and without log columns (n = 1 and 5)
    xs4 = [mp.mpf(1) / 4000 * 2 ** (mp.mpf(j) / 2) for j in range(14)]
    P = lambda k: (lambda t: t**k)
    PL = lambda k: (lambda t: t**k * mp.log(t))
    fn4: Dict[str, object] = {}
    for n in (1, 5):
        red = [(n - nu_bound(n, 0, mp.mpf(1), xv)) / lead(n, 0, 1, xv) for xv in xs4]
        for name, cols in (("powers_0..7_no_log", [P(k) for k in range(8)]),
                           ("powers_0..6_plus_x2log_x3log_x4log", [P(k) for k in range(7)] + [PL(2), PL(3), PL(4)])):
            co, res = fit(xs4, red, cols)
            fn4[f"n{n}_{name}"] = {"a0": mp.nstr(co[0], 12), "a1": mp.nstr(co[1], 12), "resid": mp.nstr(res, 3),
                                   "x2log_coef": mp.nstr(co[7], 8) if len(cols) > 8 else None}
    fn4["-80/63"] = mp.nstr(mp.mpf(-80) / 63, 12)
    OUT["FN4"] = fn4
    mp.mp.dps = 50

    # FN3: phase shifts; Seaton's limit; leading form with and without J13's factor
    fn3: Dict[str, object] = {}
    xv = mp.mpf(1) / 100
    for l in (0, 1):
        ns = [10, 14, 20, 28, 40]
        dn = [n - nu_bound(n, l, mp.mpf(1), xv) for n in ns]
        cob, _ = fit(ns, dn, [lambda t: 1, lambda t: 1 / mp.mpf(t) ** 2, lambda t: 1 / mp.mpf(t) ** 4])
        ks = [mp.mpf(1) / 20, mp.mpf(1) / 30, mp.mpf(1) / 40, mp.mpf(1) / 60]
        ph = [phase(k, l, mp.mpf(1), xv) / mp.pi for k in ks]
        coc, _ = fit(ks, ph, [lambda t: 1, lambda t: t**2, lambda t: t**4])
        fn3[f"l{l}_Seaton"] = {"delta_inf_bound": mp.nstr(cob[0], 20), "Delta0/pi_cont": mp.nstr(coc[0], 20),
                               "rel_diff": mp.nstr(abs(cob[0] - coc[0]) / abs(cob[0]), 3)}
        rows = {}
        for k in (mp.mpf(1) / 2, mp.mpf(1), mp.mpf(2)):
            E = k**2 / 2
            eta = 1 / k
            d = phase(k, l, mp.mpf(1), mp.mpf(1) / 1000) / mp.pi
            ld = lead(E, l, 1, mp.mpf(1) / 1000, bound=False)
            fac = 1 / (1 - mp.e ** (-2 * mp.pi * eta))
            rows[f"k{k}"] = {"ratio_registered(no factor)": mp.nstr(d / ld, 10),
                             "ratio_with_J13_factor": mp.nstr(d / (ld * fac), 10)}
        fn3[f"l{l}_leading_ZR=1e-3"] = rows
    OUT["FN3"] = fn3

    # second route for FN2 leading term: first-order perturbation by quadrature, n = 2, l = 0, 1
    def first_order(n, l, xv):
        r_ = sp.symbols("r", positive=True)
        from sympy.physics.hydrogen import R_nl
        Rn = sp.lambdify(r_, R_nl(n, l, r_, 1), "mpmath")
        dV = lambda r: 1 / r - (3 * xv**2 - r**2) / (2 * xv**3)
        dE = mp.quad(lambda r: dV(r) * Rn(r) ** 2 * r**2, [0, xv])
        return -dE * n**3
    OUT["FN2_route2_first_order"] = {f"n{n}_l{l}": mp.nstr(first_order(n, l, mp.mpf(1) / 1000) / lead(n, l, 1, mp.mpf(1) / 1000), 12)
                                     for n, l in ((2, 0), (5, 0), (2, 1), (5, 1))}

    print(json.dumps(OUT, indent=1, default=str))
    with open("debug/data/finite_nucleus_door.json", "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=1, default=str)


if __name__ == "__main__":
    run()
