"""Hydrogen-continuum door: tests HC1-HC4 of docs/predictions/2026-10-18_hydrogen_continuum.md.

Registered in commit fb4fdd85 before this script ran. Writes debug/data/hydrogen_continuum_door.json.

Route A: the bound-state oscillator-strength density Q_i(n) = f(i -> n) dn/dE as an analytic
function of n (Laplace transform of the confluent-hypergeometric radial function, a terminating
2F1 at integer n), validated against exact rationals at integer n, then evaluated at n = i eta.
Route B: direct radial quadrature against mpmath's regular Coulomb function, energy-normalized.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp

TRANS = [((1, 0), 1), ((2, 0), 1), ((2, 1), 0), ((2, 1), 2)]  # (n_i, l_i) -> l
OUT: Dict[str, object] = {}


# ---------------------------------------------------------------- HC1
def hc1() -> Dict[str, object]:
    n, chi, rho = sp.symbols("n chi rho", positive=True)
    f = sp.sin(n * chi) / sp.sin(chi)
    lap_s3 = sp.diff(sp.sin(chi) ** 2 * sp.diff(f, chi), chi) / sp.sin(chi) ** 2
    ok_s3 = sp.simplify(lap_s3 + (n**2 - 1) * f) == 0
    g = sp.sin(rho * chi) / sp.sinh(chi)
    lap_h3 = sp.diff(sp.sinh(chi) ** 2 * sp.diff(g, chi), chi) / sp.sinh(chi) ** 2
    ok_h3 = sp.simplify(lap_h3 + (rho**2 + 1) * g) == 0
    cont = sp.simplify((f.subs(n, sp.I * rho).subs(chi, sp.I * chi)) / g)
    cas = sp.expand((n**2 - 1).subs(n, sp.I * rho))
    en = sp.simplify((-sp.Symbol("Z") ** 2 / (2 * n**2)).subs(n, sp.I * sp.Symbol("Z") / sp.Symbol("k")))
    return {"S3_zonal_eigen_n2-1": bool(ok_s3), "H3_zonal_eigen_1+rho2": bool(ok_h3),
            "continued_S3_zonal_over_H3_zonal": str(cont), "casimir_continued": str(cas),
            "E_continued": str(en)}


# ---------------------------------------------------------------- bound radial functions
def initial_coeffs(ni: int, li: int, Z: int) -> Tuple[mp.mpf, List[Fraction]]:
    """R_i(r) = norm * sum_j a_j r^(li+j) e^(-Z r/ni); returns (norm^2 as Fraction-free mpf, a_j)."""
    k = ni - li - 1
    alpha = 2 * li + 1
    c = Fraction(2 * Z, ni)
    a = []
    for j in range(k + 1):
        lag = Fraction((-1) ** j * comb(k + alpha, k - j), factorial(j))
        a.append(lag * c ** (li + j))
    norm2 = c**3 * Fraction(factorial(ni - li - 1), 2 * ni * factorial(ni + li))
    return norm2, a


def A2(n, l: int, Z: int):
    """A_l(n)^2 = (2Z/n)^3 P_l(n)/(2n) / ((2l+1)!)^2, P_l(n) = prod_{s=-l..l}(n+s); analytic in n."""
    P = mp.mpf(1)
    for s in range(-l, l + 1):
        P *= n + s
    return (2 * Z / n) ** 3 * P / (2 * n) / mp.factorial(2 * l + 1) ** 2


def h2f1(a, b: int, c: int, z, n, ni: int):
    """2F1(a, b; c; z) by Euler's transformation, (1-z)^(c-a-b) 2F1(c-a, c-b; c; z). For every
    dipole transition tested c - b <= 0, so the second factor is a finite polynomial and the whole
    is one closed form, analytic in n; 1 - z = (n - ni)/(n + ni), principal branch."""
    mneg = c - b
    assert mneg <= 0
    poly, term = mp.mpf(1), mp.mpf(1)
    for j in range(-mneg):
        term *= (c - a + j) * (c - b + j) / mp.mpf((c + j) * (j + 1)) * z
        poly += term
    return mp.power((n - ni) / (n + ni), c - a - b) * poly


def I2_routeA(n, ni: int, li: int, l: int, Z: int):
    """I(n)^2 for I = int R_i R_{n l} r^3 dr, analytic in n (no square roots)."""
    norm2, a = initial_coeffs(ni, li, Z)
    beta = mp.mpf(Z) / ni + Z / n
    c = 2 * mp.mpf(Z) / n
    tot = 0
    for j, aj in enumerate(a):
        s = li + j + l + 4
        tot += mp.mpf(aj.numerator) / aj.denominator * mp.gamma(s) * beta ** (-s) * h2f1(l + 1 - n, s, 2 * l + 2, c / beta, n, ni)
    return mp.mpf(norm2.numerator) / norm2.denominator * A2(n, l, Z) * (2 * mp.mpf(Z) / n) ** (2 * l) * tot**2


def Q_routeA(n, ni: int, li: int, l: int, Z: int):
    Ei = -mp.mpf(Z) ** 2 / (2 * ni**2)
    En = -mp.mpf(Z) ** 2 / (2 * n**2)
    f = mp.mpf(2) / 3 * (En - Ei) * max(l, li) / (2 * li + 1) * I2_routeA(n, ni, li, l, Z)
    return f * n**3 / Z**2


def f_exact(n: int, ni: int, li: int, l: int, Z: int) -> Fraction:
    """Exact oscillator strength from explicit Laguerre polynomials (independent of the 2F1 route)."""
    norm2_i, a = initial_coeffs(ni, li, Z)
    norm2_f, b = initial_coeffs(n, l, Z)
    beta = Fraction(Z, ni) + Fraction(Z, n)
    tot = Fraction(0)
    for j, aj in enumerate(a):
        for q, bq in enumerate(b):
            p = li + j + l + q + 3  # int r^p e^(-beta r) = p!/beta^(p+1)
            tot += aj * bq * Fraction(factorial(p)) / beta ** (p + 1)
    I2 = norm2_i * norm2_f * tot * tot
    dE = Fraction(Z * Z, 2 * ni * ni) - Fraction(Z * Z, 2 * n * n)
    return Fraction(2, 3) * dE * Fraction(max(l, li), 2 * li + 1) * I2


# ---------------------------------------------------------------- route B
def dfdE_routeB(k, ni: int, li: int, l: int, Z: int):
    norm2, a = initial_coeffs(ni, li, Z)
    eta_s = -mp.mpf(Z) / k

    def integrand(r):
        poly = sum(mp.mpf(aj.numerator) / aj.denominator * r ** (li + j) for j, aj in enumerate(a))
        return poly * mp.e ** (-Z * r / ni) * mp.coulombf(l, eta_s, k * r) * r**2

    L = mp.mpf(ni) / Z
    pts = [0] + [L * t for t in (0.5, 1, 2, 4, 8, 16, 32, 64, 128, 200)]
    I = mp.sqrt(mp.mpf(norm2.numerator) / norm2.denominator) * mp.sqrt(2 / (mp.pi * k)) * mp.quad(integrand, pts)
    E = k**2 / 2
    Ei = -mp.mpf(Z) ** 2 / (2 * ni**2)
    return mp.mpf(2) / 3 * (E - Ei) * max(l, li) / (2 * li + 1) * I**2


def contact_routeB(k, Z: int):
    eta_s = -mp.mpf(Z) / k
    rho = mp.mpf(10) ** -30
    C0 = mp.coulombf(0, eta_s, rho) / rho
    return (2 / (mp.pi * k)) * (C0 * k) ** 2 / (4 * mp.pi)


# ---------------------------------------------------------------- HC2, HC3, HC4
def run() -> None:
    OUT["HC1"] = hc1()
    print("HC1", OUT["HC1"], flush=True)

    n, Z = sp.symbols("n Z", positive=True)
    bound_contact = sp.simplify((4 * Z**3 / n**3) / (4 * sp.pi) * n**3 / Z**2)
    OUT["HC2_bound_symbolic"] = str(bound_contact)

    mp.mp.dps = 40
    # route A at integer n vs exact rationals; 1s also vs the registered closed form
    valid: Dict[str, object] = {}
    for (ni, li), l in TRANS:
        worst = mp.mpf(0)
        for nn in range(ni + 1, 61):
            fe = f_exact(nn, ni, li, l, 1)
            fa = Q_routeA(mp.mpf(nn), ni, li, l, 1) / mp.mpf(nn) ** 3
            fev = mp.mpf(fe.numerator) / fe.denominator
            worst = max(worst, abs(fa - fev) / fev)
        valid[f"{ni}{'sp'[li]}->l{l}"] = mp.nstr(worst, 3)
    closed_ok = all(
        f_exact(nn, 1, 0, 1, 1) == Fraction(2**8 * nn**5 * (nn - 1) ** (2 * nn - 4), 3 * (nn + 1) ** (2 * nn + 4))
        for nn in range(2, 61))
    OUT["routeA_vs_exact_integer_n_2_60_max_rel"] = valid
    OUT["1s_registered_closed_form_exact_n2_60"] = closed_ok
    print("valid", valid, closed_ok, flush=True)

    mp.mp.dps = 30
    rows = []
    worst_univ = mp.mpf(0)
    worst_imag = mp.mpf(0)
    worst_contact = mp.mpf(0)
    for Zv in (1, 2):
        for k in (mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf(1), mp.mpf(2), mp.mpf(5)):
            eta = Zv / k
            factor = 1 / (1 - mp.e ** (-2 * mp.pi * eta))
            cb = contact_routeB(k, Zv)
            dev_c = abs(cb / (mp.mpf(Zv) / mp.pi * factor) - 1)
            worst_contact = max(worst_contact, dev_c)
            for (ni, li), l in TRANS:
                qa = Q_routeA(mp.mpc(0, eta), ni, li, l, Zv)
                worst_imag = max(worst_imag, abs(mp.im(qa)) / abs(qa))
                vb = dfdE_routeB(k, ni, li, l, Zv)
                ratio = vb / mp.re(qa)
                dev = abs(ratio / factor - 1)
                worst_univ = max(worst_univ, dev)
                rows.append({"Z": Zv, "k": mp.nstr(k, 3), "trans": f"{ni}{'sp'[li]}->l{l}",
                             "routeB": mp.nstr(vb, 20), "QA(i eta)": mp.nstr(mp.re(qa), 20),
                             "ratio/factor-1": mp.nstr(dev, 3)})
                print(rows[-1], flush=True)
    OUT["HC3_rows"] = rows
    OUT["HC3_max_|ratio/factor-1|"] = mp.nstr(worst_univ, 3)
    OUT["HC3_max_rel_imag_QA"] = mp.nstr(worst_imag, 3)
    OUT["HC2_max_|contact/(Z/pi factor)-1|"] = mp.nstr(worst_contact, 3)

    # thresholds: bound side at large n vs continuum as k -> 0, and 1s against 2^8 e^-4/3
    mp.mp.dps = 40
    thr = {}
    for (ni, li), l in TRANS:
        cell = {}
        for big in (10**2, 10**3):
            qb = Q_routeA(mp.mpf(big), ni, li, l, 1)
            eta = mp.mpf(big)
            qc = mp.re(Q_routeA(mp.mpc(0, eta), ni, li, l, 1)) / (1 - mp.e ** (-2 * mp.pi * eta))
            cell[f"bound_n={big}"] = mp.nstr(qb, 15)
            cell[f"cont_eta={big}"] = mp.nstr(qc, 15)
        thr[f"{ni}{'sp'[li]}->l{l}"] = cell
    thr["2^8 e^-4/3"] = mp.nstr(mp.mpf(2) ** 8 * mp.e ** -4 / 3, 15)
    OUT["thresholds"] = thr
    print("thr", thr, flush=True)

    # HC4: at eta = 1, 2, strip the factor and the predicted exponential; the rest should be rational
    mp.mp.dps = 50
    hc4 = {}
    for (ni, li), l in TRANS:
        for eta in (1, 2):
            k = mp.mpf(1) / eta
            vb = dfdE_routeB(k, ni, li, l, 1)
            t = vb * (1 - mp.e ** (-2 * mp.pi * eta)) / mp.e ** (-4 * eta * mp.acot(mp.mpf(eta) / ni))
            ident = mp.identify(t, tol=mp.mpf(10) ** -35)
            argG = mp.im(mp.loggamma(l + 1 + 1j * eta))
            rel = mp.pslq([t, 1, argG, mp.log(k) if eta != 1 else mp.log(3), mp.euler],
                          tol=mp.mpf(10) ** -35, maxcoeff=10**15, maxsteps=10**6)
            hc4[f"{ni}{'sp'[li]}->l{l}_eta{eta}"] = {"t": mp.nstr(t, 30), "identify": str(ident), "pslq": str(rel)}
            print(hc4[f"{ni}{'sp'[li]}->l{l}_eta{eta}"], flush=True)
    ctrl = mp.mpf(7) / 13 + 2 * mp.euler
    argG1 = mp.im(mp.loggamma(2 + 1j))
    hc4["control_7/13+2gamma"] = str(mp.pslq([ctrl, 1, argG1, mp.log(3), mp.euler],
                                            tol=mp.mpf(10) ** -35, maxcoeff=10**15, maxsteps=10**6))
    OUT["HC4"] = hc4
    with open("debug/data/hydrogen_continuum_door.json", "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print("DONE")


if __name__ == "__main__":
    run()
