"""Large-spin door: tests LS1-LS5 of docs/predictions/2026-10-18_large_spin.md.

Registered in commit fb4fdd85 before this script ran. Writes debug/data/large_spin_door.json.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp

OUT: Dict[str, object] = {}


# ---------------------------------------------------------------- LS1
def ls1() -> Dict[str, object]:
    S, n = sp.symbols("S n", positive=True)
    m = S - n
    elem = sp.sqrt((S + m) * (S - m + 1))            # S_- |S,m>, m = S - n
    scaled = elem / sp.sqrt(2 * S)
    target = sp.sqrt(n + 1) * sp.sqrt(1 - n / (2 * S))
    sym_ok = sp.simplify(scaled**2 - target**2) == 0
    series = sp.series(1 - sp.sqrt(1 - n / (2 * S)), S, sp.oo, 3).removeO()
    rate_ok = sp.simplify(series - (n / (4 * S) + n**2 / (32 * S**2))) == 0
    comm_ok = True
    for Sv in range(1, 21):
        dim = 2 * Sv + 1
        Am = sp.zeros(dim, dim)  # A_S = S_+/sqrt(2S) in the n basis
        for k in range(dim - 1):
            Am[k, k + 1] = sp.sqrt(sp.Rational((2 * Sv - k) * (k + 1), 2 * Sv))
        C = Am * Am.T - Am.T * Am
        for k in range(dim):
            if sp.simplify(C[k, k] - (1 - sp.Rational(k, Sv))) != 0:
                comm_ok = False
        if any(C[i, j] != 0 for i in range(dim) for j in range(dim) if i != j):
            comm_ok = False
    return {"matrix_element_symbolic": bool(sym_ok), "rate_n_over_4S": bool(rate_ok),
            "commutator_1_minus_n_over_S_S1_to_20": comm_ok}


# ---------------------------------------------------------------- d-matrix
def d_half_pi_sq_parts(S: int, mp_: int, m: int) -> Tuple[Fraction, int]:
    """d^S_{mp,m}(pi/2) = sign * sqrt(F) exactly; return (F, sign). Wigner's sum, integer S."""
    pref = factorial(S + mp_) * factorial(S - mp_) * factorial(S + m) * factorial(S - m)
    tot = Fraction(0)
    kmin, kmax = max(0, m - mp_), min(S + m, S - mp_)
    for k in range(kmin, kmax + 1):
        den = factorial(S + m - k) * factorial(k) * factorial(S - k - mp_) * factorial(k - m + mp_)
        tot += Fraction((-1) ** (k - m + mp_), den)
    val_sq = Fraction(pref) * tot * tot / Fraction(4**S)  # cos^a sin^b at pi/2, a + b = 2S
    sign = 1 if tot > 0 else (-1 if tot < 0 else 0)
    return val_sq, sign


def ls2() -> Dict[str, object]:
    mp.mp.dps = 50
    ok_eig = True
    worst = mp.mpf(0)
    for Sv in (3, 7, 12):
        dim = 2 * Sv + 1
        ms = list(range(Sv, -Sv - 1, -1))
        Sx = mp.zeros(dim, dim)
        for i, mv in enumerate(ms[:-1]):  # S_- : m -> m-1
            el = mp.sqrt((Sv + mv) * (Sv - mv + 1))
            Sx[i + 1, i] += el / 2
            Sx[i, i + 1] += el / 2
        for mcol in ms:
            vec = mp.matrix(dim, 1)
            for i, mrow in enumerate(ms):
                F, sg = d_half_pi_sq_parts(Sv, mrow, mcol)
                vec[i] = sg * mp.sqrt(mp.mpf(F.numerator) / F.denominator)
            r = mp.norm(Sx * vec - mcol * vec)
            worst = max(worst, r)
            if r > mp.mpf(10) ** -40:
                ok_eig = False
    return {"d_columns_are_Sx_eigenvectors_S3_7_12": ok_eig, "worst_residual": mp.nstr(worst, 5),
            "squares_rational": "by construction (exact Fraction sum)"}


# ---------------------------------------------------------------- LS3
XS = [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(5, 2), Fraction(3)]


def p_exact(N: int, x: Fraction) -> Fraction:
    """p_S(x) = sqrt(S) C(2S, S+m)/4^S at S = N^2, m = x N; exactly rational."""
    S = N * N
    mf = x * N
    assert mf.denominator == 1
    m = int(mf)
    return Fraction(N * comb(2 * S, S + m), 4**S)


def fit_coeffs(Ns: List[int], vals: List[mp.mpf], nterms: int, half: bool) -> List[mp.mpf]:
    """Least-squares fit vals ~ sum_j a_j S^(-e_j); integer exponents, or half-integers too."""
    exps = [mp.mpf(j) / 2 for j in range(nterms)] if half else [mp.mpf(j) for j in range(nterms)]
    A = mp.matrix(len(Ns), len(exps))
    b = mp.matrix(len(Ns), 1)
    for i, N in enumerate(Ns):
        S = mp.mpf(N) ** 2
        for j, e in enumerate(exps):
            A[i, j] = S ** (-e)
        b[i] = vals[i]
    sol = mp.lu_solve(A.T * A, A.T * b)
    return [sol[j] for j in range(len(exps))]


def ls3() -> Dict[str, object]:
    mp.mp.dps = 150
    Ns = list(range(40, 401, 20))
    res: Dict[str, object] = {}
    c1_ok, c2_vals = True, {}
    half_max = mp.mpf(0)
    for x in XS:
        xm = mp.mpf(x.numerator) / x.denominator
        g = mp.e ** (-xm**2) / mp.sqrt(mp.pi)
        vals = []
        for N in Ns:
            p = p_exact(N, x)
            vals.append(mp.mpf(p.numerator) / p.denominator / g - 1)
        co = fit_coeffs(Ns, vals, 9, half=False)
        c1_pred = mp.mpf(-1) / 8 + xm**2 / 2 - xm**4 / 6
        err1 = abs(co[1] - c1_pred)
        if err1 > mp.mpf(10) ** -12:
            c1_ok = False
        cohalf = fit_coeffs(Ns, vals, 14, half=True)
        half_max = max(half_max, max(abs(cohalf[j]) for j in (1, 3, 5)))
        c2_vals[str(x)] = co[2]
        res[str(x)] = {"c0": mp.nstr(co[0], 5), "c1": mp.nstr(co[1], 25), "c1_pred": mp.nstr(c1_pred, 25),
                       "c1_err": mp.nstr(err1, 3), "c2": mp.nstr(co[2], 25)}
    rat = {xs: str(mp.identify(v, tol=mp.mpf(10) ** -20)) for xs, v in c2_vals.items()}
    A = mp.matrix(len(XS), 5)
    b = mp.matrix(len(XS), 1)
    for i, x in enumerate(XS):
        xm = mp.mpf(x.numerator) / x.denominator
        for j in range(5):
            A[i, j] = xm ** (2 * j)
        b[i] = c2_vals[str(x)]
    sol = mp.lu_solve(A.T * A, A.T * b)
    resid = mp.norm(A * sol - b)
    poly = [str(mp.identify(sol[j], tol=mp.mpf(10) ** -18)) for j in range(5)]
    v = c2_vals["1"]
    rel = mp.pslq([v, 1, mp.log(2), mp.euler, mp.zeta(3)], tol=mp.mpf(10) ** -25, maxcoeff=10**4, maxsteps=10**6)
    planted = mp.mpf(37) / 11 + 5 * mp.log(2)
    rel_ctrl = mp.pslq([planted, 1, mp.log(2), mp.euler, mp.zeta(3)], tol=mp.mpf(10) ** -25,
                       maxcoeff=10**4, maxsteps=10**6)
    res["c1_all_match_1e-12"] = c1_ok
    res["max_half_integer_coeff"] = mp.nstr(half_max, 5)
    res["c2_identified"] = rat
    res["c2_poly_in_x2_coeffs"] = poly
    res["c2_poly_residual"] = mp.nstr(resid, 5)
    res["c2_x1_pslq_with_log2_gamma_zeta3"] = str(rel)
    res["planted_control_37/11+5log2"] = str(rel_ctrl)
    p0 = p_exact(400, Fraction(0))
    res["p0_N400_minus_pi^-1/2"] = mp.nstr(mp.mpf(p0.numerator) / p0.denominator - 1 / mp.sqrt(mp.pi), 5)
    return res


# ---------------------------------------------------------------- LS4
def hermite_fn(n: int, x: mp.mpf) -> mp.mpf:
    return mp.hermite(n, x) * mp.e ** (-x**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))


def ls4() -> Dict[str, object]:
    mp.mp.dps = 80
    Ns = list(range(20, 161, 10))
    out: Dict[str, object] = {}
    for n in range(5):
        for x in (Fraction(1, 2), Fraction(1), Fraction(3, 2)):
            xm = mp.mpf(x.numerator) / x.denominator
            h = hermite_fn(n, xm)
            ratios = []
            for N in Ns:
                S = N * N
                m = int(x * N)
                F, sg = d_half_pi_sq_parts(S, S - n, m)
                comp = sg * mp.sqrt(mp.mpf(F.numerator) / F.denominator) * mp.sqrt(N)  # S^(1/4) = sqrt(N)
                ratios.append((-1) ** ((S - m) % 2) * comp / h)  # convention phase (-1)^(S-m)
            signs = sorted(set(1 if r > 0 else -1 for r in ratios))
            sgn = signs[0] if len(signs) == 1 else 0
            vals = [abs(r) - 1 for r in ratios]
            co = fit_coeffs(Ns, vals, 7, half=False)
            coh = fit_coeffs(Ns, vals, 11, half=True)
            out[f"n{n}_x{x}"] = {"sign": sgn, "limit_minus_1": mp.nstr(co[0], 5), "c1": mp.nstr(co[1], 15),
                                 "half_coeffs_S^-1/2_S^-3/2": [mp.nstr(coh[j], 3) for j in (1, 3)]}
    # route 2: numerical eigendecomposition of S_x at S = 30, vs closed-form d, n = 0..4
    mp.mp.dps = 50
    Sv = 30
    dim = 2 * Sv + 1
    ms = list(range(Sv, -Sv - 1, -1))
    Sx = mp.zeros(dim, dim)
    for i, mv in enumerate(ms[:-1]):
        el = mp.sqrt((Sv + mv) * (Sv - mv + 1))
        Sx[i + 1, i] += el / 2
        Sx[i, i + 1] += el / 2
    E, Q = mp.eigsy(Sx)
    worst = mp.mpf(0)
    for col in range(dim):
        mval = int(mp.nint(E[col]))
        for nrow in range(5):
            F, sg = d_half_pi_sq_parts(Sv, Sv - nrow, mval)
            closed = sg * mp.sqrt(mp.mpf(F.numerator) / F.denominator)
            num = Q[nrow, col]
            # fix the overall eigenvector sign from the n = 0 row
            F0, s0 = d_half_pi_sq_parts(Sv, Sv, mval)
            flip = 1 if (Q[0, col] * s0) > 0 else -1
            worst = max(worst, abs(flip * num - closed))
    out["route2_eigsy_S30_vs_closed_form_max_diff"] = mp.nstr(worst, 5)
    return out


# ---------------------------------------------------------------- LS5
def ls5() -> Dict[str, object]:
    S, w2 = sp.symbols("S w2", positive=True)
    meas = (2 * S + 1) / (2 * S) * (1 + w2 / (2 * S)) ** (-2 * S - 2)
    lim_meas = sp.limit(meas, S, sp.oo)
    ok_meas = sp.simplify(lim_meas - sp.exp(-w2)) == 0
    lims = []
    for k in range(6):
        e = sp.binomial(2 * S, k) / (2 * S) ** k
        lims.append(sp.simplify(sp.limit(sp.expand_func(e), S, sp.oo) - sp.Rational(1, math.factorial(k))) == 0)
    return {"measure_contracts_to_exp(-|w|^2)": bool(ok_meas), "components_to_w^n/sqrt(n!)_n0_5": all(lims)}


if __name__ == "__main__":
    OUT["LS1"] = ls1()
    print("LS1", OUT["LS1"], flush=True)
    OUT["LS2"] = ls2()
    print("LS2", OUT["LS2"], flush=True)
    OUT["LS5"] = ls5()
    print("LS5", OUT["LS5"], flush=True)
    OUT["LS3"] = ls3()
    print("LS3", json.dumps(OUT["LS3"], indent=1, default=str), flush=True)
    OUT["LS4"] = ls4()
    print("LS4", json.dumps(OUT["LS4"], indent=1, default=str), flush=True)
    with open("debug/data/large_spin_door.json", "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print("DONE")
