"""Guards for Paper 59 sec:scope, the equation of the three-centre integral in the bond angle (records
docs/predictions/2026-10-28_nuclear_attraction_angle_operator.md (NA), 2026-10-28_two_electron_bottom_up.md (B2),
2026-10-29_two_electron_large_l.md (B2X)).  Data: tests/data/paper59_angle_2026_10_29.json, built by
debug/build_angle_guard_data_2026_10_29.py (a copy of stored results, prime 2^61 - 1; no debug/ import here).

x = cos(gamma), f(x) = sum_L F_L P_L(x); the generic sequences replace the constants of the exact split by random
residues mod p.  Claims guarded, as the paper states them:
  * NA: the one-electron integral (case N1) obeys an operator sum a_ij x^j T^i, T = (1 - x^2) d/dx, of order 8 and
    degree 31; at x* = (D1^2 + D2^2)/(2 D1 D2) its exponents are 0..5, 3/2 and 5/2; every other rational singular
    point has integer exponents only;
  * B2: the two-electron (XY|XZ) (case W1) coefficients obey a recurrence in L of order 41 and degree 926; among its
    rational characteristic roots only 17/25 and 25/17 (both mapping to x*) carry integer beta;
  * B2X: the exact coefficients decay as (17/25)^L L^-4 (local exponent approaching -4, one sign), the exact rationals
    reproduce the generic sequence the recurrence annihilates, and their ratios match an independent quadrature.
Scope: rational singular points only (irreducible higher-degree factors of the leading coefficient are not analysed
here); one prime; the bottom-up order-9 operator of record BU and the 1970-unknown bound of record AO3 are
record-backed, not guarded here.  The decay leg works at L <= 101, so it cannot exclude an arbitrarily small beta = 3
weight: fire-tested, it catches a planted (17/25)^L L^3 piece at 16% of F_100 and misses one at 1.6%.  The zero weight
rests on record B2X (L = 300, beta_eff = -4.075).
Fire-tested 2026-10-29 (debug/qa/fire_test.py), every plant fired: T built from 1 + x^2; exponent 1/2 for 3/2; a
wrong x action; the recurrence shifted by one in L; the integer-beta root set narrowed; a beta = 3 piece; an
alternating (-1875/2057)^L piece; the exact data reduced with wrong weights; the quadrature with zX and zH swapped.
"""
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import nmod, nmod_poly
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

DATA = json.loads((Path(__file__).parent / "data" / "paper59_angle_2026_10_29.json").read_text())
P = DATA["p"]


def ratrec(a):
    r0, r1, s0, s1 = P, int(a) % P, 0, 1
    B = math.isqrt(P // 2)
    while r1 > B:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    return Fraction(r1 * (1 if s1 > 0 else -1), abs(s1))


def roots(poly):
    """Rational roots (with multiplicity) of an nmod_poly; non-linear factors as ('deg', d)."""
    out = []
    for g, m in poly.factor()[1]:
        out += [ratrec(-g.coeffs()[0] / g.coeffs()[1]) if g.degree() == 1 else ("deg", g.degree())] * m
    return out


# ---- NA: the order-8 operator --------------------------------------------------------------------------------------
def _X(s):
    return [(nmod(l, P) / nmod(2 * l - 1, P) * s[l - 1] if l else nmod(0, P)) + nmod(l + 1, P) / nmod(2 * l + 3, P)
            * s[l + 1] for l in range(len(s) - 1)]


def _T(s):
    return [nmod((l + 1) * (l + 2), P) / nmod(2 * l + 3, P) * s[l + 1]
            - (nmod((l - 1) * l, P) / nmod(2 * l - 1, P) * s[l - 1] if l else nmod(0, P)) for l in range(len(s) - 1)]


def na_residual(a):
    """Rows of sum a_ij x^j T^i applied to the N1 generic sequence (coefficient space)."""
    seq = [nmod(v, P) for v in DATA["na"]["gen"]]
    keys = [tuple(k) for k in DATA["na"]["keys"]]
    R, Q = max(i for i, _ in keys), max(j for _, j in keys)
    t, cols = [seq], {}
    for i in range(R + 1):
        if i:
            t.append(_T(t[-1]))
        c = t[i]
        for j in range(Q + 1):
            if j:
                c = _X(c)
            cols[(i, j)] = c
    n = min(len(c) for c in cols.values())
    return [sum((nmod(ak, P) * cols[k][l] for k, ak in zip(keys, a)), nmod(0, P)) for l in range(n)]


def standard_form(a):
    keys = [tuple(k) for k in DATA["na"]["keys"]]
    R = max(i for i, _ in keys)
    w = nmod_poly([1, 0, -1], P)
    ops = [{0: nmod_poly([1], P)}]
    for _ in range(R):
        new = {}
        for k, b in ops[-1].items():
            new[k] = new.get(k, nmod_poly([], P)) + w * b.derivative()
            new[k + 1] = new.get(k + 1, nmod_poly([], P)) + w * b
        ops.append(new)
    c = [nmod_poly([], P) for _ in range(R + 1)]
    for (i, j), ak in zip(keys, a):
        for k, b in ops[i].items():
            c[k] += nmod_poly([0] * j + [int(ak)], P) * b
    return c


def exponents_at(c, x0):
    sh = [ck(nmod_poly([int(nmod(x0.numerator, P) / nmod(x0.denominator, P)), 1], P)) for ck in c]
    val = [next((i for i, q in enumerate(s.coeffs()) if int(q)), None) for s in sh]
    delta = min(v - k for k, v in enumerate(val) if v is not None)
    ind = nmod_poly([], P)
    for k, s in enumerate(sh):
        if val[k] is not None and val[k] - k == delta:
            fall = nmod_poly([1], P)
            for t in range(k):
                fall *= nmod_poly([-t, 1], P)
            ind += fall * nmod_poly([int(s.coeffs()[val[k]])], P)
    return sorted(roots(ind), key=lambda r: (isinstance(r, tuple), r if not isinstance(r, tuple) else 0))


def test_na_operator_annihilates_the_one_electron_sequence():
    a = DATA["na"]["a"]
    assert max(i for i, _ in DATA["na"]["keys"]) == 8 and max(j for _, j in DATA["na"]["keys"]) == 31
    res = na_residual(a)
    assert len(res) > 300 and all(int(v) == 0 for v in res)
    bad = list(a)
    bad[100] = (bad[100] + 1) % P
    assert sum(1 for v in na_residual(bad) if int(v)) > 300


def test_na_branch_point_is_x_star_with_exponents_three_halves_and_five_halves():
    c = standard_form(DATA["na"]["a"])
    D1, D2 = (Fraction(s) for s in DATA["na"]["D"])
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    assert xs == Fraction(125, 117)
    assert exponents_at(c, xs) == [0, 1, Fraction(3, 2), 2, Fraction(5, 2), 3, 4, 5]
    lead = c[-1]
    others = [r for r in roots(lead) if not isinstance(r, tuple) and r != xs]
    assert others
    for x0 in set(others):
        ex = exponents_at(c, x0)
        assert all(not isinstance(e, tuple) and e.denominator == 1 for e in ex), (x0, ex)


# ---- B2: the order-41 recurrence -----------------------------------------------------------------------------------
def _rec():
    return [nmod_poly(cf, P) for cf in DATA["b2"]["polys"]]


def b2_bad_rows(polys):
    g = [nmod(v, P) for v in DATA["b2"]["gen"]]
    m = len(polys) - 1
    return sum(1 for L in range(len(g) - m)
               if int(sum((polys[k](nmod(L, P)) * g[L + k] for k in range(m + 1)), nmod(0, P))))


def test_b2_recurrence_annihilates_the_two_electron_sequence():
    polys = _rec()
    assert len(polys) - 1 == 41 and max(q.degree() for q in polys) == 926
    assert len(DATA["b2"]["gen"]) - 41 >= 350
    assert b2_bad_rows(polys) == 0
    polys[7] = polys[7] + nmod_poly([1], P)
    assert b2_bad_rows(polys) > 350


def b2_asymptotics(polys):
    m, d = len(polys) - 1, max(q.degree() for q in polys)
    chi = nmod_poly([int(q.coeffs()[d]) if q.degree() == d else 0 for q in polys], P)
    out = {}
    for f, mult in chi.factor()[1]:
        if f.degree() != 1:
            continue
        lam = -f.coeffs()[0] / f.coeffs()[1]
        q = [sum((polys[k] * (nmod(math.comb(k, j), P) * lam ** k) for k in range(j, m + 1)), nmod_poly([], P))
             for j in range(m + 1)]
        top = max(qj.degree() - j for j, qj in enumerate(q) if qj.degree() >= 0)
        ind = nmod_poly([], P)
        for j, qj in enumerate(q):
            if qj.degree() >= 0 and qj.degree() - j == top:
                fall = nmod_poly([1], P)
                for t in range(j):
                    fall *= nmod_poly([-t, 1], P)
                ind += fall * nmod_poly([int(qj.coeffs()[qj.degree()])], P)
        out[ratrec(lam)] = roots(ind)
    return out


def test_b2_integer_beta_only_at_x_star():
    asym = b2_asymptotics(_rec())
    xs = Fraction(457, 425)
    integer_at = {lam for lam, bs in asym.items() if lam != 0
                  and any(not isinstance(b, tuple) and b.denominator == 1 for b in bs)}
    assert integer_at == {Fraction(17, 25), Fraction(25, 17)}
    assert all((lam + 1 / lam) / 2 == xs for lam in integer_at)
    assert sorted(asym[Fraction(17, 25)]) == [-6, -5, -4, 1, 3]
    for lam, bs in asym.items():
        if lam != 0 and lam not in integer_at:
            assert bs and all(not isinstance(b, tuple) and b.denominator == 2 for b in bs), (lam, bs)


# ---- B2X: the true coefficients ------------------------------------------------------------------------------------
def _const(key):
    e0, facs = key
    v = mp.exp(mp.mpf(Fraction(e0).numerator) / Fraction(e0).denominator)
    for f in facs:
        a = Fraction(f[1]) if len(f) > 1 else None
        if f[0] == "E1":
            v *= mp.e1(mp.mpf(a.numerator) / a.denominator)
        elif f[0] == "ln":
            v *= mp.log(mp.mpf(a.numerator) / a.denominator)
        elif f[0] == "g":
            v *= mp.euler
        else:
            raise AssertionError(f"unexpected constant {f}")
    return v


def exact_F(dps=700):
    w1 = DATA["w1"]
    mp.mp.dps = dps
    cv = [_const(k) for k in w1["keys"]]
    return {int(L): sum(mp.mpf(int(n)) / int(d) * cv[int(i)] for i, (n, d) in cs.items())
            for L, cs in w1["coeffs"].items()}


def test_b2x_exact_rationals_reproduce_the_annihilated_sequence():
    w1 = DATA["w1"]
    for L, cs in w1["coeffs"].items():
        s = sum((nmod(w1["weights"][int(i)], P) * nmod(int(n), P) / nmod(int(d), P) for i, (n, d) in cs.items()),
                nmod(0, P))
        assert int(s) == DATA["b2"]["gen"][int(L)], L


def test_b2x_coefficients_decay_with_power_minus_four():
    F = exact_F()
    lam = mp.mpf(17) / 25
    beta = {}
    for L in (50, 100):
        assert F[L - 1] > 0 and F[L] > 0 and F[L + 1] > 0
        beta[L] = float(mp.log(F[L + 1] / F[L] / lam) / mp.log(1 + mp.mpf(1) / L))
    assert abs(beta[100] + 4) < 0.5, beta
    assert beta[50] < beta[100] < -4, beta          # approaches -4 from below
    lo = exact_F(900)
    assert abs(lo[100] / F[100] - 1) < mp.mpf(10) ** -20      # precision is not what sets beta


def _tanh_sinh(n, h):
    j = np.arange(-n, n + 1) * h
    x = np.tanh(0.5 * np.pi * np.sinh(j))
    w = 0.5 * np.pi * np.cosh(j) / np.cosh(0.5 * np.pi * np.sinh(j)) ** 2 * h
    keep = (0.5 * (x + 1) > 0) & (0.5 * (x + 1) < 1)
    return 0.5 * (x + 1)[keep], 0.5 * w[keep]


def quad_F(D1, D2, zX, zH, lmax):
    """Independent route: the channel formula of tests/test_paper59_water_bond_angle.py (momentum space)."""
    def G(k, t, D):
        Dl = np.sqrt(k ** 2 * t * (1 - t) + t * zX ** 2 + (1 - t) * zH ** 2)
        return t * (1 - t) * np.exp(-D * Dl) * (D ** 2 * Dl ** 2 + 3 * D * Dl + 3) / Dl ** 5
    t, wt = _tanh_sinh(240, 0.0175)
    u, wu = leggauss(1600)
    u, wu = 0.5 * (u + 1), 0.5 * wu
    k = 3.0 * u / (1 - u)
    wk = wu * 3.0 / (1 - u) ** 2
    G1, G2 = G(k[:, None], t[None, :], D1), G(k[:, None], t[None, :], D2)
    return np.array([(2 * l + 1) * np.sum(wk * ((spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt)
                                          * ((spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt)) for l in range(lmax + 1)])


def test_b2x_exact_split_matches_independent_quadrature():
    F = exact_F(60)
    q = quad_F(1.7, 2.5, 1.4, 1.1, 5)
    for L in (1, 2, 3, 4):
        assert abs(float(F[L + 1] / F[L]) / (q[L + 1] / q[L]) - 1) < 1e-8, L
