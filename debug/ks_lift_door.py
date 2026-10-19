"""Pre-registered KS1-KS5: hydrogen as a constrained four-dimensional oscillator.

Record: docs/predictions/2026-10-19_ks_lift.md (committed be16fccd before this script ran).

KS1/KS2: exact integer arithmetic on Bargmann monomials (A+, A-, B+, B-).
KS3: sympy push-down of the Q = 0 oscillator states to R^3, checked against the Sturmian equation.
KS4: exact inner products in R^4 and the 1/r-weighted ones in R^3.
KS5: spectral zeta sums at 50 digits against closed forms.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp

OUT = Path(__file__).resolve().parent / "data" / "ks_lift_door.json"
Mono = Tuple[int, int, int, int]  # (a+, a-, b+, b-)
Op = Dict[Mono, Dict[Mono, Fraction]]


def level_basis(n: int) -> List[Mono]:
    """Q = 0 monomials at N = 2(n-1): a+ + b+ = a- + b- = n-1."""
    k = n - 1
    return [(ap, am, k - ap, k - am) for ap in range(k + 1) for am in range(k + 1)]


def q_zero_count(N: int) -> int:
    c = 0
    for ap in range(N + 1):
        for am in range(N + 1 - ap):
            for bp in range(N + 1 - ap - am):
                bm = N - ap - am - bp
                if (ap - am) + (bp - bm) == 0:
                    c += 1
    return c


def apply_term(mono: Mono, raise_i: int, lower_i: int) -> Tuple[Mono, int]:
    """(variable raise_i) x (d / d variable lower_i) on a monomial: returns (image, coefficient)."""
    e = list(mono)
    c = e[lower_i]
    if c == 0:
        return mono, 0
    e[lower_i] -= 1
    e[raise_i] += 1
    return tuple(e), c


def op_matrix(basis: List[Mono], terms: List[Tuple[Fraction, int, int]], diag=None) -> List[List[Fraction]]:
    idx = {m: i for i, m in enumerate(basis)}
    d = len(basis)
    M = [[Fraction(0)] * d for _ in range(d)]
    for j, m in enumerate(basis):
        if diag is not None:
            M[j][j] += diag(m)
        for coef, ri, li in terms:
            img, c = apply_term(m, ri, li)
            if c:
                M[idx[img]][j] += coef * c
    return M


def mm(A, B):
    d = len(A)
    Bt = list(zip(*B))
    return [[sum(A[i][k] * Bt[j][k] for k in range(d) if A[i][k]) for j in range(d)] for i in range(d)]


def madd(A, B, s=1):
    return [[a + s * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def scal(A, s):
    return [[s * a for a in r] for r in A]


def eye(d):
    return [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]


def is_scalar(A, s) -> bool:
    d = len(A)
    return all(A[i][j] == (s if i == j else 0) for i in range(d) for j in range(d))


AP, AM, BP, BM = 0, 1, 2, 3
H = Fraction(1, 2)


def su2_ops(basis):
    """J (on A+,B+) and J' (conjugate doublet on A-,B-), as matrices."""
    Jz = op_matrix(basis, [], diag=lambda m: H * (m[AP] - m[BP]))
    Jp = op_matrix(basis, [(Fraction(1), AP, BP)])          # A+ d/dB+
    Jm = op_matrix(basis, [(Fraction(1), BP, AP)])          # B+ d/dA+
    Kz = op_matrix(basis, [], diag=lambda m: -H * (m[AM] - m[BM]))
    Kp = op_matrix(basis, [(Fraction(-1), BM, AM)])         # -B- d/dA-
    Km = op_matrix(basis, [(Fraction(-1), AM, BM)])         # -A- d/dB-
    return (Jz, Jp, Jm), (Kz, Kp, Km)


def casimir(z, p, m):
    return madd(madd(mm(m, p), mm(z, z)), z)


def stage_ks12(nmax_count=15, nmax_minpoly=8, nmax_casimir=10):
    res = {"KS1_count": {}, "KS1_odd_levels_empty": True, "KS1_l_content": {}, "KS1_minpoly": {},
           "KS2": {}}
    for N in range(1, 2 * (nmax_count - 1) + 1, 2):
        if q_zero_count(N):
            res["KS1_odd_levels_empty"] = False
    for n in range(1, nmax_count + 1):
        cnt = q_zero_count(2 * (n - 1))
        basis = level_basis(n)
        lz = Counter((m[AP] - m[BP] - m[AM] + m[BM]) // 2 for m in basis)  # integer on Q = 0
        lcont = {l: lz[l] - lz.get(l + 1, 0) for l in range(0, n + 2)}
        ok = cnt == n * n and all(lcont[l] == (1 if l < n else 0) for l in lcont)
        res["KS1_count"][n] = cnt
        res["KS1_l_content"][n] = ok
    for n in range(1, max(nmax_minpoly, nmax_casimir) + 1):
        basis = level_basis(n)
        (Jz, Jp, Jm), (Kz, Kp, Km) = su2_ops(basis)
        if n <= nmax_minpoly:
            Lz, Lp, Lm = madd(Jz, Kz), madd(Jp, Kp), madd(Jm, Km)
            L2 = casimir(Lz, Lp, Lm)
            P = eye(len(basis))
            for l in range(n):
                P = mm(P, madd(L2, scal(eye(len(basis)), l * (l + 1)), -1))
            zero = all(v == 0 for r in P for v in r)
            # and no proper sub-product vanishes: drop the top l
            Q = eye(len(basis))
            for l in range(n - 1):
                Q = mm(Q, madd(L2, scal(eye(len(basis)), l * (l + 1)), -1))
            res["KS1_minpoly"][n] = {"product_zero": zero, "shorter_product_nonzero":
                                     (n == 1) or any(v != 0 for r in Q for v in r)}
        if n <= nmax_casimir:
            J2 = casimir(Jz, Jp, Jm)
            K2 = casimir(Kz, Kp, Km)
            C = scal(madd(J2, K2), 2)
            res["KS2"][n] = {"J2_equals_K2": J2 == K2, "casimir_is_n2m1": is_scalar(C, n * n - 1)}
        print("KS1/2 n", n, flush=True)
    return res


# ---------------- KS3 / KS4 (sympy) ----------------
z1, w1, z2, w2 = sp.symbols("z1 w1 z2 w2")  # w = conj(z), treated as independent
om = sp.symbols("omega", positive=True)
X1, X2, X3 = sp.symbols("x1 x2 x3", real=True)
RR = sp.sqrt(X1**2 + X2**2 + X3**2)


def oscillator_poly(mono: Mono) -> sp.Expr:
    """P with (a+_A)^ap (a-_A)^am (a+_B)^bp (a-_B)^bm G = P G, G = exp(-2 omega (z1 w1 + z2 w2))."""
    G = sp.exp(-2 * om * (z1 * w1 + z2 * w2))
    F = G
    ops = [(z1, w1)] * mono[AP] + [(w1, z1)] * mono[AM] + [(z2, w2)] * mono[BP] + [(w2, z2)] * mono[BM]
    for mult, dvar in ops:
        F = 4 * om * mult * F - 2 * sp.diff(F, dvar)
    return sp.expand(sp.simplify(F / G))


RS = sp.symbols("rs", positive=True)  # stands for r = |x| in pushed-down polynomials


def push_down(P: sp.Expr) -> sp.Expr:
    r_ = RS
    b11, b22 = (r_ + X3) / 2, (r_ - X3) / 2
    b12, b21 = (X1 + sp.I * X2) / 2, (X1 - sp.I * X2) / 2   # z1 w2, z2 w1
    out = 0
    for (a, b, c, d), coef in sp.Poly(P, z1, w1, z2, w2).terms():
        assert a + c == b + d, (a, b, c, d)
        k11 = min(a, b)
        a_, b_ = a - k11, b - k11
        term = b11**k11
        if a_:
            term *= b12**a_ * b22**(c)      # z1 with w2 (a_ of them), then z2 with remaining w2 = c
        else:
            term *= b21**b_ * b22**(d)      # w1 with z2 (b_), then z2 left = d with w2
        out += coef * term
    return sp.expand(out)


def lap3(f):
    return sum(sp.diff(f, v, 2) for v in (X1, X2, X3))


def stage_ks3(nmax=3):
    res = {}
    for n in range(1, nmax + 1):
        p0 = 2 * om
        Z = n * p0
        for mono in level_basis(n):
            Q = push_down(oscillator_poly(mono)).subs(RS, RR)
            f = Q * sp.exp(-p0 * RR)
            resid = -lap3(f) / 2 - Z / RR * f + p0**2 / 2 * f
            resid = sp.simplify(sp.expand(resid / sp.exp(-p0 * RR)))
            res[f"n{n}_{mono}"] = bool(resid == 0)
            print("KS3", n, mono, resid == 0, flush=True)
    return res


def gauss_moment(k: int, a) -> sp.Expr:
    if k % 2:
        return sp.Integer(0)
    m = k // 2
    return sp.gamma(sp.Rational(2 * m + 1, 2)) / a**sp.Rational(2 * m + 1, 2)


def inner4(P: sp.Expr, Pp: sp.Expr) -> sp.Expr:
    u1, u2, u3, u4 = sp.symbols("u1 u2 u3 u4", real=True)
    sub = {z1: u1 + sp.I * u2, w1: u1 - sp.I * u2, z2: u3 + sp.I * u4, w2: u3 - sp.I * u4}
    Pc = sp.conjugate(P.subs(sub))
    prod = sp.expand(Pc * Pp.subs(sub))
    a = 4 * om
    tot = 0
    for (e1, e2, e3, e4), coef in sp.Poly(prod, u1, u2, u3, u4).terms():
        tot += coef * gauss_moment(e1, a) * gauss_moment(e2, a) * gauss_moment(e3, a) * gauss_moment(e4, a)
    return sp.simplify(tot)


def inner3_weighted(Q: sp.Expr, Qp: sp.Expr) -> sp.Expr:
    rs = RS
    prod = sp.expand(sp.conjugate(Q) * Qp)
    tot = 0
    for (a, b, c, k), coef in sp.Poly(prod, X1, X2, X3, rs).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        ang = 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * \
            sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(a + b + c + 3, 2))
        p = a + b + c + k + 1          # r^(a+b+c+k) r^2 / r
        rad = sp.factorial(p) / (4 * om)**(p + 1)
        tot += coef * ang * rad
    return sp.simplify(tot)


def stage_ks4():
    states = [(1, m) for m in level_basis(1)] + [(2, m) for m in level_basis(2)]
    three = [(3, (2, 2, 0, 0)), (3, (1, 1, 1, 1)), (3, (2, 1, 0, 1))]
    states += three
    polys = {s: oscillator_poly(s[1]) for s in states}
    downs = {s: push_down(polys[s]) for s in states}
    ratios, zeros_ok, others = [], True, []
    pairs = [(s, t) for s in states for t in states if (s[0] <= 2 and t[0] <= 2) or (s == t)]
    for s, t in pairs:
        I4 = inner4(polys[s], polys[t])
        I3 = inner3_weighted(downs[s], downs[t])
        if I3 == 0 or I4 == 0:
            zeros_ok &= (I3 == 0 and I4 == 0)
        else:
            q = sp.nsimplify(sp.simplify(I4 / I3))
            ratios.append(str(q))
        print("KS4", s, t, flush=True)
    return {"distinct_ratios": sorted(set(ratios)), "n_nonzero": len(ratios), "zeros_consistent": zeros_ok}


def stage_ks5():
    mp.mp.dps = 60
    z = mp.zeta
    out = {}
    cas2 = mp.nsum(lambda n: n**2 / (n**2 - 1)**2, [2, mp.inf])
    out["casimir_s2_minus_closed"] = mp.nstr(cas2 - (mp.pi**2 / 12 + mp.mpf(1) / 16), 5)
    for s in (4, 6):
        cas = mp.nsum(lambda n: n**2 / (n**2 - 1)**s, [2, mp.inf])
        ks = mp.nsum(lambda n: n**2 * n**(-s), [1, mp.inf])
        eul = mp.nsum(lambda N: (N + 1) * (N + 2) / 2 * N**(-s), [1, mp.inf])
        en = mp.nsum(lambda N: (N + 1) * (N + 2) / 2 * (N + mp.mpf(3) / 2)**(-s), [0, mp.inf])
        en_closed = (mp.zeta(s - 2, 1.5) - mp.zeta(s, 1.5) / 4) / 2
        hur = lambda t: (2**t - 1) * z(t) - 2**t
        out[f"s{s}"] = {
            "ks_minus_zeta(s-2)": mp.nstr(ks - z(s - 2), 5),
            "euler_minus_closed": mp.nstr(eul - (z(s - 2) + 3 * z(s - 1) + 2 * z(s)) / 2, 5),
            "energy_minus_hurwitz": mp.nstr(en - en_closed, 5),
            "hurwitz_3/2_identity": mp.nstr(mp.zeta(s, 1.5) - hur(s), 5),
            "casimir_value": mp.nstr(cas, 30),
        }
        # the Casimir ladder at s = 4, 6 in closed form: partial fractions in A = 1/(n-1), B = 1/(n+1)
    return out


def run(stages):
    data = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "ks12" in stages:
        data["KS1_KS2"] = stage_ks12()
    if "ks3" in stages:
        data["KS3"] = stage_ks3()
    if "ks4" in stages:
        data["KS4"] = stage_ks4()
    if "ks5" in stages:
        data["KS5"] = stage_ks5()
    OUT.write_text(json.dumps(data, indent=1, default=str))
    print(json.dumps(data, indent=1, default=str))


if __name__ == "__main__":
    run(sys.argv[1:] or ["ks12", "ks3", "ks4", "ks5"])
