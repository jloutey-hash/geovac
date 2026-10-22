"""GT: local data (exponents, log structure / local monodromy type) at finite points and infinity,
and the exponent-sum (Fuchs) bounds for order-1 right factors and for 2-dim invariant subspaces."""
from __future__ import annotations
import itertools
from fractions import Fraction
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

# true singular points of L9 in x = t^2 (t = 1/2, 2, 5/2, 3) and t = 0, infinity
XPTS = [Fraction(1, 4), Fraction(4), Fraction(25, 4), Fraction(9)]


def at_infinity(L):
    """operator in w = 1/x with the same 'theta_t-like' derivation 2 w d/dw (theta_t = -that)."""
    D = G.op_degree(L)
    out = []
    for k, c in enumerate(L):
        cs = [int(v) for v in c.coeffs()] + [0] * (D + 1 - len(c.coeffs()))
        rev = flint.fmpz_poly(cs[::-1])
        out.append(rev * ((-1) ** k))
    return out


def exps_finite_t(L, x0):
    """t-local exponents at t = sqrt(x0) (x0 > 0): same as the h-exponents, h = x - x0"""
    return G.exps_of(G.indicial_at(L, x0))


def exps_t0(L):
    return G.exps_of(G.indicial_t0(L))


def exps_inf_growth(L):
    """mu with y ~ t^mu at infinity"""
    return G.exps_of(G.indicial_inf(L))


def log_structure(L, x0, p=H.P0, N=80):
    """for each distinct exponent (in the local variable h = x - x0; x0 = 'inf' uses w = 1/x),
    whether its Frobenius series exists (no log obstruction).  Returns list of (exponent, multiplicity, ok)."""
    if x0 == "inf":
        Lw = at_infinity(L)
        pt = Fraction(0)
        exps, rest = G.exps_of(G.indicial_at(Lw, 0))
        Luse = Lw
    else:
        pt = Fraction(x0)
        exps, rest = G.exps_of(G.indicial_at(L, pt))
        Luse = L
    out = []
    for e in sorted(set(exps)):
        mult = exps.count(e)
        try:
            G.frob_mod(G.op_mod(Luse, p), pt, e, N, p)
            ok = True
        except ZeroDivisionError:
            ok = False
        out.append((str(e), mult, ok))
    return out, rest


def monodromy_type(L, x0):
    """trivial (no logs, integer exponents), pseudo-reflection, unipotent block sizes etc. (heuristic text):
    the number of log-free Frobenius solutions = number of Jordan blocks within each integer class."""
    info, rest = log_structure(L, x0)
    classes = {}
    for e, m, ok in info:
        q = Fraction(e)
        key = q - (q.numerator // q.denominator)          # class mod Z
        classes.setdefault(key, []).append((q, m, ok))
    desc = []
    for key, lst in classes.items():
        dim = sum(m for _, m, _ in lst)
        # a Jordan block count lower bound: log-free solutions exist exactly for the 'ok' exponents (one each)
        nfree = sum(1 for _, _, ok in lst if ok)
        desc.append({"class": str(key), "dim": dim, "log_free_solutions": nfree,
                     "exponents": [(str(q), m, ok) for q, m, ok in lst]})
    return desc


# ------------------------------------------------------------------ Fuchs-type bounds
def finite_points_t():
    """true singular points of L9 in t other than 0: +-1/2, +-2, +-5/2, +-3 (as x0 values twice)"""
    return [x for x in XPTS for _ in (0, 1)]


def hyperexp_candidates(L):
    """order-1 right factors: s = t^e0 prod_a (t - a)^{e_a} P(t), P polynomial, P(a) != 0;
    deg P = mu - e0 - sum e_a must be a non-negative integer.  Returns list of admissible combos."""
    e0s = sorted(set(exps_t0(L)[0]))
    mus = sorted(set(exps_inf_growth(L)[0]))
    per = [sorted(set(exps_finite_t(L, x0)[0])) for x0 in finite_points_t()]
    hits = []
    for e0 in e0s:
        for mu in mus:
            for combo in itertools.product(*per):
                N = mu - e0 - sum(combo)
                if N >= 0 and Fraction(N).denominator == 1:
                    hits.append((e0, mu, combo, int(N)))
    return hits


def wronskian_candidates(L):
    """2-dim invariant subspaces: theta-Wronskian W = t^{f0} prod (t-a)^{f_a} P, f0 = e1+e2 (t = 0),
    f_a = e1 + e2 - 1 (a != 0), growth at infinity mu1 + mu2; deg P = growth - f0 - sum f_a >= 0."""
    def pairsums(ex, shift):
        return sorted(set(ex[i] + ex[j] - shift for i in range(len(ex)) for j in range(i + 1, len(ex))))
    e0 = exps_t0(L)[0]
    mu = exps_inf_growth(L)[0]
    f0s = pairsums(e0, 0)
    gms = pairsums(mu, 0)
    per = [pairsums(exps_finite_t(L, x0)[0], 1) for x0 in finite_points_t()]
    hits = []
    for f0 in f0s:
        for gm in gms:
            for combo in itertools.product(*per):
                N = gm - f0 - sum(combo)
                if N >= 0 and Fraction(N).denominator == 1:
                    hits.append((f0, gm, combo, int(N)))
    return hits


def hyperexp_search(L, e0, combo, N, p=H.P0):
    """mod p: is there P of degree <= N (in t) with L(t^e0 prod (t-a)^e_a P) = 0?  The ansatz is
    tested on its series at a regular point; nonzero count = number of independent solutions."""
    # work at t-level via series in x around a regular point is awkward for odd t-structure; use the
    # exact t-expansion at t = 0 instead: s = t^e0 * prod (1 - t/a)^{e_a} * P(t).
    t = sp.Symbol("t")
    pts = [sp.Rational(1, 2), sp.Rational(-1, 2), 2, -2, sp.Rational(5, 2), sp.Rational(-5, 2), 3, -3]
    # finite_points_t ordering: x0 = 1/4,1/4,4,4,25/4,25/4,9,9 -> pair with +a, -a
    NT = 4 * (N + 1) + 40 + 2 * G.op_degree(L)
    pre = [1] + [0] * (NT - 1)
    for a, ea in zip(pts, combo):
        # (1 - t/a)^{ea} mod p as a t-series
        ser = [1]
        ea_ = Fraction(ea)
        ia = G.frac_mod(Fraction(1) / Fraction(int(sp.Rational(a).p), int(sp.Rational(a).q)), p)
        eam = G.frac_mod(ea_, p)
        for n in range(1, NT):
            # coefficient: binom(ea, n) (-1/a)^n
            ser.append(ser[-1] * ((eam - (n - 1)) % p) % p * G.inv(n, p) % p * (-ia) % p)
        pre = [int(v) for v in (flint.nmod_poly(pre, p) * flint.nmod_poly(ser, p)).coeffs()][:NT]
        pre += [0] * (NT - len(pre))
    # columns: t^e0 * pre * t^j, j = 0..N ; apply L (theta_t form, coefficients in x = t^2) to t-series
    e0m = G.frac_mod(Fraction(e0), p)
    Lp = G.op_mod(L, p)
    cols = []
    for j in range(N + 1):
        s = [0] * j + pre[: NT - j]
        # L applied: sum_k c_k(t^2) theta^k ; theta^k t^(e0+n) = (e0+n)^k
        res = [0] * NT
        for k, c in enumerate(Lp):
            th = [s[n] * pow((e0m + n) % p, k, p) % p for n in range(NT)]
            cc = [int(v) for v in c.coeffs()]
            # c(t^2): insert zeros
            ct = [0] * (2 * len(cc))
            for i_, v in enumerate(cc):
                ct[2 * i_] = v
            prod = flint.nmod_poly(ct, p) * flint.nmod_poly(th, p)
            pc = [int(v) for v in prod.coeffs()][:NT]
            for n, v in enumerate(pc):
                res[n] = (res[n] + v) % p
        cols.append(res)
    valid = NT - 2 * G.op_degree(L) - 4
    rows = [[col[n] for col in cols] for n in range(valid)]
    return len(G.nullspace_mod(rows, N + 1, p))
