"""Mod-p series toolkit for FP (Q-form operators sum_i t^(2i) Q_i(theta_t)); flint nmod_mat for ranks.
Series: dict exponent -> coefficient mod p, complete below a stated top exponent."""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Tuple

import flint
import sympy as sp

u = sp.Symbol("u")
P61 = (1 << 61) - 1


def fr(x, p=P61):
    x = Fraction(x)
    return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def qcoeffs(Q, p=P61):
    out = []
    for q in Q:
        pq = sp.Poly(sp.expand(sp.sympify(q).subs(sp.Symbol("u"), u)), u)
        cs = [0] * (pq.degree() + 1) if not pq.is_zero else [0]
        for (m,), c in pq.terms():
            c = sp.Rational(c)
            cs[m] = fr(Fraction(int(c.p), int(c.q)), p)
        out.append(cs)
    return out


def ev(cs, x, p=P61):
    acc = 0
    xm = x % p
    for c in reversed(cs):
        acc = (acc * xm + c) % p
    return acc


def frob(qc, lam, NT, p=P61, free_values=None):
    """t^lam sum_n c_n t^(2n), n < NT; free slots (resonances) set to free_values.get(exp, 0).
    Returns (series, free_exps) or (None, exp) if a resonance is inconsistent."""
    c = [1]
    free = []
    for n in range(1, NT):
        s_ = 0
        for i in range(1, len(qc)):
            if n - i >= 0 and c[n - i]:
                s_ = (s_ + ev(qc[i], lam + 2 * (n - i), p) * c[n - i]) % p
        d = ev(qc[0], lam + 2 * n, p)
        if d == 0:
            if s_ != 0:
                return None, lam + 2 * n
            free.append(lam + 2 * n)
            c.append((free_values or {}).get(lam + 2 * n, 0) % p)
        else:
            c.append((-s_) * pow(d, p - 2, p) % p)
    return {lam + 2 * n: cn for n, cn in enumerate(c) if cn}, free


def apply_Q(qc, ser, top, p=P61):
    out = {}
    for i, cs in enumerate(qc):
        for e, c in ser.items():
            if e + 2 * i >= top:
                continue
            v = ev(cs, e, p)
            if v:
                out[e + 2 * i] = (out.get(e + 2 * i, 0) + v * c) % p
    return {e: c for e, c in out.items() if c}


def theta(ser, j, p=P61):
    return {e: c * pow(e % p, j, p) % p for e, c in ser.items()}


def mul(a, b, top, p=P61):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = (out.get(e1 + e2, 0) + c1 * c2) % p
    return out


def rank(rows, ncols, p=P61):
    if not rows:
        return 0
    M = flint.nmod_mat(len(rows), ncols, [x for r in rows for x in r], p)
    return M.rank()


def nullspace(rows, ncols, p=P61):
    M = flint.nmod_mat(len(rows), ncols, [x for r in rows for x in r], p)
    X, nul = M.nullspace()
    return [[int(X[i, k]) for i in range(ncols)] for k in range(nul)]


def guess_system(series_list, k, I, p=P61):
    """rows for: operators sum_{j<=k, i<=I} a_ij t^(2i) theta^j killing every (series, top) in the list."""
    labels = [(i, j) for i in range(I + 1) for j in range(k + 1)]
    rows = []
    for ser, top in series_list:
        exps = sorted(e for e in set(e + 2 * i for e in ser for i in range(I + 1)) if e < top)
        lam = min(ser)
        allexps = list(range(lam, top, 2))
        idx = {e: n for n, e in enumerate(allexps)}
        R = [[0] * len(labels) for _ in allexps]
        for e, c in ser.items():
            for col, (i, j) in enumerate(labels):
                ee = e + 2 * i
                if ee < top:
                    R[idx[ee]][col] = (R[idx[ee]][col] + c * pow(e % p, j, p)) % p
        rows += R
    return rows, labels


def min_annihilator_order(series_list, kmax, Imax, p=P61):
    """smallest k with nullity > 0 at t-degree 2*Imax; then the smallest I for it."""
    for k in range(1, kmax + 1):
        rows, labels = guess_system(series_list, k, Imax, p)
        if len(rows) < len(labels) + 10:
            return None, "too few terms"
        nl = len(labels) - rank(rows, len(labels), p)
        if nl:
            lo, hi = 0, Imax
            while lo < hi:
                mid = (lo + hi) // 2
                rr, ll = guess_system(series_list, k, mid, p)
                if len(ll) - rank(rr, len(ll), p):
                    hi = mid
                else:
                    lo = mid + 1
            return k, lo
    return None, None
