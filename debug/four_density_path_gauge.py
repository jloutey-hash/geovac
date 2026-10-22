"""FP gauge tool: an exact map A = sum_{j<r} a_j(t) theta^j sending Sol(X) into Sol(M) (M o A = 0 mod X).
Mod p search: a_j = sum_{d in [lo, hi], step 2} c_jd t^d / den(t), den = (t^2-1)^a (t^2-4)^b (t^2-9)^c;
condition M(A f) = 0 on a Frobenius solution f of X (series).  Multi-prime CRT + rational reconstruction,
then exact checks in the delta = theta/2, w = t^2 toolkit (four_density_star_ops): M o A = 0 mod X
(right division remainder zero) and the annihilator of A(Sol X) (Gaussian elimination over Q(w))."""
from __future__ import annotations

import time
from fractions import Fraction
from typing import Dict, List

import sympy as sp

import four_density_path_tools as T
import four_density_star_ops as O

u, t = T.u, sp.Symbol("t")
FACT = {"a": t**2 - 1, "b": t**2 - 4, "c": t**2 - 9}


def den_of(pw):
    d = sp.Integer(1)
    for k, e in pw.items():
        d *= FACT[k] ** e
    return sp.expand(d)


def series_inv(den, NT, p):
    cs = [T.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)), p) for c in sp.Poly(den, t).all_coeffs()[::-1]]
    out = [pow(cs[0], p - 2, p)]
    for n in range(1, NT):
        s_ = sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % p
        out.append((-s_) * out[0] % p)
    return {n: c for n, c in enumerate(out) if c}


def system(MQ, XQ, lam, r, pw, lo, hi, NT, p, free_values=None):
    qM = T.qcoeffs(MQ, p)
    qX = T.qcoeffs(XQ, p)
    f, free = T.frob(qX, lam, NT, p, free_values)
    assert f is not None, "X has no log-free solution at lam"
    ftop = lam + 2 * NT
    w = series_inv(den_of(pw), 2 * NT + 4, p)
    top = ftop + lo
    cols, labels = [], []
    for j in range(r):
        g = T.mul(w, T.theta(f, j, p), ftop, p)
        for d in range(lo, hi + 1, 2):
            sh = {e + d: c for e, c in g.items() if e + d < top}
            cols.append(T.apply_Q(qM, sh, top, p))
            labels.append((j, d))
    exps = sorted(set().union(*cols))
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return rows, labels


def nullity(MQ, XQ, lam, r, pw, lo, hi, NT=600, p=T.P61):
    rows, labels = system(MQ, XQ, lam, r, pw, lo, hi, NT, p)
    return len(labels) - T.rank(rows, len(labels), p), len(labels), len(rows)


def primes(n):
    out, q = [], (1 << 61) - 1
    while len(out) < n:
        if sp.isprime(q):
            out.append(q)
        q = sp.prevprime(q)
    return out


def ratrec(a, m):
    a %= m
    if a == 0:
        return Fraction(0)
    r0, r1, s0, s1 = m, a, 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if abs(s1) > bound or s1 == 0:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def reconstruct(MQ, XQ, lam, r, pw, lo, hi, NT=600, maxp=40, verbose=True):
    ps = primes(maxp)
    acc, mod, prev = None, 1, None
    t0 = time.time()
    for n_, p in enumerate(ps):
        rows, labels = system(MQ, XQ, lam, r, pw, lo, hi, NT, p)
        ns = T.nullspace(rows, len(labels), p)
        assert len(ns) == 1, f"nullity {len(ns)}"
        v = ns[0]
        piv = max(i for i, x in enumerate(v) if x)
        iv = pow(v[piv], p - 2, p)
        v = [x * iv % p for x in v]
        if acc is None:
            acc, mod = v, p
        else:
            acc = [(a + mod * ((b - a) * pow(mod, -1, p) % p)) % (mod * p) for a, b in zip(acc, v)]
            mod *= p
        rr = [ratrec(a, mod) for a in acc]
        ok = all(x is not None for x in rr)
        if ok and rr == prev:
            if verbose:
                print(f"   gauge stable after {n_+1} primes ({time.time()-t0:.0f}s)", flush=True)
            break
        prev = rr if ok else None
    else:
        raise RuntimeError("no stable reconstruction")
    a = [sp.Integer(0)] * r
    for c_, (j, d) in zip(rr, labels):
        a[j] += sp.Rational(c_.numerator, c_.denominator) * t**d
    den = den_of(pw)
    return [sp.cancel(x / den) for x in a]


def shiftQ(Q, s_):
    return [sp.expand(sp.sympify(q).subs(sp.Symbol("u"), u).subs(u, u + s_)) for q in Q]


def to_delta(a_even):
    """sum_j a_j(t) theta^j (a_j even in t) -> function-left delta form (theta^j = 2^j delta^j)."""
    return [O.from_theta_t([x])[0] * 2**j for j, x in enumerate(a_even)]


def exact_checks(MQ, XQ, a):
    """returns (M o A = 0 mod X, annihilator of A(Sol X) in delta form or None if A not injective).
    Odd A handled as A = t A_e with M -> t^-1 M t (theta -> theta + 1)."""
    odd = any(sp.cancel(x.subs(t, -t) + x) == 0 and x != 0 for x in a)
    if odd:
        ae = [sp.cancel(x / t) for x in a]
        Mop = O.from_Q(shiftQ(MQ, 1))
    else:
        ae = a
        Mop = O.from_Q(shiftQ(MQ, 0))
    Ae = to_delta(ae)
    Xop = O.from_Q(shiftQ(XQ, 0))
    comp = O.compose(Mop, Ae)
    _, rem = O.rdiv(comp, Xop)
    ok = all(c == 0 for c in rem)
    r = len(Xop) - 1
    S = [-c / Xop[r] for c in Xop[:r]]

    def dvec(vv):
        out = [O.d(vv[0]) + vv[-1] * S[0]]
        for k in range(1, r):
            out.append(vv[k - 1] + O.d(vv[k]) + vv[-1] * S[k])
        return out
    v0 = list(Ae) + [O.K(0)] * (r - len(Ae))
    vecs = [v0]
    for _ in range(r):
        vecs.append(dvec(vecs[-1]))
    Mx = [[vecs[k][i] for k in range(r + 1)] for i in range(r)]
    piv, rk = [], 0
    for col in range(r + 1):
        pr = next((i for i in range(rk, r) if Mx[i][col] != 0), None)
        if pr is None:
            continue
        Mx[rk], Mx[pr] = Mx[pr], Mx[rk]
        iv = 1 / Mx[rk][col]
        Mx[rk] = [x * iv for x in Mx[rk]]
        for i in range(r):
            if i != rk and Mx[i][col] != 0:
                fct = Mx[i][col]
                Mx[i] = [x - fct * y for x, y in zip(Mx[i], Mx[rk])]
        piv.append(col)
        rk += 1
    free = [c for c in range(r + 1) if c not in piv]
    if free != [r]:
        return ok, None, odd
    row = [O.K(0)] * (r + 1)
    row[r] = O.K(1)
    for i, pc in enumerate(piv):
        row[pc] = -Mx[i][r]
    if odd:   # back to the unshifted frame: R = t R_sh t^-1 = R_sh with delta -> delta - 1/2
        out = [O.K(0)] * len(row)
        for k, c in enumerate(row):
            for j in range(k + 1):
                out[j] += c * sp.binomial(k, j) * sp.Rational(-1, 2) ** (k - j)
        row = out
    return ok, row, odd
