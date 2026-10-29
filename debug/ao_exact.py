"""Exact decomposition of the bond-angle Legendre coefficients of the two-electron three-centre integral (XY|XZ)
(record AO3; numeric routes debug/ao_radial.py, debug/ao_radial2.py).  With rational D1, D2, zX, zH:

    F_L / (4 zX^3 zH^3) = sum_c q_{L,c} c,   c = e^{e0} * (product of at most two of E1(z), ln v, gamma, K(...)),

    F_L = 4 pi^2 N^2 int int g1(r1) g2(r2) r<^L / r>^(L+1),   g_i(r) = r^2 e^{-zX r} G_L^{(D_i)}(r)
(G_L as in debug/na_exact.py).  Inner integrals C(r1) = int_0^r1 g2 r^L, T(r1) = int_r1^inf g2 r^(-L-1) by the
antiderivatives Phi_n(b; r) = e^{br} P_n(r) + kappa_n E1(-b r) (exact recursions in n, summed Horner-style); outer
integrals of r^n e^{br} by J_n (debug/na_exact.Integrals) and of r^n e^{br} E1(beta r) by
  (n+1) K_n + b K_{n+1} = [r^{n+1} e^{br} E1(beta r)] + J_n(b - beta),
with the base K_{-1} = int e^{br} E1(beta r) / r dr a named constant per (b, beta, piece).  Pieces split at D1, D2;
the piece at r = 0 by Hadamard finite part throughout (the total is integrable, so the finite parts sum to it).
Field pluggable as in debug/na_exact.py.
"""
from fractions import Fraction

import flint

from na_exact import Field, Integrals, bessel_tables


def kconv(key):
    """na_exact key (e0, kind, arg) -> (e0, factors)."""
    e0, kind, arg = key
    if kind == '1':
        return (e0, ())
    if kind == 'g':
        return (e0, (('g',),))
    return (e0, ((kind, arg),))


def kmul(a, b):
    return (a[0] + b[0], tuple(sorted(a[1] + b[1], key=repr)))


def add(d, key, v):
    if key in d:
        d[key] = d[key] + v
    else:
        d[key] = v


def mono_G(L, inner, D, zH, it, kt, K):
    """G_L^{(D)}(r) = -2[(2L+1) i_L(x) k_L(y) + x i_{L+1}(x) k_L(y) - y i_L(x) k_{L+1}(y)] as {(e0, b): {n: coef}}."""
    out = {}

    def put(e0, b, n, v):
        d = out.setdefault((e0, b), {})
        add(d, n, K(-2) * v)

    if inner:
        y = zH * D
        kL = sum((K.pw(y, -j) * v for j, v in kt[L].items()), K(0))
        kL1 = sum((K.pw(y, -j) * v for j, v in kt[L + 1].items()), K(0))
        for s in (1, -1):
            for j, v in it[L][s].items():
                put(-y, s * zH, -j, K.pw(zH, -j) * v * (K(2 * L + 1) * kL - K(y) * kL1))
            for j, v in it[L + 1][s].items():
                put(-y, s * zH, 1 - j, K.pw(zH, 1 - j) * v * kL)
    else:
        x = zH * D
        for s in (1, -1):
            iL = sum((K.pw(x, -j) * v for j, v in it[L][s].items()), K(0))
            iL1 = sum((K.pw(x, -j) * v for j, v in it[L + 1][s].items()), K(0))
            for j, v in kt[L].items():
                put(s * x, -zH, -j, K.pw(zH, -j) * v * (K(2 * L + 1) * iL + K(x) * iL1))
            for j, v in kt[L + 1].items():
                put(s * x, -zH, 1 - j, -K.pw(zH, 1 - j) * v * iL)
    return out


def g_of(L, inner, D, zX, zH, it, kt, K, shift):
    """r^(2+shift) e^{-zX r} G_L as {(e0, b): {n: coef}}."""
    G = mono_G(L, inner, D, zH, it, kt, K)
    return {(e0, b - zX): {n + 2 + shift: v for n, v in poly.items()} for (e0, b), poly in G.items()}


class Mul:
    """Laurent-polynomial products through flint polynomials."""

    def __init__(self, K):
        self.K = K
        self.p = K.p

    def mul(self, a, b):
        if not a or not b:
            return {}
        oa, ob = min(a), min(b)
        if self.p is None:
            A = flint.fmpq_poly([a.get(i + oa, 0) for i in range(max(a) - oa + 1)])
            B = flint.fmpq_poly([b.get(i + ob, 0) for i in range(max(b) - ob + 1)])
        else:
            A = flint.nmod_poly([int(a.get(i + oa, 0)) for i in range(max(a) - oa + 1)], self.p)
            B = flint.nmod_poly([int(b.get(i + ob, 0)) for i in range(max(b) - ob + 1)], self.p)
        C = A * B
        cs = C.coeffs()
        o = oa + ob
        K = self.K
        if self.p is None:
            return {i + o: c for i, c in enumerate(cs) if c != 0}
        return {i + o: flint.nmod(int(c), self.p) for i, c in enumerate(cs) if int(c) != 0}


def phi_sum(c, b, K):
    """sum_n c_n Phi_n(b; r) = e^{br} P(r) + kappa E1(-b r);  returns (P dict, kappa)."""
    P = {}
    kappa = K(0)
    if not c:
        return P, kappa
    ib = K(1) / K(Fraction(b))
    pos = sorted((n for n in c if n >= 0), reverse=True)
    carry = K(0)
    if pos:
        for n in range(pos[0], -1, -1):
            a = c.get(n, K(0)) + carry
            if a != 0:
                add(P, n, a * ib)
            carry = -K(n) * a * ib
    neg = [n for n in c if n <= -2]
    carry = K(0)
    if neg:
        for n in range(min(neg), -1):
            a = c.get(n, K(0)) + carry
            if a != 0:
                add(P, n + 1, a / K(n + 1))
            carry = -K(Fraction(b)) * a / K(n + 1)
    kappa = -(c.get(-1, K(0)) + carry)
    return P, kappa


def phi_at(P, kappa, b, R, K):
    """e^{bR} P(R) + kappa E1(-bR) as {key: coef}; R = 0 finite part; R = None (infinity) -> 0."""
    out = {}
    if R is None:
        return out
    if R == 0:
        s = K(0)
        for k, v in P.items():
            if k <= 0:
                s += v * K.pw(b, -k) * K.ifact(-k)
        if s != 0:
            add(out, (Fraction(0), ()), s)
        if kappa != 0:                       # FP E1(-b eps) = -gamma - ln|b|
            add(out, (Fraction(0), (('g',),)), -kappa)
            add(out, (Fraction(0), (('ln', abs(Fraction(b))),)), -kappa)
        return out
    s = K(0)
    for k, v in P.items():
        s += v * K.pw(R, k)
    if s != 0:
        add(out, (Fraction(b) * R, ()), s)
    if kappa != 0:
        add(out, (Fraction(0), (('E1', -Fraction(b) * R),)), kappa)
    return out


class KInt:
    """K_n(b, beta; A, B) = int_A^B r^n e^{br} E1(beta r) dr as {key: coef}, cached tables."""

    def __init__(self, K, J):
        self.K = K
        self.J = J
        self.cache = {}

    def endpoint(self, b, beta, R, m):
        K = self.K
        if R is None:
            return {}
        if R == 0:
            # eps^0 coefficient of eps^m e^{b eps} E1(beta eps), E1(z) = -gamma - ln|z| - sum_k (-z)^k/(k k!)
            if m > 0:
                return {}
            s = -m
            out = {}
            c0 = K.pw(b, s) * K.ifact(s)                  # b^s / s!
            add(out, (Fraction(0), (('g',),)), -c0)
            add(out, (Fraction(0), (('ln', abs(Fraction(beta))),)), -c0)
            r = K(0)
            for k in range(1, s + 1):
                r -= K.pw(-Fraction(beta), k) / K(k) * K.ifact(k) * K.pw(b, s - k) * K.ifact(s - k)
            if r != 0:
                add(out, (Fraction(0), ()), r)
            return out
        return {(Fraction(b) * R, (('E1', Fraction(beta) * R),)): K.pw(R, m)}

    def table(self, b, beta, A, B, nmin, nmax):
        K = self.K
        key = (b, beta, A, B)
        tab = self.cache.setdefault(key, {})
        if -1 not in tab:
            tab[-1] = {(Fraction(0), (('K', Fraction(b), Fraction(beta), A, B),)): K(1)}
        bb = b - beta
        Jt = self.J.table(bb, A, B, min(nmin, -1), max(nmax, -1))
        n = max(tab)
        while n < nmax:                      # K_{n+1} = ([r^{n+1} e^{br} E1] + J_n(b-beta) - (n+1) K_n) / b
            d = {}
            ib = K(1) / K(Fraction(b))
            for kk, v in self.endpoint(b, beta, B, n + 1).items():
                add(d, kk, v * ib)
            for kk, v in self.endpoint(b, beta, A, n + 1).items():
                add(d, kk, -v * ib)
            for kk, v in Jt[n].items():
                add(d, kconv(kk), v * ib)
            for kk, v in tab[n].items():
                add(d, kk, -K(n + 1) * v * ib)
            tab[n + 1] = d
            n += 1
        n = min(tab)
        while n > nmin:                      # K_{n-1} = ([r^n e^{br} E1] + J_{n-1} - b K_n) / n
            m = n - 1
            d = {}
            inv = K(1) / K(m + 1)
            for kk, v in self.endpoint(b, beta, B, m + 1).items():
                add(d, kk, v * inv)
            for kk, v in self.endpoint(b, beta, A, m + 1).items():
                add(d, kk, -v * inv)
            for kk, v in Jt[m].items():
                add(d, kconv(kk), v * inv)
            for kk, v in tab[n].items():
                add(d, kk, -K(Fraction(b)) * v * inv)
            tab[m] = d
            n = m
        return tab


def dotJ(J, poly, b, A, B, res, pref, K):
    """res += pref * sum_n poly[n] J_n(b; A, B)  (pref: (key, coef) of a constant multiplying the whole term)."""
    if not poly:
        return
    tab = J.table(b, A, B, min(min(poly), -1), max(max(poly), -1))
    pk, pv = pref
    acc = {}
    for n, c in poly.items():
        for kk, v in tab[n].items():
            add(acc, kk, c * v)
    for kk, v in acc.items():
        add(res, kmul(pk, kconv(kk)), pv * v)


def dotK(Kt, poly, b, beta, A, B, res, pref, K):
    if not poly:
        return
    tab = Kt.table(b, beta, A, B, min(min(poly), -1), max(max(poly), -1))
    pk, pv = pref
    acc = {}
    for n, c in poly.items():
        for kk, v in tab[n].items():
            add(acc, kk, c * v)
    for kk, v in acc.items():
        add(res, kmul(pk, kk), pv * v)


def coefficients(D1, D2, zX, zH, Lmax, K, Lmin=0):
    D1, D2, zX, zH = (Fraction(v) for v in (D1, D2, zX, zH))
    it, kt = bessel_tables(K, Lmax + 2)
    J = Integrals(K)
    Kt = KInt(K, J)
    M = Mul(K)
    lo, hi = sorted([D1, D2])
    pieces = [(Fraction(0), lo), (lo, hi), (hi, None)] if lo != hi else [(Fraction(0), lo), (lo, None)]
    one = (Fraction(0), ())
    out = []
    for L in range(Lmin, Lmax + 1):
        res = {}
        info = []
        for (A, B) in pieces:
            mid = A + 1 if B is None else (A + B) / 2
            g1 = g_of(L, mid < D1, D1, zX, zH, it, kt, K, 0)
            g2C = g_of(L, mid < D2, D2, zX, zH, it, kt, K, L)          # g2 r^L
            g2T = g_of(L, mid < D2, D2, zX, zH, it, kt, K, -L - 1)     # g2 r^(-L-1)
            info.append((A, B, g1, g2C, g2T))
        # constants: full-piece integrals of g2 r^L (for C) and g2 r^(-L-1) (for T)
        fullC, fullT = [], []
        for (A, B, g1, g2C, g2T) in info:
            dC, dT = {}, {}
            for (e0, b), poly in g2C.items():
                dotJ(J, poly, b, A, B, dC, ((e0, ()), K(1)), K)
            for (e0, b), poly in g2T.items():
                dotJ(J, poly, b, A, B, dT, ((e0, ()), K(1)), K)
            fullC.append(dC)
            fullT.append(dT)
        for k, (A, B, g1, g2C, g2T) in enumerate(info):
            constC, constT = {}, {}
            for j in range(k):
                for kk, v in fullC[j].items():
                    add(constC, kk, v)
            for j in range(k + 1, len(info)):
                for kk, v in fullT[j].items():
                    add(constT, kk, v)
            # variable parts: C = sum e^{e0}[Phi(r) - Phi(A)], T = sum e^{e0}[Phi(B) - Phi(r)]
            varC, varT = [], []
            for (e0, b), poly in g2C.items():
                P, kap = phi_sum(poly, b, K)
                varC.append((e0, b, P, kap))
                for kk, v in phi_at(P, kap, b, A, K).items():
                    add(constC, kmul((e0, ()), kk), -v)
            for (e0, b), poly in g2T.items():
                P, kap = phi_sum(poly, b, K)
                varT.append((e0, b, P, kap))
                for kk, v in phi_at(P, kap, b, B, K).items():
                    add(constT, kmul((e0, ()), kk), v)
            for (e1, b1), p1 in g1.items():
                pC = {n - L - 1: v for n, v in p1.items()}            # g1 r^(-L-1)
                pT = {n + L: v for n, v in p1.items()}                # g1 r^L
                for kk, v in constC.items():
                    dotJ(J, pC, b1, A, B, res, (kmul((e1, ()), kk), v), K)
                for kk, v in constT.items():
                    dotJ(J, pT, b1, A, B, res, (kmul((e1, ()), kk), v), K)
                for (e0, b, P, kap), sgn, pp in [(x, 1, pC) for x in varC] + [(x, -1, pT) for x in varT]:
                    prod = M.mul(pp, P)
                    dotJ(J, prod, b1 + b, A, B, res, ((e1 + e0, ()), K(sgn)), K)
                    if kap != 0:
                        dotK(Kt, pp, b1, -b, A, B, res, ((e1 + e0, ()), K(sgn) * kap), K)
        out.append({k: v for k, v in res.items() if v != 0})
    return out


def value(key, mp):
    e0, facs = key
    v = mp.exp(mp.mpf(e0.numerator) / e0.denominator)
    for f in facs:
        v *= factor_value(f, mp)
    return v


_KV = {}


def factor_value(f, mp):
    q = lambda x: mp.mpf(x.numerator) / x.denominator
    if f[0] == 'g':
        return mp.euler
    if f[0] == 'ln':
        return mp.log(q(f[1]))
    if f[0] == 'E1':
        a = q(f[1])
        return mp.e1(a) if a > 0 else -mp.ei(-a)
    key = (f, mp.mp.dps)
    if key not in _KV:
        _, b, beta, A, B = f
        b, beta = q(b), q(beta)
        E = lambda r: mp.e1(beta * r) if beta > 0 else -mp.ei(-beta * r)
        hi = mp.inf if B is None else q(B)
        if A == 0:
            c = -mp.euler - mp.log(abs(beta))
            # FP int_0^B e^{br}E1(beta r)/r = int_0^B [e^{br}E1(beta r) - (c - ln r)]/r + FP int_0^B (c - ln r)/r
            if B is None:
                raise ValueError('K with A = 0 and B = inf not used')
            v = mp.quad(lambda r: (mp.exp(b * r) * E(r) - (c - mp.log(r))) / r, [0, hi])
            v += c * mp.log(hi) - mp.log(hi) ** 2 / 2
        else:
            v = mp.quad(lambda r: mp.exp(b * r) * E(r) / r, [q(A), hi])
        _KV[key] = v
    return _KV[key]
