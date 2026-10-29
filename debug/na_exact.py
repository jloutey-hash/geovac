"""Exact decomposition of the Legendre coefficients of the one-electron three-centre nuclear-attraction integral
(record NA; numeric route debug/ao_na.py).  With rational D1, D2, zX, zH:

    F_L / (2 pi N) = sum_c q_{L,c} * c,

c running over a short list of constants  e^{e0} * {1, E1(z), ln(v), gamma}  (e0, z, v rational), q_{L,c} exact.
Built bottom-up: i_L, k_L as e^{+-x} times exact polynomials in 1/x (three-term rules in L); the density multipole
rho_L = (2L+1)/2 N e^{-zX r} G_L,  G_L = -2[(2L+1) i_L(x) k_L(y) + x i_{L+1}(x) k_L(y) - y i_L(x) k_{L+1}(y)],
x = zH min(r, D1), y = zH max(r, D1);  F_L = 2 pi N int r^2 e^{-zX r} G_L(r) w_L(r) dr,  w_L = r<^L / r>^(L+1) of
(r, D2).  Radial integrals int r^n e^{b r} by integration by parts down to int e^{b r}/r = E1 differences; the piece
starting at r = 0 by Hadamard finite part (the 1/r^k terms cancel in the sum), which brings gamma and ln|b|.

The field K is pluggable: fmpq (exact, for checks) or nmod (mod a prime, for guessing).
"""
from fractions import Fraction

import flint


class Field:
    def __init__(self, p=None):
        self.p = p

    def __call__(self, q):
        q = Fraction(q)
        if self.p is None:
            return flint.fmpq(q.numerator, q.denominator)
        return flint.nmod(q.numerator, self.p) / flint.nmod(q.denominator, self.p)

    def zero(self):
        return self(0)


def bessel_tables(K, Lmax):
    """i[L] = {sigma: {j: coef}} for e^{sigma x} x^{-j};  k[L] = {j: coef} for e^{-x} x^{-j}; L = 0..Lmax."""
    h = K(Fraction(1, 2))
    i = [{1: {1: h}, -1: {1: -h}}, {1: {1: h, 2: -h}, -1: {1: h, 2: h}}]
    k = [{1: K(1)}, {1: K(1), 2: K(1)}]
    for L in range(1, Lmax):
        c = K(2 * L + 1)
        new = {}
        for s in (1, -1):
            d = dict(i[L - 1][s])
            for j, v in i[L][s].items():
                d[j + 1] = d.get(j + 1, K(0)) - c * v
            new[s] = d
        i.append(new)
        d = dict(k[L - 1])
        for j, v in k[L].items():
            d[j + 1] = d.get(j + 1, K(0)) + c * v
        k.append(d)
    return i, k


def add(d, key, v):
    d[key] = d[key] + v if key in d else v


class Integrals:
    """J_n(b; A, B) = int_A^B r^n e^{b r} dr as {key: coef}, key = (e0, kind, arg); A = 0 means finite part,
    B = None means infinity.  Cached per (b, A, B)."""

    def __init__(self, K):
        self.K = K
        self.cache = {}

    def endpoint(self, b, R, m):
        """[r^m e^{b r}] at r = R as {key: coef}; R = None -> 0; R = 0 -> finite part."""
        K = self.K
        if R is None:
            return {}
        if R == 0:
            if m > 0:
                return {}
            # coefficient of eps^0 in eps^m e^{b eps}
            f = Fraction(1)
            for t in range(1, -m + 1):
                f /= t
            return {(Fraction(0), '1', None): K(Fraction(b) ** (-m) * f)}
        return {(Fraction(b) * R, '1', None): K(Fraction(R) ** m)}

    def table(self, b, A, B, nmin, nmax):
        key = (b, A, B)
        tab = self.cache.get(key, {})
        K = self.K
        if -1 not in tab:
            J = {}
            if b == 0:
                if A == 0:
                    J = {(Fraction(0), 'ln', Fraction(B)): K(1)}
                else:
                    J = {(Fraction(0), 'ln', Fraction(B)): K(1), (Fraction(0), 'ln', Fraction(A)): K(-1)}
            else:
                if A == 0:
                    J = {(Fraction(0), 'g', None): K(-1), (Fraction(0), 'ln', abs(Fraction(b))): K(-1)}
                else:
                    J = {(Fraction(0), 'E1', -Fraction(b) * A): K(1)}
                if B is not None:
                    add(J, (Fraction(0), 'E1', -Fraction(b) * B), K(-1))
            tab[-1] = J
        n = max(tab) + 1
        while n <= nmax:                       # J_n = ([r^n e^{br}] - n J_{n-1}) / b   (b != 0)
            J = {}
            if b == 0:
                for kk, v in self.endpoint(0, B, n + 1).items():
                    add(J, kk, v / K(n + 1))
                for kk, v in self.endpoint(0, A, n + 1).items():
                    add(J, kk, -v / K(n + 1))
            else:
                ib = K(1) / K(Fraction(b))
                for kk, v in self.endpoint(b, B, n).items():
                    add(J, kk, v * ib)
                for kk, v in self.endpoint(b, A, n).items():
                    add(J, kk, -v * ib)
                for kk, v in tab[n - 1].items():
                    add(J, kk, -K(n) * v * ib)
            tab[n] = J
            n += 1
        n = min(tab) - 1
        while n >= nmin:                       # J_n = ([r^{n+1} e^{br}] - b J_{n+1}) / (n+1)
            J = {}
            inv = K(1) / K(n + 1)
            for kk, v in self.endpoint(b, B, n + 1).items():
                add(J, kk, v * inv)
            for kk, v in self.endpoint(b, A, n + 1).items():
                add(J, kk, -v * inv)
            if b != 0:
                for kk, v in tab[n + 1].items():
                    add(J, kk, -K(Fraction(b)) * v * inv)
            tab[n] = J
            n -= 1
        self.cache[key] = tab
        return tab


def coefficients(D1, D2, zX, zH, Lmax, K):
    """List over L of {constant key: q_{L,c}} with F_L / (2 pi N) = sum q c."""
    D1, D2, zX, zH = (Fraction(v) for v in (D1, D2, zX, zH))
    it, kt = bessel_tables(K, Lmax + 2)
    I = Integrals(K)
    lo, hi = sorted([D1, D2])
    pieces = [(Fraction(0), lo), (lo, hi), (hi, None)] if lo != hi else [(Fraction(0), lo), (lo, None)]
    out = []
    for L in range(Lmax + 1):
        res = {}
        for (A, B) in pieces:
            mid = A + 1 if B is None else (A + B) / 2
            inner = mid < D1                     # r < D1: x = zH r variable, y = zH D1 constant
            near = mid < D2
            # monomials {(e0, b, n): coef}:  e^{e0} r^n e^{b r}
            mono = {}

            def put(e0, b, n, v):
                add(mono, (e0, b, n), v)

            if inner:
                y = zH * D1
                kL = sum((K(y ** (-j)) * v for j, v in kt[L].items()), K(0))
                kL1 = sum((K(y ** (-j)) * v for j, v in kt[L + 1].items()), K(0))
                # (2L+1) i_L(x) kL + x i_{L+1}(x) kL - y i_L(x) kL1, times e^{-y}
                for s in (1, -1):
                    for j, v in it[L][s].items():
                        put(-y, s * zH, -j, K(zH ** (-j)) * v * (K(2 * L + 1) * kL - K(y) * kL1))
                    for j, v in it[L + 1][s].items():
                        put(-y, s * zH, 1 - j, K(zH ** (1 - j)) * v * kL)
            else:
                x = zH * D1
                for s in (1, -1):
                    iL = sum((K(x ** (-j)) * v for j, v in it[L][s].items()), K(0))
                    iL1 = sum((K(x ** (-j)) * v for j, v in it[L + 1][s].items()), K(0))
                    for j, v in kt[L].items():
                        put(s * x, -zH, -j, K(zH ** (-j)) * v * (K(2 * L + 1) * iL + K(x) * iL1))
                    for j, v in kt[L + 1].items():
                        put(s * x, -zH, 1 - j, -K(zH ** (1 - j)) * v * iL)
            # times -2 r^2 e^{-zX r} w_L(r)
            wpow = L if near else -L - 1
            wc = K(D2 ** (-L - 1)) if near else K(D2 ** L)
            for (e0, b, n), v in mono.items():
                bb, nn = b - zX, n + 2 + wpow
                tab = I.table(bb, A, B, nn, nn) if nn >= -1 else I.table(bb, A, B, nn, -1)
                for (e1, kind, arg), u in tab[nn].items():
                    add(res, (e0 + e1, kind, arg), K(-2) * wc * v * u)
        out.append({k: v for k, v in res.items() if v != 0})
    return out


def value(key, mp):
    e0, kind, arg = key
    f = mp.exp(mp.mpf(e0.numerator) / e0.denominator)
    if kind == '1':
        return f
    if kind == 'g':
        return f * mp.euler
    a = mp.mpf(arg.numerator) / arg.denominator
    if kind == 'ln':
        return f * mp.log(a)
    return f * (mp.e1(a) if a > 0 else -mp.ei(-a))
