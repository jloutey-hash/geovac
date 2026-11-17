"""Depth-2 integration-by-parts reducer for record FX (design: debug/fx_design_memo.md section 4).

Reduces definite integrals over one radial panel [A, B] (B = None for infinity; A = 0 means Hadamard finite part)

    int_A^B e^{lam r} X(r) E(r) dr,   X rational with poles at 0 and finitely many p,
    E = product of 0, 1 or 2 factors E1(c r + d)   (E1(beta r) is the factor (beta, 0)),

to a finite list of constants, exactly over a pluggable field K (na_exact.Field: fmpq, or nmod mod a prime).

Constants are keyed  CK = (e0, atoms):  e^{e0} * prod(atoms),  atoms a sorted tuple of
    ('E1', x)        E1(x), x rational (x < 0 means -Ei(-x), as na_exact.value)
    ('ln', v)        ln v, v > 0 rational          ('g',)   Euler's gamma
    ('B', lam, E, p, A, B)   the named base  int_A^B e^{lam r} E(r) / (r - p) dr  (len(E) >= 1; A = 0 = p finite part)
    any other atom is opaque to the reducer and carried through (fx_exact uses a tag atom).
The depth-0 base int e^{lam r}/(r - p) is reduced to E1 / ln / gamma keys in na_exact's convention.

Rational functions are kept in partial-fraction form PF = (poly part, {pole: [_, c_1, .., c_K]}), c_k the coefficient
of (r - pole)^{-k}.  The reduction, per bucket (lam, E, CK), with t = r - p:
    polynomial P:   e^{lam r} S antiderivative of e^{lam r} P  (S' + lam S = P; S = int P at lam = 0),
    principal T:    e^{lam r} U + C int e^{lam r}/t   (U' + lam U = T - C/t),
    int e^{lam r} X E = [e^{lam r} Y E] - int e^{lam r} Y E'  (Y = S or U),
    E1(c r + d)' = -e^{-d} e^{-c r} / (r - q), q = -d/c,  so every derivative lowers the depth by one and adds a pole.
Boundary terms are evaluated exactly; at A = 0 the finite part is the constant term of the expansion in r with every
ln r dropped, using E1(beta r) = -gamma - ln beta - ln r + Ein(beta r) and the Taylor series of E1(c r + d) at d != 0.
"""
from __future__ import annotations

from fractions import Fraction as Fr

import flint


# ---------------------------------------------------------------- field-generic polynomials
class Ring:
    def __init__(self, K):
        self.K = K
        self.p = K.p
        self.one, self.zero = K(1), K(0)

    def poly(self, coeffs):
        if self.p is None:
            return flint.fmpq_poly(list(coeffs))
        return flint.nmod_poly([int(c) for c in coeffs], self.p)

    def lin(self, c0, c1):
        return self.poly([self.K(c0), self.K(c1)])

    def coeffs(self, P):
        return list(P.coeffs()) if P.degree() >= 0 else []

    def is0(self, v):
        return v == 0


def nz(v):
    return v != 0


# ---------------------------------------------------------------- partial fractions
class PF:
    __slots__ = ('poly', 'pp')

    def __init__(self, poly, pp=None):
        self.poly = poly
        self.pp = pp or {}

    def add(self, other):
        self.poly = self.poly + other.poly
        for p, c in other.pp.items():
            mine = self.pp.get(p)
            if mine is None:
                self.pp[p] = list(c)
            else:
                if len(mine) < len(c):
                    mine.extend([c[0] * 0] * (len(c) - len(mine)))
                for k in range(1, len(c)):
                    mine[k] = mine[k] + c[k]


def _series_inv_linear_pow(R, p, n, M):
    """Coefficients j < M of (r - p)^{-n} at r = 0 (p != 0):  (-p)^{-n} binom(n+j-1, j) p^{-j}."""
    K = R.K
    s0 = K(Fr(-1) ** n / Fr(p) ** n)
    out = [s0]
    ip = K(1 / Fr(p))
    for j in range(1, M):
        out.append(out[-1] * K(n + j - 1) / K(j) * ip)
    return out


def pf_convert(R, N, I, p=None, Kp=0):
    """PF of N(r) / (r^I (r - p)^Kp); I may be negative (extra positive powers)."""
    K = R.K
    if I < 0:
        N = N.left_shift(-I)
        I = 0
    if N.degree() < 0:
        return PF(R.poly([]))
    if p is None or Kp == 0:
        cs = R.coeffs(N)
        pp = {}
        if I > 0:
            low = cs[:I] + [K(0)] * max(0, I - len(cs))
            pp[Fr(0)] = [K(0)] + [low[I - k] for k in range(1, I + 1)]
        return PF(R.poly(cs[I:]) if len(cs) > I else R.poly([]), pp)
    D = R.poly([K(1)]).left_shift(I) * (R.lin(-Fr(p), 1) ** Kp)
    Q, Rem = divmod(N, D)
    pp = {}
    if Rem.degree() >= 0:
        if I > 0:
            ser = R.poly(_series_inv_linear_pow(R, p, Kp, I))
            low = Rem.mul_low(ser, I)
            lc = R.coeffs(low) + [K(0)] * I
            pp[Fr(0)] = [K(0)] + [lc[I - k] for k in range(1, I + 1)]
        sh = Rem(R.lin(p, 1))                               # Rem(t + p)
        ser = R.poly(_series_inv_linear_pow(R, -Fr(p), I, Kp)) if I > 0 else R.poly([K(1)])
        low = sh.mul_low(ser, Kp)
        lc = R.coeffs(low) + [K(0)] * Kp
        pp[Fr(p)] = [K(0)] + [lc[Kp - k] for k in range(1, Kp + 1)]
    return PF(Q, pp)


# ---------------------------------------------------------------- constant keys
def ck_mul(c1, c2):
    return (c1[0] + c2[0], tuple(sorted(c1[1] + c2[1], key=repr)))


def ck_shift(ck, de0):
    return (ck[0] + de0, ck[1])


def add(d, k, v):
    if k in d:
        d[k] = d[k] + v
    else:
        d[k] = v


# ---------------------------------------------------------------- the reducer
class Reducer:
    """Buckets for one panel [A, B]; add terms with put(), then reduce() -> {CK: coefficient}."""

    def __init__(self, R, A, B):
        self.R, self.K = R, R.K
        self.A, self.B = Fr(A), (None if B is None else Fr(B))
        self.buckets = {}
        self.out = {}
        self.stats = dict(buckets=0, named=0)

    def put(self, lam, E, ck, pf):
        key = (Fr(lam), tuple(sorted(E)), ck)
        if key in self.buckets:
            self.buckets[key].add(pf)
        else:
            self.buckets[key] = pf

    def put_rational(self, lam, E, ck, N, I, p=None, Kp=0):
        self.put(lam, E, ck, pf_convert(self.R, N, I, p, Kp))

    def emit(self, ck, v):
        if nz(v):
            add(self.out, ck, v)

    # -- antiderivatives -------------------------------------------------------------
    def _poly_anti(self, lam, P):
        K = self.K
        if lam == 0:
            return P.integral()
        cs = self.R.coeffs(P)
        if not cs:
            return P
        il = K(1) / K(lam)
        s = [K(0)] * len(cs)
        s[-1] = cs[-1] * il
        for k in range(len(cs) - 2, -1, -1):
            s[k] = (cs[k] - K(k + 1) * s[k + 1]) * il
        return self.R.poly(s)

    def _pp_anti(self, lam, c):
        """U (list, index k -> coeff of t^-k) and C for principal part c (list with c[0] unused)."""
        K = self.K
        Kt = len(c) - 1
        U = [K(0)] * max(Kt, 1)
        if Kt == 0:
            return U, K(0)
        if lam == 0:
            for k in range(2, Kt + 1):
                U[k - 1] = c[k] / K(1 - k)
            return U, c[1]
        L = K(lam)
        if Kt >= 2:
            U[Kt - 1] = -c[Kt] / K(Kt - 1)
            for j in range(Kt - 1, 1, -1):
                U[j - 1] = (L * U[j] - c[j]) / K(j - 1)
            return U, c[1] - L * U[1]
        return U, c[1]

    # -- boundary values ---------------------------------------------------------------
    def _boundary(self, lam, E, ck, Y, pole, R, sign):
        """sign * [e^{lam r} Y(r) E(r)] at r = R.  Y: ring poly (pole None) or principal list at pole."""
        K = self.K
        if R is None:
            assert lam - sum(c for c, _ in E) < 0, ('non-decaying term at infinity', lam, E)
            return
        if R == 0:
            return self._boundary_fp(lam, E, ck, Y, pole, sign)
        x = K(R)
        if pole is None:
            val = Y(x)
        else:
            t = K(R - pole)
            it, acc, pw = K(1) / t, K(0), K(1)
            for k in range(1, len(Y)):
                pw = pw * it
                acc += Y[k] * pw
            val = acc
        if not nz(val):
            return
        atoms = tuple(('E1', c * R + d) for c, d in E)
        for c, d in E:
            assert c * R + d != 0
        self.emit(ck_mul(ck, (lam * R, atoms)), K(sign) * val)

    def _E_series(self, E, M):
        """Non-log part of E(r) at r = 0 as {(e0, atoms): coefficient list j < M}."""
        K = self.K
        out = {(Fr(0), ()): [K(1)] + [K(0)] * (M - 1)}
        for c, d in E:
            fac = {}
            if d == 0:                                        # E1(c r): -gamma - ln c + Ein(c r)
                one = [K(1)] + [K(0)] * (M - 1)
                fac[(Fr(0), (('g',),))] = [-v for v in one]
                fac[(Fr(0), (('ln', Fr(c)),))] = [-v for v in one]
                ein = [K(0)] * M
                f = Fr(1)
                for k in range(1, M):
                    f = f * k
                    ein[k] = K(Fr(-1) ** (k + 1) * Fr(c) ** k / (k * f))
                fac[(Fr(0), ())] = ein
            else:                                             # E1(d + c r) = E1(d) - e^{-d} int_0^r c e^{-cs}/(d+cs) ds
                one = [K(1)] + [K(0)] * (M - 1)
                fac[(Fr(0), (('E1', Fr(d)),))] = one
                if M > 1:
                    ex = [K(1)]
                    for j in range(1, M):
                        ex.append(ex[-1] * K(Fr(-c) / j))
                    geo = [K(Fr(c) / d)]
                    for j in range(1, M):
                        geo.append(geo[-1] * K(Fr(-c) / d))
                    pr = self.R.coeffs(self.R.poly(ex).mul_low(self.R.poly(geo), M - 1)) + [K(0)] * M
                    tau = [K(0)] + [-pr[j - 1] / K(j) for j in range(1, M)]
                    fac[(Fr(-d), ())] = tau
            new = {}
            for (e1, a1), s1 in out.items():
                for (e2, a2), s2 in fac.items():
                    pr = self.R.coeffs(self.R.poly(s1).mul_low(self.R.poly(s2), M)) + [K(0)] * M
                    key = (e1 + e2, tuple(sorted(a1 + a2, key=repr)))
                    if key in new:
                        new[key] = [x + y for x, y in zip(new[key], pr[:M])]
                    else:
                        new[key] = pr[:M]
            out = new
        return out

    def _boundary_fp(self, lam, E, ck, Y, pole, sign):
        """sign * (constant term of e^{lam r} Y(r) E(r) at r -> 0, ln r dropped)."""
        K = self.K
        if pole is None:
            cs = self.R.coeffs(Y)
            ycoef = {0: cs[0]} if cs else {}
        elif pole == 0:
            ycoef = {k: Y[k] for k in range(1, len(Y)) if nz(Y[k])}
        else:
            t = K(-pole)
            it, acc, pw = K(1) / t, K(0), K(1)
            for k in range(1, len(Y)):
                pw = pw * it
                acc += Y[k] * pw
            ycoef = {0: acc}
        if not ycoef:
            return
        M = max(ycoef) + 1
        ser = self._E_series(E, M)
        ex = [K(1)]
        for j in range(1, M):
            ex.append(ex[-1] * K(Fr(lam) / j))
        for (e0, atoms), s in ser.items():
            full = self.R.coeffs(self.R.poly(s).mul_low(self.R.poly(ex), M)) + [K(0)] * M
            v = sum((ycoef[k] * full[k] for k in ycoef), K(0))
            self.emit(ck_mul(ck, (e0, atoms)), K(sign) * v)

    # -- base integrals ------------------------------------------------------------------
    def _base(self, lam, E, ck, p, C):
        K = self.K
        A, B = self.A, self.B
        if not E:
            if lam != 0:
                if A == 0 and p == 0:
                    self.emit(ck_mul(ck, (Fr(0), (('g',),))), -C)
                    self.emit(ck_mul(ck, (Fr(0), (('ln', abs(Fr(lam))),))), -C)
                else:
                    self.emit(ck_mul(ck, (lam * p, (('E1', -lam * (A - p)),))), C)
                if B is not None:
                    self.emit(ck_mul(ck, (lam * p, (('E1', -lam * (B - p)),))), -C)
                else:
                    assert lam < 0
            else:
                assert B is not None
                self.emit(ck_mul(ck, (Fr(0), (('ln', abs(B - p)),))), C)
                if not (A == 0 and p == 0):
                    self.emit(ck_mul(ck, (Fr(0), (('ln', abs(A - p)),))), -C)
            return
        self.stats['named'] += 1
        self.emit(ck_mul(ck, (Fr(0), (('B', Fr(lam), E, Fr(p), A, B),))), C)

    # -- derivative terms --------------------------------------------------------------------
    def _deriv(self, lam, E, ck, Y, pole):
        """- int e^{lam r} Y E'  ->  buckets of depth len(E) - 1."""
        K = self.K
        for i, (c, d) in enumerate(E):
            rest = E[:i] + E[i + 1:]
            q = Fr(-d) / c
            lam2, ck2 = lam - c, ck_shift(ck, -Fr(d))
            if pole is None:                                  # S(r)/(r - q)
                Q, rem = divmod(Y, self.R.lin(-q, 1))
                cs = self.R.coeffs(rem)
                pp = {q: [K(0), cs[0]]} if cs else {}
                self.put(lam2, rest, ck2, PF(Q, pp))
            elif pole == q:
                self.put(lam2, rest, ck2, PF(self.R.poly([]), {q: [K(0), K(0)] + list(Y[1:])}))
            else:                                             # U(t)/(t + delta), t = r - pole
                delta = Fr(pole) - q
                idl = K(1) / K(delta)
                Kt = len(Y) - 1
                w = [K(0)] * (Kt + 2)
                for j in range(Kt, 0, -1):
                    w[j] = (Y[j] - w[j + 1]) * idl
                mt = K(-delta)
                imt, acc, pw = K(1) / mt, K(0), K(1)
                for k in range(1, Kt + 1):
                    pw = pw * imt
                    acc += Y[k] * pw
                self.put(lam2, rest, ck2, PF(self.R.poly([]), {Fr(pole): w[:Kt + 1], q: [K(0), acc]}))

    # -- main loop ------------------------------------------------------------------------------
    def reduce(self):
        A, B = self.A, self.B
        while self.buckets:
            key = max(self.buckets, key=lambda k: len(k[1]))
            pf = self.buckets.pop(key)
            lam, E, ck = key
            self.stats['buckets'] += 1
            if pf.poly.degree() >= 0:
                S = self._poly_anti(lam, pf.poly)
                self._boundary(lam, E, ck, S, None, B, 1)
                self._boundary(lam, E, ck, S, None, A, -1)
                if E:
                    self._deriv(lam, E, ck, S, None)
            for p, c in pf.pp.items():
                if not any(nz(v) for v in c[1:]):
                    continue
                U, C = self._pp_anti(lam, c)
                Ul = list(U)                            # index k -> t^-k (Ul[0] unused)
                if any(nz(v) for v in Ul[1:]):
                    self._boundary(lam, E, ck, Ul, p, B, 1)
                    self._boundary(lam, E, ck, Ul, p, A, -1)
                    if E:
                        self._deriv(lam, E, ck, Ul, p)
                if nz(C):
                    self._base(lam, E, ck, p, C)
        return {k: v for k, v in self.out.items() if nz(v)}


# ---------------------------------------------------------------- canonical form, reduction mod p, numerics
def _ln_primes(v):
    v = Fr(v)
    out = {}
    for n, s in ((v.numerator, 1), (v.denominator, -1)):
        if n > 1:
            for q, e in flint.fmpz(n).factor():
                out[int(q)] = out.get(int(q), 0) + s * e
    return out


def canon(res, K):
    """Expand every ln atom over primes (ln of rationals are Q-dependent); drop zeros."""
    out = {}
    for (e0, atoms), v in res.items():
        terms = [((), K(1))]
        for at in atoms:
            if at[0] == 'ln':
                lp = _ln_primes(at[1])
                terms = [(a + (('ln', Fr(q)),), w * K(e)) for a, w in terms for q, e in lp.items()]
            else:
                terms = [(a + (at,), w) for a, w in terms]
        for a, w in terms:
            add(out, (e0, tuple(sorted(a, key=repr))), v * w)
    return {k: v for k, v in out.items() if nz(v)}


def reduce_mod(res_q, p):
    out = {}
    for k, v in res_q.items():
        x = flint.nmod(int(v.p), p) / flint.nmod(int(v.q), p)
        if int(x):
            out[k] = x
    return out


def _mpq(mp, x):
    x = Fr(x)
    return mp.mpf(x.numerator) / x.denominator


def _E1(mp, x):
    return mp.e1(x) if x > 0 else -mp.ei(-x)


def base_value(mp, lam, E, p, A, B):
    lam, p = _mpq(mp, lam), _mpq(mp, p)
    E = [(_mpq(mp, c), _mpq(mp, d)) for c, d in E]

    def h(r):
        v = mp.exp(lam * r)
        for c, d in E:
            v *= _E1(mp, c * r + d)
        return v
    Bv = mp.inf if B is None else _mpq(mp, B)
    if not (A == 0 and p == 0):
        return mp.quad(lambda r: h(r) / (r - p), [_mpq(mp, A), Bv])
    # finite part at 0: h = h0 + h1 ln r + o(1)
    h0, h1 = mp.mpf(1), mp.mpf(0)
    for c, d in E:
        if d == 0:
            a0, a1 = -mp.euler - mp.log(c), mp.mpf(-1)
        else:
            a0, a1 = _E1(mp, d), mp.mpf(0)
        h0, h1 = h0 * a0, h0 * a1 + h1 * a0
    top = mp.mpf(1) if B is None else min(Bv, mp.mpf(1))
    v = mp.quad(lambda r: (h(r) - h0 - h1 * mp.log(r)) / r, [0, top])
    v += h0 * mp.log(top) + h1 * mp.log(top) ** 2 / 2
    if B is None or Bv > 1:
        v += mp.quad(lambda r: h(r) / r, [1, Bv])
    return v


def atom_value(mp, at, cache):
    if at in cache:
        return cache[at]
    kind = at[0]
    if kind == 'E1':
        v = _E1(mp, _mpq(mp, at[1]))
    elif kind == 'ln':
        v = mp.log(_mpq(mp, at[1]))
    elif kind == 'g':
        v = +mp.euler
    elif kind == 'B':
        v = base_value(mp, *at[1:])
    else:
        raise KeyError(at)
    cache[at] = v
    return v


def value(res, mp, cache=None):
    cache = {} if cache is None else cache
    s = mp.mpf(0)
    for (e0, atoms), v in res.items():
        t = mp.exp(_mpq(mp, e0)) * mp.mpf(int(v.p)) / int(v.q)
        for at in atoms:
            t *= atom_value(mp, at, cache)
        s += t
    return s
