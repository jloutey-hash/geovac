"""Shibuya-Wulfman (SW) metric between Coulomb-Sturmians of ANY (n, l, m) on
any number of centres, in any orientation (m-seam rung 3, added 2026-10-10).

Extends ``geovac/sturmian_sigma_basis.py`` (sigma basis, m = 0, two centres
on one axis, position-space prolate quadrature) to general l, m and
geometry, by an independent momentum-space route.

Conventions.  Shared scale k = 1 (distances are S = kD).  Basis functions
are L2-normalised Coulomb-Sturmians chi_{nl}(r) Y_a(rhat) (radial = hydrogenic
with Z = n k), with REAL normalised harmonics Y_a built from exact
polynomials.  The SW metric is the global inner product

    S_{ab} = <chi_a(. - X_a)| (p^2 + k^2) / (2 k^2) |chi_b(. - X_b)>,

so each centre's own block is the identity.  With chi~_{nl,a}(p) =
(-i)^l f_nl(p) Y_a(phat) and e^{i p.d} = sum_L i^L (2L+1) j_L(p d)
P_L(dhat.phat), d = X_a - X_b:

    S_ab(d) = sum_L i^{L + l_a - l_b} (2L+1) A_L(a, b; dhat) I^L(|d|),
    A_L = int Y_a Y_b P_L(dhat . phat) dOmega          (exact, Gaunt-terminated),
    I^L = int_0^oo f_a f_b (p^2 + 1)/2 j_L(p |d|) p^2 dp,

the last in CLOSED FORM by one residue at p = i (polynomial in D, 1/D times
e^{-D}, no pi), with mpmath quadrature kept as an independent referee.

f_nl is the hydrogen momentum radial (Bethe-Salpeter) at Z = n k, which for
Sturmians depends on p / k alone (shared scale).  The fixed-S kernel posing
of the one-electron many-centre problem is K = D^1/2 S D^1/2,
D = diag(Z_centre / n): eigenvalues are k, and only two-centre quantities
enter (Shibuya-Wulfman 1965; Avery & Avery 2017).
"""
from __future__ import annotations

from functools import lru_cache
from math import factorial
from typing import Dict, List, Sequence, Tuple

import mpmath as mp
import sympy as sp
from sympy import Rational, gamma

_x, _y, _z = sp.symbols('x y z', real=True)
_ux, _uy, _uz = sp.symbols('ux uy uz', real=True)

Label = Tuple[int, int, int, str]   # (n, l, m >= 0, 'c' | 's')


def _sphere_int(poly: sp.Expr) -> sp.Expr:
    total = sp.Integer(0)
    for (a, b, d), coef in sp.Poly(sp.expand(poly), _x, _y, _z).terms():
        if a % 2 or b % 2 or d % 2:
            continue
        total += coef * 2 * gamma(Rational(a + 1, 2)) * gamma(Rational(b + 1, 2)) \
            * gamma(Rational(d + 1, 2)) / gamma(Rational(a + b + d + 3, 2))
    return total


@lru_cache(maxsize=None)
def real_harmonic(l: int, m: int, kind: str) -> sp.Expr:
    """Normalised real harmonic polynomial (on the unit sphere):
    (x + i y)^m r^{l-m} (d^m P_l/dt^m)(z/r), cos or sin part."""
    t = sp.Symbol('t')
    q = sp.Poly(sp.diff(sp.legendre(l, t), t, m), t)
    rad = sp.expand(sum(coef * _z**j * (_x**2 + _y**2 + _z**2) ** ((l - m - j) // 2)
                        for (j,), coef in q.terms()))
    w = sp.expand((_x + sp.I * _y) ** m)
    p = sp.expand((sp.re(w) if kind == 'c' else sp.im(w)) * rad)
    return sp.expand(p / sp.sqrt(_sphere_int(p ** 2)))


def labels(n_max: int, l_max: int) -> List[Label]:
    out: List[Label] = []
    for n in range(1, n_max + 1):
        for l in range(min(n - 1, l_max) + 1):
            out.append((n, l, 0, 'c'))
            for m in range(1, l + 1):
                out += [(n, l, m, 'c'), (n, l, m, 's')]
    return out


@lru_cache(maxsize=None)
def _angular_fn(la: int, ma: int, ka: str, lb: int, mb: int, kb: str, L: int):
    """(2L+1) int Y_a Y_b P_L(u . rhat) dOmega as a function of u (unit)."""
    expr = _sphere_int(real_harmonic(la, ma, ka) * real_harmonic(lb, mb, kb)
                       * sp.expand(sp.legendre(L, _ux * _x + _uy * _y + _uz * _z)))
    expr = sp.expand((2 * L + 1) * expr)
    if expr == 0:
        return None
    return sp.lambdify((_ux, _uy, _uz), expr, 'mpmath')


def f_mom(n: int, l: int, p):
    """Momentum radial of the L2-normalised Sturmian (n, l) at k = 1."""
    p = mp.mpf(p)
    u = p * p
    pref = mp.sqrt(mp.mpf(2) / mp.pi * factorial(n - l - 1) / factorial(n + l)) \
        * 2 ** (2 * l + 2) * factorial(l)
    # hydrogen at Z = n: F(p) = n^{-3/2} F_H(p/n); F_H has n^2 n^l (p/n)^l / ((p/n)^2 n^2 + 1)^{l+2}
    xg = (u - 1) / (u + 1)
    return pref * mp.sqrt(n) * p ** l / (u + 1) ** (l + 2) * mp.gegenbauer(n - l - 1, l + 1, xg)


_p = sp.Symbol('p')
_D = sp.Symbol('D', positive=True)


def _f_exact(n: int, l: int) -> sp.Expr:
    pref = sp.sqrt(Rational(2) / sp.pi * Rational(factorial(n - l - 1), factorial(n + l)))         * 2 ** (2 * l + 2) * factorial(l)
    return pref * sp.sqrt(n) * _p ** l / (_p ** 2 + 1) ** (l + 2)         * sp.gegenbauer(n - l - 1, l + 1, (_p ** 2 - 1) / (_p ** 2 + 1))


def _h1(L: int, x: sp.Expr) -> sp.Expr:
    """Spherical Hankel h_L^(1)(x) (closed form)."""
    return (-sp.I) ** (L + 1) * sp.exp(sp.I * x) / x * sum(
        sp.I ** k / (sp.factorial(k) * (2 * x) ** k) * sp.factorial(L + k) / sp.factorial(L - k)
        for k in range(L + 1))


@lru_cache(maxsize=None)
def radial_I_exact(na: int, la: int, nb: int, lb: int, L: int) -> sp.Expr:
    """I^L(D) in closed form.  The integrand E(p) j_L(pD), E = f_a f_b (p^2+1)/2 p^2,
    is even and E ~ p^{la+lb+2} at 0 (L <= la+lb), so int_0^oo E j_L =
    (1/2) Re int_R E h_L^(1)(pD) dp = Re[pi i Res_{p=i} E h_L^(1)(pD)]:
    a polynomial in D, 1/D times e^{-D}; pi cancels."""
    E = sp.cancel(_f_exact(na, la) * _f_exact(nb, lb) * (_p ** 2 + 1) / 2 * _p ** 2)
    res = _residue_at_i(E, L)
    return sp.simplify(sp.re(sp.expand(sp.pi * sp.I * res)))


def _trunc_mul(a: List[sp.Expr], b: List[sp.Expr], order: int) -> List[sp.Expr]:
    out = [sp.Integer(0)] * (order + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(order + 1 - i):
            if j < len(b) and b[j] != 0:
                out[i + j] += ai * b[j]
    return [sp.expand(v) for v in out]


def _residue_at_i(E: sp.Expr, L: int) -> sp.Expr:
    """Res_{p=i} E(p) h_L^(1)(pD) for E = num / (c (p^2+1)^m), by truncated
    Taylor series in t = p - i (polynomial arithmetic only)."""
    num, den = sp.fraction(sp.together(E))
    dpoly = sp.Poly(den, _p)
    m = dpoly.degree() // 2
    c = sp.Poly(den, _p).LC()
    assert sp.expand(den - c * (_p ** 2 + 1) ** m) == 0, "denominator is not c (p^2+1)^m"
    t = sp.Symbol('t')
    order = m - 1
    nser = sp.Poly(sp.expand(num.subs(_p, sp.I + t)), t).all_coeffs()[::-1]
    nser = (nser + [sp.Integer(0)] * (order + 1))[:order + 1]
    # (2i + t)^{-m}
    pole = [(2 * sp.I) ** (-m) * sp.binomial(-m, j) * (1 / (2 * sp.I)) ** j for j in range(order + 1)]
    expo = [(sp.I * _D) ** j / sp.factorial(j) for j in range(order + 1)]
    hk = [sp.Integer(0)] * (order + 1)
    for k in range(L + 1):
        ak = sp.I ** k * sp.factorial(L + k) / (sp.factorial(k) * sp.factorial(L - k) * 2 ** k)
        ser = [sp.I ** (-(k + 1)) * sp.binomial(-(k + 1), j) * (1 / sp.I) ** j for j in range(order + 1)]
        for j in range(order + 1):
            hk[j] += ak * _D ** (-(k + 1)) * ser[j]
    acc = _trunc_mul(_trunc_mul(_trunc_mul(nser, pole, order), expo, order), hk, order)
    return (-sp.I) ** (L + 1) * sp.exp(-_D) * acc[order] / c


@lru_cache(maxsize=None)
def _radial_fn(na: int, la: int, nb: int, lb: int, L: int):
    expr = radial_I_exact(na, la, nb, lb, L)
    if L > 0:
        at0 = sp.Integer(0)          # j_L(0) = 0 for L > 0
    else:                            # j_0(0) = 1: int_0^oo E dp = Re[pi i Res_{p=i} E]
        E = sp.cancel(_f_exact(na, la) * _f_exact(nb, lb) * (_p ** 2 + 1) / 2 * _p ** 2)
        at0 = sp.simplify(sp.re(sp.expand(sp.pi * sp.I * sp.residue(E, _p, sp.I))))
    fn = sp.lambdify(_D, expr, 'mpmath')
    return fn, at0


def radial_I(na: int, la: int, nb: int, lb: int, L: int, D: str):
    """I^L(D) from the exact closed form (D as an exact decimal string)."""
    Dv = mp.mpf(D)
    fn, at0 = _radial_fn(na, la, nb, lb, L)
    return mp.mpf(sp.N(at0, mp.mp.dps + 5)) if Dv == 0 else fn(Dv)


@lru_cache(maxsize=None)
def radial_I_quad(na: int, la: int, nb: int, lb: int, L: int, D: str):
    """Independent referee: I^L(D) by mpmath quadrature of the momentum integral."""
    Dv = mp.mpf(D)
    def g(p):
        jl = mp.mpf(1) if (Dv == 0 and L == 0) else (
            mp.mpf(0) if Dv == 0 else mp.sqrt(mp.pi / (2 * p * Dv)) * mp.besselj(L + mp.mpf(1) / 2, p * Dv))
        return f_mom(na, la, p) * f_mom(nb, lb, p) * (p * p + 1) / 2 * jl * p * p
    if Dv == 0:
        return mp.quad(g, [0, 1, 4, mp.inf])
    pts = [0] + [mp.pi * j / Dv for j in range(1, 40)] + [mp.inf]
    return mp.quad(g, pts)


def _key(D) -> str:
    return mp.nstr(mp.mpf(D), mp.mp.dps)


def sw_block(ba: Sequence[Label], bb: Sequence[Label], d: Sequence) -> mp.matrix:
    """SW metric block between basis ba at X_a and bb at X_b, d = X_a - X_b."""
    d = [mp.mpf(v) for v in d]
    D = mp.sqrt(sum(v * v for v in d))
    u = [v / D for v in d] if D != 0 else [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    M = mp.matrix(len(ba), len(bb))
    for i, (na, la, ma, ka) in enumerate(ba):
        for j, (nb, lb, mb, kb) in enumerate(bb):
            val = mp.mpf(0)
            for L in range(abs(la - lb), la + lb + 1, 2):
                fn = _angular_fn(la, ma, ka, lb, mb, kb, L)
                if fn is None:
                    continue
                phase = (-1) ** ((L + la - lb) // 2)
                val += phase * fn(*u) * radial_I(na, la, nb, lb, L, _key(D))
            M[i, j] = val
    return M


def assemble(centres: Sequence[Tuple[Sequence, float, Sequence[Label]]]
             ) -> Tuple[mp.matrix, List[mp.mpf], List[Tuple[int, Label]]]:
    """Global SW metric S, diag D = Z_c / n, and (centre, label) index.
    ``centres`` = [(X, Z, labels), ...] with X in units of 1/k."""
    index = [(c, lab) for c, (_, _, labs) in enumerate(centres) for lab in labs]
    N = len(index)
    S = mp.matrix(N, N)
    offs, o = [], 0
    for _, _, labs in centres:
        offs.append(o)
        o += len(labs)
    for ca, (Xa, _, la) in enumerate(centres):
        for cb, (Xb, _, lb) in enumerate(centres):
            if cb < ca:
                continue
            blk = sw_block(la, lb, [Xa[t] - Xb[t] for t in range(3)])
            for i in range(len(la)):
                for j in range(len(lb)):
                    S[offs[ca] + i, offs[cb] + j] = blk[i, j]
                    S[offs[cb] + j, offs[ca] + i] = blk[i, j]
    Dg = [mp.mpf(centres[c][1]) / lab[0] for c, lab in index]
    return S, Dg, index
