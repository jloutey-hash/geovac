"""Helper for tests/test_paper58_exact_order.py -- order exactly 9 for the unequal-exponent exchange integrals.

Three independent pieces, none of which reads or imports debug/:

  1. `hand_table`: the derivative table of the nine functions as it is written in the independence proof (typed
     here by hand from the definitions of the nine functions, NOT produced by the helper calculus).  The test
     compares it with `_pachucki_extension.Calc.diff` identically over Q(t).
  2. Exact linear algebra for leg (D): coordinates of K, K', ..., K^(8) on the nine functions over Q(t) (Laurent
     in R), row denominators in t, and the determinant at a rational R as an INTEGER polynomial in t, obtained by
     evaluation at integer t (fraction-free Bareiss elimination on Python ints) and exact Newton interpolation.
     Pure Python; shares no code with the sympy DomainMatrix route or the mod-p route of the derivation and its
     review.
  3. mpmath versions of the nine functions for real (also irrational) exponents and the candidate log
     coefficients A_i of the proof's Step 1, with a difference operator that kills analytic functions to fifth
     order and so exposes a leftover (anything) x ln R.

Function order everywhere: _pachucki_extension.NAMES9 = f0, Ei(-(a+b)R), e^{-(a+b)R}, e^{-2aR}, e^{-2bR},
e^{2aR}Ei(-2(a+b)R), e^{2bR}Ei(-2(a+b)R), psi_a, psi_b   (indices 0..8 here, g_1..g_9 in the proof).
"""
from __future__ import annotations

from fractions import Fraction
from math import lcm
from typing import Dict, List, Sequence, Tuple

import sympy as sp

import _pachucki_extension as X


# ====================================================================== 1. the proof's derivative table
def hand_table(ctx: X.Ctx) -> List[Dict[int, Dict[int, object]]]:
    """Row i: coordinates of g_i' on the nine functions, {column: {power of R: coefficient}}.

    Typed from the definitions (c = 1/(2(a^2 - b^2)), psi_x = e^{-2xR}[Ei(2(x-y)R) - 2 Ei((x-y)R) - LC]):
        f0'    = -FM1/R,  FM1 = c (g7 - g6 - g9 + g8)           (Pachucki labels 59, 60)
        Ei(-(a+b)R)' = e^{-(a+b)R}/R
        (e^{2aR} Ei(-2(a+b)R))' = 2a (..) + e^{-2bR}/R          (and a <-> b)
        psi_a' = -2a psi_a + e^{-2bR}/R - 2 e^{-(a+b)R}/R       (and a <-> b)
    """
    a, b, one = ctx.a, ctx.b, ctx.one
    c = one / (2 * (a * a - b * b))
    return [
        {5: {-1: c}, 6: {-1: -c}, 7: {-1: -c}, 8: {-1: c}},
        {2: {-1: one}},
        {2: {0: -(a + b)}},
        {3: {0: -2 * a}},
        {4: {0: -2 * b}},
        {5: {0: 2 * a}, 4: {-1: one}},
        {6: {0: 2 * b}, 3: {-1: one}},
        {7: {0: -2 * a}, 4: {-1: one}, 2: {-1: -2 * one}},
        {8: {0: -2 * b}, 3: {-1: one}, 2: {-1: -2 * one}},
    ]


def helper_table(ctx: X.Ctx):
    """The same table from the reviewed calculus: in_nine(Calc.diff(g_i)).  None in a row = g_i' leaves the span."""
    calc = X.Calc(ctx)
    out = []
    for g, _ in X.nine_functions(ctx):
        co = X.in_nine(ctx, calc.diff(g))
        out.append(None if co is None else {j: dict(d) for j, d in enumerate(co) if d})
    return out


def fm1_on_nine(ctx: X.Ctx):
    """Coordinates of the helper's F^(-1) (label 60, Calc.FM1) on the nine functions, or None."""
    return X.in_nine(ctx, X.Calc(ctx).FM1)


# ====================================================================== 2. exact linear algebra for leg (D)
def stored_K(ctx: X.Ctx, cft: Dict[str, Dict[str, str]]) -> X.Func:
    """K(R; 1, t) as a Func over Q(t) from stored coefficient strings {name: {power of R: rational function}}."""
    assert ctx.kind == "Qt"
    K: X.Func = {}
    for (g, _), nm in zip(X.nine_functions(ctx), X.NAMES9):
        for p, s in cft.get(nm, {}).items():
            K = X.fadd(K, X.fscale(g, ctx.K.from_sympy(sp.sympify(s, locals={"t": ctx.t})), int(p)))
    return K


def derivative_rows(ctx: X.Ctx, K: X.Func, n: int = 8):
    """rows[j][i] = {power of R: coefficient} of K^(j) on g_i, j = 0..n, by the helper calculus.  Raises if a
    derivative leaves the span of the nine."""
    rows = []
    for j, d in enumerate(X.Calc(ctx).derivatives(K, n)):
        co = X.in_nine(ctx, d)
        if co is None:
            raise AssertionError(f"K^({j}) is not in the span of the nine functions")
        rows.append(co)
    return rows


def _frac(c) -> Fraction:
    return Fraction(int(c.numerator), int(c.denominator))


def linear_multiplicities(den) -> Tuple[Tuple[int, int, int], object]:
    """(multiplicities of t, t - 1, t + 1 in the polynomial `den`, what is left after dividing them out)."""
    t = den.ring.gens[0]
    mult = []
    for lin in (t, t - 1, t + 1):
        k = 0
        while True:
            q, r = den.div(lin)
            if r != 0:
                break
            den, k = q, k + 1
        mult.append(k)
    return tuple(mult), den


def denominator_audit(rows) -> Tuple[List[Tuple[int, int, int]], List[str]]:
    """Per row the largest multiplicities of t, t - 1, t + 1 over the denominators of its entries, and the list
    of every OTHER (non-constant) denominator cofactor met, as strings."""
    per_row, other = [], []
    for rw in rows:
        m = [0, 0, 0]
        for cdict in rw:
            for x in cdict.values():
                mult, rest = linear_multiplicities(x.denom)
                m = [max(u, v) for u, v in zip(m, mult)]
                if not rest.is_ground:
                    other.append(str(rest))
        per_row.append(tuple(m))
    return per_row, other


def clear_rows(ctx: X.Ctx, rows, per_row):
    """Row i times t^e0 (t-1)^e1 (t+1)^e2 (e = per_row[i]): entries {power of R: {degree in t: Fraction}}.
    Raises if an entry is still not a polynomial in t (a denominator factor other than t, t - 1, t + 1)."""
    tF = ctx.b
    out = []
    for rw, (e0, e1, e2) in zip(rows, per_row):
        D = tF ** e0 * (tF - 1) ** e1 * (tF + 1) ** e2
        crow = []
        for cdict in rw:
            ent = {}
            for p, x in cdict.items():
                y = x * D
                if not y.denom.is_ground:
                    raise AssertionError(f"entry not cleared by t, t - 1, t + 1: denominator {y.denom}")
                dc = _frac(y.denom.coeff(1))
                ent[p] = {e: _frac(cc) / dc for (e,), cc in y.numer.terms()}
            crow.append(ent)
        out.append(crow)
    return out


def integer_matrix_at(cleared, R0: Fraction):
    """The cleared 9 x 9 matrix at R = R0 as integer coefficient lists in t (ascending), each row scaled by the
    lcm of its coefficient denominators.  Returns (matrix, [row scales], degree bound of the determinant)."""
    M, scales, bound = [], [], 0
    for crow in cleared:
        polys = []
        for ent in crow:
            acc: Dict[int, Fraction] = {}
            for p, pol in ent.items():
                w = R0 ** p
                for e, cc in pol.items():
                    acc[e] = acc.get(e, Fraction(0)) + w * cc
            polys.append({e: v for e, v in acc.items() if v != 0})
        s = lcm(*[v.denominator for pol in polys for v in pol.values()]) if any(polys) else 1
        deg = max((max(pol) for pol in polys if pol), default=0)
        M.append([[int(pol.get(e, 0) * s) for e in range(max(pol) + 1)] if pol else [0] for pol in polys])
        scales.append(s)
        bound += deg
    return M, scales, bound


def _horner(cs: Sequence[int], x: int) -> int:
    acc = 0
    for c in reversed(cs):
        acc = acc * x + c
    return acc


def bareiss_det(M: List[List[int]]) -> int:
    """Determinant of an integer matrix, fraction-free (Bareiss); exact."""
    M = [r[:] for r in M]
    n, sign, prev = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1]


def fraction_det(M: List[List[Fraction]]) -> Fraction:
    """Determinant over Q by Gaussian elimination (a second, different elimination; used for the cross-check)."""
    M = [r[:] for r in M]
    n, d = len(M), Fraction(1)
    for i in range(n):
        pr = next((k for k in range(i, n) if M[k][i] != 0), None)
        if pr is None:
            return Fraction(0)
        if pr != i:
            M[i], M[pr] = M[pr], M[i]
            d = -d
        d *= M[i][i]
        for k in range(i + 1, n):
            if M[k][i] != 0:
                f = M[k][i] / M[i][i]
                M[k] = [x - f * y for x, y in zip(M[k], M[i])]
    return d


def interpolate_integer_polynomial(x0: int, ys: Sequence[int]) -> List[int]:
    """The polynomial of degree < len(ys) through (x0 + k, ys[k]), ascending integer coefficients.  Forward
    differences are exact integers; the binomial basis is expanded in Fractions.  Raises if a coefficient is not
    an integer (it must be, for the determinant of an integer polynomial matrix)."""
    n = len(ys)
    d = list(ys)
    lead = []
    for k in range(n):
        lead.append(d[0])
        d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
    # p(x) = sum_k lead[k] binom(x - x0, k);  Horner in the Newton basis, from the top
    poly = [Fraction(0)]
    for k in range(n - 1, -1, -1):
        # poly <- poly * (x - x0 - k) / (k + 1) + lead[k]
        new = [Fraction(0)] * (len(poly) + 1)
        for i, cc in enumerate(poly):
            new[i + 1] += cc
            new[i] -= cc * (x0 + k)
        poly = [v / (k + 1) for v in new]
        poly[0] += lead[k]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    if any(v.denominator != 1 for v in poly):
        raise AssertionError("interpolated determinant is not an integer polynomial")
    return [int(v) for v in poly]


def determinant_polynomial(cleared, R0: Fraction, extra: int = 6):
    """(N, scale): N(t) = det of the cleared matrix at R = R0 times `scale` (product of the integer row scales),
    as ascending integer coefficients.  Interpolated on (degree bound + 1) integer nodes and then CHECKED on
    `extra` further nodes."""
    M, scales, bound = integer_matrix_at(cleared, R0)
    n = bound + 1
    x0 = -(n // 2)

    def value(x: int) -> int:
        return bareiss_det([[_horner(cs, x) for cs in row] for row in M])

    N = interpolate_integer_polynomial(x0, [value(x0 + k) for k in range(n)])
    for x in range(x0 + n, x0 + n + extra):
        if _horner(N, x) != value(x):
            raise AssertionError("determinant polynomial fails at an extra node: degree bound too small")
    scale = 1
    for s in scales:
        scale *= s
    return N, scale


def strip_linear(N: List[int]) -> Tuple[Tuple[int, int, int], List[int]]:
    """Divide an integer polynomial (ascending) by t, t - 1, t + 1 as often as possible."""
    mult = []
    for r in (0, 1, -1):
        k = 0
        while len(N) > 1 and _horner(N, r) == 0:
            q = [0] * (len(N) - 1)                       # synthetic division by (t - r)
            acc = 0
            for i in range(len(N) - 1, 0, -1):
                acc = N[i] + acc * r
                q[i - 1] = acc
            N, k = q, k + 1
        mult.append(k)
    return tuple(mult), N


def common_factor(residuals: Sequence[List[int]]):
    """gcd over Q of the given integer polynomials as a sympy Poly in t."""
    t = sp.Symbol("t")
    g = sp.Poly(list(reversed(residuals[0])), t, domain="ZZ")
    for r in residuals[1:]:
        g = sp.gcd(g, sp.Poly(list(reversed(r)), t, domain="ZZ"))
    return g


def real_roots_above_one(g) -> int:
    """Number of distinct real roots of the sympy Poly g in the open interval (1, oo)."""
    if g.is_zero:
        raise AssertionError("zero polynomial")
    if g.degree() <= 0:
        return 0
    t = g.gens[0]
    n = int(g.count_roots(1, sp.oo))
    if g.eval(1) == 0:
        n -= 1
    assert sp.Poly(g, t).degree() == g.degree()
    return n


def rows_at(ctx: X.Ctx, rows, tv: Fraction, R0: Fraction) -> List[List[Fraction]]:
    """Rows over Q(t) (or over Q) evaluated at t = tv, R = R0."""
    return [[sum((X.qt_to_fraction(ctx, x, tv) * R0 ** p for p, x in cdict.items()), Fraction(0)) for cdict in rw]
            for rw in rows]


# ====================================================================== 3. numerics of the proof's Step 1
def nine_mp(a, b):
    """mpmath callables [g2, ..., g9] (indices 1..8 of NAMES9) and F^(-1), for real 0 < a < b (mpf; may be
    irrational).  Typed from the definitions; the test compares them with the helper's Func at rational pairs."""
    import mpmath as mp
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    g = [
        lambda R: mp.ei(-(a + b) * R),
        lambda R: mp.exp(-(a + b) * R),
        lambda R: mp.exp(-2 * a * R),
        lambda R: mp.exp(-2 * b * R),
        lambda R: mp.exp(2 * a * R) * mp.ei(-2 * (a + b) * R),
        lambda R: mp.exp(2 * b * R) * mp.ei(-2 * (a + b) * R),
        lambda R: mp.exp(-2 * a * R) * (mp.ei(2 * (a - b) * R) - 2 * mp.ei((a - b) * R) - LC),
        lambda R: mp.exp(-2 * b * R) * (mp.ei(2 * (b - a) * R) - 2 * mp.ei((b - a) * R) - LC),
    ]

    def fm1(R):
        return c * (g[5](R) - g[4](R) - g[7](R) + g[6](R))       # c (g7 - g6 - g9 + g8) in the proof's numbering

    return g, fm1


def log_coefficients(a, b):
    """The proof's A_2..A_9 (list of 8 callables, same order as nine_mp) and A_1, L."""
    import mpmath as mp
    c = 1 / (2 * (a * a - b * b))
    A = [
        lambda R: mp.mpf(1),
        lambda R: mp.mpf(0),
        lambda R: mp.mpf(0),
        lambda R: mp.mpf(0),
        lambda R: mp.exp(2 * a * R),
        lambda R: mp.exp(2 * b * R),
        lambda R: -mp.exp(-2 * a * R),
        lambda R: -mp.exp(-2 * b * R),
    ]

    def L(R):                      # log coefficient of F^(-1) = c (g7 - g6 - g9 + g8), from the A's above
        return c * (A[5](R) - A[4](R) - A[7](R) + A[6](R))

    def cinh(x):                   # Cinh(x) = int_0^x (cosh s - 1)/s ds = sum_{k>=1} x^{2k} / (2k (2k)!)
        s, term, k = mp.mpf(0), mp.mpf(1), 0
        while True:
            k += 1
            term = term * x * x / ((2 * k - 1) * (2 * k))
            add = term / (2 * k)
            s += add
            if abs(add) <= abs(s) * mp.mpf(10) ** (-mp.mp.dps - 8) or add == 0:
                return s

    def A1(R):
        return -2 * c * (cinh(2 * b * R) - cinh(2 * a * R))

    return A, A1, L


D5 = (-1, 5, -10, 10, -5, 1)       # fifth forward difference on the nodes e, 2e, ..., 6e


def fifth_difference(h, e):
    """sum_k D5[k] h((k+1) e).  O(e^5) for h analytic at 0; for h = R^m ln R (m = 0..4) it is e^m times a nonzero
    constant (0.0875, -0.060, 0.118, ... : the ln e part cancels, the ln k part does not)."""
    return sum(w * h((k + 1) * e) for k, w in enumerate(D5))


def f0_fifth_difference(fm1, e):
    """The fifth difference of f0 on e..6e from f0' = -F^(-1)/R alone: f0(ke) = f0(6e) + sum_{j=k}^{5} I_j,
    I_j = int_{je}^{(j+1)e} F^(-1)(s)/s ds; the weights sum to zero, so f0(6e) drops out."""
    import mpmath as mp
    I = [mp.quad(lambda s: fm1(s) / s, [j * e, (j + 1) * e]) for j in range(1, 6)]
    tot = mp.mpf(0)
    for k, w in enumerate(D5):               # node (k+1) e
        tot += w * sum(I[j - 1] for j in range(k + 1, 6))
    return tot
