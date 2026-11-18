"""OP52 derivation memo, lemma checks (sympy, exact).

L1/L3  the half-integer exponents at x* of the order-2 operator with solutions s B1, s B2 (s = (x* - x)^{1/2},
       B1(x*) = B2(x*) != 0 after normalization) are {1/2 + ord B1, 1/2 + ord B1 + v}, v = ord(B2/B2(x*) - B1/B1(x*))
       -- the collision resolves at the valuation of the difference, not at +1 by fiat (v = 1, 2, 3 shown).
M1 is not automatic: ONE function f = s log(x1 - x) (x1 != x*) has minimal operator of order 2 with TWO half-integer
       exponents {1/2, 3/2} at x*: its continuation around x1 adds 2 pi i s, so the (-1)-eigenspace is 2-dimensional.
       Contrast f = s r(x) + a(x), r rational: one half-integer exponent (rank one, record OP52 K6b).
"""
from __future__ import annotations

import sympy as sp

x = sp.Symbol('x')
XS = sp.Rational(5, 3)
X1 = sp.Integer(-2)


def op2(y1, y2):
    """coefficients (p0, p1) of y'' + p1 y' + p0 y = 0 with solutions y1, y2."""
    W = sp.simplify(y1 * sp.diff(y2, x) - y2 * sp.diff(y1, x))
    p1 = sp.simplify(-sp.diff(W, x) / W)
    p0 = sp.simplify(-(sp.diff(y1, x, 2) + p1 * sp.diff(y1, x)) / y1)
    for y in (y1, y2):
        assert sp.simplify(sp.diff(y, x, 2) + p1 * sp.diff(y, x) + p0 * y) == 0
    return p0, p1


def indicial_at(p0, p1, x0):
    """roots of n(n-1) + a n + b, t = x0 - x:  y'' + p1 y' + p0 y with p1 ~ -a/(x0-x)... computed via the limit."""
    t = sp.Symbol('t')
    a = sp.limit(sp.simplify(-(p1.subs(x, x0 - t)) * t), t, 0)       # d/dx = -d/dt
    b = sp.limit(sp.simplify(p0.subs(x, x0 - t) * t ** 2), t, 0)
    n = sp.Symbol('n')
    return sorted(sp.solve(n * (n - 1) + a * n + b, n))


def main():
    s = sp.sqrt(XS - x)
    rows = []
    for v in (1, 2, 3):
        B1 = (XS - x) * (1 + x)                         # ord 1
        B2 = (XS - x) * (1 + x) * (1 + (XS - x) ** v / 7)   # same lead, differs at relative order v
        p0, p1 = op2(s * B1, s * B2)
        rows.append((v, indicial_at(p0, p1, XS)))
    for v, r in rows:
        print('collision, relative order v =', v, '-> exponents at x*:', r)
        assert r == [sp.Rational(3, 2), sp.Rational(3, 2) + v]
    # M1 counterexample: one function, two branches
    f = s * sp.log(x - X1)
    g = s                                              # the continuation difference (f around x1 minus f) / (2 pi i)
    p0, p1 = op2(f, g)
    r = indicial_at(p0, p1, XS)
    print('f = s log(x - x1): exponents at x*:', r)
    assert r == [sp.Rational(1, 2), sp.Rational(3, 2)]
    print('all lemma checks pass')


if __name__ == '__main__':
    main()
