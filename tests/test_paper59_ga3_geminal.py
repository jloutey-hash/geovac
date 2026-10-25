"""Paper 59 sec:density_ring: the cusp-correct geminal f = r e^{-gamma r} on both edges at gamma = 2 brings
no new type (GA3).

Record docs/predictions/2026-10-23_geminal_shapes_g2.md (committed 037c1ea2 before computing).
Written as a separate step from the result; no debug/ import.  Data tests/data/paper59_ga3_2026_10_23.json.

Three guards.
(1) The stored order-two composition factor of Y_lin's annihilator at gamma = 2 is the image of the Legendre
    operator with z = (t^2 - 3)/(2t), twisted by x^(1/2), under the stored exact map.  The source is rebuilt
    here from the Legendre equation.  Wrong answers rejected: the gamma = 3/2 source; the real-momentum sign;
    a map perturbed by 1/1000 in one coefficient.
(2) That factor is a different operator from the exponential shape's order-two factor at gamma = 2 (not
    proportional; x-degree 36 against 30), although both map from the same Legendre operator.
(3) The four order-one factors are algebraic: each stored closed-form solution is killed by its factor, and the
    two branched at t = 4 carry the exponents -11/2 and -13/2 there, the one branched at t = 2 the exponent
    -5/2.  Wrong answers rejected: a solution offered to the wrong factor; a wrong exponent.
Not re-run here: the order-nine annihilator of Y_lin, its minimality bound, the K3 map and the quadrature on
Y_lin (record only).
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from test_paper59_c2b_elliptic import _residual

x = sp.Symbol("x")                      # plain, as in the shared residual helper
xp = sp.Symbol("x", positive=True)      # for the algebraic solutions
t = sp.Symbol("t", positive=True)
DATA = json.loads((Path(__file__).parent / "data" / "paper59_ga3_2026_10_23.json").read_text())
GAMMA = sp.Integer(2)
Z_SIGN = -1


def _poly(s, var=None):
    v = x if var is None else var
    return sp.expand(sp.sympify(str(s).replace("^", "**"), locals={"x": v}))


def _dform(theta_coeffs):
    c0, c1, c2 = [_poly(v) for v in theta_coeffs]
    return [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]          # theta = 2x d/dx


def _source(gamma=None, z_sign=None):
    gam = GAMMA if gamma is None else gamma
    zs = Z_SIGN if z_sign is None else z_sign
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    zz = (t**2 + 1 + zs * gam**2) / (2 * t)
    rt = sp.sqrt(t)
    W, W1, W2 = rt * h0, h0 / (2 * rt) + rt * h1, -h0 / (4 * t * rt) + h1 / rt + rt * h2   # w = t^(1/2) H
    zp, zpp = sp.diff(zz, t), sp.diff(zz, t, 2)
    w1 = W1 / zp
    w2 = (W2 - zpp * w1) / zp**2
    num = sp.expand(sp.numer(sp.together(sp.expand(((1 - zz**2) * w2 - 2 * zz * w1 - W / 4) / rt))))
    l0, l1, l2 = [num.coeff(h) for h in (h0, h1, h2)]
    e2, e1, e0 = sp.expand(sp.cancel(l2 / t)), sp.expand(l1), sp.expand(sp.cancel(l0 / t))
    b2 = sp.expand((4 * x * e2).subs(t, sp.sqrt(x)))
    b1 = sp.expand((2 * e2 + 2 * e1).subs(t, sp.sqrt(x)))
    b0 = sp.expand(e0.subs(t, sp.sqrt(x)))
    L = -sp.Rational(1, 2) / x                                       # y = x^(1/2) H
    tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
    den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
    return [sp.expand(sp.cancel(v * den)) for v in tw]


def _map():
    m = DATA["map"]
    return _poly(m["N0"]), _poly(m["N1"]), sp.expand(sp.sympify(m["E"]))


def _exact(target, source, N0, N1, E):
    U, V = _residual(target, source, N0, N1, E)
    return U == 0 and V == 0


def test_order_two_piece_is_the_image_of_the_legendre_operator():
    target = _dform(DATA["lin_theta"])
    N0, N1, E = _map()
    assert N1 != 0
    assert _exact(target, _source(), N0, N1, E)
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    assert not _exact(target, _source(), N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
    assert not _exact(target, _source(z_sign=+1), N0, N1, E)
    assert not _exact(target, _source(gamma=sp.Rational(3, 2)), N0, N1, E)


def test_order_two_piece_differs_from_the_exponential_shapes():
    lin = [_poly(v) for v in DATA["lin_theta"]]
    exp_ = [_poly(v) for v in DATA["exp_C2b_theta"]]
    assert sp.degree(lin[2], x) == 36 and sp.degree(exp_[2], x) == 30
    ratios = {sp.cancel(a / b) for a, b in zip(lin, exp_)}
    assert not (len(ratios) == 1 and next(iter(ratios)).is_number)


def _theta_apply(theta_coeffs, f):
    cs = [_poly(v, xp) for v in theta_coeffs]
    out, g = 0, f
    for c in cs:
        out += c * g
        g = 2 * xp * sp.diff(g, xp)
    return sp.simplify(sp.powsimp(sp.together(out / f), force=True))


def _solution(k):
    return sp.sympify(DATA["order_one_solutions"][k].replace("^", "**"), locals={"x": xp})


def test_order_one_pieces_are_algebraic_with_the_stated_exponents():
    F = DATA["order_one_factors"]
    for k in ("factor_1", "factor_2", "factor_3", "factor_4"):
        assert _theta_apply(F[k], _solution(k)) == 0, k
    # a solution offered to the wrong factor is not killed
    assert _theta_apply(F["factor_3"], _solution(k="factor_4")) != 0
    # exponents: -5/2 at x = 4 (t = 2) for factor_2; -11/2, -13/2 at x = 16 (t = 4) for factors 3 and 4
    for k, pt, e in (("factor_2", 4, sp.Rational(-5, 2)), ("factor_3", 16, sp.Rational(-11, 2)), ("factor_4", 16, sp.Rational(-13, 2))):
        lead = sp.limit(_solution(k) * (xp - pt) ** (-e), xp, pt)
        assert lead.is_finite and lead != 0, (k, lead)
        assert sp.limit(_solution(k) * (xp - pt) ** (-e + 1), xp, pt) == 0, k
