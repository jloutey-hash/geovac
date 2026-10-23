"""Paper 59 sec:density_ring: at gamma = 2 the order-two piece of the correlated triangle is the
Legendre operator of the imaginary-momentum reading.

Record docs/predictions/2026-10-22_c2b_feynman_g2.md (CF5-CF7, committed db45a5bb before computing).
Written as a separate step from the result; no debug/ import.  The order-two composition factor (theta form in
x = t^2) and the map A = (N0 + N1 d/dx)/E are loaded from tests/data/paper59_cf_g2_2026_10_22.json.  The source
is rebuilt here: the Legendre operator of degree -1/2 in z = (t^2 + 1 + Z_SIGN g^2)/(2t), for
H = t^(-1/2) w(z), moved to x and twisted by x^(1/2).

Wrong answers rejected: the gamma = 3/2 operator in place of gamma = 2; the real-momentum argument
(Z_SIGN = +1); a map perturbed by 1/1000 in one coefficient.
Not re-run here: the order-nine annihilator of Y, its composition series and the K3 map (record only).
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from test_paper59_c2b_elliptic import _residual

x = sp.Symbol("x")
t = sp.Symbol("t", positive=True)
DATA = json.loads((Path(__file__).parent / "data" / "paper59_cf_g2_2026_10_22.json").read_text())
GAMMA = sp.Integer(2)
Z_SIGN = -1


def _poly(s):
    return sp.expand(sp.sympify(s.replace("^", "**"), locals={"x": x}))


def _target():
    c0, c1, c2 = [_poly(v) for v in DATA["C2b_theta"]]
    return [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]          # theta = 2x d/dx


def _source():
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    z = (t**2 + 1 + Z_SIGN * GAMMA**2) / (2 * t)
    rt = sp.sqrt(t)
    W, W1, W2 = rt * h0, h0 / (2 * rt) + rt * h1, -h0 / (4 * t * rt) + h1 / rt + rt * h2   # w = t^(1/2) H
    zp, zpp = sp.diff(z, t), sp.diff(z, t, 2)
    w1 = W1 / zp
    w2 = (W2 - zpp * w1) / zp**2
    num = sp.expand(sp.numer(sp.together(sp.expand(((1 - z**2) * w2 - 2 * z * w1 - W / 4) / rt))))
    l0, l1, l2 = [num.coeff(h) for h in (h0, h1, h2)]
    # in x = t^2: H_t = 2t H_x, H_tt = 4x H_xx + 2 H_x; divide by t
    e2, e1, e0 = sp.expand(sp.cancel(l2 / t)), sp.expand(l1), sp.expand(sp.cancel(l0 / t))
    b2 = sp.expand((4 * x * e2).subs(t, sp.sqrt(x)))
    b1 = sp.expand((2 * e2 + 2 * e1).subs(t, sp.sqrt(x)))
    b0 = sp.expand(e0.subs(t, sp.sqrt(x)))
    # y = x^(1/2) H: H = phi y with phi'/phi = L = -1/(2x)
    L = -sp.Rational(1, 2) / x
    tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
    den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
    return [sp.expand(sp.cancel(v * den)) for v in tw]


def _map():
    return _poly(DATA["N0"]), _poly(DATA["N1"]), sp.expand(sp.sympify(DATA["E"]))


def test_order_two_piece_is_the_image_of_the_legendre_operator():
    N0, N1, E = _map()
    assert N1 != 0
    U, V = _residual(_target(), _source(), N0, N1, E)
    assert U == 0 and V == 0
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    Up, Vp = _residual(_target(), _source(), N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
    assert not (Up == 0 and Vp == 0)


def test_j_invariant_at_gamma_2():
    g2 = GAMMA**2
    z2 = (x + 1 + Z_SIGN * g2) ** 2 / (4 * x)
    j = sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2)
    assert sp.cancel(j - 16 * (x + 3) ** 6 / (x * (x - 1) ** 2 * (x - 9) ** 2)) == 0
