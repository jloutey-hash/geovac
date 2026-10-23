"""Paper 59 sec:density_ring: the cusp-correct geminal r e^{-gamma r} brings no new type at gamma = 3/2.

Record docs/predictions/2026-10-22_geminal_shapes.md (GA1, committed 5e40da9a before computing).  Written as a
separate step from the result; no debug/ import.  The order-two composition factor of Y_lin's annihilator (theta
form in x = t^2) and the map A = (N0 + N1 d/dx)/E are loaded from tests/data/paper59_ga1_2026_10_22.json.  The
source is rebuilt here: the Legendre operator of degree -1/2 in z = (t^2 + 1 - gamma^2)/(2t) for H = t^(-1/2) w(z),
moved to x and twisted by x^(1/2), as in tests/test_paper59_c2b_feynman_g2.py.

Guarded: (a) the factor is the image of that Legendre operator under the stored map, exactly, and a 1/1000
change in the map breaks it; (b) the factor is a different operator from the exponential case's C2b.
Wrong answers rejected: the gamma = 2 Legendre operator; the real-momentum argument; the exponential case's
factor C2b in place of this one.
Not re-run here: the order-nine annihilator, its composition series, the K3 map (record only).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import sympy as sp

from test_paper59_c2b_elliptic import _residual
from test_paper59_c2b_feynman_g2 import _source as _source_g2
import test_paper59_c2b_feynman_g2 as G2

x = sp.Symbol("x")
HERE = Path(__file__).parent
DATA = json.loads((HERE / "data" / "paper59_ga1_2026_10_22.json").read_text())
C2B = json.loads((HERE / "data" / "paper59_c2b_2026_10_22.json").read_text())["C2b_theta"]
GAMMA = sp.Rational(3, 2)
Z_SIGN = -1
FACTOR_KEY = "factor_theta"      # "c2b" plants the exponential case's factor


def _poly(s):
    return sp.expand(sp.sympify(s.replace("^", "**"), locals={"x": x}))


def _target():
    src = DATA[FACTOR_KEY] if FACTOR_KEY in DATA else C2B
    c0, c1, c2 = [_poly(v) for v in src]
    return [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]


def _source():
    G2.GAMMA, G2.Z_SIGN = GAMMA, Z_SIGN
    try:
        return _source_g2()
    finally:
        G2.GAMMA, G2.Z_SIGN = sp.Integer(2), -1


@pytest.mark.slow
def test_factor_is_the_image_of_the_legendre_operator():
    N0, N1, E = _poly(DATA["N0"]), _poly(DATA["N1"]), sp.expand(sp.sympify(DATA["E"]))
    assert N1 != 0
    U, V = _residual(_target(), _source(), N0, N1, E)
    assert U == 0 and V == 0
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    Up, Vp = _residual(_target(), _source(), N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
    assert not (Up == 0 and Vp == 0)


def test_factor_differs_from_the_exponential_case():
    f = [_poly(v) for v in DATA["factor_theta"]]
    c = [_poly(v) for v in C2B]
    assert sp.degree(f[2], x) != sp.degree(c[2], x) or len({sp.cancel(a / b) for a, b in zip(f, c)}) > 1
