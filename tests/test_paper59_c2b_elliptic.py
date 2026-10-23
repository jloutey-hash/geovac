"""Paper 59 sec:density_ring: the order-two piece C2b of the correlated triangle at gamma = 3/2 is elliptic.

Record docs/predictions/2026-10-22_c2b_identification.md (EC1-EC3, committed 751a89fa before computing).
Written as a separate step from the result; no debug/ import.  C2b (theta = t d/dt coefficients in
x = t^2) and the map A = (N0 + N1 d/dx)/E found by the search are loaded from
tests/data/paper59_c2b_2026_10_22.json.  Everything else is rebuilt here from P and R:

  P = 16x^2 + 152x + 25,  R = (4x - 5)(16x^2 - 616x + 25),  P^3 - R^2 = 1728 x (4x-1)^2 (4x-25)^2,
  J = P^3 / (x (4x-1)^2 (4x-25)^2),  E_x : Y^2 = X^3 - 3x^2 P X - 2x^3 R  (j = J),
  source = x^(1/2) P^(-1/4) 2F1(1/12, 5/12; 1; 1728/J), pulled back and gauged to polynomial form.

Not re-run here: the numerical monodromy (integer traces at 80 digits) and the exhaustive degree-six
J search; both are in the record.  The exact map implies the monodromy statement.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

x = sp.Symbol("x")
DATA = json.loads((Path(__file__).parent / "data" / "paper59_c2b_2026_10_22.json").read_text())
D = x * (4 * x - 1) ** 2 * (4 * x - 25) ** 2


def _poly(s):
    return sp.expand(sp.sympify(s.replace("^", "**"), locals={"x": x}))


def _P():
    return _poly(DATA["P"])


def _R():
    return _poly(DATA["R"])


def _c2b_dx():
    """C2b in d/dx form: theta = 2x d/dx, theta^2 = 4x^2 d^2 + 4x d"""
    c0, c1, c2 = [_poly(s) for s in DATA["C2b_theta"]]
    return [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]


def _source(P, half_x=True):
    """2F1(1/12, 5/12; 1; z), z = 1728/J pulled back, then y = phi w with phi = x^(1/2) P^(-1/4)"""
    a, b, c = sp.Rational(1, 12), sp.Rational(5, 12), sp.Integer(1)
    J = P**3 / D
    z = sp.cancel(1728 / J)
    zp = sp.cancel(sp.diff(z, x))
    p = sp.cancel((c - (a + b + 1) * z) / (z * (1 - z)) * zp - sp.diff(zp, x) / zp)
    q = sp.cancel(-a * b / (z * (1 - z)) * zp**2)
    L = -sp.Rational(1, 4) * sp.diff(P, x) / P + (sp.Rational(1, 2) / x if half_x else 0)
    p2 = sp.cancel(p - 2 * L)
    q2 = sp.cancel(q - L * p + L**2 - sp.diff(L, x))
    den = sp.lcm(sp.fraction(p2)[1], sp.fraction(q2)[1])
    return [sp.cancel(q2 * den), sp.cancel(p2 * den), sp.expand(den)]


def _residual(target, source, N0, N1, E):
    """target(A y) for y in Sol(source), A = (N0 + N1 d)/E, as the pair (coefficient of y, of y')"""
    b0, b1, b2 = source

    def d(f, g):
        return (sp.cancel(sp.diff(f, x) - g * b0 / b2), sp.cancel(f + sp.diff(g, x) - g * b1 / b2))

    z0 = (sp.cancel(N0 / E), sp.cancel(N1 / E))
    z1 = d(*z0)
    z2 = d(*z1)
    a0, a1, a2 = target
    return sp.cancel(a2 * z2[0] + a1 * z1[0] + a0 * z0[0]), sp.cancel(a2 * z2[1] + a1 * z1[1] + a0 * z0[1])


def _map():
    return _poly(DATA["N0"]), _poly(DATA["N1"]), _poly(DATA["E"])


def test_family_and_j_invariant():
    P, R = _P(), _R()
    assert sp.expand(P**3 - R**2 - 1728 * D) == 0
    A, B = -3 * x**2 * P, -2 * x**3 * R
    j = sp.cancel(1728 * 4 * A**3 / (4 * A**3 + 27 * B**2))
    assert sp.cancel(j - P**3 / D) == 0
    disc = sp.factor(-16 * (4 * A**3 + 27 * B**2))
    assert sp.cancel(disc / (x**7 * (4 * x - 1) ** 2 * (4 * x - 25) ** 2)).is_number
    # orders at infinity on a rational elliptic surface (weights 4, 6, 12): fibre I1
    assert (4 - sp.degree(A, x), 6 - sp.degree(B, x), 12 - sp.degree(sp.expand(4 * A**3 + 27 * B**2), x)) == (0, 0, 1)


@pytest.mark.slow
def test_c2b_is_the_image_of_the_elliptic_period_operator():
    target, source = _c2b_dx(), _source(_P())
    N0, N1, E = _map()
    assert sp.degree(N1, x) >= 0 and N1 != 0          # a genuine order-one map
    U, V = _residual(target, source, N0, N1, E)
    assert U == 0 and V == 0
    # sharp: a 1/1000 change in the largest coefficient of N0 breaks it
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    N0p = N0 + sp.Rational(1, 1000) * cs[k] * x**k
    Up, Vp = _residual(target, source, N0p, N1, E)
    assert not (Up == 0 and Vp == 0)
    # the twist matters: without the factor x^(1/2) the same map fails
    Uw, Vw = _residual(target, _source(_P(), half_x=False), N0, N1, E)
    assert not (Uw == 0 and Vw == 0)


def test_source_solution_is_a_period_of_the_family():
    mp.mp.dps = 40
    P, R = _P(), _R()
    fP, fR = sp.lambdify(x, P, "mpmath"), sp.lambdify(x, R, "mpmath")
    vals = []
    for xv in (mp.mpf("0.004"), mp.mpf("0.01"), mp.mpf("0.03")):
        a, b = -3 * xv**2 * fP(xv), -2 * xv**3 * fR(xv)
        e1, e2, e3 = sorted([mp.re(r) for r in mp.polyroots([1, 0, a, b], extraprec=200)], reverse=True)
        w2 = mp.pi / mp.agm(mp.sqrt(e1 - e3), mp.sqrt(e2 - e3))
        J = fP(xv) ** 3 / (xv * (4 * xv - 1) ** 2 * (4 * xv - 25) ** 2)
        src = mp.sqrt(xv) * fP(xv) ** mp.mpf(-0.25) * mp.hyp2f1(mp.mpf(1) / 12, mp.mpf(5) / 12, 1, 1728 / J)
        vals.append(xv * w2 / src)
    for v in vals:
        assert abs(v - mp.pi / mp.sqrt(3)) < mp.mpf(10) ** -30
