"""Paper 59 sec:density_ring: why the correlated triangle drops at gamma = 1 (equal masses, not colliding
thresholds).

Record docs/predictions/2026-10-23_gamma1_mechanism.md (MG1-MG3, committed 037c1ea2 before computing).
Written as a separate step from the result; no debug/ import.  Data tests/data/paper59_mg_2026_10_23.json.

Two guards.
(1) At gamma = 3, where the pseudo-threshold t = gamma - 1 falls on the density-pair threshold t = 2, the
    elliptic piece is still there: the stored order-two factor at gamma = 3 is the image of the Legendre operator
    with z = (t^2 - 8)/(2t) under the stored exact map; likewise at gamma = 1/2 with z = (t^2 + 3/4)/(2t).  The
    source is rebuilt here from the Legendre equation, moved to x = t^2 and twisted by x^(1/2).
    Wrong answers rejected: the gamma = 3 factor against the gamma = 1/2 source (and the reverse); the
    real-momentum sign; a map perturbed by 1/1000 in one coefficient.
(2) The kite's operator in z has order six at gamma = 1/2, 2 and 3 and order two at gamma = 1, and the rational
    roots of its leading coefficient are exactly minus the squared pair masses {4, 4 gamma^2, (1+gamma)^2,
    (1-gamma)^2} (plus apparent points that carry no exponent obstruction, none at these gamma), collapsing to
    {-4} at gamma = 1.  Wrong answers rejected: the pair-mass set of a different gamma.
Not re-run here: the order-nine annihilator of Y, its composition series, the K3 map and the quadrature on Y
(record only).
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from test_paper59_c2b_elliptic import _residual

x = sp.Symbol("x")
z = sp.Symbol("z")
u = sp.Symbol("u")
t = sp.Symbol("t", positive=True)
DATA = json.loads((Path(__file__).parent / "data" / "paper59_mg_2026_10_23.json").read_text())


def _poly(s, var=x):
    return sp.expand(sp.sympify(str(s).replace("^", "**"), locals={var.name: var}))


def _target(g):
    c0, c1, c2 = [_poly(v) for v in DATA["gammas"][g]["C2b_theta"]]
    return [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]          # theta = 2x d/dx


def _source(g, z_sign=-1):
    """Legendre operator of degree -1/2 in z = (t^2 + 1 + z_sign gamma^2)/(2t) for H = t^(-1/2) w(z),
    moved to x = t^2 and twisted by x^(1/2)."""
    gam = sp.Rational(g)
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    zz = (t**2 + 1 + z_sign * gam**2) / (2 * t)
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


def _map(g):
    m = DATA["gammas"][g]["map"]
    return _poly(m["N0"]), _poly(m["N1"]), sp.expand(sp.sympify(m["E"]))


def _exact(target, source, N0, N1, E):
    U, V = _residual(target, source, N0, N1, E)
    return U == 0 and V == 0


def test_elliptic_piece_present_at_gamma_3_and_at_gamma_one_half():
    for g in ("3", "1/2"):
        N0, N1, E = _map(g)
        assert N1 != 0
        assert _exact(_target(g), _source(g), N0, N1, E)
        cs = sp.Poly(N0, x).all_coeffs()[::-1]
        k = max(range(len(cs)), key=lambda i: abs(cs[i]))
        assert not _exact(_target(g), _source(g), N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
        assert not _exact(_target(g), _source(g, z_sign=+1), N0, N1, E)
    # the gamma = 3 factor is not the image of the gamma = 1/2 source under either map, and the reverse
    N0, N1, E = _map("3")
    assert not _exact(_target("3"), _source("1/2"), N0, N1, E)
    N0, N1, E = _map("1/2")
    assert not _exact(_target("1/2"), _source("3"), N0, N1, E)


def _kite_lead_roots(g):
    Q = [_poly(q, u) for q in DATA["gammas"][g]["kite_Q"]]
    r = max(sp.degree(q, u) for q in Q)
    assert r == DATA["gammas"][g]["kite_order"]
    lead = sum(sp.Poly(q, u).coeff_monomial(u**r) * z**i for i, q in enumerate(Q))   # coefficient of (d/dz)^r is z^r * lead
    roots = sp.roots(sp.Poly(sp.expand(lead), z), filter="Q")
    return r, {sp.Rational(k) for k in roots}


def _pair_masses_sq(g):
    gam = sp.Rational(g)
    return {-sp.Integer(4), -4 * gam**2, -((1 + gam) ** 2), -((1 - gam) ** 2)}


def test_kite_singular_points_are_minus_the_squared_pair_masses():
    for g in ("1/2", "2", "3"):
        r, roots = _kite_lead_roots(g)
        assert r == 6
        assert roots == _pair_masses_sq(g), (g, roots)
        assert roots != _pair_masses_sq("3/2")
    r, roots = _kite_lead_roots("1")
    assert r == 2
    assert _pair_masses_sq("1") == {-sp.Integer(4), sp.Integer(0)}
    # at gamma = 1 the nonzero genuine point is -4 alone; 20 and -12 are apparent (record)
    assert -sp.Integer(4) in roots
