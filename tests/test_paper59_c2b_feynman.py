"""Paper 59 sec:density_ring: C2b as a two-dimensional bubble at imaginary momentum.

Record docs/predictions/2026-10-22_c2b_feynman.md (CF1-CF4, committed d68e6fe5 before computing).
Written as a separate step from the result; no debug/ import.  Claims guarded:

  (a) twisted by t, the Legendre operator of degree -1/2 in z = (t^2 + 1 - g^2)/(2t), g = 3/2, is the EC
      period operator (the j-line 2F1(1/12, 5/12; 1; 1728/J) pulled back and gauged by x^(1/2) P^(-1/4),
      built in tests/test_paper59_c2b_elliptic.py) coefficient for coefficient, up to a scalar;
  (b) C(t) = int_0^oo K0(tx) K0(x) cosh(g x) dx solves the untwisted operator, and the registered candidate
      with exp(-g x) does not;
  (c) C(t) = pi^2 P_(-1/2)(z) / (4 sqrt t).

Wrong answers each part rejects: the real-momentum argument z = (t^2 + 1 + g^2)/(2t); a degree other than -1/2;
a missing twist; exp(-g x) for cosh(g x); a wrong prefactor.
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

from test_paper59_c2b_elliptic import _P, _source, x

t = sp.Symbol("t", positive=True)
G = sp.Rational(3, 2)
NU_NU1 = -sp.Rational(1, 4)          # nu (nu + 1) at nu = -1/2
Z_SIGN = -1                          # z = (t^2 + 1 + Z_SIGN g^2)/(2t): imaginary momentum
TWIST = True
KERNEL = "cosh"
PREF = 4                             # C = pi^2 P / (PREF sqrt t)


def _z():
    return (t**2 + 1 + Z_SIGN * G**2) / (2 * t)


def _legendre_t(power):
    """[c0, c1, c2] in d/dt for y = t^power w(z(t)), w a Legendre function: (1 - z^2) w'' - 2 z w' + NU_NU1 w = 0"""
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    p = sp.Rational(power)
    # w = t^(-p) y
    W = t ** (-p) * h0
    W1 = -p * t ** (-p - 1) * h0 + t ** (-p) * h1
    W2 = p * (p + 1) * t ** (-p - 2) * h0 - 2 * p * t ** (-p - 1) * h1 + t ** (-p) * h2
    z = _z()
    zp, zpp = sp.diff(z, t), sp.diff(z, t, 2)
    w1 = W1 / zp
    w2 = (W2 - zpp * w1) / zp**2
    num = sp.expand(sp.numer(sp.together(sp.expand(((1 - z**2) * w2 - 2 * z * w1 + NU_NU1 * W) * t**p))))
    return [num.coeff(h) for h in (h0, h1, h2)]


def _ec_source_in_t():
    """the EC period operator (in x = t^2) rewritten in d/dt: y_x = y_t/(2t), y_xx = (y_tt - y_t/t)/(4t^2)"""
    b0, b1, b2 = [sp.expand(v.subs(x, t**2)) for v in _source(_P())]
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    e = sp.expand(sp.numer(sp.together(b2 * (h2 - h1 / t) / (4 * t**2) + b1 * h1 / (2 * t) + b0 * h0)))
    return [e.coeff(h) for h in (h0, h1, h2)]


def test_twisted_legendre_is_the_ec_period_operator():
    leg = _legendre_t(sp.Rational(1, 2) if TWIST else 0)     # y = t^(1/2) w: the EC twist x^(1/2) times t^(-1/2)
    ec = _ec_source_in_t()
    # the same operator: every coefficient ratio is one common function of t
    ratios = {sp.cancel(a / b) for a, b in zip(leg, ec)}
    assert len(ratios) == 1


def _residual(tv, kernel):
    k = {"cosh": lambda y: mp.cosh(mp.mpf(3) / 2 * y), "exp": lambda y: mp.exp(-mp.mpf(3) / 2 * y)}[kernel]
    f0 = lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * k(y)
    f1 = lambda y: -y * mp.besselk(1, tv * y) * mp.besselk(0, y) * k(y)
    f2 = lambda y: y**2 * (mp.besselk(0, tv * y) + mp.besselk(1, tv * y) / (tv * y)) * mp.besselk(0, y) * k(y)
    cuts = [0, 0.5, 1, 2, 4, 8, 16, 32, mp.inf]
    vals = [mp.quad(f, cuts) for f in (f0, f1, f2)]
    cs = _legendre_t(-sp.Rational(1, 2))                       # H = t^(-1/2) w
    terms = [sp.lambdify(t, c, "mpmath")(tv) * v for c, v in zip(cs, vals)]
    return abs(sum(terms)) / max(abs(v) for v in terms)


@pytest.mark.slow
def test_cosh_integral_solves_it_and_the_exp_candidate_does_not():
    mp.mp.dps = 25
    for tv in (mp.mpf("0.8"), mp.mpf("3.2")):
        assert _residual(tv, KERNEL) < mp.mpf(10) ** -18
    assert _residual(mp.mpf("1.7"), "exp") > mp.mpf(10) ** -3


@pytest.mark.slow
def test_closed_form_is_the_legendre_function():
    mp.mp.dps = 25
    g = mp.mpf(3) / 2
    for tv in (mp.mpf("0.8"),):
        C = mp.quad(lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * mp.cosh(g * y), [0, 0.5, 1, 2, 4, 8, 16, 32, mp.inf])
        z = (tv**2 + 1 + Z_SIGN * g**2) / (2 * tv)
        rhs = mp.pi**2 / (PREF * mp.sqrt(tv)) * mp.re(mp.legenp(-0.5, 0, z, type=2 if abs(z) < 1 else 3))
        assert abs(C - rhs) / abs(C) < mp.mpf(10) ** -18
