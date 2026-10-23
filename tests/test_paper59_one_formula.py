"""Paper 59 sec:density_ring: one function, L(a, b; c) = int_0^oo K0(ax) K0(bx) cos(cx) dx
= pi^2 P_(-1/2)(z) / (4 sqrt(ab)), z = (a^2 + b^2 + c^2)/(2ab), carries both elliptic pieces.

Written as a separate step from the paper paragraph; no debug/ import.  Claims guarded:
  (a) at c = 0 it is the three-centre radial core: with a_i = c_i^(-1/2),
      (2/(pi sqrt(c1 c2))) L(a1, a2; 0) = K(m)/sqrt(cmax), m = 1 - cmin/cmax (Paper 59 eq:period);
  (b) at c = i gamma, z = cos(theta) of the triangle with sides t, 1, gamma (theta opposite gamma), the modulus
      lambda = (1 - z)/2 = sin^2(theta/2), and lambda (1 - lambda) = (area / t)^2 = -Kallen/(16 t^2);
  (c) the right angle gamma^2 = t^2 + 1 gives lambda = 1/2 and j = 1728, for every gamma; at gamma = 3/2 that is
      x = 5/4, where J1 = P^3/(x (4x-1)^2 (4x-25)^2) = 1728.
Wrong answers rejected: a real-momentum sign in z; a wrong prefactor; a wrong Heron constant; the wrong right
angle (gamma^2 = t^2 - 1).
"""
from __future__ import annotations

import mpmath as mp
import sympy as sp

Z_SIGN = 1            # z = (a^2 + b^2 + Z_SIGN c^2)/(2ab) with c real
PREF = 4              # L = pi^2 P / (PREF sqrt(ab))
HERON = 16            # HERON * area^2 = -Kallen
RIGHT = 1             # right angle: gamma^2 = t^2 + RIGHT


def test_zero_momentum_is_the_three_centre_fibre():
    mp.mp.dps = 25
    for c1, c2 in ((mp.mpf("0.21"), mp.mpf("0.09")), (mp.mpf("0.05"), mp.mpf("0.2"))):
        a, b = 1 / mp.sqrt(c1), 1 / mp.sqrt(c2)
        moment = mp.quad(lambda y: mp.besselk(0, a * y) * mp.besselk(0, b * y), [0, 0.5, 1, 2, 4, 8, mp.inf])
        L = mp.pi**2 / (PREF * mp.sqrt(a * b)) * mp.legenp(-0.5, 0, (a * a + b * b) / (2 * a * b), type=3)
        assert abs(moment - L) / L < mp.mpf(10) ** -20
        cmax, cmin = max(c1, c2), min(c1, c2)
        period = mp.ellipk(1 - cmin / cmax) / mp.sqrt(cmax)
        assert abs(2 / (mp.pi * mp.sqrt(c1 * c2)) * L - period) / period < mp.mpf(10) ** -20


def test_real_momentum_cosine_form():
    mp.mp.dps = 20
    a, b, c = mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf(2)
    lhs = mp.quad(lambda y: mp.besselk(0, a * y) * mp.besselk(0, b * y) * mp.cos(c * y), [0, 0.5, 1, 2, 4, 8, 16, mp.inf])
    z = (a * a + b * b + Z_SIGN * c * c) / (2 * a * b)
    rhs = mp.pi**2 / (PREF * mp.sqrt(a * b)) * mp.legenp(-0.5, 0, z, type=3)
    assert abs(lhs - rhs) / abs(lhs) < mp.mpf(10) ** -15


def test_imaginary_momentum_is_a_triangle():
    t, g, th = sp.symbols("t gamma theta", positive=True)
    z = (t**2 + 1 + Z_SIGN * (sp.I * g) ** 2) / (2 * t)
    # law of cosines for the side gamma opposite theta
    assert sp.simplify(z.subs(g, sp.sqrt(t**2 + 1 - 2 * t * sp.cos(th))) - sp.cos(th)) == 0
    lam = (1 - z) / 2
    assert sp.simplify(lam.subs(g, sp.sqrt(t**2 + 1 - 2 * t * sp.cos(th))) - sp.sin(th / 2) ** 2) == 0
    kallen = (t**2 - (1 + g) ** 2) * (t**2 - (1 - g) ** 2)
    s = (t + 1 + g) / 2
    area2 = s * (s - t) * (s - 1) * (s - g)                 # Heron
    assert sp.expand(HERON * area2 + kallen) == 0
    assert sp.cancel(lam * (1 - lam) - area2 / t**2) == 0


def test_right_angle_is_the_cm_fibre():
    t, g = sp.symbols("t gamma", positive=True)
    z = (t**2 + 1 - g**2) / (2 * t)
    zr = sp.simplify(z.subs(g, sp.sqrt(t**2 + RIGHT)))
    lam = sp.simplify((1 - zr) / 2)
    j = sp.simplify(64 * (zr**2 + 3) ** 3 / (1 - zr**2) ** 2)
    assert lam == sp.Rational(1, 2) and j == 1728
    x = sp.Rational(9, 4) - RIGHT                          # t^2 at gamma = 3/2
    P = 16 * x**2 + 152 * x + 25
    assert P**3 / (x * (4 * x - 1) ** 2 * (4 * x - 25) ** 2) == 1728
