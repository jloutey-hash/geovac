"""Paper 7 sec:continuum (zero energy) -- the third shape of Fock's projection.

Record: docs/predictions/2026-10-18_zero_energy.md (ZE1-ZE3, committed before the computation in
bd2435d1). Written as its own step after the result; it does not import the computation script.
Each test names the wrong answer it rejects.

  ZE1  the zero-energy radial function sqrt(r) J_{2l+1}(sqrt(8Zr)) solves the E = 0 radial equation;
       on it the Runge-Lenz commutator vanishes and A^2 = Z^2 (slow: symbolic operators), with an
       E < 0 control where Pauli's -2i E L_z survives.
  ZE2  R_nl^2 n^3/Z^2 -> 2 J_{2l+1}(sqrt(8Zr))^2 / r with corrections in 1/n^2 only; the continuum
       density tends to the same function with corrections in k^2.
  ZE3  the threshold radial density is pi-free: its Taylor coefficients are rational, and the
       contact density per unit energy is Z/pi only through the solid angle 1/(4 pi).
"""
from __future__ import annotations

import random

import mpmath as mp
import pytest
import sympy as sp


def _target(l, Z, r):
    return 2 * mp.besselj(2 * l + 1, mp.sqrt(8 * Z * r)) ** 2 / r


def _bound_density(n, l, Z, r):
    x = 2 * mp.mpf(Z) * r / n
    lognorm = 3 * mp.log(2 * mp.mpf(Z) / n) + mp.loggamma(n - l) - mp.log(2 * n) - mp.loggamma(n + l + 1)
    R = mp.exp(lognorm / 2) * x**l * mp.e ** (-x / 2) * mp.laguerre(n - l - 1, 2 * l + 1, x)
    return R**2 * mp.mpf(n) ** 3 / Z**2


def test_ze1_zero_energy_function_solves_the_radial_equation():
    """Rejects the argument sqrt(4Zr) or the order 2l in place of sqrt(8Zr) and 2l+1."""
    r, Z = sp.symbols("r Z", positive=True)
    for l in range(3):
        u = sp.sqrt(r) * sp.besselj(2 * l + 1, sp.sqrt(8 * Z * r))
        ode = -sp.diff(u, r, 2) / 2 + (l * (l + 1) / (2 * r**2) - Z / r) * u
        for rv, Zv in ((sp.Rational(3, 7), 1), (sp.Rational(5, 2), 2)):
            assert abs(sp.N(ode.subs({r: rv, Z: Zv}), 30)) < 1e-25


def test_ze2_bound_side_limit_and_no_odd_powers():
    """Rejects a limit constant other than 2 (for example sqrt 2), a different Bessel argument, and a
    1/n correction: n^2 (ratio - 1) must settle as n doubles."""
    mp.mp.dps = 40
    for Z in (1, 2):
        for l in (0, 1):
            r = mp.mpf(2)
            T = _target(l, Z, r)
            s = [(_bound_density(n, l, Z, r) / T - 1) * n**2 for n in (200, 400, 800)]
            assert abs(s[2] - s[1]) < abs(s[1] - s[0]) and abs(s[2] - s[1]) < mp.mpf("0.01") * (1 + abs(s[2]))
            assert abs(s[2]) < 50


def test_ze2_continuum_side_meets_the_same_function():
    """Rejects a continuum limit other than the bound one, and a k^1 correction."""
    mp.mp.dps = 30
    Z, l, r = 1, 0, mp.mpf(2)
    T = _target(l, Z, r)
    s = []
    for k in (mp.mpf(1) / 100, mp.mpf(1) / 200, mp.mpf(1) / 400):
        dens = (2 / (mp.pi * k)) * (mp.coulombf(l, -Z / k, k * r) / r) ** 2
        s.append((dens / T - 1) / k**2)
    assert abs(s[2] - s[1]) < mp.mpf("0.01") * abs(s[2])
    # and the bound 1/n^2 coefficient is minus Z^2 times this k^2 coefficient (one function of E)
    mp.mp.dps = 40
    b = [(_bound_density(n, l, Z, r) / T - 1) * n**2 for n in (400, 800)]
    assert abs(b[1] + Z**2 * s[2]) < mp.mpf("0.01") * abs(s[2])


def test_ze3_threshold_density_is_pi_free():
    """Rejects any pi in the radial threshold density, and a contact density other than Z/pi."""
    x = sp.symbols("x", positive=True)
    for l in range(3):
        nu = 2 * l + 1
        ser = sp.series(sp.besselj(nu, 2 * sp.sqrt(x)) ** 2 / x**nu, x, 0, 8).removeO()
        assert all(c.is_rational for c in sp.Poly(sp.expand(ser), x).all_coeffs())
    r, Z = sp.symbols("r Z", positive=True)
    contact = sp.limit(2 * sp.besselj(1, sp.sqrt(8 * Z * r)) ** 2 / r, r, 0) / (4 * sp.pi)
    assert sp.simplify(contact - Z / sp.pi) == 0


# ---------------------------------------------------------------- ZE1, operator level (slow)
X = sp.symbols("x y z", real=True)
_rr = sp.sqrt(sum(c**2 for c in X))
_EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def _p(i, f):
    return -sp.I * sp.diff(f, X[i])


def _L(i, f):
    return sum(e * X[j] * _p(k, f) for (a, j, k), e in _EPS.items() if a == i)


def _A(i, f, Z):
    pxL = sum(e * _p(j, _L(k, f)) for (a, j, k), e in _EPS.items() if a == i)
    Lxp = sum(e * _L(j, _p(k, f)) for (a, j, k), e in _EPS.items() if a == i)
    return (pxL - Lxp) / 2 - Z * X[i] / _rr * f


@pytest.mark.slow
def test_ze1_runge_lenz_contracts_to_the_euclidean_algebra():
    """Rejects a nonvanishing [A_x, A_y] or A^2 != Z^2 on a zero-energy state, and an operator
    implementation that would also kill the commutator at E < 0 (the 2p control keeps -2i E L_z)."""
    Z = 1
    pt = dict(zip(X, (sp.Rational(7, 100), -sp.Rational(11, 100), sp.Rational(5, 100))))
    # zero-energy p state, m = 1: J_3(sqrt(8r)) / r^(3/2) (x + i y), as its series (16 terms)
    g = sum((-2 * Z * _rr) ** k / (sp.factorial(k) * sp.factorial(k + 3)) for k in range(16))
    psi = g * (X[0] + sp.I * X[1])
    Ax, Ay, Az = (_A(i, psi, Z) for i in range(3))
    scale = abs(sp.N(psi.subs(pt), 40))
    comm = _A(0, Ay, Z) - _A(1, Ax, Z)
    A2 = _A(0, Ax, Z) + _A(1, Ay, Z) + _A(2, Az, Z) - Z**2 * psi
    assert abs(sp.N(comm.subs(pt), 40)) / scale < 1e-28
    assert abs(sp.N(A2.subs(pt), 40)) / scale < 1e-28
    psi2 = sp.exp(-_rr / 2) * (X[0] + sp.I * X[1])  # 2p, E = -1/8
    Bx, By = _A(0, psi2, Z), _A(1, psi2, Z)
    c2 = sp.N((_A(0, By, Z) - _A(1, Bx, Z)).subs(pt), 40)
    pred = sp.N((-2 * sp.I * sp.Rational(-1, 8) * _L(2, psi2)).subs(pt), 40)
    assert abs(c2) > 1e-3 and abs(c2 - pred) < 1e-28 * abs(pred)
