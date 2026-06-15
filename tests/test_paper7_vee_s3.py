"""
Paper 7 §V — eq:f0_s3 master density-overlap formula backing test.

Paper 7 §V presents a closed-form S^3 density-overlap formula

    F^0(a,b) = (4Z/pi) * \\int_0^\\infty Phi_a(t) Phi_b(t) dt        (eq:f0_s3)

and claims it reproduces the exact hydrogenic Slater integrals
(eq:slater_exact):

    F^0(1s,1s) = 5Z/8,  F^0(1s,2s) = 17Z/81,  F^0(2s,2s) = 77Z/512,

via the three definite integrals 5*pi/32, 17*pi/324, 77*pi/2048.

The original direct test (test_vee_s3.py) was archived in v2.7.0 when
compute_vee_s3_overlap / _phi_s_orbital were removed from the codebase,
leaving this load-bearing §V claim with NO live test (coverage gap caught in
the 2026-06-14 trunk QA re-check). This test reinstates it with pure sympy
(no removed production code): it integrates the three closed-form integrands
symbolically and checks BOTH the bare integral values and the assembled
Slater integrals against eq:slater_exact. Provenance: SYMBOLIC PROOF (exact
rational/pi-rational, no numerics).
"""
from __future__ import annotations

import sympy as sp

t, Z = sp.symbols("t Z", positive=True)


def _integral(integrand):
    return sp.simplify(sp.integrate(integrand, (t, 0, sp.oo)))


def test_f0_1s1s_exact():
    """F^0(1s,1s): integral = 5*pi/32, assembled Slater = 5Z/8."""
    I = _integral(1 / (t**2 + 1) ** 4)
    assert sp.simplify(I - sp.Rational(5, 32) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(5, 8) * Z) == 0, F0


def test_f0_1s2s_exact():
    """F^0(1s,2s): integral = 17*pi/324, assembled Slater = 17Z/81."""
    num = 1 - 12 * t**2 + 32 * t**4
    den = (t**2 + 1) ** 2 * (4 * t**2 + 1) ** 4
    I = _integral(num / den)
    assert sp.simplify(I - sp.Rational(17, 324) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(17, 81) * Z) == 0, F0


def test_f0_2s2s_exact():
    """F^0(2s,2s): integral = 77*pi/2048, assembled Slater = 77Z/512."""
    num = (1 - 12 * t**2 + 32 * t**4) ** 2
    den = (4 * t**2 + 1) ** 8
    I = _integral(num / den)
    assert sp.simplify(I - sp.Rational(77, 2048) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(77, 512) * Z) == 0, F0
