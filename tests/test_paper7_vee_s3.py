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
(no removed production code).  The integrands are NOT typed in: the
projected densities Phi_1s, Phi_2s are DERIVED here -- the radial Fourier
transform of the hydrogenic 1s and 2s densities gives eq:rho1s / eq:rho2s, and
q = 2Zt (eq:projected_density) gives eq:phi1s / eq:phi2s -- and the three
integrals are their products.  Checks BOTH the bare integral values and the
assembled Slater integrals against eq:slater_exact. Provenance: SYMBOLIC PROOF
(exact rational/pi-rational, no numerics).
"""
from __future__ import annotations

import sympy as sp

t, Z = sp.symbols("t Z", positive=True)
r, q = sp.symbols("r q", positive=True)


def _integral(integrand):
    return sp.simplify(sp.integrate(integrand, (t, 0, sp.oo)))


def _form_factor(density):
    """rho~(q) = int 4 pi r^2 rho(r) sin(qr)/(qr) dr (spherical density)."""
    return sp.simplify(sp.integrate(
        4 * sp.pi * r**2 * density * sp.sin(q * r) / (q * r), (r, 0, sp.oo)))


RHO_1S = _form_factor(Z**3 / sp.pi * sp.exp(-2 * Z * r))
RHO_2S = _form_factor(Z**3 / (32 * sp.pi) * (2 - Z * r)**2 * sp.exp(-Z * r))
PHI_1S = sp.simplify(RHO_1S.subs(q, 2 * Z * t))
PHI_2S = sp.simplify(RHO_2S.subs(q, 2 * Z * t))


def test_projected_densities_match_printed_forms():
    """eq:rho1s, eq:rho2s, eq:phi1s, eq:phi2s as printed in Paper 7."""
    assert sp.simplify(RHO_1S - 16 * Z**4 / (q**2 + 4 * Z**2)**2) == 0, RHO_1S
    assert sp.simplify(RHO_2S - Z**4 * (Z**4 - 3 * Z**2 * q**2 + 2 * q**4)
                       / (Z**2 + q**2)**4) == 0, RHO_2S
    assert sp.simplify(PHI_1S - 1 / (t**2 + 1)**2) == 0, PHI_1S
    assert sp.simplify(PHI_2S - (1 - 12 * t**2 + 32 * t**4)
                       / (4 * t**2 + 1)**4) == 0, PHI_2S


def test_f0_1s1s_exact():
    """F^0(1s,1s): integral = 5*pi/32, assembled Slater = 5Z/8."""
    I = _integral(PHI_1S * PHI_1S)
    assert sp.simplify(I - sp.Rational(5, 32) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(5, 8) * Z) == 0, F0


def test_f0_1s2s_exact():
    """F^0(1s,2s): integral = 17*pi/324, assembled Slater = 17Z/81."""
    I = _integral(PHI_1S * PHI_2S)
    assert sp.simplify(I - sp.Rational(17, 324) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(17, 81) * Z) == 0, F0


def test_f0_2s2s_exact():
    """F^0(2s,2s): integral = 77*pi/2048, assembled Slater = 77Z/512."""
    I = _integral(PHI_2S * PHI_2S)
    assert sp.simplify(I - sp.Rational(77, 2048) * sp.pi) == 0, I
    F0 = sp.simplify(sp.Rational(4) * Z / sp.pi * I)
    assert sp.simplify(F0 - sp.Rational(77, 512) * Z) == 0, F0
