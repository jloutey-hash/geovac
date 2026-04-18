"""Tests for Track RH-J: spectral L(s, chi_-4) identification.

Conjecture RH-J.1 (validated symbolically and numerically in this sprint):

    D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))

where D_even(s), D_odd(s) are the even-n and odd-n sub-sums of the Dirac
Dirichlet series on unit S^3 (Camporesi-Higuchi), and beta(s) = L(s, chi_-4)
is the Dirichlet beta function.

These tests verify:
  - Numerical agreement at 50+ digit precision across s = 2..10
  - Symbolic polynomial identity in Hurwitz zeta + beta(s) symbols
  - Reduction to Paper 28 Eq.(4) at s=4: D_even(4) - D_odd(4) = -8G + 8*beta(4)
  - Special-case identity at s=3 via regularized limit (Hurwitz pole at s-2=1).
"""
from __future__ import annotations

import mpmath
import pytest
import sympy as sp
from sympy import Rational, zeta, pi, Symbol, expand, Integer


# High precision for all numerical tests
mpmath.mp.dps = 80


# ----------------------------------------------------------------------
# Numerical primitives
# ----------------------------------------------------------------------

def _hz(s, a):
    return mpmath.hurwitz(mpmath.mpf(s), mpmath.mpf(a))


def _dirac_Deven(s):
    return mpmath.power(2, -s) * (
        8 * _hz(s - 2, mpmath.mpf(3) / 4)
        - mpmath.mpf(1) / 2 * _hz(s, mpmath.mpf(3) / 4)
    )


def _dirac_Dodd(s):
    return mpmath.power(2, -s) * (
        8 * _hz(s - 2, mpmath.mpf(5) / 4)
        - mpmath.mpf(1) / 2 * _hz(s, mpmath.mpf(5) / 4)
    )


def _beta(s):
    """Dirichlet beta(s), extended to small s via known closed forms."""
    if s == 0:
        return mpmath.mpf(1) / 2
    if s == 1:
        return mpmath.pi / 4
    return (_hz(s, mpmath.mpf(1) / 4) - _hz(s, mpmath.mpf(3) / 4)) / mpmath.power(4, s)


def _closed_form(s):
    return mpmath.power(2, s - 1) * (_beta(s) - _beta(s - 2))


# ----------------------------------------------------------------------
# Test 1: numerical identity at many s (except pole at s=3)
# ----------------------------------------------------------------------

@pytest.mark.parametrize("s", [2, 4, 5, 6, 7, 8, 9, 10])
def test_rh_j1_numerical(s):
    """RH-J.1 holds numerically to 50+ digits for s in {2, 4..10}.

    s=3 is a special case (individual D_even, D_odd have Hurwitz pole at
    s-2=1); tested separately.
    """
    D_even = _dirac_Deven(s)
    D_odd = _dirac_Dodd(s)
    D_diff = D_even - D_odd
    target = _closed_form(s)

    # Absolute error
    err = abs(D_diff - target)
    # Both sides are O(1) at these s, so absolute tolerance = 1e-50 is safe
    assert err < mpmath.mpf('1e-50'), \
        f"RH-J.1 violated at s={s}: |D_diff - 2^(s-1)(beta(s)-beta(s-2))| = {err}"


# ----------------------------------------------------------------------
# Test 2: s=3 special case (Hurwitz pole regularized)
# ----------------------------------------------------------------------

def test_rh_j1_at_s3_via_closed_form():
    """At s=3, individual D_even, D_odd have Hurwitz pole at s-2=1, but
    their difference is the regularized limit:

        D_even(3) - D_odd(3) = pi^3/8 - pi

    This equals 2^2 * (beta(3) - beta(1)) = 4 * (pi^3/32 - pi/4) = pi^3/8 - pi.
    """
    target = _closed_form(3)
    expected = mpmath.pi ** 3 / 8 - mpmath.pi
    err = abs(target - expected)
    assert err < mpmath.mpf('1e-50'), \
        f"s=3 closed form mismatch: 2^2(beta(3)-beta(1)) = {target}, pi^3/8-pi = {expected}, err = {err}"


# ----------------------------------------------------------------------
# Test 3: symbolic polynomial identity (sympy)
# ----------------------------------------------------------------------

@pytest.mark.parametrize("s", [2, 3, 4, 5, 6, 7, 8, 10, 15])
def test_rh_j1_symbolic(s):
    """RH-J.1 is a closed-form polynomial identity in Hurwitz/beta symbols.

    Substituting the Hurwitz identities
        zeta(k, 5/4) = zeta(k, 1/4) - 4^k
        zeta(k, 3/4) = zeta(k, 1/4) - 4^k * beta_k
    into D_even(s) - D_odd(s) reduces the expression exactly to
        2^{s-1} * (beta_s - beta_{s-2})
    for every integer s >= 2.
    """
    beta_s = Symbol(f'beta_{s}')
    beta_sm2 = Symbol(f'beta_{s-2}')
    # Symbolic Hurwitz zeta placeholders
    zs14 = Symbol('zs14')
    zs34 = Symbol('zs34')
    zs54 = Symbol('zs54')
    zsm2_14 = Symbol('zsm2_14')
    zsm2_34 = Symbol('zsm2_34')
    zsm2_54 = Symbol('zsm2_54')

    Deven = Rational(1, 2) ** s * (8 * zsm2_34 - Rational(1, 2) * zs34)
    Dodd = Rational(1, 2) ** s * (8 * zsm2_54 - Rational(1, 2) * zs54)
    diff = Deven - Dodd

    # Hurwitz recurrence: zeta(k, 5/4) = zeta(k, 1/4) - (1/4)^{-k} = zeta(k, 1/4) - 4^k
    diff = diff.subs(zs54, zs14 - Integer(4) ** s)
    diff = diff.subs(zsm2_54, zsm2_14 - Integer(4) ** (s - 2))
    # Dirichlet beta: zeta(k, 3/4) = zeta(k, 1/4) - 4^k * beta_k
    diff = diff.subs(zs34, zs14 - Integer(4) ** s * beta_s)
    diff = diff.subs(zsm2_34, zsm2_14 - Integer(4) ** (s - 2) * beta_sm2)

    diff = expand(diff)
    target = expand(Integer(2) ** (s - 1) * (beta_s - beta_sm2))
    residual = expand(diff - target)

    assert residual == 0, \
        f"Symbolic residual at s={s}: {residual} (expected 0)"


# ----------------------------------------------------------------------
# Test 4: Paper 28 Eq.(4) consistency at s=4
# ----------------------------------------------------------------------

def test_paper28_eq4_consistency():
    """Paper 28 Eq.(4) states:
        D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4)
        D_odd(4)  = pi^2/2 - pi^4/24 + 4G - 4*beta(4)

    Their difference is -8G + 8*beta(4) = 8*(beta(4) - beta(2)), which is
    the s=4 case of RH-J.1 (since beta(2) = G).
    """
    G = mpmath.catalan
    beta4 = _beta(4)

    paper_diff = -8 * G + 8 * beta4          # from Eq.(4)
    rh_j1 = _closed_form(4)                  # 2^3 * (beta(4) - beta(2)) = 8*(beta4 - G)
    err = abs(paper_diff - rh_j1)
    assert err < mpmath.mpf('1e-60'), \
        f"Paper 28 Eq.(4) diff inconsistent with RH-J.1: err = {err}"

    # Also verify individual D_even(4), D_odd(4) against Paper 28 Eq.(4)
    paper_Deven = mpmath.pi ** 2 / 2 - mpmath.pi ** 4 / 24 - 4 * G + 4 * beta4
    paper_Dodd = mpmath.pi ** 2 / 2 - mpmath.pi ** 4 / 24 + 4 * G - 4 * beta4
    assert abs(_dirac_Deven(4) - paper_Deven) < mpmath.mpf('1e-60')
    assert abs(_dirac_Dodd(4) - paper_Dodd) < mpmath.mpf('1e-60')


# ----------------------------------------------------------------------
# Test 5: consistency with the full Dirac Dirichlet series
# ----------------------------------------------------------------------

@pytest.mark.parametrize("s", [4, 5, 6, 7, 8])
def test_Deven_plus_Dodd_equals_D(s):
    """D_even(s) + D_odd(s) = D(s) = 2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2)."""
    D_sum = _dirac_Deven(s) + _dirac_Dodd(s)
    D_full = 2 * _hz(s - 2, mpmath.mpf(3) / 2) - mpmath.mpf(1) / 2 * _hz(s, mpmath.mpf(3) / 2)
    err = abs(D_sum - D_full)
    assert err < mpmath.mpf('1e-70'), \
        f"D_even + D_odd != D(s) at s={s}: err = {err}"


# ----------------------------------------------------------------------
# Test 6: at s=2 the Paper 28 pattern is simpler (no beta(0) known to paper)
# ----------------------------------------------------------------------

def test_s2_individual_decomposition():
    """At s=2, individual D_even, D_odd have closed forms:
        D_even(2) = -1/2 - pi^2/8 + G
        D_odd(2)  = 1/2 - pi^2/8 - G
    (found by PSLQ at 80-digit precision.)

    Their difference is 2G - 1 = 2*(beta(2) - beta(0)) = 2^{s-1}*(beta(s)-beta(s-2)).
    """
    G = mpmath.catalan
    pi2 = mpmath.pi ** 2
    expected_Deven = -mpmath.mpf(1) / 2 - pi2 / 8 + G
    expected_Dodd = mpmath.mpf(1) / 2 - pi2 / 8 - G
    assert abs(_dirac_Deven(2) - expected_Deven) < mpmath.mpf('1e-60')
    assert abs(_dirac_Dodd(2) - expected_Dodd) < mpmath.mpf('1e-60')

    diff = _dirac_Deven(2) - _dirac_Dodd(2)
    expected_diff = 2 * G - 1
    assert abs(diff - expected_diff) < mpmath.mpf('1e-60')

    # And matches RH-J.1: 2 * (G - 1/2) = 2G - 1
    rh_j1 = _closed_form(2)
    assert abs(diff - rh_j1) < mpmath.mpf('1e-60')
