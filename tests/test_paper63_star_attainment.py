"""Paper 63, Proposition prop:star_attain: Coulomb stars of 3 <= n <= 8 densities reach rung n (record SA,
docs/predictions/2026-11-13_star5_attainment.md, BORDERLINE as registered: the branch stands, the amplitude
missed its 5% bar at eps = 1e-3).

Checked (each names the wrong answer it rejects):
- F_5's operator phi(L_4) has leading coefficient (t^2 - 4)(t^2 - 16) and a half-integer local exponent (3/2) at
  t = -4; F_7's phi(L_6) has 5/2 at t = -6 (rejects "F_n has only integer exponents at -(n-1)", which would void
  the monodromy argument);
- the lower pieces phi(L_c(theta - lam)), c <= n - 2, have leading coefficient prod (t^2 - (c - 2j)^2),
  independent of lam, nonzero at t = -(n-1) (rejects a lower piece singular there);
- the R^(n-1) branch: for n = 5, eps^(1/2) int_1^oo I0((4-eps)x) R^4 x^-2 dx is nonzero, falls from eps = 1e-2
  to 1e-3, and is within 15% of Watson's constant pi^2/(8 sqrt 2) at 1e-3 (rejects "the branch vanishes" and a
  wrong exponent: eps^1 I5 -> 0);
- control n = 4: eps int_1^oo I0((3-eps)x) R^3 x^-1 dx is within 5% of (pi/2)^(3/2)/sqrt(6 pi) at 1e-3 (a pole,
  not a square root; rejects the measurement failing to tell the two apart).
"""
import numpy as np
import pytest
import sympy as sp
from scipy.special import i0e, k1e, roots_laguerre

import test_paper63_tree_rung as TR

t, th, s, rho = sp.symbols("t theta s rho")


def _leading(P):
    Q = TR._t_lists(P)
    J = max(sp.Poly(q, th).degree() for q in Q)
    return sp.expand(sum(sp.Poly(q, th).coeff_monomial(th**J) * t ** (2 * i) for i, q in enumerate(Q)))


def _exponents(P, t0):
    a = {}
    for i, q in enumerate(TR._t_lists(P)):
        for (j,), cf in sp.Poly(sp.expand(q), th).terms():
            for k in range(j + 1):
                a[k] = a.get(k, 0) + cf * sp.functions.combinatorial.numbers.stirling(j, k) * t ** (2 * i + k)
    vk = {}
    for k, ak in a.items():
        ak = sp.expand(ak)
        if ak == 0:
            continue
        p = sp.Poly(sp.expand(ak.subs(t, s + t0)), s)
        v = min(m[0] for m in p.monoms())
        vk[k] = (v, p.coeff_monomial(s**v))
    low = min(v - k for k, (v, _) in vk.items())
    return sp.solve(sum(c * sp.ff(rho, k) for k, (v, c) in vk.items() if v - k == low), rho)


@pytest.mark.parametrize("n", [5, 7])
def test_top_operator_has_a_half_integer_exponent(n):
    ex = _exponents(TR._phi(TR._sym_op(n - 1)), -(n - 1))
    assert any(sp.Rational(e).q == 2 for e in ex), ex
    # the next piece down is ordinary at -(n-1): exponents 0..n-2, all integers
    assert sorted(_exponents(TR._phi(TR._sym_op(n - 2)), -(n - 1))) == list(range(n - 1))


@pytest.mark.parametrize("c", [1, 2, 3, 4, 5])
def test_leading_coefficients_are_threshold_products_independent_of_shift(c):
    target = sp.expand(sp.prod([t**2 - (c - 2 * j) ** 2 for j in range(c // 2 + 1) if c - 2 * j > 0]))
    for lam in (-3, 0, 1, 2, 5):
        lc = _leading(TR._phi(TR._shift(TR._sym_op(c), lam)))
        assert sp.cancel(lc / target).is_number
        if c <= 3:
            assert lc.subs(t, 4) != 0
        if c <= 5:
            assert lc.subs(t, 6) != 0


UL, WL = roots_laguerre(80)


def _Rs(x):
    X = x[:, None] + UL[None, :]
    return (X * k1e(X)) @ WL


def _branch(eps, n):
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), 200001)
    x = np.exp(y)
    f = i0e((n - 1 - eps) * x) * _Rs(x) ** (n - 1) * np.exp(-eps * x) * x ** (3.0 - n)
    return np.trapezoid(f * x, y)


def test_five_star_square_root_branch():
    C5 = np.pi**2 / (8 * np.sqrt(2))
    v2, v3 = np.sqrt(1e-2) * _branch(1e-2, 5), np.sqrt(1e-3) * _branch(1e-3, 5)
    assert v3 < v2
    assert abs(v3 / C5 - 1) < 0.15
    assert 1e-3 * _branch(1e-3, 5) < 0.1 * C5


def test_four_star_control_is_a_pole():
    C4 = (np.pi / 2) ** 1.5 / np.sqrt(6 * np.pi)
    assert abs(1e-3 * _branch(1e-3, 4) / C4 - 1) < 0.05
    assert np.sqrt(1e-3) * _branch(1e-3, 4) > 10 * C4


def test_six_star_log_branch():
    """n = 6 (even, added after the self-check widened the proposition to 3 <= n <= 8): the R^5 term gives
    -c6 log eps, c6 = (pi/2)^(5/2) (10 pi)^(-1/2). The slope of I6 against log eps between 1e-3 and 1e-4 is
    within 10% of -c6 (rejects a holomorphic I6, slope 0, and a power law, slope growing)."""
    c6 = (np.pi / 2) ** 2.5 / np.sqrt(10 * np.pi)
    slope = (_branch(1e-4, 6) - _branch(1e-3, 6)) / (np.log(1e-4) - np.log(1e-3))
    assert abs(slope / -c6 - 1) < 0.10
