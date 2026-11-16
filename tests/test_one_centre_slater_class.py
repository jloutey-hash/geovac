"""Guard: the number class of one-centre two-electron integrals (inventory row 1; Paper 58 sec:census).

Written 2026-11-15 as a separate step after the fix (CLAUDE.md sec. 9).  Claims guarded, with the wrong answer each
rejects:

(1) radial R^k at unit exponent is q sqrt(s) with s square-free, s > 1 when principal quantum numbers mix:
    R^0(1s 1s; 1s 2s) = 4096 sqrt(2) / 64827.  Checked against an independent in-test quadrature to 1e-25; rejects
    the pre-fix behaviour, which returned a double-precision Fraction approximation (error ~ 7e-17).
(2) compute_rk_algebraic returns an exact Fraction only for rational values and raises otherwise; rejects a silent
    approximation.
(3) the Gaunt factor makes one-centre integrals irrational from l = 2 even with rational radial parts:
    c^1(s, p_-1) c^1(p_-1, d_0) = +-sqrt(5)/15 (sympy's Gaunt, independent of geovac).
Fire-test plants: the comment block at the end.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import pytest

from geovac.hypergeometric_slater import compute_rk_algebraic, compute_rk_exact, compute_rk_float


def _R(n: int, l: int, r):
    """Hydrogenic radial function at Z = 1, written here independently of geovac."""
    N = mp.sqrt((mp.mpf(2) / n) ** 3 * mp.factorial(n - l - 1) / (2 * n * mp.factorial(n + l)))
    x = 2 * r / n
    m, a = n - l - 1, 2 * l + 1     # generalized Laguerre L_m^a, explicit sum
    lag = sum((-1) ** j * mp.binomial(m + a, m - j) * x ** j / mp.factorial(j) for j in range(m + 1))
    return N * x ** l * mp.exp(-r / n) * lag


def _rk_quad(key, dps: int = 30):
    n1, l1, n3, l3, n2, l2, n4, l4, k = key
    with mp.workdps(dps):
        def inner(r1):
            f = lambda r2: _R(n2, l2, r2) * _R(n4, l4, r2) * r2 ** 2
            a = mp.quad(lambda r2: f(r2) * r2 ** k, [0, r1]) / r1 ** (k + 1)
            b = mp.quad(lambda r2: f(r2) / r2 ** (k + 1), [r1, mp.inf]) * r1 ** k
            return a + b
        return mp.quad(lambda r1: _R(n1, l1, r1) * _R(n3, l3, r1) * r1 ** 2 * inner(r1), [0, 2, 8, mp.inf])


def _val(q: Fraction, s: int):
    return mp.mpf(q.numerator) / q.denominator * mp.sqrt(s)


IRRATIONAL = [(1, 0, 1, 0, 1, 0, 2, 0, 0), (2, 1, 2, 1, 1, 0, 3, 0, 0),
              (1, 0, 1, 0, 1, 0, 5, 0, 0)]   # n = 5: the exact-dispatch float path
RATIONAL = [(1, 0, 2, 0, 1, 0, 2, 0, 0), (1, 0, 1, 0, 2, 1, 2, 1, 0), (1, 0, 1, 0, 2, 0, 2, 0, 0)]


def test_mixed_n_value_is_4096_sqrt2_over_64827():
    assert compute_rk_exact(1, 0, 1, 0, 1, 0, 2, 0, 0) == (Fraction(4096, 64827), 2)


@pytest.mark.slow
@pytest.mark.parametrize('key', IRRATIONAL + RATIONAL)
def test_exact_value_matches_independent_quadrature(key):  # ~2 min per case (nested 30-digit quadrature)
    q, s = compute_rk_exact(*key)
    with mp.workdps(30):
        assert abs(_val(q, s) / _rk_quad(key) - 1) < mp.mpf('1e-22'), (key, q, s)


@pytest.mark.parametrize('key', IRRATIONAL)
def test_irrational_cases_are_irrational_and_algebraic_path_refuses(key):
    q, s = compute_rk_exact(*key)
    assert s > 1 and q != 0
    with pytest.raises(ValueError):
        compute_rk_algebraic(*key)
    # the float path carries the root: agrees with the exact value to double precision
    with mp.workdps(30):
        assert abs(compute_rk_float(*key) / _val(q, s) - 1) < 1e-14


@pytest.mark.parametrize('key', RATIONAL)
def test_rational_cases_are_exact_fractions(key):
    q, s = compute_rk_exact(*key)
    assert s == 1
    assert compute_rk_algebraic(*key) == q


def test_gaunt_factor_irrational_from_l2():
    from sympy import Rational, simplify, sqrt
    from sympy.physics.wigner import wigner_3j

    def ck(k, l, m, lp, mp_):   # c^k(lm, l'm') = (-1)^m sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0)(l k l'; -m, m-m', m')
        return (-1) ** abs(m) * sqrt((2 * l + 1) * (2 * lp + 1)) * wigner_3j(l, k, lp, 0, 0, 0) *             wigner_3j(l, k, lp, -m, m - mp_, mp_)
    prod = simplify(ck(1, 0, 0, 1, -1) * ck(1, 1, -1, 2, 0))
    assert simplify(prod ** 2 - Rational(1, 45)) == 0      # (sqrt(5)/15)^2
    assert prod.is_rational is False


# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_one_centre_slater_class.py --plant-in
# geovac/hypergeometric_slater.py --plant "OLD=>NEW"); 2026-11-15:
#     float fallback restored   "    norm_q, norm_s = _sqrt_split(nsq_prod)=>    import math as _m; norm_q, norm_s = Fraction(_m.sqrt(float(nsq_prod))).limit_denominator(10**18), 1"
#     root dropped              "    return norm_q * integral, norm_s=>    return norm_q * integral, 1"
#     odd exponent kept twice   "        sq *= prime ** (e // 2)=>        sq *= prime ** ((e + 1) // 2)"
#     float path drops root     "return float(q_exact) if s_exact == 1 else float(q_exact) * math.sqrt(s_exact)=>return float(q_exact)"
#   all four fired, 2026-11-15.  Also fired (separate guard, same defect): tests/test_certified_reference_values.py
#   -k test_fast_subset_reproduces_the_stored_table with --plant-in benchmarks/certified_reference/entries_helium.py,
#   restoring a double-precision Fraction for R^k after "if Rq == 0: continue" (the stored helium values then disagree
#   at 4e-18 against the 1e-28 tolerance).
