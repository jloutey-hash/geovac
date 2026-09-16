r"""Backing test for Paper 58 sec:census g-row: the permitted-density inflation
factor R(n_max) is BOUNDED and converges to R_inf = 289777/18471 ~ 15.69.

Companion to tests/test_paper58_census.py (which pins the counts at n_max=2,3).
This test pins the CLOSED-FORM result added in the v5.12.4 investigation:

    c(mu)      = (n-|mu|)(n-|mu|+1)/2                 orbitals/center at |m|=mu
    one_center = 2 * sum_t ( sum_mu c(mu) c(t-mu) )^2  (degree-11 polynomial)
    m_rule     = 8 * one_center
    genuine    = 7 * one_center + builder
    builder    = period-2 degree-11 quasi-polynomial in n_max
    R          = genuine / builder = 1 + 7 * one_center/builder  ->  289777/18471

The enumerators here reproduce the census predicate of test_paper58_census.py
(validated at the two published anchors). builder uses the O(n^5) reformulation
(sum over (l_p,l_q,l_r,l_s), closed-form m-count) so the fit reaches large
n_max. The derived leading coefficients are one_center: 19379/1247400,
builder: 6157/831600, giving the single-valued limit.
"""
from __future__ import annotations

from fractions import Fraction

import pytest

R_INF = Fraction(289777, 18471)
LEAD_ONE_CENTER = Fraction(19379, 1247400)
LEAD_BUILDER = Fraction(6157, 831600)


# --- enumerators (census predicate of test_paper58_census.py) ---------------
def _c(mu, n_max):
    mu = abs(mu)
    return 0 if mu > n_max - 1 else (n_max - mu) * (n_max - mu + 1) // 2


def sum_d_squared(n_max):
    ms = range(-(n_max - 1), n_max)
    return sum(sum(_c(m, n_max) * _c(t - m, n_max) for m in ms) ** 2
               for t in range(-2 * (n_max - 1), 2 * (n_max - 1) + 1))


def one_center(n_max):
    return 2 * sum_d_squared(n_max)


def m_rule(n_max):
    return 16 * sum_d_squared(n_max)


def _pmcount(la, lb, M):
    """#{(ma,mb): |ma|<=la, |mb|<=lb, mb-ma=M}."""
    lo, hi = max(-la, -lb - M), min(la, lb - M)
    return max(0, hi - lo + 1)


def _start(lo, top):
    return lo if (top - lo) % 2 == 0 else lo + 1


def builder(n_max):
    """2 * (single-center quartets, m-rule, with intersecting multipole L-ranges)."""
    tot = 0
    for lp in range(n_max):
        rp = n_max - lp
        for lq in range(n_max):
            rpq = rp * (n_max - lq)
            top1, adl1, par1 = lp + lq, abs(lp - lq), (lp + lq) % 2
            for lr in range(n_max):
                rr = n_max - lr
                for ls in range(n_max):
                    top2 = lr + ls
                    if (top2 % 2) != par1:
                        continue
                    adl2 = abs(lr - ls)
                    tmin = min(top1, top2)
                    w = rpq * rr * (n_max - ls)
                    inner = 0
                    for mu in range(tmin + 1):
                        if max(_start(max(adl1, mu), top1),
                               _start(max(adl2, mu), top2)) > tmin:
                            continue
                        for M1 in ((0,) if mu == 0 else (mu, -mu)):
                            inner += _pmcount(lp, lq, M1) * _pmcount(lr, ls, -M1)
                    tot += w * inner
    return 2 * tot


def genuine(n_max):
    return 7 * one_center(n_max) + builder(n_max)


# --- tests ------------------------------------------------------------------
def test_anchors_match_published_census():
    """The closed-form enumerators reproduce test_paper58_census.py's g row."""
    assert (genuine(2), builder(2)) == (2944, 214)      # 13.8x
    assert (genuine(3), builder(3)) == (114280, 7600)   # 15.0x


def test_structural_identities():
    for nm in range(2, 11):
        oc = one_center(nm)
        assert m_rule(nm) == 8 * oc
        # genuine = m_rule - one_center + builder = 7*one_center + builder
        assert genuine(nm) == m_rule(nm) - oc + builder(nm)
        assert genuine(nm) == 7 * oc + builder(nm)


def test_inflation_monotone_increasing_and_bounded():
    Rs = [Fraction(genuine(nm), builder(nm)) for nm in range(2, 13)]
    assert all(b > a for a, b in zip(Rs, Rs[1:])), "R must strictly increase"
    assert all(R < R_INF for R in Rs), "R must stay below the limit"
    # the approach is monotone from below and the gap shrinks
    gaps = [R_INF - R for R in Rs]
    assert all(g2 < g1 for g1, g2 in zip(gaps, gaps[1:]))
    assert float(Rs[0]) == pytest.approx(13.757, abs=1e-3)   # n_max=2
    assert Rs[-1] > Fraction(1567, 100)                       # n_max=12 > 15.67


def test_limit_arithmetic_from_leading_coefficients():
    """R_inf = 1 + 7 * lead(one_center)/lead(builder)."""
    assert 1 + 7 * LEAD_ONE_CENTER / LEAD_BUILDER == R_INF
    assert float(R_INF) == pytest.approx(15.6882, abs=1e-4)


@pytest.mark.slow
def test_leading_coefficients_derived_from_enumerator():
    """Independently confirm the limit: fit degree-11 polynomials to the
    enumerator (one_center single; builder period-2), verify on held-out
    points, and check the leading coefficients give R_inf."""
    sp = pytest.importorskip("sympy")
    n = sp.Symbol("n")
    NS = list(range(2, 29))

    oc_vals = [(x, one_center(x)) for x in NS]
    oc_poly = sp.interpolate(oc_vals[:12], n)
    assert all(int(oc_poly.subs(n, x)) == v for x, v in oc_vals), "one_center deg-11"
    assert sp.LC(sp.Poly(sp.expand(oc_poly), n)) == sp.Rational(19379, 1247400)

    b_vals = [(x, builder(x)) for x in NS]
    leads = set()
    for par in (0, 1):
        sub = [(x, v) for x, v in b_vals if x % 2 == par]
        poly = sp.interpolate(sub[:12], n)
        held = sub[12:]
        assert held, "need held-out points for the builder branch"
        assert all(int(poly.subs(n, x)) == v for x, v in held), "builder branch fit"
        leads.add(sp.LC(sp.Poly(sp.expand(poly), n)))
    assert leads == {sp.Rational(6157, 831600)}, "both branches share the lead coeff"

    R_inf = 1 + 7 * sp.Rational(19379, 1247400) / sp.Rational(6157, 831600)
    assert R_inf == sp.Rational(289777, 18471)
