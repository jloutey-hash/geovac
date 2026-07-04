"""
Regression tests locking the two verified Paper 2 corrections applied in the
confidence-review-infra sprint (2026-06-01).
See debug/sprint_confidence_review_infra_memo.md.
"""
import pytest
import sympy as sp
from mpmath import mp, mpf, pi, zeta, polyroots


@pytest.fixture(autouse=True)
def _mp_dps_hermetic():
    """Restore mpmath precision after each test (the corpus-wide mp.dps
    collection-order hermeticity discipline, v4.64.0)."""
    saved = mp.dps
    yield
    mp.dps = saved


def test_paper2_circulant_roots_satisfy_hermiticity_premise():
    """Paper 2 footnote 2 premise: all three roots of s^3 - K s + 1 = 0
    are real with s_i^2 < 4K/3, which forces the circulant's b, c to be
    complex conjugates (Sec. VI D).  Added 2026-07-04 (cert run): the
    footnote's 'proven' rested on this numerically-true but previously
    untested premise."""
    mp.dps = 30
    K = pi * (42 + zeta(2) - mpf(1) / 40)
    roots = polyroots([1, 0, -K, 1])
    assert len(roots) == 3
    for r in roots:
        assert abs(r.imag) < mpf(10) ** -20      # all real
        assert r.real ** 2 < 4 * K / 3           # the Hermiticity premise
    # and with sensible margin: the closest root sits ~25% below the bound
    margin = max(float(r.real ** 2 / (4 * K / 3)) for r in roots)
    assert margin < 0.80


def test_paper2_cubic_root_value():
    """Paper 2 Table II / abstract: the cubic 1/a + a^2 = K with
    K = pi(42 + zeta(2) - 1/40) has physical root 1/alpha = 137.036011...,
    ABOVE CODATA (not the pre-fix 137.035987, which was below)."""
    mp.dps = 30
    K = pi * (42 + zeta(2) - mpf(1) / 40)
    roots = polyroots([1, 0, -K, 1])  # a^3 - K a + 1 = 0
    a = [r.real for r in roots
         if abs(r.imag) < mpf(10) ** -20 and 0 < r.real < mpf('0.1')][0]
    inv_alpha = 1 / a
    # paper's corrected printed value
    assert abs(inv_alpha - mpf('137.036011')) < mpf('1e-5')
    # matches CODATA-2018 to ~8.8e-8 (the headline)
    codata = mpf('137.035999084')
    assert abs(inv_alpha - codata) / codata < mpf('1e-7')
    # the corrected value is ABOVE CODATA (the sign of the fix)
    assert inv_alpha > codata


def test_paper2_viiid_identity_k0():
    """Paper 2 §VIII.E (S³ specificity): with the corrected lower limit k=0,
    the ratio B_formal(m)/N(m) = d collapses to numerator d(m-2)(m+d+2),
    root m=2 - matching the paper's stated quadratic. k=1 does NOT give
    this (asserted below since the 1st-cert remediation, 2026-07-03)."""
    d, m, k = sp.symbols('d m k', positive=True)
    g = sp.binomial(k + d, d) - sp.binomial(k + d - 2, d)
    absL = k * (k + d - 1)
    N = sp.summation(g, (k, 0, m))
    Bf = sp.Rational(1, 2) * sp.summation(g * absL, (k, 0, m))
    num = sp.simplify(sp.numer(sp.together(sp.simplify(Bf / N - d))))
    # roots of the paper's quadratic m^2 + d m - 2(d+2) = 0 are m=2 and m=-(d+2)
    assert sp.simplify(num.subs(m, 2)) == 0
    assert sp.simplify(num.subs(m, -(d + 2))) == 0
    # N(2) = (d+1)(d+4)/2 under the corrected k=0 lower limit
    assert sp.simplify(N.subs(m, 2) - (d + 1) * (d + 4) / 2) == 0
    # the falsifiable converse: the k=1 lower limit does NOT vanish at m=2
    N1 = sp.summation(g, (k, 1, m))
    Bf1 = sp.Rational(1, 2) * sp.summation(g * absL, (k, 1, m))
    num1 = sp.simplify(sp.numer(sp.together(sp.simplify(Bf1 / N1 - d))))
    assert sp.simplify(num1.subs(m, 2)) != 0
