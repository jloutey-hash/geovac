"""
Regression tests locking the two verified Paper 2 corrections applied in the
confidence-review-infra sprint (2026-06-01).
See debug/sprint_confidence_review_infra_memo.md.
"""
import sympy as sp
from mpmath import mp, mpf, pi, zeta, polyroots


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
    """Paper 2 §VIII.D: with the corrected lower limit k=0, the ratio
    B_formal(m)/N(m) = d collapses to numerator d(m-2)(m+d+2), root m=2 -
    matching the paper's stated quadratic. k=1 does NOT give this."""
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
