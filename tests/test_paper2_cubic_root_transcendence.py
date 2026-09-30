"""Paper 2 / Paper 18: the cubic's root alpha_K is transcendental (a one-line proof).

The proof uses two facts, checked here:
  (1) K = pi (B + F - Delta) with B = 42, F = pi^2/6, Delta = 1/40 is the
      nonconstant Q-polynomial (1679/40) pi + pi^3/6 (so it is transcendental by
      Lindemann: pi is not a root of any nonzero rational polynomial);
  (2) at the cubic's smallest positive root, K = alpha_K^2 + 1/alpha_K, so an
      algebraic alpha_K would force K algebraic.
Plus the counterexample printed in the papers, (x-1)(x-pi+1) = x^2 - pi x + (pi-1),
which vanishes at x = 1 (the earlier printed x^2 - pi x + pi(pi-1) does NOT).
Rejects: a wrong expansion of K, a root that does not satisfy the cubic, and the
retired counterexample.  Nothing here concerns the physical alpha.
"""
import mpmath as mp
import sympy as sp


def test_K_is_nonconstant_rational_polynomial_in_pi():
    pi = sp.pi
    K = pi * (42 + pi**2 / 6 - sp.Rational(1, 40))
    poly = sp.Poly(sp.expand(K), pi)
    assert poly.all_coeffs() == [sp.Rational(1, 6), 0, sp.Rational(1679, 40), 0]
    assert poly.degree() == 3


def test_root_satisfies_K_equals_alpha_sq_plus_inverse():
    mp.mp.dps = 50
    K = mp.pi * (42 + mp.pi**2 / 6 - mp.mpf(1) / 40)
    a = mp.findroot(lambda x: x**3 - K * x + 1, mp.mpf(1) / 137)
    assert a > 0 and abs(a**3 - K * a + 1) < mp.mpf(10) ** -45
    assert abs(K - (a**2 + 1 / a)) < mp.mpf(10) ** -40
    # smallest positive root (the other positive root is ~ sqrt(K))
    assert a < 1


def test_printed_counterexample_and_retired_one():
    x, pi = sp.symbols("x"), sp.pi
    good = x**2 - pi * x + (pi - 1)
    assert sp.expand(good - (x - 1) * (x - pi + 1)) == 0
    assert sp.simplify(good.subs(x, 1)) == 0
    retired = x**2 - pi * x + pi * (pi - 1)
    assert sp.simplify(retired.subs(x, 1)) != 0
