"""Backing tests for Sprint Topos-3 (Paper 57 SS open, Bohr-site remark,
two-center leg): the two-center frame meet is exactly the m-grading.
SELF-CONTAINED (mirrors debug/compute_topos3_exact_meet.py); every
load-bearing step exact (Mulliken-Ruedenberg auxiliary integrals; the
overlap is e^{-p}(U e^q + V e^{-q}) with U, V exact rationals, so
S = 0 <=> U = V = 0 by Lindemann independence).

Pins:
  1. Machinery certification: normalized <1s|1s>(R) equals the classical
     (1 + R + R^2/3) e^{-R} as an EXACT RATIONAL IDENTITY at R = 1, 2, 7/2
     (7/3, 13/3, 103/12).
  2. Classical cancellation: the prolate-coordinate integrand clears to a
     polynomial (asserted inside overlap_UV for every pair used).
  3. Two-center meet = m-grading (fast leg): at (Z2=1, R=1, n_max=3) the
     m = 1 block (d = 3) and m = 2 block have FULL exact support; ladder
     arithmetic 14 -> 9 -> 5.
  4. Full m = 0 block (d = 6) full support at (Z2=1, R=1) and (Z2=3, R=1).
"""

from fractions import Fraction
from math import comb, factorial

import sympy as sp

xi, eta = sp.symbols("xi eta")


def radial_poly_coeffs(Z: Fraction, n: int, l: int):
    lam = sp.Rational(2 * Z.numerator, Z.denominator * n)
    lag = [sp.Rational((-1) ** j * comb(n + l, n - l - 1 - j), factorial(j))
           for j in range(n - l)]
    return {l + j: c * lam ** j for j, c in enumerate(lag)}, \
        sp.Rational(Z.numerator, Z.denominator * n)


def radial_norm_sq(Z: Fraction, n: int, l: int):
    coeffs, a = radial_poly_coeffs(Z, n, l)
    return sp.nsimplify(sum(
        c1 * c2 * sp.factorial(k1 + k2 + 2) / (2 * a) ** (k1 + k2 + 3)
        for k1, c1 in coeffs.items() for k2, c2 in coeffs.items()))


def theta_norm_sq(l: int, m: int):
    return sp.Rational(2 * factorial(l + m), (2 * l + 1) * factorial(l - m))


def G_lm(l: int, m: int, x):
    t = sp.Symbol("_t")
    return sp.diff(sp.legendre(l, t), t, m).subs(t, x)


def overlap_UV(Z1, n1, l1, Z2, n2, l2, m, R):
    m = abs(m)
    Rs = sp.Rational(R.numerator, R.denominator)
    half = Rs / 2
    c1, a = radial_poly_coeffs(Z1, n1, l1)
    c2, b = radial_poly_coeffs(Z2, n2, l2)
    ra, rb = half * (xi + eta), half * (xi - eta)
    ct_a = (1 + xi * eta) / (xi + eta)
    ct_b = (xi * eta - 1) / (xi - eta)
    s2 = (xi ** 2 - 1) * (1 - eta ** 2)
    expr = (sum(c * ra ** k for k, c in c1.items())
            * sum(c * rb ** k for k, c in c2.items())
            * s2 ** m / ((xi + eta) * (xi - eta)) ** m
            * G_lm(l1, m, ct_a) * G_lm(l2, m, ct_b)
            * half ** 3 * (xi ** 2 - eta ** 2))
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    assert den == 1                              # classical cancellation
    poly = sp.Poly(sp.expand(num), xi, eta)
    p = (a + b) * Rs / 2
    q = (a - b) * Rs / 2
    max_i = max(mm[0] for mm in poly.monoms())
    max_j = max(mm[1] for mm in poly.monoms())
    a_c = {0: 1 / p}
    for i in range(1, max_i + 1):
        a_c[i] = (1 + i * a_c[i - 1]) / p
    if q != 0:
        u, v = {0: 1 / q}, {0: -1 / q}
        for j in range(1, max_j + 1):
            u[j] = sp.Rational((-1) ** j) / q + j * u[j - 1] / q
            v[j] = sp.Rational(-1) / q + j * v[j - 1] / q
        U = sum(c * a_c[i] * u[j] for (i, j), c in
                zip(poly.monoms(), poly.coeffs()))
        V = sum(c * a_c[i] * v[j] for (i, j), c in
                zip(poly.monoms(), poly.coeffs()))
        return sp.nsimplify(U), sp.nsimplify(V)
    U = sum(c * a_c[i] * sp.Rational(2, j + 1)
            for (i, j), c in zip(poly.monoms(), poly.coeffs()) if j % 2 == 0)
    return sp.nsimplify(U), None


def nonzero(U, V):
    return (U != 0) if V is None else (U != 0 or V != 0)


M_BLOCKS_NMAX3 = {0: [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)],
                  1: [(2, 1), (3, 1), (3, 2)],
                  2: [(3, 2)]}


def test_machinery_exact_crosscheck():
    expect = {Fraction(1): sp.Rational(7, 3), Fraction(2): sp.Rational(13, 3),
              Fraction(7, 2): sp.Rational(103, 12)}
    nsq = radial_norm_sq(Fraction(1), 1, 0) * theta_norm_sq(0, 0)
    for R, want in expect.items():
        U, V = overlap_UV(Fraction(1), 1, 0, Fraction(1), 1, 0, 0, R)
        assert V is None                          # equal rates -> q = 0
        assert sp.simplify(U / nsq - want) == 0   # exact rational identity


def test_two_center_meet_m1_m2_blocks_and_ladder():
    R = Fraction(1)
    for m in (1, 2):
        states = M_BLOCKS_NMAX3[m]
        for aa in states:
            for bb in states:
                U, V = overlap_UV(Fraction(1), aa[0], aa[1],
                                  Fraction(1), bb[0], bb[1], m, R)
                assert nonzero(U, V), (m, aa, bb)   # full exact support
    # ladder arithmetic at n_max = 3
    N = sum((2 * l + 1) * (3 - l) for l in range(3))
    same_center = sum(2 * l + 1 for l in range(3))
    two_center = len(range(-2, 3))
    assert (N, same_center, two_center) == (14, 9, 5)


def test_two_center_meet_m0_block_full_support():
    states = M_BLOCKS_NMAX3[0]
    for Z2 in (Fraction(1), Fraction(3)):
        for aa in states:
            for bb in states:
                U, V = overlap_UV(Fraction(1), aa[0], aa[1],
                                  Z2, bb[0], bb[1], 0, Fraction(1))
                assert nonzero(U, V), (Z2, aa, bb)
