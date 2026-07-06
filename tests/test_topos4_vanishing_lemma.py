"""Sprint Topos-4 — general-parameter vanishing lemma for inter-frame
hydrogenic overlaps (Paper 57 sec:open_bohr, Family-1 follow-on).

Pins the lemma that closes Topos-2's panel-scoped (=>) direction:

  For same-l inter-frame overlaps, S(n1,l;n2,l;Z1,Z2) factors as a positive
  prefactor times an explicit polynomial P_{n1,n2,l}(t) of degree
  D = n1+n2-2l-2 in the rate variable t = (Z1/n1 - Z2/n2)/(Z1/n1 + Z2/n2).
  S = 0  <=>  t is a root of P.  The roots are three structural families:
    (i)   t = 0                    rate coincidence   <=> |n1-n2| >= 2
    (ii)  t = (n2-n1)/(n2+n1)      same-Z orthogonality (always, n1 != n2)
    (iii) residual roots           generically irrational; the "sporadic"
                                   zeros are exactly the RATIONAL residual
                                   roots (in Z,n<7: (5,5,3)->+-1/3 and
                                   (4,5,2)->+-1/2).

Self-contained + exact (no import from debug/).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial

import sympy as sp


# ---- exact hydrogenic radial overlap (unnormalized; support-exact) -------

def _genlag(k: int, alpha: int):
    return [Fraction((-1) ** j * comb(k + alpha, k - j), factorial(j))
            for j in range(k + 1)]


def overlap(Z1: int, n1: int, Z2: int, n2: int, l: int) -> Fraction:
    """int_0^inf R_{n1 l}^{Z1} R_{n2 l}^{Z2} r^2 dr, exact Fraction."""
    lag1, lag2 = _genlag(n1 - l - 1, 2 * l + 1), _genlag(n2 - l - 1, 2 * l + 1)
    p, q = Fraction(Z1, n1), Fraction(Z2, n2)
    rate = p + q
    tot = Fraction(0)
    for j1, c1 in enumerate(lag1):
        a1 = c1 * (2 * p) ** j1
        for j2, c2 in enumerate(lag2):
            a2 = c2 * (2 * q) ** j2
            k = (l + j1) + (l + j2) + 2
            tot += a1 * a2 * Fraction(factorial(k)) / rate ** (k + 1)
    return tot


# ---- P(t) = N(1+t, 1-t), N = homogeneous-degree-D overlap numerator ------

def poly_in_t(n1: int, n2: int, l: int):
    t = sp.symbols('t', real=True)
    lag1, lag2 = _genlag(n1 - l - 1, 2 * l + 1), _genlag(n2 - l - 1, 2 * l + 1)
    D = (n1 - l - 1) + (n2 - l - 1)
    p, q = 1 + t, 1 - t
    N = sp.Integer(0)
    for j1, c1 in enumerate(lag1):
        for j2, c2 in enumerate(lag2):
            k = 2 * l + j1 + j2 + 2
            coeff = sp.Rational(c1.numerator, c1.denominator) * \
                sp.Rational(c2.numerator, c2.denominator) * \
                2 ** (j1 + j2) * sp.factorial(k)
            N += coeff * p ** j1 * q ** j2 * (p + q) ** (D - j1 - j2)
    P = sp.Poly(sp.expand(N), t)
    if P.degree() >= 1:
        P = sp.Poly(sp.primitive(P.as_expr())[1], t)
    return P, D, t


# ---- tests ---------------------------------------------------------------

def test_convention_pins():
    """The three exact Topos-2 pins (unnormalized convention control)."""
    assert overlap(1, 1, 2, 1, 0) == Fraction(2, 27)
    assert overlap(1, 1, 2, 2, 0) == Fraction(-1, 4)
    assert overlap(1, 2, 2, 2, 1) == Fraction(256, 81)


def test_boundary_case_5_5_3():
    """The Topos-2 boundary 'sporadic' zero is P_{5,5,3} = 1 - 9 t^2."""
    P, D, t = poly_in_t(5, 5, 3)
    assert D == 2
    assert sp.expand(P.as_expr() - (1 - 9 * t ** 2)) == 0
    # (Z=1,n=5) vs (Z'=2,n=5): rates 1/5, 2/5 -> t = -1/3, a root
    t_val = (Fraction(1, 5) - Fraction(2, 5)) / (Fraction(1, 5) + Fraction(2, 5))
    assert t_val == Fraction(-1, 3)
    assert P.eval(sp.Rational(-1, 3)) == 0
    assert overlap(1, 5, 2, 5, 3) == 0            # the overlap itself vanishes
    # symmetric-Jacobi (Gegenbauer) identity P_2^{(3,3)} for this D=2 case
    assert sp.expand(P.as_expr()
                     - sp.Rational(-4, 5) * sp.jacobi(2, 3, 3, t)) == 0


def _matches_symmetric_jacobi(P, D, t):
    """True iff P is proportional to a symmetric Jacobi P_D^{(a,a)}, a in 0..8."""
    for a in range(0, 9):
        J = sp.Poly(sp.expand(sp.jacobi(D, a, a, t)), t)
        if J.degree() == D and \
                sp.expand(P.as_expr() * J.LC() - P.LC() * J.as_expr()) == 0:
            return True
    return False


def test_gegenbauer_scope():
    """Honest scope of the symmetric-Jacobi identification:
    P = Gegenbauer P_2^{(l,l)} holds ONLY for the single-radial-node family
    n1=n2=l+2 (D=2); for n1=n2 with D>=4, P is even in t (frame-swap
    symmetry) but matches NO symmetric Jacobi (the claim is not general)."""
    t = sp.symbols('t', real=True)
    # (a) D=2 single-radial-node family: P IS the symmetric Jacobi P_2^{(l,l)}
    for l in range(0, 4):
        P, D, _ = poly_in_t(l + 2, l + 2, l)
        assert D == 2
        J = sp.Poly(sp.expand(sp.jacobi(2, l, l, t)), t)
        assert sp.expand(P.as_expr() * J.LC() - P.LC() * J.as_expr()) == 0
    # (b) higher-D equal-n cases: even in t, but NOT a symmetric Jacobi
    for (n, l) in [(3, 0), (4, 0), (4, 1)]:
        P, D, _ = poly_in_t(n, n, l)
        assert D >= 4
        assert sp.expand(P.as_expr() - P.as_expr().subs(t, -t)) == 0   # even
        assert not _matches_symmetric_jacobi(P, D, t)                  # not Jacobi


def test_rate_coincidence_criterion():
    """t = 0 is a root of P  <=>  |n1 - n2| >= 2  (structural, over a grid)."""
    for n1 in range(1, 7):
        for n2 in range(1, 7):
            for l in range(0, min(n1, n2)):
                P, D, t = poly_in_t(n1, n2, l)
                assert (P.eval(0) == 0) == (abs(n1 - n2) >= 2)


def test_orthogonality_root():
    """t = (n2-n1)/(n2+n1) is always a root for n1 != n2
    (= standard same-Z hydrogenic orthogonality <n1 l|n2 l> = 0)."""
    for n1 in range(1, 7):
        for n2 in range(1, 7):
            if n1 == n2:
                continue
            for l in range(0, min(n1, n2)):
                P, D, t = poly_in_t(n1, n2, l)
                assert P.eval(sp.Rational(n2 - n1, n2 + n1)) == 0
                # and the actual same-Z overlap vanishes
                assert overlap(1, n1, 1, n2, l) == 0


def test_integer_zero_census_fully_accounted():
    """Every integer zero (Z,n<7) is uniquely rate-coincidence /
    orthogonality / residual, and the residuals are exactly the
    (5,5,3) rational-residual (+-1/3) hits."""
    rate_coin = orth = 0
    residual = []
    N = 7
    for n1 in range(1, N):
        for n2 in range(n1, N):
            for l in range(0, min(n1, n2)):
                for Z1 in range(1, N):
                    for Z2 in range(1, N):
                        if overlap(Z1, n1, Z2, n2, l) != 0:
                            continue
                        p, q = Fraction(Z1, n1), Fraction(Z2, n2)
                        if p == q:
                            rate_coin += 1
                        elif Z1 == Z2:
                            orth += 1
                        else:
                            residual.append((Z1, n1, Z2, n2, l))
    assert rate_coin == 37
    assert orth == 210
    assert len(residual) == 6
    # every residual zero lives on the (5,5,3) case
    assert all((n1, n2, l) == (5, 5, 3) for (_, n1, _, n2, l) in residual)


def test_residual_is_the_only_sporadic_family():
    """Across Z,n<7, the ONLY (n1,n2,l) with a rational residual root are
    (4,5,2) [+-1/2] and (5,5,3) [+-1/3]."""
    cases = []
    for n1 in range(1, 7):
        for n2 in range(n1, 7):
            for l in range(0, min(n1, n2)):
                P, D, t = poly_in_t(n1, n2, l)
                Pt = sp.Poly(P.as_expr(), t)
                while Pt.degree() >= 1 and Pt.eval(0) == 0:
                    Pt = sp.div(Pt, sp.Poly(t, t))[0]
                if n1 != n2:
                    to = sp.Rational(n2 - n1, n2 + n1)
                    while Pt.degree() >= 1 and Pt.eval(to) == 0:
                        Pt = sp.div(Pt, sp.Poly(t - to, t))[0]
                rat = [r for r in (Pt.ground_roots() if Pt.degree() >= 1
                                   else {}) if abs(r) < 1]
                if rat:
                    cases.append((n1, n2, l, sorted(str(r) for r in rat)))
    assert cases == [(4, 5, 2, ['-1/2', '1/2']),
                     (5, 5, 3, ['-1/3', '1/3'])]
