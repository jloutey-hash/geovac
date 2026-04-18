"""Test suite for Track RH-F: Galois / number-field structure of Ihara zeros.

These tests verify concrete factorizations and Galois-group computations
against the closed forms from Paper 29.

Run: pytest tests/test_galois_ihara.py -v
"""

import pytest
import sympy as sp
from sympy.polys.numberfields import galois_group

s = sp.symbols("s")


def _n_poly_factors(expr, var):
    """Count number of polynomial-of-var factors (with multiplicity)."""
    n = 0
    for a in sp.Mul.make_args(expr):
        base, exp = a.as_base_exp()
        try:
            deg = sp.Poly(base, var).degree()
        except (sp.PolynomialError, sp.GeneratorsNeeded):
            deg = 0
        if deg > 0:
            try:
                mult = int(exp)
            except TypeError:
                mult = 1
            n += mult
    return n


class TestCyclotomicFactorizations:
    """The S^3 Coulomb Ihara zeta at max_n=3 has the cyclotomic factors
    s^2 +/- s + 1 (= Phi_6 and Phi_3). These should factor over Q(omega)."""

    def test_s2_minus_s_plus_1_is_phi_6(self):
        p = s**2 - s + 1
        assert sp.simplify(p - sp.cyclotomic_poly(6, s)) == 0

    def test_s2_plus_s_plus_1_is_phi_3(self):
        p = s**2 + s + 1
        assert sp.simplify(p - sp.cyclotomic_poly(3, s)) == 0

    def test_phi_6_splits_over_Q_omega(self):
        """s^2 - s + 1 is irreducible over Q but factors over Q(sqrt -3)."""
        p = s**2 - s + 1
        # Over Q
        fq = sp.factor(p)
        assert _n_poly_factors(fq, s) == 1
        # Over Q(sqrt -3)
        fqo = sp.factor(p, extension=[sp.sqrt(-3)])
        assert _n_poly_factors(fqo, s) == 2

    def test_phi_4_splits_over_Q_i(self):
        """s^2 + 1 = Phi_4 factors over Q(i), not over Q(sqrt 2) or Q(sqrt -3)."""
        p = s**2 + 1
        assert _n_poly_factors(sp.factor(p), s) == 1
        assert _n_poly_factors(sp.factor(p, extension=[sp.I]), s) == 2
        assert _n_poly_factors(sp.factor(p, extension=[sp.sqrt(2)]), s) == 1
        assert _n_poly_factors(sp.factor(p, extension=[sp.sqrt(-3)]), s) == 1


class TestSmallCubicsAreS3:
    """The S^3 Coulomb Ihara zeros' cubic factors 2s^3 +/- s^2 + s +/- 1 and
    the Dirac-A Rule A cubic families all have Galois group S_3 (order 6)."""

    @pytest.mark.parametrize(
        "poly_expr, expected_disc",
        [
            (2 * s**3 - s**2 + s - 1, -83),
            (2 * s**3 + s**2 + s + 1, -83),
            (2 * s**3 - 2 * s**2 + 2 * s - 1, -44),
            (2 * s**3 + 2 * s**2 + 2 * s + 1, -44),
        ],
    )
    def test_cubic_is_S3(self, poly_expr, expected_disc):
        poly = sp.Poly(poly_expr, s)
        disc = sp.discriminant(poly_expr, s)
        assert int(disc) == expected_disc
        # Non-square discriminant => non-abelian Galois => must be S_3 for cubic
        gg, _ = galois_group(poly)
        assert gg.order() == 6
        assert gg.is_solvable
        assert not gg.is_abelian

    def test_cubic_irreducible_over_Qi(self):
        """2s^3 - s^2 + s - 1 should stay irreducible over Q(i) because the
        quadratic resolvent is Q(sqrt -83), not Q(i)."""
        p = 2 * s**3 - s**2 + s - 1
        fq = sp.factor(p)
        fqi = sp.factor(p, extension=[sp.I])
        assert _n_poly_factors(fq, s) == _n_poly_factors(fqi, s) == 1


class TestS5QuadraticFactors:
    """S^5 Bargmann-Segal Ihara zetas contain abelian quadratic factors."""

    def test_2s2_plus_1_over_Q_zeta_8(self):
        """2s^2 + 1 has roots +/- i/sqrt 2; splits over Q(sqrt -2) = Q(i sqrt 2)."""
        p = 2 * s**2 + 1
        # disc = -8 => resolvent Q(sqrt -8) = Q(i sqrt 2) subset of Q(zeta_8)
        disc = sp.discriminant(p, s)
        assert int(disc) == -8
        # Over Q
        assert _n_poly_factors(sp.factor(p), s) == 1
        # Over Q(i): does it split? i * 1/sqrt2 is in Q(zeta_8), not Q(i)
        fqi = sp.factor(p, extension=[sp.I])
        # 2s^2 + 1 = 2(s - i/sqrt2)(s + i/sqrt2). i/sqrt2 requires sqrt(2).
        assert _n_poly_factors(fqi, s) == 1  # stays irreducible over Q(i) alone
        # Over Q(zeta_8) = Q(i, sqrt 2)
        fqz8 = sp.factor(p, extension=[sp.sqrt(2), sp.I])
        assert _n_poly_factors(fqz8, s) == 2


class TestS5N3P12IsSymmetric6:
    """The degree-12 polynomial P_12(s) in the S^5 N_max=3 Ihara zeta reduces
    to a degree-6 polynomial in u = s^2 whose Galois group is S_6 (order 720).
    This is the non-solvable factor of Paper 29."""

    def test_p12_even_in_s(self):
        p12 = (432 * s**12 + 666 * s**10 + 374 * s**8 + 135 * s**6
               + 47 * s**4 + 11 * s**2 + 1)
        # All odd-power coefficients are zero
        coeffs = sp.Poly(p12, s).all_coeffs()
        for i, c in enumerate(coeffs):
            power = 12 - i
            if power % 2 == 1:
                assert c == 0

    def test_p12_u_reduced_Galois_is_S6(self):
        """P_12(s) = P_6(s^2). Gal(P_6 / Q) is S_6, order 720, NON-solvable."""
        u = sp.symbols("u")
        p6 = 432 * u**6 + 666 * u**5 + 374 * u**4 + 135 * u**3 + 47 * u**2 + 11 * u + 1
        # Irreducible over Q
        assert _n_poly_factors(sp.factor(p6), u) == 1
        gg, _ = galois_group(sp.Poly(p6, u))
        assert gg.order() == 720
        assert not gg.is_solvable
        assert not gg.is_abelian


class TestDiracB_n2_saturation_polynomial:
    """Dirac-S^3 Rule B at n_max=2 saturates the Ramanujan bound. The factor
    4s^2 + 1 has roots at s = +/- i/2 exactly on |s| = 1/2 = 1/sqrt(q_max)
    for q_max = 4 — these are boundary Ihara zeros."""

    def test_4s2_plus_1_splits_over_Q_i(self):
        p = 4 * s**2 + 1
        disc = sp.discriminant(p, s)
        assert int(disc) == -16  # square-free part -1 => Q(i) resolvent
        fqi = sp.factor(p, extension=[sp.I])
        assert _n_poly_factors(fqi, s) == 2
        # Roots are +/- i/2
        roots = [sp.Rational(0) + sp.I / 2, sp.Rational(0) - sp.I / 2]
        for r in roots:
            assert p.subs(s, r) == 0
            assert abs(r) == sp.Rational(1, 2)


class TestAlgebraicityEnvelope:
    """Every nontrivial factor has integer coefficients (i.e. every Ihara
    zero is an algebraic integer multiplied by a rational)."""

    @pytest.mark.parametrize(
        "poly_expr",
        [
            s**2 - s + 1, s**2 + s + 1,
            2 * s**3 - s**2 + s - 1, 2 * s**3 + s**2 + s + 1,
            2 * s**2 + 1, 3 * s**4 + 3 * s**2 + 1,
            24 * s**6 + 21 * s**4 + s**2 - 1,
            s**2 + 1, 3 * s**2 + 1, 4 * s**2 + 1,
            2 * s**2 - s + 1, 2 * s**2 + s + 1,
            12 * s**4 + 9 * s**2 - 1,
            9 * s**2 + 1,
        ],
    )
    def test_integer_coefficients(self, poly_expr):
        poly = sp.Poly(poly_expr, s)
        for c in poly.all_coeffs():
            assert c.is_Integer
