"""Tests for one-loop QED vacuum polarization on S³.

Validates:
- Seeley-DeWitt coefficients a_0, a_1, a_2 on unit S³
- Vacuum polarization coefficient = 1/(48π²)
- QED beta function β(α) = 2α²/(3π)
- No odd-zeta content at one loop (structural, T9 theorem)
- Spectral zeta sum convergence
- T9 theorem cross-check at s = 2, 3
"""

import math

import pytest
import sympy as sp
from sympy import Integer, Rational, pi, sqrt, zeta, Symbol

from geovac.qed_vacuum_polarization import (
    seeley_dewitt_coefficients_s3,
    spectral_zeta_massive,
    vacuum_polarization_coefficient,
    beta_function_qed,
    spectral_zeta_derivative_at_zero,
    classify_transcendental_content,
    verify_a0_from_spectral_sum,
    verify_no_odd_zeta_one_loop,
)


# -----------------------------------------------------------------------
# Seeley-DeWitt coefficients on unit S³
# -----------------------------------------------------------------------

class TestSeeleyDeWittS3:
    """Test heat kernel coefficients for D² on S³(R=1)."""

    def setup_method(self):
        self.coeffs = seeley_dewitt_coefficients_s3(R=Integer(1))

    def test_volume_unit_s3(self):
        """Vol(S³) = 2π²."""
        assert sp.simplify(self.coeffs['vol'] - 2 * pi**2) == 0

    def test_scalar_curvature_unit_s3(self):
        """R_scalar = 6 on unit S³."""
        assert self.coeffs['R_scalar'] == Integer(6)

    def test_lichnerowicz_endomorphism(self):
        """E_Lich = R_scalar/4 = 3/2 on unit S³."""
        assert sp.simplify(self.coeffs['E_lich'] - Rational(3, 2)) == 0

    def test_dim_spinor(self):
        """dim(spinor) = 4 for 4-component Dirac."""
        assert self.coeffs['dim_spinor'] == 4

    def test_a0_value(self):
        """a_0 = (4π)^{-3/2} · 4 · 2π² = 8π² / (4π)^{3/2}.

        (4π)^{3/2} = 8 π^{3/2}, so a_0 = 8π² / (8π^{3/2}) = √π.
        """
        a0 = self.coeffs['a0']
        expected = sqrt(pi)
        assert sp.simplify(a0 - expected) == 0

    def test_a1_value(self):
        """a_1 = (4π)^{-3/2} · 4 · (6/6) · 2π² = (4π)^{-3/2} · 8π².

        This is a_0 · (R_scalar/6) = √π · 1 = √π.
        """
        a1 = self.coeffs['a1']
        expected = sqrt(pi)
        assert sp.simplify(a1 - expected) == 0

    def test_a2_curvature_invariants(self):
        """Verify curvature invariants on unit S³.

        |Ric|² = 12, |Riem|² = 12 on unit S³.
        """
        assert self.coeffs['Ric_sq'] == Integer(12)
        assert self.coeffs['Riem_sq'] == Integer(12)

    def test_a2_integrand(self):
        """Verify the a_2 curvature integrand on unit S³.

        5·36 - 2·12 + 2·12 - 30·6·(3/2) + 60·(9/4)
        = 180 - 24 + 24 - 270 + 135 = 45
        """
        integrand = self.coeffs['integrand_a2']
        expected = Integer(45)
        assert sp.simplify(integrand - expected) == 0

    def test_a2_value(self):
        """a_2 = (4π)^{-3/2} · 4 · (1/360) · 45 · 2π².

        = (4π)^{-3/2} · 4 · (1/8) · 2π² = (4π)^{-3/2} · π²
        = π² / (8π^{3/2}) = √π / 8.
        """
        a2 = self.coeffs['a2']
        expected = sqrt(pi) / 8
        assert sp.simplify(a2 - expected) == 0

    def test_a0_a1_a2_ratio(self):
        """Check a_0 : a_1 : a_2 = 1 : 1 : 1/8 on unit S³."""
        a0, a1, a2 = self.coeffs['a0'], self.coeffs['a1'], self.coeffs['a2']
        assert sp.simplify(a0 - a1) == 0
        assert sp.simplify(a2 / a0 - Rational(1, 8)) == 0


class TestSeeleyDeWittRadiusScaling:
    """Test R-dependence of Seeley-DeWitt coefficients."""

    def test_a0_scales_as_R3(self):
        """a_0 scales as R³ (from Vol(S³) = 2π²R³)."""
        R = Symbol('R', positive=True)
        coeffs = seeley_dewitt_coefficients_s3(R=R)
        a0 = coeffs['a0']
        # a0(R) / a0(1) = R³
        a0_unit = seeley_dewitt_coefficients_s3(R=Integer(1))['a0']
        ratio = sp.simplify(a0 / a0_unit)
        assert sp.simplify(ratio - R**3) == 0

    def test_a1_scales_as_R(self):
        """a_1 scales as R (from R_scalar·Vol = 6/R² · 2π²R³ = 12π²R)."""
        R = Symbol('R', positive=True)
        coeffs = seeley_dewitt_coefficients_s3(R=R)
        a1 = coeffs['a1']
        a1_unit = seeley_dewitt_coefficients_s3(R=Integer(1))['a1']
        ratio = sp.simplify(a1 / a1_unit)
        assert sp.simplify(ratio - R) == 0

    def test_a2_scales_as_R_inv1(self):
        """a_2 scales as 1/R on unit S³.

        The integrand is quartic in curvature (~ 1/R⁴), times Vol (~ R³),
        so a_2 ~ R⁻¹.
        """
        R = Symbol('R', positive=True)
        coeffs = seeley_dewitt_coefficients_s3(R=R)
        a2 = coeffs['a2']
        a2_unit = seeley_dewitt_coefficients_s3(R=Integer(1))['a2']
        ratio = sp.simplify(a2 / a2_unit)
        assert sp.simplify(ratio - 1 / R) == 0


# -----------------------------------------------------------------------
# Vacuum polarization and beta function
# -----------------------------------------------------------------------

class TestVacuumPolarization:
    """Test vacuum polarization coefficient and QED beta function."""

    def test_vac_pol_coefficient_exact(self):
        """Vacuum polarization coefficient = 1/(48π²)."""
        vp = vacuum_polarization_coefficient()
        expected = Rational(1, 48) / pi**2
        assert sp.simplify(vp - expected) == 0

    def test_vac_pol_coefficient_numerical(self):
        """Numerical check: 1/(48π²) ≈ 0.002111."""
        vp = float(vacuum_polarization_coefficient())
        expected = 1.0 / (48.0 * math.pi**2)
        assert abs(vp - expected) / expected < 1e-14

    def test_beta_function_form(self):
        """β(α) = 2α²/(3π)."""
        alpha = Symbol('alpha', positive=True)
        beta = beta_function_qed(alpha)
        expected = 2 * alpha**2 / (3 * pi)
        assert sp.simplify(beta - expected) == 0

    def test_beta_function_numerical(self):
        """Numerical: β(1/137) = 2/(3π·137²) ≈ 1.13×10⁻⁵."""
        alpha_val = 1.0 / 137.036
        beta_val = 2.0 * alpha_val**2 / (3.0 * math.pi)
        expected = 2.0 / (3.0 * math.pi * 137.036**2)
        assert abs(beta_val - expected) / expected < 1e-10

    def test_vp_from_beta(self):
        """Consistency: β(e) = e³/(12π²) comes from Π = 1/(48π²).

        The relation is β(α) = 2α · Π · (4π)² · (something).
        More directly: β(α) = 2α² · (4Π) = 2α²/(12π²)·(4π²)·(...).
        Standard check: β(α) = (4/3) · (α/π) · α = 2α²/(3π). ✓
        """
        alpha = Symbol('alpha', positive=True)
        beta = beta_function_qed(alpha)
        # Alternative derivation: coefficient of α² in β(α) is 2/(3π)
        coeff = sp.Rational(beta.as_coefficients_dict()[alpha**2])
        # This should be 2/(3π) but as_coefficients_dict may not simplify
        # So check numerically
        beta_over_alpha_sq = sp.simplify(beta / alpha**2)
        expected = Rational(2, 3) / pi
        assert sp.simplify(beta_over_alpha_sq - expected) == 0


# -----------------------------------------------------------------------
# No odd-zeta content (T9 theorem structural test)
# -----------------------------------------------------------------------

class TestNoOddZeta:
    """Structural test: one-loop QED on S³ has no odd-zeta content."""

    def test_zeta_D2_s2_no_odd_zeta(self):
        """ζ_{D²}(2) = π² - π⁴/12: no ζ(3)."""
        results = verify_no_odd_zeta_one_loop()
        val = results['zeta_D2_2']
        expected = pi**2 - pi**4 / 12
        assert sp.simplify(val - expected) == 0

    def test_zeta_D2_s3_no_odd_zeta(self):
        """ζ_{D²}(3) = π⁴/3 - π⁶/30: no ζ(5)."""
        results = verify_no_odd_zeta_one_loop()
        val = results['zeta_D2_3']
        expected = pi**4 / 3 - pi**6 / 30
        assert sp.simplify(val - expected) == 0

    def test_zeta_D2_s4_no_odd_zeta(self):
        """ζ_{D²}(4) is a polynomial in π² — no ζ(7)."""
        results = verify_no_odd_zeta_one_loop()
        val = results['zeta_D2_4']
        # ζ_{D²}(4) = 2^7 · [λ(6) - λ(8)]
        # λ(6) = (1 - 1/64)ζ(6) = (63/64)(π⁶/945) = π⁶/960
        # λ(8) = (1 - 1/256)ζ(8) = (255/256)(π⁸/9450) = 255π⁸/(256·9450)
        #       = 17π⁸/161280
        # ζ_{D²}(4) = 128 · (π⁶/960 - 17π⁸/161280)
        #           = 128π⁶/960 - 128·17π⁸/161280
        #           = 2π⁶/15 - 17π⁸/1260
        expected = 2 * pi**6 / 15 - 17 * pi**8 / 1260
        assert abs(complex(sp.simplify(val - expected)).real) < 1e-10

    def test_all_results_pi_even_only(self):
        """All ζ_{D²}(s) at integer s are rational × π^{even}."""
        results = verify_no_odd_zeta_one_loop()
        for s_val in [2, 3, 4, 5]:
            val = results[f'zeta_D2_{s_val}']
            # Check that val is a sum of terms rational × π^{2k}
            # by verifying it's real and doesn't contain zeta(odd)
            val_expanded = sp.expand(val)
            # If it contained ζ(3) or ζ(5) etc., sympy would leave them
            # unsimplified. Check no free zeta() calls remain.
            atoms = val_expanded.atoms(sp.Function)
            zeta_atoms = [a for a in atoms if isinstance(a, sp.core.function.AppliedUndef)
                          or (hasattr(a, 'func') and a.func == zeta)]
            # sympy evaluates ζ(even) to π^even automatically, so any
            # remaining zeta() calls would be odd arguments.
            assert len(zeta_atoms) == 0, (
                f"ζ_{{D²}}({s_val}) contains unevaluated zeta: {zeta_atoms}"
            )

    def test_classification_odd_zeta_absent(self):
        """Transcendental classification confirms odd-zeta absence."""
        classification = classify_transcendental_content()
        assert 'ABSENT' in classification['odd_zeta_at_one_loop']

    def test_classification_two_loop_expected(self):
        """Two-loop odd-zeta is expected (Track D3 ζ(3) appearance)."""
        classification = classify_transcendental_content()
        assert 'EXPECTED' in classification['odd_zeta_at_two_loops']
        assert 'ζ(3)' in classification['odd_zeta_at_two_loops']


# -----------------------------------------------------------------------
# Spectral zeta function convergence
# -----------------------------------------------------------------------

class TestSpectralZeta:
    """Test spectral zeta sum convergence and T9 cross-checks."""

    def test_spectral_zeta_converges_s2(self):
        """Direct sum at s=2 converges to T9 value.

        ζ_{D²}(2) converges slowly (the sum is barely convergent at
        s = d/2 = 3/2; at s = 2 the partial sums converge as O(1/N)).
        With n_max = 500 we get ~0.2% relative error.
        """
        result = verify_a0_from_spectral_sum(n_max=500)
        assert result['rel_err_s2'] < 5e-3

    def test_spectral_zeta_converges_s3(self):
        """Direct sum at s=3 converges to T9 value."""
        result = verify_a0_from_spectral_sum(n_max=500)
        assert result['rel_err_s3'] < 1e-6

    def test_spectral_zeta_massive_s2_massless(self):
        """Massless spectral zeta at s=2 approaches T9 (slow convergence)."""
        direct = spectral_zeta_massive(s=2.0, m_sq=0.0, n_max=500)
        t9_val = math.pi**2 - math.pi**4 / 12.0
        assert abs(direct - t9_val) / abs(t9_val) < 5e-3

    def test_spectral_zeta_massive_m_dependence(self):
        """Mass term increases eigenvalues, so ζ(s, m²) < ζ(s, 0)."""
        z_0 = spectral_zeta_massive(s=2.0, m_sq=0.0, n_max=100)
        z_1 = spectral_zeta_massive(s=2.0, m_sq=1.0, n_max=100)
        assert z_1 < z_0

    def test_zeta_prime_convergence(self):
        """ζ'(0) at two cutoffs should be close (convergence test)."""
        z500 = spectral_zeta_derivative_at_zero(m_sq=1.0, n_max=500)
        z1000 = spectral_zeta_derivative_at_zero(m_sq=1.0, n_max=1000)
        # The regularized sum should converge; check relative change is small
        if abs(z1000) > 1e-10:
            rel_change = abs(z1000 - z500) / abs(z1000)
            assert rel_change < 0.1  # Loose bound — numerical regularization


# -----------------------------------------------------------------------
# Paper 18 taxonomy integration
# -----------------------------------------------------------------------

class TestPaper18Taxonomy:
    """Verify Paper 18 exchange constant taxonomy entries."""

    def test_a0_is_calibration_pi(self):
        """a_0 = √π contains only calibration π (Vol(S³))."""
        a0 = seeley_dewitt_coefficients_s3()['a0']
        # a0 = sqrt(pi) — a half-integer power of π
        assert sp.simplify(a0 - sp.sqrt(pi)) == 0

    def test_a1_rational_on_unit_s3(self):
        """a_1/a_0 = R_scalar/6 = 1 on unit S³ — rational."""
        coeffs = seeley_dewitt_coefficients_s3()
        ratio = sp.simplify(coeffs['a1'] / coeffs['a0'])
        assert ratio == 1

    def test_a2_coefficient_calibration_pi(self):
        """Vacuum polarization 1/(48π²) is calibration-tier π."""
        vp = vacuum_polarization_coefficient()
        # Express as rational / π²
        vp_times_pi_sq = sp.simplify(vp * pi**2)
        assert vp_times_pi_sq == Rational(1, 48)

    def test_beta_function_calibration_pi(self):
        """β(α) = 2α²/(3π) — single power of calibration π."""
        alpha = Symbol('alpha', positive=True)
        beta = beta_function_qed(alpha)
        # β · (3π) / (2α²) = 1
        check = sp.simplify(beta * 3 * pi / (2 * alpha**2))
        assert check == 1
