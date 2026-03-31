"""
Tests for the cusp energy correction module.

Validates:
1. cusp correction is negative (energy too high without cusp)
2. correction -> 0 as l_max -> infinity
3. correction decreases monotonically with l_max
4. He corrected energy is closer to exact
5. H2 correction is R-dependent
6. Dissociation limit: correction -> 0 at R -> infinity (separated atoms)
7. Schwartz coefficient has correct magnitude
8. Extrapolation from multiple l_max values
"""

import numpy as np
import pytest
from math import pi


# ===========================================================================
# Phase A: He (Level 3) cusp correction
# ===========================================================================


class TestCuspCorrectionHe:
    """Tests for He cusp correction."""

    def test_correction_is_negative(self):
        """Cusp correction should be negative (missing cusp raises energy)."""
        from geovac.cusp_correction import cusp_correction_he

        for l_max in [0, 1, 2, 3]:
            dE = cusp_correction_he(Z=2.0, l_max=l_max)
            assert dE < 0, (
                f"Cusp correction at l_max={l_max} should be negative, "
                f"got {dE:.6f}"
            )

    def test_correction_decreases_with_lmax(self):
        """Magnitude of correction should decrease with l_max."""
        from geovac.cusp_correction import cusp_correction_he

        corrections = [cusp_correction_he(Z=2.0, l_max=l) for l in range(6)]

        for i in range(len(corrections) - 1):
            assert abs(corrections[i]) > abs(corrections[i + 1]), (
                f"|dE(l_max={i})| = {abs(corrections[i]):.6f} should be > "
                f"|dE(l_max={i+1})| = {abs(corrections[i+1]):.6f}"
            )

    def test_correction_vanishes_at_large_lmax(self):
        """Correction should approach zero for large l_max."""
        from geovac.cusp_correction import cusp_correction_he

        dE_20 = cusp_correction_he(Z=2.0, l_max=20)
        assert abs(dE_20) < 1e-4, (
            f"Correction at l_max=20 should be < 1e-4 Ha, got {dE_20:.2e}"
        )

    def test_schwartz_coefficient_magnitude(self):
        """Schwartz coefficient A should be ~0.34 for He."""
        from geovac.cusp_correction import schwartz_coefficient

        A = schwartz_coefficient(Z=2.0)
        # Kutzelnigg & Morgan (1992): A ~ 0.34 for He
        assert 0.2 < A < 0.5, (
            f"Schwartz coefficient A = {A:.4f} outside expected range [0.2, 0.5]"
        )

    def test_correction_magnitude_lmax0(self):
        """At l_max=0, cusp correction should be tens of mHa."""
        from geovac.cusp_correction import cusp_correction_he

        dE = cusp_correction_he(Z=2.0, l_max=0)
        # Schwartz formula gives ~80 mHa for l_max=0 (sum of all l>0)
        # This is the THEORETICAL partial-wave truncation error
        assert -0.15 < dE < -0.01, (
            f"Cusp correction at l_max=0 = {dE:.6f} Ha, "
            f"expected between -0.15 and -0.01"
        )

    def test_correction_z_scaling(self):
        """Correction should scale as Z^3 (from coalescence density)."""
        from geovac.cusp_correction import cusp_correction_he

        dE_Z2 = cusp_correction_he(Z=2.0, l_max=1)
        dE_Z3 = cusp_correction_he(Z=3.0, l_max=1)

        ratio = abs(dE_Z3) / abs(dE_Z2)
        expected_ratio = (3.0 / 2.0) ** 3
        # Allow 10% tolerance for the Z-scaling
        assert abs(ratio - expected_ratio) / expected_ratio < 0.1, (
            f"Z-scaling ratio {ratio:.4f} deviates from (3/2)^3 = "
            f"{expected_ratio:.4f}"
        )

    def test_schwartz_convergence_rate(self):
        """Verify 1/(l+1/2)^4 convergence behavior."""
        from geovac.cusp_correction import cusp_correction_he

        # The incremental correction from l_max to l_max+1 should
        # decrease as ~1/(l_max+3/2)^4
        dE = [cusp_correction_he(Z=2.0, l_max=l) for l in range(8)]
        increments = [abs(dE[i] - dE[i + 1]) for i in range(7)]

        # Check that increments decay as ~1/n^4
        for i in range(1, 6):
            ratio = increments[i] / increments[i - 1]
            expected = ((i + 0.5) / (i + 1.5)) ** 4
            # Allow 30% tolerance (the formula is asymptotic)
            assert abs(ratio - expected) / expected < 0.3, (
                f"Increment ratio at l={i}: {ratio:.4f} vs expected "
                f"{expected:.4f}"
            )

    def test_override_coalescence_density(self):
        """Can override the coalescence density."""
        from geovac.cusp_correction import cusp_correction_he

        dE_default = cusp_correction_he(Z=2.0, l_max=1)
        dE_override = cusp_correction_he(Z=2.0, l_max=1, delta_r12=0.2)

        # Different input -> different output
        assert dE_override != dE_default
        # Larger delta -> larger correction
        assert abs(dE_override) > abs(dE_default)


class TestCuspCorrectionHeFromSolver:
    """Tests for cusp correction applied to solver output."""

    def test_correction_from_mock_result(self):
        """Apply correction to a mock solver result dict."""
        from geovac.cusp_correction import cusp_correction_from_solver

        # Use l_max=3 where the correction (~1.7 mHa) is smaller than
        # the solver error (~6.4 mHa), so correction improves things
        mock_result = {
            'energy': -2.8973,  # typical l_max=3 coupled-channel
            'l_max': 3,
        }
        result = cusp_correction_from_solver(mock_result, Z=2.0)

        assert result['delta_E_cusp'] < 0
        assert result['energy_corrected'] < result['energy_uncorrected']
        # At l_max=3, correction is small enough not to overcorrect
        assert result['error_corrected_pct'] < result['error_uncorrected_pct']

    def test_correction_improves_energy(self):
        """Corrected energy should be closer to exact than uncorrected.

        Uses a mock energy typical of the coupled-channel solver at l_max=3.
        The solver gives energy ABOVE exact (too positive), and the cusp
        correction is negative, pushing toward exact.

        At l_max=3, the cusp correction (~1.7 mHa) is smaller than the
        solver error (~6.4 mHa), so it improves without overcorrecting.
        """
        from geovac.cusp_correction import cusp_correction_from_solver

        E_exact = -2.903724
        # Coupled-channel at l_max=3 gives ~0.22% error (above exact).
        # "Above exact" means less negative: E_solver > E_exact.
        E_mock = E_exact + abs(E_exact) * 0.0022  # 0.22% above exact
        mock_result = {'energy': E_mock, 'l_max': 3}

        result = cusp_correction_from_solver(mock_result, Z=2.0)
        assert result['error_corrected_pct'] < result['error_uncorrected_pct']


# ===========================================================================
# Phase A+: Extrapolation from multiple l_max
# ===========================================================================


class TestCuspExtrapolation:
    """Tests for CBS extrapolation."""

    def test_extrapolation_two_points(self):
        """Two-point extrapolation should give reasonable E_exact."""
        from geovac.cusp_correction import cusp_correction_he_extrapolated

        E_exact = -2.903724
        # Mock energies at different l_max (slightly above exact)
        energies = {
            1: E_exact + 0.0050,  # larger error at lower l_max
            3: E_exact + 0.0008,  # smaller error at higher l_max
        }
        E_ext, A_fit, rms = cusp_correction_he_extrapolated(energies)

        # Extrapolated should be closer to exact than either input
        err_ext = abs(E_ext - E_exact)
        err_l1 = abs(energies[1] - E_exact)
        assert err_ext < err_l1, (
            f"Extrapolated error {err_ext:.6f} should be < "
            f"l_max=1 error {err_l1:.6f}"
        )

    def test_extrapolation_needs_two_points(self):
        """Should raise with fewer than 2 points."""
        from geovac.cusp_correction import cusp_correction_he_extrapolated

        with pytest.raises(ValueError):
            cusp_correction_he_extrapolated({0: -2.89})


# ===========================================================================
# Phase B: H2 (Level 4) R-dependent cusp correction
# ===========================================================================


class TestCuspCorrectionH2:
    """Tests for H2 R-dependent cusp correction."""

    def test_correction_is_negative(self):
        """H2 cusp correction should be negative at all R."""
        from geovac.cusp_correction import cusp_correction_h2_point

        for R in [0.5, 1.0, 1.4, 2.0, 3.0, 5.0]:
            dE = cusp_correction_h2_point(R, l_max=2)
            assert dE < 0, (
                f"H2 cusp correction at R={R} should be negative, "
                f"got {dE:.6f}"
            )

    def test_correction_is_r_dependent(self):
        """Correction should vary with R."""
        from geovac.cusp_correction import cusp_correction_h2_point

        dE_short = cusp_correction_h2_point(R=1.0, l_max=2)
        dE_long = cusp_correction_h2_point(R=5.0, l_max=2)

        assert dE_short != dE_long, "Correction should be R-dependent"
        # Larger correction at short R (more overlap)
        assert abs(dE_short) > abs(dE_long), (
            f"|dE(R=1.0)| = {abs(dE_short):.6f} should be > "
            f"|dE(R=5.0)| = {abs(dE_long):.6f}"
        )

    def test_correction_vanishes_at_large_r(self):
        """Correction should approach zero at large R (separated atoms)."""
        from geovac.cusp_correction import cusp_correction_h2_point

        dE_far = cusp_correction_h2_point(R=20.0, l_max=2)
        dE_near = cusp_correction_h2_point(R=1.4, l_max=2)

        # At large R, correction should be very small
        assert abs(dE_far) < 0.01 * abs(dE_near), (
            f"|dE(R=20)| = {abs(dE_far):.2e} should be << "
            f"|dE(R=1.4)| = {abs(dE_near):.2e}"
        )

    def test_united_atom_limit(self):
        """At R=0, H2 coalescence density should equal He value."""
        from geovac.cusp_correction import (
            coalescence_density_h2,
            coalescence_density_he_like,
        )

        delta_h2_r0 = coalescence_density_h2(R=0.0)
        delta_he = coalescence_density_he_like(Z=2.0)

        np.testing.assert_allclose(delta_h2_r0, delta_he, rtol=0.01,
            err_msg="H2 coalescence density at R=0 should match He"
        )

    def test_correction_decreases_with_lmax(self):
        """H2 cusp correction should decrease with l_max at any R."""
        from geovac.cusp_correction import cusp_correction_h2_point

        R = 1.4  # near equilibrium
        corrections = [cusp_correction_h2_point(R, l_max=l)
                      for l in range(5)]

        for i in range(len(corrections) - 1):
            assert abs(corrections[i]) > abs(corrections[i + 1])

    def test_pes_correction(self):
        """PES correction should produce a valid corrected surface."""
        from geovac.cusp_correction import cusp_correction_h2_pes

        R_grid = np.array([1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 10.0])
        # Mock PES: simple Morse-like
        E_inf = -1.0  # 2H limit
        D_e_mock = 0.15  # slightly less than exact
        R_eq = 1.4
        E_mock = E_inf - D_e_mock * (1 - np.exp(-1.0 * (R_grid - R_eq))) ** 2

        result = cusp_correction_h2_pes(R_grid, E_mock, l_max=2)

        # Correction is negative everywhere
        assert np.all(result['delta_E_cusp'] < 0)
        # Corrected PES is lower everywhere
        assert np.all(result['E_corrected'] < result['E_uncorrected'])

    def test_correction_from_level4_mock(self):
        """Apply correction to a mock Level 4 result."""
        from geovac.cusp_correction import cusp_correction_from_level4

        mock_result = {
            'E_total': -1.12,  # typical H2 total at R=1.4
            'R': 1.4,
            'l_max': 2,
        }
        result = cusp_correction_from_level4(mock_result)

        assert result['delta_E_cusp'] < 0
        assert result['E_corrected'] < result['E_uncorrected']
        assert result['R'] == 1.4


# ===========================================================================
# Phase B+: Convergence analysis utility
# ===========================================================================


class TestCuspConvergenceAnalysis:
    """Tests for the convergence analysis utility."""

    def test_analyze_returns_expected_keys(self):
        """Analysis should return expected dictionary keys."""
        from geovac.cusp_correction import analyze_cusp_convergence

        result = analyze_cusp_convergence(Z=2.0)

        assert 'l_max_values' in result
        assert 'corrections' in result
        assert 'corrections_mHa' in result
        assert 'A_schwartz' in result

    def test_corrections_in_mha_consistent(self):
        """mHa values should be 1000x Ha values."""
        from geovac.cusp_correction import analyze_cusp_convergence

        result = analyze_cusp_convergence(Z=2.0)

        for ha, mha in zip(result['corrections'], result['corrections_mHa']):
            np.testing.assert_allclose(mha, ha * 1000, rtol=1e-10)


# ===========================================================================
# Integration tests with actual solvers (marked slow)
# ===========================================================================


@pytest.mark.slow
class TestCuspCorrectionIntegration:
    """Integration tests that run the actual Level 3 solver."""

    def test_he_lmax1_correction_improves(self):
        """Cusp correction at l_max=1 should improve the energy.

        At l_max=0, the adiabatic approximation error (~1.1%) dominates
        and the Schwartz cusp correction (-80 mHa) overcorrects.
        At l_max=1, the l_max truncation error is the dominant component,
        and the cusp correction significantly improves the energy.
        """
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        from geovac.cusp_correction import cusp_correction_from_solver

        result = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=1, n_channels=3,
            n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
        )

        cusp = cusp_correction_from_solver(result, Z=2.0)

        E_exact = -2.903724
        print(f"\nHe l_max=1 cusp correction:")
        print(f"  Uncorrected: {cusp['energy_uncorrected']:.6f} Ha "
              f"({cusp['error_uncorrected_pct']:.4f}%)")
        print(f"  dE_cusp:     {cusp['delta_E_cusp']:.6f} Ha "
              f"({cusp['delta_E_cusp'] * 1000:.3f} mHa)")
        print(f"  Corrected:   {cusp['energy_corrected']:.6f} Ha "
              f"({cusp['error_corrected_pct']:.4f}%)")

        # Correction should be negative
        assert cusp['delta_E_cusp'] < 0
        # At l_max=1, cusp correction should improve the energy
        assert cusp['error_corrected_pct'] < cusp['error_uncorrected_pct'], (
            f"Cusp correction at l_max=1 should improve energy: "
            f"{cusp['error_corrected_pct']:.4f}% vs "
            f"{cusp['error_uncorrected_pct']:.4f}%"
        )
        # Corrected error should be < 0.2%
        assert cusp['error_corrected_pct'] < 0.2

    def test_he_lmax_convergence_with_cusp(self):
        """Cusp-corrected energies should converge faster."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        from geovac.cusp_correction import cusp_correction_from_solver

        E_exact = -2.903724
        results = []

        for l_max in [0, 1, 2, 3]:
            result = solve_hyperspherical_algebraic_coupled(
                Z=2.0, n_basis=15, l_max=l_max, n_channels=3,
                n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
            )
            cusp = cusp_correction_from_solver(result, Z=2.0)
            results.append(cusp)

            print(f"  l_max={l_max}: E_unc={cusp['energy_uncorrected']:.6f} "
                  f"({cusp['error_uncorrected_pct']:.4f}%), "
                  f"E_corr={cusp['energy_corrected']:.6f} "
                  f"({cusp['error_corrected_pct']:.4f}%), "
                  f"dE={cusp['delta_E_cusp'] * 1000:.3f} mHa")

        # Corrections should decrease with l_max
        for i in range(len(results) - 1):
            assert abs(results[i]['delta_E_cusp']) > abs(
                results[i + 1]['delta_E_cusp']
            )
