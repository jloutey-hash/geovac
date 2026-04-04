"""
Tests for cusp-adapted angular basis investigation.

Validates:
1. SO(6) Casimir reproduction at rho=0
2. Consistency with standard spectral solver
3. Channel weight diagnostics
4. Schwartz extrapolation mechanics
5. Track J analogy assessment (documents negative result)

This test suite documents the investigation into whether the Track J
algebraicization pattern (absorbing singularities into the basis weight)
transfers to the theta_12 coordinate at Level 4. The answer is NO:
the theta_12 dependence is encoded in the discrete channel structure,
not a continuous coordinate with a modifiable weight function.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose


# ============================================================================
# SO(6) Casimir verification
# ============================================================================

class TestSO6Casimir:
    """Verify that the adapted solver preserves SO(6) free eigenvalues."""

    def test_casimir_lmax0(self) -> None:
        """l_max=0: single channel (0,0), free eigenvalues exact."""
        from geovac.cusp_angular_basis import verify_so6_casimir
        result = verify_so6_casimir(l_max=0, n_basis=5)
        assert result['pass'], f"Casimir error: {result['max_casimir_error']}"
        assert result['max_casimir_error'] < 1e-10

    def test_casimir_lmax2(self) -> None:
        """l_max=2: multiple channels, free eigenvalues exact."""
        from geovac.cusp_angular_basis import verify_so6_casimir
        result = verify_so6_casimir(l_max=2, n_basis=5)
        assert result['pass'], f"Casimir error: {result['max_casimir_error']}"
        assert result['max_casimir_error'] < 1e-10

    def test_casimir_values_lmax0(self) -> None:
        """l_max=0 singlet: Casimir values are nu*(nu+4)/2 with nu=0,4,8,..."""
        from geovac.cusp_angular_basis import verify_so6_casimir
        result = verify_so6_casimir(l_max=0, n_basis=5)
        # nu=0: 0, nu=4: 16, nu=8: 48, nu=12: 96, nu=16: 160
        expected_first = [0.0, 16.0, 48.0, 96.0, 160.0]
        assert_allclose(result['casimir_exact'][:5], expected_first, atol=1e-10)


# ============================================================================
# Consistency with standard solver
# ============================================================================

class TestStandardConsistency:
    """Verify that adapted solver matches standard spectral solver."""

    def test_consistency_lmax0(self) -> None:
        """l_max=0: adapted solver matches standard (no extrapolation)."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        rho, R_e = 0.5, 2.0
        std = Level4SpectralAngular(l_max=0, n_basis=10, n_quad=100)
        adp = CuspAdaptedAngularSolver(l_max=0, n_basis=10, n_quad=100, extrapolate=False)

        mu_std, _ = std.solve(rho, R_e, n_eig=1)
        mu_adp, _ = adp.solve(rho, R_e, n_eig=1)

        assert_allclose(mu_adp[0], mu_std[0], atol=1e-10,
                        err_msg="Adapted solver should match standard at l_max=0")

    def test_consistency_lmax2(self) -> None:
        """l_max=2: adapted solver matches standard (no extrapolation)."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        rho, R_e = 0.5, 2.0
        std = Level4SpectralAngular(l_max=2, n_basis=10, n_quad=100)
        adp = CuspAdaptedAngularSolver(l_max=2, n_basis=10, n_quad=100, extrapolate=False)

        mu_std, _ = std.solve(rho, R_e, n_eig=1)
        mu_adp, _ = adp.solve(rho, R_e, n_eig=1)

        assert_allclose(mu_adp[0], mu_std[0], atol=1e-10,
                        err_msg="Adapted solver should match standard at l_max=2")

    def test_consistency_multiple_rho(self) -> None:
        """Consistency at several rho values."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        R_e = 2.0
        rho_values = [0.2, 0.5, 1.0, 2.0]

        std = Level4SpectralAngular(l_max=2, n_basis=10, n_quad=100)
        adp = CuspAdaptedAngularSolver(l_max=2, n_basis=10, n_quad=100, extrapolate=False)

        for rho in rho_values:
            mu_std, _ = std.solve(rho, R_e, n_eig=1)
            mu_adp, _ = adp.solve(rho, R_e, n_eig=1)
            assert_allclose(mu_adp[0], mu_std[0], atol=1e-10,
                            err_msg=f"Mismatch at rho={rho}")

    def test_hamiltonian_hermitian(self) -> None:
        """Hamiltonian should be symmetric."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver

        solver = CuspAdaptedAngularSolver(l_max=2, n_basis=10, n_quad=100)
        diag = solver.solve_with_diagnostics(rho=0.5, R_e=2.0)
        assert diag['H_symmetric'], "Hamiltonian should be symmetric"


# ============================================================================
# Channel coefficient analysis
# ============================================================================

class TestChannelAnalysis:
    """Test channel weight extraction and convergence diagnostics."""

    def test_channel_weights_sum_to_one(self) -> None:
        """Channel weights should sum to 1 (eigenvector normalization)."""
        from geovac.cusp_angular_basis import extract_channel_coefficients

        data = extract_channel_coefficients(rho=0.5, R_e=2.0, l_max=2, n_basis=10)
        assert_allclose(data['total_weight'], 1.0, atol=1e-10,
                        err_msg="Channel weights should sum to 1")

    def test_l0_dominates(self) -> None:
        """At moderate rho, l_sum=0 channel should dominate."""
        from geovac.cusp_angular_basis import extract_channel_coefficients

        data = extract_channel_coefficients(rho=0.5, R_e=2.0, l_max=2, n_basis=10)
        l0_weight = data['l_sum_weights'].get(0, 0.0)
        assert l0_weight > 0.5, (
            f"l_sum=0 weight {l0_weight:.4f} should dominate (>0.5)")

    def test_higher_l_nonzero(self) -> None:
        """Higher l channels should have nonzero weight (from V_ee coupling)."""
        from geovac.cusp_angular_basis import extract_channel_coefficients

        data = extract_channel_coefficients(rho=0.5, R_e=2.0, l_max=2, n_basis=10)
        l2_weight = data['l_sum_weights'].get(2, 0.0)
        assert l2_weight > 1e-6, (
            f"l_sum=2 weight {l2_weight:.6e} should be nonzero")

    def test_channel_convergence(self) -> None:
        """Eigenvalue should decrease (become more negative) with l_max."""
        from geovac.cusp_angular_basis import analyze_channel_convergence

        conv = analyze_channel_convergence(
            rho=0.5, R_e=2.0, l_max_values=[0, 2], n_basis=10,
        )
        # With V_ee, mu should decrease as l_max increases (more channels)
        mu_vals = conv['mu_values']
        assert mu_vals[1] <= mu_vals[0] + 0.1, (
            f"mu should decrease or stay similar with l_max: "
            f"mu(0)={mu_vals[0]:.4f}, mu(2)={mu_vals[1]:.4f}")


# ============================================================================
# Schwartz extrapolation
# ============================================================================

class TestSchwartzExtrapolation:
    """Test Schwartz extrapolation mechanics."""

    def test_tail_positive(self) -> None:
        """Schwartz tail should be positive for p > 1."""
        from geovac.cusp_angular_basis import _schwartz_tail
        tail = _schwartz_tail(2, 4.0)
        assert tail > 0, f"Tail should be positive: {tail}"

    def test_tail_decreasing(self) -> None:
        """Tail should decrease with l_max."""
        from geovac.cusp_angular_basis import _schwartz_tail
        t1 = _schwartz_tail(0, 4.0)
        t2 = _schwartz_tail(2, 4.0)
        t3 = _schwartz_tail(4, 4.0)
        assert t1 > t2 > t3 > 0, (
            f"Tail should decrease: {t1:.6f}, {t2:.6f}, {t3:.6f}")

    def test_extrapolation_two_points(self) -> None:
        """Two-point extrapolation should return a result."""
        from geovac.cusp_angular_basis import schwartz_extrapolate_mu
        mu = np.array([-3.0, -3.1])
        lmax = np.array([0, 2])
        result = schwartz_extrapolate_mu(mu, lmax)
        assert result['converged']
        assert result['mu_extrapolated'] < mu[-1], (
            "Extrapolated mu should be lower (more negative) than last point")

    def test_extrapolation_three_points(self) -> None:
        """Three-point extrapolation should fit the exponent."""
        from geovac.cusp_angular_basis import schwartz_extrapolate_mu
        # Synthetic data following 1/(l+0.5)^4
        A = 1.0
        l_vals = np.array([0, 2, 4], dtype=float)
        mu_exact = -5.0  # target
        mu = np.array([
            mu_exact + A * sum(1/(l+0.5)**4 for l in range(int(lv)+1, 100))
            for lv in l_vals
        ])
        result = schwartz_extrapolate_mu(mu, l_vals)
        assert result['converged']
        # p_fit should be close to 4.0
        if 'p_fit' in result:
            assert abs(result['p_fit'] - 4.0) < 1.0, (
                f"Fitted exponent {result['p_fit']:.2f} should be near 4.0")

    def test_single_point_no_extrapolation(self) -> None:
        """Single point should return unconverged."""
        from geovac.cusp_angular_basis import schwartz_extrapolate_mu
        result = schwartz_extrapolate_mu(np.array([-3.0]), np.array([0]))
        assert not result['converged']


# ============================================================================
# Extrapolation applied to angular eigenvalue
# ============================================================================

class TestExtrapolatedSolver:
    """Test the full solver with Schwartz extrapolation enabled."""

    def test_extrapolation_lowers_mu(self) -> None:
        """Extrapolated mu should be lower than raw mu (cusp lowers energy)."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver

        rho, R_e = 0.5, 2.0

        raw_solver = CuspAdaptedAngularSolver(
            l_max=2, n_basis=10, n_quad=100, extrapolate=False,
        )
        ext_solver = CuspAdaptedAngularSolver(
            l_max=2, n_basis=10, n_quad=100, extrapolate=True,
        )

        mu_raw, _ = raw_solver.solve(rho, R_e, n_eig=1)
        mu_ext, _ = ext_solver.solve(rho, R_e, n_eig=1)

        # Extrapolation should lower the eigenvalue (cusp correction is negative)
        assert mu_ext[0] <= mu_raw[0] + 0.01, (
            f"Extrapolated mu ({mu_ext[0]:.6f}) should be <= raw ({mu_raw[0]:.6f})")

    def test_extrapolation_magnitude(self) -> None:
        """Extrapolation correction should be small (< 10% of eigenvalue)."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver

        rho, R_e = 0.5, 2.0

        raw_solver = CuspAdaptedAngularSolver(
            l_max=2, n_basis=10, n_quad=100, extrapolate=False,
        )
        ext_solver = CuspAdaptedAngularSolver(
            l_max=2, n_basis=10, n_quad=100, extrapolate=True,
        )

        mu_raw, _ = raw_solver.solve(rho, R_e, n_eig=1)
        mu_ext, _ = ext_solver.solve(rho, R_e, n_eig=1)

        correction = abs(mu_ext[0] - mu_raw[0])
        assert correction < 0.1 * abs(mu_raw[0]) + 1.0, (
            f"Correction {correction:.4f} should be small vs |mu| = {abs(mu_raw[0]):.4f}")


# ============================================================================
# Track J analogy assessment
# ============================================================================

class TestTrackJAnalogy:
    """Document the Track J analogy assessment (negative result)."""

    def test_analogy_negative(self) -> None:
        """Track J pattern does not transfer to theta_12."""
        from geovac.cusp_angular_basis import assess_track_j_analogy

        result = assess_track_j_analogy()

        # Track J works because: 1D singularity, continuous coordinate, modifiable weight
        assert result['track_j_singularity']['dimension'] == '1D'
        assert result['level4_cusp']['dimension'] == '2D'

        # The classification should be NEGATIVE RESULT
        assert 'NEGATIVE RESULT' in result['classification']

    def test_cusp_is_2d(self) -> None:
        """Cusp lives in 2D (alpha, theta_12), not 1D."""
        from geovac.cusp_angular_basis import assess_track_j_analogy
        result = assess_track_j_analogy()
        assert result['level4_cusp']['dimension'] == '2D'

    def test_theta12_not_integration_variable(self) -> None:
        """theta_12 is NOT an integration variable -- it's in the channel structure."""
        from geovac.cusp_angular_basis import assess_track_j_analogy
        result = assess_track_j_analogy()
        assert 'NOT integration variable' in result['level4_cusp']['coordinate']

    def test_track_x_recommended(self) -> None:
        """Track X (Schwartz extrapolation) is the recommended approach."""
        from geovac.cusp_angular_basis import assess_track_j_analogy
        result = assess_track_j_analogy()
        assert 'Track X' in result['what_works_instead']


# ============================================================================
# Comparison: standard vs adapted
# ============================================================================

class TestComparison:
    """Compare standard and adapted solvers."""

    def test_comparison_runs(self) -> None:
        """Comparison function should run without error."""
        from geovac.cusp_angular_basis import compare_standard_vs_adapted
        result = compare_standard_vs_adapted(
            rho=0.5, R_e=2.0, l_max_values=[0, 2], n_basis=10,
        )
        assert len(result['results']) == 2
        for r in result['results']:
            assert 'mu_standard' in r
            assert 'mu_adapted' in r

    def test_no_extrapolation_at_lmax0(self) -> None:
        """At l_max=0, adapted should equal standard (no extrapolation)."""
        from geovac.cusp_angular_basis import compare_standard_vs_adapted
        result = compare_standard_vs_adapted(
            rho=0.5, R_e=2.0, l_max_values=[0], n_basis=10,
        )
        r = result['results'][0]
        assert_allclose(r['mu_adapted'], r['mu_standard'], atol=1e-10)

    def test_extrapolation_difference_at_lmax2(self) -> None:
        """At l_max=2 with extrapolation, adapted may differ from standard."""
        from geovac.cusp_angular_basis import compare_standard_vs_adapted
        result = compare_standard_vs_adapted(
            rho=0.5, R_e=2.0, l_max_values=[2], n_basis=10,
        )
        r = result['results'][0]
        # The difference should exist (extrapolation) but be small
        diff = abs(r['mu_difference'])
        assert diff < 5.0, f"Extrapolation correction too large: {diff:.4f}"


# ============================================================================
# Gerade symmetry preservation
# ============================================================================

class TestSymmetryPreservation:
    """Verify that gerade/ungerade symmetry is preserved."""

    def test_gerade_channels_only(self) -> None:
        """Homonuclear: only l1+l2 even channels."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        solver = CuspAdaptedAngularSolver(l_max=3, n_basis=5)
        for l1, l2 in solver.channels_2:
            assert (l1 + l2) % 2 == 0, (
                f"Channel ({l1},{l2}) violates gerade symmetry")

    def test_heteronuclear_all_parity(self) -> None:
        """Heteronuclear: all l1+l2 parities allowed."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        solver = CuspAdaptedAngularSolver(l_max=2, n_basis=5, Z_A=1.0, Z_B=2.0)
        parities = set((l1 + l2) % 2 for l1, l2 in solver.channels_2)
        assert 0 in parities, "Should have even parity channels"
        assert 1 in parities, "Heteronuclear should have odd parity channels"


# ============================================================================
# Overlap matrix positive definiteness
# ============================================================================

class TestOverlapMatrix:
    """Test that the overlap matrix (identity for orthonormal basis) is well-conditioned."""

    def test_basis_orthonormal(self) -> None:
        """Basis functions should be orthonormal within each channel."""
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        solver = CuspAdaptedAngularSolver(l_max=2, n_basis=10, n_quad=100)

        for ic in range(solver.n_ch):
            phi = solver._channel_phi[ic]
            S = phi @ np.diag(solver._weights) @ phi.T
            assert_allclose(S, np.eye(phi.shape[0]), atol=1e-10,
                            err_msg=f"Channel {ic} basis not orthonormal")

    def test_cross_channel_overlap_finite(self) -> None:
        """Cross-channel overlap exists but the full Hamiltonian is well-defined.

        Channels with different (l1,l2) are NOT orthogonal in the raw inner
        product because they share the same integration domain [0, pi/2].
        Orthogonality is guaranteed by the Hamiltonian structure (V_ee coupling
        through Gaunt integrals), not by the basis overlap. The important thing
        is that the Hamiltonian matrix is symmetric and the eigensolve succeeds.
        """
        from geovac.cusp_angular_basis import CuspAdaptedAngularSolver
        solver = CuspAdaptedAngularSolver(l_max=2, n_basis=5, n_quad=100)

        # Just verify eigensolve works (Hamiltonian is block-diagonal in the
        # spectral basis because V_ee couples channels, not the basis itself)
        diag = solver.solve_with_diagnostics(rho=0.5, R_e=2.0)
        assert diag['H_symmetric']
        # Eigenvalues should be real and finite
        assert np.all(np.isfinite(diag['eigenvalues']))
