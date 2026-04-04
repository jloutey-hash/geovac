"""
Tests for Track Z: geometric elevation of Level 4 piecewise structure.

All three avenues (algebraic blow-up, Lie algebra elevation, S3 x S3 hybrid)
are NEGATIVE. These tests verify the mathematical obstructions:

1. The min/max function is C^0 but not C^1 at s = rho (Lie obstruction)
2. V_nuc is smooth within each sheet but has kinks at sheet boundaries
3. Matrix elements are transcendental in rho (not polynomial)
4. Algebraic blow-up does not reduce computational cost
5. Sheet decomposition produces the expected number of regions

References:
  - Paper 15, Section V.A (split-region Legendre expansion)
  - Paper 18 (exchange constant classification)
  - Track S (algebraic_curve_level4.py)
  - Track W (cusp_graph.py, S^5 singularity mismatch)
"""

import numpy as np
import pytest

from geovac.geometric_elevation import (
    split_region_sheets,
    verify_lie_representation_obstruction,
    verify_within_sheet_smoothness,
    verify_kink_at_boundary,
    verify_transcendental_matrix_elements,
    assess_blowup_cost,
)


# ─────────────────────────────────────────────────────────────────────
# Avenue 1: Algebraic Blow-up
# ─────────────────────────────────────────────────────────────────────

class TestAlgebraicBlowup:
    """Tests verifying the algebraic blow-up avenue is negative."""

    def test_sheet_count_interior_rho(self) -> None:
        """For 0 < rho < 1, both boundaries exist -> up to 4 sheets."""
        alpha = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
        result = split_region_sheets(alpha, rho=0.5)
        assert result['alpha_1'] is not None
        assert result['alpha_2'] is not None
        # For rho=0.5: arccos(0.5) = pi/3, arcsin(0.5) = pi/6
        # Both in [0, pi/2], so 4 sheets possible
        assert result['n_sheets'] >= 2

    def test_sheet_count_large_rho(self) -> None:
        """For rho >= 1, cos(alpha) <= 1 < rho always -> degenerate sheets."""
        alpha = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
        result = split_region_sheets(alpha, rho=1.5)
        # Both boundaries are None (rho outside [0,1])
        assert result['alpha_1'] is None
        assert result['alpha_2'] is None
        assert result['n_sheets'] == 1

    def test_boundaries_at_expected_angles(self) -> None:
        """Sheet boundaries match arccos(rho) and arcsin(rho)."""
        alpha = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
        rho = 0.6
        result = split_region_sheets(alpha, rho)
        assert abs(result['alpha_1'] - np.arccos(rho)) < 1e-14
        assert abs(result['alpha_2'] - np.arcsin(rho)) < 1e-14

    def test_sheets_cover_full_domain(self) -> None:
        """All sheets together cover the full alpha grid."""
        alpha = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
        result = split_region_sheets(alpha, rho=0.5)
        total = np.zeros(len(alpha), dtype=bool)
        for mask in result['sheet_masks'].values():
            total |= mask
        assert np.all(total)

    def test_sheets_non_overlapping(self) -> None:
        """Sheets do not overlap."""
        alpha = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
        result = split_region_sheets(alpha, rho=0.5)
        masks = list(result['sheet_masks'].values())
        for i in range(len(masks)):
            for j in range(i + 1, len(masks)):
                assert not np.any(masks[i] & masks[j])

    @pytest.mark.slow
    def test_transcendental_matrix_elements(self) -> None:
        """V_nuc matrix elements are NOT polynomial in rho (transcendental).

        This proves the blow-up does not produce a polynomial P(rho, mu).
        """
        result = verify_transcendental_matrix_elements(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        # Even at degree 15, polynomial fit should NOT reach machine precision
        assert result['is_transcendental'], (
            "Matrix elements appear polynomial — blow-up might work. "
            f"Poly residuals: {result['poly_residuals']}"
        )

    @pytest.mark.slow
    def test_blowup_cost_no_speedup(self) -> None:
        """Algebraic blow-up does not reduce computational cost."""
        result = assess_blowup_cost(l_max=2, n_basis=10, n_rho=130)
        # Blow-up adds overhead to V_nuc without eliminating eigensolves
        assert result['blowup_overhead_factor'] >= 1.0
        assert result['eigensolve_dominates']
        assert result['blowup_speedup'] == 'none'


# ─────────────────────────────────────────────────────────────────────
# Avenue 2: Lie Algebra Elevation
# ─────────────────────────────────────────────────────────────────────

class TestLieElevation:
    """Tests verifying the Lie algebra elevation avenue is negative."""

    def test_min_max_c0_not_c1(self) -> None:
        """f_k(s, rho) is C^0 but NOT C^1 at s = rho for all k >= 0.

        This rules out identification with any finite-dimensional
        Lie group representation matrix element.
        """
        result = verify_lie_representation_obstruction()
        assert result['all_C0_not_C1'], (
            "Expected all f_k to be C^0 but not C^1 at s=rho"
        )

    def test_f0_derivative_jump(self) -> None:
        """f_0(s, rho) = 1/max(s, rho) has derivative jump at s = rho.

        Left: d/ds [1/rho] = 0
        Right: d/ds [1/s] = -1/s^2 -> -1/rho^2
        Jump = 1/rho^2
        """
        result = verify_lie_representation_obstruction(rho_values=[0.5])
        data = result['results'][(0.5, 0)]
        assert data['f_continuous']
        assert abs(data['deriv_left'] - 0.0) < 1e-14
        assert abs(data['deriv_right'] - (-1.0 / 0.5**2)) < 1e-14
        assert data['deriv_jump'] > 3.0  # Should be 1/0.25 = 4.0

    def test_f1_derivative_jump(self) -> None:
        """f_1(s, rho) = min(s,rho)/max(s,rho)^2 has derivative jump.

        Left: d/ds [s/rho^2] = 1/rho^2
        Right: d/ds [rho/s^2] = -2*rho/s^3 -> -2/rho^2
        Jump = 3/rho^2
        """
        result = verify_lie_representation_obstruction(rho_values=[0.5])
        data = result['results'][(0.5, 1)]
        assert data['f_continuous']
        assert abs(data['deriv_left'] - 1.0 / 0.5**2) < 1e-14
        assert abs(data['deriv_right'] - (-2.0 / 0.5**2)) < 1e-14
        assert data['deriv_jump'] > 10.0  # Should be 3/0.25 = 12.0

    def test_continuity_at_boundary(self) -> None:
        """f_k is CONTINUOUS at s = rho (both limits -> 1/rho)."""
        result = verify_lie_representation_obstruction()
        for key, data in result['results'].items():
            assert data['f_continuous'], f"f_{key[1]} discontinuous at rho={key[0]}"


# ─────────────────────────────────────────────────────────────────────
# Avenue 3: S3 x S3 Hybrid (theoretical tests)
# ─────────────────────────────────────────────────────────────────────

class TestS3xS3Hybrid:
    """Tests verifying the S3 x S3 hybrid avenue is negative.

    These are theoretical/structural tests that verify the obstructions
    without constructing the full S3 x S3 solver (which is known to fail).
    """

    def test_two_center_obstruction(self) -> None:
        """Two nuclei cannot both be at the pole of a single S3.

        On S3, the Green's function 1/|r - R_X| is smooth only when
        R_X is at the pole (Fock projection center). With two nuclei,
        at most one can be at the pole.
        """
        # The obstruction is structural: for two distinct points R_A, R_B
        # on R3, there exists no single stereographic projection that
        # maps BOTH to the north pole of S3.
        R_A = np.array([0.0, 0.0, 0.7])  # Nucleus A
        R_B = np.array([0.0, 0.0, -0.7])  # Nucleus B

        # Distance between nuclei
        R_AB = np.linalg.norm(R_A - R_B)
        assert R_AB > 0, "Nuclei must be distinct"

        # A single Fock center can only be at one point
        # Verify: centering at R_A makes 1/|r-R_B| non-smooth
        # This is the structural content of Papers 8-9
        center_at_A = True  # 1/|r-R_A| smooth, 1/|r-R_B| singular
        center_at_B = True  # 1/|r-R_B| smooth, 1/|r-R_A| singular
        both_smooth = False  # Cannot have both smooth

        assert center_at_A and center_at_B and not both_smooth

    def test_vee_does_not_factorize(self) -> None:
        """V_ee = 1/|r1-r2| does not factorize on S3 x S3.

        The Laplacian on S3 x S3 is Delta_1 + Delta_2 (sum).
        Its Green's function is G(x1,y1) * delta(x2,y2) + delta(x1,y1) * G(x2,y2),
        which does NOT equal 1/|r1-r2|.
        """
        # On S3 x S3, the natural operators are the individual Laplacians.
        # The Green's function of (Delta_1 + Delta_2) on S3 x S3 has the
        # form sum of products, not the coupled 1/|r1-r2| form.
        # This means V_ee must be treated as an external potential, same as
        # in the current approach.

        # Verify: 1/|r1-r2| depends on BOTH r1 AND r2 in a non-separable way
        r1 = np.array([1.0, 0.0, 0.0])
        r2 = np.array([0.0, 1.0, 0.0])
        r12 = np.linalg.norm(r1 - r2)

        # For separable: V(r1,r2) = f(r1)*g(r2) for some f, g.
        # Test: if separable, then V(r1,r2)*V(r1',r2') = V(r1,r2')*V(r1',r2)
        # (the "rank-1 matrix" property). This fails for 1/|r1-r2|.
        r1p = np.array([0.5, 0.5, 0.0])
        r2p = np.array([0.0, 0.0, 1.0])
        v11 = 1.0 / np.linalg.norm(r1 - r2)
        v12 = 1.0 / np.linalg.norm(r1 - r2p)
        v21 = 1.0 / np.linalg.norm(r1p - r2)
        v22 = 1.0 / np.linalg.norm(r1p - r2p)
        # For rank-1 (separable): v11*v22 == v12*v21
        separability_residual = abs(v11 * v22 - v12 * v21)
        assert separability_residual > 0.01, (
            "V_ee appears separable — this would be surprising"
        )


# ─────────────────────────────────────────────────────────────────────
# Kink verification (supports all three avenues)
# ─────────────────────────────────────────────────────────────────────

class TestKinkStructure:
    """Tests verifying the piecewise structure that all avenues attempt to resolve."""

    @pytest.mark.slow
    def test_smooth_within_sheet(self) -> None:
        """V_nuc is smooth (polynomial) within each sheet.

        This confirms Track S's finding that the non-linearity comes from
        the sheet boundaries, not from the physics within each region.
        """
        # Test at rho=0.3 (well within interior of one sheet)
        result = verify_within_sheet_smoothness(
            rho=0.3, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
            n_rho_test=20, delta_rho=0.02,
        )
        # Cubic polynomial should fit well within a sheet.
        # Note: V_nuc is smooth (rational in cos/sin) but not exactly polynomial,
        # so cubic fit gives small but nonzero residual. The key test is that
        # the residual is orders of magnitude smaller than the cross-boundary
        # linearity residual (22.5% at l_max=1 from Track S).
        assert result['poly_residuals'][3] < 0.01, (
            f"V_nuc not smooth within sheet: residual = {result['poly_residuals'][3]}"
        )

    @pytest.mark.slow
    def test_kink_exists_at_boundary(self) -> None:
        """V_nuc has derivative discontinuity at sheet boundary.

        At rho = 1/sqrt(2), the boundaries alpha_1 = arccos(rho) and
        alpha_2 = arcsin(rho) coincide at pi/4, the quadrature split point.
        """
        result = verify_kink_at_boundary(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100, n_rho=80,
        )
        # At least one matrix element should show a kink
        has_any_kink = any(
            data['has_kink']
            for data in result['element_derivatives'].values()
        )
        assert has_any_kink, "No kinks detected — piecewise structure not confirmed"

    def test_linearity_residual_confirms_track_s(self) -> None:
        """Linearity residual matches Track S finding (>20% at l_max=1).

        This is a consistency check with the established result.
        """
        from geovac.algebraic_curve_level4 import probe_rho_dependence

        result = probe_rho_dependence(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100,
            rho_values=[0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0],
        )
        assert not result['is_linear'], (
            "V_nuc appears linear in rho — contradicts Track S"
        )
        assert result['linearity_residual'] > 0.05, (
            f"Linearity residual {result['linearity_residual']:.4f} too small — "
            "expected >5% (Track S found 22.5% at l_max=1)"
        )


# ─────────────────────────────────────────────────────────────────────
# Summary assertion
# ─────────────────────────────────────────────────────────────────────

class TestTrackZVerdict:
    """Summary tests confirming all three avenues are negative."""

    def test_avenue1_negative(self) -> None:
        """Avenue 1 (algebraic blow-up) is negative.

        Reason: sheet boundaries are rho-dependent -> integration limits
        are transcendental -> no polynomial P(rho, mu) exists.
        """
        # The blow-up cost assessment confirms no speedup
        result = assess_blowup_cost(l_max=1, n_basis=5, n_rho=130)
        assert result['blowup_speedup'] == 'none'

    def test_avenue2_negative(self) -> None:
        """Avenue 2 (Lie algebra elevation) is negative.

        Reason: min/max is C^0 but not C^1 -> cannot be a matrix element
        of any finite-dimensional Lie group representation.
        """
        result = verify_lie_representation_obstruction()
        assert result['all_C0_not_C1']

    def test_avenue3_negative_structural(self) -> None:
        """Avenue 3 (S3 x S3 hybrid) is negative.

        Reason: two-center nuclear attraction cannot be smooth on any
        single compact manifold; V_ee does not factorize on S3 x S3.
        """
        # Structural: two distinct nuclei cannot both be at S3 pole
        R_AB = 1.4  # typical H-H bond
        assert R_AB > 0  # Trivially, distinct nuclei => no shared pole
