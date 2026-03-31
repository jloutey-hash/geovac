"""
Tests for Track W: Cusp-Graph Theory Investigation

Validates the mathematical argument that 1/r_12 CANNOT be absorbed into
a conformally weighted graph Laplacian on the Level 4 angular space.
"""

import numpy as np
import pytest
from geovac.cusp_graph import (
    greens_function_singularity_order,
    coulomb_singularity_order,
    check_singularity_match,
    coalescence_codimension,
    fiber_green_analysis,
    conformal_singularity_preservation,
    verify_r12_singularity_at_coalescence,
    verify_s3_green_matches_coulomb,
    full_cusp_graph_analysis,
)


# ============================================================================
# Green's function singularity order
# ============================================================================

class TestGreensFunction:
    """Test Green's function singularity orders on spheres."""

    def test_s1_logarithmic(self):
        """S^1: Green's function is logarithmic, not power-law."""
        assert greens_function_singularity_order(1) == -1

    def test_s2_logarithmic(self):
        """S^2: Green's function is logarithmic, not power-law."""
        assert greens_function_singularity_order(2) == -1

    def test_s3_order_1(self):
        """S^3: Green ~ 1/d^1 — the unique sphere matching Coulomb 1/r."""
        assert greens_function_singularity_order(3) == 1

    def test_s4_order_2(self):
        """S^4: Green ~ 1/d^2."""
        assert greens_function_singularity_order(4) == 2

    def test_s5_order_3(self):
        """S^5: Green ~ 1/d^3 — the Level 4 angular space."""
        assert greens_function_singularity_order(5) == 3

    def test_general_formula(self):
        """Green on S^n has singularity 1/d^{n-2} for n >= 3."""
        for n in range(3, 20):
            assert greens_function_singularity_order(n) == n - 2


class TestSingularityMatch:
    """Test which spheres match the Coulomb singularity."""

    def test_only_s3_matches(self):
        """S^3 is the UNIQUE sphere where Green matches Coulomb 1/r."""
        for n in range(1, 20):
            result = check_singularity_match(n)
            if n == 3:
                assert result['match'] is True, f"S^{n} should match"
            else:
                assert result['match'] is False, f"S^{n} should NOT match"

    def test_s5_mismatch_is_2(self):
        """S^5 Green (1/d^3) mismatches Coulomb (1/d^1) by d^2."""
        result = check_singularity_match(5)
        assert result['mismatch'] == 2

    def test_coulomb_order_is_1(self):
        """Coulomb 1/r always has singularity order 1."""
        assert coulomb_singularity_order() == 1


# ============================================================================
# Coalescence manifold analysis
# ============================================================================

class TestCoalescence:
    """Test coalescence manifold structure."""

    def test_codimension_2(self):
        """Coalescence manifold has codimension 2 in S^5."""
        result = coalescence_codimension(angular_dim=5)
        assert result['n_coalescence_conditions'] == 2

    def test_normal_fiber_dim_2(self):
        """Normal bundle has 2D fibers."""
        result = coalescence_codimension(angular_dim=5)
        assert result['normal_fiber_dim'] == 2

    def test_coalescence_manifold_dim_3(self):
        """Coalescence manifold is 3-dimensional."""
        result = coalescence_codimension(angular_dim=5)
        assert result['coalescence_manifold_dim'] == 3

    def test_normal_fiber_does_not_match_s3(self):
        """Normal fiber dim=2 cannot contain S^3 (dim=3)."""
        result = coalescence_codimension(angular_dim=5)
        assert result['normal_fiber_matches_S3'] is False

    def test_two_conditions(self):
        """Coalescence requires exactly two conditions."""
        result = coalescence_codimension(angular_dim=5)
        assert len(result['coalescence_conditions']) == 2
        assert 'alpha = pi/4' in result['coalescence_conditions']
        assert 'theta_12 = 0' in result['coalescence_conditions']


# ============================================================================
# Fiber Green's function analysis
# ============================================================================

class TestFiber:
    """Test fiber Green's function analysis."""

    def test_no_viable_candidate(self):
        """No sphere in the 2D normal bundle can produce 1/d^1 Green function."""
        result = fiber_green_analysis()
        assert result['any_viable'] is False

    def test_need_s3_have_2d(self):
        """Need S^3 for 1/d^1 but normal bundle only has dim=2."""
        result = fiber_green_analysis()
        assert result['need_sphere_dim'] == 3
        assert result['have_fiber_dim'] == 2

    def test_s1_does_not_match(self):
        """S^1 fits in 2D fiber but has logarithmic Green function."""
        result = fiber_green_analysis()
        s1 = result['candidates']['S^1']
        assert s1['fits_in_normal_bundle'] is True
        assert s1['matches_coulomb_singularity'] is False

    def test_s2_does_not_match(self):
        """S^2 fits in 2D fiber but has logarithmic Green function."""
        result = fiber_green_analysis()
        s2 = result['candidates']['S^2']
        assert s2['fits_in_normal_bundle'] is True
        assert s2['matches_coulomb_singularity'] is False

    def test_s3_matches_but_too_big(self):
        """S^3 has the right Green function but doesn't fit in 2D fiber."""
        result = fiber_green_analysis()
        s3 = result['candidates']['S^3']
        assert s3['matches_coulomb_singularity'] is True
        assert s3['fits_in_normal_bundle'] is False
        assert s3['viable'] is False


# ============================================================================
# Conformal invariance
# ============================================================================

class TestConformal:
    """Test conformal invariance of singularity order."""

    def test_singularity_preserved(self):
        """Conformal rescaling preserves Green function singularity order."""
        result = conformal_singularity_preservation()
        assert 'preserved' in result['statement'].lower() or \
               'preserve' in result['consequence'].lower()

    def test_s5_cannot_be_fixed(self):
        """No conformal rescaling of S^5 can convert 1/d^3 to 1/d^1."""
        result = conformal_singularity_preservation()
        assert '1/d^3' in result['consequence']
        assert '1/d^1' in result['consequence']


# ============================================================================
# Numerical verification
# ============================================================================

class TestNumerical:
    """Numerical verification of singularity orders."""

    def test_r12_singularity_order_is_1(self):
        """1/r_12 has singularity order 1 near coalescence (numerically)."""
        result = verify_r12_singularity_at_coalescence(n_points=200)
        assert result['agreement'], \
            f"Measured order {result['measured_coulomb_singularity_order']:.4f}, expected 1.0"

    def test_r12_vs_s5_mismatch(self):
        """Mismatch between r_12 singularity and S^5 Green is 2."""
        result = verify_r12_singularity_at_coalescence()
        assert result['mismatch_factor'] == 2

    def test_s3_green_proportional_to_inv_d(self):
        """S^3 Green function is proportional to 1/d (numerically)."""
        result = verify_s3_green_matches_coulomb()
        assert result['numerical_green_s3_proportional_to_1_over_d']

    def test_s3_green_singularity_match(self):
        """S^3 analysis confirms match with Coulomb."""
        result = verify_s3_green_matches_coulomb()
        assert result['singularity_analysis']['match'] is True


# ============================================================================
# Full analysis integration
# ============================================================================

class TestFullAnalysis:
    """Test the complete Track W analysis."""

    def test_verdict_is_obstruction(self):
        """Full analysis yields structural obstruction."""
        result = full_cusp_graph_analysis()
        assert 'OBSTRUCTION' in result['verdict']

    def test_classification_is_dimensionality(self):
        """Obstruction is classified as dimensionality mismatch."""
        result = full_cusp_graph_analysis()
        assert 'imensionality' in result['classification']

    def test_level1_works(self):
        """Analysis confirms Level 1 (S^3) DOES work."""
        result = full_cusp_graph_analysis()
        assert result['sphere_analysis']['S^3']['match'] is True

    def test_level4_fails(self):
        """Analysis confirms Level 4 (S^5) does NOT work."""
        result = full_cusp_graph_analysis()
        assert result['sphere_analysis']['S^5']['match'] is False

    def test_no_viable_fiber(self):
        """No viable fiber decomposition exists."""
        result = full_cusp_graph_analysis()
        assert result['fiber']['any_viable'] is False

    def test_tracks_uvx_supported(self):
        """All three companion tracks are supported by this result."""
        result = full_cusp_graph_analysis()
        for track in ['Track_U', 'Track_V', 'Track_X']:
            assert 'SUPPORTED' in result['implications_tracks_UVX'][track]

    def test_paper18_classification(self):
        """Cusp is classified as EMBEDDING exchange constant."""
        result = full_cusp_graph_analysis()
        assert 'EMBEDDING' in result['paper18_classification']


# ============================================================================
# Mathematical consistency checks
# ============================================================================

class TestMathConsistency:
    """Cross-checks on the mathematical argument."""

    def test_s3_is_unique_match(self):
        """
        S^3 is the unique sphere where the Green's function matches 1/r.
        This is the foundational fact underlying the entire GeoVac framework.
        Verify it holds for a wide range of dimensions.
        """
        matches = []
        for n in range(1, 100):
            if check_singularity_match(n)['match']:
                matches.append(n)
        assert matches == [3], f"Expected only S^3 to match, got S^{matches}"

    def test_hyperspherical_r12_formula(self):
        """
        Verify r_12 = R_e * sqrt(1 - sin(2*alpha)*cos(theta_12))
        at known points.
        """
        # At alpha = pi/4, theta_12 = 0: r_12 = 0 (coalescence)
        R_e = 1.0
        r12 = R_e * np.sqrt(1 - np.sin(2 * np.pi / 4) * np.cos(0.0))
        assert abs(r12) < 1e-15

        # At alpha = pi/4, theta_12 = pi: r_12 = R_e * sqrt(2) (maximum)
        r12 = R_e * np.sqrt(1 - np.sin(np.pi / 2) * np.cos(np.pi))
        assert abs(r12 - np.sqrt(2)) < 1e-14

        # At alpha = 0: r_12 = R_e (one electron at origin)
        # sin(0) = 0, so r_12 = R_e * sqrt(1) = R_e
        r12 = R_e * np.sqrt(1 - np.sin(0) * np.cos(0.5))
        assert abs(r12 - 1.0) < 1e-14

    def test_coalescence_is_manifold_not_point(self):
        """
        The coalescence locus is a 3D manifold (codim 2 in S^5),
        NOT an isolated point. This is the structural advantage of
        hyperspherical coordinates (Paper 13, Paper 15).
        """
        result = coalescence_codimension(5)
        assert result['coalescence_manifold_dim'] == 3
        assert result['coalescence_manifold_dim'] > 0  # not a point

    def test_green_singularity_formula(self):
        """
        Verify Green's function singularity exponent formula:
        On S^n, G ~ 1/d^{n-2} for n >= 3.

        This is a standard result from Riemannian geometry.
        The fundamental solution of Delta u = delta on R^n is:
          n=1: |x|/2
          n=2: -(1/2pi) log|x|
          n>=3: c_n / |x|^{n-2}
        On S^n, the Green's function inherits the same local singularity.
        """
        # n=3 -> 1/d^1 (Coulomb!)
        assert greens_function_singularity_order(3) == 1
        # n=5 -> 1/d^3 (Level 4 angular space)
        assert greens_function_singularity_order(5) == 3
        # n=6 -> 1/d^4 (configuration space of 2 electrons in R^3)
        assert greens_function_singularity_order(6) == 4

    def test_dimension_counting(self):
        """
        Verify dimension counting for Level 4:
        - Full config space: R^6 (two electrons in R^3)
        - Hyperradius: 1 dimension
        - Angular space: 5 dimensions (= S^5)
        - Sigma reduction: Phi integrated out, but dim of S^5 remains 5
        """
        full_config_dim = 6
        hyperradial_dim = 1
        angular_dim = full_config_dim - hyperradial_dim
        assert angular_dim == 5


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
