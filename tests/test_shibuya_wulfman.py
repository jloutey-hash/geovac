"""
Tests for Shibuya-Wulfman two-center nuclear attraction integrals.

Track CD Phase 1A: validates the multipole expansion implementation
for cross-center V_ne matrix elements.
"""

import numpy as np
import pytest

from geovac.shibuya_wulfman import (
    _angular_coefficient,
    _radial_split_integral,
    compute_cross_center_vne,
    compute_cross_center_vne_element,
    convergence_study,
)
from geovac.composed_qubit import _enumerate_states


# ---------------------------------------------------------------------------
# Angular coefficient tests
# ---------------------------------------------------------------------------

class TestAngularCoefficient:
    """Verify Gaunt selection rules and angular coefficient values."""

    def test_m_diagonal(self) -> None:
        """m1 != m2 must give zero (z-axis potential, M=0)."""
        assert _angular_coefficient(1, 1, 1, 0, 1) == 0.0
        assert _angular_coefficient(1, -1, 1, 0, 1) == 0.0
        assert _angular_coefficient(0, 0, 1, 1, 1) == 0.0

    def test_parity_rule(self) -> None:
        """l1 + L + l2 must be even."""
        assert _angular_coefficient(0, 0, 0, 0, 1) == 0.0
        assert _angular_coefficient(1, 0, 0, 0, 0) == 0.0

    def test_triangle_inequality(self) -> None:
        """L must satisfy |l1-l2| <= L <= l1+l2."""
        assert _angular_coefficient(0, 0, 0, 0, 2) == 0.0
        assert _angular_coefficient(1, 0, 1, 0, 3) == 0.0

    def test_ss_L0(self) -> None:
        """For s-s (l=l'=0), L=0 coefficient should be exactly 1.0."""
        assert _angular_coefficient(0, 0, 0, 0, 0) == pytest.approx(1.0)

    def test_sp_L1(self) -> None:
        """For s-p (l=0, l'=1), only L=1 is allowed and nonzero."""
        assert _angular_coefficient(0, 0, 1, 0, 0) == 0.0
        val = _angular_coefficient(0, 0, 1, 0, 1)
        assert abs(val) > 0.1

    def test_pp_L0_and_L2(self) -> None:
        """For p-p (l=l'=1), L=0 and L=2 are allowed."""
        assert abs(_angular_coefficient(1, 0, 1, 0, 0)) > 0.1
        assert abs(_angular_coefficient(1, 0, 1, 0, 2)) > 0.1
        assert _angular_coefficient(1, 0, 1, 0, 1) == 0.0


# ---------------------------------------------------------------------------
# Analytical validation tests
# ---------------------------------------------------------------------------

class TestAnalyticalValidation:
    """Validate against known analytical results."""

    def test_1s_core_z3(self) -> None:
        """Core 1s-1s: <1s_Z=3|-1/|r-R||1s_Z=3> vs analytical formula."""
        Z_orb, Z_nuc, R = 3.0, 1.0, 3.015
        val = compute_cross_center_vne_element(
            Z_orb, 1, 0, 0, 1, 0, 0, Z_nuc, R, L_max=2,
        )
        exact = -Z_nuc * (1.0 / R) * (1.0 - (1.0 + Z_orb * R) * np.exp(-2 * Z_orb * R))
        assert val == pytest.approx(exact, rel=1e-12)

    def test_1s_valence_z1(self) -> None:
        """Valence 1s-1s: Z=1 orbitals feeling Z_nuc=1."""
        Z_orb, Z_nuc, R = 1.0, 1.0, 3.015
        val = compute_cross_center_vne_element(
            Z_orb, 1, 0, 0, 1, 0, 0, Z_nuc, R, L_max=2,
        )
        exact = -Z_nuc * (1.0 / R) * (1.0 - (1.0 + Z_orb * R) * np.exp(-2 * Z_orb * R))
        assert val == pytest.approx(exact, rel=1e-12)

    def test_h_side_feeling_li(self) -> None:
        """H-side 1s-1s: Z_orb=1 feeling Z_nuc=3."""
        Z_orb, Z_nuc, R = 1.0, 3.0, 3.015
        val = compute_cross_center_vne_element(
            Z_orb, 1, 0, 0, 1, 0, 0, Z_nuc, R, L_max=2,
        )
        exact = -Z_nuc * (1.0 / R) * (1.0 - (1.0 + Z_orb * R) * np.exp(-2 * Z_orb * R))
        assert val == pytest.approx(exact, rel=1e-12)

    def test_large_r_point_charge(self) -> None:
        """At large R, V_ne -> -Z_B/R."""
        val = compute_cross_center_vne_element(
            1.0, 1, 0, 0, 1, 0, 0, 1.0, 30.0, L_max=2, n_grid=4000,
        )
        assert val == pytest.approx(-1.0 / 30.0, rel=0.01)


# ---------------------------------------------------------------------------
# Matrix property tests
# ---------------------------------------------------------------------------

class TestMatrixProperties:
    """Verify structural properties of the V_ne matrix."""

    def test_hermiticity(self) -> None:
        """V_ne matrix must be symmetric."""
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(3.0, states, 1.0, 3.015, L_max=2, n_grid=4000)
        np.testing.assert_allclose(vne, vne.T, atol=1e-14)

    def test_m_diagonal_matrix(self) -> None:
        """All off-m-diagonal elements must be zero."""
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(1.0, states, 1.0, 3.015, L_max=2, n_grid=4000)
        for i, (_, _, m1) in enumerate(states):
            for j, (_, _, m2) in enumerate(states):
                if m1 != m2:
                    assert abs(vne[i, j]) < 1e-14

    def test_all_diagonal_negative(self) -> None:
        """Diagonal V_ne elements must be negative (attractive)."""
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(1.0, states, 1.0, 3.015, L_max=2, n_grid=4000)
        for i in range(len(states)):
            assert vne[i, i] < 0

    def test_nonzero_count_nmax2(self) -> None:
        """For n_max=2 (5 orbitals), check nonzero element count."""
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(3.0, states, 1.0, 3.015, L_max=2, n_grid=4000)
        nz = np.count_nonzero(np.abs(vne) > 1e-10)
        assert 8 <= nz <= 15


# ---------------------------------------------------------------------------
# Convergence tests
# ---------------------------------------------------------------------------

class TestConvergence:
    """Verify L_max convergence properties."""

    def test_ss_converges_at_L0(self) -> None:
        """s-s integrals need only L=0."""
        conv = convergence_study(3.0, 1, 0, 0, 1, 0, 0, 1.0, 3.015, range(4), n_grid=4000)
        vals = list(conv.values())
        for v in vals[1:]:
            assert v == pytest.approx(vals[0], abs=1e-14)

    def test_sp_converges_at_L1(self) -> None:
        """s-p integrals need only L=1."""
        conv = convergence_study(3.0, 1, 0, 0, 2, 1, 0, 1.0, 3.015, range(4), n_grid=4000)
        assert conv[0] == pytest.approx(0.0, abs=1e-14)
        for L in range(2, 4):
            assert conv[L] == pytest.approx(conv[1], abs=1e-14)

    def test_lmax2_sufficient_for_nmax2(self) -> None:
        """For n_max=2 (l_max=1), L_max=2 equals L_max=10."""
        states = _enumerate_states(2)
        vne_L2 = compute_cross_center_vne(1.0, states, 1.0, 3.015, L_max=2, n_grid=8000)
        vne_L10 = compute_cross_center_vne(1.0, states, 1.0, 3.015, L_max=10, n_grid=8000)
        np.testing.assert_allclose(vne_L2, vne_L10, atol=1e-12)


# ---------------------------------------------------------------------------
# Analytical vs grid cross-validation
# ---------------------------------------------------------------------------

class TestAnalyticalVsGrid:
    """Verify analytical integrals match grid at high resolution."""

    def test_analytical_vs_grid_nmax2(self) -> None:
        """All n_max=2 V_ne elements: analytical vs grid at n_grid=32000."""
        from geovac.shibuya_wulfman import _radial_split_integral_grid

        states = _enumerate_states(2)
        # Analytical
        vne_anal = compute_cross_center_vne(1.0, states, 1.0, 3.015, L_max=2)
        # Grid-based at high resolution
        N = len(states)
        vne_grid = np.zeros((N, N))
        unique_nl = sorted(set((n, l) for n, l, m in states))
        from geovac.shibuya_wulfman import _angular_coefficient
        rad_cache: dict = {}
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                for L in range(3):
                    if (l1 + L + l2) % 2 != 0 or abs(l1 - l2) > L or L > l1 + l2:
                        continue
                    rad_cache[(n1, l1, n2, l2, L)] = _radial_split_integral_grid(
                        1.0, n1, l1, n2, l2, L, 3.015, n_grid=32000,
                    )
        for i, (n1, l1, m1) in enumerate(states):
            for j, (n2, l2, m2) in enumerate(states):
                if m1 != m2:
                    continue
                val = 0.0
                for L in range(3):
                    if (l1 + L + l2) % 2 != 0 or abs(l1 - l2) > L or L > l1 + l2:
                        continue
                    ang = _angular_coefficient(l1, m1, l2, m2, L)
                    if abs(ang) < 1e-15:
                        continue
                    val += ang * rad_cache.get((n1, l1, n2, l2, L), 0.0)
                vne_grid[i, j] = -1.0 * val

        # Grid at 32K should agree with analytical to ~0.5%
        # (2s orbitals with nodes converge slower on uniform grids)
        mask = np.abs(vne_anal) > 1e-10
        rel_err = np.abs(vne_anal[mask] - vne_grid[mask]) / np.abs(vne_anal[mask])
        assert np.max(rel_err) < 5e-3, f"Max rel error {np.max(rel_err):.2e} > 5e-3"

    def test_analytical_machine_precision_1s1s(self) -> None:
        """1s-1s integrals at machine precision vs closed-form formula."""
        R = 3.015
        for Z_orb, Z_nuc in [(3.0, 1.0), (1.0, 1.0), (1.0, 3.0)]:
            val = compute_cross_center_vne_element(
                Z_orb, 1, 0, 0, 1, 0, 0, Z_nuc, R, L_max=2,
            )
            exact = -Z_nuc / R * (1.0 - (1.0 + Z_orb * R) * np.exp(-2 * Z_orb * R))
            assert val == pytest.approx(exact, rel=1e-12)


# ---------------------------------------------------------------------------
# Non-collinear rotation tests
# ---------------------------------------------------------------------------


class TestNonCollinearRotation:
    """Verify direction-based rotation for off-axis nuclei."""

    def test_z_axis_reproduces_existing(self) -> None:
        """direction=(0,0,1) must match nuc_parity=+1 exactly."""
        states = _enumerate_states(2)
        vne_old = compute_cross_center_vne(
            3.0, states, 1.0, 3.015, L_max=2, nuc_parity=1,
        )
        vne_new = compute_cross_center_vne(
            3.0, states, 1.0, 3.015, L_max=2, direction=(0, 0, 1),
        )
        np.testing.assert_allclose(vne_new, vne_old, atol=1e-14)

    def test_neg_z_matches_nuc_parity(self) -> None:
        """direction=(0,0,-1) must match nuc_parity=-1 exactly."""
        states = _enumerate_states(2)
        vne_parity = compute_cross_center_vne(
            3.0, states, 1.0, 3.015, L_max=2, nuc_parity=-1,
        )
        vne_dir = compute_cross_center_vne(
            3.0, states, 1.0, 3.015, L_max=2, direction=(0, 0, -1),
        )
        np.testing.assert_allclose(vne_dir, vne_parity, atol=1e-14)

    def test_ss_rotation_invariant(self) -> None:
        """s-s matrix elements must be the same for any direction."""
        states = _enumerate_states(2)
        dirs = [
            (0, 0, 1),
            (0, 0, -1),
            (1, 0, 0),
            (0, 1, 0),
            (1 / np.sqrt(2), 0, 1 / np.sqrt(2)),
            (np.sin(52.25 * np.pi / 180), 0, np.cos(52.25 * np.pi / 180)),
        ]
        vne_ref = compute_cross_center_vne(
            1.0, states, 1.0, 3.015, L_max=2, direction=dirs[0],
        )
        # s-orbital indices: 0 (1s) and 1 (2s)
        s_indices = [i for i, (n, l, m) in enumerate(states) if l == 0]
        for d in dirs[1:]:
            vne = compute_cross_center_vne(
                1.0, states, 1.0, 3.015, L_max=2, direction=d,
            )
            for i in s_indices:
                for j in s_indices:
                    np.testing.assert_allclose(
                        vne[i, j], vne_ref[i, j], atol=1e-14,
                    )

    def test_hermiticity_arbitrary_direction(self) -> None:
        """V_ne must be symmetric for any direction."""
        states = _enumerate_states(2)
        direction = (np.sin(52.25 * np.pi / 180), 0, np.cos(52.25 * np.pi / 180))
        vne = compute_cross_center_vne(
            1.0, states, 1.0, 3.015, L_max=2, direction=direction,
        )
        np.testing.assert_allclose(vne, vne.T, atol=1e-14)

    def test_90deg_px_stronger_than_pz(self) -> None:
        """For nucleus along x, p_x coupling > p_z coupling to s orbital."""
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(
            1.0, states, 1.0, 3.015, L_max=2, direction=(1, 0, 0),
        )
        # p orbital indices for n=2: m=-1(y) at 2, m=0(z) at 3, m=+1(x) at 4
        # s orbital at 0 (1s)
        assert abs(vne[0, 4]) > abs(vne[0, 3]) + 1e-10  # s-px > s-pz
        assert abs(vne[0, 3]) < 1e-14  # s-pz should be zero for x-axis

    def test_x_axis_sp_equals_z_axis_sp(self) -> None:
        """s-p coupling along nucleus direction should be same magnitude."""
        states = _enumerate_states(2)
        vne_z = compute_cross_center_vne(
            1.0, states, 1.0, 3.015, L_max=2, direction=(0, 0, 1),
        )
        vne_x = compute_cross_center_vne(
            1.0, states, 1.0, 3.015, L_max=2, direction=(1, 0, 0),
        )
        # s-pz for z-axis should equal s-px for x-axis
        np.testing.assert_allclose(
            abs(vne_z[0, 3]), abs(vne_x[0, 4]), atol=1e-14,
        )

    def test_h2o_angle_direction(self) -> None:
        """Test at the H2O half-bond-angle (52.25 deg)."""
        half_angle = 52.25 * np.pi / 180
        direction = (np.sin(half_angle), 0, np.cos(half_angle))
        states = _enumerate_states(2)
        vne = compute_cross_center_vne(
            6.0, states, 1.0, 1.809, L_max=2, direction=direction,
        )
        # Must be symmetric
        np.testing.assert_allclose(vne, vne.T, atol=1e-14)
        # All diagonal must be negative (attractive)
        for i in range(len(states)):
            assert vne[i, i] < 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
