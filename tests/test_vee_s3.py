"""
Tests for V_ee as S3 density-overlap integral (Paper 7, Section VI).

Verifies the master formula F0(a,b) = (4Z/pi) * integral_0^inf Phi_a(t)*Phi_b(t) dt
against the three known exact Slater integrals.

These tests are TOPOLOGICAL INTEGRITY CHECKS analogous to test_fock_projection.py.
They must never be bypassed or loosened.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from geovac.lattice_index import compute_vee_s3_overlap, _phi_s_orbital


# ---------------------------------------------------------------------------
# Exact reference values (Slater integrals, Z=1)
# ---------------------------------------------------------------------------
F0_1s1s_EXACT = 5.0 / 8.0          # = 0.625
F0_1s2s_EXACT = 17.0 / 81.0        # ~ 0.20988
F0_2s2s_EXACT = 77.0 / 512.0       # ~ 0.15039

# Analytical integral values for density overlaps
INT_PHI1s_SQ = 5.0 * np.pi / 32.0        # integral_0^inf 1/(1+t^2)^4 dt
INT_PHI2s_SQ = 77.0 * np.pi / 2048.0     # integral_0^inf Phi_2s(t)^2 dt

# Tolerance for topological integrity
TOL = 1e-4   # < 0.01%


class TestVeeS3Integrals:
    """Verify F0 integrals via S3 density-overlap formula."""

    def test_vee_s3_1s1s(self) -> None:
        """F0(1s,1s) = 5/8 (Z=1)."""
        result = compute_vee_s3_overlap(1, 0, 1, 0, Z=1.0)
        assert abs(result - F0_1s1s_EXACT) / F0_1s1s_EXACT < TOL, (
            f"F0(1s,1s) = {result:.8f}, expected {F0_1s1s_EXACT:.8f}"
        )

    def test_vee_s3_1s2s(self) -> None:
        """F0(1s,2s) = 17/81 (Z=1)."""
        result = compute_vee_s3_overlap(1, 0, 2, 0, Z=1.0)
        assert abs(result - F0_1s2s_EXACT) / F0_1s2s_EXACT < TOL, (
            f"F0(1s,2s) = {result:.8f}, expected {F0_1s2s_EXACT:.8f}"
        )

    def test_vee_s3_2s2s(self) -> None:
        """F0(2s,2s) = 77/512 (Z=1)."""
        result = compute_vee_s3_overlap(2, 0, 2, 0, Z=1.0)
        assert abs(result - F0_2s2s_EXACT) / F0_2s2s_EXACT < TOL, (
            f"F0(2s,2s) = {result:.8f}, expected {F0_2s2s_EXACT:.8f}"
        )

    def test_vee_s3_z_scaling(self) -> None:
        """F0(1s,1s, Z=2) = 5*Z/8 = 1.25."""
        expected = 5.0 * 2.0 / 8.0  # = 1.25
        result = compute_vee_s3_overlap(1, 0, 1, 0, Z=2.0)
        assert abs(result - expected) / expected < TOL, (
            f"F0(1s,1s, Z=2) = {result:.8f}, expected {expected:.8f}"
        )


class TestVeeS3DensityIntegrals:
    """Verify the raw density overlap integrals on S3."""

    def test_phi_1s_integral(self) -> None:
        """integral_0^inf 1/(1+t^2)^4 dt = 5*pi/32."""
        def integrand(t: float) -> float:
            return _phi_s_orbital(1, t) ** 2
        result, _ = quad(integrand, 0, np.inf, limit=200)
        assert abs(result - INT_PHI1s_SQ) < 1e-10, (
            f"integral Phi_1s^2 = {result:.12f}, expected {INT_PHI1s_SQ:.12f}"
        )

    def test_phi_2s_integral(self) -> None:
        """integral_0^inf Phi_2s(t)^2 dt = 77*pi/2048."""
        def integrand(t: float) -> float:
            return _phi_s_orbital(2, t) ** 2
        result, _ = quad(integrand, 0, np.inf, limit=200)
        assert abs(result - INT_PHI2s_SQ) < 1e-10, (
            f"integral Phi_2s^2 = {result:.12f}, expected {INT_PHI2s_SQ:.12f}"
        )


class TestVeeS3NodeVsEdge:
    """Permanently documents that edge-distance V_ee is architecturally wrong."""

    def test_vee_s3_node_not_edge(self) -> None:
        """Chordal ansatz overestimates F0(1s,2s) by more than 10x.

        The chordal distance d^2(1s,2s) = 0.4 with kappa_ee = 5Z/2 gives
        V_ee = kappa/d^2 = 6.25 (Z=1), while the exact Slater integral is
        F0(1s,2s) = 17/81 ~ 0.21. Ratio > 29x.

        This test permanently documents why edge-distance V_ee is wrong for
        cross-shell pairs. Do NOT weaken this assertion.
        """
        # S3 overlap (correct)
        f0_exact = compute_vee_s3_overlap(1, 0, 2, 0, Z=1.0)

        # Chordal ansatz (wrong for cross-shell)
        kappa_ee = 5.0 / 2.0  # Z=1
        d_sq_1s2s = 0.4        # from S3 coordinate geometry
        v_chordal = kappa_ee / d_sq_1s2s

        ratio = v_chordal / f0_exact
        assert ratio > 10.0, (
            f"Chordal/exact ratio = {ratio:.1f}x, expected > 10x. "
            f"Chordal = {v_chordal:.4f}, exact = {f0_exact:.6f}"
        )

    def test_l_nonzero_raises(self) -> None:
        """l>0 orbitals must raise NotImplementedError, not silently fallback."""
        with pytest.raises(NotImplementedError, match="l>0"):
            compute_vee_s3_overlap(1, 0, 2, 1, Z=1.0)
