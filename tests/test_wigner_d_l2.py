"""
Tests for General-l Real Spherical Harmonic Rotation Matrices (Track CM)
========================================================================

Validates the Wigner D-matrix rotation for l=0, 1, 2, and 3, using:
  - Orthonormality: R^T R = I for random Euler angles
  - Identity: rotation by (0,0) gives identity matrix
  - Known parity: rotation by (pi, 0) gives known diagonal
  - Consistency: l=1 general matches legacy _build_rotation_matrix_l1
  - Block rotation: mixed l=0,1,2 states produce orthogonal block matrix
  - Determinant: det(R) = +1 (proper rotation, no reflections)
  - Composition: R(theta1) @ R(theta2) = R(theta1+theta2) for z-axis rotations

Author: GeoVac Development Team (Track CM, v2.1.1)
Date: April 2026
"""

import numpy as np
import pytest

from geovac.shibuya_wulfman import (
    _build_block_rotation_matrix,
    _build_rotation_matrix_l1,
    _build_rotation_matrix_real_sh,
    _small_wigner_d_element,
)


# ---------------------------------------------------------------------------
# Small Wigner d-matrix element tests
# ---------------------------------------------------------------------------

class TestSmallWignerD:
    """Tests for the complex small Wigner d-matrix element."""

    def test_d0_identity(self) -> None:
        """d^0_{0,0}(beta) = 1 for all beta."""
        for beta in [0, 0.5, 1.0, np.pi]:
            assert abs(_small_wigner_d_element(0, 0, 0, beta) - 1.0) < 1e-14

    def test_d1_diagonal_zero(self) -> None:
        """d^1_{0,0}(beta) = cos(beta)."""
        for beta in [0, 0.3, np.pi / 2, np.pi]:
            expected = np.cos(beta)
            got = _small_wigner_d_element(1, 0, 0, beta)
            assert abs(got - expected) < 1e-14, f"beta={beta}: {got} vs {expected}"

    def test_d1_off_diagonal(self) -> None:
        """d^1_{1,0}(beta) = -sin(beta)/sqrt(2)."""
        for beta in [0.3, np.pi / 4, np.pi / 2]:
            expected = -np.sin(beta) / np.sqrt(2)
            got = _small_wigner_d_element(1, 1, 0, beta)
            assert abs(got - expected) < 1e-13, f"beta={beta}: {got} vs {expected}"

    def test_d2_diagonal_zero(self) -> None:
        """d^2_{0,0}(beta) = (3cos^2(beta) - 1)/2."""
        for beta in [0, 0.5, np.pi / 3, np.pi]:
            expected = (3.0 * np.cos(beta) ** 2 - 1.0) / 2.0
            got = _small_wigner_d_element(2, 0, 0, beta)
            assert abs(got - expected) < 1e-13

    def test_unitarity(self) -> None:
        """Sum_m |d^l_{m',m}|^2 = 1 for any l, m', beta."""
        rng = np.random.default_rng(42)
        for l in [1, 2, 3]:
            for m1 in range(-l, l + 1):
                beta = rng.uniform(0, np.pi)
                total = sum(
                    _small_wigner_d_element(l, m1, m2, beta) ** 2
                    for m2 in range(-l, l + 1)
                )
                assert abs(total - 1.0) < 1e-12, (
                    f"l={l}, m'={m1}, beta={beta:.3f}: sum={total}"
                )


# ---------------------------------------------------------------------------
# Real SH rotation matrix tests
# ---------------------------------------------------------------------------

class TestRealSHRotation:
    """Tests for the real spherical harmonic rotation matrix."""

    @pytest.mark.parametrize("l", [0, 1, 2, 3])
    def test_identity(self, l: int) -> None:
        """Rotation by (0, 0) is identity."""
        R = _build_rotation_matrix_real_sh(l, 0.0, 0.0)
        dim = 2 * l + 1
        assert R.shape == (dim, dim)
        np.testing.assert_allclose(R, np.eye(dim), atol=1e-14)

    @pytest.mark.parametrize("l", [0, 1, 2, 3])
    def test_orthogonality_random(self, l: int) -> None:
        """R^T R = I for 100 random (theta, phi) pairs."""
        rng = np.random.default_rng(123 + l)
        dim = 2 * l + 1
        for _ in range(100):
            theta = rng.uniform(0, np.pi)
            phi = rng.uniform(0, 2 * np.pi)
            R = _build_rotation_matrix_real_sh(l, theta, phi)
            np.testing.assert_allclose(
                R.T @ R, np.eye(dim), atol=1e-12,
                err_msg=f"l={l}, theta={theta:.4f}, phi={phi:.4f}",
            )

    @pytest.mark.parametrize("l", [0, 1, 2, 3])
    def test_determinant_positive(self, l: int) -> None:
        """det(R) = +1 (proper rotation)."""
        rng = np.random.default_rng(456 + l)
        for _ in range(20):
            theta = rng.uniform(0, np.pi)
            phi = rng.uniform(0, 2 * np.pi)
            R = _build_rotation_matrix_real_sh(l, theta, phi)
            assert abs(np.linalg.det(R) - 1.0) < 1e-12

    def test_l1_matches_legacy(self) -> None:
        """General l=1 matches legacy _build_rotation_matrix_l1."""
        rng = np.random.default_rng(789)
        for _ in range(50):
            theta = rng.uniform(0, np.pi)
            phi = rng.uniform(0, 2 * np.pi)
            direction = (
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta),
            )
            R_legacy = _build_rotation_matrix_l1(direction)
            R_general = _build_rotation_matrix_real_sh(1, theta, phi)
            np.testing.assert_allclose(
                R_general, R_legacy, atol=1e-13,
                err_msg=f"theta={theta:.4f}, phi={phi:.4f}",
            )

    def test_l2_theta_pi(self) -> None:
        """Rotation by pi around y-axis has known structure for l=2."""
        R = _build_rotation_matrix_real_sh(2, np.pi, 0.0)
        # Should be diagonal (z-axis rotation by pi in real SH)
        # Off-diagonals should be zero to machine precision
        assert np.max(np.abs(R - np.diag(np.diag(R)))) < 1e-12
        # Diagonal entries are ±1
        for d in np.diag(R):
            assert abs(abs(d) - 1.0) < 1e-12

    def test_l2_composition(self) -> None:
        """R(theta1, phi) @ R(theta2, phi) for sequential rotations."""
        # For phi=0 rotations (around y-axis), composition should be exact
        R1 = _build_rotation_matrix_real_sh(2, 0.3, 0.0)
        R2 = _build_rotation_matrix_real_sh(2, 0.5, 0.0)
        R12 = _build_rotation_matrix_real_sh(2, 0.8, 0.0)
        # R(0.8, 0) should equal R(0.5, 0) @ R(0.3, 0)
        np.testing.assert_allclose(R2 @ R1, R12, atol=1e-12)


# ---------------------------------------------------------------------------
# Block rotation matrix tests
# ---------------------------------------------------------------------------

class TestBlockRotation:
    """Tests for the block-diagonal rotation matrix builder."""

    def test_l0_only(self) -> None:
        """Pure s-orbitals: block rotation is identity."""
        states = [(1, 0, 0), (2, 0, 0), (3, 0, 0)]
        D = _build_block_rotation_matrix(states, (0.5, 0.5, 1.0 / np.sqrt(2)))
        np.testing.assert_allclose(D, np.eye(3), atol=1e-14)

    def test_mixed_l01(self) -> None:
        """Mixed l=0 and l=1: block-diagonal structure."""
        states = [
            (1, 0, 0),
            (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1),
        ]
        D = _build_block_rotation_matrix(states, (0, 0, -1))
        assert D.shape == (5, 5)
        np.testing.assert_allclose(D.T @ D, np.eye(5), atol=1e-13)

    def test_mixed_l012(self) -> None:
        """Mixed l=0, 1, 2: n_max=3 orbital set."""
        states = [
            (1, 0, 0),
            (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1),
            (3, 0, 0), (3, 1, -1), (3, 1, 0), (3, 1, 1),
            (3, 2, -2), (3, 2, -1), (3, 2, 0), (3, 2, 1), (3, 2, 2),
        ]
        rng = np.random.default_rng(111)
        for _ in range(20):
            d = rng.normal(size=3)
            d = d / np.linalg.norm(d)
            direction = tuple(d)
            D = _build_block_rotation_matrix(states, direction)
            assert D.shape == (14, 14)
            np.testing.assert_allclose(
                D.T @ D, np.eye(14), atol=1e-12,
            )
            assert abs(np.linalg.det(D) - 1.0) < 1e-11

    def test_block_diagonal_structure(self) -> None:
        """Off-block elements are zero."""
        states = [
            (1, 0, 0),                                       # block 0: [0]
            (2, 0, 0),                                       # block 1: [1]
            (2, 1, -1), (2, 1, 0), (2, 1, 1),               # block 2: [2,3,4]
            (3, 2, -2), (3, 2, -1), (3, 2, 0), (3, 2, 1), (3, 2, 2),  # block 3: [5..9]
        ]
        D = _build_block_rotation_matrix(states, (1, 0, 0))
        # s-orbital blocks should be identity (1x1)
        assert abs(D[0, 0] - 1.0) < 1e-14
        assert abs(D[1, 1] - 1.0) < 1e-14
        # Cross-block entries should be zero
        # Block 0 vs block 2
        assert np.max(np.abs(D[0, 2:5])) < 1e-14
        assert np.max(np.abs(D[2:5, 0])) < 1e-14
        # Block 2 vs block 3
        assert np.max(np.abs(D[2:5, 5:10])) < 1e-14
        assert np.max(np.abs(D[5:10, 2:5])) < 1e-14

    def test_empty_states(self) -> None:
        """Empty state list returns 0x0 identity."""
        D = _build_block_rotation_matrix([], (0, 0, 1))
        assert D.shape == (0, 0)

    def test_z_direction_is_identity(self) -> None:
        """Direction (0,0,1) gives identity rotation."""
        states = [
            (1, 0, 0),
            (2, 1, -1), (2, 1, 0), (2, 1, 1),
            (3, 2, -2), (3, 2, -1), (3, 2, 0), (3, 2, 1), (3, 2, 2),
        ]
        D = _build_block_rotation_matrix(states, (0, 0, 1))
        np.testing.assert_allclose(D, np.eye(9), atol=1e-13)
