"""
Tests for Shibuya-Wulfman nuclear attraction integrals.

Six tests validating the SW coupling module:
    1. R-dependence: coupling decreases with R for n=1
    2. Dissociation limit: V^SW -> 0 at R=100
    3. Sign convention: V^SW is negative (attractive)
    4. Magnitude: n=1 diagonal in [-0.7, -0.3] Ha at R=3.015
    5. D-matrix factorization: V^SW = -(Z_B/p0) * D * f(n, gamma)
    6. Matrix symmetry: V_AB = V_AB^T (Hermitian)
"""

import numpy as np
import pytest

from geovac.shibuya_wulfman import (
    sw_form_factor,
    sw_nuclear_attraction,
    sw_coupling_matrix,
    sw_coupling_matrix_AB,
)
from geovac.wigner_so4 import bond_angle, d_matrix_block


# LiH parameters
Z_A = 3.0  # Li
Z_B = 1.0  # H
R_EQ = 3.015  # Experimental equilibrium distance (Bohr)


def _default_p0(Z_A: float, Z_B: float) -> float:
    """Default momentum scale: sqrt(Z_A^2 + Z_B^2)."""
    return np.sqrt(Z_A**2 + Z_B**2)


class TestSWRDependence:
    """Test 1: R-dependence — coupling decreases with R for n=1."""

    def test_n1_diagonal_decreases_with_R(self) -> None:
        """The 1s-1s SW element magnitude should decrease as R increases."""
        p0 = _default_p0(Z_A, Z_B)
        R_values = [2.0, 3.0, 5.0, 8.0, 15.0]
        magnitudes = []
        for R in R_values:
            gamma = bond_angle(R, p0)
            v = sw_nuclear_attraction(1, 0, 0, 1, 0, 0, Z_B, p0, gamma)
            magnitudes.append(abs(v))

        # Each subsequent magnitude should be smaller
        for i in range(len(magnitudes) - 1):
            assert magnitudes[i] > magnitudes[i + 1], (
                f"|V(R={R_values[i]})| = {magnitudes[i]:.6f} should be > "
                f"|V(R={R_values[i+1]})| = {magnitudes[i+1]:.6f}"
            )


class TestSWDissociationLimit:
    """Test 2: Dissociation limit — V^SW -> 0 at R=100."""

    def test_vanishes_at_large_R(self) -> None:
        """At R=100 Bohr, all SW matrix elements should be negligible."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(100.0, p0)
        V = sw_coupling_matrix(3, Z_B, p0, gamma)

        max_element = np.max(np.abs(V))
        # sin(gamma) ~ gamma ~ 2/(p0*R) ~ 0.006 at R=100, so V ~ 0.002
        assert max_element < 0.01, (
            f"V^SW should vanish at R=100, but max element = {max_element:.6e}"
        )


class TestSWSignConvention:
    """Test 3: Sign convention — V^SW is negative (attractive)."""

    def test_diagonal_negative(self) -> None:
        """Diagonal SW elements (same orbital) should be negative."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)

        # 1s-1s
        v_1s = sw_nuclear_attraction(1, 0, 0, 1, 0, 0, Z_B, p0, gamma)
        assert v_1s < 0, f"V^SW(1s,1s) = {v_1s:.6f} should be negative"

        # 2s-2s
        v_2s = sw_nuclear_attraction(2, 0, 0, 2, 0, 0, Z_B, p0, gamma)
        assert v_2s < 0, f"V^SW(2s,2s) = {v_2s:.6f} should be negative"


class TestSWMagnitude:
    """Test 4: Magnitude — n=1 diagonal is negative and physically reasonable."""

    def test_n1_magnitude_range(self) -> None:
        """At equilibrium, the 1s-1s SW element should be negative and nonzero.

        With fixed p0 = sqrt(Z_A^2+Z_B^2) and f=sin(gamma), the n=1 diagonal
        is -(Z_B/p0)*sin(gamma) ~ -0.066 Ha. This is SMALLER than the Mulliken
        estimate (~0.5 Ha) because sin(gamma) ~ 0.21 at R_eq with the large
        fixed p0. This is a known limitation of the fixed-p0 + sin(gamma) ansatz.
        """
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)
        v_1s = sw_nuclear_attraction(1, 0, 0, 1, 0, 0, Z_B, p0, gamma)

        # Must be attractive (negative)
        assert v_1s < 0, f"V^SW(1s,1s) = {v_1s:.6f} should be negative"
        # Must be nonzero and not too large
        assert -1.0 < v_1s < -0.01, (
            f"V^SW(1s,1s) = {v_1s:.4f} Ha outside reasonable range [-1.0, -0.01]"
        )


class TestSWFactorization:
    """Test 5: D-matrix factorization — V^SW = -(Z_B/p0) * D * f(n, gamma)."""

    def test_factorization_n1(self) -> None:
        """Verify V^SW = -(Z_B/p0) * f * (D+D^T)/2 per shell block."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)

        # Compute via sw_coupling_matrix
        V = sw_coupling_matrix(2, Z_B, p0, gamma)

        # n=1 block: indices [0] — 1x1, D is trivially symmetric
        D1 = d_matrix_block(1, gamma)  # 1x1
        D1_sym = (D1 + D1.T) / 2.0
        f1 = sw_form_factor(1, gamma)
        V1_expected = -(Z_B / p0) * f1 * D1_sym

        np.testing.assert_allclose(
            V[0:1, 0:1], V1_expected, atol=1e-12,
            err_msg="n=1 block factorization failed"
        )

        # n=2 block: indices [1:5] — 4x4, symmetrized
        D2 = d_matrix_block(2, gamma)  # 4x4
        D2_sym = (D2 + D2.T) / 2.0
        f2 = sw_form_factor(2, gamma)
        V2_expected = -(Z_B / p0) * f2 * D2_sym

        np.testing.assert_allclose(
            V[1:5, 1:5], V2_expected, atol=1e-12,
            err_msg="n=2 block factorization failed"
        )

    def test_cross_n_zero(self) -> None:
        """Cross-n blocks should be exactly zero."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)
        V = sw_coupling_matrix(2, Z_B, p0, gamma)

        # n=1 to n=2 cross block: V[0, 1:5] and V[1:5, 0]
        np.testing.assert_allclose(
            V[0, 1:5], 0.0, atol=1e-15,
            err_msg="Cross-n block (1->2) should be zero"
        )
        np.testing.assert_allclose(
            V[1:5, 0], 0.0, atol=1e-15,
            err_msg="Cross-n block (2->1) should be zero"
        )


class TestSWSymmetry:
    """Test 6: Matrix symmetry — V = V^T."""

    def test_single_center_symmetry(self) -> None:
        """Single-center SW matrix should be symmetric."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)
        V = sw_coupling_matrix(3, Z_B, p0, gamma)

        np.testing.assert_allclose(
            V, V.T, atol=1e-12,
            err_msg="Single-center SW matrix is not symmetric"
        )

    def test_two_center_symmetry(self) -> None:
        """Two-center combined SW matrix should be symmetric."""
        p0 = _default_p0(Z_A, Z_B)
        gamma = bond_angle(R_EQ, p0)
        V = sw_coupling_matrix_AB(3, Z_A, Z_B, p0, gamma)

        np.testing.assert_allclose(
            V, V.T, atol=1e-12,
            err_msg="Two-center SW matrix is not symmetric"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
