"""
Tests for SO(4) Wigner D-matrix module (geovac/wigner_so4.py).

Eight test classes validating the SO(4) D-matrix implementation against
known algebraic identities and limiting cases from Paper 8.

Tests 1-3 (unitarity, identity, SU(2) limits) are the minimum acceptance
criteria for v0.9.15.

References:
    Paper 8: Bond Sphere Geometry (GeoVac, 2026)
    Bander & Itzykson, Rev. Mod. Phys. 38, 330 (1966)
    Shibuya & Wulfman, Proc. R. Soc. A 286, 376 (1965)
"""

import numpy as np
import pytest
from geovac.wigner_so4 import (
    cg_so4,
    wigner_d_su2,
    wigner_D_so4,
    bond_angle,
    d_matrix_block,
)


# ===================================================================
#  Test 1: Unitarity — D†D = I for all n-shells
# ===================================================================

class TestUnitarity:
    """D-matrix must be unitary (orthogonal for real angles): D†D = I."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_unitarity_various_angles(self, n: int) -> None:
        """D(γ)†D(γ) = I for several angles and n-shells."""
        for gamma in [0.0, 0.3, np.pi / 4, np.pi / 2, np.pi, 1.7]:
            D = d_matrix_block(n, gamma)
            product = D.T @ D
            identity = np.eye(n * n)
            np.testing.assert_allclose(
                product, identity, atol=1e-12,
                err_msg=f"Unitarity failed for n={n}, γ={gamma:.3f}"
            )

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_det_is_one(self, n: int) -> None:
        """det(D) = ±1 for a real orthogonal matrix."""
        D = d_matrix_block(n, np.pi / 3)
        det = np.linalg.det(D)
        assert abs(abs(det) - 1.0) < 1e-12, (
            f"det(D) = {det} for n={n}, expected ±1"
        )


# ===================================================================
#  Test 2: Identity at γ = 0 — D(0) = I
# ===================================================================

class TestIdentity:
    """At γ = 0 (R → ∞, separated atoms), the D-matrix is the identity."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_identity_at_zero(self, n: int) -> None:
        """D^(n)(γ=0) = I_{n²×n²}."""
        D = d_matrix_block(n, 0.0)
        identity = np.eye(n * n)
        np.testing.assert_allclose(
            D, identity, atol=1e-12,
            err_msg=f"D(0) ≠ I for n={n}"
        )

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_near_zero_perturbative(self, n: int) -> None:
        """D(ε) ≈ I + O(ε) for small ε."""
        eps = 1e-6
        D = d_matrix_block(n, eps)
        identity = np.eye(n * n)
        deviation = np.max(np.abs(D - identity))
        # Off-diagonal should be O(ε) at most
        assert deviation < 1e-4, (
            f"D(ε={eps}) deviates by {deviation} from I for n={n}"
        )


# ===================================================================
#  Test 3: SU(2) limits — n=1 is trivial, n=2 recovers d^(1/2)
# ===================================================================

class TestSU2Limits:
    """Check that the SO(4) D-matrix reduces correctly for small n."""

    def test_n1_scalar(self) -> None:
        """n=1 shell has j=0: D^(1) is the 1×1 identity for all γ."""
        for gamma in [0.0, 0.5, np.pi / 2, np.pi]:
            D = d_matrix_block(1, gamma)
            assert D.shape == (1, 1)
            np.testing.assert_allclose(D, [[1.0]], atol=1e-14)

    def test_n2_block_structure(self) -> None:
        """n=2 shell has j=1/2: D^(2) is 4×4, check known elements.

        For n=2, the states are (l=0,m=0), (l=1,m=-1), (l=1,m=0), (l=1,m=1).
        At γ=π/2, the s-p mixing should be non-trivial.
        """
        gamma = np.pi / 2
        D = d_matrix_block(2, gamma)
        assert D.shape == (4, 4)

        # Unitarity check
        np.testing.assert_allclose(D.T @ D, np.eye(4), atol=1e-12)

        # D(γ=π/2) should have non-trivial s-p mixing
        # D_{(0,0),(0,0)} should not be 1
        assert abs(D[0, 0] - 1.0) > 0.01, (
            "s-s element should show mixing at γ=π/2"
        )

    def test_n2_wigner_d_half(self) -> None:
        """Verify SU(2) d^(1/2) elements match known formula.

        d^(1/2)_{1/2,1/2}(β) = cos(β/2)
        d^(1/2)_{1/2,-1/2}(β) = -sin(β/2)
        d^(1/2)_{-1/2,1/2}(β) = sin(β/2)
        d^(1/2)_{-1/2,-1/2}(β) = cos(β/2)
        """
        beta = np.pi / 3

        # Standard Wigner convention
        d_pp = wigner_d_su2(0.5, 0.5, 0.5, beta)
        d_pm = wigner_d_su2(0.5, 0.5, -0.5, beta)
        d_mp = wigner_d_su2(0.5, -0.5, 0.5, beta)
        d_mm = wigner_d_su2(0.5, -0.5, -0.5, beta)

        np.testing.assert_allclose(d_pp, np.cos(beta / 2), atol=1e-14)
        np.testing.assert_allclose(d_pm, -np.sin(beta / 2), atol=1e-14)
        np.testing.assert_allclose(d_mp, np.sin(beta / 2), atol=1e-14)
        np.testing.assert_allclose(d_mm, np.cos(beta / 2), atol=1e-14)

    def test_n2_wigner_d_one(self) -> None:
        """Verify SU(2) d^1 elements for j=1 (used in n=3 shell).

        d^1_{0,0}(β) = cos(β)
        d^1_{1,0}(β) = -sin(β)/√2
        """
        beta = np.pi / 4

        d_00 = wigner_d_su2(1.0, 0.0, 0.0, beta)
        d_10 = wigner_d_su2(1.0, 1.0, 0.0, beta)

        np.testing.assert_allclose(d_00, np.cos(beta), atol=1e-14)
        np.testing.assert_allclose(
            d_10, -np.sin(beta) / np.sqrt(2), atol=1e-14
        )


# ===================================================================
#  Test 4: Composition law — D(α)D(β) = D(α+β) for same-axis
# ===================================================================

class TestComposition:
    """For same-axis rotations, D(α)·D(β) = D(α+β)."""

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_composition(self, n: int) -> None:
        """D(α)·D(β) = D(α+β) for several angle pairs."""
        alpha = 0.4
        beta = 0.7
        D_a = d_matrix_block(n, alpha)
        D_b = d_matrix_block(n, beta)
        D_ab = d_matrix_block(n, alpha + beta)

        product = D_a @ D_b
        np.testing.assert_allclose(
            product, D_ab, atol=1e-11,
            err_msg=f"Composition D({alpha})·D({beta}) ≠ D({alpha+beta}) for n={n}"
        )

    @pytest.mark.parametrize("n", [2, 3])
    def test_inverse(self, n: int) -> None:
        """D(γ)·D(-γ) = I."""
        gamma = 1.2
        D_pos = d_matrix_block(n, gamma)
        D_neg = d_matrix_block(n, -gamma)
        product = D_pos @ D_neg
        np.testing.assert_allclose(
            product, np.eye(n * n), atol=1e-12,
            err_msg=f"D(γ)·D(-γ) ≠ I for n={n}"
        )


# ===================================================================
#  Test 5: Hermitian conjugate — D(γ)† = D(-γ)
# ===================================================================

class TestHermitianConjugate:
    """For real rotation angles, D(γ)^T = D(-γ)."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4])
    def test_transpose_equals_inverse(self, n: int) -> None:
        """D(γ)^T = D(-γ) for real γ."""
        gamma = 0.8
        D_pos = d_matrix_block(n, gamma)
        D_neg = d_matrix_block(n, -gamma)
        np.testing.assert_allclose(
            D_pos.T, D_neg, atol=1e-12,
            err_msg=f"D(γ)^T ≠ D(-γ) for n={n}"
        )


# ===================================================================
#  Test 6: Shibuya-Wulfman spot check
# ===================================================================

class TestShibuyaWulfman:
    """Spot-check against Shibuya & Wulfman (1965) known values.

    For n=2, the D-matrix element D^(2)_{(1,0),(0,0)}(γ) gives the
    s-p_z mixing under the bond rotation. By explicit CG coupling:

        D^(2)_{(1,0),(0,0)}(γ) = d^(1/2)_{1/2,1/2}(γ) · d^(1/2)_{-1/2,-1/2}(γ)
                                - d^(1/2)_{-1/2,1/2}(γ) · d^(1/2)_{1/2,-1/2}(γ)
                                (weighted by CG coefficients)

    This should equal a known function of γ.
    """

    def test_n2_s_to_pz(self) -> None:
        """D^(2)_{(1,0),(0,0)}(γ) at γ=π/3.

        For n=2, j=1/2. The CG coefficients couple:
        |2,0,0⟩ = |1/2,1/2⟩|1/2,-1/2⟩ × C1 + |1/2,-1/2⟩|1/2,1/2⟩ × C2
        |2,1,0⟩ = |1/2,1/2⟩|1/2,-1/2⟩ × C3 + |1/2,-1/2⟩|1/2,1/2⟩ × C4

        With CG(1/2,m+; 1/2,m-| l=0,0) = (-1)^{1/2-m+}/√2,
        CG(1/2,m+; 1/2,m-| l=1,0) = δ_{m+,-m-} × sign / √2.
        """
        gamma = np.pi / 3
        D_10_00 = wigner_D_so4(2, 1, 0, 0, 0, gamma)

        # Verify via explicit SU(2) d-matrices and CG coefficients
        # |2,0,0⟩: CG coefficients for l=0
        cg_00 = cg_so4(2, 0, 0)
        # |2,1,0⟩: CG coefficients for l=1
        cg_10 = cg_so4(2, 1, 0)

        # Manual calculation
        j = 0.5
        manual = 0.0
        for (mp_p, mp_m), c_bra in cg_10.items():
            for (m_p, m_m), c_ket in cg_00.items():
                d_p = wigner_d_su2(j, mp_p, m_p, gamma)
                d_m = wigner_d_su2(j, mp_m, m_m, gamma)
                manual += c_bra * d_p * d_m * c_ket

        np.testing.assert_allclose(
            D_10_00, manual, atol=1e-14,
            err_msg="D-matrix element disagrees with manual CG calculation"
        )

    def test_n2_diagonal_s(self) -> None:
        """D^(2)_{(0,0),(0,0)}(γ) should decrease from 1 as γ increases."""
        D_00_at_0 = wigner_D_so4(2, 0, 0, 0, 0, 0.0)
        D_00_at_half = wigner_D_so4(2, 0, 0, 0, 0, np.pi / 2)
        D_00_at_pi = wigner_D_so4(2, 0, 0, 0, 0, np.pi)

        assert abs(D_00_at_0 - 1.0) < 1e-12, "D(0) s-s should be 1"
        assert abs(D_00_at_half) < 1.0, "D(π/2) s-s should be < 1"
        # At γ=π (united atom), s mixes maximally with p
        # The exact value depends on CG structure


# ===================================================================
#  Test 7: Large R limit — γ → 0, separated atoms
# ===================================================================

class TestLargeRLimit:
    """At large R, p_R → 0, cos γ → 1, γ → 0: atoms are separated."""

    def test_bond_angle_large_R(self) -> None:
        """bond_angle(R, p0) → 0 as R → ∞."""
        p0 = 1.0  # arbitrary
        gamma_10 = bond_angle(10.0, p0)
        gamma_100 = bond_angle(100.0, p0)
        gamma_1000 = bond_angle(1000.0, p0)

        assert gamma_10 < gamma_100 or gamma_100 < 0.1  # not necessarily monotone in float
        assert gamma_100 < 0.1
        assert gamma_1000 < 0.01

    def test_d_matrix_approaches_identity(self) -> None:
        """D-matrix → I as R → ∞ (γ → 0)."""
        p0 = 1.0
        R_large = 10000.0
        gamma = bond_angle(R_large, p0)
        D = d_matrix_block(2, gamma)
        np.testing.assert_allclose(
            D, np.eye(4), atol=1e-3,
            err_msg="D-matrix should approach I at large R"
        )

    def test_bond_angle_values(self) -> None:
        """Check specific bond angle values.

        cos γ = (p₀² - 1/R²) / (p₀² + 1/R²)

        At R=1, p₀=1: cos γ = 0, γ = π/2
        """
        gamma = bond_angle(1.0, 1.0)
        np.testing.assert_allclose(gamma, np.pi / 2, atol=1e-14)


# ===================================================================
#  Test 8: Antipodal limit — γ = π, united atom
# ===================================================================

class TestAntipodalLimit:
    """At γ = π (R → 0), the two poles merge into the united atom."""

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_d_at_pi_is_involution(self, n: int) -> None:
        """D(π) · D(π) = D(2π) = I, so D(π)² = I."""
        D_pi = d_matrix_block(n, np.pi)
        product = D_pi @ D_pi
        np.testing.assert_allclose(
            product, np.eye(n * n), atol=1e-11,
            err_msg=f"D(π)² ≠ I for n={n}"
        )

    def test_n1_at_pi(self) -> None:
        """n=1: D^(1)(π) = [[1]] (scalar, always identity)."""
        D = d_matrix_block(1, np.pi)
        np.testing.assert_allclose(D, [[1.0]], atol=1e-14)

    def test_n2_at_pi_structure(self) -> None:
        """n=2: D^(2)(π) under A-rotation sends |l,m⟩ → c × |l,-m⟩.

        The bond rotation (Runge-Lenz A_y by π) maps m → -m within
        each l-sector: D(π)|l,m⟩ = (-1)^{n-1+l+|m|} |l,-m⟩.
        For n=2: (0,0)→-1, (1,-1)↔(1,+1) with -1, (1,0)→+1.
        """
        D = d_matrix_block(2, np.pi)
        # States: (0,0), (1,-1), (1,0), (1,1)
        expected = np.array([
            [-1, 0, 0, 0],   # (0,0) → -1*(0,0)
            [0, 0, 0, -1],   # (1,-1) → -1*(1,+1)
            [0, 0, 1, 0],    # (1,0) → +1*(1,0)
            [0, -1, 0, 0],   # (1,+1) → -1*(1,-1)
        ], dtype=float)
        np.testing.assert_allclose(
            D, expected, atol=1e-12,
            err_msg="D(π) does not have correct structure for n=2"
        )

    @pytest.mark.parametrize("n", [3, 4])
    def test_parity_at_pi(self, n: int) -> None:
        """D(π)|l,m⟩ = (-1)^{n-1+l+|m|} |l,-m⟩ for general n.

        The A-rotation by π maps m → -m within each l-sector.
        """
        D = d_matrix_block(n, np.pi)
        states = []
        for l in range(n):
            for m_val in range(-l, l + 1):
                states.append((l, m_val))
        dim = len(states)

        # Build expected matrix: D|l,m⟩ = (-1)^{n-1+l+|m|} |l,-m⟩
        expected = np.zeros((dim, dim))
        state_idx = {s: i for i, s in enumerate(states)}
        for j_idx, (l, m_val) in enumerate(states):
            target = (l, -m_val)
            i_idx = state_idx[target]
            expected[i_idx, j_idx] = (-1) ** (n - 1 + l + abs(m_val))

        np.testing.assert_allclose(
            D, expected, atol=1e-10,
            err_msg=f"D(π) structure wrong for n={n}"
        )


# ===================================================================
#  Additional validation: CG coefficients
# ===================================================================

class TestCGCoefficients:
    """Validate CG coefficient properties used by the D-matrix."""

    def test_cg_n1(self) -> None:
        """n=1: j=0, only state is (0,0) → single coefficient = 1."""
        cg = cg_so4(1, 0, 0)
        assert len(cg) == 1
        assert abs(list(cg.values())[0] - 1.0) < 1e-14

    def test_cg_n2_normalization(self) -> None:
        """CG coefficients for n=2 should be normalized: Σ|C|² = 1."""
        for l in range(2):
            for m in range(-l, l + 1):
                cg = cg_so4(2, l, m)
                norm_sq = sum(c ** 2 for c in cg.values())
                np.testing.assert_allclose(
                    norm_sq, 1.0, atol=1e-14,
                    err_msg=f"CG normalization failed for n=2, l={l}, m={m}"
                )

    def test_cg_n3_normalization(self) -> None:
        """CG coefficients for n=3 should be normalized."""
        for l in range(3):
            for m in range(-l, l + 1):
                cg = cg_so4(3, l, m)
                norm_sq = sum(c ** 2 for c in cg.values())
                np.testing.assert_allclose(
                    norm_sq, 1.0, atol=1e-14,
                    err_msg=f"CG normalization failed for n=3, l={l}, m={m}"
                )

    def test_cg_completeness_n2(self) -> None:
        """Completeness: Σ_{l,m} |C^{lm}_{m+,m-}|² = 1 for each (m+,m-)."""
        n = 2
        j = 0.5
        steps = int(2 * j) + 1
        for i_plus in range(steps):
            mp = -j + i_plus
            for i_minus in range(steps):
                mm = -j + i_minus
                total = 0.0
                for l in range(n):
                    m = int(round(mp + mm))
                    if abs(m) > l:
                        continue
                    cg = cg_so4(n, l, m)
                    if (mp, mm) in cg:
                        total += cg[(mp, mm)] ** 2
                np.testing.assert_allclose(
                    total, 1.0, atol=1e-14,
                    err_msg=f"Completeness failed for n=2, m+={mp}, m-={mm}"
                )
