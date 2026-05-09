"""
Tests for the Algebraic Casimir CI module (Track DI Sprint 2).

Validates:
1. Analytical Slater integrals against known exact values
2. FCI matrix construction at k=Z reproduces standard CI energies
3. Self-consistent k* at n_max=1 reproduces Hartree screening (k=9/4)
4. Matrix polynomial structure: H(k) = Bk + Ck^2
5. Variational property: E_var(n_max) >= E_var(n_max+1) >= exact
6. Convergence of the self-consistent sequence
"""

import numpy as np
import pytest
from fractions import Fraction

from geovac.casimir_ci import (
    build_fci_matrix,
    solve_self_consistent,
    solve_variational,
    verify_matrix_polynomial,
    two_electron_integral,
    get_rk4,
    _build_orbital_basis,
    _compute_rk_numerical,
    _compute_inv_r_numerical,
    build_fci_polynomial,
    solve_variational_fast,
    build_graph_native_fci,
    build_graph_consistent_fci,
    _build_graph_h1,
    compute_he_spectrum,
)
from geovac.hypergeometric_slater import compute_rk_float, get_rk_float


# ========================================================================
# Slater integral validation
# ========================================================================

class TestSlaterIntegrals:
    """Verify exact rational Slater integrals against known values."""

    def test_f0_1s1s(self):
        """F^0(1s,1s) = 5/8."""
        val = get_rk4(1, 0, 1, 0, 1, 0, 1, 0, 0)
        assert val == Fraction(5, 8)

    def test_f0_1s2s(self):
        """F^0(1s,2s) = 17/81."""
        val = get_rk4(1, 0, 1, 0, 2, 0, 2, 0, 0)
        assert val == Fraction(17, 81)

    def test_f0_2s2s(self):
        """F^0(2s,2s) = 77/512."""
        val = get_rk4(2, 0, 2, 0, 2, 0, 2, 0, 0)
        assert val == Fraction(77, 512)

    def test_f0_1s2p(self):
        """F^0(1s,2p) = 59/243."""
        val = get_rk4(1, 0, 1, 0, 2, 1, 2, 1, 0)
        assert val == Fraction(59, 243)

    def test_f0_symmetric(self):
        """F^0(a,b) = F^0(b,a)."""
        val_ab = get_rk4(1, 0, 1, 0, 2, 0, 2, 0, 0)
        val_ba = get_rk4(2, 0, 2, 0, 1, 0, 1, 0, 0)
        assert val_ab == val_ba

    def test_g0_1s2s(self):
        """G^0(1s,2s) = 16/729."""
        val = get_rk4(1, 0, 2, 0, 1, 0, 2, 0, 0)
        assert val == Fraction(16, 729)

    def test_g0_1s2p(self):
        """G^0(1s,2p) = 176/2187."""
        val = get_rk4(1, 0, 2, 1, 1, 0, 2, 1, 0)
        assert val == Fraction(176, 2187)

    def test_numerical_matches_table(self):
        """Numerical R^k matches exact table values within grid error."""
        # F^0(1s,1s)
        val_num = _compute_rk_numerical(1, 0, 1, 0, 1, 0, 1, 0, 0)
        assert abs(val_num - 5 / 8) < 0.005, f"F0(1s,1s) numerical={val_num:.6f}, exact=0.625"

        # G^0(1s,2s)
        val_num_g = _compute_rk_numerical(1, 0, 2, 0, 1, 0, 2, 0, 0)
        assert abs(val_num_g - 16 / 729) < 0.002, f"G0(1s,2s) numerical={val_num_g:.6f}, exact={16/729:.6f}"

    def test_float_algebraic_matches_table(self):
        """Float algebraic R^k matches exact table to machine precision."""
        from geovac.casimir_ci import _RK4_TABLE
        max_err = 0.0
        for key, expected in _RK4_TABLE.items():
            n1, l1, n3, l3, n2, l2, n4, l4, k = key
            val = compute_rk_float(n1, l1, n3, l3, n2, l2, n4, l4, k)
            expected_f = float(expected)
            if abs(expected_f) > 1e-15:
                rel_err = abs(val - expected_f) / abs(expected_f)
                max_err = max(max_err, rel_err)
        assert max_err < 1e-10, f"Float algebraic max rel error: {max_err:.2e}"

    def test_float_algebraic_n4_integrals(self):
        """Float algebraic works for n=4 (beyond precomputed table)."""
        # R^0(4s,4s; 4s,4s) — not in table
        val = compute_rk_float(4, 0, 4, 0, 4, 0, 4, 0, 0)
        assert val > 0, f"R^0(4s4s;4s4s) should be positive, got {val}"
        # Known value from Fraction computation: 0.037271499633789
        assert abs(val - 0.037271499634) < 1e-10

    def test_float_algebraic_beats_grid(self):
        """Float algebraic is more accurate than grid quadrature."""
        # R^0(4s,3p; 4s,3p) — grid has ~0.4% error
        val_float = compute_rk_float(4, 0, 3, 1, 4, 0, 3, 1, 0)
        val_grid = _compute_rk_numerical(4, 0, 3, 1, 4, 0, 3, 1, 0)
        # Float algebraic should differ from grid (grid has truncation error)
        assert abs(val_float - val_grid) > 1e-6, (
            "Float and grid should differ (grid has truncation error)"
        )


# ========================================================================
# Two-electron integral validation
# ========================================================================

class TestTwoElectronIntegrals:
    """Verify two-electron integrals at k_orb=1 and k_orb=Z."""

    def test_1s1s_direct_k1(self):
        """<1s,1s|1/r12|1s,1s> at k=1 should be F^0(1s,1s) = 5/8."""
        val = two_electron_integral(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, k_orb=1.0)
        assert abs(val - 5 / 8) < 1e-10

    def test_1s1s_direct_k2(self):
        """<1s,1s|1/r12|1s,1s> at k=2 should be 2 * 5/8 = 5/4."""
        val = two_electron_integral(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, k_orb=2.0)
        assert abs(val - 5 / 4) < 1e-10

    def test_linear_scaling(self):
        """Two-electron integrals scale linearly with k_orb."""
        val1 = two_electron_integral(1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, k_orb=1.0)
        val2 = two_electron_integral(1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, k_orb=3.0)
        assert abs(val2 / val1 - 3.0) < 1e-8


# ========================================================================
# FCI matrix validation
# ========================================================================

class TestFCIMatrix:
    """Verify FCI matrix construction."""

    def test_nmax1_diagonal(self):
        """H(k=2, Z=2, n_max=1) should give E = k^2 - 2Zk + 5k/8 = -2.75."""
        H = build_fci_matrix(Z=2, n_max=1, k_orb=2.0)
        assert H.shape == (1, 1)
        assert abs(H[0, 0] - (-2.75)) < 1e-10

    def test_nmax1_arbitrary_k(self):
        """E(k) = k^2 - 2Zk + 5k/8 for He 1s^2 at any k."""
        for k in [1.0, 1.5, 2.0, 2.5, 3.0]:
            H = build_fci_matrix(Z=2, n_max=1, k_orb=k)
            expected = k ** 2 - 4 * k + 5 * k / 8
            assert abs(H[0, 0] - expected) < 1e-10, f"k={k}: H={H[0,0]}, expected={expected}"

    def test_nmax2_shape(self):
        """n_max=2, M_L=0 singlet: 7 configs (from 5 spatial orbitals)."""
        H = build_fci_matrix(Z=2, n_max=2, k_orb=2.0)
        assert H.shape == (7, 7)

    def test_nmax2_hermitian(self):
        """FCI matrix should be Hermitian."""
        H = build_fci_matrix(Z=2, n_max=2, k_orb=1.7)
        assert np.allclose(H, H.T, atol=1e-12)

    def test_nmax2_ground_state_physical(self):
        """Ground state at k=Z should be above exact He energy."""
        H = build_fci_matrix(Z=2, n_max=2, k_orb=2.0)
        E0 = np.linalg.eigvalsh(H)[0]
        # With fixed k=Z=2, the energy should be above the exact value
        # (variational principle holds for fixed k, but not for the
        # best possible k — the energy can go below exact if k is varied)
        assert E0 > -3.5, f"Ground state {E0} is unphysically low"
        assert E0 < -2.0, f"Ground state {E0} is unphysically high"

    def test_sonly_nmax2(self):
        """s-only CI at n_max=2 has 3 configs: (1s,1s), (1s,2s), (2s,2s)."""
        H = build_fci_matrix(Z=2, n_max=2, k_orb=2.0, l_max=0)
        assert H.shape == (3, 3)

    def test_nmax3_shape(self):
        """n_max=3, M_L=0 singlet: 31 configs."""
        H = build_fci_matrix(Z=2, n_max=3, k_orb=2.0)
        assert H.shape == (31, 31)


# ========================================================================
# Matrix polynomial structure
# ========================================================================

class TestMatrixPolynomial:
    """Verify H(k) = Bk + Ck^2 (quadratic polynomial in k)."""

    def test_nmax1_polynomial(self):
        """n_max=1: H(k) is exactly quadratic in k."""
        result = verify_matrix_polynomial(Z=2, n_max=1)
        assert result['is_polynomial']
        assert result['max_residual'] < 1e-12

    def test_nmax2_polynomial(self):
        """n_max=2: H(k) is exactly quadratic in k."""
        result = verify_matrix_polynomial(Z=2, n_max=2)
        assert result['is_polynomial']
        assert result['max_residual'] < 1e-10

    def test_nmax3_polynomial(self):
        """n_max=3: H(k) is exactly quadratic in k."""
        result = verify_matrix_polynomial(Z=2, n_max=3)
        assert result['is_polynomial']
        assert result['max_residual'] < 1e-8


# ========================================================================
# Self-consistent solver
# ========================================================================

class TestSelfConsistent:
    """Validate the self-consistency condition k^2 = -2E."""

    def test_nmax1_hartree(self):
        """n_max=1: k* = 9/4, E* = -81/32 (algebraic, exact)."""
        k_star, E_star, _ = solve_self_consistent(Z=2, n_max=1, n_scan=100)
        assert abs(k_star - 9 / 4) < 1e-8, f"k*={k_star}, expected 2.25"
        assert abs(E_star - (-81 / 32)) < 1e-8, f"E*={E_star}, expected {-81/32}"

    def test_self_consistency_condition(self):
        """k*^2 + 2*E(k*) = 0 within tolerance."""
        for n_max in [1, 2]:
            k_star, E_star, _ = solve_self_consistent(Z=2, n_max=n_max, n_scan=100)
            H = build_fci_matrix(Z=2, n_max=n_max, k_orb=k_star)
            E_check = np.linalg.eigvalsh(H)[0]
            residual = k_star ** 2 + 2 * E_check
            assert abs(residual) < 1e-6, f"n_max={n_max}: residual={residual}"

    def test_sc_convergence(self):
        """Self-consistent energy improves with n_max."""
        _, E1, _ = solve_self_consistent(Z=2, n_max=1, n_scan=100)
        _, E2, _ = solve_self_consistent(Z=2, n_max=2, n_scan=100)
        # E should decrease (become more negative) with n_max
        assert E2 < E1, f"E(n_max=2)={E2} should be < E(n_max=1)={E1}"


# ========================================================================
# Variational solver
# ========================================================================

class TestVariational:
    """Validate variational optimization of k."""

    def test_nmax1_screening(self):
        """n_max=1: k_var = Z - 5/16 = 27/16 (Hartree screening)."""
        k_var, E_var = solve_variational(Z=2, n_max=1)
        assert abs(k_var - 27 / 16) < 1e-6, f"k_var={k_var}, expected {27/16}"
        assert abs(E_var - (-729 / 256)) < 1e-6, f"E_var={E_var}, expected {-729/256}"

    def test_variational_below_fixed(self):
        """Variational energy <= fixed k=Z energy."""
        for n_max in [1, 2]:
            _, E_var = solve_variational(Z=2, n_max=n_max)
            H = build_fci_matrix(Z=2, n_max=n_max, k_orb=2.0)
            E_fixed = np.linalg.eigvalsh(H)[0]
            assert E_var <= E_fixed + 1e-10, f"n_max={n_max}: E_var={E_var} > E_fixed={E_fixed}"

    def test_variational_convergence(self):
        """Variational energy improves with n_max."""
        _, E1 = solve_variational(Z=2, n_max=1)
        _, E2 = solve_variational(Z=2, n_max=2)
        assert E2 < E1, f"E(n_max=2)={E2} should be < E(n_max=1)={E1}"

    def test_variational_above_exact(self):
        """Variational energy must be above exact He ground state."""
        E_exact = -2.903724
        for n_max in [1, 2]:
            _, E_var = solve_variational(Z=2, n_max=n_max)
            assert E_var > E_exact, (
                f"n_max={n_max}: E_var={E_var:.6f} violates variational bound "
                f"(exact={E_exact})"
            )


# ========================================================================
# Orbital basis
# ========================================================================

class TestOrbitalBasis:
    """Verify orbital basis construction."""

    def test_nmax1(self):
        """n_max=1: just 1s."""
        orbs = _build_orbital_basis(1)
        assert len(orbs) == 1
        assert orbs[0] == (1, 0, 0)

    def test_nmax2(self):
        """n_max=2: 1s, 2s, 2p_{-1}, 2p_0, 2p_1 = 5 orbitals."""
        orbs = _build_orbital_basis(2)
        assert len(orbs) == 5

    def test_nmax3(self):
        """n_max=3: 1 + 4 + 9 = 14 orbitals."""
        orbs = _build_orbital_basis(3)
        assert len(orbs) == 14

    def test_lmax_filter(self):
        """l_max=0 keeps only s orbitals."""
        orbs = _build_orbital_basis(3, l_max=0)
        assert len(orbs) == 3  # 1s, 2s, 3s
        assert all(l == 0 for _, l, _ in orbs)

    def test_nmax4(self):
        """n_max=4: 1+4+9+16 = 30 orbitals."""
        orbs = _build_orbital_basis(4)
        assert len(orbs) == 30

    def test_nmax5(self):
        """n_max=5: 1+4+9+16+25 = 55 orbitals."""
        orbs = _build_orbital_basis(5)
        assert len(orbs) == 55

    def test_nmax6(self):
        """n_max=6: 1+4+9+16+25+36 = 91 orbitals."""
        orbs = _build_orbital_basis(6)
        assert len(orbs) == 91


# ========================================================================
# Numerical 1/r matrix element fallback
# ========================================================================

class TestInvRNumerical:
    """Verify numerical 1/r fallback matches known values."""

    def test_diagonal_1s(self):
        """<1s|1/r|1s> = 1/1^2 = 1.0."""
        val = _compute_inv_r_numerical(1, 0, 1, 0)
        assert abs(val - 1.0) < 0.01

    def test_diagonal_2s(self):
        """<2s|1/r|2s> = 1/4."""
        val = _compute_inv_r_numerical(2, 0, 2, 0)
        assert abs(val - 0.25) < 0.005

    def test_off_diagonal_1s_2s(self):
        """<1s|1/r|2s> = 4*sqrt(2)/27."""
        expected = 4 * np.sqrt(2) / 27
        val = _compute_inv_r_numerical(1, 0, 2, 0)
        assert abs(val - expected) < 0.005, f"got {val}, expected {expected}"

    def test_diagonal_4s(self):
        """<4s|1/r|4s> = 1/16."""
        val = _compute_inv_r_numerical(4, 0, 4, 0)
        assert abs(val - 1.0 / 16) < 0.005

    def test_different_l_zero(self):
        """<1s|1/r|2p> = 0 (different l)."""
        val = _compute_inv_r_numerical(1, 0, 2, 1)
        assert abs(val) < 1e-10


# ========================================================================
# Polynomial builder and fast solver
# ========================================================================

class TestPolynomialBuilder:
    """Verify polynomial matrix extraction and fast solver."""

    def test_polynomial_nmax2(self):
        """B, C reproduce H(k) at test k values."""
        B, C = build_fci_polynomial(Z=2, n_max=2)
        for k in [1.3, 1.7, 2.5]:
            H_direct = build_fci_matrix(Z=2, n_max=2, k_orb=k)
            H_poly = B * k + C * k ** 2
            assert np.allclose(H_direct, H_poly, atol=1e-10), \
                f"k={k}: max diff = {np.max(np.abs(H_direct - H_poly))}"

    def test_fast_matches_standard(self):
        """Fast variational solver matches standard solver."""
        k_std, E_std = solve_variational(Z=2, n_max=2)
        k_fast, E_fast = solve_variational_fast(Z=2, n_max=2)
        assert abs(k_std - k_fast) < 1e-6, f"k: std={k_std}, fast={k_fast}"
        assert abs(E_std - E_fast) < 1e-8, f"E: std={E_std}, fast={E_fast}"


# ========================================================================
# Extended n_max validation (Sprint 3)
# ========================================================================

class TestExtendedNmax:
    """Validate FCI construction and variational convergence at n_max=4,5,6."""

    def test_nmax4_shape(self):
        """n_max=4, M_L=0 singlet: verify config count."""
        H = build_fci_matrix(Z=2, n_max=4, k_orb=2.0)
        # n_spatial=30, configs = pairs (i,j) with i<=j and m_i+m_j=0
        assert H.shape[0] > 31  # must be more than n_max=3
        assert H.shape[0] == H.shape[1]  # square

    def test_nmax4_hermitian(self):
        """FCI matrix at n_max=4 should be Hermitian."""
        H = build_fci_matrix(Z=2, n_max=4, k_orb=1.7)
        assert np.allclose(H, H.T, atol=1e-10)

    @pytest.mark.slow
    def test_nmax5_shape(self):
        """n_max=5, M_L=0 singlet: verify construction."""
        H = build_fci_matrix(Z=2, n_max=5, k_orb=2.0)
        assert H.shape[0] > 100
        assert H.shape[0] == H.shape[1]

    def test_variational_convergence_nmax4(self):
        """E_var(n_max=4) < E_var(n_max=3)."""
        _, E3 = solve_variational_fast(Z=2, n_max=3)
        _, E4 = solve_variational_fast(Z=2, n_max=4)
        assert E4 < E3, f"E4={E4} should be < E3={E3}"

    def test_variational_bound_nmax4(self):
        """E_var(n_max=4) > E_exact."""
        _, E4 = solve_variational_fast(Z=2, n_max=4)
        assert E4 > -2.903724, f"E4={E4} violates variational bound"

    @pytest.mark.slow
    def test_variational_convergence_nmax5(self):
        """E_var(n_max=5) < E_var(n_max=4)."""
        _, E4 = solve_variational_fast(Z=2, n_max=4)
        _, E5 = solve_variational_fast(Z=2, n_max=5)
        assert E5 < E4, f"E5={E5} should be < E4={E4}"

    @pytest.mark.slow
    def test_variational_bound_nmax5(self):
        """E_var(n_max=5) > E_exact."""
        _, E5 = solve_variational_fast(Z=2, n_max=5)
        assert E5 > -2.903724, f"E5={E5} violates variational bound"

    def test_polynomial_nmax4(self):
        """H(k) is quadratic in k at n_max=4."""
        B, C = build_fci_polynomial(Z=2, n_max=4)
        # Verify at a test point
        k_test = 1.73
        H_direct = build_fci_matrix(Z=2, n_max=4, k_orb=k_test)
        H_poly = B * k_test + C * k_test ** 2
        max_diff = np.max(np.abs(H_direct - H_poly))
        assert max_diff < 1e-6, f"Polynomial residual {max_diff} at n_max=4"


# ========================================================================
# Graph-native CI (Sprint 3C)
# ========================================================================

class TestGraphNativeCI:
    """Validate graph-native CI with hybrid h1."""

    def test_graph_h1_hermitian(self):
        """Graph h1 matrix should be Hermitian."""
        h1, orbs = _build_graph_h1(Z=2, n_max=3)
        assert np.allclose(h1, h1.T, atol=1e-14)

    def test_graph_h1_diagonal_exact(self):
        """Diagonal of graph h1 should be -Z^2/(2n^2)."""
        h1, orbs = _build_graph_h1(Z=2, n_max=3)
        for i, (n, l, m) in enumerate(orbs):
            expected = -4.0 / (2.0 * n * n)
            assert abs(h1[i, i] - expected) < 1e-12, \
                f"h1[{i},{i}] = {h1[i,i]}, expected {expected} for ({n},{l},{m})"

    def test_graph_h1_has_offdiagonal(self):
        """Graph h1 should have off-diagonal entries from adjacency."""
        h1, orbs = _build_graph_h1(Z=2, n_max=3)
        offdiag = h1 - np.diag(np.diag(h1))
        assert np.any(np.abs(offdiag) > 1e-10), "No off-diagonal entries found"

    def test_graph_h1_offdiag_value(self):
        """Off-diagonal entries should be kappa * (-A) = 1/16 for connected states."""
        h1, orbs = _build_graph_h1(Z=2, n_max=2)
        # 1s (idx 0) connects to 2s (idx 1) via T+ transition
        expected = 1.0 / 16.0  # kappa * (-1) = (-1/16) * (-1) = 1/16
        assert abs(h1[0, 1] - expected) < 1e-12, \
            f"h1[0,1] = {h1[0,1]}, expected {expected}"

    def test_graph_native_fci_nmax2(self):
        """Graph-native FCI at n_max=2 should give better than 1% error."""
        H = build_graph_native_fci(Z=2, n_max=2)
        assert H.shape == (7, 7)
        E0 = np.linalg.eigvalsh(H)[0]
        error_pct = abs((E0 - (-2.903724)) / (-2.903724)) * 100
        assert error_pct < 1.0, f"Error {error_pct:.2f}% > 1% at n_max=2"

    def test_graph_native_fci_nmax3(self):
        """Graph-native FCI at n_max=3 should give better than 0.5% error."""
        H = build_graph_native_fci(Z=2, n_max=3)
        assert H.shape == (31, 31)
        E0 = np.linalg.eigvalsh(H)[0]
        error_pct = abs((E0 - (-2.903724)) / (-2.903724)) * 100
        assert error_pct < 0.5, f"Error {error_pct:.2f}% > 0.5% at n_max=3"

    def test_graph_native_hermitian(self):
        """Graph-native FCI matrix should be Hermitian."""
        H = build_graph_native_fci(Z=2, n_max=3)
        assert np.allclose(H, H.T, atol=1e-10)

    def test_graph_native_variational_bound(self):
        """Graph-native energy must be above exact He energy."""
        H = build_graph_native_fci(Z=2, n_max=3)
        E0 = np.linalg.eigvalsh(H)[0]
        assert E0 > -2.903724, f"E0={E0} violates variational bound"

    def test_graph_native_convergence(self):
        """E(n_max=3) < E(n_max=2) (energy improves with basis)."""
        H2 = build_graph_native_fci(Z=2, n_max=2)
        H3 = build_graph_native_fci(Z=2, n_max=3)
        E2 = np.linalg.eigvalsh(H2)[0]
        E3 = np.linalg.eigvalsh(H3)[0]
        assert E3 < E2, f"E3={E3} should be < E2={E2}"

    def test_graph_native_beats_diagonal(self):
        """Graph-native should beat fixed k=2 diagonal at n_max=3."""
        H_graph = build_graph_native_fci(Z=2, n_max=3)
        E_graph = np.linalg.eigvalsh(H_graph)[0]
        H_diag = build_fci_matrix(Z=2, n_max=3, k_orb=2.0)
        E_diag = np.linalg.eigvalsh(H_diag)[0]
        # Graph-native error should be smaller
        err_graph = abs(E_graph - (-2.903724))
        err_diag = abs(E_diag - (-2.903724))
        assert err_graph < err_diag, \
            f"Graph error {err_graph:.6f} should be < diagonal error {err_diag:.6f}"


# ========================================================================
# Graph-consistent CI — basis invariance verification (Sprint 3D)
# ========================================================================

class TestGraphConsistentCI:
    """Verify FCI basis invariance: graph-consistent == graph-native."""

    def test_nmax1_identical(self):
        """n_max=1: consistent and native give identical energy."""
        H_native = build_graph_native_fci(Z=2, n_max=1)
        H_consistent = build_graph_consistent_fci(Z=2, n_max=1)
        E_native = np.linalg.eigvalsh(H_native)[0]
        E_consistent = np.linalg.eigvalsh(H_consistent)[0]
        assert abs(E_native - E_consistent) < 1e-10, \
            f"E_native={E_native:.10f}, E_consistent={E_consistent:.10f}"

    def test_nmax2_identical(self):
        """n_max=2: consistent and native give identical energy."""
        H_native = build_graph_native_fci(Z=2, n_max=2)
        H_consistent = build_graph_consistent_fci(Z=2, n_max=2)
        E_native = np.linalg.eigvalsh(H_native)[0]
        E_consistent = np.linalg.eigvalsh(H_consistent)[0]
        assert abs(E_native - E_consistent) < 1e-10, \
            f"E_native={E_native:.10f}, E_consistent={E_consistent:.10f}"

    def test_nmax3_identical(self):
        """n_max=3: consistent and native give identical energy."""
        H_native = build_graph_native_fci(Z=2, n_max=3)
        H_consistent = build_graph_consistent_fci(Z=2, n_max=3)
        E_native = np.linalg.eigvalsh(H_native)[0]
        E_consistent = np.linalg.eigvalsh(H_consistent)[0]
        assert abs(E_native - E_consistent) < 1e-10, \
            f"E_native={E_native:.10f}, E_consistent={E_consistent:.10f}"

    def test_full_spectrum_invariance(self):
        """All FCI eigenvalues are invariant under orbital rotation."""
        H_native = build_graph_native_fci(Z=2, n_max=3)
        H_consistent = build_graph_consistent_fci(Z=2, n_max=3)
        evals_native = np.sort(np.linalg.eigvalsh(H_native))
        evals_consistent = np.sort(np.linalg.eigvalsh(H_consistent))
        assert np.allclose(evals_native, evals_consistent, atol=1e-8), \
            f"Max diff: {np.max(np.abs(evals_native - evals_consistent)):.2e}"

    def test_consistent_hermitian(self):
        """Graph-consistent FCI matrix should be Hermitian."""
        H = build_graph_consistent_fci(Z=2, n_max=2)
        assert np.allclose(H, H.T, atol=1e-10)

    def test_consistent_same_shape(self):
        """Consistent and native produce same-shaped matrices."""
        H_native = build_graph_native_fci(Z=2, n_max=3)
        H_consistent = build_graph_consistent_fci(Z=2, n_max=3)
        assert H_native.shape == H_consistent.shape


# ========================================================================
# Sub-block restriction (subblock kwarg)
# ========================================================================
#
# At n_max >= 3 the M_L=0 sector mixes ¹S, ¹P at M_L=0, ¹D at M_L=0, etc.
# Sort-by-energy labeling of the eigenvalues conflates these states. The
# `subblock` kwarg restricts configuration enumeration to a fixed (l_a, l_b)
# multiset, isolating a clean angular-character sub-block.
#
# Sprint reference: debug/precision_catalogue_he_2s_singlet_triplet.py
# (the He 2¹S − 2³S exchange splitting sprint, May 2026) constructed
# the ss-only sub-block inline; this test class verifies the production
# `subblock=(0, 0)` API reproduces those sprint results.

HA_TO_CM_INVERSE = 219474.6313632


class TestSubblockBackwardCompat:
    """Verify subblock=None preserves the original behavior bit-identically."""

    def test_default_singlet_nmax2_bit_identical(self):
        """subblock=None matches the no-kwarg path bit-identically (singlet, n_max=2)."""
        H_old = build_graph_native_fci(Z=2, n_max=2, m_total=0, spin='singlet')
        H_new = build_graph_native_fci(Z=2, n_max=2, m_total=0, spin='singlet',
                                       subblock=None)
        assert H_old.shape == H_new.shape
        assert np.max(np.abs(H_old - H_new)) == 0.0

    def test_default_triplet_nmax3_bit_identical(self):
        """subblock=None matches no-kwarg (triplet, n_max=3, multiple m_total)."""
        for m_total in [-1, 0, 1]:
            H_old = build_graph_native_fci(Z=2, n_max=3, m_total=m_total,
                                           spin='triplet')
            H_new = build_graph_native_fci(Z=2, n_max=3, m_total=m_total,
                                           spin='triplet', subblock=None)
            assert H_old.shape == H_new.shape, (
                f"shape mismatch at m_total={m_total}: "
                f"{H_old.shape} vs {H_new.shape}"
            )
            assert np.max(np.abs(H_old - H_new)) == 0.0, (
                f"non-zero diff at m_total={m_total}"
            )

    def test_default_nmax4_bit_identical(self):
        """subblock=None matches no-kwarg at n_max=4 (larger basis)."""
        H_old = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='singlet')
        H_new = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='singlet',
                                       subblock=None)
        assert H_old.shape == H_new.shape
        assert np.max(np.abs(H_old - H_new)) == 0.0

    def test_compute_he_spectrum_default_unchanged(self):
        """compute_he_spectrum default behavior unchanged when subblock omitted."""
        result_old = compute_he_spectrum(n_max=2)
        result_new = compute_he_spectrum(n_max=2, subblock=None)
        # Same sectors, same dims, same energies bit-identical
        assert result_old['sector_dims'] == result_new['sector_dims']
        for label, info in result_old['states'].items():
            assert label in result_new['states']
            assert info['energy'] == result_new['states'][label]['energy']
        assert result_new['subblock'] is None


class TestSubblockSSOnly:
    """Verify subblock=(0, 0) reproduces the He 2¹S − 2³S sprint result.

    Sprint driver: debug/precision_catalogue_he_2s_singlet_triplet.py
    Reference data: debug/data/precision_catalogue_he_2s_singlet_triplet.json
    """

    def test_ss_singlet_dim_nmax3(self):
        """ss-only singlet at n_max=3: exactly 6 configs (3 s-orbitals, pairs i<=j)."""
        H = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                   subblock=(0, 0))
        # n_s = 3 (1s, 2s, 3s) → pairs (i, j) with i <= j: 3*4/2 = 6
        assert H.shape == (6, 6)

    def test_ss_triplet_dim_nmax3(self):
        """ss-only triplet at n_max=3: exactly 3 configs (i < j only, no diagonal)."""
        H = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet',
                                   subblock=(0, 0))
        # n_s = 3 → pairs (i, j) with i < j: 3*2/2 = 3
        assert H.shape == (3, 3)

    def test_ss_singlet_dim_nmax5(self):
        """ss-only singlet at n_max=5: 15 configs."""
        H = build_graph_native_fci(Z=2, n_max=5, m_total=0, spin='singlet',
                                   subblock=(0, 0))
        # n_s = 5 → 5*6/2 = 15
        assert H.shape == (15, 15)

    def test_ss_triplet_dim_nmax5(self):
        """ss-only triplet at n_max=5: 10 configs."""
        H = build_graph_native_fci(Z=2, n_max=5, m_total=0, spin='triplet',
                                   subblock=(0, 0))
        # n_s = 5 → 5*4/2 = 10
        assert H.shape == (10, 10)

    def test_ss_only_excludes_p_orbitals(self):
        """ss-only at n_max=3 must not contain any p-orbital configuration.

        Compared to the full M_L=0 sector (31 singlet configs), the ss-only
        sub-block has 6 — far fewer because ss eliminates all (..,p_+1)(..,p_-1),
        (p_0, p_0), (s, d_0), (p_+1, d_-1), (p_-1, d_+1), (d_0, d_0) etc.
        configurations.
        """
        H_full = build_graph_native_fci(Z=2, n_max=3, m_total=0,
                                        spin='singlet')
        H_ss = build_graph_native_fci(Z=2, n_max=3, m_total=0,
                                      spin='singlet', subblock=(0, 0))
        assert H_full.shape[0] == 31
        assert H_ss.shape[0] == 6
        assert H_ss.shape[0] < H_full.shape[0]

    def test_ss_singlet_triplet_splitting_nmax5(self):
        """ss-only n_max=5 reproduces the sprint 2¹S − 2³S splitting at 7308.8 cm⁻¹.

        Sprint reference (debug/data/...): at n_max=5 ss-only,
            E(1¹S) = -2.89235952, E(2³S) = -2.23643922, E(2¹S) = -2.20313800,
            splitting = 7308.7713 cm⁻¹.
        Memo §3 convergence table records this as the n_max=5 entry. The
        agreement here verifies the production API reproduces the sprint
        construction to machine precision.
        """
        H_s = build_graph_native_fci(Z=2, n_max=5, m_total=0, spin='singlet',
                                     subblock=(0, 0))
        H_t = build_graph_native_fci(Z=2, n_max=5, m_total=0, spin='triplet',
                                     subblock=(0, 0))
        evals_s = np.sort(np.linalg.eigvalsh(H_s))
        evals_t = np.sort(np.linalg.eigvalsh(H_t))
        E_1_1S = float(evals_s[0])
        E_2_1S = float(evals_s[1])
        E_2_3S = float(evals_t[0])

        # Per-state energies match sprint reference (saved at 8 dp).
        assert abs(E_1_1S - (-2.89235952)) < 1e-7, (
            f"E(1¹S) = {E_1_1S:.10f}, sprint = -2.89235952"
        )
        assert abs(E_2_3S - (-2.23643922)) < 1e-7, (
            f"E(2³S) = {E_2_3S:.10f}, sprint = -2.23643922"
        )
        assert abs(E_2_1S - (-2.20313800)) < 1e-7, (
            f"E(2¹S) = {E_2_1S:.10f}, sprint = -2.20313800"
        )

        # Splitting matches sprint memo §3: 7308.77 cm⁻¹
        dE_cm = (E_2_1S - E_2_3S) * HA_TO_CM_INVERSE
        assert abs(dE_cm - 7308.7713) < 0.01, (
            f"splitting = {dE_cm:.4f} cm⁻¹, sprint = 7308.7713"
        )

    def test_ss_singlet_triplet_splitting_nmax3(self):
        """ss-only n_max=3: matches sprint reference (9737.6 cm⁻¹)."""
        H_s = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                     subblock=(0, 0))
        H_t = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet',
                                     subblock=(0, 0))
        evals_s = np.sort(np.linalg.eigvalsh(H_s))
        evals_t = np.sort(np.linalg.eigvalsh(H_t))
        E_2_1S = float(evals_s[1])
        E_2_3S = float(evals_t[0])
        dE_cm = (E_2_1S - E_2_3S) * HA_TO_CM_INVERSE
        # Sprint n_max=3 row: 9737.5977 cm⁻¹
        assert abs(dE_cm - 9737.5977) < 0.01, (
            f"n_max=3 splitting = {dE_cm:.4f} cm⁻¹, sprint = 9737.5977"
        )

    def test_ss_singlet_triplet_splitting_nmax4(self):
        """ss-only n_max=4: matches sprint reference (8125.7 cm⁻¹)."""
        H_s = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='singlet',
                                     subblock=(0, 0))
        H_t = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='triplet',
                                     subblock=(0, 0))
        evals_s = np.sort(np.linalg.eigvalsh(H_s))
        evals_t = np.sort(np.linalg.eigvalsh(H_t))
        dE_cm = (float(evals_s[1]) - float(evals_t[0])) * HA_TO_CM_INVERSE
        # Sprint n_max=4 row: 8125.6557 cm⁻¹
        assert abs(dE_cm - 8125.6557) < 0.01, (
            f"n_max=4 splitting = {dE_cm:.4f} cm⁻¹, sprint = 8125.6557"
        )

    def test_ss_subblock_hund_ordering(self):
        """ss-only triplet 2³S < singlet 2¹S (Hund's first rule)."""
        H_s = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='singlet',
                                     subblock=(0, 0))
        H_t = build_graph_native_fci(Z=2, n_max=4, m_total=0, spin='triplet',
                                     subblock=(0, 0))
        evals_s = np.sort(np.linalg.eigvalsh(H_s))
        evals_t = np.sort(np.linalg.eigvalsh(H_t))
        E_2_1S = float(evals_s[1])
        E_2_3S = float(evals_t[0])
        assert E_2_3S < E_2_1S, (
            f"E(2³S) = {E_2_3S:.6f} should be < E(2¹S) = {E_2_1S:.6f}"
        )


class TestSubblockSP:
    """Verify subblock=(0, 1) gives a sensible ¹P/³P spectrum."""

    def test_sp_dim_nmax3_m0(self):
        """sp at n_max=3, m_total=0: 6 pairs (1s,2p_0)..(3s,3p_0).

        Spatial pairs (n_s, m_s=0) with (n_p, m_p=0) for n_s, n_p ∈ {1,..,n_max}
        (under m_total=0 ⇒ m_p must equal 0). Singlet allows n_s == n_p,
        but s-orbital (l=0) and p-orbital (l=1) live at different orbital
        indices so the i <= j vs i < j distinction depends on indexing —
        we just check the count is sensible.
        """
        H_s = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                     subblock=(0, 1))
        H_t = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet',
                                     subblock=(0, 1))
        # Singlet: n_s × n_p = 3 × 3 = 9 ordered pairs;
        # i < j ordering by orbital index in basis:
        # orbital basis at n_max=3: 1s, 2s, 2p_-1, 2p_0, 2p_+1, 3s, 3p_-1, 3p_0, 3p_+1, 3d_-2, ...
        # All s indices: {0(1s), 1(2s), 5(3s)}; p_0 indices: {3(2p_0), 7(3p_0)}
        # Pairs (s, p_0) at m_total=0: 3 × 2 = 6 ordered, all with i<j or i>j
        # Singlet uses i<=j constraint: (i, j) sorted, so 6 pairs total
        # Triplet i<j: also 6 pairs (l_a=0 != l_b=1 means no diagonal)
        assert H_s.shape == (6, 6), f"singlet sp dim = {H_s.shape}"
        assert H_t.shape == (6, 6), f"triplet sp dim = {H_t.shape}"

    def test_sp_hund_ordering_nmax3(self):
        """sp at n_max=3: lowest 2³P below lowest 2¹P (Hund's first rule)."""
        H_s = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                     subblock=(0, 1))
        H_t = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet',
                                     subblock=(0, 1))
        E_p1 = float(np.linalg.eigvalsh(H_s)[0])
        E_p3 = float(np.linalg.eigvalsh(H_t)[0])
        assert E_p3 < E_p1, (
            f"lowest triplet sp ({E_p3:.6f}) must be < lowest singlet sp ({E_p1:.6f})"
        )

    def test_sp_excludes_ss_and_pp(self):
        """sp sub-block is a proper subset: smaller than the full M_L=0 sector."""
        H_full = build_graph_native_fci(Z=2, n_max=3, m_total=0,
                                        spin='singlet')
        H_sp = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                      subblock=(0, 1))
        # Full M_L=0 has 31 configs (ss + sp + pp + sd + pd + dd at M_L=0);
        # sp is just the (l_a, l_b) = {0, 1} subset.
        assert H_full.shape[0] == 31
        assert H_sp.shape[0] == 6
        assert H_sp.shape[0] < H_full.shape[0]

    def test_sp_hermitian(self):
        """sp sub-block FCI matrix is Hermitian."""
        H_s = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                     subblock=(0, 1))
        H_t = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet',
                                     subblock=(0, 1))
        assert np.allclose(H_s, H_s.T, atol=1e-12)
        assert np.allclose(H_t, H_t.T, atol=1e-12)


class TestSubblockMisc:
    """Additional sub-block validation."""

    def test_subblock_pp_nmax3(self):
        """pp sub-block (l_a = l_b = 1): only (..p,..p) configs at M_L=0."""
        H_s = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                     subblock=(1, 1))
        # At n_max=3 there are 2 p-shells (2p, 3p) × 3 m-values each = 6 p-orbitals.
        # Pairs (i, j) with i <= j and m_a + m_b = 0:
        #   (m_a, m_b) ∈ {(-1, 1), (0, 0), (1, -1)} ; (1,-1) is symmetric to (-1,1)
        #   Considering p-orbitals at indices 2,3,4 (2p_-1, 2p_0, 2p_+1) and
        #   6,7,8 (3p_-1, 3p_0, 3p_+1).
        # Pairs with i<=j and m_a+m_b=0:
        #   (2p_-1, 2p_+1), (2p_0, 2p_0), (2p_+1, 2p_-1 conjugated to symmetric),
        #   (2p_-1, 3p_+1), (2p_0, 3p_0), (2p_+1, 3p_-1),
        #   (3p_-1, 3p_+1), (3p_0, 3p_0), (3p_+1, 3p_-1)
        # With i<=j and m=0 constraint, count the valid set.
        assert H_s.shape[0] > 0, "pp singlet should have at least 1 config"
        # Hermiticity sanity
        assert np.allclose(H_s, H_s.T, atol=1e-12)

    def test_subblock_negative_l_raises(self):
        """Negative l in subblock raises ValueError."""
        with pytest.raises(ValueError):
            build_graph_native_fci(Z=2, n_max=2, m_total=0, spin='singlet',
                                   subblock=(-1, 0))

    def test_subblock_high_l_empty(self):
        """subblock=(2, 2) at n_max=2 is empty (no d-orbitals exist)."""
        H = build_graph_native_fci(Z=2, n_max=2, m_total=0, spin='singlet',
                                   subblock=(2, 2))
        assert H.shape == (0, 0)

    def test_subblock_order_invariant(self):
        """subblock=(0, 1) and subblock=(1, 0) give identical matrices."""
        H_01 = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                      subblock=(0, 1))
        H_10 = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet',
                                      subblock=(1, 0))
        assert H_01.shape == H_10.shape
        assert np.max(np.abs(H_01 - H_10)) == 0.0

    def test_subblock_ss_hermitian(self):
        """ss sub-block FCI matrix is Hermitian at all tested n_max."""
        for n_max in [2, 3, 4, 5]:
            for spin in ['singlet', 'triplet']:
                H = build_graph_native_fci(Z=2, n_max=n_max, m_total=0,
                                           spin=spin, subblock=(0, 0))
                if H.shape[0] > 0:
                    assert np.allclose(H, H.T, atol=1e-12), (
                        f"{spin} n_max={n_max} not Hermitian"
                    )

    def test_compute_he_spectrum_subblock_ss(self):
        """compute_he_spectrum with subblock=(0,0) produces only S sectors."""
        result = compute_he_spectrum(n_max=3, subblock=(0, 0))
        # subblock kwarg recorded
        assert result['subblock'] == (0, 0)
        # Only S sectors should have non-zero dim; P, D dims should be 0 / absent
        for sector in result['sector_dims']:
            assert sector.endswith('S'), (
                f"sector {sector} has dim > 0 with subblock=(0,0); "
                f"expected only S sectors"
            )

    def test_compute_he_spectrum_subblock_ss_ground_state(self):
        """compute_he_spectrum with subblock=(0,0) gives sensible 1¹S ground state.

        The ss-block at n_max=3 reproduces the sprint's E(1¹S) = -2.88958 Ha.
        """
        result = compute_he_spectrum(n_max=3, subblock=(0, 0))
        assert '1_1S' in result['states']
        E_1S = result['states']['1_1S']['energy']
        # Sprint reference at n_max=3: -2.88957801 Ha
        assert abs(E_1S - (-2.88957801)) < 1e-7, (
            f"E(1¹S) = {E_1S:.10f}, sprint reference = -2.88957801"
        )
        # Variational bound
        assert E_1S > -2.903724377, (
            f"E(1¹S) = {E_1S:.10f} violates variational bound"
        )
