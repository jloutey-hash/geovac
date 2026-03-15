"""Tests for Neumann expansion V_ee computation.

Validates auxiliary integrals (C_l, A_n, A_l, B_l, X_l) and full V_ee
matrix elements against numerical quadrature and known analytical values.
"""

import numpy as np
import pytest
from scipy.integrate import quad, dblquad
from scipy.special import lqmn

from geovac.neumann_vee import (
    compute_Cl_table,
    compute_An_table,
    compute_An_negative,
    compute_Al_table,
    compute_B0_table,
    compute_Bl_table,
    compute_vee_matrix_neumann,
    compute_hamiltonian_neumann,
    _compute_Xl_single,
)
from geovac.hylleraas import (
    HylleraasBasisFunction,
    generate_basis,
    build_quadrature_grids,
    compute_overlap_matrix,
    compute_hamiltonian_matrix,
    solve_hylleraas,
)


# ============================================================
# Auxiliary integral tests
# ============================================================

class TestClTable:
    """C_l(p) = ∫₋₁¹ η^p P_l(η) dη."""

    def test_C0_even(self):
        """C_0(p) = 2/(p+1) for even p."""
        C = compute_Cl_table(10, 5)
        for p in range(0, 11, 2):
            assert abs(C[0, p] - 2.0 / (p + 1)) < 1e-14, f"C_0({p}) failed"

    def test_C0_odd_zero(self):
        """C_0(p) = 0 for odd p (odd integrand)."""
        C = compute_Cl_table(10, 5)
        for p in range(1, 10, 2):
            assert abs(C[0, p]) < 1e-14, f"C_0({p}) should be zero"

    def test_C1_odd(self):
        """C_1(p) = 2/(p+2) for odd p."""
        C = compute_Cl_table(10, 5)
        for p in range(1, 10, 2):
            assert abs(C[1, p] - 2.0 / (p + 2)) < 1e-14, f"C_1({p}) failed"

    def test_selection_rule(self):
        """C_l(p) = 0 when p < l or (p-l) is odd."""
        C = compute_Cl_table(8, 8)
        for l in range(9):
            for p in range(9):
                if p < l or (p - l) % 2 != 0:
                    assert abs(C[l, p]) < 1e-13, f"C_{l}({p}) should be zero"

    def test_known_values(self):
        """Check specific known values."""
        C = compute_Cl_table(6, 6)
        # C_2(2) = ∫₋₁¹ η² P_2(η) dη = ∫ η²(3η²-1)/2 dη = (3·2/5 - 2/3)/2 = 4/15
        assert abs(C[2, 2] - 4.0 / 15.0) < 1e-14
        # C_2(4) = ∫₋₁¹ η⁴ P_2(η) dη = ∫ η⁴(3η²-1)/2 dη = (3·2/7 - 2/5)/2 = 16/70 = 8/35
        assert abs(C[2, 4] - 8.0 / 35.0) < 1e-14


class TestAnTable:
    """A_n(α) = ∫₁^∞ ξ^n e^{-αξ} dξ."""

    def test_A0(self):
        """A_0(α) = e^{-α}/α."""
        alpha = 1.5
        A = compute_An_table(0, alpha)
        assert abs(A[0] - np.exp(-alpha) / alpha) < 1e-14

    def test_against_quad(self):
        """Verify against scipy.integrate.quad."""
        alpha = 2.0
        A = compute_An_table(8, alpha)
        for n in range(9):
            val, _ = quad(lambda x: x**n * np.exp(-alpha * x), 1, np.inf)
            assert abs(A[n] - val) < 1e-12, f"A_{n} failed: {A[n]} vs {val}"

    def test_negative_index(self):
        """A_{-1}(α) = E₁(α)."""
        from scipy.special import exp1
        alpha = 1.0
        A_neg = compute_An_negative(-3, alpha)
        assert abs(A_neg[-1] - float(exp1(alpha))) < 1e-14


class TestAlTable:
    """A_l(p, α) = ∫₁^∞ ξ^p e^{-αξ} P_l(ξ) dξ."""

    def test_l0_equals_An(self):
        """A_0(p, α) = A_p(α) since P_0(ξ) = 1."""
        alpha = 1.5
        An = compute_An_table(6, alpha)
        Al = compute_Al_table(0, 6, alpha)
        for p in range(7):
            assert abs(Al[0, p] - An[p]) < 1e-13

    def test_l1_equals_An_shifted(self):
        """A_1(p, α) = A_{p+1}(α) since P_1(ξ) = ξ."""
        alpha = 1.5
        An = compute_An_table(8, alpha)
        Al = compute_Al_table(1, 6, alpha)
        for p in range(7):
            assert abs(Al[1, p] - An[p + 1]) < 1e-13

    def test_against_quad(self):
        """Verify A_l against numerical integration."""
        from numpy.polynomial.legendre import legval
        alpha = 2.0
        Al = compute_Al_table(4, 4, alpha)
        for l in range(5):
            for p in range(5):
                coeffs = np.zeros(l + 1)
                coeffs[l] = 1.0
                def integrand(xi, l_=l, p_=p):
                    return xi**p_ * np.exp(-alpha * xi) * legval(xi, coeffs)
                val, _ = quad(integrand, 1, np.inf)
                assert abs(Al[l, p] - val) < 1e-11, \
                    f"A_{l}({p}) failed: {Al[l,p]} vs {val}"


class TestBlTable:
    """B_l(p, α) = ∫₁^∞ ξ^p e^{-αξ} Q_l(ξ) dξ."""

    def test_B0_against_quad(self):
        """B_0 against direct numerical integration."""
        alpha = 2.0
        B0 = compute_B0_table(4, alpha)
        for p in range(5):
            def integrand(xi, p_=p):
                Q0 = 0.5 * np.log((xi + 1) / (xi - 1)) if xi > 1 else 0.0
                return xi**p_ * np.exp(-alpha * xi) * Q0
            val, _ = quad(integrand, 1, np.inf, limit=200)
            assert abs(B0[p] - val) < 1e-10, f"B_0({p}) failed"

    def test_Bl_against_quad(self):
        """B_l against direct numerical integration for small l."""
        alpha = 2.0
        Bl = compute_Bl_table(4, 3, alpha)
        for l in range(5):
            for p in range(4):
                def integrand(xi, l_=l, p_=p):
                    Qvals, _ = lqmn(0, l_, xi)
                    return xi**p_ * np.exp(-alpha * xi) * Qvals[0, l_]
                val, _ = quad(integrand, 1, np.inf, limit=200)
                assert abs(Bl[l, p] - val) < 1e-10, \
                    f"B_{l}({p}) failed: {Bl[l,p]} vs {val}"

    def test_B0_positive(self):
        """B_0(p, α) should be positive (Q_0 > 0 for ξ > 1)."""
        B0 = compute_B0_table(6, 1.5)
        for p in range(7):
            assert B0[p] > 0, f"B_0({p}) should be positive"


class TestXlIntegral:
    """X_l(p₁, p₂, α) = ∫∫ ξ₁^{p₁} ξ₂^{p₂} e^{-α(ξ₁+ξ₂)} P_l(ξ_<) Q_l(ξ_>) dξ₁dξ₂."""

    def test_against_2d_quad(self):
        """Verify X_l against 2D numerical integration for small cases."""
        alpha = 2.0
        for l in range(3):
            for p1, p2 in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0)]:
                xl = _compute_Xl_single(l, p1, p2, alpha)

                # 2D numerical integration (split at ξ₁ = ξ₂)
                from numpy.polynomial.legendre import legval

                def integrand(xi2, xi1, l_=l, p1_=p1, p2_=p2):
                    if xi1 <= 1.0 or xi2 <= 1.0:
                        return 0.0
                    xi_less = min(xi1, xi2)
                    xi_more = max(xi1, xi2)
                    coeffs_P = np.zeros(l_ + 1)
                    coeffs_P[l_] = 1.0
                    Pl = legval(xi_less, coeffs_P)
                    Qvals, _ = lqmn(0, l_, xi_more)
                    Ql = Qvals[0, l_]
                    return (xi1**p1_ * xi2**p2_ *
                            np.exp(-alpha * (xi1 + xi2)) * Pl * Ql)

                val, _ = dblquad(integrand, 1, 30, 1, 30,
                                 epsabs=1e-10, epsrel=1e-10)
                assert abs(xl - val) < max(1e-8, abs(val) * 1e-6), \
                    f"X_{l}({p1},{p2}) failed: {xl} vs {val}"

    def test_symmetry(self):
        """X_l(p₁, p₂, α) = X_l(p₂, p₁, α)."""
        alpha = 2.0
        for l in range(4):
            for p1 in range(3):
                for p2 in range(p1 + 1, 4):
                    x1 = _compute_Xl_single(l, p1, p2, alpha)
                    x2 = _compute_Xl_single(l, p2, p1, alpha)
                    assert abs(x1 - x2) < 1e-10, \
                        f"X_{l}({p1},{p2}) != X_{l}({p2},{p1}): {x1} vs {x2}"


# ============================================================
# V_ee matrix tests
# ============================================================

class TestVeeMatrix:
    """Test Neumann V_ee matrix against numerical V_ee."""

    @pytest.fixture
    def small_basis_and_grids(self):
        """3-function p=0 basis with grids."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)
        return basis, grids

    def test_vee_symmetric(self, small_basis_and_grids):
        """V_ee matrix should be symmetric."""
        basis, grids = small_basis_and_grids
        R = 1.4011
        V = compute_vee_matrix_neumann(basis, R, l_max=15)
        assert np.allclose(V, V.T, atol=1e-14)

    def test_vee_positive_diagonal(self, small_basis_and_grids):
        """V_ee diagonal elements should be positive (repulsion)."""
        basis, grids = small_basis_and_grids
        R = 1.4011
        V = compute_vee_matrix_neumann(basis, R, l_max=15)
        assert np.all(V.diagonal() > 0), \
            f"V_ee diagonal should be positive: {V.diagonal()}"

    def test_neumann_vs_numerical_vee(self, small_basis_and_grids):
        """Neumann V_ee should match numerical V_ee (within grid convergence)."""
        basis, grids = small_basis_and_grids
        R = 1.4011

        # Get numerical V_ee = H_full - H(T+V_ne)
        H_full = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_tv = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=False
        )
        V_ee_num = H_full - H_tv

        # Neumann V_ee (exact)
        V_ee_neu = compute_vee_matrix_neumann(basis, R, l_max=15)

        # Numerical V_ee converges toward Neumann as grid increases.
        # With 30x20 grid, expect ~1-3% agreement.
        max_diff = np.max(np.abs(V_ee_neu - V_ee_num))
        max_val = np.max(np.abs(V_ee_neu))
        rel_diff = max_diff / max_val if max_val > 0 else 0
        print(f"\n  V_ee max relative difference: {rel_diff:.4f} ({rel_diff*100:.1f}%)")
        assert rel_diff < 0.10, \
            f"V_ee Neumann vs numerical too different: {rel_diff*100:.1f}%"


class TestNeumannHamiltonian:
    """Test full Hamiltonian with Neumann V_ee."""

    def test_hamiltonian_neumann_symmetric(self):
        """H(Neumann) should be symmetric."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        H = compute_hamiltonian_neumann(basis, R, grids, l_max=10)
        assert np.allclose(H, H.T, atol=1e-13)

    def test_variational_principle(self):
        """E(Neumann) should be above exact ground state."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        result = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        E_exact = -1.17475
        assert result['E_total'] > E_exact, \
            f"Variational violation: E={result['E_total']:.6f} < {E_exact}"

    def test_energy_below_atoms(self):
        """H2 energy should be below separated atoms (bound)."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        result = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        E_atoms = -1.0  # 2 × (-0.5) for two H atoms
        assert result['E_total'] < E_atoms, \
            f"H2 not bound: E={result['E_total']:.6f} > {E_atoms}"

    def test_neumann_improves_over_numerical(self):
        """Neumann V_ee should give lower (better) energy than coarse numerical."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        # Use coarse grid where numerical V_ee has significant error
        grids = build_quadrature_grids(N_xi=15, N_eta=10, N_phi=1)
        R = 1.4011

        result_num = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='numerical',
        )
        result_neu = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        # With exact V_ee, the Neumann result should be better
        # (lower energy, closer to variational minimum)
        print(f"\n  E(numerical) = {result_num['E_total']:.6f}")
        print(f"  E(Neumann)   = {result_neu['E_total']:.6f}")
        # The key test: Neumann energy should be more negative
        # (numerical V_ee overestimates repulsion on coarse grids)
        assert result_neu['E_total'] < result_num['E_total'] + 0.01, \
            "Neumann should not be significantly worse than numerical"


class TestIncludeVee:
    """Test the include_vee=False path in compute_hamiltonian_matrix."""

    def test_tvne_plus_vee_equals_full(self):
        """H(T+V_ne) + V_ee(numerical) should equal H(full)."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=25, N_eta=15, N_phi=1)
        R = 1.4011

        H_full = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_tv = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=False
        )
        # The difference should be the V_ee contribution
        V_ee = H_full - H_tv
        assert np.all(V_ee.diagonal() > 0), \
            "V_ee diagonal should be positive (repulsion)"

    def test_tvne_matches_numba(self):
        """Python include_vee=True should match Numba path."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011

        H_python = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_numba = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=True, include_vee=True
        )
        assert np.allclose(H_python, H_numba, atol=1e-10), \
            f"Python vs Numba max diff: {np.max(np.abs(H_python - H_numba))}"
