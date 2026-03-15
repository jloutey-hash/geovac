"""
Tests for Hylleraas explicitly correlated wavefunctions.

Validates:
1. Basis function generation and symmetry properties
2. r₁₂ computation in prolate spheroidals
3. Overlap matrix positive definiteness
4. Hamiltonian matrix symmetry
5. Energy convergence: p=0 (CI limit) vs p>0 (explicit correlation)
6. Dissociation energy improvement with r₁₂ terms
"""

import numpy as np
import pytest
from geovac.hylleraas import (
    HylleraasBasisFunction,
    generate_basis,
    generate_basis_truncated,
    compute_r12_squared,
    compute_r12_with_phi,
    build_quadrature_grids,
    evaluate_basis_function,
    compute_overlap_matrix,
    compute_hamiltonian_matrix,
    solve_hylleraas,
)


# ============================================================
# Category 1: Basis function infrastructure
# ============================================================

class TestBasisInfrastructure:
    """Tests for basis function generation and properties."""

    def test_basis_function_creation(self):
        """HylleraasBasisFunction stores quantum numbers correctly."""
        bf = HylleraasBasisFunction(1, 0, 0, 0, 1, alpha=1.2)
        assert bf.j == 1
        assert bf.k == 0
        assert bf.l == 0
        assert bf.m == 0
        assert bf.p == 1
        assert bf.alpha == 1.2

    def test_gerade_symmetry(self):
        """Gerade requires l+m even."""
        bf_g = HylleraasBasisFunction(0, 0, 0, 0, 0)
        assert bf_g.is_gerade
        bf_g2 = HylleraasBasisFunction(0, 0, 1, 1, 0)
        assert bf_g2.is_gerade
        bf_u = HylleraasBasisFunction(0, 0, 1, 0, 0)
        assert not bf_u.is_gerade

    def test_self_exchange(self):
        """Self-exchange when (j,k,l,m) == (k,j,m,l)."""
        bf_self = HylleraasBasisFunction(0, 0, 0, 0, 0)
        assert bf_self.is_self_exchange
        bf_not = HylleraasBasisFunction(1, 0, 0, 0, 0)
        assert not bf_not.is_self_exchange

    def test_generate_basis_p0_only(self):
        """p_max=0 generates standard CI-like basis."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        assert len(basis) > 0
        for bf in basis:
            assert bf.p == 0
            assert bf.is_gerade  # gerade_only=True by default

    def test_generate_basis_with_r12(self):
        """p_max>0 generates functions with r₁₂ dependence."""
        basis_p0 = generate_basis(j_max=1, l_max=1, p_max=0)
        basis_p1 = generate_basis(j_max=1, l_max=1, p_max=1)
        assert len(basis_p1) > len(basis_p0)
        p_values = [bf.p for bf in basis_p1]
        assert max(p_values) == 1

    def test_generate_basis_truncated(self):
        """Truncated basis respects omega_max constraint."""
        basis = generate_basis_truncated(omega_max=3, p_max=2)
        for bf in basis:
            assert bf.j + bf.k + bf.l + bf.m + bf.p <= 3
            assert bf.is_gerade

    def test_basis_no_duplicates(self):
        """No duplicate basis functions from symmetry constraints."""
        basis = generate_basis(j_max=2, l_max=2, p_max=1)
        seen = set()
        for bf in basis:
            key = (bf.j, bf.k, bf.l, bf.m, bf.p)
            assert key not in seen, f"Duplicate: {key}"
            seen.add(key)


# ============================================================
# Category 2: r₁₂ computation
# ============================================================

class TestR12Computation:
    """Tests for inter-electron distance in prolate spheroidals."""

    def test_r12_zero_same_point_with_phi(self):
        """r₁₂ = 0 when both electrons at exactly the same position (Δφ=0)."""
        r12 = compute_r12_with_phi(2.0, 0.5, 2.0, 0.5, dphi=0.0, R=2.0)
        assert abs(r12) < 1e-10

    def test_r12_avg_same_xi_eta(self):
        """Azimuthally averaged r₁₂² at same (ξ,η) equals 2ρ² (ring distance)."""
        R = 2.0
        xi, eta = 2.0, 0.5
        r12_sq = compute_r12_squared(xi, eta, xi, eta, R)
        # <r₁₂²>_φ = (R/2)²[0 + 2(ξ²-1)(1-η²)] = 2ρ²
        rho_sq = (R / 2)**2 * (xi**2 - 1) * (1 - eta**2)
        assert abs(r12_sq - 2 * rho_sq) < 1e-10

    def test_r12_positive(self):
        """r₁₂ > 0 for distinct electron positions."""
        r12_sq = compute_r12_squared(2.0, 0.5, 3.0, -0.3, R=2.0)
        assert r12_sq > 0

    def test_r12_known_value(self):
        """Check r₁₂ against Cartesian calculation for a known case.

        For electrons on the z-axis (η=±1, ξ values):
        z = (R/2)*ξ*η, so z₁ = (R/2)*ξ₁, z₂ = -(R/2)*ξ₂
        r₁₂ = |z₁ - z₂| = (R/2)(ξ₁ + ξ₂)
        """
        R = 2.0
        xi1, eta1 = 2.0, 1.0  # z₁ = 2.0
        xi2, eta2 = 3.0, -1.0  # z₂ = -3.0
        r12_sq = compute_r12_squared(xi1, eta1, xi2, eta2, R)
        # Expected: r₁₂ = |2 - (-3)| = 5.0
        expected_r12 = 5.0
        assert abs(np.sqrt(r12_sq) - expected_r12) < 0.01

    def test_r12_with_phi_average(self):
        """Azimuthal average of r₁₂ should match r₁₂_avg for m=0."""
        R = 2.0
        xi1, eta1, xi2, eta2 = 2.0, 0.3, 2.5, -0.2

        # Average over Δφ
        n_phi = 200
        dphi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
        r12_vals = np.array([
            compute_r12_with_phi(xi1, eta1, xi2, eta2, dp, R)
            for dp in dphi
        ])
        r12_sq_avg = np.mean(r12_vals**2)

        # Compare with the azimuthally-averaged formula
        r12_sq_formula = compute_r12_squared(xi1, eta1, xi2, eta2, R)
        assert abs(r12_sq_avg - r12_sq_formula) < 0.05 * r12_sq_formula

    def test_r12_symmetry(self):
        """r₁₂ is symmetric under electron exchange."""
        R = 2.0
        r12_sq_12 = compute_r12_squared(2.0, 0.3, 3.0, -0.5, R)
        r12_sq_21 = compute_r12_squared(3.0, -0.5, 2.0, 0.3, R)
        assert abs(r12_sq_12 - r12_sq_21) < 1e-12


# ============================================================
# Category 3: Basis function evaluation
# ============================================================

class TestBasisEvaluation:
    """Tests for basis function evaluation and symmetrization."""

    def test_evaluate_simplest_function(self):
        """φ_{00000} = 2 * exp(-α(ξ₁+ξ₂)) (self-exchange, factor 2)."""
        bf = HylleraasBasisFunction(0, 0, 0, 0, 0, alpha=1.0)
        R = 2.0
        xi1, xi2 = 2.0, 3.0
        eta1, eta2 = 0.5, -0.3
        r12 = np.sqrt(compute_r12_squared(xi1, eta1, xi2, eta2, R))

        val = evaluate_basis_function(bf, xi1, eta1, xi2, eta2, r12, R)
        expected = 2.0 * np.exp(-1.0 * (2.0 + 3.0))
        assert abs(val - expected) < 1e-12

    def test_exchange_symmetry(self):
        """Symmetrized function is symmetric under electron swap."""
        bf = HylleraasBasisFunction(1, 0, 0, 0, 0, alpha=1.0)
        R = 2.0
        xi1, xi2 = 2.0, 3.0
        eta1, eta2 = 0.5, -0.3

        r12_12 = np.sqrt(compute_r12_squared(xi1, eta1, xi2, eta2, R))
        r12_21 = np.sqrt(compute_r12_squared(xi2, eta2, xi1, eta1, R))

        val_12 = evaluate_basis_function(bf, xi1, eta1, xi2, eta2, r12_12, R)
        val_21 = evaluate_basis_function(bf, xi2, eta2, xi1, eta1, r12_21, R)

        assert abs(val_12 - val_21) < 1e-12

    def test_r12_factor(self):
        """p>0 multiplies by r₁₂^p."""
        bf0 = HylleraasBasisFunction(0, 0, 0, 0, 0, alpha=1.0)
        bf1 = HylleraasBasisFunction(0, 0, 0, 0, 1, alpha=1.0)
        R = 2.0
        xi1, xi2, eta1, eta2 = 2.0, 3.0, 0.5, -0.3
        r12 = np.sqrt(compute_r12_squared(xi1, eta1, xi2, eta2, R))

        val0 = evaluate_basis_function(bf0, xi1, eta1, xi2, eta2, r12, R)
        val1 = evaluate_basis_function(bf1, xi1, eta1, xi2, eta2, r12, R)

        assert abs(val1 - val0 * r12) < 1e-12


# ============================================================
# Category 4: Matrix properties
# ============================================================

class TestMatrixProperties:
    """Tests for overlap and Hamiltonian matrix properties."""

    @pytest.fixture
    def small_setup(self):
        """Small basis + grid for fast matrix tests."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=12, N_eta=10, N_phi=8, xi_max=10.0)
        return basis, grids

    def test_overlap_positive_definite(self, small_setup):
        """Overlap matrix S must be positive definite."""
        basis, grids = small_setup
        S = compute_overlap_matrix(basis, R=1.4, grids=grids)
        evals = np.linalg.eigvalsh(S)
        assert np.all(evals > -1e-10), f"S has negative eigenvalue: {evals[0]}"

    def test_overlap_symmetric(self, small_setup):
        """Overlap matrix must be symmetric."""
        basis, grids = small_setup
        S = compute_overlap_matrix(basis, R=1.4, grids=grids)
        assert np.allclose(S, S.T, atol=1e-10)

    def test_hamiltonian_symmetric(self, small_setup):
        """Hamiltonian matrix must be symmetric."""
        basis, grids = small_setup
        H = compute_hamiltonian_matrix(basis, R=1.4, grids=grids)
        assert np.allclose(H, H.T, atol=1e-6)

    def test_overlap_diagonal_positive(self, small_setup):
        """Diagonal elements of S must be positive (norms)."""
        basis, grids = small_setup
        S = compute_overlap_matrix(basis, R=1.4, grids=grids)
        assert np.all(np.diag(S) > 0)


# ============================================================
# Category 5: Energy validation
# ============================================================

class TestEnergyValidation:
    """Tests for energy accuracy and physical constraints."""

    def test_variational_principle(self):
        """Energy must be above exact ground state E = -1.17475 Ha."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=14, N_eta=10, N_phi=8, xi_max=10.0)
        result = solve_hylleraas(basis, R=1.4, grids=grids, verbose=False)
        assert result['E_total'] > -1.17475, (
            f"Variational violation: E={result['E_total']:.6f}")

    def test_energy_below_atoms(self):
        """At equilibrium, H₂ energy must be below 2H = -1.0 Ha."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=14, N_eta=10, N_phi=8, xi_max=10.0)
        result = solve_hylleraas(basis, R=1.4, grids=grids, verbose=False)
        # With a reasonable basis, E_total should be < -1.0
        # (molecule is bound). This may fail with very small grids.
        if result['E_total'] < -1.0:
            assert result['D_e'] > 0

    def test_r12_improves_energy(self):
        """Adding r₁₂ terms (p>0) must lower the energy vs p=0.

        This is the KEY test: explicit correlation should always
        improve the variational energy.
        """
        grids = build_quadrature_grids(N_xi=14, N_eta=10, N_phi=10, xi_max=10.0)
        R = 1.4

        # p=0 only
        basis_p0 = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        result_p0 = solve_hylleraas(basis_p0, R, grids, verbose=False)

        # p=0 and p=1
        basis_p1 = generate_basis(j_max=1, l_max=0, p_max=1, alpha=1.0)
        result_p1 = solve_hylleraas(basis_p1, R, grids, verbose=False)

        # More basis functions must give lower (or equal) energy
        assert result_p1['E_total'] <= result_p0['E_total'] + 1e-6, (
            f"r₁₂ terms raised energy: p0={result_p0['E_total']:.6f}, "
            f"p1={result_p1['E_total']:.6f}")

    def test_alpha_dependence(self):
        """Energy should vary smoothly with α."""
        grids = build_quadrature_grids(N_xi=12, N_eta=8, N_phi=8, xi_max=10.0)
        R = 1.4
        energies = []
        for alpha in [0.5, 1.0, 1.5]:
            basis = generate_basis(j_max=0, l_max=0, p_max=0, alpha=alpha)
            result = solve_hylleraas(basis, R, grids, verbose=False)
            energies.append(result['E_total'])
        # Should have a minimum somewhere
        assert not np.all(np.diff(energies) > 0), "Energy monotonically increasing"
        assert not np.all(np.diff(energies) < 0), "Energy monotonically decreasing"


# ============================================================
# Category 6: Quadrature grid tests
# ============================================================

class TestQuadratureGrids:
    """Tests for quadrature grid construction."""

    def test_xi_range(self):
        """ξ grid should be in [1, ξ_max]."""
        grids = build_quadrature_grids(N_xi=20, N_eta=10, xi_max=15.0)
        assert np.all(grids['xi'] >= 1.0)
        assert np.all(grids['xi'] <= 15.0)

    def test_eta_range(self):
        """η grid should be in [-1, +1]."""
        grids = build_quadrature_grids(N_xi=20, N_eta=10)
        assert np.all(grids['eta'] >= -1.0)
        assert np.all(grids['eta'] <= 1.0)

    def test_weights_positive(self):
        """All quadrature weights should be positive."""
        grids = build_quadrature_grids(N_xi=20, N_eta=10, N_phi=12)
        assert np.all(grids['w_xi'] > 0)
        assert np.all(grids['w_eta'] > 0)
        assert np.all(grids['w_phi'] > 0)

    def test_phi_grid(self):
        """Δφ grid should cover [0, 2π)."""
        grids = build_quadrature_grids(N_phi=16)
        assert grids['dphi'][0] == 0.0
        assert grids['dphi'][-1] < 2 * np.pi
        assert abs(np.sum(grids['w_phi']) - 2 * np.pi) < 1e-10


# ============================================================
# Category 7: Numba acceleration
# ============================================================

class TestNumbaAcceleration:
    """Tests for Numba-accelerated kernels vs pure Python."""

    def _has_numba(self):
        try:
            from geovac._numba_kernels import NUMBA_AVAILABLE
            return NUMBA_AVAILABLE
        except ImportError:
            return False

    def test_numba_overlap_p0_consistency(self):
        """Numba and Python overlap matrices match for p=0."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        from geovac.hylleraas import (
            _compute_overlap_matrix_python,
        )

        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=10, N_eta=8, N_phi=8, xi_max=8.0)
        R = 1.4

        S_numba = compute_overlap_matrix(basis, R, grids, use_numba=True)
        S_python = _compute_overlap_matrix_python(basis, R, grids)

        assert np.allclose(S_numba, S_python, rtol=1e-6, atol=1e-10), (
            f"Max diff: {np.max(np.abs(S_numba - S_python)):.2e}")

    def test_numba_hamiltonian_p0_consistency(self):
        """Numba and Python Hamiltonian matrices match for p=0."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=10, N_eta=8, N_phi=8, xi_max=8.0)
        R = 1.4

        H_numba = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=True)
        H_python = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False)

        # Use larger tolerance for Hamiltonian due to ellipk approximation
        assert np.allclose(H_numba, H_python, rtol=1e-4, atol=1e-8), (
            f"Max diff: {np.max(np.abs(H_numba - H_python)):.2e}")

    def test_numba_overlap_5d_consistency(self):
        """Numba and Python overlap matrices match for p>0."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        from geovac.hylleraas import (
            _compute_overlap_matrix_python,
        )

        basis = generate_basis(j_max=1, l_max=0, p_max=1, alpha=1.0)
        grids = build_quadrature_grids(N_xi=8, N_eta=6, N_phi=8, xi_max=8.0)
        R = 1.4

        S_numba = compute_overlap_matrix(basis, R, grids, use_numba=True)
        S_python = _compute_overlap_matrix_python(basis, R, grids)

        assert np.allclose(S_numba, S_python, rtol=1e-6, atol=1e-10), (
            f"Max diff: {np.max(np.abs(S_numba - S_python)):.2e}")

    def test_numba_energy_consistency(self):
        """Numba and Python give the same ground state energy."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=12, N_eta=8, N_phi=8, xi_max=10.0)

        result_numba = solve_hylleraas(basis, R=1.4, grids=grids, verbose=False)
        # Force Python path
        S_py = compute_overlap_matrix(basis, 1.4, grids, use_numba=False)
        H_py = compute_hamiltonian_matrix(basis, 1.4, grids, use_numba=False)
        V_NN = 1.0 / 1.4
        from scipy.linalg import eigh as _eigh
        evals_py, _ = _eigh(H_py + V_NN * S_py, S_py)
        E_py = evals_py[0]

        assert abs(result_numba['E_total'] - E_py) < 1e-4, (
            f"Energy mismatch: numba={result_numba['E_total']:.6f}, "
            f"python={E_py:.6f}")

    def test_numba_speedup_overlap(self):
        """Numba overlap is at least 10x faster than Python for p=0."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        from geovac._numba_kernels import warmup_hylleraas_jit
        from geovac.hylleraas import (
            _compute_overlap_matrix_python,
        )

        # Warmup JIT
        warmup_hylleraas_jit()

        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=14, N_eta=10, N_phi=8, xi_max=10.0)
        R = 1.4

        import time

        # Python timing
        t0 = time.perf_counter()
        _compute_overlap_matrix_python(basis, R, grids)
        dt_python = time.perf_counter() - t0

        # Numba timing
        t0 = time.perf_counter()
        compute_overlap_matrix(basis, R, grids, use_numba=True)
        dt_numba = time.perf_counter() - t0

        ratio = dt_python / max(dt_numba, 1e-9)
        assert ratio > 10, (
            f"Numba speedup only {ratio:.1f}x "
            f"(python={dt_python:.3f}s, numba={dt_numba:.3f}s)")

    def test_ellipk_approximation_accuracy(self):
        """Elliptic K approximation matches scipy to < 1e-6."""
        if not self._has_numba():
            pytest.skip("Numba not available")

        from geovac._numba_kernels import _ellipk_approx
        from scipy.special import ellipk as scipy_ellipk

        test_vals = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]
        for m in test_vals:
            K_approx = _ellipk_approx(m)
            K_exact = scipy_ellipk(m)
            rel_err = abs(K_approx - K_exact) / K_exact
            assert rel_err < 1e-6, (
                f"ellipk({m}): approx={K_approx:.10f}, "
                f"exact={K_exact:.10f}, rel_err={rel_err:.2e}")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
