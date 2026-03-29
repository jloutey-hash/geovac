"""
Tests for the associated Laguerre algebraic radial solver (Track J Sub-agent 2).

Validates:
1. Algebraic H and S matrices agree with generalized Laguerre quadrature reference
2. Eigenvalues of the generalized eigenvalue problem agree between methods
3. H2+ pi-state energies agree with FD solver
4. Overlap, potential, and kinetic matrices individually validated
"""
import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from scipy.linalg import eigh
from scipy.special import roots_genlaguerre, eval_genlaguerre, gammaln
from geovac.prolate_spheroidal_lattice import (
    _build_laguerre_matrices_algebraic,
    _build_laguerre_matrices_algebraic_associated,
    _associated_laguerre_moment_matrices,
    _lowered_moment_matrix,
    _stieltjes_matrix,
    ProlateSpheroidalLattice,
)


def _build_quadrature_reference_associated(
    N: int, alpha: float, m: int, A: float, a_param: float, c2: float
) -> tuple:
    """Build H and S for m!=0 using associated Laguerre basis via quadrature.

    Uses standard Gauss-Laguerre quadrature (weight e^{-x}) with explicit
    x^{alpha_L} factors in the integrand. This avoids the 1/x accuracy issues
    that arise when using generalized Laguerre quadrature with singular integrands.

    Returns separate (K, V_reg, V_cent, S) for individual component validation.
    """
    from scipy.special import roots_laguerre, eval_laguerre
    alpha_L = abs(m)
    sigma = abs(m) / 2.0
    beta = 2.0 * alpha

    # Standard Gauss-Laguerre quadrature (weight e^{-x})
    n_quad = max(5 * N, 150)
    x_q, w_q = roots_laguerre(n_quad)

    # Evaluate associated Laguerre polynomials at quadrature points
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_genlaguerre(n, alpha_L, x_q)

    xi_q = 1.0 + x_q / beta
    pref = beta ** (-(alpha_L + 1))

    # --- Overlap: S = pref * int x^{alpha_L} e^{-x} L_i L_j dx ---
    x_aL = x_q ** alpha_L
    wL = w_q[np.newaxis, :] * L_vals * x_aL[np.newaxis, :]
    S = pref * (wL @ L_vals.T)

    # --- Regular potential ---
    q_reg = A + a_param * xi_q - c2 * xi_q ** 2
    wLq = w_q[np.newaxis, :] * L_vals * x_aL[np.newaxis, :] * q_reg[np.newaxis, :]
    V_reg = pref * (wLq @ L_vals.T)

    # --- Centrifugal potential ---
    # V_cent = -m^2 beta^{1-alpha_L} int x^{alpha_L-1}/(x+2beta) e^{-x} L_i L_j dx
    weight_cent = x_q ** (alpha_L - 1) / (x_q + 2.0 * beta)
    wLc = w_q[np.newaxis, :] * L_vals * weight_cent[np.newaxis, :]
    V_cent = -m ** 2 * beta ** (1 - alpha_L) * (wLc @ L_vals.T)

    # --- Kinetic (integration by parts) ---
    # g_n(x) = (sigma+n-x/2) L_n - (n+alpha_L) L_{n-1}
    g_vals = np.zeros((N, n_quad))
    for n in range(N):
        g_vals[n] = (sigma + n - x_q / 2) * L_vals[n]
        if n >= 1:
            g_vals[n] -= (n + alpha_L) * L_vals[n - 1]

    # K = -pref * int x^{alpha_L-1}(x+2beta) e^{-x} g_i g_j dx
    weight_K = x_q ** (alpha_L - 1) * (x_q + 2.0 * beta)
    w_gK = w_q[np.newaxis, :] * g_vals * weight_K[np.newaxis, :]
    K = -pref * (w_gK @ g_vals.T)

    return K, V_reg, V_cent, S


# ============================================================
# Test 1: Matrix element agreement with quadrature
# ============================================================

class TestMatrixElementAgreement:
    """Validate algebraic matrix elements against quadrature reference."""

    @pytest.mark.parametrize("m", [1, 2, 3])
    @pytest.mark.parametrize("alpha", [0.5, 1.0, 2.5])
    def test_overlap_agreement(self, m, alpha):
        """Overlap matrix: algebraic vs quadrature."""
        N = 15
        A, a_param, c2 = 1.5, 3.0, 2.0
        _, S_alg = _build_laguerre_matrices_algebraic_associated(
            N, alpha, A, a_param, c2, m
        )
        _, _, _, S_quad = _build_quadrature_reference_associated(
            N, alpha, m, A, a_param, c2
        )
        max_err = np.max(np.abs(S_alg - S_quad))
        assert max_err < 1e-10, f"|m|={m}, alpha={alpha}: overlap max err {max_err:.2e}"

    @pytest.mark.parametrize("m", [1, 2, 3])
    @pytest.mark.parametrize("alpha", [0.5, 1.0, 2.5])
    def test_full_H_agreement(self, m, alpha):
        """Full Hamiltonian: algebraic vs quadrature."""
        N = 15
        A, a_param, c2 = 1.5, 3.0, 2.0
        H_alg, _ = _build_laguerre_matrices_algebraic_associated(
            N, alpha, A, a_param, c2, m
        )
        K_q, V_reg_q, V_cent_q, _ = _build_quadrature_reference_associated(
            N, alpha, m, A, a_param, c2
        )
        H_quad = K_q + V_reg_q + V_cent_q
        max_err = np.max(np.abs(H_alg - H_quad))
        scale = max(np.max(np.abs(H_quad)), 1.0)
        assert max_err < 1e-8 * scale, (
            f"|m|={m}, alpha={alpha}: H max err {max_err:.2e} (scale {scale:.2e})"
        )

    @pytest.mark.parametrize("m", [1, 2])
    def test_kinetic_component(self, m):
        """Kinetic matrix alone: algebraic vs quadrature."""
        N = 12
        alpha = 1.5
        A, a_param, c2 = 0.0, 0.0, 0.0  # Zero potential isolates kinetic

        H_alg, S_alg = _build_laguerre_matrices_algebraic_associated(
            N, alpha, A, a_param, c2, m
        )
        K_q, V_reg_q, V_cent_q, S_q = _build_quadrature_reference_associated(
            N, alpha, m, A, a_param, c2
        )
        # With A=a_param=c2=0: V_reg = 0, so H_alg = K_alg + V_cent_alg
        # Compare K+V_cent since they can't be separated in the algebraic version
        H_quad = K_q + V_cent_q
        max_err = np.max(np.abs(H_alg - H_quad))
        scale = max(np.max(np.abs(H_quad)), 1.0)
        assert max_err < 1e-8 * scale, (
            f"|m|={m}: K+V_cent max err {max_err:.2e} (scale {scale:.2e})"
        )


# ============================================================
# Test 2: Symmetry and positive-definiteness
# ============================================================

class TestMatrixProperties:
    """Validate structural properties of the algebraic matrices."""

    @pytest.mark.parametrize("m", [1, 2, 3])
    def test_H_symmetric(self, m):
        """H matrix must be real symmetric."""
        N = 15
        H, S = _build_laguerre_matrices_algebraic_associated(
            N, 1.5, 2.0, 3.0, 2.0, m
        )
        assert np.allclose(H, H.T, atol=1e-14), "H not symmetric"

    @pytest.mark.parametrize("m", [1, 2, 3])
    def test_S_symmetric_positive_definite(self, m):
        """Overlap must be symmetric positive definite."""
        N = 15
        _, S = _build_laguerre_matrices_algebraic_associated(
            N, 1.5, 2.0, 3.0, 2.0, m
        )
        assert np.allclose(S, S.T, atol=1e-14), "S not symmetric"
        evals = np.linalg.eigvalsh(S)
        assert np.all(evals > 0), f"S not positive definite: min eval = {evals[0]:.2e}"


# ============================================================
# Test 3: GEP eigenvalue agreement
# ============================================================

class TestEigenvalueAgreement:
    """Eigenvalues of H c = lambda S c must agree between methods."""

    @pytest.mark.parametrize("m", [1, 2])
    @pytest.mark.parametrize("alpha", [1.0, 2.0])
    def test_eigenvalue_agreement(self, m, alpha):
        """All eigenvalues of the GEP must agree between algebraic and quadrature."""
        N = 15
        A, a_param, c2 = 2.0, 4.0, 3.0

        H_alg, S_alg = _build_laguerre_matrices_algebraic_associated(
            N, alpha, A, a_param, c2, m
        )
        K_q, V_reg_q, V_cent_q, S_q = _build_quadrature_reference_associated(
            N, alpha, m, A, a_param, c2
        )
        H_quad = K_q + V_reg_q + V_cent_q

        evals_alg = np.sort(eigh(H_alg, S_alg, eigvals_only=True))
        evals_quad = np.sort(eigh(H_quad, S_q, eigvals_only=True))

        # Top eigenvalues are most relevant (root-finding uses the largest)
        for i in range(min(5, N)):
            idx = N - 1 - i
            err = abs(evals_alg[idx] - evals_quad[idx])
            scale = max(abs(evals_quad[idx]), 1.0)
            assert err < 1e-8 * scale, (
                f"Eigenvalue {i} mismatch: alg={evals_alg[idx]:.10f}, "
                f"quad={evals_quad[idx]:.10f}, err={err:.2e}"
            )


# ============================================================
# Test 4: Delegation from _build_laguerre_matrices_algebraic
# ============================================================

class TestDelegation:
    """The m=0 algebraic function should delegate to associated version for m!=0."""

    def test_m0_unchanged(self):
        """m=0 path still works (no regression)."""
        N = 15
        H, S = _build_laguerre_matrices_algebraic(N, 1.5, 2.0, 3.0, 2.0, 0)
        assert H.shape == (N, N)
        assert S.shape == (N, N)

    @pytest.mark.parametrize("m", [1, -1, 2])
    def test_delegation_works(self, m):
        """m!=0 delegates to associated builder without error."""
        N = 15
        H, S = _build_laguerre_matrices_algebraic(N, 1.5, 2.0, 3.0, 2.0, m)
        assert H.shape == (N, N)
        assert S.shape == (N, N)


# ============================================================
# Test 5: H2+ pi-state validation against FD solver
# ============================================================

class TestH2PlusPiState:
    """Validate against FD solver for H2+ pi states (|m|=1)."""

    @pytest.mark.parametrize("R", [2.0, 4.0])
    def test_pi_state_spectral_vs_fd(self, R):
        """Spectral algebraic pi-state energy agrees with FD to < 0.01 Ha."""
        # FD reference
        fd = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1, N_xi=5000, xi_max=30.0,
            m=1, n_angular=0, n_radial=0, radial_method='fd'
        )
        E_fd = fd.total_energy()

        # Spectral algebraic
        sp = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1,
            m=1, n_angular=0, n_radial=0,
            radial_method='spectral', n_basis=25,
            matrix_method='algebraic'
        )
        E_sp = sp.total_energy()

        err = abs(E_sp - E_fd)
        assert err < 0.01, (
            f"R={R}: spectral algebraic E={E_sp:.6f}, FD E={E_fd:.6f}, "
            f"diff={err:.6f}"
        )

    @pytest.mark.parametrize("R", [2.0, 4.0])
    def test_pi_state_algebraic_vs_quadrature_spectral(self, R):
        """Algebraic (associated basis) and quadrature (ordinary basis) both converge to FD.

        Note: the two spectral methods use DIFFERENT bases (associated vs ordinary
        Laguerre), so they differ at finite N. Both should agree with FD to < 0.01 Ha.
        """
        fd = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1, N_xi=5000, xi_max=30.0,
            m=1, n_angular=0, n_radial=0, radial_method='fd'
        )
        E_fd = fd.total_energy()

        sp_alg = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1,
            m=1, n_angular=0, n_radial=0,
            radial_method='spectral', n_basis=25,
            matrix_method='algebraic'
        )
        sp_quad = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1,
            m=1, n_angular=0, n_radial=0,
            radial_method='spectral', n_basis=25,
            matrix_method='quadrature'
        )
        E_alg = sp_alg.total_energy()
        E_quad = sp_quad.total_energy()

        err_alg = abs(E_alg - E_fd)
        err_quad = abs(E_quad - E_fd)
        assert err_alg < 0.01, f"R={R}: algebraic vs FD err {err_alg:.6f}"
        assert err_quad < 0.01, f"R={R}: quadrature vs FD err {err_quad:.6f}"

    def test_m2_delta_state(self):
        """|m|=2 (delta state) spectral algebraic agrees with FD."""
        R = 2.0
        fd = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1, N_xi=5000, xi_max=30.0,
            m=2, n_angular=0, n_radial=0, radial_method='fd'
        )
        E_fd = fd.total_energy()

        sp = ProlateSpheroidalLattice(
            R=R, Z_A=1, Z_B=1,
            m=2, n_angular=0, n_radial=0,
            radial_method='spectral', n_basis=25,
            matrix_method='algebraic'
        )
        E_sp = sp.total_energy()

        err = abs(E_sp - E_fd)
        assert err < 0.01, (
            f"m=2 delta: spectral E={E_sp:.6f}, FD E={E_fd:.6f}, diff={err:.6f}"
        )


# ============================================================
# Test 6: Consistency between m and -m
# ============================================================

class TestMSignSymmetry:
    """Energy should be independent of sign of m."""

    @pytest.mark.parametrize("m", [1, 2])
    def test_m_sign_symmetry(self, m):
        R = 2.0
        sp_pos = ProlateSpheroidalLattice(
            R=R, m=m, radial_method='spectral', n_basis=20,
            matrix_method='algebraic'
        )
        sp_neg = ProlateSpheroidalLattice(
            R=R, m=-m, radial_method='spectral', n_basis=20,
            matrix_method='algebraic'
        )
        E_pos = sp_pos.total_energy()
        E_neg = sp_neg.total_energy()
        assert abs(E_pos - E_neg) < 1e-12, (
            f"m=+{m} and m=-{m} give different energies: {E_pos} vs {E_neg}"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
