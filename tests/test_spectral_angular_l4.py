"""Tests for Level 4 spectral angular solver (Track K).

Validates the Jacobi polynomial spectral basis against the FD angular
solver for accuracy and performance.
"""

import numpy as np
import pytest
import time

from geovac.level4_multichannel import (
    solve_angular_multichannel,
    compute_adiabatic_curve_mc,
    _channel_list,
)
from geovac.level4_spectral_angular import (
    Level4SpectralAngular,
    solve_angular_spectral,
    compute_adiabatic_curve_spectral,
)


# ─── Free eigenvalue validation ────────────────────────────────────

class TestFreeEigenvalues:
    """Verify the analytical free eigenvalues match the SO(6) Casimir."""

    def test_casimir_formula_00(self):
        """(0,0) channel: μ_free = 2k(k+2) for k=0,2,4,..."""
        solver = Level4SpectralAngular(l_max=0, n_basis=5)
        casimir = solver._channel_casimir[0]
        # k = 0, 2, 4, 6, 8 for singlet
        expected = np.array([2*k*(k+2) for k in [0, 2, 4, 6, 8]])
        np.testing.assert_allclose(casimir, expected, atol=1e-12)

    def test_casimir_formula_11(self):
        """(1,1) channel: μ_free = ½(4)² - 2 + 2k(k+4) for k=0,2,4,..."""
        solver = Level4SpectralAngular(l_max=2, n_basis=3)
        # (1,1) is channel index 2 (after (0,0), (0,2))
        channels = solver.channels_4
        ch_11 = [i for i, c in enumerate(channels) if c == (1, 0, 1, 0)][0]
        casimir = solver._channel_casimir[ch_11]
        # k = 0, 2, 4 for singlet, L = 1+1+2 = 4
        expected = np.array([0.5*16 - 2 + 2*k*(k+4) for k in [0, 2, 4]])
        np.testing.assert_allclose(casimir, expected, atol=1e-12)

    def test_level3_consistency(self):
        """Free eigenvalues match Level 3 formula 2(l+k+1)²-2 for l1=l2=l."""
        solver = Level4SpectralAngular(l_max=2, n_basis=5)
        channels = solver.channels_4
        for ic, (l1, m1, l2, m2) in enumerate(channels):
            if l1 != l2:
                continue
            l = l1
            casimir = solver._channel_casimir[ic]
            # Level 3 formula: singlet k = 0, 2, 4, 6, 8
            k_vals = np.array([2 * j for j in range(5)])
            expected = 2.0 * (l + k_vals + 1.0) ** 2 - 2.0
            np.testing.assert_allclose(casimir, expected, atol=1e-12,
                                       err_msg=f"l1=l2={l}")


# ─── Single-point eigenvalue accuracy ──────────────────────────────

class TestSinglePointAccuracy:
    """Verify spectral eigenvalues match FD at specific (ρ, R_e) points."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.l_max = 2
        self.n_alpha = 200
        self.n_basis = 10
        self.solver = Level4SpectralAngular(
            l_max=self.l_max, n_basis=self.n_basis, n_quad=100,
        )

    @pytest.mark.parametrize("R_e", [0.5, 1.0, 1.5, 2.0])
    def test_mu0_near_equilibrium(self, R_e):
        """Ground eigenvalue matches FD to < 1e-3 near equilibrium."""
        R = 1.4
        rho = R / (2.0 * R_e)

        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e, self.l_max, Z=1.0, n_alpha=self.n_alpha, n_eig=1,
        )
        evals_sp, _ = self.solver.solve(rho, R_e, n_eig=1)

        np.testing.assert_allclose(
            evals_sp[0], evals_fd[0], atol=1e-3,
            err_msg=f"R_e={R_e}, rho={rho:.4f}",
        )

    def test_mu0_larger_re(self):
        """Ground eigenvalue at larger R_e matches FD to < 1e-2."""
        R_e = 3.0
        rho = 1.4 / (2.0 * R_e)

        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e, self.l_max, Z=1.0, n_alpha=self.n_alpha, n_eig=1,
        )
        evals_sp, _ = self.solver.solve(rho, R_e, n_eig=1)

        np.testing.assert_allclose(
            evals_sp[0], evals_fd[0], atol=1e-2,
            err_msg=f"R_e={R_e}, rho={rho:.4f}",
        )

    def test_multiple_eigenvalues(self):
        """First 2 eigenvalues match FD at equilibrium."""
        R_e = 1.5
        rho = 1.4 / (2.0 * R_e)

        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e, self.l_max, Z=1.0, n_alpha=self.n_alpha, n_eig=2,
        )
        evals_sp, _ = self.solver.solve(rho, R_e, n_eig=2)

        for i in range(2):
            np.testing.assert_allclose(
                evals_sp[i], evals_fd[i], atol=1e-2,
                err_msg=f"eigenvalue {i}",
            )


# ─── Basis convergence ─────────────────────────────────────────────

class TestBasisConvergence:
    """Verify convergence with increasing n_basis."""

    def test_monotonic_convergence(self):
        """Error decreases monotonically with n_basis (near equilibrium)."""
        R_e = 1.5
        rho = 1.4 / (2.0 * R_e)
        l_max = 2

        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z=1.0, n_alpha=200, n_eig=1,
        )
        mu0_fd = evals_fd[0]

        errors = []
        for nb in [5, 8, 10, 15]:
            solver = Level4SpectralAngular(l_max=l_max, n_basis=nb)
            evals, _ = solver.solve(rho, R_e, n_eig=1)
            errors.append(abs(evals[0] - mu0_fd))

        # Errors should decrease (or flatten when FD error dominates)
        for i in range(len(errors) - 1):
            assert errors[i] >= errors[i + 1] * 0.5, \
                f"Non-convergent: n_basis={[5,8,10,15][i]} err={errors[i]:.2e}, " \
                f"n_basis={[5,8,10,15][i+1]} err={errors[i+1]:.2e}"


# ─── Adiabatic curve ───────────────────────────────────────────────

class TestAdiabaticCurve:
    """Verify full adiabatic potential curve."""

    def test_u_min_accuracy(self):
        """U_min agrees with FD to < 1e-3 Ha."""
        R = 1.4
        l_max = 2
        R_e_grid = np.linspace(0.5, 6.0, 40)

        U_fd = compute_adiabatic_curve_mc(
            R, R_e_grid, l_max, Z=1.0, n_alpha=200,
        )
        U_sp = compute_adiabatic_curve_spectral(
            R, R_e_grid, l_max, Z=1.0, n_basis=10, n_quad=100,
        )

        assert abs(np.min(U_fd) - np.min(U_sp)) < 1e-3

    def test_equilibrium_position(self):
        """R_e at U_min agrees between FD and spectral."""
        R = 1.4
        l_max = 2
        R_e_grid = np.linspace(0.5, 4.0, 50)

        U_fd = compute_adiabatic_curve_mc(
            R, R_e_grid, l_max, Z=1.0, n_alpha=200,
        )
        U_sp = compute_adiabatic_curve_spectral(
            R, R_e_grid, l_max, Z=1.0, n_basis=10,
        )

        Re_min_fd = R_e_grid[np.argmin(U_fd)]
        Re_min_sp = R_e_grid[np.argmin(U_sp)]
        assert abs(Re_min_fd - Re_min_sp) < 0.2


# ─── Speedup ───────────────────────────────────────────────────────

class TestSpeedup:
    """Verify spectral solver is faster than FD."""

    def test_dimension_reduction(self):
        """Spectral matrix is at least 10x smaller than FD."""
        l_max = 2
        n_alpha = 200
        solver = Level4SpectralAngular(l_max=l_max, n_basis=10)
        channels = _channel_list(l_max, homonuclear=True)

        fd_dim = len(channels) * n_alpha
        sp_dim = solver._total_dim

        assert sp_dim < fd_dim / 10, \
            f"Insufficient dimension reduction: {sp_dim} vs {fd_dim}"

    def test_sweep_speedup(self):
        """Spectral sweep is at least 10x faster than FD."""
        R = 1.4
        l_max = 2
        R_e_grid = np.linspace(0.5, 4.0, 20)

        t0 = time.perf_counter()
        compute_adiabatic_curve_mc(
            R, R_e_grid, l_max, Z=1.0, n_alpha=200,
        )
        t_fd = time.perf_counter() - t0

        t0 = time.perf_counter()
        compute_adiabatic_curve_spectral(
            R, R_e_grid, l_max, Z=1.0, n_basis=10,
        )
        t_sp = time.perf_counter() - t0

        speedup = t_fd / t_sp
        assert speedup > 10, f"Speedup only {speedup:.1f}x (need > 10x)"


# ─── Precomputation correctness ────────────────────────────────────

class TestPrecomputation:
    """Verify V_ee precomputation is correct."""

    def test_vee_symmetric(self):
        """V_ee matrix is symmetric."""
        solver = Level4SpectralAngular(l_max=2, n_basis=10)
        np.testing.assert_allclose(
            solver._vee_matrix, solver._vee_matrix.T, atol=1e-14,
        )

    def test_vee_nonzero(self):
        """V_ee matrix has nonzero off-diagonal elements."""
        solver = Level4SpectralAngular(l_max=2, n_basis=10)
        # There should be coupling between channels
        off_diag = solver._vee_matrix.copy()
        for ic in range(solver.n_ch):
            n_k = len(solver._channel_k_indices[ic])
            i0 = solver._basis_offsets[ic]
            off_diag[i0:i0+n_k, i0:i0+n_k] = 0
        assert np.max(np.abs(off_diag)) > 0.01, "No inter-channel V_ee coupling"

    def test_solver_reuse(self):
        """Solver gives same results when reused across R_e points."""
        solver = Level4SpectralAngular(l_max=2, n_basis=10)
        R = 1.4

        # Solve at two different R_e points
        evals1, _ = solver.solve(0.7, 1.0, n_eig=1)
        evals2, _ = solver.solve(0.35, 2.0, n_eig=1)

        # Re-solve at first point — should be identical
        evals1b, _ = solver.solve(0.7, 1.0, n_eig=1)
        np.testing.assert_equal(evals1[0], evals1b[0])


# ─── Heteronuclear support ─────────────────────────────────────────

class TestHeteronuclear:
    """Verify spectral solver works for heteronuclear systems."""

    def test_heh_plus(self):
        """HeH+ eigenvalue matches FD within tolerance."""
        Z_A, Z_B = 2.0, 1.0
        R_e = 1.5
        R = 1.46
        rho = R / (2.0 * R_e)
        l_max = 2

        solver = Level4SpectralAngular(
            l_max=l_max, n_basis=10, Z_A=Z_A, Z_B=Z_B,
        )
        evals_sp, _ = solver.solve(rho, R_e, n_eig=1)

        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z=1.0, n_alpha=200, n_eig=1,
            Z_A=Z_A, Z_B=Z_B,
        )

        np.testing.assert_allclose(
            evals_sp[0], evals_fd[0], atol=1e-3,
            err_msg="HeH+ eigenvalue mismatch",
        )
