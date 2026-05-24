"""
ARCHIVED 2026-05-23 (Cleanup Track B): These three test classes
(TestSpectralRadialSolver, TestSpectralCoupledChannel, TestAlgebraicLaguerreMatrices)
were extracted from tests/test_hyperspherical_he.py in the redirect-before-archive
cleanup. They depend on solve_radial_spectral, solve_coupled_radial_spectral,
_build_laguerre_matrices_dirichlet, and the
solve_helium(radial_method='spectral', n_basis_radial, alpha_radial, matrix_method)
kwargs — all of which were removed from geovac/hyperspherical_radial.py in v2.7.0
(commit 8d692a0). The simplified hyperspherical_radial.py in current code is pure
FD. _build_laguerre_SK_algebraic still exists in geovac/level3_variational.py but
the rest of the test infrastructure here is dead, so the tests cannot be run.

The remaining live tests (TestGauntIntegrals, TestAngularSolver, TestAdiabaticCurve,
TestHeliumSolver, TestRadialSolver) remain in tests/test_hyperspherical_he.py.
"""

import numpy as np
import pytest

from geovac.hyperspherical_angular import solve_angular, gaunt_integral
from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
)
from geovac.hyperspherical_radial import (
    solve_radial, solve_helium,
    solve_radial_spectral, solve_coupled_radial_spectral,
    _build_laguerre_matrices_dirichlet, _build_laguerre_SK_algebraic,
)


E_EXACT = -2.903724


class TestSpectralRadialSolver:
    """Tests for the spectral Laguerre hyperradial solver (Lane 2, v2.0.9)."""

    @pytest.fixture(scope="class")
    def he_spectral(self) -> dict:
        """Spectral single-channel He solve using FD angular input."""
        return solve_helium(
            Z=2.0, l_max=0, n_alpha=200,
            N_R_angular=100, R_min=0.05, R_max=30.0,
            N_R_radial=3000, verbose=False,
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
        )

    @pytest.fixture(scope="class")
    def he_fd(self) -> dict:
        """FD single-channel He solve for comparison."""
        return solve_helium(
            Z=2.0, l_max=0, n_alpha=200,
            N_R_angular=100, R_min=0.05, R_max=30.0,
            N_R_radial=3000, verbose=False,
            radial_method='fd',
        )

    def test_spectral_fd_consistency(self, he_spectral: dict,
                                     he_fd: dict) -> None:
        """Spectral and FD agree within FD discretization error (~0.01 Ha)."""
        diff = abs(he_spectral['energy'] - he_fd['energy'])
        assert diff < 0.01, f"Spectral-FD diff {diff:.6f} Ha > 0.01 Ha"

    def test_spectral_energy_close_to_exact(self, he_spectral: dict) -> None:
        """Spectral single-channel within 0.1% of exact."""
        E = he_spectral['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.001, f"Error {err:.4f} > 0.1%"

    def test_spectral_dimension_reduction(self) -> None:
        """Spectral uses 25 basis functions vs 3000 FD grid points (120x)."""
        # n_basis_radial=25, N_R_radial=3000
        assert 3000 / 25 >= 100

    def test_spectral_convergence_with_n_basis(self) -> None:
        """Energy converges monotonically as n_basis increases."""
        energies = []
        for nb in [10, 15, 20, 25]:
            r = solve_helium(
                Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
                radial_method='spectral', n_basis_radial=nb, alpha_radial=1.5,
                verbose=False,
            )
            energies.append(r['energy'])
        # Should converge: differences shrink
        diffs = [abs(energies[i+1] - energies[i]) for i in range(len(energies)-1)]
        assert diffs[-1] < diffs[0], "Not converging with n_basis"

    def test_spectral_alpha_insensitivity(self) -> None:
        """Energy is insensitive to alpha choice over a reasonable range."""
        energies = []
        for alpha in [1.0, 1.5, 2.0, 2.5]:
            r = solve_helium(
                Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
                radial_method='spectral', n_basis_radial=25, alpha_radial=alpha,
                verbose=False,
            )
            energies.append(r['energy'])
        spread = max(energies) - min(energies)
        assert spread < 0.001, f"Alpha sensitivity: spread={spread:.6f} Ha"

    def test_spectral_wavefunction_decays(self, he_spectral: dict) -> None:
        """Wavefunction decays at large R."""
        F = he_spectral['wavefunction']
        assert np.max(np.abs(F[-100:])) / np.max(np.abs(F)) < 0.01

    def test_spectral_default_unchanged(self) -> None:
        """Default radial_method='fd' produces same result as before."""
        r = solve_helium(
            Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
            N_R_radial=2000, verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT)
        assert err < 0.001


class TestSpectralCoupledChannel:
    """Tests for spectral Laguerre in coupled-channel mode."""

    def test_spectral_coupled_fd_consistency(self) -> None:
        """Spectral and FD coupled-channel agree within 0.001 Ha."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r_fd = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='fd', verbose=False,
        )
        r_sp = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            verbose=False,
        )
        diff = abs(r_fd['energy'] - r_sp['energy'])
        assert diff < 0.001, f"Coupled-channel diff {diff:.6f} Ha > 0.001 Ha"

    def test_spectral_coupled_lmax3_accuracy(self) -> None:
        """Coupled-channel at l_max=3 gives ~0.22% error (matching FD ceiling)."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=3, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT)
        assert 0.15 < err * 100 < 0.30, \
            f"Expected 0.15-0.30% error, got {err*100:.4f}%"

    def test_spectral_coupled_speedup(self) -> None:
        """Spectral coupled-channel solver achieves >10x radial-only speedup."""
        import time
        from geovac.algebraic_angular import AlgebraicAngularSolver
        from geovac.algebraic_coupled_channel import compute_algebraic_coupling
        from scipy.interpolate import CubicSpline
        from geovac.hyperspherical_radial import (
            solve_coupled_radial,
            solve_coupled_radial_spectral,
        )

        solver = AlgebraicAngularSolver(2.0, 15, 0)
        R_grid = np.concatenate([
            np.linspace(0.1, 1.0, 67), np.linspace(1.0, 5.0, 67),
            np.linspace(5.0, 30.0, 67),
        ])
        R_grid = np.unique(R_grid)
        coupling = compute_algebraic_coupling(solver, R_grid, 3,
                                              compute_exact_dPdR=True)
        V_sp = []
        for ch in range(3):
            V_ch = effective_potential(R_grid, coupling['mu'][ch])
            V_sp.append(CubicSpline(R_grid, V_ch, extrapolate=True))
        P_sp = [[None]*3 for _ in range(3)]
        Q_sp = [[None]*3 for _ in range(3)]
        for mu in range(3):
            for nu in range(3):
                P_sp[mu][nu] = CubicSpline(R_grid, coupling['P'][mu,nu],
                                           extrapolate=True)
                Q_sp[mu][nu] = CubicSpline(R_grid, coupling['Q_exact'][mu,nu],
                                           extrapolate=True)

        t0 = time.time()
        solve_coupled_radial(V_sp, P_sp, Q_splines=Q_sp, n_channels=3,
                             R_min=0.1, R_max=30.0, N_R=3000, n_states=5,
                             sigma=-3.0, include_Q=True)
        t_fd = time.time() - t0

        t0 = time.time()
        solve_coupled_radial_spectral(V_sp, P_sp, Q_splines=Q_sp,
                                      n_channels=3, n_basis=25, alpha=1.5,
                                      R_min=0.1, n_states=5, include_Q=True)
        t_sp = time.time() - t0

        speedup = t_fd / t_sp
        assert speedup > 10, f"Speedup {speedup:.1f}x < 10x target"


class TestAlgebraicLaguerreMatrices:
    """Tests for algebraic Laguerre matrix elements (Track H, v2.0.10).

    Validates that the three-term Laguerre recurrence produces overlap S
    and kinetic K matrices matching Gauss-Laguerre quadrature to machine
    precision. Potential V remains quadrature (V_eff is transcendental).
    """

    def test_algebraic_vs_quadrature_overlap(self) -> None:
        """Algebraic overlap S matches quadrature to quadrature precision.

        The algebraic result is exact; differences are quadrature roundoff.
        M2 elements grow as ~6n^2, so absolute diffs scale with matrix size.
        Tolerance 1e-10 is ~1e-14 relative for the largest elements.
        """
        alpha = 1.5
        n_basis = 25
        V_dummy = lambda R: np.zeros_like(R)
        S_q, K_q, _, _, _, _, _, _ = _build_laguerre_matrices_dirichlet(
            V_dummy, n_basis, alpha, 0.05, matrix_method='quadrature'
        )
        S_a, K_a = _build_laguerre_SK_algebraic(n_basis, alpha)
        max_diff_S = np.max(np.abs(S_a - S_q))
        assert max_diff_S < 1e-10, f"Overlap max diff {max_diff_S:.2e} > 1e-10"

    def test_algebraic_vs_quadrature_kinetic(self) -> None:
        """Algebraic kinetic K matches quadrature to quadrature precision."""
        alpha = 1.5
        n_basis = 25
        V_dummy = lambda R: np.zeros_like(R)
        S_q, K_q, _, _, _, _, _, _ = _build_laguerre_matrices_dirichlet(
            V_dummy, n_basis, alpha, 0.05, matrix_method='quadrature'
        )
        S_a, K_a = _build_laguerre_SK_algebraic(n_basis, alpha)
        max_diff_K = np.max(np.abs(K_a - K_q))
        assert max_diff_K < 1e-10, f"Kinetic max diff {max_diff_K:.2e} > 1e-10"

    def test_algebraic_overlap_symmetry(self) -> None:
        """Algebraic overlap is symmetric."""
        S, _ = _build_laguerre_SK_algebraic(25, 1.5)
        assert np.allclose(S, S.T, atol=1e-15)

    def test_algebraic_kinetic_symmetry(self) -> None:
        """Algebraic kinetic is symmetric."""
        _, K = _build_laguerre_SK_algebraic(25, 1.5)
        assert np.allclose(K, K.T, atol=1e-15)

    def test_algebraic_overlap_positive_definite(self) -> None:
        """Algebraic overlap has all positive eigenvalues."""
        S, _ = _build_laguerre_SK_algebraic(25, 1.5)
        evals = np.linalg.eigvalsh(S)
        assert np.all(evals > 0), f"Non-positive overlap eigenvalue: {evals.min()}"

    def test_algebraic_single_channel_energy(self) -> None:
        """Algebraic single-channel energy matches quadrature to < 1e-14 Ha."""
        r_q = solve_helium(
            Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            matrix_method='quadrature', verbose=False,
        )
        r_a = solve_helium(
            Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            matrix_method='algebraic', verbose=False,
        )
        diff = abs(r_q['energy'] - r_a['energy'])
        assert diff < 1e-10, f"Energy diff {diff:.2e} Ha > 1e-10 Ha"

    def test_algebraic_coupled_channel_energy(self) -> None:
        """Algebraic coupled-channel energy matches quadrature to < 1e-10 Ha."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r_q = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            q_mode='exact', radial_method='spectral',
            n_basis_radial=25, alpha_radial=1.5,
            matrix_method='quadrature', verbose=False,
        )
        r_a = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            q_mode='exact', radial_method='spectral',
            n_basis_radial=25, alpha_radial=1.5,
            matrix_method='algebraic', verbose=False,
        )
        diff = abs(r_q['energy'] - r_a['energy'])
        assert diff < 1e-10, f"Coupled energy diff {diff:.2e} Ha > 1e-10 Ha"

    def test_algebraic_ceiling_unchanged(self) -> None:
        """Algebraic coupled at l_max=3 gives same 0.15-0.30% ceiling."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=3, n_channels=3,
            q_mode='exact', radial_method='spectral',
            n_basis_radial=25, alpha_radial=1.5,
            matrix_method='algebraic', verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT) * 100
        assert 0.15 < err < 0.30, f"Expected 0.15-0.30% error, got {err:.4f}%"

    def test_algebraic_alpha_range(self) -> None:
        """Algebraic S and K match quadrature across alpha = 0.5 to 3.0."""
        V_dummy = lambda R: np.zeros_like(R)
        for alpha in [0.5, 1.0, 2.0, 3.0]:
            S_q, K_q, _, _, _, _, _, _ = _build_laguerre_matrices_dirichlet(
                V_dummy, 20, alpha, 0.05, matrix_method='quadrature'
            )
            S_a, K_a = _build_laguerre_SK_algebraic(20, alpha)
            assert np.max(np.abs(S_a - S_q)) < 1e-10, \
                f"Overlap mismatch at alpha={alpha}"
            assert np.max(np.abs(K_a - K_q)) < 1e-10, \
                f"Kinetic mismatch at alpha={alpha}"

    def test_algebraic_n_basis_range(self) -> None:
        """Algebraic S and K match quadrature across n_basis = 10 to 30."""
        V_dummy = lambda R: np.zeros_like(R)
        for nb in [10, 15, 20, 25, 30]:
            S_q, K_q, _, _, _, _, _, _ = _build_laguerre_matrices_dirichlet(
                V_dummy, nb, 1.5, 0.05, matrix_method='quadrature'
            )
            S_a, K_a = _build_laguerre_SK_algebraic(nb, 1.5)
            assert np.max(np.abs(S_a - S_q)) < 1e-10, \
                f"Overlap mismatch at n_basis={nb}"
            assert np.max(np.abs(K_a - K_q)) < 1e-10, \
                f"Kinetic mismatch at n_basis={nb}"

    def test_algebraic_default_backward_compatible(self) -> None:
        """Default matrix_method='quadrature' is backward compatible."""
        r = solve_helium(
            Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT)
        assert err < 0.001
