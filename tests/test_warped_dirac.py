"""Tests for the warped Dirac module (G4-4a first move).

Verifies:
- 2D Cl(2, 0) gamma matrix algebra (F2 algebraic backbone)
- Discrete disk-Dirac construction (Hermitian, positive)
- S^2 Dirac spectrum (Camporesi-Higuchi)
- Constant-warp factorization at the heat-trace level (F1)
- Rank-2 spinor enhancement in continuum limit (F3 rough)
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.gravity.warped_dirac import (
    GAMMA_2D_1,
    GAMMA_2D_2,
    GAMMA_2D_5,
    I_2,
    DiscreteDirac2D,
    DiscreteDiskDirac,
    DiscreteDiskScalar,
    S2DiracSpectrum,
    VariableWarpDirac,
    WarpedDiracConstant,
    verify_F1_factorization,
    verify_F2_chirality,
    verify_F3_continuum_recovery_rough,
    verify_F4_tip_regular,
    verify_F6_riemannian_limit,
    verify_F7_factorization_loss,
    verify_gamma_algebra_2d,
)


# =================================================================
# Gamma matrix algebra (F2 backbone)
# =================================================================


class TestGammaAlgebra:
    def test_gamma_1_squared_is_identity(self):
        assert np.allclose(GAMMA_2D_1 @ GAMMA_2D_1, I_2)

    def test_gamma_2_squared_is_identity(self):
        assert np.allclose(GAMMA_2D_2 @ GAMMA_2D_2, I_2)

    def test_gamma_5_squared_is_identity(self):
        assert np.allclose(GAMMA_2D_5 @ GAMMA_2D_5, I_2)

    def test_gamma_1_gamma_2_anticommute(self):
        assert np.allclose(
            GAMMA_2D_1 @ GAMMA_2D_2 + GAMMA_2D_2 @ GAMMA_2D_1,
            np.zeros((2, 2)),
        )

    def test_gamma_5_is_minus_i_gamma_1_gamma_2(self):
        expected = -1j * (GAMMA_2D_1 @ GAMMA_2D_2)
        assert np.allclose(GAMMA_2D_5, expected)

    def test_gamma_5_anticommutes_with_gamma_1(self):
        anticom = GAMMA_2D_5 @ GAMMA_2D_1 + GAMMA_2D_1 @ GAMMA_2D_5
        assert np.allclose(anticom, np.zeros((2, 2)))

    def test_gamma_5_anticommutes_with_gamma_2(self):
        anticom = GAMMA_2D_5 @ GAMMA_2D_2 + GAMMA_2D_2 @ GAMMA_2D_5
        assert np.allclose(anticom, np.zeros((2, 2)))

    def test_verify_gamma_algebra_machine_precision(self):
        results = verify_gamma_algebra_2d()
        for name, residual in results.items():
            assert residual < 1e-14, f"{name} residual = {residual}"


# =================================================================
# Discrete disk-Dirac
# =================================================================


class TestDiscreteDiskDirac:
    def test_construction(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        assert disk.N_rho == 20
        assert disk.a == 0.5
        assert disk.N_phi == 12
        assert disk.R == 10.0
        assert np.isclose(disk.h_phi, 2 * np.pi / 12)

    def test_hilbert_dim(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        # 2 (spin) * 20 * 12 = 480
        assert disk.hilbert_dim == 480

    def test_invalid_N_rho(self):
        with pytest.raises(ValueError):
            DiscreteDiskDirac(N_rho=0, a=0.5, N_phi=12)

    def test_invalid_a(self):
        with pytest.raises(ValueError):
            DiscreteDiskDirac(N_rho=20, a=0.0, N_phi=12)

    def test_invalid_N_phi(self):
        with pytest.raises(ValueError):
            DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=1)

    def test_radial_laplacian_hermitian(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        H = disk._hermitian_radial_laplacian(m_eff=1.5)
        assert np.allclose(H, H.T)  # symmetric (real)

    def test_radial_laplacian_positive_eigenvalues(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        H = disk._hermitian_radial_laplacian(m_eff=1.5)
        evals = np.linalg.eigvalsh(H)
        assert np.all(evals > -1e-10)  # numerically nonnegative

    def test_squared_eigenvalues_count(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        # Each m_eff gives N_rho eigenvalues, doubled for rank-2 spinor
        # Total = 2 * N_rho * N_phi
        eigs = disk.squared_eigenvalues()
        assert len(eigs) == 2 * 20 * 12

    def test_squared_eigenvalues_nonnegative(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        eigs = disk.squared_eigenvalues()
        assert np.all(eigs > -1e-10)

    def test_squared_eigenvalues_doubled(self):
        """Each scalar eigenvalue appears exactly twice in spinor spectrum."""
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        eigs = disk.squared_eigenvalues()
        # Sort and check consecutive pairs are equal
        sorted_eigs = np.sort(eigs)
        # Group into pairs (since each scalar eigenvalue is doubled)
        for i in range(0, len(sorted_eigs), 2):
            assert np.isclose(sorted_eigs[i], sorted_eigs[i + 1]), (
                f"Pair at i={i}: {sorted_eigs[i]} vs {sorted_eigs[i+1]}"
            )

    def test_heat_trace_positive_decreasing(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.5, N_phi=12)
        K_small_t = disk.heat_trace(0.01)
        K_med_t = disk.heat_trace(0.1)
        K_large_t = disk.heat_trace(1.0)
        assert K_small_t > 0
        assert K_med_t > 0
        assert K_large_t > 0
        assert K_small_t > K_med_t > K_large_t


# =================================================================
# S^2 Dirac spectrum
# =================================================================


class TestS2DiracSpectrum:
    def test_construction(self):
        sphere = S2DiracSpectrum(l_max=4, r_h=2.0)
        assert sphere.l_max == 4
        assert sphere.r_h == 2.0

    def test_invalid_l_max(self):
        with pytest.raises(ValueError):
            S2DiracSpectrum(l_max=-1)

    def test_invalid_r_h(self):
        with pytest.raises(ValueError):
            S2DiracSpectrum(l_max=4, r_h=0.0)

    def test_n_modes(self):
        # sum_{n=0}^{l_max} 8(n+1)
        sphere = S2DiracSpectrum(l_max=2, r_h=1.0)
        # n=0: 8, n=1: 16, n=2: 24. Total = 48
        assert sphere.n_modes() == 48

    def test_eigenvalues_at_unit_r_h(self):
        """At r_h = 1, |lambda_n^2| = (n+1)^2."""
        sphere = S2DiracSpectrum(l_max=2, r_h=1.0)
        eigs = sphere.squared_eigenvalues()
        assert len(eigs) == 48
        # Smallest eigenvalue at n=0: 1
        assert np.isclose(np.min(eigs), 1.0)
        # Largest at n=2: 9
        assert np.isclose(np.max(eigs), 9.0)

    def test_eigenvalues_scaled_by_r_h(self):
        """At r_h, |lambda_n^2| = ((n+1)/r_h)^2."""
        sphere = S2DiracSpectrum(l_max=1, r_h=2.0)
        eigs = sphere.squared_eigenvalues()
        # n=0: (1/2)^2 = 1/4
        # n=1: (2/2)^2 = 1.0
        assert np.isclose(np.min(eigs), 0.25)
        assert np.isclose(np.max(eigs), 1.0)

    def test_multiplicities(self):
        """g_n^Dirac = 4(n+1), each sign +/-, so 8(n+1) total squared."""
        sphere = S2DiracSpectrum(l_max=2, r_h=1.0)
        eigs = sphere.squared_eigenvalues()
        # Count each unique eigenvalue
        for n in range(3):
            lam_sq = (n + 1) ** 2
            count = int(np.sum(np.isclose(eigs, lam_sq)))
            assert count == 8 * (n + 1), (
                f"n={n}: expected {8*(n+1)}, got {count}"
            )


# =================================================================
# F1: factorization at constant warp
# =================================================================


class TestF1Factorization:
    def test_F1_small_parameters(self):
        """F1 at small panel size; direct = factorized to ~ float64 sum."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        t_values = [0.05, 0.1, 0.5, 1.0]
        results = verify_F1_factorization(disk, sphere, t_values, tol=1e-10)
        assert results["all_passed"], results

    def test_F1_includes_all_t(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        t_values = [0.05, 0.5]
        results = verify_F1_factorization(disk, sphere, t_values)
        # Per-t entries present
        assert "0.05" in results
        assert "0.5" in results

    def test_F1_bit_exact_at_machine_precision(self):
        """Tightest tolerance: factorization holds to ~ float64."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=6)
        sphere = S2DiracSpectrum(l_max=2, r_h=1.5)
        t_values = [0.1, 0.5]
        results = verify_F1_factorization(disk, sphere, t_values, tol=1e-12)
        # At small panel size, summation residual should be ~ 1e-13
        for t in t_values:
            assert results[str(t)]["rel_err"] < 1e-12, results[str(t)]


# =================================================================
# F2: chirality grading
# =================================================================


class TestF2Chirality:
    def test_F2_passes(self):
        results = verify_F2_chirality()
        assert results["passed"]

    def test_F2_algebra_machine_precision(self):
        results = verify_F2_chirality()
        for name, residual in results["algebra"].items():
            assert residual < 1e-14, f"{name}: {residual}"


# =================================================================
# F3 rough: rank-2 spinor enhancement
# =================================================================


class TestF3RoughRecovery:
    def test_F3_rough_panel(self):
        """K_Dirac / K_scalar -> 2 at small t (rank-2 spinor enhancement)."""
        disk = DiscreteDiskDirac(N_rho=50, a=0.2, N_phi=24)
        # At sprint-scale N_phi = 24, modest tolerance for the rank-2 ratio.
        t_values = [0.5, 1.0]
        results = verify_F3_continuum_recovery_rough(
            disk, t_values, tol_rank_2=0.25
        )
        # All t values should have positive K and reasonable ratio
        for t in t_values:
            entry = results[str(t)]
            assert entry["K_Dirac"] > 0
            assert entry["K_scalar"] > 0
            # Ratio should be in range [1.5, 2.5]
            assert 1.5 < entry["ratio"] < 2.5, (
                f"t={t}: ratio={entry['ratio']}"
            )


# =================================================================
# WarpedDiracConstant
# =================================================================


class TestWarpedDiracConstant:
    def test_construction(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        warp = WarpedDiracConstant(disk=disk, sphere=sphere)
        assert warp.disk is disk
        assert warp.sphere is sphere

    def test_hilbert_dim(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        warp = WarpedDiracConstant(disk=disk, sphere=sphere)
        # disk: 2 * 10 * 8 = 160
        # sphere: 8 + 16 + 24 = 48
        # cigar: 160 * 48 = 7680
        assert warp.hilbert_dim == 7680

    def test_eigenvalues_count(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        warp = WarpedDiracConstant(disk=disk, sphere=sphere)
        eigs = warp.squared_eigenvalues()
        assert len(eigs) == warp.hilbert_dim

    def test_heat_trace_factorized_matches_direct(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        warp = WarpedDiracConstant(disk=disk, sphere=sphere)
        for t in [0.1, 0.5, 1.0]:
            K_fact = warp.heat_trace_factorized(t)
            K_dir = warp.heat_trace_direct(t)
            assert np.isclose(K_fact, K_dir, rtol=1e-10), (
                f"t={t}: K_fact={K_fact}, K_dir={K_dir}"
            )

    def test_heat_trace_decreasing_in_t(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        warp = WarpedDiracConstant(disk=disk, sphere=sphere)
        K_01 = warp.heat_trace_factorized(0.1)
        K_05 = warp.heat_trace_factorized(0.5)
        K_10 = warp.heat_trace_factorized(1.0)
        assert K_01 > K_05 > K_10 > 0


# =================================================================
# Explicit DiscreteDirac2D (G4-4a week 2)
# =================================================================


class TestDiscreteDirac2DConstruction:
    def test_construction(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        assert d2d.N_rho == 10
        assert d2d.a == 0.5
        assert d2d.N_phi == 8
        assert d2d.R == 5.0
        assert d2d.hilbert_dim == 160  # 2 * 10 * 8

    def test_fourier_modes_count(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        modes = d2d.fourier_modes()
        assert len(modes) == 8

    def test_fourier_modes_nonzero(self):
        """Anti-periodic phi: even smallest mode has m_eff > 0."""
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=12)
        modes = d2d.fourier_modes()
        assert all(m > 0 for m in modes), modes

    def test_L_k_hermitian(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        L = d2d.L_k(0)
        assert np.allclose(L, L.T)

    def test_sqrt_L_k_hermitian(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        sL = d2d.sqrt_L_k(0)
        assert np.allclose(sL, sL.T.conj())

    def test_sqrt_L_k_squared_recovers_L_k(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        for k_idx in range(8):
            L = d2d.L_k(k_idx)
            sL = d2d.sqrt_L_k(k_idx)
            assert np.allclose(sL @ sL, L, atol=1e-10), (
                f"k_idx={k_idx}: sqrt(L)^2 != L"
            )

    def test_dirac_block_shape(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        D_k = d2d.dirac_block(0)
        assert D_k.shape == (20, 20)

    def test_dirac_block_hermitian(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        for k_idx in range(8):
            D_k = d2d.dirac_block(k_idx)
            assert np.allclose(D_k, D_k.T.conj()), f"k_idx={k_idx}"

    def test_chirality_block_squared_is_identity(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        g5 = d2d.chirality_block()
        assert np.allclose(g5 @ g5, np.eye(20))


class TestDiscreteDirac2DOperatorLevel:
    def test_F2_operator_level_anticommutes_per_block(self):
        """Operator-level F2: {gamma^5, D_k} = 0 to machine precision."""
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=12)
        results = d2d.verify_chirality_anticommute_per_block(tol=1e-12)
        assert results["passed"], results
        assert results["max_residual"] < 1e-12

    def test_F2_hermitian_per_block(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=12)
        assert d2d.verify_hermitian_per_block(tol=1e-12)

    def test_D_squared_matches_factorized(self):
        """D^2 explicit spectrum matches DiscreteDiskDirac.squared_eigenvalues."""
        N_rho, a, N_phi = 15, 0.3, 12
        d2d = DiscreteDirac2D(N_rho=N_rho, a=a, N_phi=N_phi)
        disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_phi)
        factorized = disk.squared_eigenvalues()
        result = d2d.verify_D_squared_matches_factorized(factorized, tol=1e-10)
        assert result["passed"], result

    def test_eigenvalues_paired(self):
        """D spectrum: each positive eigenvalue has a negative partner."""
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        evals = d2d.eigenvalues()
        pos = sorted([v for v in evals if v > 1e-12])
        neg_mag = sorted([-v for v in evals if v < -1e-12])
        assert len(pos) == len(neg_mag)
        for p, n in zip(pos, neg_mag):
            assert np.isclose(p, n, atol=1e-10)

    def test_eigenvalues_count(self):
        d2d = DiscreteDirac2D(N_rho=10, a=0.5, N_phi=8)
        evals = d2d.eigenvalues()
        assert len(evals) == d2d.hilbert_dim

    def test_lowest_positive_eigenvalue_is_real_positive(self):
        d2d = DiscreteDirac2D(N_rho=50, a=0.2, N_phi=24)
        lam_min = d2d.lowest_positive_eigenvalue()
        assert lam_min > 0


class TestDiscreteDirac2DContinuumLimit:
    def test_lowest_mode_approaches_pi_over_R(self):
        """|lambda_min| -> pi/R as substrate refines (j_{1/2, 1} = pi)."""
        # Three progressively finer panels at the same R = 10
        panels = [
            (50, 0.2, 24),
            (100, 0.1, 48),
            (200, 0.05, 96),
        ]
        errs = []
        for Nr, a, Np in panels:
            d2d = DiscreteDirac2D(N_rho=Nr, a=a, N_phi=Np)
            lam_min = d2d.lowest_positive_eigenvalue()
            R = Nr * a
            target = np.pi / R
            rel_err = (lam_min - target) / target
            errs.append(abs(rel_err))
        # Errors should DECREASE as substrate refines
        assert errs[0] > errs[1] > errs[2], errs
        # Finest panel: within 1%
        assert errs[-1] < 0.01, f"Finest panel rel_err = {errs[-1]}"

    def test_lowest_mode_finest_within_one_percent(self):
        """At UV-fine panel, |lambda_min| within 1% of pi/R."""
        d2d = DiscreteDirac2D(N_rho=200, a=0.05, N_phi=96)
        lam_min = d2d.lowest_positive_eigenvalue()
        target = np.pi / 10.0
        rel_err = (lam_min - target) / target
        assert abs(rel_err) < 0.01, f"rel_err = {rel_err}"


# =================================================================
# VariableWarpDirac (G4-4b-a first move)
# =================================================================


class TestVariableWarpDiracConstruction:
    def test_smooth_tip_construction(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        assert var.disk is disk
        assert var.sphere is sphere
        assert var.r_h == 2.0
        assert len(var.warp_profile) == 10

    def test_smooth_tip_warp_at_tip(self):
        """At rho_1 = a small, r(rho_1) ~ r_h."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.05, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        # rho_1 = 0.05, r_h = 2.0
        # r(0.05) = 2 * sqrt(1 + (0.05/2)^2) = 2 * sqrt(1.000625) ~ 2.000625
        assert np.isclose(var.warp_profile[0], 2.0006249023, atol=1e-8)

    def test_smooth_tip_warp_asymptotic(self):
        """At rho >> r_h, r(rho) -> rho (Schwarzschild far field)."""
        disk = DiscreteDiskDirac(N_rho=200, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=1.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=1.0)
        # At rho_200 = 100, r ~ sqrt(1 + 100^2) ~ 100
        assert np.isclose(var.warp_profile[-1], np.sqrt(1 + 100**2), atol=1e-6)

    def test_constant_construction(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.constant(disk, sphere, r_h=2.0)
        assert np.all(var.warp_profile == 2.0)

    def test_invalid_warp_profile_length(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        with pytest.raises(ValueError):
            VariableWarpDirac(
                disk=disk, sphere=sphere,
                warp_profile=np.ones(5), r_h=2.0,
            )

    def test_invalid_warp_profile_negative(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        with pytest.raises(ValueError):
            VariableWarpDirac(
                disk=disk, sphere=sphere,
                warp_profile=np.zeros(10), r_h=2.0,
            )

    def test_invalid_r_h(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        with pytest.raises(ValueError):
            VariableWarpDirac(
                disk=disk, sphere=sphere,
                warp_profile=np.ones(10), r_h=0.0,
            )

    def test_rho_array(self):
        disk = DiscreteDiskDirac(N_rho=5, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        # rho_k = (k+1) * a for k = 0, ..., 4
        expected = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
        assert np.allclose(var.rho_array, expected)

    def test_warp_derivative_finite_and_positive(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.2, N_phi=12)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        deriv = var.warp_derivative_over_warp()
        assert np.all(np.isfinite(deriv))
        assert np.all(deriv > 0)

    def test_H_block_hermitian(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        for n in range(3):
            for k_phi in range(8):
                H = var.H_block(n, k_phi)
                assert np.allclose(H, H.T)

    def test_heat_trace_positive_decreasing(self):
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=8)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        var = VariableWarpDirac.smooth_tip(disk, sphere, r_h=2.0)
        K_01 = var.heat_trace(0.1)
        K_05 = var.heat_trace(0.5)
        K_10 = var.heat_trace(1.0)
        assert K_01 > K_05 > K_10 > 0


class TestF4TipRegular:
    def test_F4_passes_smooth_tip(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.2, N_phi=12)
        res = verify_F4_tip_regular(disk=disk, r_h=2.0)
        assert res["passed"], res

    def test_F4_finite_substrate(self):
        disk = DiscreteDiskDirac(N_rho=20, a=0.2, N_phi=12)
        res = verify_F4_tip_regular(disk=disk, r_h=2.0)
        assert res["deriv_array_finite"]

    def test_F4_tip_value_finite_and_small(self):
        """At a << r_h, r'/r(rho_1) is finite and O(a/r_h^2).

        With centered FD on the warp profile (G4-4b-d fix), the boundary
        site uses one-sided FD which has O(a) error. The structural F4
        content is regularity (no singularity), not exact analytical
        match. Check value is finite and within an order of magnitude
        of the analytical a/r_h^2.
        """
        disk = DiscreteDiskDirac(N_rho=30, a=0.1, N_phi=12)
        r_h = 2.0
        res = verify_F4_tip_regular(disk=disk, r_h=r_h)
        # a / r_h^2 = 0.025 analytical; FD boundary value can differ by ~1.5x
        assert np.isfinite(res["deriv_at_rho_1"])
        assert res["deriv_at_rho_1"] > 0
        # Within factor 3 of analytical: structural regularity check
        assert 0.025 / 3 < res["deriv_at_rho_1"] < 0.025 * 3

    def test_F4_max_value_bounded(self):
        """r'/r has a maximum (at rho = r_h, value 1/(2 r_h))."""
        disk = DiscreteDiskDirac(N_rho=100, a=0.1, N_phi=12)
        r_h = 2.0
        res = verify_F4_tip_regular(disk=disk, r_h=r_h)
        # max r'/r = 1/(2 r_h) = 0.25
        assert res["max_deriv_value"] < 0.30  # within 20% of analytical max


class TestF6RiemannianLimit:
    def test_F6_bit_exact_small_panel(self):
        """LOAD-BEARING: variable @ constant warp = constant warp bit-exact."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=6)
        sphere = S2DiracSpectrum(l_max=2, r_h=1.5)
        res = verify_F6_riemannian_limit(
            disk=disk, sphere=sphere, r_h=1.5,
            t_values=[0.05, 0.1, 0.5, 1.0], tol=1e-10,
        )
        assert res["all_passed"], res

    def test_F6_machine_precision(self):
        """Relative error ~ float64 sum noise (< 1e-13)."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=6)
        sphere = S2DiracSpectrum(l_max=2, r_h=1.5)
        res = verify_F6_riemannian_limit(
            disk=disk, sphere=sphere, r_h=1.5,
            t_values=[0.1, 0.5], tol=1e-12,
        )
        for t in ["0.1", "0.5"]:
            assert res[t]["rel_err"] < 1e-12

    def test_F6_mismatched_r_h_raises(self):
        """Verify raises when sphere.r_h != r_h."""
        disk = DiscreteDiskDirac(N_rho=10, a=0.5, N_phi=6)
        sphere = S2DiracSpectrum(l_max=2, r_h=2.0)
        with pytest.raises(ValueError):
            verify_F6_riemannian_limit(
                disk=disk, sphere=sphere, r_h=1.5, t_values=[0.1],
            )


class TestF7FactorizationLoss:
    def test_F7_positive_delta_at_variable_warp(self):
        """K_var > K_const at smooth-tip variable warp (sign + for r > r_h)."""
        disk = DiscreteDiskDirac(N_rho=20, a=0.3, N_phi=12)
        sphere = S2DiracSpectrum(l_max=3, r_h=2.0)
        res = verify_F7_factorization_loss(
            disk=disk, sphere=sphere, r_h=2.0,
            t_values=[0.1, 0.5, 1.0],
        )
        for t in ["0.1", "0.5", "1.0"]:
            assert res[t]["sign_positive"], f"t={t}: {res[t]}"

    def test_F7_monotonic_in_warp_variation(self):
        """Delta_fact grows as r_h decreases (more variation)."""
        disk = DiscreteDiskDirac(N_rho=20, a=0.3, N_phi=12)
        # Same disk, different r_h => different warp variation
        deltas = []
        for r_h in [5.0, 2.0, 1.0]:
            sphere = S2DiracSpectrum(l_max=3, r_h=r_h)
            res = verify_F7_factorization_loss(
                disk=disk, sphere=sphere, r_h=r_h, t_values=[0.5],
            )
            deltas.append(res["0.5"]["Delta_fact"])
        # Monotonic increase as r_h decreases
        assert deltas[0] < deltas[1] < deltas[2], deltas

    def test_F7_ratio_above_one(self):
        disk = DiscreteDiskDirac(N_rho=15, a=0.3, N_phi=12)
        sphere = S2DiracSpectrum(l_max=3, r_h=2.0)
        res = verify_F7_factorization_loss(
            disk=disk, sphere=sphere, r_h=2.0, t_values=[0.5],
        )
        assert res["0.5"]["ratio"] > 1.0
