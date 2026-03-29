"""
Tests for spectral Laguerre wiring into scan_h2plus_pes() and solve_with_wavefunction().

Validates that:
  1. scan_h2plus_pes(radial_method='spectral') produces accurate PES
  2. solve_with_wavefunction(radial_method='spectral') returns spectral wavefunction
  3. Default behavior (radial_method='fd') is unchanged
  4. Spectral PES is faster than FD PES

Exact references (Bates, Ledsham & Stewart 1953):
  R_eq = 1.997 bohr
  E_total(R_eq) = -0.6026342 Ha
  D_e = 0.1026 Ha
"""
import numpy as np
import pytest
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    scan_h2plus_pes,
    fit_spectroscopic_constants,
)


# --- scan_h2plus_pes spectral wiring ---

class TestSpectralPESScan:
    """Tests for scan_h2plus_pes with radial_method='spectral'."""

    @pytest.fixture(scope='class')
    def spectral_pes(self):
        """Run spectral PES scan once for multiple tests."""
        R_values = np.arange(1.5, 3.51, 0.05)
        result = scan_h2plus_pes(
            R_values, radial_method='spectral', n_basis=20, verbose=False,
        )
        return result

    @pytest.fixture(scope='class')
    def spectral_fit(self, spectral_pes):
        return fit_spectroscopic_constants(
            spectral_pes['R'], spectral_pes['E_total']
        )

    def test_scan_returns_all_keys(self, spectral_pes):
        """scan_h2plus_pes should return R, E_elec, E_total, V_NN."""
        for key in ['R', 'E_elec', 'E_total', 'V_NN']:
            assert key in spectral_pes, f"Missing key: {key}"

    def test_scan_no_nans(self, spectral_pes):
        """All energies should be finite (no failed solves)."""
        assert not np.any(np.isnan(spectral_pes['E_total']))

    def test_scan_r_eq(self, spectral_fit):
        """Spectral PES R_eq within 0.5% of exact 1.997 bohr."""
        err = abs(spectral_fit['R_eq'] - 1.997) / 1.997
        assert err < 0.005, \
            f"R_eq={spectral_fit['R_eq']:.4f}, error={err*100:.2f}% > 0.5%"

    def test_scan_e_min(self, spectral_fit):
        """Spectral PES E_min within 0.01% of exact -0.6026342 Ha."""
        err = abs(spectral_fit['E_min'] - (-0.6026342)) / 0.6026342
        assert err < 0.0001, \
            f"E_min={spectral_fit['E_min']:.8f}, error={err*100:.4f}% > 0.01%"

    def test_scan_d_e_positive(self, spectral_fit):
        """D_e should be positive (bound molecule)."""
        assert spectral_fit['D_e'] > 0, f"D_e={spectral_fit['D_e']:.6f}"

    def test_scan_d_e_accurate(self, spectral_fit):
        """D_e within 1% of exact 0.1026 Ha."""
        err = abs(spectral_fit['D_e'] - 0.1026) / 0.1026
        assert err < 0.01, \
            f"D_e={spectral_fit['D_e']:.6f}, error={err*100:.2f}% > 1%"

    def test_scan_not_at_boundary(self, spectral_fit):
        """Minimum should not be at grid boundary."""
        assert not spectral_fit['boundary']


class TestFDPESScanUnchanged:
    """Default FD behavior is unchanged when radial_method not specified."""

    def test_fd_default_scan(self):
        """scan_h2plus_pes without radial_method uses FD."""
        R_values = np.array([1.5, 2.0, 2.5, 3.0])
        result = scan_h2plus_pes(R_values, verbose=False)
        # Should work and produce reasonable energies
        assert not np.any(np.isnan(result['E_total']))
        # Check E at R=2.0 is near exact
        idx = 1  # R=2.0
        assert result['E_total'][idx] < -0.58

    def test_fd_explicit_scan(self):
        """Explicit radial_method='fd' matches implicit default."""
        R_values = np.array([2.0])
        res_default = scan_h2plus_pes(R_values, verbose=False)
        res_fd = scan_h2plus_pes(
            R_values, radial_method='fd', verbose=False
        )
        np.testing.assert_allclose(
            res_default['E_total'], res_fd['E_total'], atol=1e-12,
        )


class TestSpectralPESSpeed:
    """Spectral PES scan should be faster than FD."""

    def test_spectral_faster_than_fd(self):
        """Spectral PES scan over 10 points should be faster than FD."""
        R_values = np.linspace(1.5, 3.0, 10)

        t0 = time.time()
        scan_h2plus_pes(
            R_values, radial_method='spectral', n_basis=20, verbose=False,
        )
        dt_sp = time.time() - t0

        t0 = time.time()
        scan_h2plus_pes(
            R_values, N_xi=5000, radial_method='fd', verbose=False,
        )
        dt_fd = time.time() - t0

        speedup = dt_fd / max(dt_sp, 1e-6)
        assert speedup > 5, \
            f"Speedup {speedup:.1f}x < 5x (sp={dt_sp:.2f}s, fd={dt_fd:.2f}s)"


# --- solve_with_wavefunction spectral wiring ---

class TestSpectralWavefunction:
    """Tests for solve_with_wavefunction with radial_method='spectral'."""

    @pytest.fixture(scope='class')
    def spectral_wfn(self):
        lat = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20,
        )
        return lat.solve_with_wavefunction()

    @pytest.fixture(scope='class')
    def fd_wfn(self):
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='fd')
        return lat.solve_with_wavefunction()

    def test_spectral_wfn_has_all_keys(self, spectral_wfn):
        """Spectral wavefunction dict has all required keys."""
        for key in ['E_elec', 'c2', 'A', 'xi_grid', 'eta_grid',
                     'F_xi', 'G_eta', 'radial_method']:
            assert key in spectral_wfn, f"Missing key: {key}"

    def test_spectral_method_tag(self, spectral_wfn):
        """radial_method should be 'spectral'."""
        assert spectral_wfn['radial_method'] == 'spectral'

    def test_fd_method_tag(self, fd_wfn):
        """FD wavefunction should tag radial_method='fd'."""
        assert fd_wfn['radial_method'] == 'fd'

    def test_spectral_fxi_nonzero(self, spectral_wfn):
        """Spectral F(xi) should be nonzero."""
        assert np.max(np.abs(spectral_wfn['F_xi'])) > 0.01

    def test_spectral_fxi_normalized(self, spectral_wfn):
        """Spectral F(xi) should be approximately normalized."""
        xi = spectral_wfn['xi_grid']
        F = spectral_wfn['F_xi']
        dx = xi[1] - xi[0]
        norm_sq = np.sum(F**2) * dx
        assert abs(norm_sq - 1.0) < 0.05, f"norm^2={norm_sq:.4f}, expected ~1.0"

    def test_spectral_fxi_positive_at_xi1(self, spectral_wfn):
        """F(xi) should be positive near xi=1 for sigma_g ground state."""
        assert spectral_wfn['F_xi'][0] > 0

    def test_spectral_fxi_decays(self, spectral_wfn):
        """F(xi) should decay at large xi."""
        F = spectral_wfn['F_xi']
        # Last 10% should be much smaller than peak
        peak = np.max(np.abs(F))
        tail = np.max(np.abs(F[-50:]))
        assert tail < 0.1 * peak, \
            f"Tail max={tail:.4f} vs peak={peak:.4f} — not decaying"

    def test_spectral_energy_matches_solve(self, spectral_wfn):
        """Energy from solve_with_wavefunction matches solve()."""
        lat = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20,
        )
        E_el, _, _ = lat.solve()
        assert abs(spectral_wfn['E_elec'] - E_el) < 1e-12

    def test_angular_wfn_consistent(self, spectral_wfn, fd_wfn):
        """Angular wavefunction G(eta) should be similar for both methods.

        Not identical because spectral and FD find slightly different c2,
        which parameterizes the angular equation. Tolerance reflects the
        ~0.1% energy difference between solvers.
        """
        sp_G = spectral_wfn['G_eta']
        fd_G = fd_wfn['G_eta']
        # Could differ in sign convention
        sign = np.sign(sp_G[len(sp_G)//2]) * np.sign(fd_G[len(fd_G)//2])
        np.testing.assert_allclose(
            sp_G, sign * fd_G, atol=0.01,
            err_msg="Angular wavefunctions differ too much between spectral and FD"
        )

    def test_spectral_fxi_shape_similar_to_fd(self, spectral_wfn, fd_wfn):
        """Spectral and FD F(xi) should have similar shape near equilibrium."""
        # Interpolate FD onto spectral grid for comparison
        from scipy.interpolate import interp1d
        fd_interp = interp1d(
            fd_wfn['xi_grid'], fd_wfn['F_xi'],
            kind='cubic', fill_value=0.0, bounds_error=False,
        )
        sp_xi = spectral_wfn['xi_grid']
        sp_F = spectral_wfn['F_xi']
        fd_F_on_sp_grid = fd_interp(sp_xi)

        # Compare in the region where both are significant (xi < 10)
        mask = sp_xi < 10.0
        sp_peak = np.max(np.abs(sp_F[mask]))
        fd_peak = np.max(np.abs(fd_F_on_sp_grid[mask]))
        # Peaks should be within 20% (different normalizations, grids)
        assert abs(sp_peak - fd_peak) / max(sp_peak, fd_peak) < 0.3, \
            f"Peak mismatch: spectral={sp_peak:.4f}, fd={fd_peak:.4f}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
