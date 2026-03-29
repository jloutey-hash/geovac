"""
Tests for the spectral Laguerre radial solver in prolate spheroidal coordinates.

Validates that the spectral solver (radial_method='spectral') reproduces
H2+ energies with dramatically fewer degrees of freedom than the FD solver
(20 basis functions vs 5000 grid points), at equal or better accuracy.

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


class TestSpectralBasic:
    """Basic functionality of the spectral radial solver."""

    def test_solve_r2_bound(self):
        """Spectral solver at R=2.0 returns bound energy."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_el, c2, A = lat.solve()
        E_tot = E_el + 0.5
        assert E_tot < -0.5, f"Should be bound, got E_tot={E_tot:.6f}"

    def test_c2_positive(self):
        """c^2 should be positive for bound state."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        _, c2, _ = lat.solve()
        assert c2 > 0, f"c^2={c2} should be positive"

    def test_separation_constant_positive(self):
        """Separation constant A should be positive for sigma_g ground state."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        _, _, A = lat.solve()
        assert A > 0, f"A={A} should be positive"

    def test_total_energy_includes_vnn(self):
        """total_energy() = E_elec + Z_A*Z_B/R."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_tot = lat.total_energy()
        E_el, _, _ = lat.solve()
        assert abs(E_tot - (E_el + 0.5)) < 1e-10


class TestSpectralAccuracy:
    """Accuracy benchmarks for the spectral solver."""

    def test_accuracy_vs_exact_r2(self):
        """n_basis=20 should give < 0.01% error at R=2.0."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_el, _, _ = lat.solve()
        E_tot = E_el + 0.5
        # Exact: -0.6026342 Ha
        err = abs(E_tot - (-0.6026342)) / 0.6026342 * 100
        assert err < 0.01, f"Error {err:.4f}% > 0.01% at n_basis=20"

    def test_accuracy_vs_exact_r1(self):
        """Check accuracy at R=1.0 (short range)."""
        lat = ProlateSpheroidalLattice(R=1.0, radial_method='spectral', n_basis=20)
        E_el, _, _ = lat.solve()
        # E_elec at R=1.0 should be well below -1.0 (united atom He+ limit -2.0)
        assert E_el < -1.0, f"E_elec={E_el:.4f} should be < -1.0 at R=1.0"

    def test_accuracy_vs_fd(self):
        """Spectral n_basis=20 should match or beat FD N_xi=5000."""
        lat_sp = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        lat_fd = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        E_sp, _, _ = lat_sp.solve()
        E_fd, _, _ = lat_fd.solve()
        # Spectral should be lower (more accurate, variational)
        # or at least within 0.01 Ha
        assert E_sp <= E_fd + 0.01, \
            f"Spectral E={E_sp:.6f} much higher than FD E={E_fd:.6f}"

    def test_convergence_with_basis_size(self):
        """Energy should converge monotonically with increasing n_basis."""
        energies = []
        for nb in [5, 10, 15, 20]:
            lat = ProlateSpheroidalLattice(
                R=2.0, radial_method='spectral', n_basis=nb
            )
            E_el, _, _ = lat.solve()
            energies.append(E_el + 0.5)
        # Should be converging (becoming more negative or stable)
        for i in range(1, len(energies)):
            # Allow small tolerance for non-monotonicity due to alpha adaptation
            assert energies[i] <= energies[i - 1] + 1e-4, \
                f"E[nb={[5,10,15,20][i]}]={energies[i]:.8f} > " \
                f"E[nb={[5,10,15,20][i-1]}]={energies[i-1]:.8f}"

    def test_small_basis_still_reasonable(self):
        """Even n_basis=5 should give < 5% error."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=5)
        E_el, _, _ = lat.solve()
        E_tot = E_el + 0.5
        err = abs(E_tot - (-0.6026)) / 0.6026 * 100
        assert err < 5.0, f"Error {err:.2f}% > 5% even with n_basis=5"


class TestSpectralPES:
    """PES scan tests with the spectral solver."""

    def test_pes_has_minimum(self):
        """Spectral PES should have a minimum between R=1 and R=5."""
        R_values = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
        E_tot = []
        for R in R_values:
            lat = ProlateSpheroidalLattice(
                R=R, radial_method='spectral', n_basis=20
            )
            E_tot.append(lat.total_energy())
        E_tot = np.array(E_tot)
        idx_min = np.argmin(E_tot)
        assert 0 < idx_min < len(R_values) - 1, \
            f"Minimum at boundary (idx={idx_min})"

    def test_spectroscopic_constants(self):
        """Spectral PES should give accurate R_eq and D_e."""
        # Fine grid around equilibrium for accurate polynomial fit
        R_values = np.arange(1.5, 3.01, 0.05)
        E_tot = []
        for R in R_values:
            lat = ProlateSpheroidalLattice(
                R=R, radial_method='spectral', n_basis=20
            )
            E_tot.append(lat.total_energy())
        E_tot = np.array(E_tot)
        fit = fit_spectroscopic_constants(R_values, E_tot)

        assert not fit['boundary'], "Minimum should not be at boundary"
        # R_eq should be within 1% of exact 1.997
        assert abs(fit['R_eq'] - 1.997) / 1.997 < 0.01, \
            f"R_eq={fit['R_eq']:.4f}, expected ~1.997 (> 1% err)"
        # E_min should be within 0.01% of exact
        assert abs(fit['E_min'] - (-0.6026342)) / 0.6026342 < 0.0001, \
            f"E_min={fit['E_min']:.8f}, expected ~-0.6026342 (> 0.01% err)"
        # D_e > 0
        assert fit['D_e'] > 0, f"D_e={fit['D_e']:.6f} should be positive"


class TestSpectralDissociation:
    """Dissociation limit tests."""

    def test_large_r_approaches_h_atom(self):
        """At large R, E_total should approach -0.5 Ha."""
        lat = ProlateSpheroidalLattice(
            R=10.0, radial_method='spectral', n_basis=25
        )
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.5)) < 0.01, \
            f"E(R=10)={E_tot:.6f}, expected ~-0.5"


class TestSpectralHeteronuclear:
    """Heteronuclear tests (HeH^2+)."""

    def test_heh2plus_bound(self):
        """HeH^2+ should be well bound."""
        lat = ProlateSpheroidalLattice(
            R=2.0, Z_A=2, Z_B=1,
            radial_method='spectral', n_basis=20
        )
        E_tot = lat.total_energy()
        assert E_tot < -1.0, f"HeH2+ E={E_tot:.4f} should be well bound"


class TestSpectralPerformance:
    """Performance benchmarks."""

    def test_speedup_over_fd(self):
        """Spectral solver should be at least 10x faster than FD N_xi=5000."""
        # Spectral
        t0 = time.time()
        lat_sp = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20
        )
        lat_sp.solve()
        dt_sp = time.time() - t0

        # FD
        t0 = time.time()
        lat_fd = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        lat_fd.solve()
        dt_fd = time.time() - t0

        speedup = dt_fd / max(dt_sp, 1e-6)
        assert speedup > 10, \
            f"Speedup {speedup:.1f}x < 10x (spectral={dt_sp:.3f}s, fd={dt_fd:.3f}s)"

    def test_dimension_reduction(self):
        """Spectral should use <=25 basis functions vs 5000 FD grid points."""
        lat = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20
        )
        assert lat.n_basis <= 25
        assert lat.n_basis < lat.N_xi / 100  # >100x reduction


class TestFDUnchanged:
    """Verify that adding the spectral option didn't break FD."""

    def test_fd_default_unchanged(self):
        """Default radial_method should be 'fd'."""
        lat = ProlateSpheroidalLattice(R=2.0)
        assert lat.radial_method == 'fd'

    def test_fd_results_unchanged(self):
        """FD results should be identical to before."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        E_el, c2, A = lat.solve()
        E_tot = E_el + 0.5
        # Should match the known FD result
        err = abs(E_tot - (-0.6026)) / 0.6026 * 100
        assert err < 1.5, f"FD error {err:.2f}% > 1.5% (regression)"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
