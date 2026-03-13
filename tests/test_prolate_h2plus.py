"""
Tests for the prolate spheroidal lattice H2+ solver.

Validates that the separated-equation solver in prolate spheroidal
coordinates reproduces the exact H2+ potential energy surface
without free parameters.

Exact references:
  R_eq = 1.997 bohr
  E(R_eq) = -0.6026 Ha
  D_e = 0.1026 Ha
  E(H atom) = -0.5000 Ha
"""
import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    scan_h2plus_pes,
    fit_spectroscopic_constants,
)


class TestProlateSpheroidalBasic:
    """Basic sanity checks."""

    def test_solve_r2(self):
        """Ground state at R=2.0 returns reasonable energy."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=3000)
        E_el, c2, A = lat.solve()
        E_tot = E_el + 0.5
        assert E_tot < -0.5, "Should be bound"
        assert abs(E_tot - (-0.6026)) / 0.6026 < 0.02, \
            f"E_total={E_tot:.4f}, expected ~-0.6026 (< 2% err)"

    def test_c2_positive(self):
        """c^2 should be positive for bound state."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=2000)
        _, c2, _ = lat.solve()
        assert c2 > 0, f"c^2={c2} should be positive"

    def test_separation_constant(self):
        """Angular separation constant A should be positive for ground state."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=2000)
        _, _, A = lat.solve()
        assert A > 0, f"A={A} should be positive for sigma_g"

    def test_total_energy(self):
        """total_energy() includes V_NN."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=2000)
        E_tot = lat.total_energy()
        E_el, _, _ = lat.solve()
        assert abs(E_tot - (E_el + 0.5)) < 1e-10


class TestRadialConvergence:
    """Grid convergence of the radial solver."""

    def test_convergence_direction(self):
        """Energy should decrease (become more negative) with finer grid."""
        energies = []
        for N_xi in [500, 1000, 2000, 5000]:
            lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi)
            E_el, _, _ = lat.solve()
            energies.append(E_el + 0.5)
        # Each should be lower than the previous (converging from above)
        for i in range(1, len(energies)):
            assert energies[i] <= energies[i - 1] + 1e-6, \
                f"E[{i}]={energies[i]:.6f} > E[{i-1}]={energies[i-1]:.6f}"

    def test_5000_grid_accuracy(self):
        """N_xi=5000 should give < 1.5% error at R=2.0."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000)
        E_el, _, _ = lat.solve()
        E_tot = E_el + 0.5
        err = abs(E_tot - (-0.6026)) / 0.6026 * 100
        assert err < 1.5, f"Error {err:.2f}% > 1.5%"


class TestPES:
    """Potential energy surface tests."""

    def test_pes_has_minimum(self):
        """PES should have a minimum between R=1 and R=5."""
        R_values = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
        pes = scan_h2plus_pes(R_values, N_xi=3000, verbose=False)
        idx_min = np.argmin(pes['E_total'])
        assert idx_min > 0 and idx_min < len(R_values) - 1, \
            "Minimum should be interior"

    def test_bound_at_equilibrium(self):
        """E(R_eq) should be below the dissociation limit -0.5 Ha."""
        R_values = np.array([1.5, 2.0, 2.5, 3.0])
        pes = scan_h2plus_pes(R_values, N_xi=3000, verbose=False)
        E_min = np.min(pes['E_total'])
        assert E_min < -0.5, f"E_min={E_min:.4f} > -0.5 (not bound)"

    def test_spectroscopic_constants(self):
        """Fit should give reasonable R_eq and D_e."""
        R_values = np.arange(1.0, 5.01, 0.25)
        pes = scan_h2plus_pes(R_values, N_xi=5000, verbose=False)
        fit = fit_spectroscopic_constants(pes['R'], pes['E_total'])

        assert not fit['boundary'], "Minimum should not be at boundary"
        assert 1.5 < fit['R_eq'] < 3.0, \
            f"R_eq={fit['R_eq']:.3f} out of range [1.5, 3.0]"
        assert fit['D_e'] > 0, f"D_e={fit['D_e']:.4f} should be positive"
        assert abs(fit['E_min'] - (-0.6026)) / 0.6026 < 0.02, \
            f"E_min error > 2%"


class TestDissociationLimit:
    """Tests for the R -> infinity behavior."""

    def test_large_r_approaches_h_atom(self):
        """At large R, E_total should approach -0.5 Ha."""
        # Need large xi_max and N_xi for the extended wavefunction
        lat = ProlateSpheroidalLattice(R=10.0, N_xi=8000, xi_max=80.0)
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.5)) < 0.03, \
            f"E(R=10)={E_tot:.4f}, expected ~-0.5"


class TestHeteronuclear:
    """Test with Z_A != Z_B (e.g., HeH2+)."""

    def test_heh2plus_bound(self):
        """HeH2+ (Z_A=2, Z_B=1) should have E < -1.0 (bound relative to H+p)."""
        lat = ProlateSpheroidalLattice(R=2.0, Z_A=2, Z_B=1, N_xi=5000)
        E_tot = lat.total_energy()
        # One-electron HeH2+: E should be significantly negative
        # The united atom limit (Li2+) has E = -Z^2/2 = -4.5 Ha
        # At R=2 it should be somewhere between -4.5 and -2.0
        assert E_tot < -1.0, f"HeH2+ E={E_tot:.4f} should be well bound"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
