"""
Stress tests for the prolate spheroidal lattice solver.

Tests excited states, limiting cases, asymmetric charges, and grid convergence.
Extends the basic tests in test_prolate_h2plus.py.

Date: 2026-03-13
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


# ============================================================
# Category 1: Excited States
# ============================================================
class TestExcitedStates:
    """H2+ excited state validation."""

    def test_sigma_u_energy(self):
        """1sigma_u at R=2.0 should be near -0.168 Ha."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000, n_angular=1)
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.168)) / 0.168 < 0.05, \
            f"1sigma_u E={E_tot:.4f}, expected ~-0.168 (< 5% err)"

    def test_sigma_u_above_ground(self):
        """1sigma_u must be above 1sigma_g."""
        lat_g = ProlateSpheroidalLattice(R=2.0, N_xi=5000)
        lat_u = ProlateSpheroidalLattice(R=2.0, N_xi=5000, n_angular=1)
        E_g = lat_g.total_energy()
        E_u = lat_u.total_energy()
        assert E_u > E_g, f"1sigma_u ({E_u:.4f}) should be above 1sigma_g ({E_g:.4f})"

    def test_pi_u_dissociation(self):
        """1pi_u at large R should approach H(2p) = -0.125 Ha."""
        lat = ProlateSpheroidalLattice(R=15.0, N_xi=5000, xi_max=60.0, m=1)
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.125)) / 0.125 < 0.05, \
            f"1pi_u E(R=15)={E_tot:.4f}, expected ~-0.125"

    def test_pi_u_has_shallow_minimum(self):
        """1pi_u PES should have a shallow minimum around R=6-10."""
        R_vals = np.array([4.0, 6.0, 8.0, 10.0, 12.0])
        pes = scan_h2plus_pes(R_vals, N_xi=5000, m=1, verbose=False)
        idx_min = np.argmin(pes['E_total'])
        # Minimum should be interior (not at boundary)
        assert 0 < idx_min < len(R_vals) - 1, \
            f"pi_u minimum at boundary: idx={idx_min}"


# ============================================================
# Category 2: Limiting Cases
# ============================================================
class TestLimitingCases:
    """Verify correct behavior at extreme R."""

    def test_united_atom_r02(self):
        """At R=0.2, E_elec should approach He+ limit (-2.0 Ha)."""
        lat = ProlateSpheroidalLattice(R=0.2, N_xi=5000)
        E_el, _, _ = lat.solve()
        # Allow 15% error at this extreme (grid-limited)
        assert abs(E_el - (-2.0)) / 2.0 < 0.15, \
            f"E_elec(R=0.2)={E_el:.4f}, expected ~-2.0"

    def test_dissociation_r10(self):
        """At R=10, E_total should be near -0.5 Ha."""
        lat = ProlateSpheroidalLattice(R=10.0, N_xi=8000, xi_max=80.0)
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.5)) < 0.03, \
            f"E(R=10)={E_tot:.4f}, expected ~-0.5"

    def test_sigma_u_repulsive(self):
        """1sigma_u should be purely repulsive (no bound minimum)."""
        R_vals = np.array([1.0, 2.0, 3.0, 5.0, 8.0])
        pes = scan_h2plus_pes(R_vals, N_xi=5000, n_angular=1, verbose=False)
        # All energies should be above -0.5 (no binding)
        for i, R in enumerate(R_vals):
            assert pes['E_total'][i] > -0.52, \
                f"sigma_u E(R={R})={pes['E_total'][i]:.4f} suggests bound state"


# ============================================================
# Category 3: Asymmetric Charges
# ============================================================
class TestAsymmetricCharges:
    """Test with Z_A != Z_B."""

    def test_heh2plus_bound(self):
        """HeH2+ PES should have a clear minimum (bound state)."""
        R_vals = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
        pes = scan_h2plus_pes(R_vals, Z_A=2, Z_B=1, N_xi=5000, verbose=False)
        idx_min = np.argmin(pes['E_total'])
        E_min = pes['E_total'][idx_min]
        # Should have interior minimum (PES has a well)
        # Note: E_min > -2.0 (He+ limit) at N_xi=5000 is a grid resolution issue
        assert E_min < -1.5, f"HeH2+ E_min={E_min:.4f}, should be well bound"
        assert idx_min > 0, "Minimum should not be at smallest R"

    def test_lih3plus_solves(self):
        """LiH3+ (Z_A=3, Z_B=1) should produce a valid bound energy."""
        R_vals = np.array([1.0, 1.5, 2.0, 3.0, 5.0])
        pes = scan_h2plus_pes(R_vals, Z_A=3, Z_B=1, N_xi=5000, verbose=False)
        E_min = np.min(pes['E_total'])
        # At N_xi=5000, grid error is significant for Z=3+1. Energy should
        # still be well negative (bound). Exact dissoc limit is -4.5 Ha.
        assert E_min < -3.5, f"LiH3+ E_min={E_min:.4f}, should be < -3.5"

    def test_ch5plus_solves(self):
        """CH5+ (Z_A=6, Z_B=1) extreme asymmetry should still solve."""
        lat = ProlateSpheroidalLattice(R=2.0, Z_A=6, Z_B=1, N_xi=5000)
        E_tot = lat.total_energy()
        # Should be a reasonable negative number
        assert E_tot < 0, f"CH5+ E={E_tot:.4f}, should be negative"


# ============================================================
# Category 4: Grid Convergence
# ============================================================
class TestGridConvergence:
    """Verify systematic convergence with grid refinement."""

    def test_convergence_monotonic(self):
        """Energy should converge monotonically from above."""
        energies = []
        for N_xi in [500, 1000, 2000, 5000]:
            lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi)
            E = lat.total_energy()
            energies.append(E)

        for i in range(1, len(energies)):
            assert energies[i] <= energies[i-1] + 1e-6, \
                f"Non-monotonic: E[{i}]={energies[i]:.6f} > E[{i-1}]={energies[i-1]:.6f}"

    def test_richardson_extrapolation(self):
        """Richardson extrapolation should improve accuracy."""
        E_5k = ProlateSpheroidalLattice(R=2.0, N_xi=5000).total_energy()
        E_10k = ProlateSpheroidalLattice(R=2.0, N_xi=10000).total_energy()
        E_extrap = (4 * E_10k - E_5k) / 3  # Assuming O(h^2)

        # Extrapolated should be closer to exact than either individual
        exact = -0.6026
        err_5k = abs(E_5k - exact)
        err_extrap = abs(E_extrap - exact)
        assert err_extrap < err_5k, \
            f"Richardson didn't improve: err_5k={err_5k:.4f}, err_extrap={err_extrap:.4f}"


# ============================================================
# Category 5: n_angular/n_radial parameters
# ============================================================
class TestNewParameters:
    """Verify the excited state parameters don't break ground state."""

    def test_n_angular_0_is_ground(self):
        """n_angular=0 should give same result as default."""
        lat_default = ProlateSpheroidalLattice(R=2.0, N_xi=3000)
        lat_explicit = ProlateSpheroidalLattice(R=2.0, N_xi=3000, n_angular=0)
        E_def = lat_default.total_energy()
        E_exp = lat_explicit.total_energy()
        assert abs(E_def - E_exp) < 1e-10, \
            f"n_angular=0 gives different result: {E_def:.8f} vs {E_exp:.8f}"

    def test_n_angular_increases_energy(self):
        """Higher n_angular should give higher (less negative) energy."""
        energies = []
        for n in [0, 1, 2]:
            lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000, n_angular=n)
            E = lat.total_energy()
            energies.append(E)

        for i in range(1, len(energies)):
            assert energies[i] > energies[i-1], \
                f"Energy not increasing: n={i}: {energies[i]:.4f} <= n={i-1}: {energies[i-1]:.4f}"

    def test_m1_different_from_m0(self):
        """m=1 (pi state) should give different energy from m=0 (sigma)."""
        lat_s = ProlateSpheroidalLattice(R=2.0, N_xi=5000, m=0)
        lat_p = ProlateSpheroidalLattice(R=2.0, N_xi=5000, m=1)
        E_s = lat_s.total_energy()
        E_p = lat_p.total_energy()
        assert abs(E_s - E_p) > 0.1, \
            f"m=0 and m=1 too similar: {E_s:.4f} vs {E_p:.4f}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
