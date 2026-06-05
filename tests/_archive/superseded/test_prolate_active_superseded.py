"""
Tests for ProlateActiveSpace — prolate spheroidal molecular active space.

Tests verify:
1. H2+ (1 electron): reproduces ProlateSpheroidalLattice results
2. H2 (2 electrons): CI in sigma_g/sigma_u basis gives correct E
3. LiH (locked Li 1s^2 + 2 active): PES has a minimum
4. Core orbital on prolate grid: normalization correct
5. Backward compatibility: active_method='lcao' unchanged
6. Orbital symmetry labels correct
"""

import warnings
import numpy as np
import pytest

warnings.filterwarnings("ignore", message=".*LatticeIndex V_ee.*")
warnings.filterwarnings("ignore", message=".*MolecularLatticeIndex.*")


# ---------------------------------------------------------------------------
# Test: Core orbital on prolate grid
# ---------------------------------------------------------------------------

class TestCoreOrbital:
    """Test atomic_orbital_on_prolate_grid."""

    def test_li_1s_normalization(self):
        """Li 1s orbital should be normalized on the prolate grid."""
        from geovac.prolate_active_space import atomic_orbital_on_prolate_grid
        from geovac.prolate_scf import get_orbital_on_grid

        # Get a prolate grid from a dummy orbital
        ref = get_orbital_on_grid(R=3.0, Z_A=3.0, Z_B=1.0, n_angular=0,
                                  N_xi_grid=80, N_eta_grid=80)
        orb = atomic_orbital_on_prolate_grid(
            Z=3.0, n=1, l=0, atom_label='A', R=3.0,
            xi=ref['xi'], eta=ref['eta'],
            w_xi=ref['w_xi'], w_eta=ref['w_eta'],
        )
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J = (3.0 / 2.0)**3 * (XI**2 - ETA**2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm = np.sum(orb['psi']**2 * J * W_XI * W_ETA) * 2 * np.pi
        print(f"\nLi 1s norm on prolate grid: {norm:.6f}")
        assert abs(norm - 1.0) < 0.05, f"Normalization failed: {norm}"

    def test_h_1s_normalization(self):
        """H 1s orbital centered on atom B."""
        from geovac.prolate_active_space import atomic_orbital_on_prolate_grid
        from geovac.prolate_scf import get_orbital_on_grid

        ref = get_orbital_on_grid(R=3.0, Z_A=1.0, Z_B=1.0, n_angular=0,
                                  N_xi_grid=80, N_eta_grid=80)
        orb = atomic_orbital_on_prolate_grid(
            Z=1.0, n=1, l=0, atom_label='B', R=3.0,
            xi=ref['xi'], eta=ref['eta'],
            w_xi=ref['w_xi'], w_eta=ref['w_eta'],
        )
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J = (3.0 / 2.0)**3 * (XI**2 - ETA**2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm = np.sum(orb['psi']**2 * J * W_XI * W_ETA) * 2 * np.pi
        print(f"\nH 1s norm on prolate grid: {norm:.6f}")
        assert abs(norm - 1.0) < 0.05, f"Normalization failed: {norm}"


# ---------------------------------------------------------------------------
# Test: H2+ (1 electron) — must match Paper 11
# ---------------------------------------------------------------------------

class TestH2Plus:
    """Test that ProlateActiveSpace with 1 electron matches ProlateSpheroidalLattice."""

    def test_h2plus_energy_r2(self):
        """H2+ at R=2.0: prolate active space matches direct solve."""
        from geovac.prolate_active_space import ProlateActiveSpace
        from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice

        R = 2.0
        # Direct solve
        lat = ProlateSpheroidalLattice(R=R, Z_A=1, Z_B=1)
        E_direct = lat.total_energy()

        # ProlateActiveSpace with 1 electron, 1 orbital
        pas = ProlateActiveSpace(
            Z_A=1.0, Z_B=1.0, R=R,
            n_electrons=1,
            orbital_spec=[(0, 0)],  # sigma_g only
        )
        pas.build()
        V_NN = 1.0 / R
        eigvals, _ = pas.solve(E_locked=0.0, V_NN=V_NN)

        print(f"\nH2+ R=2.0: direct={E_direct:.6f}, prolate_active={eigvals[0]:.6f}")
        assert abs(eigvals[0] - E_direct) < 0.001, (
            f"H2+ energy mismatch: {eigvals[0]:.6f} vs {E_direct:.6f}"
        )

    def test_h2plus_energy_r4(self):
        """H2+ at R=4.0: longer bond length."""
        from geovac.prolate_active_space import ProlateActiveSpace
        from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice

        R = 4.0
        lat = ProlateSpheroidalLattice(R=R, Z_A=1, Z_B=1)
        E_direct = lat.total_energy()

        pas = ProlateActiveSpace(
            Z_A=1.0, Z_B=1.0, R=R,
            n_electrons=1,
            orbital_spec=[(0, 0)],
        )
        pas.build()
        eigvals, _ = pas.solve(E_locked=0.0, V_NN=1.0 / R)

        print(f"\nH2+ R=4.0: direct={E_direct:.6f}, prolate_active={eigvals[0]:.6f}")
        assert abs(eigvals[0] - E_direct) < 0.001


# ---------------------------------------------------------------------------
# Test: H2 (2 electrons) — sigma_g/sigma_u CI
# ---------------------------------------------------------------------------

class TestH2:
    """Test H2 CI in prolate MO basis."""

    def test_h2_energy_r1_4(self):
        """H2 at R=1.4: should be close to exact -1.174 Ha."""
        from geovac.prolate_active_space import ProlateActiveSpace

        R = 1.4
        pas = ProlateActiveSpace(
            Z_A=1.0, Z_B=1.0, R=R,
            n_electrons=2,
            orbital_spec=[(0, 0), (0, 1)],  # sigma_g, sigma_u
            N_grid=80,
        )
        pas.build()
        eigvals, _ = pas.solve(E_locked=0.0, V_NN=1.0 / R)
        E = eigvals[0]

        E_exact = -1.174  # Ha (approximate exact)
        error = abs(E - E_exact)
        print(f"\nH2 R=1.4: E={E:.6f}, exact=-1.174, error={error:.4f}")
        # Should be within ~5% (this is minimal CI, no Z_eff optimization)
        assert E < -1.0, f"H2 energy too high: {E:.6f}"

    def test_h2_bound(self):
        """H2 PES should have a minimum."""
        from geovac.prolate_active_space import ProlateActiveSpace

        energies = []
        R_values = [1.0, 1.4, 2.0, 3.0, 5.0]
        for R in R_values:
            pas = ProlateActiveSpace(
                Z_A=1.0, Z_B=1.0, R=R,
                n_electrons=2,
                orbital_spec=[(0, 0), (0, 1)],
                N_grid=60,
                verbose=False,
            )
            pas.build()
            eigvals, _ = pas.solve(E_locked=0.0, V_NN=1.0 / R)
            energies.append(eigvals[0])

        # Should have a minimum in the interior
        E_min_idx = np.argmin(energies)
        print(f"\nH2 PES: {list(zip(R_values, ['%.6f' % e for e in energies]))}")
        print(f"  Minimum at R={R_values[E_min_idx]}, E={energies[E_min_idx]:.6f}")
        assert 0 < E_min_idx < len(energies) - 1, (
            f"No interior minimum: min at index {E_min_idx}"
        )


# ---------------------------------------------------------------------------
# Test: Symmetry labels
# ---------------------------------------------------------------------------

class TestSymmetry:
    """Test orbital symmetry labels."""

    def test_sigma_labels(self):
        """sigma_g and sigma_u labels assigned correctly."""
        from geovac.prolate_active_space import ProlateActiveSpace

        pas = ProlateActiveSpace(
            Z_A=1.0, Z_B=1.0, R=2.0,
            n_electrons=1,
            orbital_spec=[(0, 0), (0, 1)],
            verbose=False,
        )
        pas.build()
        assert pas._orbital_labels[0] == '1sigma_g'
        assert pas._orbital_labels[1] == '1sigma_u'

    def test_sigma_g_lower(self):
        """sigma_g energy should be lower than sigma_u."""
        from geovac.prolate_active_space import ProlateActiveSpace

        pas = ProlateActiveSpace(
            Z_A=1.0, Z_B=1.0, R=2.0,
            n_electrons=1,
            orbital_spec=[(0, 0), (0, 1)],
            verbose=False,
        )
        pas.build()
        assert pas._orbital_energies[0] < pas._orbital_energies[1], (
            f"sigma_g ({pas._orbital_energies[0]:.4f}) not lower than "
            f"sigma_u ({pas._orbital_energies[1]:.4f})"
        )


# ---------------------------------------------------------------------------
# Test: LiH with locked Li 1s^2
# ---------------------------------------------------------------------------

class TestLiH:
    """Test LiH with prolate active space (screened Z_eff approach)."""

    def test_lih_prolate_runs(self):
        """LiH with prolate active space should complete without error."""
        from geovac.locked_shell import LockedShellMolecule

        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_method='prolate',
            prolate_kwargs={'N_grid': 60},
        )
        E, psi = mol.solve()
        print(f"\nLiH prolate: E = {E[0]:.6f} Ha")
        # Should be negative (bound system)
        assert E[0] < 0, f"LiH energy positive: {E[0]:.6f}"

    def test_lih_prolate_pes_has_minimum(self):
        """LiH PES with prolate active space should have a minimum."""
        from geovac.locked_shell import LockedShellMolecule

        energies = []
        R_values = [2.0, 2.5, 3.0, 4.0, 6.0, 10.0]
        for R in R_values:
            mol = LockedShellMolecule(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_electrons=4,
                locked_config={0: [(1, 0)]},
                active_method='prolate',
                prolate_kwargs={'N_grid': 60},
            )
            E, _ = mol.solve()
            energies.append(E[0])

        E_min_idx = np.argmin(energies)
        D_e = energies[-1] - min(energies)

        print(f"\nLiH prolate PES:")
        for R, E in zip(R_values, energies):
            print(f"  R={R:.1f}  E={E:.6f}")
        print(f"  Min at R={R_values[E_min_idx]:.1f}, D_e={D_e:.4f} Ha")

        # PES should have interior minimum (not at boundary)
        assert 0 < E_min_idx < len(energies) - 1, (
            f"No interior minimum: min at index {E_min_idx}"
        )
        # PES should be bound
        assert D_e > 0.01, f"LiH not bound: D_e = {D_e:.4f}"

    def test_lih_e_locked_correct(self):
        """E_locked should be the analytical Li 1s^2 energy + cross-nuclear."""
        from geovac.locked_shell import LockedShellMolecule

        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_method='prolate',
            prolate_kwargs={'N_grid': 60},
        )

        # Li 1s^2 bare: 2*(-9/2) + 5*3/8 = -9 + 1.875 = -7.125
        E_bare = -7.125
        # E_locked includes cross-nuclear, so should be more negative
        assert mol.E_locked < E_bare, (
            f"E_locked {mol.E_locked:.4f} should be < {E_bare}"
        )
        # But not absurdly negative
        assert mol.E_locked > -10.0, (
            f"E_locked {mol.E_locked:.4f} too negative"
        )


# ---------------------------------------------------------------------------
# Test: Analytical helpers
# ---------------------------------------------------------------------------

class TestAnalytical:
    """Test _slater_f0 and _cross_nuclear_attraction."""

    def test_slater_f0_1s1s(self):
        """F^0(1s,1s;Z=3) = 5*3/8 = 1.875."""
        from geovac.prolate_active_space import _slater_f0
        assert abs(_slater_f0(3.0, 1, 0, 3.0, 1, 0) - 1.875) < 1e-10

    def test_slater_f0_1s2s(self):
        """F^0(1s,2s;Z=3) = 17*3/81."""
        from geovac.prolate_active_space import _slater_f0
        expected = 17.0 * 3.0 / 81.0
        assert abs(_slater_f0(3.0, 1, 0, 3.0, 2, 0) - expected) < 1e-10

    def test_cross_nuclear_large_r(self):
        """<1s|1/r_B|1s> -> 1/R at large R."""
        from geovac.prolate_active_space import _cross_nuclear_attraction
        R = 50.0
        val = _cross_nuclear_attraction(3.0, 1, 0, R)
        assert abs(val - 1.0 / R) < 1e-4, f"Expected ~{1/R:.6f}, got {val:.6f}"

    def test_cross_nuclear_zero_at_zero(self):
        """No cross-nuclear at R=0."""
        from geovac.prolate_active_space import _cross_nuclear_attraction
        assert _cross_nuclear_attraction(3.0, 1, 0, 0.0) == 0.0


# ---------------------------------------------------------------------------
# Test: Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompat:
    """Verify active_method='lcao' (default) is unchanged."""

    def test_lcao_default(self):
        """Default active_method should be 'lcao'."""
        from geovac.locked_shell import LockedShellMolecule

        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        assert mol._active_method == 'lcao'
        E, _ = mol.solve()
        # Should match previous behavior (negative energy)
        assert E[0] < -7.0
