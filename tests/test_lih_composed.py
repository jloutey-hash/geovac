"""
End-to-end integration tests for the LiH composed graph pipeline.

Tests the full architecture: core screening + Level 4 valence + PES + spectroscopy.
The primary success criterion is R_eq > 2.7 bohr (proving core screening fixes
the LCAO equilibrium geometry problem where R_eq = 2.5 vs experiment 3.015).
"""

import time
import numpy as np
import pytest

from geovac.lih_composed import LiHComposedSolver, EXPT, MU_LIH_AU, HARTREE_TO_CM


@pytest.fixture(scope="module")
def solver():
    """Run the full pipeline once and share across tests."""
    s = LiHComposedSolver(
        l_max=2, n_alpha=100, zeff_mode='screened',
        use_pk=True, verbose=True,
    )

    t0 = time.time()

    # Step 1: Core
    s.solve_core()

    # Step 2: PES scan (18 R-points)
    R_grid = np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])
    s.scan_pes(R_grid=R_grid, n_Re=300)

    # Step 3: Morse fit
    s.fit_spectroscopic_constants()

    # Step 4: Nuclear lattice
    s.build_nuclear_lattice(J_max=10)

    total_time = time.time() - t0
    print(f"\n  Total pipeline time: {total_time:.1f}s")

    # Print full summary
    s.timings['total'] = total_time
    s._print_summary()

    return s


class TestCoreScreening:
    """Tests for the Li+ core solve."""

    def test_core_energy_accuracy(self, solver: LiHComposedSolver) -> None:
        """E_core should be within 1% of exact Li+ energy (-7.2799 Ha)."""
        err = abs(solver.E_core - EXPT['E_Li_plus']) / abs(EXPT['E_Li_plus'])
        print(f"  Core energy error: {err*100:.2f}%")
        assert err < 0.01, (
            f"Core energy {solver.E_core:.6f} too far from "
            f"exact {EXPT['E_Li_plus']}"
        )

    def test_z_eff_limits(self, solver: LiHComposedSolver) -> None:
        """Z_eff should be 3 at r=0 and ~1 at large r."""
        z0 = solver.core.z_eff(0.0)
        z_large = solver.core.z_eff(10.0)
        assert abs(z0 - 3.0) < 0.01, f"Z_eff(0) = {z0}, expected 3.0"
        assert abs(z_large - 1.0) < 0.1, f"Z_eff(10) = {z_large}, expected ~1.0"


class TestPES:
    """Tests for the potential energy surface."""

    def test_pes_has_minimum(self, solver: LiHComposedSolver) -> None:
        """PES should have a clear minimum in R in [2, 6] bohr."""
        R_eq = solver.pes_result['R_eq']
        assert 2.0 <= R_eq <= 6.0, (
            f"R_eq = {R_eq} outside expected range [2, 6]"
        )

    def test_r_eq_improvement_over_lcao(self, solver: LiHComposedSolver) -> None:
        """
        R_eq should be > 2.7 bohr -- the KEY test.
        LCAO gave 2.5; experiment is 3.015.
        Core screening + PK pseudopotential should push R_eq outward.
        """
        R_eq = solver.pes_result['R_eq']
        print(f"  R_eq = {R_eq:.3f} bohr"
              f" (target: > 2.7, expt: 3.015, LCAO: 2.5)")
        assert R_eq > 2.7, (
            f"R_eq = {R_eq:.3f} not improved over LCAO (2.5). "
            f"Expected > 2.7 from core screening + PK."
        )

    def test_molecule_is_bound(self, solver: LiHComposedSolver) -> None:
        """D_e should be positive (molecule is bound)."""
        D_e = solver.pes_result['D_e']
        print(f"  D_e = {D_e:.6f} Ha (expt: {EXPT['D_e']})")
        assert D_e > 0, f"D_e = {D_e}, molecule is not bound!"


class TestSpectroscopy:
    """Tests for spectroscopic constants."""

    def test_omega_e_range(self, solver: LiHComposedSolver) -> None:
        """omega_e should be 1000-2000 cm-1 (expt: 1405)."""
        omega = solver.spectro['omega_e']
        print(f"  omega_e = {omega:.1f} cm-1 (expt: {EXPT['omega_e']})")
        assert 1000 < omega < 2000, f"omega_e = {omega} outside [1000, 2000]"

    def test_B_e_range(self, solver: LiHComposedSolver) -> None:
        """B_e should be 4-12 cm-1 (expt: 7.51)."""
        B_e = solver.spectro['B_e']
        print(f"  B_e = {B_e:.2f} cm-1 (expt: {EXPT['B_e']})")
        assert 4 < B_e < 12, f"B_e = {B_e} outside [4, 12]"

    def test_d_e_positive(self, solver: LiHComposedSolver) -> None:
        """Fitted D_e from Morse should be positive."""
        D_e = solver.spectro['D_e']
        assert D_e > 0, f"Fitted D_e = {D_e}, should be positive"


class TestNuclearLattice:
    """Tests for the nuclear rovibrational spectrum."""

    def test_vmax_sufficient(self, solver: LiHComposedSolver) -> None:
        """Should have v_max > 5 bound vibrational levels."""
        v_max = solver.nuclear.vib.v_max
        print(f"  v_max = {v_max}")
        assert v_max > 5, f"v_max = {v_max}, expected > 5"

    def test_fundamental_frequency_range(self, solver: LiHComposedSolver) -> None:
        """Fundamental v01 should be 800-1800 cm-1 (expt: 1359)."""
        E0 = solver.nuclear.vib.morse_energy(0)
        E1 = solver.nuclear.vib.morse_energy(1)
        nu01 = (E1 - E0) * HARTREE_TO_CM
        print(f"  v01 = {nu01:.1f} cm-1 (expt: ~1359)")
        assert 800 < nu01 < 1800, f"v01 = {nu01} outside [800, 1800]"

    def test_rotational_structure(self, solver: LiHComposedSolver) -> None:
        """v=0 J=1 should be slightly above v=0 J=0."""
        E_00 = solver.nuclear.rovibrational_energy(0, 0)
        E_01 = solver.nuclear.rovibrational_energy(0, 1)
        gap = (E_01 - E_00) * HARTREE_TO_CM
        print(f"  J=0->1 gap = {gap:.2f} cm-1"
              f" (expect ~ 2 B_e = {2*solver.spectro['B_e']:.1f})")
        assert gap > 0, "J=1 should be above J=0"
        assert (abs(gap - 2 * solver.spectro['B_e'])
                / (2 * solver.spectro['B_e']) < 0.2)


class TestPerformance:
    """Performance tests."""

    def test_total_time(self, solver: LiHComposedSolver) -> None:
        """Entire pipeline should complete in < 15 minutes."""
        total = solver.timings.get('total', 0)
        print(f"  Total time: {total:.1f}s")
        assert total < 900, (
            f"Pipeline took {total:.0f}s, expected < 900s (15 min)"
        )
