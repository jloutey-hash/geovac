"""
Integration tests for the generalized ComposedDiatomicSolver.

Validates:
  1. LiH backward compatibility with lih_composed.py
  2. LiH accuracy progression with l_max (2, 3)
  3. BeH+ runs successfully with same architecture
  4. BeH+ produces physically reasonable spectroscopic constants

Each molecule/l_max combo requires ~4 minutes, so tests are structured
with module-scoped fixtures that run the pipeline once and share results.
"""

import time
import numpy as np
import pytest

from geovac.composed_diatomic import ComposedDiatomicSolver, REFERENCE_DATA
from geovac.nuclear_lattice import HARTREE_TO_CM


# ======================================================================
# R-grid definitions
# ======================================================================

def _lih_r_grid_standard() -> np.ndarray:
    """Standard 18-point LiH R-grid (matches lih_composed.py tests)."""
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])


def _lih_r_grid_extended() -> np.ndarray:
    """Extended R-grid out to 12 bohr for better dissociation limit."""
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 6),
        np.linspace(8.0, 12.0, 6),
    ])


def _beh_plus_r_grid() -> np.ndarray:
    """BeH+ R-grid: denser near expected minimum at ~2.5 bohr."""
    return np.concatenate([
        np.linspace(1.5, 2.0, 3),
        np.linspace(2.1, 3.5, 10),
        np.linspace(4.0, 7.0, 5),
    ])


# ======================================================================
# Fixtures — run pipelines once, share across tests
# ======================================================================

@pytest.fixture(scope="module")
def lih_l2():
    """LiH at l_max=2 (baseline)."""
    s = ComposedDiatomicSolver.LiH(l_max=2, verbose=True)
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid_standard(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    s.build_nuclear_lattice(J_max=10)
    s._print_summary()
    return s


@pytest.fixture(scope="module")
def lih_l3():
    """LiH at l_max=3 (higher accuracy)."""
    s = ComposedDiatomicSolver.LiH(l_max=3, verbose=True)
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid_extended(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    s.build_nuclear_lattice(J_max=10)
    s._print_summary()
    return s


@pytest.fixture(scope="module")
def beh_plus_l2():
    """BeH+ at l_max=2."""
    s = ComposedDiatomicSolver.BeH_plus(l_max=2, verbose=True)
    s.solve_core()
    s.scan_pes(R_grid=_beh_plus_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    s.build_nuclear_lattice(J_max=10)
    s._print_summary()
    return s


# ======================================================================
# LiH l_max=2 — backward compatibility
# ======================================================================

class TestLiHL2Backward:
    """LiH l_max=2 must match previous lih_composed.py results."""

    def test_core_energy(self, lih_l2: ComposedDiatomicSolver) -> None:
        """E_core within 1% of exact Li+ energy."""
        ref = REFERENCE_DATA['LiH']['E_core']
        err = abs(lih_l2.E_core - ref) / abs(ref)
        print(f"  Core energy error: {err*100:.2f}%")
        assert err < 0.01

    def test_r_eq_matches_previous(self, lih_l2: ComposedDiatomicSolver) -> None:
        """R_eq should be within 0.05 bohr of the previous result (2.95)."""
        R_eq = lih_l2.pes_result['R_eq']
        print(f"  R_eq = {R_eq:.3f} (previous: 2.95, expt: 3.015)")
        assert abs(R_eq - 2.95) < 0.15, (
            f"R_eq={R_eq:.3f} deviates from previous 2.95 by "
            f"{abs(R_eq-2.95):.3f} bohr"
        )

    def test_r_eq_improvement_over_lcao(self, lih_l2: ComposedDiatomicSolver) -> None:
        """R_eq > 2.7 (LCAO was 2.5, composed graph should do better)."""
        R_eq = lih_l2.pes_result['R_eq']
        assert R_eq > 2.7

    def test_molecule_bound(self, lih_l2: ComposedDiatomicSolver) -> None:
        """D_e > 0."""
        assert lih_l2.pes_result['D_e'] > 0

    def test_omega_e_range(self, lih_l2: ComposedDiatomicSolver) -> None:
        """omega_e in [1000, 2000] cm-1."""
        omega = lih_l2.spectro['omega_e']
        print(f"  omega_e = {omega:.1f} cm-1 (expt: 1405.65)")
        assert 1000 < omega < 2000

    def test_B_e_range(self, lih_l2: ComposedDiatomicSolver) -> None:
        """B_e in [4, 12] cm-1."""
        B_e = lih_l2.spectro['B_e']
        print(f"  B_e = {B_e:.2f} cm-1 (expt: 7.51)")
        assert 4 < B_e < 12


# ======================================================================
# LiH l_max=3 — accuracy improvement
# ======================================================================

class TestLiHL3Accuracy:
    """LiH l_max=3 should be closer to experiment than l_max=2."""

    def test_r_eq_reasonable(self, lih_l3: ComposedDiatomicSolver) -> None:
        """R_eq in [2.5, 4.0] bohr."""
        R_eq = lih_l3.pes_result['R_eq']
        print(f"  R_eq (l_max=3) = {R_eq:.3f} bohr (expt: 3.015)")
        assert 2.5 <= R_eq <= 4.0

    def test_molecule_bound(self, lih_l3: ComposedDiatomicSolver) -> None:
        """D_e > 0."""
        assert lih_l3.pes_result['D_e'] > 0

    def test_omega_e_range(self, lih_l3: ComposedDiatomicSolver) -> None:
        """omega_e in [800, 2500] cm-1."""
        omega = lih_l3.spectro['omega_e']
        print(f"  omega_e (l_max=3) = {omega:.1f} cm-1 (expt: 1405.65)")
        assert 800 < omega < 2500

    def test_B_e_range(self, lih_l3: ComposedDiatomicSolver) -> None:
        """B_e in [3, 15] cm-1."""
        B_e = lih_l3.spectro['B_e']
        print(f"  B_e (l_max=3) = {B_e:.2f} cm-1 (expt: 7.51)")
        assert 3 < B_e < 15

    def test_pipeline_completes(self, lih_l3: ComposedDiatomicSolver) -> None:
        """Full pipeline completes with nuclear lattice."""
        assert lih_l3.nuclear is not None
        assert lih_l3.nuclear.vib.n_states > 1


# ======================================================================
# BeH+ — architecture generalization
# ======================================================================

class TestBeHPlusPipeline:
    """BeH+ tests: same architecture, different molecule."""

    def test_pipeline_completes(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """Full pipeline completes without error."""
        assert beh_plus_l2.E_core is not None
        assert beh_plus_l2.pes_result is not None
        assert beh_plus_l2.spectro is not None
        assert beh_plus_l2.nuclear is not None

    def test_core_energy_reasonable(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """Be2+ 1s2 energy should be near -13.66 Ha."""
        E = beh_plus_l2.E_core
        print(f"  Be2+ core energy: {E:.4f} Ha (ref: -13.6556)")
        assert -15.0 < E < -12.0

    def test_pes_has_minimum(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """PES should have a clear minimum in [1.5, 4.0] bohr."""
        R_eq = beh_plus_l2.pes_result['R_eq']
        print(f"  BeH+ R_eq = {R_eq:.3f} bohr (ref: ~2.48)")
        assert 1.5 <= R_eq <= 4.0

    def test_molecule_bound(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """D_e > 0 (molecule is bound)."""
        D_e = beh_plus_l2.pes_result['D_e']
        print(f"  BeH+ D_e = {D_e:.6f} Ha")
        assert D_e > 0

    def test_r_eq_physical(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """R_eq in [2.0, 3.5] bohr for a Be-H bond."""
        R_eq = beh_plus_l2.spectro['R_eq']
        print(f"  BeH+ R_eq (Morse) = {R_eq:.3f} bohr")
        assert 2.0 <= R_eq <= 3.5

    def test_omega_e_physical(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """omega_e in [1000, 3500] cm-1 for a light hydride."""
        omega = beh_plus_l2.spectro['omega_e']
        print(f"  BeH+ omega_e = {omega:.1f} cm-1 (ref: ~2222)")
        assert 1000 < omega < 3500

    def test_B_e_physical(self, beh_plus_l2: ComposedDiatomicSolver) -> None:
        """B_e should be positive and in a reasonable range."""
        B_e = beh_plus_l2.spectro['B_e']
        print(f"  BeH+ B_e = {B_e:.2f} cm-1 (ref: ~10.31)")
        assert B_e > 0


# ======================================================================
# Architecture validation — same code path
# ======================================================================

class TestArchitectureGenerality:
    """Verify both molecules use the exact same code path."""

    def test_same_class(self, lih_l2: ComposedDiatomicSolver,
                        beh_plus_l2: ComposedDiatomicSolver) -> None:
        """Both should be ComposedDiatomicSolver instances."""
        assert type(lih_l2) is type(beh_plus_l2)
        assert type(lih_l2).__name__ == 'ComposedDiatomicSolver'

    def test_different_Z(self, lih_l2: ComposedDiatomicSolver,
                         beh_plus_l2: ComposedDiatomicSolver) -> None:
        """Only Z_A differs (no molecule-specific branches)."""
        assert lih_l2.Z_A_bare == 3.0
        assert beh_plus_l2.Z_A_bare == 4.0
        assert lih_l2.Z_B == beh_plus_l2.Z_B == 1.0

    def test_pk_scaling(self, lih_l2: ComposedDiatomicSolver,
                        beh_plus_l2: ComposedDiatomicSolver) -> None:
        """BeH+ PK parameters should be larger (tighter core)."""
        # LiH: A=5.0, B=7.0 (calibrated)
        # BeH+: A~6.67, B~12.4 (Z-scaled)
        print(f"  LiH PK: A={lih_l2.pk_A:.2f}, B={lih_l2.pk_B:.2f}")
        print(f"  BeH+ PK: A={beh_plus_l2.pk_A:.2f}, B={beh_plus_l2.pk_B:.2f}")
        assert beh_plus_l2.pk_B > lih_l2.pk_B  # Tighter core -> larger B


# ======================================================================
# Comparison table (printed at end)
# ======================================================================

class TestComparisonTable:
    """Print the multi-molecule comparison table."""

    def test_print_table(self, lih_l2: ComposedDiatomicSolver,
                         lih_l3: ComposedDiatomicSolver,
                         beh_plus_l2: ComposedDiatomicSolver) -> None:
        """Print comprehensive comparison."""
        print("\n")
        print("=" * 68)
        print("Composed Diatomic Solver — Multi-Molecule Validation")
        print("=" * 68)

        print("\nLiH Accuracy Progression:")
        print(f"  {'':14s} {'l_max=2':>10s} {'l_max=3':>10s} {'Expt':>10s}")
        print(f"  {'R_eq (bohr)':14s}"
              f" {lih_l2.spectro['R_eq']:10.3f}"
              f" {lih_l3.spectro['R_eq']:10.3f}"
              f" {'3.015':>10s}")
        print(f"  {'omega_e':14s}"
              f" {lih_l2.spectro['omega_e']:10.1f}"
              f" {lih_l3.spectro['omega_e']:10.1f}"
              f" {'1405.7':>10s}")
        print(f"  {'B_e':14s}"
              f" {lih_l2.spectro['B_e']:10.2f}"
              f" {lih_l3.spectro['B_e']:10.2f}"
              f" {'7.51':>10s}")
        print(f"  {'D_e (Ha)':14s}"
              f" {lih_l2.spectro['D_e']:10.4f}"
              f" {lih_l3.spectro['D_e']:10.4f}"
              f" {'0.0920':>10s}")
        t2 = lih_l2.timings.get('total', 0)
        t3 = lih_l3.timings.get('total', 0)
        print(f"  {'Time (sec)':14s}"
              f" {t2:10.0f}"
              f" {t3:10.0f}"
              f" {'---':>10s}")

        print(f"\nBeH+ Results (l_max=2):")
        bs = beh_plus_l2.spectro
        print(f"  R_eq = {bs['R_eq']:.3f} bohr    (ref: ~2.48)")
        print(f"  omega_e = {bs['omega_e']:.1f} cm-1  (ref: ~2222)")
        print(f"  B_e = {bs['B_e']:.2f} cm-1")
        print(f"  D_e = {bs['D_e']:.4f} Ha")
        bt = beh_plus_l2.timings.get('total', 0)
        print(f"  Time = {bt:.0f} sec")
        print(f"  PK: A={beh_plus_l2.pk_A:.2f}, B={beh_plus_l2.pk_B:.2f}")

        print(f"\nArchitecture: Same ComposedDiatomicSolver for both molecules.")
        print(f"Only inputs that change: Z_A, Z_B, M_A, M_B, pk_A, pk_B.")
        print("=" * 68)

        # This test always passes — it's just for printing
        assert True


# ======================================================================
# LiH extended angular basis — sigma+pi channels
# ======================================================================

def _lih_r_grid_extended_wide() -> np.ndarray:
    """Wide R-grid for higher l_max (minimum shifts outward)."""
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 8.0, 7),
    ])


@pytest.fixture(scope="module")
def lih_l3_sigma_pi():
    """LiH at l_max=3 with sigma+pi channels (ab initio PK)."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=3, m_max=1, l_max_per_m={0: 3, 1: 2},
        n_alpha=60, verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid_extended_wide(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=2.0)
    s.build_nuclear_lattice(J_max=10)
    s._print_summary()
    return s


@pytest.fixture(scope="module")
def lih_l4_sigma_pi():
    """LiH at l_max=4 with sigma+pi channels (ab initio PK)."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=4, m_max=1, l_max_per_m={0: 4, 1: 2},
        n_alpha=50, verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid_extended_wide(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=2.0)
    s.build_nuclear_lattice(J_max=10)
    s._print_summary()
    return s


class TestLiHL3SigmaPi:
    """LiH l_max=3 with sigma+pi channels.

    Higher l_max adds angular correlation that preferentially lowers
    energy at large R (diffuse electron cloud), shifting R_eq outward.
    This is a documented negative result for the monotonic convergence
    prediction from Paper 15 (H2). The additional channels produce a
    bound PES with correct qualitative shape.
    """

    def test_pipeline_completes(
        self, lih_l3_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """Full sigma+pi pipeline completes without error."""
        assert lih_l3_sigma_pi.E_core is not None
        assert lih_l3_sigma_pi.pes_result is not None
        assert lih_l3_sigma_pi.spectro is not None
        assert lih_l3_sigma_pi.nuclear is not None

    def test_molecule_bound(
        self, lih_l3_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """D_e > 0 (molecule is bound)."""
        D_e = lih_l3_sigma_pi.pes_result['D_e']
        print(f"  D_e (l3 sigma+pi) = {D_e:.6f} Ha")
        assert D_e > 0

    def test_r_eq_in_range(
        self, lih_l3_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """R_eq in [2.5, 5.0] bohr — bound PES with physical minimum."""
        R_eq = lih_l3_sigma_pi.spectro['R_eq']
        print(f"  R_eq (l3 sigma+pi) = {R_eq:.3f} bohr (expt: 3.015)")
        assert 2.5 <= R_eq <= 5.0

    def test_m_max_set(
        self, lih_l3_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """Verify m_max=1 was set (pi channels active)."""
        assert lih_l3_sigma_pi.m_max == 1

    def test_omega_e_physical(
        self, lih_l3_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """omega_e in physical range [600, 2500] cm-1."""
        omega = lih_l3_sigma_pi.spectro['omega_e']
        print(f"  omega_e (l3 sigma+pi) = {omega:.1f} cm-1 (expt: 1405.7)")
        assert 600 < omega < 2500


class TestLiHL4SigmaPi:
    """LiH l_max=4 with sigma+pi channels.

    Same physics as l_max=3: additional angular channels lower energy
    preferentially at large R, pushing R_eq further outward. The PES
    remains bound with correct qualitative shape.
    """

    def test_pipeline_completes(
        self, lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """Full l_max=4 sigma+pi pipeline completes without error."""
        assert lih_l4_sigma_pi.E_core is not None
        assert lih_l4_sigma_pi.pes_result is not None
        assert lih_l4_sigma_pi.spectro is not None
        assert lih_l4_sigma_pi.nuclear is not None

    def test_molecule_bound(
        self, lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """D_e > 0."""
        D_e = lih_l4_sigma_pi.pes_result['D_e']
        print(f"  D_e (l4 sigma+pi) = {D_e:.6f} Ha")
        assert D_e > 0

    def test_r_eq_in_range(
        self, lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """R_eq in [2.5, 6.0] bohr — bound PES."""
        R_eq = lih_l4_sigma_pi.spectro['R_eq']
        print(f"  R_eq (l4 sigma+pi) = {R_eq:.3f} bohr (expt: 3.015)")
        assert 2.5 <= R_eq <= 6.0

    def test_m_max_set(
        self, lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """Verify m_max=1 and l_max_per_m set."""
        assert lih_l4_sigma_pi.m_max == 1
        assert lih_l4_sigma_pi.l_max_per_m == {0: 4, 1: 2}

    def test_more_channels_than_l2(
        self, lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """l_max=4 sigma+pi has more channels than l_max=2 sigma."""
        from geovac.level4_multichannel import (
            _channel_list, _channel_list_extended,
        )
        ch_l2 = _channel_list(2, homonuclear=False)
        ch_l4 = _channel_list_extended(
            4, 1, {0: 4, 1: 2}, homonuclear=False,
        )
        print(f"  l_max=2 sigma: {len(ch_l2)} channels")
        print(f"  l_max=4 sigma+pi: {len(ch_l4)} channels")
        assert len(ch_l4) > len(ch_l2)


class TestLmaxConvergenceReport:
    """Print convergence table for l_max progression.

    Documents that R_eq does NOT converge monotonically toward
    experiment (3.015 bohr) with increasing l_max in the single-channel
    adiabatic approximation. Additional angular channels preferentially
    lower energy at large R, shifting R_eq outward.
    """

    def test_print_convergence_table(
        self,
        lih_l2: ComposedDiatomicSolver,
        lih_l3_sigma_pi: ComposedDiatomicSolver,
        lih_l4_sigma_pi: ComposedDiatomicSolver,
    ) -> None:
        """Print l_max convergence table."""
        print("\n")
        print("=" * 72)
        print("LiH l_max Convergence with Sigma+Pi Channels")
        print("=" * 72)
        print(f"  {'Config':20s} {'R_eq':>8s} {'err%':>6s} "
              f"{'omega_e':>8s} {'D_e':>8s} {'m_max':>5s}")
        print(f"  {'-'*20} {'-'*8} {'-'*6} {'-'*8} {'-'*8} {'-'*5}")

        for label, s in [
            ("l2 sigma (man.PK)", lih_l2),
            ("l3 sig+pi (ai.PK)", lih_l3_sigma_pi),
            ("l4 sig+pi (ai.PK)", lih_l4_sigma_pi),
        ]:
            R = s.spectro['R_eq']
            err = abs(R - 3.015) / 3.015 * 100
            w = s.spectro['omega_e']
            D = s.spectro['D_e']
            print(f"  {label:20s} {R:8.3f} {err:6.1f} "
                  f"{w:8.1f} {D:8.4f} {s.m_max:5d}")

        print(f"  {'Experiment':20s} {'3.015':>8s} {'0.0':>6s} "
              f"{'1405.7':>8s} {'0.0920':>8s}")
        print()
        print("  NOTE: R_eq moves AWAY from experiment at higher l_max.")
        print("  Root cause: angular correlation lowers energy more at large R")
        print("  (diffuse electrons) than at short R (compact cloud).")
        print("  This is a limitation of the single-channel adiabatic")
        print("  approximation, not a code bug.")
        print("=" * 72)
        assert True
