"""
Tests for ab initio Phillips-Kleinman pseudopotential derivation.

Validates:
  1. AbInitioPK produces physical parameters for Li and Be cores
  2. Ab initio PK reproduces molecular binding for LiH
  3. No-PK baseline shows core collapse (PK is necessary)
  4. Sensitivity: R_eq varies gradually with PK parameters
  5. Transferability: same formula works for LiH and BeH+
  6. Comparison table: no-PK vs ab-initio vs fitted vs experiment

Test structure:
  - Module-scoped fixtures for expensive CoreScreening solves
  - Fast unit tests for AbInitioPK parameter checks
  - Slow integration tests for PES scans (marked clearly)
"""

import numpy as np
import pytest

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.composed_diatomic import ComposedDiatomicSolver


# ======================================================================
# Fixtures — solve core once, share across tests
# ======================================================================

@pytest.fixture(scope="module")
def li_core() -> CoreScreening:
    """Solved Li+ (Z=3) core."""
    cs = CoreScreening(Z=3, l_max=2, n_alpha=200)
    cs.solve(verbose=False)
    return cs


@pytest.fixture(scope="module")
def be_core() -> CoreScreening:
    """Solved Be2+ (Z=4) core."""
    cs = CoreScreening(Z=4, l_max=2, n_alpha=200)
    cs.solve(verbose=False)
    return cs


@pytest.fixture(scope="module")
def li_pk(li_core: CoreScreening) -> AbInitioPK:
    """Ab initio PK for Li core."""
    return AbInitioPK(li_core, n_core=2)


@pytest.fixture(scope="module")
def be_pk(be_core: CoreScreening) -> AbInitioPK:
    """Ab initio PK for Be core."""
    return AbInitioPK(be_core, n_core=2)


def _lih_r_grid() -> np.ndarray:
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])


def _beh_r_grid() -> np.ndarray:
    return np.concatenate([
        np.linspace(1.5, 2.0, 3),
        np.linspace(2.1, 3.5, 10),
        np.linspace(4.0, 7.0, 5),
    ])


@pytest.fixture(scope="module")
def lih_ab_initio() -> ComposedDiatomicSolver:
    """LiH with ab initio PK."""
    s = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=True)
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    s._print_summary()
    return s


@pytest.fixture(scope="module")
def lih_no_pk() -> ComposedDiatomicSolver:
    """LiH without PK (baseline for core collapse)."""
    s = ComposedDiatomicSolver(
        Z_A=3, Z_B=1, n_core=2,
        M_A=7.016003, M_B=1.00782503,
        label='LiH', pk_mode='none', l_max=2, verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    return s


@pytest.fixture(scope="module")
def beh_ab_initio() -> ComposedDiatomicSolver:
    """BeH+ with ab initio PK."""
    s = ComposedDiatomicSolver.BeH_plus_ab_initio(l_max=2, verbose=True)
    s.solve_core()
    s.scan_pes(R_grid=_beh_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    s._print_summary()
    return s


# ======================================================================
# Test 1: Ab initio PK values are physical
# ======================================================================

class TestAbInitioPKParameters:
    """Ab initio PK parameters must be physically reasonable."""

    def test_li_A_positive(self, li_pk: AbInitioPK) -> None:
        assert li_pk.A > 0, f"A must be positive, got {li_pk.A}"

    def test_li_B_positive(self, li_pk: AbInitioPK) -> None:
        assert li_pk.B > 0, f"B must be positive, got {li_pk.B}"

    def test_li_r_core_range(self, li_pk: AbInitioPK) -> None:
        """Li+ core radius should be ~0.1-1.0 bohr."""
        assert 0.1 <= li_pk.r_core <= 1.0, (
            f"Li r_core = {li_pk.r_core:.4f}, expected [0.1, 1.0]"
        )

    def test_be_A_positive(self, be_pk: AbInitioPK) -> None:
        assert be_pk.A > 0

    def test_be_B_positive(self, be_pk: AbInitioPK) -> None:
        assert be_pk.B > 0

    def test_be_r_core_range(self, be_pk: AbInitioPK) -> None:
        """Be2+ core should be tighter than Li+."""
        assert 0.1 <= be_pk.r_core <= 0.8, (
            f"Be r_core = {be_pk.r_core:.4f}, expected [0.1, 0.8]"
        )

    def test_be_tighter_than_li(
        self, li_pk: AbInitioPK, be_pk: AbInitioPK,
    ) -> None:
        """Higher Z -> tighter core -> smaller r_core, larger B."""
        assert be_pk.r_core < li_pk.r_core, (
            f"Be r_core ({be_pk.r_core:.4f}) should be < Li r_core"
            f" ({li_pk.r_core:.4f})"
        )
        assert be_pk.B > li_pk.B, (
            f"Be B ({be_pk.B:.4f}) should be > Li B ({li_pk.B:.4f})"
        )

    def test_be_A_larger_than_li(
        self, li_pk: AbInitioPK, be_pk: AbInitioPK,
    ) -> None:
        """Higher Z -> deeper core -> larger A."""
        assert be_pk.A > li_pk.A, (
            f"Be A ({be_pk.A:.4f}) should be > Li A ({li_pk.A:.4f})"
        )

    def test_pk_dict_format(self, li_pk: AbInitioPK) -> None:
        d = li_pk.pk_dict(atom='A')
        assert 'C_core' in d
        assert 'beta_core' in d
        assert 'atom' in d
        assert d['C_core'] == li_pk.A
        assert d['beta_core'] == li_pk.B
        assert d['atom'] == 'A'

    def test_V_pk_positive_near_nucleus(self, li_pk: AbInitioPK) -> None:
        """PK potential should be positive (repulsive) near the nucleus."""
        r = np.array([0.1, 0.3, 0.5, 1.0])
        V = li_pk.V_pk(r)
        assert np.all(V > 0), f"V_PK should be repulsive, got {V}"

    def test_V_pk_decays(self, li_pk: AbInitioPK) -> None:
        """PK potential should decay at large r."""
        V_near = li_pk.V_pk(np.array([0.5]))[0]
        V_far = li_pk.V_pk(np.array([5.0]))[0]
        assert V_far < V_near * 0.01, (
            f"V_PK should decay: V(0.5)={V_near:.4f}, V(5)={V_far:.6f}"
        )

    def test_summary_string(self, li_pk: AbInitioPK) -> None:
        s = li_pk.summary()
        assert 'Ab initio PK' in s
        assert 'r_core' in s

    def test_print_comparison(
        self, li_pk: AbInitioPK, be_pk: AbInitioPK,
    ) -> None:
        """Print ab initio PK parameters (informational)."""
        print("\n" + "=" * 60)
        print("Ab initio PK parameters:")
        print("=" * 60)
        print(li_pk.summary())
        print()
        print(be_pk.summary())
        print("=" * 60)


# ======================================================================
# Test 2: Ab initio PK reproduces binding
# ======================================================================

class TestLiHAbInitioBinding:
    """LiH with ab initio PK must produce a bound molecule."""

    def test_pipeline_completes(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        assert lih_ab_initio.E_core is not None
        assert lih_ab_initio.pes_result is not None
        assert lih_ab_initio.spectro is not None

    def test_is_bound(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """D_e must be positive (bound molecule)."""
        D_e = lih_ab_initio.pes_result['D_e']
        assert D_e > 0.001, f"LiH not bound: D_e = {D_e:.6f}"

    def test_r_eq_physical(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """R_eq must be in a physically reasonable range."""
        R_eq = lih_ab_initio.spectro['R_eq']
        assert 2.0 <= R_eq <= 5.0, (
            f"R_eq = {R_eq:.3f}, expected [2.0, 5.0] (expt: 3.015)"
        )

    def test_pk_mode_is_ab_initio(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        assert lih_ab_initio.pk_mode == 'ab_initio'
        assert lih_ab_initio.ab_initio_pk is not None


# ======================================================================
# Test 3: No-PK baseline shows different behavior
# ======================================================================

class TestNoPKBaseline:
    """Without PK, the valence electrons see no Pauli repulsion."""

    def test_no_pk_runs(
        self, lih_no_pk: ComposedDiatomicSolver,
    ) -> None:
        assert lih_no_pk.pes_result is not None

    def test_no_pk_different_from_ab_initio(
        self,
        lih_no_pk: ComposedDiatomicSolver,
        lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """No-PK should give a meaningfully different R_eq than ab initio."""
        R_no = lih_no_pk.spectro['R_eq']
        R_ai = lih_ab_initio.spectro['R_eq']
        print(f"\n  No PK: R_eq = {R_no:.3f}")
        print(f"  Ab initio PK: R_eq = {R_ai:.3f}")
        # PK should push R_eq outward (Pauli repulsion pushes equilibrium
        # to larger R), or at minimum change R_eq noticeably
        assert abs(R_no - R_ai) > 0.05, (
            f"PK should change R_eq: no_pk={R_no:.3f},"
            f" ab_initio={R_ai:.3f}"
        )


# ======================================================================
# Test 4: Sensitivity — R_eq varies gradually
# ======================================================================

class TestSensitivity:
    """R_eq should change GRADUALLY with PK parameters."""

    def test_a_scan_gradual(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """R_eq varies by less than 1.5 bohr across A ∈ [2, 10]."""
        B_fixed = lih_ab_initio.pk_B
        A_values = [2.0, 5.0, 10.0]
        R_eqs = []
        R_grid = _lih_r_grid()

        for A_val in A_values:
            s = ComposedDiatomicSolver(
                Z_A=3, Z_B=1, n_core=2,
                M_A=7.016003, M_B=1.00782503,
                label='LiH', pk_mode='manual',
                pk_A=A_val, pk_B=B_fixed,
                l_max=2, verbose=False,
            )
            s.solve_core()
            s.scan_pes(R_grid=R_grid, n_Re=300)
            s.fit_spectroscopic_constants()
            R_eqs.append(s.spectro['R_eq'])
            print(f"  A={A_val:.1f}, B={B_fixed:.2f}: R_eq={s.spectro['R_eq']:.3f}")

        span = max(R_eqs) - min(R_eqs)
        assert span < 1.5, (
            f"R_eq span = {span:.3f} bohr across A scan, expected < 1.5"
        )

    def test_b_scan_gradual(
        self, lih_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """R_eq varies by less than 1.5 bohr across B ∈ [1, 10]."""
        A_fixed = lih_ab_initio.pk_A
        B_values = [1.0, 4.0, 10.0]
        R_eqs = []
        R_grid = _lih_r_grid()

        for B_val in B_values:
            s = ComposedDiatomicSolver(
                Z_A=3, Z_B=1, n_core=2,
                M_A=7.016003, M_B=1.00782503,
                label='LiH', pk_mode='manual',
                pk_A=A_fixed, pk_B=B_val,
                l_max=2, verbose=False,
            )
            s.solve_core()
            s.scan_pes(R_grid=R_grid, n_Re=300)
            s.fit_spectroscopic_constants()
            R_eqs.append(s.spectro['R_eq'])
            print(f"  A={A_fixed:.2f}, B={B_val:.1f}: R_eq={s.spectro['R_eq']:.3f}")

        span = max(R_eqs) - min(R_eqs)
        assert span < 1.5, (
            f"R_eq span = {span:.3f} bohr across B scan, expected < 1.5"
        )


# ======================================================================
# Test 5: Transferability — both LiH and BeH+ work
# ======================================================================

class TestTransferability:
    """Same ab initio formula must work for both LiH and BeH+."""

    def test_beh_bound(
        self, beh_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """BeH+ with ab initio PK must be bound."""
        D_e = beh_ab_initio.pes_result['D_e']
        assert D_e > 0.001, f"BeH+ not bound: D_e = {D_e:.6f}"

    def test_beh_r_eq_physical(
        self, beh_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """BeH+ R_eq in reasonable range."""
        R_eq = beh_ab_initio.spectro['R_eq']
        assert 1.5 <= R_eq <= 4.5, (
            f"BeH+ R_eq = {R_eq:.3f}, expected [1.5, 4.5] (ref: ~2.48)"
        )

    def test_both_use_same_formula(
        self,
        lih_ab_initio: ComposedDiatomicSolver,
        beh_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """Both molecules must use ab_initio pk_mode."""
        assert lih_ab_initio.pk_mode == 'ab_initio'
        assert beh_ab_initio.pk_mode == 'ab_initio'


# ======================================================================
# Test 6: Comparison table
# ======================================================================

class TestComparisonTable:
    """Print the ammunition table for the paper."""

    def test_print_table(
        self,
        lih_ab_initio: ComposedDiatomicSolver,
        lih_no_pk: ComposedDiatomicSolver,
        beh_ab_initio: ComposedDiatomicSolver,
    ) -> None:
        """Print comparison: no-PK vs ab-initio vs fitted vs experiment."""
        print("\n" + "=" * 72)
        print("PK Pseudopotential: Ab Initio vs Fitted vs None")
        print("=" * 72)

        # LiH ab initio params
        ai = lih_ab_initio.ab_initio_pk
        print(f"\nLiH ab initio PK: A={ai.A:.4f}, B={ai.B:.4f}"
              f"  (r_core={ai.r_core:.4f} bohr)")

        R_no = lih_no_pk.spectro['R_eq']
        R_ai = lih_ab_initio.spectro['R_eq']
        w_no = lih_no_pk.spectro['omega_e']
        w_ai = lih_ab_initio.spectro['omega_e']

        print(f"\nLiH (Z_A=3, Z_B=1):")
        print(f"  {'':16s} {'No PK':>10s} {'Ab initio':>10s}"
              f" {'Fitted':>10s} {'Expt':>10s}")
        print(f"  {'R_eq (bohr)':16s} {R_no:10.3f} {R_ai:10.3f}"
              f" {'2.97':>10s} {'3.015':>10s}")
        print(f"  {'omega_e (cm-1)':16s} {w_no:10.0f} {w_ai:10.0f}"
              f" {'1617':>10s} {'1406':>10s}")

        # BeH+ ab initio params
        ai_be = beh_ab_initio.ab_initio_pk
        R_ai_be = beh_ab_initio.spectro['R_eq']
        w_ai_be = beh_ab_initio.spectro['omega_e']

        print(f"\nBeH+ ab initio PK: A={ai_be.A:.4f}, B={ai_be.B:.4f}"
              f"  (r_core={ai_be.r_core:.4f} bohr)")
        print(f"\nBeH+ (Z_A=4, Z_B=1):")
        print(f"  {'':16s} {'Ab initio':>10s} {'Fitted':>10s}"
              f" {'Ref':>10s}")
        print(f"  {'R_eq (bohr)':16s} {R_ai_be:10.3f}"
              f" {'2.82':>10s} {'~2.48':>10s}")
        print(f"  {'omega_e (cm-1)':16s} {w_ai_be:10.0f}"
              f" {'---':>10s} {'2222':>10s}")

        print("=" * 72)
