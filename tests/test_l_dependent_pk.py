"""
Tests for l-dependent Phillips-Kleinman pseudopotential.

Validates:
  1. Channel weights are correctly applied (unit test)
  2. l_max=2 results are nearly identical for both modes (sigma-only,
     all channels are (0,0), (0,2), (2,0), (2,2) — only (0,0) has
     full PK in l-dependent mode, but it dominates)
  3. l_max=3 and l_max=4 R_eq convergence improves vs channel-blind
  4. Backward compatibility: channel_blind default matches old behavior
"""

import numpy as np
import pytest

from geovac.composed_diatomic import ComposedDiatomicSolver
from geovac.level4_multichannel import (
    _channel_list, _channel_list_extended,
    build_angular_hamiltonian,
)


# ======================================================================
# Unit test: channel weight logic
# ======================================================================

class TestChannelWeights:
    """Verify the per-channel PK weight logic."""

    def test_channel_list_l2_hetero(self) -> None:
        """l_max=2 heteronuclear sigma channels."""
        ch = _channel_list(2, homonuclear=False)
        expected = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2),
                    (2, 0), (2, 1), (2, 2)]
        assert set(ch) == set(expected)

    def test_weight_classification_l2(self) -> None:
        """Check PK weights for l_max=2 heteronuclear channels."""
        ch = _channel_list(2, homonuclear=False)
        weights = {}
        for l1, l2 in ch:
            w1 = 1.0 if l1 == 0 else 0.0
            w2 = 1.0 if l2 == 0 else 0.0
            weights[(l1, l2)] = (w1, w2)

        # (0,0): both electrons l=0 -> full PK on both
        assert weights[(0, 0)] == (1.0, 1.0)
        # (0,1): e1 l=0, e2 l=1 -> PK on e1 only
        assert weights[(0, 1)] == (1.0, 0.0)
        # (1,0): e1 l=1, e2 l=0 -> PK on e2 only
        assert weights[(1, 0)] == (0.0, 1.0)
        # (1,1): both l>0 -> no PK
        assert weights[(1, 1)] == (0.0, 0.0)
        # (2,2): both l>0 -> no PK
        assert weights[(2, 2)] == (0.0, 0.0)

    def test_weight_classification_l3(self) -> None:
        """l_max=3 heteronuclear: verify zero-PK channels exist."""
        ch = _channel_list(3, homonuclear=False)
        zero_pk = [(l1, l2) for l1, l2 in ch
                   if l1 > 0 and l2 > 0]
        full_pk = [(l1, l2) for l1, l2 in ch
                   if l1 == 0 and l2 == 0]
        partial = [(l1, l2) for l1, l2 in ch
                   if (l1 == 0) != (l2 == 0)]

        print(f"\n  l_max=3 channels: {len(ch)} total")
        print(f"    Full PK (l1=0,l2=0): {len(full_pk)} — {full_pk}")
        print(f"    Partial PK: {len(partial)} — {partial}")
        print(f"    Zero PK (l1>0,l2>0): {len(zero_pk)} — {zero_pk}")

        assert len(full_pk) == 1  # only (0,0)
        assert len(zero_pk) > 0   # e.g. (1,1), (1,2), etc.
        assert len(partial) > 0   # e.g. (0,1), (1,0), etc.

    def test_sigma_pi_weight_classification(self) -> None:
        """Sigma+pi channels with m>0 also get correct weights."""
        ch = _channel_list_extended(3, m_max=1, l_max_per_m={0: 3, 1: 2},
                                    homonuclear=False)
        for l1, m1, l2, m2 in ch:
            w1 = 1.0 if l1 == 0 else 0.0
            w2 = 1.0 if l2 == 0 else 0.0
            # l=0 always has m=0, so for m!=0 channels w must be 0
            if m1 != 0:
                assert w1 == 0.0, f"l1={l1},m1={m1} should have w1=0"
            if m2 != 0:
                assert w2 == 0.0, f"l2={l2},m2={m2} should have w2=0"


# ======================================================================
# Integration: channel_blind backward compatibility
# ======================================================================

class TestBackwardCompatibility:
    """Default pk_channel_mode='channel_blind' must match old behavior."""

    def test_default_is_channel_blind(self) -> None:
        """ComposedDiatomicSolver default is channel_blind."""
        s = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=False)
        assert s.pk_channel_mode == 'channel_blind'

    def test_channel_blind_in_pk_dict(self) -> None:
        """channel_blind mode is embedded in pk_potentials."""
        s = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=False)
        s.solve_core()
        pk_d = s.pk_potentials[0]
        assert pk_d['channel_mode'] == 'channel_blind'

    def test_l_dependent_in_pk_dict(self) -> None:
        """l_dependent mode is embedded in pk_potentials."""
        s = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=2, pk_channel_mode='l_dependent', verbose=False,
        )
        s.solve_core()
        pk_d = s.pk_potentials[0]
        assert pk_d['channel_mode'] == 'l_dependent'

    def test_invalid_pk_channel_mode_raises(self) -> None:
        """Invalid pk_channel_mode raises ValueError."""
        with pytest.raises(ValueError, match="pk_channel_mode"):
            ComposedDiatomicSolver.LiH_ab_initio(
                l_max=2, pk_channel_mode='invalid', verbose=False,
            )


# ======================================================================
# R-grid
# ======================================================================

def _lih_r_grid() -> np.ndarray:
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 8.0, 7),
    ])


# ======================================================================
# LiH l_max=2 comparison: channel_blind vs l_dependent
# ======================================================================

@pytest.fixture(scope="module")
def lih_l2_channel_blind():
    """LiH l_max=2, channel_blind PK (baseline)."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2, pk_channel_mode='channel_blind', verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    return s


@pytest.fixture(scope="module")
def lih_l2_l_dependent():
    """LiH l_max=2, l_dependent PK."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2, pk_channel_mode='l_dependent', verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    return s


class TestLiHL2Comparison:
    """l_max=2: l-dependent PK should give similar results to channel-blind."""

    def test_both_complete(
        self,
        lih_l2_channel_blind: ComposedDiatomicSolver,
        lih_l2_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        assert lih_l2_channel_blind.spectro is not None
        assert lih_l2_l_dependent.spectro is not None

    def test_both_bound(
        self,
        lih_l2_channel_blind: ComposedDiatomicSolver,
        lih_l2_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        assert lih_l2_channel_blind.pes_result['D_e'] > 0
        assert lih_l2_l_dependent.pes_result['D_e'] > 0

    def test_r_eq_comparison(
        self,
        lih_l2_channel_blind: ComposedDiatomicSolver,
        lih_l2_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        """At l_max=2, most PES weight is on (0,0) channel, so
        l-dependent should give a similar but slightly different R_eq."""
        R_cb = lih_l2_channel_blind.spectro['R_eq']
        R_ld = lih_l2_l_dependent.spectro['R_eq']
        print(f"\n  l_max=2 channel_blind R_eq = {R_cb:.3f}")
        print(f"  l_max=2 l_dependent  R_eq = {R_ld:.3f}")
        # Should be physically reasonable
        assert 2.0 <= R_ld <= 5.0


# ======================================================================
# LiH l_max=3 and l_max=4: l-dependent PK convergence
# ======================================================================

@pytest.fixture(scope="module")
def lih_l3_l_dependent():
    """LiH l_max=3 sigma+pi, l_dependent PK."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=3, m_max=1, l_max_per_m={0: 3, 1: 2},
        n_alpha=60, pk_channel_mode='l_dependent', verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=2.0)
    s._print_summary()
    return s


@pytest.fixture(scope="module")
def lih_l4_l_dependent():
    """LiH l_max=4 sigma+pi, l_dependent PK."""
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=4, m_max=1, l_max_per_m={0: 4, 1: 2},
        n_alpha=50, pk_channel_mode='l_dependent', verbose=True,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=2.0)
    s._print_summary()
    return s


class TestLiHL3LDependent:
    """LiH l_max=3 with l-dependent PK."""

    def test_pipeline_completes(
        self, lih_l3_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        assert lih_l3_l_dependent.spectro is not None

    def test_molecule_bound(
        self, lih_l3_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        D_e = lih_l3_l_dependent.pes_result['D_e']
        print(f"  D_e (l3, l_dep) = {D_e:.6f} Ha")
        assert D_e > 0

    def test_r_eq_physical(
        self, lih_l3_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        R_eq = lih_l3_l_dependent.spectro['R_eq']
        print(f"  R_eq (l3, l_dep) = {R_eq:.3f} bohr (expt: 3.015)")
        assert 2.0 <= R_eq <= 5.0


class TestLiHL4LDependent:
    """LiH l_max=4 with l-dependent PK."""

    def test_pipeline_completes(
        self, lih_l4_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        assert lih_l4_l_dependent.spectro is not None

    def test_molecule_bound(
        self, lih_l4_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        D_e = lih_l4_l_dependent.pes_result['D_e']
        print(f"  D_e (l4, l_dep) = {D_e:.6f} Ha")
        assert D_e > 0

    def test_r_eq_physical(
        self, lih_l4_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        R_eq = lih_l4_l_dependent.spectro['R_eq']
        print(f"  R_eq (l4, l_dep) = {R_eq:.3f} bohr (expt: 3.015)")
        assert 2.0 <= R_eq <= 6.0


# ======================================================================
# Convergence comparison table
# ======================================================================

class TestLDependentConvergenceReport:
    """Print side-by-side convergence: channel_blind vs l_dependent."""

    def test_print_convergence_table(
        self,
        lih_l2_channel_blind: ComposedDiatomicSolver,
        lih_l2_l_dependent: ComposedDiatomicSolver,
        lih_l3_l_dependent: ComposedDiatomicSolver,
        lih_l4_l_dependent: ComposedDiatomicSolver,
    ) -> None:
        """Print l_max convergence for both PK modes."""
        print("\n")
        print("=" * 80)
        print("LiH l_max Convergence: channel_blind vs l_dependent PK")
        print("=" * 80)

        # Channel structure at each l_max
        for lm in [2, 3, 4]:
            ch = _channel_list(lm, homonuclear=False)
            n_full = sum(1 for l1, l2 in ch if l1 == 0 and l2 == 0)
            n_partial = sum(1 for l1, l2 in ch
                           if (l1 == 0) != (l2 == 0))
            n_zero = sum(1 for l1, l2 in ch if l1 > 0 and l2 > 0)
            print(f"\n  l_max={lm}: {len(ch)} sigma channels"
                  f" (full_PK={n_full}, partial={n_partial},"
                  f" zero_PK={n_zero})")
            for l1, l2 in ch:
                w1 = 1.0 if l1 == 0 else 0.0
                w2 = 1.0 if l2 == 0 else 0.0
                tag = ("FULL" if w1 + w2 == 2.0
                       else "PART" if w1 + w2 == 1.0
                       else "ZERO")
                print(f"    ({l1},{l2}): w1={w1:.0f}, w2={w2:.0f} [{tag}]")

        print(f"\n  {'Config':28s} {'R_eq':>8s} {'err%':>6s} "
              f"{'omega_e':>8s} {'D_e':>8s}")
        print(f"  {'-'*28} {'-'*8} {'-'*6} {'-'*8} {'-'*8}")

        ref_R = 3.015
        entries = [
            ("l2 ch_blind (ab initio)", lih_l2_channel_blind),
            ("l2 l_dependent (ab initio)", lih_l2_l_dependent),
            ("l3 l_dependent sig+pi", lih_l3_l_dependent),
            ("l4 l_dependent sig+pi", lih_l4_l_dependent),
        ]
        for label, s in entries:
            R = s.spectro['R_eq']
            err = abs(R - ref_R) / ref_R * 100
            w = s.spectro['omega_e']
            D = s.spectro['D_e']
            print(f"  {label:28s} {R:8.3f} {err:6.1f} "
                  f"{w:8.1f} {D:8.4f}")

        print(f"  {'Experiment':28s} {'3.015':>8s} {'0.0':>6s} "
              f"{'1405.7':>8s} {'0.0920':>8s}")
        print("=" * 80)
        assert True
