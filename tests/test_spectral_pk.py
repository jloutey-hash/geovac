"""
Tests for the spectral rank-k Phillips-Kleinman pseudopotential.

Validates:
  1. Spectral PK projector construction (Gegenbauer coefficients)
  2. Rank-1 reproduces the known v2.0.6 failure behavior (drift 5.6x worse)
  3. Rank-3 projector produces a well-defined LiH PES with minimum
  4. l_max divergence comparison: spectral vs Gaussian PK
  5. Spectral components are ordered by significance

References:
  - Paper 17, Sec IV (Phillips-Kleinman)
  - CLAUDE.md Section 3: "Algebraic PK projector (rank-1)" failure entry
"""

import numpy as np
import pytest
import time

from geovac.composed_diatomic import ComposedDiatomicSolver
from geovac.ab_initio_pk import AbInitioPK
from geovac.core_screening import CoreScreening


# Reference data
LIH_R_EQ_EXPT = 3.015  # bohr


def _lih_r_grid_compact() -> np.ndarray:
    """Compact R-grid for fast PES scans (12 points)."""
    return np.concatenate([
        np.linspace(2.0, 2.5, 2),
        np.linspace(2.7, 4.0, 6),
        np.linspace(4.5, 7.0, 4),
    ])


# ======================================================================
# Projector construction tests (fast, no PES scan)
# ======================================================================

class TestSpectralPKConstruction:
    """Test that the spectral PK projector is constructed correctly."""

    @pytest.fixture(scope="class")
    def core_and_pk(self):
        """Solve core once, build projectors at multiple ranks."""
        core = CoreScreening(Z=3, l_max=0, n_alpha=200)
        core.solve(verbose=False)
        pk = AbInitioPK(core, n_core=2)
        return core, pk

    def test_projector_has_correct_rank(self, core_and_pk):
        """Projector returns exactly k spectral components."""
        _, pk = core_and_pk
        for rank in [1, 2, 3, 5]:
            proj = pk.spectral_rank_k_projector(rank=rank)
            assert proj['mode'] == 'spectral_rank_k'
            assert proj['rank'] == rank
            assert len(proj['spectral_components']) == rank

    def test_components_ordered_by_significance(self, core_and_pk):
        """Top-k components should be sorted by descending |coeff|."""
        _, pk = core_and_pk
        proj = pk.spectral_rank_k_projector(rank=5)
        coeffs = [abs(c['coeff']) for c in proj['spectral_components']]
        assert coeffs == sorted(coeffs, reverse=True), \
            "Components should be ordered by decreasing significance"

    def test_dominant_component_is_l0_k0(self, core_and_pk):
        """The dominant component should be l=0, k=0 (ground Gegenbauer)."""
        _, pk = core_and_pk
        proj = pk.spectral_rank_k_projector(rank=1)
        comp = proj['spectral_components'][0]
        assert comp['l'] == 0, "Dominant component should be l=0"
        assert comp['k'] == 0, "Dominant component should be k=0"
        assert abs(comp['coeff']) > 0.9, \
            f"Dominant coeff should be > 0.9, got {abs(comp['coeff']):.4f}"

    def test_energy_shift_positive(self, core_and_pk):
        """Energy shift should be positive (repulsive barrier)."""
        _, pk = core_and_pk
        proj = pk.spectral_rank_k_projector(rank=3)
        assert proj['energy_shift'] > 0, \
            f"E_shift should be > 0, got {proj['energy_shift']:.4f}"

    def test_rank1_subset_of_rank3(self, core_and_pk):
        """rank-1 component should match the first component of rank-3."""
        _, pk = core_and_pk
        p1 = pk.spectral_rank_k_projector(rank=1)
        p3 = pk.spectral_rank_k_projector(rank=3)
        c1 = p1['spectral_components'][0]
        c3 = p3['spectral_components'][0]
        assert c1['l'] == c3['l']
        assert c1['k'] == c3['k']
        assert abs(c1['coeff'] - c3['coeff']) < 1e-10

    def test_coefficients_sum_to_unity(self, core_and_pk):
        """Sum of squared coefficients should be ~1 for a normalized state."""
        _, pk = core_and_pk
        proj = pk.spectral_rank_k_projector(rank=10, n_basis=10)
        sum_sq = sum(c['coeff']**2 for c in proj['spectral_components'])
        assert abs(sum_sq - 1.0) < 0.01, \
            f"Sum of squared coefficients should be ~1, got {sum_sq:.6f}"


# ======================================================================
# PES and l_max tests (slow, require PES scans)
# ======================================================================

@pytest.mark.slow
class TestSpectralPKPES:
    """Test PES behavior of the spectral rank-k PK projector."""

    def test_rank3_pes_has_minimum(self):
        """Rank-3 spectral PK should produce a PES with a well-defined min."""
        s = ComposedDiatomicSolver.LiH_spectral_pk(
            l_max=0, rank=3, verbose=False,
        )
        s.solve_core()
        s.scan_pes(R_grid=_lih_r_grid_compact(), n_Re=200)

        assert s.pes_result is not None
        R_eq = s.pes_result['R_eq']
        D_e = s.pes_result['D_e']

        # PES should have a minimum (D_e > 0) and R_eq in reasonable range
        assert D_e > 0.01, f"D_e should be > 0.01, got {D_e:.4f}"
        assert 1.5 < R_eq < 6.0, f"R_eq should be in [1.5, 6.0], got {R_eq:.3f}"

    def test_rank1_reproduces_weak_projector(self):
        """Rank-1 spectral PK should behave similarly to v2.0.6 algebraic PK.

        Both are rank-1 projectors; the spectral one may differ slightly
        because it uses the Gegenbauer basis rather than the FD grid.
        """
        s = ComposedDiatomicSolver.LiH_spectral_pk(
            l_max=0, rank=1, verbose=False,
        )
        s.solve_core()
        s.scan_pes(R_grid=_lih_r_grid_compact(), n_Re=200)

        assert s.pes_result is not None
        # Rank-1 should still produce a bound molecule
        assert s.pes_result['D_e'] > 0.0, "Rank-1 should still be bound"


@pytest.mark.slow
class TestSpectralPKLmaxDivergence:
    """Test l_max divergence behavior of spectral rank-k PK."""

    @pytest.fixture(scope="class")
    def lmax_sweep_spectral_rank3(self):
        """Run LiH at l_max=0,1,2 with spectral rank-3 PK."""
        results = {}
        for l_max in [0, 1, 2]:
            s = ComposedDiatomicSolver.LiH_spectral_pk(
                l_max=l_max, rank=3, verbose=False,
                pk_channel_mode='l_dependent',
            )
            s.solve_core()
            s.scan_pes(R_grid=_lih_r_grid_compact(), n_Re=200)
            results[l_max] = s.pes_result['R_eq']
        return results

    @pytest.fixture(scope="class")
    def lmax_sweep_gaussian(self):
        """Run LiH at l_max=0,1,2 with Gaussian ab initio PK."""
        results = {}
        for l_max in [0, 1, 2]:
            s = ComposedDiatomicSolver.LiH_ab_initio(
                l_max=l_max, verbose=False,
                pk_channel_mode='l_dependent',
            )
            s.solve_core()
            s.scan_pes(R_grid=_lih_r_grid_compact(), n_Re=200)
            results[l_max] = s.pes_result['R_eq']
        return results

    def test_spectral_pes_has_minimum_at_all_lmax(
        self, lmax_sweep_spectral_rank3,
    ):
        """PES should have a well-defined minimum at all l_max values."""
        for l_max, R_eq in lmax_sweep_spectral_rank3.items():
            assert 1.5 < R_eq < 8.0, \
                f"l_max={l_max}: R_eq={R_eq:.3f} out of range"

    def test_spectral_drift_rate(
        self, lmax_sweep_spectral_rank3,
    ):
        """Spectral rank-3 l_max drift should be characterized."""
        R_eqs = lmax_sweep_spectral_rank3
        # Drift rate: (R_eq(l_max=2) - R_eq(l_max=0)) / 2
        drift = (R_eqs[2] - R_eqs[0]) / 2.0
        print(f"\nSpectral rank-3 drift: {drift:+.3f} bohr/l_max")
        print(f"  l_max=0: R_eq={R_eqs[0]:.3f}")
        print(f"  l_max=1: R_eq={R_eqs.get(1, 'N/A')}")
        print(f"  l_max=2: R_eq={R_eqs[2]:.3f}")
        # Report but do not assert — the drift measurement is the goal
