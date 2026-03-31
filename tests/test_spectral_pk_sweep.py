"""
Spectral rank-k PK vs Gaussian PK: l_max divergence comparison.

This test module runs the full l_max sweep (0, 1, 2, 3) for both
spectral rank-k PK and Gaussian ab initio PK, measuring R_eq drift.

All tests are marked slow since each l_max requires a PES scan.
"""

import numpy as np
import pytest
import time

from geovac.composed_diatomic import ComposedDiatomicSolver


LIH_R_EQ_EXPT = 3.015  # bohr


def _r_grid() -> np.ndarray:
    """Compact R-grid for PES scans."""
    return np.concatenate([
        np.linspace(2.0, 2.5, 2),
        np.linspace(2.7, 4.0, 6),
        np.linspace(4.5, 7.0, 4),
    ])


def _run_lih(pk_mode: str, l_max: int, **kwargs) -> float:
    """Run LiH and return R_eq."""
    if pk_mode.startswith('spectral_rank_'):
        s = ComposedDiatomicSolver.LiH_spectral_pk(
            l_max=l_max,
            rank=int(pk_mode.split('_')[-1]),
            verbose=False,
            pk_channel_mode='l_dependent',
            **kwargs,
        )
    elif pk_mode == 'ab_initio':
        s = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=l_max, verbose=False,
            pk_channel_mode='l_dependent',
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown pk_mode: {pk_mode}")
    s.solve_core()
    s.scan_pes(R_grid=_r_grid(), n_Re=200)
    return s.pes_result['R_eq']


# ======================================================================
# Rank sensitivity at l_max=0 (validates rank-1 reproduces known failure)
# ======================================================================

@pytest.mark.slow
class TestRankSensitivity:
    """Test rank sensitivity at l_max=0."""

    @pytest.fixture(scope="class")
    def rank_results(self):
        """Run rank 1, 2, 3 at l_max=0."""
        results = {}
        for rank in [1, 2, 3]:
            R_eq = _run_lih(f'spectral_rank_{rank}', l_max=0)
            results[rank] = R_eq
        return results

    def test_all_ranks_produce_minimum(self, rank_results):
        """All ranks should produce a PES with a bound state."""
        for rank, R_eq in rank_results.items():
            assert 1.5 < R_eq < 6.0, \
                f"rank={rank}: R_eq={R_eq:.3f} out of range"

    def test_rank_sensitivity_report(self, rank_results):
        """Report rank sensitivity (informational)."""
        print("\n--- Rank sensitivity at l_max=0 ---")
        for rank, R_eq in sorted(rank_results.items()):
            err = abs(R_eq - LIH_R_EQ_EXPT) / LIH_R_EQ_EXPT * 100
            print(f"  rank={rank}: R_eq={R_eq:.3f} bohr  err={err:.1f}%")


# ======================================================================
# l_max divergence sweep
# ======================================================================

@pytest.mark.slow
class TestLmaxDivergence:
    """Full l_max sweep comparing spectral rank-3 vs Gaussian PK."""

    @pytest.fixture(scope="class")
    def spectral_sweep(self):
        """Run spectral rank-3 at l_max=0,1,2,3."""
        results = {}
        for l_max in [0, 1, 2, 3]:
            t0 = time.time()
            R_eq = _run_lih('spectral_rank_3', l_max=l_max)
            dt = time.time() - t0
            results[l_max] = {'R_eq': R_eq, 'time': dt}
        return results

    @pytest.fixture(scope="class")
    def gaussian_sweep(self):
        """Run Gaussian ab initio at l_max=0,1,2,3."""
        results = {}
        for l_max in [0, 1, 2, 3]:
            t0 = time.time()
            R_eq = _run_lih('ab_initio', l_max=l_max)
            dt = time.time() - t0
            results[l_max] = {'R_eq': R_eq, 'time': dt}
        return results

    def test_spectral_pes_has_minimum_at_all_lmax(self, spectral_sweep):
        """Spectral PES should have a minimum at all l_max."""
        for l_max, data in spectral_sweep.items():
            assert 1.5 < data['R_eq'] < 8.0, \
                f"l_max={l_max}: R_eq={data['R_eq']:.3f} out of range"

    def test_gaussian_pes_has_minimum_at_all_lmax(self, gaussian_sweep):
        """Gaussian PES should have a minimum at all l_max."""
        for l_max, data in gaussian_sweep.items():
            assert 1.5 < data['R_eq'] < 8.0, \
                f"l_max={l_max}: R_eq={data['R_eq']:.3f} out of range"

    def test_divergence_comparison_report(
        self, spectral_sweep, gaussian_sweep,
    ):
        """Report the divergence comparison (primary result)."""
        print("\n" + "=" * 72)
        print("LiH l_max Divergence: Spectral rank-3 vs Gaussian ab initio")
        print("=" * 72)
        print(f"{'l_max':>5s}  {'Spec R_eq':>10s}  {'Gauss R_eq':>10s}"
              f"  {'Spec err%':>10s}  {'Gauss err%':>10s}"
              f"  {'Spec t':>7s}  {'Gauss t':>7s}")
        for l_max in [0, 1, 2, 3]:
            sr = spectral_sweep[l_max]['R_eq']
            gr = gaussian_sweep[l_max]['R_eq']
            se = abs(sr - LIH_R_EQ_EXPT) / LIH_R_EQ_EXPT * 100
            ge = abs(gr - LIH_R_EQ_EXPT) / LIH_R_EQ_EXPT * 100
            st = spectral_sweep[l_max]['time']
            gt = gaussian_sweep[l_max]['time']
            print(f"{l_max:5d}  {sr:10.3f}  {gr:10.3f}"
                  f"  {se:10.1f}  {ge:10.1f}"
                  f"  {st:7.1f}s  {gt:7.1f}s")

        spec_drift = (
            (spectral_sweep[3]['R_eq'] - spectral_sweep[0]['R_eq']) / 3.0
        )
        gauss_drift = (
            (gaussian_sweep[3]['R_eq'] - gaussian_sweep[0]['R_eq']) / 3.0
        )
        print(f"\nDrift rate (bohr/l_max):")
        print(f"  Spectral rank-3:    {spec_drift:+.3f}")
        print(f"  Gaussian ab initio: {gauss_drift:+.3f}")
        if abs(gauss_drift) > 1e-6:
            ratio = spec_drift / gauss_drift
            print(f"  Ratio (spectral/gaussian): {ratio:.2f}x")

        # Store for assertions
        self._spec_drift = spec_drift
        self._gauss_drift = gauss_drift

    def test_spectral_l2_competitive(self, spectral_sweep):
        """Spectral rank-3 at l_max=2 should have R_eq error < 15%."""
        R_eq = spectral_sweep[2]['R_eq']
        err = abs(R_eq - LIH_R_EQ_EXPT) / LIH_R_EQ_EXPT * 100
        assert err < 15.0, \
            f"l_max=2 R_eq error {err:.1f}% exceeds 15% threshold"
