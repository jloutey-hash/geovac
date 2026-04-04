"""
Tests for Track X2 cusp correction applied to 2D variational Level 4 H2 solver.

Validates:
1. 2D solver at l_max=2 sigma-only reproduces Paper 15 Table III (~79% D_e)
2. Cusp correction is negative at all R
3. Cusp correction is R-dependent (larger at short R)
4. Corrected D_e > uncorrected D_e
5. Cusp correction vanishes at dissociation

All tests here run the actual Level 4 2D solver and are marked @pytest.mark.slow.
"""

import numpy as np
import pytest
import time

# All tests in this module are slow (2D solver)
pytestmark = pytest.mark.slow

# Exact D_e for H2 from Kolos & Wolniewicz
D_E_EXACT = 0.17447  # Ha

# R grid for PES scan
R_GRID = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0, 10.0])


def _run_2d_pes_scan(
    R_grid: np.ndarray,
    l_max: int = 2,
    m_max: int = 0,
    n_alpha: int = 60,
    n_Re: int = 120,
) -> dict:
    """Run Level 4 2D solver at multiple R points."""
    from geovac.level4_multichannel import solve_level4_h2_multichannel

    E_totals = []
    per_point = []
    t0 = time.time()

    for R in R_grid:
        t_start = time.time()
        result = solve_level4_h2_multichannel(
            R=float(R),
            l_max=l_max,
            m_max=m_max,
            n_coupled=-1,
            n_alpha=n_alpha,
            n_Re=n_Re,
            verbose=False,
        )
        t_end = time.time()

        E_totals.append(result['E_total'])
        per_point.append({
            'R': float(R),
            'E_total': float(result['E_total']),
            'D_e': float(result['D_e']),
            'D_e_pct': float(result['D_e_pct']),
            'l_max': l_max,
            'wall_time': t_end - t_start,
        })

    total_time = time.time() - t0

    return {
        'R_grid': R_grid,
        'E_totals': np.array(E_totals),
        'per_point': per_point,
        'l_max': l_max,
        'total_time': total_time,
    }


def _apply_cusp(scan: dict) -> dict:
    """Apply cusp correction to PES scan results."""
    from geovac.cusp_correction import (
        cusp_correction_h2_pes,
        cusp_correction_h2_point,
        coalescence_density_h2,
    )

    R_grid = scan['R_grid']
    E_pes = scan['E_totals']
    l_max = scan['l_max']

    pes_result = cusp_correction_h2_pes(R_grid, E_pes, l_max=l_max, E_atoms=-1.0)

    per_point_cusp = []
    for pp in scan['per_point']:
        dE = cusp_correction_h2_point(pp['R'], l_max=l_max)
        delta = coalescence_density_h2(pp['R'])
        per_point_cusp.append({
            'R': pp['R'],
            'E_uncorrected': pp['E_total'],
            'dE_cusp': float(dE),
            'E_corrected': pp['E_total'] + dE,
            'delta_r12': float(delta),
        })

    pes_result['per_point_cusp'] = per_point_cusp
    return pes_result


# ============================================================================
# Shared fixture: run the 2D PES scan once, reuse across tests
# ============================================================================

@pytest.fixture(scope='module')
def pes_scan_2d_lmax2():
    """Run the Level 4 2D PES scan at l_max=2 sigma-only once for all tests."""
    scan = _run_2d_pes_scan(R_GRID, l_max=2, m_max=0, n_alpha=60, n_Re=120)
    pes_result = _apply_cusp(scan)
    return {
        'scan': scan,
        'pes_result': pes_result,
    }


# ============================================================================
# Tests
# ============================================================================


class TestCuspCorrection2D:
    """Validate Track X2 cusp correction on 2D variational Level 4 H2 PES."""

    def test_2d_baseline_consistency(self, pes_scan_2d_lmax2):
        """2D solver at l_max=2 sigma-only should give D_e ~ 79.5% +/- 2%
        (Paper 15 Table III)."""
        pes_result = pes_scan_2d_lmax2['pes_result']
        D_e_pct = pes_result['D_e_pct_uncorrected']

        print(f"\n  2D l_max=2 sigma-only D_e: {D_e_pct:.1f}% "
              f"(Paper 15 Table III: ~79.5%)")

        # Paper 15 Table III: l_max=2, m_max=0 gives ~79.5%
        # Allow +/-2% for grid resolution differences
        assert 77.0 < D_e_pct < 82.0, (
            f"D_e = {D_e_pct:.1f}%, expected ~79.5% +/- 2%"
        )

    def test_2d_cusp_correction_negative(self, pes_scan_2d_lmax2):
        """Cusp correction must be negative at all R."""
        pes_result = pes_scan_2d_lmax2['pes_result']

        for pp in pes_result['per_point_cusp']:
            assert pp['dE_cusp'] < 0, (
                f"Cusp correction at R={pp['R']} should be negative, "
                f"got {pp['dE_cusp']:.6f} Ha"
            )

    def test_2d_cusp_correction_r_dependent(self, pes_scan_2d_lmax2):
        """Cusp correction should be larger at short R and smaller at large R."""
        pes_result = pes_scan_2d_lmax2['pes_result']
        per_point = pes_result['per_point_cusp']

        dE_short = None
        dE_long = None
        for pp in per_point:
            if abs(pp['R'] - 1.0) < 0.01:
                dE_short = pp['dE_cusp']
            if abs(pp['R'] - 10.0) < 0.01:
                dE_long = pp['dE_cusp']

        assert dE_short is not None and dE_long is not None
        assert abs(dE_short) > abs(dE_long), (
            f"|dE(R=1.0)| = {abs(dE_short):.6f} should be > "
            f"|dE(R=10.0)| = {abs(dE_long):.6f}"
        )

    def test_2d_corrected_de_improves(self, pes_scan_2d_lmax2):
        """Corrected D_e should be greater than uncorrected D_e."""
        pes_result = pes_scan_2d_lmax2['pes_result']

        D_e_unc = pes_result['D_e_uncorrected']
        D_e_corr = pes_result['D_e_corrected']

        print(f"\n  D_e uncorrected: {D_e_unc:.6f} Ha "
              f"({pes_result['D_e_pct_uncorrected']:.1f}%)")
        print(f"  D_e corrected:   {D_e_corr:.6f} Ha "
              f"({pes_result['D_e_pct_corrected']:.1f}%)")
        print(f"  Improvement:     "
              f"{pes_result['D_e_pct_corrected'] - pes_result['D_e_pct_uncorrected']:+.1f} pp")

        assert D_e_corr > D_e_unc, (
            f"Corrected D_e ({D_e_corr:.6f}) should be > "
            f"uncorrected D_e ({D_e_unc:.6f})"
        )

    def test_2d_cusp_vanishes_at_dissociation(self, pes_scan_2d_lmax2):
        """Cusp correction should become very small at large R."""
        pes_result = pes_scan_2d_lmax2['pes_result']
        per_point = pes_result['per_point_cusp']

        dE_10 = None
        dE_1p4 = None
        for pp in per_point:
            if abs(pp['R'] - 10.0) < 0.01:
                dE_10 = pp['dE_cusp']
            if abs(pp['R'] - 1.4) < 0.01:
                dE_1p4 = pp['dE_cusp']

        assert dE_10 is not None and dE_1p4 is not None
        ratio = abs(dE_10) / abs(dE_1p4)
        assert ratio < 0.05, (
            f"|dE(R=10)| / |dE(R=1.4)| = {ratio:.4f}, expected < 0.05"
        )


# ============================================================================
# Reference: l_max=4 sigma+pi headline results (too expensive for CI)
# ============================================================================
#
# These results were computed in the Track X2 diagnostic run.
# Run manually with:
#   python -c "from tests.test_cusp_correction_2d import _run_2d_pes_scan, _apply_cusp; ..."
#
# Results (n_alpha=60, n_Re=120, ~180-290s per point):
#
# | Solver    | l_max | Channels | D_e uncorr | dE_cusp(R_eq) | D_e corr |
# |-----------|-------|----------|------------|---------------|----------|
# | Adiabatic | 2     | sigma    | 89.7%      | -1.74 mHa     | 90.7%    |
# | 2D        | 2     | sigma    | 79.8%      | -1.74 mHa     | 80.8%    |
# | 2D        | 4     | sigma+pi | 94.1%      | -0.39 mHa     | 94.3%    |
#
# The cusp correction is solver-independent (same l_max -> same Schwartz
# formula -> same correction at each R). The correction magnitude at
# l_max=4 is much smaller than l_max=2 because the partial-wave tail
# starts at l=5 instead of l=3: correction scales as ~1/(l_max+3/2)^4.
#
# The 95.1% projection in Paper 15 was based on the l_max=2 differential
# correction (1.74 mHa), incorrectly applied to the l_max=4 baseline.
# The correct l_max=4 correction is only 0.39 mHa at R_eq, giving 94.3%.
#
# Full per-point data saved in debug/track_x2_2d/pes_scan_2d_lmax4.json
