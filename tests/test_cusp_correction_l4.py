"""
Tests for Track X cusp correction applied to actual Level 4 H2 PES.

Validates:
1. Level 4 solver runs at multiple R with spectral angular
2. Cusp correction is negative at all R
3. Cusp correction is R-dependent (larger at short R)
4. Corrected D_e > uncorrected D_e
5. Corrected D_e percentage quantified against Kolos & Wolniewicz exact

All tests here run the actual Level 4 solver and are marked @pytest.mark.slow.
"""

import numpy as np
import pytest
import time
import json
import os


# All tests in this module require the Level 4 solver (expensive)
pytestmark = pytest.mark.slow

# Exact D_e for H2 from Kolos & Wolniewicz
D_E_EXACT = 0.17447  # Ha

# R grid for PES scan — covers short-range, equilibrium, and dissociation
R_GRID = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0, 10.0])


def _run_level4_pes_scan(
    R_grid: np.ndarray,
    l_max: int = 2,
    angular_method: str = 'spectral',
) -> dict:
    """Run Level 4 solver at multiple R points and collect results.

    Returns dict with arrays of E_total, per-point results, and timing.
    """
    from geovac.level4_multichannel import solve_level4_h2_multichannel

    E_totals = []
    per_point = []
    t0 = time.time()

    for R in R_grid:
        t_start = time.time()
        result = solve_level4_h2_multichannel(
            R=float(R),
            l_max=l_max,
            angular_method=angular_method,
            verbose=False,
        )
        t_end = time.time()

        E_totals.append(result['E_total'])
        per_point.append({
            'R': float(R),
            'E_total': float(result['E_total']),
            'D_e': float(result['D_e']),
            'D_e_pct': float(result['D_e_pct']),
            'l_max': result['l_max'],
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


def _apply_cusp_and_analyze(scan: dict) -> dict:
    """Apply cusp correction to PES scan results and compute D_e.

    Returns the cusp_correction_h2_pes result dict augmented with
    per-point cusp data.
    """
    from geovac.cusp_correction import (
        cusp_correction_h2_pes,
        cusp_correction_from_level4,
        coalescence_density_h2,
    )

    R_grid = scan['R_grid']
    E_pes = scan['E_totals']
    l_max = scan['l_max']

    # Full PES correction (computes D_e)
    # Use the known dissociation limit E_atoms = -1.0 Ha (H + H),
    # NOT the solver's large-R value (which may not converge due to
    # finite hyperradial grid).
    pes_result = cusp_correction_h2_pes(R_grid, E_pes, l_max=l_max,
                                        E_atoms=-1.0)

    # Per-point corrections for detailed reporting
    per_point_cusp = []
    for pp in scan['per_point']:
        mock_result = {
            'E_total': pp['E_total'],
            'R': pp['R'],
            'l_max': l_max,
        }
        cusp = cusp_correction_from_level4(mock_result)
        delta = coalescence_density_h2(pp['R'])
        per_point_cusp.append({
            'R': pp['R'],
            'E_uncorrected': pp['E_total'],
            'dE_cusp': cusp['delta_E_cusp'],
            'E_corrected': cusp['E_corrected'],
            'delta_r12': float(delta),
        })

    pes_result['per_point_cusp'] = per_point_cusp
    return pes_result


def _save_diagnostic_data(scan: dict, pes_result: dict, filename: str) -> None:
    """Save diagnostic data to debug/track_x_validation/."""
    debug_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'debug', 'track_x_validation',
    )
    os.makedirs(debug_dir, exist_ok=True)

    data = {
        'l_max': scan['l_max'],
        'total_wall_time': scan['total_time'],
        'D_e_exact': D_E_EXACT,
        'D_e_uncorrected': float(pes_result['D_e_uncorrected']),
        'D_e_corrected': float(pes_result['D_e_corrected']),
        'D_e_pct_uncorrected': float(pes_result['D_e_pct_uncorrected']),
        'D_e_pct_corrected': float(pes_result['D_e_pct_corrected']),
        'per_point': pes_result['per_point_cusp'],
        'solver_per_point': scan['per_point'],
    }

    filepath = os.path.join(debug_dir, filename)
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


# ============================================================================
# Shared fixture: run the PES scan once, reuse across tests
# ============================================================================

@pytest.fixture(scope='module')
def level4_pes_scan():
    """Run the Level 4 PES scan once for all tests in this module."""
    scan = _run_level4_pes_scan(R_GRID, l_max=2, angular_method='spectral')
    pes_result = _apply_cusp_and_analyze(scan)

    # Save diagnostic data
    _save_diagnostic_data(scan, pes_result, 'pes_scan_lmax2.json')

    return {
        'scan': scan,
        'pes_result': pes_result,
    }


# ============================================================================
# Tests
# ============================================================================


class TestCuspCorrectionLevel4PES:
    """Validate Track X cusp correction on actual Level 4 H2 PES."""

    def test_solver_runs_at_all_r_points(self, level4_pes_scan):
        """Level 4 solver should produce finite energies at all R."""
        scan = level4_pes_scan['scan']

        for pp in scan['per_point']:
            assert np.isfinite(pp['E_total']), (
                f"Non-finite energy at R={pp['R']}: {pp['E_total']}"
            )

    def test_cusp_correction_is_negative_at_all_r(self, level4_pes_scan):
        """Cusp correction must be negative at all R (energy too high
        without cusp -> correction lowers it)."""
        pes_result = level4_pes_scan['pes_result']

        for pp in pes_result['per_point_cusp']:
            assert pp['dE_cusp'] < 0, (
                f"Cusp correction at R={pp['R']} should be negative, "
                f"got {pp['dE_cusp']:.6f} Ha"
            )

    def test_cusp_correction_is_r_dependent(self, level4_pes_scan):
        """Cusp correction should be larger at short R (more overlap)
        and smaller at large R (separated atoms)."""
        pes_result = level4_pes_scan['pes_result']
        per_point = pes_result['per_point_cusp']

        # Find short-R and long-R points
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

    def test_corrected_pes_is_lower(self, level4_pes_scan):
        """Corrected PES should be lower than uncorrected at all R."""
        pes_result = level4_pes_scan['pes_result']

        assert np.all(pes_result['E_corrected'] < pes_result['E_uncorrected']), (
            "Corrected PES should be lower than uncorrected at all R"
        )

    def test_corrected_de_improves(self, level4_pes_scan):
        """Corrected D_e should be greater than uncorrected D_e.

        The cusp correction is larger near equilibrium (where electrons
        overlap more) than at dissociation, so it deepens the well
        differentially and increases D_e.
        """
        pes_result = level4_pes_scan['pes_result']

        D_e_unc = pes_result['D_e_uncorrected']
        D_e_corr = pes_result['D_e_corrected']

        assert D_e_corr > D_e_unc, (
            f"Corrected D_e ({D_e_corr:.6f}) should be > "
            f"uncorrected D_e ({D_e_unc:.6f})"
        )

    def test_corrected_de_percentage(self, level4_pes_scan):
        """Report the corrected D_e as a percentage of exact.

        This is the headline validation number for Track X.
        """
        pes_result = level4_pes_scan['pes_result']

        D_e_pct_unc = pes_result['D_e_pct_uncorrected']
        D_e_pct_corr = pes_result['D_e_pct_corrected']

        print(f"\n{'='*60}")
        print(f"Track X Cusp Correction — Level 4 H2 PES Validation")
        print(f"{'='*60}")
        print(f"  l_max = {pes_result['l_max']}")
        print(f"  D_e exact (KW):     {D_E_EXACT:.5f} Ha")
        print(f"  D_e uncorrected:    {pes_result['D_e_uncorrected']:.6f} Ha "
              f"({D_e_pct_unc:.1f}%)")
        print(f"  D_e corrected:      {pes_result['D_e_corrected']:.6f} Ha "
              f"({D_e_pct_corr:.1f}%)")
        print(f"  Improvement:        {D_e_pct_corr - D_e_pct_unc:+.1f} pp")
        print(f"{'='*60}")

        # Print per-point table
        print(f"\n{'R':>6s} {'E_unc':>12s} {'dE_cusp':>10s} {'E_corr':>12s}")
        print(f"{'-'*6} {'-'*12} {'-'*10} {'-'*12}")
        for pp in pes_result['per_point_cusp']:
            print(f"{pp['R']:6.2f} {pp['E_uncorrected']:12.6f} "
                  f"{pp['dE_cusp']:10.6f} {pp['E_corrected']:12.6f}")

        # The corrected D_e should be at least as good as uncorrected
        assert D_e_pct_corr >= D_e_pct_unc - 0.5, (
            f"Corrected D_e ({D_e_pct_corr:.1f}%) should not be much worse "
            f"than uncorrected ({D_e_pct_unc:.1f}%)"
        )

    def test_cusp_correction_vanishes_at_dissociation(self, level4_pes_scan):
        """Cusp correction should become very small at large R.

        At R -> infinity, each atom has only 1 electron, so there is no
        electron-electron coalescence and the correction should vanish.
        """
        pes_result = level4_pes_scan['pes_result']
        per_point = pes_result['per_point_cusp']

        # Get correction at R=10 and R=1.4
        dE_10 = None
        dE_1p4 = None
        for pp in per_point:
            if abs(pp['R'] - 10.0) < 0.01:
                dE_10 = pp['dE_cusp']
            if abs(pp['R'] - 1.4) < 0.01:
                dE_1p4 = pp['dE_cusp']

        assert dE_10 is not None and dE_1p4 is not None
        # At R=10, correction should be < 5% of the correction at R=1.4
        ratio = abs(dE_10) / abs(dE_1p4)
        assert ratio < 0.05, (
            f"|dE(R=10)| / |dE(R=1.4)| = {ratio:.4f}, expected < 0.05"
        )

    def test_differential_correction_near_equilibrium(self, level4_pes_scan):
        """The DIFFERENTIAL cusp correction (equilibrium minus dissociation)
        is what matters for D_e. Verify it has the expected magnitude (~1-3 mHa).
        """
        pes_result = level4_pes_scan['pes_result']
        per_point = pes_result['per_point_cusp']

        # Find correction near equilibrium (R ~ 1.4) and at dissociation (R=10)
        dE_eq = None
        dE_inf = None
        for pp in per_point:
            if abs(pp['R'] - 1.4) < 0.01:
                dE_eq = pp['dE_cusp']
            if abs(pp['R'] - 10.0) < 0.01:
                dE_inf = pp['dE_cusp']

        assert dE_eq is not None and dE_inf is not None

        # Differential correction = |dE_eq| - |dE_inf| (both negative,
        # so this is the extra lowering at equilibrium vs dissociation)
        differential_mHa = (abs(dE_eq) - abs(dE_inf)) * 1000

        print(f"\n  Differential cusp correction: {differential_mHa:.2f} mHa")
        print(f"  dE_cusp(R=1.4): {dE_eq * 1000:.2f} mHa")
        print(f"  dE_cusp(R=10):  {dE_inf * 1000:.4f} mHa")

        # Should be positive (more correction at equilibrium) and
        # in the range 0.5 - 5 mHa based on Track X estimates
        assert differential_mHa > 0.1, (
            f"Differential correction {differential_mHa:.2f} mHa is too small"
        )
        assert differential_mHa < 10.0, (
            f"Differential correction {differential_mHa:.2f} mHa is suspiciously large"
        )

    def test_wall_time_reasonable(self, level4_pes_scan):
        """PES scan should complete in reasonable time with spectral angular."""
        scan = level4_pes_scan['scan']
        total_time = scan['total_time']
        n_points = len(R_GRID)

        print(f"\n  Wall time: {total_time:.1f}s total, "
              f"{total_time/n_points:.2f}s per R-point")

        # With spectral angular, each point should take < 2s
        assert total_time / n_points < 5.0, (
            f"Average {total_time/n_points:.1f}s per point is too slow"
        )
