"""
Tests for cusp correction wiring in the composed-geometry pipeline (Track AB).

Validates that the Schwartz partial-wave cusp correction (Track X) is correctly
integrated into ComposedDiatomicSolver.scan_pes() as a per-fiber post-processing
correction to the electronic energy.

Test structure:
  - Fast unit tests (no PES scan): parameter wiring, standalone consistency
  - Slow integration tests: full LiH PES scan with/without cusp correction

References:
  - geovac/cusp_correction.py (Track X, v2.0.14)
  - geovac/composed_diatomic.py (Level 5 pipeline)
"""

import numpy as np
import pytest

from geovac.composed_diatomic import ComposedDiatomicSolver
from geovac.cusp_correction import cusp_correction_h2_point


# ======================================================================
# Fast unit tests — no PES scan required
# ======================================================================

class TestCuspCorrectionWiring:
    """Verify parameter plumbing without running expensive PES scans."""

    def test_default_cusp_off(self):
        """cusp_correction defaults to False."""
        solver = ComposedDiatomicSolver.LiH_ab_initio(l_max=2)
        assert solver.cusp_correction is False

    def test_cusp_on_via_kwarg(self):
        """cusp_correction=True propagates through class methods."""
        solver = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=2, cusp_correction=True)
        assert solver.cusp_correction is True

    def test_cusp_on_beh_plus(self):
        """cusp_correction propagates to BeH+ constructor."""
        solver = ComposedDiatomicSolver.BeH_plus_ab_initio(
            l_max=2, cusp_correction=True)
        assert solver.cusp_correction is True

    def test_cusp_on_manual_pk(self):
        """cusp_correction works with manual PK mode."""
        solver = ComposedDiatomicSolver.LiH(
            l_max=2, cusp_correction=True)
        assert solver.cusp_correction is True

    def test_standalone_cusp_values(self):
        """Verify standalone cusp_correction_h2_point gives expected range."""
        # At l_max=2, R~3 bohr, correction should be ~-1 to -2 mHa
        dc = cusp_correction_h2_point(R=3.0, l_max=2)
        assert dc < 0, "Cusp correction must be negative"
        assert abs(dc) < 0.01, "Correction should be < 10 mHa"
        assert abs(dc) > 0.0001, "Correction should be > 0.1 mHa"

    def test_cusp_vanishes_at_large_R(self):
        """Cusp correction vanishes at dissociation (large R)."""
        dc_close = cusp_correction_h2_point(R=3.0, l_max=2)
        dc_far = cusp_correction_h2_point(R=20.0, l_max=2)
        assert abs(dc_far) < abs(dc_close), (
            "Correction should be smaller at large R")

    def test_cusp_scales_with_lmax(self):
        """Cusp correction decreases with increasing l_max."""
        dc_l2 = cusp_correction_h2_point(R=3.0, l_max=2)
        dc_l4 = cusp_correction_h2_point(R=3.0, l_max=4)
        assert abs(dc_l4) < abs(dc_l2), (
            "Higher l_max should give smaller correction")


# ======================================================================
# Slow integration tests — full PES scans
# ======================================================================

def _short_r_grid():
    """Minimal R-grid for fast integration testing (7 points)."""
    return np.array([2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0])


@pytest.fixture(scope="module")
def lih_no_cusp():
    """LiH l_max=2, ab initio PK, NO cusp correction."""
    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2, cusp_correction=False, verbose=False)
    solver.solve_core()
    result = solver.scan_pes(R_grid=_short_r_grid(), n_Re=300)
    return solver, result


@pytest.fixture(scope="module")
def lih_with_cusp():
    """LiH l_max=2, ab initio PK, WITH cusp correction."""
    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2, cusp_correction=True, verbose=False)
    solver.solve_core()
    result = solver.scan_pes(R_grid=_short_r_grid(), n_Re=300)
    return solver, result


@pytest.mark.slow
class TestLiHCuspCorrection:
    """Integration tests for LiH with cusp correction."""

    def test_no_cusp_delta_zero(self, lih_no_cusp):
        """Without cusp correction, delta_cusp is identically zero."""
        _, result = lih_no_cusp
        assert np.all(result['delta_cusp'] == 0.0)

    def test_cusp_negative_at_all_R(self, lih_with_cusp):
        """Cusp correction is negative at every R-point."""
        _, result = lih_with_cusp
        dc = result['delta_cusp']
        valid = ~np.isnan(result['E_composed'])
        assert np.all(dc[valid] < 0), (
            f"delta_cusp should be < 0 everywhere, got {dc[valid]}")

    def test_cusp_lowers_energy(self, lih_no_cusp, lih_with_cusp):
        """E_composed(cusp=True) < E_composed(cusp=False) at all R."""
        _, res_off = lih_no_cusp
        _, res_on = lih_with_cusp
        E_off = res_off['E_composed']
        E_on = res_on['E_composed']
        valid = ~np.isnan(E_off) & ~np.isnan(E_on)
        assert np.all(E_on[valid] < E_off[valid]), (
            "Cusp-corrected energy should be lower at every R")

    def test_standalone_consistency(self, lih_with_cusp):
        """Pipeline cusp matches standalone cusp_correction_h2_point."""
        solver, result = lih_with_cusp
        R_grid = result['R']
        dc_pipeline = result['delta_cusp']
        for i, R in enumerate(R_grid):
            if np.isnan(result['E_composed'][i]):
                continue
            dc_standalone = cusp_correction_h2_point(R, l_max=solver.l_max)
            assert abs(dc_pipeline[i] - dc_standalone) < 1e-14, (
                f"Mismatch at R={R}: pipeline={dc_pipeline[i]}, "
                f"standalone={dc_standalone}")

    def test_magnitude_lmax2(self, lih_with_cusp):
        """At l_max=2, |delta_cusp| near R_eq should be O(0.1-1 mHa).

        For LiH, the valence bond pair has Z_A_eff=1, Z_B=1 (like H2),
        but R_eq~3.2 bohr is larger than H2's 1.4 bohr, so the
        coalescence density (and thus the correction) is smaller.
        Expected ~0.1-0.5 mHa near R_eq.
        """
        _, result = lih_with_cusp
        # Find point nearest R_eq (around 3.0 bohr for LiH)
        R_grid = result['R']
        i_eq = np.argmin(np.abs(R_grid - 3.0))
        dc = abs(result['delta_cusp'][i_eq])
        assert 0.00005 < dc < 0.005, (
            f"|delta_cusp| at R~3 should be 0.05-5 mHa, got {dc*1000:.3f} mHa")

    def test_uncorrected_fields_present(self, lih_with_cusp):
        """When cusp is on, uncorrected comparison fields are stored."""
        _, result = lih_with_cusp
        assert 'E_composed_uncorrected' in result
        assert 'R_eq_uncorrected' in result
        assert 'E_min_uncorrected' in result
        assert 'D_e_uncorrected' in result

    def test_uncorrected_fields_absent_when_off(self, lih_no_cusp):
        """When cusp is off, uncorrected fields are NOT stored."""
        _, result = lih_no_cusp
        assert 'E_composed_uncorrected' not in result

    def test_cusp_improves_or_preserves_D_e(self, lih_no_cusp, lih_with_cusp):
        """Cusp correction should not worsen D_e significantly.

        The correction is largest at short R and smallest at large R,
        so it should deepen the well (improve D_e) or leave it nearly
        unchanged.
        """
        _, res_off = lih_no_cusp
        _, res_on = lih_with_cusp
        D_e_off = res_off['D_e']
        D_e_on = res_on['D_e']
        # D_e should increase or stay within 10% (generous margin for
        # R-dependent correction shape)
        assert D_e_on >= D_e_off * 0.90, (
            f"D_e should not worsen by >10%: off={D_e_off:.6f}, "
            f"on={D_e_on:.6f}")

    def test_benchmark_report(self, lih_no_cusp, lih_with_cusp):
        """Print benchmark comparison (always passes, for reporting)."""
        _, res_off = lih_no_cusp
        _, res_on = lih_with_cusp
        ref_R_eq = 3.015
        ref_D_e = 0.0920

        R_eq_off = res_off['R_eq']
        R_eq_on = res_on['R_eq']
        D_e_off = res_off['D_e']
        D_e_on = res_on['D_e']
        E_min_off = res_off['E_min']
        E_min_on = res_on['E_min']

        R_eq_err_off = abs(R_eq_off - ref_R_eq) / ref_R_eq * 100
        R_eq_err_on = abs(R_eq_on - ref_R_eq) / ref_R_eq * 100

        dc = res_on['delta_cusp']
        R_grid = res_on['R']
        i_eq = np.argmin(res_on['E_composed'][~np.isnan(res_on['E_composed'])])

        print("\n" + "=" * 60)
        print("LiH l_max=2 Cusp Correction Benchmark")
        print("=" * 60)
        print(f"  R_eq:  {R_eq_off:.4f} -> {R_eq_on:.4f} bohr "
              f"(ref {ref_R_eq}, err {R_eq_err_off:.2f}% -> {R_eq_err_on:.2f}%)")
        print(f"  E_min: {E_min_off:.6f} -> {E_min_on:.6f} Ha "
              f"(delta {(E_min_on-E_min_off)*1000:+.3f} mHa)")
        print(f"  D_e:   {D_e_off:.6f} -> {D_e_on:.6f} Ha "
              f"(ref {ref_D_e})")
        print(f"  Cusp at R_eq: {dc[i_eq]*1000:.4f} mHa")
        print(f"  Cusp range: {dc.min()*1000:.4f} to {dc.max()*1000:.4f} mHa")
        print("=" * 60)
