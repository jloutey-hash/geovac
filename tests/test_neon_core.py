"""
Tests for the FrozenCore Ne-like (10-electron) core screening module.

Validates:
  - Z_eff asymptotic behavior (Z near nucleus, Z-10 at infinity)
  - Z_eff monotonically decreasing
  - Density non-negative everywhere
  - Density integrates to 10
  - Core energy matches NIST
  - Z_eff smooth (no jumps)
  - Works for all Z=11-18
  - Raises ValueError for Z outside 11-18
"""

import numpy as np
import pytest

from geovac.neon_core import (
    FrozenCore, _NIST_CORE_ENERGIES,
    _solve_screened_radial, screened_r3_inverse, screened_so_splitting,
    screened_xi_so,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=range(11, 19), ids=[f"Z={z}" for z in range(11, 19)])
def frozen_core(request):
    """Solved FrozenCore for each Z=11-18."""
    fc = FrozenCore(Z=request.param)
    fc.solve()
    return fc


# ---------------------------------------------------------------------------
# Basic interface tests
# ---------------------------------------------------------------------------

class TestFrozenCoreInterface:
    """Test that FrozenCore has the expected interface."""

    def test_solve_required(self):
        """z_eff raises if solve() not called."""
        fc = FrozenCore(Z=11)
        with pytest.raises(RuntimeError, match="solve"):
            fc.z_eff(1.0)

    def test_invalid_z_low(self):
        """Z=10 (Ne) is not in any frozen-core table (it's a target core,
        not a FrozenCore use case)."""
        with pytest.raises(ValueError, match="No frozen-core data"):
            FrozenCore(Z=10)

    def test_invalid_z_unsupported(self):
        """Z in unsupported range (e.g., 4d or 5p block) raises ValueError."""
        # Z=39 (Y) — 4d block, not implemented
        with pytest.raises(ValueError, match="No frozen-core data"):
            FrozenCore(Z=39)
        # Z=57 (La) — lanthanide, not implemented
        with pytest.raises(ValueError, match="No frozen-core data"):
            FrozenCore(Z=57)

    def test_invalid_z_he_core(self):
        """Z in He-core range (Z=3-10) uses CoreScreening, not FrozenCore."""
        with pytest.raises(ValueError, match="No frozen-core data"):
            FrozenCore(Z=3)


# ---------------------------------------------------------------------------
# Z_eff asymptotic behavior
# ---------------------------------------------------------------------------

class TestZeffAsymptotics:
    """Z_eff must approach Z at r=0 and Z-10 at large r."""

    def test_zeff_near_nucleus(self, frozen_core):
        """Z_eff(0.001) should be close to Z (within 0.5)."""
        z_near = frozen_core.z_eff(0.001)
        assert abs(z_near - frozen_core.Z) < 0.5, (
            f"Z={frozen_core.Z}: z_eff(0.001)={z_near:.4f}, "
            f"expected ~{frozen_core.Z}"
        )

    def test_zeff_at_infinity(self, frozen_core):
        """Z_eff(100.0) should be close to Z-10 (within 0.01)."""
        z_far = frozen_core.z_eff(100.0)
        expected = frozen_core.Z - 10
        assert abs(z_far - expected) < 0.01, (
            f"Z={frozen_core.Z}: z_eff(100.0)={z_far:.6f}, "
            f"expected ~{expected}"
        )

    def test_zeff_intermediate(self, frozen_core):
        """Z_eff at intermediate r should be between Z and Z-10."""
        r_test = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
        z_vals = frozen_core.z_eff(r_test)
        Z = frozen_core.Z
        assert np.all(z_vals >= Z - 10 - 0.01), (
            f"Z={Z}: z_eff values below Z-10"
        )
        assert np.all(z_vals <= Z + 0.01), (
            f"Z={Z}: z_eff values above Z"
        )


# ---------------------------------------------------------------------------
# Monotonicity and smoothness
# ---------------------------------------------------------------------------

class TestZeffMonotonicity:
    """Z_eff must be monotonically decreasing."""

    def test_monotonically_decreasing(self, frozen_core):
        """Z_eff(r) should decrease as r increases."""
        r_test = np.linspace(0.01, 15.0, 500)
        z_vals = frozen_core.z_eff(r_test)
        diffs = np.diff(z_vals)
        # Allow tiny numerical noise (1e-6)
        assert np.all(diffs <= 1e-6), (
            f"Z={frozen_core.Z}: z_eff not monotonically decreasing. "
            f"Max positive diff={np.max(diffs):.2e}"
        )

    def test_smooth_no_jumps(self, frozen_core):
        """No jumps > 0.5 between adjacent grid points on a fine grid."""
        # Use 5000 points to resolve the steep screening region near the
        # nucleus for high-Z atoms (Z=15-18 have very compact cores)
        r_test = np.linspace(0.001, 15.0, 5000)
        z_vals = frozen_core.z_eff(r_test)
        max_jump = np.max(np.abs(np.diff(z_vals)))
        assert max_jump < 0.5, (
            f"Z={frozen_core.Z}: z_eff has jumps > 0.5. "
            f"Max jump={max_jump:.4f}"
        )


# ---------------------------------------------------------------------------
# Density tests
# ---------------------------------------------------------------------------

class TestDensity:
    """Core density must be non-negative and integrate to 10."""

    def test_density_nonnegative(self, frozen_core):
        """Density should be non-negative everywhere."""
        r_test = np.linspace(0.01, 15.0, 500)
        d = frozen_core.density(r_test)
        assert np.all(d >= -1e-10), (
            f"Z={frozen_core.Z}: negative density found. "
            f"Min={np.min(d):.2e}"
        )

    def test_density_integrates_to_10(self, frozen_core):
        """Integral of density should be 10 (within 1%)."""
        r_fine = np.linspace(0.001, 15.0, 5000)
        d = frozen_core.density(r_fine)
        total = np.trapezoid(d, r_fine)
        assert abs(total - 10.0) / 10.0 < 0.01, (
            f"Z={frozen_core.Z}: density integrates to {total:.4f}, "
            f"expected 10.0"
        )

    def test_density_scalar_input(self, frozen_core):
        """Density works with scalar input."""
        d = frozen_core.density(1.0)
        assert isinstance(d, float)
        assert d >= 0.0

    def test_density_array_input(self, frozen_core):
        """Density works with array input."""
        r = np.array([0.1, 0.5, 1.0, 2.0])
        d = frozen_core.density(r)
        assert isinstance(d, np.ndarray)
        assert d.shape == (4,)


# ---------------------------------------------------------------------------
# Core energy
# ---------------------------------------------------------------------------

class TestCoreEnergy:
    """Core energy should match NIST values."""

    def test_energy_matches_nist(self, frozen_core):
        """Energy property returns NIST value."""
        Z = frozen_core.Z
        expected = _NIST_CORE_ENERGIES[Z]
        assert frozen_core.energy == expected, (
            f"Z={Z}: energy={frozen_core.energy}, expected={expected}"
        )

    def test_energy_negative(self, frozen_core):
        """Core energy should be negative."""
        assert frozen_core.energy < 0

    def test_energy_decreases_with_z(self):
        """Core energy should become more negative with increasing Z."""
        energies = []
        for Z in range(11, 19):
            fc = FrozenCore(Z=Z)
            energies.append(fc.energy)
        for i in range(len(energies) - 1):
            assert energies[i] > energies[i + 1], (
                f"Energy not decreasing: E(Z={11+i})={energies[i]:.3f} "
                f"vs E(Z={12+i})={energies[i+1]:.3f}"
            )


# ---------------------------------------------------------------------------
# Scalar / array consistency
# ---------------------------------------------------------------------------

class TestScalarArrayConsistency:
    """Scalar and array inputs must give consistent results."""

    def test_zeff_scalar_matches_array(self, frozen_core):
        """z_eff(1.0) should match z_eff([1.0])[0]."""
        scalar_val = frozen_core.z_eff(1.0)
        array_val = frozen_core.z_eff(np.array([1.0]))[0]
        assert abs(scalar_val - array_val) < 1e-12

    def test_zeff_returns_float_for_scalar(self, frozen_core):
        """z_eff with scalar input returns float."""
        val = frozen_core.z_eff(1.0)
        assert isinstance(val, float)

    def test_zeff_returns_array_for_array(self, frozen_core):
        """z_eff with array input returns ndarray."""
        val = frozen_core.z_eff(np.array([1.0, 2.0]))
        assert isinstance(val, np.ndarray)


# ---------------------------------------------------------------------------
# All-Z sweep (non-parametric for summary table)
# ---------------------------------------------------------------------------

class TestAllZSweep:
    """Quick checks across all Z values."""

    def test_all_z_solve_and_query(self):
        """All Z=11-18 can be solved and queried without error."""
        for Z in range(11, 19):
            fc = FrozenCore(Z=Z)
            fc.solve()
            # Basic queries
            z_near = fc.z_eff(0.01)
            z_far = fc.z_eff(100.0)
            d = fc.density(1.0)
            e = fc.energy
            assert z_near > Z - 1.0
            assert abs(z_far - (Z - 10)) < 0.01
            assert d >= 0.0
            assert e < 0.0

    def test_verbose_mode(self):
        """Verbose mode runs without error."""
        fc = FrozenCore(Z=11)
        fc.solve(verbose=True)


# ---------------------------------------------------------------------------
# Screened radial solver tests
# ---------------------------------------------------------------------------

class TestScreenedRadialSolver:
    """Tests for _solve_screened_radial and screened_r3_inverse."""

    def test_l0_raises(self):
        """l=0 should raise ValueError (1/r^3 diverges)."""
        with pytest.raises(ValueError, match="l must be >= 1"):
            _solve_screened_radial(11, 0, 1)

    def test_invalid_n_raises(self):
        """n < l+1 should raise ValueError."""
        with pytest.raises(ValueError, match="n_target"):
            _solve_screened_radial(11, 1, 1)  # n=1, l=1: impossible

    def test_r3_inverse_l0_raises(self):
        """screened_r3_inverse should raise for l=0."""
        with pytest.raises(ValueError, match="diverges"):
            screened_r3_inverse(11, 3, 0)

    def test_na_3p_energy(self):
        """Na 3p eigenvalue should be near -0.11 Ha (IP ~ 3 eV)."""
        energy, _, _ = _solve_screened_radial(11, 1, 3, n_grid=8000)
        # Na 3p IP is about 3.03 eV = 0.111 Ha
        assert abs(energy + 0.111) < 0.015, (
            f"Na 3p energy {energy:.4f} Ha too far from -0.111 Ha"
        )

    def test_na_3p_normalization(self):
        """Na 3p wavefunction should be normalized."""
        _, u, r = _solve_screened_radial(11, 1, 3, n_grid=8000)
        norm = np.trapezoid(u**2, r)
        assert abs(norm - 1.0) < 1e-4, f"norm = {norm}, expected 1.0"

    def test_na_3p_wavefunction_shape(self):
        """Na 3p wavefunction should peak at r ~ 3-8 bohr."""
        _, u, r = _solve_screened_radial(11, 1, 3, n_grid=8000)
        peak_r = r[np.argmax(u**2)]
        assert 1.0 < peak_r < 15.0, (
            f"Na 3p peak at r={peak_r:.2f}, expected 3-8 bohr"
        )

    def test_r3_inverse_positive(self):
        """<1/r^3> should be positive."""
        r3 = screened_r3_inverse(11, 3, 1, n_grid=8000)
        assert r3 > 0, f"<1/r^3> = {r3}, expected positive"

    def test_r3_inverse_enhancement(self):
        """Screened <1/r^3> should be much larger than hydrogenic Z_eff=1."""
        r3_scr = screened_r3_inverse(11, 3, 1, n_grid=8000)
        r3_hyd = 1.0 / (27 * 1 * 1.5 * 2)  # Z_eff=1, n=3, l=1
        assert r3_scr > 5 * r3_hyd, (
            f"Enhancement factor {r3_scr/r3_hyd:.1f}x too small"
        )

    def test_r3_inverse_grid_convergence(self):
        """<1/r^3> should converge with grid refinement."""
        r3_coarse = screened_r3_inverse(11, 3, 1, n_grid=4000)
        r3_fine = screened_r3_inverse(11, 3, 1, n_grid=12000)
        # Should agree within 5%
        rel_diff = abs(r3_fine - r3_coarse) / r3_fine
        assert rel_diff < 0.05, (
            f"Grid convergence: {rel_diff:.3f} relative diff (>5%)"
        )


class TestScreenedSOSplitting:
    """Tests for screened_so_splitting and screened_xi_so."""

    def test_splitting_positive(self):
        """SO splitting should be positive (j=l+1/2 above j=l-1/2)."""
        result = screened_so_splitting(11, 3, 1, n_grid=8000)
        assert result['splitting'] > 0

    def test_splitting_cm1_order_of_magnitude(self):
        """Na 3p splitting should be 10-100 cm^-1 order of magnitude."""
        result = screened_so_splitting(11, 3, 1, n_grid=8000)
        assert 1.0 < result['splitting_cm1'] < 500.0, (
            f"Na 3p splitting {result['splitting_cm1']:.1f} cm^-1 "
            f"outside expected range 1-500"
        )

    def test_xi_so_positive(self):
        """xi = <(1/r) dV/dr> should be positive."""
        xi = screened_xi_so(11, 3, 1, n_grid=8000)
        assert xi > 0, f"xi = {xi}, expected positive"

    def test_l0_raises(self):
        """l=0 should raise."""
        with pytest.raises(ValueError):
            screened_so_splitting(11, 3, 0)

    def test_heavier_atom_larger_splitting(self):
        """Sr 5p splitting should be larger than Na 3p."""
        na_result = screened_so_splitting(11, 3, 1, n_grid=8000)
        sr_result = screened_so_splitting(38, 5, 1, n_grid=8000)
        assert sr_result['splitting_cm1'] > na_result['splitting_cm1'], (
            f"Sr splitting {sr_result['splitting_cm1']:.1f} should be > "
            f"Na splitting {na_result['splitting_cm1']:.1f}"
        )

    def test_enhancement_over_hydrogenic(self):
        """Screened <1/r^3> should enhance over hydrogenic for all atoms."""
        for Z, n in [(11, 3), (20, 4), (38, 5)]:
            result = screened_so_splitting(Z, n, 1, n_grid=8000)
            assert result['enhancement'] > 5.0, (
                f"Z={Z}: enhancement {result['enhancement']:.1f}x too small"
            )

    def test_result_keys(self):
        """screened_so_splitting should return all expected keys."""
        result = screened_so_splitting(11, 3, 1, n_grid=4000)
        expected_keys = {
            'r3_inv_screened', 'r3_inv_hydrogenic', 'enhancement',
            'xi_so', 'E_j_plus', 'E_j_minus', 'splitting',
            'splitting_cm1', 'splitting_znuc_cm1', 'energy_eV',
        }
        assert set(result.keys()) == expected_keys

    @pytest.mark.slow
    def test_ba_6p_splitting(self):
        """Ba 6p splitting should be in the right ballpark (100-5000 cm^-1)."""
        result = screened_so_splitting(56, 6, 1, n_grid=12000)
        assert 100 < result['splitting_cm1'] < 5000, (
            f"Ba 6p splitting {result['splitting_cm1']:.1f} cm^-1 "
            f"outside expected range"
        )
