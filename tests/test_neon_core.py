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

from geovac.neon_core import FrozenCore, _NIST_CORE_ENERGIES


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
