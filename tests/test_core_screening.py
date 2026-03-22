"""
Tests for core electron screening function.

Validates that Z_eff(r) from the hyperspherical two-electron wavefunction
has the correct asymptotic behavior, monotonicity, normalization, and
physical values for He (Z=2) and Li+ (Z=3).
"""

import numpy as np
import pytest

from geovac.core_screening import CoreScreening, compute_core_density, compute_z_eff


# Shared parameters for faster tests (lower resolution, still physically correct)
TEST_PARAMS = dict(l_max=2, n_alpha=100, n_radial=400, N_R_angular=120, N_R_radial=1200)


@pytest.fixture(scope="module")
def he_screening():
    """Solve He (Z=2) core screening once for all He tests."""
    cs = CoreScreening(Z=2, **TEST_PARAMS)
    cs.solve(verbose=False)
    return cs


@pytest.fixture(scope="module")
def li_screening():
    """Solve Li+ (Z=3) core screening once for all Li tests."""
    cs = CoreScreening(Z=3, **TEST_PARAMS)
    cs.solve(verbose=False)
    return cs


# --- Normalization ---

class TestNormalization:
    def test_he_density_integrates_to_2(self, he_screening):
        """The density of He must integrate to 2 (two electrons)."""
        total = np.trapezoid(he_screening.n_r, he_screening.r_grid)
        assert abs(total - 2.0) < 0.05, f"He density integral = {total}, expected 2.0"

    def test_li_density_integrates_to_2(self, li_screening):
        """The density of Li+ core must integrate to 2 (two electrons)."""
        total = np.trapezoid(li_screening.n_r, li_screening.r_grid)
        assert abs(total - 2.0) < 0.05, f"Li+ density integral = {total}, expected 2.0"


# --- Asymptotic limits ---

class TestAsymptotics:
    def test_he_z_eff_near_nucleus(self, he_screening):
        """Z_eff(r -> 0) should approach Z = 2 (no screening at nucleus)."""
        z_eff_0 = he_screening.z_eff(0.01)
        assert z_eff_0 > 1.8, f"He Z_eff(0.01) = {z_eff_0}, expected > 1.8"

    def test_li_z_eff_near_nucleus(self, li_screening):
        """Z_eff(r -> 0) should approach Z = 3 (no screening at nucleus)."""
        z_eff_0 = li_screening.z_eff(0.01)
        assert z_eff_0 > 2.8, f"Li+ Z_eff(0.01) = {z_eff_0}, expected > 2.8"

    def test_he_z_eff_far_away(self, he_screening):
        """Z_eff(r -> inf) should approach 0 for neutral He (full screening)."""
        z_eff_inf = he_screening.z_eff(15.0)
        assert z_eff_inf < 0.2, f"He Z_eff(15) = {z_eff_inf}, expected ~0"

    def test_li_z_eff_far_away(self, li_screening):
        """Z_eff(r -> inf) should approach 1 for Li+ (valence sees +1)."""
        z_eff_inf = li_screening.z_eff(15.0)
        assert abs(z_eff_inf - 1.0) < 0.2, f"Li+ Z_eff(15) = {z_eff_inf}, expected ~1.0"


# --- Monotonicity ---

class TestMonotonicity:
    def test_he_z_eff_monotonically_decreasing(self, he_screening):
        """Z_eff must decrease monotonically (more screening at larger r)."""
        r_test = np.linspace(0.05, 15.0, 200)
        z_eff = he_screening.z_eff(r_test)
        diffs = np.diff(z_eff)
        # Allow tiny positive diffs from interpolation noise
        assert np.all(diffs < 0.01), "He Z_eff is not monotonically decreasing"

    def test_li_z_eff_monotonically_decreasing(self, li_screening):
        """Z_eff must decrease monotonically (more screening at larger r)."""
        r_test = np.linspace(0.05, 15.0, 200)
        z_eff = li_screening.z_eff(r_test)
        diffs = np.diff(z_eff)
        assert np.all(diffs < 0.01), "Li+ Z_eff is not monotonically decreasing"


# --- Smoothness ---

class TestSmoothness:
    def test_he_z_eff_smooth(self, he_screening):
        """Z_eff should be smooth (no discontinuities)."""
        r_test = np.linspace(0.1, 10.0, 500)
        z_eff = he_screening.z_eff(r_test)
        # Check that max |jump| between adjacent points is small
        max_jump = np.max(np.abs(np.diff(z_eff)))
        dr = r_test[1] - r_test[0]
        # Max derivative ~Z (at steepest), so max_jump ~ Z * dr
        assert max_jump < 0.1, f"He Z_eff max jump = {max_jump}, too large"

    def test_li_z_eff_smooth(self, li_screening):
        """Z_eff should be smooth (no discontinuities)."""
        r_test = np.linspace(0.1, 10.0, 500)
        z_eff = li_screening.z_eff(r_test)
        max_jump = np.max(np.abs(np.diff(z_eff)))
        assert max_jump < 0.1, f"Li+ Z_eff max jump = {max_jump}, too large"


# --- Physical values ---

class TestPhysicalValues:
    def test_he_z_eff_at_1bohr(self, he_screening):
        """He Z_eff at 1 bohr should be between 0 and 2."""
        z_eff_1 = he_screening.z_eff(1.0)
        assert 0.0 < z_eff_1 < 2.0, f"He Z_eff(1.0) = {z_eff_1}, out of range"

    def test_li_z_eff_at_1bohr(self, li_screening):
        """Li+ Z_eff at 1 bohr should be between 1 and 3."""
        z_eff_1 = li_screening.z_eff(1.0)
        assert 1.0 < z_eff_1 < 3.0, f"Li+ Z_eff(1.0) = {z_eff_1}, out of range"

    def test_he_energy_reasonable(self, he_screening):
        """He ground-state energy should be close to exact -2.9037 Ha."""
        E = he_screening.energy
        err = abs(E - (-2.9037)) / 2.9037 * 100
        assert err < 1.0, f"He energy = {E:.4f}, error = {err:.2f}%"

    def test_li_energy_reasonable(self, li_screening):
        """Li+ ground-state energy should be close to exact -7.2799 Ha."""
        E = li_screening.energy
        err = abs(E - (-7.2799)) / 7.2799 * 100
        assert err < 1.0, f"Li+ energy = {E:.4f}, error = {err:.2f}%"


# --- Standalone function API ---

class TestStandaloneFunctions:
    def test_compute_core_density_returns_tuple(self):
        """compute_core_density should return (r_grid, n_r) tuple."""
        r, n = compute_core_density(Z=2, l_max=1, n_alpha=50, n_radial=100,
                                    N_R_angular=60, N_R_radial=600)
        assert r.shape == (100,)
        assert n.shape == (100,)
        assert np.all(n >= 0)

    def test_compute_z_eff_returns_array(self):
        """compute_z_eff should return Z_eff on the given r_grid."""
        r_test = np.array([0.1, 1.0, 5.0])
        z_eff = compute_z_eff(Z=2, r_grid=r_test, l_max=1, n_alpha=50,
                              n_radial=100, N_R_angular=60, N_R_radial=600)
        assert z_eff.shape == (3,)
        assert z_eff[0] > z_eff[1] > z_eff[2]


# --- CoreScreening class API ---

class TestCoreScreeningClass:
    def test_unsolved_raises(self):
        """Accessing results before solve() should raise RuntimeError."""
        cs = CoreScreening(Z=2)
        with pytest.raises(RuntimeError):
            cs.z_eff(1.0)

    def test_scalar_input(self, he_screening):
        """z_eff with scalar input should return float."""
        result = he_screening.z_eff(1.0)
        assert isinstance(result, float)

    def test_array_input(self, he_screening):
        """z_eff with array input should return ndarray."""
        result = he_screening.z_eff(np.array([0.5, 1.0, 2.0]))
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)

    def test_density_positive(self, he_screening):
        """density() should return non-negative values."""
        r_test = np.linspace(0.1, 10.0, 50)
        rho = he_screening.density(r_test)
        assert np.all(rho >= 0)


# --- Summary printout ---

def test_print_summary(he_screening, li_screening):
    """Print Z_eff summary for manual physics verification."""
    r_points = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0])

    print("\n" + "=" * 60)
    print("Core Screening Summary")
    print("=" * 60)

    for label, cs in [("He (Z=2)", he_screening), ("Li+ (Z=3)", li_screening)]:
        z_eff = cs.z_eff(r_points)
        print(f"\n{label}:")
        print(f"  Energy = {cs.energy:.6f} Ha")
        print(f"  {'r (bohr)':<12} {'Z_eff(r)':<12}")
        print(f"  {'-'*24}")
        for r, z in zip(r_points, z_eff):
            print(f"  {r:<12.2f} {z:<12.4f}")

    total_he = np.trapezoid(he_screening.n_r, he_screening.r_grid)
    total_li = np.trapezoid(li_screening.n_r, li_screening.r_grid)
    print(f"\nDensity normalization: He = {total_he:.4f}, Li+ = {total_li:.4f}")
    print("=" * 60)
