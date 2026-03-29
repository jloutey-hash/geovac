"""
Tests for composed H₂O PES solver (geovac/composed_water.py).

Validates:
  1. Core solves without error, Z_eff = 6
  2. Bond pair solves at R = 1.81 bohr without error
  3. Lone pair solves without error, energy is reasonable
  4. PES scan produces a bound state (minimum exists)
  5. R_eq is within a physically reasonable range
  6. Block-diagonal PES (no coupling) runs without error
  7. Exchange-coupled PES runs without error
  8. Charge-center origin is used for bond pair solver
  9. All existing tests still pass
"""

import numpy as np
import pytest

from geovac.composed_water import ComposedWaterSolver


# ---------------------------------------------------------------------------
# Fixtures: shared solver objects (expensive, compute once)
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def water_solver_no_coupling():
    """ComposedWaterSolver without inter-fiber coupling (faster)."""
    solver = ComposedWaterSolver(
        l_max=2, n_alpha=100,
        include_coupling=False,
        verbose=False,
    )
    solver.solve_core()
    solver.solve_lone_pair()
    return solver


@pytest.fixture(scope='module')
def water_core_result(water_solver_no_coupling):
    """Just the solved core object."""
    return water_solver_no_coupling


@pytest.fixture(scope='module')
def water_pes_no_coupling(water_solver_no_coupling):
    """PES scan without coupling (coarse grid for speed)."""
    R_grid = np.array([1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0])
    return water_solver_no_coupling.scan_pes(R_grid=R_grid, n_Re=200)


# ---------------------------------------------------------------------------
# 1. Core solves correctly
# ---------------------------------------------------------------------------

class TestCore:
    def test_core_solves(self, water_core_result):
        """Core energy should be finite and negative."""
        assert water_core_result.E_core is not None
        assert water_core_result.E_core < 0

    def test_z_eff(self, water_core_result):
        """Z_eff = Z_O - n_core = 8 - 2 = 6."""
        assert water_core_result.Z_eff == 6.0

    def test_pk_derived(self, water_core_result):
        """Ab initio PK should be derived from core."""
        assert water_core_result.ab_initio_pk is not None
        assert water_core_result.pk_A > 0
        assert water_core_result.pk_B > 0


# ---------------------------------------------------------------------------
# 2. Bond pair solves at R = 1.81
# ---------------------------------------------------------------------------

class TestBondPair:
    def test_bond_solves(self, water_core_result):
        """Bond pair at R_eq should return finite energy."""
        E_bond = water_core_result._solve_bond_at_R(R=1.81, n_Re=200)
        assert np.isfinite(E_bond)
        assert E_bond < 0  # bound state

    def test_bond_returns_full_result(self, water_core_result):
        """Full result dict should contain expected keys."""
        result = water_core_result._solve_bond_at_R(
            R=1.81, n_Re=200, return_full=True)
        assert 'E_elec' in result
        assert 'wavefunction' in result
        assert np.isfinite(result['E_elec'])


# ---------------------------------------------------------------------------
# 3. Lone pair solves
# ---------------------------------------------------------------------------

class TestLonePair:
    def test_lone_pair_solves(self, water_core_result):
        """Lone pair energy should be finite and negative."""
        assert water_core_result.E_lone_pair is not None
        assert water_core_result.E_lone_pair < 0

    def test_lone_pair_channel_data(self, water_core_result):
        """Channel data should be extracted for coupling."""
        cd = water_core_result._lone_pair_channel_data
        assert cd is not None
        assert 'channels' in cd
        assert 'ch_weights' in cd
        assert len(cd['channels']) > 0

    def test_lone_pair_energy_reasonable(self, water_core_result):
        """Lone pair energy for Z_eff=6 should be more negative than He (Z=2)."""
        # He ground state ~ -2.90 Ha. Z_eff=6 should be much lower.
        assert water_core_result.E_lone_pair < -20.0


# ---------------------------------------------------------------------------
# 4. PES scan produces bound state
# ---------------------------------------------------------------------------

class TestPESScan:
    def test_pes_has_minimum(self, water_pes_no_coupling):
        """PES should have a well-defined minimum (bound state)."""
        E_valid = np.array(water_pes_no_coupling['E_valid'])
        R_valid = np.array(water_pes_no_coupling['R_valid'])
        # Minimum should not be at the endpoints
        i_min = np.argmin(E_valid)
        assert i_min > 0, "Minimum at left edge — PES not bound"
        assert i_min < len(E_valid) - 1, "Minimum at right edge — PES not bound"

    def test_d_e_positive(self, water_pes_no_coupling):
        """Dissociation energy should be positive (bound)."""
        assert water_pes_no_coupling['D_e'] > 0

    def test_r_eq_physical(self, water_pes_no_coupling):
        """R_eq should be in the range 1.0–3.0 bohr."""
        R_eq = water_pes_no_coupling['R_eq']
        assert 1.0 < R_eq < 3.0, f"R_eq={R_eq} outside physical range"


# ---------------------------------------------------------------------------
# 5. R_eq in reasonable range
# ---------------------------------------------------------------------------

class TestREquilibrium:
    def test_r_eq_within_range(self, water_pes_no_coupling):
        """R_eq should be within 1.0–4.0 bohr (expt: 1.81, expect ~1.4 uncoupled with charge_center)."""
        R_eq = water_pes_no_coupling['R_eq']
        assert 1.0 <= R_eq <= 4.0, f"R_eq={R_eq:.3f} outside 1.0–4.0 range"


# ---------------------------------------------------------------------------
# 6. Block-diagonal PES (no coupling) runs without error
# ---------------------------------------------------------------------------

class TestNoCoupling:
    def test_no_coupling_runs(self, water_pes_no_coupling):
        """PES scan without coupling should complete without error."""
        assert water_pes_no_coupling is not None
        assert len(water_pes_no_coupling['R_valid']) >= 3


# ---------------------------------------------------------------------------
# 7. Exchange-coupled PES runs without error
# ---------------------------------------------------------------------------

class TestWithCoupling:
    @pytest.mark.slow
    def test_coupled_pes_runs(self):
        """PES scan with inter-fiber coupling should complete."""
        solver = ComposedWaterSolver(
            l_max=2, n_alpha=100,
            include_coupling=True,
            verbose=False,
        )
        solver.solve_core()
        solver.solve_lone_pair()
        R_grid = np.array([2.0, 2.5, 3.0, 3.5])
        pes = solver.scan_pes(R_grid=R_grid, n_Re=200)
        assert pes is not None
        assert len(pes['R_valid']) >= 3

    @pytest.mark.slow
    def test_coupling_affects_energy(self):
        """Exchange coupling should modify total energy."""
        solver_off = ComposedWaterSolver(
            l_max=2, n_alpha=100,
            include_coupling=False,
            verbose=False,
        )
        solver_off.solve_core()
        solver_off.solve_lone_pair()

        solver_on = ComposedWaterSolver(
            l_max=2, n_alpha=100,
            include_coupling=True,
            verbose=False,
        )
        solver_on.solve_core()
        solver_on.solve_lone_pair()

        R_grid = np.array([2.0, 2.5, 3.0])
        pes_off = solver_off.scan_pes(R_grid=R_grid, n_Re=200)
        pes_on = solver_on.scan_pes(R_grid=R_grid, n_Re=200)

        E_off = np.array(pes_off['E_valid'])
        E_on = np.array(pes_on['E_valid'])

        # Energies should differ when coupling is on
        assert not np.allclose(E_off, E_on, atol=1e-6), \
            "Coupling had no effect on energies"


# ---------------------------------------------------------------------------
# 8. Nuclear repulsion geometry
# ---------------------------------------------------------------------------

class TestNuclearRepulsion:
    def test_v_nn_at_known_geometry(self):
        """V_NN at R=1.809, θ=104.5° should match analytic formula."""
        solver = ComposedWaterSolver(verbose=False)
        R = 1.809
        theta = 1.8238
        d_HH = 2.0 * R * np.sin(theta / 2.0)
        # V_NN = 2*Z_O*Z_H/R + Z_H^2/d_HH
        V_expected = 2.0 * 8.0 * 1.0 / R + 1.0 / d_HH
        V_computed = solver._nuclear_repulsion(R)
        np.testing.assert_allclose(V_computed, V_expected, rtol=1e-6)

    def test_v_nn_decreases_with_R(self):
        """Nuclear repulsion should decrease as R increases."""
        solver = ComposedWaterSolver(verbose=False)
        V1 = solver._nuclear_repulsion(1.5)
        V2 = solver._nuclear_repulsion(2.5)
        assert V1 > V2


# ---------------------------------------------------------------------------
# 9. Constructor validation
# ---------------------------------------------------------------------------

class TestConstructor:
    def test_default_angles(self):
        """Default angles should be set correctly."""
        solver = ComposedWaterSolver(verbose=False)
        assert abs(solver.theta_HOH - 1.8238) < 0.001
        assert abs(solver.Z_eff - 6.0) < 1e-10

    def test_custom_angles(self):
        """Custom angles should be accepted."""
        solver = ComposedWaterSolver(
            theta_HOH=np.pi/2,
            theta_bond_lone=2.0,
            theta_lone_lone=1.8,
            verbose=False,
        )
        assert abs(solver.theta_HOH - np.pi/2) < 1e-10
        assert abs(solver.theta_bond_lone - 2.0) < 1e-10
        assert abs(solver.theta_lone_lone - 1.8) < 1e-10

    def test_invalid_n_core(self):
        """n_core != 2 should raise ValueError."""
        with pytest.raises(ValueError):
            ComposedWaterSolver(n_core=4, verbose=False)

    def test_invalid_pk_mode(self):
        """Invalid pk_mode should raise ValueError."""
        with pytest.raises(ValueError):
            ComposedWaterSolver(pk_mode='invalid', verbose=False)


# ---------------------------------------------------------------------------
# 10. Charge-center origin
# ---------------------------------------------------------------------------

class TestChargeCenter:
    def test_bond_uses_charge_center(self, water_core_result):
        """Bond pair solver should use origin='charge_center'."""
        result = water_core_result._solve_bond_at_R(
            R=1.81, n_Re=200, return_full=True)
        assert result.get('origin') == 'charge_center', \
            f"Expected origin='charge_center', got '{result.get('origin')}'."

    def test_charge_center_z0_formula(self, water_core_result):
        """z0 should equal R*(Z_A-Z_B)/(2*(Z_A+Z_B)) for charge_center."""
        R = 1.81
        result = water_core_result._solve_bond_at_R(
            R=R, n_Re=200, return_full=True)
        Z_A = water_core_result.Z_eff  # 6
        Z_B = water_core_result.Z_H    # 1
        z0_expected = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
        np.testing.assert_allclose(
            result['z0'], z0_expected, rtol=1e-10,
            err_msg="z0 formula mismatch")
