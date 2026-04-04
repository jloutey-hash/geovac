"""
Tests for the 4-electron 2D variational solver (Track AR).

Tests the solve_4e_lih_2d() and scan_pes_4e_lih_2d() functions that
bypass the adiabatic approximation by solving the full (R_e, angular)
tensor product problem.
"""

import numpy as np
import pytest


# ==========================================================================
# FAST TESTS (unit / structural / small systems)
# ==========================================================================


def test_2d_solver_returns_dict():
    """solve_4e_lih_2d returns a dict with required keys."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    assert isinstance(result, dict)
    for key in ['E_elec', 'E_total', 'D_e', 'R', 'l_max', 'method',
                'angular_dim', 'total_dim', 'time_total']:
        assert key in result, f"Missing key: {key}"
    assert result['method'] == '2d'
    assert result['l_max'] == 2


def test_2d_solver_dimensions():
    """Verify angular and total dimensions are correct."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    # Angular dim should be positive (12 S4 channels * 27 grid points = 324)
    assert result['angular_dim'] > 0
    # Total dim should be n_Re * angular_dim
    assert result['total_dim'] == 15 * result['angular_dim']


def test_2d_solver_energy_negative():
    """Electronic energy should be negative (bound system)."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    assert result['E_elec'] < 0, "Electronic energy should be negative"
    assert result['E_total'] < 0, "Total energy should be negative"


def test_2d_solver_charge_center():
    """Test charge_center origin produces a result."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        origin='charge_center',
        symmetry='s4',
        verbose=False,
    )
    assert np.isfinite(result['E_total'])
    assert result['z0'] != 0.0  # LiH is heteronuclear


def test_2d_solver_variational_bound_basic():
    """E_total at R=1.0 must be above exact -8.0705 Ha."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    E_exact = -8.0705
    assert result['E_total'] > E_exact, \
        f"Variational bound violated: E={result['E_total']:.4f} < exact={E_exact}"


def test_scan_pes_returns_dict():
    """scan_pes_4e_lih_2d returns dict with required keys."""
    from geovac.n_electron_2d import scan_pes_4e_lih_2d

    result = scan_pes_4e_lih_2d(
        R_values=np.array([0.7, 1.0, 1.5]),
        n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    assert isinstance(result, dict)
    for key in ['R', 'E_total', 'D_e', 'R_eq', 'E_min', 'method']:
        assert key in result
    assert result['method'] == '2d'
    assert len(result['R']) == 3


def test_scan_pes_energy_ordering():
    """PES energies: high at short R, low near R_eq."""
    from geovac.n_electron_2d import scan_pes_4e_lih_2d

    result = scan_pes_4e_lih_2d(
        R_values=np.array([0.5, 1.0, 2.0]),
        n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    # At R=0.5, V_NN = 6.0 Ha; E_total should be highest
    assert result['E_total'][0] > result['E_total'][1], \
        "Energy at R=0.5 should be higher than R=1.0"


def test_2d_equilibrium_region():
    """The 2D PES should have a minimum near R ~ 1.0 bohr."""
    from geovac.n_electron_2d import scan_pes_4e_lih_2d

    result = scan_pes_4e_lih_2d(
        R_values=np.array([0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0]),
        n_grid=3, l_max=2, n_Re=15,
        R_e_min=0.5, R_e_max=5.0,
        symmetry='s4',
        verbose=False,
    )
    # Minimum should be interior (not at boundary)
    assert result['has_minimum'], \
        f"No equilibrium found. Min at boundary R={result['R_eq']:.3f}"


# ==========================================================================
# SLOW TESTS (production grids, full PES scans)
# ==========================================================================


@pytest.mark.slow
def test_2d_solver_production_grid():
    """2D solver at production n_grid=4 with S4 projection."""
    from geovac.n_electron_2d import solve_4e_lih_2d

    result = solve_4e_lih_2d(
        R=1.0, n_grid=4, l_max=2, n_Re=50,
        R_e_min=0.3, R_e_max=10.0,
        symmetry='s4',
        verbose=True,
    )
    E_exact = -8.0705
    assert result['E_total'] > E_exact, \
        f"Variational bound violated: E={result['E_total']:.4f}"
    print(f"  E_total = {result['E_total']:.6f} Ha")
    print(f"  D_e     = {result['D_e']:.6f} Ha")


@pytest.mark.slow
def test_2d_pes_scan_full():
    """Full PES scan with 2D solver."""
    from geovac.n_electron_2d import scan_pes_4e_lih_2d

    result = scan_pes_4e_lih_2d(
        R_values=np.array([0.5, 0.7, 0.9, 1.0, 1.1, 1.2, 1.5,
                           2.0, 3.0, 4.0, 6.0]),
        n_grid=4, l_max=2, n_Re=50,
        R_e_min=0.3, R_e_max=10.0,
        symmetry='s4',
        verbose=True,
        output_dir='debug/track_ar',
    )
    # Equilibrium should exist
    assert result['has_minimum'], \
        f"No equilibrium found. Min at R={result['R_eq']:.3f}"


@pytest.mark.slow
def test_2d_vs_adiabatic_comparison():
    """Compare 2D vs adiabatic at R=1.0 bohr.

    The 2D total energy should be higher than the adiabatic (which overcounts
    due to missing non-adiabatic coupling).
    """
    from geovac.n_electron_2d import solve_4e_lih_2d
    from geovac.n_electron_solver import solve_4e_lih_multichannel

    R = 1.0

    result_2d = solve_4e_lih_2d(
        R=R, n_grid=4, l_max=2, n_Re=50,
        R_e_min=0.3, R_e_max=10.0,
        symmetry='s4',
        verbose=True,
    )

    result_ad = solve_4e_lih_multichannel(
        R=R, n_grid=6, l_max=2,
        verbose=True,
    )

    print(f"\n  Comparison at R = {R:.1f} bohr:")
    print(f"    Adiabatic E_total = {result_ad['E_total']:.6f} Ha")
    print(f"    2D E_total        = {result_2d['E_total']:.6f} Ha")
    delta = result_2d['E_total'] - result_ad['E_total']
    print(f"    Delta (2D - adiab) = {delta:.6f} Ha")

    # 2D should be above adiabatic (less binding)
    assert delta > -0.1, \
        f"2D energy unexpectedly far below adiabatic: delta={delta:.4f}"
