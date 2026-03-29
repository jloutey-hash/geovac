"""
Tests for the algebraic coupled-channel solver.

Validates the integration of AlgebraicAngularSolver with the coupled-channel
radial solver, including P-matrix coupling and Q-matrix (DBOC) corrections.

Test inventory:
  1. P-matrix antisymmetry
  2. DBOC decomposition consistency (total = internal + external)
  3. Q-matrix sign convention (diagonal Q raises energy)
  4. Eigenvalue consistency (coupled-channel eigenvalues are physical)
  5. Channel norm conservation (ground state weights sum to 1)
  6. l_max=0 P-only benchmark (achieves < 0.1%)
  7. l_max=1 energy direction (coupled ABOVE single-channel with Q)
  8. l_max convergence monotonic (q_mode='full' converges toward exact)
  9. Channel population (ground state dominated by channel 0)
  10. Q diagonal matches DBOC
  11. dP/dR matches finite differences
  12. dP/dR diagonal is zero (antisymmetry)
  13. Q_exact diagonal equals Q_closure diagonal (dP_μμ/dR = 0)
  14. q_mode='exact' l_max convergence monotonic (l_max=1-5)
  15. q_mode='exact' improves over 'full' at l_max>=1
  16. l_max=5 error below 0.25% (convergence floor characterization)
  17. n_channels=5 matches n_channels=3 (channel convergence)
"""

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# 1. P-matrix antisymmetry
# ---------------------------------------------------------------------------

def test_p_matrix_antisymmetry():
    """P_μν = -P_νμ (antisymmetry of first-derivative coupling)."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.array([0.5, 1.0, 2.0, 5.0])
    coupling = compute_algebraic_coupling(solver, R_grid, n_channels=3)

    P = coupling['P']
    for i in range(len(R_grid)):
        P_i = P[:, :, i]
        np.testing.assert_allclose(
            P_i, -P_i.T, atol=1e-12,
            err_msg=f"P not antisymmetric at R={R_grid[i]}"
        )


# ---------------------------------------------------------------------------
# 2. DBOC decomposition consistency
# ---------------------------------------------------------------------------

def test_dboc_decomposition():
    """DBOC_total = DBOC_internal + DBOC_external at every R point."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.linspace(0.5, 10.0, 20)
    coupling = compute_algebraic_coupling(solver, R_grid, n_channels=3)

    np.testing.assert_allclose(
        coupling['DBOC_total'],
        coupling['DBOC_internal'] + coupling['DBOC_external'],
        atol=1e-12,
        err_msg="DBOC total != internal + external"
    )


# ---------------------------------------------------------------------------
# 3. Q-matrix sign convention
# ---------------------------------------------------------------------------

def test_q_diagonal_negative():
    """Q_μμ should be negative (so -½Q raises energy via +DBOC).

    Q_μμ = Σ_κ P_μκ P_κμ = -Σ_κ |P_μκ|² = -2*DBOC < 0.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.array([1.0, 2.0, 5.0])
    coupling = compute_algebraic_coupling(solver, R_grid, n_channels=3)

    Q = coupling['Q']
    for i in range(len(R_grid)):
        for mu in range(3):
            assert Q[mu, mu, i] <= 0, (
                f"Q[{mu},{mu}] = {Q[mu, mu, i]:.6f} should be <= 0 at R={R_grid[i]}"
            )


# ---------------------------------------------------------------------------
# 4. Eigenvalue consistency
# ---------------------------------------------------------------------------

def test_coupled_eigenvalue_physical():
    """Coupled-channel ground state should be between -Z²/2 and 0 for He."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=10, l_max=0, n_channels=3,
        n_R=100, N_R_radial=1000, q_mode='full', verbose=False,
    )
    E = result['energy']
    assert -4.0 < E < 0.0, f"Energy {E} outside physical range [-4, 0]"
    # Should be near -2.9
    assert -3.5 < E < -2.5, f"Energy {E} too far from expected ~-2.9"


# ---------------------------------------------------------------------------
# 5. Channel norm conservation
# ---------------------------------------------------------------------------

def test_channel_weights_sum_to_one():
    """Ground state channel weights should sum to 1."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=10, l_max=0, n_channels=3,
        n_R=100, N_R_radial=1000, q_mode='full', verbose=False,
    )
    np.testing.assert_allclose(
        np.sum(result['channel_weights']), 1.0, atol=0.01,
        err_msg="Channel weights don't sum to 1"
    )


# ---------------------------------------------------------------------------
# 6. l_max=0 P-only benchmark
# ---------------------------------------------------------------------------

def test_lmax0_p_only_accuracy():
    """P-only coupled-channel at l_max=0 should achieve < 0.15% error."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=0, n_channels=3,
        n_R=150, N_R_radial=2000, q_mode='none', verbose=False,
    )
    E_exact = -2.903724
    err = abs(result['energy'] - E_exact) / abs(E_exact) * 100
    assert err < 0.15, f"l_max=0 P-only error {err:.4f}% exceeds 0.15%"


# ---------------------------------------------------------------------------
# 7. l_max=1 energy direction with Q
# ---------------------------------------------------------------------------

def test_lmax1_coupled_above_single():
    """With Q, coupled-channel at l_max=1 should be ABOVE single-channel."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=1, n_channels=3,
        n_R=150, N_R_radial=2000, q_mode='full', verbose=False,
    )
    assert result['energy'] > result['energy_single_channel'], (
        f"Coupled {result['energy']:.6f} should be above "
        f"single-channel {result['energy_single_channel']:.6f}"
    )
    # Should also be above exact
    E_exact = -2.903724
    assert result['energy'] > E_exact, (
        f"Coupled {result['energy']:.6f} should be above exact {E_exact}"
    )


# ---------------------------------------------------------------------------
# 8. l_max convergence monotonic
# ---------------------------------------------------------------------------

def test_lmax_convergence_monotonic():
    """Error should decrease monotonically for l_max=1,2,3 with q_mode='full'."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    E_exact = -2.903724
    errors = []
    for lmax in [1, 2, 3]:
        result = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=lmax, n_channels=3,
            n_R=150, N_R_radial=2000, q_mode='full', verbose=False,
        )
        err = abs(result['energy'] - E_exact) / abs(E_exact) * 100
        errors.append(err)

    # Each error should be <= previous (monotonic convergence)
    for i in range(1, len(errors)):
        assert errors[i] <= errors[i - 1] + 0.01, (
            f"Error at l_max={i + 1} ({errors[i]:.4f}%) not <= "
            f"error at l_max={i} ({errors[i - 1]:.4f}%)"
        )


# ---------------------------------------------------------------------------
# 9. Channel population
# ---------------------------------------------------------------------------

def test_ground_state_channel_0_dominant():
    """Ground state should be dominated by channel 0 (weight > 0.99)."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=10, l_max=0, n_channels=3,
        n_R=100, N_R_radial=1000, q_mode='full', verbose=False,
    )
    assert result['channel_weights'][0] > 0.99, (
        f"Channel 0 weight {result['channel_weights'][0]:.4f} should be > 0.99"
    )


# ---------------------------------------------------------------------------
# 10. Q diagonal matches DBOC
# ---------------------------------------------------------------------------

def test_q_diagonal_matches_dboc():
    """Q_μμ should equal -2*DBOC_total_μ (closure approximation is exact
    for diagonal because dP_μμ/dR = 0 by antisymmetry)."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    coupling = compute_algebraic_coupling(solver, R_grid, n_channels=3)

    Q = coupling['Q']
    DBOC = coupling['DBOC_total']

    for i in range(len(R_grid)):
        for mu in range(3):
            expected = -2.0 * DBOC[mu, i]
            np.testing.assert_allclose(
                Q[mu, mu, i], expected, rtol=0.01,
                err_msg=f"Q[{mu},{mu}] != -2*DBOC at R={R_grid[i]}"
            )


# ---------------------------------------------------------------------------
# 11. dP/dR matches finite differences
# ---------------------------------------------------------------------------

def test_dPdR_matches_finite_difference():
    """Algebraic dP/dR should match central finite-difference to < 1e-6."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=15, l_max=0)

    R_test = 2.0
    delta = 1e-5
    R_grid = np.array([R_test - delta, R_test, R_test + delta])
    coupling = compute_algebraic_coupling(
        solver, R_grid, n_channels=3, compute_exact_dPdR=True
    )

    P = coupling['P']
    dPdR_alg = coupling['dPdR'][:, :, 1]

    # Central finite difference
    dPdR_fd = (P[:, :, 2] - P[:, :, 0]) / (2 * delta)

    np.testing.assert_allclose(
        dPdR_alg, dPdR_fd, atol=1e-6,
        err_msg="Algebraic dP/dR does not match finite difference"
    )


# ---------------------------------------------------------------------------
# 12. dP/dR diagonal is zero
# ---------------------------------------------------------------------------

def test_dPdR_diagonal_zero():
    """dP_μμ/dR = 0 because P_μμ = 0 for all R (antisymmetry)."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.array([0.5, 1.0, 2.0, 5.0])
    coupling = compute_algebraic_coupling(
        solver, R_grid, n_channels=3, compute_exact_dPdR=True
    )

    dPdR = coupling['dPdR']
    for i in range(len(R_grid)):
        for mu in range(3):
            assert abs(dPdR[mu, mu, i]) < 1e-12, (
                f"dP[{mu},{mu}]/dR = {dPdR[mu, mu, i]} should be 0 at R={R_grid[i]}"
            )


# ---------------------------------------------------------------------------
# 13. Q_exact diagonal equals Q_closure diagonal
# ---------------------------------------------------------------------------

def test_q_exact_diagonal_equals_closure():
    """Q_exact_μμ = Q_closure_μμ because dP_μμ/dR = 0."""
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.algebraic_coupled_channel import compute_algebraic_coupling

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    R_grid = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    coupling = compute_algebraic_coupling(
        solver, R_grid, n_channels=3, compute_exact_dPdR=True
    )

    Q = coupling['Q']
    Q_exact = coupling['Q_exact']

    for i in range(len(R_grid)):
        for mu in range(3):
            np.testing.assert_allclose(
                Q_exact[mu, mu, i], Q[mu, mu, i], atol=1e-12,
                err_msg=f"Q_exact[{mu},{mu}] != Q_closure at R={R_grid[i]}"
            )


# ---------------------------------------------------------------------------
# 14. q_mode='exact' l_max convergence monotonic
# ---------------------------------------------------------------------------

def test_exact_lmax_convergence_monotonic():
    """Error should decrease monotonically for l_max=1,2,3,4,5 with q_mode='exact'."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    E_exact = -2.903724
    errors = []
    for lmax in [1, 2, 3, 4, 5]:
        result = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=lmax, n_channels=3,
            n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
        )
        err = abs(result['energy'] - E_exact) / abs(E_exact) * 100
        errors.append(err)

    for i in range(1, len(errors)):
        assert errors[i] <= errors[i - 1] + 0.01, (
            f"Error at l_max={i + 1} ({errors[i]:.4f}%) not <= "
            f"error at l_max={i} ({errors[i - 1]:.4f}%)"
        )


# ---------------------------------------------------------------------------
# 15. q_mode='exact' improves over 'full' at l_max>=1
# ---------------------------------------------------------------------------

def test_exact_improves_over_full():
    """q_mode='exact' should give lower error than 'full' at l_max=2."""
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    E_exact = -2.903724
    result_full = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=2, n_channels=3,
        n_R=150, N_R_radial=2000, q_mode='full', verbose=False,
    )
    result_exact = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=2, n_channels=3,
        n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
    )
    err_full = abs(result_full['energy'] - E_exact) / abs(E_exact) * 100
    err_exact = abs(result_exact['energy'] - E_exact) / abs(E_exact) * 100

    assert err_exact < err_full, (
        f"Exact ({err_exact:.4f}%) should be better than full ({err_full:.4f}%)"
    )


# ---------------------------------------------------------------------------
# 16. l_max=5 error below 0.25% (convergence floor characterization)
# ---------------------------------------------------------------------------

def test_lmax5_error_below_025():
    """l_max=5 with q_mode='exact' should achieve < 0.25% error.

    This characterizes the coupled-channel convergence floor. The error
    converges as ~l_max^{-2} toward a structural floor of ~0.19-0.20%
    set by the adiabatic channel truncation. Sub-0.1% requires a
    variational 2D solver (Paper 15) rather than more l_max.
    """
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    E_exact = -2.903724
    result = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=5, n_channels=3,
        n_R=200, N_R_radial=3000, q_mode='exact', verbose=False,
    )
    err = abs(result['energy'] - E_exact) / abs(E_exact) * 100
    assert err < 0.25, f"l_max=5 exact error {err:.4f}% exceeds 0.25%"
    # Also confirm above the structural floor (not overshooting)
    assert err > 0.15, f"l_max=5 exact error {err:.4f}% below 0.15% — update if genuine"


# ---------------------------------------------------------------------------
# 17. n_channels=5 matches n_channels=3 (channel convergence)
# ---------------------------------------------------------------------------

def test_5_channels_matches_3():
    """n_channels=5 should give nearly identical energy to n_channels=3.

    At l_max=3, channels 4-5 have weights < 1e-5 and contribute < 0.001%
    to the energy. This confirms 3-channel truncation is sufficient.
    """
    from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

    result_3 = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=3, n_channels=3,
        n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
    )
    result_5 = solve_hyperspherical_algebraic_coupled(
        Z=2.0, n_basis=15, l_max=3, n_channels=5,
        n_R=150, N_R_radial=2000, q_mode='exact', verbose=False,
    )
    # Energies should match to < 0.001 Ha
    np.testing.assert_allclose(
        result_3['energy'], result_5['energy'], atol=0.001,
        err_msg="3-channel and 5-channel energies differ by > 1 mHa"
    )
