"""Tests for 4-electron coupled-channel solver (Track AO).

Tests the coupled-channel radial solver for the N-electron mol-frame
hyperspherical system. Validates P-matrix computation, DBOC diagnostics,
and the sanity check that 1-channel reproduces adiabatic.

Track AO result is a NEGATIVE RESULT: coupled-channel does not fix the
overbinding because FD P-matrices are too noisy and eigenvalue gaps
are too small for reliable non-adiabatic coupling.
"""

import numpy as np
import pytest
from geovac.n_electron_solver import (
    compute_multichannel_angular_sweep_4e,
    compute_p_matrix_fd_4e,
    solve_4e_lih_coupled,
    CENTRIFUGAL_4E,
)


# ---------------------------------------------------------------------------
# P-matrix diagnostics
# ---------------------------------------------------------------------------

class TestPMatrixDiagnostics:
    """Test P-matrix computation from finite differences."""

    @pytest.fixture(scope='class')
    def p_data_lmax2(self):
        """Compute P-matrix at l_max=2 for LiH at R=1.0."""
        R = 1.0
        R_e_grid = np.concatenate([
            np.linspace(0.5, 1.5, 8),
            np.linspace(1.5, 4.0, 10),
            np.linspace(4.0, 8.0, 5),
        ])
        R_e_grid = np.unique(R_e_grid)

        angular_data = compute_multichannel_angular_sweep_4e(
            R=R, R_e_grid=R_e_grid,
            Z_A=3.0, Z_B=1.0, z0=0.0,
            n_grid=6, l_max=2, n_channels=3,
            symmetry='s4', verbose=False,
        )
        p_data = compute_p_matrix_fd_4e(angular_data, n_channels=3)
        return angular_data, p_data

    def test_p_matrix_shape(self, p_data_lmax2):
        """P-matrix has correct shape."""
        _, p_data = p_data_lmax2
        P = p_data['P']
        assert P.shape[0] == 3  # n_channels
        assert P.shape[1] == 3
        assert P.shape[2] > 10  # n_Re points

    def test_p_diagonal_zero(self, p_data_lmax2):
        """Diagonal P-matrix elements should be ~0 (norm conservation)."""
        _, p_data = p_data_lmax2
        P = p_data['P']
        n_Re = P.shape[2]
        # Check interior points (central difference gives exact 0 on diagonal)
        for i in range(1, n_Re - 1):
            for mu in range(3):
                assert abs(P[mu, mu, i]) < 1e-10, \
                    f"P[{mu},{mu},{i}] = {P[mu, mu, i]:.2e} (expected ~0)"

    def test_p_antisymmetric_noisy(self, p_data_lmax2):
        """FD P-matrix at Level 4N has poor antisymmetry.

        This is part of the Track AO negative result: the FD P-matrix is too
        noisy for reliable coupled-channel corrections. The eigenvalue gaps are
        O(0.01-0.05), amplifying FD errors. We verify that SOME antisymmetry
        exists (not completely random), but accept the high violation rate
        as a documented structural finding.
        """
        _, p_data = p_data_lmax2
        P = p_data['P']
        n_Re = P.shape[2]
        P_max = np.max(np.abs(P))
        atol = 0.3 * P_max  # Generous tolerance for FD noise
        n_violations = 0
        n_checked = 0
        for i in range(2, n_Re - 2):
            for nu in range(3):
                for mu in range(nu + 1, 3):
                    if max(abs(P[nu, mu, i]), abs(P[mu, nu, i])) > 0.05 * P_max:
                        diff = abs(P[nu, mu, i] + P[mu, nu, i])
                        n_checked += 1
                        if diff > atol:
                            n_violations += 1
        # High violation rate is expected due to small eigenvalue gaps
        if n_checked > 0:
            frac = n_violations / n_checked
            assert frac < 0.8, \
                f"P should show SOME antisymmetry: " \
                f"{n_violations}/{n_checked} ({frac*100:.0f}%) violations"

    def test_p_frob_norm_finite(self, p_data_lmax2):
        """Frobenius norm of P should be finite and positive."""
        _, p_data = p_data_lmax2
        assert np.all(np.isfinite(p_data['P_frob_norm']))
        assert np.max(p_data['P_frob_norm']) > 0

    def test_dboc_positive(self, p_data_lmax2):
        """DBOC should be non-negative (it is a sum of squares)."""
        _, p_data = p_data_lmax2
        DBOC = p_data['DBOC']
        assert np.all(DBOC >= -1e-15)

    def test_p_much_larger_than_level3(self, p_data_lmax2):
        """P-matrix at Level 4N should be much larger than Level 3.

        This documents the structural finding: eigenvalue gaps at Level 4N
        are O(0.01-0.05), making P ~ V/gap very large. At Level 3, gaps are
        O(1) and P peaks at ~0.3-0.5. The 10-30x amplification is the root
        cause of the coupled-channel failure.
        """
        _, p_data = p_data_lmax2
        # Level 3 P peak is ~0.3-0.5
        # Level 4N should be > 1.0 (we observe ~3.3)
        assert np.max(p_data['P_frob_norm']) > 1.0, \
            "P-matrix should be significantly larger than Level 3"

    def test_dboc_much_larger_than_level3(self, p_data_lmax2):
        """DBOC at Level 4N should be much larger than Level 3.

        Level 3 DBOC peaks at ~0.05-0.15 Ha. Level 4N DBOC peaks at ~5 Ha.
        """
        _, p_data = p_data_lmax2
        assert np.max(p_data['DBOC'][0]) > 1.0, \
            "DBOC should be significantly larger than Level 3"


# ---------------------------------------------------------------------------
# Coupled-channel solver
# ---------------------------------------------------------------------------

class TestCoupledChannelSolver:
    """Test the coupled-channel radial solver.

    Track AO negative result: coupled-channel overshoots below adiabatic
    at most R-values due to large FD P-matrices with small eigenvalue gaps.
    Tests validate the infrastructure works correctly even though the
    physics result is negative.
    """

    @pytest.mark.slow
    def test_single_channel_reproduces_adiabatic(self):
        """1-channel coupled should reproduce adiabatic result."""
        R = 1.0
        result = solve_4e_lih_coupled(
            R, n_grid=6, l_max=2, n_channels=1,
            q_mode='none', verbose=False,
        )
        diff = abs(result['E_total'] - result['E_total_single'])
        assert diff < 0.01, \
            f"1-channel coupled ({result['E_total']:.6f}) != " \
            f"adiabatic ({result['E_total_single']:.6f}), diff={diff:.6f}"

    @pytest.mark.slow
    def test_channel_weights_sum_to_one(self):
        """Channel weights should sum to 1."""
        R = 1.0
        result = solve_4e_lih_coupled(
            R, n_grid=6, l_max=2, n_channels=3,
            q_mode='none', verbose=False,
        )
        assert abs(np.sum(result['channel_weights']) - 1.0) < 0.01

    @pytest.mark.slow
    def test_coupled_energy_finite(self):
        """Coupled-channel energy should be finite."""
        R = 1.0
        result = solve_4e_lih_coupled(
            R, n_grid=6, l_max=2, n_channels=3,
            q_mode='diagonal', verbose=False,
        )
        assert np.isfinite(result['E_total'])
        assert np.isfinite(result['E_total_single'])

    @pytest.mark.slow
    def test_negative_result_documented(self):
        """Coupled-channel overshoots below adiabatic (negative result).

        At 3+ channels, the P-coupling without adequate Q causes the
        energy to overshoot below the adiabatic energy (wrong direction).
        This is the structural negative result of Track AO.
        """
        R = 1.0
        result_none = solve_4e_lih_coupled(
            R, n_grid=6, l_max=2, n_channels=3,
            q_mode='none', verbose=False,
        )
        # With P-coupling only (no Q), energy should go DOWN (overshoot)
        delta = result_none['E_total'] - result_none['E_total_single']
        # We observe delta ~ -0.24 Ha at R=1.0 (overshoot below adiabatic)
        assert delta < 0.1, \
            f"Expected overshoot or near-zero, got delta = {delta:+.6f}"
