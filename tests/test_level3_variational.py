"""
Tests for the Level 3 (He) 2D variational solver on S^5.

Track DI, Sprint 1: Break the adiabatic 0.19-0.20% floor by treating
hyperradius R and hyperangle alpha simultaneously in a tensor-product
spectral basis.

Key results validated:
  - 2D solver is variational (energy above exact)
  - Breaks the adiabatic floor at l_max >= 2
  - l_max convergence is monotonic
  - Radial basis converges quickly (n_R >= 20 sufficient)
  - Cusp correction at l_max=4 achieves < 0.01% error
"""

import numpy as np
import pytest

from geovac.level3_variational import (
    solve_he_variational_2d,
    _build_radial_operator_matrices,
)

E_EXACT = -2.903724377034119598  # Pekeris (Ha)


class TestRadialOperators:
    """Validate the radial operator matrices."""

    def test_overlap_positive_definite(self) -> None:
        S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(20, 1.5, 0.05)
        eigs = np.linalg.eigvalsh(S_R)
        assert np.all(eigs > 0), "Overlap matrix must be positive definite"

    def test_kinetic_positive_semidefinite(self) -> None:
        S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(20, 1.5, 0.05)
        eigs = np.linalg.eigvalsh(K_R)
        assert np.all(eigs >= -1e-12), "Kinetic matrix must be positive semidefinite"

    def test_inv_R_positive_definite(self) -> None:
        S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(20, 1.5, 0.05)
        eigs = np.linalg.eigvalsh(M_inv_R)
        assert np.all(eigs > 0), "1/R operator matrix must be positive definite"

    def test_inv_R2_positive_definite(self) -> None:
        S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(20, 1.5, 0.05)
        eigs = np.linalg.eigvalsh(M_inv_R2)
        assert np.all(eigs > 0), "1/R^2 operator matrix must be positive definite"

    def test_matrices_symmetric(self) -> None:
        S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(20, 1.5, 0.05)
        for name, M in [("S_R", S_R), ("K_R", K_R), ("M_inv_R", M_inv_R),
                         ("M_inv_R2", M_inv_R2)]:
            assert np.allclose(M, M.T, atol=1e-14), f"{name} must be symmetric"


class TestVariational2D:
    """Core validation of the 2D variational solver."""

    def test_variational_bound(self) -> None:
        """Energy must be above exact (variational principle)."""
        res = solve_he_variational_2d(
            n_basis_R=20, n_basis_alpha=15, l_max=3, alpha_R=1.5,
        )
        assert res['energies'][0] > E_EXACT, \
            f"Variational violation: {res['energies'][0]:.8f} < {E_EXACT:.8f}"

    def test_lmax0_basic(self) -> None:
        """l_max=0 should give ~0.9% error (s-wave only)."""
        res = solve_he_variational_2d(
            n_basis_R=25, n_basis_alpha=15, l_max=0, alpha_R=1.5,
        )
        assert res['error_pct'] < 1.0, f"l_max=0 error too large: {res['error_pct']:.2f}%"
        assert res['error_pct'] > 0.5, f"l_max=0 error suspiciously small: {res['error_pct']:.2f}%"

    def test_lmax_monotonic_convergence(self) -> None:
        """Energy must decrease monotonically with l_max."""
        energies = []
        for lm in range(5):
            res = solve_he_variational_2d(
                n_basis_R=20, n_basis_alpha=10, l_max=lm, alpha_R=1.5,
            )
            energies.append(res['energies'][0])

        for i in range(1, len(energies)):
            assert energies[i] < energies[i - 1], \
                f"Energy not decreasing: l_max={i} ({energies[i]:.6f}) >= l_max={i-1} ({energies[i-1]:.6f})"

    def test_breaks_adiabatic_floor(self) -> None:
        """2D solver at l_max=5 must beat the adiabatic 0.19-0.20% floor.

        The adiabatic coupled-channel solver (v2.0.8) converges to a
        structural floor of 0.19-0.20% because it neglects non-adiabatic
        coupling. The 2D solver eliminates this approximation.
        """
        res = solve_he_variational_2d(
            n_basis_R=25, n_basis_alpha=20, l_max=5, alpha_R=2.0,
        )
        assert res['error_pct'] < 0.10, \
            f"Failed to break adiabatic floor: {res['error_pct']:.4f}% (target < 0.10%)"

    def test_radial_convergence(self) -> None:
        """Energy should be stable for n_R >= 20."""
        E_prev = None
        for n_R in [20, 25, 30]:
            res = solve_he_variational_2d(
                n_basis_R=n_R, n_basis_alpha=15, l_max=2, alpha_R=1.5,
            )
            if E_prev is not None:
                assert abs(res['energies'][0] - E_prev) < 1e-5, \
                    f"Radial not converged at n_R={n_R}"
            E_prev = res['energies'][0]

    def test_angular_convergence(self) -> None:
        """Increasing angular basis should lower energy."""
        E_prev = None
        for n_a in [10, 20, 30]:
            res = solve_he_variational_2d(
                n_basis_R=25, n_basis_alpha=n_a, l_max=0, alpha_R=1.5,
            )
            if E_prev is not None:
                assert res['energies'][0] <= E_prev + 1e-10, \
                    f"Angular convergence violated at n_alpha={n_a}"
            E_prev = res['energies'][0]

    def test_dimension_consistency(self) -> None:
        """Check that reported dimensions are correct."""
        res = solve_he_variational_2d(
            n_basis_R=15, n_basis_alpha=10, l_max=3, alpha_R=1.5,
        )
        assert res['dim_R'] == 15
        assert res['dim_alpha'] == 10 * 4  # 4 channels (l=0,1,2,3)
        assert res['dim_total'] == 15 * 40

    def test_multiple_states(self) -> None:
        """Requesting multiple states should give ordered eigenvalues."""
        res = solve_he_variational_2d(
            n_basis_R=20, n_basis_alpha=10, l_max=0, alpha_R=1.5,
            n_states=3,
        )
        assert len(res['energies']) == 3
        assert res['energies'][0] < res['energies'][1] < res['energies'][2]


class TestCuspCorrection:
    """Test cusp-corrected accuracy targets."""

    def test_cusp_correction_sub_01pct(self) -> None:
        """l_max=4 + cusp correction should achieve < 0.01% error.

        This is the headline result of Track DI Phase 1.
        """
        from geovac.cusp_correction import cusp_correction_he

        res = solve_he_variational_2d(
            n_basis_R=25, n_basis_alpha=40, l_max=4, alpha_R=2.0,
            n_quad_alpha=150,
        )
        dE_cusp = cusp_correction_he(Z=2.0, l_max=4)
        E_corrected = res['energies'][0] + dE_cusp
        error_pct = abs((E_corrected - E_EXACT) / E_EXACT) * 100.0

        assert error_pct < 0.01, \
            f"Cusp-corrected error {error_pct:.5f}% exceeds 0.01% target"

    def test_raw_lmax5_sub_005pct(self) -> None:
        """Raw 2D solver at l_max=5 with large basis should be < 0.05%."""
        res = solve_he_variational_2d(
            n_basis_R=25, n_basis_alpha=30, l_max=5, alpha_R=2.0,
        )
        assert res['error_pct'] < 0.05, \
            f"Raw l_max=5 error {res['error_pct']:.4f}% exceeds 0.05% target"


@pytest.mark.slow
class TestPrecision:
    """Extended precision tests (computationally expensive)."""

    def test_precision_solver(self) -> None:
        """Full precision solver pipeline."""
        from geovac.level3_variational import solve_he_precision

        result = solve_he_precision(
            n_basis_R=25, n_basis_alpha=40,
            l_max_range=list(range(6)),
            alpha_R=2.0,
        )

        assert result['best_error_pct'] < 0.01, \
            f"Best error {result['best_error_pct']:.5f}% exceeds 0.01%"

        # Raw energies should be monotonically decreasing with l_max
        lmax_vals = sorted(result['raw_energies'].keys())
        for i in range(1, len(lmax_vals)):
            E_curr = result['raw_energies'][lmax_vals[i]]
            E_prev = result['raw_energies'][lmax_vals[i - 1]]
            assert E_curr < E_prev, "Non-monotonic l_max convergence"

    def test_lmax7_convergence(self) -> None:
        """l_max=7 with large angular basis: raw < 0.03%."""
        res = solve_he_variational_2d(
            n_basis_R=25, n_basis_alpha=40, l_max=7, alpha_R=2.0,
            n_quad_alpha=150,
        )
        assert res['error_pct'] < 0.03, \
            f"l_max=7 error {res['error_pct']:.4f}% exceeds 0.03%"
