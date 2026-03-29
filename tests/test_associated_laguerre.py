"""
Tests for associated Laguerre moment matrices, lowered moment matrix,
and Stieltjes integral matrix (Track J infrastructure).

Validates:
1. Band structure and symmetry of all matrices
2. Agreement with Gauss-Laguerre quadrature reference
3. Reduction to ordinary Laguerre case at alpha=0
4. Stability diagnostics across matrix dimensions
"""
import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from scipy.special import roots_genlaguerre, eval_genlaguerre, gammaln
from geovac.prolate_spheroidal_lattice import (
    _laguerre_moment_matrices,
    _associated_laguerre_moment_matrices,
    _lowered_moment_matrix,
    _stieltjes_matrix,
)


def _quadrature_moment(N: int, alpha: int, k: int, n_quad: int = 200) -> np.ndarray:
    """Compute moment matrix M_k via Gauss-Laguerre quadrature (reference)."""
    x_q, w_q = roots_genlaguerre(n_quad, alpha)
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_genlaguerre(n, alpha, x_q)
    xk = x_q ** k
    wL = w_q[np.newaxis, :] * L_vals * xk[np.newaxis, :]
    return wL @ L_vals.T


def _quadrature_lowered(N: int, alpha: int, n_quad: int = 200) -> np.ndarray:
    """Compute lowered moment matrix via quadrature with weight x^(alpha-1) e^{-x}."""
    x_q, w_q = roots_genlaguerre(n_quad, alpha - 1)
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_genlaguerre(n, alpha, x_q)
    wL = w_q[np.newaxis, :] * L_vals
    return wL @ L_vals.T


def _quadrature_stieltjes(N: int, alpha: int, a: float, n_quad: int = 300) -> np.ndarray:
    """Compute Stieltjes matrix via quadrature."""
    x_q, w_q = roots_genlaguerre(n_quad, alpha)
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_genlaguerre(n, alpha, x_q)
    inv_xa = 1.0 / (x_q + a)
    wL = w_q[np.newaxis, :] * L_vals * inv_xa[np.newaxis, :]
    return wL @ L_vals.T


# ============================================================
# Phase 1: Associated Laguerre Moment Matrices
# ============================================================

class TestAssociatedMomentMatrices:
    """Validate M0, M1, M2 for associated Laguerre polynomials."""

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_m0_diagonal(self, alpha: int):
        """M0 should be diagonal with h_n entries."""
        N = 10
        M0, _, _ = _associated_laguerre_moment_matrices(N, alpha)
        off_diag = M0 - np.diag(np.diag(M0))
        assert np.max(np.abs(off_diag)) < 1e-15

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_m0_values(self, alpha: int):
        """M0 diagonal should equal h_n = Gamma(n+alpha+1)/n!."""
        N = 10
        M0, _, _ = _associated_laguerre_moment_matrices(N, alpha)
        for n in range(N):
            h_n = np.exp(gammaln(n + alpha + 1) - gammaln(n + 1))
            assert abs(M0[n, n] - h_n) < 1e-12 * h_n

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_m1_tridiagonal(self, alpha: int):
        """M1 should be tridiagonal."""
        N = 10
        _, M1, _ = _associated_laguerre_moment_matrices(N, alpha)
        for i in range(N):
            for j in range(N):
                if abs(i - j) > 1:
                    assert abs(M1[i, j]) < 1e-15

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_m2_pentadiagonal(self, alpha: int):
        """M2 should be pentadiagonal."""
        N = 10
        _, _, M2 = _associated_laguerre_moment_matrices(N, alpha)
        for i in range(N):
            for j in range(N):
                if abs(i - j) > 2:
                    assert abs(M2[i, j]) < 1e-15

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_symmetry(self, alpha: int):
        """All moment matrices should be symmetric."""
        N = 15
        M0, M1, M2 = _associated_laguerre_moment_matrices(N, alpha)
        assert np.max(np.abs(M0 - M0.T)) < 1e-14
        assert np.max(np.abs(M1 - M1.T)) < 1e-14
        assert np.max(np.abs(M2 - M2.T)) < 1e-14

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_vs_quadrature_m0(self, alpha: int):
        """M0 should match quadrature reference."""
        N = 15
        M0, _, _ = _associated_laguerre_moment_matrices(N, alpha)
        M0_ref = _quadrature_moment(N, alpha, 0)
        rel_err = np.max(np.abs(M0 - M0_ref)) / np.max(np.abs(M0_ref))
        assert rel_err < 1e-12, f"M0 rel error {rel_err:.2e}"

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_vs_quadrature_m1(self, alpha: int):
        """M1 should match quadrature reference."""
        N = 15
        _, M1, _ = _associated_laguerre_moment_matrices(N, alpha)
        M1_ref = _quadrature_moment(N, alpha, 1)
        rel_err = np.max(np.abs(M1 - M1_ref)) / np.max(np.abs(M1_ref))
        assert rel_err < 1e-12, f"M1 rel error {rel_err:.2e}"

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_vs_quadrature_m2(self, alpha: int):
        """M2 should match quadrature reference."""
        N = 15
        _, _, M2 = _associated_laguerre_moment_matrices(N, alpha)
        M2_ref = _quadrature_moment(N, alpha, 2)
        rel_err = np.max(np.abs(M2 - M2_ref)) / np.max(np.abs(M2_ref))
        assert rel_err < 1e-12, f"M2 rel error {rel_err:.2e}"

    def test_alpha_zero_matches_ordinary(self):
        """At alpha=0, associated matrices should equal ordinary Laguerre."""
        N = 15
        M0, M1, M2 = _associated_laguerre_moment_matrices(N, 0)
        M0_ord, M1_ord, M2_ord = _laguerre_moment_matrices(N)
        assert np.max(np.abs(M0 - M0_ord)) < 1e-14
        assert np.max(np.abs(M1 - M1_ord)) < 1e-14
        assert np.max(np.abs(M2 - M2_ord)) < 1e-14


# ============================================================
# Phase 2: Lowered Moment Matrix
# ============================================================

class TestLoweredMoment:
    """Validate M_{-1} lowered moment matrix."""

    def test_alpha_zero_raises(self):
        """alpha=0 should raise ValueError (integral diverges)."""
        with pytest.raises(ValueError):
            _lowered_moment_matrix(10, 0)

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_symmetry(self, alpha: int):
        """M_{-1} should be symmetric."""
        N = 15
        M = _lowered_moment_matrix(N, alpha)
        assert np.max(np.abs(M - M.T)) < 1e-14

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_positive_definite(self, alpha: int):
        """M_{-1} should be positive definite."""
        N = 15
        M = _lowered_moment_matrix(N, alpha)
        evals = np.linalg.eigvalsh(M)
        assert np.all(evals > 0), f"Negative eigenvalue found: {evals.min()}"

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_vs_quadrature(self, alpha: int):
        """M_{-1} should match quadrature reference."""
        N = 15
        M = _lowered_moment_matrix(N, alpha)
        M_ref = _quadrature_lowered(N, alpha)
        rel_err = np.max(np.abs(M - M_ref)) / np.max(np.abs(M_ref))
        assert rel_err < 1e-12, f"M_{{-1}} rel error {rel_err:.2e}"

    @pytest.mark.parametrize("alpha", [1, 2, 3])
    def test_structure(self, alpha: int):
        """M_{-1}[i,j] = sum_{p=0}^{min(i,j)} h_p^(alpha-1) — monotonically increasing along rows."""
        N = 10
        M = _lowered_moment_matrix(N, alpha)
        # Each row i: M[i,0] <= M[i,1] <= ... <= M[i,i] = M[i,i+1] = ... = M[i,N-1]
        for i in range(N):
            for j in range(min(i, N - 1)):
                assert M[i, j] <= M[i, j + 1] + 1e-14
            # Constant for j >= i
            for j in range(i, N - 1):
                assert abs(M[i, j] - M[i, j + 1]) < 1e-12


# ============================================================
# Phase 3: Stieltjes Integral Matrix
# ============================================================

class TestStieltjesMatrix:
    """Validate Stieltjes integral matrix J."""

    @pytest.mark.parametrize("alpha,a", [
        (1, 1.0), (1, 4.0), (1, 10.0),
        (2, 1.0), (2, 4.0), (2, 10.0),
        (3, 1.0), (3, 4.0), (3, 10.0),
    ])
    def test_symmetry(self, alpha: int, a: float):
        """J should be symmetric."""
        N = 15
        J = _stieltjes_matrix(N, alpha, a)
        asym = np.max(np.abs(J - J.T))
        assert asym < 1e-10, f"Asymmetry {asym:.2e} at alpha={alpha}, a={a}"

    @pytest.mark.parametrize("alpha,a", [
        (1, 1.0), (1, 4.0), (1, 10.0),
        (2, 1.0), (2, 4.0), (2, 10.0),
        (3, 1.0), (3, 4.0), (3, 10.0),
    ])
    def test_positive_definite(self, alpha: int, a: float):
        """J should be positive definite."""
        N = 15
        J = _stieltjes_matrix(N, alpha, a)
        evals = np.linalg.eigvalsh(J)
        assert np.all(evals > 0), f"Negative eigenvalue {evals.min():.6e}"

    @pytest.mark.parametrize("alpha,a", [
        (1, 4.0), (1, 10.0),
        (2, 4.0), (2, 10.0),
        (3, 4.0), (3, 10.0),
    ])
    def test_vs_quadrature_large_a(self, alpha: int, a: float):
        """J should match quadrature to machine precision for a >= 4."""
        N = 20
        J = _stieltjes_matrix(N, alpha, a)
        J_ref = _quadrature_stieltjes(N, alpha, a)
        max_elem = np.max(np.abs(J_ref))
        rel_err = np.max(np.abs(J - J_ref)) / max_elem
        assert rel_err < 1e-10, (
            f"Stieltjes rel error {rel_err:.2e} at alpha={alpha}, a={a}"
        )

    @pytest.mark.parametrize("alpha,a", [
        (1, 1.0), (2, 1.0), (3, 1.0),
    ])
    def test_vs_quadrature_small_a(self, alpha: int, a: float):
        """For small a, forward recurrence limits accuracy to ~1e-8."""
        N = 15
        J = _stieltjes_matrix(N, alpha, a)
        J_ref = _quadrature_stieltjes(N, alpha, a)
        max_elem = np.max(np.abs(J_ref))
        rel_err = np.max(np.abs(J - J_ref)) / max_elem
        assert rel_err < 1e-7, (
            f"Stieltjes rel error {rel_err:.2e} at alpha={alpha}, a={a}"
        )

    def test_a_negative_raises(self):
        """Negative shift parameter should raise ValueError."""
        with pytest.raises(ValueError):
            _stieltjes_matrix(10, 1, -1.0)

    def test_a_zero_raises(self):
        """Zero shift parameter should raise ValueError."""
        with pytest.raises(ValueError):
            _stieltjes_matrix(10, 1, 0.0)

    @pytest.mark.parametrize("alpha", [0, 1, 2])
    def test_large_a_limit(self, alpha: int):
        """For large a, J ~ M0/a - M1/a^2 (two-term expansion of 1/(x+a))."""
        N = 10
        a = 1000.0
        J = _stieltjes_matrix(N, alpha, a)
        M0, M1, _ = _associated_laguerre_moment_matrices(N, alpha)
        for n in range(N):
            # Two-term expansion: J[n,n] ~ M0[n,n]/a - M1[n,n]/a^2
            expected = M0[n, n] / a - M1[n, n] / (a * a)
            assert abs(J[n, n] - expected) / max(abs(expected), 1e-30) < 0.01, (
                f"Large-a limit failed at n={n}, alpha={alpha}"
            )


# ============================================================
# Stability Diagnostics
# ============================================================

class TestStabilityDiagnostics:
    """Track precision as a function of matrix dimension N."""

    def test_moment_stability(self):
        """Moment matrices should maintain precision up to N=30."""
        alpha = 2
        for N in [5, 10, 15, 20, 25, 30]:
            M0, M1, M2 = _associated_laguerre_moment_matrices(N, alpha)
            M0_ref = _quadrature_moment(N, alpha, 0, n_quad=300)
            M1_ref = _quadrature_moment(N, alpha, 1, n_quad=300)
            M2_ref = _quadrature_moment(N, alpha, 2, n_quad=300)
            err0 = np.max(np.abs(M0 - M0_ref)) / np.max(np.abs(M0_ref))
            err1 = np.max(np.abs(M1 - M1_ref)) / np.max(np.abs(M1_ref))
            err2 = np.max(np.abs(M2 - M2_ref)) / np.max(np.abs(M2_ref))
            assert err0 < 1e-11, f"M0 precision loss at N={N}: {err0:.2e}"
            assert err1 < 1e-11, f"M1 precision loss at N={N}: {err1:.2e}"
            assert err2 < 1e-11, f"M2 precision loss at N={N}: {err2:.2e}"

    def test_lowered_stability(self):
        """Lowered moment matrix should maintain precision up to N=30."""
        alpha = 2
        for N in [5, 10, 15, 20, 25, 30]:
            M = _lowered_moment_matrix(N, alpha)
            M_ref = _quadrature_lowered(N, alpha, n_quad=300)
            rel_err = np.max(np.abs(M - M_ref)) / np.max(np.abs(M_ref))
            assert rel_err < 1e-11, f"M_{{-1}} precision loss at N={N}: {rel_err:.2e}"

    @pytest.mark.parametrize("alpha,a", [(1, 4.0), (2, 4.0), (3, 4.0)])
    def test_stieltjes_stability_production(self, alpha: int, a: float):
        """Stieltjes at production parameters (a >= 4) should be stable to N=30."""
        for N in [5, 10, 15, 20, 25, 30]:
            J = _stieltjes_matrix(N, alpha, a)
            J_ref = _quadrature_stieltjes(N, alpha, a, n_quad=300)
            max_elem = np.max(np.abs(J_ref))
            rel_err = np.max(np.abs(J - J_ref)) / max_elem
            assert rel_err < 1e-10, (
                f"Stieltjes precision loss at N={N}, alpha={alpha}, a={a}: {rel_err:.2e}"
            )


# ============================================================
# Alpha=0 Stieltjes (ordinary Laguerre)
# ============================================================

class TestStieltjesAlphaZero:
    """Stieltjes matrix should work for alpha=0 too."""

    @pytest.mark.parametrize("a", [1.0, 4.0, 10.0])
    def test_vs_quadrature(self, a: float):
        """Alpha=0 Stieltjes should match quadrature."""
        N = 15
        J = _stieltjes_matrix(N, 0, a)
        J_ref = _quadrature_stieltjes(N, 0, a)
        max_elem = np.max(np.abs(J_ref))
        rel_err = np.max(np.abs(J - J_ref)) / max_elem
        # Threshold depends on a: tighter for large a (backward segment dominates)
        tol = 1e-8 if a < 2 else 1e-10
        assert rel_err < tol, f"Alpha=0 Stieltjes rel error {rel_err:.2e} at a={a}"
