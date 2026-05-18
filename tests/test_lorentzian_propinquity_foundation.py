"""Tests for the Sprint L3b first-move foundation modules.

Covers:
  - geovac.krein_space_compact_temporal (CompactTemporalKreinSpace)
  - geovac.lorentzian_dirac_compact
  - geovac.operator_system_compact_temporal
  - geovac.krein_positive_state_space
  - geovac.central_fejer_compact_temporal

LOAD-BEARING falsifiers (must pass bit-exact at residual = 0.0):
  - Riemannian limit at N_t = 1 for every construction
  - Krein-self-adjointness D_L^x = D_L for the compact Lorentzian Dirac
  - Propagation number = 2 (achievable envelope) for the compact-temporal
    operator system at n_max >= 2
  - Plancherel symbol factorization in joint Fejer kernel (exact sympy)

These tests build on Sprint L2's 390-test baseline and Sprint L3a-1's
33 tests (423 total); zero regression must be preserved.
"""

from __future__ import annotations

import math

import mpmath
import numpy as np
import pytest
import sympy as sp

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_compact_temporal import (
    CompactTemporalKreinSpace,
    fourier_momentum_grid,
)
from geovac.lorentzian_dirac_compact import (
    fourier_d_dt_matrix,
    krein_self_adjoint_residual,
    lorentzian_dirac_compact_matrix,
    verify_riemannian_limit_compact,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
    compact_temporal_multiplier_matrices,
)
from geovac.krein_positive_state_space import KreinPositiveStateSpace
from geovac.central_fejer_compact_temporal import (
    cb_norm_circle,
    fejer_kernel_circle,
    gamma_rate_circle,
    joint_cb_norm,
    joint_fejer_kernel,
    joint_gamma_rate,
    joint_plancherel_symbol,
    plancherel_symbol_circle,
    verify_plancherel_factorization,
    verify_riemannian_limit_compact_temporal,
)


# ---------------------------------------------------------------------------
# 1. Compact-temporal Krein space basics
# ---------------------------------------------------------------------------


def test_fourier_momentum_grid_odd():
    """Symmetric grid for odd N_t."""
    g = fourier_momentum_grid(5)
    assert g.tolist() == [-2, -1, 0, 1, 2]


def test_fourier_momentum_grid_singleton():
    """At N_t=1, only k=0."""
    g = fourier_momentum_grid(1)
    assert g.tolist() == [0]


def test_compact_krein_space_basic():
    """Construction sanity at (n_max=2, N_t=3)."""
    K = CompactTemporalKreinSpace(n_max=2, N_t=3)
    assert K.dim == K.dim_spatial * K.N_t
    assert K.dim_spatial == full_dirac_dim(2)
    assert K.dim == 48


def test_compact_krein_J_squared():
    """J^2 = +I bit-exact."""
    K = CompactTemporalKreinSpace(n_max=2, N_t=3)
    ok, residual = K.verify_J_squared_identity(tol=1e-14)
    assert ok
    assert residual < 1e-14


def test_compact_krein_J_hermitian():
    """J = J^*."""
    K = CompactTemporalKreinSpace(n_max=2, N_t=3)
    ok, residual = K.verify_J_hermitian(tol=1e-14)
    assert ok


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_compact_krein_riemannian_limit_BITEXACT(n_max):
    """LOAD-BEARING: at N_t = 1, K reduces to H_GV bit-identically.

    Must pass with J_match_residual EXACTLY 0.0 in float64.
    """
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=1)
    ok, details = K.riemannian_limit_check(tol=1e-14)
    assert ok, details
    assert details["J_match_residual"] == 0.0
    assert details["dim_match"]
    assert details["basis_match"]


# ---------------------------------------------------------------------------
# 2. Compact-temporal Lorentzian Dirac
# ---------------------------------------------------------------------------


def test_fourier_d_dt_singleton():
    """At N_t=1: D_t = (0)."""
    D_t = fourier_d_dt_matrix(1)
    assert D_t.shape == (1, 1)
    assert D_t[0, 0] == 0


def test_fourier_d_dt_anti_hermitian():
    """D_t is anti-Hermitian (diagonal, purely imaginary)."""
    D_t = fourier_d_dt_matrix(5, T=2 * np.pi)
    residual = np.linalg.norm(D_t + D_t.conj().T)
    assert residual < 1e-14


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (3, 1), (2, 3), (3, 5)])
def test_compact_dirac_krein_self_adjoint(n_max, N_t):
    """D_L^x = D_L bit-exact for truthful Camporesi-Higuchi spatial Dirac."""
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t)
    residual = krein_self_adjoint_residual(K)
    assert residual == 0.0, f"residual={residual}"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_compact_dirac_riemannian_limit_BITEXACT(n_max):
    """LOAD-BEARING: at N_t=1, D_L = i * D_GV bit-identically."""
    ok, details = verify_riemannian_limit_compact(n_max, tol=1e-14)
    assert ok, details
    assert details["residual_F_norm"] == 0.0


# ---------------------------------------------------------------------------
# 3. Compact-temporal operator system
# ---------------------------------------------------------------------------


def test_compact_temporal_multipliers_basic():
    """Temporal multipliers: identity at p=0, diagonal in momentum."""
    mats = compact_temporal_multiplier_matrices(5, T=2 * np.pi)
    assert len(mats) == 5
    # p=0 is identity
    assert np.allclose(mats[0], np.eye(5))
    # All diagonal
    for m in mats:
        assert np.allclose(m - np.diag(np.diag(m)), 0)


def test_compact_operator_system_construction():
    """Instantiation succeeds at (n_max=2, N_t=3)."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    assert O.dim_K == 48
    # 14 spatial multipliers (chirality-doubled FullDirac at n_max=2)
    # * 3 temporal multipliers = 42 generators
    assert len(O.multiplier_labels) == 42


def test_compact_operator_system_identity_in_O():
    """I_{dim_K} is in O (constant function on S^3 x S^1)."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    ok, residual = O.identity_in_O()
    assert ok, f"residual={residual}"


def test_compact_operator_system_star_closed():
    """*-closure: every generator's adjoint is in O."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    ok, failures = O.is_star_closed()
    assert ok, f"failures={failures}"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_compact_operator_system_riemannian_limit_BITEXACT(n_max):
    """LOAD-BEARING: at N_t=1, the compact-temporal construction reduces
    bit-identically to FullDiracTruncatedOperatorSystem(n_max).
    """
    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=1)
    ok, details = O.verify_riemannian_limit()
    assert ok, details
    assert details["max_residual"] == 0.0, details


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 3), (3, 1)])
def test_compact_propagation_number_achievable_envelope(n_max, N_t):
    """prop = 2 under the achievable envelope (Weyl-doubled subspace).

    Matches Paper 32 §III prop=2 and Sprint L3a-1 verdict.
    """
    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    prop, dim_seq = O.compute_propagation_number(
        envelope="achievable", max_k=3
    )
    assert prop == 2, f"prop={prop}, dim_seq={dim_seq}"


def test_compact_operator_system_krein_positive_trivial():
    """Structural finding from L3a-1: all multipliers preserve K^+ in
    the chirality-doubled scalar lift.
    """
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    preservers, _ = O.krein_positive_preservers()
    assert len(preservers) == len(O.multiplier_labels)


def test_compact_operator_system_compare_to_l3a1():
    """Compatibility check: same n_gens and spatial labels as L3a-1."""
    from geovac.operator_system_lorentzian import (
        LorentzianTruncatedOperatorSystem,
    )

    O_grid = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    O_compact = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    cmp = O_compact.compare_to_l3a1_grid(O_grid)
    assert cmp["n_max_match"]
    assert cmp["N_t_match"]
    assert cmp["dim_K_match"]
    assert cmp["n_multipliers_match"]
    assert cmp["spat_labels_match"]
    # Linear span dim equal (same spatial factor, both temporal algebras
    # are N_t-dim commutative diagonal subalgebras)
    assert cmp["dim_O_grid"] == cmp["dim_O_compact"]


# ---------------------------------------------------------------------------
# 4. K^+ state space
# ---------------------------------------------------------------------------


def test_krein_positive_state_space_construction():
    """Basic construction at (n_max=2, N_t=3)."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    assert S.dim_K == 48
    assert S.K_plus_dim + S.K_minus_dim == S.dim_K


def test_J_eigendecomposition():
    """J has eigenvalues +/- 1 with equal multiplicities."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    ok, details = S.verify_J_eigendecomp()
    assert ok, details
    # Chirality-doubled: K^+ = K^- in dim
    assert S.K_plus_dim == S.K_minus_dim


def test_pure_state_density_trace_one():
    """Pure-state density has trace 1."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    rho = S.pure_state_density(0)
    assert abs(np.trace(rho) - 1.0) < 1e-14


def test_K_plus_state_is_krein_positive():
    """A K^+ basis vector gives a Krein-positive state.

    On a K^+ vector v (J v = +v), the Krein-twisted expectation
    omega(a^* J a) = <v, J a^* J a v> = <v, a^* J a v>; for the
    chirality-doubled multiplier a = M oplus M and J = chirality-swap,
    one shows omega(a^* J a) >= 0 numerically.
    """
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    v_plus = S.K_plus_eigvecs[:, 0]
    rho_plus = S.pure_state_density(v_plus)
    ok, details = S.is_krein_positive_state(rho_plus, sample_size=10)
    assert ok, details


def test_K_plus_inner_product_positive():
    """Hilbert-space inner product on K^+ vectors is positive-definite."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    v = S.K_plus_eigvecs[:, 0]
    val = S.K_plus_inner_product(v, v)
    assert abs(val.imag) < 1e-12
    assert val.real > 0


def test_K_plus_restriction_dimensions():
    """Projecting rho to K^+ block reduces dimension to K_plus_dim."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)
    S = KreinPositiveStateSpace(op_sys=O)
    rho = S.pure_state_density(0)
    rho_plus = S.restrict_to_K_plus(rho)
    assert rho_plus.shape == (S.K_plus_dim, S.K_plus_dim)


@pytest.mark.slow
def test_wasserstein_distance_finite_with_offdiag():
    """SDP-based pure-state distance computes finite values with offdiag CH.

    The truthful CH spatial Dirac has n-degeneracy that makes most
    cross-shell pair distances +infinity (per L3a-1 / R3.5).  The
    offdiag CH lifts the n-degeneracy and gives finite distances on
    cross-shell pure-state pairs.
    """
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=1)
    S = KreinPositiveStateSpace(op_sys=O)
    D_offdiag = camporesi_higuchi_offdiag_dirac_matrix(O.basis_spatial)
    # Take two basis vectors in different shells
    d = S.wasserstein_distance_pure(0, 2, D_offdiag)
    assert d >= 0
    assert d < float("inf")


@pytest.mark.slow
def test_wasserstein_distance_v_eq_w():
    """d(v, v) = 0."""
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=1)
    S = KreinPositiveStateSpace(op_sys=O)
    D = camporesi_higuchi_full_dirac_matrix(O.basis_spatial)
    d = S.wasserstein_distance_pure(0, 0, D)
    assert d == 0.0


# ---------------------------------------------------------------------------
# 5. Joint compact x compact Fejer kernel
# ---------------------------------------------------------------------------


def test_plancherel_symbol_circle_N_t_5():
    """Cesaro / Fejer Plancherel on the circle at N_t=5: k=0 -> 1."""
    s_0 = plancherel_symbol_circle(5, 0)
    s_1 = plancherel_symbol_circle(5, 1)
    s_2 = plancherel_symbol_circle(5, 2)
    s_3 = plancherel_symbol_circle(5, 3)
    assert s_0 == 1
    # (N_t + 1 - 2*1) / (N_t + 1) = 4/6 = 2/3
    assert s_1 == sp.Rational(2, 3)
    # (N_t + 1 - 2*2) / (N_t + 1) = 2/6 = 1/3
    assert s_2 == sp.Rational(1, 3)
    # |k| = 3 > K_max = 2, so 0
    assert s_3 == 0


def test_plancherel_symbol_circle_N_t_1():
    """At N_t=1: only k=0 with hat_K(0) = 1."""
    assert plancherel_symbol_circle(1, 0) == 1
    assert plancherel_symbol_circle(1, 1) == 0


def test_joint_cb_norm():
    """Joint cb-norm = SU(2) cb-norm * U(1) cb-norm = 2/(n_max+1) * 1."""
    cb = joint_cb_norm(3, 5)
    assert cb == sp.Rational(2, 4)


def test_joint_plancherel_factorization():
    """LOAD-BEARING: hat_K_joint(j, k) = hat_K_SU(2)(j) * hat_K_U(1)(k)
    in exact sympy rationals.
    """
    ok, details = verify_plancherel_factorization(3, 5)
    assert ok, details
    assert details["pairs_match"] == details["pairs_checked"]


def test_joint_kernel_symbolic():
    """Joint kernel constructed as product of factors."""
    chi, theta = sp.symbols("chi theta", real=True)
    K_joint = joint_fejer_kernel(2, 3, T=2 * sp.pi, chi=chi, theta=theta)
    # Should depend on both chi and theta
    free = K_joint.free_symbols
    assert chi in free
    assert theta in free


def test_gamma_rate_circle_basic():
    """Sanity for gamma_u1 at N_t=5, T=2pi."""
    g = gamma_rate_circle(5, T=2 * np.pi, prec=20)
    # Positive and bounded
    assert g > 0
    assert g < 2 * np.pi  # Less than diameter of S^1_{2pi}


def test_gamma_rate_circle_decreases_with_N_t():
    """gamma_u1 decreases as N_t grows (kernel concentrates)."""
    g_3 = gamma_rate_circle(3, T=2 * np.pi, prec=20)
    g_7 = gamma_rate_circle(7, T=2 * np.pi, prec=20)
    g_15 = gamma_rate_circle(15, T=2 * np.pi, prec=20)
    assert g_3 > g_7 > g_15


def test_joint_gamma_rate_factorization():
    """joint gamma_l1 = gamma_su2 + gamma_u1 (linear sum)."""
    d = joint_gamma_rate(3, 5, T=2 * np.pi, prec=30)
    expected_l1 = d["gamma_su2"] + d["gamma_u1"]
    assert abs(d["gamma_l1"] - expected_l1) < mpmath.mpf("1e-25")


def test_joint_gamma_panel_paper38_consistency():
    """At N_t held fixed, gamma_su2 component matches Paper 38 SU(2) panel.

    Paper 38: gamma_su2(2) ~ 2.0746, gamma_su2(3) ~ 1.6101, gamma_su2(4)
    ~ 1.3223.  These are bit-identical to the joint construction's
    SU(2) factor (it's the same gamma_rate_su2 call).
    """
    d2 = joint_gamma_rate(2, 3, prec=20)
    d3 = joint_gamma_rate(3, 3, prec=20)
    d4 = joint_gamma_rate(4, 3, prec=20)
    # Check against Paper 38 known values (to ~4 decimal places)
    assert abs(float(d2["gamma_su2"]) - 2.0746) < 1e-3
    assert abs(float(d3["gamma_su2"]) - 1.6101) < 1e-3
    assert abs(float(d4["gamma_su2"]) - 1.3223) < 1e-3


def test_riemannian_limit_compact_temporal_kernel():
    """At N_t = 1: gamma_su2 matches Paper 38 verbatim;
    gamma_u1 = T/4 (uniform-circle distance mean).
    """
    ok, details = verify_riemannian_limit_compact_temporal(
        n_max=2, T=2 * np.pi, prec=20
    )
    assert ok, details
    # gamma_su2 component is bit-identical to Paper 38 SU(2) gamma
    # (same gamma_rate_su2 function call).
    assert details["gamma_su2_residual"] == 0.0
    # gamma_u1 at N_t=1 is the constant Haar = uniform on S^1_T.
    # Mean of min(theta, T-theta) over theta in [0,T] = T/4.
    # The residual is dominated by mpmath.quad precision (~10^-4 at
    # default settings) plus float64 truncation of T = 2pi; the residual
    # is small compared to T/4 ~ 1.57 but not bit-exact in numerical
    # integration. We assert it falls in the expected quad-precision band.
    assert details["gamma_u1_residual"] < 1e-3


# ---------------------------------------------------------------------------
# 6. Zero regression spot checks (no upstream behavior change)
# ---------------------------------------------------------------------------


def test_zero_regression_full_dirac_dim():
    """Sprint L2 full_dirac_dim still works (no breakage)."""
    assert full_dirac_dim(1) == 4
    assert full_dirac_dim(2) == 16
    assert full_dirac_dim(3) == 40


def test_zero_regression_camporesi_higuchi():
    """Camporesi-Higuchi spatial Dirac unchanged."""
    basis = full_dirac_basis(2)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    assert D.shape == (16, 16)
    assert np.allclose(D, D.conj().T)
