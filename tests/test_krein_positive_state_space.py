"""Tests for `geovac.krein_positive_state_space` (Sprint L3b first move).

These tests verify:

  (a) Instantiation of `KreinPositiveStateSpace` on a compact-temporal
      truncated operator system at small panel cells (n_max in {1, 2},
      N_t in {1, 3}).

  (b) J eigendecomposition: chirality doubling
        dim(+1) = dim(-1) = dim_K / 2
      exact at every tested cell (J is a 0/+/-1 permutation; eigenvalues
      are exactly +/- 1 to floating-point precision).

  (c) Projectors P_+ and P_- are orthogonal idempotents with
      P_+ + P_- = I, P_+ P_- = 0.

  (d) Pure-state densities from K^+ eigenvectors are Hermitian, trace-1,
      rank-1, and satisfy the structural Krein-positivity trace marker
      Tr(rho J) = +1 (K^+ pure states have J-positive expectation).

  (e) K^- pure states have Tr(rho J) = -1 (the "antimatter" half), which
      is the structural distinguishing trace marker.

  (f) The existing module's `is_krein_positive_state` check passes on
      BOTH K^+ and K^- pure-state densities (this matches the L3a-1
      "trivial Krein-positivity at operator-multiplier level" finding;
      every multiplier commutes with J, so omega(M^* J M) = Tr(rho M^* M)
      >= 0 holds for any density matrix). The state-space-level
      distinguishing observable is the J-trace marker Tr(rho J), NOT the
      generator-level inequality.

  (g) Wasserstein-Kantorovich SDP infrastructure at small scale:
        d(idx, idx) = 0 (sanity, no SDP needed -- direct identity check)
        d(idx_v, idx_w) finite at (n_max=2, N_t=1) with offdiag CH spatial
                          Dirac
      symmetry d(v, w) = d(w, v) by SDP feasible-set symmetry argument
      (Hermitian feasible region is invariant under x -> -x).

  (h) Riemannian limit at N_t = 1 (LOAD-BEARING analog): the K^+ state
      space dimension matches the spatial-only (chirality-doubled
      FullDirac) setup; specifically dim_K = dim_spatial,
      dim_K_plus = dim_spatial / 2.

If any of (b), (h) fails at N_t = 1, STOP and escalate (the load-bearing
structure has changed). The other tests can be re-examined under panel
context.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_dim,
)
from geovac.krein_positive_state_space import KreinPositiveStateSpace
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def op_sys_2_1():
    """(n_max=2, N_t=1) -- smallest non-trivial cell, dim_K = 16."""
    return CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=1)


@pytest.fixture
def op_sys_2_3():
    """(n_max=2, N_t=3) -- dim_K = 48."""
    return CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=3)


@pytest.fixture
def op_sys_1_1():
    """(n_max=1, N_t=1) -- minimum cell, dim_K = 4."""
    return CompactTemporalTruncatedOperatorSystem(n_max=1, N_t=1)


@pytest.fixture
def state_space_2_1(op_sys_2_1):
    return KreinPositiveStateSpace(op_sys_2_1)


@pytest.fixture
def state_space_2_3(op_sys_2_3):
    return KreinPositiveStateSpace(op_sys_2_3)


@pytest.fixture
def state_space_1_1(op_sys_1_1):
    return KreinPositiveStateSpace(op_sys_1_1)


# ---------------------------------------------------------------------------
# (a) Instantiation
# ---------------------------------------------------------------------------


def test_instantiate_small_panel(op_sys_1_1):
    ks = KreinPositiveStateSpace(op_sys_1_1)
    assert ks.dim_K == op_sys_1_1.dim_K
    assert ks.J.shape == (ks.dim_K, ks.dim_K)


def test_instantiate_2_1(state_space_2_1):
    ks = state_space_2_1
    assert ks.dim_K == 16  # full_dirac_dim(2) * 1


def test_instantiate_2_3(state_space_2_3):
    ks = state_space_2_3
    assert ks.dim_K == 48  # full_dirac_dim(2) * 3


def test_repr_format(state_space_2_1):
    s = repr(state_space_2_1)
    assert "KreinPositiveStateSpace" in s
    assert "dim_K=16" in s


# ---------------------------------------------------------------------------
# (b) J eigendecomposition: chirality doubling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (2, 3), (1, 3)])
def test_chirality_doubling_exact(n_max, N_t):
    """dim(+1) = dim(-1) = dim_K / 2 exactly."""
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    ks = KreinPositiveStateSpace(op_sys)
    assert ks.K_plus_dim + ks.K_minus_dim == ks.dim_K
    assert ks.K_plus_dim == ks.dim_K // 2
    assert ks.K_minus_dim == ks.dim_K // 2


def test_J_eigenvalues_pm1(state_space_2_1):
    """J eigenvalues are exactly +/- 1."""
    ks = state_space_2_1
    ok, det = ks.verify_J_eigendecomp()
    assert ok
    assert det["eigenvalues_pm1"]
    # Should be very close to +/- 1 numerically
    assert abs(det["min_eig"] - (-1.0)) < 1e-12
    assert abs(det["max_eig"] - (+1.0)) < 1e-12


def test_J_eigenvectors_shapes(state_space_2_1):
    """K^+ and K^- eigenvector blocks have the right shapes."""
    ks = state_space_2_1
    assert ks.K_plus_eigvecs.shape == (ks.dim_K, ks.K_plus_dim)
    assert ks.K_minus_eigvecs.shape == (ks.dim_K, ks.K_minus_dim)


def test_J_eigvecs_orthonormal(state_space_2_1):
    """eigh returns orthonormal eigenvectors."""
    ks = state_space_2_1
    # K+ block orthonormal
    U_plus = ks.K_plus_eigvecs
    I_plus = U_plus.conj().T @ U_plus
    assert np.allclose(I_plus, np.eye(ks.K_plus_dim), atol=1e-12)
    # K- block orthonormal
    U_minus = ks.K_minus_eigvecs
    I_minus = U_minus.conj().T @ U_minus
    assert np.allclose(I_minus, np.eye(ks.K_minus_dim), atol=1e-12)


def test_J_eigvecs_are_eigenvectors(state_space_2_1):
    """J U_plus = +U_plus, J U_minus = -U_minus."""
    ks = state_space_2_1
    JU_plus = ks.J @ ks.K_plus_eigvecs
    JU_minus = ks.J @ ks.K_minus_eigvecs
    assert np.allclose(JU_plus, +ks.K_plus_eigvecs, atol=1e-12)
    assert np.allclose(JU_minus, -ks.K_minus_eigvecs, atol=1e-12)


# ---------------------------------------------------------------------------
# (c) Projectors P_+, P_-: orthogonal idempotents summing to I
# ---------------------------------------------------------------------------


def test_projectors_idempotent(state_space_2_1):
    ks = state_space_2_1
    P_plus = ks.P_plus
    P_minus = ks.P_minus
    # P_+^2 = P_+
    assert np.allclose(P_plus @ P_plus, P_plus, atol=1e-12)
    assert np.allclose(P_minus @ P_minus, P_minus, atol=1e-12)


def test_projectors_orthogonal(state_space_2_1):
    """P_+ P_- = 0 (orthogonal eigenspaces)."""
    ks = state_space_2_1
    P_plus = ks.P_plus
    P_minus = ks.P_minus
    assert np.allclose(P_plus @ P_minus, 0, atol=1e-12)


def test_projectors_sum_to_identity(state_space_2_1):
    """P_+ + P_- = I_dim_K."""
    ks = state_space_2_1
    P_plus = ks.P_plus
    P_minus = ks.P_minus
    I_K = np.eye(ks.dim_K, dtype=np.complex128)
    assert np.allclose(P_plus + P_minus, I_K, atol=1e-12)


def test_projector_traces_equal_dims(state_space_2_3):
    """Tr(P_+) = dim_K_plus, Tr(P_-) = dim_K_minus."""
    ks = state_space_2_3
    tr_plus = np.trace(ks.P_plus).real
    tr_minus = np.trace(ks.P_minus).real
    assert abs(tr_plus - ks.K_plus_dim) < 1e-10
    assert abs(tr_minus - ks.K_minus_dim) < 1e-10


# ---------------------------------------------------------------------------
# (d) Pure-state densities from K^+ eigenvectors
# ---------------------------------------------------------------------------


def test_pure_state_density_hermitian(state_space_2_1):
    ks = state_space_2_1
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    assert np.allclose(rho, rho.conj().T, atol=1e-12)


def test_pure_state_density_trace_one(state_space_2_1):
    ks = state_space_2_1
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    assert abs(np.trace(rho).real - 1.0) < 1e-12


def test_pure_state_density_rank_one(state_space_2_1):
    """A rank-1 projector has spectrum {1, 0, ..., 0}."""
    ks = state_space_2_1
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    eigs = np.linalg.eigvalsh(rho)
    # One eigenvalue = 1, rest = 0
    assert abs(eigs.max() - 1.0) < 1e-10
    # Sum of all except the top = 0
    assert abs(eigs[:-1].sum()) < 1e-10


def test_pure_state_K_plus_J_trace(state_space_2_1):
    """K^+ pure-state structural marker: Tr(rho J) = +1.

    For v in K^+ (J v = +v): rho = |v><v|, so
    Tr(rho J) = <v| J |v> = <v| v> = 1.
    """
    ks = state_space_2_1
    for idx in range(min(4, ks.K_plus_dim)):
        v = ks.K_plus_eigvecs[:, idx]
        rho = ks.pure_state_density(v)
        tr_J = np.trace(rho @ ks.J).real
        assert abs(tr_J - 1.0) < 1e-10, (
            f"K+ idx={idx}: Tr(rho J) = {tr_J}, expected +1"
        )


def test_pure_state_K_minus_J_trace(state_space_2_1):
    """K^- pure-state structural marker: Tr(rho J) = -1.

    For v in K^- (J v = -v): Tr(rho J) = <v| J |v> = -<v|v> = -1.
    This is the distinguishing trace observable between K^+ and K^-
    pure states.
    """
    ks = state_space_2_1
    for idx in range(min(4, ks.K_minus_dim)):
        v = ks.K_minus_eigvecs[:, idx]
        rho = ks.pure_state_density(v)
        tr_J = np.trace(rho @ ks.J).real
        assert abs(tr_J - (-1.0)) < 1e-10, (
            f"K- idx={idx}: Tr(rho J) = {tr_J}, expected -1"
        )


# ---------------------------------------------------------------------------
# (f) Krein-positivity check: passes on BOTH K^+ and K^- pure states
# ---------------------------------------------------------------------------


def test_K_plus_pure_state_krein_positive(state_space_2_1):
    """A K^+ pure state passes the generator-level Krein-positivity check.

    Per the L3a-1 'trivial Krein-positivity at operator-multiplier level'
    finding (every multiplier commutes with J because J^2 = I and the
    multipliers are block-diagonal in chirality), the check passes for
    any density matrix.  We verify the explicit pass for a K^+ pure state.
    """
    ks = state_space_2_1
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    ok, det = ks.is_krein_positive_state(rho)
    assert ok
    assert det["n_violated"] == 0
    # All values are real and non-negative
    assert det["min_real_value"] >= -1e-10
    assert det["max_imag_value"] < 1e-10


def test_K_minus_pure_state_krein_positive_too(state_space_2_1):
    """K^- pure state also passes the generator-level inequality.

    This is the EXPECTED behavior matching the module docstring's
    'trivial Krein-positivity at operator-multiplier level' note.  The
    state-space-level distinguishing observable is Tr(rho J) (see
    test_pure_state_K_minus_J_trace), NOT the operator-multiplier
    inequality.
    """
    ks = state_space_2_1
    v = ks.K_minus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    ok, det = ks.is_krein_positive_state(rho)
    # Both K+ and K- pass at the multiplier level (L3a-1 finding)
    assert ok
    assert det["n_violated"] == 0


def test_krein_positive_check_handles_sample_size(state_space_2_3):
    """sample_size limits the number of generators checked."""
    ks = state_space_2_3
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    ok_5, det_5 = ks.is_krein_positive_state(rho, sample_size=5)
    assert det_5["n_total"] == 5
    ok_all, det_all = ks.is_krein_positive_state(rho)
    assert det_all["n_total"] == len(ks.op_sys.multiplier_matrices)


# ---------------------------------------------------------------------------
# (g) Wasserstein-Kantorovich SDP distance
# ---------------------------------------------------------------------------


def test_wasserstein_distance_diagonal(state_space_2_1):
    """d(idx, idx) = 0 exactly (no SDP solve needed)."""
    ks = state_space_2_1
    D = np.eye(ks.dim_K, dtype=np.complex128)
    d = ks.wasserstein_distance_pure(0, 0, D)
    assert d == 0.0


def test_wasserstein_distance_symmetric(state_space_2_1):
    """d(v, w) = d(w, v) by SDP feasible-set symmetry.

    Use offdiag CH spatial Dirac to get a finite, non-trivial distance.
    At N_t = 1 the Krein Dirac reduces to the spatial Dirac.
    """
    pytest.importorskip("cvxpy")
    ks = state_space_2_1
    # Build offdiag CH spatial Dirac on the FullDirac basis
    D_spat = camporesi_higuchi_offdiag_dirac_matrix(ks.op_sys.basis_spatial)
    # At N_t = 1, D = D_spat (kron with I_1)
    D = D_spat
    # Pick two indices with a chance of finite distance
    d_01 = ks.wasserstein_distance_pure(0, 5, D)
    d_10 = ks.wasserstein_distance_pure(5, 0, D)
    # Both should be non-negative; equal up to SDP tolerance
    assert d_01 >= 0.0
    assert d_10 >= 0.0
    if np.isfinite(d_01) and np.isfinite(d_10):
        assert abs(d_01 - d_10) < 1e-4, (
            f"d(0,5) = {d_01}, d(5,0) = {d_10}, not symmetric"
        )


def test_wasserstein_returns_finite_offdiag_CH(state_space_2_1):
    """At (n_max=2, N_t=1) with offdiag CH spatial Dirac, the Wasserstein
    distance produces some finite value (not all +inf).

    The offdiag CH Dirac is the SDP-bounding device from R3.5 / L3a-1;
    it has native cross-chirality off-diagonal couplings that break the
    n-degeneracy obstruction of the truthful CH.
    """
    pytest.importorskip("cvxpy")
    ks = state_space_2_1
    D_spat = camporesi_higuchi_offdiag_dirac_matrix(ks.op_sys.basis_spatial)
    D = D_spat  # at N_t=1
    # Sample a few pairs; expect at least one to be finite (and >= 0)
    finite_count = 0
    for idx_v, idx_w in [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2)]:
        d = ks.wasserstein_distance_pure(idx_v, idx_w, D)
        assert d >= 0.0
        if np.isfinite(d):
            finite_count += 1
    assert finite_count > 0, (
        "No finite Wasserstein distances under offdiag CH at (n_max=2, N_t=1)"
    )


def test_wasserstein_invalid_index_raises(state_space_2_1):
    """Out-of-range index raises ValueError."""
    ks = state_space_2_1
    D = np.eye(ks.dim_K, dtype=np.complex128)
    with pytest.raises(ValueError):
        ks.wasserstein_distance_pure(ks.dim_K, 0, D)
    with pytest.raises(ValueError):
        ks.wasserstein_distance_pure(0, -1, D)


def test_wasserstein_invalid_D_shape_raises(state_space_2_1):
    """D with wrong shape raises ValueError."""
    ks = state_space_2_1
    D_wrong = np.eye(5, dtype=np.complex128)
    with pytest.raises(ValueError):
        ks.wasserstein_distance_pure(0, 1, D_wrong)


# ---------------------------------------------------------------------------
# (h) Riemannian limit at N_t = 1 (LOAD-BEARING)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_N_t_1(n_max):
    """LOAD-BEARING: at N_t = 1, dim_K = full_dirac_dim(n_max) and
    dim_K_plus = full_dirac_dim(n_max) / 2.

    This is the structural analog of L2-B's Riemannian-limit check at
    the state-space level: the K^+ block reduces to the chirality-+
    half of the spatial FullDirac sector.
    """
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=1)
    ks = KreinPositiveStateSpace(op_sys)
    expected_dim_spatial = full_dirac_dim(n_max)
    assert ks.dim_K == expected_dim_spatial
    assert ks.K_plus_dim == expected_dim_spatial // 2
    assert ks.K_minus_dim == expected_dim_spatial // 2


def test_riemannian_limit_chirality_doubling_strict(state_space_2_1):
    """At N_t = 1, dim_K_plus = dim_K_minus is forced by chirality
    doubling (FullDirac sector is built as Weyl + anti-Weyl with the
    same multiplicity per level)."""
    ks = state_space_2_1
    assert ks.K_plus_dim == ks.K_minus_dim
    # Each is exactly half
    assert 2 * ks.K_plus_dim == ks.dim_K


# ---------------------------------------------------------------------------
# (e -- bonus) restrict_to_K_plus produces the K^+ block in K^+ basis
# ---------------------------------------------------------------------------


def test_restrict_to_K_plus_shape(state_space_2_1):
    """restrict_to_K_plus returns a (K_plus_dim, K_plus_dim) matrix."""
    ks = state_space_2_1
    # Use an arbitrary density (mix of K^+ and K^- to make the projection
    # non-trivial)
    v_plus = ks.K_plus_eigvecs[:, 0]
    v_minus = ks.K_minus_eigvecs[:, 0]
    mix = (v_plus + v_minus) / np.sqrt(2)
    rho_mix = ks.pure_state_density(mix)
    rho_plus_block = ks.restrict_to_K_plus(rho_mix)
    assert rho_plus_block.shape == (ks.K_plus_dim, ks.K_plus_dim)


def test_restrict_to_K_plus_pure_K_plus_idempotent(state_space_2_1):
    """For a pure K^+ state in K^+, the restriction has trace 1 in the
    K^+ block (full mass preserved)."""
    ks = state_space_2_1
    v = ks.K_plus_eigvecs[:, 0]
    rho = ks.pure_state_density(v)
    rho_plus_block = ks.restrict_to_K_plus(rho)
    # In the K^+ basis, the projection picks up only the (idx, idx) coordinate
    assert abs(np.trace(rho_plus_block).real - 1.0) < 1e-10
