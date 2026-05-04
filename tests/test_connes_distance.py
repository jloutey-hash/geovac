"""Tests for geovac.connes_distance: Connes distance d_D(phi_v, phi_w) on
the GeoVac truncated operator system O_{n_max}.

Verifies metric structure (positivity, symmetry, triangle inequality) and
sanity properties of the SDP-based computation. The absolute numerical
values are placeholder-convention dependent; the tested invariants are
not.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.operator_system import HyperLabel, TruncatedOperatorSystem

# Skip the entire module if cvxpy is not installed (it is an optional
# dependency for the round-3 metric layer).
cvxpy = pytest.importorskip("cvxpy")

from geovac.connes_distance import (  # noqa: E402  (import after skip)
    DiracProxy,
    _kernel_of_commutator_in_O,
    compute_connes_distance,
    compute_distance_matrix,
    default_dirac,
    fock_graph_distance,
    graph_distance_matrix,
)


# ---------------------------------------------------------------------------
# Trivial cases
# ---------------------------------------------------------------------------


def test_distance_at_n_max_1_zero():
    """At n_max = 1 there is one basis vector; d(v, v) = 0 trivially."""
    op = TruncatedOperatorSystem(1)
    d = compute_connes_distance(op, 0, 0)
    assert d == 0.0


def test_distance_diagonal_zero_n_max_2():
    """d(v, v) = 0 for every basis vector at n_max = 2."""
    op = TruncatedOperatorSystem(2)
    for v_idx in range(op.dim_H):
        d = compute_connes_distance(op, v_idx, v_idx)
        assert d == 0.0


# ---------------------------------------------------------------------------
# Symmetry: d(v, w) = d(w, v)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "v_idx,w_idx",
    [(0, 1), (0, 3), (1, 3), (0, 2)],
)
def test_distance_symmetry_n_max_2(v_idx, w_idx):
    """d(v, w) = d(w, v) at n_max = 2 to SDP precision."""
    op = TruncatedOperatorSystem(2)
    d_vw = compute_connes_distance(op, v_idx, w_idx)
    d_wv = compute_connes_distance(op, w_idx, v_idx)
    if np.isinf(d_vw) and np.isinf(d_wv):
        return
    # Both finite: agreement to relative SDP precision.
    rel = abs(d_vw - d_wv) / max(abs(d_vw), abs(d_wv), 1e-10)
    assert rel < 1e-3, f"d(v,w)={d_vw}, d(w,v)={d_wv}, rel={rel}"


# ---------------------------------------------------------------------------
# Positivity: d(v, w) >= 0
# ---------------------------------------------------------------------------


def test_distances_nonnegative_n_max_2():
    """Every Connes distance is >= 0 (positive-semidefinite metric)."""
    op = TruncatedOperatorSystem(2)
    for i in range(op.dim_H):
        for j in range(i + 1, op.dim_H):
            d = compute_connes_distance(op, i, j)
            assert d >= 0.0, f"negative distance d({i},{j})={d}"


# ---------------------------------------------------------------------------
# Triangle inequality
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "u_idx,v_idx,w_idx",
    [(0, 1, 3), (0, 3, 4), (1, 3, 0)],
)
def test_triangle_inequality_n_max_2(u_idx, v_idx, w_idx):
    """d(u, w) <= d(u, v) + d(v, w) for finite distances."""
    op = TruncatedOperatorSystem(2)
    d_uw = compute_connes_distance(op, u_idx, w_idx)
    d_uv = compute_connes_distance(op, u_idx, v_idx)
    d_vw = compute_connes_distance(op, v_idx, w_idx)
    if not (np.isfinite(d_uw) and np.isfinite(d_uv) and np.isfinite(d_vw)):
        # Triangle inequality is vacuous if any side is +infinity.
        return
    # Allow 1e-3 SDP slack
    assert d_uw <= d_uv + d_vw + 1e-3 * max(d_uv + d_vw, 1.0), (
        f"triangle violation: d(u,w)={d_uw}, d(u,v)+d(v,w)={d_uv+d_vw}"
    )


# ---------------------------------------------------------------------------
# Distance matrix at n_max = 2
# ---------------------------------------------------------------------------


def test_distance_matrix_shape_and_zero_diagonal():
    """compute_distance_matrix returns (N, N) symmetric with zero diagonal."""
    op = TruncatedOperatorSystem(2)
    M = compute_distance_matrix(op)
    assert M.shape == (op.dim_H, op.dim_H)
    # Symmetry to SDP precision.
    diff = M - M.T
    finite = np.isfinite(M) & np.isfinite(M.T)
    if finite.any():
        assert (np.abs(diff)[finite] < 1e-3 * (np.abs(M)[finite] + 1.0)).all()
    # Diagonal exactly zero.
    np.testing.assert_array_equal(np.diag(M), np.zeros(op.dim_H))


# ---------------------------------------------------------------------------
# Dirac proxies and the role of degeneracy
# ---------------------------------------------------------------------------


def test_shell_scalar_dirac_collapses_metric():
    """With the diagonal Camporesi-Higuchi-style Dirac the metric is
    +infinity on (almost) every pair, exhibiting the eigenvalue-degeneracy
    obstruction. This is the truthful diagnostic."""
    op = TruncatedOperatorSystem(2)
    D = DiracProxy(n_max=2, mode="shell_scalar").matrix(op.basis)
    # Pair (0, 1) is in different shells but the kernel still contains
    # diagonal-by-shell multipliers, so the SDP is unbounded -> +inf.
    d = compute_connes_distance(op, 0, 1, D=D)
    assert np.isinf(d) or d > 1e6


def test_offdiag_dirac_kernel_is_identity_only():
    """With the off-diagonal Dirac proxy ker([D, .]) cap O = C * I."""
    op = TruncatedOperatorSystem(2)
    D = default_dirac(op)
    ker = _kernel_of_commutator_in_O(op, D)
    # Single basis element, proportional to identity.
    assert len(ker) == 1
    K = ker[0]
    # K is proportional to I:
    K_normalized = K / np.linalg.norm(K)
    I_normalized = np.eye(op.dim_H) / np.linalg.norm(np.eye(op.dim_H))
    overlap = abs(np.trace(K_normalized.conj().T @ I_normalized))
    assert overlap > 0.999, f"kernel basis overlap with I = {overlap}"


# ---------------------------------------------------------------------------
# Fock-graph distance baseline
# ---------------------------------------------------------------------------


def test_graph_distance_self_zero():
    """d_graph(v, v) = 0 for every node."""
    for n in range(1, 4):
        op = TruncatedOperatorSystem(n)
        for b in op.basis:
            assert fock_graph_distance(b, b) == 0


def test_graph_distance_adjacent_shells():
    """d_graph(|1,0,0>, |2,0,0>) = 1 (Delta n = 1, Delta l = 0, Delta m = 0)."""
    v = HyperLabel(n=1, l=0, m=0)
    w = HyperLabel(n=2, l=0, m=0)
    assert fock_graph_distance(v, w) == 1


def test_graph_distance_matrix_n_max_2():
    """d_graph matrix has zero diagonal and is symmetric."""
    op = TruncatedOperatorSystem(2)
    G = graph_distance_matrix(op)
    assert G.shape == (op.dim_H, op.dim_H)
    np.testing.assert_array_equal(G, G.T)
    np.testing.assert_array_equal(np.diag(G), np.zeros(op.dim_H, dtype=int))


# ---------------------------------------------------------------------------
# SO(3) m-symmetry: d(|n,l,-m>, |n,l,+m>) = 0 by operator-system content
# ---------------------------------------------------------------------------


def test_m_reflection_symmetry_zero_distance_n_max_2():
    """Pure-state distance d(|2,1,-1>, |2,1,+1>) = 0 because every multiplier
    in O has identical diagonal entries at the m=-1 and m=+1 indices (SO(3)
    symmetry of the operator system construction). This is structural,
    independent of the Dirac choice."""
    op = TruncatedOperatorSystem(2)
    # Find indices of |2,1,-1> and |2,1,+1>.
    idx_m1 = op.basis.index(HyperLabel(n=2, l=1, m=-1))
    idx_p1 = op.basis.index(HyperLabel(n=2, l=1, m=+1))
    d = compute_connes_distance(op, idx_m1, idx_p1)
    # SDP precision: |d| < a few mHz of the structural zero.
    assert d < 1e-2, f"d(m=-1, m=+1) = {d}, expected ~0 by SO(3) symmetry"
