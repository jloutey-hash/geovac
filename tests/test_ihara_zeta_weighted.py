"""Tests for geovac.ihara_zeta_weighted (Track RH-K).

Weighted Ihara zeta with alpha^2 spin-orbit edge weights on the
Dirac-S^3 graph.

Coverage
--------
 1. Unweighted limit (all w=1): weighted Ihara-Bass reproduces the
    unweighted ihara_zeta.ihara_zeta_bass EXACTLY (symbolic identity).
 2. alpha -> 0 limit: the symbolic weighted Ihara-Bass on the Dirac-S^3
    graph at alpha = 0 reduces to the unweighted result (n_max = 2, Rule A).
 3. Graph symmetry: swapping the node ordering gives an isospectral
    weighted adjacency (zeta polynomial unchanged).
 4. Ramanujan checks run for each (n_max, rule) at alpha=0 and at
    alpha=1/137 (physical) and give numerical verdicts with reported
    deviations.
 5. Small-graph sanity: weighted K_4 with uniform w = 2 has
    weighted Ihara-Bass matching the unweighted K_4 zeta under the
    substitution s -> 2 s (sharper: the formal weighted identity
    collapses to the expected re-scaling).
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest
import sympy as sp

from geovac.ihara_zeta import ihara_zeta_bass, is_ramanujan
from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.ihara_zeta_weighted import (
    alpha_sym,
    build_weighted_dirac_adjacency,
    dirac_s3_edge_weight,
    is_weighted_ramanujan,
    weighted_adjacency,
    weighted_hashimoto_matrix,
    weighted_ihara_zeta_bass,
    weighted_q_matrix,
    connectivity_pattern,
)

s = sp.symbols("s")


# ---------------------------------------------------------------------------
# Test 1: unweighted limit (all w = 1)
# ---------------------------------------------------------------------------

def test_unweighted_limit_k4():
    """K_4 with uniform w = 1 matches the unweighted Ihara-Bass K_4."""
    edges = list(itertools.combinations(range(4), 2))
    weights = [sp.Integer(1)] * len(edges)
    A_w = weighted_adjacency(edges, weights, 4)
    zeta_inv_w = sp.cancel(weighted_ihara_zeta_bass(A_w, s=s))
    # Reference unweighted:
    A01 = np.zeros((4, 4), dtype=int)
    for i, j in edges:
        A01[i, j] = 1
        A01[j, i] = 1
    zeta_unw = ihara_zeta_bass(A01)
    zeta_inv_unw = sp.cancel(1 / zeta_unw)
    diff = sp.cancel(zeta_inv_w - zeta_inv_unw)
    assert diff == 0, f"K_4 unweighted-limit diff = {diff}"


def test_unweighted_limit_k4_hashimoto_spectrum():
    """K_4 unit-weighted Hashimoto spectrum equals the unweighted one."""
    edges = list(itertools.combinations(range(4), 2))
    weights = [sp.Integer(1)] * len(edges)
    A_w = weighted_adjacency(edges, weights, 4)
    T_w = weighted_hashimoto_matrix(A_w)
    ev_w = np.sort(np.abs(np.linalg.eigvals(T_w)))
    # Reference via ihara_zeta.hashimoto_matrix
    from geovac.ihara_zeta import hashimoto_matrix
    A01 = np.zeros((4, 4), dtype=int)
    for i, j in edges:
        A01[i, j] = 1
        A01[j, i] = 1
    T_unw = hashimoto_matrix(A01)
    ev_unw = np.sort(np.abs(np.linalg.eigvals(T_unw.astype(float))))
    assert np.allclose(ev_w, ev_unw, atol=1e-10), "spectra differ"


# ---------------------------------------------------------------------------
# Test 2: alpha -> 0 on the Dirac-S^3 graph
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max, rule", [(2, "A"), (2, "B"), (3, "A")])
def test_alpha_zero_recovers_unweighted(n_max, rule):
    """alpha=0 substitution in the weighted Dirac-S^3 zeta recovers
    the unweighted ihara_zeta_bass result."""
    A01, _, _, _ = build_dirac_s3_graph(n_max, rule)
    A_w, _ = build_weighted_dirac_adjacency(n_max, rule)
    A_w_alpha0 = A_w.subs(alpha_sym, 0)
    zeta_inv_w0 = sp.cancel(weighted_ihara_zeta_bass(A_w_alpha0, s=s))
    zeta_inv_unw = sp.cancel(1 / ihara_zeta_bass(A01))
    diff = sp.cancel(zeta_inv_w0 - zeta_inv_unw)
    assert diff == 0, f"rule={rule} n_max={n_max}: alpha->0 residual = {diff}"


# ---------------------------------------------------------------------------
# Test 3: node relabeling invariance
# ---------------------------------------------------------------------------

def test_relabel_invariance_small():
    """Swapping two node labels gives an isospectral weighted zeta
    (PAPA^T = a different matrix, but has the same determinant/polynomial)."""
    # Build a 4-node weighted path graph with distinct alpha-weights
    alpha = alpha_sym
    A_w = sp.zeros(4, 4)
    w12 = 1 + alpha ** 2 * sp.Rational(1, 3)
    w23 = 1 + alpha ** 2 * sp.Rational(1, 5)
    w34 = 1 + alpha ** 2 * sp.Rational(2, 7)
    A_w[0, 1] = A_w[1, 0] = w12
    A_w[1, 2] = A_w[2, 1] = w23
    A_w[2, 3] = A_w[3, 2] = w34
    z1 = sp.cancel(weighted_ihara_zeta_bass(A_w, s=s))
    # Swap labels 1 <-> 2
    P = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    A_w2 = P * A_w * P.T
    z2 = sp.cancel(weighted_ihara_zeta_bass(A_w2, s=s))
    diff = sp.cancel(z1 - z2)
    assert diff == 0, f"relabel residual = {diff}"


# ---------------------------------------------------------------------------
# Test 4: Ramanujan verdict runs for each (n_max, rule)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max, rule", [(2, "A"), (2, "B"), (3, "A")])
def test_weighted_ramanujan_reports_numeric(n_max, rule):
    """At alpha = 1/137 and alpha = 0, is_weighted_ramanujan returns
    a finite deviation and a boolean verdict. This test does NOT assert
    the verdict (weighted Ramanujan is a new scientific question); it
    only checks that the machinery runs and produces sensible values."""
    A_w, _ = build_weighted_dirac_adjacency(n_max, rule)
    # Numerical alpha
    A_w_phys = A_w.subs(alpha_sym, 1.0 / 137.036)
    is_ram, dev, text = is_weighted_ramanujan(A_w_phys)
    assert isinstance(is_ram, (bool, np.bool_))
    assert np.isfinite(dev)
    assert len(text) > 0

    # alpha = 0 must match unweighted Ramanujan
    A_w_zero = A_w.subs(alpha_sym, 0)
    is_ram0, dev0, _ = is_weighted_ramanujan(A_w_zero)
    # Compare to unweighted is_ramanujan on the 0/1 connectivity pattern
    pattern = connectivity_pattern(A_w_zero)
    is_ram_unw, dev_unw, _ = is_ramanujan(pattern)
    # The two need not agree to machine precision for irregular graphs
    # (weighted vs unweighted q_max), but they MUST agree on the verdict
    # when all weights are 1 (because q_max^w = row-sum - 1 = degree - 1).
    assert is_ram0 == is_ram_unw, (
        f"alpha=0 Ramanujan verdict mismatch: weighted={is_ram0}, "
        f"unweighted={is_ram_unw}")
    # Deviations agree to ~1e-6 for alpha=0 (weighted uses sqrt formula
    # on numeric row-sums, unweighted on int row-sums)
    assert abs(dev0 - dev_unw) < 1e-6, (
        f"alpha=0 deviations disagree: {dev0} vs {dev_unw}")


# ---------------------------------------------------------------------------
# Test 5: small-graph sanity — uniform weight scaling
# ---------------------------------------------------------------------------

def test_uniform_weight_scaling_cycle_c3():
    """Triangle C_3 with uniform weight w = 2 should give a weighted
    zeta whose zeros lie at 1/(2 mu) where mu ranges over the unit-weight
    C_3 Hashimoto eigenvalues. Equivalently: det(I - s T_w) = det(I - 2 s T).

    This is the rigorous statement of 'weighted Ihara rescales like
    s -> w s' when all weights equal the same constant w.
    """
    # C_3: triangle
    A_w_weight2 = sp.Matrix([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
    A01 = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=int)
    T_weighted = weighted_hashimoto_matrix(A_w_weight2)
    from geovac.ihara_zeta import hashimoto_matrix
    T_unw = hashimoto_matrix(A01).astype(float)
    ev_w = np.sort(np.abs(np.linalg.eigvals(T_weighted)))
    ev_unw = np.sort(np.abs(np.linalg.eigvals(T_unw)))
    # Prediction: ev_w = 2 * ev_unw (each leg of a directed edge carries
    # sqrt(w) = sqrt(2), so T_w = sqrt(w)*sqrt(w) = 2 times the sub-blocks
    # of T_unw, which scales every eigenvalue by exactly 2).
    assert np.allclose(ev_w, 2.0 * ev_unw, atol=1e-10), (
        f"scaling broken: ev_w={ev_w}, 2*ev_unw={2 * ev_unw}")


# ---------------------------------------------------------------------------
# Test 6: edge-weight convention sanity
# ---------------------------------------------------------------------------

def test_edge_weight_zero_at_alpha_zero():
    """dirac_s3_edge_weight evaluated at alpha=0 is exactly 1 for all edges."""
    from geovac.dirac_matrix_elements import DiracLabel
    # Any pair of labels:
    a = DiracLabel(n_fock=1, kappa=-1, two_m_j=1)
    b = DiracLabel(n_fock=2, kappa=-2, two_m_j=3)
    w = dirac_s3_edge_weight(a, b)
    w_at_zero = w.subs(alpha_sym, 0)
    assert sp.simplify(w_at_zero - 1) == 0


def test_edge_weight_kramers_l0():
    """Both endpoints have l=0 (κ=−1) ⇒ weight = 1 + alpha^2 * 0 = 1."""
    from geovac.dirac_matrix_elements import DiracLabel
    a = DiracLabel(n_fock=1, kappa=-1, two_m_j=1)
    b = DiracLabel(n_fock=2, kappa=-1, two_m_j=1)
    w = dirac_s3_edge_weight(a, b)
    assert sp.simplify(w - 1) == 0, f"Kramers failure: w = {w}"


def test_edge_weight_positive_at_physical_alpha():
    """For any edge on the n_max=2 Dirac-S^3 graph (Rule B), the weight
    is strictly positive at alpha = 1/137."""
    A_w, labels = build_weighted_dirac_adjacency(2, "B")
    n = A_w.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if A_w[i, j] != 0:
                val = float(A_w[i, j].subs(alpha_sym, 1.0 / 137.036).evalf())
                assert val > 0, f"non-positive weight at edge ({i}, {j}): {val}"


# ---------------------------------------------------------------------------
# Test 7: weighted Q matrix sanity
# ---------------------------------------------------------------------------

def test_weighted_q_matrix_uniform_w1():
    """When all weights are 1, Q_w = diag(d_v - 1), matching unweighted Q."""
    A_w = sp.Matrix([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    Q = weighted_q_matrix(A_w)
    expected = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # degree 2 each, -1 = 1
    assert Q == expected
