"""Entropy-per-pair-state decomposition for Paper 27 sec:cusp (result 3).

Backs the co-location claim that the 2026-08-28 cert demoted to an inference:
"the pair-state with maximal V_ee contribution is also the pair-state on which
the two-body-generated entropy concentrates."  A decomposition now exists.

Method.  In the singlet pair-state basis the V_ee graph's NODES and the CI
configurations are the same objects, so entropy can be attributed to graph
EDGES directly: zero one coupling V[ref, J] out of the dominant reference
node, re-solve, and measure the entropy lost.  The attribution is exact (no
perturbative step) and its control is sharp -- an edge not incident on the
reference must cost no entropy at all.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from debug.energy_graph_exploration import (  # noqa: E402
    build_vee_matrix, build_h1_pair_matrix,
)
from debug.archive.chemistry_qc_arc.entanglement_molecular import (  # noqa: E402
    build_1rdm_from_singlet_ci, compute_entanglement_entropy,
)

Z_HE = 2


def _system(n_max):
    V, pairs, orbitals = build_vee_matrix(n_max, m_total=0)
    h1 = build_h1_pair_matrix(n_max, Z_HE, m_total=0)
    if isinstance(h1, tuple):
        h1 = h1[0]
    return V, h1, pairs, orbitals


def _entropy(V, h1, pairs, n_orb):
    w, U = np.linalg.eigh(h1 + V)
    rho = build_1rdm_from_singlet_ci(U[:, 0], pairs, n_orb)
    S, _ = compute_entanglement_entropy(rho)
    return w[0], U[:, 0], S


def _decompose(n_max):
    V, h1, pairs, orbitals = _system(n_max)
    n_orb = len(orbitals)
    _, c0, S_full = _entropy(V, h1, pairs, n_orb)
    ref = int(np.argmax(np.abs(c0)))

    rows = []
    for J in range(len(pairs)):
        if J == ref or abs(V[ref, J]) < 1e-14:
            continue
        Vc = V.copy()
        Vc[ref, J] = Vc[J, ref] = 0.0
        dS = S_full - _entropy(Vc, h1, pairs, n_orb)[2]
        rows.append((dS, abs(V[ref, J]), J))
    rows.sort(reverse=True)
    return V, h1, pairs, orbitals, ref, S_full, rows


def _hottest_edge(V):
    off = V - np.diag(np.diag(V))
    return set(np.unravel_index(np.argmax(np.abs(off)), off.shape))


def test_paper27_entropy_concentrates_on_the_cusp_edge():
    """The V_ee-hottest edge IS the entropy-dominant edge, and it carries a
    large share of the total -- the co-location, measured."""
    V, h1, pairs, orbitals, ref, S_full, rows = _decompose(3)

    assert S_full > 1e-4, f"no entropy to attribute: {S_full}"
    # the reference is the (1s,1s) pair-state
    assert orbitals[pairs[ref][0]][:2] == (1, 0) and orbitals[pairs[ref][1]][:2] == (1, 0)

    top_dS, _, top_J = rows[0]
    assert {ref, top_J} == _hottest_edge(V), (
        "entropy-dominant edge is not the V_ee-hottest edge")
    # and it is dominant by a clear margin, not a photo finish
    assert top_dS / S_full > 0.35, f"top edge carries only {top_dS/S_full:.1%}"
    assert top_dS > 1.3 * rows[1][0], "top edge not clearly ahead of the runner-up"


def test_paper27_entropy_locus_control_and_grading():
    """Two things that make the co-location non-trivial: (i) a LARGE V_ee edge
    that does not touch the reference costs exactly zero entropy -- so it is
    mass at the cusp vertex, not V_ee mass anywhere; (ii) across all edges out
    of the reference the entropy ordering tracks the V_ee ordering."""
    V, h1, pairs, orbitals, ref, S_full, rows = _decompose(3)
    n_orb = len(orbitals)

    # (i) control
    cands = [(a, b) for a in range(len(pairs)) for b in range(a + 1, len(pairs))
             if ref not in (a, b) and abs(V[a, b]) > 1e-6]
    a, b = max(cands, key=lambda ab: abs(V[ab[0], ab[1]]))
    assert abs(V[a, b]) > 0.3 * rows[0][1], "control edge should be comparably large"
    Vc = V.copy()
    Vc[a, b] = Vc[b, a] = 0.0
    dS_ctrl = S_full - _entropy(Vc, h1, pairs, n_orb)[2]
    assert abs(dS_ctrl) < 1e-12 * max(1.0, S_full), (
        f"non-reference edge moved the entropy by {dS_ctrl:.2e}")

    # (ii) graded correspondence (Spearman on the edges out of the reference)
    dS = np.array([r[0] for r in rows])
    vv = np.array([r[1] for r in rows])
    rx = np.argsort(np.argsort(dS))
    ry = np.argsort(np.argsort(vv))
    rho = float(np.corrcoef(rx, ry)[0, 1])
    assert rho > 0.8, f"entropy/V_ee edge ordering only weakly related: rho={rho:.2f}"


@pytest.mark.slow
def test_paper27_entropy_locus_basis_robust():
    """The co-location is not an n_max=3 artifact: same edge at n_max=4, with
    a comparable share of the total entropy."""
    V, h1, pairs, orbitals, ref, S_full, rows = _decompose(4)
    top_dS, _, top_J = rows[0]
    assert {ref, top_J} == _hottest_edge(V)
    assert top_dS / S_full > 0.35, f"n_max=4 top edge carries {top_dS/S_full:.1%}"
