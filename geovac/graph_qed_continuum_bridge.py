"""
Graph-native to continuum QED bridge: identifying the projection exchange constant.
====================================================================================

This module systematically connects the graph-native QED (finite rational
matrices, pi-free) to the continuum spectral-sum QED (infinite mode sums,
transcendental results) and identifies WHERE transcendentals enter.

**HEADLINE RESULT** (found during implementation):

At ANY finite truncation, BOTH the graph VP and the continuum VP are
rational numbers (no pi).  The graph VP Tr(Pi_graph) is strictly richer
than the continuum VP Pi_continuum_truncated because the graph CG
projection from Dirac (n,kappa,m_j) states to Fock (n,l,m) states opens
coupling channels that the SO(4) channel count W classifies as zero.
At n_max_fock=2 (n_max_ch=1), the graph has Tr(Pi) = 224/75 while
the continuum VP is exactly ZERO (the only allowed triple (1,1,1)
has W=0).

The projection exchange constant at finite truncation (when both are
nonzero, i.e. n_max_fock >= 3) IS RATIONAL.  Transcendental content
(pi) enters ONLY via:
  (1) The infinite-sum regularization (Hurwitz zeta -> pi^{even})
  (2) Vol(S^3) = 2*pi^2 normalization for dimensionful quantities

This is the CALIBRATION tier of Paper 18's exchange constant taxonomy.

Architecture
------------
(A) Truncated continuum VP at matching n_max.

(B) Projection exchange constant extraction at n_max_fock >= 3.

(C) Per-edge graph VP decomposition and comparison with continuum triples.

(D) Convergence analysis: graph VP vs continuum VP as n_max increases.

(E) Self-energy bridge: structural zero at n_ext=0.

References
----------
- GeoVac graph_qed_vertex.py (graph-native VP)
- GeoVac graph_qed_propagator.py (graph electron propagator)
- GeoVac graph_qed_photon.py (graph photon propagator, L1)
- GeoVac qed_vacuum_polarization.py (continuum VP, Seeley-DeWitt)
- GeoVac qed_self_energy.py (continuum self-energy, structural zero)
- GeoVac qed_vertex.py (continuum two-loop, vertex selection rules)
- GeoVac Paper 18 (exchange constant taxonomy)
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Rational, pi, sqrt, simplify

__all__ = [
    # Part A
    "continuum_vp_mode_contribution",
    "continuum_vp_truncated",
    # Part B
    "compute_projection_constant",
    "classify_projection_constant",
    # Part C
    "graph_vp_edge_decomposition",
    "continuum_vp_mode_decomposition",
    "mode_ratio_analysis",
    # Part D
    "convergence_table",
    # Part E
    "graph_self_energy_n0",
    "self_energy_structural_zero_bridge",
    # Driver
    "run_bridge_analysis",
]


# ---------------------------------------------------------------------------
# Continuum spectrum helpers (CH convention: n >= 0)
# ---------------------------------------------------------------------------

def _lambda_n_ch(n: int) -> Rational:
    """Absolute Dirac eigenvalue |lambda_n| = n + 3/2 (CH convention)."""
    return Rational(2 * n + 3, 2)


def _g_n_dirac(n: int) -> int:
    """Full Dirac degeneracy g_n = 2(n+1)(n+2)."""
    return 2 * (n + 1) * (n + 2)


def _mu_q(q: int) -> int:
    """Hodge-1 Laplacian eigenvalue mu_q = q(q+2), q >= 1."""
    return q * (q + 2)


def _d_q_transverse(q: int) -> int:
    """Transverse photon degeneracy d_q^T = q(q+2)."""
    return q * (q + 2)


def _vertex_allowed(n1: int, n2: int, q: int) -> bool:
    """SO(4) vertex selection rule: triangle + parity + q >= 1."""
    if q < 1:
        return False
    if q < abs(n1 - n2) or q > n1 + n2:
        return False
    if (n1 + n2 + q) % 2 == 0:
        return False
    return True


def _so4_channel_count(n1: int, n2: int, q: int) -> int:
    """Count SO(4) vector harmonic channels (0, 1, or 2).

    Uses the SU(2)_L x SU(2)_R triangle check consistent with
    qed_vertex.so4_channel_count.
    """
    if not _vertex_allowed(n1, n2, q):
        return 0

    j1_L = Fraction(n1 + 1, 2)
    j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2)
    j2_R = Fraction(n2 + 1, 2)

    count = 0
    # Component A: ((q+1)/2, (q-1)/2)
    jg_L_A = Fraction(q + 1, 2)
    jg_R_A = Fraction(q - 1, 2)
    if (jg_R_A >= 0
            and abs(j1_L - jg_L_A) <= j2_L <= j1_L + jg_L_A
            and abs(j1_R - jg_R_A) <= j2_R <= j1_R + jg_R_A):
        count += 1

    # Component B: ((q-1)/2, (q+1)/2)
    jg_L_B = Fraction(q - 1, 2)
    jg_R_B = Fraction(q + 1, 2)
    if (jg_L_B >= 0
            and abs(j1_L - jg_L_B) <= j2_L <= j1_L + jg_L_B
            and abs(j1_R - jg_R_B) <= j2_R <= j1_R + jg_R_B):
        count += 1

    return count


# ---------------------------------------------------------------------------
# Convention bridge
# ---------------------------------------------------------------------------

def _n_ch_from_n_fock(n_fock: int) -> int:
    return n_fock - 1

def _n_fock_from_n_ch(n_ch: int) -> int:
    return n_ch + 1

def _n_max_ch_from_n_max_fock(n_max_fock: int) -> int:
    return n_max_fock - 1


# =========================================================================
# Part A: Truncated continuum VP
# =========================================================================

def continuum_vp_mode_contribution(
    n1_ch: int,
    n2_ch: int,
    q: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> Rational:
    """Single (n1, n2, q) continuum VP mode contribution.

    Pi_cont(n1,n2,q) = W(n1,n2,q) * g(n1) * g(n2) * d_T(q)
                       / (|lam(n1)|^{2*s_e} * |lam(n2)|^{2*s_e} * mu(q)^{s_gamma})

    This is RATIONAL at every mode (no pi). The transcendental content
    comes from summing to infinity, not from individual modes.
    """
    W = _so4_channel_count(n1_ch, n2_ch, q)
    if W == 0:
        return Rational(0)

    g1 = _g_n_dirac(n1_ch)
    g2 = _g_n_dirac(n2_ch)
    d_T = _d_q_transverse(q)
    lam1 = _lambda_n_ch(n1_ch)
    lam2 = _lambda_n_ch(n2_ch)
    mu = _mu_q(q)

    return Rational(W * g1 * g2 * d_T) / (lam1 ** (2 * s_e) * lam2 ** (2 * s_e) * mu ** s_gamma)


def continuum_vp_truncated(
    n_max_ch: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> Rational:
    """Continuum VP spectral sum truncated to n_max_ch modes."""
    total = Rational(0)
    for n1 in range(n_max_ch + 1):
        for n2 in range(n_max_ch + 1):
            q_hi = n1 + n2
            for q in range(max(1, abs(n1 - n2)), q_hi + 1):
                total += continuum_vp_mode_contribution(n1, n2, q, s_e, s_gamma)
    return total


# =========================================================================
# Part B: Projection exchange constant
# =========================================================================

def compute_projection_constant(
    n_max_fock: int = 3,
    graph_trace: Optional[Rational] = None,
) -> Dict[str, Any]:
    """Projection exchange constant C = Pi_continuum_trunc / Tr(Pi_graph).

    At finite truncation, BOTH quantities are rational, so C is rational.
    Pi enters only in the infinite-sum limit.

    NOTE: At n_max_fock=2 the continuum VP is zero (the only allowed
    triple (1,1,1) has SO(4) channel count W=0), so the ratio is zero.
    Use n_max_fock >= 3 for a meaningful projection constant.
    """
    n_max_ch = _n_max_ch_from_n_max_fock(n_max_fock)

    if graph_trace is None:
        from geovac.graph_qed_vertex import compute_vacuum_polarization
        vp = compute_vacuum_polarization(n_max=n_max_fock, t=Rational(0), exact=True)
        graph_trace = sp.nsimplify(vp['Pi'].trace(), rational=True)

    continuum_trunc = continuum_vp_truncated(n_max_ch)

    C = None
    if graph_trace is not None and graph_trace != 0:
        C = simplify(continuum_trunc / graph_trace)

    is_rational = isinstance(C, (sp.Rational, sp.Integer)) if C is not None else None

    return {
        'n_max_fock': n_max_fock,
        'n_max_ch': n_max_ch,
        'graph_trace': str(graph_trace),
        'graph_trace_float': float(graph_trace) if graph_trace is not None else None,
        'continuum_truncated': str(continuum_trunc),
        'continuum_truncated_float': float(continuum_trunc),
        'projection_constant': str(C) if C is not None else None,
        'projection_constant_float': float(C) if C is not None else None,
        'is_rational': is_rational,
        'continuum_is_zero': continuum_trunc == 0,
        'graph_is_zero': graph_trace == 0 if graph_trace is not None else True,
    }


def classify_projection_constant(n_max_fock: int = 3) -> Dict[str, str]:
    """Classify the projection exchange constant in Paper 18's taxonomy.

    The central finding: at finite truncation, the projection constant
    is RATIONAL (no pi).  The graph VP is richer than the continuum VP
    because the CG projection opens channels that W=0 in the SO(4)
    channel counting.

    Pi enters only in:
    1. Infinite-sum regularization (Hurwitz zeta -> pi^{even})
    2. Vol(S^3) = 2*pi^2 normalization for dimensionful quantities
    """
    result = compute_projection_constant(n_max_fock)
    classification = {}

    if result['continuum_is_zero']:
        classification['status'] = (
            'CONTINUUM VP IS ZERO at this truncation. The only allowed '
            'triples have SO(4) channel count W=0. The graph VP is '
            'nonzero because the CG projection opens additional channels.'
        )
        classification['paper_18_tier'] = 'N/A (zero continuum)'
    elif result['is_rational']:
        classification['status'] = (
            'RATIONAL: at finite truncation, both graph and continuum '
            'VPs are rational. The projection constant is rational with '
            'no pi content.'
        )
        classification['pi_entry_point'] = (
            'Pi enters ONLY when: (1) the continuum sum is extended to '
            'infinity and regularized via Hurwitz zeta (producing pi^{even} '
            'from Bernoulli numbers), or (2) the VP is normalized by '
            'Vol(S^3) = 2*pi^2 for dimensionful physical quantities.'
        )
        classification['paper_18_tier'] = 'CALIBRATION'
        classification['mechanism'] = (
            'The pi is from the Riemannian volume form of S^3. '
            'The graph is pi-free at every finite truncation. '
            'Transcendental content is injected exclusively by the '
            'continuum limit (n_max -> infinity) and/or manifold measure.'
        )
    else:
        classification['status'] = f'NOT RATIONAL: {result["projection_constant"]}'
        classification['paper_18_tier'] = 'UNKNOWN'

    classification['projection_constant_value'] = result['projection_constant']
    classification['graph_trace'] = result['graph_trace']
    classification['continuum_truncated'] = result['continuum_truncated']

    return classification


# =========================================================================
# Part C: Mode-by-mode decomposition
# =========================================================================

def graph_vp_edge_decomposition(
    n_max_fock: int = 2,
) -> Dict[str, Any]:
    """Decompose graph VP trace by Fock edge.

    Each edge e of the Fock graph contributes Pi[e,e] = Tr(V_e^T G V_e G).
    At t=0 with diagonal G, this is a rational number.

    The per-edge numerator is universally 32 (= 2^5) at n_max_fock=2,3.
    The denominator encodes the Dirac eigenvalue products from the
    two electron propagators.
    """
    from geovac.graph_qed_vertex import (
        build_projection_matrix,
        build_vertex_tensor,
        vertex_tensor_to_matrices,
    )
    from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
    from geovac.graph_qed_photon import build_fock_graph

    P, dirac_labels, fock_states = build_projection_matrix(n_max_fock)
    fock_data = build_fock_graph(n_max_fock)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max_fock, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    op = DiracGraphOperator(n_max=n_max_fock, t=Rational(0))
    G_e, _ = electron_propagator(op, exact=True)

    edge_data = []
    total = Rational(0)
    for e_idx, (v1, v2) in enumerate(fock_data.edges):
        s1 = fock_data.states[v1]
        s2 = fock_data.states[v2]
        Ve = V_mats[e_idx]
        Pi_ee = sp.nsimplify((Ve.T * G_e * Ve * G_e).trace(), rational=True)
        total += Pi_ee

        # Classify the edge by shell information
        n1_fock, l1, m1 = s1
        n2_fock, l2, m2 = s2
        edge_type = 'inter-shell' if n1_fock != n2_fock else 'intra-shell'

        edge_data.append({
            'edge_idx': e_idx,
            'fock_node_1': list(s1),
            'fock_node_2': list(s2),
            'edge_type': edge_type,
            'n_fock_1': n1_fock,
            'n_fock_2': n2_fock,
            'Pi_ee': str(Pi_ee),
            'Pi_ee_float': float(Pi_ee),
            'numerator': int(sp.Rational(Pi_ee).p),
            'denominator': int(sp.Rational(Pi_ee).q),
        })

    return {
        'n_max_fock': n_max_fock,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'edges': edge_data,
        'total_trace': str(total),
        'total_trace_float': float(total),
    }


def _enumerate_allowed_triples(n_max_ch: int) -> List[Tuple[int, int, int]]:
    """Enumerate all (n1, n2, q) triples allowed by vertex selection rule."""
    triples = []
    for n1 in range(n_max_ch + 1):
        for n2 in range(n_max_ch + 1):
            for q in range(max(1, abs(n1 - n2)), n1 + n2 + 1):
                if _vertex_allowed(n1, n2, q):
                    triples.append((n1, n2, q))
    return triples


def continuum_vp_mode_decomposition(
    n_max_ch: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> Dict[str, Any]:
    """Decompose continuum VP by (n1, n2, q) mode triples."""
    triples = _enumerate_allowed_triples(n_max_ch)

    contributions = []
    total = Rational(0)
    for n1, n2, q in triples:
        W = _so4_channel_count(n1, n2, q)
        contrib = continuum_vp_mode_contribution(n1, n2, q, s_e, s_gamma)
        total += contrib
        contributions.append({
            'n1_ch': n1,
            'n2_ch': n2,
            'q': q,
            'W': W,
            'allowed': True,
            'contribution': str(contrib),
            'contribution_float': float(contrib),
        })

    return {
        'n_max_ch': n_max_ch,
        'n_allowed_triples': len(triples),
        'n_nonzero_W': sum(1 for c in contributions if c['W'] > 0),
        'n_zero_W_but_allowed': sum(1 for c in contributions if c['W'] == 0),
        'contributions': contributions,
        'total': str(total),
        'total_float': float(total),
    }


def mode_ratio_analysis(n_max_fock: int = 3) -> Dict[str, Any]:
    """Compare graph VP total with continuum VP total.

    The graph and continuum use DIFFERENT decompositions:
    - Graph: photon edge index (E_fock edges)
    - Continuum: (n1, n2, q) mode triples

    These are different labelings. The meaningful comparison is at the
    TOTAL (trace) level.

    Key finding: at n_max_fock=2, the graph VP = 224/75 while the
    continuum VP = 0.  This means the graph is strictly richer.
    """
    n_max_ch = _n_max_ch_from_n_max_fock(n_max_fock)

    from geovac.graph_qed_vertex import compute_vacuum_polarization
    vp = compute_vacuum_polarization(n_max=n_max_fock, t=Rational(0), exact=True)
    graph_trace = sp.nsimplify(vp['Pi'].trace(), rational=True)

    continuum_total = continuum_vp_truncated(n_max_ch)

    ratio = None
    if graph_trace != 0:
        ratio = simplify(continuum_total / graph_trace)

    # Also compute the "graph excess" = graph - continuum
    excess = simplify(graph_trace - continuum_total) if graph_trace is not None else None

    cont_decomp = continuum_vp_mode_decomposition(n_max_ch)

    return {
        'n_max_fock': n_max_fock,
        'n_max_ch': n_max_ch,
        'graph_trace': str(graph_trace),
        'graph_trace_float': float(graph_trace),
        'continuum_total': str(continuum_total),
        'continuum_total_float': float(continuum_total),
        'ratio_cont_over_graph': str(ratio) if ratio is not None else None,
        'ratio_float': float(ratio) if ratio is not None else None,
        'ratio_is_rational': isinstance(ratio, (sp.Rational, sp.Integer)) if ratio is not None else None,
        'graph_excess': str(excess) if excess is not None else None,
        'graph_excess_float': float(excess) if excess is not None else None,
        'n_continuum_triples_with_nonzero_W': cont_decomp['n_nonzero_W'],
        'n_continuum_triples_allowed_but_W0': cont_decomp['n_zero_W_but_allowed'],
        'n_graph_edges': vp['E_fock'],
    }


# =========================================================================
# Part D: Convergence
# =========================================================================

def _continuum_vp_full_estimate(
    n_max_sum: int = 200,
    s_e: int = 2,
    s_gamma: int = 1,
) -> float:
    """Estimate the full (untruncated) continuum VP sum numerically."""
    total = 0.0
    for n1 in range(n_max_sum + 1):
        g1 = 2.0 * (n1 + 1) * (n1 + 2)
        lam1 = n1 + 1.5
        for n2 in range(n_max_sum + 1):
            g2 = 2.0 * (n2 + 1) * (n2 + 2)
            lam2 = n2 + 1.5
            for q in range(max(1, abs(n1 - n2)), n1 + n2 + 1):
                if (n1 + n2 + q) % 2 == 0:
                    continue
                W = _so4_channel_count(n1, n2, q)
                if W == 0:
                    continue
                d_T = q * (q + 2.0)
                mu = q * (q + 2.0)
                total += (W * g1 * g2 * d_T
                          / (lam1 ** (2 * s_e) * lam2 ** (2 * s_e)
                             * mu ** s_gamma))
    return total


def convergence_table(
    n_max_values: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """Build a convergence table for graph VP and continuum VP vs n_max.

    Shows how the ratio C = Pi_cont / Tr(Pi_graph) evolves with truncation,
    and what fraction of the full continuum VP is captured.
    """
    if n_max_values is None:
        n_max_values = [1, 2, 3]

    full_estimate = _continuum_vp_full_estimate(n_max_sum=50)

    rows = []
    for n_max_fock in n_max_values:
        n_max_ch = _n_max_ch_from_n_max_fock(n_max_fock)

        # Graph VP
        from geovac.graph_qed_vertex import compute_vacuum_polarization
        if n_max_fock <= 3:
            vp = compute_vacuum_polarization(n_max=n_max_fock, t=Rational(0), exact=True)
            g_trace_exact = sp.nsimplify(vp['Pi'].trace(), rational=True)
            g_trace = float(g_trace_exact)
        else:
            vp = compute_vacuum_polarization(n_max=n_max_fock, t=Rational(0), exact=False)
            g_trace = float(np.trace(vp['Pi']))
            g_trace_exact = None

        # Continuum VP truncated
        if n_max_fock <= 4:
            c_trunc_exact = continuum_vp_truncated(n_max_ch)
            c_trunc = float(c_trunc_exact)
        else:
            c_trunc = _continuum_vp_full_estimate(n_max_sum=n_max_ch)
            c_trunc_exact = None

        ratio = c_trunc / g_trace if g_trace != 0 else None
        frac = c_trunc / full_estimate if full_estimate != 0 else None

        rows.append({
            'n_max_fock': n_max_fock,
            'n_max_ch': n_max_ch,
            'graph_trace': g_trace,
            'graph_trace_exact': str(g_trace_exact) if g_trace_exact is not None else None,
            'continuum_truncated': c_trunc,
            'continuum_truncated_exact': str(c_trunc_exact) if c_trunc_exact is not None else None,
            'projection_constant': ratio,
            'fraction_of_full_continuum': frac,
        })

    return {
        'full_continuum_estimate_n50': full_estimate,
        'rows': rows,
    }


# =========================================================================
# Part E: Self-energy bridge
# =========================================================================

def graph_self_energy_n0(n_max_fock: int = 2) -> Dict[str, Any]:
    """Check ground-state vertex couplings in the graph.

    The continuum self-energy has Sigma(n_ext=0) = 0 by vertex parity:
    at n_CH=0, the selection rule requires n_int + q odd AND q = n_int
    (triangle), giving 2*n_int = odd, which is impossible.

    The GRAPH has ground-state vertex couplings because the CG projection
    maps Dirac ground states onto Fock nodes that ARE connected by edges.
    """
    from geovac.graph_qed_vertex import (
        build_projection_matrix,
        build_vertex_tensor,
    )
    from geovac.graph_qed_photon import build_fock_graph

    P, dirac_labels, fock_states = build_projection_matrix(n_max_fock)
    fock_data = build_fock_graph(n_max_fock)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max_fock, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )

    ground_indices = [
        i for i, lab in enumerate(dirac_labels) if lab.n_fock == 1
    ]

    ground_couplings = [
        {'a': a, 'b': b, 'e': e, 'val': str(val)}
        for a, b, e, val in entries
        if a in ground_indices
    ]

    return {
        'n_max_fock': n_max_fock,
        'n_ground_states': len(ground_indices),
        'ground_state_labels': [
            f"(n={dirac_labels[i].n_fock}, k={dirac_labels[i].kappa}, "
            f"2mj={dirac_labels[i].two_m_j})"
            for i in ground_indices
        ],
        'n_ground_vertex_couplings': len(ground_couplings),
        'has_ground_couplings': len(ground_couplings) > 0,
        'continuum_structural_zero': True,
        'continuum_mechanism': (
            'At n_CH=0: n_int + q must be odd AND q = n_int (triangle). '
            '=> 2*n_int = odd => impossible. Zero allowed triples.'
        ),
    }


def self_energy_structural_zero_bridge(n_max_fock: int = 2) -> Dict[str, Any]:
    """Compute graph self-energy for ground-state Dirac labels.

    Sigma_graph[a] = sum_{e1,e2,b} V_{e1}[a,b] G_e[b,b] G_gamma[e1,e2] V_{e2}[b,a]
    """
    from geovac.graph_qed_vertex import (
        build_projection_matrix,
        build_vertex_tensor,
        vertex_tensor_to_matrices,
    )
    from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
    from geovac.graph_qed_photon import build_fock_graph, compute_photon_propagator

    P, dirac_labels, fock_states = build_projection_matrix(n_max_fock)
    fock_data = build_fock_graph(n_max_fock)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max_fock, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    op = DiracGraphOperator(n_max=n_max_fock, t=Rational(0))
    G_e, _ = electron_propagator(op, exact=True)

    photon = compute_photon_propagator(n_max_fock, exact=True)
    G_gamma = photon.G_gamma

    ground_indices = [
        i for i, lab in enumerate(dirac_labels) if lab.n_fock == 1
    ]

    ground_self_energies = {}
    for a in ground_indices:
        sigma_a = Rational(0)
        for b in range(N_dirac):
            g_b = G_e[b, b]
            for e1 in range(E_fock):
                v_ab = V_mats[e1][a, b]
                if v_ab == 0:
                    continue
                for e2 in range(E_fock):
                    v_ba = V_mats[e2][b, a]
                    if v_ba == 0:
                        continue
                    g_gg = G_gamma[e1, e2] if G_gamma is not None else (
                        Rational(1) if e1 == e2 else Rational(0)
                    )
                    sigma_a += v_ab * g_b * g_gg * v_ba

        sigma_a = sp.nsimplify(sigma_a, rational=True)
        lab = dirac_labels[a]
        key = f"(n={lab.n_fock}, k={lab.kappa}, 2mj={lab.two_m_j})"
        ground_self_energies[key] = {
            'value': str(sigma_a),
            'value_float': float(sigma_a),
            'is_zero': sigma_a == 0,
        }

    all_zero = all(v['is_zero'] for v in ground_self_energies.values())

    return {
        'n_max_fock': n_max_fock,
        'ground_state_self_energies': ground_self_energies,
        'all_ground_state_zero': all_zero,
        'continuum_structural_zero': True,
        'structural_zero_preserved': all_zero,
        'interpretation': (
            'PRESERVED: graph self-energy vanishes at ground state, '
            'same structural zero as continuum.'
            if all_zero else
            'BROKEN: graph self-energy is NONZERO at ground state. '
            'The CG projection from Dirac to Fock basis opens '
            'couplings that the continuum vertex parity selection '
            'rule forbids. The structural zero is a continuum '
            'selection rule, not preserved by the finite graph.'
        ),
    }


# =========================================================================
# Driver
# =========================================================================

def run_bridge_analysis(
    n_max_fock: int = 2,
    save_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run the complete graph-to-continuum bridge analysis."""
    results: Dict[str, Any] = {
        'module': 'geovac.graph_qed_continuum_bridge',
        'n_max_fock': n_max_fock,
    }

    n_max_ch = _n_max_ch_from_n_max_fock(n_max_fock)

    # Part A
    cont_trunc = continuum_vp_truncated(n_max_ch)
    results['part_a'] = {
        'n_max_ch': n_max_ch,
        'continuum_vp_truncated': str(cont_trunc),
        'continuum_vp_truncated_float': float(cont_trunc),
    }

    # Part B: projection constant at requested n_max AND at n_max=3
    results['part_b_at_requested'] = compute_projection_constant(n_max_fock)
    results['part_b_classification'] = classify_projection_constant(n_max_fock)
    if n_max_fock < 3:
        results['part_b_at_nmax3'] = compute_projection_constant(3)
        results['part_b_classification_nmax3'] = classify_projection_constant(3)

    # Part C: edge decomposition + mode comparison
    results['part_c_graph_edges'] = graph_vp_edge_decomposition(n_max_fock)
    results['part_c_continuum_modes'] = continuum_vp_mode_decomposition(n_max_ch)
    results['part_c_ratio'] = mode_ratio_analysis(
        max(n_max_fock, 3)  # need n_max >= 3 for nonzero continuum
    )

    # Part D: convergence
    nmax_list = sorted(set([1, 2, 3] + ([n_max_fock] if n_max_fock <= 3 else [])))
    results['part_d'] = convergence_table(n_max_values=nmax_list)

    # Part E: self-energy bridge
    results['part_e_couplings'] = graph_self_energy_n0(n_max_fock)
    results['part_e_self_energy'] = self_energy_structural_zero_bridge(n_max_fock)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open("w") as f:
            json.dump(results, f, indent=2, default=str)

    return results


if __name__ == "__main__":
    output = Path(__file__).parent.parent / "debug" / "data" / "gn7_continuum_bridge.json"
    results = run_bridge_analysis(n_max_fock=2, save_path=output)
    print(f"Saved: {output}")

    print("\n=== KEY FINDINGS ===")
    pc = results['part_b_at_requested']
    print(f"Graph VP trace (n_max_fock=2): {pc['graph_trace']}")
    print(f"Continuum VP (n_max_ch=1):     {pc['continuum_truncated']}")
    print(f"Projection constant:           {pc['projection_constant']}")

    if 'part_b_at_nmax3' in results:
        pc3 = results['part_b_at_nmax3']
        print(f"\nGraph VP trace (n_max_fock=3): {pc3['graph_trace']}")
        print(f"Continuum VP (n_max_ch=2):     {pc3['continuum_truncated']}")
        print(f"Projection constant:           {pc3['projection_constant']}")
        print(f"Is rational:                   {pc3['is_rational']}")

    cls = results.get('part_b_classification_nmax3', results['part_b_classification'])
    print(f"\nPaper 18 tier: {cls.get('paper_18_tier', 'N/A')}")

    se = results['part_e_self_energy']
    print(f"\nSelf-energy structural zero preserved: {se['structural_zero_preserved']}")
    for k, v in se['ground_state_self_energies'].items():
        print(f"  {k}: Sigma = {v['value']} ({v['value_float']:.6f})")
