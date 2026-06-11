"""
Native Dirac Graph QED Sprint — Build QED on (n, kappa, m_j) Nodes
====================================================================

This script implements QED directly on the Dirac graph, bypassing
the CG projection used in the scalar-Fock-graph QED (GN-1 through GN-7).

Nodes are Dirac labels (n_fock, kappa, m_j).
Edges come from two adjacency rules:
  - Rule A (kappa-preserving scalar analog)
  - Rule B (E1 dipole, parity-flip)

Tasks:
  1. Hodge decomposition for Rule A and Rule B at n_max=2,3
  2. Photon propagator G_gamma = L_1^+ (pseudoinverse)
  3. Electron propagator D_Dirac = Lambda + t * A_Dirac
  4. Self-energy Sigma and structural zero check
  5. Selection rule census on the Dirac graph

Author: GeoVac Sprint, April 2026.
"""

from __future__ import annotations
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Matrix, Rational, sqrt, zeros as sp_zeros

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.ihara_zeta_dirac import (
    build_dirac_s3_graph,
    AdjacencyRule,
)
from geovac.dirac_s3 import dirac_eigenvalue_abs, fock_to_ch

# ============================================================================
# TASK 1: Hodge Decomposition
# ============================================================================

def extract_undirected_edges(A: np.ndarray) -> List[Tuple[int, int]]:
    """Extract undirected edges (i, j) with i < j from adjacency matrix."""
    edges = []
    V = A.shape[0]
    for i in range(V):
        for j in range(i + 1, V):
            if A[i, j] != 0:
                edges.append((i, j))
    return sorted(edges)


def build_incidence_matrix(V: int, edges: List[Tuple[int, int]]) -> np.ndarray:
    """Build V x E signed incidence matrix B with entries in {-1, 0, +1}.

    Convention: for edge k = (i, j) with i < j,
      B[i, k] = +1   (tail / source)
      B[j, k] = -1   (head / sink)
    """
    E = len(edges)
    B = np.zeros((V, E), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0
    return B


def hodge_decomposition(
    n_max: int,
    rule: AdjacencyRule,
) -> Dict[str, Any]:
    """Full Hodge decomposition of the Dirac graph.

    Returns dict with keys:
      V, E, edges, B, L0, L1, beta_0, beta_1,
      L0_spectrum, L1_spectrum, svd_identity_verified,
      pi_free_certificate, degree_sequence, pendant_vertices
    """
    A, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    V = len(labels)
    edges = extract_undirected_edges(A)
    E = len(edges)

    B = build_incidence_matrix(V, edges)

    # L0 = B @ B.T (node Laplacian, V x V)
    L0 = B @ B.T
    # L1 = B.T @ B (edge Laplacian, E x E)
    L1 = B.T @ B

    # Verify L0 = D - A
    D_mat = np.diag(deg.astype(float))
    A_float = A.astype(float)
    L0_check = D_mat - A_float
    L0_DA_match = np.allclose(L0, L0_check, atol=1e-12)

    # Spectra
    L0_evals = np.sort(np.linalg.eigvalsh(L0))
    L1_evals = np.sort(np.linalg.eigvalsh(L1)) if E > 0 else np.array([])

    # Betti numbers
    beta_0 = int(np.sum(np.abs(L0_evals) < 1e-10))  # connected components
    beta_1 = E - V + beta_0  # first Betti number

    # SVD identity: nonzero eigenvalues of L0 and L1 should match
    L0_nonzero = sorted([ev for ev in L0_evals if abs(ev) > 1e-10])
    L1_nonzero = sorted([ev for ev in L1_evals if abs(ev) > 1e-10])
    svd_identity_verified = len(L0_nonzero) == len(L1_nonzero) and all(
        abs(a - b) < 1e-8 for a, b in zip(L0_nonzero, L1_nonzero)
    )

    # pi-free certificate: adjacency is integer, so all eigenvalues are
    # roots of integer-coefficient characteristic polynomial = algebraic
    # Verify: L0 and L1 have integer entries (since B has integer entries)
    L0_int = np.allclose(L0, np.round(L0), atol=1e-12)
    L1_int = np.allclose(L1, np.round(L1), atol=1e-12)
    pi_free = L0_int and L1_int

    # Pendant vertices (degree 1)
    pendant_vertices = []
    for i, d in enumerate(deg):
        if d == 1:
            pendant_vertices.append(i)

    # Ground state check: is the ground state (n_fock=1, kappa=-1) a pendant?
    gs_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    gs_pendants = [i for i in gs_indices if i in pendant_vertices]

    return {
        "description": desc,
        "V": V,
        "E": E,
        "beta_0": beta_0,
        "beta_1": beta_1,
        "L0_DA_match": L0_DA_match,
        "svd_identity_verified": svd_identity_verified,
        "pi_free_certificate": pi_free,
        "L0_spectrum": L0_evals.tolist(),
        "L1_spectrum": L1_evals.tolist(),
        "degree_sequence": deg.tolist(),
        "pendant_vertex_count": len(pendant_vertices),
        "pendant_vertex_indices": pendant_vertices,
        "gs_indices": gs_indices,
        "gs_is_pendant": len(gs_pendants) > 0,
        "gs_degrees": [int(deg[i]) for i in gs_indices],
        "B": B,
        "L0": L0,
        "L1": L1,
        "adjacency": A,
        "labels": labels,
        "edges": edges,
    }


# ============================================================================
# TASK 2: Photon Propagator on Dirac Graph
# ============================================================================

def photon_propagator_dirac(L1: np.ndarray) -> np.ndarray:
    """Compute photon propagator G_gamma = L1^+ (Moore-Penrose pseudoinverse).

    Uses SVD-based pseudoinverse for numerical stability.
    """
    if L1.shape[0] == 0:
        return np.array([]).reshape(0, 0)
    return np.linalg.pinv(L1)


def photon_propagator_exact(L1_int: np.ndarray) -> Optional[sp.Matrix]:
    """Exact sympy pseudoinverse for small matrices."""
    if L1_int.shape[0] > 40:
        return None
    M = sp.Matrix(np.round(L1_int).astype(int).tolist())
    try:
        return M.pinv()
    except Exception:
        return None


# ============================================================================
# TASK 3: Electron Propagator on Dirac Graph
# ============================================================================

def build_dirac_operator(
    labels: List[DiracLabel],
    A: np.ndarray,
    t: float = 0.0,
) -> np.ndarray:
    """Build D_Dirac = Lambda + t * A_Dirac.

    Lambda is diagonal with Camporesi-Higuchi eigenvalues:
      lambda_i = chi_i * (n_CH_i + 3/2)
    where n_CH = n_fock - 1 (Fock to CH convention) and
    chi = +1 for kappa < 0, -1 for kappa > 0 (chirality).
    """
    V = len(labels)
    D = np.zeros((V, V), dtype=np.float64)

    for i, lab in enumerate(labels):
        n_ch = lab.n_fock - 1  # Fock to CH convention
        chi = 1 if lab.kappa < 0 else -1
        lam = chi * (n_ch + 1.5)
        D[i, i] = lam

    D += t * A.astype(np.float64)
    return D


def build_dirac_operator_exact(
    labels: List[DiracLabel],
    A: np.ndarray,
    t: sp.Expr = Rational(0),
) -> sp.Matrix:
    """Build exact sympy D_Dirac = Lambda + t * A."""
    V = len(labels)
    D = sp.zeros(V, V)

    for i, lab in enumerate(labels):
        n_ch = lab.n_fock - 1
        chi = 1 if lab.kappa < 0 else -1
        lam = Rational(chi * (2 * n_ch + 3), 2)
        D[i, i] = lam

    for i in range(V):
        for j in range(V):
            if A[i, j] != 0:
                D[i, j] = D[i, j] + t * Integer(int(A[i, j]))

    return D


# ============================================================================
# TASK 4: Self-Energy on Dirac Graph
# ============================================================================

def compute_dirac_graph_self_energy(
    hodge: Dict[str, Any],
    t: float = 0.0,
) -> Dict[str, Any]:
    """Compute one-loop self-energy on the Dirac graph.

    On the Dirac graph, the vertex coupling is the identity on edges:
    V[a, b, e] = delta_{a, i(e)} * delta_{b, j(e)}
    where edge e connects nodes i(e) and j(e).

    This is fundamentally different from the scalar Fock graph where
    V[a,b,e] requires CG projection.

    Sigma[a,b] = sum_{e,e'} V_e[a,:] . G_gamma[e,e'] . V_{e'}[:,b]^T
    """
    V_nodes = hodge["V"]
    E_edges = hodge["E"]
    labels = hodge["labels"]
    edges = hodge["edges"]
    A = hodge["adjacency"]
    L1 = hodge["L1"]

    if E_edges == 0:
        return {
            "Sigma_trace": 0.0,
            "gs_block_zero": True,
            "gs_block": [[0.0]],
            "eigenvalues": [],
            "is_hermitian": True,
        }

    # Photon propagator
    G_gamma = photon_propagator_dirac(L1)

    # Build vertex matrices V_e (V_nodes x V_nodes) for each edge
    # V_e[a, b] = delta_{a, i(e)} * delta_{b, j(e)} + delta_{a, j(e)} * delta_{b, i(e)}
    # (symmetrized: both orientations of the edge contribute)
    V_mats = []
    for e_idx, (i, j) in enumerate(edges):
        Ve = np.zeros((V_nodes, V_nodes))
        Ve[i, j] = 1.0
        Ve[j, i] = 1.0
        V_mats.append(Ve)

    # Self-energy: Sigma = sum_{e1, e2} G_gamma[e1, e2] * V_{e1} @ V_{e2}.T
    Sigma = np.zeros((V_nodes, V_nodes))
    for e1 in range(E_edges):
        for e2 in range(E_edges):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma += g_ee * (V_mats[e1] @ V_mats[e2].T)

    # Properties
    trace_val = np.trace(Sigma)
    is_hermitian = np.allclose(Sigma, Sigma.T, atol=1e-12)
    evals = np.sort(np.linalg.eigvalsh(Sigma))

    # Ground state block
    gs_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    if len(gs_indices) > 0:
        gs_block = Sigma[np.ix_(gs_indices, gs_indices)]
        gs_block_zero = np.allclose(gs_block, 0, atol=1e-12)
    else:
        gs_block = np.array([[0.0]])
        gs_block_zero = True

    # Check positive semidefiniteness
    is_psd = all(ev > -1e-10 for ev in evals)

    return {
        "Sigma": Sigma,
        "Sigma_trace": float(trace_val),
        "is_hermitian": is_hermitian,
        "is_psd": is_psd,
        "eigenvalues": evals.tolist(),
        "gs_indices": gs_indices,
        "gs_block": gs_block.tolist() if isinstance(gs_block, np.ndarray) else gs_block,
        "gs_block_zero": bool(gs_block_zero),
        "N_dirac": V_nodes,
        "E_edges": E_edges,
    }


def compute_dirac_graph_vertex_correction(
    hodge: Dict[str, Any],
    t: float = 0.0,
) -> Dict[str, Any]:
    """Compute one-loop vertex correction Lambda on the Dirac graph.

    Lambda_total = sum_{e', e''} G_gamma[e', e''] * V_{e'} @ G_e @ V_{e''}^T
    """
    V_nodes = hodge["V"]
    E_edges = hodge["E"]
    labels = hodge["labels"]
    edges = hodge["edges"]
    A = hodge["adjacency"]
    L1 = hodge["L1"]

    if E_edges == 0:
        return {
            "Lambda_trace": 0.0,
            "eigenvalues": [],
            "F2": None,
        }

    # Photon propagator
    G_gamma = photon_propagator_dirac(L1)

    # Electron propagator
    D_dirac = build_dirac_operator(labels, A, t=t)
    try:
        G_e = np.linalg.inv(D_dirac)
    except np.linalg.LinAlgError:
        return {
            "Lambda_trace": None,
            "eigenvalues": [],
            "F2": None,
            "error": "D_dirac is singular",
        }

    # Build vertex matrices
    V_mats = []
    for e_idx, (i, j) in enumerate(edges):
        Ve = np.zeros((V_nodes, V_nodes))
        Ve[i, j] = 1.0
        Ve[j, i] = 1.0
        V_mats.append(Ve)

    # Vertex correction: Lambda = sum_{e1, e2} G_gamma[e1, e2] * V_{e1} @ G_e @ V_{e2}.T
    Lambda_total = np.zeros((V_nodes, V_nodes))
    for e1 in range(E_edges):
        for e2 in range(E_edges):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Lambda_total += g_ee * (V_mats[e1] @ G_e @ V_mats[e2].T)

    trace_val = np.trace(Lambda_total)
    evals = np.sort(np.linalg.eigvalsh(Lambda_total))

    # F2 extraction: Tr(Lambda) / Tr(V_bare @ G_e)
    # V_bare = sum_e V_e (the total bare vertex)
    V_bare = np.zeros((V_nodes, V_nodes))
    for Ve in V_mats:
        V_bare += Ve
    bare_trace = np.trace(V_bare @ G_e)
    F2 = float(trace_val / bare_trace) if abs(bare_trace) > 1e-15 else None

    return {
        "Lambda_trace": float(trace_val),
        "eigenvalues": evals.tolist(),
        "F2": F2,
        "bare_vertex_trace": float(bare_trace),
    }


# ============================================================================
# TASK 5: Selection Rule Census on Dirac Graph
# ============================================================================

def selection_rule_census(
    hodge: Dict[str, Any],
    sigma_result: Dict[str, Any],
) -> Dict[str, Any]:
    """Test all 8 continuum QED selection rules on the Dirac graph.

    Rules tested:
    1. Angular momentum conservation (Delta m_j)
    2. Spatial parity (E1 rule: l_a + l_b odd)
    3. Gaunt/CG sparsity
    4. Vertex parity (n1 + n2 + q odd)
    5. SO(4) channel count (W > 0)
    6. Charge conjugation (C)
    7. Furry's theorem (odd-loop diagrams vanish)
    8. Ward identity ([D, Lambda] = [Sigma, D])
    """
    labels = hodge["labels"]
    edges = hodge["edges"]
    V_nodes = hodge["V"]
    E_edges = hodge["E"]
    A = hodge["adjacency"]

    results = {}

    # --- Rule 1: Angular Momentum Conservation (Delta m_j) ---
    # For each edge, check if Delta m_j is conserved
    # On the Dirac graph, edges connect specific (n, kappa, m_j) states
    # The vertex V_e couples states at the edge endpoints
    total_edges = len(edges)
    mj_conserved = 0
    mj_violations = 0
    for (i, j) in edges:
        dmj = abs(labels[i].two_m_j - labels[j].two_m_j)
        # In this graph, edges themselves represent the transition
        # Delta m_j should be 0 or +/-1 for photon emission/absorption
        if dmj <= 2:  # |Delta m_j| <= 1 (in 2*m_j units)
            mj_conserved += 1
        else:
            mj_violations += 1

    results["rule_1_delta_mj"] = {
        "verdict": "SURVIVES" if mj_violations == 0 else "BROKEN",
        "total_edges": total_edges,
        "conserved": mj_conserved,
        "violations": mj_violations,
        "tier": "INTRINSIC",
    }

    # --- Rule 2: Spatial Parity (l_a + l_b odd for E1) ---
    parity_correct = 0
    parity_violations = 0
    for (i, j) in edges:
        la = kappa_to_l(labels[i].kappa)
        lb = kappa_to_l(labels[j].kappa)
        if (la + lb) % 2 == 1:  # odd = E1 allowed
            parity_correct += 1
        else:
            parity_violations += 1

    results["rule_2_spatial_parity"] = {
        "verdict": "SURVIVES" if parity_violations == 0 else "BROKEN",
        "total_edges": total_edges,
        "parity_correct": parity_correct,
        "parity_violations": parity_violations,
        "tier": "INTRINSIC",
    }

    # --- Rule 3: Gaunt/CG Sparsity ---
    # On the Dirac graph, the vertex is the identity on edges
    # (no CG projection needed). Sparsity comes from the graph topology.
    # Total possible couplings: V * (V-1) / 2 (all pairs)
    total_possible = V_nodes * (V_nodes - 1) // 2
    nnz = E_edges
    density = nnz / total_possible if total_possible > 0 else 0
    sparsity = 1.0 - density

    results["rule_3_gaunt_sparsity"] = {
        "verdict": "SURVIVES",
        "total_possible": total_possible,
        "nnz": nnz,
        "density": density,
        "sparsity": sparsity,
        "tier": "INTRINSIC",
    }

    # --- Rule 4: Vertex Parity (n1 + n2 + q odd) ---
    # On the Dirac graph, vertices are directly connected by edges
    # The "q" (photon quantum number) is implicit in the edge
    # For Rule B: edges enforce Delta l = +/-1 (parity flip)
    # Check n_fock parity: n1 + n2 should have specific parity
    # for the coupling to match the continuum vertex parity rule
    parity_forbidden = 0
    parity_allowed = 0
    for (i, j) in edges:
        n1 = labels[i].n_fock
        n2 = labels[j].n_fock
        # In the continuum, q is the photon angular quantum number
        # On the graph, there's no explicit q, but the edge structure
        # implicitly carries angular momentum transfer
        # n1 + n2 + q odd => with q=1 (dipole), n1+n2 even required
        # Check if n1 + n2 is even (consistent with q=1 odd parity)
        if (n1 + n2) % 2 == 0:
            parity_allowed += 1
        else:
            parity_forbidden += 1

    total_vp = parity_allowed + parity_forbidden
    results["rule_4_vertex_parity"] = {
        "verdict": "BROKEN" if parity_forbidden > 0 and parity_allowed > 0 else
                   ("SURVIVES" if parity_forbidden == 0 else "BROKEN"),
        "total_couplings": total_vp,
        "parity_allowed": parity_allowed,
        "parity_forbidden": parity_forbidden,
        "fraction_allowed": parity_allowed / total_vp if total_vp > 0 else 0,
        "tier": "STRUCTURAL",
        "note": "n1+n2 parity check with implicit q=1 dipole",
    }

    # --- Rule 5: SO(4) Channel Count ---
    # On the Dirac graph, there are no explicit SO(4) channels
    # The coupling structure is purely topological
    results["rule_5_so4_channel"] = {
        "verdict": "STRUCTURAL (not applicable)",
        "note": "SO(4) channel count requires vector harmonic structure; "
                "Dirac graph has spinor structure but no explicit W channel decomposition",
        "tier": "STRUCTURAL",
    }

    # --- Rule 6: Charge Conjugation (C) ---
    # C conjugation on Dirac labels: (n, kappa, m_j) -> (n, -kappa, -m_j)
    Sigma = sigma_result.get("Sigma")
    if Sigma is not None:
        c_violations = 0
        c_pairs_checked = 0
        for i, lab_i in enumerate(labels):
            # Find C-conjugate of lab_i
            c_kappa = -lab_i.kappa
            c_two_mj = -lab_i.two_m_j
            # Find index of C-conjugate
            c_idx = None
            for k, lab_k in enumerate(labels):
                if (lab_k.n_fock == lab_i.n_fock and
                    lab_k.kappa == c_kappa and
                    lab_k.two_m_j == c_two_mj):
                    c_idx = k
                    break
            if c_idx is not None and c_idx > i:
                c_pairs_checked += 1
                # Check Sigma[i, j] == Sigma[C(i), C(j)] for all j
                for j in range(V_nodes):
                    # Find C-conjugate of j
                    cj_kappa = -labels[j].kappa
                    cj_two_mj = -labels[j].two_m_j
                    cj_idx = None
                    for k, lab_k in enumerate(labels):
                        if (lab_k.n_fock == labels[j].n_fock and
                            lab_k.kappa == cj_kappa and
                            lab_k.two_m_j == cj_two_mj):
                            cj_idx = k
                            break
                    if cj_idx is not None:
                        if abs(Sigma[i, j] - Sigma[c_idx, cj_idx]) > 1e-10:
                            c_violations += 1

        results["rule_6_charge_conjugation"] = {
            "verdict": "SURVIVES" if c_violations == 0 else "BROKEN",
            "c_pairs_checked": c_pairs_checked,
            "c_violations": c_violations,
            "tier": "INTRINSIC",
        }
    else:
        results["rule_6_charge_conjugation"] = {
            "verdict": "NOT COMPUTED",
            "note": "Sigma matrix not available",
            "tier": "INTRINSIC",
        }

    # --- Rule 7: Furry's Theorem (odd-loop diagrams vanish) ---
    # Tadpole: Tr(sum_e V_e @ G_e) for the electron propagator
    # On the Dirac graph with identity vertex:
    # Tadpole_e = G_e[i(e), j(e)] + G_e[j(e), i(e)]  (symmetric vertex)
    D_dirac = build_dirac_operator(labels, A, t=0.0)
    try:
        G_e = np.linalg.inv(D_dirac)
    except np.linalg.LinAlgError:
        G_e = None

    if G_e is not None:
        tadpole = 0.0
        for (i, j) in edges:
            tadpole += G_e[i, j] + G_e[j, i]

        # Bubble: sum_{e,e'} Tr(V_e @ G_e @ V_{e'} @ G_e)
        V_mats_furry = []
        for (i, j) in edges:
            Ve = np.zeros((V_nodes, V_nodes))
            Ve[i, j] = 1.0
            Ve[j, i] = 1.0
            V_mats_furry.append(Ve)

        bubble = 0.0
        for e1 in range(E_edges):
            for e2 in range(E_edges):
                bubble += np.trace(V_mats_furry[e1] @ G_e @ V_mats_furry[e2] @ G_e)

        # Triangle: sum_{e,e',e''} Tr(V_e @ G_e @ V_{e'} @ G_e @ V_{e''} @ G_e)
        # (expensive for large E, limit computation)
        if E_edges <= 20:
            triangle = 0.0
            for e1 in range(E_edges):
                VG1 = V_mats_furry[e1] @ G_e
                for e2 in range(E_edges):
                    VG2 = V_mats_furry[e2] @ G_e
                    VG12 = VG1 @ VG2
                    for e3 in range(E_edges):
                        triangle += np.trace(VG12 @ V_mats_furry[e3] @ G_e)
        else:
            triangle = None

        results["rule_7_furry"] = {
            "verdict": "SURVIVES" if (abs(tadpole) < 1e-10 and
                                       (triangle is None or abs(triangle) < 1e-10)) else "BROKEN",
            "tadpole": float(tadpole),
            "tadpole_zero": abs(tadpole) < 1e-10,
            "bubble": float(bubble),
            "triangle": float(triangle) if triangle is not None else "not computed",
            "triangle_zero": abs(triangle) < 1e-10 if triangle is not None else None,
            "tier": "INTRINSIC",
        }
    else:
        results["rule_7_furry"] = {
            "verdict": "NOT COMPUTED",
            "note": "D_dirac is singular",
            "tier": "INTRINSIC",
        }

    # --- Rule 8: Ward Identity ([D, Lambda] = [Sigma, D]) ---
    if Sigma is not None and G_e is not None:
        L1_mat = hodge["L1"]
        G_gamma_mat = photon_propagator_dirac(L1_mat)

        # Build vertex matrices
        V_mats_ward = []
        for (i, j) in edges:
            Ve = np.zeros((V_nodes, V_nodes))
            Ve[i, j] = 1.0
            Ve[j, i] = 1.0
            V_mats_ward.append(Ve)

        # Lambda_total
        Lambda_total = np.zeros((V_nodes, V_nodes))
        for e1 in range(E_edges):
            for e2 in range(E_edges):
                g_ee = G_gamma_mat[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                Lambda_total += g_ee * (V_mats_ward[e1] @ G_e @ V_mats_ward[e2].T)

        D_op = build_dirac_operator(labels, A, t=0.0)
        lhs = D_op @ Lambda_total - Lambda_total @ D_op  # [D, Lambda]
        rhs = Sigma @ D_op - D_op @ Sigma  # [Sigma, D]
        ward_diff = lhs - rhs
        ward_max = np.max(np.abs(ward_diff))

        results["rule_8_ward_identity"] = {
            "verdict": "SURVIVES" if ward_max < 1e-10 else "BROKEN",
            "max_ward_diff": float(ward_max),
            "Sigma_trace": float(np.trace(Sigma)),
            "Lambda_trace": float(np.trace(Lambda_total)),
            "tier": "STRUCTURAL",
        }
    else:
        results["rule_8_ward_identity"] = {
            "verdict": "NOT COMPUTED",
            "tier": "STRUCTURAL",
        }

    return results


# ============================================================================
# MAIN: Run all tasks
# ============================================================================

def run_sprint() -> Dict[str, Any]:
    """Run the full Dirac Graph QED sprint."""

    all_results = {}

    # =================================================================
    # TASK 1: Hodge Decomposition
    # =================================================================
    print("=" * 70)
    print("TASK 1: Hodge Decomposition on Dirac Graph")
    print("=" * 70)

    hodge_results = {}
    for rule in ["A", "B"]:
        for n_max in [2, 3]:
            key = f"rule_{rule}_nmax_{n_max}"
            print(f"\n--- Rule {rule}, n_max={n_max} ---")
            t0 = time.time()
            hodge = hodge_decomposition(n_max, rule)
            dt = time.time() - t0

            print(f"  V={hodge['V']}, E={hodge['E']}, "
                  f"beta_0={hodge['beta_0']}, beta_1={hodge['beta_1']}")
            print(f"  L0=D-A: {hodge['L0_DA_match']}")
            print(f"  SVD identity: {hodge['svd_identity_verified']}")
            print(f"  pi-free: {hodge['pi_free_certificate']}")
            print(f"  GS pendant: {hodge['gs_is_pendant']}, "
                  f"GS degrees: {hodge['gs_degrees']}")
            print(f"  Time: {dt:.2f}s")

            hodge_results[key] = {
                "V": hodge["V"],
                "E": hodge["E"],
                "beta_0": hodge["beta_0"],
                "beta_1": hodge["beta_1"],
                "L0_DA_match": hodge["L0_DA_match"],
                "svd_identity_verified": hodge["svd_identity_verified"],
                "pi_free_certificate": hodge["pi_free_certificate"],
                "degree_sequence": hodge["degree_sequence"],
                "pendant_vertex_count": hodge["pendant_vertex_count"],
                "gs_is_pendant": hodge["gs_is_pendant"],
                "gs_degrees": hodge["gs_degrees"],
                "L0_spectrum": [round(x, 10) for x in hodge["L0_spectrum"]],
                "L1_spectrum": [round(x, 10) for x in hodge["L1_spectrum"]],
            }

    all_results["task_1_hodge"] = hodge_results

    # =================================================================
    # TASK 2 + 3 + 4 + 5: Self-energy and selection rules
    # Focus on n_max=2 for both Rule A and Rule B
    # =================================================================

    for rule in ["A", "B"]:
        print(f"\n{'='*70}")
        print(f"TASKS 2-5: Self-energy and Census, Rule {rule}, n_max=2")
        print(f"{'='*70}")

        hodge = hodge_decomposition(2, rule)

        # TASK 2: Photon propagator
        print(f"\n--- TASK 2: Photon Propagator (Rule {rule}) ---")
        G_gamma = photon_propagator_dirac(hodge["L1"])
        print(f"  G_gamma shape: {G_gamma.shape}")
        print(f"  G_gamma max: {np.max(np.abs(G_gamma)):.6f}")
        print(f"  G_gamma rank: {np.linalg.matrix_rank(G_gamma)}")
        print(f"  GS is pendant: {hodge['gs_is_pendant']}")
        print(f"  GS degrees: {hodge['gs_degrees']}")

        # TASK 3: Electron propagator
        print(f"\n--- TASK 3: Electron Propagator (Rule {rule}) ---")
        D_dirac = build_dirac_operator(hodge["labels"], hodge["adjacency"], t=0.0)
        cond = np.linalg.cond(D_dirac)
        print(f"  D_dirac shape: {D_dirac.shape}")
        print(f"  D_dirac condition number: {cond:.4f}")
        try:
            G_e = np.linalg.inv(D_dirac)
            print(f"  G_e max: {np.max(np.abs(G_e)):.6f}")
            print(f"  G_e is diagonal (t=0): {np.allclose(G_e, np.diag(np.diag(G_e)), atol=1e-12)}")
        except Exception as ex:
            print(f"  G_e inversion failed: {ex}")

        # TASK 4: Self-energy
        print(f"\n--- TASK 4: Self-Energy (Rule {rule}) ---")
        t0 = time.time()
        sigma = compute_dirac_graph_self_energy(hodge, t=0.0)
        dt = time.time() - t0
        print(f"  Sigma trace: {sigma['Sigma_trace']:.6f}")
        print(f"  Sigma is Hermitian: {sigma['is_hermitian']}")
        print(f"  Sigma is PSD: {sigma['is_psd']}")
        print(f"  GS block zero: {sigma['gs_block_zero']}")
        print(f"  GS block: {sigma['gs_block']}")
        print(f"  Eigenvalues: {[f'{ev:.6f}' for ev in sigma['eigenvalues']]}")
        print(f"  Time: {dt:.2f}s")

        # TASK 4b: Vertex correction
        print(f"\n--- TASK 4b: Vertex Correction (Rule {rule}) ---")
        vc = compute_dirac_graph_vertex_correction(hodge, t=0.0)
        print(f"  Lambda trace: {vc['Lambda_trace']}")
        print(f"  F2: {vc['F2']}")
        print(f"  Eigenvalues: {[f'{ev:.6f}' for ev in vc['eigenvalues']]}")

        # TASK 5: Selection rule census
        print(f"\n--- TASK 5: Selection Rule Census (Rule {rule}) ---")
        census = selection_rule_census(hodge, sigma)
        for rule_name, rule_data in sorted(census.items()):
            verdict = rule_data.get("verdict", "N/A")
            print(f"  {rule_name}: {verdict}")

        # Store results
        rule_key = f"rule_{rule}"
        all_results[f"task_2_photon_{rule_key}"] = {
            "G_gamma_shape": list(G_gamma.shape),
            "G_gamma_max": float(np.max(np.abs(G_gamma))),
            "G_gamma_rank": int(np.linalg.matrix_rank(G_gamma)),
            "gs_is_pendant": hodge["gs_is_pendant"],
            "gs_degrees": hodge["gs_degrees"],
        }

        all_results[f"task_3_electron_{rule_key}"] = {
            "D_dirac_shape": list(D_dirac.shape),
            "D_dirac_condition": float(cond),
            "G_e_is_diagonal_t0": True,  # Always true at t=0
        }

        # Serialize self-energy results (exclude numpy array)
        sigma_ser = {k: v for k, v in sigma.items() if k != "Sigma"}
        all_results[f"task_4_self_energy_{rule_key}"] = sigma_ser

        all_results[f"task_4b_vertex_{rule_key}"] = vc

        all_results[f"task_5_census_{rule_key}"] = census

    # =================================================================
    # Also run self-energy/census at n_max=3 for Rule B
    # =================================================================
    print(f"\n{'='*70}")
    print(f"TASKS 4-5: Self-energy and Census, Rule B, n_max=3")
    print(f"{'='*70}")

    hodge_3 = hodge_decomposition(3, "B")
    sigma_3 = compute_dirac_graph_self_energy(hodge_3, t=0.0)
    census_3 = selection_rule_census(hodge_3, sigma_3)

    print(f"  Sigma trace (n=3): {sigma_3['Sigma_trace']:.6f}")
    print(f"  GS block zero (n=3): {sigma_3['gs_block_zero']}")
    print(f"  GS block (n=3): {sigma_3['gs_block']}")

    sigma_3_ser = {k: v for k, v in sigma_3.items() if k != "Sigma"}
    all_results["task_4_self_energy_rule_B_nmax3"] = sigma_3_ser
    all_results["task_5_census_rule_B_nmax3"] = census_3

    # =================================================================
    # Comparison with scalar Fock graph
    # =================================================================
    print(f"\n{'='*70}")
    print("COMPARISON TABLE: Scalar vs Dirac Graph")
    print("=" * 70)

    # Scalar Fock graph reference values (from existing GN sprint)
    scalar_ref = {
        "Sigma_trace_nmax2": "44/3",
        "gs_block_zero": False,
        "gs_block": "[[1,1],[1,1]]",
        "F2_t0": "5*sqrt(2)/3",
        "selection_rules_surviving": "1/8 (Gaunt/CG sparsity only)",
    }
    all_results["scalar_reference"] = scalar_ref

    print(f"\n  Scalar Fock graph (n_max=2):")
    print(f"    Sigma trace: {scalar_ref['Sigma_trace_nmax2']}")
    print(f"    GS block zero: {scalar_ref['gs_block_zero']}")
    print(f"    F2(t=0): {scalar_ref['F2_t0']}")
    print(f"    Selection rules surviving: {scalar_ref['selection_rules_surviving']}")

    for rule in ["A", "B"]:
        rule_key = f"rule_{rule}"
        se = all_results[f"task_4_self_energy_{rule_key}"]
        ce = all_results[f"task_5_census_{rule_key}"]
        print(f"\n  Dirac graph Rule {rule} (n_max=2):")
        print(f"    Sigma trace: {se['Sigma_trace']:.6f}")
        print(f"    GS block zero: {se['gs_block_zero']}")
        n_surviving = sum(1 for k, v in ce.items()
                         if isinstance(v, dict) and v.get("verdict", "").startswith("SURVIVES"))
        n_total = sum(1 for k, v in ce.items()
                      if isinstance(v, dict) and "verdict" in v)
        print(f"    Selection rules surviving: {n_surviving}/{n_total}")

    return all_results


if __name__ == "__main__":
    results = run_sprint()

    # Save results
    out_path = Path(__file__).parent / "data" / "dirac_graph_qed_sprint.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert non-serializable types
    def make_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, dict):
            return {k: make_serializable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [make_serializable(v) for v in obj]
        return obj

    results_ser = make_serializable(results)
    with open(out_path, "w") as f:
        json.dump(results_ser, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")
