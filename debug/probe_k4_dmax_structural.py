"""
Probe: Is κ = -1/d_max² structural across sphere dimensions?

S³ (Coulomb): κ = -1/16 = -1/4², d_max = dim(R⁴) = 4
S⁵ (HO):     predicted κ_S5 = -1/36 = -1/6², d_max = dim(R⁶) = 6

This script tests the hypothesis by:
1. Building the S³ Fock graph and verifying κ = -1/16 maps eigenvalues to energies
2. Building the S⁵ Bargmann-Segal graph and testing κ_S5 = -1/36
3. Testing S¹ and S² analogs
4. Comparing alternative explanations: -1/(2*N_init)², -1/(2*l_max+2)², etc.

Output: debug/data/probe_k4_dmax_structural.json
"""

from __future__ import annotations

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import json
import numpy as np
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple


# ============================================================================
# S³ Coulomb graph: verify κ = -1/16
# ============================================================================

def build_s3_fock_graph_laplacian(n_max: int) -> Tuple[np.ndarray, List[Tuple[int, int, int]]]:
    """
    Build the S³ Fock graph (binary adjacency, unit weights).

    Nodes: (n, l, m) with n = 1..n_max, l = 0..n-1, m = -l..l
    Edges:
      - Angular: (n,l,m) <-> (n,l,m±1)  [L± transitions]
      - Radial:  (n,l,m) <-> (n±1,l,m)  [T± transitions, same l,m]

    Returns L = D - A (unnormalized graph Laplacian) and the node list.
    """
    nodes = []
    idx_map = {}
    i = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                nodes.append((n, l, m))
                idx_map[(n, l, m)] = i
                i += 1

    V = len(nodes)
    A = np.zeros((V, V))

    for n, l, m in nodes:
        i = idx_map[(n, l, m)]
        # Angular: m -> m+1
        if m < l and (n, l, m+1) in idx_map:
            j = idx_map[(n, l, m+1)]
            A[i, j] = 1.0
            A[j, i] = 1.0
        # Radial: n -> n+1 (same l, m)
        if n < n_max and (n+1, l, m) in idx_map:
            j = idx_map[(n+1, l, m)]
            A[i, j] = 1.0
            A[j, i] = 1.0

    D = np.diag(A.sum(axis=1))
    L = D - A
    return L, nodes


def test_s3_kappa(n_max: int = 4) -> Dict:
    """
    Verify that κ = -1/16 applied to the S³ graph Laplacian
    reproduces the hydrogen spectrum E_n = -1/(2n²).

    The S³ Laplace-Beltrami eigenvalues are λ_n = n²-1.
    The graph Laplacian L = D - A, when we compute κ*L, should give
    eigenvalues that map to -1/(2n²) for each n-shell.

    Actually, the graph Laplacian eigenvalues within each (l,m) sector
    are a tridiagonal problem. The key relationship is:
    κ * (graph eigenvalue) ≈ -1/(2n²)

    More precisely: the S³ conformal equivalence says
    κ * (n² - 1) = -1/(2n²) is NOT exact; instead
    κ is determined by requiring the RIGHT spectrum in the combined
    H = κ*(D-A) + V formulation.
    """
    L, nodes = build_s3_fock_graph_laplacian(n_max)
    evals = np.sort(np.linalg.eigvalsh(L))

    # Exact hydrogen energies with degeneracy
    exact_energies = []
    for n in range(1, n_max + 1):
        for _ in range(n * n):  # n² degeneracy
            exact_energies.append(-0.5 / (n * n))
    exact_energies = np.sort(exact_energies)

    # The graph Laplacian eigenvalues by (l,m) sector
    # Within each (l,m) sector, the nodes are (l+1, l, m), (l+2, l, m), ...
    # connected by T± radial transitions — a path graph of length n_max - l.
    # The degree matrix D has entries 2 (interior) or 1 (boundary).
    # So L restricted to each (l,m) sector is tridiagonal.

    # The full L eigenvalues:
    # For a path graph of length p, eigenvalues are 2 - 2*cos(k*π/(p+1)) for k=1..p
    # but the GeoVac graph also has angular edges within each n-shell.

    # Let's just compute κ*L eigenvalues and compare to exact
    kappa = -1.0/16.0
    H_graph = kappa * L
    evals_graph = np.sort(np.linalg.eigvalsh(H_graph))

    # Compare
    results = {
        "n_max": n_max,
        "V": len(nodes),
        "kappa": -1/16,
        "d_max": 4,
        "kappa_equals_neg_inv_dmax_sq": True,
        "graph_laplacian_evals": evals.tolist(),
        "kappa_L_evals": evals_graph.tolist(),
        "exact_energies": exact_energies,
    }

    # Check if κ*L eigenvalues are close to -1/(2n²)
    # (They won't be exact because the graph Laplacian alone doesn't reproduce
    # the full hydrogen Hamiltonian — you need the diagonal V term too.
    # κ maps the KINETIC part.)

    return results


# ============================================================================
# S⁵ Bargmann-Segal graph
# ============================================================================

def test_s5_kappa(N_max: int = 4) -> Dict:
    """
    Test whether κ_S5 = -1/36 = -1/d_max² maps the S⁵ graph Laplacian
    to the HO spectrum.

    Key difference from S³: the S⁵ graph has NON-UNIFORM edge weights
    (squared dipole matrix elements), while the S³ graph has binary weights.

    The S⁵ Laplace-Beltrami eigenvalues are N(N+4).
    The HO energies are ℏω(N + 3/2).

    For S³: eigenvalues n²-1, energies -1/(2n²).
    The relationship is: κ*(n²-1) ≈ ... but actually κ is the scale
    that maps the graph Laplacian to the Fock momentum-space kinetic energy.

    For S⁵: what κ would map graph Laplacian eigenvalues to HO energies?
    """
    from geovac.nuclear.bargmann_graph import build_bargmann_graph

    g = build_bargmann_graph(N_max)

    # Get the adjacency matrix (with non-uniform rational weights)
    A = g.adjacency_dense()

    # Graph Laplacian L = D - A
    L = g.graph_laplacian_dense()

    # Eigenvalues of L
    evals_L = np.sort(np.linalg.eigvalsh(L))

    # The diagonal (HO Hamiltonian)
    H_diag = g.hamiltonian_dense(hw=1.0)
    diag_vals = np.sort(np.diag(H_diag))

    # S⁵ Laplace-Beltrami eigenvalues: N(N+4) with degeneracy (N+1)(N+2)/2
    s5_LB_evals = []
    for N in range(N_max + 1):
        deg = (N + 1) * (N + 2) // 2
        for _ in range(deg):
            s5_LB_evals.append(N * (N + 4))
    s5_LB_evals = np.sort(s5_LB_evals).astype(float)

    # ---- Test 1: Does κ_S5 = -1/36 map L eigenvalues to anything meaningful?
    kappa_pred = -1.0 / 36.0
    kappa_L_pred = kappa_pred * evals_L

    # ---- Test 2: What κ_S5 would be needed to map L evals to HO energies?
    # If κ*λ_L = ℏω(N+3/2), then κ = (N+3/2)/λ_L
    # But L has a zero eigenvalue, so this is singular.
    # Instead, look at non-zero eigenvalues.

    # ---- Test 3: Does L have eigenvalues proportional to N(N+4)?
    # I.e., does L = c * (Laplace-Beltrami on S⁵)?

    # ---- Test 4: What if we look at D + A instead of D - A?
    # Or -A (just the adjacency)?
    evals_A = np.sort(np.linalg.eigvalsh(A))

    # ---- Test 5: Ratio analysis
    # For S³: the per-shell average of L eigenvalues should be ~n²-1
    # For S⁵: the per-shell average of L eigenvalues should be ~N(N+4)
    shell_L_avgs = []
    idx = 0
    for N in range(N_max + 1):
        deg = (N + 1) * (N + 2) // 2
        # The nodes are grouped by N-shell. But L mixes across shells
        # (edges connect N <-> N+1), so the eigenvalues don't decompose
        # by shell. We need to look at the full spectrum.
        idx += deg

    # ---- Test 6: For S³, the key insight from Paper 7 is that the
    # graph Laplacian on S³ with binary weights has eigenvalues that
    # BLOCK-DECOMPOSE by (l,m) sector (because T± only connects same l,m).
    # Within each (l,m) sector of size (n_max - l), the tridiagonal L
    # has known eigenvalues.
    #
    # For S⁵, the graph also has a block structure: edges only connect
    # N to N±1 (bipartite-like). Let's check the block structure.

    # Degree distribution
    degree = A.sum(axis=1)
    nodes = g.nodes
    degree_by_N = {}
    for i, (N, l, m) in enumerate(nodes):
        if N not in degree_by_N:
            degree_by_N[N] = []
        degree_by_N[N].append(degree[i])

    # ---- Test 7: Can we find ANY κ that maps L eigenvalues to
    # a linear function of N (i.e., the HO spectrum)?
    # If L eigenvalues are {λ_i}, and we want κ*λ_i = c*N_i + d,
    # that requires L eigenvalues to be linear in N — but L eigenvalues
    # are NOT grouped by N (L is not block-diagonal in N).

    # However, the S³ graph Laplacian is NOT block-diagonal in n either.
    # The key is that within each (l,m) sector, it IS block-diagonal,
    # and the eigenvalues of the (l,m)-sector tridiagonal match n²-1.

    # ---- Test 8: For S⁵, look at the (l,m_l)-restricted subgraph
    # Edges connect (N,l,m) to (N+1,l±1,m') — so they CHANGE l.
    # Unlike S³ where T± preserves l, the S⁵ dipole transitions change l.
    # This means L does NOT decompose into (l,m) sectors.

    # ---- Key structural difference:
    # S³: edges preserve (l,m), change n. Laplacian block-diags by (l,m).
    # S⁵: edges change (N,l,m) simultaneously. No simple block decomposition.

    # ---- Test 9: Fit κ by least-squares on the eigenvalue ratios.
    # For non-zero eigenvalues, find κ such that κ*λ_L best approximates
    # the S⁵ LB eigenvalues N(N+4) or HO energies N+3/2.

    # Remove zero eigenvalues for ratio analysis
    nonzero_mask_L = np.abs(evals_L) > 1e-10
    evals_L_nz = evals_L[nonzero_mask_L]

    # Test: ratio of consecutive non-zero eigenvalues
    if len(evals_L_nz) > 1:
        ratios = evals_L_nz[1:] / evals_L_nz[:-1]
    else:
        ratios = np.array([])

    # ---- Test 10: Check if L eigenvalues match S⁵ LB eigenvalues
    # (up to a constant scale)
    # S⁵ LB: 0, 5, 12, 21, 32, 45 for N=0,1,2,3,4,5
    # with degeneracies 1, 3, 6, 10, 15, 21

    # Sort both and try to match
    n_zero_L = np.sum(np.abs(evals_L) < 1e-10)
    n_zero_LB = np.sum(np.abs(s5_LB_evals) < 1e-10)

    # Best-fit scale: κ_fit = Σ(LB * L) / Σ(L²) for non-zero pairs
    if len(evals_L_nz) > 0 and len(s5_LB_evals) > 0:
        # Both should have 1 zero eigenvalue. Align non-zero parts.
        s5_nz = s5_LB_evals[s5_LB_evals > 1e-10]
        # sizes may differ; take the shorter
        n_match = min(len(evals_L_nz), len(s5_nz))
        if n_match > 0:
            kappa_fit_LB = np.dot(s5_nz[:n_match], evals_L_nz[:n_match]) / np.dot(evals_L_nz[:n_match], evals_L_nz[:n_match])
            residual_LB = np.sqrt(np.mean((s5_nz[:n_match] - kappa_fit_LB * evals_L_nz[:n_match])**2))
        else:
            kappa_fit_LB = None
            residual_LB = None
    else:
        kappa_fit_LB = None
        residual_LB = None

    # Best-fit scale: κ_fit = Σ(HO * L) / Σ(L²) for HO energies
    ho_evals = diag_vals  # already sorted
    if len(evals_L_nz) > 0:
        ho_nz = ho_evals[ho_evals > ho_evals[0] + 1e-10]  # skip GS degeneracy
        # Actually HO energies are N+3/2, all positive. Just use all.
        n_match_ho = min(len(evals_L), len(ho_evals))
        kappa_fit_HO = np.dot(ho_evals[:n_match_ho], evals_L[:n_match_ho]) / np.dot(evals_L[:n_match_ho], evals_L[:n_match_ho])
        residual_HO = np.sqrt(np.mean((ho_evals[:n_match_ho] - kappa_fit_HO * evals_L[:n_match_ho])**2))
    else:
        kappa_fit_HO = None
        residual_HO = None

    results = {
        "N_max": N_max,
        "V": g.n_nodes,
        "E": len(g.adjacency),
        "d_max_S5": 6,
        "kappa_predicted": -1/36,
        "graph_laplacian_evals": sorted(evals_L.tolist()),
        "n_zero_eigenvalues_L": int(n_zero_L),
        "adjacency_evals": sorted(evals_A.tolist()),
        "s5_LB_evals": sorted(s5_LB_evals.tolist()),
        "ho_energies": sorted(ho_evals.tolist()),
        "degree_by_N": {str(k): [float(x) for x in v] for k, v in degree_by_N.items()},
        "kappa_fit_to_LB": float(kappa_fit_LB) if kappa_fit_LB is not None else None,
        "residual_LB": float(residual_LB) if residual_LB is not None else None,
        "kappa_fit_to_HO": float(kappa_fit_HO) if kappa_fit_HO is not None else None,
        "residual_HO": float(residual_HO) if residual_HO is not None else None,
        "edge_weights_uniform": False,  # S⁵ has non-uniform weights
    }

    return results


# ============================================================================
# S¹ path/cycle graph
# ============================================================================

def test_s1_kappa() -> Dict:
    """
    S¹ = circle, ambient dimension 2.
    Predicted: κ_S1 = -1/4 = -1/2².

    The "graph" for S¹: a cycle graph C_p or path graph P_p.
    The Laplacian eigenvalues of a cycle C_p are 2 - 2*cos(2πk/p) for k=0..p-1.
    The LB eigenvalues of S¹ are m² for m ∈ Z.

    For a path graph P_p (nodes 0..p-1):
    Eigenvalues of L are 2 - 2*cos(πk/p) for k=0..p-1.

    On S¹: LB eigenvalue = m², spectrum on unit circle is E_m = m² for m=0,±1,±2,...
    The "free particle on a circle" has energies m²ℏ²/(2mR²).

    There's no natural "Fock projection" for S¹ — S¹ is too simple.
    But we can check: does -1/4 times path graph eigenvalues give m²?
    """
    results_s1 = {}

    for p in [3, 5, 7, 10]:
        # Path graph P_p
        A = np.zeros((p, p))
        for i in range(p - 1):
            A[i, i+1] = 1.0
            A[i+1, i] = 1.0
        D = np.diag(A.sum(axis=1))
        L = D - A
        evals = np.sort(np.linalg.eigvalsh(L))

        # Cycle graph C_p
        A_cyc = np.zeros((p, p))
        for i in range(p):
            A_cyc[i, (i+1) % p] = 1.0
            A_cyc[(i+1) % p, i] = 1.0
        D_cyc = np.diag(A_cyc.sum(axis=1))
        L_cyc = D_cyc - A_cyc
        evals_cyc = np.sort(np.linalg.eigvalsh(L_cyc))

        # S¹ LB eigenvalues (first p): 0, 1, 1, 4, 4, 9, 9, ...
        s1_lb = []
        for m in range(p):
            if m == 0:
                s1_lb.append(0)
            else:
                s1_lb.append(m * m)
                s1_lb.append(m * m)
        s1_lb = sorted(s1_lb[:p])

        # Test kappa = -1/4
        kappa_s1 = -1.0 / 4.0
        kappa_L_path = kappa_s1 * evals
        kappa_L_cycle = kappa_s1 * evals_cyc

        results_s1[f"p={p}"] = {
            "path_L_evals": evals.tolist(),
            "cycle_L_evals": evals_cyc.tolist(),
            "kappa_path_evals": kappa_L_path.tolist(),
            "kappa_cycle_evals": kappa_L_cycle.tolist(),
            "s1_LB_first_p": s1_lb,
        }

    return {"s1_tests": results_s1, "d_max_S1": 2, "kappa_predicted": -1/4}


# ============================================================================
# S² sphere graph
# ============================================================================

def test_s2_kappa() -> Dict:
    """
    S² lives in R³, d_max = 3. Predicted κ_S2 = -1/9 = -1/3².

    S² LB eigenvalues: l(l+1) with degeneracy 2l+1.

    There's no canonical "Fock graph" for S². But we can construct
    the natural graph on spherical harmonic labels (l, m):
    - Nodes: (l, m) for l = 0..l_max, m = -l..l
    - Edges: (l,m) <-> (l,m±1) [angular, within same l]

    This is the ANGULAR part of the S³ graph (no radial transitions).
    Each l-shell is a path graph of length 2l+1.
    """
    results_s2 = {}

    for l_max in [2, 3, 4]:
        nodes = []
        idx_map = {}
        i = 0
        for l in range(l_max + 1):
            for m in range(-l, l + 1):
                nodes.append((l, m))
                idx_map[(l, m)] = i
                i += 1

        V = len(nodes)
        A = np.zeros((V, V))

        for l, m in nodes:
            i_node = idx_map[(l, m)]
            # m -> m+1 within same l
            if m < l and (l, m+1) in idx_map:
                j_node = idx_map[(l, m+1)]
                A[i_node, j_node] = 1.0
                A[j_node, i_node] = 1.0

        D = np.diag(A.sum(axis=1))
        L = D - A
        evals = np.sort(np.linalg.eigvalsh(L))

        # S² LB eigenvalues: l(l+1) with degeneracy 2l+1
        s2_lb = []
        for l in range(l_max + 1):
            for _ in range(2*l + 1):
                s2_lb.append(l * (l + 1))
        s2_lb = sorted(s2_lb)

        kappa_s2 = -1.0 / 9.0
        kappa_L = kappa_s2 * evals

        results_s2[f"l_max={l_max}"] = {
            "V": V,
            "L_evals": evals.tolist(),
            "kappa_L_evals": kappa_L.tolist(),
            "s2_LB_evals": s2_lb,
        }

    return {"s2_tests": results_s2, "d_max_S2": 3, "kappa_predicted": -1/9}


# ============================================================================
# Main analysis: structural test of κ = -1/d_max²
# ============================================================================

def analyze_s3_graph_eigenvalues(n_max: int = 4) -> Dict:
    """
    Detailed analysis of S³ graph Laplacian eigenvalues vs n²-1.

    The S³ Fock graph with binary weights, restricted to a single (l,m) sector,
    is a path graph of length (n_max - l). The tridiagonal Laplacian of
    a path graph P_p has eigenvalues 2 - 2*cos(πk/(p+1)) for k=1..p.

    Paper 7 proves that the FULL graph Laplacian (D-A) on the S³ Fock graph
    has eigenvalues that, within each (l,m) sector, are exactly
    the eigenvalues of the path graph. The n²-1 connection comes through
    the conformal equivalence, not through a direct eigenvalue match.

    The actual relationship is:
    H = κ*(D-A) gives the kinetic energy operator
    The full Hamiltonian is H_full = κ*(D-A) + V where V is the diagonal
    And the eigenvalues of H_full are -Z²/(2n²).

    So κ is the scale factor for the KINETIC part only.
    The question is: what determines κ = -1/16?
    """
    L, nodes = build_s3_fock_graph_laplacian(n_max)
    V_count = len(nodes)

    # Block-decompose by (l,m)
    sectors = {}
    for i, (n, l, m) in enumerate(nodes):
        key = (l, m)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(i)

    sector_evals = {}
    for (l, m), indices in sectors.items():
        if len(indices) == 0:
            continue
        L_block = L[np.ix_(indices, indices)]
        ev = np.sort(np.linalg.eigvalsh(L_block))
        sector_evals[f"({l},{m})"] = {
            "size": len(indices),
            "eigenvalues": ev.tolist(),
            # Path graph of length p has evals 2 - 2*cos(πk/(p+1)) for k=1..p
            "path_graph_evals": [2 - 2*np.cos(np.pi*k/(len(indices)+1))
                                 for k in range(1, len(indices)+1)],
        }

    # Check: do sector eigenvalues match path graph eigenvalues?
    match_path = True
    max_diff = 0.0
    for key, data in sector_evals.items():
        diff = max(abs(a - b) for a, b in zip(data["eigenvalues"], data["path_graph_evals"]))
        max_diff = max(max_diff, diff)
        if diff > 1e-10:
            match_path = False

    return {
        "n_max": n_max,
        "V": V_count,
        "sectors_match_path_graph": match_path,
        "max_sector_path_diff": max_diff,
        "sector_eigenvalues": sector_evals,
    }


def analyze_s5_graph_eigenvalues(N_max: int = 3) -> Dict:
    """
    Detailed analysis of S⁵ Bargmann-Segal graph Laplacian eigenvalues.

    Key structural difference from S³:
    - S³ graph has BINARY edge weights and decomposes by (l,m) sector
    - S⁵ graph has NON-UNIFORM edge weights (squared dipole matrix elements)
      and does NOT decompose by (l,m) sector (edges change l)

    This means the S⁵ graph Laplacian is a fundamentally different object.
    """
    from geovac.nuclear.bargmann_graph import build_bargmann_graph

    g = build_bargmann_graph(N_max)
    A = g.adjacency_dense()
    L = g.graph_laplacian_dense()

    evals_L = np.sort(np.linalg.eigvalsh(L))

    # S⁵ LB eigenvalues for comparison
    s5_lb = []
    for N in range(N_max + 1):
        deg = (N + 1) * (N + 2) // 2
        for _ in range(deg):
            s5_lb.append(float(N * (N + 4)))
    s5_lb = np.sort(s5_lb)

    # Check: do L eigenvalues have any proportionality to N(N+4)?
    # Count zero eigenvalues
    n_zero = int(np.sum(np.abs(evals_L) < 1e-10))
    n_zero_lb = int(np.sum(np.abs(s5_lb) < 1e-10))

    # Non-zero eigenvalues
    evals_nz = evals_L[np.abs(evals_L) > 1e-10]
    s5_lb_nz = s5_lb[s5_lb > 1e-10]

    # Try least-squares fit: L_eval ≈ c * N(N+4) + d
    # If L eigenvalues are proportional to S⁵ LB, then c should be constant
    # and d should be 0.

    # Since eigenvalues don't decompose by shell, let's try a different approach:
    # Check if L eigenvalues, sorted, are proportional to sorted S⁵ LB eigenvalues.
    n_match = min(len(evals_nz), len(s5_lb_nz))
    if n_match > 0:
        # Least-squares fit: evals_nz ≈ c * s5_lb_nz
        c_fit = np.dot(evals_nz[:n_match], s5_lb_nz[:n_match]) / np.dot(s5_lb_nz[:n_match], s5_lb_nz[:n_match])
        residual = np.sqrt(np.mean((evals_nz[:n_match] - c_fit * s5_lb_nz[:n_match])**2))
        max_rel_err = np.max(np.abs(evals_nz[:n_match] - c_fit * s5_lb_nz[:n_match]) / (np.abs(evals_nz[:n_match]) + 1e-15))
    else:
        c_fit = None
        residual = None
        max_rel_err = None

    # Also check: adjacency eigenvalues
    evals_A = np.sort(np.linalg.eigvalsh(A))

    # Key question: is there ANY scale κ such that κ*L has eigenvalues
    # matching N(N+4) or (N+3/2)?
    # If c_fit exists and residual is small, then κ = 1/c_fit would work.

    # Also compute per-shell degree sums
    nodes = g.nodes
    degree = A.sum(axis=1)
    shell_avg_degree = {}
    for i, (N, l, m) in enumerate(nodes):
        if N not in shell_avg_degree:
            shell_avg_degree[N] = []
        shell_avg_degree[N].append(degree[i])

    avg_degrees = {N: np.mean(degs) for N, degs in shell_avg_degree.items()}

    return {
        "N_max": N_max,
        "V": g.n_nodes,
        "E": len(g.adjacency),
        "n_zero_L": n_zero,
        "n_zero_LB": n_zero_lb,
        "L_evals_sorted": evals_L.tolist(),
        "s5_LB_evals_sorted": s5_lb.tolist(),
        "A_evals_sorted": evals_A.tolist(),
        "scale_fit_L_to_LB": float(c_fit) if c_fit is not None else None,
        "residual_L_to_LB": float(residual) if residual is not None else None,
        "max_rel_error_L_to_LB": float(max_rel_err) if max_rel_err is not None else None,
        "avg_degree_by_shell": {str(k): float(v) for k, v in avg_degrees.items()},
        "edge_weights_non_uniform": True,
    }


# ============================================================================
# The definitive test: κ from the S³ conformal map
# ============================================================================

def definitive_kappa_test() -> Dict:
    """
    The definitive analysis of what κ = -1/16 actually IS.

    From Paper 7: The S³ graph Laplacian eigenvalues in each (l,m) sector
    are path-graph eigenvalues. The conformal map from S³ to R³ via
    stereographic projection introduces the scale factor.

    The key equation from Paper 0 / Paper 7:
    κ = -1/d_max² where d_max is the "packing dimension" from the
    packing construction.

    For S³ in R⁴: d_max = 4, κ = -1/16
    For S⁵ in R⁶: d_max = 6, predicted κ_S5 = -1/36

    But there are alternative explanations:
    - κ = -1/(2*N_init)² where N_init = 2 (Paper 0 initial count)
    - κ = -(d₀/2)⁴ = -(1/2)⁴ = -1/16 (packing diameter d₀ = 1)
    - κ = -1/(2*l_max+2)² with l_max = 1

    To distinguish: we need to check the S⁵ case.
    If κ_S5 = -1/36 = -1/6²: confirms d_max = dim(ambient R^d).
    If κ_S5 = -1/16: would suggest κ is universal (not d_max-dependent).
    If κ_S5 = -1/64 = -1/8²: would suggest 2*d = 2*dim(S^d).
    """

    # ---- S³ analysis ----
    # The Fock projection gives: E_n = -Z²/(2n²)
    # The S³ LB eigenvalue for the n-th shell is: λ_n = n² - 1
    # The graph Laplacian diagonal (from D matrix) for a vertex in shell n
    # has degree contributions.
    #
    # The key: κ*(n²-1) does NOT directly equal -1/(2n²).
    # Instead, κ is the coefficient of the KINETIC operator κ*(D-A).
    # The full Hamiltonian is κ*(D-A) + V where V_ii = -Z/n²_i.
    # The eigenvalues of this full H are -Z²/(2n²).
    #
    # What PAPER 0 derives: κ = -1/(d_max)² comes from the packing
    # construction. d_max = 4 for 3D (the ambient dimension of S³ ⊂ R⁴).
    # This is proven from geometric packing, not fitted.

    # ---- S⁵ analysis ----
    # For the 3D HO on S⁵:
    # The S⁵ LB eigenvalue for the N-th shell is: λ_N = N(N+4)
    # The HO energy is: E_N = ℏω(N + 3/2)
    # The S⁵ packing dimension would be d_max = 6 (ambient dimension of S⁵ ⊂ R⁶).
    #
    # But there's a KEY structural difference:
    # For S³, the spectrum n²-1 is QUADRATIC in n.
    # The Coulomb energy -1/(2n²) is INVERSE QUADRATIC in n.
    # These are related by the Fock conformal map p₀ = √(-2E) = Z/n.
    #
    # For S⁵, the spectrum N(N+4) is QUADRATIC in N.
    # The HO energy N+3/2 is LINEAR in N.
    # Paper 24 says the HO Hamiltonian is DIAGONAL (not from D-A).
    # The graph adjacency encodes dipole transitions, NOT the spectrum.

    # ---- The fundamental asymmetry ----
    # S³ (Coulomb): (D-A) computes the spectrum. Graph Laplacian IS the kinetic operator.
    # S⁵ (HO): (D-A) does NOT compute the spectrum. Hamiltonian is diagonal.
    #
    # This means there's NO natural "κ_S5" in the same sense as κ_S3.
    # The S⁵ graph encodes transition structure, not energy eigenvalues.

    # Let's quantify this by checking if (D-A) on S⁵ has any spectral
    # relationship to N(N+4) or N+3/2.

    from geovac.nuclear.bargmann_graph import build_bargmann_graph

    results_by_Nmax = {}

    for N_max in [2, 3, 4]:
        g = build_bargmann_graph(N_max)
        A = g.adjacency_dense()
        L = g.graph_laplacian_dense()
        evals_L = np.sort(np.linalg.eigvalsh(L))

        # HO energies (diagonal)
        H_diag = g.hamiltonian_dense(hw=1.0)
        ho_evals = np.sort(np.diag(H_diag))

        # S⁵ LB eigenvalues
        s5_lb = []
        for N in range(N_max + 1):
            deg = (N + 1) * (N + 2) // 2
            for _ in range(deg):
                s5_lb.append(float(N * (N + 4)))
        s5_lb = np.sort(s5_lb)

        # Test: κ_pred * L eigenvalues vs HO energies
        kappa_pred = -1.0 / 36.0
        kappa_L = kappa_pred * evals_L

        # Test: κ_pred * L eigenvalues vs S⁵ LB eigenvalues
        # Since one L eigenvalue is 0, remove it for comparison
        evals_L_nz = evals_L[evals_L > 1e-10]
        s5_lb_nz = s5_lb[s5_lb > 1e-10]

        # Compute the BEST κ that maps L evals to S⁵ LB
        if len(evals_L_nz) > 0 and len(s5_lb_nz) > 0:
            n_match = min(len(evals_L_nz), len(s5_lb_nz))
            kappa_best_LB = np.dot(s5_lb_nz[:n_match], evals_L_nz[:n_match]) / \
                            np.dot(evals_L_nz[:n_match], evals_L_nz[:n_match])
            res_LB = np.sqrt(np.mean((s5_lb_nz[:n_match] - kappa_best_LB * evals_L_nz[:n_match])**2))
            rel_res_LB = res_LB / np.mean(s5_lb_nz[:n_match])
        else:
            kappa_best_LB = None
            res_LB = None
            rel_res_LB = None

        # Compute the BEST κ that maps L evals to HO energies
        # Remove the zero from L
        # Actually HO energies include the ground state 3/2, so include all
        n_match_ho = min(len(evals_L), len(ho_evals))
        if n_match_ho > 0 and np.dot(evals_L[:n_match_ho], evals_L[:n_match_ho]) > 1e-15:
            kappa_best_HO = np.dot(ho_evals[:n_match_ho], evals_L[:n_match_ho]) / \
                            np.dot(evals_L[:n_match_ho], evals_L[:n_match_ho])
            res_HO = np.sqrt(np.mean((ho_evals[:n_match_ho] - kappa_best_HO * evals_L[:n_match_ho])**2))
            rel_res_HO = res_HO / np.mean(ho_evals[:n_match_ho])
        else:
            kappa_best_HO = None
            res_HO = None
            rel_res_HO = None

        results_by_Nmax[f"N_max={N_max}"] = {
            "V": g.n_nodes,
            "E": len(g.adjacency),
            "L_evals": evals_L.tolist(),
            "HO_evals": ho_evals.tolist(),
            "S5_LB_evals": s5_lb.tolist(),
            "kappa_pred_L_evals": kappa_L.tolist(),
            "kappa_best_fit_to_LB": float(kappa_best_LB) if kappa_best_LB is not None else None,
            "residual_to_LB": float(res_LB) if res_LB is not None else None,
            "rel_residual_to_LB": float(rel_res_LB) if rel_res_LB is not None else None,
            "kappa_best_fit_to_HO": float(kappa_best_HO) if kappa_best_HO is not None else None,
            "residual_to_HO": float(res_HO) if res_HO is not None else None,
            "rel_residual_to_HO": float(rel_res_HO) if rel_res_HO is not None else None,
        }

    return results_by_Nmax


# ============================================================================
# Alternative decompositions of -1/16
# ============================================================================

def alternative_decompositions() -> Dict:
    """
    All known ways to write -1/16:

    (a) -1/d_max² where d_max = 4 = dim(R⁴) = ambient dimension of S³
    (b) -1/(2*N_init)² where N_init = 2 (Paper 0 initial state count)
    (c) -(d₀/2)⁴ where d₀ = 1 (packing diameter) — i.e. -(1/2)⁴
    (d) -1/(2*(l_max+1))² with l_max = 1 (Paper 2 cutoff)
    (e) -1/(2*dim(S³))² = -1/(2*3)² = -1/36... NO, that gives -1/36.
    (f) -(1/2)^(d_max) = -(1/2)^4 = -1/16
    (g) -1/(n_max_Paper0)^(d_max) where n_max=2, d_max=2: -1/4... no.

    Actually: d_max = 4 appears in Paper 0 as the dimension where the
    packing construction terminates. Paper 0's d_max is NOT literally
    dim(R⁴) — it's determined by the packing axiom "pack Planck's
    constant into disks."

    The coincidence is: d_max_packing = 4 = dim(ambient space of S³).
    The question is whether this is a deep structural relationship.
    """
    decompositions = {
        "kappa": -1/16,
        "decompositions": {
            "(a) -1/d_max^2 (d_max=4=dim(R^4))": {
                "formula": "-1/4^2",
                "value": -1/16,
                "match": True,
                "d_max_meaning": "packing dimension = ambient dimension of S^3"
            },
            "(b) -1/(2*N_init)^2 (N_init=2)": {
                "formula": "-1/(2*2)^2",
                "value": -1/16,
                "match": True,
                "meaning": "N_init = Paper 0 initial state count (2 states)"
            },
            "(c) -(1/2)^4": {
                "formula": "-(1/2)^4",
                "value": -1/16,
                "match": True,
                "meaning": "Half^(packing_dim) or (d_0/2)^d_max"
            },
            "(d) -1/(2*(l_max+1))^2 (l_max=1)": {
                "formula": "-1/(2*2)^2",
                "value": -1/16,
                "match": True,
                "meaning": "l_max=1 is Paper 2's angular cutoff — same as (b)"
            },
            "(e) -1/(dim(S^3)+1)^2 = -1/4^2": {
                "formula": "-1/4^2",
                "value": -1/16,
                "match": True,
                "meaning": "dim(S^d) + 1 = dim(R^{d+1}) = d_max — same as (a)"
            },
        },
        "note": "All decompositions reduce to either d_max=4 or N_init=2. "
                "The coincidence d_max = 2*N_init = dim(R^4) = dim(S^3)+1 "
                "makes it impossible to distinguish these from the S^3 case alone. "
                "The S^5 test is the discriminant: "
                "d_max(S^5) = 6, N_init = 2, dim(S^5)+1 = 6."
    }

    # S⁵ predictions under each hypothesis
    decompositions["s5_predictions"] = {
        "(a) d_max=6=dim(R^6)": {"kappa_S5": -1/36, "formula": "-1/6^2"},
        "(b) N_init=2 universal": {"kappa_S5": -1/16, "formula": "-1/(2*2)^2 (universal)"},
        "(c) (1/2)^6": {"kappa_S5": -1/64, "formula": "-(1/2)^6"},
        "(d) l_max=1 universal": {"kappa_S5": -1/16, "formula": "-1/4^2 (universal)"},
    }

    return decompositions


# ============================================================================
# Master analysis
# ============================================================================

def main():
    out_path = Path(__file__).parent / "data" / "probe_k4_dmax_structural.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Probe: Is kappa = -1/d_max^2 structural across sphere dimensions?")
    print("=" * 70)

    # 1. S³ sector decomposition
    print("\n--- S³ Coulomb graph sector analysis ---")
    s3_sectors = analyze_s3_graph_eigenvalues(n_max=4)
    print(f"  Sectors match path graph: {s3_sectors['sectors_match_path_graph']}")
    print(f"  Max diff from path graph: {s3_sectors['max_sector_path_diff']:.2e}")

    # 2. S⁵ Bargmann-Segal graph analysis
    print("\n--- S⁵ Bargmann-Segal graph analysis ---")
    s5_analysis = analyze_s5_graph_eigenvalues(N_max=3)
    print(f"  V = {s5_analysis['V']}, E = {s5_analysis['E']}")
    print(f"  Number of zero L eigenvalues: {s5_analysis['n_zero_L']}")
    print(f"  Number of zero S⁵ LB eigenvalues: {s5_analysis['n_zero_LB']}")
    print(f"  Scale fit L → S⁵ LB: {s5_analysis['scale_fit_L_to_LB']}")
    print(f"  Residual L → S⁵ LB: {s5_analysis['residual_L_to_LB']}")
    print(f"  Max rel error L → S⁵ LB: {s5_analysis['max_rel_error_L_to_LB']}")

    print(f"\n  L eigenvalues: {[f'{x:.4f}' for x in s5_analysis['L_evals_sorted']]}")
    print(f"  S⁵ LB evals:   {[f'{x:.1f}' for x in s5_analysis['s5_LB_evals_sorted']]}")

    # 3. Definitive kappa test across N_max
    print("\n--- Definitive kappa test (S^5) ---")
    definitive = definitive_kappa_test()
    for key, data in definitive.items():
        print(f"\n  {key}:")
        print(f"    V={data['V']}, E={data['E']}")
        print(f"    kappa_best_fit to S^5 LB: {data['kappa_best_fit_to_LB']}")
        print(f"    Rel residual to S^5 LB: {data['rel_residual_to_LB']}")
        print(f"    kappa_best_fit to HO: {data['kappa_best_fit_to_HO']}")
        print(f"    Rel residual to HO: {data['rel_residual_to_HO']}")

    # 4. Alternative decompositions
    print("\n--- Alternative decompositions of -1/16 ---")
    alts = alternative_decompositions()
    for name, data in alts["decompositions"].items():
        print(f"  {name}: {data['formula']} = {data['value']}")
    print(f"\n  S⁵ predictions:")
    for name, data in alts["s5_predictions"].items():
        print(f"    {name}: κ_S5 = {data['kappa_S5']}")

    # 5. Key structural finding
    print("\n" + "=" * 70)
    print("STRUCTURAL ANALYSIS")
    print("=" * 70)

    # The critical finding: does the S⁵ graph Laplacian L have eigenvalues
    # proportional to the S⁵ LB eigenvalues N(N+4)?

    # Let's check directly with exact edge weights at N_max=3
    from geovac.nuclear.bargmann_graph import build_bargmann_graph
    g = build_bargmann_graph(3)
    A = g.adjacency_dense()
    L = g.graph_laplacian_dense()
    evals_L = np.sort(np.linalg.eigvalsh(L))

    # S⁵ LB eigenvalues at N_max=3: 0(x1), 5(x3), 12(x6), 21(x10)
    s5_lb_3 = [0.0] + [5.0]*3 + [12.0]*6 + [21.0]*10
    s5_lb_3 = np.sort(s5_lb_3)

    print(f"\nS⁵ graph at N_max=3 (V={g.n_nodes}, E={len(g.adjacency)}):")
    print(f"  L eigenvalues:      {[f'{x:.6f}' for x in evals_L]}")
    print(f"  S⁵ LB eigenvalues:  {[f'{x:.1f}' for x in s5_lb_3]}")

    # Check: how many distinct L eigenvalues? How many distinct LB eigenvalues?
    unique_L = np.unique(np.round(evals_L, 8))
    unique_LB = np.unique(s5_lb_3)
    print(f"\n  Distinct L eigenvalues:  {len(unique_L)}")
    print(f"  Distinct S⁵ LB eigenvalues: {len(unique_LB)}")
    print(f"  Distinct L values: {[f'{x:.6f}' for x in unique_L]}")
    print(f"  Distinct LB values: {[f'{x:.1f}' for x in unique_LB]}")

    # Check degeneracies of L
    print(f"\n  L eigenvalue degeneracies:")
    for val in unique_L:
        count = np.sum(np.abs(evals_L - val) < 1e-6)
        print(f"    λ = {val:.6f}: degeneracy {count}")

    # If L eigenvalues DON'T have the same degeneracy structure as S⁵ LB,
    # then L ≠ c * (Laplace-Beltrami on S⁵).

    # Check the NUMBER of distinct eigenvalues
    # S⁵ LB at N_max=3 has 4 distinct eigenvalues (0, 5, 12, 21)
    # with degeneracies (1, 3, 6, 10)

    # CRITICAL CHECK: does the graph Laplacian have the SAME number of
    # distinct eigenvalues as the LB operator?
    same_distinct_count = len(unique_L) == len(unique_LB)

    print(f"\n  Same number of distinct eigenvalues: {same_distinct_count}")

    # If not, then L is NOT proportional to LB, and κ = -1/d_max² cannot
    # work for S⁵ in the same way it works for S³.

    # ---- The S³ comparison ----
    # For S³ at n_max=3 (14 nodes):
    L_s3, nodes_s3 = build_s3_fock_graph_laplacian(3)
    evals_s3 = np.sort(np.linalg.eigvalsh(L_s3))
    unique_s3 = np.unique(np.round(evals_s3, 8))

    # S³ LB eigenvalues: n²-1 with degeneracy n²
    # At n_max=3: 0(x1), 3(x4), 8(x9)
    s3_lb = [0.0] + [3.0]*4 + [8.0]*9
    s3_lb = np.sort(s3_lb)
    unique_s3_lb = np.unique(s3_lb)

    print(f"\nS³ graph at n_max=3 (V={len(nodes_s3)}):")
    print(f"  L eigenvalues:      {[f'{x:.6f}' for x in evals_s3]}")
    print(f"  S³ LB eigenvalues:  {[f'{x:.1f}' for x in s3_lb]}")
    print(f"  Distinct L values: {[f'{x:.6f}' for x in unique_s3]}")
    print(f"  Distinct LB values: {[f'{x:.1f}' for x in unique_s3_lb]}")

    # S³ L eigenvalue degeneracies
    print(f"\n  S³ L eigenvalue degeneracies:")
    for val in unique_s3:
        count = np.sum(np.abs(evals_s3 - val) < 1e-6)
        print(f"    λ = {val:.6f}: degeneracy {count}")

    # Check: do S³ L eigenvalues match S³ LB eigenvalues?
    s3_L_matches_LB = (len(unique_s3) == len(unique_s3_lb))
    if s3_L_matches_LB:
        # Check proportionality
        ratios_s3 = unique_s3[1:] / unique_s3_lb[1:]
        s3_proportional = np.allclose(ratios_s3, ratios_s3[0], rtol=1e-6)
        s3_scale = ratios_s3[0] if s3_proportional else None
    else:
        s3_proportional = False
        s3_scale = None

    print(f"\n  S³ L proportional to S³ LB: {s3_proportional}")
    if s3_scale is not None:
        print(f"  Scale factor: {s3_scale:.6f}")
        print(f"  κ would be: {-1.0/(s3_scale * 2):.6f}")  # E = -Z²/(2n²), so need more care

    # ---- VERDICT ----
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)

    # The key structural finding from Paper 24 (quoted in CLAUDE.md):
    # "The graph adjacency encodes dipole transitions only, NOT the spectrum
    # (in contrast to the Coulomb S³ case where (D-A) computes the spectrum)"
    #
    # This means there IS NO natural κ_S5 in the same sense as κ_S3.
    # κ = -1/16 maps D-A to the kinetic energy on S³ because the Fock
    # conformal map makes the Laplacian encode the spectrum.
    # On S⁵, the Hamiltonian is diagonal (not from D-A), so there's no
    # analogous κ.

    # But let's still check: does the L spectrum have ANY structure matching LB?

    s5_L_matches_LB_structure = same_distinct_count
    if same_distinct_count:
        # Check degeneracy match
        s5_L_degs = [int(np.sum(np.abs(evals_L - val) < 1e-6)) for val in unique_L]
        s5_LB_degs = [int(np.sum(np.abs(s5_lb_3 - val) < 1e-6)) for val in unique_LB]
        degs_match = (s5_L_degs == s5_LB_degs)
    else:
        degs_match = False
        s5_L_degs = [int(np.sum(np.abs(evals_L - val) < 1e-6)) for val in unique_L]
        s5_LB_degs = [int(np.sum(np.abs(s5_lb_3 - val) < 1e-6)) for val in unique_LB]

    # Compute what κ_S5 would need to be if we FORCED L → LB
    if s5_L_matches_LB_structure and len(unique_L) > 1 and len(unique_LB) > 1:
        forced_ratios = unique_LB[1:] / unique_L[1:]
        forced_kappa = 1.0 / forced_ratios[0]  # L * forced_kappa ≈ LB
    else:
        forced_ratios = None
        forced_kappa = None

    verdict = {
        "s3_L_matches_LB_degeneracy": s3_proportional,
        "s3_L_scale_to_LB": float(s3_scale) if s3_scale is not None else None,
        "s5_L_distinct_count": len(unique_L),
        "s5_LB_distinct_count": len(unique_LB),
        "s5_L_matches_LB_distinct_count": same_distinct_count,
        "s5_L_degeneracies": s5_L_degs,
        "s5_LB_degeneracies": s5_LB_degs,
        "s5_degeneracies_match": degs_match,
        "structural_verdict": None,  # filled below
        "explanation": None,
    }

    if not same_distinct_count or not degs_match:
        verdict["structural_verdict"] = "NEGATIVE — the hypothesis κ = -1/d_max² does NOT transfer to S⁵"
        verdict["explanation"] = (
            f"The S⁵ Bargmann-Segal graph Laplacian L has {len(unique_L)} distinct eigenvalues "
            f"with degeneracies {s5_L_degs}, while the S⁵ Laplace-Beltrami operator has "
            f"{len(unique_LB)} distinct eigenvalues N(N+4) with degeneracies {s5_LB_degs}. "
            f"The degeneracy structures do not match, so L is NOT proportional to the "
            f"Laplace-Beltrami operator on S⁵. This is consistent with Paper 24's finding: "
            f"the S⁵ graph adjacency encodes dipole transitions, not the spectrum. "
            f"The S³ case is special because the Fock conformal map makes D-A compute the "
            f"spectrum directly. No such conformal map exists for S⁵ (Paper 23's rigidity theorem). "
            f"Therefore κ = -1/d_max² is a COULOMB-SPECIFIC relationship tied to the Fock "
            f"projection, not a universal structural law across sphere dimensions."
        )
    else:
        verdict["structural_verdict"] = "NEEDS FURTHER ANALYSIS — degeneracy structures match"

    print(f"\n  {verdict['structural_verdict']}")
    print(f"\n  {verdict['explanation']}")

    # ---- Save results ----
    all_results = {
        "probe": "κ = -1/d_max² structural test",
        "date": "2026-04-19",
        "s3_sector_analysis": {
            "sectors_match_path_graph": s3_sectors["sectors_match_path_graph"],
            "max_diff": s3_sectors["max_sector_path_diff"],
        },
        "s3_vs_LB": {
            "n_max": 3,
            "L_distinct_evals": unique_s3.tolist(),
            "LB_distinct_evals": unique_s3_lb.tolist(),
            "proportional": s3_proportional,
            "scale": float(s3_scale) if s3_scale is not None else None,
        },
        "s5_analysis": {
            "N_max": 3,
            "V": g.n_nodes,
            "E": len(g.adjacency),
            "L_evals": evals_L.tolist(),
            "L_distinct_evals": unique_L.tolist(),
            "L_degeneracies": s5_L_degs,
            "LB_distinct_evals": unique_LB.tolist(),
            "LB_degeneracies": s5_LB_degs,
            "A_evals": np.sort(np.linalg.eigvalsh(A)).tolist(),
        },
        "definitive_test": definitive,
        "alternative_decompositions": alts,
        "verdict": verdict,
    }

    with out_path.open("w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print(f"\nWrote: {out_path}")

    return all_results


def supplementary_s3_check() -> Dict:
    """
    The CORRECT test for S^3: check that H = kappa*(D-A) + V gives
    the hydrogen spectrum -Z^2/(2n^2).

    V_ii = -Z/n^2 is the node weight (potential energy).
    kappa = -1/16 is the kinetic scale.
    The eigenvalues of H should be E_n = -Z^2/(2n^2) with degeneracy n^2.
    """
    results = {}
    for n_max in [2, 3, 4, 5]:
        L, nodes = build_s3_fock_graph_laplacian(n_max)
        V_count = len(nodes)

        # Potential energy diagonal: V_ii = -Z/n^2 (Z=1)
        V = np.zeros(V_count)
        for i, (n, l, m) in enumerate(nodes):
            V[i] = -1.0 / (n * n)

        # H = kappa * L + V
        kappa = -1.0 / 16.0
        H = kappa * L + np.diag(V)
        evals_H = np.sort(np.linalg.eigvalsh(H))

        # Exact hydrogen energies
        exact = []
        for n in range(1, n_max + 1):
            for _ in range(n * n):
                exact.append(-0.5 / (n * n))
        exact = np.sort(exact)

        # Error
        max_err = np.max(np.abs(evals_H - exact))
        rms_err = np.sqrt(np.mean((evals_H - exact)**2))

        results[f"n_max={n_max}"] = {
            "V": V_count,
            "H_evals": evals_H.tolist(),
            "exact_evals": exact,
            "max_abs_error": float(max_err),
            "rms_error": float(rms_err),
            "kappa_used": kappa,
        }

        # Also try WITHOUT the diagonal V (just kappa*L)
        H_kinetic = kappa * L
        evals_kinetic = np.sort(np.linalg.eigvalsh(H_kinetic))

        results[f"n_max={n_max}"]["kinetic_only_evals"] = evals_kinetic.tolist()

    return results


def supplementary_s5_check() -> Dict:
    """
    Try to find a kappa for S^5 such that kappa*(D-A) + diagonal(HO)
    gives the correct S^5 LB eigenvalues N(N+4), or vice versa.

    Also try: kappa*(D-A) + alpha*diagonal = target_spectrum.
    """
    from geovac.nuclear.bargmann_graph import build_bargmann_graph

    results = {}
    for N_max in [2, 3, 4]:
        g = build_bargmann_graph(N_max)
        A = g.adjacency_dense()
        L = g.graph_laplacian_dense()
        H_diag = g.hamiltonian_dense(hw=1.0)
        V = g.n_nodes

        # Target 1: HO energies (N + 3/2) with correct degeneracies
        ho_target = np.sort(np.diag(H_diag))

        # Target 2: S^5 LB eigenvalues N(N+4)
        lb_target = []
        for N in range(N_max + 1):
            deg = (N + 1) * (N + 2) // 2
            for _ in range(deg):
                lb_target.append(float(N * (N + 4)))
        lb_target = np.sort(lb_target)

        # Test various kappa values
        test_kappas = [-1/36, -1/16, -1/64, -1/4, 1/36, 1/16, 0.0]
        kappa_results = {}

        for kappa in test_kappas:
            # H = kappa * L + H_diag
            H = kappa * L + H_diag
            evals_H = np.sort(np.linalg.eigvalsh(H))

            # Compare to HO target
            err_ho = np.max(np.abs(evals_H - ho_target))

            # Compare to LB target
            err_lb = np.max(np.abs(evals_H - lb_target))

            kappa_results[f"kappa={kappa:.6f}"] = {
                "H_evals": evals_H.tolist(),
                "max_err_vs_HO": float(err_ho),
                "max_err_vs_LB": float(err_lb),
            }

        # Now try H = kappa * L + alpha * H_diag
        # For the S^3 case: H = kappa*L + V, where V = -1/n^2 diagonal
        # and the eigenvalues are -1/(2n^2).
        # Note that V is NOT the spectrum -- V_ii = -1/n_i^2 is the POTENTIAL.
        # The spectrum -1/(2n^2) emerges from diagonalizing H = kappa*L + V.

        # For S^5: the "potential" analog is unclear.
        # The diagonal IS the spectrum (N+3/2).
        # There's no separate potential.

        # Try: does kappa*L alone have any relation to N(N+4)?
        evals_L = np.sort(np.linalg.eigvalsh(L))

        # For each kappa, compute kappa*L eigenvalues and check vs LB
        best_kappa_for_LB = None
        best_err_LB = 1e10
        for trial_kappa in np.linspace(-100, 100, 10001):
            trial_evals = trial_kappa * evals_L
            err = np.max(np.abs(trial_evals - lb_target))
            if err < best_err_LB:
                best_err_LB = err
                best_kappa_for_LB = trial_kappa

        results[f"N_max={N_max}"] = {
            "V": V,
            "kappa_sweep_results": kappa_results,
            "best_kappa_for_LB_spectrum": float(best_kappa_for_LB),
            "best_err_LB_spectrum": float(best_err_LB),
            "L_evals": evals_L.tolist(),
            "ho_target": ho_target.tolist(),
            "lb_target": lb_target.tolist(),
        }

    return results


def critical_reanalysis() -> Dict:
    """
    CRITICAL REANALYSIS: kappa = -1/16 is NOT kappa*(D-A) -> energies.

    The actual role of kappa in GeoVac:
    In casimir_ci._build_graph_h1:
      h1[i,i] = -Z^2/(2n^2)  (exact diagonal)
      h1[i,j] = kappa * (-A[i,j]) = +(1/16)*A[i,j]  (off-diagonal coupling)

    So kappa determines the INTER-SHELL COUPLING STRENGTH.
    The diagonal already has the exact spectrum.
    kappa scales how much the graph adjacency mixes different n-shells.

    For S^3: the adjacency A connects (n,l,m) <-> (n+1,l,m) [radial]
    and (n,l,m) <-> (n,l,m+1) [angular]. kappa = 1/16 is the coupling.

    For S^5: the adjacency A has squared dipole matrix elements as weights.
    These weights are already physical (from the r operator).
    The "coupling" is built into the weights, not into a separate kappa.

    Key test: for S^3, verify that the ADJACENCY EIGENVALUE SCALE
    determines kappa. The graph Laplacian eigenvalues of the path graph
    within each (l,m) sector are 2 - 2*cos(pi*k/(p+1)).
    The maximum eigenvalue of a path graph of length p is ~ 4 (for large p).
    The maximum eigenvalue of the S^3 adjacency is 2*degree_max.
    """
    # S^3 graph: check adjacency eigenvalues
    results = {}
    for n_max in [2, 3, 4, 5, 6]:
        L, nodes = build_s3_fock_graph_laplacian(n_max)
        V_count = len(nodes)

        # Extract adjacency from L (L = D - A, so A = D - L)
        D = np.diag(np.diag(L))
        A = D - L

        evals_A = np.sort(np.linalg.eigvalsh(A))

        # S^3 LB eigenvalues: n^2 - 1
        s3_lb = []
        for n in range(1, n_max + 1):
            for _ in range(n * n):
                s3_lb.append(n * n - 1)
        s3_lb = np.sort(s3_lb).astype(float)

        # Check: adjacency eigenvalue spread
        max_A = evals_A[-1]
        min_A = evals_A[0]
        spread_A = max_A - min_A

        # Check: max degree
        max_degree = np.max(np.diag(D))

        # The graph-native h1 is:
        # h1 = diag(-Z^2/(2n^2)) + (1/16) * A
        # = diag(exact_energies) + (1/16) * A
        #
        # If we define the perturbation parameter epsilon = 1/16,
        # then at epsilon=0, h1 = exact diagonal.
        # At epsilon=1/16, the off-diagonal coupling mixes shells.
        #
        # The EXACT hydrogenic h1 in the (n,l,m) basis is diagonal:
        # h1_exact[i,i] = -Z^2/(2n_i^2)
        # h1_exact[i,j] = 0 for i != j (eigenbasis of H)
        #
        # But the GRAPH-NATIVE h1 is NOT diagonal in the energy eigenbasis.
        # The graph adjacency A couples different n-shells.
        # kappa = -1/16 gives the RIGHT AMOUNT of coupling.
        #
        # Question: is 1/16 related to 1/d_max^2 = 1/16 for S^3?

        # Check: what eigenvalues does h1 = diag(-1/(2n^2)) + (1/16)*A give?
        h1 = np.zeros((V_count, V_count))
        for i, (n, l, m) in enumerate(nodes):
            h1[i, i] = -0.5 / (n * n)
        h1 += (1.0/16.0) * A

        evals_h1 = np.sort(np.linalg.eigvalsh(h1))

        # Exact energies
        exact = []
        for n in range(1, n_max + 1):
            for _ in range(n * n):
                exact.append(-0.5 / (n * n))
        exact = np.sort(exact)

        # Check match
        max_err = np.max(np.abs(evals_h1 - exact))

        # Now check: does the graph eigenbasis h1 from casimir_ci match?
        # h1 has diagonal = exact energies, off-diagonal = (1/16)*A
        # The eigenvalues of this h1 should be close to exact for small coupling.

        results[f"n_max={n_max}"] = {
            "V": V_count,
            "A_max_eval": float(max_A),
            "A_min_eval": float(min_A),
            "A_spread": float(spread_A),
            "max_degree": float(max_degree),
            "h1_max_error": float(max_err),
            "h1_eigenvalues": evals_h1.tolist(),
            "exact_energies": exact,
        }

    # Now for S^5: the Bargmann graph already has weighted adjacency.
    # The "coupling" IS the squared dipole matrix element.
    # There's no separate kappa.
    # The diagonal IS the HO spectrum. Adding ANY A contribution changes it.
    #
    # This confirms: kappa = -1/16 is specific to the S^3 graph construction
    # where the adjacency has binary weights and needs to be rescaled.

    return results


if __name__ == "__main__":
    results = main()

    # Supplementary checks
    print("\n" + "=" * 70)
    print("SUPPLEMENTARY: S^3 H = kappa*L + V check")
    print("=" * 70)
    s3_supp = supplementary_s3_check()
    for key, data in s3_supp.items():
        print(f"\n  {key}: max_abs_error = {data['max_abs_error']:.2e}, rms = {data['rms_error']:.2e}")
        if data['max_abs_error'] < 1e-10:
            print(f"    EXACT MATCH: kappa={data['kappa_used']} with V=-1/n^2 reproduces -1/(2n^2)")
        else:
            print(f"    NOT exact. H eigenvalues:")
            for i, (h, e) in enumerate(zip(data['H_evals'], data['exact_evals'])):
                if abs(h - e) > 1e-10:
                    print(f"      [{i}] H={h:.8f}  exact={e:.8f}  diff={h-e:.2e}")

    print("\n" + "=" * 70)
    print("SUPPLEMENTARY: S^5 kappa sweep")
    print("=" * 70)
    s5_supp = supplementary_s5_check()
    for key, data in s5_supp.items():
        print(f"\n  {key}:")
        print(f"    Best kappa for L->LB spectrum: {data['best_kappa_for_LB_spectrum']:.6f}")
        print(f"    Best max error: {data['best_err_LB_spectrum']:.4f}")
        for k2, d2 in data['kappa_sweep_results'].items():
            print(f"    {k2}: err vs HO={d2['max_err_vs_HO']:.4f}, err vs LB={d2['max_err_vs_LB']:.4f}")

    print("\n" + "=" * 70)
    print("CRITICAL REANALYSIS: kappa as off-diagonal adjacency scale")
    print("=" * 70)
    crit = critical_reanalysis()
    for key, data in crit.items():
        print(f"\n  {key} (V={data['V']}):")
        print(f"    Adjacency eigenvalue range: [{data['A_min_eval']:.4f}, {data['A_max_eval']:.4f}]")
        print(f"    Adjacency spread: {data['A_spread']:.4f}")
        print(f"    Max degree: {data['max_degree']:.0f}")
        print(f"    h1 = diag(-1/(2n^2)) + (1/16)*A max error: {data['h1_max_error']:.6e}")

    # Update saved results
    results["supplementary_s3"] = s3_supp
    results["supplementary_s5"] = s5_supp
    results["critical_reanalysis"] = crit

    out_path = Path(__file__).parent / "data" / "probe_k4_dmax_structural.json"
    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nUpdated: {out_path}")
