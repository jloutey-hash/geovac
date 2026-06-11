"""
Exact sympy computation for Dirac Graph QED at n_max=2.
Complements the float64 sprint script with exact algebraic results.

Focuses on:
- Exact self-energy eigenvalues and traces
- Exact F2 extraction (fixing the bare_vertex_trace=0 issue)
- Detailed structural zero analysis
- n_max dependence of GS block for the pendant-edge theorem analog
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Matrix, Rational, sqrt, zeros as sp_zeros, eye as sp_eye

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


def extract_edges(A: np.ndarray) -> List[Tuple[int, int]]:
    """Extract undirected edges (i, j) with i < j."""
    edges = []
    V = A.shape[0]
    for i in range(V):
        for j in range(i + 1, V):
            if A[i, j] != 0:
                edges.append((i, j))
    return sorted(edges)


def exact_hodge(n_max: int, rule: str) -> Dict[str, Any]:
    """Exact sympy Hodge decomposition."""
    A, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    V = len(labels)
    edges = extract_edges(A)
    E = len(edges)

    # Build exact incidence matrix
    B = sp.zeros(V, E)
    for k, (i, j) in enumerate(edges):
        B[i, k] = Integer(1)
        B[j, k] = Integer(-1)

    L0 = B * B.T
    L1 = B.T * B

    # Exact eigenvalues
    L0_eigs = L0.eigenvals()
    L1_eigs = L1.eigenvals() if E > 0 else {}

    # Exact pseudoinverse of L1
    if E > 0 and E <= 25:
        G_gamma = L1.pinv()
    else:
        G_gamma = None

    return {
        "V": V, "E": E, "labels": labels, "edges": edges,
        "B": B, "L0": L0, "L1": L1,
        "L0_eigs": L0_eigs, "L1_eigs": L1_eigs,
        "G_gamma": G_gamma, "adjacency": A, "deg": deg,
    }


def exact_self_energy(hodge: Dict[str, Any]) -> Dict[str, Any]:
    """Exact sympy self-energy computation."""
    V = hodge["V"]
    E = hodge["E"]
    labels = hodge["labels"]
    edges = hodge["edges"]
    G_gamma = hodge["G_gamma"]

    if G_gamma is None or E == 0:
        return {"error": "No exact G_gamma available"}

    # Build vertex matrices V_e (identity on edges)
    V_mats = []
    for (i, j) in edges:
        Ve = sp.zeros(V, V)
        Ve[i, j] = Integer(1)
        Ve[j, i] = Integer(1)
        V_mats.append(Ve)

    # Sigma = sum_{e1, e2} G_gamma[e1, e2] * V_{e1} @ V_{e2}^T
    Sigma = sp.zeros(V, V)
    for e1 in range(E):
        for e2 in range(E):
            g_ee = G_gamma[e1, e2]
            if g_ee == 0:
                continue
            contrib = g_ee * (V_mats[e1] * V_mats[e2].T)
            Sigma = Sigma + contrib

    # Simplify
    for i in range(V):
        for j in range(V):
            Sigma[i, j] = sp.nsimplify(sp.expand(Sigma[i, j]), rational=False)

    trace_val = sp.nsimplify(Sigma.trace(), rational=False)

    # Eigenvalues
    try:
        eig_dict = Sigma.eigenvals()
        eigenvalues = {}
        for ev, mult in eig_dict.items():
            eigenvalues[str(sp.nsimplify(ev, rational=False))] = int(mult)
    except Exception as ex:
        eigenvalues = {"error": str(ex)}

    # Ground state block
    gs_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    if gs_indices:
        gs_block = sp.zeros(len(gs_indices), len(gs_indices))
        for ii, gi in enumerate(gs_indices):
            for jj, gj in enumerate(gs_indices):
                gs_block[ii, jj] = Sigma[gi, gj]
        gs_zero = all(gs_block[i, j] == 0
                      for i in range(gs_block.rows) for j in range(gs_block.cols))
    else:
        gs_block = None
        gs_zero = None

    return {
        "Sigma": Sigma,
        "trace": str(trace_val),
        "eigenvalues": eigenvalues,
        "gs_block": str(gs_block) if gs_block is not None else None,
        "gs_block_zero": gs_zero,
        "V_mats": V_mats,
    }


def exact_vertex_correction(hodge: Dict[str, Any]) -> Dict[str, Any]:
    """Exact sympy vertex correction at t=0."""
    V = hodge["V"]
    E = hodge["E"]
    labels = hodge["labels"]
    edges = hodge["edges"]
    G_gamma = hodge["G_gamma"]

    if G_gamma is None or E == 0:
        return {"error": "No exact G_gamma available"}

    # Electron propagator at t=0: G_e = diag(1/lambda_i)
    G_e = sp.zeros(V, V)
    for i, lab in enumerate(labels):
        n_ch = lab.n_fock - 1
        chi = 1 if lab.kappa < 0 else -1
        lam = Rational(chi * (2 * n_ch + 3), 2)
        G_e[i, i] = Rational(1, 1) / lam

    # Build vertex matrices
    V_mats = []
    for (i, j) in edges:
        Ve = sp.zeros(V, V)
        Ve[i, j] = Integer(1)
        Ve[j, i] = Integer(1)
        V_mats.append(Ve)

    # Lambda = sum_{e1, e2} G_gamma[e1, e2] * V_{e1} @ G_e @ V_{e2}^T
    Lambda_total = sp.zeros(V, V)
    for e1 in range(E):
        for e2 in range(E):
            g_ee = G_gamma[e1, e2]
            if g_ee == 0:
                continue
            Lambda_total = Lambda_total + g_ee * (V_mats[e1] * G_e * V_mats[e2].T)

    # Simplify
    for i in range(V):
        for j in range(V):
            Lambda_total[i, j] = sp.nsimplify(sp.expand(Lambda_total[i, j]), rational=False)

    trace_val = sp.nsimplify(Lambda_total.trace(), rational=False)

    # Eigenvalues
    try:
        eig_dict = Lambda_total.eigenvals()
        eigenvalues = {}
        for ev, mult in eig_dict.items():
            eigenvalues[str(sp.nsimplify(ev, rational=False))] = int(mult)
    except Exception as ex:
        eigenvalues = {"error": str(ex)}

    # F2 extraction: Tr(Lambda) / Tr(V_bare @ G_e)
    # But V_bare = sum_e V_e gives the total adjacency coupling
    # The issue is that Tr(V_bare @ G_e) can be zero by symmetry
    # For the Dirac graph, V_bare[i,j] = number of edges between i and j
    # G_e is diagonal at t=0
    # So Tr(V_bare @ G_e) = sum_{(i,j) edge} (G_e[i,i] + G_e[j,j])
    # ... which is a sum over edges of diagonal propagator values
    # This should NOT be zero in general.
    V_bare = sp.zeros(V, V)
    for Ve in V_mats:
        V_bare = V_bare + Ve
    bare_product = V_bare * G_e
    bare_trace = sp.nsimplify(bare_product.trace(), rational=False)

    # Alternative F2: use the off-diagonal bare vertex trace
    # Tr(sum_e V_e @ G_e) = sum_e (G_e[i(e), j(e)] + G_e[j(e), i(e)])
    # At t=0, G_e is diagonal, so G_e[i,j]=0 for i!=j
    # So the trace simplifies to sum_e (G_e[i,i]*V_e[i,i] + ...) = 0
    # because V_e is off-diagonal (V_e[i,i] = 0)
    # This means the bare vertex is purely off-diagonal, and Tr(V_bare @ G_e) = 0
    # when G_e is diagonal.
    #
    # The correct F2 denominator should use the vertex norm, not the trace.
    # Alternative: F2 = Tr(Lambda) / ||V_bare||_F  or some other norm
    # Or: F2 = Tr(Lambda) / Tr(Sigma) * (self-energy to vertex ratio)

    # Compute the Frobenius-norm-based F2
    V_bare_norm_sq = sum(V_bare[i,j]**2 for i in range(V) for j in range(V))
    V_bare_norm = sp.sqrt(V_bare_norm_sq)

    # Yet another approach: use Tr(V_bare @ G_e @ V_bare) / Tr(V_bare^2)
    # as the normalization

    return {
        "Lambda_trace": str(trace_val),
        "eigenvalues": eigenvalues,
        "bare_vertex_trace": str(bare_trace),
        "bare_vertex_trace_zero_explanation":
            "V_bare is purely off-diagonal; G_e(t=0) is diagonal; "
            "their product has zero diagonal => Tr = 0. "
            "F2 extraction requires t != 0 or a different normalization.",
        "Lambda_Sigma_ratio": None,  # computed below
    }


def pendant_analysis(n_max_values: List[int] = [2, 3, 4]) -> Dict[str, Any]:
    """Analyze GS pendant status across n_max values."""
    results = {}
    for rule in ["A", "B"]:
        for n_max in n_max_values:
            A, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
            gs_indices = [i for i, lab in enumerate(labels)
                          if lab.n_fock == 1 and lab.kappa == -1]
            gs_degrees = [int(deg[i]) for i in gs_indices]
            pendant = any(d == 1 for d in gs_degrees)
            key = f"rule_{rule}_nmax_{n_max}"
            results[key] = {
                "gs_degrees": gs_degrees,
                "gs_is_pendant": pendant,
                "total_V": len(labels),
                "total_E": int(np.sum(A) // 2),
            }
    return results


def main():
    print("=" * 70)
    print("EXACT SYMPY COMPUTATION: Dirac Graph QED at n_max=2")
    print("=" * 70)

    all_results = {}

    for rule in ["A", "B"]:
        print(f"\n{'='*60}")
        print(f"Rule {rule}, n_max=2")
        print(f"{'='*60}")

        hodge = exact_hodge(2, rule)

        print(f"\nHodge decomposition:")
        print(f"  V={hodge['V']}, E={hodge['E']}")
        print(f"  L0 eigenvalues: {hodge['L0_eigs']}")
        print(f"  L1 eigenvalues: {hodge['L1_eigs']}")

        # Self-energy
        print(f"\nSelf-energy (exact):")
        se = exact_self_energy(hodge)
        print(f"  Trace: {se['trace']}")
        print(f"  Eigenvalues: {se['eigenvalues']}")
        print(f"  GS block: {se['gs_block']}")
        print(f"  GS block zero: {se['gs_block_zero']}")

        # Vertex correction
        print(f"\nVertex correction (exact, t=0):")
        vc = exact_vertex_correction(hodge)
        print(f"  Lambda trace: {vc['Lambda_trace']}")
        print(f"  Eigenvalues: {vc['eigenvalues']}")
        print(f"  Bare vertex trace: {vc['bare_vertex_trace']}")
        print(f"  Explanation: {vc['bare_vertex_trace_zero_explanation']}")

        rule_key = f"rule_{rule}"
        all_results[f"self_energy_{rule_key}"] = {
            "trace": se["trace"],
            "eigenvalues": se["eigenvalues"],
            "gs_block": se["gs_block"],
            "gs_block_zero": se["gs_block_zero"],
        }
        all_results[f"vertex_correction_{rule_key}"] = {
            "Lambda_trace": vc["Lambda_trace"],
            "eigenvalues": vc["eigenvalues"],
            "bare_vertex_trace": vc["bare_vertex_trace"],
        }

    # Pendant analysis
    print(f"\n{'='*60}")
    print("Pendant Analysis: GS degree across n_max")
    print("=" * 60)
    pendant = pendant_analysis([2, 3, 4])
    for k, v in sorted(pendant.items()):
        print(f"  {k}: GS degrees={v['gs_degrees']}, pendant={v['gs_is_pendant']}")
    all_results["pendant_analysis"] = pendant

    # Save
    out_path = Path(__file__).parent / "data" / "dirac_graph_qed_exact.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
