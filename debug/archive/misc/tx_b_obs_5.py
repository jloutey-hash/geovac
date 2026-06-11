"""TX-B Observable 5: Static Coulomb-like potential V_ab on Hopf S^3 graph
from Green's function sum at fixed nodes.

POSITIVE PREDICTION (no pi).

The discrete Green's function on the Fock S^3 graph for the SHIFTED Laplacian
(L + epsilon I) at two fixed nodes a, b is

    G_ab(epsilon) = sum_{n: lambda_n + epsilon != 0} psi_n(a) psi_n*(b)
                                                    / (lambda_n + epsilon)

where {(lambda_n, psi_n)} are the graph Laplacian eigenpairs.  We use a small
positive epsilon to regularize the zero mode (or, equivalently, the
pseudoinverse/Moore-Penrose Green's function with the constant mode projected
out).

This is the discrete analog of the static Coulomb-like potential V(theta_1,
theta_2) = G(theta_1, theta_2) on S^3, evaluated at two fixed angular positions
encoded as graph nodes.

PREDICTION: The values are algebraic (rational + sqrt extensions of integers),
NO pi.  The graph is finite, all matrix entries are integer/rational, all
eigenvalues are algebraic, all eigenvector entries are algebraic, the
Green's function sum is a finite sum of algebraic terms.

Mechanism per Paper 35: there is NO continuous integration.  The "spectral
sum" sum_n psi(a) psi(b)/lambda is a discrete sum, not an integral.  Pi
should NOT appear.

We verify by:
 1. Building the Fock S^3 graph at n_max = 2, 3 using geovac.lattice.
 2. Computing the regularized Green's function matrix (L + epsilon I)^{-1}
    in exact rational arithmetic via sympy.
 3. Reading off entries between specified node pairs and checking
    transcendental content.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
import numpy as np

OUT = Path(__file__).parent / "data" / "tx_b_obs_5.json"


def build_fock_s3_adjacency_sympy(n_max: int):
    """Build the Fock S^3 graph adjacency with raising/lowering operators
    in EXACT sympy Rational arithmetic (binary edge weights = 1).

    Nodes are (n, l, m) with n in 1..n_max, l in 0..n-1, m in -l..l.
    Edges (binary, undirected, no self-loops):
      - Angular: (n,l,m) <-> (n,l,m+/-1) if |m+/-1| <= l
      - Radial:  (n,l,m) <-> (n+/-1,l,m) if l <= (n+/-1)-1 (i.e., l < n_other)
    """
    states = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    N = len(states)
    idx = {s: i for i, s in enumerate(states)}

    A = sp.zeros(N, N)
    edges = 0
    for n, l, m in states:
        i = idx[(n, l, m)]
        # Angular raising/lowering
        for dm in (-1, +1):
            if abs(m + dm) <= l:
                if (n, l, m + dm) in idx:
                    j = idx[(n, l, m + dm)]
                    if A[i, j] == 0:
                        A[i, j] = 1
                        A[j, i] = 1
                        edges += 1
        # Radial raising/lowering
        for dn in (-1, +1):
            n2 = n + dn
            if 1 <= n2 <= n_max and l < n2:
                if (n2, l, m) in idx:
                    j = idx[(n2, l, m)]
                    if A[i, j] == 0:
                        A[i, j] = 1
                        A[j, i] = 1
                        edges += 1
    return states, A, edges


def graph_laplacian(A):
    """L = D - A in sympy."""
    deg = [sum(A[i, :]) for i in range(A.rows)]
    D = sp.diag(*deg)
    return D - A


def regularized_greens_function(L, epsilon):
    """Compute (L + epsilon I)^{-1} in sympy exact arithmetic.

    With epsilon = 1 (rational), the Laplacian shift is invertible because
    L is positive semidefinite (kernel = constants), and (L + I)^{-1} is
    well-defined.

    Returns the full N x N inverse matrix.
    """
    N = L.rows
    I = sp.eye(N)
    M = L + epsilon * I
    return M.inv()


def expr_contains_pi(expr) -> bool:
    """Check if a sympy expression contains pi anywhere."""
    return sp.pi in expr.atoms() or any(sp.pi in a.atoms() for a in expr.atoms()
                                         if hasattr(a, 'atoms'))


def main():
    results_per_nmax = []

    for n_max in [2, 3]:
        states, A, n_edges = build_fock_s3_adjacency_sympy(n_max)
        N = len(states)
        L = graph_laplacian(A)

        # Regularized Green's function with epsilon = 1 (exact rational)
        G = regularized_greens_function(L, sp.Rational(1, 1))

        # Sample several node pairs and check transcendental content
        # Strategy: pick nodes (1,0,0), (n_max,0,0), (n_max,1,0), etc.
        sample_pairs = []
        for i, sa in enumerate(states):
            for j, sb in enumerate(states):
                if i <= j:
                    sample_pairs.append((i, j, sa, sb))

        # Check ALL entries
        all_entries_rational = True
        any_entry_has_pi = False
        sample_values = []
        for (i, j, sa, sb) in sample_pairs[:20]:    # log first 20 for record
            val = sp.simplify(G[i, j])
            is_rat = isinstance(val, sp.Rational)
            has_pi = expr_contains_pi(val)
            if not is_rat:
                # Could be in algebraic-extension ring (e.g., contain sqrt)
                # Check whether it's an algebraic number
                try:
                    is_algebraic = val.is_algebraic
                except AttributeError:
                    is_algebraic = None
            else:
                is_algebraic = True
            if has_pi:
                any_entry_has_pi = True
            sample_values.append({
                "node_a": str(sa),
                "node_b": str(sb),
                "G_ab": str(val),
                "is_rational": bool(is_rat),
                "contains_pi": bool(has_pi),
            })

        # Full check: scan entire matrix for pi
        for i in range(N):
            for j in range(N):
                if expr_contains_pi(G[i, j]):
                    any_entry_has_pi = True
                    break
            if any_entry_has_pi:
                break

        # Also check whether the Laplacian eigenvalues themselves are algebraic
        # (they always are -- they're roots of a polynomial with integer
        # coefficients), but check for pi.
        eigenvals = list(L.eigenvals().keys())
        any_eig_has_pi = any(expr_contains_pi(ev) for ev in eigenvals)

        results_per_nmax.append({
            "n_max": n_max,
            "N_nodes": N,
            "n_edges": n_edges,
            "graph_eigenvalue_pi_check": "no pi in any L eigenvalue",
            "any_eigenvalue_contains_pi": bool(any_eig_has_pi),
            "first_5_eigenvalues": [str(ev) for ev in eigenvals[:5]],
            "regularization_epsilon": "1",
            "any_G_entry_contains_pi": bool(any_entry_has_pi),
            "sample_G_entries_first_20": sample_values,
        })

    # Determine final verdict
    pi_predicted = False
    pi_observed = any(r["any_G_entry_contains_pi"] or r["any_eigenvalue_contains_pi"]
                      for r in results_per_nmax)

    out = {
        "observable_id": 5,
        "observable_name": "Static Coulomb-like potential on Hopf S^3 graph (regularized Green's function)",
        "predicted_class": "rational + algebraic-extension (sqrt of integers); no pi",
        "pi_predicted": pi_predicted,
        "computed_class": ("algebraic, no pi" if not pi_observed else "contains pi"),
        "pi_observed": pi_observed,
        "match_with_prediction": (pi_predicted == pi_observed),
        "construction_note": ("Regularized graph Green's function (L + I)^{-1} in "
                               "exact sympy rational arithmetic.  No continuous "
                               "integration.  No temporal compactification."),
        "results_per_nmax": results_per_nmax,
        "verdict": ("CONFIRMS Paper 35 Prediction 1: NO pi appears because "
                     "the discrete Green's function involves no continuous integration."
                     if not pi_observed
                     else "REFUTES Paper 35 Prediction 1: pi appears in a discrete sum, "
                          "violating the prediction"),
        "paper_34_projections": [
            "Fock conformal projection (Paper 7 S^3 graph)",
            "Discrete Green's function evaluation (no continuous integration)"
        ],
        "temporal_integration_present": False,
        "specific_integration": "None -- regularized matrix inverse over a finite "
                                "graph, equivalent to a discrete sum over eigenpairs",
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps({
        "observable": "Hopf S^3 graph regularized Green's function",
        "predicted_pi": pi_predicted,
        "observed_pi": pi_observed,
        "n_max_tested": [r["n_max"] for r in results_per_nmax],
        "match": (pi_predicted == pi_observed),
    }, indent=2))


if __name__ == "__main__":
    main()
