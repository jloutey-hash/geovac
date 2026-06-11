"""Driver script for GN-1: Fock graph Hodge decomposition investigation.

Generates debug/data/gn1_fock_graph_hodge.json with full summary data
for the S^3 Fock scalar graph at n_max=1,2,3,4.

Key findings documented:
- β₀ = n_max (one connected component per angular momentum sector l=0,...,n_max-1)
- β₁ = E - V + β₀ (loop count)
- SVD identity: nonzero spectra of L₀ and L₁ are identical
- π-free certificate: all eigenvalues are algebraic integers (no π)
- But NOT all rational: l≥2 sectors produce √5-containing eigenvalues at n_max=3
- Hodge-1 gap: μₙ^{Hodge-1} = n(n+2) while λₙ^{scalar} = n²-1, gap = 2n+1 (Ricci shift)
"""
from __future__ import annotations

import json
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import sympy as sp

from geovac.fock_graph_hodge import FockGraphHodge


def serialize_spectrum(spectrum):
    """Convert sympy eigenvalue list to JSON-serializable format."""
    result = []
    for ev in spectrum:
        ev_expr = sp.sympify(ev)
        entry = {
            "sympy_str": str(ev_expr),
            "float_value": float(ev_expr),
            "is_integer": ev_expr.is_integer,
            "is_rational": ev_expr.is_rational,
            "is_algebraic": True,  # structural: integer matrix → algebraic eigenvalues
            "has_pi": bool(ev_expr.has(sp.pi)),
            "has_sqrt": bool(ev_expr.has(sp.sqrt)),
        }
        result.append(entry)
    return result


def analyze_hodge(n_max: int) -> dict:
    """Build FockGraphHodge for n_max and collect all analysis results."""
    print(f"  Analyzing n_max={n_max}...")
    hodge = FockGraphHodge(n_max)

    # Basic topology
    data = {
        "n_max": n_max,
        "n_nodes": hodge.n_nodes,
        "n_edges": hodge.n_edges,
        "betti_0": hodge.betti_0,
        "betti_1": hodge.betti_1,
        "euler_characteristic": hodge.n_nodes - hodge.n_edges,
    }

    # Verification results
    data["checks"] = {
        "L0_equals_D_minus_A": hodge.verify_L0_equals_D_minus_A(),
        "svd_identity": hodge.verify_svd_identity(),
        "betti_formula": hodge.verify_betti_formula(),
        "pi_free": hodge.verify_pi_free(),
        "algebraic_integer": hodge.verify_algebraic_integer(),
    }

    # Spectra (exact sympy)
    if n_max <= 3:
        L0_spec = hodge.node_laplacian_spectrum()
        L1_spec = hodge.edge_laplacian_spectrum()

        data["L0_spectrum"] = serialize_spectrum(L0_spec)
        data["L1_spectrum"] = serialize_spectrum(L1_spec)

        # Count zeros and irrationals
        data["L0_zero_count"] = sum(1 for ev in L0_spec if float(sp.sympify(ev)) < 1e-10)
        data["L0_irrational_count"] = sum(
            1 for ev in L0_spec
            if sp.sympify(ev).is_rational is False
        )
        data["L1_zero_count"] = sum(1 for ev in L1_spec if float(sp.sympify(ev)) < 1e-10)

        # Sorted float spectra
        data["L0_spectrum_floats_sorted"] = sorted(float(sp.sympify(ev)) for ev in L0_spec)
        data["L1_spectrum_floats_sorted"] = sorted(float(sp.sympify(ev)) for ev in L1_spec)

    else:
        # numpy only for larger n_max
        L0_spec_np = np.sort(hodge.node_laplacian_spectrum_numpy())
        L1_spec_np = np.sort(hodge.edge_laplacian_spectrum_numpy())

        data["L0_spectrum_floats_sorted"] = L0_spec_np.tolist()
        data["L1_spectrum_floats_sorted"] = L1_spec_np.tolist()
        data["L0_zero_count"] = int(np.sum(L0_spec_np < 1e-9))
        data["L1_zero_count"] = int(np.sum(L1_spec_np < 1e-9))

    # Hodge-1 comparison
    hodge1_comparison = hodge.compare_with_hodge1_spectrum()
    data["hodge1_comparison"] = [
        {
            "n": r["n"],
            "lambda_scalar": r["lambda_scalar"],  # n²-1 (graph edge Laplacian)
            "mu_hodge1": r["mu_hodge1"],          # n(n+2) (continuum Hodge-1)
            "gap_ricci": r["gap_ricci"],          # Ricci shift = 2n+1
        }
        for r in hodge1_comparison
    ]

    return data


def main():
    output_path = os.path.join(
        os.path.dirname(__file__), "data", "gn1_fock_graph_hodge.json"
    )

    print("GN-1: Fock Graph Hodge Decomposition Investigation")
    print("=" * 55)

    results = {
        "sprint": "GN-1",
        "description": (
            "Hodge decomposition of the S^3 Fock scalar graph: "
            "signed incidence B, node Laplacian L0=BB^T, edge Laplacian L1=B^TB. "
            "Verifies β₀=n_max (angular sectors), SVD identity, π-free certificate. "
            "Documents the Bochner-Weitzenböck Ricci shift gap=2n+1 between the "
            "graph edge spectrum λₙ=n²-1 and the continuum Hodge-1 spectrum μₙ=n(n+2)."
        ),
        "key_findings": [
            "β₀ = n_max: the Fock graph has n_max connected components (one per l-sector)",
            "T± transitions preserve l, L± preserve (n,l): no transition changes l",
            "β₁ = E - V + n_max loops (closed angular orbits at each level)",
            "SVD identity: nonzero spectra of L₀ and L₁ are identical",
            "π-free: all eigenvalues are algebraic integers, no transcendentals",
            "NOT all rational: l=2 sector at n_max=3 produces √5-containing eigenvalues",
            "Hodge-1 gap = 2n+1 (Ricci shift +2 is continuum curvature, absent from graph)",
        ],
        "n_max_results": {},
    }

    for n_max in [1, 2, 3, 4]:
        data = analyze_hodge(n_max)
        results["n_max_results"][str(n_max)] = data

        print(f"\n  n_max={n_max}: V={data['n_nodes']}, E={data['n_edges']}, "
              f"b0={data['betti_0']}, b1={data['betti_1']}")
        print(f"    Checks: {data['checks']}")
        if "L0_irrational_count" in data:
            print(f"    L0 zeros={data['L0_zero_count']}, irrationals={data['L0_irrational_count']}")

    # Tabulate Hodge-1 comparison at n_max=3
    print("\n  Hodge-1 Comparison (n_max=3 levels):")
    print(f"  {'n':>4} | {'lambda_n':>10} | {'mu_hodge1':>12} | {'gap=2n+1':>10}")
    print("  " + "-" * 44)
    for r in results["n_max_results"]["3"]["hodge1_comparison"]:
        print(f"  {r['n']:>4} | {r['lambda_scalar']:>10} | {r['mu_hodge1']:>12} | {r['gap_ricci']:>10}")

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to {output_path}")
    print("\nAll checks passed!" if all(
        all(v for v in results["n_max_results"][str(n)]["checks"].values())
        for n in [1, 2, 3, 4]
    ) else "WARNING: some checks failed")


if __name__ == "__main__":
    main()
