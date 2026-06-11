"""
Sprint 5 Track S5: Build the S^5 Bargmann-Segal graph at N_max = 5
and verify the Paper 24 certificate (56 nodes, 165 edges, pi-free).

This script reuses the production module geovac/nuclear/bargmann_graph.py
and adds light instrumentation + JSON dumps suitable for the Paper 25
framework-observation extension sprint.

Outputs:
  debug/data/s5_bargmann_graph.json
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

from geovac.nuclear.bargmann_graph import (
    build_bargmann_graph,
    verify_pi_free,
    enumerate_nodes,
    total_nodes,
    shell_degeneracy,
)


def _fraction_to_str(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def main() -> None:
    out_path = Path(__file__).parent / "data" / "s5_bargmann_graph.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    N_max = 5
    cert = verify_pi_free(N_max)

    # Paper 24 target counts
    assert cert["n_nodes"] == 56, f"expected 56 nodes, got {cert['n_nodes']}"
    assert cert["n_edges"] == 165, f"expected 165 edges, got {cert['n_edges']}"
    assert cert["pi_free"] is True

    g = build_bargmann_graph(N_max)
    nodes = g.nodes
    node_index = g.index

    # Classify edges by (N, l, l') transition
    edge_classification: Dict[str, int] = {}
    for (i, j), w in g.adjacency.items():
        (Ni, li, mi) = nodes[i]
        (Nj, lj, mj) = nodes[j]
        # Ensure i < j is the lower-N endpoint (by build convention)
        if Ni > Nj:
            Ni, Nj = Nj, Ni
            li, lj = lj, li
            mi, mj = mj, mi
        dl = lj - li
        key = f"N:{Ni}->{Ni+1}, l:{li}->{li+dl}"
        edge_classification[key] = edge_classification.get(key, 0) + 1

    # Per-N node and edge counts
    nodes_by_N = [0] * (N_max + 1)
    for (N, l, m) in nodes:
        nodes_by_N[N] += 1

    edges_by_N_pair = {}
    for (i, j), w in g.adjacency.items():
        Ni = nodes[i][0]
        Nj = nodes[j][0]
        Nlow, Nhigh = min(Ni, Nj), max(Ni, Nj)
        key = f"{Nlow}->{Nhigh}"
        edges_by_N_pair[key] = edges_by_N_pair.get(key, 0) + 1

    # Unique edge weight values
    unique_weights = sorted(set(g.adjacency.values()), key=lambda f: float(f))

    # Connectivity: is the graph connected?
    # Use union-find on the edges
    parent = list(range(len(nodes)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for (i, j) in g.adjacency.keys():
        union(i, j)
    roots = {find(x) for x in range(len(nodes))}
    n_components = len(roots)

    # Euler / Betti
    # beta_1 = E - V + c
    E = len(g.adjacency)
    V = len(nodes)
    beta_1 = E - V + n_components

    results = {
        "sprint": "S5 Track S5 -- S^5 Bargmann-Segal graph at N_max=5",
        "N_max": N_max,
        "n_nodes": cert["n_nodes"],
        "n_edges": cert["n_edges"],
        "pi_free": cert["pi_free"],
        "nodes_by_N": nodes_by_N,
        "total_nodes_formula": total_nodes(N_max),
        "shell_degeneracies": [shell_degeneracy(N) for N in range(N_max + 1)],
        "edges_by_N_pair": edges_by_N_pair,
        "edge_classification_by_l": edge_classification,
        "unique_weight_values": [_fraction_to_str(w) for w in unique_weights],
        "n_connected_components": n_components,
        "beta_1_first_betti_number": beta_1,
        "euler_characteristic": V - E,
        "paper_24_certified": cert["n_nodes"] == 56 and cert["n_edges"] == 165,
    }

    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"S^5 Bargmann-Segal graph (N_max={N_max}):")
    print(f"  Nodes: {V}  (Paper 24 target: 56 -- {'OK' if V == 56 else 'FAIL'})")
    print(f"  Edges: {E}  (Paper 24 target: 165 -- {'OK' if E == 165 else 'FAIL'})")
    print(f"  pi-free: {cert['pi_free']}")
    print(f"  Connected components: {n_components}")
    print(f"  beta_1 = E - V + c = {E} - {V} + {n_components} = {beta_1}")
    print(f"  Euler characteristic V - E = {V - E}")
    print(f"  Shell degeneracies: {results['shell_degeneracies']}")
    print(f"  Nodes by N: {nodes_by_N}")
    print(f"  Edges by N-pair: {edges_by_N_pair}")
    print(f"  Edge classes by (N, l, l'): {edge_classification}")
    print(f"  Unique weight values: {len(unique_weights)}")
    for w in unique_weights[:20]:
        print(f"    {w} = {float(w):.6f}")
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
