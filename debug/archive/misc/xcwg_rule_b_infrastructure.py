"""
XCWG Rule B infrastructure survey.
===================================

Survey the Dirac-Rule-B graph (Paper 29 §RH-C) at n_max = 2, 3, 4 for
the forthcoming U(1) Wilson lattice gauge construction. Tabulates:

  §1 Edge set definition (explicit formula)
  §2 Graph properties: V, E, c, beta_1, degree, diameter, path length
  §3 Cross-l edge fraction vs scalar Fock graph
  §4 Plaquette census at n_max = 3 (primitive non-backtracking 4-walks)
  §5 Effective dimension indicators (plaquette / edge ratios)
  §6 Recommendation: is Rule B a substrate for 3D Wilson gauge?

Output: debug/data/xcwg_rule_b_infrastructure.json
        debug/xcwg_rule_b_infrastructure_memo.md (separate writeup)

This is a *survey* script; no production code modified.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from geovac.dirac_matrix_elements import kappa_to_l
from geovac.fock_graph_hodge import FockGraphHodge
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def count_components(A: np.ndarray) -> int:
    """BFS component count."""
    n = A.shape[0]
    if n == 0:
        return 0
    seen = [False] * n
    comps = 0
    for s in range(n):
        if seen[s]:
            continue
        comps += 1
        stack = [s]
        seen[s] = True
        while stack:
            u = stack.pop()
            for v in np.nonzero(A[u])[0]:
                if not seen[int(v)]:
                    seen[int(v)] = True
                    stack.append(int(v))
    return comps


def bfs_distances(A: np.ndarray, src: int) -> np.ndarray:
    """BFS shortest path from src (single-source). Returns -1 for unreachable."""
    n = A.shape[0]
    d = np.full(n, -1, dtype=int)
    d[src] = 0
    queue = [src]
    while queue:
        u = queue.pop(0)
        for v in np.nonzero(A[u])[0]:
            v = int(v)
            if d[v] < 0:
                d[v] = d[u] + 1
                queue.append(v)
    return d


def all_pairs_metrics(A: np.ndarray) -> Tuple[int, float, int]:
    """Return (diameter, avg path length, num reachable pairs) over largest CC.

    For the disconnected scalar Fock graph we compute these on the GIANT
    component; for connected Rule B we compute them globally.
    """
    n = A.shape[0]
    if n == 0:
        return 0, 0.0, 0
    # Compute global diameter / avg-path on connected pairs only
    total_dist = 0
    pair_count = 0
    diameter = 0
    for s in range(n):
        d = bfs_distances(A, s)
        for t in range(n):
            if t == s:
                continue
            if d[t] >= 0:
                total_dist += d[t]
                pair_count += 1
                if d[t] > diameter:
                    diameter = int(d[t])
    avg = (total_dist / pair_count) if pair_count > 0 else 0.0
    return diameter, avg, pair_count // 2  # undirected: halve


def edge_list(A: np.ndarray) -> List[Tuple[int, int]]:
    """Undirected edge list (i < j)."""
    n = A.shape[0]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] != 0:
                edges.append((i, j))
    return edges


def cross_l_fraction(A: np.ndarray, l_per_node: List[int]) -> Tuple[int, int, float]:
    """Return (E_cross, E_total, fraction) where cross-l means la != lb."""
    edges = edge_list(A)
    E_total = len(edges)
    E_cross = sum(1 for (i, j) in edges if l_per_node[i] != l_per_node[j])
    frac = (E_cross / E_total) if E_total > 0 else 0.0
    return E_cross, E_total, frac


# ---------------------------------------------------------------------------
# Plaquette enumeration
# ---------------------------------------------------------------------------

def primitive_4_plaquettes(A: np.ndarray) -> List[Tuple[int, int, int, int]]:
    """Enumerate primitive non-backtracking closed walks of length 4.

    A 4-plaquette is a closed walk v0 -> v1 -> v2 -> v3 -> v0 with:
      - v_{i+1} adjacent to v_i (graph edge)
      - non-backtracking: v_{i+1} != v_{i-1}  (i.e. no immediate reverse)
      - distinct vertices (primitive): all four vertices distinct
      - closes: v3 adjacent to v0

    We return one representative per orbit under cyclic rotation +
    reflection (dihedral group D_4 of order 8) by sorting and using the
    canonical lex-smallest rotation. Duplicates are then deduped.

    Note: this enumerates the actual 4-cycles (C_4 subgraphs).
    """
    n = A.shape[0]
    seen = set()
    cycles = []
    for v0 in range(n):
        neigh0 = [int(j) for j in np.nonzero(A[v0])[0]]
        for v1 in neigh0:
            if v1 == v0:
                continue
            neigh1 = [int(j) for j in np.nonzero(A[v1])[0] if int(j) != v0]
            for v2 in neigh1:
                if v2 == v0 or v2 == v1:
                    continue
                neigh2 = [int(j) for j in np.nonzero(A[v2])[0]
                          if int(j) != v1 and int(j) != v0]
                # require v2-v3 edge AND v3-v0 edge AND v3 != v1
                for v3 in neigh2:
                    if v3 == v0 or v3 == v1 or v3 == v2:
                        continue
                    # Closing edge v3 -> v0
                    if A[v3, v0] == 0:
                        continue
                    # Canonical form: rotations and reflections of (v0,v1,v2,v3)
                    cyc = (v0, v1, v2, v3)
                    rotations = [cyc, cyc[1:] + cyc[:1],
                                 cyc[2:] + cyc[:2], cyc[3:] + cyc[:3]]
                    reflected = [tuple(reversed(r)) for r in rotations]
                    canon = min(rotations + reflected)
                    if canon not in seen:
                        seen.add(canon)
                        cycles.append(canon)
    return cycles


def classify_plaquettes_by_l(
    cycles: List[Tuple[int, int, int, int]], l_per_node: List[int]
) -> Dict[str, int]:
    """Classify 4-plaquettes as within-shell or cross-shell.

    Within-shell: all 4 vertices have same l.
    Cross-shell:  at least 2 distinct l-values among the 4 vertices.
    Further split by the multiset of l-values.
    """
    within = 0
    cross = 0
    l_signature_counts: Counter = Counter()
    for cyc in cycles:
        ls = tuple(sorted(l_per_node[v] for v in cyc))
        l_signature_counts[ls] += 1
        if len(set(ls)) == 1:
            within += 1
        else:
            cross += 1
    return {
        "within_shell": within,
        "cross_shell": cross,
        "total": within + cross,
        "by_l_signature": {str(k): v for k, v in l_signature_counts.items()},
    }


def shortest_cycle_length(A: np.ndarray) -> int:
    """Girth: length of shortest cycle. Returns 0 if forest (no cycles).

    BFS-based girth computation: for each vertex, BFS until we find a
    non-tree edge closing a cycle.
    """
    n = A.shape[0]
    if n == 0:
        return 0
    best = float("inf")
    for src in range(n):
        # BFS with parent tracking
        d = np.full(n, -1, dtype=int)
        parent = np.full(n, -1, dtype=int)
        d[src] = 0
        queue = [src]
        while queue:
            u = queue.pop(0)
            for v in np.nonzero(A[u])[0]:
                v = int(v)
                if d[v] < 0:
                    d[v] = d[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif parent[u] != v:
                    # Non-tree edge: cycle of length d[u] + d[v] + 1
                    cyc_len = d[u] + d[v] + 1
                    if cyc_len < best:
                        best = cyc_len
        if best == 3:
            break  # cannot get smaller (no self-loops)
    return int(best) if best != float("inf") else 0


# ---------------------------------------------------------------------------
# Scalar Fock graph utilities
# ---------------------------------------------------------------------------

def scalar_fock_adjacency_and_labels(n_max: int) -> Tuple[np.ndarray, List[Tuple[int, int, int]]]:
    """Build scalar Fock graph adjacency matrix and (n, l, m) labels."""
    fgh = FockGraphHodge(n_max)
    A_sparse = fgh._lattice.adjacency
    A_dense = (A_sparse.toarray() != 0).astype(int)
    np.fill_diagonal(A_dense, 0)
    labels = list(fgh.states)
    return A_dense, labels


# ---------------------------------------------------------------------------
# Main survey
# ---------------------------------------------------------------------------

def survey_one(n_max: int) -> Dict:
    """Survey both Rule B Dirac graph and scalar Fock graph at given n_max."""
    out: Dict = {"n_max": n_max}

    # --- Rule B Dirac graph
    A_B, dirac_labels, deg_B, desc_B = build_dirac_s3_graph(n_max, "B")
    l_per_node_B = [kappa_to_l(lab.kappa) for lab in dirac_labels]
    V_B = int(A_B.shape[0])
    E_B = int(A_B.sum()) // 2
    c_B = count_components(A_B)
    beta1_B = E_B - V_B + c_B

    Ec_B, Et_B, frac_B = cross_l_fraction(A_B, l_per_node_B)
    diam_B, avg_path_B, pairs_B = all_pairs_metrics(A_B)
    deg_counts_B = Counter(int(d) for d in deg_B)

    girth_B = shortest_cycle_length(A_B)

    # Degree distribution by l-shell
    deg_by_l_B: Dict[int, List[int]] = defaultdict(list)
    for i, lab in enumerate(dirac_labels):
        deg_by_l_B[kappa_to_l(lab.kappa)].append(int(deg_B[i]))

    rule_b_data = {
        "V": V_B,
        "E": E_B,
        "c": c_B,
        "beta_1": beta1_B,
        "diameter": diam_B,
        "avg_path_length": avg_path_B,
        "girth": girth_B,
        "degree_min": int(deg_B.min()) if V_B > 0 else 0,
        "degree_max": int(deg_B.max()) if V_B > 0 else 0,
        "degree_avg": float(deg_B.mean()) if V_B > 0 else 0.0,
        "degree_distribution": dict(sorted(deg_counts_B.items())),
        "cross_l_edges": Ec_B,
        "total_edges": Et_B,
        "cross_l_fraction": frac_B,
        "deg_min_by_l": {l: min(dl) for l, dl in deg_by_l_B.items()},
        "deg_max_by_l": {l: max(dl) for l, dl in deg_by_l_B.items()},
        "deg_avg_by_l": {l: float(np.mean(dl)) for l, dl in deg_by_l_B.items()},
        "labels": [(lab.n_fock, lab.kappa, lab.two_m_j) for lab in dirac_labels],
        "l_per_node": l_per_node_B,
    }
    out["rule_b"] = rule_b_data

    # --- Scalar Fock graph
    A_S, scalar_labels = scalar_fock_adjacency_and_labels(n_max)
    l_per_node_S = [nlm[1] for nlm in scalar_labels]
    V_S = int(A_S.shape[0])
    E_S = int(A_S.sum()) // 2
    c_S = count_components(A_S)
    beta1_S = E_S - V_S + c_S

    Ec_S, Et_S, frac_S = cross_l_fraction(A_S, l_per_node_S)
    diam_S, avg_path_S, pairs_S = all_pairs_metrics(A_S)
    deg_S = A_S.sum(axis=1)
    deg_counts_S = Counter(int(d) for d in deg_S)
    girth_S = shortest_cycle_length(A_S)

    scalar_data = {
        "V": V_S,
        "E": E_S,
        "c": c_S,
        "beta_1": beta1_S,
        "diameter": diam_S,
        "avg_path_length": avg_path_S,
        "girth": girth_S,
        "degree_min": int(deg_S.min()) if V_S > 0 else 0,
        "degree_max": int(deg_S.max()) if V_S > 0 else 0,
        "degree_avg": float(deg_S.mean()) if V_S > 0 else 0.0,
        "degree_distribution": dict(sorted(deg_counts_S.items())),
        "cross_l_edges": Ec_S,
        "total_edges": Et_S,
        "cross_l_fraction": frac_S,
        "labels": [list(nlm) for nlm in scalar_labels],
        "l_per_node": l_per_node_S,
    }
    out["scalar"] = scalar_data

    # --- Plaquettes at this n_max (4-cycles)
    # Compute for both graphs (n_max=2 scalar has no cycles)
    plaq_B = primitive_4_plaquettes(A_B)
    plaq_S = primitive_4_plaquettes(A_S)
    out["rule_b"]["num_4_plaquettes"] = len(plaq_B)
    out["rule_b"]["plaquette_classification"] = classify_plaquettes_by_l(
        plaq_B, l_per_node_B
    )
    out["scalar"]["num_4_plaquettes"] = len(plaq_S)
    out["scalar"]["plaquette_classification"] = classify_plaquettes_by_l(
        plaq_S, l_per_node_S
    )

    # Effective-dimension indicators
    out["rule_b"]["plaquettes_per_edge"] = (
        len(plaq_B) / E_B if E_B > 0 else 0.0
    )
    out["rule_b"]["plaquettes_per_vertex"] = (
        len(plaq_B) / V_B if V_B > 0 else 0.0
    )
    out["scalar"]["plaquettes_per_edge"] = (
        len(plaq_S) / E_S if E_S > 0 else 0.0
    )
    out["scalar"]["plaquettes_per_vertex"] = (
        len(plaq_S) / V_S if V_S > 0 else 0.0
    )

    return out


def main():
    print("XCWG Rule B Infrastructure Survey")
    print("=" * 70)

    results: Dict[str, Dict] = {}
    for n_max in [2, 3, 4]:
        print(f"\nSurveying n_max = {n_max}...")
        r = survey_one(n_max)
        results[f"n_max_{n_max}"] = r

        rb = r["rule_b"]
        sc = r["scalar"]
        print(f"  Rule B: V={rb['V']}, E={rb['E']}, c={rb['c']}, "
              f"beta1={rb['beta_1']}, girth={rb['girth']}")
        print(f"          cross-l: {rb['cross_l_edges']}/{rb['total_edges']} "
              f"= {rb['cross_l_fraction']:.4f}")
        print(f"          diam={rb['diameter']}, avg_path={rb['avg_path_length']:.3f}")
        print(f"          plaquettes (4-cycles): {rb['num_4_plaquettes']}")
        print(f"          within-shell={rb['plaquette_classification']['within_shell']}, "
              f"cross-shell={rb['plaquette_classification']['cross_shell']}")
        print(f"  Scalar: V={sc['V']}, E={sc['E']}, c={sc['c']}, "
              f"beta1={sc['beta_1']}, girth={sc['girth']}")
        print(f"          cross-l: {sc['cross_l_edges']}/{sc['total_edges']} "
              f"= {sc['cross_l_fraction']:.4f}")
        print(f"          plaquettes (4-cycles): {sc['num_4_plaquettes']}")

    # Save to JSON
    outpath = Path("debug/data/xcwg_rule_b_infrastructure.json")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {outpath}")

    # Print key numbers for memo
    print("\n" + "=" * 70)
    print("KEY NUMBERS")
    print("=" * 70)
    for nm in [2, 3, 4]:
        r = results[f"n_max_{nm}"]["rule_b"]
        s = results[f"n_max_{nm}"]["scalar"]
        print(f"\nn_max = {nm}")
        print(f"  Rule B cross-l fraction = {r['cross_l_fraction']:.4f}")
        print(f"  Rule B plaquette/edge   = {r['plaquettes_per_edge']:.4f}")
        if r["num_4_plaquettes"] > 0:
            css = r["plaquette_classification"]
            print(f"  Rule B cross-shell plaq fraction = "
                  f"{css['cross_shell'] / css['total']:.4f}")
        print(f"  Scalar cross-l fraction = {s['cross_l_fraction']:.4f}")
        print(f"  Scalar plaquette/edge   = {s['plaquettes_per_edge']:.4f}")


if __name__ == "__main__":
    main()
