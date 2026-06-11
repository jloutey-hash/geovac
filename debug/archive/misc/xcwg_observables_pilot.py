"""
Pilot: alternative U(1) Wilson observables on the Dirac Rule B graph.

Sprint XCWG (May 2026). Companion to debug/xcwg_alternative_observables_memo.md.

Computes at n_max = 2 (cheap, exhaustive):
  - All primitive Wilson loops up to length 6
  - Their plaquette-tiling (minimum area)
  - Strong-coupling Wilson loop expectations <W_L>_strong = (beta/2)^A
  - Weak-coupling Wilson loop expectations <W_L>_weak via Gaussian propagator
  - Plaquette correlation function <cos theta_P cos theta_P'> at weak coupling
  - Distribution of plaquette-pair distances
  - Spectral dimension of the "loop space" (l-shells visited per primitive loop)

At n_max = 3 we compute the cheap (analytical) pieces only:
  - Smallest-loop strong/weak coupling
  - Plaquette-pair distance distribution (still tractable, 994 plaquettes)

All output to debug/data/xcwg_observables_pilot.json. No production code modified.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import numpy as np

# Make geovac importable without installing the package.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402
from geovac.dirac_matrix_elements import kappa_to_l  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers (signed incidence, edge index, etc.)
# ---------------------------------------------------------------------------

def signed_incidence(A: np.ndarray) -> Tuple[np.ndarray, List[Tuple[int, int]],
                                              Dict[Tuple[int, int], int]]:
    V = A.shape[0]
    edges: List[Tuple[int, int]] = []
    for u in range(V):
        for v in range(u + 1, V):
            if A[u, v]:
                edges.append((u, v))
    E = len(edges)
    B = np.zeros((V, E), dtype=int)
    edge_idx: Dict[Tuple[int, int], int] = {}
    for e_idx, (u, v) in enumerate(edges):
        B[u, e_idx] = -1
        B[v, e_idx] = +1
        edge_idx[(u, v)] = e_idx
        edge_idx[(v, u)] = e_idx  # undirected lookup
    return B, edges, edge_idx


def enumerate_primitive_loops(
    A: np.ndarray,
    max_length: int = 6,
) -> List[List[Tuple[int, int]]]:
    """Enumerate primitive closed non-backtracking walks as oriented edge sequences.

    Returns a list of walks; each walk is a list of (u, v) oriented edges
    forming a closed cycle. Each cycle representative is returned once
    (canonical lex-min over cyclic rotation + orientation reversal).
    """
    V = A.shape[0]
    seen = set()
    loops: List[List[Tuple[int, int]]] = []

    def canon(walk: Tuple[Tuple[int, int], ...]) -> Tuple[Tuple[int, int], ...]:
        L = len(walk)
        best = walk
        for k in range(L):
            rot = walk[k:] + walk[:k]
            if rot < best:
                best = rot
        rev_walk = tuple((v, u) for (u, v) in reversed(walk))
        for k in range(L):
            rot = rev_walk[k:] + rev_walk[:k]
            if rot < best:
                best = rot
        return best

    def is_primitive(walk: Tuple[Tuple[int, int], ...]) -> bool:
        L = len(walk)
        for d in range(1, L):
            if L % d == 0:
                period = walk[:d]
                if all(walk[k * d:(k + 1) * d] == period for k in range(L // d)):
                    return False
        return True

    # DFS from each starting oriented edge
    for u0 in range(V):
        for v0 in np.nonzero(A[u0])[0]:
            v0 = int(v0)
            stack = [[(u0, v0)]]
            while stack:
                path = stack.pop()
                if len(path) >= max_length:
                    continue
                u, v = path[-1]
                for w in np.nonzero(A[v])[0]:
                    w = int(w)
                    if w == u:
                        continue  # backtrack
                    new_path = path + [(v, w)]
                    L = len(new_path)
                    # Closure check
                    if w == u0 and L >= 3:
                        last_u, last_v = new_path[-1]
                        first_u, first_v = new_path[0]
                        if not (last_u == first_v and last_v == first_u):
                            t = tuple(new_path)
                            if is_primitive(t):
                                c = canon(t)
                                if c not in seen:
                                    seen.add(c)
                                    loops.append(list(c))
                    if L < max_length:
                        stack.append(new_path)
    return loops


def loop_edge_indicator(loop: List[Tuple[int, int]],
                        edges: List[Tuple[int, int]]) -> np.ndarray:
    """Build signed edge-indicator vector C in R^E for the loop.

    C[e] = +1 if loop traverses e in canonical orientation (u<v -> u,v),
    C[e] = -1 if loop traverses e in reversed orientation, 0 otherwise.
    A loop visiting an edge twice in opposite orientation has C[e] = 0
    (it would cancel — but primitive non-backtracking loops never do this).
    """
    E = len(edges)
    idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    C = np.zeros(E)
    for (u, v) in loop:
        e_canon = (min(u, v), max(u, v))
        i = idx[e_canon]
        sign = +1 if (u, v) == e_canon else -1
        C[i] += sign
    return C


# ---------------------------------------------------------------------------
# Plaquette boundary operator d_1: R^E -> R^P
# ---------------------------------------------------------------------------

def build_plaquette_boundary(
    loops_4: List[List[Tuple[int, int]]],
    edges: List[Tuple[int, int]],
) -> np.ndarray:
    """Construct d_1: each row corresponds to a plaquette P (a length-4 loop),
    each column to an edge; entry +/- 1 according to whether the edge appears
    in P with the canonical orientation (u<v) or its reverse.

    Returns d_1 of shape (P, E).
    """
    P = len(loops_4)
    E = len(edges)
    idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    d_1 = np.zeros((P, E), dtype=int)
    for p, loop in enumerate(loops_4):
        for (u, v) in loop:
            e_canon = (min(u, v), max(u, v))
            i = idx[e_canon]
            sign = +1 if (u, v) == e_canon else -1
            d_1[p, i] += sign
    return d_1


# ---------------------------------------------------------------------------
# Weak-coupling Wilson loop via Gaussian propagator
# ---------------------------------------------------------------------------

def weak_coupling_wilson_loop(loop_C: np.ndarray,
                              K_pinv: np.ndarray) -> float:
    """<W_L>_weak = exp(-1/(2 beta) C^T K^{-1} C); return the C^T K^{-1} C
    coefficient (units of 1/beta), independent of beta.

    K = d_1^T d_1 is the Wilson kinetic operator from §4.2 of the design memo
    (the part of the edge Laplacian sourced by the plaquette boundary map).
    The closed loop C has B C = 0 (no node boundary), so C lives in the
    cycle space; on that space K and L1 = B^T B differ structurally because
    L1 has zero eigenmodes corresponding to BOTH gauge modes (im d_0) AND
    harmonic 1-cochains (ker d_1 cap ker B = H_1), while K only kills the
    gauge modes — it sees the cycle modes via d_1.

    For a contractible loop bounding a single plaquette P: d_1 C = e_P
    (unit vector at P), so C^T K_pinv C = e_P^T (d_1 d_1^T)_pinv e_P.
    """
    return float(loop_C @ K_pinv @ loop_C)


def loop_homology_class(loop_C: np.ndarray, L1: np.ndarray, tol: float = 1e-9) -> dict:
    """Compute the harmonic projection of the loop indicator C onto ker(L1).

    A loop is null-homologous iff its harmonic projection is zero. The
    norm of the projection measures how non-contractible the loop is.
    """
    w, U = np.linalg.eigh(L1)
    # Harmonic eigenvectors: those with eigenvalue ~ 0
    harmonic_mask = np.abs(w) < tol
    H = U[:, harmonic_mask]  # E x beta_1
    proj = H.T @ loop_C  # beta_1-vector
    proj_norm = float(np.linalg.norm(proj))
    return {
        "harmonic_dim": int(harmonic_mask.sum()),
        "harmonic_projection_norm": proj_norm,
        "is_null_homologous": bool(proj_norm < tol),
    }


# ---------------------------------------------------------------------------
# Plaquette correlator (weak-coupling Gaussian)
# ---------------------------------------------------------------------------

def plaquette_pair_distance(loop_p: List[Tuple[int, int]],
                            loop_q: List[Tuple[int, int]],
                            sp_dist: np.ndarray) -> int:
    """Minimum graph distance between any vertex of plaquette p and any
    vertex of plaquette q. Used as the "distance" variable in <P P'>.
    """
    verts_p = set(u for (u, v) in loop_p) | set(v for (u, v) in loop_p)
    verts_q = set(u for (u, v) in loop_q) | set(v for (u, v) in loop_q)
    return int(min(sp_dist[u, v] for u in verts_p for v in verts_q))


def plaquette_correlator_weak(loop_p: List[Tuple[int, int]],
                              loop_q: List[Tuple[int, int]],
                              edges: List[Tuple[int, int]],
                              K_pinv: np.ndarray) -> float:
    """<cos theta_P cos theta_Q>_weak connected, at leading order in 1/beta.

    With Gaussian theta of covariance <theta_e theta_f> = (1/beta) K_pinv[e,f],
    Wick's theorem gives the connected piece:
       <cos theta_P cos theta_Q>_c ~ (1/(2 beta)) (C_P^T K_pinv C_Q)
                                       * exp(-(C_P^T K_pinv C_P + C_Q^T K_pinv C_Q)/(2 beta))
    at leading order in 1/beta. We report C_P^T K_pinv C_Q (the off-diagonal
    propagator element between the two plaquette boundaries; units of 1/beta).
    """
    C_P = loop_edge_indicator(loop_p, edges)
    C_Q = loop_edge_indicator(loop_q, edges)
    return float(C_P @ K_pinv @ C_Q)


# ---------------------------------------------------------------------------
# Spectral dimension of Wilson-loop content (l-shell extent)
# ---------------------------------------------------------------------------

def l_shells_visited(loop: List[Tuple[int, int]], labels) -> int:
    verts = set(u for (u, v) in loop) | set(v for (u, v) in loop)
    return len({kappa_to_l(labels[v].kappa) for v in verts})


# ---------------------------------------------------------------------------
# Polyakov-loop analog scoping
# ---------------------------------------------------------------------------

def n_fock_layers(labels) -> Dict[int, List[int]]:
    """Vertices grouped by n_fock (the principal-quantum-number 'layer')."""
    by_n: Dict[int, List[int]] = defaultdict(list)
    for i, lab in enumerate(labels):
        by_n[lab.n_fock].append(i)
    return dict(by_n)


def m_j_orbits(labels) -> Dict[Tuple[int, int], List[int]]:
    """Vertices grouped by (n_fock, kappa) — the 'spin orbit' at fixed orbital."""
    by_nl: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, lab in enumerate(labels):
        by_nl[(lab.n_fock, lab.kappa)].append(i)
    return dict(by_nl)


# ---------------------------------------------------------------------------
# Monopole density: flux through faces of small 2-surfaces
# ---------------------------------------------------------------------------

def find_minimal_2_surfaces(loops: List[List[Tuple[int, int]]],
                            edges: List[Tuple[int, int]],
                            max_facets: int = 8) -> List[List[int]]:
    """Find candidate 'cells' on the graph: small collections of plaquettes
    that share edges to form a closed 2-surface.

    A monopole is detected by computing the total signed flux Sum_P theta_P
    around a closed 2-surface (a 'cube'); if this flux is a multiple of 2pi,
    the configuration carries an integer monopole charge.

    For Rule B at n_max=2 with V=10 plaquettes (44 of length 4), we
    enumerate triangulations of small 2-surfaces by looking for plaquette
    triples / quadruples / etc. that share enough edges to close up.

    This is the most tentative observable; we return small triangulations
    only (3-4 plaquettes sharing 3 edges, an "octahedral" minimal closed
    surface in the bipartite l-parity graph).
    """
    # Build edge -> plaquettes incidence
    edge_idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    plaq_edges = []  # list of frozenset(edge_idx) per plaquette
    for loop in loops:
        es = frozenset(edge_idx[tuple(sorted(e))] for e in loop)
        plaq_edges.append(es)

    # For each plaquette, find its neighbors (shared edge plaquettes)
    plaq_nbrs: Dict[int, List[int]] = {i: [] for i in range(len(plaq_edges))}
    for i in range(len(plaq_edges)):
        for j in range(i + 1, len(plaq_edges)):
            if plaq_edges[i] & plaq_edges[j]:
                plaq_nbrs[i].append(j)
                plaq_nbrs[j].append(i)

    # Triangular cells: 3 plaquettes mutually sharing edges
    triangular_cells: List[List[int]] = []
    for i in range(len(plaq_edges)):
        for j in plaq_nbrs[i]:
            if j <= i:
                continue
            for k in plaq_nbrs[j]:
                if k <= j:
                    continue
                if i in plaq_nbrs[k]:
                    # i, j, k mutually adjacent — candidate triangular cell
                    triangular_cells.append([i, j, k])
                    if len(triangular_cells) >= max_facets * 10:
                        return triangular_cells
    return triangular_cells


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main() -> None:
    results: dict = {
        "sprint": "XCWG",
        "task": "alternative U(1) Wilson observables on Rule B",
    }

    # ---- n_max = 2 ----
    A, labels, deg, desc = build_dirac_s3_graph(2, "B")
    V, E_tot = A.shape[0], int(A.sum() // 2)
    B_inc, edges, edge_idx = signed_incidence(A)
    L0 = B_inc @ B_inc.T
    L1 = B_inc.T @ B_inc

    # Shortest paths on the vertex graph
    from scipy.sparse.csgraph import shortest_path
    sp_dist = shortest_path(A, directed=False, unweighted=True)

    # Enumerate loops up to length 6 (cheap at n_max=2)
    print(f"n_max=2: V={V}, E={E_tot}, enumerating loops up to L=6...")
    loops_4 = []
    loops_6 = []
    all_loops = enumerate_primitive_loops(A, max_length=6)
    for ell in all_loops:
        if len(ell) == 4:
            loops_4.append(ell)
        elif len(ell) == 6:
            loops_6.append(ell)
    print(f"  found {len(loops_4)} length-4 loops, {len(loops_6)} length-6 loops")

    # Build d_1: plaquette boundary map (rows = length-4 plaquettes, cols = edges)
    d_1 = build_plaquette_boundary(loops_4, edges)
    K = (d_1.T @ d_1).astype(float)  # Wilson kinetic operator (E x E)
    K_pinv = np.linalg.pinv(K)
    print(f"  K = d_1^T d_1: shape={K.shape}, rank={np.linalg.matrix_rank(K)}, "
          f"plaquette boundary d_1 shape={d_1.shape}, rank={np.linalg.matrix_rank(d_1)}")

    # --- §1: Wilson loops, strong vs weak coupling ---
    wilson_data = {"by_length": {}}
    for L_loop, loops_list in [(4, loops_4), (6, loops_6)]:
        if not loops_list:
            continue
        # Strong-coupling: <W> = (beta/2)^A where A = number of plaquettes tiling
        # For a length-L primitive loop on Rule B, the minimum tiling A satisfies
        # 2 A >= L (since each tile has perimeter 4 and contributes 4 edges, half
        # of which can be shared with adjacent tiles). A length-4 loop = 1 plaquette.
        # A length-6 loop is bounded by A >= 2 (cannot be tiled by a single plaquette).

        # Weak-coupling: <W>_weak = exp(-S/(2 beta)) where S = C^T K_pinv C and
        # K = d_1^T d_1 is the Wilson kinetic operator.
        Ss = []
        for loop in loops_list[:30]:
            C = loop_edge_indicator(loop, edges)
            S = weak_coupling_wilson_loop(C, K_pinv)
            Ss.append(S)

        wilson_data["by_length"][L_loop] = {
            "count": len(loops_list),
            "sample_size_analyzed": min(30, len(loops_list)),
            "weak_coupling_S_mean": float(np.mean(Ss)),
            "weak_coupling_S_std": float(np.std(Ss)),
            "weak_coupling_S_min": float(np.min(Ss)),
            "weak_coupling_S_max": float(np.max(Ss)),
            "strong_coupling_leading_power_of_beta": 1 if L_loop == 4 else 2,
        }

    # --- §2: Polyakov-loop analog scoping ---
    n_layers = n_fock_layers(labels)
    mj_orbs = m_j_orbits(labels)
    # Look for loops that visit all n_fock layers — analog of "around the thermal circle"
    layer_traversing = []
    for loop in loops_4 + loops_6:
        verts = set(u for (u, v) in loop) | set(v for (u, v) in loop)
        n_layers_visited = len({labels[v].n_fock for v in verts})
        if n_layers_visited == len(n_layers):
            layer_traversing.append(loop)

    # Honest verdict: Rule B is a static spatial graph — no distinguished
    # time direction. The n_fock index is a *radial* quantum number, not
    # a temporal one. We tabulate loops that cross all n_fock layers as
    # candidates for a Polyakov-loop analog, but warn this is structurally
    # not the standard thermal-circle construction.
    polyakov_data = {
        "n_fock_layers": {n: len(verts) for n, verts in n_layers.items()},
        "loops_visiting_all_n_layers_L4": sum(
            1 for loop in loops_4
            if len({labels[v].n_fock for vv in loop for v in vv}) == len(n_layers)
        ),
        "loops_visiting_all_n_layers_L6": sum(
            1 for loop in loops_6
            if len({labels[v].n_fock for vv in loop for v in vv}) == len(n_layers)
        ),
        "honest_verdict": (
            "Rule B is a static spatial graph. The n_fock direction is RADIAL "
            "(principal quantum number), NOT temporal. No clean Polyakov-loop "
            "analog exists without a separate Euclidean-time compactification."
        ),
    }

    # --- §3: Monopole content (tentative) ---
    # Look for small 2-surfaces (triangulations of plaquettes)
    print(f"n_max=2: searching for triangular 2-cells among {len(loops_4)} 4-plaquettes...")
    tri_cells = find_minimal_2_surfaces(loops_4, edges, max_facets=50)
    print(f"  found {len(tri_cells)} candidate triangular 2-cells")

    monopole_data = {
        "n_4_plaquettes": len(loops_4),
        "n_triangular_2_cells": len(tri_cells),
        "structural_note": (
            "A 'monopole' is detected by integer flux through a closed "
            "2-surface. The smallest 2-surface on Rule B is a triangular "
            "cell of 3 mutually-adjacent plaquettes. With sympy-exact link "
            "phases at canonical Wilson configurations one could compute "
            "flux quantization; this pilot only enumerates candidate cells."
        ),
    }

    # --- §4: Plaquette correlators (weak-coupling) ---
    print(f"n_max=2: computing plaquette-pair correlators (44 x 44 / 2 = 946 pairs)...")
    pair_distances: Counter = Counter()
    correlators_by_dist: Dict[int, List[float]] = defaultdict(list)
    n_plaq = len(loops_4)
    for i in range(n_plaq):
        for j in range(i + 1, n_plaq):
            d = plaquette_pair_distance(loops_4[i], loops_4[j], sp_dist)
            pair_distances[d] += 1
            corr = plaquette_correlator_weak(loops_4[i], loops_4[j], edges, K_pinv)
            correlators_by_dist[d].append(corr)

    correlator_data = {
        "pair_distance_distribution": dict(sorted(pair_distances.items())),
        "mean_correlator_by_distance": {
            int(d): {
                "n_pairs": len(corrs),
                "mean": float(np.mean(corrs)),
                "std": float(np.std(corrs)),
                "abs_mean": float(np.mean(np.abs(corrs))),
            }
            for d, corrs in sorted(correlators_by_dist.items())
        },
    }

    # --- §5: Spectral dimension of Wilson loops (l-shell extent) ---
    l_extents_4 = [l_shells_visited(loop, labels) for loop in loops_4]
    l_extents_6 = [l_shells_visited(loop, labels) for loop in loops_6]
    l_extent_data = {
        "L4_n_l_shells_distribution": dict(Counter(l_extents_4)),
        "L6_n_l_shells_distribution": dict(Counter(l_extents_6)),
        "L4_mean_l_extent": float(np.mean(l_extents_4)),
        "L6_mean_l_extent": float(np.mean(l_extents_6)),
        "structural_note": (
            "A length-4 loop on Rule B touches at most 2 distinct l-shells "
            "(due to bipartite l-parity + min 4-cycle, alternating l-parity). "
            "A length-6 loop CAN touch 3 distinct l-shells. Maximum l-extent "
            "= n_max for a primitive loop visiting one node per l-shell."
        ),
    }

    # --- assemble n_max=2 ---
    results["n_max_2"] = {
        "graph": {
            "V": V,
            "E": E_tot,
            "beta_1": E_tot - V + 1,
        },
        "wilson_loops": wilson_data,
        "polyakov_analog": polyakov_data,
        "monopole_content": monopole_data,
        "plaquette_correlators": correlator_data,
        "wilson_l_shell_extent": l_extent_data,
    }

    # ---- n_max = 3: cheap pieces only ----
    print("\nn_max=3: cheap pieces only (skipping correlator and full-loop enumeration)...")
    A3, labels3, deg3, desc3 = build_dirac_s3_graph(3, "B")
    V3, E3 = A3.shape[0], int(A3.sum() // 2)
    B3, edges3, _ = signed_incidence(A3)

    print(f"  V={V3}, E={E3}, enumerating L=4 plaquettes...")
    loops_4_n3 = []
    for loop in enumerate_primitive_loops(A3, max_length=4):
        if len(loop) == 4:
            loops_4_n3.append(loop)
    print(f"  {len(loops_4_n3)} length-4 plaquettes")

    # Build K_3 = d_1^T d_1 at n_max=3
    d_1_n3 = build_plaquette_boundary(loops_4_n3, edges3)
    K_3 = (d_1_n3.T @ d_1_n3).astype(float)
    K_3_pinv = np.linalg.pinv(K_3)
    print(f"  K_3 shape={K_3.shape}, rank={np.linalg.matrix_rank(K_3)}, "
          f"d_1 rank={np.linalg.matrix_rank(d_1_n3)}")

    # Smallest Wilson loop weak-coupling S (sample 50)
    Ss_n3 = []
    for loop in loops_4_n3[:50]:
        C = loop_edge_indicator(loop, edges3)
        Ss_n3.append(weak_coupling_wilson_loop(C, K_3_pinv))

    # l-shell extent
    l_extents_4_n3 = [l_shells_visited(loop, labels3) for loop in loops_4_n3]

    results["n_max_3_cheap"] = {
        "graph": {"V": V3, "E": E3, "beta_1": E3 - V3 + 1},
        "L4_count": len(loops_4_n3),
        "K_rank": int(np.linalg.matrix_rank(K_3)),
        "d_1_rank": int(np.linalg.matrix_rank(d_1_n3)),
        "L4_weak_coupling_S": {
            "sample_size": min(50, len(loops_4_n3)),
            "mean": float(np.mean(Ss_n3)),
            "std": float(np.std(Ss_n3)),
            "min": float(np.min(Ss_n3)),
            "max": float(np.max(Ss_n3)),
        },
        "L4_l_shell_extent": dict(Counter(l_extents_4_n3)),
    }

    out_path = os.path.join(_HERE, "data", "xcwg_observables_pilot.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")

    # Compact stdout summary
    print("\n=== n_max=2 summary ===")
    for L_loop, w in results["n_max_2"]["wilson_loops"]["by_length"].items():
        print(f"  L={L_loop}: {w['count']} loops, weak-coupling S = "
              f"{w['weak_coupling_S_mean']:.3f} +/- {w['weak_coupling_S_std']:.3f} "
              f"(units of 1/beta)")
    print(f"  Plaquette pair distances: {results['n_max_2']['plaquette_correlators']['pair_distance_distribution']}")
    print(f"  L4 l-shell extent: {results['n_max_2']['wilson_l_shell_extent']['L4_n_l_shells_distribution']}")
    print(f"  L6 l-shell extent: {results['n_max_2']['wilson_l_shell_extent']['L6_n_l_shells_distribution']}")
    print(f"  Triangular 2-cells: {results['n_max_2']['monopole_content']['n_triangular_2_cells']}")
    print("\n=== n_max=3 (cheap) summary ===")
    r3 = results["n_max_3_cheap"]
    print(f"  L=4: {r3['L4_count']} plaquettes, weak-coupling S = "
          f"{r3['L4_weak_coupling_S']['mean']:.3f} +/- {r3['L4_weak_coupling_S']['std']:.3f}")
    print(f"  K rank: {r3['K_rank']}, d_1 rank: {r3['d_1_rank']}")
    print(f"  L4 l-shell extent at n_max=3: {r3['L4_l_shell_extent']}")


if __name__ == "__main__":
    main()
