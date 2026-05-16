"""
Pilot computation for U(1) Wilson lattice gauge on Dirac Rule B graph.

Sprint XCWG (May 2026). Companion to debug/xcwg_u1_wilson_rule_b_design_memo.md.

Constructs explicitly at n_max = 2 and n_max = 3:
  - Rule B adjacency (Paper 29 RH-C)
  - Oriented edges, signed incidence B
  - Node Laplacian L_0 = B B^T
  - Edge Laplacian L_1 = B^T B
  - Numerical Hodge identity verification (nonzero spec(L_0) = nonzero spec(L_1))
  - Plaquette counts at lengths 4, 6, 8 by enumerating primitive non-backtracking
    closed walks
  - Weak-coupling kinetic-term sanity check at n_max = 2

All output to debug/data/xcwg_u1_wilson_rule_b_pilot.json (no production code
modifications, per sprint constraints).
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np

# Make geovac importable without installing the package.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402


# ---------------------------------------------------------------------------
# Signed incidence B and Hodge Laplacians L_0, L_1
# ---------------------------------------------------------------------------

def signed_incidence(adjacency: np.ndarray) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """Build signed incidence matrix B of shape (V, E_undir).

    Convention: for an undirected edge e = {u, v} with u < v, B[u, e] = -1
    and B[v, e] = +1. Then L_0 = B B^T is the standard graph Laplacian (deg
    on diagonal, -1 for adjacent pairs) and L_1 = B^T B is the edge
    Laplacian. We list undirected edges in lex order (u, v) with u < v.
    """
    A = np.asarray(adjacency)
    V = A.shape[0]
    edges: List[Tuple[int, int]] = []
    for u in range(V):
        for v in range(u + 1, V):
            if A[u, v]:
                edges.append((u, v))
    E = len(edges)
    B = np.zeros((V, E), dtype=int)
    for e_idx, (u, v) in enumerate(edges):
        B[u, e_idx] = -1
        B[v, e_idx] = 1
    return B, edges


def hodge_identity_check(L0: np.ndarray, L1: np.ndarray, tol: float = 1e-9) -> dict:
    """Verify nonzero spec(L_0) = nonzero spec(L_1) numerically.

    Returns a dict with the two sorted nonzero spectra, the largest absolute
    discrepancy on the shared portion, the kernel dimensions, and a PASS/FAIL
    flag.
    """
    w0 = np.sort(np.linalg.eigvalsh(L0).real)
    w1 = np.sort(np.linalg.eigvalsh(L1).real)
    nz0 = w0[w0 > tol]
    nz1 = w1[w1 > tol]
    kd0 = int((w0 <= tol).sum())
    kd1 = int((w1 <= tol).sum())
    # Hodge identity: |E| - |V| + c = beta_1 of L_1, |V| - c = rank(L_0).
    if len(nz0) == len(nz1):
        max_diff = float(np.max(np.abs(nz0 - nz1)))
        passed = max_diff < tol * max(1.0, float(np.max(np.abs(nz0))))
    else:
        max_diff = float("inf")
        passed = False
    return {
        "kernel_dim_L0": kd0,
        "kernel_dim_L1": kd1,
        "nonzero_count_L0": int(len(nz0)),
        "nonzero_count_L1": int(len(nz1)),
        "shared_max_diff": max_diff,
        "hodge_identity_passed": bool(passed),
        "nonzero_spec_L0": [round(float(x), 9) for x in nz0],
        "nonzero_spec_L1": [round(float(x), 9) for x in nz1],
    }


# ---------------------------------------------------------------------------
# Plaquette enumeration (primitive closed non-backtracking walks)
# ---------------------------------------------------------------------------

def enumerate_plaquettes(
    adjacency: np.ndarray,
    max_length: int = 8,
) -> Dict[int, int]:
    """Count primitive closed non-backtracking walks by length.

    A closed NB walk of length L on the (oriented) edge set is a sequence
    e_1 ... e_L with head(e_i) = tail(e_{i+1}) and e_{i+1} != reverse(e_i),
    plus head(e_L) = tail(e_1) and e_1 != reverse(e_L). "Primitive" means
    the walk is not a proper power of a shorter walk. Two walks differing
    by cyclic rotation are identified. We also identify a walk with its
    orientation-reversal (since the Wilson action picks up Re tr U_P which
    is invariant under reversal up to conjugation).

    Returns
    -------
    counts : dict {length: count}
        The number of distinct primitive plaquettes of each length up to
        max_length.
    """
    A = np.asarray(adjacency)
    V = A.shape[0]

    # Build oriented edges
    oriented: List[Tuple[int, int]] = []
    idx_of: Dict[Tuple[int, int], int] = {}
    for u in range(V):
        for v in np.nonzero(A[u])[0]:
            oriented.append((u, int(v)))
            idx_of[(u, int(v))] = len(oriented) - 1

    def rev(e_idx: int) -> int:
        u, v = oriented[e_idx]
        return idx_of[(v, u)]

    seen: set = set()
    counts: Counter = Counter()

    def is_primitive(walk: Tuple[int, ...]) -> bool:
        L = len(walk)
        for d in range(1, L):
            if L % d == 0:
                period = walk[:d]
                if all(walk[k * d:(k + 1) * d] == period for k in range(L // d)):
                    return False
        return True

    def canonical(walk: Tuple[int, ...]) -> Tuple[int, ...]:
        # Up to cyclic rotation and orientation reversal, pick the lex-min.
        L = len(walk)
        best = walk
        for k in range(L):
            rot = walk[k:] + walk[:k]
            if rot < best:
                best = rot
        rev_walk = tuple(rev(e) for e in reversed(walk))
        for k in range(L):
            rot = rev_walk[k:] + rev_walk[:k]
            if rot < best:
                best = rot
        return best

    # DFS from each starting oriented edge.
    for start in range(len(oriented)):
        stack: List[List[int]] = [[start]]
        while stack:
            path = stack.pop()
            if len(path) >= max_length:
                continue
            last = path[-1]
            u, v = oriented[last]
            for w in np.nonzero(A[v])[0]:
                w = int(w)
                nxt = idx_of[(v, w)]
                # No backtracking
                if nxt == rev(last):
                    continue
                new_path = path + [nxt]
                L = len(new_path)
                # Check closure
                first = new_path[0]
                fu, fv = oriented[first]
                if w == fu and L >= 3:
                    # closure NB: first edge must not be reverse of last new edge
                    if nxt != rev(first):
                        t = tuple(new_path)
                        if is_primitive(t):
                            can = canonical(t)
                            if can not in seen:
                                seen.add(can)
                                counts[L] += 1
                if L < max_length:
                    stack.append(new_path)
    return dict(sorted(counts.items()))


# ---------------------------------------------------------------------------
# Weak-coupling kinetic-term verification
# ---------------------------------------------------------------------------

def weak_coupling_kinetic_test(
    L1: np.ndarray,
    plaquettes_by_length: Dict[int, int],
) -> dict:
    """Sanity check that the quadratic expansion of S_W has the form a^T L_1 a.

    For U(1) Wilson, S_W = beta * sum_P (1 - cos(theta_P)) where theta_P =
    sum_{e in P} theta_e (oriented sum). To second order in theta_e:
        S_W ~ (beta/2) sum_P theta_P^2 = (beta/2) theta^T (d_1^T d_1) theta,
    where d_1 is the plaquette-to-edge boundary map. For a triangle-free
    graph with all primitive cycles of equal length L_0, the operator
    d_1^T d_1 is *not* identical to L_1 = B^T B; in general L_1 contains
    both the exact part d_0 d_0^* (from "longitudinal" 1-cochains) and the
    co-exact part d_1^T d_1 (from plaquette boundaries), so the strict
    statement is that L_1 governs the propagation of the full edge mode
    space, with d_1^T d_1 controlling the projected gauge-invariant part.

    We report L_1 dimensions and the number of plaquettes of each length;
    a full d_1 construction is in §5 of the memo.
    """
    return {
        "L1_dim": int(L1.shape[0]),
        "plaquettes_by_length": plaquettes_by_length,
        "smallest_plaquette_length": (
            min(plaquettes_by_length) if plaquettes_by_length else None
        ),
    }


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main() -> None:
    results: dict = {"sprint": "XCWG", "graph_family": "Dirac Rule B (Paper 29 RH-C)"}

    for n_max in [2, 3]:
        A, labels, deg, desc = build_dirac_s3_graph(n_max, "B")
        V = int(A.shape[0])
        E = int(A.sum()) // 2
        B_inc, edges = signed_incidence(A)
        L0 = B_inc @ B_inc.T
        L1 = B_inc.T @ B_inc

        hodge = hodge_identity_check(L0, L1)

        # Plaquette enumeration: at n_max=3 this can be expensive; cap length 8.
        max_len = 8 if n_max == 2 else 6
        plaqs = enumerate_plaquettes(A, max_length=max_len)
        weak = weak_coupling_kinetic_test(L1, plaqs)

        # Per-kappa block sizes (Rule B mixes kappa, but the kappa partition
        # is still meaningful for diagnosing structure).
        kappa_counts = Counter(lab.kappa for lab in labels)

        results[f"n_max_{n_max}"] = {
            "description": desc,
            "V": V,
            "E": E,
            "beta_1_from_incidence": int(L1.shape[0] - V + 1),  # c = 1 for Rule B beyond n_max=1
            "degree_sequence": deg.tolist(),
            "degree_distribution": dict(Counter(deg.tolist())),
            "kappa_block_sizes": {int(k): int(c) for k, c in sorted(kappa_counts.items())},
            "hodge": hodge,
            "plaquettes": weak,
            "L1_eigenvalues_sample": [round(float(x), 6)
                                       for x in np.sort(np.linalg.eigvalsh(L1).real)[:20]],
        }

    out_path = os.path.join(_HERE, "data", "xcwg_u1_wilson_rule_b_pilot.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Wrote {out_path}")

    # Brief stdout summary
    for k in ["n_max_2", "n_max_3"]:
        r = results[k]
        print(f"\n{k}:")
        print(f"  V={r['V']}, E={r['E']}, beta_1={r['beta_1_from_incidence']}")
        print(f"  kappa blocks: {r['kappa_block_sizes']}")
        print(f"  Hodge identity passed: {r['hodge']['hodge_identity_passed']}")
        print(f"  kernel dims (L0, L1): "
              f"({r['hodge']['kernel_dim_L0']}, {r['hodge']['kernel_dim_L1']})")
        print(f"  Plaquettes by length: {r['plaquettes']['plaquettes_by_length']}")


if __name__ == "__main__":
    main()
