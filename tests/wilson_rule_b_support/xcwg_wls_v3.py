"""XCWG Wilson-loop scaling v3: bias-corrected + simple-cycle option.

v3 fixes from v1:
  (a) Randomized DFS start-vertex order (v2) to avoid DFS-from-v0=0 bias.
  (b) Optional simple-cycle filter: a "simple cycle" is a closed walk with
      no vertex repetition (every vertex appears at most once before closure).
      Composite walks (figure-eights, theta-graphs) are then EXCLUDED.

Why simple cycles matter:
  In lattice gauge theory, "Wilson loops" are primitive simple cycles. A
  figure-eight (two simple cycles joined at a vertex) has Wilson expectation
  <W> that factorizes: <W_C1 W_C2> ≈ <W_C1><W_C2>. The action S of a
  figure-eight is essentially the sum of S of its component simple cycles.
  Including figure-eights in the average mixes higher-L Wilson statistics
  with lower-L statistics; for the area-law/perimeter-law diagnostic on
  primitive Wilson loops, ONLY simple cycles are the canonical objects.

Sprint outputs both the unrestricted (all primitive non-backtracking) and
the simple-cycle-only means for every (n_max, L).
"""
from __future__ import annotations
import os, sys, time, json
from typing import Dict, List, Iterable, Tuple
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from xcwg_wilson_loop_scaling import (
    signed_incidence, adjacency_list, _canonical_walk, _is_primitive,
    vertex_walk_to_edge_indicator, build_K_from_L4_plaquettes,
    fit_scaling_exponent,
)

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)


def enumerate_primitive_closed_walks_streaming(
    adj: List[List[int]],
    L: int,
    vertex_order: List[int] = None,
    simple_only: bool = False,
) -> Iterable[Tuple[int, ...]]:
    """Generator of canonical primitive closed non-backtracking walks.

    simple_only: if True, restrict to simple cycles (no repeated vertex
        in the walk before closure).
    """
    V = len(adj)
    seen = set()
    order = vertex_order if vertex_order is not None else list(range(V))
    for v0 in order:
        stack = [(v0,)]
        while stack:
            path = stack.pop()
            k = len(path)
            cur = path[-1]
            prev = path[-2] if k >= 2 else -1
            if k == L:
                if v0 in adj[cur]:
                    if k < 2:
                        continue
                    if v0 == path[-2]:
                        continue
                    if path[-1] == path[1]:
                        continue
                    closed = tuple(path)
                    if simple_only and len(set(closed)) != L:
                        continue
                    if _is_primitive(closed):
                        canon = _canonical_walk(closed + (v0,))
                        if canon not in seen:
                            seen.add(canon)
                            yield canon
                continue
            for w in adj[cur]:
                if w == prev:
                    continue
                if simple_only and w in path:
                    continue  # avoid revisiting vertices
                stack.append(path + (w,))


def measure_S(adj, edges, edge_idx, K_pinv, L, V, cap, simple_only=False, seed=42):
    rng = np.random.default_rng(seed)
    order = rng.permutation(V).tolist()
    n_edges = len(edges)
    t0 = time.time()
    sum_S = 0.0; sum_S2 = 0.0; count = 0
    S_min = +np.inf; S_max = -np.inf
    capped = False
    running = []
    for walk in enumerate_primitive_closed_walks_streaming(
        adj, L, vertex_order=order, simple_only=simple_only):
        C = vertex_walk_to_edge_indicator(walk, edge_idx, n_edges)
        S = float(C @ K_pinv @ C)
        sum_S += S; sum_S2 += S * S
        if S < S_min: S_min = S
        if S > S_max: S_max = S
        count += 1
        if count % max(1, cap // 6) == 0:
            running.append({"n": count, "mean": sum_S / count})
        if count >= cap:
            capped = True
            break
    dt = time.time() - t0
    mean = sum_S / count if count else float("nan")
    var = max(0.0, sum_S2 / count - mean * mean) if count else 0.0
    return {
        "L": L, "count": count, "mean": float(mean), "std": float(np.sqrt(var)),
        "min": float(S_min) if count else float("nan"),
        "max": float(S_max) if count else float("nan"),
        "walltime_sec": dt, "capped": capped, "running_means": running,
        "simple_only": simple_only,
    }


def run_n_max(n_max, target_lengths, caps_full, caps_simple):
    print(f"\n=== n_max = {n_max} ===")
    A, _, _, _ = build_dirac_s3_graph(n_max, 'B')
    V = A.shape[0]; E_tot = int(A.sum() // 2)
    print(f"  V={V}, E={E_tot}, beta_1={E_tot-V+1}")
    B_inc, edges, edge_idx = signed_incidence(A)
    adj = adjacency_list(A)
    plaqs = list(enumerate_primitive_closed_walks_streaming(adj, 4))
    d_1, K = build_K_from_L4_plaquettes(plaqs, edge_idx, E_tot)
    K_pinv = np.linalg.pinv(K)
    print(f"  {len(plaqs)} L=4 plaquettes; K rank = {np.linalg.matrix_rank(K)}")

    by_L_full = {}
    by_L_simple = {}
    for L in target_lengths:
        print(f"  L={L} FULL (cap={caps_full[L]})...")
        r_full = measure_S(adj, edges, edge_idx, K_pinv, L, V, caps_full[L], simple_only=False)
        print(f"    count={r_full['count']}, <S>={r_full['mean']:.5f}+/- {r_full['std']:.5f}"
              f"{', CAPPED' if r_full['capped'] else ''}")
        by_L_full[L] = r_full
        print(f"  L={L} SIMPLE-only (cap={caps_simple[L]})...")
        r_simple = measure_S(adj, edges, edge_idx, K_pinv, L, V, caps_simple[L], simple_only=True)
        print(f"    count={r_simple['count']}, <S>={r_simple['mean']:.5f}+/- {r_simple['std']:.5f}"
              f"{', CAPPED' if r_simple['capped'] else ''}")
        by_L_simple[L] = r_simple

    L_used = [L for L in target_lengths if by_L_full[L]["count"] > 0]
    fit_full = fit_scaling_exponent(L_used, [by_L_full[L]["mean"] for L in L_used])
    L_used_s = [L for L in target_lengths if by_L_simple[L]["count"] > 0]
    fit_simple = fit_scaling_exponent(L_used_s, [by_L_simple[L]["mean"] for L in L_used_s])
    print(f"  FULL fit alpha={fit_full['alpha']:.4f}, R^2={fit_full['r_squared']:.4f}, "
          f"verdict={fit_full['verdict']}")
    print(f"  SIMPLE fit alpha={fit_simple['alpha']:.4f}, R^2={fit_simple['r_squared']:.4f}, "
          f"verdict={fit_simple['verdict']}")
    return {
        "n_max": n_max,
        "graph": {"V": V, "E": E_tot, "beta_1": E_tot-V+1},
        "K_rank": int(np.linalg.matrix_rank(K)),
        "plaquettes": len(plaqs),
        "wilson_action_full": by_L_full,
        "wilson_action_simple": by_L_simple,
        "scaling_fit_full": fit_full,
        "scaling_fit_simple": fit_simple,
    }


def main():
    out = {"sprint": "XCWG Wilson-loop scaling v3 (bias-corrected + simple-cycle)",
           "date": "2026-05-15",
           "notes": ("Randomized DFS start-vertex order + simple-cycle filter."
                     " Both unrestricted (FULL) and simple-cycle-only (SIMPLE) means"
                     " reported for each n_max, L.")}

    # Strategy: full enumeration where possible, randomized cap where not.
    # Simple-cycle filter cuts walk counts substantially, so its cap can be lower.
    out["n_max_2"] = run_n_max(
        2, [4, 6, 8],
        caps_full={4: 10_000_000, 6: 10_000_000, 8: 10_000_000},
        caps_simple={4: 10_000_000, 6: 10_000_000, 8: 10_000_000},
    )
    out["n_max_3"] = run_n_max(
        3, [4, 6, 8],
        caps_full={4: 10_000_000, 6: 10_000_000, 8: 1_000_000},
        caps_simple={4: 10_000_000, 6: 10_000_000, 8: 1_000_000},
    )
    out["n_max_4"] = run_n_max(
        4, [4, 6, 8],
        caps_full={4: 10_000_000, 6: 10_000_000, 8: 500_000},
        caps_simple={4: 10_000_000, 6: 1_000_000, 8: 500_000},
    )
    out["n_max_5"] = run_n_max(
        5, [4, 6],
        caps_full={4: 500_000, 6: 500_000},
        caps_simple={4: 500_000, 6: 500_000},
    )
    out["n_max_5"]["L8_skipped"] = "Compute-intractable at L=8 on n_max=5"

    # Cross summary
    alphas = {}
    for k in ["n_max_2", "n_max_3", "n_max_4", "n_max_5"]:
        alphas[k] = {
            "alpha_full": out[k]["scaling_fit_full"]["alpha"],
            "alpha_simple": out[k]["scaling_fit_simple"]["alpha"],
            "r2_full": out[k]["scaling_fit_full"]["r_squared"],
            "r2_simple": out[k]["scaling_fit_simple"]["r_squared"],
            "verdict_full": out[k]["scaling_fit_full"]["verdict"],
            "verdict_simple": out[k]["scaling_fit_simple"]["verdict"],
        }
    out["cross_n_max_alpha"] = alphas

    a_full = [v["alpha_full"] for v in alphas.values()]
    a_simple = [v["alpha_simple"] for v in alphas.values()]
    out["witness_assessment"] = {
        "alpha_full_trajectory": a_full,
        "alpha_simple_trajectory": a_simple,
        "max_alpha_full": float(max(a_full)),
        "max_alpha_simple": float(max(a_simple)),
        "comment_full": (
            "FULL (all primitive closed walks): "
            + ("area-law (alpha~2)" if all(abs(a-2)<0.4 for a in a_full)
               else "perimeter-law (alpha~1)" if all(abs(a-1)<0.4 for a in a_full)
               else "INTERMEDIATE")
        ),
        "comment_simple": (
            "SIMPLE (vertex-disjoint cycles only): "
            + ("area-law (alpha~2)" if all(abs(a-2)<0.4 for a in a_simple)
               else "perimeter-law (alpha~1)" if all(abs(a-1)<0.4 for a in a_simple)
               else "INTERMEDIATE")
        ),
    }

    out_path = os.path.join(_HERE, 'data', 'xcwg_wilson_loop_scaling.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # Print table
    print("\n" + "="*78)
    print(f"{'n_max':<7}{'alpha_FULL':<14}{'verdict':<14}{'alpha_SIMPLE':<14}{'verdict':<14}")
    print("="*78)
    for k, v in alphas.items():
        nm = k.replace("n_max_", "")
        print(f"{nm:<7}{v['alpha_full']:<14.4f}{v['verdict_full']:<14}"
              f"{v['alpha_simple']:<14.4f}{v['verdict_simple']:<14}")
    print()
    print(f"FULL-mode comment: {out['witness_assessment']['comment_full']}")
    print(f"SIMPLE-mode comment: {out['witness_assessment']['comment_simple']}")


if __name__ == "__main__":
    main()
