"""
Sprint ST-SU3 plaquette tabulation and Wilson-loop Monte Carlo.

Produces:
  - Plaquette counts at length L = 4, 6, 8 for N_max = 2, 3, 4.
  - First-Betti numbers β_1 = E - V + c.
  - Wilson loop expectation values at N_max = 2 (smallest non-trivial)
    via Metropolis Monte Carlo at beta in {0.5, 1.0, 2.0, 5.0}.

Output: debug/data/st_su3_plaquettes.json
        debug/data/st_su3_wilson_loops.json
"""

import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp_sparse

from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.su3_wilson_s5 import (
    OrientedEdge,
    bargmann_adjacency_dense,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    monte_carlo_wilson_expectation,
)


def graph_stats(N_max: int):
    """Compute V, E, c, beta_1 for the Bargmann graph at N_max."""
    g = build_bargmann_graph(N_max)
    A = bargmann_adjacency_dense(N_max)
    V = A.shape[0]
    E = int(A.sum() // 2)
    # Connected components via BFS
    from collections import deque
    visited = [False] * V
    c = 0
    for s in range(V):
        if visited[s]:
            continue
        c += 1
        q = deque([s])
        visited[s] = True
        while q:
            u = q.popleft()
            for v in np.nonzero(A[u])[0]:
                if not visited[v]:
                    visited[v] = True
                    q.append(int(v))
    beta_1 = E - V + c
    return {"V": V, "E": E, "c": c, "beta_1": beta_1}


def plaquette_count_by_length(
    N_max: int, max_length: int = 8
):
    """Enumerate plaquettes up to max_length, group by length."""
    A = bargmann_adjacency_dense(N_max)
    plaqs = enumerate_plaquettes(A, max_length=max_length, both_orientations=False)
    counts = {}
    for L in range(2, max_length + 1):
        counts[L] = sum(1 for P in plaqs if len(P) == L)
    return counts, len(plaqs)


def plaquette_table(N_max_list, max_length: int = 8):
    """Build a table of {N_max: stats + plaquette counts by length}."""
    table = {}
    for N_max in N_max_list:
        t0 = time.time()
        stats = graph_stats(N_max)
        counts, total = plaquette_count_by_length(N_max, max_length)
        elapsed = time.time() - t0
        table[N_max] = {
            **stats,
            "plaquette_counts": {f"L={L}": v for L, v in counts.items() if v > 0},
            "total_plaquettes_up_to_L_{}".format(max_length): total,
            "elapsed_sec": round(elapsed, 2),
        }
    return table


def wilson_loop_run(
    N_max: int,
    beta_list,
    n_samples: int = 800,
    n_thermalize: int = 200,
    seed: int = 12,
    max_length: int = 4,
):
    """Run Metropolis Monte Carlo for the Wilson loop on a primitive
    plaquette, at several beta values."""
    A = bargmann_adjacency_dense(N_max)
    oriented, _ = enumerate_oriented_edges(A)
    plaqs = enumerate_plaquettes(A, max_length=max_length, both_orientations=False)
    if len(plaqs) == 0:
        return {"N_max": N_max, "skipped": "no plaquettes"}

    forward_edges = [
        OrientedEdge(e.source, e.target)
        for e in oriented if e.source < e.target
    ]
    test_loop = plaqs[0]  # first primitive plaquette

    out = {
        "N_max": N_max,
        "n_links": len(forward_edges),
        "n_plaquettes": len(plaqs),
        "loop_length": len(test_loop),
        "n_samples": n_samples,
        "n_thermalize": n_thermalize,
        "seed": seed,
        "update_mode": "local",
        "beta_results": [],
    }
    for beta in beta_list:
        t0 = time.time()
        mean, stderr = monte_carlo_wilson_expectation(
            test_loop, plaqs, forward_edges,
            beta=beta, n_samples=n_samples, n_thermalize=n_thermalize, seed=seed,
            update="local",
        )
        elapsed = time.time() - t0
        out["beta_results"].append({
            "beta": beta,
            "W_mean": mean,
            "W_stderr": stderr,
            "elapsed_sec": round(elapsed, 1),
        })
        print(
            f"  beta={beta:5.2f}  <W>={mean:+.4f} +- {stderr:.4f}  ({elapsed:.1f}s)"
        )
    return out


if __name__ == "__main__":
    print("=" * 60)
    print("Sprint ST-SU3: Plaquette tabulation")
    print("=" * 60)
    table = plaquette_table([1, 2, 3, 4], max_length=8)
    for N_max, info in table.items():
        print(f"\nN_max = {N_max}:")
        for k, v in info.items():
            print(f"  {k} = {v}")

    out_path = Path("debug/data/st_su3_plaquettes.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({"sprint": "ST-SU3", "table": table}, f, indent=2)
    print(f"\nSaved: {out_path}")

    print("\n" + "=" * 60)
    print("Wilson-loop Monte Carlo at N_max = 2 (smallest non-trivial)")
    print("Local updates with auto-tuned scale ~ 1/sqrt(beta)")
    print("=" * 60)
    wl_n2 = wilson_loop_run(
        N_max=2,
        beta_list=[0.5, 1.0, 2.0, 5.0, 10.0],
        n_samples=600,
        n_thermalize=200,
        seed=12,
    )

    print("\n" + "=" * 60)
    print("Wilson-loop Monte Carlo at N_max = 3")
    print("=" * 60)
    wl_n3 = wilson_loop_run(
        N_max=3,
        beta_list=[0.5, 1.0, 2.0, 5.0],
        n_samples=400,
        n_thermalize=200,
        seed=12,
    )

    out_path2 = Path("debug/data/st_su3_wilson_loops.json")
    with open(out_path2, "w") as f:
        json.dump(
            {"sprint": "ST-SU3", "N_max=2": wl_n2, "N_max=3": wl_n3},
            f, indent=2,
        )
    print(f"\nSaved: {out_path2}")
