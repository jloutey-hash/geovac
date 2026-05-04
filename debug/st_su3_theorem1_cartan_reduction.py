"""
Sprint ST-SU3 Theorem 1: SU(3) Wilson Cartan-torus reduction.

Verifies at machine precision on the Bargmann-Segal graph at N_max = 2, 3
that the SU(3) Wilson action restricted to the maximal Cartan torus
T = U(1) x U(1) reduces to the U(1) x U(1) Wilson action with character
(1/3)(cos a_P + cos b_P + cos(a_P + b_P)).

Output: debug/data/st_su3_theorem1_cartan.json
"""

import json
from pathlib import Path

import numpy as np

from geovac.su3_wilson_s5 import (
    bargmann_adjacency_dense,
    cartan_links_from_phases,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    u1xu1_action_from_su3,
    wilson_action,
)


def cartan_reduction_test(N_max: int, seed: int = 7, n_trials: int = 5):
    """Random Cartan phases on every forward edge; verify reduction."""
    A = bargmann_adjacency_dense(N_max)
    oriented, _ = enumerate_oriented_edges(A)
    plaqs = enumerate_plaquettes(A, max_length=8, both_orientations=False)
    forward = [(e.source, e.target) for e in oriented if e.source < e.target]

    rng = np.random.default_rng(seed)
    diffs = []
    for trial in range(n_trials):
        phases_a = {k: float(rng.uniform(-np.pi, np.pi)) for k in forward}
        phases_b = {k: float(rng.uniform(-np.pi, np.pi)) for k in forward}
        links = cartan_links_from_phases(phases_a, phases_b)
        beta = float(rng.uniform(0.1, 5.0))
        S_su3 = wilson_action(plaqs, links, beta)
        S_u1xu1 = u1xu1_action_from_su3(plaqs, phases_a, phases_b, beta)
        diff = abs(S_su3 - S_u1xu1)
        diffs.append({
            "trial": trial,
            "beta": beta,
            "S_SU3": S_su3,
            "S_U1xU1": S_u1xu1,
            "diff": diff,
        })

    return {
        "N_max": N_max,
        "n_plaquettes": len(plaqs),
        "n_forward_edges": len(forward),
        "n_trials": n_trials,
        "max_diff": max(d["diff"] for d in diffs) if diffs else 0.0,
        "trials": diffs,
        "verdict": (
            "Cartan reduction holds at machine precision"
            if (diffs and max(d["diff"] for d in diffs) < 1e-10)
            else "POSSIBLE DEVIATION (diff exceeds 1e-10)"
        ),
    }


if __name__ == "__main__":
    out = {
        "sprint": "ST-SU3",
        "theorem": "Theorem 1: Cartan torus reduction SU(3) -> U(1) x U(1)",
        "Bargmann_S5_graph": True,
        "results": {},
    }
    for N_max in [2, 3]:
        print(f"--- Bargmann N_max = {N_max} ---")
        r = cartan_reduction_test(N_max)
        out["results"][f"N_max={N_max}"] = r
        print(f"  plaquettes: {r['n_plaquettes']}")
        print(f"  forward edges: {r['n_forward_edges']}")
        print(f"  max |S_SU3 - S_U1xU1|: {r['max_diff']:.2e}")
        print(f"  verdict: {r['verdict']}")

    out["overall_verdict"] = (
        "POSITIVE: Cartan torus reduction is exact at machine precision"
        if all(r["max_diff"] < 1e-10 for r in out["results"].values())
        else "MIXED: see per-N_max diff"
    )
    print()
    print(f"OVERALL: {out['overall_verdict']}")

    out_path = Path("debug/data/st_su3_theorem1_cartan.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {out_path}")
