"""Diagnostic: investigate the n_max=4 L=6 enumeration drift.

The running <S> decreased from 0.086 (25k walks) to 0.073 (315k walks).
This is a 15% drift. Are early walks biased toward higher S?

Hypothesis: DFS from v0=0 first explores walks that revisit shell-0
vertices, which have higher local connectivity, contributing to higher
S per walk.
"""
import os, sys, time, json
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from xcwg_wilson_loop_scaling import (
    signed_incidence, adjacency_list, enumerate_primitive_closed_walks,
    vertex_walk_to_edge_indicator, build_K_from_L4_plaquettes,
)

A, _, _, _ = build_dirac_s3_graph(4, 'B')
B_inc, edges, edge_idx = signed_incidence(A)
adj = adjacency_list(A)
plaqs = list(enumerate_primitive_closed_walks(adj, 4))
d_1, K = build_K_from_L4_plaquettes(plaqs, edge_idx, len(edges))
K_pinv = np.linalg.pinv(K)
E_tot = len(edges)

# Compute all 315k L=6 walks and their S values; report running mean at
# decimated points
print("Enumerating all 315k L=6 walks at n_max=4...")
t0 = time.time()
Ss = []
for i, w in enumerate(enumerate_primitive_closed_walks(adj, 6)):
    C = vertex_walk_to_edge_indicator(w, edge_idx, E_tot)
    Ss.append(float(C @ K_pinv @ C))
print(f"  {len(Ss)} walks in {time.time()-t0:.1f}s")

# Running means
Ss = np.array(Ss)
print("\nRunning mean vs walk-index:")
for frac in [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0]:
    idx = int(frac * len(Ss))
    if idx == 0:
        idx = 1
    print(f"  {frac*100:5.1f}% ({idx} walks): mean = {Ss[:idx].mean():.6f}, std = {Ss[:idx].std():.6f}")

# Shuffle and check
rng = np.random.default_rng(42)
Ss_shuffled = Ss.copy()
rng.shuffle(Ss_shuffled)
print("\nShuffled running mean:")
for frac in [0.05, 0.25, 0.50, 1.0]:
    idx = int(frac * len(Ss_shuffled))
    print(f"  {frac*100:5.1f}% ({idx} walks): mean = {Ss_shuffled[:idx].mean():.6f}")

# Now also check n_max=5 L=6: should have similar issue
A5, _, _, _ = build_dirac_s3_graph(5, 'B')
B5, edges5, edge_idx5 = signed_incidence(A5)
adj5 = adjacency_list(A5)
plaqs5 = list(enumerate_primitive_closed_walks(adj5, 4))
d1_5, K5 = build_K_from_L4_plaquettes(plaqs5, edge_idx5, len(edges5))
K5_pinv = np.linalg.pinv(K5)
print("\nEnumerating first 400k L=6 walks at n_max=5...")
Ss5 = []
t0 = time.time()
for i, w in enumerate(enumerate_primitive_closed_walks(adj5, 6)):
    C = vertex_walk_to_edge_indicator(w, edge_idx5, len(edges5))
    Ss5.append(float(C @ K5_pinv @ C))
    if i >= 400_000:
        break
print(f"  {len(Ss5)} walks in {time.time()-t0:.1f}s")
Ss5 = np.array(Ss5)
print("Running mean for n_max=5 L=6:")
for n in [25000, 50000, 100000, 200000, 300000, 400000]:
    if n <= len(Ss5):
        print(f"  {n:6d}: mean = {Ss5[:n].mean():.6f}")

# Save
out = {
    "n_max_4_L6": {
        "running_means": {str(int(frac*len(Ss))): float(Ss[:int(frac*len(Ss))].mean())
                          for frac in [0.05, 0.10, 0.25, 0.50, 0.75, 1.0]},
        "shuffled_running_means": {str(int(frac*len(Ss_shuffled))): float(Ss_shuffled[:int(frac*len(Ss_shuffled))].mean())
                                    for frac in [0.05, 0.25, 0.50, 1.0]},
        "final_mean_full": float(Ss.mean()),
        "drift_first10_vs_last10": {
            "first_10pct": float(Ss[:len(Ss)//10].mean()),
            "last_10pct": float(Ss[9*len(Ss)//10:].mean()),
            "relative_drift": float((Ss[:len(Ss)//10].mean() - Ss[9*len(Ss)//10:].mean()) / Ss.mean()),
        },
    },
    "n_max_5_L6": {
        "running_means": {str(n): float(Ss5[:n].mean()) for n in [25000, 50000, 100000, 200000, 300000, 400000] if n <= len(Ss5)},
        "n_enumerated": len(Ss5),
    },
}
with open(os.path.join(_HERE, 'data', 'xcwg_wls_diagnostic.json'), 'w') as f:
    json.dump(out, f, indent=2)
print("\nWrote debug/data/xcwg_wls_diagnostic.json")
