"""Diagnostic: which L=6 walks does new enumeration find that pilot doesn't?

Pilot finds 144 L=6 loops at n_max=2. New code finds 568. Factor ~4.
Are the extras figure-eight loops (vertex repetition)?
"""
import os
import sys
from collections import Counter

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import numpy as np
from geovac.ihara_zeta_dirac import build_dirac_s3_graph

from xcwg_wilson_loop_scaling import (
    signed_incidence, adjacency_list, enumerate_primitive_closed_walks,
)

A, labels, deg, desc = build_dirac_s3_graph(2, "B")
adj = adjacency_list(A)

# Count L=6 walks via new method, classify by vertex-repetition
new_walks = list(enumerate_primitive_closed_walks(adj, 6))
print(f"New method: {len(new_walks)} L=6 canonical walks")

# Classify by unique vertex count
unique_v_counts = Counter()
for walk in new_walks:
    n_unique = len(set(walk))
    unique_v_counts[n_unique] += 1
print(f"Distribution of unique vertices per walk (L=6 walks have 6 vertex slots):")
for k, c in sorted(unique_v_counts.items()):
    print(f"  {k} unique vertices: {c} walks")

# So walks with 6 unique vertices = "simple cycles" (Hamiltonian-like)
# Walks with <6 unique = vertex repetition (figure-eights, theta-graphs, etc.)
print()
print(f"  simple cycles (6 unique): {unique_v_counts.get(6, 0)}")
print(f"  walks with vertex repetition: {sum(c for k, c in unique_v_counts.items() if k < 6)}")

# Now count walks via pilot method (oriented edges)
from xcwg_observables_pilot import enumerate_primitive_loops
pilot_walks = enumerate_primitive_loops(A, max_length=6)
pilot_6 = [w for w in pilot_walks if len(w) == 6]
print()
print(f"Pilot method: {len(pilot_6)} L=6 walks")

# Classify pilot walks too
pilot_unique_counts = Counter()
for w in pilot_6:
    vs = set()
    for (u, v) in w:
        vs.add(u); vs.add(v)
    pilot_unique_counts[len(vs)] += 1
print("Pilot distribution:")
for k, c in sorted(pilot_unique_counts.items()):
    print(f"  {k} unique vertices: {c} walks")

# Diagnose: pilot's filter at line 119-122 might reject closure-backtrack cases differently
# Let me check the simplest non-overlapping case
print()
print("First 3 new-method walks with 6 unique vertices:")
for w in new_walks[:5]:
    if len(set(w)) == 6:
        print(f"  walk={w}")
