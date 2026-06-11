"""Detailed comparison: pilot 144 vs new 568 L=6 walks at n_max=2.

Determine whether the pilot is rejecting figure-eights via some path/canon logic,
or whether the new method has a bug.
"""
import os, sys
from collections import Counter
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)

import numpy as np
from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from xcwg_wilson_loop_scaling import adjacency_list, enumerate_primitive_closed_walks

A, labels, deg, desc = build_dirac_s3_graph(2, "B")
adj = adjacency_list(A)
V = A.shape[0]

# Take ONE figure-eight walk from new method
new_walks = list(enumerate_primitive_closed_walks(adj, 6))
figure_eights = [w for w in new_walks if len(set(w)) == 5]
print(f"Figure-eight examples (5 unique vertices, 6 walk slots):")
for w in figure_eights[:3]:
    print(f"  walk={w}, unique={sorted(set(w))}")

# Take one and convert to oriented-edge form for pilot comparison
example = figure_eights[0]
L = len(example)
print(f"\nExample: {example}")
oriented_edges = []
for i in range(L):
    u, v = example[i], example[(i+1) % L]
    oriented_edges.append((u, v))
print(f"  oriented edges: {oriented_edges}")
print(f"  primitive in EDGE sense? {all(oriented_edges[i:i+1] != oriented_edges[0:1] for i in range(1,L))}")

# Check non-backtracking
print("  non-backtracking check:")
for i in range(L):
    u, v = oriented_edges[i]
    nu, nv = oriented_edges[(i+1) % L]
    if v != nu:
        print(f"    edge {i}->{i+1}: {(u,v)} -> {(nu,nv)} BROKEN (v={v} != nu={nu})")
    if u == nv:
        print(f"    edge {i}->{i+1}: {(u,v)} -> {(nu,nv)} BACKTRACK (u={u} == nv={nv})")
print()

# Now check whether pilot's enumeration would find this exact edge walk
# Pilot starts from each (u0, v0) oriented edge and DFS-extends...
# Let's actually run the pilot enumeration and look for THIS walk

from xcwg_observables_pilot import enumerate_primitive_loops
pilot_walks = enumerate_primitive_loops(A, max_length=6)
pilot_6 = [w for w in pilot_walks if len(w) == 6]
print(f"Pilot found {len(pilot_6)} L=6 walks")

# Look for the figure-eight walk in pilot output
# The pilot canon'd the walk, so we need to look at canonical form
# Pilot canon uses lex-min over rotations of oriented edges + reversal
def pilot_canon(walk_tuple):
    L = len(walk_tuple)
    best = walk_tuple
    for k in range(L):
        rot = walk_tuple[k:] + walk_tuple[:k]
        if rot < best:
            best = rot
    rev_walk = tuple((v, u) for (u, v) in reversed(walk_tuple))
    for k in range(L):
        rot = rev_walk[k:] + rev_walk[:k]
        if rot < best:
            best = rot
    return best

# Convert our figure-eight to oriented-edge form, then canon
oe_tuple = tuple(oriented_edges)
canon_oe = pilot_canon(oe_tuple)
print(f"\nFigure-eight as canonical-oriented-edge: {canon_oe}")
print(f"  is in pilot output? {canon_oe in pilot_6}")

# Now sanity-check: do any of the pilot walks have vertex repetition?
print()
print("Pilot walks by unique-vertex count:")
pilot_uv = Counter()
for w in pilot_6:
    vs = set()
    for (u, v) in w:
        vs.add(u); vs.add(v)
    pilot_uv[len(vs)] += 1
for k in sorted(pilot_uv):
    print(f"  {k} unique vertices: {pilot_uv[k]}")
