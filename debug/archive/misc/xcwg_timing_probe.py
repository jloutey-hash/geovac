"""Timing probe: how long does L=8 enumeration take at n_max=3?"""
import os, sys, time
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)

import numpy as np
from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from xcwg_wilson_loop_scaling import (
    signed_incidence, adjacency_list, enumerate_primitive_closed_walks,
)

for n_max in [3]:
    A, labels, deg, desc = build_dirac_s3_graph(n_max, 'B')
    V = A.shape[0]
    E = int(A.sum() // 2)
    print(f"n_max={n_max}: V={V}, E={E}")
    adj = adjacency_list(A)
    for L in [4, 6, 8]:
        t0 = time.time()
        count = 0
        for w in enumerate_primitive_closed_walks(adj, L):
            count += 1
            if count >= 50_000:
                break
        dt = time.time() - t0
        print(f"  L={L}: {count} canonical walks in {dt:.1f}s ({count/dt:.0f}/s)")
