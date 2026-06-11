"""Minimal import-sanity probe for Sprint L3b-2f-beta.1.

Test whether the geovac module imports cleanly. If this hangs, the
environmental issue from L3b-2f-alpha is still present.
"""
from __future__ import annotations
import sys
import time

print("[start] python interpreter ready", flush=True)
t0 = time.time()

print("[import 1/3] importing numpy...", flush=True)
import numpy as np
print(f"  numpy {np.__version__} imported in {time.time() - t0:.2f}s", flush=True)

t1 = time.time()
print("[import 2/3] importing geovac.krein_space_compact_temporal...", flush=True)
from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
print(f"  imported in {time.time() - t1:.2f}s", flush=True)

t2 = time.time()
print("[import 3/3] importing geovac.lorentzian_dirac_compact + operator_system...", flush=True)
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import CompactTemporalTruncatedOperatorSystem
print(f"  imported in {time.time() - t2:.2f}s", flush=True)

print(f"[done] all imports OK in {time.time() - t0:.2f}s total", flush=True)
