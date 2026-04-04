"""
Track AD Part 1 (minimal): PK-free diagnostic — one l_max at a time.
Run: python pk_free_minimal.py <pk_mode> <l_max>
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.composed_diatomic import ComposedDiatomicSolver

pk_mode = sys.argv[1] if len(sys.argv) > 1 else 'none'
l_max = int(sys.argv[2]) if len(sys.argv) > 2 else 2

R_EQ_REF = 3.015

R_grid = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 10),
    np.linspace(4.5, 7.0, 5),
])

print(f"pk_mode={pk_mode}, l_max={l_max}, adiabatic")
print(flush=True)

t0 = time.time()
solver = ComposedDiatomicSolver.LiH_ab_initio(
    l_max=l_max,
    pk_mode=pk_mode,
    level4_method='adiabatic',
    pk_channel_mode='l_dependent',
    verbose=True,
)
solver.solve_core()
print(f"Core solved in {time.time()-t0:.1f}s", flush=True)

pes = solver.scan_pes(R_grid=R_grid)

R_eq = pes['R_eq']
err_pct = (R_eq - R_EQ_REF) / R_EQ_REF * 100
elapsed = time.time() - t0

print(f"\nRESULT: pk_mode={pk_mode}, l_max={l_max}")
print(f"  R_eq = {R_eq:.4f} bohr")
print(f"  err  = {err_pct:+.1f}%")
print(f"  D_e  = {pes['D_e']:.6f} Ha")
print(f"  E_min = {pes['E_min']:.6f} Ha")
print(f"  time = {elapsed:.0f}s")
