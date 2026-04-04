"""
Track AD Part 2: Test R_ref candidates on LiH at l_max=2,3.

Tests top 3 candidates: r_half_screen, r_pk_width, r_avg.
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.composed_diatomic import ComposedDiatomicSolver

R_EQ_REF = 3.015

R_grid = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 10),
    np.linspace(4.5, 7.0, 5),
])

# Get candidate values first
rref_method = sys.argv[1] if len(sys.argv) > 1 else 'r_half_screen'
l_max = int(sys.argv[2]) if len(sys.argv) > 2 else 2

print(f"Testing R_ref method='{rref_method}', l_max={l_max}")
print(flush=True)

t0 = time.time()
solver = ComposedDiatomicSolver.LiH_ab_initio(
    l_max=l_max,
    pk_rref=rref_method,
    pk_channel_mode='l_dependent',
    verbose=True,
)
solver.solve_core()
print(f"  R_ref = {solver.pk_rref:.4f} bohr")

pes = solver.scan_pes(R_grid=R_grid)

R_eq = pes['R_eq']
err_pct = (R_eq - R_EQ_REF) / R_EQ_REF * 100
elapsed = time.time() - t0

print(f"\nRESULT: rref={rref_method} ({solver.pk_rref:.4f}), l_max={l_max}")
print(f"  R_eq = {R_eq:.4f} bohr")
print(f"  err  = {err_pct:+.1f}%")
print(f"  D_e  = {pes['D_e']:.6f} Ha")
print(f"  E_min = {pes['E_min']:.6f} Ha")
print(f"  time = {elapsed:.0f}s")
