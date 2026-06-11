"""Quick test: verify 2D solver works at a single R point for LiH."""
import sys
sys.path.insert(0, '.')

from geovac.composed_diatomic import ComposedDiatomicSolver

# Test that the parameter is accepted and solver runs
solver = ComposedDiatomicSolver.LiH_ab_initio(
    l_max=2,
    level4_method='variational_2d',
    pk_channel_mode='l_dependent',
    verbose=True,
)
print(f"level4_method = {solver.level4_method}")
solver.solve_core()

# Single-point test at equilibrium
import numpy as np
E = solver._solve_valence_at_R(3.015, n_Re=150)
print(f"\nE_elec at R=3.015: {E:.6f} Ha")
print("SUCCESS: 2D solver accepts PK and runs in composed pipeline")
