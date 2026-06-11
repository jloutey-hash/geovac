"""Lane 3: LiH l_max=2 comparison — adiabatic vs variational_2d."""
import sys, time
import numpy as np
sys.path.insert(0, '.')
from geovac.composed_diatomic import ComposedDiatomicSolver

R_EQ_EXPT = 3.015
R_grid = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 14),
    np.linspace(4.5, 7.0, 5),
])

for method in ['adiabatic', 'variational_2d']:
    print(f"\n{'='*60}")
    print(f"  l_max=2, method={method}")
    print(f"{'='*60}")
    t0 = time.time()
    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2,
        level4_method=method,
        pk_channel_mode='l_dependent',
        verbose=True,
    )
    solver.solve_core()
    pes = solver.scan_pes(R_grid=R_grid, n_Re=200)
    wall = time.time() - t0
    R_eq = pes['R_eq']
    D_e = pes['D_e']
    err = abs(R_eq - R_EQ_EXPT) / R_EQ_EXPT * 100
    print(f"\n  RESULT: R_eq={R_eq:.3f} ({err:.1f}%), D_e={D_e:.6f}, wall={wall:.1f}s")
