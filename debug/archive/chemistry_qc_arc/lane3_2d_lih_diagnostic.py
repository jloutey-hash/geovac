"""
Lane 3 diagnostic: Test LiH composed geometry with 2D variational solver.

Compares adiabatic vs variational_2d at l_max=2,3,4.
Key question: does R_eq diverge with l_max (adiabatic) or stabilize (2D)?
"""
import sys
import time
import numpy as np

sys.path.insert(0, '.')
from geovac.composed_diatomic import ComposedDiatomicSolver

# Reference: R_eq = 3.015 bohr (expt)
R_EQ_EXPT = 3.015

# Denser grid near equilibrium for better R_eq resolution
R_grid = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 14),
    np.linspace(4.5, 7.0, 5),
])

results = {}

for l_max in [2, 3, 4]:
    for method in ['adiabatic', 'variational_2d']:
        label = f"l_max={l_max}, {method}"
        print(f"\n{'='*60}")
        print(f"  {label}")
        print(f"{'='*60}")

        t0 = time.time()
        try:
            solver = ComposedDiatomicSolver.LiH_ab_initio(
                l_max=l_max,
                level4_method=method,
                pk_channel_mode='l_dependent',
                verbose=True,
            )
            solver.solve_core()
            pes = solver.scan_pes(R_grid=R_grid, n_Re=200)

            R_eq = pes['R_eq']
            D_e = pes['D_e']
            err_pct = abs(R_eq - R_EQ_EXPT) / R_EQ_EXPT * 100
            wall = time.time() - t0

            results[(l_max, method)] = {
                'R_eq': R_eq, 'D_e': D_e,
                'err_pct': err_pct, 'wall_s': wall,
            }
            print(f"\n  RESULT: R_eq={R_eq:.3f} ({err_pct:.1f}% error), "
                  f"D_e={D_e:.6f} Ha, wall={wall:.1f}s")
        except Exception as e:
            wall = time.time() - t0
            results[(l_max, method)] = {'error': str(e), 'wall_s': wall}
            print(f"\n  FAILED: {e} (wall={wall:.1f}s)")

# Summary table
print(f"\n\n{'='*70}")
print(f"  SUMMARY: LiH 2D vs Adiabatic")
print(f"{'='*70}")
print(f"  {'l_max':>5s}  {'method':>15s}  {'R_eq':>7s}  {'err%':>6s}  {'D_e':>10s}  {'wall':>6s}")
print(f"  {'-'*5}  {'-'*15}  {'-'*7}  {'-'*6}  {'-'*10}  {'-'*6}")
for (l, m), r in sorted(results.items()):
    if 'error' in r:
        print(f"  {l:5d}  {m:>15s}  {'FAIL':>7s}  {'':>6s}  {'':>10s}  {r['wall_s']:6.1f}")
    else:
        print(f"  {l:5d}  {m:>15s}  {r['R_eq']:7.3f}  {r['err_pct']:6.1f}  "
              f"{r['D_e']:10.6f}  {r['wall_s']:6.1f}")

# Check l_max convergence for each method
print(f"\n  l_max convergence (R_eq drift):")
for method in ['adiabatic', 'variational_2d']:
    r_vals = []
    for l in [2, 3, 4]:
        r = results.get((l, method))
        if r and 'R_eq' in r:
            r_vals.append(r['R_eq'])
    if len(r_vals) >= 2:
        drift = r_vals[-1] - r_vals[0]
        print(f"    {method:>15s}: R_eq drift l2->l4 = {drift:+.3f} bohr "
              f"({'+' if drift > 0 else ''}{drift/(4-2):.3f}/l_max)")
