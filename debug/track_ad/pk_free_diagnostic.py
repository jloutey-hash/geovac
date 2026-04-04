"""
Track AD Part 1: PK-free diagnostic for LiH l_max drift.

Runs LiH composed solver at l_max=1,2,3,4 with pk_mode='none'
and both level4_method='adiabatic' and 'variational_2d'.
Compares to pk_mode='ab_initio' baseline.
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.composed_diatomic import ComposedDiatomicSolver

R_EQ_REF = 3.015  # LiH experimental R_eq

# Dense R grid near equilibrium for accurate R_eq
R_grid = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 10),
    np.linspace(4.5, 7.0, 5),
])

results = []

for pk_mode in ['none', 'ab_initio']:
    for level4_method in ['adiabatic', 'variational_2d']:
        for l_max in [1, 2, 3]:
            label = f"pk={pk_mode}, method={level4_method}, l_max={l_max}"
            print(f"\n{'='*60}")
            print(f"  {label}")
            print(f"{'='*60}")

            t0 = time.time()
            try:
                solver = ComposedDiatomicSolver.LiH_ab_initio(
                    l_max=l_max,
                    pk_mode=pk_mode,
                    level4_method=level4_method,
                    pk_channel_mode='l_dependent',
                    verbose=True,
                )
                solver.solve_core()
                pes = solver.scan_pes(R_grid=R_grid)

                R_eq = pes['R_eq']
                E_min = pes['E_min']
                D_e = pes['D_e']
                err_pct = (R_eq - R_EQ_REF) / R_EQ_REF * 100

                elapsed = time.time() - t0
                results.append({
                    'pk_mode': pk_mode,
                    'method': level4_method,
                    'l_max': l_max,
                    'R_eq': R_eq,
                    'err_pct': err_pct,
                    'E_min': E_min,
                    'D_e': D_e,
                    'time': elapsed,
                })
                print(f"\n  >>> R_eq={R_eq:.3f}, err={err_pct:+.1f}%, "
                      f"D_e={D_e:.4f}, time={elapsed:.0f}s")

            except Exception as e:
                elapsed = time.time() - t0
                print(f"\n  >>> FAILED after {elapsed:.0f}s: {e}")
                results.append({
                    'pk_mode': pk_mode,
                    'method': level4_method,
                    'l_max': l_max,
                    'R_eq': np.nan,
                    'err_pct': np.nan,
                    'E_min': np.nan,
                    'D_e': np.nan,
                    'time': elapsed,
                })

# Summary table
print(f"\n\n{'='*80}")
print("SUMMARY: PK-free vs ab_initio PK diagnostic")
print(f"{'='*80}")
print(f"{'pk_mode':>12s} {'method':>16s} {'l_max':>5s} {'R_eq':>8s} "
      f"{'err%':>8s} {'D_e':>8s} {'time':>6s}")
print("-" * 80)
for r in results:
    print(f"{r['pk_mode']:>12s} {r['method']:>16s} {r['l_max']:>5d} "
          f"{r['R_eq']:>8.3f} {r['err_pct']:>+8.1f} {r['D_e']:>8.4f} "
          f"{r['time']:>6.0f}s")

# Compute drift rates
print(f"\n\nDRIFT ANALYSIS:")
print("-" * 60)
for pk_mode in ['none', 'ab_initio']:
    for method in ['adiabatic', 'variational_2d']:
        subset = [r for r in results if r['pk_mode'] == pk_mode
                  and r['method'] == method and not np.isnan(r['R_eq'])]
        if len(subset) >= 2:
            l_vals = [r['l_max'] for r in subset]
            R_vals = [r['R_eq'] for r in subset]
            # Linear regression for drift rate
            if len(l_vals) >= 2:
                slope = np.polyfit(l_vals, R_vals, 1)[0]
                print(f"  pk={pk_mode:>12s}, method={method:>16s}: "
                      f"drift = {slope:+.3f} bohr/l_max")

# Save results
outfile = os.path.join(os.path.dirname(__file__), 'pk_free_results.txt')
with open(outfile, 'w') as f:
    f.write("pk_mode,method,l_max,R_eq,err_pct,E_min,D_e,time\n")
    for r in results:
        f.write(f"{r['pk_mode']},{r['method']},{r['l_max']},"
                f"{r['R_eq']:.4f},{r['err_pct']:.2f},"
                f"{r['E_min']:.6f},{r['D_e']:.6f},{r['time']:.1f}\n")
print(f"\nResults saved to {outfile}")
