"""
Track AD Part 2: Sweep R_ref values to find optimal and compare candidates.

Tests R_ref = 0.5, 1.0, 2.0, 3.0, 5.0 at l_max=2 to map sensitivity.
Also tests 'auto' (r_half_screen) for comparison.
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

l_max = 2
results = []

# Test several R_ref values + baseline (no R-dep PK)
rref_values = [
    ('baseline', None),
    ('auto', 'auto'),
    ('0.5', 0.5),
    ('1.0', 1.0),
    ('2.0', 2.0),
    ('3.0', 3.0),
    ('5.0', 5.0),
    ('10.0', 10.0),
]

for label, rref in rref_values:
    print(f"\n{'='*60}")
    print(f"  R_ref={label}, l_max={l_max}")
    print(f"{'='*60}")

    t0 = time.time()
    try:
        solver = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=l_max,
            pk_rref=rref,
            pk_channel_mode='l_dependent',
            verbose=True,
        )
        solver.solve_core()
        rref_val = solver.pk_rref
        pes = solver.scan_pes(R_grid=R_grid)

        R_eq = pes['R_eq']
        err_pct = (R_eq - R_EQ_REF) / R_EQ_REF * 100
        elapsed = time.time() - t0
        results.append({
            'label': label,
            'R_ref': rref_val,
            'R_eq': R_eq,
            'err_pct': err_pct,
            'D_e': pes['D_e'],
            'time': elapsed,
        })
        print(f"\n  >>> R_eq={R_eq:.3f}, err={err_pct:+.1f}%")

    except Exception as e:
        import traceback
        traceback.print_exc()
        elapsed = time.time() - t0
        results.append({
            'label': label,
            'R_ref': None,
            'R_eq': np.nan,
            'err_pct': np.nan,
            'D_e': np.nan,
            'time': elapsed,
        })

# Summary
print(f"\n\n{'='*70}")
print(f"R_ref SWEEP at l_max={l_max}")
print(f"{'='*70}")
print(f"{'label':>10s} {'R_ref':>8s} {'R_eq':>8s} {'err%':>8s} "
      f"{'D_e':>8s} {'time':>6s}")
print("-" * 70)
for r in results:
    rref_str = f"{r['R_ref']:.4f}" if r['R_ref'] is not None else "---"
    print(f"{r['label']:>10s} {rref_str:>8s} "
          f"{r['R_eq']:>8.3f} {r['err_pct']:>+8.1f} {r['D_e']:>8.4f} "
          f"{r['time']:>6.0f}s")

outfile = os.path.join(os.path.dirname(__file__), 'rref_sweep_results.txt')
with open(outfile, 'w') as f:
    f.write("label,R_ref,l_max,R_eq,err_pct,D_e,time\n")
    for r in results:
        rref_str = f"{r['R_ref']:.4f}" if r['R_ref'] is not None else ""
        f.write(f"{r['label']},{rref_str},{l_max},"
                f"{r['R_eq']:.4f},{r['err_pct']:.2f},"
                f"{r['D_e']:.6f},{r['time']:.1f}\n")
print(f"\nSaved to {outfile}")
