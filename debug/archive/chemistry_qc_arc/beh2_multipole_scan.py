"""
BeH2 PES scan comparing all interbond coupling modes.

Modes tested:
  1. block_diagonal (no interbond coupling)
  2. exchange (S·F⁰ factorized)
  3. direct_exchange (channel-resolved, no S·F⁰ factorization)
  4. multipole_exchange k_max=0 (should match direct_exchange)
  5. multipole_exchange k_max=1 (+ dipole)
  6. multipole_exchange k_max=2 (+ quadrupole)

Reports R_eq for each mode and comparison to reference (2.507 bohr).
"""

import json
import time
import numpy as np
from pathlib import Path

from geovac.composed_triatomic import ComposedTriatomicSolver

# Reference
R_EQ_REF = 2.507  # bohr

R_grid = np.arange(1.5, 5.1, 0.2)  # coarser grid for speed

modes = [
    ('block_diagonal', dict(include_interbond=False)),
    ('exchange', dict(interbond_mode='exchange')),
    ('direct_exchange', dict(interbond_mode='direct_exchange')),
    ('multipole_k0', dict(interbond_mode='multipole_exchange', k_max=0)),
    ('multipole_k1', dict(interbond_mode='multipole_exchange', k_max=1)),
    ('multipole_k2', dict(interbond_mode='multipole_exchange', k_max=2)),
]

results = {}

for name, kwargs in modes:
    print(f"\n{'='*64}")
    print(f"  Mode: {name}")
    print(f"{'='*64}")

    t0 = time.time()
    solver = ComposedTriatomicSolver.BeH2(l_max=2, verbose=True, **kwargs)
    result = solver.run_all(R_grid=R_grid, n_Re=300)
    elapsed = time.time() - t0

    R_eq = result['spectro']['R_eq']
    D_e = result['spectro']['D_e']
    err_pct = abs(R_eq - R_EQ_REF) / R_EQ_REF * 100

    results[name] = {
        'R_eq': R_eq,
        'D_e': D_e,
        'E_min': result['spectro']['E_min'],
        'R_eq_error_pct': err_pct,
        'elapsed_s': elapsed,
        'R': result['pes']['R'],
        'E_total': result['pes']['E_total'],
        'V_inter': result['pes']['V_inter'],
    }

    print(f"\n  >> {name}: R_eq = {R_eq:.4f} bohr, "
          f"error = {err_pct:.1f}%, D_e = {D_e:.6f} Ha, "
          f"time = {elapsed:.0f}s")

# Summary table
print(f"\n\n{'='*72}")
print(f"  SUMMARY: BeH2 R_eq comparison across interbond modes")
print(f"{'='*72}")
print(f"  {'Mode':<22s} {'R_eq':>8s} {'Error':>8s} {'D_e':>10s} {'Time':>8s}")
print(f"  {'-'*22} {'-'*8} {'-'*8} {'-'*10} {'-'*8}")
print(f"  {'Reference':<22s} {R_EQ_REF:8.4f} {'---':>8s} {'0.147':>10s} {'---':>8s}")
for name, r in results.items():
    print(f"  {name:<22s} {r['R_eq']:8.4f} {r['R_eq_error_pct']:7.1f}% "
          f"{r['D_e']:10.6f} {r['elapsed_s']:7.0f}s")

# Check multipole_k0 matches direct_exchange
if 'direct_exchange' in results and 'multipole_k0' in results:
    r_de = results['direct_exchange']['R_eq']
    r_mk0 = results['multipole_k0']['R_eq']
    match = abs(r_de - r_mk0) < 0.01
    print(f"\n  multipole_k0 matches direct_exchange: "
          f"{'YES' if match else 'NO'} "
          f"(delta = {abs(r_de - r_mk0):.4f})")

# Fitted target is 12% error
target = 0.12 * R_EQ_REF
print(f"\n  Target (fitted): {R_EQ_REF * (1 - 0.12):.3f} - "
      f"{R_EQ_REF * (1 + 0.12):.3f} bohr (12% error band)")

# Save data
out = {k: {kk: vv for kk, vv in v.items()
            if kk not in ('R', 'E_total', 'V_inter')}
       for k, v in results.items()}
out_file = Path('debug/data/beh2_multipole_exchange.json')
out_file.parent.mkdir(parents=True, exist_ok=True)
with open(out_file, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n  Data saved to {out_file}")
