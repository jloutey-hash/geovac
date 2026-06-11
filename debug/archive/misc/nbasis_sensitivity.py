"""
Test whether the ~0.21% error floor in the coupled-channel He solver
comes from n_basis truncation in the Gegenbauer spectral basis.

Sweeps n_basis = 10, 15, 20, 25, 30 at fixed l_max=3, n_channels=3,
q_mode='exact'.  Also tests N_R_radial=5000 at n_basis=15 to check
radial grid convergence.

Output: formatted table + debug/data/nbasis_sensitivity.json
"""

import json
import time
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

E_EXACT = -2.903724

# Common parameters
common = dict(
    Z=2.0,
    l_max=3,
    n_channels=3,
    n_R=200,
    N_R_radial=3000,
    q_mode='exact',
    verbose=True,
)

results = []

# --- n_basis sweep ---
for nb in [10, 15, 20, 25, 30]:
    print(f"\n{'='*60}")
    print(f"  n_basis = {nb}")
    print(f"{'='*60}")
    t0 = time.time()
    res = solve_hyperspherical_algebraic_coupled(n_basis=nb, **common)
    wall = time.time() - t0
    E = res['energy']
    err = abs((E - E_EXACT) / E_EXACT) * 100
    results.append({
        'n_basis': nb,
        'N_R_radial': 3000,
        'energy': float(E),
        'error_pct': float(err),
        'wall_time_s': round(wall, 2),
        'label': f'n_basis={nb}',
    })
    print(f"  E = {E:.6f} Ha, error = {err:.4f}%, time = {wall:.1f}s")

# --- radial grid convergence check ---
print(f"\n{'='*60}")
print(f"  n_basis=15, N_R_radial=5000 (radial grid check)")
print(f"{'='*60}")
t0 = time.time()
res = solve_hyperspherical_algebraic_coupled(
    n_basis=15, N_R_radial=5000, Z=2.0, l_max=3, n_channels=3,
    n_R=200, q_mode='exact', verbose=True,
)
wall = time.time() - t0
E = res['energy']
err = abs((E - E_EXACT) / E_EXACT) * 100
results.append({
    'n_basis': 15,
    'N_R_radial': 5000,
    'energy': float(E),
    'error_pct': float(err),
    'wall_time_s': round(wall, 2),
    'label': 'n_basis=15, N_R=5000',
})
print(f"  E = {E:.6f} Ha, error = {err:.4f}%, time = {wall:.1f}s")

# --- Print formatted table ---
print(f"\n\n{'='*72}")
print(f"  n_basis Sensitivity for He Coupled-Channel (l_max=3, n_ch=3, exact Q)")
print(f"{'='*72}")
print(f"  {'Label':<25s} {'Energy (Ha)':>14s} {'Error %':>10s} {'Time (s)':>10s}")
print(f"  {'-'*25} {'-'*14} {'-'*10} {'-'*10}")
for r in results:
    print(f"  {r['label']:<25s} {r['energy']:>14.6f} {r['error_pct']:>10.4f} {r['wall_time_s']:>10.1f}")
print(f"  {'-'*25} {'-'*14} {'-'*10} {'-'*10}")
print(f"  Exact He:               {E_EXACT:>14.6f}")

# --- Save JSON ---
out_path = os.path.join(os.path.dirname(__file__), 'data', 'nbasis_sensitivity.json')
with open(out_path, 'w') as f:
    json.dump({'exact_energy': E_EXACT, 'results': results}, f, indent=2)
print(f"\nSaved to {out_path}")
