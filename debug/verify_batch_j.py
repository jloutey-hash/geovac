"""Verify batch cross-atom J matches original quad-based computation."""
import time
import sys
import os
sys.path.insert(0, '.')

import numpy as np
from geovac.lattice_index import (
    compute_cross_atom_J, compute_cross_atom_J_batch, _form_factor_nl,
    _form_factor_nl_grid
)

# First verify form factor grid matches scalar version
print("=== Form Factor Grid vs Scalar ===")
q_grid = np.linspace(0.01, 50, 200)
test_cases = [(1, 0, 5.0), (2, 0, 5.0), (2, 1, 5.0), (3, 0, 5.0),
              (3, 1, 5.0), (3, 2, 5.0), (1, 0, 11.0), (3, 2, 17.0)]

for n, l, Z in test_cases:
    grid_vals = _form_factor_nl_grid(n, l, Z, q_grid)
    scalar_vals = np.array([_form_factor_nl(n, l, Z, q) for q in q_grid])
    max_err = np.max(np.abs(grid_vals - scalar_vals))
    rel_err = np.max(np.abs(grid_vals - scalar_vals) / (np.abs(scalar_vals) + 1e-15))
    print(f"  (n={n},l={l},Z={Z:.0f}): max_abs_err={max_err:.2e}, max_rel_err={rel_err:.2e}")

# Now verify J integrals
print("\n=== Cross-Atom J: Batch vs Original ===")
# Clear caches to force fresh computation
cache_dir = os.path.join('geovac', 'cache')
for f in os.listdir(cache_dir):
    if f.startswith('cross_atom_J_') and 'R3.0150' in f:
        os.remove(os.path.join(cache_dir, f))

pairs = [
    (1, 0, 1, 0), (1, 0, 2, 0), (1, 0, 2, 1), (1, 0, 3, 0),
    (1, 0, 3, 1), (1, 0, 3, 2),
    (2, 0, 1, 0), (2, 0, 2, 0), (2, 0, 2, 1),
    (2, 1, 1, 0), (2, 1, 2, 0), (2, 1, 2, 1), (2, 1, 3, 2),
    (3, 0, 1, 0), (3, 0, 3, 2),
    (3, 1, 1, 0), (3, 1, 3, 1),
    (3, 2, 1, 0), (3, 2, 2, 1), (3, 2, 3, 2),
]

R = 3.015
ZA = 3.0
ZB = 1.0

# Batch computation
t0 = time.perf_counter()
batch_results = compute_cross_atom_J_batch(pairs, R, ZA, ZB)
t_batch = time.perf_counter() - t0
print(f"Batch time: {t_batch:.3f}s")

# Clear caches again for fair comparison
for f in os.listdir(cache_dir):
    if f.startswith('cross_atom_J_') and 'R3.0150' in f:
        os.remove(os.path.join(cache_dir, f))

# Original quad computation
t0 = time.perf_counter()
quad_results = {}
for na, la, nb, lb in pairs:
    quad_results[(na, la, nb, lb)] = compute_cross_atom_J(na, la, nb, lb, R, ZA, ZB)
t_quad = time.perf_counter() - t0
print(f"Quad time: {t_quad:.3f}s")
print(f"Speedup: {t_quad/t_batch:.1f}x")

print(f"\n{'Pair':<16} {'Batch':>12} {'Quad':>12} {'Abs Err':>12} {'Rel Err':>12}")
print("-" * 68)
max_rel = 0
for key in sorted(pairs):
    b = batch_results[key]
    q = quad_results[key]
    abs_err = abs(b - q)
    rel_err = abs_err / (abs(q) + 1e-15)
    max_rel = max(max_rel, rel_err)
    print(f"  {str(key):<14} {b:12.8f} {q:12.8f} {abs_err:12.2e} {rel_err:12.2e}")

print(f"\nMax relative error: {max_rel:.2e}")
if max_rel < 1e-4:
    print("PASS: Batch agrees with quad to < 0.01%")
else:
    print("WARN: Relative error exceeds 0.01%")

# Now test with higher Z (Na/Cl scales)
print("\n=== High-Z Test (Z=11, Z=17) ===")
for f in os.listdir(cache_dir):
    if f.startswith('cross_atom_J_') and 'R4.0000' in f:
        os.remove(os.path.join(cache_dir, f))

pairs_hZ = [
    (1, 0, 1, 0), (2, 1, 2, 1), (3, 2, 3, 2),
    (1, 0, 3, 2), (3, 2, 1, 0),
]

t0 = time.perf_counter()
batch_hZ = compute_cross_atom_J_batch(pairs_hZ, 4.0, 11.0, 17.0)
t_batch_hZ = time.perf_counter() - t0

for f in os.listdir(cache_dir):
    if f.startswith('cross_atom_J_') and 'R4.0000' in f:
        os.remove(os.path.join(cache_dir, f))

t0 = time.perf_counter()
quad_hZ = {}
for na, la, nb, lb in pairs_hZ:
    quad_hZ[(na, la, nb, lb)] = compute_cross_atom_J(na, la, nb, lb, 4.0, 11.0, 17.0)
t_quad_hZ = time.perf_counter() - t0

print(f"Batch: {t_batch_hZ:.3f}s, Quad: {t_quad_hZ:.3f}s, Speedup: {t_quad_hZ/t_batch_hZ:.1f}x")
print(f"\n{'Pair':<16} {'Batch':>12} {'Quad':>12} {'Rel Err':>12}")
for key in sorted(pairs_hZ):
    b = batch_hZ[key]
    q = quad_hZ[key]
    rel_err = abs(b - q) / (abs(q) + 1e-15)
    print(f"  {str(key):<14} {b:12.8f} {q:12.8f} {rel_err:12.2e}")
