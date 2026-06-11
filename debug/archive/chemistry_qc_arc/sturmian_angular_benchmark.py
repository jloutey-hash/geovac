"""
Benchmark script for Sturmian angular basis investigation.

The key insight: the Sturmian basis is constructed in a LARGE free basis
(n_construct=50), then truncated to n_basis eigenvectors.  These n_basis
vectors span a better-adapted subspace than the first n_basis free functions.

Runs Steps 1-4:
  1. Monopole decomposition
  2. R0 optimization
  3. l_max=0 benchmarks (Sturmian vs free basis at matched effective dimension)
  4. l_max=1 benchmarks
"""

import numpy as np
import json
import time
from math import pi

N_CONSTRUCT = 50  # Large free basis for Sturmian construction

# Step 1: Monopole decomposition
print("=" * 70)
print("STEP 1: Monopole Decomposition")
print("=" * 70)

from geovac.algebraic_angular_sturmian import SturmianAngularSolver

for l_max in [0, 1, 2]:
    solver = SturmianAngularSolver(
        Z=2.0, n_basis=15, l_max=l_max, R0=1.0, n_construct=N_CONSTRUCT,
    )
    decomp = solver.monopole_decomposition()
    offdiag = solver.off_diagonal_analysis()

    print(f"\nl_max = {l_max} (n_construct={N_CONSTRUCT}, n_basis=15):")
    print(f"  V_nuc Frobenius:      {decomp['V_nuc_frobenius']:.4f}")
    print(f"  V_ee  Frobenius:      {decomp['V_ee_frobenius']:.4f}")
    print(f"  V_total Frobenius:    {decomp['V_total_frobenius']:.4f}")
    print(f"  Monopole fraction:    {decomp['monopole_fraction']:.1%}")
    print(f"  Remainder fraction:   {decomp['remainder_fraction']:.1%}")
    print(f"  Projection error:     {decomp['projection_error']:.2e}")
    print(f"  Mono off-diag (proj): {offdiag['monopole_offdiag_norm']:.6f}")
    print(f"  V_ee off-diag (proj): {offdiag['remainder_offdiag_norm']:.6f}")

# Step 2: R0 optimization
print("\n" + "=" * 70)
print("STEP 2: R0 Optimization (n_basis=10, n_construct=50, l_max=0)")
print("=" * 70)

from geovac.algebraic_angular_sturmian import optimize_R0

R0_result = optimize_R0(
    Z=2.0, n_basis=10, l_max=0, n_construct=N_CONSTRUCT,
    R0_values=[0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0],
    n_R=150, N_R_radial=1500,
)
print(f"\nBest R0 = {R0_result['best_R0']:.1f}, "
      f"E = {R0_result['best_energy']:.6f} Ha")

best_R0 = R0_result['best_R0']

# Step 3: l_max=0 benchmarks
print("\n" + "=" * 70)
print("STEP 3: l_max=0 Benchmarks — Sturmian (truncated) vs Free Basis")
print("=" * 70)

from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian
from geovac.algebraic_angular import solve_hyperspherical_algebraic

E_exact = -2.903724

# Compare: n_basis Sturmian functions (from n_construct=50) vs n_basis free functions
basis_sizes = [3, 5, 8, 10, 15, 20, 25]

print(f"\nSturmian: n_construct={N_CONSTRUCT}, R0={best_R0}")
print(f"Free: standard sin(2n*alpha) basis")
print(f"\n{'n_basis':>8s}  {'E_Sturm':>12s}  {'err_S(%)':>10s}  "
      f"{'E_Free':>12s}  {'err_F(%)':>10s}  {'improvement':>12s}")

results_sturm = {}
results_free = {}

for nb in basis_sizes:
    res_s = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=nb, l_max=0, R0=best_R0, n_construct=N_CONSTRUCT,
        n_R=200, R_max=30.0, N_R_radial=2000, verbose=False,
    )
    res_f = solve_hyperspherical_algebraic(
        Z=2.0, n_basis=nb, l_max=0,
        n_R=200, R_max=30.0, N_R_radial=2000, verbose=False,
    )

    err_s = abs(res_s['energy'] - E_exact) / abs(E_exact) * 100
    err_f = abs(res_f['energy'] - E_exact) / abs(E_exact) * 100

    if err_f > 1e-10:
        improvement = f"{err_f / err_s:.2f}x" if err_s > 1e-10 else "inf"
    else:
        improvement = "N/A"

    results_sturm[nb] = {'energy': res_s['energy'], 'error': err_s}
    results_free[nb] = {'energy': res_f['energy'], 'error': err_f}

    print(f"  {nb:6d}  {res_s['energy']:12.6f}  {err_s:10.4f}  "
          f"{res_f['energy']:12.6f}  {err_f:10.4f}  {improvement:>12s}")

# Step 3b: n_construct sensitivity
print("\n--- n_construct sensitivity (n_basis=10, R0=" + f"{best_R0}) ---")
for nc in [15, 20, 30, 40, 50, 60, 80]:
    res = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=10, l_max=0, R0=best_R0, n_construct=nc,
        n_R=200, R_max=30.0, N_R_radial=2000, verbose=False,
    )
    err = abs(res['energy'] - E_exact) / abs(E_exact) * 100
    print(f"  n_construct={nc:3d}: E={res['energy']:.6f}, err={err:.4f}%")

# Step 4: l_max comparison
print("\n" + "=" * 70)
print("STEP 4: l_max Comparison (n_basis=15)")
print("=" * 70)

# Optimize R0 for l_max=1
print("\nOptimizing R0 for l_max=1...")
R0_l1 = optimize_R0(
    Z=2.0, n_basis=10, l_max=1, n_construct=N_CONSTRUCT,
    R0_values=[0.5, 1.0, 1.5, 2.0, 3.0, 5.0],
    n_R=150, N_R_radial=1500,
)
best_R0_l1 = R0_l1['best_R0']
print(f"Best R0 for l_max=1: {best_R0_l1}")

print(f"\n{'l_max':>6s}  {'E_Sturm':>12s}  {'err_S(%)':>10s}  "
      f"{'E_Free':>12s}  {'err_F(%)':>10s}  {'R0':>6s}")

results_lmax = {}
for l_max in [0, 1, 2, 3]:
    R0_use = best_R0 if l_max == 0 else best_R0_l1

    res_s = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=15, l_max=l_max, R0=R0_use, n_construct=N_CONSTRUCT,
        n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
    )
    res_f = solve_hyperspherical_algebraic(
        Z=2.0, n_basis=15, l_max=l_max,
        n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
    )

    err_s = abs(res_s['energy'] - E_exact) / abs(E_exact) * 100
    err_f = abs(res_f['energy'] - E_exact) / abs(E_exact) * 100

    results_lmax[l_max] = {
        'E_sturm': res_s['energy'], 'err_sturm': err_s,
        'E_free': res_f['energy'], 'err_free': err_f,
        'R0': R0_use,
    }

    print(f"  {l_max:5d}  {res_s['energy']:12.6f}  {err_s:10.4f}  "
          f"{res_f['energy']:12.6f}  {err_f:10.4f}  {R0_use:6.1f}")

# Step 5: mu(R) convergence at fixed R
print("\n" + "=" * 70)
print("STEP 5: mu(R) Convergence at Fixed R (l_max=0)")
print("=" * 70)

from geovac.algebraic_angular import AlgebraicAngularSolver as FreeSolver

# Reference: large free basis
ref = FreeSolver(Z=2.0, n_basis=80, l_max=0, n_quad=200)

for R in [1.0, 2.0, 5.0, 10.0]:
    evals_ref, _ = ref.solve(R, n_channels=1)
    mu_ref = evals_ref[0]

    print(f"\nR = {R:.1f} bohr (ref mu = {mu_ref:.8f}):")
    print(f"  {'n_basis':>8s}  {'mu_Sturm':>14s}  {'err_S':>12s}  "
          f"{'mu_Free':>14s}  {'err_F':>12s}")

    for nb in [3, 5, 8, 10, 15]:
        s = SturmianAngularSolver(
            Z=2.0, n_basis=nb, l_max=0, R0=best_R0, n_construct=N_CONSTRUCT,
        )
        f = FreeSolver(Z=2.0, n_basis=nb, l_max=0)

        evals_s, _ = s.solve(R, n_channels=1)
        evals_f, _ = f.solve(R, n_channels=1)

        err_s = abs(evals_s[0] - mu_ref)
        err_f = abs(evals_f[0] - mu_ref)

        print(f"  {nb:8d}  {evals_s[0]:14.8f}  {err_s:12.2e}  "
              f"{evals_f[0]:14.8f}  {err_f:12.2e}")

# Save all results
all_results = {
    'R0_optimization': {
        'R0_values': R0_result['R0_values'],
        'energies': R0_result['energies'],
        'best_R0': R0_result['best_R0'],
    },
    'l_max_0_benchmark': {
        'sturmian': {str(k): v for k, v in results_sturm.items()},
        'free': {str(k): v for k, v in results_free.items()},
    },
    'l_max_comparison': {str(k): v for k, v in results_lmax.items()},
    'n_construct': N_CONSTRUCT,
}

with open('debug/data/sturmian_angular_benchmark.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\n\nResults saved to debug/data/sturmian_angular_benchmark.json")
