"""Convergence study: Neumann V_ee vs numerical V_ee for H2.

Compares energy convergence with exact (Neumann) vs grid-limited (numerical)
electron-electron repulsion at various basis sizes.
"""

import time
import sys
import numpy as np

from geovac.hylleraas import (
    generate_basis, build_quadrature_grids, solve_hylleraas, optimize_alpha,
)

R = 1.4011
E_exact = -1.17475
D_e_exact = 0.1745

grids = build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)

print("=" * 70)
print("NEUMANN V_ee CONVERGENCE STUDY")
print(f"R = {R}, E_exact = {E_exact}, D_e_exact = {D_e_exact}")
print("=" * 70)
sys.stdout.flush()

configs = [
    (0, 0, 0, "j=0,l=0"),
    (1, 0, 0, "j=1,l=0"),
    (1, 1, 0, "j=1,l=1"),
    (2, 1, 0, "j=2,l=1"),
    (2, 2, 0, "j=2,l=2"),
    (3, 2, 0, "j=3,l=2"),
    (3, 3, 0, "j=3,l=3"),
]

results = []

for j_max, l_max, p_max, label in configs:
    def basis_gen(alpha, j_max=j_max, l_max=l_max, p_max=p_max):
        return generate_basis(j_max=j_max, l_max=l_max, p_max=p_max, alpha=alpha)

    n_bf = len(basis_gen(1.0))
    print(f"\n--- {label} ({n_bf} bf) ---")
    sys.stdout.flush()

    # Numerical V_ee
    t0 = time.time()
    try:
        opt_num = optimize_alpha(
            basis_gen, R, grids, alpha_range=(0.5, 2.5),
            verbose=False, vee_method='numerical',
        )
        E_num = opt_num['E_opt']
        alpha_num = opt_num['alpha_opt']
        dt_num = time.time() - t0
        D_e_num = -1.0 - E_num
        pct_num = D_e_num / D_e_exact * 100
        print(f"  Numerical: E={E_num:.6f}, D_e={pct_num:.1f}%, "
              f"alpha={alpha_num:.3f}, {dt_num:.0f}s")
    except Exception as e:
        E_num = None
        dt_num = time.time() - t0
        print(f"  Numerical FAILED: {e}")
    sys.stdout.flush()

    # Neumann V_ee
    t0 = time.time()
    try:
        opt_neu = optimize_alpha(
            basis_gen, R, grids, alpha_range=(0.5, 2.5),
            verbose=False, vee_method='neumann', l_max_neumann=20,
        )
        E_neu = opt_neu['E_opt']
        alpha_neu = opt_neu['alpha_opt']
        dt_neu = time.time() - t0
        D_e_neu = -1.0 - E_neu
        pct_neu = D_e_neu / D_e_exact * 100
        print(f"  Neumann:   E={E_neu:.6f}, D_e={pct_neu:.1f}%, "
              f"alpha={alpha_neu:.3f}, {dt_neu:.0f}s")
    except Exception as e:
        E_neu = None
        dt_neu = time.time() - t0
        print(f"  Neumann FAILED: {e}")
    sys.stdout.flush()

    if E_num is not None and E_neu is not None:
        dE = E_neu - E_num
        print(f"  dE = {dE:.6f} Ha")
        results.append((label, n_bf, E_num, pct_num, dt_num,
                         E_neu, pct_neu, dt_neu, alpha_num, alpha_neu))

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"{'Config':<12} {'N':>4} {'E_num':>10} {'D_e%':>6} {'t':>5} "
      f"{'E_neu':>10} {'D_e%':>6} {'t':>5} {'dE':>10}")
print("-" * 80)
for r in results:
    label, n_bf, E_num, pct_num, dt_num, E_neu, pct_neu, dt_neu, a_num, a_neu = r
    dE = E_neu - E_num
    print(f"{label:<12} {n_bf:>4} {E_num:>10.6f} {pct_num:>5.1f}% {dt_num:>4.0f}s "
          f"{E_neu:>10.6f} {pct_neu:>5.1f}% {dt_neu:>4.0f}s {dE:>10.6f}")

print(f"\nExact: E = {E_exact}, D_e = {D_e_exact} Ha")
sys.stdout.flush()
