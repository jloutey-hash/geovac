"""
Hylleraas convergence study for H2.

Phase 9 of the Prolate Spheroidal stress tests.

Tests systematic expansion of the James-Coolidge basis:
1. p_max=0 (CI limit — no explicit correlation)
2. p_max=1 (linear r₁₂ — cusp condition)
3. p_max=2 (quadratic r₁₂)

Tracks D_e convergence toward exact 0.1745 Ha.

Reference:
    James & Coolidge (1933): E = -1.17447 Ha with 13 terms
    Exact: E = -1.17475 Ha, D_e = 0.1745 Ha
"""

import numpy as np
import time
import sys
sys.path.insert(0, '.')

from geovac.hylleraas import (
    generate_basis,
    generate_basis_truncated,
    build_quadrature_grids,
    solve_hylleraas,
    optimize_alpha,
)


def run_p0_convergence():
    """Phase 9a: p=0 convergence (CI limit, fast 4D integration)."""
    print("=" * 70)
    print("PHASE 9a: p=0 CONVERGENCE (CI LIMIT)")
    print("=" * 70)
    print()

    R = 1.4
    D_e_exact = 0.1745

    # Medium grid — sufficient for p=0
    grids = build_quadrature_grids(N_xi=18, N_eta=14, N_phi=10, xi_max=12.0)

    results = []

    # Scan alpha for best value with a medium basis
    print("--- Alpha optimization (j_max=2, l_max=1) ---")
    best_alpha = 0.9
    best_E = 0.0
    for alpha in np.arange(0.5, 1.6, 0.1):
        basis = generate_basis(j_max=2, l_max=1, p_max=0, alpha=alpha)
        r = solve_hylleraas(basis, R, grids, verbose=False)
        if r['E_total'] < best_E:
            best_E = r['E_total']
            best_alpha = alpha
        print(f"  alpha={alpha:.1f}: E={r['E_total']:.6f}, "
              f"D_e={r['D_e']:.4f} ({r['D_e_pct']:.1f}%)")
    print(f"  Best alpha = {best_alpha:.1f}")
    print()

    # Convergence with basis size
    print("--- Basis size convergence (alpha={:.1f}) ---".format(best_alpha))
    print(f"{'j_max':>5s} {'l_max':>5s} {'N_bf':>5s} {'E_total':>12s} "
          f"{'D_e':>8s} {'%exact':>8s} {'Time':>8s}")
    print("-" * 60)

    for j_max, l_max in [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2),
                          (3, 1), (3, 2)]:
        basis = generate_basis(j_max=j_max, l_max=l_max, p_max=0,
                               alpha=best_alpha)
        t0 = time.time()
        r = solve_hylleraas(basis, R, grids, verbose=False)
        dt = time.time() - t0
        pct = r['D_e'] / D_e_exact * 100
        print(f"{j_max:5d} {l_max:5d} {len(basis):5d} {r['E_total']:12.6f} "
              f"{r['D_e']:8.4f} {pct:7.1f}% {dt:7.1f}s")
        results.append({
            'j_max': j_max, 'l_max': l_max, 'p_max': 0,
            'n_basis': len(basis), 'E_total': r['E_total'],
            'D_e': r['D_e'], 'D_e_pct': pct, 'time': dt,
            'alpha': best_alpha,
        })

    print()
    return results, best_alpha


def run_p1_convergence(best_alpha: float):
    """Phase 9b: p=1 convergence (linear r₁₂, 5D integration)."""
    print("=" * 70)
    print("PHASE 9b: p=0,1 CONVERGENCE (LINEAR r12)")
    print("=" * 70)
    print()
    print("NOTE: 5D integration is slow. Using minimal grids.")
    print()

    R = 1.4
    D_e_exact = 0.1745

    # Smaller grid for 5D integration
    grids = build_quadrature_grids(N_xi=10, N_eta=8, N_phi=8, xi_max=10.0)

    results = []

    print(f"{'j_max':>5s} {'l_max':>5s} {'N_bf':>5s} {'E_total':>12s} "
          f"{'D_e':>8s} {'%exact':>8s} {'Time':>8s}")
    print("-" * 60)

    for j_max, l_max in [(0, 0), (1, 0), (1, 1)]:
        basis = generate_basis(j_max=j_max, l_max=l_max, p_max=1,
                               alpha=best_alpha)
        t0 = time.time()
        r = solve_hylleraas(basis, R, grids, verbose=False)
        dt = time.time() - t0
        pct = r['D_e'] / D_e_exact * 100
        print(f"{j_max:5d} {l_max:5d} {len(basis):5d} {r['E_total']:12.6f} "
              f"{r['D_e']:8.4f} {pct:7.1f}% {dt:7.1f}s")
        results.append({
            'j_max': j_max, 'l_max': l_max, 'p_max': 1,
            'n_basis': len(basis), 'E_total': r['E_total'],
            'D_e': r['D_e'], 'D_e_pct': pct, 'time': dt,
            'alpha': best_alpha,
        })

    print()
    return results


def run_comparison():
    """Phase 9c: Direct p=0 vs p=1 comparison at same basis."""
    print("=" * 70)
    print("PHASE 9c: p=0 vs p=1 COMPARISON")
    print("=" * 70)
    print()

    R = 1.4
    D_e_exact = 0.1745
    grids = build_quadrature_grids(N_xi=10, N_eta=8, N_phi=8, xi_max=10.0)

    print(f"{'Basis':>20s} {'p_max':>5s} {'N_bf':>5s} {'E_total':>12s} "
          f"{'D_e':>8s} {'%exact':>8s}")
    print("-" * 65)

    for label, j_max, l_max in [("j1,l0", 1, 0), ("j1,l1", 1, 1)]:
        for p_max in [0, 1]:
            basis = generate_basis(j_max=j_max, l_max=l_max, p_max=p_max,
                                   alpha=0.9)
            r = solve_hylleraas(basis, R, grids, verbose=False)
            pct = r['D_e'] / D_e_exact * 100
            print(f"{label:>20s} {p_max:5d} {len(basis):5d} "
                  f"{r['E_total']:12.6f} {r['D_e']:8.4f} {pct:7.1f}%")
        print()

    print()


def main():
    t_start = time.time()

    print()
    print("*" * 70)
    print("  HYLLERAAS CONVERGENCE STUDY — H2 MOLECULE")
    print("  Phase 9 of Prolate Spheroidal Stress Tests")
    print("*" * 70)
    print(f"  Reference: D_e(exact) = 0.1745 Ha, E(exact) = -1.17475 Ha")
    print(f"  James-Coolidge (1933): E = -1.17447 Ha (13 terms)")
    print()

    # Phase 9a: p=0 (fast)
    results_p0, best_alpha = run_p0_convergence()

    # Phase 9c: comparison (medium speed)
    run_comparison()

    # Phase 9b: p=1 (slow)
    results_p1 = run_p1_convergence(best_alpha)

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    all_results = results_p0 + results_p1
    print(f"{'p_max':>5s} {'j_max':>5s} {'l_max':>5s} {'N_bf':>5s} "
          f"{'D_e':>8s} {'%exact':>8s}")
    print("-" * 45)
    for r in all_results:
        print(f"{r['p_max']:5d} {r['j_max']:5d} {r['l_max']:5d} "
              f"{r['n_basis']:5d} {r['D_e']:8.4f} {r['D_e_pct']:7.1f}%")

    dt_total = time.time() - t_start
    print(f"\nTotal time: {dt_total:.0f}s ({dt_total/60:.1f} min)")

    # Save results
    with open('debug/data/hylleraas_convergence.txt', 'w') as f:
        f.write("# Hylleraas Convergence Study\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"# R = 1.4 bohr\n")
        f.write(f"# D_e(exact) = 0.1745 Ha\n\n")
        f.write(f"{'p_max':>5s} {'j_max':>5s} {'l_max':>5s} {'N_bf':>5s} "
                f"{'E_total':>12s} {'D_e':>10s} {'pct_exact':>10s} "
                f"{'time_s':>8s} {'alpha':>6s}\n")
        for r in all_results:
            f.write(f"{r['p_max']:5d} {r['j_max']:5d} {r['l_max']:5d} "
                    f"{r['n_basis']:5d} {r['E_total']:12.6f} "
                    f"{r['D_e']:10.6f} {r['D_e_pct']:10.1f} "
                    f"{r['time']:8.1f} {r['alpha']:6.2f}\n")
    print("\nResults saved to debug/data/hylleraas_convergence.txt")


if __name__ == '__main__':
    main()
