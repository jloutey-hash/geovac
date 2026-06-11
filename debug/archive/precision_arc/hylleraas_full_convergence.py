"""
Full Hylleraas convergence study for H2.

Target: Reproduce James & Coolidge (1933) result of D_e = 0.1745 Ha (99.7%)
using the prolate spheroidal Hylleraas implementation with Numba acceleration.

Systematic expansion:
  1. Vary j_max from 0 to 4
  2. Vary l_max from 0 to 3
  3. Vary p_max from 0 to 4
  4. Optimize alpha at each level
  5. Record: N_terms, E_total, D_e, % of exact, wall time

Also uses generate_basis_truncated(omega_max) for comparison.

Key questions:
  - What's the smallest basis that achieves 95%? 99%?
  - How does optimal alpha vary with basis size?
  - Which quantum numbers (j,k,l,m,p) contribute most?
  - Is there a "diminishing returns" point?

Reference values:
  James & Coolidge (1933): E = -1.17447 Ha with 13 terms
  Kolos & Wolniewicz:      E = -1.17475 Ha (exact non-relativistic)
  D_e(exact) = 0.17475 Ha, we use 0.1745 as in the codebase
"""

import numpy as np
import time
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.hylleraas import (
    generate_basis,
    generate_basis_truncated,
    build_quadrature_grids,
    solve_hylleraas,
    optimize_alpha,
)

# Try to warm up Numba
try:
    from geovac._numba_kernels import NUMBA_AVAILABLE, warmup_hylleraas_jit
    if NUMBA_AVAILABLE:
        print("Warming up Numba JIT kernels...")
        t0 = time.time()
        warmup_hylleraas_jit()
        print(f"  JIT warmup: {time.time() - t0:.1f}s")
    else:
        print("WARNING: Numba not available, using pure Python (will be slow)")
except ImportError:
    NUMBA_AVAILABLE = False
    print("WARNING: Numba not available")


# ============================================================
# Configuration
# ============================================================

R_EQ = 1.4011  # Equilibrium bond length (bohr)
D_E_EXACT = 0.1745  # Exact D_e (Ha) for normalization
E_EXACT = -1.17475  # Exact total energy (Ha)
ALPHA_RANGE = (0.7, 1.3)  # Narrowed from (0.5, 2.5) based on initial runs


def get_grids(has_r12: bool, n_basis: int) -> dict:
    """Choose quadrature grids based on basis complexity.

    Larger bases need finer grids for numerical stability,
    but we balance accuracy vs. compute cost.
    """
    if not has_r12:
        # p=0: 4D integral, relatively cheap
        if n_basis <= 5:
            return build_quadrature_grids(N_xi=20, N_eta=16, N_phi=8, xi_max=12.0)
        elif n_basis <= 15:
            return build_quadrature_grids(N_xi=24, N_eta=18, N_phi=8, xi_max=12.0)
        else:
            return build_quadrature_grids(N_xi=28, N_eta=20, N_phi=8, xi_max=14.0)
    else:
        # p>0: 5D integral with FD kinetic, expensive
        if n_basis <= 5:
            return build_quadrature_grids(N_xi=16, N_eta=12, N_phi=10, xi_max=12.0)
        elif n_basis <= 15:
            return build_quadrature_grids(N_xi=18, N_eta=14, N_phi=12, xi_max=12.0)
        elif n_basis <= 30:
            return build_quadrature_grids(N_xi=20, N_eta=14, N_phi=12, xi_max=12.0)
        else:
            return build_quadrature_grids(N_xi=22, N_eta=16, N_phi=14, xi_max=14.0)


def run_single_point(
    j_max: int, l_max: int, p_max: int, label: str = "",
    use_truncated: bool = False, omega_max: int = 0,
) -> dict:
    """Run a single Hylleraas calculation with alpha optimization.

    Returns dict with all results including timing.
    """
    # Generate basis
    if use_truncated:
        def gen(alpha):
            return generate_basis_truncated(omega_max, p_max, alpha)
        basis_test = gen(1.0)
        tag = f"omega={omega_max},p={p_max}"
    else:
        def gen(alpha):
            return generate_basis(j_max, l_max, p_max, alpha)
        basis_test = gen(1.0)
        tag = f"j={j_max},l={l_max},p={p_max}"

    n_basis = len(basis_test)
    has_r12 = p_max > 0
    grids = get_grids(has_r12, n_basis)

    print(f"\n{'='*60}")
    print(f"  {label}: {tag}  ({n_basis} terms)")
    print(f"  Grid: {grids['N_xi']}x{grids['N_eta']}x{grids['N_phi']}")
    print(f"{'='*60}")

    t_start = time.time()

    try:
        opt_result = optimize_alpha(
            gen, R_EQ, grids,
            alpha_range=ALPHA_RANGE,
            verbose=True,
        )
        result = opt_result['result']
        alpha_opt = opt_result['alpha_opt']
    except Exception as e:
        print(f"  ERROR: {e}")
        return {
            'label': label, 'tag': tag, 'n_basis': n_basis,
            'j_max': j_max, 'l_max': l_max, 'p_max': p_max,
            'E_total': np.nan, 'D_e': np.nan, 'D_e_pct': np.nan,
            'alpha_opt': np.nan, 'wall_time': time.time() - t_start,
            'error': str(e),
        }

    wall_time = time.time() - t_start

    # Extract coefficient analysis
    coeffs = result['coeffs']
    basis_final = result['basis']
    coeff_info = []
    for i, (bf, c) in enumerate(zip(basis_final, coeffs)):
        coeff_info.append({
            'idx': i, 'j': bf.j, 'k': bf.k, 'l': bf.l, 'm': bf.m,
            'p': bf.p, 'coeff': c, 'abs_coeff': abs(c),
        })
    coeff_info.sort(key=lambda x: -x['abs_coeff'])

    print(f"\n  === RESULT ===")
    print(f"  E_total = {result['E_total']:.8f} Ha")
    print(f"  D_e     = {result['D_e']:.6f} Ha  ({result['D_e_pct']:.1f}% of exact)")
    print(f"  alpha   = {alpha_opt:.4f}")
    print(f"  Time    = {wall_time:.1f}s")
    print(f"\n  Top 5 coefficients:")
    for ci in coeff_info[:5]:
        print(f"    [{ci['idx']}] (j,k,l,m,p)=({ci['j']},{ci['k']},{ci['l']},"
              f"{ci['m']},{ci['p']})  c={ci['coeff']:+.6f}")

    return {
        'label': label, 'tag': tag, 'n_basis': n_basis,
        'j_max': j_max, 'l_max': l_max, 'p_max': p_max,
        'E_total': result['E_total'],
        'D_e': result['D_e'],
        'D_e_pct': result['D_e_pct'],
        'alpha_opt': alpha_opt,
        'wall_time': wall_time,
        'coeffs': coeff_info,
        'S_cond': np.linalg.cond(result['S']),
    }


# ============================================================
# Main convergence study
# ============================================================

def main():
    print("=" * 70)
    print("  HYLLERAAS FULL CONVERGENCE STUDY FOR H2")
    print(f"  R = {R_EQ} bohr, D_e(exact) = {D_E_EXACT} Ha")
    print(f"  Numba: {'YES' if NUMBA_AVAILABLE else 'NO'}")
    print("=" * 70)

    all_results = []
    t_global = time.time()

    # --------------------------------------------------------
    # Part 1: Systematic expansion with generate_basis
    # --------------------------------------------------------
    print("\n" + "#" * 70)
    print("# PART 1: Systematic basis expansion (generate_basis)")
    print("#" * 70)

    # Ordered from smallest to largest — trimmed for tractability
    cases = [
        # (j_max, l_max, p_max, label)
        # --- p=0: CI limit (fast, 4D) ---
        (0, 0, 0, "A1"),   # 1 term: simplest
        (1, 0, 0, "A2"),   # 3 terms
        (2, 0, 0, "A3"),   # 6 terms
        (1, 1, 0, "A4"),   # 6 terms, with eta
        (2, 1, 0, "A5"),   # 12 terms
        (2, 2, 0, "A6"),   # 27 terms
        # --- p=1: first r12 terms (5D, slower) ---
        (1, 0, 1, "B1"),   # minimal + r12
        (1, 1, 1, "B2"),   # + eta + r12
        (2, 1, 1, "B3"),   # medium + r12
        # --- p=2: more correlation ---
        (1, 0, 2, "C1"),   # small + r12^2
        (1, 1, 2, "C2"),   # + eta + r12^2
        (2, 1, 2, "C3"),   # medium + r12^2
    ]

    for j_max, l_max, p_max, label in cases:
        # Quick size check
        basis_test = generate_basis(j_max, l_max, p_max, alpha=1.0)
        n_test = len(basis_test)

        # Skip if too large (>40 terms with p>0 would take hours)
        if n_test > 40 and p_max > 0:
            print(f"\n  Skipping {label}: {n_test} terms with p>0 (too slow)")
            continue

        result = run_single_point(j_max, l_max, p_max, label=label)
        all_results.append(result)

        # Early termination if we hit 99%
        if result['D_e_pct'] >= 99.0:
            print(f"\n  *** 99% D_e ACHIEVED with {label}! ***")

    # --------------------------------------------------------
    # Part 2: omega_max truncation (James-Coolidge style)
    # --------------------------------------------------------
    print("\n" + "#" * 70)
    print("# PART 2: Total-power truncation (generate_basis_truncated)")
    print("#" * 70)

    omega_cases = [
        # (omega_max, p_max, label)
        (2, 1, "T1"),
        (3, 2, "T2"),
        (4, 2, "T3"),
        (5, 2, "T5"),
    ]

    for omega_max, p_max, label in omega_cases:
        basis_test = generate_basis_truncated(omega_max, p_max, alpha=1.0)
        n_test = len(basis_test)

        if n_test > 80 and p_max > 0:
            print(f"\n  Skipping {label}: omega={omega_max},p={p_max} "
                  f"({n_test} terms, too slow)")
            continue

        result = run_single_point(
            0, 0, p_max, label=label,
            use_truncated=True, omega_max=omega_max,
        )
        all_results.append(result)

    # --------------------------------------------------------
    # Summary table
    # --------------------------------------------------------
    total_time = time.time() - t_global

    print("\n" + "=" * 90)
    print("  CONVERGENCE SUMMARY TABLE")
    print("=" * 90)
    print(f"{'Label':<6} {'Basis':<22} {'N':>4} {'E_total':>12} "
          f"{'D_e':>9} {'%exact':>7} {'alpha':>6} {'Time':>8} {'S_cond':>10}")
    print("-" * 90)

    for r in all_results:
        if np.isnan(r.get('E_total', np.nan)):
            print(f"{r['label']:<6} {r['tag']:<22} {r['n_basis']:>4} "
                  f"{'ERROR':>12} {'':>9} {'':>7} {'':>6} "
                  f"{r['wall_time']:>7.0f}s")
            continue

        s_cond = r.get('S_cond', np.nan)
        print(f"{r['label']:<6} {r['tag']:<22} {r['n_basis']:>4} "
              f"{r['E_total']:>12.8f} {r['D_e']:>9.6f} "
              f"{r['D_e_pct']:>6.1f}% "
              f"{r['alpha_opt']:>6.3f} {r['wall_time']:>7.0f}s "
              f"{s_cond:>10.1e}")

    print("-" * 90)
    print(f"  Total wall time: {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"  Exact: E = {E_EXACT:.5f} Ha, D_e = {D_E_EXACT:.4f} Ha")

    # --------------------------------------------------------
    # Alpha trajectory
    # --------------------------------------------------------
    print("\n" + "=" * 60)
    print("  ALPHA TRAJECTORY")
    print("=" * 60)
    for r in all_results:
        if not np.isnan(r.get('alpha_opt', np.nan)):
            print(f"  {r['label']:<6} N={r['n_basis']:>3}  "
                  f"alpha={r['alpha_opt']:.4f}  D_e%={r['D_e_pct']:.1f}")

    # --------------------------------------------------------
    # Best result analysis
    # --------------------------------------------------------
    valid = [r for r in all_results if not np.isnan(r.get('E_total', np.nan))]
    if valid:
        best = min(valid, key=lambda r: r['E_total'])
        print(f"\n{'='*60}")
        print(f"  BEST RESULT: {best['label']} ({best['tag']})")
        print(f"  E_total = {best['E_total']:.8f} Ha")
        print(f"  D_e     = {best['D_e']:.6f} Ha ({best['D_e_pct']:.1f}% of exact)")
        print(f"  alpha   = {best['alpha_opt']:.4f}")
        print(f"  N_terms = {best['n_basis']}")
        print(f"  Time    = {best['wall_time']:.0f}s")

        if 'coeffs' in best:
            print(f"\n  Coefficient analysis (top 10):")
            for ci in best['coeffs'][:10]:
                print(f"    (j,k,l,m,p)=({ci['j']},{ci['k']},{ci['l']},"
                      f"{ci['m']},{ci['p']})  c={ci['coeff']:+.8f}  "
                      f"|c|={ci['abs_coeff']:.8f}")

    # --------------------------------------------------------
    # Milestones
    # --------------------------------------------------------
    print(f"\n{'='*60}")
    print("  MILESTONES")
    print("=" * 60)
    for threshold in [50, 75, 85, 90, 95, 97, 99]:
        hits = [r for r in valid if r['D_e_pct'] >= threshold]
        if hits:
            first = min(hits, key=lambda r: r['n_basis'])
            print(f"  {threshold}% D_e: {first['label']} ({first['tag']}, "
                  f"N={first['n_basis']}, {first['wall_time']:.0f}s)")
        else:
            print(f"  {threshold}% D_e: NOT REACHED")

    # --------------------------------------------------------
    # Save results to file
    # --------------------------------------------------------
    outfile = os.path.join(os.path.dirname(__file__),
                           'data', 'hylleraas_convergence.txt')
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, 'w') as f:
        f.write("# Hylleraas convergence study for H2\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"# R = {R_EQ} bohr\n")
        f.write(f"# Numba: {NUMBA_AVAILABLE}\n\n")
        f.write(f"{'Label':<6} {'Basis':<22} {'N':>4} {'E_total':>14} "
                f"{'D_e':>10} {'pct':>7} {'alpha':>7} {'time_s':>8}\n")
        for r in all_results:
            if np.isnan(r.get('E_total', np.nan)):
                continue
            f.write(f"{r['label']:<6} {r['tag']:<22} {r['n_basis']:>4} "
                    f"{r['E_total']:>14.8f} {r['D_e']:>10.6f} "
                    f"{r['D_e_pct']:>6.1f}% {r['alpha_opt']:>7.4f} "
                    f"{r['wall_time']:>7.0f}\n")
    print(f"\n  Results saved to {outfile}")

    # --------------------------------------------------------
    # Convergence plot
    # --------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: D_e vs N_terms
        ax = axes[0]
        ns = [r['n_basis'] for r in valid]
        des = [r['D_e_pct'] for r in valid]
        labels_plot = [r['label'] for r in valid]

        # Color by p_max
        colors = []
        for r in valid:
            if r['p_max'] == 0:
                colors.append('blue')
            elif r['p_max'] == 1:
                colors.append('green')
            elif r['p_max'] == 2:
                colors.append('orange')
            else:
                colors.append('red')

        ax.scatter(ns, des, c=colors, s=80, zorder=5, edgecolors='k')
        for n, d, lab in zip(ns, des, labels_plot):
            ax.annotate(lab, (n, d), fontsize=7, ha='left',
                        xytext=(4, 4), textcoords='offset points')
        ax.axhline(100, color='gray', linestyle='--', alpha=0.5, label='Exact')
        ax.axhline(99, color='gray', linestyle=':', alpha=0.5, label='99%')
        ax.set_xlabel('Number of basis functions')
        ax.set_ylabel('D_e (% of exact)')
        ax.set_title('Hylleraas Convergence: D_e vs Basis Size')
        ax.legend(['Exact', '99%', 'p=0', 'p=1', 'p=2', 'p≥3'])
        ax.grid(True, alpha=0.3)

        # Right: alpha trajectory
        ax = axes[1]
        ax.scatter(ns, [r['alpha_opt'] for r in valid],
                   c=colors, s=80, zorder=5, edgecolors='k')
        for n, a, lab in zip(ns, [r['alpha_opt'] for r in valid], labels_plot):
            ax.annotate(lab, (n, a), fontsize=7, ha='left',
                        xytext=(4, 4), textcoords='offset points')
        ax.set_xlabel('Number of basis functions')
        ax.set_ylabel('Optimal alpha')
        ax.set_title('Optimal Alpha vs Basis Size')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plotfile = os.path.join(os.path.dirname(__file__),
                                'plots', 'hylleraas_convergence.png')
        os.makedirs(os.path.dirname(plotfile), exist_ok=True)
        plt.savefig(plotfile, dpi=150)
        print(f"  Plot saved to {plotfile}")
        plt.close()
    except ImportError:
        print("  matplotlib not available, skipping plot")

    print("\n  CONVERGENCE STUDY COMPLETE")
    return all_results


if __name__ == '__main__':
    main()
