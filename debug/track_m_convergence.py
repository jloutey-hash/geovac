"""Track M — Prolate CI l_max Convergence Study.

Tests whether Paper 12's prolate spheroidal CI for H₂ converges algebraically
with increasing angular basis (l_max) or saturates at a ceiling.

Key question for Paper 18: Is Level 4's transcendental μ(ρ,R) a coordinate
artifact (Scenario A) or fundamentally necessary (Scenario B)?

Strategy:
  - Fix j_max=3 (radial saturation demonstrated at j=2→3: only 0.2% change)
  - Scan l_max from 0 to highest feasible
  - Use Neumann V_ee (exact, algebraic — no quadrature error in V_ee)
  - Fix alpha=1.0 (near-optimal; full optimization adds 10× cost)
  - Also run alpha-optimized at key points for validation

Existing data (from debug/data/neumann_convergence.txt):
  j=0,l=0:  1 bf, 61.8%    j=2,l=2: 27 bf, 92.2%
  j=1,l=0:  3 bf, 66.1%    j=3,l=2: 46 bf, 92.4%
  j=1,l=1:  6 bf, 74.3%    j=3,l=3: 72 bf, 92.4%
  j=2,l=1: 12 bf, 75.9%

Date: 2026-03-29
"""

import numpy as np
import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.hylleraas import (
    generate_basis,
    build_quadrature_grids,
    solve_hylleraas,
    optimize_alpha,
)

# ============================================================
# Configuration
# ============================================================

R_EQ = 1.4011       # Equilibrium bond length (bohr)
D_E_EXACT = 0.1745  # Exact D_e (Ha)
E_EXACT = -1.17475  # Exact total energy (Ha)
J_MAX = 3           # Fixed radial truncation

# Grid for T+V_ne quadrature (V_ee is algebraic via Neumann)
GRIDS = build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)


def run_single(j_max: int, l_max: int, alpha: float = 1.0,
               grids: dict = None) -> dict:
    """Run a single Hylleraas solve with Neumann V_ee."""
    if grids is None:
        grids = GRIDS
    basis = generate_basis(j_max=j_max, l_max=l_max, p_max=0, alpha=alpha)
    n_bf = len(basis)

    t0 = time.time()
    result = solve_hylleraas(
        basis, R_EQ, grids, verbose=False,
        vee_method='neumann', l_max_neumann=20,
    )
    dt = time.time() - t0

    D_e = result['D_e']
    D_e_pct = result['D_e_pct']
    D_e_error = 100.0 - D_e_pct  # Error in %

    return {
        'j_max': j_max,
        'l_max': l_max,
        'n_bf': n_bf,
        'alpha': alpha,
        'E_total': result['E_total'],
        'D_e': D_e,
        'D_e_pct': D_e_pct,
        'D_e_error': D_e_error,
        'time': dt,
    }


def run_alpha_optimized(j_max: int, l_max: int, grids: dict = None) -> dict:
    """Run with alpha optimization for validation."""
    if grids is None:
        grids = GRIDS

    def gen(alpha):
        return generate_basis(j_max=j_max, l_max=l_max, p_max=0, alpha=alpha)

    n_bf = len(gen(1.0))
    t0 = time.time()
    opt = optimize_alpha(
        gen, R_EQ, grids, alpha_range=(0.7, 1.3),
        verbose=False, vee_method='neumann', l_max_neumann=20,
    )
    dt = time.time() - t0

    r = opt['result']
    return {
        'j_max': j_max,
        'l_max': l_max,
        'n_bf': n_bf,
        'alpha': opt['alpha_opt'],
        'E_total': r['E_total'],
        'D_e': r['D_e'],
        'D_e_pct': r['D_e_pct'],
        'D_e_error': 100.0 - r['D_e_pct'],
        'time': dt,
    }


# ============================================================
# Phase 1: l_max scan at fixed alpha=1.0, j_max=3
# ============================================================

def phase1_lmax_scan():
    """Scan l_max with fixed j_max=3, alpha=1.0."""
    print("=" * 70)
    print("PHASE 1: l_max SCAN (j_max=3, alpha=1.0, Neumann V_ee)")
    print(f"R = {R_EQ}, D_e_exact = {D_E_EXACT}, E_exact = {E_EXACT}")
    print("=" * 70)
    print()

    print(f"{'l_max':>5s} {'N_bf':>5s} {'E_total':>12s} {'D_e':>10s} "
          f"{'D_e%':>8s} {'Error%':>8s} {'Time':>8s}")
    print("-" * 65)

    results = []
    for l_max in range(8):  # 0 through 7
        basis = generate_basis(j_max=J_MAX, l_max=l_max, p_max=0, alpha=1.0)
        n_bf = len(basis)

        # Estimate time; skip if > 60 minutes based on scaling
        if results and results[-1]['time'] > 0:
            ratio = (n_bf / results[-1]['n_bf']) ** 2.5
            est_time = results[-1]['time'] * ratio
            if est_time > 3600:
                print(f"{l_max:5d} {n_bf:5d} {'SKIPPED':>12s} "
                      f"{'':>10s} {'':>8s} {'':>8s} "
                      f"~{est_time:.0f}s est")
                continue

        r = run_single(J_MAX, l_max)
        results.append(r)
        print(f"{r['l_max']:5d} {r['n_bf']:5d} {r['E_total']:12.6f} "
              f"{r['D_e']:10.6f} {r['D_e_pct']:7.1f}% "
              f"{r['D_e_error']:7.2f}% {r['time']:7.1f}s")
        sys.stdout.flush()

    return results


# ============================================================
# Phase 2: Alpha-optimized validation at key l_max values
# ============================================================

def phase2_alpha_validation(phase1_results):
    """Run alpha-optimized at key points to check if alpha=1.0 is adequate."""
    print()
    print("=" * 70)
    print("PHASE 2: ALPHA OPTIMIZATION AT KEY l_max VALUES")
    print("=" * 70)
    print()

    # Only validate at l_max values where phase 1 ran
    l_max_vals = [r['l_max'] for r in phase1_results]
    # Pick a few key values
    key_vals = [v for v in [0, 2, 3, 4] if v in l_max_vals]

    print(f"{'l_max':>5s} {'N_bf':>5s} {'alpha':>7s} {'E_total':>12s} "
          f"{'D_e%':>8s} {'Error%':>8s} {'Time':>8s} {'vs_a=1':>8s}")
    print("-" * 70)

    results = []
    for l_max in key_vals:
        r = run_alpha_optimized(J_MAX, l_max)
        # Find the alpha=1.0 result for comparison
        r_fixed = next((x for x in phase1_results if x['l_max'] == l_max), None)
        delta = r['D_e_pct'] - r_fixed['D_e_pct'] if r_fixed else 0.0
        results.append(r)
        print(f"{r['l_max']:5d} {r['n_bf']:5d} {r['alpha']:7.3f} "
              f"{r['E_total']:12.6f} {r['D_e_pct']:7.1f}% "
              f"{r['D_e_error']:7.2f}% {r['time']:7.1f}s "
              f"{delta:+7.2f}%")
        sys.stdout.flush()

    return results


# ============================================================
# Phase 3: Convergence analysis
# ============================================================

def phase3_analysis(results):
    """Fit convergence model and determine scenario."""
    print()
    print("=" * 70)
    print("PHASE 3: CONVERGENCE ANALYSIS")
    print("=" * 70)
    print()

    # Extract l_max >= 1 data (l=0 is qualitatively different)
    data = [(r['l_max'], r['D_e_pct'], r['D_e_error'])
            for r in results if r['l_max'] >= 1]
    l_vals = np.array([d[0] for d in data])
    pct_vals = np.array([d[1] for d in data])
    err_vals = np.array([d[2] for d in data])

    # Check for convergence
    if len(data) < 3:
        print("Insufficient data for convergence analysis.")
        return

    # 1. Successive differences
    print("Successive D_e% improvements:")
    for i in range(1, len(data)):
        delta = pct_vals[i] - pct_vals[i-1]
        ratio = delta / (pct_vals[i-1] - pct_vals[i-2]) if i >= 2 and abs(pct_vals[i-1] - pct_vals[i-2]) > 1e-10 else float('nan')
        print(f"  l_max {l_vals[i-1]} -> {l_vals[i]}: "
              f"+{delta:.3f}% (ratio: {ratio:.3f})")

    # 2. Log-log fit for algebraic convergence: error ~ C * l^{-p}
    mask = err_vals > 0.01  # Only fit positive errors
    if np.sum(mask) >= 2:
        log_l = np.log(l_vals[mask])
        log_err = np.log(err_vals[mask])
        # Linear fit in log-log space
        coeffs = np.polyfit(log_l, log_err, 1)
        p_fit = -coeffs[0]  # Power law exponent
        C_fit = np.exp(coeffs[1])
        residuals = log_err - np.polyval(coeffs, log_l)
        r_squared = 1 - np.var(residuals) / np.var(log_err) if np.var(log_err) > 0 else 0

        print(f"\nPower-law fit: D_e_error ~ {C_fit:.1f} * l_max^(-{p_fit:.2f})")
        print(f"  R² = {r_squared:.4f}")
        print(f"  Extrapolated error at l_max=10: {C_fit * 10**(-p_fit):.2f}%")
        print(f"  Extrapolated error at l_max=20: {C_fit * 20**(-p_fit):.2f}%")
        print(f"  Extrapolated error at l_max=50: {C_fit * 50**(-p_fit):.2f}%")

    # 3. Check for saturation (consecutive changes < 0.05%)
    last_changes = [pct_vals[i] - pct_vals[i-1] for i in range(1, len(data))]
    if len(last_changes) >= 2:
        recent = last_changes[-2:]
        if all(abs(c) < 0.05 for c in recent):
            plateau_val = pct_vals[-1]
            print(f"\n*** SATURATION DETECTED ***")
            print(f"  D_e% has plateaued at {plateau_val:.1f}%")
            print(f"  Last two changes: {recent[0]:+.3f}%, {recent[1]:+.3f}%")
            print(f"  Remaining gap: {100 - plateau_val:.1f}% = basis completeness limit")
            print(f"\n  VERDICT: Scenario B — prolate CI saturates at ~{plateau_val:.0f}%.")
            print(f"  The mol-frame hyperspherical reparameterization (Paper 15)")
            print(f"  is fundamentally necessary to access the remaining {100 - plateau_val:.1f}%.")
        elif all(abs(c) < 0.5 for c in recent) and p_fit < 1 if 'p_fit' in dir() else False:
            print(f"\n  VERDICT: Scenario C — algebraic but very slow (p={p_fit:.2f} < 1).")
        else:
            print(f"\n  VERDICT: Scenario A — still converging.")

    # 4. Compare with Hylleraas partial-wave convergence rate
    print(f"\n  Reference: Kutzelnigg (1985) partial-wave convergence for 1/r₁₂ cusp:")
    print(f"  - Natural parity: O(l^(-4)) for energy increments")
    print(f"  - Unnatural parity: O(l^(-2)) for energy increments")
    if 'p_fit' in dir():
        print(f"  - This study: O(l^(-{p_fit:.1f})) — ", end="")
        if p_fit > 3:
            print("consistent with natural parity (fast convergence)")
        elif p_fit > 1:
            print("between natural and unnatural parity")
        else:
            print("slower than expected — suggests structural limitation")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("Track M — Prolate CI l_max Convergence Study")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print()

    # Phase 1: l_max scan
    results1 = phase1_lmax_scan()

    # Phase 2: alpha validation (only at fast l_max values)
    results2 = phase2_alpha_validation(results1)

    # Phase 3: convergence analysis
    phase3_analysis(results1)

    # Save results
    outfile = os.path.join(os.path.dirname(__file__), 'data',
                           'track_m_convergence.txt')
    with open(outfile, 'w') as f:
        f.write("Track M — Prolate CI l_max Convergence Study\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"R = {R_EQ}, j_max = {J_MAX}, alpha = 1.0\n")
        f.write(f"Neumann V_ee (l_max_neumann=20), grid 30x20\n\n")

        f.write(f"{'l_max':>5s} {'N_bf':>5s} {'E_total':>12s} {'D_e':>10s} "
                f"{'D_e%':>8s} {'Error%':>8s} {'Time':>8s}\n")
        f.write("-" * 65 + "\n")
        for r in results1:
            f.write(f"{r['l_max']:5d} {r['n_bf']:5d} {r['E_total']:12.6f} "
                    f"{r['D_e']:10.6f} {r['D_e_pct']:7.1f}% "
                    f"{r['D_e_error']:7.2f}% {r['time']:7.1f}s\n")

        if results2:
            f.write("\nAlpha-optimized validation:\n")
            for r in results2:
                f.write(f"  l={r['l_max']}: alpha={r['alpha']:.3f}, "
                        f"D_e={r['D_e_pct']:.1f}%, E={r['E_total']:.6f}\n")

    print(f"\nResults saved to {outfile}")
