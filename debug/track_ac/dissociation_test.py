"""Track AC Part 2: Dissociation limit convergence.

Runs Level 4 solver at large R to check convergence to E_atoms = -1.0 Ha.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import json
import time
from geovac.level4_multichannel import solve_level4_h2_multichannel

E_ATOMS = -1.0  # Ha (exact: 2 × H ground state)

R_VALUES = [5.0, 10.0, 15.0, 20.0, 30.0, 50.0]


def run_dissociation(l_max=4, m_max=1):
    """Run PES at large R to check dissociation limit."""
    l_max_per_m = {0: l_max, 1: min(l_max, 2)} if m_max > 0 else None

    results = []
    for R in R_VALUES:
        t0 = time.time()
        res = solve_level4_h2_multichannel(
            R=R, l_max=l_max, m_max=m_max,
            l_max_per_m=l_max_per_m,
            n_coupled=1,  # adiabatic
            n_alpha=200, n_Re=400,
            angular_method='spectral', n_basis_angular=10,
            n_quad_angular=100, verbose=False
        )
        dt = time.time() - t0
        E = res['E_total']
        diff = E - E_ATOMS
        results.append({
            'R': float(R),
            'E_total': float(E),
            'E_minus_Eatoms': float(diff),
            'time': float(dt),
        })
        print(f"  R={R:5.1f}: E={E:.6f}, E-E_atoms={diff:+.6f} Ha ({dt:.1f}s)", flush=True)

    return results


if __name__ == '__main__':
    print("=" * 70)
    print("PART 2: Dissociation limit convergence (l_max=4, adiabatic)")
    print("=" * 70)

    print("\nSigma+pi (m_max=1, pi frozen at 2):")
    results_sp = run_dissociation(l_max=4, m_max=1)

    print("\nSigma-only (m_max=0):")
    results_s = run_dissociation(l_max=4, m_max=0)

    # Summary table
    print("\n" + "=" * 70)
    print("Dissociation Table (l_max=4, adiabatic)")
    print("=" * 70)
    print(f"{'R':>6} | {'E(s+p)':>10} | {'E-E_at':>10} | {'E(s)':>10} | {'E-E_at':>10}")
    print("-" * 70)
    for sp, s in zip(results_sp, results_s):
        print(f"{sp['R']:>6.1f} | {sp['E_total']:>10.6f} | "
              f"{sp['E_minus_Eatoms']:>+10.6f} | {s['E_total']:>10.6f} | "
              f"{s['E_minus_Eatoms']:>+10.6f}")

    # Also check l_max=2 for comparison
    print("\n\nSigma+pi l_max=2:")
    results_2 = run_dissociation(l_max=2, m_max=1)

    all_results = {
        'sigma_pi_lmax4': results_sp,
        'sigma_only_lmax4': results_s,
        'sigma_pi_lmax2': results_2,
    }

    with open('debug/track_ac/dissociation_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\nSaved to debug/track_ac/dissociation_results.json")
