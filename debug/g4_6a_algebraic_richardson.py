"""G4-6a algebraic approach: Richardson extrapolation for A coefficient.

Strategy: Since B = 1/6 is KNOWN exactly (verified at 0.001%), we can extract
A = 1/(24*pi) by measuring tip(t) at SMALL t where A/t dominates, subtracting
the known B, and Richardson-extrapolating the radial discretization error (a -> 0).

This is the SAME method that G4-5a used for B (measure at large t, Richardson
in a), but applied at small t for A.

The key question: does tip(t_fixed) converge under radial refinement at fixed
small t? If yes, Richardson extrapolation gives A without multi-axis sweeps.
"""

import numpy as np
import json
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral, DiscreteWedgeDiracSpectral
)


def compute_tip_at_t(a, N_rho, N_phi, t_values, alpha_eps=0.1):
    """Compute tip(t) = dK/dalpha|_{alpha=1} - K_disk for given substrate."""
    alpha_plus = 1.0 + alpha_eps
    alpha_minus = 1.0 - alpha_eps

    disk = DiscreteDiskDiracSpectral(N_rho, a, N_phi)
    wedge_plus = DiscreteWedgeDiracSpectral(N_rho, a, N_phi, alpha_plus)
    wedge_minus = DiscreteWedgeDiracSpectral(N_rho, a, N_phi, alpha_minus)

    tips = []
    for t in t_values:
        K_disk = disk.heat_trace(t)
        K_wp = wedge_plus.heat_trace(t)
        K_wm = wedge_minus.heat_trace(t)
        dK_dalpha = (K_wp - K_wm) / (2 * alpha_eps)
        tip = dK_dalpha - K_disk
        tips.append(tip)

    return np.array(tips)


def richardson_extrapolate(values, a_values, order=2):
    """Richardson extrapolation assuming error ~ a^order."""
    if len(values) < 2:
        return values[-1]
    # Use last two values
    v1, v2 = values[-2], values[-1]
    a1, a2 = a_values[-2], a_values[-1]
    ratio = (a1 / a2) ** order
    return (ratio * v2 - v1) / (ratio - 1)


def main():
    print("=" * 72)
    print("G4-6a ALGEBRAIC: Richardson extrapolation for A = 1/(24*pi)")
    print("=" * 72)

    A_continuum = 1.0 / (24 * np.pi)
    B_continuum = 1.0 / 6.0

    # Substrate parameters: radial refinement panel
    # Hold R = N_rho * a = 10 fixed
    R = 10.0
    N_phi = 120  # Fixed azimuthal count
    alpha_eps = 0.1

    substrates = [
        (0.100, 100),
        (0.050, 200),
        (0.025, 400),
        (0.0125, 800),
    ]

    # t-values spanning from UV (where A/t dominates) to intermediate
    # Choosing FIXED t-values (not scaled with a) so we can Richardson in a
    t_test = np.array([0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])

    print(f"\n  A_continuum = 1/(24*pi) = {A_continuum:.8f}")
    print(f"  B_continuum = 1/6 = {B_continuum:.8f}")
    print(f"  R = {R}, N_phi = {N_phi}, eps = {alpha_eps}")
    print(f"  Substrates: {[(a, N) for a, N in substrates]}")

    # Compute tip at each substrate and each t
    all_tips = {}
    a_values = []

    for a, N_rho in substrates:
        print(f"\n  Computing substrate a={a}, N_rho={N_rho}...", end=" ", flush=True)
        tips = compute_tip_at_t(a, N_rho, N_phi, t_test, alpha_eps)
        all_tips[a] = tips
        a_values.append(a)
        print(f"done. tip(0.1)={tips[2]:.6f}, tip(1.0)={tips[5]:.6f}")

    a_arr = np.array(a_values)

    # Extract A_est at each (a, t): A_est = t * (tip(t) - B_continuum)
    print("\n" + "=" * 72)
    print("  A EXTRACTION: A_est(a, t) = t * (tip(t) - 1/6)")
    print("=" * 72)
    print(f"\n  {'t':>8} | ", end="")
    for a in a_arr:
        print(f"  a={a:.4f}", end="")
    print(f" | {'Richardson':>12}  {'Recovery':>10}")
    print("  " + "-" * 90)

    results = {}
    for i, t in enumerate(t_test):
        A_ests = []
        for a in a_arr:
            tip_val = all_tips[a][i]
            A_est = t * (tip_val - B_continuum)
            A_ests.append(A_est)

        # Richardson extrapolation (order 1 and 2)
        A_rich_1 = richardson_extrapolate(A_ests, a_arr, order=1)
        A_rich_2 = richardson_extrapolate(A_ests, a_arr, order=2)

        recovery_1 = A_rich_1 / A_continuum * 100
        recovery_2 = A_rich_2 / A_continuum * 100

        print(f"  {t:8.4f} | ", end="")
        for A_e in A_ests:
            print(f"  {A_e:+.5f}", end="")
        print(f" | R1:{A_rich_1:+.6f}  {recovery_1:+.1f}%")
        print(f"  {'':>8} | {'':>{10*len(a_arr)}} | R2:{A_rich_2:+.6f}  {recovery_2:+.1f}%")

        results[float(t)] = {
            'A_ests': [float(x) for x in A_ests],
            'A_richardson_order1': float(A_rich_1),
            'A_richardson_order2': float(A_rich_2),
            'recovery_order1_pct': float(recovery_1),
            'recovery_order2_pct': float(recovery_2),
        }

    # Also check: convergence of tip(t) itself under a-refinement
    print("\n" + "=" * 72)
    print("  CONVERGENCE DIAGNOSTIC: tip(t) at fixed t, varying a")
    print("=" * 72)
    print(f"\n  {'t':>8} | ", end="")
    for a in a_arr:
        print(f"  a={a:.4f}", end="")
    print(f" | {'continuum':>10}  {'ratio fine/coarse':>18}")
    print("  " + "-" * 90)

    for i, t in enumerate(t_test):
        tip_cont = A_continuum / t + B_continuum
        print(f"  {t:8.4f} | ", end="")
        for a in a_arr:
            print(f"  {all_tips[a][i]:+.5f}", end="")
        # Ratio of finest to coarsest
        ratio = all_tips[a_arr[-1]][i] / all_tips[a_arr[0]][i] if all_tips[a_arr[0]][i] != 0 else 0
        print(f" | {tip_cont:10.5f}  {ratio:.4f}")

    # Summary
    print("\n" + "=" * 72)
    print("  SUMMARY")
    print("=" * 72)
    # Find best Richardson recovery
    best_t = None
    best_recovery = 0
    for t_val, res in results.items():
        r2 = abs(res['recovery_order2_pct'] - 100)
        if best_t is None or r2 < best_recovery:
            best_recovery = r2
            best_t = t_val
    print(f"\n  Best Richardson-2 recovery: t = {best_t}, "
          f"A_est = {results[best_t]['A_richardson_order2']:.8f}, "
          f"recovery = {results[best_t]['recovery_order2_pct']:.1f}%")
    print(f"  Target: A_cont = {A_continuum:.8f}")

    # Save
    output = {
        'parameters': {
            'R': R, 'N_phi': N_phi, 'alpha_eps': alpha_eps,
            'substrates': [(a, N) for a, N in substrates],
            'A_continuum': A_continuum,
            'B_continuum': B_continuum,
        },
        't_values': [float(t) for t in t_test],
        'tips_by_a': {str(a): [float(x) for x in tips] for a, tips in all_tips.items()},
        'A_extraction': results,
    }
    with open('debug/data/g4_6a_algebraic_richardson.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Data saved to debug/data/g4_6a_algebraic_richardson.json")


if __name__ == '__main__':
    main()
