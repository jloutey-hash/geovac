"""G4-6a isotropic refinement: scale BOTH N_rho and N_phi as 1/a.

The single-axis radial sweep (g4_6a_algebraic_richardson.py) showed that
tip(t) converges to B (not A/t + B) because N_phi is fixed. The A coefficient
requires azimuthal UV content.

This diagnostic tests ISOTROPIC refinement: N_rho = R/a, N_phi = C/a, so
both resolutions improve together. This is the physical continuum limit.

If tip(t) under isotropic refinement converges to A/t + B at intermediate t,
then Richardson extrapolation in the isotropic parameter gives A cleanly.
"""

import numpy as np
import json
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral, DiscreteWedgeDiracSpectral
)


def compute_tip(a, N_rho, N_phi, t_values, alpha_eps=0.1):
    """Compute tip(t) = dK/dalpha|_{alpha=1} - K_disk."""
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


def main():
    print("=" * 72)
    print("G4-6a ISOTROPIC: Richardson with joint (N_rho, N_phi) refinement")
    print("=" * 72)

    A_cont = 1.0 / (24 * np.pi)
    B_cont = 1.0 / 6.0

    R = 10.0
    alpha_eps = 0.1

    # Isotropic refinement: N_phi = N_phi_0 / a (proportional to N_rho)
    # Base: a=0.1, N_rho=100, N_phi=100 (1 mode per radial site)
    substrates = [
        (0.200, 50, 50),
        (0.100, 100, 100),
        (0.050, 200, 200),
        (0.025, 400, 400),
    ]

    # Fixed t-values for A extraction
    t_test = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])

    print(f"\n  A_cont = {A_cont:.8f}")
    print(f"  B_cont = {B_cont:.8f}")
    print(f"  R = {R}, eps = {alpha_eps}")
    print(f"  ISOTROPIC substrates (a, N_rho, N_phi):")
    for a, Nr, Np in substrates:
        print(f"    a={a:.3f}, N_rho={Nr}, N_phi={Np}, dim={2*Nr*Np}")

    all_tips = {}
    a_values = []
    for a, N_rho, N_phi in substrates:
        print(f"\n  Computing a={a}, N_rho={N_rho}, N_phi={N_phi}...",
              end=" ", flush=True)
        tips = compute_tip(a, N_rho, N_phi, t_test, alpha_eps)
        all_tips[a] = tips
        a_values.append(a)
        print(f"tip(0.5)={tips[2]:.6f}, tip(1.0)={tips[3]:.6f}")

    a_arr = np.array(a_values)

    # Convergence diagnostic
    print("\n" + "=" * 72)
    print("  CONVERGENCE: tip(t) vs continuum under isotropic refinement")
    print("=" * 72)
    print(f"\n  {'t':>6} | ", end="")
    for a in a_arr:
        print(f" a={a:.3f}", end="")
    print(f" | {'cont':>8} | {'fine/cont':>10}")
    print("  " + "-" * 80)

    for i, t in enumerate(t_test):
        tip_cont = A_cont / t + B_cont
        print(f"  {t:6.3f} | ", end="")
        for a in a_arr:
            print(f" {all_tips[a][i]:+.4f}", end="")
        recovery = all_tips[a_arr[-1]][i] / tip_cont if tip_cont != 0 else 0
        print(f" | {tip_cont:8.5f} | {recovery*100:8.1f}%")

    # A extraction with Richardson
    print("\n" + "=" * 72)
    print("  A EXTRACTION: A_est = t * (tip(t) - B_cont)")
    print("=" * 72)
    print(f"\n  {'t':>6} | ", end="")
    for a in a_arr:
        print(f" a={a:.3f}", end="")
    print(f" | {'Rich-1':>10} {'%':>6} | {'Rich-2':>10} {'%':>6}")
    print("  " + "-" * 90)

    for i, t in enumerate(t_test):
        A_ests = [t * (all_tips[a][i] - B_cont) for a in a_arr]

        # Richardson order 1 (linear convergence)
        v1, v2 = A_ests[-2], A_ests[-1]
        a1, a2 = a_arr[-2], a_arr[-1]
        r1 = (a1/a2)
        A_rich1 = (r1 * v2 - v1) / (r1 - 1)
        rec1 = A_rich1 / A_cont * 100

        # Richardson order 2 (quadratic convergence)
        r2 = (a1/a2)**2
        A_rich2 = (r2 * v2 - v1) / (r2 - 1)
        rec2 = A_rich2 / A_cont * 100

        print(f"  {t:6.3f} | ", end="")
        for A_e in A_ests:
            print(f" {A_e:+.4f}", end="")
        print(f" | {A_rich1:+.7f} {rec1:+.1f}% | {A_rich2:+.7f} {rec2:+.1f}%")

    # Convergence ORDER diagnostic
    print("\n" + "=" * 72)
    print("  CONVERGENCE ORDER: |tip(t) - tip_cont| vs a")
    print("=" * 72)
    print(f"\n  {'t':>6} | {'log(err) ratios (consecutive pairs)':>50} | {'est. order':>10}")
    print("  " + "-" * 80)

    for i, t in enumerate(t_test):
        tip_cont = A_cont / t + B_cont
        errs = [abs(all_tips[a][i] - tip_cont) for a in a_arr]
        ratios = []
        for j in range(1, len(errs)):
            if errs[j] > 0 and errs[j-1] > 0:
                r = np.log(errs[j-1]/errs[j]) / np.log(a_arr[j-1]/a_arr[j])
                ratios.append(r)
            else:
                ratios.append(float('nan'))
        print(f"  {t:6.3f} | ", end="")
        for r in ratios:
            print(f" {r:8.3f}", end="")
        mean_order = np.nanmean(ratios) if ratios else 0
        print(f" | {mean_order:8.3f}")

    # Save
    output = {
        'method': 'isotropic_refinement',
        'parameters': {'R': R, 'alpha_eps': alpha_eps},
        'substrates': [(a, Nr, Np) for a, Nr, Np in substrates],
        'A_continuum': A_cont,
        'B_continuum': B_cont,
        't_values': [float(t) for t in t_test],
        'tips_by_a': {str(a): [float(x) for x in v] for a, v in all_tips.items()},
    }
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/g4_6a_isotropic_richardson.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Saved to debug/data/g4_6a_isotropic_richardson.json")


if __name__ == '__main__':
    main()
