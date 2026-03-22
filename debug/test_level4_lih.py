"""
Level 4 LiH diagnostic: locked Li 1s² core in hyperspherical coordinates.

Compares three core models:
  'bare'    — Z_eff = Z (no screening), diagnostic baseline
  'zeff'    — Z_eff = Z - n_core (integer screening)
  'hartree' — Z_eff + Core Hartree penetration correction

Success criteria:
  - R_eq shifts toward 3.015 bohr (expt)
  - D_e improves toward 2.515 eV (expt)
  - PES has a proper minimum (bound state)

Date: 2026-03-20
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.level4_multichannel import solve_level4_lih


R_eq_expt = 3.015   # bohr
D_e_expt_eV = 2.515  # eV
D_e_expt_Ha = D_e_expt_eV / 27.2114


def test_single_point() -> None:
    """Quick single-point comparison of core models at R=3.015."""
    print("=" * 70)
    print("1. Single-point comparison at R = 3.015 bohr")
    print("=" * 70)

    results = {}
    for model in ['zeff', 'hartree']:
        print(f"\n--- core_model = '{model}' ---")
        try:
            r = solve_level4_lih(
                R=3.015,
                core_model=model,
                l_max=4,
                n_alpha=200,
                n_Re=400,
                verbose=True,
            )
            results[model] = r
        except Exception as e:
            print(f"  FAILED: {e}")
            results[model] = None

    print("\n" + "=" * 70)
    print("Summary:")
    print(f"{'Model':<10} {'E_total (Ha)':>14} {'D_e (eV)':>10} {'Bound?':>8}")
    print("-" * 50)
    for model in ['zeff', 'hartree']:
        r = results.get(model)
        if r is None:
            print(f"{model:<10} {'FAILED':>14}")
            continue
        bound = "YES" if r['D_e'] > 0 else "NO"
        print(f"{model:<10} {r['E_total']:>14.6f} {r['D_e_eV']:>10.3f} {bound:>8}")
    print(f"{'expt':<10} {'':>14} {D_e_expt_eV:>10.3f}")


def test_pes_scan() -> None:
    """PES scan with the zeff model (the safest baseline)."""
    print("\n" + "=" * 70)
    print("2. PES scan (zeff model)")
    print("=" * 70)

    R_values = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
    energies = []

    for R in R_values:
        print(f"\n--- R = {R:.1f} bohr ---")
        try:
            r = solve_level4_lih(
                R=R,
                core_model='zeff',
                l_max=4,
                n_alpha=150,
                n_Re=300,
                verbose=True,
            )
            energies.append(r['E_total'])
        except Exception as e:
            print(f"  FAILED: {e}")
            energies.append(np.nan)

    print("\n" + "=" * 70)
    print("PES Summary (zeff):")
    print(f"{'R (bohr)':>10} {'E_total (Ha)':>14}")
    print("-" * 28)
    for R, E in zip(R_values, energies):
        print(f"{R:>10.2f} {E:>14.6f}")

    # Find minimum
    valid = [(R, E) for R, E in zip(R_values, energies) if not np.isnan(E)]
    if valid:
        R_min, E_min = min(valid, key=lambda x: x[1])
        E_inf = valid[-1][1]
        D_e_Ha = E_inf - E_min
        D_e_eV = D_e_Ha * 27.2114
        print(f"\nR_eq ~ {R_min:.1f} bohr (expt: {R_eq_expt})")
        print(f"D_e  ~ {D_e_eV:.3f} eV (expt: {D_e_expt_eV})")


def test_energy_decomposition() -> None:
    """Show energy decomposition at R_eq."""
    print("\n" + "=" * 70)
    print("3. Energy decomposition at R = 3.015 bohr (zeff)")
    print("=" * 70)

    r = solve_level4_lih(
        R=3.015,
        core_model='zeff',
        l_max=4,
        verbose=True,
    )

    print(f"\n  E_core_isolated = {r['E_core_isolated']:.6f} Ha")
    print(f"  V_cross_nuc     = {r['V_cross_nuc']:.6f} Ha")
    print(f"  E_locked        = {r['E_locked']:.6f} Ha")
    print(f"  E_elec(active)  = {r['E_elec']:.6f} Ha")
    print(f"  V_NN            = {r['V_NN']:.6f} Ha")
    print(f"  E_total         = {r['E_total']:.6f} Ha")
    print(f"  Z_eff_A         = {r['Z_eff_A']:.2f}")
    print(f"  Z_eff_B         = {r['Z_eff_B']:.2f}")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', choices=['single', 'pes', 'decomp', 'all'],
                        default='single')
    args = parser.parse_args()

    if args.test in ('single', 'all'):
        test_single_point()
    if args.test in ('decomp', 'all'):
        test_energy_decomposition()
    if args.test in ('pes', 'all'):
        test_pes_scan()
