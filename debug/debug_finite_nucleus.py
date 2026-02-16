"""
Finite Nucleus Experiment for Li+ (Z=3) Overbinding Correction
===============================================================

Hypothesis: The ~2% overbinding in Li+ (E=-7.42 vs target -7.28)
is caused by the point nucleus approximation.

The n=1 node weight W=-Z/1²=-3.0 is too strongly binding.
A finite nucleus softens this to W=-Z/R_eff, reducing overbinding.

This script sweeps R_eff to find the value that gives the correct
Li+ energy, then tests if it transfers to Be2+.

Date: February 15, 2026
"""

import numpy as np
import time
import sys
sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian, CALIBRATED_KINETIC_SCALE

# Reference values (NIST)
TARGETS = {
    'He':   -2.903724,
    'Li+':  -7.279913,
    'Be2+': -13.655566,
}


def run_lithium_sweep(r_eff_values: list, max_n: int = 7, verbose: bool = True) -> dict:
    """
    Sweep R_eff for Li+ with Split Scaling v1 + Finite Nucleus.

    Parameters:
    -----------
    r_eff_values : list of float
        Core radius values to test
    max_n : int
        Basis set size (default: 7)
    verbose : bool
        Print progress

    Returns:
    --------
    dict: {R_eff: {'energy', 'error_pct', 'time'}}
    """
    Z_target = 3
    Z_ref = 2
    E_target = TARGETS['Li+']

    kinetic_Z_scale = (Z_target / Z_ref)**2
    scaled_kinetic = CALIBRATED_KINETIC_SCALE * kinetic_Z_scale

    results = {}

    if verbose:
        print(f"\n{'='*70}")
        print(f"FINITE NUCLEUS SWEEP: Li+ (Z=3)")
        print(f"{'='*70}")
        print(f"  Target energy:  {E_target:.6f} Ha")
        print(f"  Kinetic scale:  {scaled_kinetic:.6f}")
        print(f"  max_n:          {max_n}")
        print(f"  R_eff values:   {r_eff_values}")
        print(f"{'='*70}\n")

    for r_eff in r_eff_values:
        t0 = time.time()

        # Build with finite nucleus
        mol = MoleculeHamiltonian(
            nuclei=[(0.0, 0.0, 0.0)],
            nuclear_charges=[Z_target],
            max_n=max_n,
            kinetic_scale=scaled_kinetic,
            nuclear_model='finite' if r_eff > 1.0 else 'point',
            core_radius=r_eff
        )

        # Apply split scaling
        mol.apply_isoelectronic_scaling(Z_ref=Z_ref, Z_target=Z_target)
        mol._build_molecular_hamiltonian()

        # Optimize Z_eff
        result = mol.optimize_effective_charge(
            method='full_ci', n_points=15, z_range=(0.7, 1.0)
        )
        mol.set_effective_charges(result['z_eff_optimal'])

        # Compute energy
        energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        t1 = time.time()

        E_computed = energies[0]
        error_pct = 100 * (E_computed - E_target) / abs(E_target)
        abs_error_pct = abs(error_pct)

        results[r_eff] = {
            'energy': E_computed,
            'error_pct': error_pct,
            'abs_error_pct': abs_error_pct,
            'z_eff': result['z_eff_optimal'][0],
            'time': t1 - t0
        }

        if verbose:
            sign = "overbinding" if error_pct < 0 else "underbinding"
            print(f"  R_eff={r_eff:.3f}  E={E_computed:.6f}  "
                  f"err={error_pct:+.3f}% ({sign})  "
                  f"Z_eff={result['z_eff_optimal'][0]:.3f}  "
                  f"t={t1-t0:.1f}s")

    return results


def run_beryllium_check(r_eff: float, max_n: int = 8, verbose: bool = True) -> dict:
    """
    Test optimal R_eff on Be2+ to check transferability.
    """
    Z_target = 4
    Z_ref = 2
    E_target = TARGETS['Be2+']

    kinetic_Z_scale = (Z_target / Z_ref)**2
    scaled_kinetic = CALIBRATED_KINETIC_SCALE * kinetic_Z_scale

    if verbose:
        print(f"\n{'='*70}")
        print(f"TRANSFERABILITY CHECK: Be2+ (Z=4) with R_eff={r_eff:.3f}")
        print(f"{'='*70}")
        print(f"  Target energy: {E_target:.6f} Ha")

    t0 = time.time()

    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0)],
        nuclear_charges=[Z_target],
        max_n=max_n,
        kinetic_scale=scaled_kinetic,
        nuclear_model='finite' if r_eff > 1.0 else 'point',
        core_radius=r_eff
    )

    mol.apply_isoelectronic_scaling(Z_ref=Z_ref, Z_target=Z_target)
    mol._build_molecular_hamiltonian()

    result = mol.optimize_effective_charge(
        method='full_ci', n_points=15, z_range=(0.7, 1.0)
    )
    mol.set_effective_charges(result['z_eff_optimal'])

    energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
    t1 = time.time()

    E_computed = energies[0]
    error_pct = 100 * (E_computed - E_target) / abs(E_target)

    if verbose:
        sign = "overbinding" if error_pct < 0 else "underbinding"
        print(f"\n  Energy:     {E_computed:.6f} Ha")
        print(f"  Target:     {E_target:.6f} Ha")
        print(f"  Error:      {error_pct:+.3f}% ({sign})")
        print(f"  Z_eff:      {result['z_eff_optimal'][0]:.3f}")
        print(f"  Time:       {t1-t0:.1f}s")

    return {
        'energy': E_computed,
        'error_pct': error_pct,
        'abs_error_pct': abs(error_pct),
        'z_eff': result['z_eff_optimal'][0],
        'time': t1 - t0
    }


def run_helium_check(r_eff: float, max_n: int = 5, verbose: bool = True) -> dict:
    """
    Verify R_eff doesn't break Helium (Z=2, the reference system).
    He uses no split scaling, so finite nucleus directly affects it.
    """
    E_target = TARGETS['He']

    if verbose:
        print(f"\n{'='*70}")
        print(f"REGRESSION CHECK: He (Z=2) with R_eff={r_eff:.3f}")
        print(f"{'='*70}")
        print(f"  Target energy: {E_target:.6f} Ha")

    t0 = time.time()

    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0)],
        nuclear_charges=[2],
        max_n=max_n,
        nuclear_model='finite' if r_eff > 1.0 else 'point',
        core_radius=r_eff
    )

    result = mol.optimize_effective_charge(
        method='full_ci', n_points=12, z_range=(0.7, 1.0)
    )
    mol.set_effective_charges(result['z_eff_optimal'])

    energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
    t1 = time.time()

    E_computed = energies[0]
    error_pct = 100 * (E_computed - E_target) / abs(E_target)

    if verbose:
        sign = "overbinding" if error_pct < 0 else "underbinding"
        print(f"\n  Energy:     {E_computed:.6f} Ha")
        print(f"  Target:     {E_target:.6f} Ha")
        print(f"  Error:      {error_pct:+.3f}% ({sign})")
        print(f"  Z_eff:      {result['z_eff_optimal'][0]:.3f}")
        print(f"  Time:       {t1-t0:.1f}s")

    return {
        'energy': E_computed,
        'error_pct': error_pct,
        'abs_error_pct': abs(error_pct),
        'z_eff': result['z_eff_optimal'][0],
        'time': t1 - t0
    }


if __name__ == '__main__':
    print("=" * 70)
    print("FINITE NUCLEUS EXPERIMENT")
    print("Hypothesis: Point nucleus singularity causes ~2% overbinding")
    print("Method: Core softening via R_eff parameter")
    print("=" * 70)

    # Phase 1: Coarse sweep on Li+
    print("\n\n--- PHASE 1: COARSE SWEEP (Li+) ---")
    r_eff_coarse = [1.0, 1.05, 1.10, 1.15, 1.20, 1.30, 1.50, 2.00]
    results_coarse = run_lithium_sweep(r_eff_coarse, max_n=7)

    # Find the R_eff closest to zero error
    best_r = min(results_coarse, key=lambda r: results_coarse[r]['abs_error_pct'])
    best_err = results_coarse[best_r]['abs_error_pct']
    print(f"\n  >>> Best coarse R_eff = {best_r:.3f} (|error| = {best_err:.3f}%)")

    # Phase 2: Fine sweep around best value
    print("\n\n--- PHASE 2: FINE SWEEP (Li+) ---")
    delta = 0.05
    r_eff_fine = np.linspace(max(1.0, best_r - delta), best_r + delta, 8).tolist()
    results_fine = run_lithium_sweep(r_eff_fine, max_n=7)

    best_r_fine = min(results_fine, key=lambda r: results_fine[r]['abs_error_pct'])
    best_err_fine = results_fine[best_r_fine]['abs_error_pct']
    print(f"\n  >>> Best fine R_eff = {best_r_fine:.4f} (|error| = {best_err_fine:.3f}%)")

    # Phase 3: Transferability check
    print("\n\n--- PHASE 3: TRANSFERABILITY ---")
    be_result = run_beryllium_check(best_r_fine, max_n=8)
    he_result = run_helium_check(best_r_fine, max_n=5)

    # Summary
    print("\n\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"\n  Optimal R_eff: {best_r_fine:.4f}")
    print(f"\n  System    | Point Nucleus   | Finite Nucleus (R_eff={best_r_fine:.3f})")
    print(f"  ---------+------------------+-------------------------------------------")

    # Point nucleus reference (R_eff=1.0)
    li_point = results_coarse.get(1.0, results_fine.get(1.0))
    print(f"  Li+      | {li_point['error_pct']:+.3f}%          | {results_fine[best_r_fine]['error_pct']:+.3f}%")
    print(f"  Be2+     | (from v0.4.1)     | {be_result['error_pct']:+.3f}%")
    print(f"  He       | (baseline)        | {he_result['error_pct']:+.3f}%")

    print(f"\n  Hypothesis {'CONFIRMED' if best_err_fine < 1.0 else 'PARTIAL'}: "
          f"Finite nucleus reduces Li+ error to {best_err_fine:.3f}%")

    if be_result['abs_error_pct'] < 2.0:
        print(f"  Transferability: CONFIRMED (Be2+ error = {be_result['abs_error_pct']:.3f}%)")
    else:
        print(f"  Transferability: PARTIAL (Be2+ error = {be_result['abs_error_pct']:.3f}%)")

    print(f"\n  He regression: {'PASS' if he_result['abs_error_pct'] < 5.0 else 'FAIL'} "
          f"({he_result['abs_error_pct']:.3f}%)")
