"""
Resolution Limit Analysis - Basis Set Convergence Study
=========================================================

Tests whether the remaining ~2-3% error in isoelectronic series
is due to finite basis set size (max_n) rather than fundamental
physics issues.

Hypothesis: Error → 0 as max_n → ∞

Method: Split Scaling v1 with varying max_n

Author: GeoVac Development Team
Date: February 15, 2026
"""

import numpy as np
import sys
import os
import time

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac import MoleculeHamiltonian, CALIBRATED_KINETIC_SCALE

# Reference values (NIST)
E_LITHIUM_ION = -7.27991  # Li+ ground state (Ha)
E_BERYLLIUM_DICATION = -13.65556  # Be2+ ground state (Ha)


def test_lithium_convergence(max_n_values, verbose=True):
    """
    Test Li+ energy convergence with increasing basis size.

    Parameters:
    -----------
    max_n_values : list of int
        Basis sizes to test (e.g., [6, 8, 10, 12, 15])
    verbose : bool
        Print detailed results

    Returns:
    --------
    results : dict
        Keys: max_n, values: (energy, error_pct, time)
    """
    Z_target = 3
    Z_ref = 2

    # Split scaling factors
    kinetic_Z_scale = (Z_target / Z_ref)**2
    potential_Z_scale = Z_target / Z_ref

    results = {}

    if verbose:
        print("\n" + "="*80)
        print(" "*20 + "LITHIUM ION (Li+) - BASIS SET CONVERGENCE")
        print("="*80)
        print("\nMethod: Split Scaling v1 (T,V_nuc ~ Z², V_ee ~ Z)")
        print(f"Target: {E_LITHIUM_ION:.6f} Ha")
        print("\n" + "-"*80)
        print(f"{'max_n':<10} {'States':<10} {'Energy (Ha)':<15} {'Error (%)':<12} {'Time (s)':<10}")
        print("-"*80)

    for max_n in max_n_values:
        t0 = time.time()

        # Build with TARGET Z (correct V_ee scaling)
        base_kinetic_scale = CALIBRATED_KINETIC_SCALE
        scaled_kinetic = base_kinetic_scale * kinetic_Z_scale

        mol = MoleculeHamiltonian(
            nuclei=[(0.0, 0.0, 0.0)],
            nuclear_charges=[Z_target],
            max_n=max_n,
            kinetic_scale=scaled_kinetic
        )

        # Apply split scaling for isoelectronic series
        mol.apply_isoelectronic_scaling(Z_ref=Z_ref, Z_target=Z_target)
        mol._build_molecular_hamiltonian()

        # Optimize Z_eff
        result = mol.optimize_effective_charge(
            method='full_ci',
            n_points=15,
            z_range=(0.7, 1.0)
        )
        mol.set_effective_charges(result['z_eff_optimal'])

        # Compute energy
        energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        E = energies[0]

        t1 = time.time()

        # Calculate error
        error_pct = 100 * abs(E - E_LITHIUM_ION) / abs(E_LITHIUM_ION)
        n_states = mol.lattices[0].num_states

        results[max_n] = {
            'energy': E,
            'error_pct': error_pct,
            'time': t1 - t0,
            'n_states': n_states
        }

        if verbose:
            print(f"{max_n:<10} {n_states:<10} {E:<15.6f} {error_pct:<12.2f} {t1-t0:<10.2f}")

    if verbose:
        print("-"*80)
        print(f"\nConvergence trend:")
        errors = [results[n]['error_pct'] for n in max_n_values]
        if len(errors) >= 2:
            improvement = errors[0] - errors[-1]
            print(f"  Error reduction: {improvement:.2f} percentage points")
            print(f"  From {errors[0]:.2f}% (max_n={max_n_values[0]}) → {errors[-1]:.2f}% (max_n={max_n_values[-1]})")

            # Extrapolate to infinite basis
            if len(max_n_values) >= 3:
                # Fit: E(n) = E_∞ + A/n^α
                # Simple estimate: use last two points
                n1, n2 = max_n_values[-2], max_n_values[-1]
                E1, E2 = results[n1]['energy'], results[n2]['energy']

                # Assume α=2 (typical for basis set convergence)
                A = (E1 - E2) / (1/n1**2 - 1/n2**2)
                E_inf = E2 - A/n2**2
                error_inf = 100 * abs(E_inf - E_LITHIUM_ION) / abs(E_LITHIUM_ION)

                print(f"\n  Extrapolated E(∞): {E_inf:.6f} Ha")
                print(f"  Extrapolated error: {error_inf:.2f}%")

        print("="*80)

    return results


def test_beryllium_convergence(max_n_values, verbose=True):
    """
    Test Be2+ energy convergence with increasing basis size.

    Parameters:
    -----------
    max_n_values : list of int
        Basis sizes to test (e.g., [6, 8, 10, 12, 15])
    verbose : bool
        Print detailed results

    Returns:
    --------
    results : dict
        Keys: max_n, values: (energy, error_pct, time)
    """
    Z_target = 4
    Z_ref = 2

    # Split scaling factors
    kinetic_Z_scale = (Z_target / Z_ref)**2
    potential_Z_scale = Z_target / Z_ref

    results = {}

    if verbose:
        print("\n" + "="*80)
        print(" "*15 + "BERYLLIUM DICATION (Be2+) - BASIS SET CONVERGENCE")
        print("="*80)
        print("\nMethod: Split Scaling v1 (T,V_nuc ~ Z², V_ee ~ Z)")
        print(f"Target: {E_BERYLLIUM_DICATION:.6f} Ha")
        print("\n" + "-"*80)
        print(f"{'max_n':<10} {'States':<10} {'Energy (Ha)':<15} {'Error (%)':<12} {'Time (s)':<10}")
        print("-"*80)

    for max_n in max_n_values:
        t0 = time.time()

        # Build with TARGET Z (correct V_ee scaling)
        base_kinetic_scale = CALIBRATED_KINETIC_SCALE
        scaled_kinetic = base_kinetic_scale * kinetic_Z_scale

        mol = MoleculeHamiltonian(
            nuclei=[(0.0, 0.0, 0.0)],
            nuclear_charges=[Z_target],
            max_n=max_n,
            kinetic_scale=scaled_kinetic
        )

        # Apply split scaling for isoelectronic series
        mol.apply_isoelectronic_scaling(Z_ref=Z_ref, Z_target=Z_target)
        mol._build_molecular_hamiltonian()

        # Optimize Z_eff
        result = mol.optimize_effective_charge(
            method='full_ci',
            n_points=15,
            z_range=(0.7, 1.0)
        )
        mol.set_effective_charges(result['z_eff_optimal'])

        # Compute energy
        energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        E = energies[0]

        t1 = time.time()

        # Calculate error
        error_pct = 100 * abs(E - E_BERYLLIUM_DICATION) / abs(E_BERYLLIUM_DICATION)
        n_states = mol.lattices[0].num_states

        results[max_n] = {
            'energy': E,
            'error_pct': error_pct,
            'time': t1 - t0,
            'n_states': n_states
        }

        if verbose:
            print(f"{max_n:<10} {n_states:<10} {E:<15.6f} {error_pct:<12.2f} {t1-t0:<10.2f}")

    if verbose:
        print("-"*80)
        print(f"\nConvergence trend:")
        errors = [results[n]['error_pct'] for n in max_n_values]
        if len(errors) >= 2:
            improvement = errors[0] - errors[-1]
            print(f"  Error reduction: {improvement:.2f} percentage points")
            print(f"  From {errors[0]:.2f}% (max_n={max_n_values[0]}) → {errors[-1]:.2f}% (max_n={max_n_values[-1]})")

            # Extrapolate to infinite basis
            if len(max_n_values) >= 3:
                n1, n2 = max_n_values[-2], max_n_values[-1]
                E1, E2 = results[n1]['energy'], results[n2]['energy']

                A = (E1 - E2) / (1/n1**2 - 1/n2**2)
                E_inf = E2 - A/n2**2
                error_inf = 100 * abs(E_inf - E_BERYLLIUM_DICATION) / abs(E_BERYLLIUM_DICATION)

                print(f"\n  Extrapolated E(∞): {E_inf:.6f} Ha")
                print(f"  Extrapolated error: {error_inf:.2f}%")

        print("="*80)

    return results


def main():
    """Run complete resolution limit analysis."""
    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + "RESOLUTION LIMIT ANALYSIS - v0.4.1".center(78) + "#")
    print("#" + "  Basis Set Convergence for Isoelectronic Series".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80)
    print("\nGoal: Determine if remaining error is due to finite basis size")
    print("Method: Split Scaling v1 with increasing max_n")

    # Test Li+ convergence
    max_n_values_li = [6, 8, 10, 12]
    results_li = test_lithium_convergence(max_n_values_li, verbose=True)

    # Test Be2+ convergence
    max_n_values_be = [6, 8, 10, 12]
    results_be = test_beryllium_convergence(max_n_values_be, verbose=True)

    # Final summary
    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + "SUMMARY".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80)

    print("\nLi+ Convergence:")
    for n in max_n_values_li:
        r = results_li[n]
        status = "OK" if r['error_pct'] < 1.0 else "  "
        print(f"  max_n={n:2d}: {r['error_pct']:5.2f}% error  {status}")

    print("\nBe2+ Convergence:")
    for n in max_n_values_be:
        r = results_be[n]
        status = "OK" if r['error_pct'] < 1.0 else "  "
        print(f"  max_n={n:2d}: {r['error_pct']:5.2f}% error  {status}")

    # Check if we reached <1% error
    li_best = min(results_li[n]['error_pct'] for n in max_n_values_li)
    be_best = min(results_be[n]['error_pct'] for n in max_n_values_be)

    print("\n" + "="*80)
    if li_best < 1.0 and be_best < 1.0:
        print("OK SUCCESS: <1% ERROR ACHIEVED!")
        print("The remaining error WAS due to finite basis size.")
        print("Split Scaling + High Resolution = Universal <1% Engine!")
    else:
        print(f"Best errors: Li+ {li_best:.2f}%, Be2+ {be_best:.2f}%")
        if li_best < 2.0 and be_best < 3.0:
            print("GOOD: Chemical accuracy achieved (~2-3% error)")
            print("Further improvement requires max_n > 12 or relativistic corrections")
        else:
            print("Further investigation needed for convergence")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
