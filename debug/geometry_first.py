"""
Geometry First: Lattice Torsion Experiment
==========================================

The nucleus is a topological defect. The connection between the
innermost electron shell (n=1) and the singularity has TORSION.

This torsion modifies the METRIC (adjacency/kinetic), not the
POTENTIAL. The potential W = -Z/n^2 remains pure topology.

    A_ij -> A_ij * (1 - gamma)    for edges touching n=1 core

Question: Can we correct the ~2% overbinding in Li+ by modifying
kinetic topology instead of potential energy?

Date: February 15, 2026
"""

import numpy as np
import time
import sys
import io

sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian, CALIBRATED_KINETIC_SCALE

# Reference values (NIST)
TARGETS = {
    'He':   -2.903724,
    'Li+':  -7.279913,
    'Be2+': -13.655566,
}


def _build_and_solve(Z_target: int, Z_ref: int, max_n: int,
                     torsion: float, use_split_scaling: bool = True) -> dict:
    """
    Build and solve an isoelectronic system with given torsion.

    Returns dict with energy, error_pct, z_eff, time.
    """
    if use_split_scaling:
        kinetic_Z_scale = (Z_target / Z_ref)**2
        scaled_kinetic = CALIBRATED_KINETIC_SCALE * kinetic_Z_scale
    else:
        scaled_kinetic = CALIBRATED_KINETIC_SCALE

    t0 = time.time()

    # Suppress build output
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()

    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0)],
        nuclear_charges=[Z_target],
        max_n=max_n,
        kinetic_scale=scaled_kinetic,
        lattice_torsion=torsion
    )

    if use_split_scaling:
        mol.apply_isoelectronic_scaling(Z_ref=Z_ref, Z_target=Z_target)
        mol._build_molecular_hamiltonian()

    # Optimize Z_eff
    result = mol.optimize_effective_charge(
        method='full_ci', n_points=15, z_range=(0.7, 1.0)
    )
    mol.set_effective_charges(result['z_eff_optimal'])

    # Compute energy
    energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')

    sys.stdout = old_stdout
    t1 = time.time()

    return {
        'energy': energies[0],
        'z_eff': result['z_eff_optimal'][0],
        'time': t1 - t0
    }


def sweep_torsion(system: str, torsion_values: list,
                  verbose: bool = True) -> dict:
    """
    Sweep lattice torsion for a given system.

    Parameters:
    -----------
    system : str
        'Li+', 'Be2+', or 'He'
    torsion_values : list of float
        Torsion values to test
    """
    config = {
        'Li+':  {'Z_target': 3, 'Z_ref': 2, 'max_n': 7,
                  'split_scaling': True},
        'Be2+': {'Z_target': 4, 'Z_ref': 2, 'max_n': 8,
                  'split_scaling': True},
        'He':   {'Z_target': 2, 'Z_ref': 2, 'max_n': 5,
                  'split_scaling': False},
    }

    if system not in config:
        raise ValueError(f"Unknown system: {system}")

    cfg = config[system]
    E_target = TARGETS[system]

    if verbose:
        print(f"\n{'='*70}")
        print(f"TORSION SWEEP: {system} (Z={cfg['Z_target']})")
        print(f"{'='*70}")
        print(f"  Target:     {E_target:.6f} Ha")
        print(f"  max_n:      {cfg['max_n']}")
        print(f"  Split scaling: {cfg['split_scaling']}")
        print(f"  Torsion values: {len(torsion_values)} points")
        print(f"{'='*70}")
        print(f"\n  {'gamma':>8}  {'Energy':>12}  {'Error':>10}  "
              f"{'Direction':>14}  {'Z_eff':>6}  {'Time':>6}")
        print(f"  {'-'*8}  {'-'*12}  {'-'*10}  {'-'*14}  {'-'*6}  {'-'*6}")

    results = {}
    for gamma in torsion_values:
        r = _build_and_solve(
            Z_target=cfg['Z_target'],
            Z_ref=cfg['Z_ref'],
            max_n=cfg['max_n'],
            torsion=gamma,
            use_split_scaling=cfg['split_scaling']
        )

        error_pct = 100 * (r['energy'] - E_target) / abs(E_target)
        direction = "overbinding" if error_pct < 0 else "underbinding"

        results[gamma] = {
            'energy': r['energy'],
            'error_pct': error_pct,
            'abs_error_pct': abs(error_pct),
            'z_eff': r['z_eff'],
            'time': r['time']
        }

        if verbose:
            print(f"  {gamma:8.5f}  {r['energy']:12.6f}  "
                  f"{error_pct:+9.3f}%  {direction:>14}  "
                  f"{r['z_eff']:.3f}  {r['time']:.1f}s")

    return results


def find_zero_crossing(results: dict) -> float:
    """
    Interpolate to find the torsion value where error = 0.
    """
    gammas = sorted(results.keys())
    errors = [results[g]['error_pct'] for g in gammas]

    # Find sign change
    for i in range(len(errors) - 1):
        if errors[i] * errors[i+1] <= 0:
            # Linear interpolation
            g0, g1 = gammas[i], gammas[i+1]
            e0, e1 = errors[i], errors[i+1]
            gamma_zero = g0 - e0 * (g1 - g0) / (e1 - e0)
            return gamma_zero

    # No crossing found - return best
    best_idx = min(range(len(errors)), key=lambda i: abs(errors[i]))
    return gammas[best_idx]


if __name__ == '__main__':
    print("=" * 70)
    print("GEOMETRY FIRST: LATTICE TORSION EXPERIMENT")
    print("=" * 70)
    print("\nHypothesis: The nuclear defect has torsion that deforms")
    print("the local metric, reducing hopping to the n=1 core.")
    print("This is a KINETIC correction, not a POTENTIAL correction.")
    print("The potential W = -Z/n^2 remains PURE TOPOLOGY.")

    # ========================================================
    # Phase 1: Coarse sweep on Li+
    # ========================================================
    print("\n\n" + "=" * 70)
    print("PHASE 1: COARSE SWEEP (Li+)")
    print("=" * 70)

    gamma_coarse = [0.0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.10]
    results_li_coarse = sweep_torsion('Li+', gamma_coarse)

    # Find zero crossing
    gamma_zero = find_zero_crossing(results_li_coarse)
    print(f"\n  >>> Interpolated zero-crossing: gamma = {gamma_zero:.6f}")

    # ========================================================
    # Phase 2: Fine sweep around zero crossing
    # ========================================================
    print("\n\n" + "=" * 70)
    print("PHASE 2: FINE SWEEP (Li+)")
    print("=" * 70)

    delta = 0.003
    gamma_fine = np.linspace(
        max(0.0, gamma_zero - delta),
        gamma_zero + delta,
        8
    ).tolist()
    results_li_fine = sweep_torsion('Li+', gamma_fine)

    gamma_best = min(results_li_fine,
                     key=lambda g: results_li_fine[g]['abs_error_pct'])
    best_err = results_li_fine[gamma_best]['abs_error_pct']
    print(f"\n  >>> Best torsion: gamma = {gamma_best:.6f}")
    print(f"  >>> Li+ error:   {results_li_fine[gamma_best]['error_pct']:+.4f}%")

    # ========================================================
    # Phase 3: Transferability
    # ========================================================
    print("\n\n" + "=" * 70)
    print(f"PHASE 3: TRANSFERABILITY (gamma = {gamma_best:.6f})")
    print("=" * 70)

    # Be2+ with same torsion
    results_be = sweep_torsion('Be2+', [0.0, gamma_best])

    # He with same torsion (regression check)
    results_he = sweep_torsion('He', [0.0, gamma_best])

    # ========================================================
    # Final Summary
    # ========================================================
    print("\n\n" + "=" * 70)
    print("FINAL SUMMARY: GEOMETRY FIRST")
    print("=" * 70)

    li_base = results_li_coarse[0.0]
    li_best = results_li_fine[gamma_best]
    be_base = results_be[0.0]
    be_best = results_be[gamma_best]
    he_base = results_he[0.0]
    he_best = results_he[gamma_best]

    print(f"\n  Optimal torsion: gamma = {gamma_best:.6f}")
    print(f"\n  {'System':<8} | {'gamma=0 (flat)':>16} | "
          f"{'gamma={0:.5f}'.format(gamma_best):>18} | {'Change':>10}")
    print(f"  {'-'*8}-+-{'-'*16}-+-{'-'*18}-+-{'-'*10}")
    print(f"  {'Li+':<8} | {li_base['error_pct']:>+13.3f}%   | "
          f"{li_best['error_pct']:>+15.3f}%   | "
          f"{li_best['abs_error_pct'] - li_base['abs_error_pct']:>+9.3f}")
    print(f"  {'Be2+':<8} | {be_base['error_pct']:>+13.3f}%   | "
          f"{be_best['error_pct']:>+15.3f}%   | "
          f"{be_best['abs_error_pct'] - be_base['abs_error_pct']:>+9.3f}")
    print(f"  {'He':<8} | {he_base['error_pct']:>+13.3f}%   | "
          f"{he_best['error_pct']:>+15.3f}%   | "
          f"{he_best['abs_error_pct'] - he_base['abs_error_pct']:>+9.3f}")

    print(f"\n  Key question: Is the torsion UNIVERSAL?")
    print(f"  - Li+  improvement: "
          f"{li_base['abs_error_pct']:.3f}% -> {li_best['abs_error_pct']:.3f}%")
    print(f"  - Be2+ improvement: "
          f"{be_base['abs_error_pct']:.3f}% -> {be_best['abs_error_pct']:.3f}%")
    print(f"  - He   regression:  "
          f"{he_base['abs_error_pct']:.3f}% -> {he_best['abs_error_pct']:.3f}%")

    if li_best['abs_error_pct'] < 1.0:
        print(f"\n  VERDICT: Lattice torsion CORRECTS the overbinding!")
        print(f"           The nuclear defect IS a metric deformation.")
    else:
        print(f"\n  VERDICT: Torsion helps but does not fully correct.")
        print(f"           May need Z-dependent torsion or higher-order terms.")

    print(f"\n{'='*70}")
    print(f"'The potential is pure. The metric has torsion.'")
    print(f"{'='*70}")
