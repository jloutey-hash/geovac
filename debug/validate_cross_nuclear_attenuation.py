"""
Cross-Nuclear Attenuation Scan
===============================
Scan attenuation functions and alpha parameters for H2 and LiH.

For each (function, alpha) pair, compute PES and extract R_eq, D_e, k.
Goal: find a universal (f, alpha) that works for both molecules.

Date: 2026-03-12
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
from typing import Dict, List, Tuple, Optional
import time

from geovac.lattice_index import MolecularLatticeIndex, LatticeIndex


# Reference data
H2_REF = {'R_eq': 1.401, 'D_e': 0.1745, 'k': 0.369}
LIH_REF = {'R_eq': 3.015, 'D_e': 0.092, 'k': 0.0659}


def compute_pes(
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    R_values: np.ndarray,
    attenuation: Optional[str] = None,
    alpha: float = 1.0,
    fci_method: str = 'auto',
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute PES E(R) for a given attenuation setting. Returns (R, E)."""
    energies = []
    for R in R_values:
        mol = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons,
            vee_method='slater_full',
            fci_method=fci_method,
            cross_nuclear_attenuation=attenuation,
            cross_nuclear_alpha=alpha,
        )
        E = mol.compute_ground_state(n_states=1)[0][0]
        energies.append(E)
    return R_values, np.array(energies)


def fit_pes(R: np.ndarray, E: np.ndarray, E_sep: float) -> Dict[str, float]:
    """Extract R_eq, D_e, k from a PES."""
    idx_min = np.argmin(E)

    # Check if minimum is at boundary (no true minimum)
    if idx_min == 0 or idx_min == len(R) - 1:
        # Try parabolic fit anyway but flag it
        lo = max(0, idx_min - 2)
        hi = min(len(R), idx_min + 3)
        if hi - lo < 3:
            return {'R_eq': R[idx_min], 'D_e': E_sep - E[idx_min],
                    'k': float('nan'), 'boundary': True}

    # Fit parabola near minimum
    lo = max(0, idx_min - 2)
    hi = min(len(R), idx_min + 3)
    if hi - lo < 3:
        lo = max(0, hi - 3)

    coeffs = np.polyfit(R[lo:hi], E[lo:hi], 2)
    a, b, c = coeffs
    R_eq = -b / (2 * a)
    E_min = a * R_eq**2 + b * R_eq + c
    k = 2 * a

    return {
        'R_eq': R_eq,
        'D_e': E_sep - E_min,
        'k': k,
        'boundary': (idx_min == 0 or idx_min == len(R) - 1),
    }


def scan_attenuation(
    name: str,
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    n_electrons_A: int,
    n_electrons_B: int,
    R_values: np.ndarray,
    ref: Dict[str, float],
    fci_method: str = 'auto',
) -> List[Dict]:
    """Scan all attenuation functions and alpha values."""

    print(f"\n{'='*80}")
    print(f"  {name} Cross-Nuclear Attenuation Scan")
    print(f"{'='*80}")

    # Isolated atom energies
    atom_A = LatticeIndex(
        n_electrons=n_electrons_A, max_n=nmax_A,
        nuclear_charge=Z_A, vee_method='slater_full',
        h1_method='exact', fci_method='auto',
    )
    E_A = atom_A.compute_ground_state(n_states=1)[0][0]

    atom_B = LatticeIndex(
        n_electrons=n_electrons_B, max_n=nmax_B,
        nuclear_charge=Z_B, vee_method='slater_full',
        h1_method='exact', fci_method='auto',
    )
    E_B = atom_B.compute_ground_state(n_states=1)[0][0]
    E_sep = E_A + E_B
    print(f"  E_sep = {E_sep:.6f} Ha")

    # Attenuation configs to scan
    configs = [
        (None, 0.0, 'bare'),
        ('linear', 0.5, 'linear(0.5)'),
        ('linear', 1.0, 'linear(1.0)'),
        ('linear', 2.0, 'linear(2.0)'),
        ('linear', 5.0, 'linear(5.0)'),
        ('exp', 0.5, 'exp(0.5)'),
        ('exp', 1.0, 'exp(1.0)'),
        ('exp', 2.0, 'exp(2.0)'),
        ('exp', 5.0, 'exp(5.0)'),
        ('pade', 0.5, 'pade(0.5)'),
        ('pade', 1.0, 'pade(1.0)'),
        ('pade', 2.0, 'pade(2.0)'),
        ('pade', 5.0, 'pade(5.0)'),
        ('quadratic', 0.0, 'quadratic'),
    ]

    print(f"\n  {'Config':>16s} {'R_eq':>8s} {'D_e':>8s} {'k':>8s} "
          f"{'R_err%':>8s} {'D_err%':>8s} {'k_err%':>8s} {'Notes':>10s}")
    print(f"  {'-'*76}")

    results = []
    t0 = time.time()

    for atten, alpha, label in configs:
        _, E_arr = compute_pes(
            Z_A, Z_B, nmax_A, nmax_B,
            n_electrons, R_values,
            attenuation=atten, alpha=alpha,
            fci_method=fci_method,
        )

        fit = fit_pes(R_values, E_arr, E_sep)

        R_err = abs(fit['R_eq'] - ref['R_eq']) / ref['R_eq'] * 100
        D_err = abs(fit['D_e'] - ref['D_e']) / ref['D_e'] * 100 if fit['D_e'] > 0 else float('inf')
        k_err = abs(fit['k'] - ref['k']) / ref['k'] * 100 if not np.isnan(fit['k']) else float('inf')

        notes = ''
        if fit['boundary']:
            notes = 'BOUNDARY'
        if fit['k'] < 0:
            notes = 'CONCAVE'
        if fit['D_e'] < 0:
            notes = 'UNBOUND'

        print(f"  {label:>16s} {fit['R_eq']:8.3f} {fit['D_e']:8.4f} {fit['k']:8.4f} "
              f"{R_err:8.1f} {D_err:8.1f} {k_err:8.1f} {notes:>10s}")

        results.append({
            'label': label,
            'attenuation': atten,
            'alpha': alpha,
            'R_eq': fit['R_eq'],
            'D_e': fit['D_e'],
            'k': fit['k'],
            'R_err': R_err,
            'D_err': D_err,
            'k_err': k_err,
            'boundary': fit['boundary'],
            'E_arr': E_arr,
        })

    elapsed = time.time() - t0
    print(f"\n  Scan completed in {elapsed:.1f}s")

    # Find best config
    valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
    if valid:
        best = min(valid, key=lambda r: r['R_err'])
        print(f"\n  Best (by R_eq): {best['label']} — R_eq={best['R_eq']:.3f} "
              f"({best['R_err']:.1f}%), D_e={best['D_e']:.4f} ({best['D_err']:.1f}%), "
              f"k={best['k']:.4f} ({best['k_err']:.1f}%)")
    else:
        print(f"\n  No valid minimum found for any configuration!")

    return results


def save_results(
    results_h2: List[Dict], results_lih: List[Dict],
    R_h2: np.ndarray, R_lih: np.ndarray,
) -> None:
    """Save results to data files."""
    for name, results, R_vals in [('h2', results_h2, R_h2), ('lih', results_lih, R_lih)]:
        filepath = f'debug/data/cross_nuclear_attenuation_{name}.txt'
        with open(filepath, 'w') as f:
            f.write(f"# {name.upper()} Cross-Nuclear Attenuation Results\n")
            f.write(f"# {'Config':>16s} {'R_eq':>8s} {'D_e':>8s} {'k':>8s} "
                    f"{'R_err%':>8s} {'D_err%':>8s} {'k_err%':>8s}\n")
            for r in results:
                f.write(f"  {r['label']:>16s} {r['R_eq']:8.3f} {r['D_e']:8.4f} "
                        f"{r['k']:8.4f} {r['R_err']:8.1f} {r['D_err']:8.1f} "
                        f"{r['k_err']:8.1f}\n")
            # Also save PES for each config
            f.write(f"\n# PES data: R  E(bare)  E(best)\n")
            bare = results[0]['E_arr']
            best_valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
            if best_valid:
                best = min(best_valid, key=lambda r: r['R_err'])
                for i, R in enumerate(R_vals):
                    f.write(f"  {R:8.4f} {bare[i]:14.8f} {best['E_arr'][i]:14.8f}\n")
        print(f"  Saved to {filepath}")


def write_report(
    results_h2: List[Dict], results_lih: List[Dict],
) -> None:
    """Write analysis report."""
    filepath = 'debug/CROSS_NUCLEAR_ATTENUATION.md'

    # Find configs that work for both
    h2_valid = {r['label']: r for r in results_h2
                if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    lih_valid = {r['label']: r for r in results_lih
                 if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}

    common = set(h2_valid.keys()) & set(lih_valid.keys())

    with open(filepath, 'w') as f:
        f.write("# Cross-Nuclear Attenuation Analysis\n\n")
        f.write("**Date:** 2026-03-12\n\n")

        f.write("## Problem\n")
        f.write("Cross-nuclear attraction (V_cross) is the dominant driver of PES shape errors.\n")
        f.write("It grows too strong at short R because the graph framework uses unperturbed\n")
        f.write("atomic densities with no orthogonalization penalty.\n\n")

        f.write("## Method\n")
        f.write("Attenuate cross-nuclear attraction: V_cross(a) *= f(S_a)\n")
        f.write("where S_a = sum_b S(a,b)^2 is total squared overlap.\n\n")

        f.write("## H2 Results (nmax=3, expt: R_eq=1.401, D_e=0.1745, k=0.369)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |\n")
        f.write(f"|--------|------|-----|---|--------|--------|--------|\n")
        for r in results_h2:
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f}% | {r['D_err']:.1f}% | "
                    f"{r['k_err']:.1f}% |\n")

        f.write(f"\n## LiH Results (nmax=3, expt: R_eq=3.015, D_e=0.092, k=0.0659)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |\n")
        f.write(f"|--------|------|-----|---|--------|--------|--------|\n")
        for r in results_lih:
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f}% | {r['D_err']:.1f}% | "
                    f"{r['k_err']:.1f}% |\n")

        f.write(f"\n## Universal Configs (valid minimum for both molecules)\n\n")
        if common:
            f.write(f"| Config | H2 R_err% | H2 D_err% | LiH R_err% | LiH D_err% |\n")
            f.write(f"|--------|-----------|-----------|------------|------------|\n")
            for label in sorted(common):
                h2 = h2_valid[label]
                lih = lih_valid[label]
                f.write(f"| {label} | {h2['R_err']:.1f}% | {h2['D_err']:.1f}% | "
                        f"{lih['R_err']:.1f}% | {lih['D_err']:.1f}% |\n")

            # Best universal
            best_label = min(common, key=lambda l: h2_valid[l]['R_err'] + lih_valid[l]['R_err'])
            f.write(f"\n**Best universal:** {best_label}\n")
            f.write(f"- H2: R_eq={h2_valid[best_label]['R_eq']:.3f} ({h2_valid[best_label]['R_err']:.1f}%), "
                    f"D_e={h2_valid[best_label]['D_e']:.4f} ({h2_valid[best_label]['D_err']:.1f}%)\n")
            f.write(f"- LiH: R_eq={lih_valid[best_label]['R_eq']:.3f} ({lih_valid[best_label]['R_err']:.1f}%), "
                    f"D_e={lih_valid[best_label]['D_e']:.4f} ({lih_valid[best_label]['D_err']:.1f}%)\n")
        else:
            f.write("No configuration produced valid minima for both molecules.\n")

        f.write(f"\n## Conclusion\n\n")
        if common:
            f.write("Overlap-attenuated cross-nuclear attraction can produce correct PES shapes.\n")
            f.write("The attenuation is physically motivated (Pauli orthogonalization screening).\n")
        else:
            f.write("No universal attenuation found. The cross-nuclear issue may require\n")
            f.write("a different approach (e.g., full Löwdin orthogonalization of the basis).\n")

    print(f"  Report saved to {filepath}")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("Cross-Nuclear Attenuation Scan")
    print("=" * 80)

    # --- H2 (fast: nmax=3, 1540 SDs per point) ---
    R_h2 = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])

    results_h2 = scan_attenuation(
        name='H2',
        Z_A=1, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=2,
        n_electrons_A=1,
        n_electrons_B=1,
        R_values=R_h2,
        ref=H2_REF,
    )

    # --- LiH (slow: nmax=3, 367k SDs per point) ---
    R_lih = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0])

    results_lih = scan_attenuation(
        name='LiH',
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=4,
        n_electrons_A=3,
        n_electrons_B=1,
        R_values=R_lih,
        ref=LIH_REF,
    )

    # --- Save and report ---
    save_results(results_h2, results_lih, R_h2, R_lih)
    write_report(results_h2, results_lih)

    # --- Final summary ---
    print(f"\n{'='*80}")
    print("  SUMMARY")
    print(f"{'='*80}")

    h2_valid = [r for r in results_h2 if not r['boundary'] and r['k'] > 0]
    lih_valid = [r for r in results_lih if not r['boundary'] and r['k'] > 0]

    h2_labels = {r['label'] for r in h2_valid}
    lih_labels = {r['label'] for r in lih_valid}
    universal = h2_labels & lih_labels

    if universal:
        print(f"\n  Universal configs (valid for both): {sorted(universal)}")
        for label in sorted(universal):
            h2 = next(r for r in results_h2 if r['label'] == label)
            lih = next(r for r in results_lih if r['label'] == label)
            score = h2['R_err'] + lih['R_err']
            print(f"    {label:>16s}: H2 R_err={h2['R_err']:.1f}%, "
                  f"LiH R_err={lih['R_err']:.1f}%, combined={score:.1f}%")
    else:
        print(f"\n  No universal config found.")
        if h2_valid:
            best_h2 = min(h2_valid, key=lambda r: r['R_err'])
            print(f"  Best H2: {best_h2['label']} (R_err={best_h2['R_err']:.1f}%)")
        if lih_valid:
            best_lih = min(lih_valid, key=lambda r: r['R_err'])
            print(f"  Best LiH: {best_lih['label']} (R_err={best_lih['R_err']:.1f}%)")
