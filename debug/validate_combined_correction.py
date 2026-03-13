"""
Combined Correction: Attenuation + Fock-Weighted Kinetic
=========================================================
Scan the 2D parameter space (a0, lambda) with:
  - pade attenuation, beta=2.0 (Fock-style n-dependent)
  - Fock-weighted kinetic correction (diagonal)

Attenuation shifts R_eq; kinetic tunes curvature (k).

Date: 2026-03-12
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
from typing import Dict, List, Tuple, Optional
import time
import sys

from geovac.lattice_index import MolecularLatticeIndex, LatticeIndex


# Reference data
H2_REF = {'R_eq': 1.401, 'D_e': 0.1745, 'k': 0.369}
LIH_REF = {'R_eq': 3.015, 'D_e': 0.092, 'k': 0.0659}


def compute_pes(
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    R_values: np.ndarray,
    a0: float,
    lam: float,
    fci_method: str = 'auto',
) -> np.ndarray:
    """Compute PES with combined attenuation + kinetic correction."""
    energies = []
    for R in R_values:
        mol = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons,
            vee_method='slater_full',
            fci_method=fci_method,
            cross_nuclear_attenuation='pade',
            cross_nuclear_alpha=a0,
            cross_nuclear_beta=2.0,
            t_corr_fock_weighted=True,
            t_corr_lambda=lam,
        )
        E = mol.compute_ground_state(n_states=1)[0][0]
        energies.append(E)
    return np.array(energies)


def fit_pes(R: np.ndarray, E: np.ndarray, E_sep: float) -> Dict[str, float]:
    """Extract R_eq, D_e, k from PES via parabolic fit near minimum."""
    idx_min = np.argmin(E)
    boundary = (idx_min == 0 or idx_min == len(R) - 1)

    lo = max(0, idx_min - 2)
    hi = min(len(R), idx_min + 3)
    if hi - lo < 3:
        lo = max(0, hi - 3)
    if hi - lo < 3:
        return {'R_eq': R[idx_min], 'D_e': E_sep - E[idx_min],
                'k': float('nan'), 'boundary': True}

    coeffs = np.polyfit(R[lo:hi], E[lo:hi], 2)
    a, b, c = coeffs
    R_eq = -b / (2 * a)
    E_min = a * R_eq**2 + b * R_eq + c
    k = 2 * a

    return {
        'R_eq': R_eq,
        'D_e': E_sep - E_min,
        'k': k,
        'boundary': boundary,
    }


def scan_molecule(
    name: str,
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    n_electrons_A: int,
    n_electrons_B: int,
    R_values: np.ndarray,
    ref: Dict[str, float],
    a0_values: List[float],
    lam_values: List[float],
    fci_method: str = 'auto',
) -> List[Dict]:
    """Scan 2D (a0, lambda) grid for one molecule."""

    print(f"\n{'='*80}")
    print(f"  {name} Combined Correction Scan (pade, beta=2, Fock-weighted kinetic)")
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

    print(f"\n  {'Config':>18s} {'R_eq':>7s} {'D_e':>7s} {'k':>7s} "
          f"{'R%':>6s} {'D%':>6s} {'k%':>6s} {'Notes':>8s}")
    print(f"  {'-'*68}")

    results = []
    t0 = time.time()

    for a0 in a0_values:
        for lam in lam_values:
            label = f"a0={a0:.1f},L={lam:+.2f}"

            E_arr = compute_pes(
                Z_A, Z_B, nmax_A, nmax_B,
                n_electrons, R_values,
                a0=a0, lam=lam,
                fci_method=fci_method,
            )

            fit = fit_pes(R_values, E_arr, E_sep)

            R_err = abs(fit['R_eq'] - ref['R_eq']) / ref['R_eq'] * 100
            D_err = abs(fit['D_e'] - ref['D_e']) / ref['D_e'] * 100 if fit['D_e'] > 0 else float('inf')
            k_err = abs(fit['k'] - ref['k']) / ref['k'] * 100 if not np.isnan(fit['k']) else float('inf')

            notes = ''
            if fit['boundary']:
                notes = 'BNDRY'
            if fit['k'] < 0:
                notes = 'CONCAV'
            if fit['D_e'] < 0:
                notes = 'UNBND'

            print(f"  {label:>18s} {fit['R_eq']:7.3f} {fit['D_e']:7.4f} {fit['k']:7.4f} "
                  f"{R_err:6.1f} {D_err:6.1f} {k_err:6.1f} {notes:>8s}")
            sys.stdout.flush()

            results.append({
                'label': label,
                'a0': a0,
                'lam': lam,
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

    # Best by combined score (R + k error) among valid configs
    valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
    if valid:
        best = min(valid, key=lambda r: r['R_err'] + r['k_err'])
        print(f"  Best (R+k score): {best['label']} — R_eq={best['R_eq']:.3f} "
              f"({best['R_err']:.1f}%), D_e={best['D_e']:.4f} ({best['D_err']:.1f}%), "
              f"k={best['k']:.4f} ({best['k_err']:.1f}%)")

        # Also best by R_eq alone
        best_R = min(valid, key=lambda r: r['R_err'])
        print(f"  Best (R_eq only): {best_R['label']} — R_eq={best_R['R_eq']:.3f} "
              f"({best_R['R_err']:.1f}%), k={best_R['k']:.4f} ({best_R['k_err']:.1f}%)")
    else:
        print(f"  No valid minimum found!")

    return results


def save_results(
    results: List[Dict], name: str, R_values: np.ndarray,
) -> None:
    """Save results to data file."""
    filepath = f'debug/data/combined_correction_{name}.txt'
    with open(filepath, 'w') as f:
        f.write(f"# {name.upper()} Combined Correction (pade beta=2 + Fock-weighted kinetic)\n")
        f.write(f"# {'Config':>18s} {'R_eq':>8s} {'D_e':>8s} {'k':>8s} "
                f"{'R_err%':>8s} {'D_err%':>8s} {'k_err%':>8s}\n")
        for r in results:
            f.write(f"  {r['label']:>18s} {r['R_eq']:8.3f} {r['D_e']:8.4f} "
                    f"{r['k']:8.4f} {r['R_err']:8.1f} {r['D_err']:8.1f} "
                    f"{r['k_err']:8.1f}\n")
    print(f"  Saved to {filepath}")


def write_report(
    results_h2: List[Dict], results_lih: List[Dict],
) -> None:
    """Write combined analysis report."""
    filepath = 'debug/COMBINED_CORRECTION.md'

    h2_valid = {r['label']: r for r in results_h2
                if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    lih_valid = {r['label']: r for r in results_lih
                 if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}

    with open(filepath, 'w') as f:
        f.write("# Combined Correction: Attenuation + Fock-Weighted Kinetic\n\n")
        f.write("**Date:** 2026-03-12\n\n")

        f.write("## Method\n\n")
        f.write("```\n")
        f.write("V_cross(a) *= 1 / (1 + alpha_eff(a) * S_a)    [attenuation]\n")
        f.write("  alpha_eff(a) = a0 * (n/Z)^2                  [Fock-style]\n")
        f.write("h1_diag(a) += lambda * S2w(a)                  [kinetic correction]\n")
        f.write("  S2w = sum_b S(a,b)^2 * (Z_A/n_a)^2 * (Z_B/n_b)^2  [Fock-weighted]\n")
        f.write("```\n\n")

        # H2 table
        f.write("## H2 Results (nmax=3, expt: R_eq=1.401, D_e=0.1745, k=0.369)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R% | D% | k% |\n")
        f.write(f"|--------|------|-----|---|----|----|----|\n")
        for r in results_h2:
            notes = ""
            if r['D_e'] < 0:
                notes = " UNBND"
            elif r['boundary']:
                notes = " BNDRY"
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f} | {r['D_err']:.1f} | "
                    f"{r['k_err']:.1f}{notes} |\n")

        # H2 best per-molecule
        if h2_valid:
            best_h2 = min(h2_valid.values(), key=lambda r: r['R_err'] + r['k_err'])
            f.write(f"\n**Best H2:** {best_h2['label']} — "
                    f"R_eq={best_h2['R_eq']:.3f} ({best_h2['R_err']:.1f}%), "
                    f"D_e={best_h2['D_e']:.4f} ({best_h2['D_err']:.1f}%), "
                    f"k={best_h2['k']:.4f} ({best_h2['k_err']:.1f}%)\n")

        # LiH table
        f.write(f"\n## LiH Results (nmax=3, expt: R_eq=3.015, D_e=0.092, k=0.0659)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R% | D% | k% |\n")
        f.write(f"|--------|------|-----|---|----|----|----|\n")
        for r in results_lih:
            notes = ""
            if r['D_e'] < 0:
                notes = " UNBND"
            elif r['boundary']:
                notes = " BNDRY"
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f} | {r['D_err']:.1f} | "
                    f"{r['k_err']:.1f}{notes} |\n")

        # LiH best per-molecule
        if lih_valid:
            best_lih = min(lih_valid.values(), key=lambda r: r['R_err'] + r['k_err'])
            f.write(f"\n**Best LiH:** {best_lih['label']} — "
                    f"R_eq={best_lih['R_eq']:.3f} ({best_lih['R_err']:.1f}%), "
                    f"D_e={best_lih['D_e']:.4f} ({best_lih['D_err']:.1f}%), "
                    f"k={best_lih['k']:.4f} ({best_lih['k_err']:.1f}%)\n")

        # Per-molecule fit quality
        f.write(f"\n## Per-Molecule Fit Quality\n\n")
        f.write("Can each molecule be individually fit to <20% error on all metrics?\n\n")

        for mol_name, mol_valid, ref in [
            ('H2', h2_valid, H2_REF), ('LiH', lih_valid, LIH_REF)
        ]:
            good = [r for r in mol_valid.values()
                    if r['R_err'] < 20 and r['D_err'] < 30 and r['k_err'] < 30]
            if good:
                best_g = min(good, key=lambda r: r['R_err'] + r['D_err'] + r['k_err'])
                f.write(f"**{mol_name}**: YES — {best_g['label']} achieves "
                        f"R={best_g['R_err']:.1f}%, D={best_g['D_err']:.1f}%, "
                        f"k={best_g['k_err']:.1f}%\n\n")
            else:
                # Find closest
                if mol_valid:
                    closest = min(mol_valid.values(),
                                  key=lambda r: r['R_err'] + r['D_err'] + r['k_err'])
                    f.write(f"**{mol_name}**: NO — closest is {closest['label']} with "
                            f"R={closest['R_err']:.1f}%, D={closest['D_err']:.1f}%, "
                            f"k={closest['k_err']:.1f}%\n\n")
                else:
                    f.write(f"**{mol_name}**: NO valid configs found.\n\n")

        # Universality test
        f.write(f"## Universality Test\n\n")
        common = set(h2_valid.keys()) & set(lih_valid.keys())
        if common:
            # For each common config, compute combined score
            scored = []
            for label in common:
                h2 = h2_valid[label]
                lih = lih_valid[label]
                score = h2['R_err'] + h2['k_err'] + lih['R_err'] + lih['k_err']
                scored.append((label, h2, lih, score))
            scored.sort(key=lambda x: x[3])

            f.write(f"| Config | H2 R% | H2 k% | LiH R% | LiH k% | Score |\n")
            f.write(f"|--------|-------|-------|--------|--------|-------|\n")
            for label, h2, lih, score in scored[:10]:
                f.write(f"| {label} | {h2['R_err']:.1f} | {h2['k_err']:.1f} | "
                        f"{lih['R_err']:.1f} | {lih['k_err']:.1f} | {score:.1f} |\n")

            best_label, best_h2, best_lih, best_score = scored[0]
            f.write(f"\n**Best universal:** {best_label} (score={best_score:.1f})\n")

            # Check strict criteria
            h2_ok = best_h2['R_err'] < 10 and best_h2['k_err'] < 30
            lih_ok = best_lih['R_err'] < 30 and best_lih['D_err'] < 50
            if h2_ok and lih_ok:
                f.write(f"\nUNIVERSAL SUCCESS: meets criteria for both molecules.\n")
            else:
                f.write(f"\nNot universal — best config still has trade-offs.\n")
        else:
            f.write("No configs produced valid minima for both molecules.\n")

        # Key question answers
        f.write(f"\n## Key Questions\n\n")
        f.write(f"### 1. Can each molecule be individually fit?\n\n")
        for mol_name, mol_valid in [('H2', h2_valid), ('LiH', lih_valid)]:
            good = [r for r in mol_valid.values()
                    if r['R_err'] < 20 and r['D_err'] < 30 and r['k_err'] < 30]
            f.write(f"- **{mol_name}**: {'YES' if good else 'NO'} "
                    f"({len(good)} configs meet <20%R/<30%D/<30%k)\n")

        f.write(f"\n### 2. Is there universal (a0, lambda)?\n\n")
        if common and scored:
            best_label, bh2, blih, bscore = scored[0]
            f.write(f"Best universal: {best_label}\n")
            f.write(f"- H2:  R={bh2['R_err']:.1f}%, k={bh2['k_err']:.1f}%\n")
            f.write(f"- LiH: R={blih['R_err']:.1f}%, k={blih['k_err']:.1f}%\n")
        else:
            f.write("No overlap in parameter space.\n")

        f.write(f"\n### 3. What is the molecule-dependence?\n\n")
        if h2_valid and lih_valid:
            best_h2 = min(h2_valid.values(), key=lambda r: r['R_err'] + r['k_err'])
            best_lih = min(lih_valid.values(), key=lambda r: r['R_err'] + r['k_err'])
            f.write(f"- H2 prefers:  a0={best_h2['a0']}, lambda={best_h2['lam']}\n")
            f.write(f"- LiH prefers: a0={best_lih['a0']}, lambda={best_lih['lam']}\n")
            if abs(best_h2['lam'] - best_lih['lam']) > 0.01:
                f.write(f"- Lambda differs by {abs(best_h2['lam'] - best_lih['lam']):.2f}\n")
                f.write(f"- H2 needs lambda={best_h2['lam']:+.2f} (")
                f.write("soften" if best_h2['lam'] < 0 else "stiffen")
                f.write(f"), LiH needs lambda={best_lih['lam']:+.2f} (")
                f.write("soften" if best_lih['lam'] < 0 else "stiffen")
                f.write(")\n")

    print(f"  Report saved to {filepath}")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("Combined Correction Scan: Attenuation + Fock-Weighted Kinetic")
    print("=" * 80)
    sys.stdout.flush()

    # Parameter grid
    a0_values = [0.2, 0.3, 0.4]
    lam_values = [-0.10, -0.05, 0.0, +0.05, +0.10]

    # --- H2 (fast: 1,540 SDs per point) ---
    R_h2 = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])

    results_h2 = scan_molecule(
        name='H2',
        Z_A=1, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=2,
        n_electrons_A=1,
        n_electrons_B=1,
        R_values=R_h2,
        ref=H2_REF,
        a0_values=a0_values,
        lam_values=lam_values,
    )

    # --- LiH (367k SDs per point, Numba direct CI) ---
    R_lih = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0])

    results_lih = scan_molecule(
        name='LiH',
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=4,
        n_electrons_A=3,
        n_electrons_B=1,
        R_values=R_lih,
        ref=LIH_REF,
        a0_values=a0_values,
        lam_values=lam_values,
    )

    # --- Save and report ---
    save_results(results_h2, 'h2', R_h2)
    save_results(results_lih, 'lih', R_lih)
    write_report(results_h2, results_lih)

    # --- Final summary ---
    print(f"\n{'='*80}")
    print("  FINAL SUMMARY")
    print(f"{'='*80}")

    for mol_name, results, ref in [
        ('H2', results_h2, H2_REF), ('LiH', results_lih, LIH_REF)
    ]:
        valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
        if valid:
            best = min(valid, key=lambda r: r['R_err'] + r['k_err'])
            good = [r for r in valid if r['R_err'] < 20 and r['D_err'] < 30 and r['k_err'] < 30]
            print(f"\n  {mol_name}: {len(good)} configs meet <20%R/<30%D/<30%k")
            print(f"    Best (R+k): {best['label']} — "
                  f"R={best['R_eq']:.3f}({best['R_err']:.1f}%) "
                  f"D_e={best['D_e']:.4f}({best['D_err']:.1f}%) "
                  f"k={best['k']:.4f}({best['k_err']:.1f}%)")

    # Universality
    h2_valid = {r['label']: r for r in results_h2
                if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    lih_valid = {r['label']: r for r in results_lih
                 if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    common = set(h2_valid.keys()) & set(lih_valid.keys())
    if common:
        scored = [(l, h2_valid[l]['R_err'] + h2_valid[l]['k_err'] +
                      lih_valid[l]['R_err'] + lih_valid[l]['k_err'])
                  for l in common]
        scored.sort(key=lambda x: x[1])
        print(f"\n  Universal (top 3):")
        for label, score in scored[:3]:
            h2 = h2_valid[label]
            lih = lih_valid[label]
            print(f"    {label:>18s}: H2 R={h2['R_err']:.1f}%,k={h2['k_err']:.1f}% | "
                  f"LiH R={lih['R_err']:.1f}%,k={lih['k_err']:.1f}% | "
                  f"score={score:.1f}")
    sys.stdout.flush()
