"""
N-Dependent Cross-Nuclear Attenuation Scan
============================================
Core orbitals (low n/Z) should be attenuated less; valence (high n/Z) more.

alpha_eff(a) = alpha0 * (n_a / Z_a)^beta

Parameters:
    alpha0 in {0.3, 0.5, 0.7, 1.0}
    beta   in {0.5, 1.0, 2.0}
    Also test Fock-style: beta=2 (i.e. alpha_eff = alpha0 * (n/Z)^2)

f = pade (best from uniform scan)

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

# Success criteria
H2_CRITERIA = {'R_eq_pct': 15.0, 'k_pct': 20.0}
LIH_CRITERIA = {'R_eq_pct': 30.0, 'D_e_pct': 50.0}


def compute_pes(
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    R_values: np.ndarray,
    alpha0: float,
    beta: Optional[float],
    fci_method: str = 'auto',
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute PES E(R) with n-dependent pade attenuation."""
    energies = []
    for R in R_values:
        mol = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons,
            vee_method='slater_full',
            fci_method=fci_method,
            cross_nuclear_attenuation='pade',
            cross_nuclear_alpha=alpha0,
            cross_nuclear_beta=beta,
        )
        E = mol.compute_ground_state(n_states=1)[0][0]
        energies.append(E)
    return R_values, np.array(energies)


def fit_pes(R: np.ndarray, E: np.ndarray, E_sep: float) -> Dict[str, float]:
    """Extract R_eq, D_e, k from a PES via parabolic fit near minimum."""
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
    configs: List[Tuple[float, Optional[float], str]],
    fci_method: str = 'auto',
) -> List[Dict]:
    """Scan all (alpha0, beta) configs for one molecule."""

    print(f"\n{'='*80}")
    print(f"  {name} N-Dependent Attenuation Scan (pade)")
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

    # Show alpha_eff values for LiH to verify n-dependence
    if name == 'LiH':
        print(f"\n  Alpha_eff values (alpha0=0.5 examples):")
        print(f"  {'Orbital':>10s} {'n/Z':>6s}", end='')
        for beta_ex in [0.5, 1.0, 2.0]:
            print(f" {'b='+str(beta_ex):>8s}", end='')
        print()
        for orb_name, n, Z in [('Li 1s', 1, 3), ('Li 2s', 2, 3), ('Li 3s', 3, 3),
                                ('H 1s', 1, 1), ('H 2s', 2, 1), ('H 3s', 3, 1)]:
            nZ = n / Z
            print(f"  {orb_name:>10s} {nZ:6.3f}", end='')
            for beta_ex in [0.5, 1.0, 2.0]:
                a_eff = 0.5 * nZ ** beta_ex
                print(f" {a_eff:8.4f}", end='')
            print()

    print(f"\n  {'Config':>20s} {'R_eq':>8s} {'D_e':>8s} {'k':>8s} "
          f"{'R_err%':>8s} {'D_err%':>8s} {'k_err%':>8s} {'Notes':>10s}")
    print(f"  {'-'*82}")

    results = []
    t0 = time.time()

    for alpha0, beta, label in configs:
        _, E_arr = compute_pes(
            Z_A, Z_B, nmax_A, nmax_B,
            n_electrons, R_values,
            alpha0=alpha0, beta=beta,
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

        print(f"  {label:>20s} {fit['R_eq']:8.3f} {fit['D_e']:8.4f} {fit['k']:8.4f} "
              f"{R_err:8.1f} {D_err:8.1f} {k_err:8.1f} {notes:>10s}")
        sys.stdout.flush()

        results.append({
            'label': label,
            'alpha0': alpha0,
            'beta': beta,
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

    # Best by R_eq error among valid configs
    valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
    if valid:
        best = min(valid, key=lambda r: r['R_err'])
        print(f"  Best (by R_eq): {best['label']} — R_eq={best['R_eq']:.3f} "
              f"({best['R_err']:.1f}%), D_e={best['D_e']:.4f} ({best['D_err']:.1f}%), "
              f"k={best['k']:.4f} ({best['k_err']:.1f}%)")
    else:
        print(f"  No valid minimum found!")

    return results


def save_results(
    results: List[Dict], name: str, R_values: np.ndarray,
) -> None:
    """Save results to data file."""
    filepath = f'debug/data/n_dependent_attenuation_{name}.txt'
    with open(filepath, 'w') as f:
        f.write(f"# {name.upper()} N-Dependent Cross-Nuclear Attenuation (pade)\n")
        f.write(f"# {'Config':>20s} {'R_eq':>8s} {'D_e':>8s} {'k':>8s} "
                f"{'R_err%':>8s} {'D_err%':>8s} {'k_err%':>8s}\n")
        for r in results:
            f.write(f"  {r['label']:>20s} {r['R_eq']:8.3f} {r['D_e']:8.4f} "
                    f"{r['k']:8.4f} {r['R_err']:8.1f} {r['D_err']:8.1f} "
                    f"{r['k_err']:8.1f}\n")

        # PES data for bare and best
        f.write(f"\n# PES data: R  E(bare)  E(best)\n")
        bare = results[0]['E_arr']
        valid = [r for r in results if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0]
        if valid:
            best = min(valid, key=lambda r: r['R_err'])
            for i, R in enumerate(R_values):
                f.write(f"  {R:8.4f} {bare[i]:14.8f} {best['E_arr'][i]:14.8f}\n")

    print(f"  Saved to {filepath}")


def write_report(
    results_h2: List[Dict], results_lih: List[Dict],
) -> None:
    """Write analysis report."""
    filepath = 'debug/N_DEPENDENT_ATTENUATION.md'

    # Find valid configs
    h2_valid = {r['label']: r for r in results_h2
                if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    lih_valid = {r['label']: r for r in results_lih
                 if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    common = set(h2_valid.keys()) & set(lih_valid.keys())

    with open(filepath, 'w') as f:
        f.write("# N-Dependent Cross-Nuclear Attenuation\n\n")
        f.write("**Date:** 2026-03-12\n\n")

        f.write("## Motivation\n\n")
        f.write("Uniform pade(0.5) gives nearly exact H2 force constant (6.9% err) but\n")
        f.write("overcorrects R_eq. For LiH, it creates a valid minimum at R=4.68 (expt 3.02).\n")
        f.write("Core orbitals should be attenuated less than valence.\n\n")

        f.write("## Method\n\n")
        f.write("```\n")
        f.write("alpha_eff(a) = alpha0 * (n_a / Z_a)^beta\n")
        f.write("V_cross(a) *= 1 / (1 + alpha_eff(a) * S_a)\n")
        f.write("```\n\n")
        f.write("- Core (Li 1s, n/Z=0.33): alpha_eff small -> minimal attenuation\n")
        f.write("- Valence (H 1s, n/Z=1.0): alpha_eff = alpha0 -> full attenuation\n\n")

        f.write("## H2 Results (nmax=3, expt: R_eq=1.401, D_e=0.1745, k=0.369)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |\n")
        f.write(f"|--------|------|-----|---|--------|--------|--------|\n")
        for r in results_h2:
            notes = ""
            if r['D_e'] < 0:
                notes = " UNBOUND"
            elif r['boundary']:
                notes = " BOUNDARY"
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f}% | {r['D_err']:.1f}% | "
                    f"{r['k_err']:.1f}%{notes} |\n")

        f.write(f"\n## LiH Results (nmax=3, expt: R_eq=3.015, D_e=0.092, k=0.0659)\n\n")
        f.write(f"| Config | R_eq | D_e | k | R_err% | D_err% | k_err% |\n")
        f.write(f"|--------|------|-----|---|--------|--------|--------|\n")
        for r in results_lih:
            notes = ""
            if r['D_e'] < 0:
                notes = " UNBOUND"
            elif r['boundary']:
                notes = " BOUNDARY"
            f.write(f"| {r['label']} | {r['R_eq']:.3f} | {r['D_e']:.4f} | "
                    f"{r['k']:.4f} | {r['R_err']:.1f}% | {r['D_err']:.1f}% | "
                    f"{r['k_err']:.1f}%{notes} |\n")

        # Universal configs
        f.write(f"\n## Universal Configs (valid for both molecules)\n\n")
        if common:
            f.write(f"| Config | H2 R_eq (err%) | H2 k (err%) | LiH R_eq (err%) | LiH D_e (err%) |\n")
            f.write(f"|--------|---------------|------------|----------------|---------------|\n")
            sorted_common = sorted(common,
                key=lambda l: h2_valid[l]['R_err'] + lih_valid[l]['R_err'])
            for label in sorted_common:
                h2 = h2_valid[label]
                lih = lih_valid[label]
                f.write(f"| {label} | {h2['R_eq']:.3f} ({h2['R_err']:.1f}%) | "
                        f"{h2['k']:.4f} ({h2['k_err']:.1f}%) | "
                        f"{lih['R_eq']:.3f} ({lih['R_err']:.1f}%) | "
                        f"{lih['D_e']:.4f} ({lih['D_err']:.1f}%) |\n")

            best_label = sorted_common[0]
            f.write(f"\n**Best universal:** {best_label}\n")
            h2 = h2_valid[best_label]
            lih = lih_valid[best_label]
            f.write(f"- H2: R_eq={h2['R_eq']:.3f} ({h2['R_err']:.1f}%), "
                    f"D_e={h2['D_e']:.4f} ({h2['D_err']:.1f}%), "
                    f"k={h2['k']:.4f} ({h2['k_err']:.1f}%)\n")
            f.write(f"- LiH: R_eq={lih['R_eq']:.3f} ({lih['R_err']:.1f}%), "
                    f"D_e={lih['D_e']:.4f} ({lih['D_err']:.1f}%), "
                    f"k={lih['k']:.4f} ({lih['k_err']:.1f}%)\n")

            # Check success criteria
            f.write(f"\n## Success Criteria Check\n\n")
            h2_pass_R = h2['R_err'] <= H2_CRITERIA['R_eq_pct']
            h2_pass_k = h2['k_err'] <= H2_CRITERIA['k_pct']
            lih_pass_R = lih['R_err'] <= LIH_CRITERIA['R_eq_pct']
            lih_pass_D = lih['D_err'] <= LIH_CRITERIA['D_e_pct']
            f.write(f"| Criterion | Target | Achieved | Pass? |\n")
            f.write(f"|-----------|--------|----------|-------|\n")
            f.write(f"| H2 R_eq | <{H2_CRITERIA['R_eq_pct']:.0f}% | {h2['R_err']:.1f}% | {'YES' if h2_pass_R else 'NO'} |\n")
            f.write(f"| H2 k | <{H2_CRITERIA['k_pct']:.0f}% | {h2['k_err']:.1f}% | {'YES' if h2_pass_k else 'NO'} |\n")
            f.write(f"| LiH R_eq | <{LIH_CRITERIA['R_eq_pct']:.0f}% | {lih['R_err']:.1f}% | {'YES' if lih_pass_R else 'NO'} |\n")
            f.write(f"| LiH D_e | <{LIH_CRITERIA['D_e_pct']:.0f}% | {lih['D_err']:.1f}% | {'YES' if lih_pass_D else 'NO'} |\n")
        else:
            f.write("No configuration produced valid minima for both molecules.\n")

        # Comparison with uniform
        f.write(f"\n## Comparison with Uniform pade(0.5)\n\n")
        f.write("| Metric | Uniform pade(0.5) | Best n-dependent |\n")
        f.write("|--------|-------------------|------------------|\n")
        uniform_h2 = {'R_eq': 1.728, 'R_err': 23.3, 'D_e': 0.1078, 'D_err': 38.2, 'k': 0.3435, 'k_err': 6.9}
        uniform_lih = {'R_eq': 4.675, 'R_err': 55.1, 'D_e': 0.1341, 'D_err': 45.8, 'k': 0.0218, 'k_err': 66.9}
        if common:
            best_label = sorted_common[0]
            h2b = h2_valid[best_label]
            lihb = lih_valid[best_label]
            f.write(f"| H2 R_eq | {uniform_h2['R_eq']:.3f} ({uniform_h2['R_err']:.1f}%) | "
                    f"{h2b['R_eq']:.3f} ({h2b['R_err']:.1f}%) |\n")
            f.write(f"| H2 k | {uniform_h2['k']:.4f} ({uniform_h2['k_err']:.1f}%) | "
                    f"{h2b['k']:.4f} ({h2b['k_err']:.1f}%) |\n")
            f.write(f"| LiH R_eq | {uniform_lih['R_eq']:.3f} ({uniform_lih['R_err']:.1f}%) | "
                    f"{lihb['R_eq']:.3f} ({lihb['R_err']:.1f}%) |\n")
            f.write(f"| LiH D_e | {uniform_lih['D_e']:.4f} ({uniform_lih['D_err']:.1f}%) | "
                    f"{lihb['D_e']:.4f} ({lihb['D_err']:.1f}%) |\n")
        else:
            f.write("| (no valid n-dependent config) | — | — |\n")

        f.write(f"\n## Conclusion\n\n")
        if common:
            all_pass = (h2_pass_R and h2_pass_k and lih_pass_R and lih_pass_D)
            if all_pass:
                f.write(f"**SUCCESS**: {best_label} meets all criteria.\n")
            else:
                passed = sum([h2_pass_R, h2_pass_k, lih_pass_R, lih_pass_D])
                f.write(f"N-dependent attenuation passes {passed}/4 criteria with {best_label}.\n")
        else:
            f.write("N-dependent attenuation did not produce universal valid minima.\n")

    print(f"  Report saved to {filepath}")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("N-Dependent Cross-Nuclear Attenuation Scan")
    print("=" * 80)
    sys.stdout.flush()

    # Configs: (alpha0, beta, label)
    # beta=None means uniform (baseline comparison)
    configs = [
        # Uniform baseline (beta=None)
        (0.5, None, 'uniform(0.5)'),
        # alpha0=0.3
        (0.3, 0.5, 'a0=0.3,b=0.5'),
        (0.3, 1.0, 'a0=0.3,b=1.0'),
        (0.3, 2.0, 'a0=0.3,b=2.0'),
        # alpha0=0.5
        (0.5, 0.5, 'a0=0.5,b=0.5'),
        (0.5, 1.0, 'a0=0.5,b=1.0'),
        (0.5, 2.0, 'a0=0.5,b=2.0'),
        # alpha0=0.7
        (0.7, 0.5, 'a0=0.7,b=0.5'),
        (0.7, 1.0, 'a0=0.7,b=1.0'),
        (0.7, 2.0, 'a0=0.7,b=2.0'),
        # alpha0=1.0
        (1.0, 0.5, 'a0=1.0,b=0.5'),
        (1.0, 1.0, 'a0=1.0,b=1.0'),
        (1.0, 2.0, 'a0=1.0,b=2.0'),
    ]

    # --- H2 (homonuclear: n/Z always integer, beta has limited effect) ---
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
        configs=configs,
    )

    # --- LiH (heteronuclear: this is where n-dependence matters) ---
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
        configs=configs,
    )

    # --- Save and report ---
    save_results(results_h2, 'h2', R_h2)
    save_results(results_lih, 'lih', R_lih)
    write_report(results_h2, results_lih)

    # --- Final summary ---
    print(f"\n{'='*80}")
    print("  SUMMARY: N-Dependent vs Uniform")
    print(f"{'='*80}")

    h2_valid = {r['label']: r for r in results_h2
                if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    lih_valid = {r['label']: r for r in results_lih
                 if not r['boundary'] and r['k'] > 0 and r['D_e'] > 0}
    common = set(h2_valid.keys()) & set(lih_valid.keys())

    if common:
        print(f"\n  Universal configs: {len(common)}")
        sorted_common = sorted(common,
            key=lambda l: h2_valid[l]['R_err'] + lih_valid[l]['R_err'])
        for label in sorted_common[:5]:
            h2 = h2_valid[label]
            lih = lih_valid[label]
            score = h2['R_err'] + lih['R_err']
            print(f"    {label:>20s}: H2 R={h2['R_eq']:.3f}({h2['R_err']:.1f}%) "
                  f"k={h2['k']:.4f}({h2['k_err']:.1f}%) | "
                  f"LiH R={lih['R_eq']:.3f}({lih['R_err']:.1f}%) "
                  f"D_e={lih['D_e']:.4f}({lih['D_err']:.1f}%) | "
                  f"score={score:.1f}%")
    else:
        print("\n  No universal config found.")

    print(f"\n  Uniform pade(0.5) reference:")
    print(f"    H2:  R=1.728(23.3%) k=0.3435(6.9%)")
    print(f"    LiH: R=4.675(55.1%) D_e=0.1341(45.8%)")
    sys.stdout.flush()
