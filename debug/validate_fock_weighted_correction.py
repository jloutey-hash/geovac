"""
Fock-weighted overlap kinetic correction scan for LiH.

The uniform overlap correction failed because diffuse 2s/3s orbitals maintain
large overlaps at all R. The Fock weighting w = (Z_A/n_a)^2 * (Z_B/n_b)^2
suppresses these, concentrating the kinetic repulsion on compact core orbitals.

For LiH (Z_A=3, Z_B=1):
  Li 1s <-> H 1s: w = 9*1 = 9
  Li 1s <-> H 2s: w = 9*0.25 = 2.25
  Li 2s <-> H 1s: w = 2.25*1 = 2.25
  Li 2s <-> H 2s: w = 2.25*0.25 = 0.5625
  Li 3s <-> H 3s: w = 1*0.111 = 0.111

This gives ~80x discrimination between core-core and diffuse-diffuse overlaps.

Date: 2026-03-11
Version: v0.9.37
"""
import os
import sys
import time
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
    compute_bsse_correction,
)
from geovac.lattice_index import compute_overlap_element

NMAX = 3
LAMBDA_VALUES = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
R_SCAN = [1.5, 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0]

TARGET_R_EQ = 3.015
TARGET_D_CP = 0.092


def print_weighted_overlap_diagnostic() -> None:
    """Print weighted vs unweighted overlap sums at R=1.5 and R=5.0."""
    Z_A, Z_B = 3, 1
    print("\n" + "=" * 70)
    print("FOCK-WEIGHTED OVERLAP DIAGNOSTIC")
    print("=" * 70)

    for R in [1.5, 3.0, 5.0]:
        print(f"\n  R = {R} bohr:")
        print(f"  {'Orbital':>10s}  {'sum S^2 (unif)':>14s}  {'sum S^2*w (Fock)':>16s}  {'ratio unif':>10s}  {'ratio Fock':>10s}")
        print(f"  {'-'*70}")

        # Store values for ratio computation
        unif_vals = {}
        fock_vals = {}

        # Li orbitals
        for n_a in range(1, NMAX + 1):
            label = f"Li {n_a}s"
            sum_unif = 0.0
            sum_fock = 0.0
            for n_b in range(1, NMAX + 1):
                s = compute_overlap_element(n_a, 0, n_b, 0, Z_A, Z_B, R)
                w = (Z_A / n_a)**2 * (Z_B / n_b)**2
                sum_unif += s * s
                sum_fock += s * s * w
            unif_vals[label] = sum_unif
            fock_vals[label] = sum_fock

        # H orbitals
        for n_b in range(1, NMAX + 1):
            label = f"H  {n_b}s"
            sum_unif = 0.0
            sum_fock = 0.0
            for n_a in range(1, NMAX + 1):
                s = compute_overlap_element(n_a, 0, n_b, 0, Z_A, Z_B, R)
                w = (Z_A / n_a)**2 * (Z_B / n_b)**2
                sum_unif += s * s
                sum_fock += s * s * w
            unif_vals[label] = sum_unif
            fock_vals[label] = sum_fock

        # Print with ratios relative to R=1.5 (computed on second pass)
        for label in unif_vals:
            print(f"  {label:>10s}  {unif_vals[label]:14.6f}  {fock_vals[label]:16.6f}")

    # Now print R=1.5/R=5.0 ratios
    print(f"\n  R-discrimination (R=1.5 / R=5.0 ratio):")
    print(f"  {'Orbital':>10s}  {'uniform ratio':>14s}  {'Fock ratio':>14s}")
    print(f"  {'-'*45}")

    for n_a in range(1, NMAX + 1):
        label = f"Li {n_a}s"
        s15_u, s50_u, s15_f, s50_f = 0.0, 0.0, 0.0, 0.0
        for n_b in range(1, NMAX + 1):
            s15 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 1.5)
            s50 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 5.0)
            w = (3 / n_a)**2 * (1 / n_b)**2
            s15_u += s15**2
            s50_u += s50**2
            s15_f += s15**2 * w
            s50_f += s50**2 * w
        r_u = s15_u / max(s50_u, 1e-15)
        r_f = s15_f / max(s50_f, 1e-15)
        print(f"  {label:>10s}  {r_u:14.1f}  {r_f:14.1f}")

    for n_b in range(1, NMAX + 1):
        label = f"H  {n_b}s"
        s15_u, s50_u, s15_f, s50_f = 0.0, 0.0, 0.0, 0.0
        for n_a in range(1, NMAX + 1):
            s15 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 1.5)
            s50 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 5.0)
            w = (3 / n_a)**2 * (1 / n_b)**2
            s15_u += s15**2
            s50_u += s50**2
            s15_f += s15**2 * w
            s50_f += s50**2 * w
        r_u = s15_u / max(s50_u, 1e-15)
        r_f = s15_f / max(s50_f, 1e-15)
        print(f"  {label:>10s}  {r_u:14.1f}  {r_f:14.1f}")

    # Total
    tot15_u, tot50_u, tot15_f, tot50_f = 0.0, 0.0, 0.0, 0.0
    for n_a in range(1, NMAX + 1):
        for n_b in range(1, NMAX + 1):
            s15 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 1.5)
            s50 = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, 5.0)
            w = (3 / n_a)**2 * (1 / n_b)**2
            tot15_u += s15**2
            tot50_u += s50**2
            tot15_f += s15**2 * w
            tot50_f += s50**2 * w
    print(f"  {'TOTAL':>10s}  {tot15_u/tot50_u:14.1f}  {tot15_f/tot50_f:14.1f}")
    print()
    sys.stdout.flush()


def main() -> None:
    print("=" * 70)
    print("FOCK-WEIGHTED LiH OVERLAP KINETIC CORRECTION SCAN")
    print(f"Lambda values: {LAMBDA_VALUES}")
    print(f"R scan: {R_SCAN}")
    print("=" * 70)
    sys.stdout.flush()

    t_start = time.time()

    # Diagnostic: overlap R-discrimination
    print_weighted_overlap_diagnostic()

    # Atomic references
    print("Atomic references...")
    sys.stdout.flush()
    li = LatticeIndex(n_electrons=3, max_n=NMAX, nuclear_charge=3,
                      vee_method='slater_full', h1_method='exact', fci_method='auto')
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(n_electrons=1, max_n=NMAX, nuclear_charge=1,
                     vee_method='slater_full', h1_method='exact', fci_method='auto')
    E_h = h.compute_ground_state(n_states=1)[0][0]
    E_sep = E_li + E_h
    print(f"  E(Li)={E_li:.6f}, E(H)={E_h:.6f}, E_sep={E_sep:.6f}")
    sys.stdout.flush()

    # BSSE
    print("\nBSSE...")
    sys.stdout.flush()
    bsse_data = compute_bsse_correction(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=3.015, n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
    )
    bsse = bsse_data['BSSE']
    print(f"  BSSE = {bsse:.6f}")
    sys.stdout.flush()

    # Scan
    all_results: Dict[float, List[Dict[str, float]]] = {}

    for lam in LAMBDA_VALUES:
        print(f"\n{'='*50}")
        print(f"  lambda = {lam}  (Fock-weighted)")
        print(f"{'='*50}")
        sys.stdout.flush()
        results = []
        for R in R_SCAN:
            t0 = time.time()
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
                R=R, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                cross_nuclear_method='exact',
                cross_atom_vee=True,
                t_corr_lambda=lam,
                t_corr_fock_weighted=True,
            )
            eigvals, _ = mol.compute_ground_state(n_states=1)
            E_mol = eigvals[0]
            D_raw = E_sep - E_mol
            D_cp = D_raw + bsse
            dt = time.time() - t0
            print(f"  R={R:.3f}: E={E_mol:.5f}, D_cp={D_cp:.4f} ({dt:.0f}s)")
            sys.stdout.flush()
            results.append({'R': R, 'E_mol': E_mol, 'D_raw': D_raw, 'D_cp': D_cp})
        all_results[lam] = results

        # Check for equilibrium
        D_arr = np.array([r['D_cp'] for r in results])
        idx_max = np.argmax(D_arr)
        if 0 < idx_max < len(D_arr) - 1:
            # Parabolic interpolation
            R_arr = np.array([r['R'] for r in results])
            R0, R1, R2 = R_arr[idx_max-1], R_arr[idx_max], R_arr[idx_max+1]
            D0, D1, D2 = D_arr[idx_max-1], D_arr[idx_max], D_arr[idx_max+1]
            denom = 2.0 * ((R1-R0)*(D1-D2) - (R1-R2)*(D1-D0))
            if abs(denom) > 1e-15:
                R_eq = R1 - ((R1-R0)**2*(D1-D2) - (R1-R2)**2*(D1-D0)) / denom
            else:
                R_eq = R1
            print(f"  => EQUILIBRIUM near R={R_eq:.3f}, D_cp_max={D_arr[idx_max]:.4f}")
        else:
            if D_arr[0] >= D_arr[-1]:
                print(f"  => NO EQUILIBRIUM (D_cp monotonically decreasing)")
            else:
                print(f"  => NO EQUILIBRIUM (D_cp increasing toward large R)")
        sys.stdout.flush()

    # Generate outputs
    print("\n\nGenerating output files...")
    sys.stdout.flush()

    # Data file
    datapath = os.path.join(PROJECT_ROOT, "debug", "data",
                            "lih_fock_weighted_correction.txt")
    os.makedirs(os.path.dirname(datapath), exist_ok=True)
    with open(datapath, 'w', encoding='utf-8') as f:
        f.write("LiH Fock-Weighted Kinetic Correction Scan -- v0.9.37, 2026-03-11\n")
        f.write(f"nmax={NMAX}, cross_nuclear=exact, cross_atom_vee=True\n")
        f.write(f"E(Li)={E_li:.6f}, E(H)={E_h:.6f}, E_sep={E_sep:.6f}\n")
        f.write(f"BSSE={bsse:.6f}\n")
        f.write(f"Weighting: Fock w=(Z_A/n_a)^2*(Z_B/n_b)^2\n\n")

        for lam, results in all_results.items():
            D_arr = np.array([r['D_cp'] for r in results])
            idx_max = np.argmax(D_arr)
            has_eq = 0 < idx_max < len(D_arr) - 1
            status = f"EQ near R={results[idx_max]['R']:.1f}" if has_eq else "NO EQ"
            f.write(f"lambda = {lam}  [{status}]\n")
            f.write(f"  {'R':>7s}  {'E_mol':>10s}  {'D_raw':>8s}  {'D_cp':>8s}\n")
            f.write(f"  {'-'*42}\n")
            for r in results:
                f.write(f"  {r['R']:7.3f}  {r['E_mol']:10.5f}  {r['D_raw']:8.4f}  {r['D_cp']:8.4f}\n")
            f.write("\n")

        f.write("\nSUMMARY\n")
        f.write(f"{'lambda':>8s}  {'R_eq':>8s}  {'D_cp':>8s}  {'status':>20s}\n")
        f.write("-" * 50 + "\n")
        for lam, results in all_results.items():
            D_arr = np.array([r['D_cp'] for r in results])
            R_arr = np.array([r['R'] for r in results])
            idx_max = np.argmax(D_arr)
            has_eq = 0 < idx_max < len(D_arr) - 1
            if has_eq:
                f.write(f"{lam:8.3f}  {R_arr[idx_max]:8.3f}  {D_arr[idx_max]:8.4f}  {'EQUILIBRIUM':>20s}\n")
            else:
                f.write(f"{lam:8.3f}  {'---':>8s}  {'---':>8s}  {'NO EQUILIBRIUM':>20s}\n")

    print(f"Data saved to {datapath}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(all_results)))

        for (lam, results), color in zip(all_results.items(), colors):
            R_arr = [r['R'] for r in results]
            D_cp = [r['D_cp'] for r in results]
            D_arr_np = np.array(D_cp)
            idx_max = np.argmax(D_arr_np)
            has_eq = 0 < idx_max < len(D_arr_np) - 1
            label = f'lam={lam}'
            if has_eq:
                label += f' (eq~{R_arr[idx_max]:.1f})'
            ax1.plot(R_arr, D_cp, 'o-', color=color, label=label, markersize=4)
            if has_eq:
                ax1.plot(R_arr[idx_max], D_cp[idx_max], 's', color=color, markersize=8)

        ax1.axhline(TARGET_D_CP, color='red', linestyle='--', alpha=0.5, label=f'expt D_cp={TARGET_D_CP}')
        ax1.axvline(TARGET_R_EQ, color='red', linestyle=':', alpha=0.5, label=f'expt R_eq={TARGET_R_EQ}')
        ax1.axhline(0.0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
        ax1.set_xlabel('R (bohr)')
        ax1.set_ylabel('D_cp (Ha)')
        ax1.set_title('Fock-Weighted: D_cp(R)')
        ax1.legend(fontsize=7, loc='upper right')
        ax1.grid(True, alpha=0.3)

        # Compare Fock vs uniform for best lambda
        # Find best Fock lambda (one with equilibrium closest to experiment)
        best_lam = None
        for lam, results in all_results.items():
            D_arr_np = np.array([r['D_cp'] for r in results])
            idx = np.argmax(D_arr_np)
            if 0 < idx < len(D_arr_np) - 1:
                if best_lam is None:
                    best_lam = lam
        if best_lam is not None and best_lam in all_results:
            R_fock = [r['R'] for r in all_results[best_lam]]
            D_fock = [r['D_cp'] for r in all_results[best_lam]]
            ax2.plot(R_fock, D_fock, 'o-', color='blue', label=f'Fock lam={best_lam}', linewidth=2)

        # baseline
        R_base = [r['R'] for r in all_results[0.0]]
        D_base = [r['D_cp'] for r in all_results[0.0]]
        ax2.plot(R_base, D_base, 'o-', color='gray', label='baseline (lam=0)', linewidth=1)

        ax2.axhline(TARGET_D_CP, color='red', linestyle='--', alpha=0.5, label=f'expt D_cp')
        ax2.axvline(TARGET_R_EQ, color='red', linestyle=':', alpha=0.5, label=f'expt R_eq')
        ax2.axhline(0.0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
        ax2.set_xlabel('R (bohr)')
        ax2.set_ylabel('D_cp (Ha)')
        ax2.set_title('Best Fock vs Baseline')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        fig.suptitle('LiH Fock-Weighted Kinetic Correction (nmax=3)', fontsize=13)
        fig.tight_layout()
        plotpath = os.path.join(PROJECT_ROOT, "debug", "plots",
                                "lih_fock_weighted_correction.png")
        fig.savefig(plotpath, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved to {plotpath}")
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot.")

    dt = time.time() - t_start
    print(f"\nTotal time: {dt:.0f}s ({dt/60:.1f} min)")
    print("DONE")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
