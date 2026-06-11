"""
Fast overlap kinetic correction scan — reduced grid.

Tests lambda = [0.0, 0.5, 2.0, 5.0, 10.0] at R = [1.5, 2.5, 3.0, 4.0, 5.0]
to quickly determine if ANY lambda creates an equilibrium.

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

NMAX = 3
LAMBDA_VALUES = [0.0, 0.5, 2.0, 5.0, 10.0, 20.0]
R_SCAN = [1.5, 2.5, 3.0, 4.0, 5.0]

TARGET_R_EQ = 3.015
TARGET_D_CP = 0.092


def main() -> None:
    print("=" * 70)
    print("FAST LiH OVERLAP KINETIC CORRECTION SCAN")
    print(f"Lambda values: {LAMBDA_VALUES}")
    print(f"R scan: {R_SCAN}")
    print("=" * 70)
    sys.stdout.flush()

    t_start = time.time()

    # Atomic references
    print("\nAtomic references...")
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
        print(f"\n--- lambda = {lam} ---")
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
            )
            eigvals, _ = mol.compute_ground_state(n_states=1)
            E_mol = eigvals[0]
            D_raw = E_sep - E_mol
            D_cp = D_raw + bsse
            dt = time.time() - t0
            print(f"  R={R:.1f}: E={E_mol:.5f}, D_cp={D_cp:.4f} ({dt:.0f}s)")
            sys.stdout.flush()
            results.append({'R': R, 'E_mol': E_mol, 'D_raw': D_raw, 'D_cp': D_cp})
        all_results[lam] = results

        # Check for equilibrium
        D_arr = np.array([r['D_cp'] for r in results])
        idx_max = np.argmax(D_arr)
        if 0 < idx_max < len(D_arr) - 1:
            print(f"  => EQUILIBRIUM near R={results[idx_max]['R']}, D_cp={D_arr[idx_max]:.4f}")
        else:
            if D_arr[0] >= D_arr[-1]:
                print(f"  => NO EQUILIBRIUM (D_cp monotonically decreasing)")
            else:
                print(f"  => NO EQUILIBRIUM (D_cp increasing toward large R)")
        sys.stdout.flush()

    # Also print overlap magnitudes for reference
    print("\n\n--- OVERLAP MAGNITUDES at each R ---")
    from geovac.hamiltonian import compute_overlap_element
    for R in R_SCAN:
        sum_s2 = 0.0
        pairs = []
        for n_a in range(1, NMAX + 1):
            for n_b in range(1, NMAX + 1):
                s = compute_overlap_element(n_a, 0, n_b, 0, 3, 1, R)
                sum_s2 += s * s
                if abs(s) > 0.01:
                    pairs.append(f"({n_a}s,{n_b}s)={s:.3f}")
        print(f"  R={R:.1f}: sum S^2 = {sum_s2:.4f}  [{', '.join(pairs)}]")
    sys.stdout.flush()

    # Summary
    print("\n\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for lam, results in all_results.items():
        D_arr = [r['D_cp'] for r in results]
        R_arr = [r['R'] for r in results]
        idx_max = np.argmax(D_arr)
        has_eq = 0 < idx_max < len(D_arr) - 1
        if has_eq:
            print(f"  lam={lam:6.1f}: EQUILIBRIUM at R~{R_arr[idx_max]:.1f}, D_cp={D_arr[idx_max]:.4f}")
        else:
            print(f"  lam={lam:6.1f}: NO EQ, D_cp range [{min(D_arr):.4f}, {max(D_arr):.4f}]")

    dt = time.time() - t_start
    print(f"\nTotal time: {dt:.0f}s ({dt/60:.1f} min)")
    print("DONE")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
