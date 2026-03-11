"""
Diagnostic: Atom-dependent Sturmian p0 for LiH (v0.9.22).

Each atom uses its own self-consistent p0 from isolated-atom FCI:
  p0_A = sqrt(-2 * E_Li),  p0_B = sqrt(-2 * E_H) = 1.0

Three diagnostics at R = 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0:
  1. Dissociation limit: D_e_CP <= 0.001 Ha at R=6.0?
  2. Bound minimum: PES has minimum with D_e_CP > 0?
  3. BSSE: different from -0.115 Ha baseline?
"""

import warnings
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.lattice_index import (
    MolecularLatticeIndex, compute_bsse_correction, compute_atomic_p0
)
from geovac.wigner_so4 import bond_angle


def run_diagnostics():
    """Run full LiH atomic Sturmian diagnostics."""
    output_lines = []

    def log(msg):
        print(msg)
        output_lines.append(msg)

    log("=" * 70)
    log("LiH Atomic Sturmian Diagnostics (v0.9.22)")
    log("=" * 70)

    # --- Atomic p0 values ---
    p0_A = compute_atomic_p0(3, 3)
    p0_B = compute_atomic_p0(1, 3)
    p0_AB = np.sqrt(p0_A * p0_B)

    log(f"\nAtomic p0 values:")
    log(f"  Li (Z=3, nmax=3): p0_A = {p0_A:.6f}")
    log(f"  H  (Z=1, nmax=3): p0_B = {p0_B:.6f}")
    log(f"  Geometric mean:   p0_AB = {p0_AB:.6f}")

    # --- Numerical checks ---
    log(f"\nNumerical checks:")
    eps_Li_1s = p0_A**2 / 2.0 - 3.0 * p0_A / 1.0
    eps_H_1s = p0_B**2 / 2.0 - 1.0 * p0_B / 1.0
    log(f"  Li 1s diagonal: {eps_Li_1s:.6f} Ha (hydrogenic: -4.500 Ha)")
    log(f"  H  1s diagonal: {eps_H_1s:.6f} Ha (hydrogenic: -0.500 Ha)")

    R_eq = 3.015
    cross_Li = -(1.0 / p0_A)  # -Z_B/p0_A
    cross_H_uncapped = -(3.0 / p0_B)  # -Z_A/p0_B
    cross_H_capped = -(3.0 / R_eq)    # -Z_A/R
    log(f"  Li cross-nuc (1s, Z_B=1): {cross_Li:.6f} Ha")
    log(f"  H  cross-nuc (1s, Z_A=3): uncapped={cross_H_uncapped:.6f}, "
        f"capped={cross_H_capped:.6f} Ha")

    gamma_A = bond_angle(R_eq, p0_A)
    gamma_B = bond_angle(R_eq, p0_B)
    log(f"\nBond angles at R={R_eq}:")
    log(f"  gamma_A = {gamma_A:.4f} rad ({np.degrees(gamma_A):.1f} deg)")
    log(f"  gamma_B = {gamma_B:.4f} rad ({np.degrees(gamma_B):.1f} deg)")
    log(f"  gamma_B > gamma_A: {gamma_B > gamma_A}")

    # --- PES scan ---
    R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
    E_atoms = -7.892086  # Li + H from nmax=3 FCI

    log(f"\nPES scan (nmax=3, E_atoms={E_atoms:.6f} Ha):")
    log(f"{'R (bohr)':>10} {'E_mol (Ha)':>12} {'D_e_raw':>10} {'Status':>8}")
    log("-" * 46)

    results = []
    for R in R_values:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian='atomic',
            )
            eigvals, _ = mol.compute_ground_state(n_states=1)
            E_mol = eigvals[0]

        D_e_raw = E_atoms - E_mol
        status = "bound" if E_mol < E_atoms else "unbound"
        results.append((R, E_mol, D_e_raw, status))
        log(f"{R:10.3f} {E_mol:12.6f} {D_e_raw:10.4f} {status:>8}")

    # --- BSSE correction ---
    log(f"\nBSSE correction (R={R_eq}):")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bsse = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=R_eq, n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='auto',
        )
    log(f"  E_A_own   = {bsse['E_A_own']:.6f} Ha")
    log(f"  E_B_own   = {bsse['E_B_own']:.6f} Ha")
    log(f"  E_A_ghost = {bsse['E_A_ghost']:.6f} Ha")
    log(f"  E_B_ghost = {bsse['E_B_ghost']:.6f} Ha")
    log(f"  BSSE_A    = {bsse['BSSE_A']:.6f} Ha")
    log(f"  BSSE_B    = {bsse['BSSE_B']:.6f} Ha")
    log(f"  BSSE      = {bsse['BSSE']:.6f} Ha (baseline: -0.115 Ha)")
    BSSE = bsse['BSSE']

    # --- CP-corrected PES ---
    log(f"\nCP-corrected PES:")
    log(f"{'R (bohr)':>10} {'E_mol':>12} {'D_e_raw':>10} {'D_e_CP':>10} {'Status':>8}")
    log("-" * 56)

    for R, E_mol, D_e_raw, status in results:
        D_e_CP = D_e_raw + BSSE  # BSSE is negative, so this reduces D_e
        cp_status = "bound" if D_e_CP > 0 else "unbound"
        log(f"{R:10.3f} {E_mol:12.6f} {D_e_raw:10.4f} {D_e_CP:10.4f} {cp_status:>8}")

    # --- Summary ---
    log(f"\n{'='*70}")
    log("SUMMARY")
    log(f"{'='*70}")

    # Find minimum energy
    min_idx = np.argmin([r[1] for r in results])
    R_min, E_min, D_e_min, _ = results[min_idx]
    D_e_CP_min = D_e_min + BSSE

    log(f"  Minimum E_mol: {E_min:.6f} Ha at R={R_min:.3f} bohr")
    log(f"  D_e_raw at minimum: {D_e_min:.4f} Ha")
    log(f"  D_e_CP at minimum:  {D_e_CP_min:.4f} Ha")
    log(f"  Expt: D_e = 0.092 Ha at R_eq = 3.015 bohr")

    # Diagnostic 1: dissociation limit
    D_e_6 = results[-1][2] + BSSE
    log(f"\nDiagnostic 1 - Dissociation limit (R=6.0):")
    log(f"  D_e_CP(R=6.0) = {D_e_6:.4f} Ha")
    log(f"  |D_e_CP| <= 0.001: {'PASS' if abs(D_e_6) <= 0.001 else 'FAIL'}")

    # Diagnostic 2: bound minimum
    log(f"\nDiagnostic 2 - Bound minimum:")
    log(f"  D_e_CP at R_min = {D_e_CP_min:.4f} Ha")
    log(f"  D_e_CP > 0: {'PASS' if D_e_CP_min > 0 else 'FAIL'}")

    # Diagnostic 3: BSSE comparison
    log(f"\nDiagnostic 3 - BSSE:")
    log(f"  BSSE = {BSSE:.6f} Ha")
    log(f"  Baseline = -0.115 Ha")
    bsse_changed = abs(BSSE - (-0.115)) > 0.001
    log(f"  BSSE changed from baseline: {'YES' if bsse_changed else 'NO'}")

    # Comparison with baselines
    log(f"\nComparison with baselines:")
    log(f"  Hybrid (v0.9.18) D_e_CP: 0.143 Ha")
    log(f"  Atomic Sturmian D_e_CP:   {D_e_CP_min:.4f} Ha")
    log(f"  Experiment D_e:           0.092 Ha")

    # Success criterion
    in_range = 0.05 <= D_e_CP_min <= 0.20
    log(f"\nSuccess criterion: D_e_CP in [0.05, 0.20] Ha")
    log(f"  D_e_CP = {D_e_CP_min:.4f} Ha")
    log(f"  {'FULL SUCCESS' if in_range else 'PARTIAL SUCCESS (see diagnosis)'}")

    # Save results
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'atomic_sturmian_lih_results.txt'
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(output_lines))
    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    run_diagnostics()
