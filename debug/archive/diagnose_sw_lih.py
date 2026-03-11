"""
Diagnostic script for Shibuya-Wulfman integration (v0.9.17).

Compares baseline (Mulliken, use_dmatrix=False) with SW coupling
(use_dmatrix=True, now using SW integrals instead of geometric-mean).

Three diagnostics:
    1. Dissociation limit: D_e_CP -> 0 at R=6.0
    2. Bound minimum: max D_e_CP > 0 in PES scan
    3. BSSE comparison: same BSSE (ghost atoms bypass D-matrix)
"""

import sys
import warnings
import numpy as np

sys.path.insert(0, '.')
warnings.filterwarnings('ignore')

from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
)


def compute_atoms(Z_A: int, Z_B: int, nmax: int) -> tuple:
    """Compute isolated atomic energies."""
    li_A = LatticeIndex(
        n_electrons=Z_A, max_n=nmax, nuclear_charge=Z_A,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_A = li_A.compute_ground_state(n_states=1)[0][0]

    li_B = LatticeIndex(
        n_electrons=Z_B, max_n=nmax, nuclear_charge=Z_B,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_B = li_B.compute_ground_state(n_states=1)[0][0]
    return E_A, E_B


def run_lih(R: float, nmax: int, use_dmatrix: bool) -> float:
    """Run LiH at given R and return E_mol."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        use_dmatrix=use_dmatrix,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return float(eigvals[0])


def run_bsse(R: float, nmax: int) -> dict:
    """Compute BSSE at given R (ghost atoms bypass D-matrix)."""
    return compute_bsse_correction(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
    )


if __name__ == '__main__':
    nmax = 3
    R_eq = 3.015

    print("=" * 70)
    print("PART 1: BASELINE (Mulliken, use_dmatrix=False)")
    print("=" * 70)

    E_Li, E_H = compute_atoms(3, 1, nmax)
    E_atoms = E_Li + E_H
    print(f"\nE_Li = {E_Li:.6f} Ha")
    print(f"E_H  = {E_H:.6f} Ha")
    print(f"E_Li + E_H = {E_atoms:.6f} Ha")

    print(f"\n--- LiH at R={R_eq}, Mulliken ---")
    E_mol_base = run_lih(R_eq, nmax, use_dmatrix=False)
    D_e_raw_base = E_atoms - E_mol_base
    print(f"E_mol = {E_mol_base:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_base:.6f} Ha")

    print(f"\n--- BSSE at R={R_eq} ---")
    bsse_data = run_bsse(R_eq, nmax)
    BSSE = bsse_data['BSSE']
    D_e_CP_base = D_e_raw_base + BSSE
    print(f"BSSE = {BSSE:.6f} Ha")
    print(f"D_e_CP = {D_e_CP_base:.6f} Ha")

    print("\n" + "=" * 70)
    print("PART 2: SHIBUYA-WULFMAN (use_dmatrix=True)")
    print("=" * 70)

    print(f"\n--- LiH at R={R_eq}, SW ---")
    E_mol_sw = run_lih(R_eq, nmax, use_dmatrix=True)
    D_e_raw_sw = E_atoms - E_mol_sw
    print(f"E_mol = {E_mol_sw:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_sw:.6f} Ha")
    D_e_CP_sw = D_e_raw_sw + BSSE
    print(f"D_e_CP = {D_e_CP_sw:.6f} Ha")

    print("\n" + "=" * 70)
    print("PART 3: SW PES SCAN")
    print("=" * 70)

    R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 6.0]
    print(f"\n{'R':>8} {'E_mol':>12} {'D_e_raw':>12} {'D_e_CP':>12}")
    print("-" * 48)

    pes = []
    for R in R_values:
        E = run_lih(R, nmax, use_dmatrix=True)
        De_raw = E_atoms - E
        De_cp = De_raw + BSSE
        pes.append((R, E, De_raw, De_cp))
        print(f"{R:8.3f} {E:12.6f} {De_raw:12.6f} {De_cp:12.6f}")

    print("\n" + "=" * 70)
    print("PART 4: DIAGNOSTICS")
    print("=" * 70)

    # Diag 1: Dissociation limit
    De_cp_6 = pes[-1][3]
    d1_pass = abs(De_cp_6) < 0.01
    print(f"\n1. Dissociation (R=6.0): D_e_CP = {De_cp_6:.6f} Ha "
          f"{'PASS' if d1_pass else 'FAIL'} (|D_e_CP| < 0.01)")

    # Diag 2: Bound minimum
    De_cps = [p[3] for p in pes]
    max_De = max(De_cps)
    max_idx = np.argmax(De_cps)
    d2_pass = max_De > 0
    print(f"2. Bound minimum: max D_e_CP = {max_De:.6f} Ha "
          f"at R={pes[max_idx][0]:.3f} "
          f"{'PASS' if d2_pass else 'FAIL'}")

    # Diag 3: Coupling magnitude
    from geovac.wigner_so4 import bond_angle
    from geovac.shibuya_wulfman import sw_form_factor
    p0 = np.sqrt(3**2 + 1**2)
    gamma = bond_angle(R_eq, p0)
    Z_eff = (3 + 1) / 2.0
    kappa_sw = Z_eff / p0 * sw_form_factor(1, gamma)
    print(f"3. SW coupling: kappa_SW(1s) = {kappa_sw:.6f} Ha "
          f"(Z_eff/p0={Z_eff/p0:.4f}, sin(gamma)={np.sin(gamma):.4f})")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Mulliken vs SW at R=3.015")
    print("=" * 70)
    print(f"{'':>20} {'Mulliken':>12} {'SW':>12} {'Expt':>12}")
    print(f"{'E_mol (Ha)':>20} {E_mol_base:12.6f} {E_mol_sw:12.6f}")
    print(f"{'D_e_raw (Ha)':>20} {D_e_raw_base:12.6f} {D_e_raw_sw:12.6f}")
    print(f"{'D_e_CP (Ha)':>20} {D_e_CP_base:12.6f} {D_e_CP_sw:12.6f} {'0.092':>12}")
    print(f"{'BSSE (Ha)':>20} {BSSE:12.6f} {BSSE:12.6f}")

    # Save results
    outfile = 'debug/data/sw_lih_results.txt'
    with open(outfile, 'w') as f:
        f.write("Shibuya-Wulfman LiH Diagnostic Results (v0.9.17)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"nmax = {nmax}\n")
        f.write(f"p0 = {p0:.6f}\n")
        f.write(f"gamma(R_eq) = {gamma:.6f} rad ({np.degrees(gamma):.2f} deg)\n")
        f.write(f"Z_eff = {Z_eff:.1f}\n")
        f.write(f"f(1,gamma) = sin(gamma) = {np.sin(gamma):.6f}\n")
        f.write(f"kappa_SW(1s) = {kappa_sw:.6f} Ha\n\n")
        f.write("Atomic energies:\n")
        f.write(f"  E_Li = {E_Li:.6f} Ha\n")
        f.write(f"  E_H  = {E_H:.6f} Ha\n")
        f.write(f"  E_atoms = {E_atoms:.6f} Ha\n\n")
        f.write(f"BSSE = {BSSE:.6f} Ha\n\n")
        f.write("Summary at R=3.015:\n")
        f.write(f"  Mulliken: E={E_mol_base:.6f}, D_e_CP={D_e_CP_base:.6f} Ha\n")
        f.write(f"  SW:       E={E_mol_sw:.6f}, D_e_CP={D_e_CP_sw:.6f} Ha\n")
        f.write(f"  Expt:     D_e = 0.092 Ha\n\n")
        f.write("PES scan:\n")
        f.write(f"{'R':>8} {'E_mol':>12} {'D_e_raw':>12} {'D_e_CP':>12}\n")
        for r, e, dr, dc in pes:
            f.write(f"{r:8.3f} {e:12.6f} {dr:12.6f} {dc:12.6f}\n")
        f.write(f"\nDiagnostics:\n")
        f.write(f"  1. Dissociation: {'PASS' if d1_pass else 'FAIL'} "
                f"(D_e_CP(R=6) = {De_cp_6:.6f})\n")
        f.write(f"  2. Bound minimum: {'PASS' if d2_pass else 'FAIL'} "
                f"(max D_e_CP = {max_De:.6f} at R={pes[max_idx][0]:.3f})\n")
    print(f"\nResults saved to {outfile}")
