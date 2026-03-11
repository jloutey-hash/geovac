"""
Diagnostic script for hybrid cross-atom architecture (v0.9.18).

Hybrid = Mulliken diagonal cross-nuclear + SW/D-matrix off-diagonal coupling.
Replaces bridge mechanism for off-diagonal while preserving Fourier cross-nuclear
attraction on diagonal.

Three diagnostics:
    1. Dissociation limit: D_e_CP -> 0 at R=6.0
    2. Bound minimum: max D_e_CP > 0 in PES scan
    3. BSSE comparison: same BSSE as Mulliken (ghost atoms bypass D-matrix)
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


def run_lih(R: float, nmax: int, use_dmatrix) -> float:
    """Run LiH at given R and return E_mol."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        use_dmatrix=use_dmatrix,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return float(eigvals[0])


def run_lih_with_h1_diag(R: float, nmax: int, use_dmatrix) -> tuple:
    """Run LiH and return (E_mol, h1_diag, H1_spatial_dense)."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        use_dmatrix=use_dmatrix,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    H1_dense = mol._H1_spatial.toarray()
    return float(eigvals[0]), mol._h1_diag.copy(), H1_dense, mol._n_spatial_A


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
    print("v0.9.18 HYBRID CROSS-ATOM ARCHITECTURE DIAGNOSTIC")
    print("Mulliken diagonal + SW/D-matrix off-diagonal")
    print("=" * 70)

    # --- Atomic energies ---
    E_Li, E_H = compute_atoms(3, 1, nmax)
    E_atoms = E_Li + E_H
    print(f"\nE_Li = {E_Li:.6f} Ha")
    print(f"E_H  = {E_H:.6f} Ha")
    print(f"E_Li + E_H = {E_atoms:.6f} Ha")

    # --- BSSE (same for all methods — ghost atoms don't use D-matrix) ---
    print(f"\n--- BSSE at R={R_eq} ---")
    bsse_data = run_bsse(R_eq, nmax)
    BSSE = bsse_data['BSSE']
    print(f"BSSE = {BSSE:.6f} Ha")

    # --- Part 1: Baseline (Mulliken, standard path) ---
    print("\n" + "=" * 70)
    print("PART 1: BASELINE (Mulliken, use_dmatrix=False)")
    print("=" * 70)

    E_mol_base, h1_diag_base, H1_base, nA = run_lih_with_h1_diag(
        R_eq, nmax, use_dmatrix=False)
    D_e_raw_base = E_atoms - E_mol_base
    D_e_CP_base = D_e_raw_base + BSSE
    print(f"E_mol = {E_mol_base:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_base:.6f} Ha")
    print(f"D_e_CP = {D_e_CP_base:.6f} Ha")

    # Print key H1 matrix elements
    print(f"\n  H1 diagonal (baseline):")
    print(f"    Li 1s: {h1_diag_base[0]:.6f}")
    print(f"    H 1s:  {h1_diag_base[nA]:.6f}")
    print(f"  H1 off-diagonal (1s_A, 1s_B): {H1_base[0, nA]:.6f}")

    # --- Part 2: Hybrid ---
    print("\n" + "=" * 70)
    print("PART 2: HYBRID (Mulliken diagonal + SW off-diagonal)")
    print("=" * 70)

    E_mol_hyb, h1_diag_hyb, H1_hyb, _ = run_lih_with_h1_diag(
        R_eq, nmax, use_dmatrix='hybrid')
    D_e_raw_hyb = E_atoms - E_mol_hyb
    D_e_CP_hyb = D_e_raw_hyb + BSSE
    print(f"E_mol = {E_mol_hyb:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_hyb:.6f} Ha")
    print(f"D_e_CP = {D_e_CP_hyb:.6f} Ha")

    # Print key H1 matrix elements
    print(f"\n  H1 diagonal (hybrid):")
    print(f"    Li 1s: {h1_diag_hyb[0]:.6f}")
    print(f"    H 1s:  {h1_diag_hyb[nA]:.6f}")
    print(f"  H1 off-diagonal (1s_A, 1s_B): {H1_hyb[0, nA]:.6f}")

    # Verify diagonal matches
    diag_match = np.allclose(h1_diag_base, h1_diag_hyb, atol=1e-12)
    print(f"\n  Diagonal match with baseline: {'YES' if diag_match else 'NO'}")
    if not diag_match:
        diff = np.max(np.abs(h1_diag_base - h1_diag_hyb))
        print(f"  Max diagonal difference: {diff:.2e}")

    # Compare off-diagonals
    off_diag_base = H1_base[0, nA]  # 1s_A -> 1s_B
    off_diag_hyb = H1_hyb[0, nA]
    print(f"\n  Off-diag comparison (1s_A, 1s_B):")
    print(f"    Baseline (bridges): {off_diag_base:.6f}")
    print(f"    Hybrid (SW):        {off_diag_hyb:.6f}")
    print(f"    Ratio:              {off_diag_hyb/off_diag_base:.4f}" if abs(off_diag_base) > 1e-10 else "    Ratio: N/A")

    # --- Part 3: PES scan ---
    print("\n" + "=" * 70)
    print("PART 3: HYBRID PES SCAN")
    print("=" * 70)

    R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
    print(f"\n{'R':>8} {'E_mol':>12} {'D_e_raw':>12} {'D_e_CP':>12}")
    print("-" * 48)

    pes = []
    for R in R_values:
        E = run_lih(R, nmax, use_dmatrix='hybrid')
        De_raw = E_atoms - E
        De_cp = De_raw + BSSE
        pes.append((R, E, De_raw, De_cp))
        print(f"{R:8.3f} {E:12.6f} {De_raw:12.6f} {De_cp:12.6f}")

    # --- Part 4: Diagnostics ---
    print("\n" + "=" * 70)
    print("PART 4: DIAGNOSTICS")
    print("=" * 70)

    # Diag 1: Dissociation limit
    De_cp_6 = [p[3] for p in pes if p[0] == 6.0][0]
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

    # Diag 3: BSSE comparison
    print(f"3. BSSE = {BSSE:.6f} Ha (same as Mulliken baseline)")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY: Baseline vs Hybrid at R=3.015")
    print("=" * 70)
    print(f"{'':>20} {'Baseline':>12} {'Hybrid':>12} {'Expt':>12}")
    print(f"{'E_mol (Ha)':>20} {E_mol_base:12.6f} {E_mol_hyb:12.6f}")
    print(f"{'D_e_raw (Ha)':>20} {D_e_raw_base:12.6f} {D_e_raw_hyb:12.6f}")
    print(f"{'D_e_CP (Ha)':>20} {D_e_CP_base:12.6f} {D_e_CP_hyb:12.6f} {'0.092':>12}")
    print(f"{'BSSE (Ha)':>20} {BSSE:12.6f} {BSSE:12.6f}")

    delta = E_mol_hyb - E_mol_base
    print(f"\nHybrid - Baseline = {delta:.6f} Ha "
          f"({'lower' if delta < 0 else 'higher'} energy)")

    # Is SW coupling additive or subtractive relative to bridges?
    if delta < 0:
        print("SW off-diagonal LOWERS energy relative to bridges (additive)")
    else:
        print("SW off-diagonal RAISES energy relative to bridges (subtractive)")

    # --- Save results ---
    outfile = 'debug/data/hybrid_lih_results.txt'
    with open(outfile, 'w') as f:
        f.write("Hybrid Cross-Atom Architecture LiH Diagnostic (v0.9.18)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"nmax = {nmax}\n")
        f.write(f"Architecture: Mulliken diagonal + SW/D-matrix off-diagonal\n\n")
        f.write("Atomic energies:\n")
        f.write(f"  E_Li = {E_Li:.6f} Ha\n")
        f.write(f"  E_H  = {E_H:.6f} Ha\n")
        f.write(f"  E_atoms = {E_atoms:.6f} Ha\n\n")
        f.write(f"BSSE = {BSSE:.6f} Ha\n\n")
        f.write("Summary at R=3.015:\n")
        f.write(f"  Baseline: E={E_mol_base:.6f}, D_e_CP={D_e_CP_base:.6f} Ha\n")
        f.write(f"  Hybrid:   E={E_mol_hyb:.6f}, D_e_CP={D_e_CP_hyb:.6f} Ha\n")
        f.write(f"  Expt:     D_e = 0.092 Ha\n\n")
        f.write("H1 matrix elements at R=3.015:\n")
        f.write(f"  Diagonal Li 1s:  base={h1_diag_base[0]:.6f}, "
                f"hybrid={h1_diag_hyb[0]:.6f}\n")
        f.write(f"  Diagonal H 1s:   base={h1_diag_base[nA]:.6f}, "
                f"hybrid={h1_diag_hyb[nA]:.6f}\n")
        f.write(f"  Off-diag 1sA-1sB: base={off_diag_base:.6f}, "
                f"hybrid={off_diag_hyb:.6f}\n\n")
        f.write("PES scan (hybrid):\n")
        f.write(f"{'R':>8} {'E_mol':>12} {'D_e_raw':>12} {'D_e_CP':>12}\n")
        for r, e, dr, dc in pes:
            f.write(f"{r:8.3f} {e:12.6f} {dr:12.6f} {dc:12.6f}\n")
        f.write(f"\nDiagnostics:\n")
        f.write(f"  1. Dissociation: {'PASS' if d1_pass else 'FAIL'} "
                f"(D_e_CP(R=6) = {De_cp_6:.6f})\n")
        f.write(f"  2. Bound minimum: {'PASS' if d2_pass else 'FAIL'} "
                f"(max D_e_CP = {max_De:.6f} at R={pes[max_idx][0]:.3f})\n")
        f.write(f"  3. BSSE = {BSSE:.6f} Ha (same as baseline)\n")
    print(f"\nResults saved to {outfile}")
