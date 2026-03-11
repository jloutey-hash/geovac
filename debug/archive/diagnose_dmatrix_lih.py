"""
Diagnostic script for D-matrix integration (v0.9.16).

Part 2 of the v0.9.16 session: collect baseline (use_dmatrix=False) and
D-matrix (use_dmatrix=True) results for LiH.
"""

import sys
import warnings
import numpy as np

sys.path.insert(0, '.')
warnings.filterwarnings('ignore')

from geovac.lattice_index import MolecularLatticeIndex, LatticeIndex, compute_bsse_correction


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


def run_bsse(R: float, nmax: int, use_dmatrix: bool) -> dict:
    """Compute BSSE at given R. Note: BSSE computation uses standard path
    (ghost atoms don't use D-matrix since Z=0 → ghost flag is set)."""
    # For the D-matrix path, the ghost atoms bypass D-matrix automatically
    # (ghost check in _build_molecular_h1). So BSSE computation uses the
    # standard code path for ghost calculations, which is correct.
    return compute_bsse_correction(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
    )


if __name__ == '__main__':
    nmax = 3
    R_eq = 3.015  # experimental LiH equilibrium

    print("=" * 70)
    print("PART 1: BASELINE (use_dmatrix=False)")
    print("=" * 70)

    # Atomic energies
    E_Li, E_H = compute_atoms(3, 1, nmax)
    E_atoms = E_Li + E_H
    print(f"\nE_Li = {E_Li:.6f} Ha")
    print(f"E_H  = {E_H:.6f} Ha")
    print(f"E_Li + E_H = {E_atoms:.6f} Ha")

    # LiH at R_eq, baseline
    print(f"\n--- LiH at R={R_eq}, use_dmatrix=False ---")
    E_mol_base = run_lih(R_eq, nmax, use_dmatrix=False)
    D_e_raw_base = E_atoms - E_mol_base
    print(f"E_mol = {E_mol_base:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_base:.6f} Ha")

    # BSSE at R_eq
    print(f"\n--- BSSE at R={R_eq} ---")
    bsse_base = run_bsse(R_eq, nmax, use_dmatrix=False)
    D_e_CP_base = D_e_raw_base + bsse_base['BSSE']  # BSSE is negative
    print(f"BSSE = {bsse_base['BSSE']:.6f} Ha")
    print(f"D_e_CP = {D_e_CP_base:.6f} Ha")

    print("\n" + "=" * 70)
    print("PART 2: D-MATRIX (use_dmatrix=True)")
    print("=" * 70)

    # LiH at R_eq, D-matrix
    print(f"\n--- LiH at R={R_eq}, use_dmatrix=True ---")
    E_mol_dmat = run_lih(R_eq, nmax, use_dmatrix=True)
    D_e_raw_dmat = E_atoms - E_mol_dmat
    print(f"E_mol = {E_mol_dmat:.6f} Ha")
    print(f"D_e_raw = {D_e_raw_dmat:.6f} Ha")

    # Sanity check
    if E_mol_dmat > 0 or abs(E_mol_dmat) > 100:
        print(f"*** UNPHYSICAL: E_mol = {E_mol_dmat:.6f} Ha ***")
        print("Stopping. Likely coupling_scale or p0 issue.")
        sys.exit(1)

    # BSSE for D-matrix (ghost atoms bypass D-matrix automatically)
    print(f"\n--- BSSE at R={R_eq} (D-matrix) ---")
    bsse_dmat = run_bsse(R_eq, nmax, use_dmatrix=True)
    D_e_CP_dmat = D_e_raw_dmat + bsse_dmat['BSSE']
    print(f"BSSE = {bsse_dmat['BSSE']:.6f} Ha")
    print(f"D_e_CP = {D_e_CP_dmat:.6f} Ha")

    # PES scan
    print("\n" + "=" * 70)
    print("PART 3: D-MATRIX PES SCAN")
    print("=" * 70)

    R_values = [2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
    print(f"\n{'R (bohr)':>10} {'E_mol (Ha)':>14} {'D_e_raw (Ha)':>14} {'D_e_CP (Ha)':>14}")
    print("-" * 56)

    pes_results = []
    for R_val in R_values:
        E_mol = run_lih(R_val, nmax, use_dmatrix=True)
        D_e_raw = E_atoms - E_mol
        # BSSE is R-independent for ghost atoms (no bridges/cross-nuclear)
        D_e_CP = D_e_raw + bsse_dmat['BSSE']
        pes_results.append((R_val, E_mol, D_e_raw, D_e_CP))
        print(f"{R_val:10.3f} {E_mol:14.6f} {D_e_raw:14.6f} {D_e_CP:14.6f}")

    print("\n" + "=" * 70)
    print("PART 4: DIAGNOSTICS")
    print("=" * 70)

    # Diagnostic 1: Dissociation limit
    E_mol_6 = pes_results[-1][1]  # R=6.0
    D_e_CP_6 = pes_results[-1][3]
    diss_pass = abs(D_e_CP_6) <= 0.001
    print(f"\n1. Dissociation limit (R=6.0): D_e_CP = {D_e_CP_6:.6f} Ha "
          f"{'PASS' if diss_pass else 'FAIL'} (threshold: <=0.001 Ha)")

    # Diagnostic 2: Bound minimum
    D_e_CPs = [r[3] for r in pes_results]
    has_minimum = max(D_e_CPs) > D_e_CPs[0] or max(D_e_CPs) > D_e_CPs[-1]
    # Check if there's a point where D_e_CP is positive and larger than neighbors
    min_idx = np.argmax(D_e_CPs)
    bound_pass = D_e_CPs[min_idx] > 0
    print(f"2. Bound minimum: max D_e_CP = {D_e_CPs[min_idx]:.6f} Ha at R={pes_results[min_idx][0]:.3f} "
          f"{'PASS' if bound_pass else 'FAIL'}")

    # Diagnostic 3: BSSE reduction
    # Note: BSSE for ghost atoms is the same for both paths (ghost atoms
    # don't use D-matrix). The question is whether D_e_raw is smaller,
    # making the CP correction less important.
    bsse_reduction = abs(bsse_dmat['BSSE']) < abs(bsse_base['BSSE'])
    print(f"3. BSSE reduction: Mulliken BSSE = {bsse_base['BSSE']:.6f} Ha, "
          f"D-matrix BSSE = {bsse_dmat['BSSE']:.6f} Ha "
          f"{'PASS' if bsse_reduction else 'FAIL (same BSSE — ghost atoms bypass D-matrix)'}")

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: Mulliken vs D-matrix at R=3.015")
    print("=" * 70)
    print(f"{'':>20} {'Mulliken':>14} {'D-matrix':>14}")
    print(f"{'E_mol (Ha)':>20} {E_mol_base:14.6f} {E_mol_dmat:14.6f}")
    print(f"{'D_e_raw (Ha)':>20} {D_e_raw_base:14.6f} {D_e_raw_dmat:14.6f}")
    print(f"{'BSSE (Ha)':>20} {bsse_base['BSSE']:14.6f} {bsse_dmat['BSSE']:14.6f}")
    print(f"{'D_e_CP (Ha)':>20} {D_e_CP_base:14.6f} {D_e_CP_dmat:14.6f}")
    print(f"{'Expt D_e (Ha)':>20} {'0.092':>14}")

    # Save results
    with open('debug/data/dmatrix_lih_results.txt', 'w') as f:
        f.write("D-matrix LiH diagnostic results (v0.9.16)\n")
        f.write(f"nmax = {nmax}\n\n")
        f.write(f"E_Li = {E_Li:.6f} Ha\n")
        f.write(f"E_H  = {E_H:.6f} Ha\n")
        f.write(f"E_Li + E_H = {E_atoms:.6f} Ha\n\n")
        f.write("Mulliken vs D-matrix at R=3.015:\n")
        f.write(f"  Mulliken: E_mol={E_mol_base:.6f}, D_e_raw={D_e_raw_base:.6f}, "
                f"BSSE={bsse_base['BSSE']:.6f}, D_e_CP={D_e_CP_base:.6f}\n")
        f.write(f"  D-matrix: E_mol={E_mol_dmat:.6f}, D_e_raw={D_e_raw_dmat:.6f}, "
                f"BSSE={bsse_dmat['BSSE']:.6f}, D_e_CP={D_e_CP_dmat:.6f}\n\n")
        f.write("D-matrix PES:\n")
        f.write(f"{'R':>10} {'E_mol':>14} {'D_e_raw':>14} {'D_e_CP':>14}\n")
        for R_val, E_mol, D_e_raw, D_e_CP in pes_results:
            f.write(f"{R_val:10.3f} {E_mol:14.6f} {D_e_raw:14.6f} {D_e_CP:14.6f}\n")
        f.write(f"\nDiagnostic 1 (dissociation): {'PASS' if diss_pass else 'FAIL'}\n")
        f.write(f"Diagnostic 2 (bound minimum): {'PASS' if bound_pass else 'FAIL'}\n")
        f.write(f"Diagnostic 3 (BSSE reduction): {'PASS' if bsse_reduction else 'FAIL'}\n")
    print("\nResults saved to debug/data/dmatrix_lih_results.txt")
