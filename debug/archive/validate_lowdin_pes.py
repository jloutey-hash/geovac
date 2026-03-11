"""
LiH PES scan: orthogonalized vs non-orthogonalized molecular basis.

Compares MolecularLatticeIndex with orthogonalize=True/False at nmax=3
across the same geometry points as Table III in paper_geovac_fci.tex.

Reports: E_mol, D_e_raw, D_e_CP, and condition number of overlap matrix S.

Date: 2026-03-08
Version: v0.9.12 (Lowdin orthogonalization)
"""

import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac.lattice_index import (
    MolecularLatticeIndex,
    LatticeIndex,
    compute_bsse_correction,
)

# --- Constants ---
R_POINTS = [2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0, 8.0]
Z_A, Z_B = 3, 1
NMAX = 3
N_ELECTRONS = 4
D_E_EXPT = 0.0924  # Ha

OUTPUT_FILE = os.path.join(os.path.dirname(__file__),
                           "data", "lih_lowdin_pes.txt")


def get_atomic_energies(nmax: int) -> tuple:
    """Compute separated atom energies (own basis)."""
    li = LatticeIndex(
        n_electrons=3, max_n=nmax, nuclear_charge=Z_A,
        vee_method='slater_full', h1_method='exact',
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=nmax, nuclear_charge=Z_B,
        vee_method='slater_full', h1_method='exact',
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    return E_li, E_h


def run_single_R(
    R: float, nmax: int, orthogonalize: bool, E_sep: float
) -> dict:
    """Run LiH at one R value, return energies and timing."""
    t0 = time.perf_counter()

    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=N_ELECTRONS,
        n_bridges=10, vee_method='slater_full',
        fci_method='auto', orthogonalize=orthogonalize,
    )

    # Get overlap matrix condition number (before orthogonalization)
    S = mol._compute_overlap_matrix()
    eigvals_S = np.linalg.eigvalsh(S)
    cond_S = eigvals_S[-1] / max(eigvals_S[0], 1e-16)
    min_eigval_S = eigvals_S[0]

    eigvals, _ = mol.compute_ground_state(n_states=1)
    E_mol = eigvals[0]
    t_mol = time.perf_counter() - t0

    D_e_raw = E_sep - E_mol

    # BSSE / counterpoise
    t_bsse_start = time.perf_counter()
    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax, R=R,
        n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
        orthogonalize=orthogonalize,
    )
    t_bsse = time.perf_counter() - t_bsse_start

    E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
    D_e_cp = E_ghost_sum - E_mol

    return {
        'R': R,
        'E_mol': E_mol,
        'D_e_raw': D_e_raw,
        'D_e_cp': D_e_cp,
        'BSSE': bsse['BSSE'],
        'cond_S': cond_S,
        'min_eigval_S': min_eigval_S,
        't_mol': t_mol,
        't_bsse': t_bsse,
    }


def main() -> None:
    log: list = []

    header = "LiH PES: Lowdin Orthogonalization Comparison (v0.9.12)"
    log.append(header)
    log.append(f"Date: 2026-03-08")
    log.append(f"nmax={NMAX}, Z_A={Z_A} (Li), Z_B={Z_B} (H), Ne={N_ELECTRONS}")
    log.append(f"D_e(expt) = {D_E_EXPT} Ha")
    log.append("")
    print(header)

    # --- Atomic reference ---
    print("\n--- Atomic reference energies ---")
    E_li, E_h = get_atomic_energies(NMAX)
    E_sep = E_li + E_h
    log.append(f"E(Li, nmax={NMAX}) = {E_li:.6f} Ha")
    log.append(f"E(H,  nmax={NMAX}) = {E_h:.6f} Ha")
    log.append(f"E_sep = {E_sep:.6f} Ha")
    log.append("")
    print(f"E(Li) = {E_li:.6f}, E(H) = {E_h:.6f}, E_sep = {E_sep:.6f}")

    # --- PES scan: both methods ---
    results_no: list = []
    results_orth: list = []

    for R in R_POINTS:
        print(f"\n{'='*60}")
        print(f"R = {R:.3f} bohr")
        print(f"{'='*60}")

        print("\n  [no orthogonalization]")
        r_no = run_single_R(R, NMAX, orthogonalize=False, E_sep=E_sep)
        results_no.append(r_no)
        print(f"  E_mol={r_no['E_mol']:.6f}  D_e_raw={r_no['D_e_raw']:.4f}  "
              f"D_e_CP={r_no['D_e_cp']:.4f}  BSSE={r_no['BSSE']:.4f}  "
              f"cond(S)={r_no['cond_S']:.1f}  t={r_no['t_mol']:.1f}s")

        print("\n  [orthogonalize=True]")
        r_orth = run_single_R(R, NMAX, orthogonalize=True, E_sep=E_sep)
        results_orth.append(r_orth)
        print(f"  E_mol={r_orth['E_mol']:.6f}  D_e_raw={r_orth['D_e_raw']:.4f}  "
              f"D_e_CP={r_orth['D_e_cp']:.4f}  BSSE={r_orth['BSSE']:.4f}  "
              f"cond(S)={r_orth['cond_S']:.1f}  t={r_orth['t_mol']:.1f}s")

    # --- Summary tables ---
    log.append("=" * 80)
    log.append("NON-ORTHOGONALIZED BASIS (orthogonalize=False)")
    log.append("=" * 80)
    log.append(f"{'R':>6} {'E_mol':>12} {'D_e_raw':>10} {'D_e_CP':>10} "
               f"{'BSSE':>10} {'cond(S)':>10} {'min_eig(S)':>12} {'t(s)':>8}")
    for r in results_no:
        log.append(f"{r['R']:>6.3f} {r['E_mol']:>12.6f} {r['D_e_raw']:>10.4f} "
                   f"{r['D_e_cp']:>10.4f} {r['BSSE']:>10.4f} "
                   f"{r['cond_S']:>10.1f} {r['min_eigval_S']:>12.6f} "
                   f"{r['t_mol']:>8.1f}")

    log.append("")
    log.append("=" * 80)
    log.append("LOWDIN ORTHOGONALIZED BASIS (orthogonalize=True)")
    log.append("=" * 80)
    log.append(f"{'R':>6} {'E_mol':>12} {'D_e_raw':>10} {'D_e_CP':>10} "
               f"{'BSSE':>10} {'cond(S)':>10} {'min_eig(S)':>12} {'t(s)':>8}")
    for r in results_orth:
        log.append(f"{r['R']:>6.3f} {r['E_mol']:>12.6f} {r['D_e_raw']:>10.4f} "
                   f"{r['D_e_cp']:>10.4f} {r['BSSE']:>10.4f} "
                   f"{r['cond_S']:>10.1f} {r['min_eigval_S']:>12.6f} "
                   f"{r['t_mol']:>8.1f}")

    # --- Comparison table ---
    log.append("")
    log.append("=" * 80)
    log.append("COMPARISON: Delta(orthogonalized - non-orthogonalized)")
    log.append("=" * 80)
    log.append(f"{'R':>6} {'dE_mol':>12} {'dD_e_raw':>10} {'dD_e_CP':>10} "
               f"{'dBSSE':>10} {'|BSSE| reduction':>18}")
    for r_no, r_orth in zip(results_no, results_orth):
        dE = r_orth['E_mol'] - r_no['E_mol']
        dDraw = r_orth['D_e_raw'] - r_no['D_e_raw']
        dDcp = r_orth['D_e_cp'] - r_no['D_e_cp']
        dBSSE = r_orth['BSSE'] - r_no['BSSE']
        bsse_reduction = (1.0 - abs(r_orth['BSSE']) / max(abs(r_no['BSSE']), 1e-16)) * 100
        log.append(f"{r_no['R']:>6.3f} {dE:>12.6f} {dDraw:>10.4f} "
                   f"{dDcp:>10.4f} {dBSSE:>10.4f} {bsse_reduction:>17.1f}%")

    # --- Print summary ---
    print("\n" + "\n".join(log[-len(R_POINTS)-5:]))

    # --- Write output ---
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nOutput saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
