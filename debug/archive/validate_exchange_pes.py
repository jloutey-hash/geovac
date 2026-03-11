"""
LiH PES scan: cross-atom V_ee with exchange (v0.9.13).

Computes LiH energy at 8 geometry points with cross-atom J and
Mulliken exchange K for s-s orbital pairs.

Reports: E_mol, D_e_raw, D_e_CP, BSSE at each R.

Date: 2026-03-09
Version: v0.9.13 (s-s cross-atom J+K)
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

# v0.9.13 baseline (s-s cross-atom J+K only, CP-corrected)
BASELINE_CP = {
    2.5: 0.1708,
    3.0: 0.1111,
    3.015: 0.1098,
    3.5: 0.0762,
    4.0: 0.0534,
    5.0: 0.0227,
    6.0: 0.0002,
    8.0: -0.0001,
}
BASELINE_E_MOL = {
    2.5: -8.177808,
    3.0: -8.118061,
    3.015: -8.116739,
    3.5: -8.083173,
}

OUTPUT_FILE = os.path.join(os.path.dirname(__file__),
                           "data", "lih_exchange_pes.txt")


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


def run_single_R(R: float, nmax: int, E_sep: float) -> dict:
    """Run LiH at one R value, return energies and timing."""
    t0 = time.perf_counter()

    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=N_ELECTRONS,
        n_bridges=10, vee_method='slater_full',
        fci_method='auto',
    )

    # Count exchange ERIs
    nA = mol._n_spatial_A
    n_exchange = 0
    for (a, b, c, d) in mol._eri:
        a_on_A = a < nA
        b_on_A = b < nA
        if a_on_A != b_on_A and a == d and b == c:
            n_exchange += 1

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
        'n_exchange': n_exchange,
        't_mol': t_mol,
        't_bsse': t_bsse,
    }


def main() -> None:
    log: list = []

    header = "LiH PES: Cross-Atom V_ee v0.9.13 (s-s J+K)"
    log.append(header)
    log.append("Date: 2026-03-09")
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

    # --- PES scan ---
    results: list = []

    for R in R_POINTS:
        print(f"\n{'='*60}")
        print(f"R = {R:.3f} bohr")
        print(f"{'='*60}")

        r = run_single_R(R, NMAX, E_sep)
        results.append(r)

        # Compare to baseline
        bl_cp = BASELINE_CP.get(R, None)
        delta_cp = f"  delta_CP={r['D_e_cp'] - bl_cp:.4f}" if bl_cp else ""

        print(f"  E_mol={r['E_mol']:.6f}  D_e_raw={r['D_e_raw']:.4f}  "
              f"D_e_CP={r['D_e_cp']:.4f}  BSSE={r['BSSE']:.4f}  "
              f"K_entries={r['n_exchange']}  t={r['t_mol']:.1f}s{delta_cp}")

    # --- Summary table ---
    log.append("=" * 90)
    log.append("v0.9.13 CROSS-ATOM V_ee (s-s J+K)")
    log.append("=" * 90)
    log.append(f"{'R':>6} {'E_mol':>12} {'D_e_raw':>10} {'D_e_CP':>10} "
               f"{'BSSE':>10} {'K_entries':>10} {'t(s)':>8}")
    for r in results:
        log.append(f"{r['R']:>6.3f} {r['E_mol']:>12.6f} {r['D_e_raw']:>10.4f} "
                   f"{r['D_e_cp']:>10.4f} {r['BSSE']:>10.4f} "
                   f"{r['n_exchange']:>10} {r['t_mol']:>8.1f}")

    # --- Comparison to baseline ---
    log.append("")
    log.append("=" * 90)
    log.append("COMPARISON TO v0.9.13 BASELINE (s-s only)")
    log.append("=" * 90)
    log.append(f"{'R':>6} {'D_e_CP(new)':>12} {'D_e_CP(old)':>12} "
               f"{'delta':>10} {'% change':>10}")
    for r in results:
        bl = BASELINE_CP.get(r['R'], None)
        if bl is not None:
            delta = r['D_e_cp'] - bl
            pct = delta / bl * 100 if abs(bl) > 1e-10 else 0.0
            log.append(f"{r['R']:>6.3f} {r['D_e_cp']:>12.4f} {bl:>12.4f} "
                       f"{delta:>10.4f} {pct:>9.1f}%")

    # --- Key metrics ---
    r_eq = results[2]  # R=3.015
    log.append("")
    log.append(f"D_e_CP at R_eq = {r_eq['D_e_cp']:.4f} Ha "
               f"(expt: {D_E_EXPT} Ha, error: "
               f"{abs(r_eq['D_e_cp'] - D_E_EXPT)/D_E_EXPT*100:.1f}%)")

    # --- Print summary ---
    print("\n" + "\n".join(log[-10:]))

    # --- Write output ---
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nOutput saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
