"""
LiH FCI convergence: nmax=3 baseline re-validation and nmax=4 extension.

Validates that the v0.9.10 normalization fix (_phi_s_orbital_general for n>=3)
does not regress the nmax=3 result, then extends to nmax=4 per atom.

Output: debug/data/lih_nmax4_convergence.txt

Date: 2026-03-08
"""

import os
import sys
import time
import warnings
from math import comb

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)

# Ensure geovac is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac.lattice_index import (
    MolecularLatticeIndex,
    LatticeIndex,
    compute_bsse_correction,
)

R_EQ = 3.015       # experimental equilibrium, Bohr
D_E_EXPT = 0.0924  # experimental binding energy, Ha
Z_A, Z_B = 3, 1    # Li, H
N_ELECTRONS = 4
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "data", "lih_nmax4_convergence.txt")


def run_lih(nmax: int, log_lines: list) -> dict:
    """Run LiH FCI at given nmax, returning timing and energy data."""
    header = f"===== nmax = {nmax} per atom ====="
    print(f"\n{header}")
    log_lines.append(f"\n{header}")

    # Predict basis size
    n_spatial_per_atom = sum(n**2 for n in range(1, nmax + 1))
    n_spatial = 2 * n_spatial_per_atom
    n_sp = 2 * n_spatial
    n_sd = comb(n_sp, N_ELECTRONS)
    msg = f"Predicted: {n_spatial} spatial, {n_sp} spin-orb, {n_sd:,} SDs"
    print(msg)
    log_lines.append(msg)

    # --- Molecular energy ---
    t_start = time.perf_counter()
    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
        R=R_EQ, n_electrons=N_ELECTRONS,
        n_bridges=10, vee_method='slater_full',
        fci_method='auto',
    )
    t_build = time.perf_counter() - t_start

    t_solve_start = time.perf_counter()
    eigvals, eigvecs = mol.compute_ground_state(n_states=1)
    t_solve = time.perf_counter() - t_solve_start
    t_total = time.perf_counter() - t_start

    E_mol = eigvals[0]
    n_sd_actual = mol.n_sd

    # --- Separated atom energies (own basis) ---
    t_sep_start = time.perf_counter()
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
    t_sep = time.perf_counter() - t_sep_start

    E_sep = E_li + E_h
    D_e_raw = E_sep - E_mol

    # --- BSSE / counterpoise ---
    t_bsse_start = time.perf_counter()
    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax, R=R_EQ,
        n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
    )
    t_bsse = time.perf_counter() - t_bsse_start

    E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
    D_e_cp = E_ghost_sum - E_mol
    err_pct = abs(D_e_cp - D_E_EXPT) / D_E_EXPT * 100

    # --- Report ---
    results = {
        'nmax': nmax,
        'n_sd': n_sd_actual,
        'E_mol': E_mol,
        'E_li': E_li,
        'E_h': E_h,
        'D_e_raw': D_e_raw,
        'BSSE': bsse['BSSE'],
        'D_e_cp': D_e_cp,
        'err_pct': err_pct,
        't_build': t_build,
        't_solve': t_solve,
        't_total': t_total,
        't_bsse': t_bsse,
    }

    lines = [
        f"NSD            = {n_sd_actual:,}",
        f"E_mol          = {E_mol:.6f} Ha",
        f"E(Li)          = {E_li:.6f} Ha",
        f"E(H)           = {E_h:.6f} Ha",
        f"E(Li)+E(H)     = {E_sep:.6f} Ha",
        f"D_e (raw)      = {D_e_raw:.4f} Ha",
        f"BSSE           = {bsse['BSSE']:.4f} Ha  (Li: {bsse['BSSE_A']:.4f}, H: {bsse['BSSE_B']:.4f})",
        f"D_e (CP)       = {D_e_cp:.4f} Ha",
        f"D_e (expt)     = {D_E_EXPT:.4f} Ha",
        f"Error vs expt  = {err_pct:.1f}%",
        f"",
        f"Timing:",
        f"  Build (init)   = {t_build:.1f}s",
        f"  Solve (eigen)  = {t_solve:.1f}s",
        f"  Total (mol)    = {t_total:.1f}s",
        f"  BSSE (2 ghost) = {t_bsse:.1f}s",
    ]

    if t_solve > t_build:
        lines.append(f"  Bottleneck: eigensolver ({t_solve:.0f}s > {t_build:.0f}s build)")
    else:
        lines.append(f"  Bottleneck: assembly ({t_build:.0f}s > {t_solve:.0f}s solve)")

    for line in lines:
        print(line)
        log_lines.append(line)

    return results


def main() -> None:
    log_lines: list = []
    header = "LiH FCI Convergence Study — v0.9.10 (cross-atom J fix)"
    log_lines.append(header)
    log_lines.append(f"Date: 2026-03-08")
    log_lines.append(f"R = {R_EQ} Bohr, Z_A={Z_A} (Li), Z_B={Z_B} (H), Ne={N_ELECTRONS}")
    log_lines.append(f"D_e(expt) = {D_E_EXPT} Ha")
    print(header)

    all_results = []

    # --- nmax=2 (fast baseline) ---
    try:
        r2 = run_lih(2, log_lines)
        all_results.append(r2)
    except Exception as e:
        msg = f"nmax=2 FAILED: {e}"
        print(msg)
        log_lines.append(msg)

    # --- nmax=3 (re-validation under normalization fix) ---
    try:
        r3 = run_lih(3, log_lines)
        all_results.append(r3)
    except Exception as e:
        msg = f"nmax=3 FAILED: {e}"
        print(msg)
        log_lines.append(msg)

    # --- nmax=4 ---
    log_lines.append("\n--- nmax=4 attempt ---")
    log_lines.append("WARNING: 8.2M SDs expected. This may take >600s.")
    print("\n--- nmax=4: 8.2M SDs, may be very slow ---")
    try:
        r4 = run_lih(4, log_lines)
        all_results.append(r4)
    except MemoryError as e:
        msg = f"nmax=4 OUT OF MEMORY: {e}"
        print(msg)
        log_lines.append(msg)
    except Exception as e:
        msg = f"nmax=4 FAILED after partial run: {e}"
        print(msg)
        log_lines.append(msg)

    # --- Summary table ---
    log_lines.append("\n===== CONVERGENCE SUMMARY =====")
    log_lines.append(f"{'nmax':>5} {'NSD':>12} {'D_e_CP':>10} {'err%':>8} {'time(s)':>10}")
    for r in all_results:
        log_lines.append(
            f"{r['nmax']:>5} {r['n_sd']:>12,} {r['D_e_cp']:>10.4f} "
            f"{r['err_pct']:>7.1f}% {r['t_total']:>10.1f}"
        )
    log_lines.append(f"{'expt':>5} {'':>12} {D_E_EXPT:>10.4f}")

    summary = "\n".join(log_lines[-len(all_results)-2:])
    print(f"\n{summary}")

    # --- Write output ---
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nOutput saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
