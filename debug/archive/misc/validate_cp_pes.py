"""
Counterpoise-Corrected PES for H2 and LiH
==========================================
Tests whether the BSSE-free PES has correct equilibrium geometry and force constant.

Method (Boys-Bernardi Counterpoise):
    E_CP(R) = E_AB(R) - BSSE(R)

where BSSE(R) = [E_A^mol(R) - E_A] + [E_B^mol(R) - E_B]
    E_A^mol(R) = atom A in molecular basis (ghost B at distance R)
    E_A        = atom A in its own basis (isolated)

Success criterion:
    If k_CP ~ k_Morse and R_min^CP ~ R_eq => graph structure is sound, BSSE is the issue.
    If k_CP is still wrong => fundamental limitation in graph structure.

Date: 2026-03-12
"""

import warnings
warnings.filterwarnings('ignore')  # Suppress all warnings (scipy IntegrationWarning, etc.)

import numpy as np
from typing import Dict, List, Tuple
import time

from geovac.lattice_index import (
    MolecularLatticeIndex,
    LatticeIndex,
    compute_bsse_correction,
)


# ============================================================
# Experimental / reference data
# ============================================================
# H2: R_eq = 1.401 bohr, D_e = 0.1745 Ha, omega_e = 4401 cm^-1
# Morse: k = 2 * D_e * a^2, a = omega_e * sqrt(mu / (2*D_e))
# For H2: mu = 0.5 amu = 918.076 a.u., omega_e = 0.02005 Ha
# k_H2 = mu * omega_e^2 = 918.076 * 0.02005^2 = 0.369 Ha/bohr^2
H2_REF = {
    'R_eq': 1.401,       # bohr
    'D_e': 0.1745,       # Ha
    'k': 0.369,          # Ha/bohr^2 (harmonic force constant)
}

# LiH: R_eq = 3.015 bohr, D_e = 0.092 Ha, omega_e = 1406 cm^-1
# mu = m_Li * m_H / (m_Li + m_H) = 7*1/(7+1) = 0.875 amu = 1605.6 a.u.
# omega_e = 1406 cm^-1 = 0.006407 Ha
# k = mu * omega_e^2 = 1605.6 * 0.006407^2 = 0.0659 Ha/bohr^2
LIH_REF = {
    'R_eq': 3.015,       # bohr
    'D_e': 0.092,        # Ha
    'k': 0.0659,         # Ha/bohr^2
}


def compute_pes_point(
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    R: float,
    n_electrons: int,
    n_electrons_A: int,
    n_electrons_B: int,
    fci_method: str = 'auto',
) -> Dict[str, float]:
    """Compute raw and CP-corrected energy at a single R."""

    # --- Molecular FCI ---
    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax_A, nmax_B=nmax_B,
        R=R, n_electrons=n_electrons,
        vee_method='slater_full',
        fci_method=fci_method,
    )
    E_mol = mol.compute_ground_state(n_states=1)[0][0]

    # --- BSSE correction at this R ---
    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax_A, nmax_B=nmax_B,
        R=R,
        n_electrons_A=n_electrons_A,
        n_electrons_B=n_electrons_B,
        vee_method='slater_full',
        fci_method=fci_method,
    )

    E_CP = E_mol - bsse['BSSE']  # Remove artificial BSSE lowering

    return {
        'R': R,
        'E_mol': E_mol,
        'E_CP': E_CP,
        'BSSE': bsse['BSSE'],
        'BSSE_A': bsse['BSSE_A'],
        'BSSE_B': bsse['BSSE_B'],
        'E_A_ghost': bsse['E_A_ghost'],
        'E_B_ghost': bsse['E_B_ghost'],
        'E_A_own': bsse['E_A_own'],
        'E_B_own': bsse['E_B_own'],
    }


def fit_parabola(R_vals: np.ndarray, E_vals: np.ndarray) -> Tuple[float, float, float]:
    """
    Fit E(R) = a*(R - R0)^2 + E0 near the minimum.
    Returns (R0, E0, k) where k = 2*a is the force constant.
    """
    # Use 2nd-order polynomial: E = a*R^2 + b*R + c
    coeffs = np.polyfit(R_vals, E_vals, 2)
    a, b, c = coeffs
    R0 = -b / (2 * a)
    E0 = a * R0**2 + b * R0 + c
    k = 2 * a  # Force constant (d^2E/dR^2)
    return R0, E0, k


def scan_molecule(
    name: str,
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    n_electrons_A: int,
    n_electrons_B: int,
    R_values: np.ndarray,
    ref: Dict[str, float],
    fci_method: str = 'auto',
) -> Dict:
    """Full PES scan with CP correction for one molecule."""

    print(f"\n{'='*70}")
    print(f"  {name}  (nmax_A={nmax_A}, nmax_B={nmax_B})")
    print(f"{'='*70}")

    # --- Isolated atom energies ---
    print("Computing isolated atoms...")
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
    print(f"  E({name[0]}) = {E_A:.6f} Ha")
    print(f"  E({name[-1]}) = {E_B:.6f} Ha")
    print(f"  E_sep     = {E_sep:.6f} Ha")

    # --- PES scan ---
    print(f"\nScanning {len(R_values)} R points...")
    print(f"{'R':>8s} {'E_mol':>12s} {'E_CP':>12s} {'BSSE':>10s} "
          f"{'D_e_raw':>10s} {'D_e_CP':>10s}")
    print("-" * 70)

    results = []
    t0 = time.time()

    for i, R in enumerate(R_values):
        pt = compute_pes_point(
            Z_A, Z_B, nmax_A, nmax_B, R,
            n_electrons, n_electrons_A, n_electrons_B,
            fci_method=fci_method,
        )
        results.append(pt)

        D_e_raw = E_sep - pt['E_mol']
        D_e_CP = E_sep - pt['E_CP']

        print(f"{R:8.3f} {pt['E_mol']:12.6f} {pt['E_CP']:12.6f} "
              f"{pt['BSSE']:10.6f} {D_e_raw:10.6f} {D_e_CP:10.6f}")

    elapsed = time.time() - t0
    print(f"\nScan completed in {elapsed:.1f}s")

    # --- Extract arrays ---
    R_arr = np.array([r['R'] for r in results])
    E_mol_arr = np.array([r['E_mol'] for r in results])
    E_CP_arr = np.array([r['E_CP'] for r in results])
    BSSE_arr = np.array([r['BSSE'] for r in results])

    # --- Find minima ---
    idx_raw = np.argmin(E_mol_arr)
    idx_cp = np.argmin(E_CP_arr)

    # --- Fit parabola near minimum (use 5 points around min) ---
    def fit_near_min(idx: int, E_arr: np.ndarray) -> Tuple[float, float, float]:
        lo = max(0, idx - 2)
        hi = min(len(R_arr), idx + 3)
        if hi - lo < 3:
            lo = max(0, hi - 3)
        return fit_parabola(R_arr[lo:hi], E_arr[lo:hi])

    R0_raw, E0_raw, k_raw = fit_near_min(idx_raw, E_mol_arr)
    R0_cp, E0_cp, k_cp = fit_near_min(idx_cp, E_CP_arr)

    D_e_raw = E_sep - E0_raw
    D_e_cp = E_sep - E0_cp

    # --- Report ---
    print(f"\n--- {name} Results ---")
    print(f"{'':20s} {'Raw':>12s} {'CP-corr':>12s} {'Expt':>12s}")
    print(f"{'R_eq (bohr)':20s} {R0_raw:12.4f} {R0_cp:12.4f} {ref['R_eq']:12.4f}")
    print(f"{'D_e (Ha)':20s} {D_e_raw:12.6f} {D_e_cp:12.6f} {ref['D_e']:12.6f}")
    print(f"{'k (Ha/bohr^2)':20s} {k_raw:12.6f} {k_cp:12.6f} {ref['k']:12.6f}")
    print()
    print(f"  R_eq error (raw):  {abs(R0_raw - ref['R_eq'])/ref['R_eq']*100:.1f}%")
    print(f"  R_eq error (CP):   {abs(R0_cp - ref['R_eq'])/ref['R_eq']*100:.1f}%")
    print(f"  D_e error (raw):   {abs(D_e_raw - ref['D_e'])/ref['D_e']*100:.1f}%")
    print(f"  D_e error (CP):    {abs(D_e_cp - ref['D_e'])/ref['D_e']*100:.1f}%")
    print(f"  k error (raw):     {abs(k_raw - ref['k'])/ref['k']*100:.1f}%")
    print(f"  k error (CP):      {abs(k_cp - ref['k'])/ref['k']*100:.1f}%")

    # --- BSSE R-dependence ---
    print(f"\n--- BSSE R-dependence ---")
    print(f"  BSSE range: [{BSSE_arr.min():.6f}, {BSSE_arr.max():.6f}] Ha")
    print(f"  BSSE variation: {BSSE_arr.max() - BSSE_arr.min():.6f} Ha")
    if BSSE_arr.max() - BSSE_arr.min() > 0.01:
        print(f"  ** BSSE is R-dependent — CP correction reshapes PES **")
    else:
        print(f"  BSSE is approximately R-independent — CP is a rigid shift")

    return {
        'name': name,
        'R': R_arr,
        'E_mol': E_mol_arr,
        'E_CP': E_CP_arr,
        'BSSE': BSSE_arr,
        'R0_raw': R0_raw, 'R0_cp': R0_cp,
        'D_e_raw': D_e_raw, 'D_e_cp': D_e_cp,
        'k_raw': k_raw, 'k_cp': k_cp,
        'E_sep': E_sep,
        'ref': ref,
        'results': results,
    }


def save_results(data: Dict, filepath: str) -> None:
    """Save PES data to text file."""
    with open(filepath, 'w') as f:
        f.write(f"# {data['name']} Counterpoise-Corrected PES\n")
        f.write(f"# E_sep = {data['E_sep']:.8f} Ha\n")
        f.write(f"# R0_raw = {data['R0_raw']:.4f}, R0_CP = {data['R0_cp']:.4f}, "
                f"R_eq(expt) = {data['ref']['R_eq']:.4f} bohr\n")
        f.write(f"# D_e_raw = {data['D_e_raw']:.6f}, D_e_CP = {data['D_e_cp']:.6f}, "
                f"D_e(expt) = {data['ref']['D_e']:.6f} Ha\n")
        f.write(f"# k_raw = {data['k_raw']:.6f}, k_CP = {data['k_cp']:.6f}, "
                f"k(expt) = {data['ref']['k']:.6f} Ha/bohr^2\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_mol':>14s} {'E_CP':>14s} {'BSSE':>14s} "
                f"{'D_e_raw':>14s} {'D_e_CP':>14s}\n")

        for i in range(len(data['R'])):
            R = data['R'][i]
            E_mol = data['E_mol'][i]
            E_CP = data['E_CP'][i]
            bsse = data['BSSE'][i]
            D_raw = data['E_sep'] - E_mol
            D_cp = data['E_sep'] - E_CP
            f.write(f"  {R:8.4f} {E_mol:14.8f} {E_CP:14.8f} {bsse:14.8f} "
                    f"{D_raw:14.8f} {D_cp:14.8f}\n")

    print(f"  Saved to {filepath}")


def print_verdict(h2_data: Dict, lih_data: Dict) -> None:
    """Print final diagnostic verdict."""
    print(f"\n{'='*70}")
    print("  VERDICT: Is the graph structure fundamentally sound?")
    print(f"{'='*70}\n")

    for data in [h2_data, lih_data]:
        name = data['name']
        ref = data['ref']

        R_err = abs(data['R0_cp'] - ref['R_eq']) / ref['R_eq'] * 100
        D_err = abs(data['D_e_cp'] - ref['D_e']) / ref['D_e'] * 100
        k_err = abs(data['k_cp'] - ref['k']) / ref['k'] * 100

        print(f"  {name}:")
        print(f"    R_eq:  {data['R0_cp']:.3f} vs {ref['R_eq']:.3f} bohr "
              f"({R_err:.1f}% error)")
        print(f"    D_e:   {data['D_e_cp']:.4f} vs {ref['D_e']:.4f} Ha "
              f"({D_err:.1f}% error)")
        print(f"    k:     {data['k_cp']:.4f} vs {ref['k']:.4f} Ha/bohr^2 "
              f"({k_err:.1f}% error)")

        if R_err < 10 and k_err < 30:
            print(f"    => Graph structure OK. Residual error from basis truncation.\n")
        elif R_err < 10:
            print(f"    => R_eq OK but k wrong. PES shape needs work.\n")
        else:
            print(f"    => R_eq WRONG. Fundamental structural limitation.\n")

    # Overall
    h2_Rok = abs(h2_data['R0_cp'] - H2_REF['R_eq']) / H2_REF['R_eq'] < 0.10
    lih_Rok = abs(lih_data['R0_cp'] - LIH_REF['R_eq']) / LIH_REF['R_eq'] < 0.10

    if h2_Rok and lih_Rok:
        print("  CONCLUSION: Both molecules have correct R_eq after CP correction.")
        print("  The graph structure is fundamentally sound. BSSE was the problem.")
    elif h2_Rok and not lih_Rok:
        print("  CONCLUSION: H2 is correct but LiH R_eq remains shifted.")
        print("  Graph structure works for homonuclear; heteronuclear has a")
        print("  structural limitation (likely cross-nuclear attraction model).")
    else:
        print("  CONCLUSION: Neither molecule has correct R_eq after CP correction.")
        print("  The graph structure itself has fundamental limitations beyond BSSE.")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("Counterpoise-Corrected PES Validation")
    print("=" * 70)

    # --- H2 PES (nmax=3: fast, 9 spatial orbs, C(18,2)=153 SDs) ---
    R_h2 = np.array([
        0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0,
    ])

    h2_data = scan_molecule(
        name='H2',
        Z_A=1, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=2,
        n_electrons_A=1,
        n_electrons_B=1,
        R_values=R_h2,
        ref=H2_REF,
        fci_method='auto',
    )

    # --- LiH PES (nmax=3: 28 spatial, 56 spin-orbs, ~367k SDs) ---
    R_lih = np.array([
        2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0,
    ])

    lih_data = scan_molecule(
        name='LiH',
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=4,
        n_electrons_A=3,
        n_electrons_B=1,
        R_values=R_lih,
        ref=LIH_REF,
        fci_method='auto',
    )

    # --- Save data ---
    save_results(h2_data, 'debug/data/cp_pes_h2.txt')
    save_results(lih_data, 'debug/data/cp_pes_lih.txt')

    # --- Verdict ---
    print_verdict(h2_data, lih_data)
