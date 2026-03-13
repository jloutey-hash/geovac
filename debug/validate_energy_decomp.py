"""
Energy Decomposition by Hamiltonian Term
=========================================
Decompose E(R) into individual terms to identify which has incorrect R-dependence.

Terms:
    T       : Intra-atom kinetic (graph Laplacian off-diagonal, intra-atom)
    V_nA    : Electron-nucleus A attraction (-Z_A^2/(2n^2) for A orbitals)
    V_nB    : Electron-nucleus B attraction (-Z_B^2/(2n^2) for B orbitals)
    V_cross_A: Cross-nuclear (A electrons attracted to nucleus B)
    V_cross_B: Cross-nuclear (B electrons attracted to nucleus A)
    V_bridge: Bridge hopping (inter-atom off-diagonal kinetic)
    V_ee    : Electron-electron repulsion (residual)
    V_NN    : Nuclear-nuclear repulsion (Z_A*Z_B/R)

Date: 2026-03-12
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
from typing import Dict, List
import time

from geovac.lattice_index import MolecularLatticeIndex, LatticeIndex


def decompose_at_R(
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    R: float,
    n_electrons: int,
    fci_method: str = 'auto',
) -> Dict[str, float]:
    """Compute energy decomposition at a single R."""
    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax_A, nmax_B=nmax_B,
        R=R, n_electrons=n_electrons,
        vee_method='slater_full',
        fci_method=fci_method,
    )
    eigvals, eigvecs = mol.compute_ground_state(n_states=1)
    E_total = eigvals[0]
    civec = eigvecs[:, 0]

    decomp = mol.decompose_energy(civec, E_total)
    decomp['R'] = R
    return decomp


def scan_decomposition(
    name: str,
    Z_A: int, Z_B: int,
    nmax_A: int, nmax_B: int,
    n_electrons: int,
    n_electrons_A: int,
    n_electrons_B: int,
    R_values: np.ndarray,
    fci_method: str = 'auto',
) -> List[Dict[str, float]]:
    """Full decomposition scan."""

    print(f"\n{'='*80}")
    print(f"  {name} Energy Decomposition  (nmax_A={nmax_A}, nmax_B={nmax_B})")
    print(f"{'='*80}")

    # Isolated atom energies for reference
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
    print(f"  E_sep = {E_sep:.6f} Ha  (E_A={E_A:.6f}, E_B={E_B:.6f})")

    # Column headers
    terms = ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B', 'V_bridge', 'V_ee', 'V_NN', 'E_total']
    hdr = f"{'R':>6s}"
    for t in terms:
        hdr += f" {t:>10s}"
    print(f"\n{hdr}")
    print("-" * (6 + 11 * len(terms)))

    results = []
    t0 = time.time()

    for R in R_values:
        decomp = decompose_at_R(
            Z_A, Z_B, nmax_A, nmax_B, R, n_electrons, fci_method,
        )
        results.append(decomp)

        row = f"{R:6.2f}"
        for t in terms:
            row += f" {decomp[t]:10.5f}"
        print(row)

    elapsed = time.time() - t0
    print(f"\nScan completed in {elapsed:.1f}s")

    # Analysis: compute dE/dR for each term
    print(f"\n--- R-dependence analysis (change from R_max to R_min) ---")
    R_arr = np.array([d['R'] for d in results])
    i_min, i_max = 0, len(results) - 1

    print(f"  {'Term':>12s} {'R={:.1f}'.format(R_arr[i_min]):>10s} "
          f"{'R={:.1f}'.format(R_arr[i_max]):>10s} {'Delta':>10s} {'Drives':>10s}")
    print(f"  {'-'*54}")

    for t in terms:
        v_short = results[i_min][t]
        v_long = results[i_max][t]
        delta = v_short - v_long
        # Negative delta = term becomes more negative at short R = attractive
        drives = "attract" if delta < -0.01 else ("repel" if delta > 0.01 else "~flat")
        print(f"  {t:>12s} {v_short:10.5f} {v_long:10.5f} {delta:+10.5f} {drives:>10s}")

    # Identify the dominant attractive and repulsive terms
    print(f"\n--- Dominant R-dependent terms ---")
    deltas = {}
    for t in terms:
        if t == 'E_total':
            continue
        deltas[t] = results[i_min][t] - results[i_max][t]

    sorted_terms = sorted(deltas.items(), key=lambda x: x[1])
    print(f"  Most attractive (drives R inward):")
    for t, d in sorted_terms[:3]:
        if d < -0.01:
            print(f"    {t}: {d:+.5f} Ha")

    print(f"  Most repulsive (drives R outward):")
    for t, d in reversed(sorted_terms):
        if d > 0.01:
            print(f"    {t}: {d:+.5f} Ha")

    # Check: does the sum balance correctly?
    net = sum(deltas.values())
    E_delta = results[i_min]['E_total'] - results[i_max]['E_total']
    print(f"\n  Net delta (sum of terms): {net:+.5f} Ha")
    print(f"  E_total delta:            {E_delta:+.5f} Ha")
    print(f"  Consistency check:        {abs(net - E_delta):.2e} Ha")

    return results


def save_decomposition(results: List[Dict], name: str, filepath: str) -> None:
    """Save decomposition data."""
    terms = ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B', 'V_bridge', 'V_ee', 'V_NN', 'E_total']
    with open(filepath, 'w') as f:
        f.write(f"# {name} Energy Decomposition\n")
        hdr = f"# {'R':>8s}"
        for t in terms:
            hdr += f" {t:>14s}"
        f.write(hdr + "\n")
        for d in results:
            row = f"  {d['R']:8.4f}"
            for t in terms:
                row += f" {d[t]:14.8f}"
            f.write(row + "\n")
    print(f"  Saved to {filepath}")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("Energy Decomposition by Hamiltonian Term")
    print("=" * 80)

    # --- H2 (nmax=3, homonuclear, fast) ---
    R_h2 = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])

    h2_results = scan_decomposition(
        name='H2',
        Z_A=1, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=2,
        n_electrons_A=1,
        n_electrons_B=1,
        R_values=R_h2,
        fci_method='auto',
    )

    # --- LiH (nmax=3, heteronuclear, ~367k SDs) ---
    R_lih = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0])

    lih_results = scan_decomposition(
        name='LiH',
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        n_electrons=4,
        n_electrons_A=3,
        n_electrons_B=1,
        R_values=R_lih,
        fci_method='auto',
    )

    # --- Save ---
    save_decomposition(h2_results, 'H2', 'debug/data/energy_decomp_h2.txt')
    save_decomposition(lih_results, 'LiH', 'debug/data/energy_decomp_lih.txt')

    # --- Cross-molecule comparison ---
    print(f"\n{'='*80}")
    print("  DIAGNOSIS: Which term causes the PES shape error?")
    print(f"{'='*80}")

    for name, results in [('H2', h2_results), ('LiH', lih_results)]:
        print(f"\n  {name}:")
        # Find R near equilibrium vs large R
        R_arr = np.array([d['R'] for d in results])
        E_arr = np.array([d['E_total'] for d in results])
        i_eq = np.argmin(E_arr)

        print(f"    E_min at R = {R_arr[i_eq]:.2f} bohr (E = {E_arr[i_eq]:.6f} Ha)")

        # For each term, show % of total binding energy it contributes
        d_eq = results[i_eq]
        d_far = results[-1]
        D_e = d_far['E_total'] - d_eq['E_total']

        if D_e > 0:
            print(f"    D_e = {D_e:.6f} Ha (binding energy)")
            terms_no_E = ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B',
                          'V_bridge', 'V_ee', 'V_NN']
            for t in terms_no_E:
                contrib = d_far[t] - d_eq[t]
                pct = contrib / D_e * 100
                print(f"      {t:>12s}: {contrib:+.5f} Ha ({pct:+.1f}% of D_e)")
        else:
            print(f"    No minimum in range — monotonically bound")
