"""
Lithium Hydride (LiH) - Heteronuclear Molecular Bond
=====================================================

First test of the Geometric Torsion engine on a heteronuclear molecule.

System:
    Li (Z=3) at origin  +  H (Z=1) at 3.015 Bohr
    4 electrons total (3 from Li, 1 from H)
    Bond length: 3.015 Bohr = 1.596 Angstrom

Physics:
    Li nucleus is a topological defect: gamma = 0.25 * (3 - 2) = 0.25
    H nucleus is flat: Z=1 < Z_ref=2, no torsion.

    For heteronuclear molecules we use apply_molecular_torsion() which
    applies ONLY Law 3 (per-atom torsion) without global kinetic or
    potential rescaling. Each atom keeps its natural -Z/n^2 potential.

    4-electron energy: E = 2*E_0 + 2*E_1 (fill lowest 2 orbitals)

Target:
    Total energy: -8.07 Ha (Hartree-Fock reference)

Date: February 15, 2026
"""

import numpy as np
import time
import sys

sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian, CALIBRATED_KINETIC_SCALE


def run_lih(bond_length: float = 3.015,
            max_n: int = 5,
            verbose: bool = True) -> dict:
    """
    Simulate Lithium Hydride at a given bond length.

    Uses mean-field single-particle orbitals with 4-electron filling:
    E_total = 2*E_0 + 2*E_1 (2 electrons per orbital, 2 orbitals).

    Parameters:
    -----------
    bond_length : float
        Li-H distance in Bohr (default: 3.015, equilibrium)
    max_n : int
        Basis set size (default: 5)
    verbose : bool
        Print details

    Returns:
    --------
    dict with energy, orbital energies, timing
    """
    Z_Li = 3
    Z_H = 1
    E_target = -8.07  # Hartree-Fock reference

    if verbose:
        print(f"\n{'='*70}")
        print(f"LITHIUM HYDRIDE (LiH) - HETERONUCLEAR BOND")
        print(f"{'='*70}")
        print(f"\n  Configuration:")
        print(f"    Li (Z=3) at (0, 0, 0)")
        print(f"    H  (Z=1) at ({bond_length}, 0, 0)")
        print(f"    Bond length: {bond_length:.3f} Bohr ({bond_length * 0.529177:.3f} Angstrom)")
        print(f"    Electrons:   4 (3 from Li + 1 from H)")
        print(f"\n  Method:")
        print(f"    Torsion:  Li gamma=0.25 (defect), H gamma=0.00 (flat)")
        print(f"    Solver:   Mean-field + 4-electron filling (2*E0 + 2*E1)")
        print(f"    max_n:    {max_n}")
        print(f"    Target:   {E_target:.2f} Ha")

    t0 = time.time()

    # Build with default kinetic scale and actual nuclear charges
    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0), (bond_length, 0.0, 0.0)],
        nuclear_charges=[Z_Li, Z_H],
        max_n=max_n
    )

    # Apply per-atom torsion only (no global kinetic/potential rescaling)
    mol.apply_molecular_torsion()

    t_build = time.time()

    # Solve for single-particle orbitals
    energies, wf = mol.compute_ground_state(n_states=6, method='mean_field')
    t_solve = time.time()

    # 4-electron filling: 2 per orbital, fill lowest 2
    E_total = 2 * energies[0] + 2 * energies[1]
    error_pct = 100 * (E_total - E_target) / abs(E_target)

    if verbose:
        print(f"\n{'='*70}")
        print(f"  RESULTS")
        print(f"{'='*70}")
        print(f"\n  Single-particle orbitals:")
        labels = ['1s(Li)', '2s(Li)', '2p(Li)', '2p(Li)', '1s(H)', '3s(Li)']
        for i in range(min(6, len(energies))):
            occ = '**' if i < 2 else '  '
            lbl = labels[i] if i < len(labels) else '...'
            print(f"    {occ} E_{i} = {energies[i]:.6f} Ha  ({lbl})")
        print(f"\n  4-electron energy: 2*E_0 + 2*E_1 = {E_total:.6f} Ha")
        print(f"  Target:            {E_target:.2f} Ha")
        print(f"  Error:             {error_pct:+.2f}%")
        print(f"\n  Orbital gap (HOMO-LUMO): {energies[2] - energies[1]:.4f} Ha "
              f"= {(energies[2] - energies[1]) * 27.2114:.2f} eV")
        print(f"  States: {mol.n_total_states} single-particle")
        print(f"  Build: {(t_build - t0)*1000:.0f} ms, Solve: {(t_solve - t_build)*1000:.0f} ms")

    return {
        'energy': E_total,
        'target': E_target,
        'error_pct': error_pct,
        'bond_length': bond_length,
        'orbital_energies': energies,
        'n_states': mol.n_total_states,
        'time': t_solve - t0
    }


def scan_bond_length(max_n: int = 5,
                     verbose: bool = True) -> dict:
    """
    Scan LiH bond length to find the potential energy curve.
    """
    distances = np.arange(1.5, 8.0, 0.5)

    if verbose:
        print(f"\n{'='*70}")
        print(f"LiH BOND LENGTH SCAN")
        print(f"{'='*70}")
        print(f"  max_n: {max_n}")
        print(f"  Distances: {distances[0]:.1f} to {distances[-1]:.1f} Bohr")
        print(f"\n  {'R (Bohr)':>10}  {'E_total (Ha)':>12}  {'Note':>15}")
        print(f"  {'-'*10}  {'-'*12}  {'-'*15}")

    results = {}
    E_min = 0.0
    R_min = 0.0

    for R in distances:
        r = run_lih(bond_length=R, max_n=max_n, verbose=False)
        results[R] = r

        note = ""
        if r['energy'] < E_min:
            E_min = r['energy']
            R_min = R
            note = "<-- minimum"

        if verbose:
            print(f"  {R:10.3f}  {r['energy']:12.6f}  {note:>15}")

    if verbose:
        print(f"\n  Equilibrium:    R = {R_min:.1f} Bohr, E = {E_min:.6f} Ha")
        print(f"  Experimental:   R = 3.015 Bohr")
        # Compute binding energy relative to separated atoms
        E_inf = results[max(distances)]['energy']
        binding = E_min - E_inf
        print(f"  Binding energy: {binding:.4f} Ha = {binding * 27.2114:.2f} eV")

    return results


def compare_torsion(bond_length: float = 3.015,
                    max_n: int = 5,
                    verbose: bool = True) -> None:
    """
    Compare LiH with and without torsion on Li.
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"TORSION COMPARISON: LiH at R = {bond_length} Bohr")
        print(f"{'='*70}")

    # Without torsion
    mol_flat = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0), (bond_length, 0.0, 0.0)],
        nuclear_charges=[3, 1],
        max_n=max_n
    )
    e_flat, _ = mol_flat.compute_ground_state(n_states=4, method='mean_field')
    E_flat = 2 * e_flat[0] + 2 * e_flat[1]

    # With torsion on Li
    mol_torsion = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0), (bond_length, 0.0, 0.0)],
        nuclear_charges=[3, 1],
        max_n=max_n
    )
    mol_torsion.apply_molecular_torsion()
    e_torsion, _ = mol_torsion.compute_ground_state(n_states=4, method='mean_field')
    E_torsion = 2 * e_torsion[0] + 2 * e_torsion[1]

    E_target = -8.07

    if verbose:
        print(f"\n  {'':>20}  {'Flat (gamma=0)':>16}  {'Torsion (Li 0.25)':>18}")
        print(f"  {'-'*20}  {'-'*16}  {'-'*18}")
        print(f"  {'E_0 (1s core)':>20}  {e_flat[0]:>16.6f}  {e_torsion[0]:>18.6f}")
        print(f"  {'E_1 (2s bond)':>20}  {e_flat[1]:>16.6f}  {e_torsion[1]:>18.6f}")
        print(f"  {'E_2 (LUMO)':>20}  {e_flat[2]:>16.6f}  {e_torsion[2]:>18.6f}")
        print(f"  {'4e total':>20}  {E_flat:>16.6f}  {E_torsion:>18.6f}")
        print(f"  {'Target':>20}  {E_target:>16.2f}  {E_target:>18.2f}")
        err_flat = 100 * (E_flat - E_target) / abs(E_target)
        err_torsion = 100 * (E_torsion - E_target) / abs(E_target)
        print(f"  {'Error':>20}  {err_flat:>+15.2f}%  {err_torsion:>+17.2f}%")
        print(f"\n  Torsion correction: {E_torsion - E_flat:+.4f} Ha "
              f"= {(E_torsion - E_flat) * 27.2114:+.2f} eV")


if __name__ == '__main__':
    print("=" * 70)
    print("GeoVac: Lithium Hydride (LiH) Demo")
    print("First heteronuclear test of the Geometric Torsion engine")
    print("=" * 70)

    # Torsion comparison at equilibrium
    compare_torsion(bond_length=3.015, max_n=5)

    # Main calculation
    result = run_lih(bond_length=3.015, max_n=5, verbose=True)

    # Bond length scan
    scan_bond_length(max_n=5, verbose=True)

    print(f"\n{'='*70}")
    print(f"'A twisted lithium bonds to a flat hydrogen.'")
    print(f"{'='*70}")
