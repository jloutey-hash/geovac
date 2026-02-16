"""
Lithium Hydride (LiH) - Dynamic Bridge Geometry Optimization
=============================================================

v0.5.0: Distance-dependent bridges enable molecular geometry optimization.

System:
    Li (Z=3) at origin  +  H (Z=1) at variable distance R
    4 electrons total (3 from Li, 1 from H)
    Experimental bond length: 3.015 Bohr = 1.596 Angstrom

Physics:
    Bridge weight decays exponentially with distance:
        W_bridge = A * exp(-lambda * R)

    Total energy = E_electronic + V_NN (nuclear repulsion)

    Competition:
        - Nuclear repulsion (Z_Li * Z_H / R) pushes atoms apart
        - Graph tunneling (exp(-lambda * R)) pulls atoms together
        - Equilibrium: minimum in E(R) = The Bond

    Li nucleus is a topological defect: gamma = 0.25 * (3 - 2) = 0.25
    H nucleus is flat: Z=1 < Z_ref=2, no torsion.

    4-electron energy: E = 2*E_0 + 2*E_1 (fill lowest 2 orbitals)

Date: February 15, 2026
"""

import numpy as np
import time
import sys

sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian, CALIBRATED_KINETIC_SCALE


def run_lih(bond_length: float = 3.015,
            max_n: int = 5,
            bridge_decay_rate: float = 1.0,
            bridge_amplitude: float = 1.0,
            verbose: bool = True) -> dict:
    """
    Simulate Lithium Hydride at a given bond length with dynamic bridges.

    Uses mean-field single-particle orbitals with 4-electron filling:
    E_total = 2*E_0 + 2*E_1 + V_NN (nuclear repulsion)

    Parameters:
    -----------
    bond_length : float
        Li-H distance in Bohr (default: 3.015, equilibrium)
    max_n : int
        Basis set size (default: 5)
    bridge_decay_rate : float
        Exponential decay rate lambda (1/Bohr) for bridge weights (default: 1.0)
    bridge_amplitude : float
        Pre-exponential factor A for bridge weights (default: 1.0)
    verbose : bool
        Print details

    Returns:
    --------
    dict with energy, orbital energies, bridge weight, timing
    """
    Z_Li = 3
    Z_H = 1

    t0 = time.time()

    # Build molecule with distance-dependent bridges
    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0), (bond_length, 0.0, 0.0)],
        nuclear_charges=[Z_Li, Z_H],
        max_n=max_n,
        bridge_amplitude=bridge_amplitude,
        bridge_decay_rate=bridge_decay_rate,
    )

    # Apply per-atom torsion only (no global kinetic/potential rescaling)
    mol.apply_molecular_torsion()

    t_build = time.time()

    # Solve for single-particle orbitals
    energies, wf = mol.compute_ground_state(n_states=6, method='mean_field')
    t_solve = time.time()

    # 4-electron filling: 2 per orbital, fill lowest 2
    E_electronic = 2 * energies[0] + 2 * energies[1]

    # Nuclear repulsion: V_NN = Z_Li * Z_H / R
    V_NN = mol.compute_nuclear_repulsion()

    # Total energy = electronic + nuclear repulsion
    E_total = E_electronic + V_NN

    # Bridge weight for this distance
    bridge_weight = mol.bridge_info[0]['bridge_weight'] if mol.bridge_info else 1.0

    if verbose:
        print(f"\n{'='*70}")
        print(f"LITHIUM HYDRIDE (LiH) - DYNAMIC BRIDGE BOND")
        print(f"{'='*70}")
        print(f"\n  Configuration:")
        print(f"    Li (Z=3) at (0, 0, 0)")
        print(f"    H  (Z=1) at ({bond_length}, 0, 0)")
        print(f"    Bond length: {bond_length:.3f} Bohr ({bond_length * 0.529177:.3f} Angstrom)")
        print(f"    Electrons:   4 (3 from Li + 1 from H)")
        print(f"\n  Bridge Physics:")
        print(f"    Decay rate:  lambda = {bridge_decay_rate:.2f} / Bohr")
        print(f"    Amplitude:   A = {bridge_amplitude:.2f}")
        print(f"    Weight:      W = {bridge_weight:.6f}  (A * exp(-lambda * R))")
        print(f"\n  Energy Decomposition:")
        print(f"    E_electronic: {E_electronic:.6f} Ha")
        print(f"    V_NN:         {V_NN:+.6f} Ha  (Z_Li * Z_H / R)")
        print(f"    E_total:      {E_total:.6f} Ha")
        print(f"\n  Orbitals:")
        labels = ['1s(Li)', '2s(Li)', '2p(Li)', '2p(Li)', '1s(H)', '3s(Li)']
        for i in range(min(6, len(energies))):
            occ = '**' if i < 2 else '  '
            lbl = labels[i] if i < len(labels) else '...'
            print(f"    {occ} E_{i} = {energies[i]:.6f} Ha  ({lbl})")
        print(f"\n  States: {mol.n_total_states} single-particle")
        print(f"  Build: {(t_build - t0)*1000:.0f} ms, Solve: {(t_solve - t_build)*1000:.0f} ms")

    return {
        'energy': E_total,
        'E_electronic': E_electronic,
        'V_NN': V_NN,
        'bridge_weight': bridge_weight,
        'bond_length': bond_length,
        'orbital_energies': energies,
        'n_states': mol.n_total_states,
        'time': t_solve - t0
    }


def scan_bond_length(max_n: int = 5,
                     bridge_decay_rate: float = 1.0,
                     bridge_amplitude: float = 1.0,
                     verbose: bool = True) -> dict:
    """
    Scan LiH bond length to find equilibrium geometry.

    The energy curve E(R) = E_electronic(R) + Z_Li*Z_H/R should show:
    - Too close: nuclear repulsion dominates -> energy rises
    - Too far: bridges break (W->0) -> energy rises
    - Just right: a minimum appears (The Bond)
    """
    distances = np.arange(2.0, 6.5, 0.25)

    if verbose:
        print(f"\n{'='*70}")
        print(f"LiH BOND LENGTH SCAN (Dynamic Bridges)")
        print(f"{'='*70}")
        print(f"  max_n:  {max_n}")
        print(f"  lambda: {bridge_decay_rate:.2f} / Bohr")
        print(f"  A:      {bridge_amplitude:.2f}")
        print(f"  Range:  {distances[0]:.1f} to {distances[-1]:.1f} Bohr")
        print(f"\n  {'R':>6}  {'W_bridge':>10}  {'E_elec':>10}  {'V_NN':>8}"
              f"  {'E_total':>10}  {'Note':>12}")
        print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*10}  {'-'*12}")

    results = {}
    E_min = 0.0
    R_min = 0.0

    for R in distances:
        r = run_lih(
            bond_length=R, max_n=max_n,
            bridge_decay_rate=bridge_decay_rate,
            bridge_amplitude=bridge_amplitude,
            verbose=False,
        )
        results[R] = r

        note = ""
        if r['energy'] < E_min:
            E_min = r['energy']
            R_min = R
            note = "<-- min"

        if verbose:
            print(f"  {R:6.2f}  {r['bridge_weight']:10.6f}  {r['E_electronic']:10.4f}"
                  f"  {r['V_NN']:8.4f}  {r['energy']:10.4f}  {note:>12}")

    if verbose:
        print(f"\n  {'='*70}")
        print(f"  Equilibrium:    R = {R_min:.2f} Bohr ({R_min * 0.529177:.3f} Angstrom)")
        print(f"  Experimental:   R = 3.015 Bohr (1.596 Angstrom)")
        print(f"  E_min:          {E_min:.6f} Ha")

        # Binding energy relative to separated atoms
        E_inf = results[max(distances)]['energy']
        binding = E_min - E_inf
        if binding < 0:
            print(f"  Binding energy: {binding:.4f} Ha = {binding * 27.2114:.2f} eV")
            print(f"  --> BOUND STATE FOUND!")
        else:
            print(f"  Binding energy: {binding:.4f} Ha (no minimum found)")

    return results


if __name__ == '__main__':
    print("=" * 70)
    print("GeoVac v0.5.0: LiH Dynamic Bridge Geometry Optimization")
    print("Finding the bond from Nuclear Repulsion (1/R) vs")
    print("Graph Tunneling (exp(-lambda*R))")
    print("=" * 70)

    # Calibrated bridge parameters for LiH
    # A=8.5, lambda=0.2 places the PES minimum near R~2.75 Bohr
    # with binding energy ~0.10 Ha (~2.7 eV)
    A_cal = 8.5
    lam_cal = 0.2

    # Single-point at experimental geometry
    result = run_lih(
        bond_length=3.015, max_n=5,
        bridge_amplitude=A_cal, bridge_decay_rate=lam_cal,
        verbose=True,
    )

    # Bond length scan: find the equilibrium
    scan_bond_length(
        max_n=5,
        bridge_amplitude=A_cal, bridge_decay_rate=lam_cal,
        verbose=True,
    )

    print(f"\n{'='*70}")
    print(f"'The bond emerges from the competition between")
    print(f" repulsion and tunneling on the graph.'")
    print(f"{'='*70}")
