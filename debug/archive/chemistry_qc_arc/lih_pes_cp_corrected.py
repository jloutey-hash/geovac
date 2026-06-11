"""
CP-corrected PES for LiH via Boys-Bernardi counterpoise.

Computes E_mol(R), E_ghost_sum(R), and CP-corrected binding energy D_e(R)
at several internuclear distances.

Output: debug/data/lih_pes_cp_corrected.txt
"""
import warnings
import numpy as np

warnings.filterwarnings('ignore')

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
    compute_bsse_correction,
)

NMAX = 3
VEE = 'slater_full'
FCI = 'auto'

# Own-basis atomic energies (R-independent)
li_own = LatticeIndex(
    n_electrons=3, max_n=NMAX, nuclear_charge=3,
    vee_method=VEE, h1_method='exact', fci_method=FCI,
)
E_li_own = li_own.compute_ground_state(n_states=1)[0][0]

h_own = LatticeIndex(
    n_electrons=1, max_n=NMAX, nuclear_charge=1,
    vee_method=VEE, h1_method='exact', fci_method=FCI,
)
E_h_own = h_own.compute_ground_state(n_states=1)[0][0]

E_sep_own = E_li_own + E_h_own
print(f"E(Li, own) = {E_li_own:.6f} Ha")
print(f"E(H,  own) = {E_h_own:.6f} Ha")
print(f"E_sep(own) = {E_sep_own:.6f} Ha")
print()

R_values = [2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0, 8.0]
header = f"{'R':>6s}  {'E_mol':>10s}  {'E_ghost':>10s}  {'D_e_CP':>8s}  {'D_e_raw':>8s}  {'BSSE':>8s}"
print(header)
print("-" * len(header))

for R in R_values:
    # Molecular energy
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=4,
        vee_method=VEE, fci_method=FCI,
    )
    E_mol = mol.compute_ground_state(n_states=1)[0][0]

    # Ghost energies
    li_g = MolecularLatticeIndex(
        Z_A=3, Z_B=0, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=3,
        vee_method=VEE, fci_method=FCI,
    )
    E_li_g = li_g.compute_ground_state(n_states=1)[0][0]

    h_g = MolecularLatticeIndex(
        Z_A=0, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=1,
        vee_method=VEE, fci_method=FCI,
    )
    E_h_g = h_g.compute_ground_state(n_states=1)[0][0]

    E_ghost = E_li_g + E_h_g
    D_e_cp = E_ghost - E_mol
    D_e_raw = E_sep_own - E_mol
    BSSE = E_ghost - E_sep_own

    print(f"{R:6.3f}  {E_mol:10.5f}  {E_ghost:10.5f}  {D_e_cp:8.4f}  {D_e_raw:8.4f}  {BSSE:8.4f}")
