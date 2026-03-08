"""
LiH Potential Energy Surface — GeoVac MolecularLatticeIndex FCI

Sweeps R from 2.0 to 8.0 Bohr at nmax=3.
Same-atom V_ee approximation (cross-atom ERIs deferred).

Experimental references:
    R_eq = 3.015 Bohr, E(LiH exact NR) = -8.0705 Ha
    D_e = 0.0924 Ha (2.515 eV)
"""

import sys
import numpy as np
import time

sys.path.insert(0, '.')
from geovac.lattice_index import MolecularLatticeIndex, LatticeIndex

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# PES sweep parameters
R_values = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
NMAX = 3
N_BRIDGES = 20
results = []

print(f"LiH PES Sweep: nmax={NMAX}, n_bridges={N_BRIDGES}, vee=slater_full")
print(f"R values: {R_values}")
print("=" * 70)

for R in R_values:
    t0 = time.time()
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=4,
        n_bridges=N_BRIDGES, vee_method='slater_full',
        fci_method='auto',
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    E = eigvals[0]
    dt = time.time() - t0
    results.append((R, E, mol.V_NN, E - mol.V_NN))
    print(f"R={R:.2f}: E_total={E:.6f} Ha  E_elec={E-mol.V_NN:.6f}  "
          f"V_NN={mol.V_NN:.4f}  time={dt:.1f}s")
    sys.stdout.flush()

# Compute atomic reference energies at same nmax
print("\n" + "=" * 70)
print("Atomic reference energies (nmax=3):")

li = LatticeIndex(n_electrons=3, max_n=NMAX, nuclear_charge=3,
                  vee_method='slater_full', h1_method='exact')
E_li, _ = li.compute_ground_state(n_states=1)
print(f"  E(Li, nmax={NMAX}) = {E_li[0]:.6f} Ha")

h = LatticeIndex(n_electrons=1, max_n=NMAX, nuclear_charge=1,
                 vee_method='slater_full', h1_method='exact')
E_h, _ = h.compute_ground_state(n_states=1)
print(f"  E(H, nmax={NMAX}) = {E_h[0]:.6f} Ha")
E_sep = E_li[0] + E_h[0]
print(f"  E(Li+H) = {E_sep:.6f} Ha")

# Analysis
R_arr = np.array([r for r, e, v, ee in results])
E_arr = np.array([e for r, e, v, ee in results])
idx = np.argmin(E_arr)

print("\n" + "=" * 70)
print("RESULTS SUMMARY")
print("=" * 70)
print(f"Equilibrium: R_eq = {R_arr[idx]:.2f} Bohr, E = {E_arr[idx]:.6f} Ha")
print(f"Experimental: R_eq = 3.015 Bohr, E = -8.0705 Ha")
err_R = abs(R_arr[idx] - 3.015) / 3.015 * 100
err_E = abs((E_arr[idx] - (-8.0705)) / (-8.0705)) * 100
print(f"R error: {err_R:.1f}%")
print(f"E error: {err_E:.2f}%")

D_e = E_arr[-1] - E_arr[idx]
D_e_from_atoms = E_sep - E_arr[idx]
print(f"\nBinding energy (from PES endpoints): D_e = {D_e:.6f} Ha ({D_e*27.2114:.3f} eV)")
print(f"Binding energy (from atomic): D_e = {D_e_from_atoms:.6f} Ha ({D_e_from_atoms*27.2114:.3f} eV)")
print(f"Experimental: D_e = 0.0924 Ha (2.515 eV)")

# Save results
with open('debug/data/lih_pes.txt', 'w') as f:
    f.write(f"# LiH PES: GeoVac MolecularLatticeIndex FCI\n")
    f.write(f"# Z_A=3 (Li), Z_B=1 (H), nmax={NMAX}, n_bridges={N_BRIDGES}\n")
    f.write(f"# vee_method=slater_full, same-atom V_ee approximation\n")
    f.write(f"# Equilibrium: R_eq={R_arr[idx]:.2f} Bohr, E={E_arr[idx]:.6f} Ha\n")
    f.write(f"# Expt: R_eq=3.015 Bohr, E=-8.0705 Ha\n")
    f.write(f"# E(Li,nmax={NMAX})={E_li[0]:.6f}, E(H,nmax={NMAX})={E_h[0]:.6f}\n")
    f.write(f"#\n")
    f.write(f"# R(Bohr)  E_total(Ha)  V_NN(Ha)  E_elec(Ha)\n")
    for R, E, V, Ee in results:
        f.write(f"{R:.2f}  {E:.6f}  {V:.6f}  {Ee:.6f}\n")
print("\nSaved to debug/data/lih_pes.txt")
