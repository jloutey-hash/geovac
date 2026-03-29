"""PES validation for algebraic m!=0 radial solver: H2+ 1pi_u state (m=1).

Scans R from 1.0 to 10.0 bohr, comparing matrix_method='algebraic' vs 'quadrature'
with radial_method='spectral', n_basis=20.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice

R_values = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]

results = []
for R in R_values:
    sp_alg = ProlateSpheroidalLattice(
        R=R, Z_A=1, Z_B=1,
        m=1, n_angular=0, n_radial=0,
        radial_method='spectral', n_basis=20,
        matrix_method='algebraic'
    )
    E_alg = sp_alg.total_energy()

    sp_quad = ProlateSpheroidalLattice(
        R=R, Z_A=1, Z_B=1,
        m=1, n_angular=0, n_radial=0,
        radial_method='spectral', n_basis=20,
        matrix_method='quadrature'
    )
    E_quad = sp_quad.total_energy()

    disc = E_alg - E_quad
    results.append((R, E_alg, E_quad, disc))
    print(f"R={R:5.1f}  E_alg={E_alg:12.8f}  E_quad={E_quad:12.8f}  disc={disc:+12.2e}")

# Find R_eq for each method
R_arr = np.array([r[0] for r in results])
E_alg_arr = np.array([r[1] for r in results])
E_quad_arr = np.array([r[2] for r in results])

idx_alg = np.argmin(E_alg_arr)
idx_quad = np.argmin(E_quad_arr)

max_disc = np.max(np.abs([r[3] for r in results]))

print("\n" + "="*70)
print(f"R_eq (algebraic):   {R_arr[idx_alg]:.1f} bohr  (E = {E_alg_arr[idx_alg]:.8f} Ha)")
print(f"R_eq (quadrature):  {R_arr[idx_quad]:.1f} bohr  (E = {E_quad_arr[idx_quad]:.8f} Ha)")
print(f"Max |discrepancy|:  {max_disc:.2e} Ha")
print("="*70)

# Save to file
outpath = os.path.join(os.path.dirname(__file__), 'data', 'track_j_pes_validation.txt')
with open(outpath, 'w') as f:
    f.write("H2+ 1pi_u (m=1) PES Validation: algebraic vs quadrature spectral Laguerre\n")
    f.write(f"radial_method='spectral', n_basis=20, n_angular=0, n_radial=0\n")
    f.write("="*70 + "\n\n")
    f.write(f"{'R (bohr)':>10s}  {'E_algebraic (Ha)':>18s}  {'E_quadrature (Ha)':>18s}  {'Discrepancy (Ha)':>18s}\n")
    f.write("-"*70 + "\n")
    for R, E_a, E_q, d in results:
        f.write(f"{R:10.1f}  {E_a:18.10f}  {E_q:18.10f}  {d:+18.2e}\n")
    f.write("\n")
    f.write(f"R_eq (algebraic):   {R_arr[idx_alg]:.1f} bohr  (E = {E_alg_arr[idx_alg]:.10f} Ha)\n")
    f.write(f"R_eq (quadrature):  {R_arr[idx_quad]:.1f} bohr  (E = {E_quad_arr[idx_quad]:.10f} Ha)\n")
    f.write(f"Max |discrepancy|:  {max_disc:.2e} Ha\n")

print(f"\nResults saved to {outpath}")
