"""Track AF: 2D + l-dep PK + cusp at l_max=2 only."""
import sys
sys.path.insert(0, r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric')
import os
import numpy as np
from geovac.composed_diatomic import ComposedDiatomicSolver

R_grid = np.array([
    2.0, 2.3, 2.5, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
    3.5, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0,
])
OUT = r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric\debug\track_af'

solver = ComposedDiatomicSolver.LiH_ab_initio(
    l_max=2,
    pk_channel_mode='l_dependent',
    level4_method='variational_2d',
    cusp_correction=True,
    verbose=True,
)
result = solver.run_all(R_grid=R_grid, n_Re=300)
spectro = result['spectro']
pes = result['pes']
R_eq_err = abs(spectro['R_eq'] - 3.015) / 3.015 * 100
print(f"\n=== l_max=2+cusp RESULT ===")
print(f"R_eq = {spectro['R_eq']:.4f} bohr (err {R_eq_err:.1f}%)")
print(f"E_min = {spectro['E_min']:.6f} Ha")
print(f"D_e = {spectro['D_e']:.6f} Ha")

# Grid-based R_eq
imin = np.argmin(pes['E_composed'][~np.isnan(pes['E_composed'])])
R_grid_valid = pes['R'][~np.isnan(pes['E_composed'])]
print(f"Grid R_eq = {R_grid_valid[imin]:.1f} bohr")

# Cusp correction size at R_eq
print(f"Cusp correction at grid min: {pes['delta_cusp'][imin]*1000:.3f} mHa")
print(f"Uncorrected R_eq: {pes.get('R_eq_uncorrected', 'N/A')}")

with open(os.path.join(OUT, 'pes_2d_ldep_lmax2_cusp.txt'), 'w') as f:
    f.write(f"# R_eq = {spectro['R_eq']:.4f}, err = {R_eq_err:.1f}%\n")
    f.write(f"# E_min = {spectro['E_min']:.6f}, D_e = {spectro['D_e']:.6f}\n")
    for i in range(len(pes['R'])):
        f.write(f"{pes['R'][i]:8.4f}  {pes['E_composed'][i]:14.8f}  "
                f"{pes['delta_cusp'][i]*1000:+8.3f}mHa  {pes['wall_times'][i]:8.2f}\n")
