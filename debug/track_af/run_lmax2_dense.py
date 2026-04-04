"""Track AF: Dense PES at l_max=2 for 2D + l-dependent PK."""
import sys
sys.path.insert(0, r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric')

import numpy as np
from geovac.composed_diatomic import ComposedDiatomicSolver

R_grid = np.array([
    2.0, 2.3, 2.5, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
    3.5, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0,
])

solver = ComposedDiatomicSolver.LiH_ab_initio(
    l_max=2,
    pk_channel_mode='l_dependent',
    level4_method='variational_2d',
    verbose=True,
)
result = solver.run_all(R_grid=R_grid, n_Re=300)
spectro = result['spectro']
print(f"\n=== l_max=2 RESULT ===")
print(f"R_eq = {spectro['R_eq']:.4f} bohr (err {abs(spectro['R_eq']-3.015)/3.015*100:.1f}%)")
print(f"E_min = {spectro['E_min']:.6f} Ha")
print(f"D_e = {spectro['D_e']:.6f} Ha")
print(f"omega_e = {spectro['omega_e']:.1f} cm-1")

# Save
import os
OUT = r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric\debug\track_af'
pes = result['pes']
with open(os.path.join(OUT, 'pes_2d_ldep_lmax2.txt'), 'w') as f:
    f.write(f"# R_eq = {spectro['R_eq']:.4f}, err = {abs(spectro['R_eq']-3.015)/3.015*100:.1f}%\n")
    f.write(f"# E_min = {spectro['E_min']:.6f}, D_e = {spectro['D_e']:.6f}\n")
    f.write(f"# omega_e = {spectro['omega_e']:.1f} cm-1\n")
    f.write(f"# avg time/pt = {pes['wall_times'].mean():.1f}s\n")
    for i in range(len(pes['R'])):
        f.write(f"{pes['R'][i]:8.4f}  {pes['E_composed'][i]:14.8f}  {pes['wall_times'][i]:8.2f}\n")
