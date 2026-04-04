"""Quick test: verify 2D + l-dependent PK wiring works at l_max=2."""
import sys
import os
sys.path.insert(0, r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric')

import numpy as np
from geovac.composed_diatomic import ComposedDiatomicSolver

def test_2d_ldep_pk_lmax2():
    """Run 2D + l-dep PK at l_max=2 on sparse grid."""
    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2,
        pk_channel_mode='l_dependent',
        level4_method='variational_2d',
        verbose=True,
    )
    solver.solve_core()

    # Sparse grid to check wiring
    R_grid = np.array([2.5, 3.0, 3.5, 4.0, 5.0])
    pes = solver.scan_pes(R_grid=R_grid)
    spectro = solver.fit_spectroscopic_constants()

    R_eq = spectro['R_eq']
    print(f"\n  R_eq = {R_eq:.4f} bohr")
    print(f"  R_eq err = {abs(R_eq - 3.015)/3.015*100:.1f}%")
    print(f"  D_e = {spectro['D_e']:.6f} Ha")

    # Basic sanity: R_eq should be between 2.0 and 5.0
    assert 2.0 < R_eq < 5.0, f"R_eq={R_eq} out of range"
    assert spectro['D_e'] > 0, f"D_e={spectro['D_e']} should be positive"

if __name__ == '__main__':
    test_2d_ldep_pk_lmax2()
