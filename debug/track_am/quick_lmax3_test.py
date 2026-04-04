"""Quick l_max=3 test at key R values with sparse R_e grid."""
import sys
sys.path.insert(0, 'c:/Users/jlout/OneDrive/Desktop/Project_Geometric')

import numpy as np
import time
from geovac.n_electron_solver import (
    solve_angular_4e_multichannel, CENTRIFUGAL_4E
)
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal


def fast_solve(R, l_max=3, n_grid=4, n_Re_pts=12):
    Z_A, Z_B = 3.0, 1.0
    R_e_angular = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 13.0])[:n_Re_pts]
    mu_vals = np.zeros(len(R_e_angular))
    t0 = time.time()
    for i, R_e in enumerate(R_e_angular):
        rho_A = (R / 2) / R_e
        rho_B = (R / 2) / R_e
        evals, _ = solve_angular_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
            s4_projection=True, symmetry='s4', n_states=1)
        mu_vals[i] = evals[0]
        print(f'  Re={R_e:.1f}: mu={mu_vals[i]:.2f} ({time.time()-t0:.0f}s)', flush=True)

    U = (mu_vals + CENTRIFUGAL_4E) / R_e_angular**2
    U_spline = CubicSpline(R_e_angular, U, extrapolate=True)
    R_e_min, R_e_max, n_Re = 0.3, 12.0, 200
    h = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_r = R_e_min + (np.arange(n_Re) + 1) * h
    V = U_spline(R_e_r)
    diag = np.ones(n_Re) / h**2 + V
    off = -0.5 * np.ones(n_Re - 1) / h**2
    evals_r, _ = eigh_tridiagonal(diag, off, select='i', select_range=(0, 0))
    E_elec = evals_r[0]
    E_total = E_elec + Z_A * Z_B / R
    return E_total, time.time() - t0


E_atoms = -7.4781 + (-0.5)  # Li + H

for R in [0.7, 1.0, 1.5, 2.0, 3.0]:
    print(f'\n=== R={R} bohr ===', flush=True)
    E, t = fast_solve(R, l_max=3, n_grid=4, n_Re_pts=12)
    D_e = E_atoms - E
    print(f'R={R}: E_total={E:.6f}, D_e={D_e:.6f}, time={t:.0f}s', flush=True)
