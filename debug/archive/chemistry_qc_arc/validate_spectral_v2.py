"""Quick check: spectral accuracy at different n_basis for tails."""
import sys
sys.path.insert(0, '.')

import time
import numpy as np
from geovac.level4_multichannel import compute_adiabatic_curve_mc
from geovac.level4_spectral_angular import compute_adiabatic_curve_spectral

R = 1.4
l_max = 2
R_e_grid = np.linspace(0.5, 8.0, 50)

# FD reference
U_fd = compute_adiabatic_curve_mc(R, R_e_grid, l_max, Z=1.0, n_alpha=200)
i_min = np.argmin(U_fd)

for nb in [10, 15, 20]:
    t0 = time.perf_counter()
    U_sp = compute_adiabatic_curve_spectral(
        R, R_e_grid, l_max, Z=1.0, n_basis=nb, n_quad=100,
    )
    t_sp = time.perf_counter() - t0

    delta = np.abs(U_fd - U_sp)
    # Near-equilibrium error (central 60% of grid)
    i_lo, i_hi = len(R_e_grid)//5, 4*len(R_e_grid)//5
    delta_core = delta[i_lo:i_hi]

    print(f"n_basis={nb:2d}  dim={5*nb:3d}  time={t_sp*1000:.0f}ms  "
          f"max_err={np.max(delta):.2e}  core_err={np.max(delta_core):.2e}  "
          f"U_min_err={abs(U_fd[i_min]-U_sp[i_min]):.2e}")
