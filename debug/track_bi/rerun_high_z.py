"""Rerun high-Z cases with adjusted grid parameters to check if error is numerical or structural."""

import sys
sys.path.insert(0, r"c:\Users\jlout\OneDrive\Desktop\Project_Geometric")

import numpy as np
from geovac.hyperspherical_radial import solve_helium

EXACT_ENERGIES = {
    6: -32.4062,
    7: -44.7814,
    9: -75.5316,
}

# The issue might be R_max too small. For high Z, the wavefunction is compact
# so R_max can be smaller, BUT the adiabatic potential needs enough R range
# to approach the asymptotic limit. Also try more angular grid points.

configs = [
    # Z, R_max, N_R_radial, n_alpha, label
    (6, 10.0, 3000, 150, "Z6_Rmax10_nr3k_na150"),
    (6, 15.0, 3000, 150, "Z6_Rmax15_nr3k_na150"),
    (7, 10.0, 3000, 150, "Z7_Rmax10_nr3k_na150"),
    (7, 15.0, 3000, 150, "Z7_Rmax15_nr3k_na150"),
    (9, 10.0, 3000, 150, "Z9_Rmax10_nr3k_na150"),
    (9, 15.0, 3000, 150, "Z9_Rmax15_nr3k_na150"),
    (9, 20.0, 3000, 150, "Z9_Rmax20_nr3k_na150"),
]

for Z, R_max, N_R, n_alpha, label in configs:
    print(f"\n--- {label} ---")
    try:
        result = solve_helium(Z=Z, l_max=2, verbose=True, R_max=R_max,
                              N_R_radial=N_R, n_alpha=n_alpha)
        E = result['energy']
        E_exact = EXACT_ENERGIES[Z]
        err_pct = abs((E - E_exact) / E_exact) * 100
        print(f"  E = {E:.6f}, E_exact = {E_exact:.4f}, err = {err_pct:.4f}%")
    except Exception as e:
        print(f"  FAILED: {e}")
