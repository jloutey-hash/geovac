"""
Track AD Part 2: Compute ab initio R_ref candidates from core screening.

Natural R_ref candidates:
(a) Core radius: <r> of 1s core orbital
(b) PK Gaussian width: r_PK = 1/sqrt(B) from AbInitioPK
(c) Screening length: distance where Z_eff transitions

All derived from CoreScreening with zero fitting.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK

# Solve core for Li (Z=3)
print("Solving Li 1s^2 core...")
core = CoreScreening(Z=3, l_max=2, n_alpha=200)
core.solve(verbose=True)

print(f"\nCore energy: {core.energy:.6f} Ha")

# Derive PK parameters
pk = AbInitioPK(core, n_core=2)
print(f"\n{pk.summary()}")

# Candidate R_ref values
print(f"\n{'='*60}")
print("R_ref CANDIDATES (all from core screening, zero fitting)")
print(f"{'='*60}")

# (a) Core radius variants
r_grid = core.r_grid
n_r = core.n_r

# <r> = integral r * n(r) dr / integral n(r) dr
r_avg = np.trapezoid(r_grid * n_r, r_grid) / np.trapezoid(n_r, r_grid)
print(f"\n(a1) <r> (mean core radius): {r_avg:.4f} bohr")

# <r^2>^{1/2} = RMS radius
r2_avg = np.trapezoid(r_grid**2 * n_r, r_grid) / np.trapezoid(n_r, r_grid)
r_rms = np.sqrt(r2_avg)
print(f"(a2) sqrt(<r^2>) (RMS core radius): {r_rms:.4f} bohr")

# Median radius (where N_core = 1)
N_core = core.N_core
from scipy.interpolate import CubicSpline
N_spline = CubicSpline(r_grid, N_core)
from scipy.optimize import brentq
idx = np.searchsorted(N_core, 1.0)
r_median = brentq(lambda x: N_spline(x) - 1.0, r_grid[max(idx-1, 1)], r_grid[min(idx+1, len(r_grid)-1)])
print(f"(a3) r_median (N_core=1): {r_median:.4f} bohr")

# Peak of r^2 * n(r)
radial_density = r_grid**2 * n_r
r_peak = r_grid[np.argmax(radial_density)]
print(f"(a4) r_peak (max of r^2*n(r)): {r_peak:.4f} bohr")

# (b) From PK Gaussian width parameter B
print(f"\n(b) PK Gaussian radii from B candidates:")
for name, info in pk.B_candidates.items():
    r_pk = 1.0 / np.sqrt(info['B'])
    print(f"  r_PK({name}) = 1/sqrt(B) = {r_pk:.4f} bohr  (B={info['B']:.4f})")

# Selected B (inv2)
r_pk_selected = 1.0 / np.sqrt(pk.B)
print(f"  r_PK(selected={pk.B_method}) = {r_pk_selected:.4f} bohr")

# (c) Screening length: where Z_eff transitions from Z to Z-2
# dZ_eff/dr is the negative density; find where it's steepest
r_test = np.linspace(0.05, 5.0, 500)
z_eff_vals = core.z_eff(r_test)
dz_dr = np.gradient(z_eff_vals, r_test)
i_steepest = np.argmin(dz_dr)  # most negative derivative
r_steepest = r_test[i_steepest]
print(f"\n(c1) r_steepest (max |dZ_eff/dr|): {r_steepest:.4f} bohr")

# Distance where Z_eff = Z - 1 (half screening)
z_half = core.Z - 1.0
from scipy.optimize import brentq as bq
z_spline = CubicSpline(r_test, z_eff_vals - z_half)
# Find zero crossing
sign_changes = np.where(np.diff(np.sign(z_eff_vals - z_half)))[0]
if len(sign_changes) > 0:
    r_half_screen = bq(lambda x: CubicSpline(r_test, z_eff_vals)(x) - z_half,
                        r_test[sign_changes[0]], r_test[sign_changes[0]+1])
    print(f"(c2) r_half_screen (Z_eff = Z-1): {r_half_screen:.4f} bohr")
else:
    r_half_screen = r_avg
    print(f"(c2) r_half_screen: NOT FOUND (Z_eff never reaches Z-1)")

# r where Z_eff = Z - 1.5 (75% screening)
z_75 = core.Z - 1.5
sign_changes_75 = np.where(np.diff(np.sign(z_eff_vals - z_75)))[0]
if len(sign_changes_75) > 0:
    r_75_screen = bq(lambda x: CubicSpline(r_test, z_eff_vals)(x) - z_75,
                      r_test[sign_changes_75[0]], r_test[sign_changes_75[0]+1])
    print(f"(c3) r_75_screen (Z_eff = Z-1.5): {r_75_screen:.4f} bohr")
else:
    r_75_screen = None
    print(f"(c3) r_75_screen: NOT FOUND")

# Summary table
print(f"\n\n{'='*60}")
print("SUMMARY OF R_ref CANDIDATES")
print(f"{'='*60}")
candidates = {
    'r_avg': r_avg,
    'r_rms': r_rms,
    'r_median': r_median,
    'r_peak': r_peak,
    'r_pk_inv2': r_pk_selected,
    'r_steepest': r_steepest,
    'r_half_screen': r_half_screen,
}
for name, val in candidates.items():
    print(f"  {name:>20s} = {val:.4f} bohr")

# Also compute the empirical R_ref that gives ~2.0% at l_max=4
# w_PK(R) = delta_{l,0} * min(cap, R/R_ref)
# From CLAUDE.md: empirical R_ref achieves 2.0% at l_max=4
print(f"\n\nNote: The w_PK(R) = delta_{{l,0}} * min(cap, R/R_ref) form")
print(f"with empirical R_ref achieves 2.0% R_eq at l_max=4.")
print(f"These candidates should be tested at l_max=2,3,4 to find")
print(f"which (if any) achieve <3% without fitting.")

# Save
outfile = os.path.join(os.path.dirname(__file__), 'rref_candidates.txt')
with open(outfile, 'w') as f:
    f.write("name,value_bohr\n")
    for name, val in candidates.items():
        f.write(f"{name},{val:.6f}\n")
print(f"\nSaved to {outfile}")
