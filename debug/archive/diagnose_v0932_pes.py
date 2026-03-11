"""v0.9.32 PES scan and full diagnostics for LiH MO Sturmian FCI."""
import sys
import numpy as np
sys.path.insert(0, '.')
from geovac.mo_fci import MOSturmianFCI

R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
E_atoms = -7.892  # Li + H
E_exact = -8.071  # exact LiH at R_eq

results = []
for R in R_values:
    print(f"\n--- R = {R:.3f} bohr ---")
    fci = MOSturmianFCI(Z_A=3.0, Z_B=1.0, R=R, nmax=3, n_electrons=4, n_grid=40)
    E, p0 = fci.solve(damping=0.5, max_iter=20, verbose=False)
    D_e = E_atoms - E
    print(f"  E = {E:.6f} Ha, p0* = {p0:.4f}, D_e = {D_e:.4f} Ha")
    results.append((R, E, p0, D_e))

print("\n" + "=" * 70)
print("PES SUMMARY")
print("=" * 70)
print(f"{'R':>6}  {'E_mol':>12}  {'p0*':>8}  {'D_e_raw':>10}  {'Bound?':>8}")
for R, E, p0, De in results:
    bound = "YES" if De > 0 else "NO"
    print(f"{R:6.3f}  {E:12.6f}  {p0:8.4f}  {De:10.4f}  {bound:>8}")

# Find minimum
energies = [r[1] for r in results]
min_idx = np.argmin(energies)
R_eq = results[min_idx][0]
E_min = results[min_idx][1]
print(f"\nR_eq = {R_eq:.3f} bohr (expt 3.015)")
print(f"E_min = {E_min:.6f} Ha (exact -8.071)")

# Dissociation check
E_R6 = results[-1][1]
print(f"\nDissociation check:")
print(f"  E(R=6) = {E_R6:.6f}, E_atoms = {E_atoms:.3f}")
print(f"  |E(R=6) - E_atoms| = {abs(E_R6 - E_atoms):.4f} Ha (threshold 0.1)")
if abs(E_R6 - E_atoms) < 0.1:
    print("  PASS")
else:
    print("  FAIL")

# Bound minimum?
if np.any(np.array([r[3] for r in results]) > 0):
    print("\nBound minimum: YES")
else:
    print("\nBound minimum: NO")

# R_eq in [2.5, 3.5]?
if 2.5 <= R_eq <= 3.5:
    print(f"R_eq in [2.5, 3.5]: YES ({R_eq:.3f})")
else:
    print(f"R_eq in [2.5, 3.5]: NO ({R_eq:.3f})")
