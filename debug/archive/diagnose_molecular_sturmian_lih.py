"""
Diagnostic: Molecular Sturmian β_k for LiH (v0.9.28).

Computes β_k table, self-consistency, PES, and counterpoise correction.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.molecular_sturmian import compute_beta_k
from geovac.lattice_index import MolecularLatticeIndex, compute_bsse_correction

Z_A, Z_B = 3, 1  # Li, H
nmax = 3
n_electrons = 4
p0_init = np.sqrt(10.0)

output_lines = []
def log(msg: str = "") -> None:
    print(msg)
    output_lines.append(msg)


# =====================================================================
# Part 1: β_k table at R=3.015, p0=√10
# =====================================================================
log("=" * 70)
log("MOLECULAR STURMIAN β_k TABLE (v0.9.28)")
log("=" * 70)
R = 3.015
p0 = p0_init

log(f"\nR = {R} bohr, p0 = {p0:.6f}, Z_A = {Z_A}, Z_B = {Z_B}")
log(f"{'n':>3} {'l':>3} {'beta_A':>12} {'eps_A':>12} {'beta_B':>12} {'eps_B':>12}")
log("-" * 60)

for n in range(1, nmax + 1):
    for l_val in range(n):
        bk_A = compute_beta_k(n, l_val, float(Z_A), float(Z_B), R, p0)
        bk_B = compute_beta_k(n, l_val, float(Z_B), float(Z_A), R, p0)
        eps_A = p0**2 / 2.0 - bk_A * float(Z_A) * p0 / n
        eps_B = p0**2 / 2.0 - bk_B * float(Z_B) * p0 / n
        log(f"{n:3d} {l_val:3d} {bk_A:12.6f} {eps_A:12.6f} {bk_B:12.6f} {eps_B:12.6f}")


# =====================================================================
# Part 2: Checks 1-4
# =====================================================================
log("\n" + "=" * 70)
log("NUMERICAL CHECKS")
log("=" * 70)

# Check 1: Li 1s
bk = compute_beta_k(1, 0, float(Z_A), float(Z_B), R, p0)
eps = p0**2 / 2.0 - bk * float(Z_A) * p0
log(f"\nCheck 1 (Li 1s):  β = {bk:.6f}, ε = {eps:.6f} Ha")

# Check 2: H 1s (KEY)
bk_H = compute_beta_k(1, 0, float(Z_B), float(Z_A), R, p0)
eps_H = p0**2 / 2.0 - bk_H * float(Z_B) * p0
eps_atomic_H = p0**2 / 2.0 - float(Z_B) * p0
log(f"Check 2 (H 1s):  β = {bk_H:.6f}, ε_mol = {eps_H:.6f} Ha, ε_atm = {eps_atomic_H:.6f} Ha")
log(f"  >>> H 1s diagonal is {'NEGATIVE (bound)' if eps_H < 0 else 'POSITIVE (still unbound)'}")
log(f"  >>> Improvement: {eps_atomic_H - eps_H:.4f} Ha")

# Check 3: R=100 backward compat
worst_err = 0.0
for n in range(1, 4):
    for l_val in range(n):
        bk_A = compute_beta_k(n, l_val, float(Z_A), float(Z_B), 100.0, p0)
        bk_B = compute_beta_k(n, l_val, float(Z_B), float(Z_A), 100.0, p0)
        worst_err = max(worst_err, abs(bk_A - 1.0), abs(bk_B - 1.0))
log(f"Check 3 (R=100): worst β deviation = {worst_err:.6f} (target < 0.05)")

# Check 4: all β ≥ 1
all_ge_1 = True
for n in range(1, 4):
    for l_val in range(n):
        if compute_beta_k(n, l_val, float(Z_A), float(Z_B), R, p0) < 1.0:
            all_ge_1 = False
        if compute_beta_k(n, l_val, float(Z_B), float(Z_A), R, p0) < 1.0:
            all_ge_1 = False
log(f"Check 4 (β≥1):  {'PASS' if all_ge_1 else 'FAIL'}")


# =====================================================================
# Part 3: Self-consistency at R=3.015
# =====================================================================
log("\n" + "=" * 70)
log("SELF-CONSISTENCY AT R=3.015")
log("=" * 70)

mol = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
    R=R, n_electrons=n_electrons, n_bridges=1,
    vee_method='slater', fci_method='auto',
    use_sturmian='molecular',
)
p0_star, E_mol, n_iter, converged = mol.solve_sturmian_p0(
    R=R, damping=0.5, max_iter=30, tol=1e-5
)
E_atoms = -7.892  # Li(-7.432) + H(-0.5) ≈ -7.932 (FCI); we use -7.892 from v0.9.21

log(f"\np0* = {p0_star:.6f}")
log(f"E_mol = {E_mol:.6f} Ha")
log(f"Converged = {converged} in {n_iter} iterations")
log(f"E_atoms = {E_atoms:.3f} Ha")
log(f"D_e(raw) = {E_atoms - E_mol:.6f} Ha")
log(f">>> LiH is {'BOUND' if E_mol < E_atoms else 'UNBOUND'}")


# =====================================================================
# Part 4: PES scan
# =====================================================================
log("\n" + "=" * 70)
log("POTENTIAL ENERGY SURFACE")
log("=" * 70)

R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
log(f"\n{'R':>6} {'E_mol':>12} {'D_e_raw':>10} {'D_e_CP':>10}")
log("-" * 42)

for R_val in R_values:
    try:
        mol_r = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
            R=R_val, n_electrons=n_electrons, n_bridges=1,
            vee_method='slater', fci_method='auto',
            use_sturmian='molecular',
        )
        p0_r, E_r, n_it, conv = mol_r.solve_sturmian_p0(
            R=R_val, damping=0.5, max_iter=30, tol=1e-5
        )

        D_e_raw = E_atoms - E_r

        # Counterpoise (ghost atoms bypass Sturmian — same BSSE as standard)
        try:
            bsse_dict = compute_bsse_correction(
                Z_A=Z_A, Z_B=Z_B, nmax_A=nmax, nmax_B=nmax,
                R=R_val, n_electrons_A=Z_A, n_electrons_B=Z_B,
                vee_method='slater', fci_method='auto',
            )
            bsse = bsse_dict['BSSE']
            D_e_CP = D_e_raw - bsse
        except Exception as e2:
            bsse = float('nan')
            D_e_CP = float('nan')

        log(f"{R_val:6.3f} {E_r:12.6f} {D_e_raw:10.6f} {D_e_CP:10.6f}")
    except Exception as e:
        log(f"{R_val:6.3f} {'ERROR':>12} -- {str(e)[:40]}")


# Save results
outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "data", "molecular_sturmian_lih_results.txt")
os.makedirs(os.path.dirname(outpath), exist_ok=True)
with open(outpath, "w") as f:
    f.write("\n".join(output_lines))
print(f"\nResults saved to {outpath}")
