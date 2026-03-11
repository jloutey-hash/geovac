"""v0.9.32 Diagnostic: Löwdin vs Canonical Orthogonalization + Exact J integrals.

Runs Part 1 (overlap diagnosis) and Part 2 (J integral validation) diagnostics
for LiH at R=3.015. Output goes to debug/data/mo_fci_lih_results.txt.
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from geovac.mo_fci import MOSturmianFCI

# ============================================================================
# Part 1: Löwdin Diagnosis
# ============================================================================
print("=" * 70)
print("PART 1: LÖWDIN DIAGNOSIS")
print("=" * 70)

fci = MOSturmianFCI(Z_A=3.0, Z_B=1.0, R=3.015, nmax=3, n_electrons=4, n_grid=40)

# Use the v0.9.31 converged p0
p0 = np.sqrt(-2.0 * (-10.07))  # ≈ 4.489
print(f"p0 = {p0:.4f}")

diag = fci.diagnose_orthogonalization(p0, R=3.015)

# Print overlap matrix
S = diag['S']
n = S.shape[0]
print(f"\nOverlap matrix S ({n}x{n}):")
print("     ", end="")
for j in range(n):
    print(f"  MO{j:2d}  ", end="")
print()
for i in range(n):
    print(f"MO{i:2d} ", end="")
    for j in range(n):
        print(f" {S[i,j]:7.4f}", end="")
    print()

# Eigenvalues
print(f"\nS eigenvalues (ascending):")
for i, ev in enumerate(diag['S_evals']):
    flag = " <-- NEAR-LD" if ev < 0.01 else ""
    print(f"  s_{i} = {ev:.6f}{flag}")

# H1 diagonal before/after
print(f"\nH1 diagonal BEFORE orthogonalization (raw Sturmian identity):")
for i, val in enumerate(diag['H1_raw_diag']):
    print(f"  MO{i:2d}: {val:10.6f} Ha")

print(f"\nH1 diagonal AFTER Löwdin orthogonalization:")
for i, val in enumerate(diag['H1_lowdin_diag']):
    flag = " <-- POSITIVE!" if val > 0 else ""
    print(f"  MO{i:2d}: {val:10.6f} Ha{flag}")

print(f"\nH1 diagonal AFTER canonical orthogonalization (n_dropped={diag['n_dropped']}):")
for i, val in enumerate(diag['H1_canonical_diag']):
    flag = " <-- POSITIVE!" if val > 0 else ""
    print(f"  MO{i:2d}: {val:10.6f} Ha{flag}")

# Top off-diagonal overlaps
print(f"\nTop 5 off-diagonal overlap elements:")
for i, j, sij in diag['top_offdiag']:
    mo_i = diag['mo_results'][i]
    mo_j = diag['mo_results'][j]
    print(f"  S[{i},{j}] = {sij:.6f}  "
          f"(MO n={mo_i[0]} m={mo_i[1]} n_sph={mo_i[2]} <-> "
          f"MO n={mo_j[0]} m={mo_j[1]} n_sph={mo_j[2]})")

# Key check: is H1[1,1] negative after canonical?
h1_can = diag['H1_canonical_diag']
if len(h1_can) > 1:
    print(f"\n*** KEY CHECK: H1[1,1] = {h1_can[1]:.6f} Ha ***")
    if h1_can[1] < 0:
        print("    PASS: 2-sigma orbital has NEGATIVE H1 diagonal after canonical orth")
    else:
        print("    FAIL: 2-sigma orbital still POSITIVE")

# ============================================================================
# Part 2: Exact J Integrals
# ============================================================================
print("\n" + "=" * 70)
print("PART 2: EXACT J INTEGRALS (Poisson solve)")
print("=" * 70)

from geovac.molecular_sturmian import compute_molecular_sturmian_betas, compute_exact_j_integrals

mo_results = compute_molecular_sturmian_betas(
    Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
)
print(f"Computing exact J integrals at p0={p0:.4f}, R=3.015...")
J = compute_exact_j_integrals(mo_results, 3.0, 1.0, 3.015, p0, n_grid=30)

n_mo = J.shape[0]
print(f"\nJ integral matrix ({n_mo}x{n_mo}):")
print("     ", end="")
for j in range(n_mo):
    print(f"  MO{j:2d}  ", end="")
print()
for i in range(n_mo):
    print(f"MO{i:2d} ", end="")
    for j in range(n_mo):
        print(f" {J[i,j]:7.4f}", end="")
    print()

# Validation: J_00 should be close to 5*Z_Li/8 = 1.875 Ha
J_00_exact = 5.0 * 3.0 / 8.0
print(f"\nJ_00 = {J[0,0]:.4f} Ha (atomic Li 1s self-repulsion = {J_00_exact:.4f} Ha)")
err_j00 = abs(J[0,0] - J_00_exact) / J_00_exact * 100
print(f"Error: {err_j00:.1f}%  (tolerance: 30%)")
if err_j00 < 30:
    print("PASS: J_00 within 30% of atomic value")
else:
    print("FAIL: J_00 outside 30% tolerance — check Poisson solver")

# Ordering: J_00 > J_11
if n_mo > 1:
    print(f"\nJ_00 = {J[0,0]:.4f}, J_11 = {J[1,1]:.4f}")
    if J[0, 0] > J[1, 1]:
        print("PASS: J_00 > J_11 (core > bonding self-repulsion)")
    else:
        print("FAIL: J_00 <= J_11")

# ============================================================================
# Part 3: FCI with canonical orth + exact J
# ============================================================================
print("\n" + "=" * 70)
print("PART 3: FCI AT R=3.015")
print("=" * 70)

fci2 = MOSturmianFCI(Z_A=3.0, Z_B=1.0, R=3.015, nmax=3, n_electrons=4, n_grid=40)
E_mol, p0_star = fci2.solve(damping=0.5, max_iter=20, verbose=True)

E_atoms = -7.892
E_exact = -8.071
print(f"\nE_mol = {E_mol:.6f} Ha")
print(f"E_exact = {E_exact:.3f} Ha")
print(f"E_atoms = {E_atoms:.3f} Ha")
print(f"Error vs exact: {abs(E_mol - E_exact)/abs(E_exact)*100:.1f}%")
if E_mol < E_atoms:
    print(f"D_e_raw = {E_atoms - E_mol:.4f} Ha (expt 0.092)")
    print("BOUND: YES")
else:
    print("UNBOUND")

print("\nDone.")
