"""
LiH Energy Error Diagnostic (v0.9.13).

Systematic decomposition of D_e(CP) = 0.110 Ha vs expt 0.092 Ha (19% error).
Identifies which approximation contributes the 18 mHa overbinding.

Date: 2026-03-09
"""

import os
import sys
import copy
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from scipy.integrate import quad
from geovac.lattice_index import (
    MolecularLatticeIndex,
    LatticeIndex,
    compute_bsse_correction,
)

R = 3.015
NMAX = 3
Z_A, Z_B = 3, 1
N_ELECTRONS = 4

print("=" * 70)
print("LiH ENERGY ERROR DIAGNOSTIC — v0.9.13")
print(f"R = {R} bohr, nmax = {NMAX}, Z_A = {Z_A}, Z_B = {Z_B}")
print("=" * 70)

# ======================================================================
# 1. Isolated atom energies
# ======================================================================
print("\n--- 1. ISOLATED ATOM ENERGIES ---")

li = LatticeIndex(
    n_electrons=3, max_n=NMAX, nuclear_charge=Z_A,
    vee_method='slater_full', h1_method='exact',
)
E_li = li.compute_ground_state(n_states=1)[0][0]

h = LatticeIndex(
    n_electrons=1, max_n=NMAX, nuclear_charge=Z_B,
    vee_method='slater_full', h1_method='exact',
)
E_h = h.compute_ground_state(n_states=1)[0][0]

E_sep = E_li + E_h
print(f"E(Li, nmax={NMAX}) = {E_li:.6f} Ha")
print(f"E(H,  nmax={NMAX}) = {E_h:.6f} Ha")
print(f"E_sep = E(Li) + E(H) = {E_sep:.6f} Ha")

# ======================================================================
# 2a. Full v0.9.13 baseline
# ======================================================================
print("\n--- 2a. FULL v0.9.13 CALCULATION ---")

mol_full = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
    R=R, n_electrons=N_ELECTRONS,
    n_bridges=10, vee_method='slater_full',
    fci_method='auto',
)
eigvals_full, _ = mol_full.compute_ground_state(n_states=1)
E_full = eigvals_full[0]
D_e_raw_full = E_sep - E_full
print(f"E_mol (full) = {E_full:.6f} Ha")
print(f"D_e_raw = {D_e_raw_full:.4f} Ha")

# CP correction
bsse = compute_bsse_correction(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX, R=R,
    n_electrons_A=3, n_electrons_B=1,
    vee_method='slater_full', fci_method='auto',
)
E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
D_e_cp_full = E_ghost_sum - E_full
print(f"BSSE = {bsse['BSSE']:.4f} Ha")
print(f"D_e_CP = {D_e_cp_full:.4f} Ha")
print(f"E_A_ghost = {bsse['E_A_ghost']:.6f},  E_B_ghost = {bsse['E_B_ghost']:.6f}")

# ======================================================================
# 2b. Zero ALL cross-atom V_ee (J=0, K=0)
# ======================================================================
print("\n--- 2b. ZERO ALL CROSS-ATOM V_ee ---")

mol_noVee = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
    R=R, n_electrons=N_ELECTRONS,
    n_bridges=10, vee_method='slater_full',
    fci_method='auto',
)

# Zero out cross-atom ERIs
nA = mol_noVee._n_spatial_A
keys_to_remove = []
for key in mol_noVee._eri:
    a, b, c, d = key
    a_on_A = a < nA
    b_on_A = b < nA
    # Cross-atom: one index on A, one on B
    if a_on_A != b_on_A or (len(set([a < nA, b < nA, c < nA, d < nA])) > 1):
        keys_to_remove.append(key)
for key in keys_to_remove:
    del mol_noVee._eri[key]

# Rebuild dense ERI array and re-solve
n_sp = mol_noVee._n_spatial
mol_noVee._eri_dense = np.zeros((n_sp, n_sp, n_sp, n_sp))
for (a, b, c, d), val in mol_noVee._eri.items():
    mol_noVee._eri_dense[a, b, c, d] = val
eigvals_noVee, _ = mol_noVee.compute_ground_state(n_states=1)
E_noVee = eigvals_noVee[0]
print(f"E_mol (no cross-Vee) = {E_noVee:.6f} Ha")
print(f"D_e_raw = {E_sep - E_noVee:.4f} Ha")
print(f"delta E vs full = {E_noVee - E_full:.4f} Ha (positive = less bound)")

# ======================================================================
# 2c. Zero cross-nuclear: H on Li orbitals only
# ======================================================================
print("\n--- 2c. ZERO CROSS-NUCLEAR: H nucleus on Li orbitals ---")

mol_noXnucA = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
    R=R, n_electrons=N_ELECTRONS,
    n_bridges=10, vee_method='slater_full',
    fci_method='auto',
)
# Report cross-nuclear values before zeroing
print("  Cross-nuclear values (before zeroing):")
for i, (ni, li, mi) in enumerate(mol_noXnucA._li_A.lattice.states):
    if li == 0:
        v_orig = mol_noXnucA._h1_diag[i] - (-float(Z_A)**2 / (2*ni**2))
        print(f"    Li orbital (n={ni},l={li}): V_cross = {v_orig:.6f} Ha")
for j, (nj, lj, mj) in enumerate(mol_noXnucA._li_B.lattice.states):
    if lj == 0:
        v_orig = mol_noXnucA._h1_diag[nA + j] - (-float(Z_B)**2 / (2*nj**2))
        print(f"    H  orbital (n={nj},l={lj}): V_cross = {v_orig:.6f} Ha")

# Zero H-nucleus attraction on Li orbitals
for i, (ni, li, mi) in enumerate(mol_noXnucA._li_A.lattice.states):
    if li == 0:
        v_cross = mol_noXnucA._h1_diag[i] - (-float(Z_A)**2 / (2*ni**2))
        mol_noXnucA._h1_diag[i] -= v_cross
# Rebuild H1 sparse matrix
from scipy.sparse import diags
H1_offdiag = mol_noXnucA.kinetic_scale * (-mol_noXnucA._adjacency_combined)
mol_noXnucA._H1_spatial = (diags(mol_noXnucA._h1_diag) + H1_offdiag).tocsr()
eigvals_noXnucA, _ = mol_noXnucA.compute_ground_state(n_states=1)
E_noXnucA = eigvals_noXnucA[0]
print(f"E_mol (no H->Li cross-nuc) = {E_noXnucA:.6f} Ha")
print(f"delta E vs full = {E_noXnucA - E_full:.4f} Ha")

# ======================================================================
# 2d. Zero ALL cross-nuclear attraction
# ======================================================================
print("\n--- 2d. ZERO ALL CROSS-NUCLEAR ---")

mol_noXnuc = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
    R=R, n_electrons=N_ELECTRONS,
    n_bridges=10, vee_method='slater_full',
    fci_method='auto',
)
# Zero all cross-nuclear
for i, (ni, li, mi) in enumerate(mol_noXnuc._li_A.lattice.states):
    if li == 0:
        v_cross = mol_noXnuc._h1_diag[i] - (-float(Z_A)**2 / (2*ni**2))
        mol_noXnuc._h1_diag[i] -= v_cross
for j, (nj, lj, mj) in enumerate(mol_noXnuc._li_B.lattice.states):
    if lj == 0:
        v_cross = mol_noXnuc._h1_diag[nA + j] - (-float(Z_B)**2 / (2*nj**2))
        mol_noXnuc._h1_diag[nA + j] -= v_cross
mol_noXnuc._H1_spatial = (diags(mol_noXnuc._h1_diag) +
    mol_noXnuc.kinetic_scale * (-mol_noXnuc._adjacency_combined)).tocsr()
eigvals_noXnuc, _ = mol_noXnuc.compute_ground_state(n_states=1)
E_noXnuc = eigvals_noXnuc[0]
print(f"E_mol (no cross-nuclear) = {E_noXnuc:.6f} Ha")
print(f"D_e_raw = {E_sep - E_noXnuc:.4f} Ha")
print(f"delta E vs full = {E_noXnuc - E_full:.4f} Ha")

# ======================================================================
# 2e. Zero bridge/off-diagonal connections
# ======================================================================
print("\n--- 2e. ZERO BRIDGE CONNECTIONS ---")

mol_noBridge = MolecularLatticeIndex(
    Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
    R=R, n_electrons=N_ELECTRONS,
    n_bridges=0,  # Zero bridges directly
    vee_method='slater_full',
    fci_method='auto',
)
eigvals_noBridge, _ = mol_noBridge.compute_ground_state(n_states=1)
E_noBridge = eigvals_noBridge[0]
print(f"E_mol (no bridges) = {E_noBridge:.6f} Ha")
print(f"D_e_raw = {E_sep - E_noBridge:.4f} Ha")
print(f"delta E vs full = {E_noBridge - E_full:.4f} Ha")

# ======================================================================
# 3. Verify _fourier_cross_attraction against analytical formula
# ======================================================================
print("\n--- 3. VERIFY _fourier_cross_attraction ---")

# Analytical: <1s_A(Z)| -Z_B/|r-R_B| |1s_A(Z)> for s-orbital
# = -(Z_B/R)(1 - e^{-2ZR}(1 + ZR))   [exact for 1s]
# Wait — the standard formula for <1s(Z)| 1/|r-R| |1s(Z)> uses Z as scaling:
# <1s| -Z_B/|r-R| |1s> = -Z_B × [ (1/R) - exp(-2ZR)(1/R + Z) ]
# Let me verify from first principles using Gauss's law result.

def analytical_1s_cross_nuclear(Z_self: float, Z_other: float, R: float) -> float:
    """
    Exact nuclear attraction integral for 1s orbital:
    <1s(Z)| -Z_other/|r-R| |1s(Z)>

    Using electrostatic potential of |R_{1s}(r)|² = 4Z³ exp(-2Zr):
    V(R) = -Z_other × [(1/R) ∫₀^R 4Z³ e^{-2Zr} r² dr + ∫_R^∞ 4Z³ e^{-2Zr} r dr]

    Analytical inner integral: ∫₀^R 4Z³ e^{-2Zr} r² dr
     = 1 - e^{-2ZR}(1 + 2ZR + 2Z²R²)

    Analytical outer integral: ∫_R^∞ 4Z³ e^{-2Zr} r dr
     = e^{-2ZR}(2Z²R + 2Z) / (2Z)²  ... let me compute properly
    """
    # Just compute analytically from the known result:
    # Potential of 1s density = (1/r)(1 - e^{-2Zr}(1 + 2Zr)) + Z*e^{-2Zr}
    # Wait, let me use the standard electrostatic potential:
    # Φ(R) = (1/R)(1 - e^{-2ZR}(1 + 2ZR + 2Z²R²)) + (2Z²R + 2Z)e^{-2ZR}/(4Z²)

    # Easier: use the known formula. For H-like 1s with nuclear charge Z:
    # V(R) = (1/R)(1 - (1 + ZR) e^{-2ZR}) ... NO that's wrong dimension.

    # Let me just compute from first principles using exact integrals.
    # |R_1s|² = 4Z³ exp(-2Zr)
    # inner = ∫₀^R 4Z³ e^{-2Zr} r² dr
    # Using integration by parts or standard formula:
    # ∫₀^R r² e^{-ar} dr = -(1/a)R²e^{-aR} - (2/a²)Re^{-aR} - (2/a³)e^{-aR} + 2/a³
    # with a = 2Z:
    a = 2.0 * Z_self
    eaR = np.exp(-a * R)
    inner = 4 * Z_self**3 * (
        2.0/a**3 - eaR * (R**2/a + 2*R/a**2 + 2/a**3)
    )
    # outer = ∫_R^∞ 4Z³ e^{-2Zr} r dr
    # ∫_R^∞ r e^{-ar} dr = eaR * (R/a + 1/a²)
    outer = 4 * Z_self**3 * eaR * (R/a + 1.0/a**2)

    return -Z_other * (inner / R + outer)


# Test 1: H 1s (Z=1) with Z_B=1 at R=1.4
v_computed = MolecularLatticeIndex._fourier_cross_attraction(1, 0, 1, 1, 1.4)
v_exact = analytical_1s_cross_nuclear(1.0, 1.0, 1.4)
print(f"H 1s, Z_B=1, R=1.4:")
print(f"  _fourier_cross_attraction = {v_computed:.8f} Ha")
print(f"  analytical               = {v_exact:.8f} Ha")
print(f"  difference               = {abs(v_computed - v_exact):.2e} Ha")

# Test 2: Li 1s (Z=3) with Z_B=1 at R=3.015
v_computed2 = MolecularLatticeIndex._fourier_cross_attraction(1, 0, 3, 1, R)
v_exact2 = analytical_1s_cross_nuclear(3.0, 1.0, R)
print(f"\nLi 1s, Z_B=1, R={R}:")
print(f"  _fourier_cross_attraction = {v_computed2:.8f} Ha")
print(f"  analytical               = {v_exact2:.8f} Ha")
print(f"  difference               = {abs(v_computed2 - v_exact2):.2e} Ha")

# Test 3: H 1s (Z=1) with Z_B=3 at R=3.015 (H attracted to Li nucleus)
v_computed3 = MolecularLatticeIndex._fourier_cross_attraction(1, 0, 1, 3, R)
v_exact3 = analytical_1s_cross_nuclear(1.0, 3.0, R)
print(f"\nH 1s, Z_B=3, R={R}:")
print(f"  _fourier_cross_attraction = {v_computed3:.8f} Ha")
print(f"  analytical               = {v_exact3:.8f} Ha")
print(f"  difference               = {abs(v_computed3 - v_exact3):.2e} Ha")

# Test 4: Point-charge limit (Li 1s at large R)
R_large = 100.0
v_large = MolecularLatticeIndex._fourier_cross_attraction(1, 0, 3, 1, R_large)
v_point = -1.0 / R_large
print(f"\nLi 1s, Z_B=1, R={R_large} (point charge limit):")
print(f"  _fourier_cross_attraction = {v_large:.8f} Ha")
print(f"  -Z_B/R                   = {v_point:.8f} Ha")
print(f"  difference               = {abs(v_large - v_point):.2e} Ha")

# ======================================================================
# 4. Full ERI count breakdown
# ======================================================================
print("\n--- 4. ERI COUNT BREAKDOWN ---")

nA = mol_full._n_spatial_A
n_same_A, n_same_B, n_cross_J, n_cross_K, n_other = 0, 0, 0, 0, 0

for (a, b, c, d) in mol_full._eri:
    a_A = a < nA
    b_A = b < nA
    c_A = c < nA
    d_A = d < nA

    if a_A and b_A and c_A and d_A:
        n_same_A += 1
    elif not a_A and not b_A and not c_A and not d_A:
        n_same_B += 1
    else:
        # Cross-atom: check if Coulomb (a,b,a,b) or exchange (a,b,b,a)
        if a == c and b == d:
            n_cross_J += 1
        elif a == d and b == c:
            n_cross_K += 1
        else:
            n_other += 1

print(f"Same-atom Li ERIs:  {n_same_A}")
print(f"Same-atom H ERIs:   {n_same_B}")
print(f"Cross-atom J:       {n_cross_J}")
print(f"Cross-atom K:       {n_cross_K}")
print(f"Other cross-atom:   {n_other}")
print(f"Total ERI entries:  {len(mol_full._eri)}")

# ======================================================================
# 5. Summary table
# ======================================================================
print("\n" + "=" * 70)
print("SUMMARY: ENERGY DECOMPOSITION at R=3.015")
print("=" * 70)
print(f"{'Variant':<35} {'E_mol':>12} {'D_e_raw':>10} {'delta':>10}")
print("-" * 70)
print(f"{'Full v0.9.13':<35} {E_full:>12.6f} {E_sep - E_full:>10.4f} {'---':>10}")
print(f"{'No cross-atom V_ee (J=K=0)':<35} {E_noVee:>12.6f} {E_sep - E_noVee:>10.4f} {E_noVee - E_full:>10.4f}")
print(f"{'No H->Li cross-nuclear':<35} {E_noXnucA:>12.6f} {E_sep - E_noXnucA:>10.4f} {E_noXnucA - E_full:>10.4f}")
print(f"{'No cross-nuclear (both)':<35} {E_noXnuc:>12.6f} {E_sep - E_noXnuc:>10.4f} {E_noXnuc - E_full:>10.4f}")
print(f"{'No bridges':<35} {E_noBridge:>12.6f} {E_sep - E_noBridge:>10.4f} {E_noBridge - E_full:>10.4f}")
print("-" * 70)
print(f"{'E_sep = E(Li) + E(H)':<35} {E_sep:>12.6f}")
print(f"{'D_e_CP (full)':<35} {D_e_cp_full:>12.4f}")
print(f"{'D_e (expt)':<35} {'0.0924':>12}")
print(f"{'Overbinding':<35} {D_e_cp_full - 0.0924:>12.4f}")
print("=" * 70)
print("\nNote: delta = E_variant - E_full (positive = less bound / weaker)")
