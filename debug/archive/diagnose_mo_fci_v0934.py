"""Diagnostics for v0.9.34: Dual-p0 MO Sturmian FCI.

Runs numerical checks 1-4, combined S eigenvalues, H1_HH diagonal,
dual-p0 FCI at R=3.015, PES scan, and eleven-configuration comparison.

Output: debug/data/mo_fci_lih_results.txt

Performance note: uses n_grid=20 for speed. Full production runs use n_grid=40.
"""

from __future__ import annotations

import sys
import time
import numpy as np

sys.path.insert(0, '.')
from geovac.molecular_sturmian import (
    compute_molecular_sturmian_betas,
    compute_h1_matrix,
    compute_cross_set_integrals,
)
from geovac.mo_fci import MOSturmianFCI

Z_A, Z_B = 3.0, 1.0
R_eq = 3.015
E_exact = -8.0706  # exact LiH at R=3.015
E_atoms_exact = -7.892  # sum of exact Li + H atomic energies
N_GRID = 20  # reduced for speed (production: 40)

out_lines: list[str] = []


def log(msg: str = '') -> None:
    print(msg)
    out_lines.append(msg)


log("=" * 70)
log("v0.9.34 Dual-p0 MO Sturmian FCI — LiH Diagnostics")
log(f"n_grid={N_GRID} (reduced for speed)")
log("=" * 70)

# ==========================================================================
# Check 1: Cross-set overlap consistency at equal p0
# ==========================================================================
log("\n--- CHECK 1: Cross-set overlap at equal p0 ---")
p0_test = 2.0
mo_test = compute_molecular_sturmian_betas(
    Z_A=Z_A, Z_B=Z_B, R=R_eq, p0=p0_test, nmax=2,
)
_, _, _, S_within, _ = compute_h1_matrix(
    mo_test, Z_A, Z_B, R_eq, p0_test, n_grid=N_GRID,
    orth_method='none_raw',
)
O_cross, _ = compute_cross_set_integrals(
    mo_test, mo_test, Z_A, Z_B, R_eq, p0_test, p0_test, n_grid=N_GRID,
)
n_test = min(S_within.shape[0], O_cross.shape[0])
log(f"  p0={p0_test}, nmax=2, n_mo={n_test}")
for k in range(n_test):
    if abs(S_within[k, k]) < 1e-6:
        continue
    rel_err = abs(O_cross[k, k] - S_within[k, k]) / abs(S_within[k, k])
    log(f"  k={k}: S_within={S_within[k,k]:.6f}, O_cross={O_cross[k,k]:.6f}, "
        f"err={rel_err*100:.2f}%")
check1_pass = all(
    abs(O_cross[k, k] - S_within[k, k]) / max(abs(S_within[k, k]), 1e-10) < 0.05
    for k in range(n_test) if abs(S_within[k, k]) > 1e-6
)
log(f"  {'PASS' if check1_pass else 'FAIL'}")

# ==========================================================================
# Check 2: H-scale 1sigma H1 diagonal
# ==========================================================================
log("\n--- CHECK 2: H-scale H1 diagonal at p0_H=1.0 ---")
p0_H = 1.0
mo_H = compute_molecular_sturmian_betas(
    Z_A=Z_A, Z_B=Z_B, R=R_eq, p0=p0_H, nmax=2,
)
_, _, _, _, H1_HH_raw = compute_h1_matrix(
    mo_H, Z_A, Z_B, R_eq, p0_H, n_grid=N_GRID, orth_method='none_raw',
)
valid_H = [mo for mo in mo_H if np.isfinite(mo[4])]
n_h = len(valid_H)
log(f"  p0_H={p0_H}, nmax_H=2, n_mo_H={n_h}")
all_negative = True
for k in range(n_h):
    val = H1_HH_raw[k, k]
    label = f"(n={valid_H[k][0]},m={valid_H[k][1]},ns={valid_H[k][2]},nr={valid_H[k][3]})"
    log(f"  H1_HH[{k},{k}] = {val:.6f} Ha  {label}")
    if val >= 0:
        all_negative = False
log(f"  1sigma' H1 = {H1_HH_raw[0,0]:.4f} Ha (expect ~-0.5)")
pct_from_05 = abs(H1_HH_raw[0, 0] - (-0.5)) / 0.5 * 100
log(f"  Distance from -0.5: {pct_from_05:.1f}% (tolerance 20%)")
log(f"  All negative? {'PASS' if all_negative else 'FAIL'}")

# ==========================================================================
# Check 3: Combined overlap eigenvalues
# ==========================================================================
log("\n--- CHECK 3: Combined S eigenvalues ---")
p0_Li = np.sqrt(2.0 * 7.392)
p0_H = 1.0
log(f"  p0_Li={p0_Li:.4f}, p0_H={p0_H:.4f}, R={R_eq}")

mo_Li = compute_molecular_sturmian_betas(
    Z_A=Z_A, Z_B=Z_B, R=R_eq, p0=p0_Li, nmax=3,
)
mo_H_check = compute_molecular_sturmian_betas(
    Z_A=Z_A, Z_B=Z_B, R=R_eq, p0=p0_H, nmax=2,
)
_, _, _, S_LL, _ = compute_h1_matrix(
    mo_Li, Z_A, Z_B, R_eq, p0_Li, n_grid=N_GRID, orth_method='none_raw',
)
_, _, _, S_HH_check, _ = compute_h1_matrix(
    mo_H_check, Z_A, Z_B, R_eq, p0_H, n_grid=N_GRID, orth_method='none_raw',
)
S_LH_check, _ = compute_cross_set_integrals(
    mo_Li, mo_H_check, Z_A, Z_B, R_eq, p0_Li, p0_H, n_grid=N_GRID,
)

nL = S_LL.shape[0]
nH = S_HH_check.shape[0]
S_comb = np.zeros((nL + nH, nL + nH))
S_comb[:nL, :nL] = S_LL
S_comb[nL:, nL:] = S_HH_check
S_comb[:nL, nL:] = S_LH_check
S_comb[nL:, :nL] = S_LH_check.T

s_evals = np.sort(np.linalg.eigvalsh(S_comb))
log(f"  Combined basis: {nL} Li + {nH} H = {nL+nH} total")
log(f"  S eigenvalues (sorted):")
for i, ev in enumerate(s_evals):
    log(f"    [{i:2d}] {ev:.6f}")
log(f"  Smallest: {s_evals[0]:.6f}")
log(f"  Near-LD (< 0.01)? {'YES — reduce nmax_H' if s_evals[0] < 0.01 else 'NO — OK'}")

# ==========================================================================
# Check 4: H-scale MO H1 diagonal should be negative (from check 2)
# ==========================================================================
log("\n--- CHECK 4: H-scale H1 negativity ---")
log(f"  All H-scale H1 diagonals negative: {'PASS' if all_negative else 'FAIL'}")

# ==========================================================================
# Dual-p0 FCI at R=3.015
# ==========================================================================
log("\n" + "=" * 70)
log("DUAL-P0 FCI AT R=3.015")
log("=" * 70)

t0 = time.perf_counter()
fci = MOSturmianFCI(
    Z_A=Z_A, Z_B=Z_B, R=R_eq, nmax=3, n_electrons=4,
    n_grid=N_GRID, dual_p0=True, nmax_H=2,
)
E_mol, p0_Li_star, p0_H_star = fci.solve_dual(
    p0_Li_init=np.sqrt(2.0 * 7.392),
    p0_H_init=1.0,
    damping=0.3,
    max_iter=15,
    tol=1e-4,
    verbose=True,
)
t_fci = time.perf_counter() - t0

log(f"\n  p0_Li* = {p0_Li_star:.4f}")
log(f"  p0_H*  = {p0_H_star:.4f}")
log(f"  E_mol  = {E_mol:.6f} Ha")
log(f"  E_exact = {E_exact:.4f} Ha")
err_pct = abs(E_mol - E_exact) / abs(E_exact) * 100
log(f"  Error  = {err_pct:.2f}%")
D_e_raw = E_atoms_exact - E_mol
log(f"  D_e_raw = {D_e_raw:.4f} Ha (expt 0.092 Ha)")
log(f"  Bound?  {'YES' if E_mol < E_atoms_exact else 'NO'}")
log(f"  Time   = {t_fci:.1f}s")
log(f"  n_dropped = {fci._n_dropped}")

# ==========================================================================
# PES scan (reduced set for speed)
# ==========================================================================
log("\n" + "=" * 70)
log("PES SCAN")
log("=" * 70)

R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
pes_results = []

for R in R_values:
    t1 = time.perf_counter()
    fci_pes = MOSturmianFCI(
        Z_A=Z_A, Z_B=Z_B, R=R, nmax=3, n_electrons=4,
        n_grid=N_GRID, dual_p0=True, nmax_H=2,
    )
    E, pLi, pH = fci_pes.solve_dual(
        p0_Li_init=np.sqrt(2.0 * 7.392),
        p0_H_init=1.0,
        damping=0.3, max_iter=10, tol=1e-3,
        verbose=False,
    )
    dt = time.perf_counter() - t1
    pes_results.append((R, E, pLi, pH))
    log(f"  R={R:.3f}  E={E:.6f} Ha  p0_Li={pLi:.4f}  p0_H={pH:.4f}  "
        f"t={dt:.0f}s")

# Find R_eq from PES
E_vals = [r[1] for r in pes_results]
R_vals = [r[0] for r in pes_results]
idx_min = np.argmin(E_vals)
log(f"\n  PES minimum at R={R_vals[idx_min]:.3f}, E={E_vals[idx_min]:.6f} Ha")
log(f"  Dissociation: E(R=6.0) = {E_vals[-1]:.6f} Ha, "
    f"E_atoms = {E_atoms_exact:.4f} Ha, "
    f"|delta| = {abs(E_vals[-1] - E_atoms_exact):.4f} Ha")

# ==========================================================================
# Eleven-configuration comparison at R=3.015
# ==========================================================================
log("\n" + "=" * 70)
log("ELEVEN-CONFIGURATION COMPARISON (R=3.015)")
log("=" * 70)
log(f"{'Config':<35s} {'E (Ha)':>10s} {'err%':>8s} {'D_e':>8s} {'R_eq':>6s}")
log("-" * 70)

configs = [
    ("v0.9.8  atom-cent graph", -8.097, 0.33, 0.198, "~2.0"),
    ("v0.9.9  CP-corrected", -8.097, 0.33, 0.083, "~2.5"),
    ("v0.9.11 CP PES", -8.097, 0.33, 0.093, "~2.5"),
    ("v0.9.12 Lowdin ERI", None, None, None, "FAIL"),
    ("v0.9.13 Mulliken K", -8.097, 0.33, 0.110, "~2.5"),
    ("v0.9.18 hybrid D-mat", -8.097, 0.33, 0.143, "<2.0"),
    ("v0.9.31 MO Sturm (single p0)", -10.07, 24.8, None, "<2.0"),
    ("v0.9.32 canon+exactJ", -7.131, 11.6, None, "~3.5"),
    ("v0.9.33 exact K", -6.796, 15.8, None, "6.0"),
]

for name, E, err, De, Req in configs:
    if E is not None:
        log(f"  {name:<35s} {E:>10.4f} {err:>7.1f}% "
            f"{'  ' + f'{De:.3f}' if De is not None else '  unbnd':>8s} "
            f"{Req:>6s}")
    else:
        log(f"  {name:<35s} {'---':>10s} {'---':>8s} {'---':>8s} {Req:>6s}")

# v0.9.34 dual-p0 entry
D_e_934 = E_atoms_exact - E_mol if E_mol < E_atoms_exact else None
R_eq_934 = f"{R_vals[idx_min]:.1f}" if len(R_vals) > 0 else "?"
log(f"  {'v0.9.34 dual-p0':<35s} {E_mol:>10.4f} {err_pct:>7.1f}% "
    f"{'  ' + f'{D_e_934:.3f}' if D_e_934 is not None else '  unbnd':>8s} "
    f"{R_eq_934:>6s}")

log(f"\n  Exact LiH: -8.0706 Ha, D_e=0.092 Ha, R_eq=3.015 Bohr")

# ==========================================================================
# Standard diagnostics
# ==========================================================================
log("\n" + "=" * 70)
log("STANDARD DIAGNOSTICS")
log("=" * 70)

# 1. Dissociation
E_R6 = E_vals[-1] if len(E_vals) > 0 else None
if E_R6 is not None:
    log(f"  Dissociation: E(R=6) = {E_R6:.4f}, E_atoms = {E_atoms_exact:.4f}, "
        f"|delta| = {abs(E_R6 - E_atoms_exact):.4f}")
    log(f"  {'PASS' if abs(E_R6 - E_atoms_exact) < 0.5 else 'FAIL'} "
        f"(tolerance 0.5 Ha)")

# 2. Bound minimum
log(f"  Bound minimum: E_min = {min(E_vals):.4f} < E_atoms = {E_atoms_exact:.4f}? "
    f"{'YES' if min(E_vals) < E_atoms_exact else 'NO'}")

# 3. R_eq in [2.5, 3.5]
log(f"  R_eq = {R_vals[idx_min]:.3f} in [2.5, 3.5]? "
    f"{'YES' if 2.5 <= R_vals[idx_min] <= 3.5 else 'NO'}")


# ==========================================================================
# Save results
# ==========================================================================
with open('debug/data/mo_fci_lih_results.txt', 'w') as f:
    f.write('\n'.join(out_lines))

log(f"\nResults saved to debug/data/mo_fci_lih_results.txt")
log(f"Total wall time: {time.perf_counter() - t0:.0f}s")
