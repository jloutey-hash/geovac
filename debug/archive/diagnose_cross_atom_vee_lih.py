"""Diagnostics for v0.9.35: Cross-Atom V_ee Extension to All l.

Compares all-l vs s-only vs disabled cross-atom V_ee:
  1. J_cross table for all (n_A, l_A, n_B, l_B) pairs at R=3.015
  2. dE(R) correction: E(all_l) - E(s_only) at several R values
  3. CP-corrected PES with all-l V_ee
  4. Term decomposition comparison at R=2.5 and R=3.015
  5. Twelve-configuration comparison table

Output: debug/data/cross_atom_vee_lih_v0935.txt
"""

from __future__ import annotations

import sys
import time
import numpy as np

sys.path.insert(0, '.')
from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
    compute_cross_atom_J,
)

Z_A, Z_B = 3, 1
R_eq = 3.015
E_exact = -8.0706  # exact LiH at R=3.015
E_atoms_exact = -7.892  # sum of exact Li + H atomic energies
NMAX = 3

out_lines: list[str] = []


def log(msg: str = '') -> None:
    print(msg)
    out_lines.append(msg)


def build_mol(R: float, cross_atom_vee: object = True,
              fci_method: str = 'matrix') -> MolecularLatticeIndex:
    return MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=4,
        n_bridges=10, vee_method='slater_full',
        fci_method=fci_method,
        cross_atom_vee=cross_atom_vee,
    )


t0_global = time.perf_counter()

# ==========================================================================
# Section 1: J_cross table for all (n, l) pairs at R=3.015
# ==========================================================================
log("=" * 70)
log("v0.9.35 Cross-Atom V_ee Extension -- LiH Diagnostics")
log(f"nmax={NMAX}")
log("=" * 70)

log("\n--- SECTION 1: J_cross Table at R=3.015 ---")
log(f"{'(nA,lA)':>8s}  {'(nB,lB)':>8s}  {'J_cross (Ha)':>12s}")
log("-" * 35)

# All unique (n,l) pairs up to nmax=3
nl_pairs = []
for n in range(1, NMAX + 1):
    for l in range(n):
        nl_pairs.append((n, l))

for na, la in nl_pairs:
    for nb, lb in nl_pairs:
        j = compute_cross_atom_J(na, la, nb, lb, R=R_eq, ZA=float(Z_A), ZB=float(Z_B))
        log(f"  ({na},{la})     ({nb},{lb})     {j:12.6f}")

# ==========================================================================
# Section 2: dE(R) correction -- all_l vs s_only vs off
# ==========================================================================
log("\n" + "=" * 70)
log("SECTION 2: dE(R) = E(all_l) - E(s_only) and E(s_only) - E(off)")
log("=" * 70)

R_compare = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
log(f"\n{'R':>6s}  {'E(all_l)':>10s}  {'E(s_only)':>10s}  {'E(off)':>10s}  "
    f"{'dE(all-s)':>10s}  {'dE(s-off)':>10s}")
log("-" * 70)

pes_all: list[tuple[float, float]] = []
pes_s_only: list[tuple[float, float]] = []

for R in R_compare:
    t1 = time.perf_counter()
    mol_all = build_mol(R, cross_atom_vee=True)
    E_all = mol_all.compute_ground_state(n_states=1)[0][0]

    mol_s = build_mol(R, cross_atom_vee='s_only')
    E_s = mol_s.compute_ground_state(n_states=1)[0][0]

    mol_off = build_mol(R, cross_atom_vee=False)
    E_off = mol_off.compute_ground_state(n_states=1)[0][0]

    dt = time.perf_counter() - t1
    delta_all_s = E_all - E_s
    delta_s_off = E_s - E_off

    pes_all.append((R, E_all))
    pes_s_only.append((R, E_s))

    log(f"  {R:5.3f}  {E_all:10.6f}  {E_s:10.6f}  {E_off:10.6f}  "
        f"{delta_all_s:10.6f}  {delta_s_off:10.6f}  ({dt:.0f}s)")

# ==========================================================================
# Section 3: CP-corrected PES with all-l V_ee
# ==========================================================================
log("\n" + "=" * 70)
log("SECTION 3: CP-Corrected PES (all-l V_ee)")
log("=" * 70)

# First get atomic energies at same nmax
li = LatticeIndex(n_electrons=3, max_n=NMAX, nuclear_charge=3,
                  vee_method='slater_full', h1_method='exact')
E_li = li.compute_ground_state(n_states=1)[0][0]

h = LatticeIndex(n_electrons=1, max_n=NMAX, nuclear_charge=1,
                 vee_method='slater_full', h1_method='exact')
E_h = h.compute_ground_state(n_states=1)[0][0]
E_atoms_nmax = E_li + E_h

log(f"\n  E(Li, nmax={NMAX}) = {E_li:.6f} Ha")
log(f"  E(H, nmax={NMAX})  = {E_h:.6f} Ha")
log(f"  E_atoms(nmax)      = {E_atoms_nmax:.6f} Ha")

log(f"\n{'R':>6s}  {'E_mol':>10s}  {'BSSE':>10s}  {'D_e_raw':>10s}  "
    f"{'D_e_CP':>10s}  {'err%':>6s}")
log("-" * 60)

cp_results: list[dict] = []

for R in R_compare:
    t1 = time.perf_counter()
    mol = build_mol(R, cross_atom_vee=True)
    E_mol = mol.compute_ground_state(n_states=1)[0][0]

    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX, R=R,
        n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='matrix',
    )
    E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
    D_e_raw = E_atoms_nmax - E_mol
    D_e_CP = E_ghost_sum - E_mol

    err_pct = abs(E_mol - E_exact) / abs(E_exact) * 100 if R == R_eq else 0.0
    dt = time.perf_counter() - t1

    cp_results.append({
        'R': R, 'E_mol': E_mol, 'BSSE': bsse['BSSE'],
        'D_e_raw': D_e_raw, 'D_e_CP': D_e_CP,
    })

    log(f"  {R:5.3f}  {E_mol:10.6f}  {bsse['BSSE']:10.6f}  "
        f"{D_e_raw:10.6f}  {D_e_CP:10.6f}  "
        f"{'  ' + f'{err_pct:.1f}%' if R == R_eq else '':>6s}  ({dt:.0f}s)")

# Find R_eq from PES
E_cp_vals = [r['D_e_CP'] for r in cp_results]
R_vals = [r['R'] for r in cp_results]
idx_max = np.argmax(E_cp_vals)
log(f"\n  PES minimum (max D_e_CP) at R={R_vals[idx_max]:.3f}, "
    f"D_e_CP={E_cp_vals[idx_max]:.6f} Ha")
log(f"  Experimental: R_eq=3.015, D_e=0.092 Ha")

# Also find minimum energy point
E_mol_vals = [r['E_mol'] for r in cp_results]
idx_min_E = np.argmin(E_mol_vals)
log(f"  Raw energy minimum at R={R_vals[idx_min_E]:.3f}, "
    f"E={E_mol_vals[idx_min_E]:.6f} Ha")

# ==========================================================================
# Section 4: Term decomposition at R=2.5 and R=3.015
# ==========================================================================
log("\n" + "=" * 70)
log("SECTION 4: Term Decomposition (all-l vs s-only)")
log("=" * 70)

for R_decomp in [2.5, 3.015]:
    log(f"\n  --- R = {R_decomp:.3f} ---")

    for mode, label in [(True, 'all_l'), ('s_only', 's_only')]:
        mol = build_mol(R_decomp, cross_atom_vee=mode)
        eigvals, eigvecs = mol.compute_ground_state(n_states=1)
        E_total = eigvals[0]
        civec = eigvecs[:, 0]
        decomp = mol.decompose_energy(civec, E_total)

        log(f"  [{label:>6s}] E_total={decomp['E_total']:.6f}")
        log(f"          T={decomp['T']:.6f}  V_nA={decomp['V_nA']:.6f}  "
            f"V_nB={decomp['V_nB']:.6f}")
        log(f"          V_crossA={decomp['V_cross_A']:.6f}  "
            f"V_crossB={decomp['V_cross_B']:.6f}  "
            f"V_bridge={decomp['V_bridge']:.6f}")
        log(f"          V_ee={decomp['V_ee']:.6f}  V_NN={decomp['V_NN']:.6f}")

# ==========================================================================
# Section 5: Standard diagnostics
# ==========================================================================
log("\n" + "=" * 70)
log("SECTION 5: Standard Diagnostics")
log("=" * 70)

# D_e_CP at R=3.015
cp_eq = next(r for r in cp_results if r['R'] == R_eq)
log(f"\n  D_e_CP(R=3.015) = {cp_eq['D_e_CP']:.4f} Ha "
    f"(target [0.08, 0.12], expt 0.092)")
in_target = 0.08 <= cp_eq['D_e_CP'] <= 0.12
log(f"  In target range: {'YES' if in_target else 'NO'}")

# R_eq
R_eq_found = R_vals[idx_max]
log(f"  R_eq = {R_eq_found:.3f} (target [2.7, 3.3], expt 3.015)")
r_in_target = 2.7 <= R_eq_found <= 3.3
log(f"  In target range: {'YES' if r_in_target else 'NO'}")

# Dissociation
cp_r6 = next(r for r in cp_results if r['R'] == 6.0)
log(f"  |D_e_CP(R=6)| = {abs(cp_r6['D_e_CP']):.4f} Ha (target <= 0.005)")
diss_ok = abs(cp_r6['D_e_CP']) <= 0.005
log(f"  Dissociation OK: {'YES' if diss_ok else 'NO'}")

# ==========================================================================
# Section 6: Twelve-configuration comparison
# ==========================================================================
log("\n" + "=" * 70)
log("TWELVE-CONFIGURATION COMPARISON (R=3.015)")
log("=" * 70)
log(f"{'Config':<38s} {'E (Ha)':>10s} {'err%':>8s} {'D_e':>8s} {'R_eq':>6s}")
log("-" * 74)

configs = [
    ("v0.9.8  atom-cent graph", -8.097, 0.33, 0.198, "~2.0"),
    ("v0.9.9  CP-corrected", -8.097, 0.33, 0.083, "~2.5"),
    ("v0.9.11 CP PES (s-only)", -8.097, 0.33, 0.093, "~2.5"),
    ("v0.9.12 Lowdin ERI", None, None, None, "FAIL"),
    ("v0.9.13 Mulliken K", -8.097, 0.33, 0.110, "~2.5"),
    ("v0.9.18 hybrid D-mat", -8.097, 0.33, 0.143, "<2.0"),
    ("v0.9.31 MO Sturm (single p0)", -10.07, 24.8, None, "<2.0"),
    ("v0.9.32 canon+exactJ", -7.131, 11.6, None, "~3.5"),
    ("v0.9.33 exact K", -6.796, 15.8, None, "6.0"),
    ("v0.9.34 dual-p0", -10.1, 25.0, None, "n/a"),
]

for name, E, err, De, Req in configs:
    if E is not None:
        log(f"  {name:<38s} {E:>10.4f} {err:>7.1f}% "
            f"{'  ' + f'{De:.3f}' if De is not None else '  unbnd':>8s} "
            f"{Req:>6s}")
    else:
        log(f"  {name:<38s} {'---':>10s} {'---':>8s} {'---':>8s} {Req:>6s}")

# v0.9.35 entry
E_935 = cp_eq['E_mol']
err_935 = abs(E_935 - E_exact) / abs(E_exact) * 100
D_e_935 = cp_eq['D_e_CP']
R_eq_935 = f"{R_eq_found:.1f}"
log(f"  {'v0.9.35 all-l cross-atom V_ee':<38s} {E_935:>10.4f} {err_935:>7.1f}% "
    f"{'  ' + f'{D_e_935:.3f}':>8s} {R_eq_935:>6s}")

log(f"\n  Exact LiH: -8.0706 Ha, D_e=0.092 Ha, R_eq=3.015 Bohr")

# ==========================================================================
# Save
# ==========================================================================
t_total = time.perf_counter() - t0_global
log(f"\nTotal wall time: {t_total:.0f}s")

with open('debug/data/cross_atom_vee_lih_v0935.txt', 'w') as f:
    f.write('\n'.join(out_lines))

log(f"Results saved to debug/data/cross_atom_vee_lih_v0935.txt")
