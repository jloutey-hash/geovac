"""
Sturmian Basis LiH Diagnostic Script (v0.9.20)
===============================================

Implements Paper 9 Sturmian basis for LiH and runs three diagnostics:
  1. Dissociation limit: D_e_CP ≤ 0.001 Ha at R=6.0?
  2. Bound minimum: PES has minimum with D_e_CP > 0?
  3. BSSE: different from Mulliken baseline 0.115 Ha?

Also reports: self-consistent p0, iteration count, D_e_CP vs baselines.

PES scan uses p0 converged at R=3.015 (fixed for all R points).
"""

import sys
import warnings
import time
import numpy as np

# Suppress all construction noise
warnings.filterwarnings("ignore")

sys.path.insert(0, '.')

from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
)


# ---- Configuration ----
NMAX = 3
R_EQ = 3.015
R_SCAN = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
EXACT_LIH = -8.0705    # Ha (exact NR)
EXACT_DE = 0.0924      # Ha (experimental binding)
HYBRID_DE_CP = 0.143    # Ha (v0.9.18 baseline)


def compute_atomic_energies(nmax: int) -> float:
    """E(Li) + E(H) at given nmax."""
    li = LatticeIndex(
        n_electrons=3, max_n=nmax, nuclear_charge=3,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=nmax, nuclear_charge=1,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    return E_li + E_h, E_li, E_h


def run_sturmian_self_consistency():
    """Run Sturmian self-consistency loop for LiH at R_eq."""
    print("=" * 70)
    print("STURMIAN SELF-CONSISTENCY LOOP")
    print("=" * 70)

    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R_EQ, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        use_sturmian=True,
    )

    t0 = time.time()
    p0_final, E_mol, n_iter, converged = mol.solve_sturmian_p0(
        R=R_EQ, tol=1e-6, max_iter=50,
    )
    elapsed = time.time() - t0

    print(f"\n--- Self-consistency result ---")
    print(f"  converged:    {converged}")
    print(f"  p0_final:     {p0_final:.6f}")
    print(f"  E_mol:        {E_mol:.6f} Ha")
    print(f"  n_iterations: {n_iter}")
    print(f"  wall time:    {elapsed:.1f}s")

    return mol, p0_final, E_mol, n_iter, converged


def run_pes_scan(p0_fixed: float):
    """PES scan with fixed p0 from self-consistency at R_eq."""
    print("\n" + "=" * 70)
    print(f"PES SCAN (fixed p0={p0_fixed:.6f} from R={R_EQ})")
    print("=" * 70)

    results = []
    for R in R_SCAN:
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
            R=R, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True, sturmian_p0=p0_fixed,
        )
        eigvals, _ = mol.compute_ground_state(n_states=1)
        E = eigvals[0]
        results.append((R, E))
        print(f"  R={R:.3f}  E={E:.6f} Ha")

    return results


def run_bsse_sturmian(p0_fixed: float):
    """BSSE with Sturmian basis at R_eq."""
    print("\n" + "=" * 70)
    print("BSSE (STURMIAN)")
    print("=" * 70)

    # Own-basis energies (same as non-Sturmian — single-atom has no Sturmian mode)
    li = LatticeIndex(
        n_electrons=3, max_n=NMAX, nuclear_charge=3,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=NMAX, nuclear_charge=1,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    # Ghost-basis energies: use standard (non-Sturmian) ghost since
    # Sturmian mode doesn't apply to ghost atoms
    mol_A_ghost = MolecularLatticeIndex(
        Z_A=3, Z_B=0, nmax_A=NMAX, nmax_B=NMAX,
        R=R_EQ, n_electrons=3,
        vee_method='slater_full', fci_method='auto',
    )
    E_A_ghost = mol_A_ghost.compute_ground_state(n_states=1)[0][0]

    mol_B_ghost = MolecularLatticeIndex(
        Z_A=0, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R_EQ, n_electrons=1,
        vee_method='slater_full', fci_method='auto',
    )
    E_B_ghost = mol_B_ghost.compute_ground_state(n_states=1)[0][0]

    BSSE_A = E_A_ghost - E_li
    BSSE_B = E_B_ghost - E_h
    BSSE = BSSE_A + BSSE_B

    print(f"  E(Li, own):   {E_li:.6f} Ha")
    print(f"  E(H, own):    {E_h:.6f} Ha")
    print(f"  E(Li, ghost): {E_A_ghost:.6f} Ha")
    print(f"  E(H, ghost):  {E_B_ghost:.6f} Ha")
    print(f"  BSSE_Li:      {BSSE_A:.6f} Ha")
    print(f"  BSSE_H:       {BSSE_B:.6f} Ha")
    print(f"  BSSE_total:   {BSSE:.6f} Ha")

    return {
        'E_li': E_li, 'E_h': E_h,
        'E_A_ghost': E_A_ghost, 'E_B_ghost': E_B_ghost,
        'BSSE_A': BSSE_A, 'BSSE_B': BSSE_B, 'BSSE': BSSE,
    }


def main():
    output_lines = []
    def log(msg: str = ""):
        print(msg)
        output_lines.append(msg)

    log("Sturmian Basis LiH Diagnostic (v0.9.20)")
    log(f"nmax={NMAX}, R_eq={R_EQ} bohr")
    log("")

    # 1. Atomic reference energies
    E_sep, E_li, E_h = compute_atomic_energies(NMAX)
    log(f"Atomic reference: E(Li)={E_li:.6f}, E(H)={E_h:.6f}, "
        f"E_sep={E_sep:.6f} Ha")

    # 2. Self-consistency loop
    mol, p0_final, E_mol, n_iter, converged = run_sturmian_self_consistency()

    log(f"\nSelf-consistency: converged={converged}, p0={p0_final:.6f}, "
        f"E_mol={E_mol:.6f} Ha, n_iter={n_iter}")

    # 3. PES scan with fixed p0
    pes = run_pes_scan(p0_final)

    log("\nPES (fixed p0):")
    log(f"  {'R':>6}  {'E (Ha)':>12}")
    log(f"  {'---':>6}  {'---':>12}")
    for R, E in pes:
        log(f"  {R:6.3f}  {E:12.6f}")

    # 4. BSSE
    bsse = run_bsse_sturmian(p0_final)

    # 5. Diagnostics
    log("\n" + "=" * 70)
    log("DIAGNOSTICS")
    log("=" * 70)

    # Find minimum of PES
    R_vals = [r for r, _ in pes]
    E_vals = [e for _, e in pes]
    E_min = min(E_vals)
    R_min = R_vals[E_vals.index(E_min)]

    # D_e raw (no CP)
    E_at_6 = dict(pes).get(6.0, E_sep)
    D_e_raw = E_at_6 - E_min

    # D_e CP-corrected
    E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
    D_e_cp = E_ghost_sum - E_min

    log(f"\n1. Dissociation limit (R=6.0):")
    log(f"   E(R=6.0) = {E_at_6:.6f} Ha")
    log(f"   D_e_raw(R=6.0 - R_min) = {D_e_raw:.6f} Ha")

    log(f"\n2. Bound minimum:")
    log(f"   R_min = {R_min:.3f} bohr (expt: {R_EQ})")
    log(f"   E_min = {E_min:.6f} Ha (exact: {EXACT_LIH})")
    log(f"   D_e_raw = {D_e_raw:.6f} Ha")
    log(f"   D_e_CP  = {D_e_cp:.6f} Ha (expt: {EXACT_DE})")

    log(f"\n3. BSSE:")
    log(f"   BSSE = {bsse['BSSE']:.6f} Ha (Mulliken baseline: -0.115)")

    log(f"\n4. Comparison with baselines:")
    log(f"   D_e_CP(Sturmian): {D_e_cp:.4f} Ha")
    log(f"   D_e_CP(hybrid):   {HYBRID_DE_CP:.4f} Ha")
    log(f"   D_e_CP(expt):     {EXACT_DE:.4f} Ha")

    log(f"\n5. Self-consistent p0:")
    log(f"   p0_init = {np.sqrt(5.0):.6f}")
    log(f"   p0_final = {p0_final:.6f}")
    log(f"   n_iterations = {n_iter}")
    log(f"   converged = {converged}")

    # Convergence/divergence diagnosis
    if not converged:
        if E_mol >= 0:
            log(f"\n   DIAGNOSIS: E_mol became positive (unbound state).")
            log(f"   The Sturmian basis at p0={p0_final:.4f} is too compact/diffuse.")
        elif p0_final > 10:
            log(f"\n   DIAGNOSIS: p0 growing without bound (E_mol too negative).")
            log(f"   Diffuse Sturmian basis over-stabilizes the molecule.")
        else:
            log(f"\n   DIAGNOSIS: Oscillating or slow convergence.")
            log(f"   May need damping or larger nmax.")

    # Save results
    outpath = "debug/data/sturmian_lih_results.txt"
    with open(outpath, 'w') as f:
        f.write('\n'.join(output_lines))
    print(f"\nResults saved to {outpath}")


if __name__ == '__main__':
    main()
