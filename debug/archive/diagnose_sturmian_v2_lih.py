"""
Sturmian v2 LiH Diagnostic Script (v0.9.21)
============================================

v0.9.21 replaces the frozen Fourier cross-nuclear with the p0-dependent
D-matrix formula (Paper 9 Eq. 23), resolving the Fourier-Sturmian
inconsistency that caused v0.9.20 divergence.

Diagnostics:
  1. Self-consistency: bare (damping=1.0), then damped (0.5, 0.3) if needed
  2. PES scan at fixed converged p0
  3. BSSE analysis
  4. Three standard diagnostics (dissociation, bound minimum, BSSE)
"""

import sys
import warnings
import time
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, '.')

from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex,
)


# ---- Configuration ----
NMAX = 3
R_EQ = 3.015
R_SCAN = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
EXACT_LIH = -8.0705
EXACT_DE = 0.0924
HYBRID_DE_CP = 0.143
V020_BSSE = -0.115


def compute_atomic_energies(nmax: int) -> tuple:
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


def run_self_consistency(damping: float = 1.0, max_iter: int = 50) -> tuple:
    """Run Sturmian self-consistency loop with given damping."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R_EQ, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        use_sturmian=True,
    )

    t0 = time.time()
    p0_final, E_mol, n_iter, converged = mol.solve_sturmian_p0(
        R=R_EQ, tol=1e-6, max_iter=max_iter, damping=damping,
    )
    elapsed = time.time() - t0

    return mol, p0_final, E_mol, n_iter, converged, elapsed


def run_pes_scan(p0_fixed: float) -> list:
    """PES scan with fixed p0 from self-consistency at R_eq."""
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
    return results


def run_bsse_sturmian() -> dict:
    """BSSE with standard (non-Sturmian) ghost basis at R_eq."""
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

    log("Sturmian v2 LiH Diagnostic (v0.9.21)")
    log(f"nmax={NMAX}, R_eq={R_EQ} bohr")
    log(f"Fix: D-matrix cross-nuclear replaces frozen Fourier")
    log("")

    # 1. Atomic reference energies
    E_sep, E_li, E_h = compute_atomic_energies(NMAX)
    log(f"Atomic reference: E(Li)={E_li:.6f}, E(H)={E_h:.6f}, "
        f"E_sep={E_sep:.6f} Ha")

    # 2. Self-consistency: try bare first, then damped
    log("\n" + "=" * 70)
    log("SELF-CONSISTENCY LOOP")
    log("=" * 70)

    best_result = None

    for damping in [1.0, 0.5, 0.3]:
        log(f"\n--- damping={damping:.1f} ---")
        try:
            mol, p0, E_mol, n_iter, converged, elapsed = \
                run_self_consistency(damping=damping, max_iter=50)
            log(f"  converged:    {converged}")
            log(f"  p0_final:     {p0:.6f}")
            log(f"  E_mol:        {E_mol:.6f} Ha")
            log(f"  n_iterations: {n_iter}")
            log(f"  wall time:    {elapsed:.1f}s")

            if best_result is None or converged:
                best_result = (mol, p0, E_mol, n_iter, converged, damping)

            if converged:
                log(f"  SUCCESS: converged at damping={damping:.1f}")
                break
        except Exception as e:
            log(f"  ERROR: {e}")

    mol, p0_final, E_mol, n_iter, converged, damping_used = best_result
    log(f"\nBest result: damping={damping_used:.1f}, p0={p0_final:.6f}, "
        f"E_mol={E_mol:.6f}, converged={converged}")

    # 3. PES scan with fixed p0
    log("\n" + "=" * 70)
    log(f"PES SCAN (fixed p0={p0_final:.6f})")
    log("=" * 70)

    pes = run_pes_scan(p0_final)
    log(f"\n  {'R':>6}  {'E (Ha)':>12}")
    log(f"  {'---':>6}  {'---':>12}")
    for R, E in pes:
        log(f"  {R:6.3f}  {E:12.6f}")

    # 4. BSSE
    log("\n" + "=" * 70)
    log("BSSE (standard ghost basis)")
    log("=" * 70)

    bsse = run_bsse_sturmian()
    log(f"  E(Li, own):   {bsse['E_li']:.6f} Ha")
    log(f"  E(H, own):    {bsse['E_h']:.6f} Ha")
    log(f"  E(Li, ghost): {bsse['E_A_ghost']:.6f} Ha")
    log(f"  E(H, ghost):  {bsse['E_B_ghost']:.6f} Ha")
    log(f"  BSSE_Li:      {bsse['BSSE_A']:.6f} Ha")
    log(f"  BSSE_H:       {bsse['BSSE_B']:.6f} Ha")
    log(f"  BSSE_total:   {bsse['BSSE']:.6f} Ha")
    log(f"  v0.9.20 BSSE: {V020_BSSE:.6f} Ha")
    log(f"  Changed:      {abs(bsse['BSSE'] - V020_BSSE) > 0.001}")

    # 5. Diagnostics
    log("\n" + "=" * 70)
    log("DIAGNOSTICS")
    log("=" * 70)

    R_vals = [r for r, _ in pes]
    E_vals = [e for _, e in pes]
    E_min = min(E_vals)
    R_min = R_vals[E_vals.index(E_min)]

    E_at_6 = dict(pes).get(6.0, E_sep)
    D_e_raw = E_at_6 - E_min

    E_ghost_sum = bsse['E_A_ghost'] + bsse['E_B_ghost']
    D_e_cp = E_ghost_sum - E_min

    log(f"\n1. Dissociation limit (R=6.0):")
    log(f"   E(R=6.0) = {E_at_6:.6f} Ha")
    log(f"   D_e_raw(R=6.0 - R_min) = {D_e_raw:.6f} Ha")
    log(f"   PASS: {D_e_raw > 0.001}")

    log(f"\n2. Bound minimum:")
    log(f"   R_min = {R_min:.3f} bohr (expt: {R_EQ})")
    log(f"   E_min = {E_min:.6f} Ha (exact: {EXACT_LIH})")
    log(f"   D_e_raw = {D_e_raw:.6f} Ha")
    log(f"   D_e_CP  = {D_e_cp:.6f} Ha (expt: {EXACT_DE})")

    log(f"\n3. BSSE:")
    log(f"   BSSE = {bsse['BSSE']:.6f} Ha (v0.9.20: {V020_BSSE})")
    log(f"   Changed from v0.9.20: {abs(bsse['BSSE'] - V020_BSSE) > 0.001}")

    log(f"\n4. Comparison with baselines:")
    log(f"   D_e_CP(Sturmian v2): {D_e_cp:.4f} Ha")
    log(f"   D_e_CP(hybrid):      {HYBRID_DE_CP:.4f} Ha")
    log(f"   D_e_CP(expt):        {EXACT_DE:.4f} Ha")

    log(f"\n5. Self-consistent p0:")
    log(f"   p0_init = {np.sqrt(5.0):.6f}")
    log(f"   p0_final = {p0_final:.6f}")
    log(f"   damping = {damping_used:.1f}")
    log(f"   n_iterations = {n_iter}")
    log(f"   converged = {converged}")

    # Save results
    outpath = "debug/data/sturmian_v2_lih_results.txt"
    with open(outpath, 'w') as f:
        f.write('\n'.join(output_lines))
    print(f"\nResults saved to {outpath}")


if __name__ == '__main__':
    main()
