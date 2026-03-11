"""
Diagnostic: SO(4) CG Cross-Center V_ee for LiH Sturmian (v0.9.27)

Tests whether SO(4) Clebsch-Gordan two-electron integrals can rescue the
single-p0 Sturmian from being unbound.

Expected: NO — the H 1s one-electron deficit (+1.85 Ha at p0*=3.168) is
~15x the physical LiH correlation energy (~0.12 Ha).

This script:
  1. Prints cross-center ERI table (D-matrix-weighted J values)
  2. Runs 3-point PES at R = 2.0, 3.015, 6.0 with use_so4_vee=True
  3. Gap analysis: E_mol vs E_atoms, cross-center V_ee vs H 1s deficit
  4. Conclusion statement on molecular Sturmian necessity
"""

import warnings
import sys
import time
import numpy as np

sys.path.insert(0, '.')
from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex,
)
from geovac.wigner_so4 import bond_angle, d_matrix_block


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NMAX = 3
Z_A, Z_B = 3, 1
N_ELECTRONS = 4
R_VALUES = [2.0, 3.015, 6.0]
R_EQ_EXPT = 3.015

outpath = 'debug/data/so4_vee_lih_results.txt'


def compute_atomic_energies() -> dict:
    """Compute isolated Li and H energies at nmax=3."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = LatticeIndex(
            n_electrons=3, max_n=NMAX, nuclear_charge=Z_A,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        )
        E_Li = li.compute_ground_state(n_states=1)[0][0]

        h = LatticeIndex(
            n_electrons=1, max_n=NMAX, nuclear_charge=Z_B,
            vee_method='slater_full', h1_method='exact', fci_method='matrix',
        )
        E_H = h.compute_ground_state(n_states=1)[0][0]

    return {'E_Li': E_Li, 'E_H': E_H, 'E_atoms': E_Li + E_H}


def main():
    t_start = time.perf_counter()
    lines = []

    def log(s: str = ''):
        print(s)
        lines.append(s)

    log("=" * 80)
    log("LiH SO(4) CG Cross-Center V_ee Diagnostic (v0.9.27)")
    log("nmax=%d, Z_A=%d, Z_B=%d, N_e=%d" % (NMAX, Z_A, Z_B, N_ELECTRONS))
    log("=" * 80)

    # Atomic energies
    atoms = compute_atomic_energies()
    log("\nAtomic energies (nmax=%d):" % NMAX)
    log("  E(Li) = %.6f Ha" % atoms['E_Li'])
    log("  E(H)  = %.6f Ha" % atoms['E_H'])
    log("  E(Li+H) = %.6f Ha" % atoms['E_atoms'])

    # -----------------------------------------------------------------------
    # Cross-center ERI table at R=3.015 with converged p0
    # -----------------------------------------------------------------------
    log("\n" + "=" * 80)
    log("Cross-center ERI table at R=%.3f" % R_EQ_EXPT)
    log("J(A_{n_a l_a m_a}, B_{n_b l_b m_b}) = sum D^(n_b) * F0_Sturmian")
    log("=" * 80)

    # First run self-consistency without SO4 V_ee to get converged p0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_ref = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
            R=R_EQ_EXPT, n_electrons=N_ELECTRONS,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True,
        )
        p0_ref, E_ref, _, _ = mol_ref.solve_sturmian_p0(
            R=R_EQ_EXPT, damping=0.5, max_iter=50,
        )
    log("\nReference (no cross-center V_ee): p0*=%.4f, E=%.6f Ha" % (p0_ref, E_ref))

    # Build SO4 cross-center V_ee at converged p0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_eri = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
            R=R_EQ_EXPT, n_electrons=N_ELECTRONS,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True, sturmian_p0=p0_ref,
            use_so4_vee=True,
        )
        mol_eri._build_so4_cross_vee(p0_ref)

    gamma = bond_angle(R_EQ_EXPT, p0_ref)
    log("\ngamma(R=%.3f, p0=%.4f) = %.6f rad (%.2f deg)" % (
        R_EQ_EXPT, p0_ref, gamma, np.degrees(gamma)))

    # Print ERI table
    nA = mol_eri._n_spatial_A
    states_A = mol_eri._li_A.lattice.states
    states_B = mol_eri._li_B.lattice.states

    log("\n  %-12s  %-12s  %12s  %12s" % (
        'A (n,l,m)', 'B (n,l,m)', 'J_cross', 'D_elem'))
    log("  " + "-" * 55)

    for i_a, (n_a, l_a, m_a) in enumerate(states_A):
        for i_b, (n_b, l_b, m_b) in enumerate(states_B):
            j_b = i_b + nA
            j_val = mol_eri._eri.get((i_a, j_b, i_a, j_b), 0.0)
            if abs(j_val) > 1e-10:
                # Get dominant D-matrix element
                D_n = d_matrix_block(n_b, gamma)
                def lm_idx(l, m):
                    idx = 0
                    for ll in range(l):
                        idx += 2 * ll + 1
                    return idx + m + l
                col = lm_idx(l_b, m_b)
                d_elem = D_n[col, col]  # diagonal D-matrix element
                log("  (%d,%d,%+d)      (%d,%d,%+d)      %12.6f  %12.6f" % (
                    n_a, l_a, m_a, n_b, l_b, m_b, j_val, d_elem))

    # F0_Sturmian reference values
    log("\nF0_Sturmian reference (p0=%.4f):" % p0_ref)
    log("  F0(1s,1s) = p0 * 5/8 = %.6f Ha" % (p0_ref * 5.0 / 8.0))
    log("  F0(1s,2s) = p0 * 17/81 = %.6f Ha" % (p0_ref * 17.0 / 81.0))

    # Same-center reference
    j_same_1s = mol_eri._eri.get((0, 0, 0, 0), 0.0)
    j_cross_1s = mol_eri._eri.get((0, nA, 0, nA), 0.0)
    log("\n  J_same(1s_A, 1s_A) = %.6f Ha" % j_same_1s)
    log("  J_cross(1s_A, 1s_B) = %.6f Ha" % j_cross_1s)
    log("  Ratio cross/same = %.4f" % (j_cross_1s / j_same_1s if j_same_1s else 0))

    # -----------------------------------------------------------------------
    # 3-point PES with SO(4) CG V_ee
    # -----------------------------------------------------------------------
    log("\n" + "=" * 80)
    log("3-point PES with SO(4) CG cross-center V_ee")
    log("use_sturmian=True, use_so4_vee=True, damping=0.5")
    log("=" * 80)

    pes_results = {}
    for R in R_VALUES:
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
                R=R, n_electrons=N_ELECTRONS,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True,
                use_so4_vee=True,
            )
            p0_f, E_mol, n_iter, converged = mol.solve_sturmian_p0(
                R=R, damping=0.5, max_iter=50,
            )
        dt = time.perf_counter() - t0
        pes_results[R] = {
            'p0': p0_f, 'E_mol': E_mol, 'n_iter': n_iter,
            'converged': converged, 'time': dt,
        }
        log("  R=%.3f: p0=%.4f, E=%.6f Ha, %d iter, %s (%.1fs)" % (
            R, p0_f, E_mol, n_iter,
            'CONVERGED' if converged else 'NOT CONVERGED', dt))

    # PES summary table
    log("\n  %8s  %8s  %12s  %10s  %10s  %8s" % (
        'R', 'p0', 'E_mol', 'D_e_raw', 'E-E_atoms', 'bound?'))
    log("  %s  %s  %s  %s  %s  %s" % (
        '-'*8, '-'*8, '-'*12, '-'*10, '-'*10, '-'*8))
    for R in R_VALUES:
        res = pes_results[R]
        D_e_raw = atoms['E_atoms'] - res['E_mol']
        gap = res['E_mol'] - atoms['E_atoms']
        bound = 'YES' if res['E_mol'] < atoms['E_atoms'] else 'NO'
        log("  %8.3f  %8.4f  %12.6f  %10.6f  %10.6f  %8s" % (
            R, res['p0'], res['E_mol'], D_e_raw, gap, bound))

    # -----------------------------------------------------------------------
    # Reference: v0.9.21 (no cross-center V_ee)
    # -----------------------------------------------------------------------
    log("\n" + "=" * 80)
    log("Reference: v0.9.21 Sturmian (no cross-center V_ee)")
    log("=" * 80)

    ref_results = {}
    for R in R_VALUES:
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
                R=R, n_electrons=N_ELECTRONS,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True,
                use_so4_vee=False,
            )
            p0_f, E_mol, n_iter, converged = mol.solve_sturmian_p0(
                R=R, damping=0.5, max_iter=50,
            )
        dt = time.perf_counter() - t0
        ref_results[R] = {
            'p0': p0_f, 'E_mol': E_mol, 'n_iter': n_iter,
            'converged': converged, 'time': dt,
        }
        log("  R=%.3f: p0=%.4f, E=%.6f Ha, %d iter, %s (%.1fs)" % (
            R, p0_f, E_mol, n_iter,
            'CONVERGED' if converged else 'NOT CONVERGED', dt))

    log("\n  %8s  %12s  %12s  %12s" % (
        'R', 'E(no V_ee)', 'E(SO4 V_ee)', 'Delta'))
    log("  %s  %s  %s  %s" % ('-'*8, '-'*12, '-'*12, '-'*12))
    for R in R_VALUES:
        E_ref = ref_results[R]['E_mol']
        E_so4 = pes_results[R]['E_mol']
        delta = E_so4 - E_ref
        log("  %8.3f  %12.6f  %12.6f  %12.6f" % (R, E_ref, E_so4, delta))

    # -----------------------------------------------------------------------
    # Gap analysis
    # -----------------------------------------------------------------------
    log("\n" + "=" * 80)
    log("Gap Analysis")
    log("=" * 80)

    # H 1s one-electron deficit at p0*
    for label, results in [('SO4 V_ee', pes_results), ('no V_ee', ref_results)]:
        R_eq = R_EQ_EXPT
        if R_eq in results:
            p0_star = results[R_eq]['p0']
            h1s_deficit = p0_star**2 / 2.0 - 1.0 * p0_star / 1.0 - (-0.5)
            E_mol = results[R_eq]['E_mol']
            gap = E_mol - atoms['E_atoms']
            log("\n  %s at R=%.3f:" % (label, R_eq))
            log("    p0* = %.4f" % p0_star)
            log("    H 1s deficit = eps_Sturmian(1s,p0*) - eps_exact(1s)")
            log("                 = (%.4f^2/2 - %.4f) - (-0.5)" % (p0_star, p0_star))
            log("                 = %.4f Ha" % h1s_deficit)
            log("    E_mol - E_atoms = %.4f Ha (>0 means UNBOUND)" % gap)
            log("    Cross-center V_ee at 1s = %.6f Ha" % (
                j_cross_1s if label == 'SO4 V_ee' else 0.0))

    # Physical correlation energy for reference
    D_e_phys = 0.0924  # Ha
    log("\n  Physical LiH D_e = %.4f Ha" % D_e_phys)
    log("  Physical correlation energy ~ 0.12 Ha")

    p0_star = pes_results[R_EQ_EXPT]['p0']
    h1s_deficit = p0_star**2 / 2.0 - p0_star + 0.5
    log("\n  H 1s deficit / D_e = %.1f x" % (h1s_deficit / D_e_phys))
    log("  H 1s deficit / correlation = %.1f x" % (h1s_deficit / 0.12))

    # Total cross-center V_ee contribution
    total_cross_vee = 0.0
    for i_a in range(nA):
        for i_b in range(len(states_B)):
            j_b = i_b + nA
            val = mol_eri._eri.get((i_a, j_b, i_a, j_b), 0.0)
            total_cross_vee += val
    log("\n  Total cross-center V_ee (sum of J) = %.6f Ha" % total_cross_vee)
    log("  H 1s deficit = %.6f Ha" % h1s_deficit)
    log("  Ratio V_ee / deficit = %.4f" % (
        total_cross_vee / h1s_deficit if h1s_deficit else 0))

    # -----------------------------------------------------------------------
    # Conclusion
    # -----------------------------------------------------------------------
    log("\n" + "=" * 80)
    log("CONCLUSION")
    log("=" * 80)

    E_eq_so4 = pes_results[R_EQ_EXPT]['E_mol']
    bound_so4 = E_eq_so4 < atoms['E_atoms']
    E_eq_ref = ref_results[R_EQ_EXPT]['E_mol']

    log("\n  E_mol(SO4 V_ee, R=3.015)  = %.6f Ha" % E_eq_so4)
    log("  E_mol(no V_ee, R=3.015)   = %.6f Ha" % E_eq_ref)
    log("  E_atoms                   = %.6f Ha" % atoms['E_atoms'])
    log("  Bound with SO4 V_ee?      %s" % ('YES' if bound_so4 else 'NO'))

    if not bound_so4:
        log("\n  NEGATIVE RESULT: SO(4) CG cross-center V_ee does NOT rescue")
        log("  the single-p0 Sturmian basis from the one-electron deficit.")
        log("  The H 1s deficit (+%.2f Ha) is %.0fx larger than the cross-center" % (
            h1s_deficit, h1s_deficit / max(abs(total_cross_vee), 0.001)))
        log("  V_ee correction (~%.3f Ha). The single-S3 assumption (shared p0)" % (
            abs(total_cross_vee)))
        log("  is fundamentally incompatible with heteronuclear molecules.")
        log("\n  CONCLUSION: Molecular Sturmians (atom-dependent p0 or adaptive")
        log("  energy-shell) are required for correct heteronuclear LiH binding.")
    else:
        log("\n  POSITIVE RESULT: SO(4) CG V_ee produces a bound molecule.")
        log("  This is unexpected — investigate further.")

    elapsed = time.perf_counter() - t_start
    log("\nTotal time: %.1fs" % elapsed)

    # Write results
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print("\nResults written to %s" % outpath)


if __name__ == '__main__':
    main()
