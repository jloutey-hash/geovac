"""Diagnostic script for MO Sturmian FCI on LiH (v0.9.31).

Checks 1-4, H1 matrix analysis, PES scan, eight-configuration comparison.
"""

import numpy as np
import time
from geovac.molecular_sturmian import (
    compute_molecular_sturmian_betas,
    compute_h1_matrix,
    compute_eri_matrix,
)
from geovac.mo_fci import MOSturmianFCI


def run_checks(Z_A: float = 3.0, Z_B: float = 1.0, R: float = 3.015) -> dict:
    """Run numerical checks 1-4."""
    p0 = np.sqrt(10.0)
    mo_results = compute_molecular_sturmian_betas(Z_A, Z_B, R, p0, nmax=3)
    H1, U = compute_h1_matrix(mo_results, Z_A, Z_B, R, p0, n_grid=40)

    results = {}

    # Check 1: H1[0,0] within 20% of -4.5 Ha
    h00 = H1[0, 0]
    err1 = abs(h00 - (-4.5)) / 4.5 * 100
    results['check1'] = {'H1_00': h00, 'error%': err1, 'pass': err1 < 20}

    # Check 2: H1[1,1] more negative than H1[0,0]?
    h11 = H1[1, 1]
    results['check2'] = {'H1_00': h00, 'H1_11': h11,
                         'more_negative': h11 < h00}

    # Check 3: H1 symmetric to 1e-6
    asym = np.max(np.abs(H1 - H1.T))
    results['check3'] = {'max_asymmetry': asym, 'pass': asym < 1e-6}

    # Check 4: V_NN
    V_NN = Z_A * Z_B / R
    V_NN_exact = 3.0 / 3.015
    results['check4'] = {'V_NN': V_NN, 'V_NN_exact': V_NN_exact,
                         'pass': abs(V_NN - V_NN_exact) < 1e-6}

    # H1 diagonal for all MOs
    results['H1_diag'] = {}
    for i, (n, m, ns, nr, beta) in enumerate(mo_results):
        results['H1_diag'][i] = {
            'n': n, 'm': m, 'n_sph': ns, 'n_rad': nr,
            'beta': beta, 'H1_ii': H1[i, i]
        }

    return results


def run_pes(R_values: list = None) -> list:
    """Run PES scan at multiple R values."""
    if R_values is None:
        R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]

    pes_results = []
    for R in R_values:
        t0 = time.perf_counter()
        fci = MOSturmianFCI(
            Z_A=3.0, Z_B=1.0, R=R, nmax=3, n_electrons=4, n_grid=35
        )
        E, p0 = fci.solve(damping=0.5, max_iter=15, verbose=False)
        dt = time.perf_counter() - t0

        V_NN = 3.0 / R
        E_atoms = -7.892
        D_e_raw = -(E - E_atoms) if E < 0 else float('nan')

        pes_results.append({
            'R': R, 'E_mol': E, 'p0': p0, 'V_NN': V_NN,
            'D_e_raw': D_e_raw, 'time': dt,
        })
        print(f"  R={R:.3f}: E={E:.6f} Ha, p0={p0:.4f}, D_e={D_e_raw:.4f}, t={dt:.1f}s")

    return pes_results


def write_results(checks: dict, pes: list, output_path: str) -> None:
    """Write all results to file."""
    with open(output_path, 'w') as f:
        f.write("MO Sturmian FCI -- LiH Diagnostics (v0.9.31)\n")
        f.write("=" * 60 + "\n\n")

        # Checks
        f.write("=== Numerical Checks ===\n\n")
        c1 = checks['check1']
        f.write(f"Check 1: H1[0,0] = {c1['H1_00']:.4f} Ha "
                f"(expect ~-4.5, error {c1['error%']:.1f}%) -- "
                f"{'PASS' if c1['pass'] else 'FAIL'}\n")

        c2 = checks['check2']
        f.write(f"Check 2: H1[0,0] = {c2['H1_00']:.4f}, "
                f"H1[1,1] = {c2['H1_11']:.4f} -- "
                f"{'2sigma MORE negative' if c2['more_negative'] else '2sigma LESS negative (expected: 1sigma=core, 2sigma=bonding)'}\n")

        c3 = checks['check3']
        f.write(f"Check 3: H1 symmetry max|H1-H1.T| = {c3['max_asymmetry']:.2e} -- "
                f"{'PASS' if c3['pass'] else 'FAIL'}\n")

        c4 = checks['check4']
        f.write(f"Check 4: V_NN = {c4['V_NN']:.6f} Ha "
                f"(exact {c4['V_NN_exact']:.6f}) -- "
                f"{'PASS' if c4['pass'] else 'FAIL'}\n")

        # H1 diagonal
        f.write("\n=== H1 Diagonal (Löwdin-orthogonalized) ===\n\n")
        f.write(f"{'MO':>3} {'n':>3} {'m':>3} {'n_sph':>5} {'n_rad':>5} "
                f"{'beta':>8} {'H1_ii':>10}\n")
        f.write("-" * 50 + "\n")
        for i, info in checks['H1_diag'].items():
            f.write(f"{i:3d} {info['n']:3d} {info['m']:3d} {info['n_sph']:5d} "
                    f"{info['n_rad']:5d} {info['beta']:8.4f} "
                    f"{info['H1_ii']:+10.4f} Ha\n")

        # PES
        f.write("\n=== PES Scan ===\n\n")
        f.write(f"{'R':>6} {'E_mol':>12} {'p0*':>8} {'V_NN':>8} "
                f"{'D_e_raw':>10} {'time':>6}\n")
        f.write("-" * 60 + "\n")
        for p in pes:
            f.write(f"{p['R']:6.3f} {p['E_mol']:+12.6f} {p['p0']:8.4f} "
                    f"{p['V_NN']:8.4f} {p['D_e_raw']:+10.4f} "
                    f"{p['time']:6.1f}s\n")

        # Eight-configuration comparison
        f.write("\n=== Eight-Configuration Comparison (R=3.015) ===\n\n")
        eq_result = next((p for p in pes if abs(p['R'] - 3.015) < 0.01), None)
        E_eq = eq_result['E_mol'] if eq_result else float('nan')
        D_eq = eq_result['D_e_raw'] if eq_result else float('nan')
        BSSE = -0.115  # approximate, from v0.9.11
        D_cp = D_eq - abs(BSSE) if np.isfinite(D_eq) else float('nan')

        configs = [
            ("v0.9.11 baseline (LCAO)", "+0.093", "~2.5"),
            ("v0.9.18 hybrid (SW+Fourier)", "+0.143", "<2.0"),
            ("v0.9.21 single-p0", "UNBOUND", "N/A"),
            ("v0.9.28 atomic beta_k", "UNBOUND", "N/A"),
            ("v0.9.29 MO unmapped", "UNBOUND", "N/A"),
            ("v0.9.30 MO projected", "UNPHYSICAL (-42.5 Ha)", "N/A"),
            (f"v0.9.31 MO FCI",
             f"{D_cp:+.3f}" if np.isfinite(D_cp) else "N/A",
             "report"),
            ("Experiment", "+0.092", "3.015"),
        ]

        f.write(f"{'Config':<35} {'D_e_CP':>12} {'R_eq':>8}\n")
        f.write("-" * 60 + "\n")
        for name, d_e, r_eq in configs:
            f.write(f"{name:<35} {d_e:>12} {r_eq:>8}\n")

        f.write(f"\nE_mol(R=3.015) = {E_eq:.6f} Ha\n")
        f.write(f"D_e_raw = {D_eq:.4f} Ha\n")
        f.write(f"BSSE ~ {BSSE:.3f} Ha (from v0.9.11, approximate)\n")
        f.write(f"D_e_CP ~ {D_cp:.4f} Ha (expt: 0.092 Ha)\n")

        # Find R_eq (minimum E in PES)
        if pes:
            min_p = min(pes, key=lambda x: x['E_mol'])
            f.write(f"R_eq ~ {min_p['R']:.3f} bohr (E_min = {min_p['E_mol']:.6f} Ha)\n")


if __name__ == '__main__':
    print("Running MO Sturmian FCI diagnostics for LiH...\n")

    print("=== Checks 1-4 ===")
    checks = run_checks()
    for k in ['check1', 'check2', 'check3', 'check4']:
        print(f"  {k}: {checks[k]}")

    print("\n=== PES Scan ===")
    pes = run_pes()

    output_path = 'debug/data/mo_fci_lih_results.txt'
    write_results(checks, pes, output_path)
    print(f"\nResults written to {output_path}")
