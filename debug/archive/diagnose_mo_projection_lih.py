"""v0.9.30 diagnostic: MO-projected Sturmian LiH FCI.

Runs self-consistent p0 loop with use_sturmian='molecular',
then PES scan at R = 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0.
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.lattice_index import MolecularLatticeIndex, compute_bsse_correction


def run_single_R(R: float, nmax: int = 3, damping: float = 0.5) -> dict:
    """Run MO-projected Sturmian FCI at a single R."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method='slater_full',
        use_sturmian='molecular',
    )

    p0, E_mol, n_iter, converged = mol.solve_sturmian_p0(
        R=R, damping=damping, max_iter=30, tol=1e-5
    )

    return {
        'R': R, 'p0': p0, 'E_mol': E_mol,
        'n_iter': n_iter, 'converged': converged,
    }


def main() -> None:
    # First: single-point at R=3.015
    print("=" * 70)
    print("v0.9.30: MO-Projected Sturmian LiH FCI")
    print("=" * 70)

    result_eq = run_single_R(3.015)
    E_eq = result_eq['E_mol']
    print(f"\nR=3.015: E_mol = {E_eq:.6f} Ha, "
          f"p0 = {result_eq['p0']:.6f}, "
          f"converged={result_eq['converged']}")

    # Atomic reference energies (from standard FCI)
    # Li: -7.392 Ha (nmax=3, slater_full)
    # H:  -0.500 Ha (exact)
    E_atoms = -7.892  # Li + H separated

    D_e_raw = E_atoms - E_eq
    print(f"D_e (raw) = {D_e_raw:.6f} Ha")

    # PES scan
    R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
    results = []
    for R in R_values:
        if abs(R - 3.015) < 0.01:
            results.append(result_eq)
            continue
        try:
            res = run_single_R(R)
            results.append(res)
        except Exception as e:
            print(f"R={R}: FAILED — {e}")
            results.append({'R': R, 'E_mol': float('nan'),
                            'converged': False, 'p0': float('nan'),
                            'n_iter': 0})

    # Print PES table
    print("\n" + "=" * 70)
    print("PES TABLE")
    print("=" * 70)
    print(f"{'R':>6s}  {'E_mol':>10s}  {'D_e_raw':>10s}  {'p0':>8s}  {'conv':>5s}")
    for res in results:
        R = res['R']
        E = res['E_mol']
        D_e = E_atoms - E if np.isfinite(E) else float('nan')
        print(f"{R:6.3f}  {E:10.6f}  {D_e:10.6f}  "
              f"{res['p0']:8.4f}  {'Y' if res['converged'] else 'N':>5s}")

    # Diagnostics
    print("\n" + "=" * 70)
    print("SEVEN-CONFIGURATION COMPARISON at R=3.015")
    print("=" * 70)
    print(f"{'Version':>15s}  {'D_e_CP':>10s}  {'R_eq':>8s}  {'Notes':>30s}")
    print(f"{'v0.9.11':>15s}  {'+0.093':>10s}  {'~2.5':>8s}  {'Baseline (1.0% err)':>30s}")
    print(f"{'v0.9.18':>15s}  {'+0.143':>10s}  {'<2.0':>8s}  {'Hybrid SW':>30s}")
    print(f"{'v0.9.21':>15s}  {'UNBOUND':>10s}  {'---':>8s}  {'Single-p0':>30s}")
    print(f"{'v0.9.28':>15s}  {'UNBOUND':>10s}  {'---':>8s}  {'Atomic beta_k':>30s}")
    print(f"{'v0.9.29':>15s}  {'UNBOUND':>10s}  {'---':>8s}  {'MO betas unmapped':>30s}")

    D_e_v30 = E_atoms - E_eq if np.isfinite(E_eq) else float('nan')
    bound_str = f"{D_e_v30:+.3f}" if np.isfinite(D_e_v30) and D_e_v30 > 0 else "UNBOUND"
    print(f"{'v0.9.30':>15s}  {bound_str:>10s}  {'report':>8s}  {'MO projected':>30s}")
    print(f"{'Experiment':>15s}  {'+0.092':>10s}  {'3.015':>8s}  {'':>30s}")

    # Three standard diagnostics
    print("\n" + "=" * 70)
    print("THREE STANDARD DIAGNOSTICS")
    print("=" * 70)

    # 1. Dissociation: |D_e_CP(R=6)| <= 0.001?
    r6_res = [r for r in results if abs(r['R'] - 6.0) < 0.01]
    if r6_res and np.isfinite(r6_res[0]['E_mol']):
        D_e_r6 = E_atoms - r6_res[0]['E_mol']
        print(f"1. Dissociation: D_e(R=6) = {D_e_r6:.4f} Ha "
              f"{'PASS' if abs(D_e_r6) <= 0.001 else 'FAIL'}")
    else:
        print("1. Dissociation: R=6 not available")

    # 2. Bound minimum with D_e > 0?
    if np.isfinite(D_e_v30):
        print(f"2. Bound: D_e = {D_e_v30:.4f} Ha "
              f"{'PASS' if D_e_v30 > 0 else 'FAIL'}")
    else:
        print("2. Bound: E_mol not finite — FAIL")

    # 3. R_eq in [2.5, 3.5]?
    E_vals = [(r['R'], r['E_mol']) for r in results if np.isfinite(r['E_mol'])]
    if len(E_vals) >= 3:
        R_min = min(E_vals, key=lambda x: x[1])[0]
        print(f"3. R_eq = {R_min:.3f} bohr "
              f"{'PASS' if 2.5 <= R_min <= 3.5 else 'FAIL'}")


if __name__ == '__main__':
    main()
