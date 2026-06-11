"""
Validate the prolate spheroidal lattice on H2+.

Success criterion: R_eq ~ 2.0 bohr, E ~ -0.603 Ha, D_e ~ 0.103 Ha
WITHOUT any fitted parameters.

Exact H2+: R_eq = 1.997 bohr, E = -0.6026 Ha, D_e = 0.1026 Ha
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    scan_h2plus_pes,
    fit_spectroscopic_constants,
)

EXACT_R_EQ = 1.997
EXACT_E = -0.6026
EXACT_DE = 0.1026


def main():
    print("=" * 72)
    print("  Prolate Spheroidal Lattice -- H2+ Validation")
    print("  Exact: R_eq=1.997, E=-0.6026 Ha, D_e=0.1026 Ha")
    print("=" * 72)

    # =================================================================
    # Test 1: Single point at R=2.0
    # =================================================================
    print("\n--- Test 1: Single point at R=2.0 bohr ---")
    lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000, xi_max=25.0)
    E_el, c2, A = lat.solve()
    V_NN = 0.5
    print(f"  c^2 = {c2:.6f}")
    print(f"  A = {A:.6f}")
    print(f"  E_elec = {E_el:.6f}")
    print(f"  E_total = {E_el + V_NN:.6f} (exact ~ -0.6026)")

    # =================================================================
    # Test 2: Radial grid convergence
    # =================================================================
    print("\n--- Test 2: Radial grid convergence at R=2.0 ---")
    print(f"  {'N_xi':>8s}  {'c^2':>10s}  {'E_total':>10s}  {'err%':>7s}")
    for N_xi in [500, 1000, 2000, 5000, 10000]:
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi, xi_max=25.0)
        E_el, c2, _ = lat.solve()
        E_tot = E_el + 0.5
        err = abs(E_tot - EXACT_E) / abs(EXACT_E) * 100
        print(f"  {N_xi:>8d}  {c2:10.6f}  {E_tot:10.6f}  {err:7.3f}")

    # =================================================================
    # Test 3: Fine PES near minimum
    # =================================================================
    print("\n--- Test 3: Fine PES scan (N_xi=8000) ---")
    R_fine = np.arange(1.0, 6.01, 0.1)
    pes = scan_h2plus_pes(R_fine, N_xi=8000, xi_max=25.0, verbose=False)
    fit = fit_spectroscopic_constants(pes['R'], pes['E_total'])

    R_err = abs(fit['R_eq'] - EXACT_R_EQ) / EXACT_R_EQ * 100
    E_err = abs(fit['E_min'] - EXACT_E) / abs(EXACT_E) * 100
    D_err = abs(fit['D_e'] - EXACT_DE) / EXACT_DE * 100

    print(f"  R_eq = {fit['R_eq']:.4f} bohr  (exact {EXACT_R_EQ}, err {R_err:.2f}%)")
    print(f"  E_min = {fit['E_min']:.6f} Ha  (exact {EXACT_E}, err {E_err:.2f}%)")
    print(f"  D_e = {fit['D_e']:.6f} Ha  (exact {EXACT_DE}, err {D_err:.2f}%)")
    print(f"  k = {fit['k']:.4f}")

    # Print PES table
    print(f"\n  {'R':>6s}  {'E_total':>10s}")
    for R, E_tot in zip(pes['R'], pes['E_total']):
        if np.isnan(E_tot):
            continue
        marker = " <--" if abs(R - fit['R_eq']) < 0.15 else ""
        print(f"  {R:6.2f}  {E_tot:10.6f}{marker}")

    # =================================================================
    # Test 4: Dissociation limit
    # =================================================================
    print("\n--- Test 4: Dissociation limit (should -> -0.5000 Ha) ---")
    for R in [10.0, 20.0, 50.0]:
        lat = ProlateSpheroidalLattice(R=R, N_xi=8000, xi_max=max(25.0, R * 3))
        E_el, _, _ = lat.solve()
        E_tot = E_el + 1.0 / R
        print(f"  R={R:5.1f}: E_total={E_tot:.6f}")

    # =================================================================
    # Summary
    # =================================================================
    print("\n" + "=" * 72)
    print(f"  SUMMARY (zero free parameters)")
    print(f"  R_eq:  {fit['R_eq']:.4f} vs {EXACT_R_EQ:.3f}  ({R_err:.2f}%)")
    print(f"  E_min: {fit['E_min']:.6f} vs {EXACT_E:.4f}  ({E_err:.2f}%)")
    print(f"  D_e:   {fit['D_e']:.6f} vs {EXACT_DE:.4f}  ({D_err:.2f}%)")

    PASS = E_err < 1.0
    print(f"\n  {'PASS' if PASS else 'FAIL'}: E_min err {'<' if E_err < 1.0 else '>'} 1%")
    print("=" * 72)

    # Save data
    outdir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, 'prolate_h2plus_pes.txt')
    with open(outfile, 'w') as f:
        f.write("# H2+ Prolate Spheroidal Lattice PES (separated solver)\n")
        f.write(f"# N_xi=8000, xi_max=25.0\n")
        f.write(f"# R_eq={fit['R_eq']:.4f} E_min={fit['E_min']:.6f} "
                f"D_e={fit['D_e']:.6f}\n")
        f.write(f"# {'R':>8s} {'E_elec':>12s} {'E_total':>12s}\n")
        for R, E_el, E_tot in zip(pes['R'], pes['E_elec'], pes['E_total']):
            if not np.isnan(E_tot):
                f.write(f"  {R:8.3f} {E_el:12.6f} {E_tot:12.6f}\n")
    print(f"\n  Data saved to {outfile}")


if __name__ == '__main__':
    main()
