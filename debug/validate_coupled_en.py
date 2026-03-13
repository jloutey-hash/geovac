"""
Coupled Electron-Nuclear Lattice Validation for LiH
=====================================================

Computes:
1. Morse <R>_v expectation values for LiH vibrational states
2. LiH LCAO-FCI electronic energy at each R(v)
3. E_total(v) = E_elec(R(v)) + E_nuc(v)
4. Lambda derivation from force constant matching

Answers the key question:
    Does the coupled electron-nuclear energy E_total(v) have a
    minimum at v=0 corresponding to R_eq ~ 3.015 bohr?

Output:
    debug/data/coupled_en_lih.txt  — tabulated results
    debug/EQUILIBRIUM_FROM_COUPLING.md — analysis (manual)

Date: March 2026
Version: v0.9.38+
"""

import sys
import os
import time
import numpy as np

# Ensure package is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.coupled_en_lattice import (
    compute_morse_parameters,
    compute_R_expectation_table,
    morse_force_constant,
    morse_expectation_R,
    morse_classical_turning_points,
    numerical_force_constant,
    derive_lambda_from_force_constant,
    coupled_energy_table,
    find_equilibrium_v,
)
from geovac.nuclear_lattice import DIATOMIC_CONSTANTS, HARTREE_TO_CM, AMU_TO_ME

# ======================================================================
# Configuration
# ======================================================================

MOLECULE = 'LiH'
Z_A, Z_B = 3, 1
N_ELECTRONS = 4
NMAX = 3

# How many vibrational states to compute FCI for
# (full computation for low v, skip high v where R is very large)
V_MAX_COMPUTE = 10

# Additional R scan points near equilibrium for force constant fitting
R_SCAN_EXTRA = [1.5, 2.0, 2.5, 2.8, 3.0, 3.015, 3.2, 3.5, 4.0, 5.0, 6.0]

# Lambda values to test
LAMBDA_VALUES = [0.0, 0.02, 0.05, 0.10]


def compute_separated_atoms(nmax: int) -> tuple:
    """Compute Li + H separated atom energies."""
    from geovac.lattice_index import LatticeIndex

    print("Computing separated atom energies...")
    t0 = time.time()

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

    dt = time.time() - t0
    print(f"  Li: {E_li:.6f} Ha, H: {E_h:.6f} Ha, sum: {E_li + E_h:.6f} Ha  ({dt:.1f}s)")
    return E_li, E_h


def compute_lih_fci(R: float, nmax: int,
                    t_corr_lambda: float = 0.0,
                    t_corr_fock_weighted: bool = False) -> float:
    """Compute LiH FCI energy at bond length R."""
    from geovac.lattice_index import MolecularLatticeIndex

    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax, nmax_B=nmax,
        R=R,
        n_electrons=N_ELECTRONS,
        vee_method='slater_full',
        fci_method='auto',
        cross_nuclear_method='exact',
        cross_atom_vee=True,
        t_corr_lambda=t_corr_lambda,
        t_corr_fock_weighted=t_corr_fock_weighted,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return eigvals[0]


def main() -> None:
    print("=" * 70)
    print("Coupled Electron-Nuclear Lattice Validation: LiH")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Morse parameters and <R>_v table
    # ------------------------------------------------------------------
    print("\n--- Step 1: Morse Parameters ---")
    params = compute_morse_parameters(MOLECULE)
    print(f"  r_e     = {params['r_e']:.4f} bohr")
    print(f"  D_e     = {params['D_e']:.4f} Ha ({params['D_e']*27.2114:.2f} eV)")
    print(f"  a       = {params['a']:.4f} bohr^-1")
    print(f"  lambda  = {params['lam']:.2f}")
    print(f"  j       = {params['j']:.2f}")
    print(f"  v_max   = {params['v_max']}")
    print(f"  k_Morse = {params['k_morse']:.6f} Ha/bohr^2")

    # Compute <R>_v for all states
    table_all = compute_R_expectation_table(MOLECULE)
    print(f"\n  <R>_v expectation values:")
    print(f"  {'v':>3s}  {'<R>_v':>10s}  {'E_nuc(v)':>12s}  {'r_inner':>10s}  {'r_outer':>10s}")
    omega_h = params['omega_e_hartree']
    omegax_h = params['omega_e_xe_hartree']
    for row in table_all:
        v = int(row[0])
        R_v = row[1]
        E_nuc = row[2]
        r_in, r_out = morse_classical_turning_points(
            v, params['r_e'], params['a'], params['D_e'], E_nuc)
        print(f"  {v:3d}  {R_v:10.4f}  {E_nuc:12.6f}  {r_in:10.4f}  {r_out:10.4f}")

    # ------------------------------------------------------------------
    # Step 2: Separated atom energies
    # ------------------------------------------------------------------
    print("\n--- Step 2: Separated Atoms ---")
    E_li, E_h = compute_separated_atoms(NMAX)
    E_sep = E_li + E_h

    # ------------------------------------------------------------------
    # Step 3: R scan — bare PES (lambda=0) + Fock-weighted PES
    # ------------------------------------------------------------------
    # Collect all R values: from <R>_v table + extra scan points
    R_from_v = [row[1] for row in table_all[:V_MAX_COMPUTE + 1]
                if np.isfinite(row[1])]
    R_all = sorted(set(R_from_v + R_SCAN_EXTRA))

    print(f"\n--- Step 3: Electronic PES at {len(R_all)} R-points ---")
    print(f"  R range: [{min(R_all):.3f}, {max(R_all):.3f}] bohr")

    # Store results: {lambda: {R: E_elec}}
    E_elec_by_lambda: dict = {lam: {} for lam in LAMBDA_VALUES}

    for i, R in enumerate(R_all):
        print(f"  [{i+1}/{len(R_all)}] R = {R:.4f} bohr ... ", end="", flush=True)
        t0 = time.time()

        for lam in LAMBDA_VALUES:
            fock_weighted = (lam > 0)
            E = compute_lih_fci(R, NMAX, t_corr_lambda=lam,
                                t_corr_fock_weighted=fock_weighted)
            E_elec_by_lambda[lam][R] = E

        dt = time.time() - t0
        E0 = E_elec_by_lambda[0.0][R]
        print(f"E(lam=0) = {E0:.6f} Ha  ({dt:.1f}s)")

    # ------------------------------------------------------------------
    # Step 4: Coupled E_total(v) for each lambda
    # ------------------------------------------------------------------
    print("\n--- Step 4: Coupled E_total(v) ---")

    v_table = table_all[:V_MAX_COMPUTE + 1]

    for lam in LAMBDA_VALUES:
        print(f"\n  lambda = {lam:.3f}:")
        ctable = coupled_energy_table(v_table, E_elec_by_lambda[lam], E_sep)

        print(f"  {'v':>3s}  {'R(v)':>8s}  {'E_nuc':>10s}  {'E_elec':>10s}  "
              f"{'E_total':>10s}  {'D_eff':>10s}")
        for row in ctable:
            v, Rv, Enuc, Eelec, Etot, Deff = row
            print(f"  {int(v):3d}  {Rv:8.4f}  {Enuc:10.6f}  {Eelec:10.6f}  "
                  f"{Etot:10.6f}  {Deff:10.6f}")

        v_eq, E_min = find_equilibrium_v(ctable[:, 4])
        R_eq = ctable[v_eq, 1]
        print(f"  -> Minimum at v={v_eq}, R={R_eq:.4f} bohr, E={E_min:.6f} Ha")

    # ------------------------------------------------------------------
    # Step 5: Force constant matching → derive lambda
    # ------------------------------------------------------------------
    print("\n--- Step 5: Force Constant Matching ---")
    k_morse = params['k_morse']
    R_eq_expt = params['r_e']

    R_arr = np.array(R_all)
    E_bare_arr = np.array([E_elec_by_lambda[0.0][R] for R in R_all])

    # Compute force constant of bare PES
    k_bare = numerical_force_constant(R_arr, E_bare_arr, R_eq_expt)
    print(f"  k_Morse (target)   = {k_morse:.6f} Ha/bohr^2")
    print(f"  k_bare (lambda=0)  = {k_bare:.6f} Ha/bohr^2")

    # For each non-zero lambda, compute correction contribution
    for lam in LAMBDA_VALUES:
        if lam == 0:
            continue
        E_lam_arr = np.array([E_elec_by_lambda[lam][R] for R in R_all])
        E_corr_per_lam = (E_lam_arr - E_bare_arr) / lam

        result = derive_lambda_from_force_constant(
            k_morse, R_arr, E_bare_arr, E_corr_per_lam, R_eq_expt)

        print(f"\n  Using correction shape from lambda={lam:.3f}:")
        print(f"    k_correction/lam = {result['k_correction_per_lambda']:.6f} Ha/bohr^2")
        print(f"    lambda_derived   = {result['lambda_derived']:.6f}")
        print(f"    k_total (check)  = {result['k_total']:.6f} Ha/bohr^2")

    # ------------------------------------------------------------------
    # Step 6: Alternative derivation — vibrational frequency matching
    # ------------------------------------------------------------------
    print("\n--- Step 6: Vibrational Frequency Matching ---")
    mu_me = params['mu_me']
    omega_target = params['omega_e_hartree']
    print(f"  omega_e (target) = {omega_target:.6f} Ha = {params['omega_e_hartree'] * HARTREE_TO_CM:.2f} cm^-1")

    for lam in LAMBDA_VALUES:
        if lam == 0:
            continue
        E_lam_arr = np.array([E_elec_by_lambda[lam][R] for R in R_all])
        k_lam = numerical_force_constant(R_arr, E_lam_arr, R_eq_expt)
        if k_lam > 0:
            omega_lam = np.sqrt(k_lam / mu_me)
            print(f"  lambda={lam:.3f}: k={k_lam:.6f}, "
                  f"omega={omega_lam:.6f} Ha ({omega_lam * HARTREE_TO_CM:.1f} cm^-1)")
        else:
            print(f"  lambda={lam:.3f}: k={k_lam:.6f} (no minimum)")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    outdir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, 'coupled_en_lih.txt')

    with open(outfile, 'w') as f:
        f.write("# Coupled Electron-Nuclear LiH Validation\n")
        f.write(f"# nmax={NMAX}, molecule={MOLECULE}\n")
        f.write(f"# E_Li={E_li:.8f}, E_H={E_h:.8f}, E_sep={E_sep:.8f}\n")
        f.write(f"# k_Morse={k_morse:.8f} Ha/bohr^2\n")
        f.write(f"# Morse: r_e={params['r_e']}, a={params['a']:.6f}, "
                f"lam={params['lam']:.4f}, j={params['j']:.4f}\n")
        f.write("#\n")

        # Morse expectation values
        f.write("# <R>_v table:\n")
        f.write("# v  <R>_v  E_nuc  r_inner  r_outer\n")
        for row in table_all:
            v = int(row[0])
            E_nuc = row[2]
            r_in, r_out = morse_classical_turning_points(
                v, params['r_e'], params['a'], params['D_e'], E_nuc)
            f.write(f"{v:3d}  {row[1]:10.6f}  {E_nuc:12.8f}  "
                    f"{r_in:10.6f}  {r_out:10.6f}\n")

        f.write("#\n# PES data:\n")
        f.write("# R  " + "  ".join(f"E(lam={l})" for l in LAMBDA_VALUES) + "\n")
        for R in R_all:
            line = f"{R:8.4f}"
            for lam in LAMBDA_VALUES:
                line += f"  {E_elec_by_lambda[lam][R]:12.8f}"
            f.write(line + "\n")

    print(f"\n  Results saved to {outfile}")
    print("\n" + "=" * 70)
    print("Done.")


if __name__ == '__main__':
    main()
