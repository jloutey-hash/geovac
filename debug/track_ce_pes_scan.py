#!/usr/bin/env python
"""Track CE: Balanced coupled LiH PES scan at n_max=2.

Computes FCI energies at standard R-points, verifies against existing
analytical data, and computes D_e from dissociation limit.

Nuclear repulsion note: The MolecularSpec's nuclear_repulsion_constant
includes E_core and V_cross (core-H nuclear attraction), which are
appropriate for the composed 2-electron valence-only treatment. For
4-electron FCI where core electrons are explicit AND cross-center V_ne
is included in the balanced Hamiltonian, we must use V_NN = Z_A*Z_B/R
only, to avoid double-counting.
"""
import json
import time
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.composed_qubit import lih_spec

# LiH constants
Z_A = 3  # Li
Z_B = 1  # H


def compute_balanced_fci(R_val: float, verbose: bool = True) -> dict:
    """Build balanced Hamiltonian and compute 4e FCI at given R.

    Uses V_NN = Z_A*Z_B/R as nuclear repulsion (not the spec's constant
    which includes E_core and V_cross that are already captured by
    the explicit core block and cross-center V_ne).
    """
    t0 = time.perf_counter()

    spec = lih_spec(max_n_core=2, max_n_val=2, R=R_val, include_pk=False)
    ham = build_balanced_hamiltonian(spec, R=R_val, verbose=False)

    # Override nuclear_repulsion for 4e FCI: use V_NN only
    V_NN = Z_A * Z_B / R_val
    ham['nuclear_repulsion'] = V_NN

    fci = coupled_fci_energy(ham, n_electrons=4, verbose=verbose)

    dt = time.perf_counter() - t0
    return {
        'E_fci': fci['E_coupled'],
        'V_NN': V_NN,
        'time_s': dt,
        'eigenvalues': fci['eigenvalues'],
    }


def main():
    # R-points matching existing data
    R_points = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0]
    # Extra points for dissociation assessment
    R_dissoc = [10.0, 15.0, 20.0, 50.0, 100.0]

    results_all = {}
    energies_standard = []

    print("=" * 70)
    print("Track CE: Balanced Coupled LiH PES Scan (n_max=2)")
    print("=" * 70)
    print(f"Nuclear repulsion: V_NN = Z_A*Z_B/R (4e FCI, no E_core/V_cross)")

    # ------------------------------------------------------------------
    # 1. Compute FCI at standard R-points
    # ------------------------------------------------------------------
    print("\n--- Standard PES scan ---")
    for R_val in R_points:
        print(f"\nR = {R_val:.3f} bohr:")
        result = compute_balanced_fci(R_val)
        E = result['E_fci']
        energies_standard.append(E)
        results_all[R_val] = result
        print(f"  E = {E:.10f} Ha  (V_NN={result['V_NN']:.6f}, {result['time_s']:.1f}s)")

    # ------------------------------------------------------------------
    # 2. Verify against existing analytical data
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Verification against existing analytical data")
    print("=" * 70)

    ref_path = os.path.join(os.path.dirname(__file__),
                            'data', 'balanced_coupled_lih_nmax2_analytical.json')
    with open(ref_path, 'r') as f:
        ref = json.load(f)

    E_ref = ref['E_analytical']
    max_diff = 0.0
    print(f"{'R (bohr)':>10s} {'E (this)':>16s} {'E (ref)':>16s} {'diff (mHa)':>12s}")
    print("-" * 58)
    for i, R_val in enumerate(R_points):
        diff = (energies_standard[i] - E_ref[i]) * 1000  # mHa
        max_diff = max(max_diff, abs(diff))
        print(f"{R_val:10.3f} {energies_standard[i]:16.10f} {E_ref[i]:16.10f} {diff:12.6f}")

    if max_diff < 0.1:  # 0.1 mHa tolerance (reference has limited precision)
        print(f"\nVERIFICATION PASSED: max diff = {max_diff:.6f} mHa (< 0.1 mHa)")
    elif max_diff < 1.0:
        print(f"\nVERIFICATION MARGINAL: max diff = {max_diff:.6f} mHa (< 1.0 mHa)")
    else:
        print(f"\nVERIFICATION FAILED: max diff = {max_diff:.6f} mHa (>= 1.0 mHa)")

    # ------------------------------------------------------------------
    # 3. Dissociation limit
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Dissociation limit assessment")
    print("=" * 70)

    energies_dissoc = {}
    for R_val in R_dissoc:
        print(f"\nR = {R_val:.1f} bohr:")
        result = compute_balanced_fci(R_val)
        E = result['E_fci']
        energies_dissoc[R_val] = E
        results_all[R_val] = result
        print(f"  E = {E:.10f} Ha  (V_NN={result['V_NN']:.6f}, {result['time_s']:.1f}s)")

    # Dissociation assessment
    E_5 = energies_standard[-1]  # R=5.0
    E_100 = energies_dissoc[100.0]

    # Electronic energy = E_total - V_NN. At R->inf, V_NN->0 and E_elec
    # converges to the sum of isolated atom energies.
    print(f"\n--- Convergence to dissociation limit ---")
    print(f"{'R (bohr)':>10s} {'E_total (Ha)':>16s} {'V_NN (Ha)':>10s} {'E_elec (Ha)':>14s}")
    print("-" * 55)
    all_dissoc_R = [5.0, 10.0, 15.0, 20.0, 50.0, 100.0]
    for R_val in all_dissoc_R:
        if R_val == 5.0:
            E_val = E_5
        else:
            E_val = energies_dissoc[R_val]
        V_NN = Z_A * Z_B / R_val
        E_elec = E_val - V_NN
        print(f"{R_val:10.1f} {E_val:16.10f} {V_NN:10.6f} {E_elec:14.10f}")

    # E_elec at R=100 is the best estimate of separated atoms
    E_elec_100 = E_100 - Z_A * Z_B / 100.0
    E_elec_50 = energies_dissoc[50.0] - Z_A * Z_B / 50.0
    print(f"\nE_elec(R=50) = {E_elec_50:.10f} Ha")
    print(f"E_elec(R=100) = {E_elec_100:.10f} Ha")
    print(f"E_elec(R=50) - E_elec(R=100) = {(E_elec_50 - E_elec_100)*1000:.3f} mHa")

    # ------------------------------------------------------------------
    # 4. Compute D_e
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Dissociation energy D_e")
    print("=" * 70)

    # Find equilibrium from standard points
    E_arr = np.array(energies_standard)
    R_arr = np.array(R_points)
    i_min = np.argmin(E_arr)
    E_eq = E_arr[i_min]
    R_eq_approx = R_arr[i_min]

    # Use parabolic fit around minimum for better R_eq
    if 0 < i_min < len(R_arr) - 1:
        R_3 = R_arr[i_min - 1:i_min + 2]
        E_3 = E_arr[i_min - 1:i_min + 2]
        coeffs = np.polyfit(R_3, E_3, 2)
        R_eq_fit = -coeffs[1] / (2 * coeffs[0])
        E_eq_fit = np.polyval(coeffs, R_eq_fit)
    else:
        R_eq_fit = R_eq_approx
        E_eq_fit = E_eq

    # D_e using electronic energies (subtract V_NN from both sides)
    E_eq_elec = E_eq_fit - Z_A * Z_B / R_eq_fit
    E_inf_elec = E_elec_100  # best dissociation limit
    D_e_true = E_inf_elec - E_eq_elec  # positive if bound

    # D_e using R=5 as practical limit (matches reference convention)
    D_e_R5 = E_5 - E_eq_fit  # E_5 > E_eq (less negative) => positive if bound
    # Actually: D_e = E(atoms) - E(molecule). Since E(R=5) > E(R_eq):
    # D_e_R5 = E(R=5) - E(R_eq_fit) should be positive for bound molecule
    # But the ref file has D_e = 0.161 = E_eq_parabolic - E(R=5) ... no wait
    # Ref: R_eq_analytical: 3.227, E_eq: -7.933, D_e: 0.161
    # E(R=5) = -7.772, so D_e = |E_eq - E(R=5)| = |-7.933 - (-7.772)| = 0.161
    D_e_R5 = abs(E_eq_fit - E_5)

    D_e_ref = ref.get('D_e', None)

    print(f"E(R_eq grid ~ {R_eq_approx:.3f}) = {E_eq:.10f} Ha")
    print(f"E(R_eq fit ~ {R_eq_fit:.4f})  = {E_eq_fit:.10f} Ha (parabolic)")
    print(f"E_elec(R_eq fit)             = {E_eq_elec:.10f} Ha")
    print(f"E_elec(R=100, dissoc limit)  = {E_inf_elec:.10f} Ha")
    print(f"\nD_e (electronic, R=100):  {D_e_true:.6f} Ha = {D_e_true * 27.2114:.4f} eV")
    print(f"D_e (practical, E_eq - E(R=5)): {D_e_R5:.6f} Ha = {D_e_R5 * 27.2114:.4f} eV")
    if D_e_true > 0:
        print(f"BOUND: D_e > 0")
    else:
        print(f"UNBOUND: D_e <= 0")

    if D_e_ref is not None:
        print(f"\nD_e from reference file: {D_e_ref:.6f} Ha (E_eq_parabolic - E(R=5))")

    # Exact LiH values for comparison
    D_e_exact = 0.0924  # Ha (experimental)
    R_eq_exact = 3.015   # bohr
    E_exact = -8.0705    # Ha (near-exact)
    print(f"\nExact reference: D_e = {D_e_exact} Ha, R_eq = {R_eq_exact} bohr, E = {E_exact} Ha")

    E_err_pct = abs((E_eq_fit - E_exact) / E_exact * 100)
    R_err_pct = abs((R_eq_fit - R_eq_exact) / R_eq_exact * 100)
    D_err_pct_true = abs((D_e_true - D_e_exact) / D_e_exact * 100)
    D_err_pct_R5 = abs((D_e_R5 - D_e_exact) / D_e_exact * 100)
    print(f"Errors: |E(R_eq)| = {E_err_pct:.2f}%, |R_eq| = {R_err_pct:.1f}%")
    print(f"        |D_e(R=100)| = {D_err_pct_true:.1f}%, |D_e(R=5)| = {D_err_pct_R5:.1f}%")

    # ------------------------------------------------------------------
    # 5. Full PES table
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Full PES Table")
    print("=" * 70)
    all_R = sorted(results_all.keys())
    print(f"{'R (bohr)':>10s} {'E (Ha)':>16s} {'V_NN (Ha)':>10s}")
    print("-" * 40)
    for R_val in all_R:
        E = results_all[R_val]['E_fci']
        V_NN = results_all[R_val]['V_NN']
        print(f"{R_val:10.3f} {E:16.10f} {V_NN:10.6f}")

    # ------------------------------------------------------------------
    # 6. Save results
    # ------------------------------------------------------------------
    output = {
        'n_max': 2,
        'Q': 30,
        'n_electrons': 4,
        'nuclear_repulsion_convention': 'V_NN = Z_A*Z_B/R only (4e FCI)',
        'R_points_standard': R_points,
        'E_standard': energies_standard,
        'R_points_dissoc': R_dissoc,
        'E_dissoc': [energies_dissoc[R] for R in R_dissoc],
        'R_eq_grid': float(R_eq_approx),
        'R_eq_fit': float(R_eq_fit),
        'E_eq': float(E_eq),
        'E_eq_fit': float(E_eq_fit),
        'E_elec_R100': float(E_inf_elec),
        'D_e_true_Ha': float(D_e_true),
        'D_e_true_eV': float(D_e_true * 27.2114),
        'D_e_R5_Ha': float(D_e_R5),
        'D_e_R5_eV': float(D_e_R5 * 27.2114),
        'D_e_exact_Ha': D_e_exact,
        'E_exact_Ha': E_exact,
        'R_eq_exact': R_eq_exact,
        'E_error_pct': float(E_err_pct),
        'R_eq_error_pct': float(R_err_pct),
        'D_e_true_error_pct': float(D_err_pct_true),
        'D_e_R5_error_pct': float(D_err_pct_R5),
        'verification_max_diff_mHa': float(max_diff),
        'verification_passed': max_diff < 0.1,
        'E_elec_R50_minus_R100_mHa': float((E_elec_50 - E_elec_100) * 1000),
    }

    out_path = os.path.join(os.path.dirname(__file__),
                            'data', 'track_ce_lih_pes_nmax2.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
