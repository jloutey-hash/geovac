#!/usr/bin/env python
"""Sprint 7: Balanced coupled FCI for second-row molecules NaH and MgH2.

Follows the same methodology as Track CE (LiH balanced FCI):
- Build balanced Hamiltonian at each R-point
- Override nuclear_repulsion with V_NN = sum(Z_A*Z_B/R_AB) for FCI
- Sector-restricted FCI in the (N_up, N_down) sector
- Compute PES, R_eq, D_e, and error vs experiment

NaH: [Ne] frozen core, 1 bond block, 2 valence electrons, Q=20
MgH2: [Ne] frozen core, 2 bond blocks, 4 valence electrons, Q=40

Nuclear repulsion convention: For FCI where core electrons are frozen
and cross-center V_ne is included in the balanced Hamiltonian, we use
V_NN = sum(Z_i * Z_j / R_ij) over all NUCLEAR pairs (full charges).
This avoids double-counting with E_core and V_cross already in the spec.
"""

import json
import os
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec, mgh2_spec


# ===========================================================================
# Molecular constants
# ===========================================================================

# NaH: Z_Na=11, Z_H=1, R_eq=3.566 bohr (Huber & Herzberg)
# D_e = 0.073 Ha (experimental, ~47 kcal/mol)
# Near-exact total energy ~ -162.389 Ha (CCSD(T)/CBS literature)
NAH_Z_A = 11.0
NAH_Z_B = 1.0
NAH_R_EQ_EXACT = 3.566  # bohr
NAH_DE_EXACT = 0.073    # Ha

# MgH2: Z_Mg=12, Z_H=1, linear H-Mg-H, R_eq=3.261 bohr
# D_e(MgH2 -> Mg + H2) ~ 0.069 Ha, D_e(MgH2 -> Mg + 2H) ~ 0.238 Ha
# For PES comparison: atomization energy (Mg + 2H)
MGH2_Z_MG = 12.0
MGH2_Z_H = 1.0
MGH2_R_EQ_EXACT = 3.261  # bohr
MGH2_DE_EXACT = 0.069    # Ha (-> Mg + H2, approximate)


def compute_nah_fci(R_val: float, verbose: bool = False) -> dict:
    """Build balanced NaH Hamiltonian and compute 2e FCI."""
    t0 = time.perf_counter()

    spec = nah_spec(R=R_val, include_pk=False, max_n=2)
    ham = build_balanced_hamiltonian(spec, R=R_val, verbose=False)

    # Override nuclear repulsion: V_NN = Z_Na * Z_H / R
    V_NN = NAH_Z_A * NAH_Z_B / R_val
    ham['nuclear_repulsion'] = V_NN

    fci = coupled_fci_energy(ham, n_electrons=2, verbose=verbose)

    dt = time.perf_counter() - t0
    return {
        'R': R_val,
        'E_fci': fci['E_coupled'],
        'V_NN': V_NN,
        'N_pauli': ham['N_pauli'],
        'time_s': dt,
    }


def compute_mgh2_fci(R_val: float, verbose: bool = False) -> dict:
    """Build balanced MgH2 Hamiltonian and compute 4e FCI.

    Linear H-Mg-H: Mg at origin, H at +/-R on z-axis.
    V_NN = Z_Mg*Z_H/R + Z_Mg*Z_H/R + Z_H*Z_H/(2R)
         = 2*12/R + 1/(2R) = 24.5/R
    """
    t0 = time.perf_counter()

    spec = mgh2_spec(R=R_val, include_pk=False, max_n=2)
    ham = build_balanced_hamiltonian(spec, R=R_val, verbose=False)

    # Override nuclear repulsion: all nuclear pairs
    V_NN = (MGH2_Z_MG * MGH2_Z_H / R_val +    # Mg-H1
            MGH2_Z_MG * MGH2_Z_H / R_val +    # Mg-H2
            MGH2_Z_H * MGH2_Z_H / (2 * R_val))  # H1-H2
    ham['nuclear_repulsion'] = V_NN

    fci = coupled_fci_energy(ham, n_electrons=4, verbose=verbose)

    dt = time.perf_counter() - t0
    return {
        'R': R_val,
        'E_fci': fci['E_coupled'],
        'V_NN': V_NN,
        'N_pauli': ham['N_pauli'],
        'time_s': dt,
    }


def analyze_pes(R_points, energies, R_eq_exact, D_e_exact, Z_A, Z_B,
                label, n_H=1):
    """Analyze PES: find R_eq, D_e, compute errors."""
    E_arr = np.array(energies)
    R_arr = np.array(R_points)

    # Find minimum
    i_min = np.argmin(E_arr)
    E_eq = E_arr[i_min]
    R_eq_grid = R_arr[i_min]

    # Parabolic fit around minimum for better R_eq
    R_eq_fit = R_eq_grid
    E_eq_fit = E_eq
    if 0 < i_min < len(R_arr) - 1:
        R_3 = R_arr[i_min - 1:i_min + 2]
        E_3 = E_arr[i_min - 1:i_min + 2]
        try:
            coeffs = np.polyfit(R_3, E_3, 2)
            if coeffs[0] > 0:  # upward parabola
                R_eq_fit = -coeffs[1] / (2 * coeffs[0])
                E_eq_fit = np.polyval(coeffs, R_eq_fit)
        except Exception:
            pass

    # D_e: use largest-R point as dissociation limit
    E_dissoc = E_arr[-1]
    R_dissoc = R_arr[-1]
    D_e = E_dissoc - E_eq_fit  # positive if bound

    # Check if PES is monotonically decreasing (overattraction)
    monotone_decrease = all(E_arr[i+1] <= E_arr[i] for i in range(len(E_arr)-1))
    has_minimum = i_min > 0 and i_min < len(R_arr) - 1

    # Errors
    R_eq_err_pct = abs((R_eq_fit - R_eq_exact) / R_eq_exact) * 100
    D_e_err_pct = abs((D_e - D_e_exact) / D_e_exact) * 100 if D_e_exact > 0 else float('inf')

    result = {
        'label': label,
        'R_eq_grid': float(R_eq_grid),
        'R_eq_fit': float(R_eq_fit),
        'E_eq': float(E_eq),
        'E_eq_fit': float(E_eq_fit),
        'E_dissoc': float(E_dissoc),
        'R_dissoc': float(R_dissoc),
        'D_e': float(D_e),
        'D_e_exact': D_e_exact,
        'R_eq_exact': R_eq_exact,
        'R_eq_err_pct': float(R_eq_err_pct),
        'D_e_err_pct': float(D_e_err_pct),
        'has_minimum': has_minimum,
        'monotone_decrease': monotone_decrease,
        'bound': D_e > 0 and has_minimum,
    }

    return result


def print_pes_table(R_points, energies, V_NNs, label):
    """Print PES table."""
    print(f"\n{'R (bohr)':>10s} {'E_FCI (Ha)':>16s} {'V_NN (Ha)':>12s} "
          f"{'E_elec (Ha)':>14s}")
    print("-" * 56)
    E_min = min(energies)
    for R, E, V in zip(R_points, energies, V_NNs):
        E_elec = E - V
        marker = " *" if abs(E - E_min) < 1e-10 else ""
        print(f"{R:10.3f} {E:16.10f} {V:12.6f} {E_elec:14.10f}{marker}")


def main():
    all_results = {}
    t_start = time.perf_counter()

    # ==================================================================
    # 1. NaH balanced FCI
    # ==================================================================
    print("=" * 70)
    print("NaH Balanced Coupled FCI (n_max=2, Q=20, 2 valence electrons)")
    print("=" * 70)
    print(f"V_NN = Z_Na*Z_H/R = {NAH_Z_A}*{NAH_Z_B}/R")
    print(f"Experimental: R_eq = {NAH_R_EQ_EXACT} bohr, D_e = {NAH_DE_EXACT} Ha")

    nah_R_points = [2.0, 2.5, 3.0, 3.3, 3.566, 3.8, 4.0, 4.5, 5.0,
                    6.0, 8.0, 10.0, 15.0, 20.0]
    nah_energies = []
    nah_vnns = []
    nah_times = []

    for R_val in nah_R_points:
        print(f"  R = {R_val:.3f} bohr ... ", end="", flush=True)
        try:
            result = compute_nah_fci(R_val)
            nah_energies.append(result['E_fci'])
            nah_vnns.append(result['V_NN'])
            nah_times.append(result['time_s'])
            print(f"E = {result['E_fci']:.10f} Ha  ({result['time_s']:.1f}s)")
        except Exception as e:
            print(f"FAILED: {e}")
            traceback.print_exc()
            nah_energies.append(float('nan'))
            nah_vnns.append(NAH_Z_A * NAH_Z_B / R_val)
            nah_times.append(0.0)

    # PES table
    print_pes_table(nah_R_points, nah_energies, nah_vnns, "NaH")

    # Analysis
    # Use main PES points (up to R=20) for analysis
    valid = [i for i in range(len(nah_energies))
             if not np.isnan(nah_energies[i])]
    nah_R_valid = [nah_R_points[i] for i in valid]
    nah_E_valid = [nah_energies[i] for i in valid]

    nah_analysis = analyze_pes(
        nah_R_valid, nah_E_valid,
        NAH_R_EQ_EXACT, NAH_DE_EXACT,
        NAH_Z_A, NAH_Z_B, "NaH"
    )

    print(f"\n--- NaH Analysis ---")
    print(f"Has minimum: {nah_analysis['has_minimum']}")
    print(f"Monotone decrease: {nah_analysis['monotone_decrease']}")
    if nah_analysis['has_minimum']:
        print(f"R_eq (grid): {nah_analysis['R_eq_grid']:.3f} bohr")
        print(f"R_eq (fit):  {nah_analysis['R_eq_fit']:.4f} bohr "
              f"(exact: {NAH_R_EQ_EXACT}, err: {nah_analysis['R_eq_err_pct']:.1f}%)")
        print(f"E_eq (fit):  {nah_analysis['E_eq_fit']:.10f} Ha")
        print(f"D_e:         {nah_analysis['D_e']:.6f} Ha "
              f"(exact: {NAH_DE_EXACT}, err: {nah_analysis['D_e_err_pct']:.1f}%)")
        print(f"Bound: {nah_analysis['bound']}")
    else:
        print(f"NO EQUILIBRIUM FOUND (overattraction)")
        print(f"Energy at R_eq_exact: {nah_energies[nah_R_points.index(NAH_R_EQ_EXACT)] if NAH_R_EQ_EXACT in nah_R_points else 'N/A'}")
        print(f"PES minimum at R = {nah_analysis['R_eq_grid']:.3f} bohr "
              f"(smallest R scanned)")

    all_results['NaH'] = {
        'n_max': 2,
        'Q': 20,
        'n_electrons': 2,
        'R_points': nah_R_points,
        'E_fci': nah_energies,
        'V_NN': nah_vnns,
        'times_s': nah_times,
        'analysis': nah_analysis,
    }

    # ==================================================================
    # 2. MgH2 balanced FCI
    # ==================================================================
    print("\n\n" + "=" * 70)
    print("MgH2 Balanced Coupled FCI (n_max=2, Q=40, 4 valence electrons)")
    print("=" * 70)
    print(f"V_NN = 2*Z_Mg*Z_H/R + Z_H*Z_H/(2R) = 24.5/R")
    print(f"Experimental: R_eq = {MGH2_R_EQ_EXACT} bohr, D_e = {MGH2_DE_EXACT} Ha")

    mgh2_R_points = [2.0, 2.5, 2.8, 3.0, 3.261, 3.5, 3.8, 4.0, 4.5, 5.0,
                     6.0, 8.0, 10.0, 15.0]
    mgh2_energies = []
    mgh2_vnns = []
    mgh2_times = []

    for R_val in mgh2_R_points:
        print(f"  R = {R_val:.3f} bohr ... ", end="", flush=True)
        try:
            result = compute_mgh2_fci(R_val)
            mgh2_energies.append(result['E_fci'])
            mgh2_vnns.append(result['V_NN'])
            mgh2_times.append(result['time_s'])
            print(f"E = {result['E_fci']:.10f} Ha  ({result['time_s']:.1f}s)")
        except Exception as e:
            print(f"FAILED: {e}")
            traceback.print_exc()
            mgh2_energies.append(float('nan'))
            mgh2_vnns.append(24.5 / R_val)
            mgh2_times.append(0.0)

    # PES table
    mgh2_vnns_actual = [24.5 / R for R in mgh2_R_points]
    print_pes_table(mgh2_R_points, mgh2_energies, mgh2_vnns, "MgH2")

    # Analysis
    valid = [i for i in range(len(mgh2_energies))
             if not np.isnan(mgh2_energies[i])]
    mgh2_R_valid = [mgh2_R_points[i] for i in valid]
    mgh2_E_valid = [mgh2_energies[i] for i in valid]

    mgh2_analysis = analyze_pes(
        mgh2_R_valid, mgh2_E_valid,
        MGH2_R_EQ_EXACT, MGH2_DE_EXACT,
        MGH2_Z_MG, MGH2_Z_H, "MgH2", n_H=2
    )

    print(f"\n--- MgH2 Analysis ---")
    print(f"Has minimum: {mgh2_analysis['has_minimum']}")
    print(f"Monotone decrease: {mgh2_analysis['monotone_decrease']}")
    if mgh2_analysis['has_minimum']:
        print(f"R_eq (grid): {mgh2_analysis['R_eq_grid']:.3f} bohr")
        print(f"R_eq (fit):  {mgh2_analysis['R_eq_fit']:.4f} bohr "
              f"(exact: {MGH2_R_EQ_EXACT}, err: {mgh2_analysis['R_eq_err_pct']:.1f}%)")
        print(f"E_eq (fit):  {mgh2_analysis['E_eq_fit']:.10f} Ha")
        print(f"D_e:         {mgh2_analysis['D_e']:.6f} Ha "
              f"(exact: {MGH2_DE_EXACT}, err: {mgh2_analysis['D_e_err_pct']:.1f}%)")
        print(f"Bound: {mgh2_analysis['bound']}")
    else:
        print(f"NO EQUILIBRIUM FOUND (overattraction)")
        print(f"PES minimum at R = {mgh2_analysis['R_eq_grid']:.3f} bohr "
              f"(smallest R scanned)")

    all_results['MgH2'] = {
        'n_max': 2,
        'Q': 40,
        'n_electrons': 4,
        'R_points': mgh2_R_points,
        'E_fci': mgh2_energies,
        'V_NN': mgh2_vnns,
        'times_s': mgh2_times,
        'analysis': mgh2_analysis,
    }

    # ==================================================================
    # 3. Comparison summary
    # ==================================================================
    print("\n\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)

    print(f"\n{'Molecule':>10s} {'Q':>4s} {'N_el':>5s} {'N_pauli':>8s} "
          f"{'Has min':>8s} {'R_eq fit':>10s} {'R_eq err':>9s} "
          f"{'D_e (Ha)':>10s} {'D_e err':>9s}")
    print("-" * 85)

    # LiH reference (from Track CE)
    print(f"{'LiH (ref)':>10s} {'30':>4s} {'4':>5s} {'878':>8s} "
          f"{'Yes':>8s} {'3.227':>10s} {'7.0%':>9s} "
          f"{'0.037':>10s} {'60%':>9s}")

    for mol in ['NaH', 'MgH2']:
        d = all_results[mol]
        a = d['analysis']
        Q = d['Q']
        n_el = d['n_electrons']
        N_pauli = '?'
        # Get Pauli count from first valid point
        for i, E in enumerate(d['E_fci']):
            if not np.isnan(E):
                break

        has_min = 'Yes' if a['has_minimum'] else 'No'
        R_eq_str = f"{a['R_eq_fit']:.3f}" if a['has_minimum'] else "N/A"
        R_err_str = f"{a['R_eq_err_pct']:.1f}%" if a['has_minimum'] else "N/A"
        D_e_str = f"{a['D_e']:.4f}" if a['bound'] else "unbound"
        D_err_str = f"{a['D_e_err_pct']:.1f}%" if a['bound'] else "N/A"

        print(f"{mol:>10s} {Q:>4d} {n_el:>5d} {'':>8s} "
              f"{has_min:>8s} {R_eq_str:>10s} {R_err_str:>9s} "
              f"{D_e_str:>10s} {D_err_str:>9s}")

    t_total = time.perf_counter() - t_start
    print(f"\nTotal wall time: {t_total:.1f}s")

    # ==================================================================
    # 4. Save results
    # ==================================================================
    # Convert any numpy types for JSON serialization
    def sanitize(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, (np.bool_,)):
            return bool(obj)
        elif isinstance(obj, bool):
            return obj
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: sanitize(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [sanitize(v) for v in obj]
        return obj

    out_path = os.path.join(os.path.dirname(__file__),
                            'data', 'sprint7_balanced_second_row.json')
    with open(out_path, 'w') as f:
        json.dump(sanitize(all_results), f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
