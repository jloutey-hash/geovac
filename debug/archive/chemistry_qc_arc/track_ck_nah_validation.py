#!/usr/bin/env python
"""Track CK: NaH Validation and Second-Row Convergence

Builds balanced coupled Hamiltonians for NaH and other second-row molecules,
computes FCI energies (2-electron for NaH, sector-restricted for others),
and compares resource metrics.

NaH has a frozen [Ne] core (10 electrons) with 2 valence electrons encoded.
The balanced coupled framework adds cross-center V_ne and cross-block ERIs.
"""
import json
import time
import sys
import os
import traceback

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.composed_qubit import nah_spec, mgh2_spec, hcl_spec, sih4_spec
from geovac.measurement_grouping import qwc_groups


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------
R_EQ_NAH = 3.566  # bohr (NIST)
D_E_NAH = 0.073   # Ha (experimental)


def build_and_measure(spec, R, nuclei=None, verbose=False):
    """Build balanced Hamiltonian and return resource metrics."""
    result = build_balanced_hamiltonian(spec, R=R, nuclei=nuclei, verbose=verbose)
    qubit_op = result['qubit_op']
    n_pauli = result['N_pauli']
    one_norm = result['one_norm']
    n_qwc = result['n_qwc_balanced']
    Q = result['Q']
    return result, {
        'Q': Q,
        'N_pauli': n_pauli,
        'one_norm': round(one_norm, 4),
        'n_qwc': n_qwc,
    }


def run_fci_2e(result, verbose=False):
    """Run 2-electron FCI on balanced Hamiltonian.

    For NaH with frozen [Ne] core, only 2 valence electrons are encoded.
    The nuclear_repulsion in the spec already includes E_core and V_cross,
    which is correct for the 2-electron frozen-core treatment.
    """
    return coupled_fci_energy(result, n_electrons=2, verbose=verbose)


# =========================================================================
# 1. NaH balanced coupled at n_max=2 and n_max=3
# =========================================================================
def section_1_nah_basic():
    """NaH resource metrics and FCI at R_eq."""
    print("=" * 70)
    print("Section 1: NaH Balanced Coupled at R_eq = {:.3f} bohr".format(R_EQ_NAH))
    print("=" * 70)

    results_sec1 = {}

    for max_n in [2, 3]:
        print(f"\n--- n_max = {max_n} ---")
        t0 = time.perf_counter()

        try:
            spec = nah_spec(max_n_val=max_n, R=R_EQ_NAH)
            result, metrics = build_and_measure(spec, R=R_EQ_NAH, verbose=False)
            dt_build = time.perf_counter() - t0

            print(f"  Q = {metrics['Q']}")
            print(f"  Pauli terms = {metrics['N_pauli']}")
            print(f"  1-norm = {metrics['one_norm']:.4f} Ha")
            print(f"  QWC groups = {metrics['n_qwc']}")
            print(f"  Build time = {dt_build:.1f}s")

            # FCI
            print(f"  Running 2-electron FCI...")
            t1 = time.perf_counter()
            fci = run_fci_2e(result, verbose=False)
            dt_fci = time.perf_counter() - t1
            E_fci = fci['E_coupled']
            print(f"  E_FCI = {E_fci:.6f} Ha")
            print(f"  FCI time = {dt_fci:.1f}s")
            print(f"  Nuclear repulsion (includes E_core) = {result['nuclear_repulsion']:.4f} Ha")

            # Diagnostic: valence-only energy
            E_valence = E_fci - result['nuclear_repulsion']
            print(f"  Valence electronic energy = {E_valence:.6f} Ha")

            results_sec1[f'nmax{max_n}'] = {
                **metrics,
                'E_fci': round(E_fci, 6),
                'E_valence': round(E_valence, 6),
                'nuclear_repulsion': round(result['nuclear_repulsion'], 4),
                'build_time_s': round(dt_build, 1),
                'fci_time_s': round(dt_fci, 1),
                'N_pauli_composed': result['N_pauli_composed'],
                'one_norm_composed': round(result['one_norm_composed'], 4),
                'pauli_ratio': round(result['pauli_ratio'], 2),
            }
        except NotImplementedError as e:
            dt_build = time.perf_counter() - t0
            print(f"  BLOCKED: {e}")
            results_sec1[f'nmax{max_n}'] = {
                'error': str(e),
                'note': 'n_max=3 requires l=2 rotation (Ivanic-Ruedenberg) not yet implemented',
            }

    return results_sec1


# =========================================================================
# 2. NaH PES scan at n_max=2
# =========================================================================
def section_2_nah_pes():
    """NaH PES scan at n_max=2."""
    print("\n" + "=" * 70)
    print("Section 2: NaH PES Scan (balanced coupled, n_max=2)")
    print("=" * 70)

    R_points = [2.5, 3.0, 3.566, 4.0, 5.0, 6.0, 8.0]
    pes_data = []

    print(f"\n{'R (bohr)':>10}  {'E_total (Ha)':>14}  {'Pauli':>7}  {'1-norm':>10}  {'QWC':>5}  {'time':>6}")
    print("-" * 65)

    for R_val in R_points:
        t0 = time.perf_counter()
        spec = nah_spec(max_n_val=2, R=R_val)
        result, metrics = build_and_measure(spec, R=R_val, verbose=False)

        fci = run_fci_2e(result, verbose=False)
        E_fci = fci['E_coupled']
        dt = time.perf_counter() - t0

        print(f"{R_val:10.3f}  {E_fci:14.6f}  {metrics['N_pauli']:7d}  "
              f"{metrics['one_norm']:10.4f}  {metrics['n_qwc']:5d}  {dt:6.1f}s")

        pes_data.append({
            'R': R_val,
            'E_fci': round(E_fci, 6),
            **metrics,
            'time_s': round(dt, 1),
        })

    # Find R_eq from minimum
    energies = [d['E_fci'] for d in pes_data]
    R_vals = [d['R'] for d in pes_data]
    idx_min = np.argmin(energies)
    E_min = energies[idx_min]
    R_min = R_vals[idx_min]

    # Parabolic fit around minimum for better R_eq estimate
    if 0 < idx_min < len(energies) - 1:
        R_fit = R_vals[idx_min - 1:idx_min + 2]
        E_fit = energies[idx_min - 1:idx_min + 2]
        coeffs = np.polyfit(R_fit, E_fit, 2)
        R_eq_fit = -coeffs[1] / (2 * coeffs[0])
        E_eq_fit = np.polyval(coeffs, R_eq_fit)
    else:
        R_eq_fit = R_min
        E_eq_fit = E_min

    # Dissociation energy estimate: D_e = E(large R) - E(R_eq)
    E_dissoc = energies[-1]  # R=8.0, should be near dissociation
    D_e_computed = E_dissoc - E_eq_fit

    R_eq_err = abs(R_eq_fit - R_EQ_NAH) / R_EQ_NAH * 100

    print(f"\n  R_eq (grid min) = {R_min:.3f} bohr")
    print(f"  R_eq (parabolic fit) = {R_eq_fit:.3f} bohr")
    print(f"  R_eq error = {R_eq_err:.1f}% (ref: {R_EQ_NAH:.3f} bohr)")
    print(f"  E(R_eq) = {E_eq_fit:.6f} Ha")
    print(f"  E(R=8.0) = {E_dissoc:.6f} Ha")
    print(f"  D_e = {D_e_computed:.6f} Ha (ref: {D_E_NAH:.3f} Ha)")

    pes_summary = {
        'R_points': R_vals,
        'energies': [round(e, 6) for e in energies],
        'R_eq_grid': round(R_min, 3),
        'R_eq_fit': round(R_eq_fit, 3),
        'R_eq_ref': R_EQ_NAH,
        'R_eq_error_pct': round(R_eq_err, 1),
        'E_eq_fit': round(E_eq_fit, 6),
        'D_e_computed': round(D_e_computed, 6),
        'D_e_ref': D_E_NAH,
        'pes_data': pes_data,
    }
    return pes_summary


# =========================================================================
# 3. Resource comparison across second-row molecules
# =========================================================================
def section_3_resource_comparison():
    """Build balanced Hamiltonians for NaH, MgH2, HCl, SiH4 at n_max=2."""
    print("\n" + "=" * 70)
    print("Section 3: Second-Row Resource Comparison (balanced coupled, n_max=2)")
    print("=" * 70)

    molecules = []

    # NaH
    try:
        spec = nah_spec(max_n_val=2)
        result, metrics = build_and_measure(spec, R=R_EQ_NAH, verbose=False)
        molecules.append({'name': 'NaH', 'n_electrons_encoded': 2, **metrics})
        print(f"  NaH:  Q={metrics['Q']:3d}  Pauli={metrics['N_pauli']:6d}  "
              f"1-norm={metrics['one_norm']:10.2f}  QWC={metrics['n_qwc']:4d}")
    except Exception as e:
        print(f"  NaH: FAILED - {e}")
        traceback.print_exc()

    # MgH2
    try:
        spec = mgh2_spec(max_n_val=2)
        result, metrics = build_and_measure(spec, R=3.268, verbose=False)
        molecules.append({'name': 'MgH2', 'n_electrons_encoded': 4, **metrics})
        print(f"  MgH2: Q={metrics['Q']:3d}  Pauli={metrics['N_pauli']:6d}  "
              f"1-norm={metrics['one_norm']:10.2f}  QWC={metrics['n_qwc']:4d}")
    except Exception as e:
        print(f"  MgH2: FAILED - {e}")
        traceback.print_exc()

    # HCl
    try:
        spec = hcl_spec(max_n_val=2)
        result, metrics = build_and_measure(spec, R=2.409, verbose=False)
        molecules.append({'name': 'HCl', 'n_electrons_encoded': 8, **metrics})
        print(f"  HCl:  Q={metrics['Q']:3d}  Pauli={metrics['N_pauli']:6d}  "
              f"1-norm={metrics['one_norm']:10.2f}  QWC={metrics['n_qwc']:4d}")
    except Exception as e:
        print(f"  HCl: FAILED - {e}")
        traceback.print_exc()

    # SiH4
    try:
        spec = sih4_spec(max_n_val=2)
        result, metrics = build_and_measure(spec, R=2.798, verbose=False)
        molecules.append({'name': 'SiH4', 'n_electrons_encoded': 8, **metrics})
        print(f"  SiH4: Q={metrics['Q']:3d}  Pauli={metrics['N_pauli']:6d}  "
              f"1-norm={metrics['one_norm']:10.2f}  QWC={metrics['n_qwc']:4d}")
    except Exception as e:
        print(f"  SiH4: FAILED - {e}")
        traceback.print_exc()

    # Scaling fit (Pauli vs Q)
    if len(molecules) >= 2:
        Q_arr = np.array([m['Q'] for m in molecules], dtype=float)
        P_arr = np.array([m['N_pauli'] for m in molecules], dtype=float)
        log_Q = np.log(Q_arr)
        log_P = np.log(P_arr)
        slope, intercept = np.polyfit(log_Q, log_P, 1)
        print(f"\n  Pauli vs Q scaling exponent: {slope:.2f}")
        print(f"  (log-log fit across {len(molecules)} molecules)")

        for m in molecules:
            m['log_Q'] = round(np.log(m['Q']), 4)
            m['log_Pauli'] = round(np.log(m['N_pauli']), 4)

        scaling = {
            'exponent': round(slope, 2),
            'intercept': round(intercept, 2),
            'n_molecules': len(molecules),
        }
    else:
        scaling = {'exponent': None, 'n_molecules': len(molecules)}

    return {'molecules': molecules, 'scaling': scaling}


# =========================================================================
# Main
# =========================================================================
def main():
    all_results = {}

    # Section 1: NaH basic
    try:
        all_results['section1_nah_basic'] = section_1_nah_basic()
    except Exception as e:
        print(f"\nSection 1 FAILED: {e}")
        traceback.print_exc()
        all_results['section1_nah_basic'] = {'error': str(e)}

    # Section 2: NaH PES
    try:
        all_results['section2_nah_pes'] = section_2_nah_pes()
    except Exception as e:
        print(f"\nSection 2 FAILED: {e}")
        traceback.print_exc()
        all_results['section2_nah_pes'] = {'error': str(e)}

    # Section 3: Resource comparison
    try:
        all_results['section3_resource_comparison'] = section_3_resource_comparison()
    except Exception as e:
        print(f"\nSection 3 FAILED: {e}")
        traceback.print_exc()
        all_results['section3_resource_comparison'] = {'error': str(e)}

    # Save results
    out_path = os.path.join(os.path.dirname(__file__), 'data', 'track_ck_nah_validation.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
