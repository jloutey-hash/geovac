"""
Track DF Sprint 5: Heterogeneous Nested LiH Investigation
==========================================================

Scans Z_val to find optimal valence orbital exponent, then runs FCI PES
at the optimal Z_val.  Produces comparison table against all architectures.

Usage: python debug/track_df_heterogeneous.py
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / 'tests'))

from test_nested_hyperspherical import (
    _build_nested_lih_heterogeneous,
    _build_nested_lih,
    _one_norm,
)
from geovac.coupled_composition import coupled_fci_energy
from geovac.composed_qubit import build_composed_hamiltonian, lih_spec


def fci_energy(result, n_electrons=4):
    """Compute FCI ground-state energy."""
    fci = coupled_fci_energy(result, n_electrons)
    return fci['E_coupled']


def scan_z_val(R=3.015, max_n=2, z_vals=None):
    """Scan Z_val at fixed R, return energies."""
    if z_vals is None:
        z_vals = [0.5, 0.8, 1.0, 1.2, 1.5, 1.7, 2.0, 2.5]

    print(f"\n{'='*60}")
    print(f"Z_val SCAN at R = {R} bohr, max_n = {max_n}")
    print(f"{'='*60}")
    print(f"{'Z_val':>8} {'E_FCI (Ha)':>12} {'N_Pauli':>10} {'Q':>6} "
          f"{'S_min':>8} {'1-norm_ni':>10}")
    print(f"{'-'*60}")

    results = []
    for z_val in z_vals:
        t0 = time.perf_counter()
        res = _build_nested_lih_heterogeneous(
            max_n=max_n, R=R, Z_core=3.0, Z_val=z_val,
        )
        E = fci_energy(res)
        _, lnorm_ni = _one_norm(res['qubit_op'])
        s_min = float(min(res['S_eigenvalues']))
        dt = time.perf_counter() - t0

        print(f"{z_val:8.2f} {E:12.6f} {res['N_pauli']:10d} {res['Q']:6d} "
              f"{s_min:8.4f} {lnorm_ni:10.2f}   ({dt:.1f}s)")

        results.append({
            'Z_val': z_val,
            'E_FCI': E,
            'N_pauli': res['N_pauli'],
            'Q': res['Q'],
            'S_min_eigenvalue': s_min,
            'one_norm_ni': lnorm_ni,
            'wall_time_s': dt,
        })

    # Find optimal Z_val
    best = min(results, key=lambda r: r['E_FCI'])
    print(f"\nOptimal Z_val = {best['Z_val']:.2f} "
          f"(E = {best['E_FCI']:.6f} Ha)")

    return results, best['Z_val']


def run_pes(Z_val_opt, max_n=2, R_values=None):
    """Run FCI PES at optimal Z_val."""
    if R_values is None:
        R_values = [1.5, 2.0, 2.5, 3.015, 4.0, 5.0, 6.0, 8.0]

    print(f"\n{'='*60}")
    print(f"PES SCAN at Z_val = {Z_val_opt:.2f}, max_n = {max_n}")
    print(f"{'='*60}")
    print(f"{'R (bohr)':>10} {'E_FCI (Ha)':>12} {'N_Pauli':>10} "
          f"{'1-norm_ni':>10}")
    print(f"{'-'*46}")

    pes_data = []
    for R in R_values:
        t0 = time.perf_counter()
        res = _build_nested_lih_heterogeneous(
            max_n=max_n, R=R, Z_core=3.0, Z_val=Z_val_opt,
        )
        E = fci_energy(res)
        _, lnorm_ni = _one_norm(res['qubit_op'])
        dt = time.perf_counter() - t0

        print(f"{R:10.3f} {E:12.6f} {res['N_pauli']:10d} "
              f"{lnorm_ni:10.2f}   ({dt:.1f}s)")

        pes_data.append({
            'R': R,
            'E_FCI': E,
            'N_pauli': res['N_pauli'],
            'one_norm_ni': lnorm_ni,
        })

    # Extract PES metrics
    energies = [d['E_FCI'] for d in pes_data]
    R_vals = [d['R'] for d in pes_data]
    i_min = np.argmin(energies)
    E_min = energies[i_min]
    R_min = R_vals[i_min]
    E_inf = energies[-1]
    D_e = E_inf - E_min

    print(f"\nR_eq = {R_min:.3f} bohr (exact: 3.015)")
    print(f"R_eq error = {abs(R_min - 3.015) / 3.015 * 100:.1f}%")
    print(f"E_min = {E_min:.6f} Ha")
    print(f"D_e = {D_e:.4f} Ha (exact: 0.092)")
    print(f"D_e error = {abs(D_e - 0.092) / 0.092 * 100:.1f}%")

    return pes_data, R_min, D_e


def comparison_table(pes_data, Z_val_opt, max_n=2):
    """Print comparison table against all architectures."""
    # Get heterogeneous metrics at R=3.015
    R_eq_data = next((d for d in pes_data if abs(d['R'] - 3.015) < 0.01), None)
    energies = [d['E_FCI'] for d in pes_data]
    R_vals = [d['R'] for d in pes_data]
    i_min = np.argmin(energies)
    R_min = R_vals[i_min]
    E_inf = energies[-1]
    D_e = E_inf - energies[i_min]

    # Build uniform nested at R=3.015 for comparison
    res_uni = _build_nested_lih(max_n=2, R=3.015)
    from openfermion import get_sparse_operator
    sH = get_sparse_operator(res_uni['qubit_op'], n_qubits=10).toarray()
    E_uni = float(np.linalg.eigvalsh(sH)[0])
    _, ln_uni = _one_norm(res_uni['qubit_op'])

    # Build composed+PK at R=3.015
    spec_c = lih_spec(max_n_core=2, max_n_val=2)
    res_c = build_composed_hamiltonian(spec_c, pk_in_hamiltonian=True)
    _, ln_c = _one_norm(res_c['qubit_op'])

    # Heterogeneous at R=3.015
    res_h = _build_nested_lih_heterogeneous(
        max_n=max_n, R=3.015, Z_core=3.0, Z_val=Z_val_opt,
    )
    _, ln_h = _one_norm(res_h['qubit_op'])

    print(f"\n{'='*72}")
    print(f"COMPARISON TABLE — LiH at max_n={max_n}")
    print(f"{'='*72}")
    print(f"{'Metric':<20} {'Uniform':>10} {'Hetero':>10} "
          f"{'Composed+PK':>12} {'Balanced':>10}")
    print(f"{'-'*72}")
    print(f"{'Q':<20} {res_uni['Q']:>10} {res_h['Q']:>10} "
          f"{res_c['Q']:>12} {30:>10}")
    print(f"{'N_Pauli':<20} {res_uni['N_pauli']:>10} {res_h['N_pauli']:>10} "
          f"{res_c['N_pauli']:>12} {878:>10}")
    print(f"{'1-norm_ni (Ha)':<20} {ln_uni:>10.2f} {ln_h:>10.2f} "
          f"{ln_c:>12.2f} {74.1:>10.1f}")
    print(f"{'R_eq error (%)':<20} {'33.7':>10} "
          f"{abs(R_min-3.015)/3.015*100:>10.1f} "
          f"{'~5':>12} {'7.0':>10}")
    print(f"{'D_e error (%)':<20} {'4.7':>10} "
          f"{abs(D_e-0.092)/0.092*100:>10.1f} "
          f"{'--':>12} {'--':>10}")
    print(f"{'Z_val':<20} {'3.0':>10} {Z_val_opt:>10.2f} "
          f"{'1.0':>12} {'N/A':>10}")
    print(f"{'='*72}")


def main():
    print("Track DF Sprint 5: Heterogeneous Nested LiH Investigation")
    print("=" * 60)

    t_start = time.perf_counter()

    # Phase 1: Z_val scan
    z_scan_results, Z_val_opt = scan_z_val()

    # Phase 2: PES at optimal Z_val
    pes_data, R_min, D_e = run_pes(Z_val_opt)

    # Phase 3: Comparison table
    comparison_table(pes_data, Z_val_opt)

    total_time = time.perf_counter() - t_start

    # Phase 4: Save data
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    outpath = data_dir / 'track_df_heterogeneous.json'

    save_data = {
        'track': 'DF',
        'sprint': 5,
        'description': 'Heterogeneous nested LiH investigation',
        'Z_val_optimal': Z_val_opt,
        'R_eq_bohr': R_min,
        'R_eq_error_pct': abs(R_min - 3.015) / 3.015 * 100,
        'D_e_Ha': D_e,
        'D_e_error_pct': abs(D_e - 0.092) / 0.092 * 100,
        'z_scan': z_scan_results,
        'pes': pes_data,
        'total_wall_time_s': total_time,
    }

    with open(outpath, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nData saved to {outpath}")
    print(f"Total wall time: {total_time:.1f}s")


if __name__ == '__main__':
    main()
