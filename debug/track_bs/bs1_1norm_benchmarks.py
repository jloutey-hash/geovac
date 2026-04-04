"""
Track BS-1: 1-Norm Benchmarks for Composed Qubit Hamiltonians
=============================================================

Computes total and electronic-only 1-norms for H2O, LiH, and He
at multiple basis sizes. PK fraction quantifies the classical-
partitioning advantage (Track BF).

Author: GeoVac Development Team
Date: April 2026
"""

import json
import sys
import time
from pathlib import Path

# Ensure project root is on path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from geovac.trotter_bounds import pauli_1norm


def compute_composed_1norms(builder_fn, label: str, builder_kwargs_with_pk: dict,
                            builder_kwargs_no_pk: dict, verbose: bool = True) -> dict:
    """Compute total and electronic-only 1-norms for a composed system."""
    t0 = time.perf_counter()

    # With PK
    if verbose:
        print(f"\n{'='*60}")
        print(f"  {label} — WITH PK")
        print(f"{'='*60}")
    result_pk = builder_fn(**builder_kwargs_with_pk, verbose=verbose)
    norm_total = pauli_1norm(result_pk['qubit_op'])
    Q = result_pk['Q']
    N_pauli = result_pk['N_pauli']

    # Without PK
    if verbose:
        print(f"\n{'='*60}")
        print(f"  {label} — WITHOUT PK")
        print(f"{'='*60}")
    result_no_pk = builder_fn(**builder_kwargs_no_pk, verbose=verbose)
    norm_elec = pauli_1norm(result_no_pk['qubit_op'])
    Q_check = result_no_pk['Q']
    assert Q == Q_check, f"Qubit count mismatch: {Q} vs {Q_check}"

    norm_pk = norm_total - norm_elec
    pk_frac = (norm_pk / norm_total * 100) if norm_total > 0 else 0.0

    elapsed = time.perf_counter() - t0

    row = {
        'system': label,
        'Q': Q,
        'N_pauli_total': N_pauli,
        'N_pauli_elec': result_no_pk['N_pauli'],
        '1_norm_total_Ha': round(norm_total, 4),
        '1_norm_electronic_Ha': round(norm_elec, 4),
        '1_norm_PK_Ha': round(norm_pk, 4),
        'PK_fraction_pct': round(pk_frac, 2),
        'wall_time_s': round(elapsed, 1),
    }

    if verbose:
        print(f"\n  RESULT: {label}")
        print(f"    Q = {Q}, N_pauli(total) = {N_pauli}, N_pauli(elec) = {result_no_pk['N_pauli']}")
        print(f"    1-norm total    = {norm_total:.4f} Ha")
        print(f"    1-norm elec     = {norm_elec:.4f} Ha")
        print(f"    1-norm PK       = {norm_pk:.4f} Ha")
        print(f"    PK fraction     = {pk_frac:.2f}%")
        print(f"    Wall time       = {elapsed:.1f}s")

    return row


def compute_he_1norm(max_n: int, verbose: bool = True) -> dict:
    """Compute 1-norm for He (no PK partitioning — single geometry)."""
    from geovac.vqe_benchmark import build_geovac_he

    t0 = time.perf_counter()
    if verbose:
        print(f"\n{'='*60}")
        print(f"  He — max_n={max_n}")
        print(f"{'='*60}")

    spo, of_op, n_qubits, exact_e = build_geovac_he(max_n=max_n)
    norm = pauli_1norm(of_op)
    N_pauli = len(of_op.terms)
    elapsed = time.perf_counter() - t0

    row = {
        'system': f'He (n_max={max_n})',
        'Q': n_qubits,
        'N_pauli_total': N_pauli,
        'N_pauli_elec': N_pauli,
        '1_norm_total_Ha': round(norm, 4),
        '1_norm_electronic_Ha': round(norm, 4),
        '1_norm_PK_Ha': 0.0,
        'PK_fraction_pct': 0.0,
        'wall_time_s': round(elapsed, 1),
    }

    if verbose:
        print(f"    Q = {n_qubits}, N_pauli = {N_pauli}")
        print(f"    1-norm = {norm:.4f} Ha")
        print(f"    Wall time = {elapsed:.1f}s")

    return row


def main():
    from geovac.composed_qubit import build_composed_lih, build_composed_h2o

    results = []

    # ----------------------------------------------------------------
    # He at n_max=1 and n_max=2
    # ----------------------------------------------------------------
    for n in [1, 2]:
        row = compute_he_1norm(max_n=n)
        results.append(row)

    # ----------------------------------------------------------------
    # LiH at n_max=1 (max_n_core=1, max_n_val=1) and n_max=2
    # ----------------------------------------------------------------
    for n in [1, 2]:
        row = compute_composed_1norms(
            build_composed_lih,
            f'LiH (n_max={n})',
            builder_kwargs_with_pk={'max_n_core': n, 'max_n_val': n,
                                     'pk_in_hamiltonian': True},
            builder_kwargs_no_pk={'max_n_core': n, 'max_n_val': n,
                                   'pk_in_hamiltonian': False},
        )
        results.append(row)

    # ----------------------------------------------------------------
    # H2O at n_max=2
    # ----------------------------------------------------------------
    row = compute_composed_1norms(
        build_composed_h2o,
        'H2O (n_max=2)',
        builder_kwargs_with_pk={'max_n_core': 2, 'max_n_val': 2,
                                 'pk_in_hamiltonian': True},
        builder_kwargs_no_pk={'max_n_core': 2, 'max_n_val': 2,
                               'pk_in_hamiltonian': False},
    )
    results.append(row)

    # ----------------------------------------------------------------
    # Summary table
    # ----------------------------------------------------------------
    print("\n" + "=" * 100)
    print("  1-NORM BENCHMARK SUMMARY")
    print("=" * 100)
    header = (f"{'System':<20} {'Q':>4} {'N_pauli':>8} "
              f"{'1-norm_tot':>12} {'1-norm_elec':>12} {'1-norm_PK':>10} {'PK%':>7}")
    print(header)
    print("-" * 100)
    for r in results:
        line = (f"{r['system']:<20} {r['Q']:>4} {r['N_pauli_total']:>8} "
                f"{r['1_norm_total_Ha']:>12.4f} {r['1_norm_electronic_Ha']:>12.4f} "
                f"{r['1_norm_PK_Ha']:>10.4f} {r['PK_fraction_pct']:>6.2f}%")
        print(line)
    print("=" * 100)

    # ----------------------------------------------------------------
    # Save results
    # ----------------------------------------------------------------
    out_path = Path(__file__).resolve().parent / 'bs1_results.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
