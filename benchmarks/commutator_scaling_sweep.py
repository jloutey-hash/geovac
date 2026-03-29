#!/usr/bin/env python
"""
Commutator-Based Trotter Bound Scaling Sweep
=============================================

Computes commutator-based (tight) Trotter error bounds for GeoVac He
at nmax=2,3,4 (and 5 if feasible).  Compares against the standard
1-norm bound and fits power-law scaling exponents.

Outputs:
  - Formatted tables to stdout
  - debug/data/commutator_bounds.json (raw data)

Usage:
    python benchmarks/commutator_scaling_sweep.py

Author: GeoVac Development Team
Date: March 2026
"""

import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

# Suppress LatticeIndex info messages
warnings.filterwarnings('ignore', category=UserWarning)

# Ensure geovac is importable
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder
from geovac.trotter_bounds import (
    pauli_1norm,
    pauli_commutator_bound,
    analyze_commutator_cost,
)


def _timer() -> float:
    return time.perf_counter()


def collect_commutator_data(
    Z: int,
    n_electrons: int,
    max_n: int,
    label: str,
) -> Dict[str, Any]:
    """Build GeoVac system and compute both Trotter bounds."""
    print(f"  {label} ... ", end="", flush=True)
    t0 = _timer()

    li = LatticeIndex(
        n_electrons=n_electrons,
        max_n=max_n,
        nuclear_charge=Z,
        vee_method='slater_full',
        h1_method='hybrid',
        fci_method='matrix',
    )
    enc = JordanWignerEncoder(li)
    qop = enc.build_qubit_operator()
    n_qubits = li.n_sp
    n_terms = len(qop.terms)

    t_build = _timer() - t0
    print(f"built Q={n_qubits}, N={n_terms:,} ({t_build:.1f}s) ... ", end="", flush=True)

    # Compute commutator bound
    t1 = _timer()
    result = analyze_commutator_cost(qop, epsilons=[1e-3, 1e-6], time=1.0)
    t_comm = _timer() - t1
    print(f"comm bound ({t_comm:.1f}s)")

    base = result['base']
    lam = pauli_1norm(qop)

    return {
        'label': label,
        'max_n': max_n,
        'n_qubits': n_qubits,
        'n_terms': n_terms,
        'one_norm': lam,
        'comm_bound': base['comm_bound'],
        'onenorm_bound': base['onenorm_bound'],
        'tightening_ratio': base['tightening_ratio'],
        'n_anticommuting_pairs': base['n_anticommuting_pairs'],
        'n_total_pairs': base['n_total_pairs'],
        'anticommuting_fraction': base['anticommuting_fraction'],
        'step_comparison': result['step_comparison'],
        'build_time_s': t_build,
        'comm_time_s': t_comm,
    }


def fit_power_law(Q_vals: List[float], metric_vals: List[float]) -> float:
    """Fit metric ~ Q^alpha, return alpha."""
    if len(Q_vals) < 2:
        return float('nan')
    log_Q = np.log(Q_vals)
    log_m = np.log(metric_vals)
    alpha, _ = np.polyfit(log_Q, log_m, 1)
    return float(alpha)


def main() -> None:
    print("=" * 70)
    print("Commutator-Based Trotter Bound Scaling Sweep")
    print("=" * 70)

    systems = [
        (2, 2, 2, "He nmax=2"),
        (2, 2, 3, "He nmax=3"),
        (2, 2, 4, "He nmax=4"),
    ]

    results: List[Dict[str, Any]] = []

    print("\nCollecting data:")
    for Z, n_e, nmax, label in systems:
        try:
            data = collect_commutator_data(Z, n_e, nmax, label)
            results.append(data)
        except Exception as e:
            print(f"FAILED: {e}")

    # nmax=5 skipped: N=227,338 terms -> N*N ~ 52 GB, infeasible.
    # Three data points (nmax=2,3,4) suffice for scaling fit.
    print("\n  He nmax=5: SKIPPED (N~227k -> N*N ~ 52 GB, extrapolate from fit)")

    # Print summary table
    print("\n" + "=" * 90)
    print(f"{'System':<14} {'Q':>4} {'N_terms':>10} {'lam':>10} "
          f"{'Anticomm%':>10} {'Tight/L^2':>10} {'r_comm':>8} {'r_1norm':>8}")
    print("-" * 90)
    for r in results:
        eps3 = [s for s in r['step_comparison'] if s['epsilon'] == 1e-3][0]
        print(f"{r['label']:<14} {r['n_qubits']:>4} {r['n_terms']:>10,} "
              f"{r['one_norm']:>10.4f} "
              f"{r['anticommuting_fraction']*100:>9.2f}% "
              f"{r['tightening_ratio']:>10.6f} "
              f"{eps3['comm_steps']:>8,} {eps3['onenorm_steps']:>8,}")
    print("=" * 90)

    # Fit scaling exponents
    if len(results) >= 2:
        Q_vals = [r['n_qubits'] for r in results]

        # 1-norm Trotter steps at eps=1e-3
        onenorm_steps = [
            [s for s in r['step_comparison'] if s['epsilon'] == 1e-3][0]['onenorm_steps']
            for r in results
        ]
        comm_steps = [
            [s for s in r['step_comparison'] if s['epsilon'] == 1e-3][0]['comm_steps']
            for r in results
        ]

        # Filter out zero comm_steps for log fit
        valid_comm = [(q, s) for q, s in zip(Q_vals, comm_steps) if s > 0]
        if len(valid_comm) >= 2:
            alpha_comm = fit_power_law(
                [x[0] for x in valid_comm],
                [x[1] for x in valid_comm],
            )
        else:
            alpha_comm = float('nan')

        alpha_1norm = fit_power_law(Q_vals, onenorm_steps)

        print(f"\nScaling exponents (r ~ Q^a, eps=1e-3, t=1):")
        print(f"  1-norm bound:      a = {alpha_1norm:.2f}")
        print(f"  Commutator bound:  a = {alpha_comm:.2f}")

        # Anticommuting fraction scaling
        fracs = [r['anticommuting_fraction'] for r in results]
        alpha_frac = fit_power_law(Q_vals, fracs)
        print(f"\nAnticommuting fraction scaling:")
        print(f"  f_anticomm ~ Q^{alpha_frac:.2f}")

    # Save to JSON
    output_dir = project_root / "debug" / "data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "commutator_bounds.json"

    # Convert for JSON serialization
    json_results = []
    for r in results:
        jr = dict(r)
        # Ensure all values are JSON-serializable
        for k, v in jr.items():
            if isinstance(v, np.integer):
                jr[k] = int(v)
            elif isinstance(v, np.floating):
                jr[k] = float(v)
        json_results.append(jr)

    output_data = {
        'description': 'Commutator-based Trotter bound scaling sweep for GeoVac He',
        'date': '2026-03-25',
        'systems': json_results,
    }
    if len(results) >= 2:
        output_data['scaling'] = {
            'alpha_1norm': alpha_1norm,
            'alpha_comm': alpha_comm,
            'alpha_anticomm_fraction': alpha_frac,
        }

    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    main()
