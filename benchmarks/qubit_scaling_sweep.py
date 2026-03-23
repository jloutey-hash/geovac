#!/usr/bin/env python
"""
Qubit Scaling Sweep — Comprehensive Quantum Simulation Cost Metrics
====================================================================

Systematically collects Pauli term counts, QWC measurement groups,
1-norms, and Trotter step bounds across GeoVac lattice sizes and
Gaussian reference systems.  Fits power-law scaling exponents.

Outputs:
  - Formatted tables to stdout
  - benchmarks/qubit_scaling_data.json   (raw numbers for figure generation)
  - benchmarks/qubit_scaling_summary.txt (tables as plain text)

Usage:
    python benchmarks/qubit_scaling_sweep.py

Author: GeoVac Development Team
Date: March 2026
"""

import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

# Suppress LatticeIndex info messages
warnings.filterwarnings('ignore', category=UserWarning)

# Ensure geovac is importable
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import (
    JordanWignerEncoder,
    fit_pauli_scaling,
)
from geovac.measurement_grouping import analyze_measurement_cost
from geovac.trotter_bounds import analyze_trotter_cost, pauli_1norm
from geovac.gaussian_reference import h2_sto3g, he_sto3g, build_qubit_hamiltonian

# Timeout for JW transform (seconds)
JW_TIMEOUT = 600

# Skip QWC grouping above this many Pauli terms (O(N^2) too slow)
QWC_MAX_TERMS = 50000


def _timer() -> float:
    """Return current time in seconds."""
    return time.perf_counter()


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_geovac_system(
    label: str,
    Z: int,
    n_electrons: int,
    max_n: int,
) -> Optional[Dict[str, Any]]:
    """
    Build a GeoVac LatticeIndex, encode, and collect all metrics.

    Returns None if the system times out.
    """
    print(f"  {label} ... ", end="", flush=True)
    t0 = _timer()

    # Build lattice (always use slater_full so _eri dict is available
    # for JordanWignerEncoder, even for 1-electron systems)
    li = LatticeIndex(
        n_electrons=n_electrons,
        max_n=max_n,
        nuclear_charge=Z,
        vee_method='slater_full',
        h1_method='hybrid',
        fci_method='matrix',
    )
    t_lattice = _timer() - t0

    # Encode
    t1 = _timer()
    enc = JordanWignerEncoder(li)
    analysis = enc.analyze()
    t_encode = _timer() - t1

    if t_encode > JW_TIMEOUT:
        print(f"TIMEOUT ({t_encode:.0f}s > {JW_TIMEOUT}s)")
        return None

    qop = enc.build_qubit_operator()

    # Measurement grouping (skip if too many terms — O(N^2) is prohibitive)
    t2 = _timer()
    if analysis.n_pauli_terms <= QWC_MAX_TERMS:
        meas = analyze_measurement_cost(qop)
        n_qwc = meas.n_qwc_groups
        grp_ratio = n_qwc / max(1, meas.n_pauli_terms)
    else:
        print(f"(skipping QWC: {analysis.n_pauli_terms:,} terms > {QWC_MAX_TERMS:,}) ", end="")
        n_qwc = None
        grp_ratio = None
    t_meas = _timer() - t2

    # Trotter bounds (cheap — just iterates coefficients)
    t3 = _timer()
    trotter = analyze_trotter_cost(qop)
    t_trotter = _timer() - t3

    total = _timer() - t0
    print(f"done ({total:.1f}s)")

    return {
        'label': label,
        'type': 'geovac',
        'Z': Z,
        'n_electrons': n_electrons,
        'max_n': max_n,
        'Q': analysis.n_qubits,
        'n_pauli_terms': analysis.n_pauli_terms,
        'max_pauli_weight': analysis.max_pauli_weight,
        'eri_density': analysis.eri_density,
        'n_qwc_groups': n_qwc,
        'groups_per_term': grp_ratio,
        'one_norm': trotter.one_norm,
        'lambda_per_Q': trotter.one_norm / max(1, analysis.n_qubits),
        'max_coefficient': trotter.max_coefficient,
        'trotter_steps_eps3': trotter.trotter_steps_eps3,
        'trotter_steps_eps6': trotter.trotter_steps_eps6,
        'time_lattice_s': round(t_lattice, 2),
        'time_encode_s': round(t_encode, 2),
        'time_meas_s': round(t_meas, 2),
        'time_total_s': round(total, 2),
    }


def collect_gaussian_system(
    label: str,
    system_fn,
) -> Dict[str, Any]:
    """Build a Gaussian reference system and collect all metrics."""
    print(f"  {label} ... ", end="", flush=True)
    t0 = _timer()

    sys_data = system_fn()
    _, qop, n_pauli = build_qubit_hamiltonian(sys_data)

    n_qubits = 2 * sys_data['n_spatial']
    n_spatial = sys_data['n_spatial']

    # ERI density: for Gaussian, it's the fraction of nonzero eri entries
    eri = sys_data['eri']
    n_eri_nz = int(np.count_nonzero(np.abs(eri) > 1e-12))
    n_eri_possible = n_spatial ** 4
    eri_density = n_eri_nz / max(1, n_eri_possible)

    # Max Pauli weight
    max_w = 0
    for term in qop.terms:
        if len(term) > max_w:
            max_w = len(term)

    # Measurement grouping
    meas = analyze_measurement_cost(qop)

    # Trotter bounds
    trotter = analyze_trotter_cost(qop)

    total = _timer() - t0
    print(f"done ({total:.1f}s)")

    return {
        'label': label,
        'type': 'gaussian',
        'basis': sys_data['description'],
        'Q': n_qubits,
        'n_pauli_terms': n_pauli,
        'max_pauli_weight': max_w,
        'eri_density': eri_density,
        'n_qwc_groups': meas.n_qwc_groups,
        'groups_per_term': meas.n_qwc_groups / max(1, n_pauli),
        'one_norm': trotter.one_norm,
        'lambda_per_Q': trotter.one_norm / max(1, n_qubits),
        'max_coefficient': trotter.max_coefficient,
        'trotter_steps_eps3': trotter.trotter_steps_eps3,
        'trotter_steps_eps6': trotter.trotter_steps_eps6,
        'literature_energy': sys_data['literature_energy'],
        'time_total_s': round(total, 2),
    }


# ---------------------------------------------------------------------------
# Scaling law fits
# ---------------------------------------------------------------------------

def fit_scaling_laws(
    results: List[Dict[str, Any]],
) -> Dict[str, Dict[str, float]]:
    """
    Fit power-law exponents for each metric vs Q.

    Returns dict mapping metric name to {exponent, prefactor, r_squared}.
    """
    fits = {}

    metrics = [
        ('pauli_terms', 'n_pauli_terms'),
        ('qwc_groups', 'n_qwc_groups'),
        ('one_norm', 'one_norm'),
        ('trotter_eps3', 'trotter_steps_eps3'),
        ('trotter_eps6', 'trotter_steps_eps6'),
    ]

    for name, key in metrics:
        # Filter out results where the metric is None (e.g. skipped QWC)
        valid = [(r['Q'], r[key]) for r in results if r[key] is not None]
        if len(valid) < 3:
            continue

        q_arr = np.array([v[0] for v in valid], dtype=float)
        vals = np.array([v[1] for v in valid], dtype=float)

        # Skip if any values are zero or negative (can't take log)
        if np.any(vals <= 0):
            continue

        exponent, prefactor = fit_pauli_scaling(q_arr, vals)

        # Compute R^2
        log_q = np.log(q_arr)
        log_v = np.log(vals)
        log_pred = exponent * log_q + np.log(prefactor)
        ss_res = np.sum((log_v - log_pred) ** 2)
        ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        fits[name] = {
            'exponent': round(exponent, 3),
            'prefactor': round(prefactor, 4),
            'r_squared': round(r_squared, 4),
        }

    return fits


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def format_master_table(results: List[Dict[str, Any]]) -> str:
    """Format the master results table."""
    header = (
        f"{'System':<22} {'Q':>4} {'Pauli':>9} {'ERI%':>7} "
        f"{'QWC_grp':>8} {'grp/trm':>8} {'lambda':>9} "
        f"{'lam/Q':>7} {'r(1e-3)':>9} {'r(1e-6)':>11} {'Time':>7}"
    )
    sep = "-" * len(header)

    rows = []
    for r in results:
        qwc_str = f"{r['n_qwc_groups']:>8,}" if r['n_qwc_groups'] is not None else f"{'--':>8}"
        gpt_str = f"{r['groups_per_term']:>8.3f}" if r['groups_per_term'] is not None else f"{'--':>8}"
        rows.append(
            f"{r['label']:<22} "
            f"{r['Q']:>4} "
            f"{r['n_pauli_terms']:>9,} "
            f"{r['eri_density']:>6.1%} "
            f"{qwc_str} "
            f"{gpt_str} "
            f"{r['one_norm']:>9.2f} "
            f"{r['lambda_per_Q']:>7.3f} "
            f"{r['trotter_steps_eps3']:>9,} "
            f"{r['trotter_steps_eps6']:>11,} "
            f"{r['time_total_s']:>6.1f}s"
        )

    return "\n".join(["", "MASTER TABLE: Quantum Simulation Cost Metrics",
                       "=" * len(header), header, sep] + rows + [sep, ""])


def format_scaling_table(
    fits: Dict[str, Dict[str, float]],
    series_label: str,
) -> str:
    """Format the scaling exponents table."""
    nice_names = {
        'pauli_terms': 'Pauli terms',
        'qwc_groups': 'QWC groups',
        'one_norm': '1-norm (lambda)',
        'trotter_eps3': 'Trotter steps (1e-3)',
        'trotter_eps6': 'Trotter steps (1e-6)',
    }

    header = f"{'Metric':<25} {'Exponent':>10} {'R^2':>8}"
    sep = "-" * 45

    rows = []
    for key, data in fits.items():
        name = nice_names.get(key, key)
        rows.append(
            f"{name:<25} {data['exponent']:>10.3f} {data['r_squared']:>8.4f}"
        )

    return "\n".join([
        "",
        f"SCALING EXPONENTS: {series_label} (metric ~ Q^alpha)",
        "=" * 45, header, sep,
    ] + rows + [sep, ""])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the full scaling sweep."""
    print("=" * 70)
    print("  GeoVac Qubit Scaling Sweep")
    print("  Comprehensive quantum simulation cost metrics")
    print("=" * 70)

    all_results: List[Dict[str, Any]] = []

    # Part 1: GeoVac He (2 electrons)
    print("\nPart 1: GeoVac He (Z=2, 2 electrons)")
    print("-" * 40)
    he_results = []
    for nmax in [2, 3, 4, 5]:
        result = collect_geovac_system(
            label=f"GeoVac He n={nmax}",
            Z=2, n_electrons=2, max_n=nmax,
        )
        if result is not None:
            he_results.append(result)
            all_results.append(result)

    # Part 2: GeoVac H (1 electron, control)
    print("\nPart 2: GeoVac H (Z=1, 1 electron)")
    print("-" * 40)
    h_results = []
    for nmax in [2, 3, 4, 5]:
        result = collect_geovac_system(
            label=f"GeoVac H n={nmax}",
            Z=1, n_electrons=1, max_n=nmax,
        )
        if result is not None:
            h_results.append(result)
            all_results.append(result)

    # Part 3: Gaussian baselines
    print("\nPart 3: Gaussian baselines")
    print("-" * 40)
    gaussian_results = []
    for label, fn in [("Gauss He STO-3G", he_sto3g), ("Gauss H2 STO-3G", h2_sto3g)]:
        result = collect_gaussian_system(label, fn)
        gaussian_results.append(result)
        all_results.append(result)

    # Part 4: Scaling law fits
    print("\nPart 4: Scaling law fits")
    print("-" * 40)

    he_fits = {}
    if len(he_results) >= 3:
        he_fits = fit_scaling_laws(he_results)
        print(f"  GeoVac He: {len(he_results)} data points")
        for name, data in he_fits.items():
            print(f"    {name}: exponent={data['exponent']:.3f}, R²={data['r_squared']:.4f}")
    else:
        print("  GeoVac He: insufficient data points for fit")

    h_fits = {}
    if len(h_results) >= 3:
        h_fits = fit_scaling_laws(h_results)
        print(f"  GeoVac H:  {len(h_results)} data points")
        for name, data in h_fits.items():
            print(f"    {name}: exponent={data['exponent']:.3f}, R²={data['r_squared']:.4f}")
    else:
        print("  GeoVac H: insufficient data points for fit")

    # Part 5: Output
    print("\n" + "=" * 70)
    print("  RESULTS")
    print("=" * 70)

    master = format_master_table(all_results)
    print(master)

    scaling_tables = []
    if he_fits:
        t = format_scaling_table(he_fits, "GeoVac He (2e)")
        print(t)
        scaling_tables.append(t)
    if h_fits:
        t = format_scaling_table(h_fits, "GeoVac H (1e)")
        print(t)
        scaling_tables.append(t)

    # Save JSON
    out_dir = Path(__file__).resolve().parent
    json_path = out_dir / "qubit_scaling_data.json"

    json_data = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'systems': all_results,
        'scaling_fits': {
            'geovac_he': he_fits,
            'geovac_h': h_fits,
        },
    }

    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"\nJSON saved: {json_path}")

    # Save summary text
    txt_path = out_dir / "qubit_scaling_summary.txt"
    with open(txt_path, 'w') as f:
        f.write("GeoVac Qubit Scaling Sweep\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(master)
        for t in scaling_tables:
            f.write(t)
        f.write("\n")
    print(f"Summary saved: {txt_path}")


if __name__ == '__main__':
    main()
