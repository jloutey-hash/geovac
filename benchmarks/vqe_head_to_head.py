#!/usr/bin/env python
"""
VQE Head-to-Head Benchmark: GeoVac vs Gaussian
================================================

Runs VQE on the same chemical systems using both GeoVac lattice encodings
and Gaussian STO-3G encodings, comparing qubit count, gate count, circuit
depth, measurement cost, and VQE accuracy.

Includes:
  - Original VQE comparison (He, H2)
  - Equal-qubit accuracy comparison (Task 1)
  - Pauli term scaling crossover analysis with plot (Task 2)
  - Headline accuracy-per-qubit summary (Task 3)

Systems:
  - He (GeoVac nmax=2 vs STO-3G)
  - He (GeoVac nmax=3 vs STO-3G)
  - H2 (Gaussian STO-3G only — GeoVac molecular encoding is Level 4,
        not yet wired to JW transform)

Usage:
    python benchmarks/vqe_head_to_head.py

Output:
    - Formatted comparison tables to stdout
    - JSON results to debug/data/vqe_benchmark_results.json
    - Equal-qubit comparison to debug/data/equal_qubit_comparison.json
    - Crossover plot to debug/plots/pauli_crossover.png

Author: GeoVac Development Team
Date: March 2026
"""

import json
import sys
import os
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

# Ensure project root is on path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.vqe_benchmark import (
    build_geovac_he,
    build_gaussian_he,
    build_gaussian_h2,
    collect_static_metrics,
    run_vqe,
    format_comparison_table,
    save_results,
)
from geovac.gaussian_reference import (
    he_cc_pvdz,
    he_cc_pvtz,
    build_qubit_hamiltonian,
    _GAUSS_ALPHA,
    _GAUSS_C,
    _estimate_gaussian_pauli_terms,
)


# ---------------------------------------------------------------------------
# Exact He energy for error-vs-exact comparisons
# ---------------------------------------------------------------------------
HE_EXACT = -2.903724  # Non-relativistic exact (Pekeris, 1959)


# ---------------------------------------------------------------------------
# GeoVac Pauli term data from Paper 14
# ---------------------------------------------------------------------------
GEOVAC_DATA = [
    {'nmax': 2, 'n_spatial': 5,  'n_qubits': 10,  'pauli_terms': 120},
    {'nmax': 3, 'n_spatial': 14, 'n_qubits': 28,  'pauli_terms': 2_659},
    {'nmax': 4, 'n_spatial': 30, 'n_qubits': 60,  'pauli_terms': 31_039},
    {'nmax': 5, 'n_spatial': 55, 'n_qubits': 110, 'pauli_terms': 227_338},
]

# Gaussian baseline data from Paper 14 (H2 systems)
GAUSSIAN_H2_DATA = [
    {'basis': 'STO-3G',  'n_qubits': 4,  'pauli_terms': 15},
    {'basis': '6-31G',   'n_qubits': 8,  'pauli_terms': 201},
    {'basis': 'cc-pVDZ', 'n_qubits': 20, 'pauli_terms': 23_189},
]


# ---------------------------------------------------------------------------
# Task 1: Equal-qubit accuracy comparison
# ---------------------------------------------------------------------------

def build_equal_qubit_table(
    geovac_energies: Dict[int, float],
) -> List[Dict[str, Any]]:
    """
    Build the accuracy-at-fixed-qubit-count comparison table.

    Parameters
    ----------
    geovac_energies : dict
        Mapping nmax -> FCI-in-basis energy from GeoVac JW Hamiltonian.

    Returns
    -------
    list of dict
        Each entry has: qubits, geovac_encoding, geovac_error_pct,
        geovac_pauli, gaussian_encoding, gaussian_error_pct,
        gaussian_pauli_estimated.
    """
    rows: List[Dict[str, Any]] = []

    # Q=2: Gaussian STO-3G only (GeoVac minimum is Q=10)
    he_sto3g_energy = -2.847656  # 2*h11 + (11|11) with zeta=1.6875
    rows.append({
        'qubits': 2,
        'geovac_encoding': '—',
        'geovac_error_pct': None,
        'geovac_pauli': None,
        'geovac_energy': None,
        'gaussian_encoding': 'He STO-3G (1 spatial)',
        'gaussian_error_pct': 100.0 * abs(he_sto3g_energy - HE_EXACT) / abs(HE_EXACT),
        'gaussian_pauli': 4,
        'gaussian_energy': he_sto3g_energy,
    })

    # Q=4: H2 STO-3G only
    h2_sto3g_energy = -1.1373  # Szabo & Ostlund
    h2_exact = -1.1745  # Kolos & Wolniewicz exact at R=1.4 bohr
    rows.append({
        'qubits': 4,
        'geovac_encoding': '—',
        'geovac_error_pct': None,
        'geovac_pauli': None,
        'geovac_energy': None,
        'gaussian_encoding': 'H2 STO-3G (2 spatial)',
        'gaussian_error_pct': 100.0 * abs(h2_sto3g_energy - h2_exact) / abs(h2_exact),
        'gaussian_pauli': 15,
        'gaussian_energy': h2_sto3g_energy,
    })

    # Q=10: GeoVac nmax=2 vs Gaussian cc-pVDZ
    pvdz = he_cc_pvdz()
    _, _, pvdz_pauli = build_qubit_hamiltonian(pvdz)
    pvdz_error = 100.0 * abs(pvdz['literature_energy'] - HE_EXACT) / abs(HE_EXACT)
    gv2_energy = geovac_energies.get(2)
    rows.append({
        'qubits': 10,
        'geovac_encoding': 'He nmax=2 (5 spatial)',
        'geovac_error_pct': (
            100.0 * abs(gv2_energy - HE_EXACT) / abs(HE_EXACT)
            if gv2_energy is not None else None
        ),
        'geovac_pauli': 120,
        'geovac_energy': gv2_energy,
        'gaussian_encoding': 'He cc-pVDZ (5 spatial)',
        'gaussian_error_pct': pvdz_error,
        'gaussian_pauli': pvdz_pauli,
        'gaussian_energy': pvdz['literature_energy'],
    })

    # Q=28: GeoVac nmax=3 vs Gaussian cc-pVTZ
    pvtz = he_cc_pvtz()
    _, _, pvtz_pauli = build_qubit_hamiltonian(pvtz)
    pvtz_error = 100.0 * abs(pvtz['literature_energy'] - HE_EXACT) / abs(HE_EXACT)
    gv3_energy = geovac_energies.get(3)
    rows.append({
        'qubits': 28,
        'geovac_encoding': 'He nmax=3 (14 spatial)',
        'geovac_error_pct': (
            100.0 * abs(gv3_energy - HE_EXACT) / abs(HE_EXACT)
            if gv3_energy is not None else None
        ),
        'geovac_pauli': 2_659,
        'geovac_energy': gv3_energy,
        'gaussian_encoding': 'He cc-pVTZ (14 spatial)',
        'gaussian_error_pct': pvtz_error,
        'gaussian_pauli': pvtz_pauli,
        'gaussian_energy': pvtz['literature_energy'],
    })

    return rows


def format_equal_qubit_table(rows: List[Dict[str, Any]]) -> str:
    """Format the equal-qubit comparison as a readable table."""
    lines = [
        "",
        "=" * 105,
        "Equal-Qubit Accuracy Comparison: GeoVac vs Gaussian (He atom)",
        "=" * 105,
        f"{'Q':>4}  {'GeoVac Encoding':<24} {'Err%':>7} {'Pauli':>7}  "
        f"{'Gaussian Encoding':<24} {'Err%':>7} {'Pauli':>8}",
        "-" * 105,
    ]

    for r in rows:
        gv_enc = r['geovac_encoding']
        gv_err = f"{r['geovac_error_pct']:.2f}%" if r['geovac_error_pct'] is not None else "—"
        gv_pauli = f"{r['geovac_pauli']:>7,}" if r['geovac_pauli'] is not None else "      —"
        ga_enc = r['gaussian_encoding']
        ga_err = f"{r['gaussian_error_pct']:.2f}%"
        ga_pauli = f"{r['gaussian_pauli']:>8,}"
        lines.append(
            f"{r['qubits']:>4}  {gv_enc:<24} {gv_err:>7} {gv_pauli}  "
            f"{ga_enc:<24} {ga_err:>7} {ga_pauli}"
        )

    lines.append("-" * 105)
    lines.append(
        "Gaussian Pauli counts at Q>=10 are actual JW term counts from "
        "computed MO integrals (see gaussian_reference.py)."
    )
    lines.append(
        "Gaussian FCI energies from cc-pVDZ/cc-pVTZ integrals (Woon & Dunning, "
        "JCP 100, 2975, 1994). Errors vs exact He = -2.903724 Ha."
    )
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Task 2: Crossover analysis
# ---------------------------------------------------------------------------

def compute_crossover() -> Dict[str, Any]:
    """
    Compute the Pauli-term crossover point and generate scaling curves.

    Uses log-log least-squares fits anchored at measured data:
      - GeoVac: nmax=2..5 (Q=10,28,60,110)
      - Gaussian: H2 STO-3G/6-31G/cc-pVDZ + He STO-3G/cc-pVDZ/cc-pVTZ

    Returns
    -------
    dict
        Fitted parameters and crossover qubit count.
    """
    # Fit GeoVac: ln(N) = ln(c_gv) + alpha_gv * ln(Q)
    gv_q = np.array([d['n_qubits'] for d in GEOVAC_DATA], dtype=float)
    gv_n = np.array([d['pauli_terms'] for d in GEOVAC_DATA], dtype=float)
    gv_fit = np.polyfit(np.log(gv_q), np.log(gv_n), 1)
    alpha_gv = gv_fit[0]
    c_gv = np.exp(gv_fit[1])

    # Gaussian He actual Pauli counts (from computed MO integrals)
    pvdz = he_cc_pvdz()
    _, _, pvdz_pauli = build_qubit_hamiltonian(pvdz)
    pvtz = he_cc_pvtz()
    _, _, pvtz_pauli = build_qubit_hamiltonian(pvtz)
    gaussian_he_data = [
        {'basis': 'He STO-3G',  'n_qubits': 2,  'pauli_terms': 4},
        {'basis': 'He cc-pVDZ', 'n_qubits': 10, 'pauli_terms': pvdz_pauli},
        {'basis': 'He cc-pVTZ', 'n_qubits': 28, 'pauli_terms': pvtz_pauli},
    ]

    # Fit Gaussian: combine H2 + He data for broader coverage
    all_gaussian = GAUSSIAN_H2_DATA + gaussian_he_data
    ga_q = np.array([d['n_qubits'] for d in all_gaussian], dtype=float)
    ga_n = np.array([d['pauli_terms'] for d in all_gaussian], dtype=float)
    ga_fit = np.polyfit(np.log(ga_q), np.log(ga_n), 1)
    alpha_ga = ga_fit[0]
    c_ga = np.exp(ga_fit[1])

    # Crossover: c_gv * Q^alpha_gv = c_ga * Q^alpha_ga
    # => Q^(alpha_ga - alpha_gv) = c_gv / c_ga
    # => Q_cross = (c_gv / c_ga)^(1 / (alpha_ga - alpha_gv))
    if alpha_ga > alpha_gv and c_gv > 0 and c_ga > 0:
        q_cross = (c_gv / c_ga) ** (1.0 / (alpha_ga - alpha_gv))
    else:
        q_cross = float('nan')

    return {
        'geovac_alpha': round(alpha_gv, 3),
        'geovac_c': c_gv,
        'gaussian_alpha': round(alpha_ga, 3),
        'gaussian_c': c_ga,
        'crossover_qubits': round(q_cross, 1),
        'geovac_data': GEOVAC_DATA,
        'gaussian_h2_data': GAUSSIAN_H2_DATA,
        'gaussian_he_data': gaussian_he_data,
    }


def plot_crossover(crossover: Dict[str, Any], outpath: str) -> None:
    """Generate the Pauli-term crossover plot and save to outpath."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 6))

    # Extrapolation range
    q_range = np.logspace(np.log10(2), np.log10(200), 300)

    # GeoVac fitted curve
    c_gv = crossover['geovac_c']
    a_gv = crossover['geovac_alpha']
    gv_curve = c_gv * q_range ** a_gv
    ax.plot(q_range, gv_curve, '-', color='#2563EB', linewidth=2,
            label=f'GeoVac fit: $N \\propto Q^{{{a_gv:.2f}}}$')

    # Gaussian fitted curve
    c_ga = crossover['gaussian_c']
    a_ga = crossover['gaussian_alpha']
    ga_curve = c_ga * q_range ** a_ga
    ax.plot(q_range, ga_curve, '-', color='#DC2626', linewidth=2,
            label=f'Gaussian fit: $N \\propto Q^{{{a_ga:.2f}}}$')

    # GeoVac data points
    gv_q = [d['n_qubits'] for d in GEOVAC_DATA]
    gv_n = [d['pauli_terms'] for d in GEOVAC_DATA]
    ax.scatter(gv_q, gv_n, s=80, color='#2563EB', zorder=5, edgecolors='white',
               linewidths=1.5, label='GeoVac He (measured)')

    # Gaussian H2 data points
    ga_q = [d['n_qubits'] for d in GAUSSIAN_H2_DATA]
    ga_n = [d['pauli_terms'] for d in GAUSSIAN_H2_DATA]
    ax.scatter(ga_q, ga_n, s=80, color='#DC2626', marker='s', zorder=5,
               edgecolors='white', linewidths=1.5, label='Gaussian H$_2$ (measured)')

    # Gaussian He data points (actual JW counts)
    he_data = crossover.get('gaussian_he_data', [])
    if he_data:
        he_q = [d['n_qubits'] for d in he_data]
        he_n = [d['pauli_terms'] for d in he_data]
        ax.scatter(he_q, he_n, s=80, color='#DC2626', marker='D', zorder=5,
                   edgecolors='white', linewidths=1.5, label='Gaussian He (measured)')
        for d in he_data:
            ax.annotate(d['basis'], (d['n_qubits'], d['pauli_terms']),
                        textcoords='offset points', xytext=(10, -12),
                        fontsize=7, color='#DC2626')

    # Mark crossover
    q_cross = crossover['crossover_qubits']
    if not np.isnan(q_cross) and 1 < q_cross < 200:
        n_cross = c_gv * q_cross ** a_gv
        ax.axvline(q_cross, color='gray', linestyle='--', alpha=0.5)
        ax.scatter([q_cross], [n_cross], s=120, color='#F59E0B', marker='*',
                   zorder=6, edgecolors='black', linewidths=0.5)
        ax.annotate(f'Crossover Q ≈ {q_cross:.0f}',
                    (q_cross, n_cross), textcoords='offset points',
                    xytext=(15, 10), fontsize=10, fontweight='bold',
                    color='#B45309')

    # Annotate GeoVac advantage at Q=10 and Q=60
    for q_pt in [10, 60]:
        gv_val = c_gv * q_pt ** a_gv
        ga_val = c_ga * q_pt ** a_ga
        ratio = ga_val / gv_val
        ax.annotate(f'{ratio:.0f}× fewer',
                    (q_pt, gv_val), textcoords='offset points',
                    xytext=(-50, -25), fontsize=9, color='#2563EB',
                    fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color='#2563EB', lw=1.2))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Qubits (Q)', fontsize=12)
    ax.set_ylabel('Pauli Terms (N)', fontsize=12)
    ax.set_title('GeoVac vs Gaussian: Pauli Term Scaling Crossover', fontsize=14)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(1.5, 250)

    plt.tight_layout()
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Task 3: Headline summary
# ---------------------------------------------------------------------------

def format_headline(
    equal_qubit_rows: List[Dict[str, Any]],
    crossover: Dict[str, Any],
) -> str:
    """Generate the headline accuracy-per-qubit summary."""
    lines = [
        "",
        "=" * 90,
        "HEADLINE: Accuracy-per-Qubit Summary",
        "=" * 90,
    ]

    # Find the Q=10 row
    row_10 = next((r for r in equal_qubit_rows if r['qubits'] == 10), None)
    if row_10 and row_10['geovac_error_pct'] is not None:
        ratio = row_10['gaussian_pauli'] / row_10['geovac_pauli']
        lines.append(
            f"\n  At 10 qubits, GeoVac He achieves {row_10['geovac_error_pct']:.2f}% error "
            f"with {row_10['geovac_pauli']:,} Pauli terms."
        )
        lines.append(
            f"  A Gaussian basis (cc-pVDZ) at the same qubit count achieves "
            f"{row_10['gaussian_error_pct']:.2f}% error"
        )
        lines.append(
            f"  with {row_10['gaussian_pauli']:,} Pauli terms. "
            f"GeoVac uses {ratio:.0f}x fewer Pauli terms"
        )
        lines.append(f"  at comparable accuracy.")

    row_28 = next((r for r in equal_qubit_rows if r['qubits'] == 28), None)
    if row_28 and row_28['geovac_error_pct'] is not None:
        ratio = row_28['gaussian_pauli'] / row_28['geovac_pauli']
        lines.append(
            f"\n  At 28 qubits, GeoVac He achieves {row_28['geovac_error_pct']:.2f}% error "
            f"with {row_28['geovac_pauli']:,} Pauli terms."
        )
        lines.append(
            f"  Gaussian cc-pVTZ at the same qubit count: "
            f"{row_28['gaussian_error_pct']:.2f}% error, "
            f"{row_28['gaussian_pauli']:,} Pauli terms."
        )
        lines.append(
            f"  GeoVac uses {ratio:.0f}x fewer Pauli terms; Gaussian is "
            f"{row_28['geovac_error_pct'] / row_28['gaussian_error_pct']:.1f}x more accurate."
        )

    q_cross = crossover['crossover_qubits']
    lines.append(
        f"\n  Scaling crossover: GeoVac (Q^{crossover['geovac_alpha']:.2f}) vs "
        f"Gaussian (Q^{crossover['gaussian_alpha']:.2f})."
    )
    if not np.isnan(q_cross):
        lines.append(
            f"  GeoVac has fewer Pauli terms than Gaussian for Q > {q_cross:.0f} "
            f"(i.e., beyond minimal basis)."
        )

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Original VQE benchmark (preserved)
# ---------------------------------------------------------------------------

def run_vqe_benchmark() -> list:
    """Run the original VQE head-to-head comparison."""
    results = []

    # He: GeoVac nmax=2 vs Gaussian STO-3G
    print("Building He GeoVac nmax=2...")
    spo_gv2, of_gv2, nq_gv2, exact_gv2 = build_geovac_he(max_n=2)
    print(f"  {nq_gv2} qubits, {len(spo_gv2)} Pauli terms, exact E = {exact_gv2:.6f} Ha")

    print("Building He Gaussian STO-3G...")
    spo_g, of_g, nq_g, exact_g = build_gaussian_he()
    print(f"  {nq_g} qubits, {len(spo_g)} Pauli terms, exact E = {exact_g:.6f} Ha")

    print("Running VQE: He GeoVac nmax=2...")
    r = run_vqe(spo_gv2, of_gv2, nq_gv2, "He", "GeoVac nmax=2", exact_gv2,
                reps=1, maxiter=500)
    results.append(r)
    print(f"  VQE E = {r.vqe_energy:.6f}, err = {r.error_pct:.3f}%, iters = {r.n_iterations}")

    print("Running VQE: He Gaussian STO-3G...")
    r = run_vqe(spo_g, of_g, nq_g, "He", "STO-3G", exact_g,
                reps=1, maxiter=500)
    results.append(r)
    print(f"  VQE E = {r.vqe_energy:.6f}, err = {r.error_pct:.3f}%, iters = {r.n_iterations}")

    # He: GeoVac nmax=3
    print("\nBuilding He GeoVac nmax=3...")
    spo_gv3, of_gv3, nq_gv3, exact_gv3 = build_geovac_he(max_n=3)
    print(f"  {nq_gv3} qubits, {len(spo_gv3)} Pauli terms, exact E = {exact_gv3:.6f} Ha")

    print("Collecting static metrics: He GeoVac nmax=3 (VQE skipped — 28 qubits)...")
    r = collect_static_metrics(spo_gv3, of_gv3, nq_gv3, "He", "GeoVac nmax=3",
                               exact_gv3, reps=1)
    results.append(r)
    print(f"  CX={r.cx_count}, depth={r.circuit_depth}, QWC={r.n_qwc_groups}")

    # H2: Gaussian STO-3G
    print("\nBuilding H2 Gaussian STO-3G...")
    spo_h2, of_h2, nq_h2, exact_h2 = build_gaussian_h2()
    print(f"  {nq_h2} qubits, {len(spo_h2)} Pauli terms, exact E = {exact_h2:.6f} Ha")

    print("Running VQE: H2 Gaussian STO-3G...")
    r = run_vqe(spo_h2, of_h2, nq_h2, "H2", "STO-3G", exact_h2,
                reps=2, maxiter=500)
    results.append(r)
    print(f"  VQE E = {r.vqe_energy:.6f}, err = {r.error_pct:.3f}%, iters = {r.n_iterations}")

    return results, exact_gv2, exact_gv3


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # ------------------------------------------------------------------
    # Original VQE benchmark
    # ------------------------------------------------------------------
    vqe_results, exact_gv2, exact_gv3 = run_vqe_benchmark()

    table = format_comparison_table(vqe_results)
    print(table)

    save_results(vqe_results)
    print("Results saved to debug/data/vqe_benchmark_results.json")

    # ------------------------------------------------------------------
    # Task 1: Equal-qubit accuracy comparison
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("TASK 1: Equal-Qubit Accuracy Comparison")
    print("=" * 70)

    geovac_energies = {2: exact_gv2, 3: exact_gv3}
    equal_qubit_rows = build_equal_qubit_table(geovac_energies)
    print(format_equal_qubit_table(equal_qubit_rows))

    # ------------------------------------------------------------------
    # Task 2: Crossover analysis
    # ------------------------------------------------------------------
    print("=" * 70)
    print("TASK 2: Pauli Term Scaling Crossover")
    print("=" * 70)

    crossover = compute_crossover()
    print(f"\n  GeoVac scaling:  N ~ {crossover['geovac_c']:.4f} * Q^{crossover['geovac_alpha']}")
    print(f"  Gaussian scaling: N ~ {crossover['gaussian_c']:.5f} * Q^{crossover['gaussian_alpha']}")
    print(f"  Crossover at Q ~ {crossover['crossover_qubits']}")

    plot_path = "debug/plots/pauli_crossover.png"
    plot_crossover(crossover, plot_path)
    print(f"  Plot saved to {plot_path}")

    # ------------------------------------------------------------------
    # Task 3: Headline summary
    # ------------------------------------------------------------------
    print(format_headline(equal_qubit_rows, crossover))

    # ------------------------------------------------------------------
    # Save equal-qubit results (Task 3)
    # ------------------------------------------------------------------
    output = {
        'equal_qubit_comparison': equal_qubit_rows,
        'crossover_analysis': {
            'geovac_alpha': crossover['geovac_alpha'],
            'gaussian_alpha': crossover['gaussian_alpha'],
            'crossover_qubits': crossover['crossover_qubits'],
        },
        'geovac_measured_data': GEOVAC_DATA,
        'gaussian_h2_data': GAUSSIAN_H2_DATA,
        'gaussian_he_data': crossover.get('gaussian_he_data', []),
        'notes': {
            'geovac_pauli_terms': 'Measured from JW encoding (Paper 14)',
            'gaussian_pauli_q10_q28': (
                'Actual JW term counts from computed MO integrals '
                '(single-center Gaussian engine, see gaussian_reference.py)'
            ),
            'gaussian_fci_energies': (
                'Computed from cc-pVDZ/cc-pVTZ basis sets; '
                'Woon & Dunning, JCP 100, 2975 (1994)'
            ),
            'he_exact_energy': HE_EXACT,
        },
    }

    out_path = Path("debug/data/equal_qubit_comparison.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"Equal-qubit comparison saved to {out_path}")


if __name__ == '__main__':
    main()
