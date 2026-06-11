"""
BeH2 Z_eff partitioning scan.

Runs PES scans for each Z_eff partition scheme and compares results.
Saves comparison data, best-scheme PES, and overlay plot.

Usage:
    python debug/beh2_partition_scan.py
"""

import json
import time
import numpy as np
from pathlib import Path

from geovac.composed_triatomic import ComposedTriatomicSolver

# Experimental reference
REF = {
    'R_eq': 2.507,       # bohr
    'omega_e': 2345.0,   # cm-1
    'D_e': 0.147,        # Ha
}

SCHEMES = ['full', 'equal', 'valence_scaled', 'sqrt']
R_GRID = np.arange(1.5, 5.1, 0.1)
L_MAX = 2
N_RE = 300

OUTPUT_DIR = Path(__file__).parent / 'data'
PLOT_DIR = Path(__file__).parent / 'plots'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def pct_error(computed: float, reference: float) -> float:
    return abs(computed - reference) / abs(reference) * 100


def run_scheme(scheme: str) -> dict:
    """Run full PES scan for one partition scheme."""
    print(f"\n{'='*64}")
    print(f"  SCHEME: '{scheme}'")
    print(f"{'='*64}")

    solver = ComposedTriatomicSolver.BeH2(
        l_max=L_MAX, z_eff_partition=scheme, verbose=True)

    t0 = time.time()
    solver.solve_core()
    solver.scan_pes(R_grid=R_GRID, n_Re=N_RE)
    solver.fit_spectroscopic_constants()
    wall = time.time() - t0

    s = solver.spectro
    pes = solver.pes_result

    # Check PES shape
    E_valid = np.array(pes['E_valid'])
    i_min = np.argmin(E_valid)
    is_bound = pes['D_e'] > 0
    has_wall = i_min > 0
    approaches_dissoc = E_valid[-1] > E_valid[i_min]

    result = {
        'scheme': scheme,
        'Z_eff_per_bond': solver.Z_eff,
        'Z_eff_full': solver.Z_eff_full,
        'pk_A': solver.pk_A,
        'pk_B': solver.pk_B,
        'pk_scaled': solver._scale_pk,
        'R_eq': s['R_eq'],
        'D_e': s['D_e'],
        'omega_e': s['omega_e_sym'],
        'E_min': s['E_min'],
        'R_eq_error_pct': pct_error(s['R_eq'], REF['R_eq']),
        'D_e_error_pct': pct_error(s['D_e'], REF['D_e']),
        'omega_e_error_pct': pct_error(s['omega_e_sym'], REF['omega_e']),
        'is_bound': is_bound,
        'has_repulsive_wall': has_wall,
        'approaches_dissoc': approaches_dissoc,
        'wall_time_s': wall,
        'pes': {
            'R': pes['R'],
            'E_total': pes['E_total'],
            'R_valid': pes['R_valid'],
            'E_valid': pes['E_valid'],
        },
    }
    return result


def print_comparison_table(results: list) -> None:
    """Print side-by-side comparison of all schemes."""
    print(f"\n{'='*80}")
    print("BeH2 Z_eff Partition Comparison (l_max=2)")
    print(f"{'='*80}")
    print(f"Reference: R_eq={REF['R_eq']} bohr, D_e={REF['D_e']} Ha,"
          f" omega_e={REF['omega_e']} cm-1")
    print(f"{'='*80}")

    hdr = (f"{'Scheme':<18s} {'Z_eff':>6s} {'R_eq':>7s} {'err%':>6s}"
           f" {'D_e':>8s} {'err%':>8s} {'omega_e':>8s} {'err%':>6s}"
           f" {'Bound':>5s} {'Wall':>5s} {'Diss':>5s}")
    print(hdr)
    print("-" * 80)

    for r in results:
        row = (f"{r['scheme']:<18s} {r['Z_eff_per_bond']:6.3f}"
               f" {r['R_eq']:7.3f} {r['R_eq_error_pct']:6.1f}"
               f" {r['D_e']:8.4f} {r['D_e_error_pct']:8.1f}"
               f" {r['omega_e']:8.1f} {r['omega_e_error_pct']:6.1f}"
               f" {'Y' if r['is_bound'] else 'N':>5s}"
               f" {'Y' if r['has_repulsive_wall'] else 'N':>5s}"
               f" {'Y' if r['approaches_dissoc'] else 'N':>5s}")
        print(row)

    print(f"{'='*80}")

    # Identify best R_eq (among schemes with valid PES)
    valid = [r for r in results if r['is_bound'] and r['has_repulsive_wall']]
    if valid:
        best = min(valid, key=lambda r: r['R_eq_error_pct'])
        print(f"\nBest R_eq: '{best['scheme']}' ({best['R_eq']:.3f} bohr,"
              f" {best['R_eq_error_pct']:.1f}% error)")


def make_overlay_plot(results: list) -> None:
    """Create PES overlay plot with all schemes."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    colors = {'full': 'blue', 'equal': 'red',
              'valence_scaled': 'green', 'sqrt': 'purple'}

    for r in results:
        R = np.array(r['pes']['R_valid'])
        E = np.array(r['pes']['E_valid'])
        c = colors.get(r['scheme'], 'black')
        lbl = (f"{r['scheme']} (Z_eff={r['Z_eff_per_bond']:.2f},"
               f" R_eq={r['R_eq']:.3f})")

        ax1.plot(R, E, '-o', color=c, markersize=3, label=lbl)

        # Relative to each scheme's dissociation
        E_rel = E - E[-1]
        ax2.plot(R, E_rel, '-o', color=c, markersize=3, label=lbl)

    # Experimental R_eq
    ax1.axvline(x=REF['R_eq'], color='orange', linestyle=':', linewidth=2,
                label=f'R_eq (expt) = {REF["R_eq"]}')
    ax2.axvline(x=REF['R_eq'], color='orange', linestyle=':', linewidth=2,
                label=f'R_eq (expt) = {REF["R_eq"]}')

    ax1.set_xlabel('R (bohr)')
    ax1.set_ylabel('E (Ha)')
    ax1.set_title('BeH2 PES: Z_eff Partition Comparison')
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('R (bohr)')
    ax2.set_ylabel('E - E_dissoc (Ha)')
    ax2.set_title('Binding Curves (relative to dissociation)')
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plot_path = PLOT_DIR / 'beh2_partition_comparison.png'
    plt.savefig(str(plot_path), dpi=150)
    plt.close(fig)
    print(f"\nPlot saved: {plot_path}")


def main() -> None:
    t_start = time.time()

    results = []
    for scheme in SCHEMES:
        result = run_scheme(scheme)
        results.append(result)

    print_comparison_table(results)

    # Save comparison table
    comparison_path = OUTPUT_DIR / 'beh2_zeff_partition_comparison.json'
    # Strip full PES arrays for the comparison file (keep it compact)
    comparison_data = []
    for r in results:
        entry = {k: v for k, v in r.items() if k != 'pes'}
        comparison_data.append(entry)
    with open(comparison_path, 'w') as f:
        json.dump(comparison_data, f, indent=2)
    print(f"Comparison saved: {comparison_path}")

    # Identify best scheme (lowest R_eq error with valid PES)
    valid = [r for r in results if r['is_bound'] and r['has_repulsive_wall']]
    if valid:
        best = min(valid, key=lambda r: r['R_eq_error_pct'])
        best_scheme = best['scheme']
        print(f"\nBest scheme: '{best_scheme}'")

        # Save best-scheme PES
        pes_path = OUTPUT_DIR / 'beh2_partitioned_pes.json'
        pes_data = {
            'scheme': best_scheme,
            'Z_eff_per_bond': best['Z_eff_per_bond'],
            'pk_A': best['pk_A'],
            'pk_B': best['pk_B'],
            'R': best['pes']['R'],
            'E_total': best['pes']['E_total'],
        }
        with open(pes_path, 'w') as f:
            json.dump(pes_data, f, indent=2)
        print(f"Best PES saved: {pes_path}")

        # Save best-scheme spectroscopic constants
        spectro_path = OUTPUT_DIR / 'beh2_partitioned_spectro.json'
        spectro_data = {
            'scheme': best_scheme,
            'Z_eff_per_bond': best['Z_eff_per_bond'],
            'computed': {
                'R_eq': best['R_eq'],
                'D_e': best['D_e'],
                'omega_e': best['omega_e'],
            },
            'reference': REF,
            'errors': {
                'R_eq_pct': best['R_eq_error_pct'],
                'D_e_pct': best['D_e_error_pct'],
                'omega_e_pct': best['omega_e_error_pct'],
            },
        }
        with open(spectro_path, 'w') as f:
            json.dump(spectro_data, f, indent=2)
        print(f"Best spectro saved: {spectro_path}")

    # Make overlay plot
    make_overlay_plot(results)

    print(f"\nTotal wall time: {time.time() - t_start:.0f}s")


if __name__ == '__main__':
    main()
