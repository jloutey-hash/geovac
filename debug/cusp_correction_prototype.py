#!/usr/bin/env python3
"""
Cusp Correction Prototype — DBBSC for GeoVac
=============================================

Runs the density-based basis-set correction on He and He-like ions.

Parts:
1. He at n_max = 2..7: uncorrected vs corrected energy
2. Coalescence density scaling
3. He-like ions Z=2..10 at n_max=4
4. Comparison with known 2D variational cusp correction
5. Plots

Author: GeoVac Development Team
Date: April 2026
"""

import json
import os
import sys
import time

import numpy as np

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from geovac.cusp_correction_dbbsc import (
    he_cusp_correction_from_ci,
    coalescence_density_from_ci,
    schwartz_correction,
    hurwitz_correction,
    effective_l_max,
    _radial_at_origin,
)


# Exact energies for He-like ions (NIST / Drake)
EXACT_ENERGIES = {
    1: -0.527751,    # H-
    2: -2.903724,    # He
    3: -7.279913,    # Li+
    4: -13.655566,   # Be2+
    5: -22.030972,   # B3+
    6: -32.406247,   # C4+
    8: -59.156570,   # O6+
    10: -93.906665,  # Ne8+
}


def part1_he_convergence():
    """Part 1: He cusp correction convergence with n_max."""
    print("=" * 70)
    print("PART 1: He (Z=2) Cusp Correction Convergence")
    print("=" * 70)

    Z = 2
    results = []

    for n_max in range(2, 8):
        print(f"\n--- n_max = {n_max} ---")
        t0 = time.time()
        res = he_cusp_correction_from_ci(Z, n_max, method='hurwitz')
        elapsed = time.time() - t0

        E_exact = EXACT_ENERGIES[Z]
        error_uncorrected = abs((res['E_ci'] - E_exact) / E_exact) * 100
        error_schwartz = abs((res['E_ci'] + res['delta_E_schwartz'] - E_exact) / E_exact) * 100
        error_hurwitz = abs((res['E_ci'] + res['delta_E_hurwitz'] - E_exact) / E_exact) * 100
        error_hurwitz_s3 = abs((res['E_ci'] + res['delta_E_hurwitz_s3'] - E_exact) / E_exact) * 100

        entry = {
            'n_max': n_max,
            'l_max': res['l_max'],
            'n_configs': res['n_configs'],
            'E_ci': res['E_ci'],
            'coalescence': res['coalescence'],
            'delta_E_schwartz': res['delta_E_schwartz'],
            'delta_E_hurwitz': res['delta_E_hurwitz'],
            'delta_E_hurwitz_s3': res['delta_E_hurwitz_s3'],
            'E_corrected_schwartz': res['E_ci'] + res['delta_E_schwartz'],
            'E_corrected_hurwitz': res['E_ci'] + res['delta_E_hurwitz'],
            'E_corrected_hurwitz_s3': res['E_ci'] + res['delta_E_hurwitz_s3'],
            'error_uncorrected_pct': error_uncorrected,
            'error_schwartz_pct': error_schwartz,
            'error_hurwitz_pct': error_hurwitz,
            'error_hurwitz_s3_pct': error_hurwitz_s3,
            'elapsed_s': elapsed,
        }
        results.append(entry)

        print(f"  n_configs = {res['n_configs']}, l_max = {res['l_max']}")
        print(f"  E_CI          = {res['E_ci']:.8f} Ha")
        print(f"  <d3(r12)>     = {res['coalescence']:.6f} bohr^-3")
        print(f"  dE_Schwartz   = {res['delta_E_schwartz']:.8f} Ha")
        print(f"  dE_Hurwitz    = {res['delta_E_hurwitz']:.8f} Ha")
        print(f"  dE_Hurwitz_s3 = {res['delta_E_hurwitz_s3']:.8f} Ha")
        print(f"  E_exact       = {E_exact:.8f} Ha")
        print(f"  Error (raw)   = {error_uncorrected:.4f}%")
        print(f"  Error (Schw.) = {error_schwartz:.4f}%")
        print(f"  Error (H,s4)  = {error_hurwitz:.4f}%")
        print(f"  Error (H,s3)  = {error_hurwitz_s3:.4f}%")
        print(f"  Time: {elapsed:.1f}s")

    print("\n\nSummary Table:")
    print(f"{'n_max':>5} {'l_max':>5} {'configs':>7} {'E_CI':>12} {'coal':>10} "
          f"{'err_raw%':>9} {'err_S%':>9} {'err_H4%':>9} {'err_H3%':>9}")
    print("-" * 90)
    for r in results:
        print(f"{r['n_max']:>5d} {r['l_max']:>5d} {r['n_configs']:>7d} "
              f"{r['E_ci']:>12.6f} {r['coalescence']:>10.6f} "
              f"{r['error_uncorrected_pct']:>9.4f} {r['error_schwartz_pct']:>9.4f} "
              f"{r['error_hurwitz_pct']:>9.4f} {r['error_hurwitz_s3_pct']:>9.4f}")

    return results


def part2_z_scaling():
    """Part 2: He-like ions Z=2..10 at n_max=4."""
    print("\n" + "=" * 70)
    print("PART 2: He-like Ions Z=2..10 at n_max=4")
    print("=" * 70)

    n_max = 4
    z_values = [2, 3, 4, 5, 6, 8, 10]
    results = []

    for Z in z_values:
        print(f"\n--- Z = {Z} ---")
        t0 = time.time()
        res = he_cusp_correction_from_ci(Z, n_max, method='hurwitz')
        elapsed = time.time() - t0

        E_exact = EXACT_ENERGIES.get(Z)
        if E_exact is not None:
            error_raw = abs((res['E_ci'] - E_exact) / E_exact) * 100
            error_hurwitz = abs((res['E_ci'] + res['delta_E_hurwitz'] - E_exact) / E_exact) * 100
        else:
            error_raw = None
            error_hurwitz = None

        entry = {
            'Z': Z,
            'n_max': n_max,
            'E_ci': res['E_ci'],
            'coalescence': res['coalescence'],
            'delta_E_hurwitz': res['delta_E_hurwitz'],
            'E_corrected': res['E_ci'] + res['delta_E_hurwitz'],
            'E_exact': E_exact,
            'error_raw_pct': error_raw,
            'error_hurwitz_pct': error_hurwitz,
            'elapsed_s': elapsed,
        }
        results.append(entry)

        print(f"  E_CI       = {res['E_ci']:.8f} Ha")
        print(f"  <d3(r12)>  = {res['coalescence']:.6f} bohr^-3")
        print(f"  dE_Hurwitz = {res['delta_E_hurwitz']:.8f} Ha")
        if E_exact:
            print(f"  Error raw  = {error_raw:.4f}%")
            print(f"  Error corr = {error_hurwitz:.4f}%")

    # Check Z-scaling of coalescence density
    print("\n\nCoalescence density Z-scaling:")
    print(f"{'Z':>4} {'coal':>12} {'coal/Z^3':>12} {'coal/Z^6':>14}")
    print("-" * 50)
    for r in results:
        Z = r['Z']
        coal = r['coalescence']
        print(f"{Z:>4d} {coal:>12.6f} {coal/Z**3:>12.8f} {coal/Z**6:>14.10f}")

    return results


def part3_comparison():
    """Part 3: Compare with known results."""
    print("\n" + "=" * 70)
    print("PART 3: Comparison with Known Results")
    print("=" * 70)

    # Known results from CLAUDE.md:
    # - Adiabatic coupled-channel floor: 0.19-0.20% at n_max=7
    # - 2D variational cusp correction: 0.004% at l_max=4
    # - Graph-native CI: 0.19% at n_max=7

    print("\nKnown GeoVac results:")
    print("  Adiabatic floor:           0.19-0.20%")
    print("  2D var cusp (l_max=4):     0.004%")
    print("  Graph-native CI n_max=7:   0.19%")
    print("  Graph-native CI n_max=9:   0.201%")

    # Our DBBSC correction at n_max=7
    res = he_cusp_correction_from_ci(2, 7, method='hurwitz')
    E_exact = EXACT_ENERGIES[2]
    err_raw = abs((res['E_ci'] - E_exact) / E_exact) * 100
    err_hurwitz = abs((res['E_ci'] + res['delta_E_hurwitz'] - E_exact) / E_exact) * 100

    print(f"\nDBBSC results (n_max=7):")
    print(f"  Raw error:      {err_raw:.4f}%")
    print(f"  Hurwitz corr:   {err_hurwitz:.4f}%")
    print(f"  Improvement:    {err_raw/err_hurwitz:.1f}x" if err_hurwitz > 0 else "")

    return {
        'n_max_7_raw': err_raw,
        'n_max_7_hurwitz': err_hurwitz,
    }


def generate_plots(he_results, z_results):
    """Generate diagnostic plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot_dir = os.path.join(project_root, 'debug', 'plots')

    # Plot 1: Convergence with n_max
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    n_max_vals = [r['n_max'] for r in he_results]
    err_raw = [r['error_uncorrected_pct'] for r in he_results]
    err_schwartz = [r['error_schwartz_pct'] for r in he_results]
    err_hurwitz = [r['error_hurwitz_pct'] for r in he_results]

    ax1.semilogy(n_max_vals, err_raw, 'ko-', label='Uncorrected', linewidth=2, markersize=8)
    ax1.semilogy(n_max_vals, err_schwartz, 'bs-', label='Schwartz correction', linewidth=2, markersize=8)
    ax1.semilogy(n_max_vals, err_hurwitz, 'r^-', label='Hurwitz correction', linewidth=2, markersize=8)
    ax1.axhline(y=0.004, color='green', linestyle='--', alpha=0.7, label='2D var cusp (l_max=4)')
    ax1.axhline(y=0.19, color='gray', linestyle=':', alpha=0.7, label='Adiabatic floor')
    ax1.set_xlabel('n_max', fontsize=12)
    ax1.set_ylabel('Error (%)', fontsize=12)
    ax1.set_title('He Ground State: DBBSC Cusp Correction', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(n_max_vals)

    # Plot correction magnitude
    coal_vals = [r['coalescence'] for r in he_results]
    dE_h = [abs(r['delta_E_hurwitz']) for r in he_results]

    ax2.semilogy(n_max_vals, dE_h, 'r^-', label='|dE_Hurwitz|', linewidth=2, markersize=8)
    ax2.set_xlabel('n_max', fontsize=12)
    ax2.set_ylabel('|dE| (Ha)', fontsize=12)
    ax2.set_title('Cusp Correction Magnitude', fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(n_max_vals)

    ax2_twin = ax2.twinx()
    ax2_twin.plot(n_max_vals, coal_vals, 'gD--', label='<d3(r12)>', linewidth=1.5, markersize=6)
    ax2_twin.set_ylabel('Coalescence (bohr^-3)', fontsize=12, color='green')
    ax2_twin.tick_params(axis='y', labelcolor='green')
    ax2_twin.legend(loc='upper left', fontsize=10)

    plt.tight_layout()
    path1 = os.path.join(plot_dir, 'cusp_correction_convergence.png')
    plt.savefig(path1, dpi=150)
    print(f"\nSaved: {path1}")
    plt.close()

    # Plot 2: Z scaling
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    z_vals = [r['Z'] for r in z_results]
    coal_z = [r['coalescence'] for r in z_results]
    err_raw_z = [r['error_raw_pct'] for r in z_results if r['error_raw_pct'] is not None]
    err_corr_z = [r['error_hurwitz_pct'] for r in z_results if r['error_hurwitz_pct'] is not None]
    z_with_exact = [r['Z'] for r in z_results if r['error_raw_pct'] is not None]

    ax1.loglog(z_vals, coal_z, 'ko-', linewidth=2, markersize=8)
    # Fit power law
    log_z = np.log(z_vals)
    log_coal = np.log(coal_z)
    slope, intercept = np.polyfit(log_z, log_coal, 1)
    z_fit = np.linspace(min(z_vals), max(z_vals), 100)
    ax1.loglog(z_fit, np.exp(intercept) * z_fit ** slope, 'r--',
               label=f'Z^{slope:.2f} fit', linewidth=1.5)
    ax1.set_xlabel('Z', fontsize=12)
    ax1.set_ylabel('<d3(r12)> (bohr^-3)', fontsize=12)
    ax1.set_title(f'Coalescence Density vs Z (n_max=4)', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)

    ax2.plot(z_with_exact, err_raw_z, 'ko-', label='Uncorrected', linewidth=2, markersize=8)
    ax2.plot(z_with_exact, err_corr_z, 'r^-', label='Hurwitz correction', linewidth=2, markersize=8)
    ax2.set_xlabel('Z', fontsize=12)
    ax2.set_ylabel('Error (%)', fontsize=12)
    ax2.set_title('Error vs Z (n_max=4)', fontsize=14)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    path2 = os.path.join(plot_dir, 'cusp_correction_vs_Z.png')
    plt.savefig(path2, dpi=150)
    print(f"Saved: {path2}")
    plt.close()


def main():
    """Run all parts of the cusp correction prototype."""
    print("DBBSC Cusp Correction Prototype for GeoVac")
    print("=" * 70)

    # Part 1: He convergence
    he_results = part1_he_convergence()

    # Part 2: Z scaling
    z_results = part2_z_scaling()

    # Part 3: Comparison
    comparison = part3_comparison()

    # Generate plots
    generate_plots(he_results, z_results)

    # Save results to JSON
    all_results = {
        'he_convergence': he_results,
        'z_scaling': z_results,
        'comparison': comparison,
        'metadata': {
            'description': 'DBBSC cusp correction for GeoVac graph-native CI',
            'date': '2026-04-13',
            'method': 'Schwartz/Hurwitz partial-wave extrapolation',
            'note': 'Coalescence density computed from singlet CI wavefunction',
        },
    }

    data_dir = os.path.join(project_root, 'debug', 'data')
    json_path = os.path.join(data_dir, 'cusp_correction_results.json')
    with open(json_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved results: {json_path}")


if __name__ == '__main__':
    main()
