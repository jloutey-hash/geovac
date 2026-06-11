#!/usr/bin/env python3
"""
debug_publication_figures.py
============================
Generate three publication-quality figures for arXiv submission.

Figure 1: Energy convergence vs max_n for He and Li
Figure 2: Scaling (assembly time and NNZ vs max_n)
Figure 3: GeoVac vs PySCF comparison bar chart

Author: GeoVac Development Team, March 2026
"""

import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

PLOT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'debug', 'plots')
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'debug', 'data')
os.makedirs(PLOT_DIR, exist_ok=True)


def load_scaling_data() -> list:
    """Load results from scaling study."""
    path = os.path.join(DATA_DIR, 'slater_full_scaling.json')
    with open(path) as f:
        return json.load(f)


def figure1_convergence(data: list) -> None:
    """Figure 1: Energy convergence vs max_n for He and Li."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    HE_EXACT = -2.9037
    HE_HF = -2.8617
    LI_EXACT = -7.4781
    LI_HF = -7.4327

    # --- He panel ---
    he_full_hybrid = [r for r in data
                      if 'He' in r['label'] and 'hybrid' in r['label']
                      and 'slater_full' in r['label']]
    he_full_exact = [r for r in data
                     if 'He' in r['label'] and 'exact' in r['label']
                     and 'slater_full' in r['label']]
    he_f0 = [r for r in data
             if 'He' in r['label'] and 'slater+hybrid' in r['label']]

    if he_full_hybrid:
        ns = [r['max_n'] for r in he_full_hybrid]
        errs = [r['pct_err'] for r in he_full_hybrid]
        ax1.plot(ns, errs, 'o-', color='#2196F3', lw=2, ms=8,
                 label='slater_full + hybrid h1')

    if he_full_exact:
        ns = [r['max_n'] for r in he_full_exact]
        errs = [r['pct_err'] for r in he_full_exact]
        ax1.plot(ns, errs, 's--', color='#FF9800', lw=1.5, ms=6,
                 label='slater_full + exact h1')

    if he_f0:
        ns = [r['max_n'] for r in he_f0]
        errs = [r['pct_err'] for r in he_f0]
        ax1.plot(ns, errs, '^--', color='#9E9E9E', lw=1.5, ms=6,
                 label='slater F0 + hybrid h1')

    # Reference lines
    he_hf_err = abs(HE_HF - HE_EXACT) / abs(HE_EXACT) * 100
    ax1.axhline(he_hf_err, color='red', ls=':', lw=1, alpha=0.7,
                label=f'HF limit ({he_hf_err:.1f}%)')
    ax1.axhline(0, color='green', ls=':', lw=1, alpha=0.5)

    ax1.set_xlabel('max_n (basis size)', fontsize=12)
    ax1.set_ylabel('Error vs exact (%)', fontsize=12)
    ax1.set_title('Helium (Z=2, 2e)', fontsize=14)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_ylim(bottom=-0.2)
    ax1.grid(True, alpha=0.3)

    # --- Li panel ---
    li_full_hybrid = [r for r in data
                      if 'Li' in r['label'] and 'hybrid' in r['label']
                      and 'slater_full' in r['label']]
    li_full_exact = [r for r in data
                     if 'Li' in r['label'] and 'exact' in r['label']
                     and 'slater_full' in r['label']]

    if li_full_exact:
        ns = [r['max_n'] for r in li_full_exact]
        errs = [r['pct_err'] for r in li_full_exact]
        ax2.plot(ns, errs, 'o-', color='#2196F3', lw=2, ms=8,
                 label='slater_full + exact h1')

    if li_full_hybrid:
        ns = [r['max_n'] for r in li_full_hybrid]
        errs = [r['pct_err'] for r in li_full_hybrid]
        ax2.plot(ns, errs, 's--', color='#FF9800', lw=1.5, ms=6,
                 label='slater_full + hybrid h1 (non-monotonic)')

    # Reference lines
    li_hf_err = abs(LI_HF - LI_EXACT) / abs(LI_EXACT) * 100
    ax2.axhline(li_hf_err, color='red', ls=':', lw=1, alpha=0.7,
                label=f'HF limit ({li_hf_err:.1f}%)')
    ax2.axhline(0, color='green', ls=':', lw=1, alpha=0.5)

    ax2.set_xlabel('max_n (basis size)', fontsize=12)
    ax2.set_ylabel('Error vs exact (%)', fontsize=12)
    ax2.set_title('Lithium (Z=3, 3e)', fontsize=14)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.set_ylim(bottom=-0.2)
    ax2.grid(True, alpha=0.3)

    fig.suptitle('GeoVac FCI Energy Convergence with Full Slater Integrals',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()

    out_path = os.path.join(PLOT_DIR, 'energy_convergence_slater_full.png')
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Figure 1 saved: {out_path}")
    plt.close()


def figure2_scaling(data: list) -> None:
    """Figure 2: Scaling — assembly time and NNZ vs max_n."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    he_full = [r for r in data
               if 'He' in r['label'] and 'hybrid' in r['label']
               and 'slater_full' in r['label']]
    li_full = [r for r in data
               if 'Li' in r['label'] and 'hybrid' in r['label']
               and 'slater_full' in r['label']]

    # --- NNZ panel ---
    if he_full:
        ns = [r['max_n'] for r in he_full]
        nnz = [r['nnz'] for r in he_full]
        ax1.plot(ns, nnz, 'o-', color='#2196F3', lw=2, ms=8, label='He (2e)')

    if li_full:
        ns = [r['max_n'] for r in li_full]
        nnz = [r['nnz'] for r in li_full]
        ax1.plot(ns, nnz, 's-', color='#E91E63', lw=2, ms=8, label='Li (3e)')

    ax1.set_xlabel('max_n', fontsize=12)
    ax1.set_ylabel('Non-zero elements (NNZ)', fontsize=12)
    ax1.set_title('Hamiltonian Sparsity', fontsize=14)
    ax1.set_yscale('log')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # --- Time panel ---
    if he_full:
        ns = [r['max_n'] for r in he_full]
        times = [r['t_assembly'] for r in he_full]
        ax2.plot(ns, times, 'o-', color='#2196F3', lw=2, ms=8, label='He (2e)')

    if li_full:
        ns = [r['max_n'] for r in li_full]
        times = [r['t_assembly'] for r in li_full]
        ax2.plot(ns, times, 's-', color='#E91E63', lw=2, ms=8, label='Li (3e)')

    ax2.set_xlabel('max_n', fontsize=12)
    ax2.set_ylabel('Assembly time (s)', fontsize=12)
    ax2.set_title('Construction Time', fontsize=14)
    ax2.set_yscale('log')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    fig.suptitle('GeoVac FCI Scaling (slater_full)',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()

    out_path = os.path.join(PLOT_DIR, 'scaling_slater_full.png')
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Figure 2 saved: {out_path}")
    plt.close()


def figure3_comparison(data: list) -> None:
    """Figure 3: GeoVac vs PySCF comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Data from pyscf_comparison_ci.txt + our slater_full results
    systems = ['H (1e)', 'He (2e)', 'Li (3e)']

    # GeoVac best results
    # H: max_n=30 -> 0.57%
    # He: pick best from data
    he_best = min(
        [r for r in data if 'He' in r['label'] and 'slater_full' in r['label']],
        key=lambda r: r['pct_err'],
        default=None
    )
    # Li: pick best from data
    li_best = min(
        [r for r in data if 'Li' in r['label'] and 'slater_full' in r['label']],
        key=lambda r: r['pct_err'],
        default=None
    )

    geovac_errs = [0.57]  # H max_n=30
    geovac_labels = ['max_n=30']
    if he_best:
        geovac_errs.append(he_best['pct_err'])
        geovac_labels.append(f'max_n={he_best["max_n"]}')
    else:
        geovac_errs.append(2.08)
        geovac_labels.append('max_n=4')
    if li_best:
        geovac_errs.append(li_best['pct_err'])
        geovac_labels.append(f'max_n={li_best["max_n"]}')
    else:
        geovac_errs.append(1.15)
        geovac_labels.append('max_n=3')

    # PySCF STO-3G
    pyscf_sto3g = [6.68, 3.30, None]  # Li not available

    # PySCF cc-pVQZ
    pyscf_vqz = [0.01, 1.45, None]  # Li not available

    x = np.arange(len(systems))
    width = 0.25

    bars1 = ax.bar(x - width, geovac_errs, width, label='GeoVac (slater_full)',
                   color='#2196F3', edgecolor='black', linewidth=0.5)

    # STO-3G
    sto3g_vals = [v if v is not None else 0 for v in pyscf_sto3g]
    bars2 = ax.bar(x, sto3g_vals, width, label='PySCF STO-3G',
                   color='#FF9800', edgecolor='black', linewidth=0.5)

    # cc-pVQZ
    vqz_vals = [v if v is not None else 0 for v in pyscf_vqz]
    bars3 = ax.bar(x + width, vqz_vals, width, label='PySCF cc-pVQZ',
                   color='#4CAF50', edgecolor='black', linewidth=0.5)

    # Add value labels
    for bars, vals in [(bars1, geovac_errs), (bars2, pyscf_sto3g), (bars3, pyscf_vqz)]:
        for bar, val in zip(bars, vals):
            if val is not None and val > 0:
                ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                        f'{val:.2f}%', ha='center', va='bottom', fontsize=9)
            elif val is None:
                ax.text(bar.get_x() + bar.get_width()/2., 0.05,
                        'N/A', ha='center', va='bottom', fontsize=8, color='gray')

    ax.set_ylabel('Error vs exact energy (%)', fontsize=12)
    ax.set_title('GeoVac vs PySCF: Multi-electron Atoms',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(systems, fontsize=12)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, max(geovac_errs + [v for v in pyscf_sto3g if v]) * 1.3)

    # Add note
    ax.text(0.98, 0.95,
            'GeoVac: O(V) sparse graph\nPySCF: O(N^4) integrals',
            transform=ax.transAxes, fontsize=9, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()

    out_path = os.path.join(PLOT_DIR, 'geovac_vs_pyscf.png')
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Figure 3 saved: {out_path}")
    plt.close()


def main() -> None:
    data = load_scaling_data()
    print(f"Loaded {len(data)} results\n")

    figure1_convergence(data)
    figure2_scaling(data)
    figure3_comparison(data)

    print("\n[DONE] All figures saved to debug/plots/")


if __name__ == "__main__":
    main()
