"""
Generate publication-quality figures for Paper 14.

Figures:
  1. pauli_scaling.png  — log-log Pauli term count vs qubits
  2. eri_density.png    — ERI density vs spatial orbitals

Run from repo root:
    python papers/core/paper_14_figures/generate_figures.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Colorblind-friendly palette
BLUE = '#0077BB'
RED = '#CC3311'
GRAY = '#BBBBBB'

# Output directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def figure_1_pauli_scaling() -> None:
    """Figure 1: Log-log Pauli term count vs number of qubits."""

    # --- Data ---
    gauss_q = np.array([4, 8, 20])
    gauss_p = np.array([15, 201, 23189])
    gauss_labels = ['STO-3G', '6-31G', 'cc-pVDZ']

    geov_q = np.array([10, 28, 60, 110])
    geov_p = np.array([120, 2659, 31039, 227338])
    geov_labels = [r'$n_{\max}$=2', r'$n_{\max}$=3',
                   r'$n_{\max}$=4', r'$n_{\max}$=5']

    # --- Power-law fits ---
    def fit(x: np.ndarray, y: np.ndarray):
        b, log_a = np.polyfit(np.log(x.astype(float)),
                              np.log(y.astype(float)), 1)
        return b, np.exp(log_a)

    gauss_exp, gauss_a = fit(gauss_q, gauss_p)
    geov_exp, geov_a = fit(geov_q, geov_p)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(3.4, 2.8))  # single-column width

    # Reference lines
    q_ref = np.linspace(3, 150, 200)
    ax.loglog(q_ref, 0.8 * q_ref ** 2, color=GRAY, ls=':', lw=0.8, zorder=1)
    ax.loglog(q_ref, 0.015 * q_ref ** 4, color=GRAY, ls=':', lw=0.8, zorder=1)
    ax.text(100, 0.8 * 100**2 * 1.4, r'$\sim Q^2$',
            fontsize=7, color=GRAY, ha='center')
    ax.text(22, 0.015 * 22**4 * 2.0, r'$\sim Q^4$',
            fontsize=7, color=GRAY, ha='center')

    # Fit lines (extended)
    q_fit = np.linspace(3, 150, 200)
    ax.loglog(q_fit, gauss_a * q_fit ** gauss_exp,
              color=RED, ls='--', lw=1.0, alpha=0.4, zorder=2)
    ax.loglog(q_fit, geov_a * q_fit ** geov_exp,
              color=BLUE, ls='--', lw=1.0, alpha=0.4, zorder=2)

    # Data points
    ax.loglog(gauss_q, gauss_p, 's', color=RED, ms=6, mew=1.2,
              mfc='white', zorder=4,
              label=fr'Gaussian $O(Q^{{{gauss_exp:.2f}}})$')
    ax.loglog(geov_q, geov_p, 'o', color=BLUE, ms=6, mew=1.2,
              mfc='white', zorder=4,
              label=fr'GeoVac $O(Q^{{{geov_exp:.2f}}})$')

    # Annotations
    offset_gauss = [(8, -12), (8, -12), (8, 5)]
    for i, lbl in enumerate(gauss_labels):
        ax.annotate(lbl, (gauss_q[i], gauss_p[i]),
                    textcoords='offset points', xytext=offset_gauss[i],
                    fontsize=6.5, color=RED, ha='left')

    offset_geov = [(6, 5), (6, 5), (-8, -14), (-8, -14)]
    ha_geov = ['left', 'left', 'right', 'right']
    for i, lbl in enumerate(geov_labels):
        ax.annotate(lbl, (geov_q[i], geov_p[i]),
                    textcoords='offset points', xytext=offset_geov[i],
                    fontsize=6.5, color=BLUE, ha=ha_geov[i])

    ax.set_xlabel('Number of qubits $Q$', fontsize=9)
    ax.set_ylabel('Pauli term count', fontsize=9)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc='upper left', framealpha=0.9,
              edgecolor='none', handletextpad=0.4)
    ax.set_xlim(3, 160)
    ax.set_ylim(8, 5e5)
    ax.grid(True, which='major', alpha=0.15, lw=0.5)
    ax.grid(True, which='minor', alpha=0.07, lw=0.3)

    fig.tight_layout(pad=0.4)
    out = os.path.join(SCRIPT_DIR, 'pauli_scaling.png')
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


def figure_2_eri_density() -> None:
    """Figure 2: ERI density vs spatial orbitals M."""

    # --- Data ---
    gauss_m = np.array([2, 4, 10])
    gauss_d = np.array([50.0, 50.0, 100.0])

    geov_m = np.array([5, 14, 30, 55])
    geov_d = np.array([10.4, 3.9, 2.1, 1.3])

    # --- Reference curve 1/M^2 ---
    m_ref = np.linspace(4, 60, 200)
    # Normalize: pass through geov point at M=5 (10.4%)
    c = 10.4 * 5**2
    ref_curve = c / m_ref**2

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(3.4, 2.8))

    # Reference 1/M^2 curve
    ax.semilogy(m_ref, ref_curve, color=GRAY, ls='--', lw=1.0, zorder=1,
                label=r'$\sim 1/M^2$')

    # Gaussian
    ax.semilogy(gauss_m, gauss_d, 's', color=RED, ms=6, mew=1.2,
                mfc='white', zorder=4, label='Gaussian basis')

    # GeoVac
    ax.semilogy(geov_m, geov_d, 'o', color=BLUE, ms=6, mew=1.2,
                mfc='white', zorder=4, label='GeoVac lattice')

    # Annotations — Gaussian
    gauss_labels = ['STO-3G', '6-31G', 'cc-pVDZ']
    g_offsets = [(-6, -12), (6, -12), (-6, 6)]
    g_ha = ['right', 'left', 'right']
    for i, lbl in enumerate(gauss_labels):
        ax.annotate(lbl, (gauss_m[i], gauss_d[i]),
                    textcoords='offset points', xytext=g_offsets[i],
                    fontsize=6.5, color=RED, ha=g_ha[i])

    # Annotations — GeoVac
    geov_labels = [r'$n_{\max}$=2', r'$n_{\max}$=3',
                   r'$n_{\max}$=4', r'$n_{\max}$=5']
    v_offsets = [(6, -10), (6, -10), (6, 5), (6, 5)]
    for i, lbl in enumerate(geov_labels):
        ax.annotate(lbl, (geov_m[i], geov_d[i]),
                    textcoords='offset points', xytext=v_offsets[i],
                    fontsize=6.5, color=BLUE, ha='left')

    ax.set_xlabel('Spatial orbitals $M$', fontsize=9)
    ax.set_ylabel('ERI density (%)', fontsize=9)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc='upper right', framealpha=0.9,
              edgecolor='none', handletextpad=0.4)
    ax.set_xlim(0, 60)
    ax.set_ylim(0.5, 200)
    ax.grid(True, which='major', alpha=0.15, lw=0.5)
    ax.grid(True, which='minor', alpha=0.07, lw=0.3)

    fig.tight_layout(pad=0.4)
    out = os.path.join(SCRIPT_DIR, 'eri_density.png')
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
    })
    figure_1_pauli_scaling()
    figure_2_eri_density()
    print('Done.')
