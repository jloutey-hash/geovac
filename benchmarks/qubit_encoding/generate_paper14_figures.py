"""
Generate Paper 14 figures (pauli_scaling.png and eri_density.png).

Uses the canonical Table I and Table II data from Paper 14
(paper_14_qubit_encoding.tex). The point of this script is figure
regeneration that matches the paper exactly — the underlying
benchmarks are in benchmarks/qubit_encoding/pauli_comparison.py
(GeoVac side) and in published Gaussian references (Gaussian side).

Output: papers/group4_quantum_computing/paper_14_figures/
  - pauli_scaling.png
  - eri_density.png
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Paper 14 Table I: Pauli term counts
# ----------------------------------------------------------------------

GAUSSIAN_POINTS = [
    # (label, M, Q, Pauli, ERI_density_pct)
    ('He STO-3G',     1,   2,       4, 100.0),
    ('H$_2$ STO-3G',  2,   4,      15,  50.0),
    ('6-31G$^\\dagger$', 4, 8,    201,  50.0),
    ('He cc-pVDZ',    5,  10,     156, 100.0),
    ('cc-pVDZ$^\\dagger$', 10, 20, 23189, 100.0),
    ('He cc-pVTZ',   14,  28,   21607, 100.0),
]

GEOVAC_POINTS = [
    # (label, M, Q, Pauli, ERI_density_pct)
    (r'$n_{\max}{=}2$',  5,  10,     120, 10.4),
    (r'$n_{\max}{=}3$', 14,  28,    2659,  3.9),
    (r'$n_{\max}{=}4$', 30,  60,   31039,  2.1),
    (r'$n_{\max}{=}5$', 55, 110,  227338,  1.3),
]


def fit_power_law(x, y):
    log_x = np.log(np.asarray(x, dtype=float))
    log_y = np.log(np.asarray(y, dtype=float))
    slope, intercept = np.polyfit(log_x, log_y, 1)
    return slope, np.exp(intercept)


def make_pauli_scaling(out_path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    gauss_q = np.array([p[2] for p in GAUSSIAN_POINTS])
    gauss_p = np.array([p[3] for p in GAUSSIAN_POINTS])
    geov_q  = np.array([p[2] for p in GEOVAC_POINTS])
    geov_p  = np.array([p[3] for p in GEOVAC_POINTS])

    geov_alpha, geov_a = fit_power_law(geov_q, geov_p)

    ax.loglog(gauss_q, gauss_p, 'rs', markersize=10,
              label='Gaussian (published / synthetic)')
    ax.loglog(geov_q, geov_p, 'bo', markersize=10,
              label=fr'GeoVac (lattice): $O(Q^{{{geov_alpha:.2f}}})$')

    q_fit = np.geomspace(geov_q.min() * 0.7, geov_q.max() * 1.4, 200)
    ax.loglog(q_fit, geov_a * q_fit ** geov_alpha, 'b--', alpha=0.4)

    # Reference scaling lines (Q^2, Q^4)
    q_ref = np.geomspace(2, 130, 100)
    ax.loglog(q_ref, 5 * q_ref ** 2, 'k:', alpha=0.25)
    ax.text(60, 5 * 60**2 * 1.6, r'$\sim Q^2$', fontsize=10, alpha=0.4)
    ax.loglog(q_ref, 0.3 * q_ref ** 4, 'k:', alpha=0.25)
    ax.text(15, 0.3 * 15**4 * 1.6, r'$\sim Q^4$', fontsize=10, alpha=0.4)

    # Annotate
    for label, M, Q, P, _ in GAUSSIAN_POINTS:
        ax.annotate(label, (Q, P), textcoords='offset points',
                    xytext=(7, 5), fontsize=8, color='darkred')
    for label, M, Q, P, _ in GEOVAC_POINTS:
        ax.annotate(label, (Q, P), textcoords='offset points',
                    xytext=(7, -11), fontsize=8, color='navy')

    ax.set_xlabel(r'Number of qubits $Q$', fontsize=13)
    ax.set_ylabel('Pauli string terms (Jordan--Wigner)', fontsize=13)
    ax.set_title('Pauli term count: GeoVac lattice vs Gaussian basis sets',
                 fontsize=13)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=11, loc='lower right')

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'pauli_scaling.png written: {out_path}')
    print(f'  GeoVac power-law exponent (Table I data): {geov_alpha:.3f}')


def make_eri_density(out_path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    gauss_M = np.array([p[1] for p in GAUSSIAN_POINTS], dtype=float)
    gauss_d = np.array([p[4] for p in GAUSSIAN_POINTS], dtype=float)
    geov_M  = np.array([p[1] for p in GEOVAC_POINTS], dtype=float)
    geov_d  = np.array([p[4] for p in GEOVAC_POINTS], dtype=float)

    ax.semilogy(gauss_M, gauss_d, 'rs', markersize=10,
                label='Gaussian (published / synthetic)')
    ax.semilogy(geov_M, geov_d, 'bo-', markersize=10, linewidth=1.5,
                label='GeoVac (lattice)')

    # 1/M^2 reference normalized to the n_max=2 point (M=5, d=10.4%)
    M_ref = np.linspace(geov_M.min(), geov_M.max(), 100)
    d_ref = 10.4 * (5.0 / M_ref) ** 2
    ax.semilogy(M_ref, d_ref, 'k--', alpha=0.4, label=r'$\sim 1/M^2$ reference')

    # Annotate
    for label, M, _, _, d in GAUSSIAN_POINTS:
        ax.annotate(label, (M, d), textcoords='offset points',
                    xytext=(7, 5), fontsize=8, color='darkred')
    for label, M, _, _, d in GEOVAC_POINTS:
        ax.annotate(label, (M, d), textcoords='offset points',
                    xytext=(7, -11), fontsize=8, color='navy')

    ax.set_xlabel(r'Number of spatial orbitals $M$', fontsize=13)
    ax.set_ylabel('ERI tensor density (\\%)', fontsize=13)
    ax.set_title('ERI density: angular selection rules give $1/M^2$ decay',
                 fontsize=13)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=11, loc='upper right')
    ax.set_xlim(0, 60)
    ax.set_ylim(0.5, 200)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'eri_density.png written: {out_path}')


def main():
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..', '..',
        'papers', 'group4_quantum_computing', 'paper_14_figures',
    )
    out_dir = os.path.normpath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    make_pauli_scaling(os.path.join(out_dir, 'pauli_scaling.png'))
    make_eri_density(os.path.join(out_dir, 'eri_density.png'))


if __name__ == '__main__':
    main()
