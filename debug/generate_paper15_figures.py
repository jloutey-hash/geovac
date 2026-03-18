"""
Generate publication-quality figures for Paper 15.
Output: papers/core/paper_15_figures/{level4_coordinates,convergence_lmax,adiabatic_potential}.png
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc, Circle
from matplotlib.lines import Line2D
import os
import sys

# Output directory
OUT = os.path.join(os.path.dirname(__file__), '..', 'papers', 'core', 'paper_15_figures')
os.makedirs(OUT, exist_ok=True)

# Colors (colorblind-friendly)
BLUE = '#0077BB'
RED = '#CC3311'
GRAY = '#BBBBBB'
DARK = '#333333'
GREEN = '#009988'

plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
})


def figure1_coordinates():
    """Schematic of Level 4 coordinate system."""
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 3.0))
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-1.0, 3.2)
    ax.set_aspect('equal')
    ax.axis('off')

    # Nuclei positions
    R = 1.4
    xA, yA = -R/2, 0
    xB, yB = R/2, 0

    # Draw internuclear axis
    ax.plot([xA - 0.5, xB + 0.5], [0, 0], '-', color=GRAY, lw=0.8, zorder=0)

    # Draw nuclei
    for (x, y, label) in [(xA, yA, '$A$'), (xB, yB, '$B$')]:
        circle = Circle((x, y), 0.08, fc=DARK, ec=DARK, zorder=5)
        ax.add_patch(circle)
        ax.text(x, y - 0.25, label, ha='center', va='top', fontsize=10, fontweight='bold')

    # R label
    ax.annotate('', xy=(xB, -0.45), xytext=(xA, -0.45),
                arrowprops=dict(arrowstyle='<->', color=DARK, lw=1.0))
    ax.text(0, -0.60, '$R$', ha='center', va='top', fontsize=11)

    # Midpoint
    ax.plot(0, 0, 'o', color=GRAY, ms=3, zorder=3)

    # Electron 1
    e1_x, e1_y = 1.2, 1.8
    r1 = np.sqrt(e1_x**2 + e1_y**2)
    ax.plot(e1_x, e1_y, 'o', color=BLUE, ms=6, zorder=5)
    ax.text(e1_x + 0.15, e1_y + 0.1, '$e_1$', color=BLUE, fontsize=10, fontweight='bold')

    # Electron 2
    e2_x, e2_y = -0.6, 2.4
    r2 = np.sqrt(e2_x**2 + e2_y**2)
    ax.plot(e2_x, e2_y, 'o', color=RED, ms=6, zorder=5)
    ax.text(e2_x - 0.30, e2_y + 0.1, '$e_2$', color=RED, fontsize=10, fontweight='bold')

    # r1 vector (midpoint to e1)
    ax.annotate('', xy=(e1_x, e1_y), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.2))
    ax.text(0.75, 0.75, '$r_1$', color=BLUE, fontsize=10, ha='center')

    # r2 vector (midpoint to e2)
    ax.annotate('', xy=(e2_x, e2_y), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=RED, lw=1.2))
    ax.text(-0.45, 1.05, '$r_2$', color=RED, fontsize=10, ha='center')

    # Hyperradius arc
    theta1 = np.arctan2(e1_y, e1_x)
    theta2 = np.arctan2(e2_y, e2_x)
    Re = np.sqrt(r1**2 + r2**2)
    # Draw a small arc segment at radius ~ Re * 0.3 to indicate Re
    arc_r = 0.45 * Re
    arc_angles = np.linspace(0, np.pi/2, 40)
    ax.plot(arc_r * np.cos(arc_angles), arc_r * np.sin(arc_angles),
            '--', color=GREEN, lw=1.0, zorder=2)
    ax.text(0.25, arc_r * 0.65, r'$R_e = \sqrt{r_1^2 + r_2^2}$',
            color=GREEN, fontsize=8, ha='left', style='italic')

    # Alpha angle arc between r1 and r2 directions (in the r1-r2 plane)
    # Show alpha as angle from r1 direction to r2 direction
    alpha_arc = np.linspace(theta1, theta2, 30)
    alpha_r = 0.55
    ax.plot(alpha_r * np.cos(alpha_arc), alpha_r * np.sin(alpha_arc),
            '-', color=DARK, lw=1.0)
    mid_a = (theta1 + theta2) / 2
    ax.text(alpha_r * np.cos(mid_a) + 0.15, alpha_r * np.sin(mid_a),
            r'$\alpha$', fontsize=10, ha='left', color=DARK)

    # Theta1 angle (from z-axis to r1)
    th1_arc = np.linspace(np.pi/2, theta1, 20)
    th1_r = 0.5
    ax.plot(th1_r * np.cos(th1_arc), th1_r * np.sin(th1_arc),
            '-', color=BLUE, lw=0.8, alpha=0.7)
    # Label at edge of arc
    ax.text(th1_r * 1.15, th1_r * 0.85, r'$\theta_1$', fontsize=9, color=BLUE)

    # Theta2 angle (from z-axis to r2)
    th2_arc = np.linspace(np.pi/2, theta2, 20)
    th2_r = 0.65
    ax.plot(th2_r * np.cos(th2_arc), th2_r * np.sin(th2_arc),
            '-', color=RED, lw=0.8, alpha=0.7)
    ax.text(-0.25, th2_r * 1.1, r'$\theta_2$', fontsize=9, color=RED)

    # Cusp locus: dashed line at alpha = pi/4 (equal distance from midpoint)
    cusp_len = 2.8
    cusp_angle = np.pi / 4
    ax.plot([0, cusp_len * np.cos(cusp_angle)], [0, cusp_len * np.sin(cusp_angle)],
            ':', color=GRAY, lw=1.2, zorder=1)
    ax.text(cusp_len * np.cos(cusp_angle) * 0.72, cusp_len * np.sin(cusp_angle) * 0.72 + 0.2,
            r'cusp: $\alpha = \pi/4$', fontsize=8, color='#666666',
            ha='center', rotation=45)

    # z-axis label
    ax.text(xB + 0.65, -0.05, '$z$', fontsize=10, ha='left', color=GRAY)

    fig.savefig(os.path.join(OUT, 'level4_coordinates.png'), dpi=300)
    plt.close(fig)
    print("  [OK] level4_coordinates.png")


def figure2_convergence():
    """Bar chart of D_e% vs l_max."""
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.8))

    lmax = [0, 1, 2, 3, 4]
    De_pct = [30.8, 37.3, 87.8, 88.5, 95.5]

    bars = ax.bar(lmax, De_pct, width=0.6, color=BLUE, edgecolor='white', linewidth=0.5, zorder=3)

    # Horizontal reference lines
    ax.axhline(y=92.4, color=RED, ls='--', lw=1.2, zorder=2)
    ax.text(4.35, 92.4, 'Paper 12\n(92.4%)', fontsize=7.5, color=RED,
            va='center', ha='left')

    ax.axhline(y=100, color=GRAY, ls='--', lw=1.0, zorder=2)
    ax.text(4.35, 100, 'Exact', fontsize=7.5, color='#666666',
            va='center', ha='left')

    # Annotate the quadrupole jump
    ax.annotate('quadrupole\nopens',
                xy=(2, 87.8), xytext=(0.3, 75),
                fontsize=7.5, color=DARK, ha='center',
                arrowprops=dict(arrowstyle='->', color=DARK, lw=0.8))

    # Value labels on bars — offset those near the Paper 12 line
    for i, (x, y) in enumerate(zip(lmax, De_pct)):
        if 85 < y < 94:
            # Near Paper 12 line — put label inside bar, lower
            ax.text(x, 78, f'{y:.1f}%', ha='center', va='top', fontsize=7.5,
                    fontweight='bold', color=DARK)
        else:
            ax.text(x, y + 1.5, f'{y:.1f}%', ha='center', va='bottom', fontsize=7.5,
                    fontweight='bold', color=DARK)

    ax.set_xlabel(r'$l_{\max}$')
    ax.set_ylabel(r'$D_e / D_e^{\mathrm{exact}}$ (%)')
    ax.set_xticks(lmax)
    ax.set_ylim(0, 115)
    ax.set_xlim(-0.7, 5.5)

    # Clean grid
    ax.yaxis.grid(True, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.savefig(os.path.join(OUT, 'convergence_lmax.png'), dpi=300)
    plt.close(fig)
    print("  [OK] convergence_lmax.png")


def figure3_adiabatic():
    """Adiabatic potential curves U(R_e) at R=1.4 bohr."""
    # Import the solver
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from geovac.level4_multichannel import solve_angular_multichannel

    R = 1.4
    Z = 1.0

    # R_e grid
    Re_vals = np.concatenate([
        np.linspace(0.3, 1.0, 15),
        np.linspace(1.0, 3.0, 30),
        np.linspace(3.0, 8.0, 20),
    ])

    # Compute for l_max=0, 2, 4
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.8))

    for lm, color, label in [(0, GRAY, r'$l_{\max}=0$'),
                               (2, RED, r'$l_{\max}=2$'),
                               (4, BLUE, r'$l_{\max}=4$')]:
        U_vals = []
        for Re in Re_vals:
            rho = R / (2 * Re)
            try:
                mu, _, _, _ = solve_angular_multichannel(
                    rho, Re, l_max=lm, n_alpha=80
                )
                U = (mu[0] + 15.0/8.0) / Re**2
            except Exception:
                U = np.nan
            U_vals.append(U)

        U_arr = np.array(U_vals)
        # Add V_NN
        U_total = U_arr + Z**2 / R
        ax.plot(Re_vals, U_total, '-', color=color, lw=1.3, label=label, zorder=3)

    # Reference: 2*E(H) = -1.0
    ax.axhline(y=-1.0, color=GRAY, ls=':', lw=0.8, zorder=1)
    ax.text(7.5, -0.97, r'$2E_{\mathrm{H}}$', fontsize=7.5, color='#666666', ha='right')

    ax.set_xlabel(r'$R_e$ (bohr)')
    ax.set_ylabel(r'$U(R_e) + V_{NN}$ (Ha)')
    ax.set_xlim(0.3, 8.0)
    ax.set_ylim(-3.5, 0.5)
    ax.legend(loc='lower right', framealpha=0.9, edgecolor='none')

    ax.yaxis.grid(True, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.savefig(os.path.join(OUT, 'adiabatic_potential.png'), dpi=300)
    plt.close(fig)
    print("  [OK] adiabatic_potential.png")


if __name__ == '__main__':
    print("Generating Paper 15 figures...")
    figure1_coordinates()
    figure2_convergence()
    figure3_adiabatic()
    print("Done. Figures saved to papers/core/paper_15_figures/")
