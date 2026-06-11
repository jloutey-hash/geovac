"""Generate figures for Paper 12: Algebraic Two-Electron Integrals.

Produces:
    debug/plots/paper12_convergence.png — Neumann vs numerical D_e% convergence
    debug/plots/paper12_gap_analysis.png — Energy decomposition and gap diagnosis
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Ensure output directory exists
os.makedirs("debug/plots", exist_ok=True)

# ============================================================
# Data from convergence study (debug/data/neumann_convergence.txt)
# ============================================================

configs = [
    "j=0,l=0", "j=1,l=0", "j=1,l=1",
    "j=2,l=1", "j=2,l=2", "j=3,l=2", "j=3,l=3",
]
N_bf = np.array([1, 3, 6, 12, 27, 46, 72])
E_num = np.array([-1.073192, -1.083438, -1.103485, -1.107960, -1.138836, -1.139329, -1.139792])
E_neu = np.array([-1.107907, -1.115282, -1.129606, -1.132373, -1.160961, -1.161162, -1.161304])

D_e_exact = 0.17475
E_exact = -1.17475

D_e_pct_num = (-1.0 - E_num) / D_e_exact * 100
D_e_pct_neu = (-1.0 - E_neu) / D_e_exact * 100

# r12 comparison data
r12_N = np.array([6, 18])
r12_pct = np.array([85.7, 86.8])
r12_labels = [r"$r_{12}^1$, 60×40", r"$r_{12}^{0-2}$, 20×14"]

# ============================================================
# Figure 1: Convergence curves
# ============================================================

fig, ax = plt.subplots(1, 1, figsize=(7, 5))

ax.plot(N_bf, D_e_pct_num, 'o-', color='#d62728', linewidth=2, markersize=8,
        label='Numerical $V_{ee}$ (grid 30×20)', zorder=5)
ax.plot(N_bf, D_e_pct_neu, 's-', color='#1f77b4', linewidth=2, markersize=8,
        label='Neumann $V_{ee}$ (algebraic)', zorder=5)

# r12 comparison points
ax.scatter(r12_N, r12_pct, marker='^', s=100, color='#2ca02c',
           zorder=6, label=r'$r_{12}^p$ basis (numerical $V_{ee}$)')
for i, lbl in enumerate(r12_labels):
    ax.annotate(lbl, (r12_N[i], r12_pct[i]),
                textcoords="offset points", xytext=(10, -5),
                fontsize=8, color='#2ca02c')

# Exact line
ax.axhline(y=100, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Exact $D_e$')

# Plateau annotations
ax.axhline(y=92.4, color='#1f77b4', linestyle=':', linewidth=1, alpha=0.5)
ax.axhline(y=80.0, color='#d62728', linestyle=':', linewidth=1, alpha=0.5)

ax.annotate('Neumann plateau: 92.4%', xy=(60, 92.4), xytext=(50, 96),
            fontsize=9, color='#1f77b4',
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.2))
ax.annotate('Numerical plateau: ~80%', xy=(60, 80.1), xytext=(50, 84.5),
            fontsize=9, color='#d62728',
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.2))

# Gap annotation
ax.annotate('', xy=(74, 92.4), xytext=(74, 100),
            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
ax.text(76, 96.2, '7.6%\ngap', fontsize=9, color='gray', ha='left')

ax.set_xlabel('Number of basis functions $N$', fontsize=12)
ax.set_ylabel('$D_e$ recovered (%)', fontsize=12)
ax.set_title(r'H$_2$ Dissociation Energy: Neumann vs Numerical $V_{ee}$' + '\n'
             r'$R = 1.4011$ bohr', fontsize=13)
ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
ax.set_xlim(-2, 80)
ax.set_ylim(35, 105)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("debug/plots/paper12_convergence.png", dpi=200, bbox_inches='tight')
print("Saved debug/plots/paper12_convergence.png")
plt.close()


# ============================================================
# Figure 2: Gap analysis — improvement and dE
# ============================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel (a): Percentage point improvement
improvement = D_e_pct_neu - D_e_pct_num
ax1.bar(range(len(configs)), improvement, color='#1f77b4', alpha=0.8, edgecolor='#1a5276')
ax1.set_xticks(range(len(configs)))
ax1.set_xticklabels(configs, rotation=45, ha='right', fontsize=9)
ax1.set_ylabel('Improvement (percentage points)', fontsize=11)
ax1.set_title('(a) Neumann advantage over numerical', fontsize=12)
ax1.axhline(y=12.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(5.5, 13, 'avg ≈ 12.5 pp', fontsize=9, color='gray', ha='right')
ax1.grid(True, alpha=0.3, axis='y')

for i, (imp, n) in enumerate(zip(improvement, N_bf)):
    ax1.text(i, imp + 0.3, f'{imp:.1f}', ha='center', fontsize=8, fontweight='bold')

# Panel (b): Energy difference
dE = E_neu - E_num
ax2.plot(N_bf, -dE * 1000, 'o-', color='#e67e22', linewidth=2, markersize=8)
ax2.set_xlabel('Number of basis functions $N$', fontsize=11)
ax2.set_ylabel(r'$|E_{\mathrm{num}} - E_{\mathrm{Neu}}|$ (mHa)', fontsize=11)
ax2.set_title(r'(b) Energy difference: numerical $-$ Neumann', fontsize=12)
ax2.grid(True, alpha=0.3)

for i, (n, de) in enumerate(zip(N_bf, dE)):
    ax2.annotate(f'{-de*1000:.1f}', (n, -de*1000),
                 textcoords="offset points", xytext=(8, 5), fontsize=8)

plt.tight_layout()
plt.savefig("debug/plots/paper12_gap_analysis.png", dpi=200, bbox_inches='tight')
print("Saved debug/plots/paper12_gap_analysis.png")
plt.close()


# ============================================================
# Figure 3: Grid convergence (numerical V_ee limitation)
# ============================================================

fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

# Grid convergence data for B1 (6 terms)
grids = ['12×8', '20×14', '30×20', '40×26', '60×40']
grid_N = [96, 280, 600, 1040, 2400]  # approximate grid points
grid_pct = [66.0, 76.2, 81.4, 83.7, 85.7]
grid_time = [1, 6, 45, 83, 447]

ax.semilogx(grid_N, grid_pct, 'o-', color='#d62728', linewidth=2, markersize=8,
            label='Numerical $V_{ee}$ (6 bf)')
ax.axhline(y=61.8, color='#1f77b4', linestyle='-', linewidth=2, alpha=0.7,
           label='Neumann $V_{ee}$ (6 bf, exact)')

for i, (n, p, g) in enumerate(zip(grid_N, grid_pct, grids)):
    offset = (10, 8) if i < 3 else (10, -12)
    ax.annotate(f'{g}\n{p:.1f}%', (n, p),
                textcoords="offset points", xytext=offset, fontsize=8)

ax.set_xlabel(r'Grid points ($N_\xi \times N_\eta$)', fontsize=11)
ax.set_ylabel('$D_e$ recovered (%)', fontsize=11)
ax.set_title('Grid convergence of numerical $V_{ee}$\n(6 basis functions, $\\alpha = 1.0$)',
             fontsize=12)
ax.legend(loc='lower right', fontsize=9)
ax.set_ylim(55, 90)
ax.grid(True, alpha=0.3)

ax.annotate('Grid-limited:\n~1.3 pp per\ngrid doubling',
            xy=(800, 82), fontsize=9, color='#d62728',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fce4e4', alpha=0.8))

plt.tight_layout()
plt.savefig("debug/plots/paper12_grid_convergence.png", dpi=200, bbox_inches='tight')
print("Saved debug/plots/paper12_grid_convergence.png")
plt.close()

print("\nAll Paper 12 figures generated successfully.")
