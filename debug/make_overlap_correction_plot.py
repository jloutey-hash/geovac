"""Generate overlap correction scan plot from completed data."""
import os
import sys
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Data from scan results
data = {
    0.0:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [0.4642, 0.2803, 0.1706, 0.1118, 0.0852, 0.0724, 0.0416]},
    0.05: {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [0.3941, 0.2298, 0.1329, 0.0819, 0.0743, 0.0636, 0.0359]},
    0.1:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [0.3240, 0.1794, 0.0952, 0.0681, 0.0635, 0.0547, 0.0301]},
    0.2:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [0.1840, 0.0788, 0.0317, 0.0405, 0.0418, 0.0371, 0.0188]},
    0.5:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [-0.2354, -0.1433, -0.0802, -0.0421, -0.0229, -0.0153, -0.0149]},
    1.0:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [-0.4408, -0.2991, -0.2449, -0.1786, -0.1296, -0.1016, -0.0699]},
    2.0:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [-0.6448, -0.3933, -0.2930, -0.2356, -0.1995, -0.1843, -0.1778]},
    5.0:  {'R': [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
           'D_cp': [-1.2398, -0.6597, -0.4252, -0.3185, -0.2647, -0.2406, -0.2335]},
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

colors = plt.cm.viridis(np.linspace(0, 0.9, len(data)))

for (lam, d), color in zip(data.items(), colors):
    label = f'lam={lam}'
    ax1.plot(d['R'], d['D_cp'], 'o-', color=color, label=label, markersize=4)

# Mark experiment
ax1.axhline(0.092, color='red', linestyle='--', alpha=0.6, label='expt D_cp=0.092')
ax1.axvline(3.015, color='red', linestyle=':', alpha=0.6, label='expt R_eq=3.015')
ax1.axhline(0.0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)

ax1.set_xlabel('R (bohr)')
ax1.set_ylabel('D_cp (Ha)')
ax1.set_title('CP-Corrected Binding Energy')
ax1.legend(fontsize=7, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_ylim(-0.7, 0.5)

# Zoom on the lambda=0.2 near-equilibrium
ax2.plot(data[0.0]['R'], data[0.0]['D_cp'], 'o-', color='blue', label='lam=0.0 (baseline)', markersize=5)
ax2.plot(data[0.1]['R'], data[0.1]['D_cp'], 'o-', color='green', label='lam=0.1', markersize=5)
ax2.plot(data[0.2]['R'], data[0.2]['D_cp'], 'o-', color='orange', label='lam=0.2 (near-eq)', markersize=5, linewidth=2)
ax2.plot(data[0.5]['R'], data[0.5]['D_cp'], 'o-', color='red', label='lam=0.5 (unbound)', markersize=5)

ax2.axhline(0.092, color='red', linestyle='--', alpha=0.5, label='expt D_cp')
ax2.axvline(3.015, color='red', linestyle=':', alpha=0.5, label='expt R_eq')
ax2.axhline(0.0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)

ax2.set_xlabel('R (bohr)')
ax2.set_ylabel('D_cp (Ha)')
ax2.set_title('Zoom: lam=0.2 near-equilibrium')
ax2.legend(fontsize=7, loc='lower right')
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-0.30, 0.35)

fig.suptitle('LiH Overlap Kinetic Correction: lambda scan (nmax=3)', fontsize=13)
fig.tight_layout()

outpath = os.path.join(PROJECT_ROOT, "debug", "plots", "lih_overlap_correction.png")
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Plot saved to {outpath}")
