"""
Post-process Task 10 results captured from the main run.
Saves JSON and generates all plots.
"""

import json
import numpy as np
from pathlib import Path

R_EXP = 3.015

# === Captured results from the run ===
results = {
    "task": "Task 10: Self-Consistent Projected PK Iteration",
    "parameters": {
        "Z_A_eff": 1.0, "Z_B": 1.0, "Z_A_bare": 3.0,
        "PK_A": 6.93, "PK_B": 6.80,
        "damping": 0.5, "max_iter": 20, "tol": 1e-4,
        "n_alpha": 100, "n_Re": 25,
    },

    # Part 2: Fixed geometry convergence
    "fixed_geometry": {
        "lmax2_R3.015_Re1.5": {
            "l_max": 2, "R": 3.015, "R_e": 1.5,
            "mu0_ldependent": -3.566118,
            "mu0_selfconsistent": -3.599076,
            "p1": 0.9598, "p2": 0.9598,
            "converged": True, "n_iter": 10,
            "history": [
                {"iter": 0, "mu0": -3.566118, "p1": 0.9592, "p2": 0.9592, "delta": 0.040768},
                {"iter": 1, "mu0": -3.582847, "p1": 0.9595, "p2": 0.9595, "delta": 0.020089},
                {"iter": 2, "mu0": -3.591102, "p1": 0.9597, "p2": 0.9597, "delta": 0.009900},
                {"iter": 3, "mu0": -3.595174, "p1": 0.9597, "p2": 0.9597, "delta": 0.004879},
                {"iter": 4, "mu0": -3.597181, "p1": 0.9598, "p2": 0.9598, "delta": 0.002405},
                {"iter": 5, "mu0": -3.598171, "p1": 0.9598, "p2": 0.9598, "delta": 0.001185},
                {"iter": 6, "mu0": -3.598658, "p1": 0.9598, "p2": 0.9598, "delta": 0.000584},
                {"iter": 7, "mu0": -3.598899, "p1": 0.9598, "p2": 0.9598, "delta": 0.000288},
                {"iter": 8, "mu0": -3.599017, "p1": 0.9598, "p2": 0.9598, "delta": 0.000142},
                {"iter": 9, "mu0": -3.599076, "p1": 0.9598, "p2": 0.9598, "delta": 0.000070},
            ],
        },
        "lmax2_R5.0_Re2.5": {
            "l_max": 2, "R": 5.0, "R_e": 2.5,
            "mu0_ldependent": -7.639460,
            "mu0_selfconsistent": -7.654165,
            "p1": 0.9045, "p2": 0.9045,
            "converged": True, "n_iter": 11,
            "history": [
                {"iter": 0, "mu0": -7.639460, "p1": 0.9039, "p2": 0.9039, "delta": 0.096106},
                {"iter": 1, "mu0": -7.646852, "p1": 0.9042, "p2": 0.9042, "delta": 0.047772},
                {"iter": 2, "mu0": -7.650534, "p1": 0.9043, "p2": 0.9043, "delta": 0.023747},
                {"iter": 3, "mu0": -7.652367, "p1": 0.9044, "p2": 0.9044, "delta": 0.011804},
                {"iter": 4, "mu0": -7.653278, "p1": 0.9044, "p2": 0.9044, "delta": 0.005868},
                {"iter": 5, "mu0": -7.653731, "p1": 0.9044, "p2": 0.9044, "delta": 0.002917},
                {"iter": 6, "mu0": -7.653956, "p1": 0.9044, "p2": 0.9044, "delta": 0.001450},
                {"iter": 7, "mu0": -7.654068, "p1": 0.9044, "p2": 0.9044, "delta": 0.000721},
                {"iter": 8, "mu0": -7.654124, "p1": 0.9044, "p2": 0.9044, "delta": 0.000358},
                {"iter": 9, "mu0": -7.654151, "p1": 0.9045, "p2": 0.9045, "delta": 0.000178},
                {"iter": 10, "mu0": -7.654165, "p1": 0.9045, "p2": 0.9045, "delta": 0.000089},
            ],
        },
        "lmax4_R3.015_Re1.5": {
            "l_max": 4, "R": 3.015, "R_e": 1.5,
            "mu0_ldependent": -3.627408,
            "mu0_selfconsistent": -3.666172,
            "p1": 0.9522, "p2": 0.9522,
            "converged": True, "n_iter": 10,
            "history": [
                {"iter": 0, "mu0": -3.627408, "p1": 0.9514, "p2": 0.9514, "delta": 0.048610},
                {"iter": 1, "mu0": -3.647133, "p1": 0.9518, "p2": 0.9518, "delta": 0.023885},
                {"iter": 2, "mu0": -3.656842, "p1": 0.9520, "p2": 0.9520, "delta": 0.011738},
                {"iter": 3, "mu0": -3.661619, "p1": 0.9521, "p2": 0.9521, "delta": 0.005769},
                {"iter": 4, "mu0": -3.663967, "p1": 0.9522, "p2": 0.9522, "delta": 0.002835},
                {"iter": 5, "mu0": -3.665121, "p1": 0.9522, "p2": 0.9522, "delta": 0.001394},
                {"iter": 6, "mu0": -3.665689, "p1": 0.9522, "p2": 0.9522, "delta": 0.000685},
                {"iter": 7, "mu0": -3.665968, "p1": 0.9522, "p2": 0.9522, "delta": 0.000337},
                {"iter": 8, "mu0": -3.666105, "p1": 0.9522, "p2": 0.9522, "delta": 0.000165},
                {"iter": 9, "mu0": -3.666172, "p1": 0.9522, "p2": 0.9522, "delta": 0.000081},
            ],
        },
        "lmax4_R5.0_Re2.5": {
            "l_max": 4, "R": 5.0, "R_e": 2.5,
            "mu0_ldependent": -7.911764,
            "mu0_selfconsistent": -7.930474,
            "p1": 0.8729, "p2": 0.8729,
            "converged": True, "n_iter": 12,
            "history": [
                {"iter": 0, "mu0": -7.911764, "p1": 0.8719, "p2": 0.8719, "delta": 0.128090},
                {"iter": 1, "mu0": -7.921175, "p1": 0.8724, "p2": 0.8724, "delta": 0.063570},
                {"iter": 2, "mu0": -7.925859, "p1": 0.8726, "p2": 0.8726, "delta": 0.031549},
                {"iter": 3, "mu0": -7.928187, "p1": 0.8727, "p2": 0.8727, "delta": 0.015658},
                {"iter": 4, "mu0": -7.929343, "p1": 0.8728, "p2": 0.8728, "delta": 0.007771},
                {"iter": 5, "mu0": -7.929917, "p1": 0.8728, "p2": 0.8728, "delta": 0.003857},
                {"iter": 6, "mu0": -7.930202, "p1": 0.8728, "p2": 0.8728, "delta": 0.001914},
                {"iter": 7, "mu0": -7.930343, "p1": 0.8728, "p2": 0.8728, "delta": 0.000950},
                {"iter": 8, "mu0": -7.930414, "p1": 0.8728, "p2": 0.8728, "delta": 0.000471},
                {"iter": 9, "mu0": -7.930448, "p1": 0.8729, "p2": 0.8729, "delta": 0.000234},
                {"iter": 10, "mu0": -7.930466, "p1": 0.8729, "p2": 0.8729, "delta": 0.000116},
                {"iter": 11, "mu0": -7.930474, "p1": 0.8729, "p2": 0.8729, "delta": 0.000058},
            ],
        },
    },

    # Part 2b: PK weights vs R
    "pk_vs_R": {
        "lmax2": [
            {"R": 2.0, "R_e_opt": 2.20, "p1": 0.8813, "p2": 0.8813},
            {"R": 2.5, "R_e_opt": 2.47, "p1": 0.8359, "p2": 0.8359},
            {"R": 3.0, "R_e_opt": 2.47, "p1": 0.7844, "p2": 0.7844},
            {"R": 3.5, "R_e_opt": 2.78, "p1": 0.7424, "p2": 0.7424},
            {"R": 4.0, "R_e_opt": 3.13, "p1": 0.6959, "p2": 0.6959},
            {"R": 4.5, "R_e_opt": 3.13, "p1": 0.6995, "p2": 0.6995},
            {"R": 5.0, "R_e_opt": 3.52, "p1": 0.6494, "p2": 0.6494},
            {"R": 5.5, "R_e_opt": 3.96, "p1": 0.5974, "p2": 0.5974},
            {"R": 6.5, "R_e_opt": 4.45, "p1": 0.5664, "p2": 0.5664},
        ],
        "lmax4": [
            {"R": 2.0, "R_e_opt": 2.20, "p1": 0.8494, "p2": 0.8494},
            {"R": 2.5, "R_e_opt": 2.20, "p1": 0.7415, "p2": 0.7415},
            {"R": 3.0, "R_e_opt": 2.47, "p1": 0.6615, "p2": 0.6615},
            {"R": 3.5, "R_e_opt": 2.78, "p1": 0.5782, "p2": 0.5782},
            {"R": 4.0, "R_e_opt": 3.13, "p1": 0.4970, "p2": 0.4970},
            {"R": 4.5, "R_e_opt": 3.13, "p1": 0.5035, "p2": 0.5035},
            {"R": 5.0, "R_e_opt": 3.52, "p1": 0.4256, "p2": 0.4256},
            {"R": 5.5, "R_e_opt": 3.96, "p1": 0.3567, "p2": 0.3567},
            {"R": 6.5, "R_e_opt": 4.45, "p1": 0.3278, "p2": 0.3278},
        ],
    },

    # Part 3: PES sweep results
    "pes_sweep": {
        "2": {
            "l_dependent": {"R_eq": 2.803, "E_min": -1.220212, "time": 13.9},
            "selfconsistent": {"R_eq": 2.878, "E_min": -1.256245, "time": 173.9},
        },
        "3": {
            "l_dependent": {"R_eq": 3.015, "E_min": -1.272553, "time": 30.2},
            "selfconsistent": {"R_eq": 3.026, "E_min": -1.318572, "time": 389.7},
        },
        "4": {
            "l_dependent": {"R_eq": 3.423, "E_min": -1.342805, "time": 65.0},
            "selfconsistent": {"R_eq": 3.192, "E_min": -1.387480, "time": 1027.0},
        },
    },

    # Part 4: R-dependent PK
    "r_dependent": {
        "cap1.5": {
            "2": {"R_eq": 2.837, "E_min": -1.220669, "time": 13.8},
            "3": {"R_eq": 2.975, "E_min": -1.272719, "time": 30.1},
            "4": {"R_eq": 3.075, "E_min": -1.334423, "time": 64.1},
        },
        "cap2.0": {
            "2": {"R_eq": 2.837, "E_min": -1.220669, "time": 13.6},
            "3": {"R_eq": 2.975, "E_min": -1.272719, "time": 30.2},
            "4": {"R_eq": 3.075, "E_min": -1.334423, "time": 64.2},
        },
    },
}

# Save JSON
outfile = Path(__file__).parent / 'data' / 'selfconsistent_pk.json'
outfile.parent.mkdir(exist_ok=True)
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Results saved to {outfile}")

# === Summary Table ===
print()
print("=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(f"{'Mode':<25s} {'l_max':>5s} {'R_eq':>8s} {'dR_eq':>8s} "
      f"{'% error':>8s} {'E_min':>12s}")
print("-" * 80)

for lm in ["2", "3", "4"]:
    for mode, label in [("l_dependent", "l-dependent"),
                         ("selfconsistent", "self-consistent")]:
        d = results["pes_sweep"][lm][mode]
        dR = d["R_eq"] - R_EXP
        pct = abs(dR) / R_EXP * 100
        print(f"{label:<25s} {lm:>5s} {d['R_eq']:>8.3f} "
              f"{dR:>+8.3f} {pct:>8.1f} {d['E_min']:>12.6f}")

for cap in ["cap1.5", "cap2.0"]:
    for lm in ["2", "3", "4"]:
        d = results["r_dependent"][cap][lm]
        dR = d["R_eq"] - R_EXP
        pct = abs(dR) / R_EXP * 100
        label = f"R-dep ({cap})"
        print(f"{label:<25s} {lm:>5s} {d['R_eq']:>8.3f} "
              f"{dR:>+8.3f} {pct:>8.1f} {d['E_min']:>12.6f}")

print(f"{'Experiment':<25s} {'':>5s} {R_EXP:>8.3f} "
      f"{0.0:>+8.3f} {0.0:>8.1f}")

# === R_eq drift analysis ===
print()
print("=" * 80)
print("R_eq DRIFT ANALYSIS (slope in bohr per l_max unit)")
print("=" * 80)

for mode, label in [("l_dependent", "l-dependent"),
                     ("selfconsistent", "self-consistent")]:
    lm_arr = np.array([2, 3, 4], dtype=float)
    req_arr = np.array([results["pes_sweep"][str(lm)][mode]["R_eq"]
                         for lm in [2, 3, 4]])
    slope = np.polyfit(lm_arr, req_arr, 1)[0]
    print(f"  {label}: slope = {slope:+.3f} bohr/l_max")

for cap in ["cap1.5", "cap2.0"]:
    lm_arr = np.array([2, 3, 4], dtype=float)
    req_arr = np.array([results["r_dependent"][cap][str(lm)]["R_eq"]
                         for lm in [2, 3, 4]])
    slope = np.polyfit(lm_arr, req_arr, 1)[0]
    print(f"  R-dep ({cap}): slope = {slope:+.3f} bohr/l_max")

# === Generate plots ===
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plot_dir = Path(__file__).parent / 'plots'
plot_dir.mkdir(exist_ok=True)

# --- Plot A: Iteration convergence ---
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Self-Consistent PK Iteration Convergence', fontsize=14)

cases = [
    ("lmax2_R3.015_Re1.5", "l_max=2, R=3.015, R_e=1.5", 0, 0),
    ("lmax2_R5.0_Re2.5", "l_max=2, R=5.0, R_e=2.5", 0, 1),
    ("lmax4_R3.015_Re1.5", "l_max=4, R=3.015, R_e=1.5", 1, 0),
    ("lmax4_R5.0_Re2.5", "l_max=4, R=5.0, R_e=2.5", 1, 1),
]

for key, title, row, col in cases:
    ax = axes[row, col]
    d = results["fixed_geometry"][key]
    hist = d["history"]
    iters = [h["iter"] for h in hist]
    mu0s = [h["mu0"] for h in hist]
    p1s = [h["p1"] for h in hist]
    p2s = [h["p2"] for h in hist]

    ax2 = ax.twinx()
    l1, = ax.plot(iters, mu0s, 'b-o', markersize=4, label='mu_0')
    l2, = ax2.plot(iters, p1s, 'r--s', markersize=3, label='p1')
    l3, = ax2.plot(iters, p2s, 'g--^', markersize=3, label='p2')

    ax.axhline(d["mu0_ldependent"], color='b', linestyle=':', alpha=0.5,
               label=f'l-dep mu0={d["mu0_ldependent"]:.3f}')

    ax.set_xlabel('Iteration')
    ax.set_ylabel('mu_0', color='b')
    ax2.set_ylabel('p1, p2', color='r')
    ax.set_title(title)

    lines = [l1, l2, l3]
    labels_l = ['mu_0', 'p1', 'p2']
    ax.legend(lines, labels_l, loc='center right', fontsize=8)

plt.tight_layout()
plt.savefig(plot_dir / 'selfconsistent_pk_convergence.png', dpi=150)
plt.close()
print(f"\nPlot A saved: {plot_dir}/selfconsistent_pk_convergence.png")

# --- Plot B: PK weights vs R ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Converged Self-Consistent PK Weights vs R', fontsize=14)

for idx, (lm_key, title) in enumerate([("lmax2", "l_max=2"),
                                          ("lmax4", "l_max=4")]):
    ax = axes[idx]
    data = results["pk_vs_R"][lm_key]
    Rs = [d["R"] for d in data]
    p1s = [d["p1"] for d in data]

    ax.plot(Rs, p1s, 'ro-', label='p1 = p2 (l=0 content)', linewidth=2)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5,
               label='l-dependent (w=1)')
    ax.axvline(R_EXP, color='k', linestyle=':', alpha=0.5,
               label=f'R_exp={R_EXP}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('Converged PK weight')
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.15)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(plot_dir / 'selfconsistent_pk_weights_vs_R.png', dpi=150)
plt.close()
print(f"Plot B saved: {plot_dir}/selfconsistent_pk_weights_vs_R.png")

# --- Plot C: R_eq vs l_max comparison ---
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_title('R_eq vs l_max: PK Mode Comparison (LiH)', fontsize=14)

lm_arr = [2, 3, 4]

# l-dependent
req_ldep = [results["pes_sweep"][str(lm)]["l_dependent"]["R_eq"]
            for lm in lm_arr]
ax.plot(lm_arr, req_ldep, 'b-o', label='l-dependent (w=1)',
        linewidth=2, markersize=10)

# Self-consistent
req_sc = [results["pes_sweep"][str(lm)]["selfconsistent"]["R_eq"]
          for lm in lm_arr]
ax.plot(lm_arr, req_sc, 'r-s', label='self-consistent',
        linewidth=2, markersize=10)

# R-dependent
for cap, color, marker in [("cap1.5", "green", "D"), ("cap2.0", "purple", "^")]:
    req_rdep = [results["r_dependent"][cap][str(lm)]["R_eq"]
                for lm in lm_arr]
    ax.plot(lm_arr, req_rdep, f'{color}', marker=marker,
            label=f'R-dependent ({cap})', linewidth=2, markersize=10,
            linestyle='--')

ax.axhline(R_EXP, color='k', linestyle=':', linewidth=2,
            label=f'Experiment ({R_EXP})')

# Also show Task 9 real solver results for context
task9_lm = [2, 3, 4, 5]
task9_req = [3.282, 3.471, 3.729, 3.970]
ax.plot(task9_lm, task9_req, 'k--x', label='Task 9 (real solver, l-dep)',
        linewidth=1.5, markersize=8, alpha=0.6)

ax.set_xlabel('l_max', fontsize=12)
ax.set_ylabel('R_eq (bohr)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xticks([2, 3, 4, 5])

plt.tight_layout()
plt.savefig(plot_dir / 'selfconsistent_pk_req_vs_lmax.png', dpi=150)
plt.close()
print(f"Plot C saved: {plot_dir}/selfconsistent_pk_req_vs_lmax.png")

# --- Plot D: R_eq drift slopes bar chart ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_title('R_eq Drift Rate by PK Mode', fontsize=14)

modes = ['l-dependent', 'self-consistent', 'R-dep (cap1.5)', 'R-dep (cap2.0)',
         'Task 9 (real)']
slopes = []

for mode_key in ["l_dependent", "selfconsistent"]:
    lm = np.array([2, 3, 4], dtype=float)
    req = np.array([results["pes_sweep"][str(int(l))][mode_key]["R_eq"]
                     for l in lm])
    slopes.append(np.polyfit(lm, req, 1)[0])

for cap in ["cap1.5", "cap2.0"]:
    lm = np.array([2, 3, 4], dtype=float)
    req = np.array([results["r_dependent"][cap][str(int(l))]["R_eq"]
                     for l in lm])
    slopes.append(np.polyfit(lm, req, 1)[0])

# Task 9 slope
lm9 = np.array([2, 3, 4, 5], dtype=float)
req9 = np.array(task9_req)
slopes.append(np.polyfit(lm9, req9, 1)[0])

colors = ['blue', 'red', 'green', 'purple', 'gray']
bars = ax.bar(modes, slopes, color=colors, alpha=0.7, edgecolor='black')

ax.set_ylabel('dR_eq / d(l_max) (bohr per unit)', fontsize=12)
ax.axhline(0, color='k', linewidth=0.5)
ax.grid(True, alpha=0.3, axis='y')

for bar, val in zip(bars, slopes):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f'{val:+.3f}', ha='center', fontsize=10)

plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.savefig(plot_dir / 'selfconsistent_pk_drift_rates.png', dpi=150)
plt.close()
print(f"Plot D saved: {plot_dir}/selfconsistent_pk_drift_rates.png")

print("\nAll plots generated successfully.")
