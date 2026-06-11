"""
Task 9: PK Saturation in Real Composed Geometry Solver.
Runs the full ComposedDiatomicSolver for LiH at l_max = 2..5
with l-dependent PK and compares to diagnostic model.
"""

import sys
import time
import json
import numpy as np
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.composed_diatomic import ComposedDiatomicSolver, _v_cross_nuc_1s
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel

R_EXP = 3.015
R_GRID = np.array([2.0, 2.4, 2.7, 3.0, 3.2, 3.5, 3.8, 4.2, 4.7, 5.5, 6.5])

# ===================================================================
# Solve core once
# ===================================================================
print("=" * 76)
print("Task 9: PK Saturation in Real Composed Solver")
print("=" * 76)
print("Solving core (once)...")
t0 = time.time()
core = CoreScreening(Z=3, l_max=2, n_alpha=200)
core.solve(verbose=True)
E_core = core.energy
pk = AbInitioPK(core, n_core=2)
pk_A_val = pk.A
pk_B_val = pk.B
print(f"E_core = {E_core:.6f} Ha")
print(f"PK: A={pk_A_val:.4f}, B={pk_B_val:.4f}")
print(f"Core solve: {time.time()-t0:.1f}s")

Z_A_bare = 3.0
Z_B = 1.0
Z_A_eff = 1.0  # Z_A - n_core


def scan_pes(l_max, pk_channel_mode, n_alpha=80):
    """Scan PES at given l_max with specified PK mode."""
    pk_dict = pk.pk_dict(atom='A')
    pk_dict['channel_mode'] = pk_channel_mode
    pk_pots = [pk_dict]

    E_elec_arr = []
    R_valid = []
    for R in R_GRID:
        try:
            result = solve_level4_h2_multichannel(
                R=R, Z_A=Z_A_eff, Z_B=Z_B,
                l_max=l_max, n_alpha=n_alpha, n_Re=200,
                verbose=False, m_max=0,
                pk_potentials=pk_pots,
            )
            E_elec_arr.append(result['E_elec'])
            R_valid.append(R)
        except Exception as e:
            print(f"  R={R:.2f} FAILED: {e}")

    R_valid = np.array(R_valid)
    E_elec = np.array(E_elec_arr)
    V_NN = Z_A_bare * Z_B / R_valid
    V_cross = np.array([_v_cross_nuc_1s(Z_A_bare, 2, Z_B, R) for R in R_valid])
    E_composed = E_core + V_cross + E_elec + V_NN

    # Find R_eq by cubic interpolation
    i_min = np.argmin(E_composed)
    if 0 < i_min < len(R_valid) - 1:
        i_lo = max(0, i_min - 2)
        i_hi = min(len(R_valid), i_min + 3)
        cs = CubicSpline(R_valid[i_lo:i_hi], E_composed[i_lo:i_hi])
        R_fine = np.linspace(R_valid[i_lo], R_valid[i_hi-1], 500)
        E_fine = cs(R_fine)
        i_f = np.argmin(E_fine)
        R_eq = float(R_fine[i_f])
        E_min = float(E_fine[i_f])
    else:
        R_eq = float(R_valid[i_min])
        E_min = float(E_composed[i_min])

    return R_eq, E_min, R_valid.tolist(), E_composed.tolist(), E_elec.tolist()


# ===================================================================
# l-dependent PK sweep
# ===================================================================
print()
print("=" * 76)
print("l-dependent PK sweep")
print("=" * 76)

results = {}
for l_max in [2, 3, 4, 5]:
    t0 = time.time()
    print(f"\nl_max={l_max}...", flush=True)
    R_eq, E_min, R_arr, E_arr, Ee_arr = scan_pes(l_max, 'l_dependent', n_alpha=80)
    elapsed = time.time() - t0
    err = 100 * abs(R_eq - R_EXP) / R_EXP
    results[l_max] = {
        'R_eq': R_eq, 'E_min': E_min, 'time': elapsed,
        'R': R_arr, 'E': E_arr, 'E_elec': Ee_arr,
    }
    print(f"  R_eq = {R_eq:.4f} bohr ({err:.1f}% error), "
          f"E_min = {E_min:.6f} Ha, {elapsed:.1f}s")

    if elapsed > 600:
        print("  ** Aborting: too slow **")
        break


# ===================================================================
# Channel-blind comparison
# ===================================================================
print()
print("=" * 76)
print("Channel-blind PK comparison")
print("=" * 76)

results_blind = {}
for l_max in [2, 3, 4]:
    t0 = time.time()
    print(f"\nl_max={l_max} (blind)...", flush=True)
    R_eq, E_min, _, _, _ = scan_pes(l_max, 'channel_blind', n_alpha=80)
    elapsed = time.time() - t0
    err = 100 * abs(R_eq - R_EXP) / R_EXP
    results_blind[l_max] = {'R_eq': R_eq, 'E_min': E_min, 'time': elapsed}
    print(f"  R_eq = {R_eq:.4f} bohr ({err:.1f}% error), {elapsed:.1f}s")


# ===================================================================
# Summary tables
# ===================================================================
print()
print("=" * 76)
print("SUMMARY: l-dependent PK")
print("=" * 76)

lm_list = sorted(results.keys())
r_eqs = [results[lm]['R_eq'] for lm in lm_list]
e_mins = [results[lm]['E_min'] for lm in lm_list]

print(f"{'l_max':>5}  {'R_eq':>10}  {'dR_eq':>8}  {'%err':>7}  "
      f"{'E_min':>12}  {'Time':>8}")
print("-" * 60)
for i, lm in enumerate(lm_list):
    dr = r_eqs[i] - r_eqs[i-1] if i > 0 else 0.0
    err = 100 * abs(r_eqs[i] - R_EXP) / R_EXP
    print(f"{lm:>5}  {r_eqs[i]:>10.4f}  {dr:>+8.4f}  {err:>7.1f}%  "
          f"{e_mins[i]:>12.6f}  {results[lm]['time']:>8.1f}s")

# Channel-blind table
print()
print(f"{'l_max':>5}  {'blind R_eq':>12}  {'l-dep R_eq':>12}  {'diff':>8}")
print("-" * 50)
for lm in sorted(results_blind.keys()):
    rb = results_blind[lm]['R_eq']
    rl = results.get(lm, {}).get('R_eq', float('nan'))
    print(f"{lm:>5}  {rb:>12.4f}  {rl:>12.4f}  {rb-rl:>+8.4f}")

# Diagnostic comparison
print()
print("=" * 76)
print("COMPARISON: Real vs Diagnostic Model")
print("=" * 76)
diag = {2: 3.71, 3: 4.13, 4: 4.33, 5: 4.45}
print(f"{'l_max':>5}  {'Real':>10}  {'Diag':>10}  {'Diff':>8}  "
      f"{'Real%':>8}  {'Diag%':>8}")
print("-" * 55)
for lm in lm_list:
    rr = results[lm]['R_eq']
    rd = diag.get(lm, float('nan'))
    er = 100 * abs(rr - R_EXP) / R_EXP
    ed = 100 * abs(rd - R_EXP) / R_EXP
    print(f"{lm:>5}  {rr:>10.4f}  {rd:>10.2f}  {rr-rd:>+8.4f}  "
          f"{er:>8.1f}%  {ed:>8.1f}%")


# ===================================================================
# V_NN verification
# ===================================================================
print()
print("=" * 76)
print("V_NN VERIFICATION")
print("=" * 76)
print(f"Real solver uses V_NN = Z_A_bare * Z_B / R = {Z_A_bare} / R")
print(f"Diagnostic used V_NN = Z_A_eff * Z_B / R = 2.69 / R (WRONG)")
print(f"At R=3.015: real V_NN = {Z_A_bare/3.015:.4f}, "
      f"diag V_NN = {2.69/3.015:.4f}")
print(f"Difference: {(Z_A_bare-2.69)/3.015:.4f} Ha = extra repulsion")
v_cross = _v_cross_nuc_1s(Z_A_bare, 2, Z_B, 3.015)
print(f"V_cross_nuc at R=3.015: {v_cross:.4f} Ha "
      f"(core attraction, not in diagnostic)")
print(f"Net V_NN correction: {(Z_A_bare-2.69)/3.015 + v_cross:.4f} Ha")


# ===================================================================
# Fit analysis
# ===================================================================
if len(lm_list) >= 3:
    x = np.array(lm_list, dtype=float)
    y = np.array(r_eqs)
    ss_tot = np.sum((y - np.mean(y))**2)

    def linear(x, a, b):
        return a + b * x

    def exp_sat(x, Ri, c, b):
        return Ri - c * np.exp(-b * x)

    popt_l, _ = curve_fit(linear, x, y)
    r2_l = 1 - np.sum((y - linear(x, *popt_l))**2) / ss_tot

    R_inf = None
    r2_e = -1
    popt_e = None
    try:
        popt_e, _ = curve_fit(
            exp_sat, x, y,
            p0=[y[-1] + 0.5, y[-1] - y[0] + 0.3, 0.3],
            maxfev=10000,
        )
        r2_e = 1 - np.sum((y - exp_sat(x, *popt_e))**2) / ss_tot
        R_inf = popt_e[0]
    except Exception as e:
        print(f"Exp fit failed: {e}")

    print()
    print("=" * 76)
    print("FIT ANALYSIS")
    print("=" * 76)
    print(f"Linear: slope={popt_l[1]:.4f} bohr/l_max, R2={r2_l:.6f}")
    if R_inf is not None:
        print(f"Exp saturation: R_inf={R_inf:.4f} bohr "
              f"({100*abs(R_inf-R_EXP)/R_EXP:.1f}% from expt), R2={r2_e:.6f}")

    drifts = [r_eqs[i] - r_eqs[i-1] for i in range(1, len(r_eqs))]
    print(f"Drift increments: {[f'{d:+.4f}' for d in drifts]}")


# ===================================================================
# Bottom line
# ===================================================================
print()
print("=" * 76)
print("BOTTOM LINE")
print("=" * 76)

if len(r_eqs) >= 2:
    last_dr = abs(r_eqs[-1] - r_eqs[-2])
    first_dr = abs(r_eqs[1] - r_eqs[0])
    if last_dr < first_dr * 0.5:
        print(f"Drift DIMINISHING: {first_dr:.4f} -> {last_dr:.4f} "
              f"(ratio {last_dr/first_dr:.2f})")
    else:
        print(f"Drift NOT clearly diminishing: {first_dr:.4f} -> {last_dr:.4f} "
              f"(ratio {last_dr/first_dr:.2f})")

err_final = 100 * abs(r_eqs[-1] - R_EXP) / R_EXP
print(f"Best R_eq = {r_eqs[-1]:.4f} at l_max={lm_list[-1]} "
      f"({err_final:.1f}% error)")

if R_inf is not None and r2_e > r2_l and r2_e > 0.9:
    print(f"R_eq SATURATES to R_inf = {R_inf:.3f} bohr "
          f"({100*abs(R_inf-R_EXP)/R_EXP:.1f}% from experiment)")
elif R_inf is not None:
    print(f"Exponential fit R_inf = {R_inf:.3f} but R2={r2_e:.4f} "
          f"(linear R2={r2_l:.4f})")
else:
    print("Could not determine saturation behavior")


# ===================================================================
# Save data
# ===================================================================
data_dir = Path(__file__).parent / 'data'
data_dir.mkdir(exist_ok=True)

output = {
    'task': 'Task 9: PK Saturation Real Solver',
    'R_exp': R_EXP,
    'E_core': E_core,
    'pk_A': pk_A_val,
    'pk_B': pk_B_val,
    'Z_A_eff': Z_A_eff,
    'Z_A_bare': Z_A_bare,
    'l_dependent': {str(k): v for k, v in results.items()},
    'channel_blind': {str(k): v for k, v in results_blind.items()},
    'diagnostic_r_eqs': diag,
    'summary': {
        'l_max_values': lm_list,
        'R_eq_values': r_eqs,
        'E_min_values': e_mins,
    },
}

if len(lm_list) >= 3:
    output['fits'] = {
        'linear': {
            'params': popt_l.tolist(),
            'R2': r2_l,
        },
        'exponential': {
            'params': popt_e.tolist() if popt_e is not None else None,
            'R2': r2_e,
            'R_inf': R_inf,
        },
    }

json_path = data_dir / 'pk_saturation_real.json'
with open(json_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nData saved to {json_path}")


# ===================================================================
# Plots
# ===================================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot_dir = Path(__file__).parent / 'plots'
    plot_dir.mkdir(exist_ok=True)

    # --- Plot A: R_eq vs l_max (real + diagnostic + experiment) ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    ax.plot(lm_list, r_eqs, 'o-', color='#1f77b4', linewidth=2.5,
            markersize=10, label='Real solver (l-dependent PK)', zorder=5)

    diag_lm = sorted(diag.keys())
    diag_r = [diag[lm] for lm in diag_lm]
    ax.plot(diag_lm, diag_r, 's--', color='#ff7f0e', linewidth=2,
            markersize=8, alpha=0.7, label='Diagnostic (const Z_eff=2.69)')

    bl_lm = sorted(results_blind.keys())
    bl_r = [results_blind[lm]['R_eq'] for lm in bl_lm]
    ax.plot(bl_lm, bl_r, 'D-.', color='#d62728', linewidth=2,
            markersize=8, alpha=0.7, label='Channel-blind PK')

    ax.axhline(R_EXP, color='gray', linestyle='--', linewidth=1.5,
               alpha=0.7, label=f'Experiment ({R_EXP} bohr)')

    if R_inf is not None and r2_e > 0.9:
        x_fine = np.linspace(1.5, 8, 100)
        ax.plot(x_fine, exp_sat(x_fine, *popt_e), ':', color='#2ca02c',
                linewidth=2, alpha=0.7,
                label=f'Exp fit: $R_\\infty$={R_inf:.3f}')
        ax.axhline(R_inf, color='#2ca02c', linestyle=':', alpha=0.3)

    ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
    ax.set_ylabel('$R_{\\rm eq}$ (bohr)', fontsize=14)
    ax.set_title('LiH $R_{\\rm eq}$ vs $l_{\\rm max}$: '
                 'Real Composed Solver vs Diagnostic', fontsize=13)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(range(2, 6))
    fig.tight_layout()
    fig.savefig(plot_dir / 'pk_saturation_real_vs_diag.png', dpi=150)
    plt.close(fig)
    print("Plot A saved to debug/plots/pk_saturation_real_vs_diag.png")

    # --- Plot B: PES curves at each l_max ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for i, lm in enumerate(lm_list):
        res = results[lm]
        ax.plot(res['R'], res['E'], '-o', color=colors[i % 5],
                markersize=4, linewidth=1.5, label=f'$l_{{\\rm max}}={lm}$')

    ax.axvline(R_EXP, color='gray', linestyle='--', alpha=0.5,
               label='Expt $R_{\\rm eq}$')
    ax.set_xlabel('$R$ (bohr)', fontsize=14)
    ax.set_ylabel('$E_{\\rm composed}$ (Ha)', fontsize=14)
    ax.set_title('LiH PES: Real Composed Solver (l-dependent PK)',
                 fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_dir / 'pk_saturation_real_pes.png', dpi=150)
    plt.close(fig)
    print("Plot B saved to debug/plots/pk_saturation_real_pes.png")

    # --- Plot C: Drift increments ---
    if len(lm_list) >= 3:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        drifts = [r_eqs[i] - r_eqs[i-1] for i in range(1, len(r_eqs))]
        drift_lm = lm_list[1:]

        ax1.bar(drift_lm, [abs(d) for d in drifts],
                color='#1f77b4', alpha=0.7, width=0.6)
        ax1.set_xlabel('$l_{\\rm max}$', fontsize=12)
        ax1.set_ylabel('$|\\Delta R_{\\rm eq}|$ (bohr)', fontsize=12)
        ax1.set_title('Drift Increments', fontsize=13)
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(drift_lm)

        # E_min vs l_max
        ax2.plot(lm_list, e_mins, 'o-', color='#1f77b4',
                 linewidth=2, markersize=10)
        ax2.set_xlabel('$l_{\\rm max}$', fontsize=12)
        ax2.set_ylabel('$E_{\\rm min}$ (Ha)', fontsize=12)
        ax2.set_title('Minimum Energy vs $l_{\\rm max}$', fontsize=13)
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(lm_list)

        fig.tight_layout()
        fig.savefig(plot_dir / 'pk_saturation_real_drifts.png', dpi=150)
        plt.close(fig)
        print("Plot C saved to debug/plots/pk_saturation_real_drifts.png")

except ImportError:
    print("matplotlib not available -- skipping plots.")

print("\nDone.")
