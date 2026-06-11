"""
Task 9: PK Saturation Test in the Real Composed Geometry Solver.

Runs the full ComposedDiatomicSolver.LiH_ab_initio() pipeline at
l_max = 2, 3, 4, 5, 6 (and higher if feasible) with l-dependent PK.

Key differences from the simplified diagnostic (pk_saturation_test.py):
  - Radially-varying Z_eff(r) from solved hyperspherical core
  - V_NN uses bare nuclear charges: Z_A * Z_B / R = 3.0 / R
  - V_cross_nuc: core-electron attraction to nucleus B
  - E_core from exact 2-electron hyperspherical solve
  - Variational 2D solver (not just adiabatic approximation)

References:
  Paper 17 Sec III.B (PK pseudopotential), Sec V.A (l_max divergence)
  Paper 15 (Level 4 solver)
"""

import sys
import json
import time
import numpy as np
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.composed_diatomic import ComposedDiatomicSolver

# === Parameters ===
R_EXP = 3.015      # Experimental R_eq (bohr)
R_GRID = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.5, 12),
    np.linspace(5.0, 7.0, 4),
])

# l_max values to test
L_MAX_VALUES = [2, 3, 4, 5, 6]

# Also test channel-blind for comparison at l_max=2,3,4
L_MAX_BLIND = [2, 3, 4]


def find_r_eq_interpolated(R_vals: np.ndarray, E_vals: np.ndarray) -> tuple:
    """Find R_eq by cubic interpolation near the minimum."""
    valid = ~np.isnan(E_vals)
    R_v = R_vals[valid]
    E_v = E_vals[valid]

    i_min = np.argmin(E_v)
    if i_min == 0 or i_min == len(R_v) - 1:
        return float(R_v[i_min]), float(E_v[i_min]), False

    i_lo = max(0, i_min - 2)
    i_hi = min(len(R_v), i_min + 3)
    R_sub = R_v[i_lo:i_hi]
    E_sub = E_v[i_lo:i_hi]

    if len(R_sub) < 3:
        return float(R_v[i_min]), float(E_v[i_min]), True

    cs = CubicSpline(R_sub, E_sub)
    R_fine = np.linspace(R_sub[0], R_sub[-1], 500)
    E_fine = cs(R_fine)
    i_fine = np.argmin(E_fine)
    return float(R_fine[i_fine]), float(E_fine[i_fine]), True


def exp_saturation(x, R_inf, c, beta):
    """Saturating exponential: R_inf - c * exp(-beta * x)."""
    return R_inf - c * np.exp(-beta * np.asarray(x))


def linear_model(x, a, b):
    """Linear: a + b * x."""
    return a + b * np.asarray(x)


def run_composed_pes(l_max: int, pk_channel_mode: str = 'l_dependent',
                     n_alpha: int = 100) -> dict:
    """
    Run the full ComposedDiatomicSolver for LiH at a given l_max.

    Returns dict with R_eq, E_min, PES data, PK parameters, timing.
    """
    t0 = time.time()

    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=l_max,
        pk_channel_mode=pk_channel_mode,
        n_alpha=n_alpha,
        verbose=True,
    )

    # Step 1: Solve core (same for all l_max)
    E_core = solver.solve_core()

    # Step 2: Scan PES
    pes = solver.scan_pes(R_grid=R_GRID, n_Re=300)

    elapsed = time.time() - t0

    # Interpolated R_eq
    R_eq_interp, E_min_interp, success = find_r_eq_interpolated(
        pes['R_valid'], pes['E_valid'])

    return {
        'l_max': l_max,
        'pk_channel_mode': pk_channel_mode,
        'n_alpha': n_alpha,
        'R_eq_grid': float(pes['R_eq']),
        'R_eq_interp': R_eq_interp,
        'E_min': E_min_interp,
        'D_e': float(pes['D_e']),
        'E_core': E_core,
        'pk_A': solver.pk_A,
        'pk_B': solver.pk_B,
        'Z_A_eff': solver.Z_A_eff,
        'Z_A_bare': solver.Z_A_bare,
        'time_s': elapsed,
        'R_eq_found': success,
        'R': pes['R_valid'].tolist(),
        'E_composed': pes['E_valid'].tolist(),
        'E_elec': pes['E_elec'][~np.isnan(pes['E_composed'])].tolist(),
        'V_NN_bare': pes['V_NN_bare'][~np.isnan(pes['E_composed'])].tolist(),
        'V_cross_nuc': pes['V_cross_nuc'][~np.isnan(pes['E_composed'])].tolist(),
    }


def main():
    print("=" * 76)
    print("Task 9: PK Saturation in Real Composed Geometry Solver")
    print("=" * 76)
    print(f"Experimental R_eq = {R_EXP} bohr")
    print(f"l_max sweep: {L_MAX_VALUES}")
    print(f"R grid: {len(R_GRID)} points [{R_GRID[0]:.1f}, {R_GRID[-1]:.1f}] bohr")
    print()

    # ===================================================================
    # Main sweep: l-dependent PK at l_max = 2..6
    # ===================================================================
    results_ldep = {}
    for l_max in L_MAX_VALUES:
        print(f"\n{'='*76}")
        print(f"l_max = {l_max}, pk_channel_mode = l_dependent")
        print(f"{'='*76}")

        try:
            res = run_composed_pes(l_max, pk_channel_mode='l_dependent')
            results_ldep[l_max] = res
            print(f"\n>>> R_eq = {res['R_eq_interp']:.4f} bohr "
                  f"({100*abs(res['R_eq_interp'] - R_EXP)/R_EXP:.1f}% error), "
                  f"E_min = {res['E_min']:.6f} Ha, "
                  f"time = {res['time_s']:.1f}s")
        except Exception as e:
            print(f"\n>>> FAILED at l_max={l_max}: {e}")
            import traceback
            traceback.print_exc()

        # Safety: abort if too slow
        if l_max in results_ldep and results_ldep[l_max]['time_s'] > 600:
            print(f"\n** l_max={l_max} took {results_ldep[l_max]['time_s']:.0f}s "
                  f"(>600s). Stopping sweep. **")
            break

    # ===================================================================
    # Comparison: channel-blind PK at l_max = 2, 3, 4
    # ===================================================================
    results_blind = {}
    for l_max in L_MAX_BLIND:
        if l_max > max(results_ldep.keys(), default=2):
            break
        print(f"\n{'='*76}")
        print(f"l_max = {l_max}, pk_channel_mode = channel_blind")
        print(f"{'='*76}")

        try:
            res = run_composed_pes(l_max, pk_channel_mode='channel_blind')
            results_blind[l_max] = res
            print(f"\n>>> R_eq = {res['R_eq_interp']:.4f} bohr "
                  f"({100*abs(res['R_eq_interp'] - R_EXP)/R_EXP:.1f}% error), "
                  f"time = {res['time_s']:.1f}s")
        except Exception as e:
            print(f"\n>>> channel_blind FAILED at l_max={l_max}: {e}")

    # ===================================================================
    # Analysis
    # ===================================================================
    print("\n\n" + "=" * 76)
    print("RESULTS TABLE: l-dependent PK")
    print("=" * 76)

    lm_list = sorted(results_ldep.keys())
    r_eqs = [results_ldep[lm]['R_eq_interp'] for lm in lm_list]
    e_mins = [results_ldep[lm]['E_min'] for lm in lm_list]

    print(f"{'l_max':>5}  {'R_eq':>10}  {'dR_eq':>8}  {'%err':>7}  "
          f"{'E_min':>12}  {'D_e':>10}  {'Time':>8}")
    print("-" * 76)

    for i, lm in enumerate(lm_list):
        res = results_ldep[lm]
        dr = r_eqs[i] - r_eqs[i-1] if i > 0 else 0.0
        err = 100 * abs(r_eqs[i] - R_EXP) / R_EXP
        print(f"{lm:>5}  {r_eqs[i]:>10.4f}  {dr:>+8.4f}  {err:>7.1f}%  "
              f"{e_mins[i]:>12.6f}  {res['D_e']:>10.6f}  "
              f"{res['time_s']:>8.1f}s")

    # Channel-blind comparison
    if results_blind:
        print(f"\n{'='*76}")
        print("COMPARISON: channel-blind PK")
        print(f"{'='*76}")
        print(f"{'l_max':>5}  {'R_eq(blind)':>12}  {'R_eq(l-dep)':>12}  "
              f"{'diff':>8}  {'blind %err':>10}")
        for lm in sorted(results_blind.keys()):
            rb = results_blind[lm]['R_eq_interp']
            rl = results_ldep.get(lm, {}).get('R_eq_interp', float('nan'))
            err_b = 100 * abs(rb - R_EXP) / R_EXP
            print(f"{lm:>5}  {rb:>12.4f}  {rl:>12.4f}  "
                  f"{rb - rl:>+8.4f}  {err_b:>10.1f}%")

    # ===================================================================
    # Fit analysis
    # ===================================================================
    if len(lm_list) >= 3:
        print(f"\n{'='*76}")
        print("FIT ANALYSIS")
        print(f"{'='*76}")

        x = np.array(lm_list, dtype=float)
        y = np.array(r_eqs)
        ss_tot = np.sum((y - np.mean(y))**2)

        # Linear
        try:
            popt_lin, _ = curve_fit(linear_model, x, y)
            ss_res = np.sum((y - linear_model(x, *popt_lin))**2)
            r2_lin = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"Linear: R_eq = {popt_lin[0]:.4f} + {popt_lin[1]:.4f} * l_max"
                  f"  (R²={r2_lin:.6f})")
        except Exception:
            popt_lin = None
            r2_lin = -1

        # Saturating exponential
        try:
            p0 = [r_eqs[-1] + 0.3, r_eqs[-1] - r_eqs[0] + 0.3, 0.4]
            popt_exp, _ = curve_fit(exp_saturation, x, y, p0=p0, maxfev=10000)
            ss_res = np.sum((y - exp_saturation(x, *popt_exp))**2)
            r2_exp = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            R_inf = popt_exp[0]
            print(f"Exp. saturation: R_inf = {R_inf:.4f} bohr  "
                  f"({100*abs(R_inf - R_EXP)/R_EXP:.1f}% from expt)"
                  f"  (R²={r2_exp:.6f})")
        except Exception as e:
            print(f"Exp. fit failed: {e}")
            popt_exp = None
            r2_exp = -1
            R_inf = None

        # Drift increments
        drifts = [r_eqs[i] - r_eqs[i-1] for i in range(1, len(r_eqs))]
        print(f"\nDrift increments:")
        for i, dr in enumerate(drifts):
            print(f"  l_max {lm_list[i]}->{lm_list[i+1]}: dR = {dr:+.4f} bohr")

    # ===================================================================
    # Diagnostic comparison table
    # ===================================================================
    # Values from pk_saturation_test.py (simplified model)
    diagnostic_r_eqs = {2: 3.71, 3: 4.13, 4: 4.33, 5: 4.45,
                        6: 4.50, 7: 4.54, 8: 4.56}
    print(f"\n{'='*76}")
    print("COMPARISON: Real Solver vs Diagnostic Model")
    print(f"{'='*76}")
    print(f"{'l_max':>5}  {'Real R_eq':>10}  {'Diag R_eq':>10}  "
          f"{'Diff':>8}  {'Real %err':>10}  {'Diag %err':>10}")
    for lm in lm_list:
        rr = results_ldep[lm]['R_eq_interp']
        rd = diagnostic_r_eqs.get(lm, float('nan'))
        err_r = 100 * abs(rr - R_EXP) / R_EXP
        err_d = 100 * abs(rd - R_EXP) / R_EXP if not np.isnan(rd) else float('nan')
        print(f"{lm:>5}  {rr:>10.4f}  {rd:>10.2f}  "
              f"{rr - rd:>+8.4f}  {err_r:>10.1f}%  {err_d:>10.1f}%")

    # ===================================================================
    # V_NN verification
    # ===================================================================
    print(f"\n{'='*76}")
    print("V_NN VERIFICATION")
    print(f"{'='*76}")
    if lm_list:
        res0 = results_ldep[lm_list[0]]
        print(f"Z_A_bare = {res0['Z_A_bare']}")
        print(f"Z_A_eff  = {res0['Z_A_eff']}")
        print(f"V_NN uses bare charges: Z_A_bare * Z_B / R = "
              f"{res0['Z_A_bare']} * 1.0 / R = {res0['Z_A_bare']}/R")
        print(f"Diagnostic used: Z_A_eff * Z_B / R = 2.69/R (WRONG)")
        print(f"Ratio: {res0['Z_A_bare']}/2.69 = {res0['Z_A_bare']/2.69:.3f}")
        print(f"V_NN at R=3.015: real = {res0['Z_A_bare']/3.015:.4f} Ha, "
              f"diagnostic = {2.69/3.015:.4f} Ha, "
              f"diff = {(res0['Z_A_bare'] - 2.69)/3.015:.4f} Ha")
        # V_cross_nuc
        if res0['V_cross_nuc']:
            idx_3 = min(range(len(res0['R'])),
                        key=lambda i: abs(res0['R'][i] - 3.015))
            print(f"V_cross_nuc at R~3.015: {res0['V_cross_nuc'][idx_3]:.4f} Ha")
            print(f"  (core attraction to nucleus B — not in diagnostic)")

    # ===================================================================
    # Bottom line
    # ===================================================================
    print(f"\n{'='*76}")
    print("BOTTOM LINE")
    print(f"{'='*76}")

    if len(lm_list) >= 3:
        last_drift = abs(r_eqs[-1] - r_eqs[-2]) if len(r_eqs) >= 2 else float('inf')
        first_drift = abs(r_eqs[1] - r_eqs[0]) if len(r_eqs) >= 2 else float('inf')
        diminishing = last_drift < first_drift * 0.7

        if popt_exp is not None and r2_exp > r2_lin and r2_exp > 0.95:
            print(f"R_eq SATURATES to R_inf = {R_inf:.3f} bohr "
                  f"({100*abs(R_inf - R_EXP)/R_EXP:.1f}% from experiment).")
        elif diminishing:
            print(f"Drift increments diminishing ({first_drift:.4f} -> "
                  f"{last_drift:.4f}), suggesting saturation.")
            print(f"Last computed R_eq = {r_eqs[-1]:.4f} bohr at l_max={lm_list[-1]}")
        else:
            print(f"R_eq trend unclear. Drifts: {first_drift:.4f} -> {last_drift:.4f}")
            print(f"Last R_eq = {r_eqs[-1]:.4f} at l_max={lm_list[-1]}")

        err_last = 100 * abs(r_eqs[-1] - R_EXP) / R_EXP
        print(f"Current best R_eq error: {err_last:.1f}% at l_max={lm_list[-1]}")

    # ===================================================================
    # Save data
    # ===================================================================
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)

    output = {
        'task': 'Task 9: PK Saturation in Real Composed Solver',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'R_exp': R_EXP,
        'l_dependent_results': {str(k): v for k, v in results_ldep.items()},
        'channel_blind_results': {str(k): v for k, v in results_blind.items()},
        'diagnostic_r_eqs': {str(k): v for k, v in diagnostic_r_eqs.items()},
        'summary': {
            'l_max_values': lm_list,
            'R_eq_values': r_eqs if lm_list else [],
            'E_min_values': e_mins if lm_list else [],
        },
    }

    if len(lm_list) >= 3:
        output['fits'] = {
            'linear': {
                'params': popt_lin.tolist() if popt_lin is not None else None,
                'R2': r2_lin,
            },
            'exponential': {
                'params': popt_exp.tolist() if popt_exp is not None else None,
                'R2': r2_exp,
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

        if lm_list:
            # --- Plot A: R_eq vs l_max (real + diagnostic + experiment) ---
            fig, ax = plt.subplots(1, 1, figsize=(10, 7))

            ax.plot(lm_list, r_eqs, 'o-', color='#1f77b4', linewidth=2.5,
                    markersize=10, label='Real solver (l-dependent PK)', zorder=5)

            # Diagnostic values
            diag_lm = sorted([lm for lm in diagnostic_r_eqs if lm <= max(lm_list) + 2])
            diag_r = [diagnostic_r_eqs[lm] for lm in diag_lm]
            ax.plot(diag_lm, diag_r, 's--', color='#ff7f0e', linewidth=2,
                    markersize=8, alpha=0.7, label='Diagnostic model (const Z_eff)')

            # Channel-blind
            if results_blind:
                bl_lm = sorted(results_blind.keys())
                bl_r = [results_blind[lm]['R_eq_interp'] for lm in bl_lm]
                ax.plot(bl_lm, bl_r, 'D-.', color='#d62728', linewidth=2,
                        markersize=8, alpha=0.7, label='Channel-blind PK')

            ax.axhline(R_EXP, color='gray', linestyle='--', linewidth=1.5,
                        alpha=0.7, label=f'Experiment ({R_EXP} bohr)')

            # Saturation fit
            if popt_exp is not None and r2_exp > 0.9:
                x_fine = np.linspace(1.5, max(lm_list) + 2, 100)
                ax.plot(x_fine, exp_saturation(x_fine, *popt_exp), ':',
                        color='#2ca02c', linewidth=2, alpha=0.7,
                        label=f'Exp. fit: R$_\\infty$={R_inf:.3f} bohr')
                ax.axhline(R_inf, color='#2ca02c', linestyle=':', alpha=0.3)

            ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
            ax.set_ylabel('$R_{\\rm eq}$ (bohr)', fontsize=14)
            ax.set_title('LiH $R_{\\rm eq}$ vs $l_{\\rm max}$ — '
                         'Real Composed Solver vs Diagnostic', fontsize=13)
            ax.legend(fontsize=10, loc='best')
            ax.grid(True, alpha=0.3)
            ax.set_xticks(range(2, max(max(lm_list), max(diag_lm)) + 1))
            fig.tight_layout()
            fig.savefig(plot_dir / 'pk_saturation_real_vs_diag.png', dpi=150)
            plt.close(fig)
            print(f"Plot A saved to debug/plots/pk_saturation_real_vs_diag.png")

            # --- Plot B: PES curves at each l_max ---
            fig, ax = plt.subplots(1, 1, figsize=(10, 7))
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                      '#8c564b', '#e377c2']
            for i, lm in enumerate(lm_list):
                res = results_ldep[lm]
                ax.plot(res['R'], res['E_composed'], '-o',
                        color=colors[i % len(colors)], markersize=4,
                        linewidth=1.5, label=f'$l_{{\\rm max}}={lm}$')

            ax.axvline(R_EXP, color='gray', linestyle='--', alpha=0.5,
                       label=f'Expt R_eq')
            ax.set_xlabel('$R$ (bohr)', fontsize=14)
            ax.set_ylabel('$E_{\\rm composed}$ (Ha)', fontsize=14)
            ax.set_title('LiH PES: Real Composed Solver (l-dependent PK)',
                         fontsize=13)
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(plot_dir / 'pk_saturation_real_pes.png', dpi=150)
            plt.close(fig)
            print(f"Plot B saved to debug/plots/pk_saturation_real_pes.png")

            # --- Plot C: Energy decomposition at R_eq ---
            if len(lm_list) >= 2:
                fig, axes = plt.subplots(1, 2, figsize=(14, 6))

                # Left: E_elec at ~R_eq for each l_max
                idx_eq = [min(range(len(results_ldep[lm]['R'])),
                              key=lambda i: abs(results_ldep[lm]['R'][i] - 3.0))
                          for lm in lm_list]
                e_elec_eq = [results_ldep[lm]['E_elec'][idx_eq[i]]
                             for i, lm in enumerate(lm_list)]
                v_cross_eq = [results_ldep[lm]['V_cross_nuc'][idx_eq[i]]
                              for i, lm in enumerate(lm_list)]

                axes[0].plot(lm_list, e_elec_eq, 'o-', color='#1f77b4',
                             label='$E_{\\rm elec}$')
                axes[0].plot(lm_list, v_cross_eq, 's-', color='#d62728',
                             label='$V_{\\rm cross}$')
                axes[0].set_xlabel('$l_{\\rm max}$')
                axes[0].set_ylabel('Energy (Ha)')
                axes[0].set_title(f'Energy components at R$\\approx$3.0 bohr')
                axes[0].legend()
                axes[0].grid(True, alpha=0.3)

                # Right: drift increments
                drifts = [r_eqs[i] - r_eqs[i-1] for i in range(1, len(r_eqs))]
                drift_lm = lm_list[1:]
                if drifts:
                    axes[1].bar(drift_lm, [abs(d) for d in drifts],
                                color='#1f77b4', alpha=0.7)
                    axes[1].set_xlabel('$l_{\\rm max}$')
                    axes[1].set_ylabel('$|\\Delta R_{\\rm eq}|$ (bohr)')
                    axes[1].set_title('Drift Increments')
                    axes[1].grid(True, alpha=0.3)

                fig.tight_layout()
                fig.savefig(plot_dir / 'pk_saturation_real_decomp.png', dpi=150)
                plt.close(fig)
                print(f"Plot C saved to debug/plots/pk_saturation_real_decomp.png")

    except ImportError:
        print("matplotlib not available -- skipping plots.")

    print("\nDone.")


if __name__ == '__main__':
    main()
