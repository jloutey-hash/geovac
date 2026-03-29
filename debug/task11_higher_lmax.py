"""
Task 11 Part A: R-Dependent PK at Higher l_max.

Tests whether the R-dependent PK improvement (Task 10) holds at l_max = 5, 6, 7.
Also re-runs l-dependent PK at l_max = 5, 6 for consistent comparison.

Uses the same angular eigenvalue infrastructure as Task 10 (selfconsistent_pk.py).

Key parameters (matching Task 10):
  - R-dependent PK: w(R) = delta_{l,0} * min(cap, R/R_ref), cap=1.5, R_ref=3.015
  - R grid: same 9 points as Task 10
  - Adiabatic approximation: min_{R_e} U(R_e)
"""

import sys
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.level4_multichannel import (
    _channel_list,
    build_angular_hamiltonian,
    compute_pk_pseudopotential,
)

# === Physical parameters (match Tasks 9–10) ===
Z_A_EFF = 1.0
Z_B = 1.0
Z_A_BARE = 3.0
PK_A = 6.93
PK_B_PARAM = 6.80
R_EXP = 3.015

# Grid parameters
N_ALPHA = 100
N_RE = 25
R_E_GRID = np.logspace(np.log10(0.3), np.log10(5.0), N_RE)

# PES scan grid (same as Task 10)
R_VALUES = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5])


def get_channels(l_max: int) -> list:
    return _channel_list(l_max, homonuclear=(Z_A_EFF == Z_B))


def solve_angular_ldependent(
    R: float, R_e: float, l_max: int, n_alpha: int = N_ALPHA,
) -> float:
    """Standard l-dependent PK (w=1.0 for l=0). Returns mu0."""
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A_EFF - Z_B) / (2.0 * (Z_A_EFF + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_PARAM,
        'atom': 'A',
        'channel_mode': 'l_dependent',
    }]

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A_EFF, Z_B=Z_B, z0=z0,
        pk_potentials=pk_pots,
    )
    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def solve_angular_r_dependent_pk(
    R: float, R_e: float, l_max: int, R_ref: float = R_EXP,
    pk_cap: float = 1.5, n_alpha: int = N_ALPHA,
) -> float:
    """
    R-dependent PK: w(R) = delta_{l,0} * min(cap, R/R_ref).
    """
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A_EFF - Z_B) / (2.0 * (Z_A_EFF + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    channels = get_channels(l_max)
    n_ch = len(channels)

    rho_A = (R / 2.0 - z0) / R_e
    rho_B = (R / 2.0 + z0) / R_e
    pk_pots_raw = [{'C_core': PK_A, 'beta_core': PK_B_PARAM, 'atom': 'A'}]
    V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
        alpha, rho_A, R_e, pk_pots_raw, rho_B=rho_B,
    )

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A_EFF, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )

    scale = min(pk_cap, max(1.0, R / R_ref))

    for ic, (l1, l2) in enumerate(channels):
        w1 = scale if l1 == 0 else 0.0
        w2 = scale if l2 == 0 else 0.0
        if w1 == 0.0 and w2 == 0.0:
            continue
        for i in range(n_alpha):
            ii = ic * n_alpha + i
            H[ii, ii] += R_e * w1 * V_pk_e1[i]
            H[ii, ii] += R_e * w2 * V_pk_e2[i]

    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def compute_pes(
    l_max: int,
    R_grid: np.ndarray,
    Re_grid: np.ndarray,
    mode: str = 'l_dependent',
    R_ref: float = R_EXP,
    pk_cap: float = 1.5,
) -> dict:
    """Compute PES curve with adiabatic approximation."""
    n_R = len(R_grid)
    E_elec = np.zeros(n_R)

    for ir, R in enumerate(R_grid):
        U_vals = np.zeros(len(Re_grid))
        for ire, R_e in enumerate(Re_grid):
            if mode == 'l_dependent':
                mu0 = solve_angular_ldependent(R, R_e, l_max)
            elif mode == 'r_dependent':
                mu0 = solve_angular_r_dependent_pk(
                    R, R_e, l_max, R_ref=R_ref, pk_cap=pk_cap)
            else:
                raise ValueError(f"Unknown mode: {mode}")
            U_vals[ire] = (mu0 + 15.0 / 8.0) / (R_e ** 2)
        E_elec[ir] = np.min(U_vals)

    V_nn = Z_A_EFF * Z_B / R_grid
    E_total = E_elec + V_nn

    return {
        'R': R_grid.tolist(),
        'E_total': E_total.tolist(),
        'E_elec': E_elec.tolist(),
    }


def find_r_eq(R_vals: np.ndarray, E_vals: np.ndarray) -> tuple:
    """Find R_eq by cubic interpolation near the minimum."""
    i_min = np.argmin(E_vals)
    if i_min == 0 or i_min == len(R_vals) - 1:
        return float(R_vals[i_min]), float(E_vals[i_min]), False

    i_lo = max(0, i_min - 2)
    i_hi = min(len(R_vals), i_min + 3)
    R_sub = R_vals[i_lo:i_hi]
    E_sub = E_vals[i_lo:i_hi]

    if len(R_sub) < 3:
        return float(R_vals[i_min]), float(E_vals[i_min]), True

    cs = CubicSpline(R_sub, E_sub)
    R_fine = np.linspace(R_sub[0], R_sub[-1], 500)
    E_fine = cs(R_fine)
    i_fine_min = np.argmin(E_fine)
    return float(R_fine[i_fine_min]), float(E_fine[i_fine_min]), True


def main():
    print("=" * 76)
    print("Task 11 Part A: R-Dependent PK at Higher l_max")
    print("=" * 76)
    print(f"R-dependent PK: w(R) = delta_{{l,0}} * min(1.5, R/{R_EXP})")
    print(f"Experimental R_eq = {R_EXP} bohr")
    print()

    # Collect all results from Tasks 9-10 and new computations
    # Task 10 data (l_max 2-4 for l-dependent, self-consistent, r-dependent)
    all_results = {
        'l_dependent': {},
        'r_dependent': {},
        'selfconsistent': {},
    }

    # Pre-fill Task 10 results
    task10_ldep = {
        2: {'R_eq': 2.803, 'E_min': -1.220212, 'time': 13.9},
        3: {'R_eq': 3.015, 'E_min': -1.272553, 'time': 30.2},
        4: {'R_eq': 3.423, 'E_min': -1.342805, 'time': 65.0},
    }
    task10_rdep = {
        2: {'R_eq': 2.837, 'E_min': -1.220669, 'time': 13.8},
        3: {'R_eq': 2.975, 'E_min': -1.272719, 'time': 30.1},
        4: {'R_eq': 3.075, 'E_min': -1.334423, 'time': 64.1},
    }
    task10_sc = {
        2: {'R_eq': 2.878, 'E_min': -1.256245, 'time': 173.9},
        3: {'R_eq': 3.026, 'E_min': -1.318572, 'time': 389.7},
        4: {'R_eq': 3.192, 'E_min': -1.387480, 'time': 1027.0},
    }

    all_results['l_dependent'] = dict(task10_ldep)
    all_results['r_dependent'] = dict(task10_rdep)
    all_results['selfconsistent'] = dict(task10_sc)

    # l_max values to compute fresh
    l_max_new = [5, 6, 7]
    max_time_per_point = 1800  # abort if single l_max takes > 30 min

    # ===================================================================
    # Run l-dependent PK at l_max = 5, 6, 7
    # ===================================================================
    for l_max in l_max_new:
        n_ch = len(get_channels(l_max))
        print(f"\n{'='*76}")
        print(f"l-dependent PK: l_max = {l_max} ({n_ch} channels)")
        print(f"{'='*76}")

        t0 = time.time()
        try:
            pes = compute_pes(l_max, R_VALUES, R_E_GRID, mode='l_dependent')
            elapsed = time.time() - t0
            R_eq, E_min, ok = find_r_eq(
                np.array(pes['R']), np.array(pes['E_total']))
            err = 100 * abs(R_eq - R_EXP) / R_EXP

            all_results['l_dependent'][l_max] = {
                'R_eq': R_eq, 'E_min': E_min, 'time': elapsed,
                'pes': pes,
            }
            print(f"  R_eq = {R_eq:.4f} bohr ({err:.1f}% error)")
            print(f"  E_min = {E_min:.6f} Ha, time = {elapsed:.1f}s")

            if elapsed > max_time_per_point:
                print(f"\n  ** Took {elapsed:.0f}s (> {max_time_per_point}s). "
                      f"Stopping l-dependent sweep. **")
                break
        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()
            break

    # ===================================================================
    # Run R-dependent PK at l_max = 5, 6, 7
    # ===================================================================
    for l_max in l_max_new:
        # Only run if we managed l-dependent at this l_max (similar cost)
        if l_max not in all_results['l_dependent']:
            print(f"\n  Skipping R-dependent at l_max={l_max} "
                  f"(l-dependent didn't complete)")
            break

        n_ch = len(get_channels(l_max))
        print(f"\n{'='*76}")
        print(f"R-dependent PK (cap=1.5): l_max = {l_max} ({n_ch} channels)")
        print(f"{'='*76}")

        t0 = time.time()
        try:
            pes = compute_pes(l_max, R_VALUES, R_E_GRID, mode='r_dependent',
                              pk_cap=1.5)
            elapsed = time.time() - t0
            R_eq, E_min, ok = find_r_eq(
                np.array(pes['R']), np.array(pes['E_total']))
            err = 100 * abs(R_eq - R_EXP) / R_EXP

            all_results['r_dependent'][l_max] = {
                'R_eq': R_eq, 'E_min': E_min, 'time': elapsed,
                'pes': pes,
            }
            print(f"  R_eq = {R_eq:.4f} bohr ({err:.1f}% error)")
            print(f"  E_min = {E_min:.6f} Ha, time = {elapsed:.1f}s")

            if elapsed > max_time_per_point:
                print(f"\n  ** Took {elapsed:.0f}s. Stopping R-dependent sweep. **")
                break
        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()
            break

    # ===================================================================
    # Summary Table
    # ===================================================================
    print("\n\n" + "=" * 76)
    print("COMPREHENSIVE SUMMARY TABLE")
    print("=" * 76)
    print(f"{'Mode':<25s} {'l_max':>5s} {'n_ch':>5s} {'R_eq':>8s} "
          f"{'dR_eq':>8s} {'% err':>7s} {'E_min':>12s} {'Time':>8s}")
    print("-" * 80)

    for mode_name, mode_key in [('l-dependent', 'l_dependent'),
                                 ('R-dependent (cap=1.5)', 'r_dependent'),
                                 ('self-consistent', 'selfconsistent')]:
        lm_list = sorted(all_results[mode_key].keys())
        prev_req = None
        for lm in lm_list:
            data = all_results[mode_key][lm]
            R_eq = data['R_eq']
            dR = R_eq - prev_req if prev_req is not None else 0.0
            err = 100 * abs(R_eq - R_EXP) / R_EXP
            n_ch = len(get_channels(lm))
            t = data.get('time', 0)
            print(f"{mode_name:<25s} {lm:>5d} {n_ch:>5d} {R_eq:>8.3f} "
                  f"{dR:>+8.3f} {err:>7.1f} {data['E_min']:>12.6f} "
                  f"{t:>7.1f}s")
            prev_req = R_eq
        print()

    print(f"{'Experiment':<25s} {'':>5s} {'':>5s} {R_EXP:>8.3f}")

    # ===================================================================
    # Drift rate analysis
    # ===================================================================
    print("\n" + "=" * 76)
    print("DRIFT RATE ANALYSIS")
    print("=" * 76)

    for mode_name, mode_key in [('l-dependent', 'l_dependent'),
                                 ('R-dependent (cap=1.5)', 'r_dependent'),
                                 ('self-consistent', 'selfconsistent')]:
        lm_list = sorted(all_results[mode_key].keys())
        if len(lm_list) >= 2:
            reqs = [all_results[mode_key][lm]['R_eq'] for lm in lm_list]
            slope = np.polyfit(lm_list, reqs, 1)[0]
            drifts = [reqs[i] - reqs[i-1] for i in range(1, len(reqs))]
            print(f"\n  {mode_name}:")
            print(f"    Overall slope: {slope:+.4f} bohr/l_max")
            print(f"    Drift increments:", end="")
            for i, dr in enumerate(drifts):
                print(f"  {lm_list[i]}->{lm_list[i+1]}: {dr:+.4f}", end="")
            print()

    # ===================================================================
    # Save results
    # ===================================================================
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)

    # Strip PES data for JSON (keep it compact)
    save_data = {}
    for mode_key in ['l_dependent', 'r_dependent', 'selfconsistent']:
        save_data[mode_key] = {}
        for lm, data in all_results[mode_key].items():
            entry = {k: v for k, v in data.items() if k != 'pes'}
            save_data[mode_key][str(lm)] = entry
            if 'pes' in data:
                entry['pes'] = data['pes']

    output = {
        'task': 'Task 11 Part A: R-Dependent PK at Higher l_max',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'R_exp': R_EXP,
        'pk_cap': 1.5,
        'R_ref': R_EXP,
        'results': save_data,
    }

    json_path = data_dir / 'task11_higher_lmax.json'
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

        # --- Plot 1: R_eq vs l_max for all modes ---
        fig, ax = plt.subplots(figsize=(10, 7))

        mode_styles = {
            'l_dependent': ('o-', '#1f77b4', 'l-dependent'),
            'r_dependent': ('D--', '#2ca02c', 'R-dependent (cap=1.5)'),
            'selfconsistent': ('s-.', '#d62728', 'self-consistent'),
        }

        for mode_key, (style, color, label) in mode_styles.items():
            lm_list = sorted(all_results[mode_key].keys())
            reqs = [all_results[mode_key][lm]['R_eq'] for lm in lm_list]
            ax.plot(lm_list, reqs, style, color=color, label=label,
                    linewidth=2, markersize=8)

        ax.axhline(R_EXP, color='gray', linestyle=':', linewidth=2,
                    alpha=0.7, label=f'Experiment ({R_EXP} bohr)')
        ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
        ax.set_ylabel('$R_{\\rm eq}$ (bohr)', fontsize=14)
        ax.set_title('Task 11: LiH $R_{\\rm eq}$ vs $l_{\\rm max}$ — '
                     'All PK Modes', fontsize=13)
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        max_lm = max(max(all_results[k].keys()) for k in all_results if all_results[k])
        ax.set_xticks(range(2, max_lm + 1))
        fig.tight_layout()
        fig.savefig(plot_dir / 'task11_req_vs_lmax.png', dpi=150)
        plt.close(fig)
        print(f"Plot 1 saved: debug/plots/task11_req_vs_lmax.png")

        # --- Plot 2: Drift increments (log scale) ---
        fig, ax = plt.subplots(figsize=(10, 7))

        for mode_key, (style, color, label) in mode_styles.items():
            lm_list = sorted(all_results[mode_key].keys())
            reqs = [all_results[mode_key][lm]['R_eq'] for lm in lm_list]
            if len(reqs) >= 2:
                drifts = [abs(reqs[i] - reqs[i-1]) for i in range(1, len(reqs))]
                drift_lm = lm_list[1:]
                ax.semilogy(drift_lm, drifts, style.replace('-', '-'),
                            color=color, label=label,
                            linewidth=2, markersize=8)

        ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
        ax.set_ylabel('$|\\Delta R_{\\rm eq}|$ (bohr)', fontsize=14)
        ax.set_title('Task 11: Drift Increments $|\\Delta R_{\\rm eq}|$ '
                     'vs $l_{\\rm max}$', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, which='both')
        fig.tight_layout()
        fig.savefig(plot_dir / 'task11_drift_increments.png', dpi=150)
        plt.close(fig)
        print(f"Plot 2 saved: debug/plots/task11_drift_increments.png")

        # --- Plot 3: PES curves for R-dependent PK ---
        fig, ax = plt.subplots(figsize=(10, 7))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                  '#8c564b']
        lm_plot = [lm for lm in sorted(all_results['r_dependent'].keys())
                   if 'pes' in all_results['r_dependent'][lm]]
        # Also include l_max 2, 4 even without pes data
        for i, lm in enumerate(sorted(all_results['r_dependent'].keys())):
            data = all_results['r_dependent'][lm]
            if 'pes' in data:
                ax.plot(data['pes']['R'], data['pes']['E_total'],
                        '-o', color=colors[i % len(colors)], markersize=4,
                        linewidth=1.5, label=f'$l_{{\\rm max}}={lm}$')

        ax.axvline(R_EXP, color='gray', linestyle='--', alpha=0.5,
                   label=f'Expt $R_{{eq}}$')
        ax.set_xlabel('$R$ (bohr)', fontsize=14)
        ax.set_ylabel('$E_{\\rm total}$ (Ha)', fontsize=14)
        ax.set_title('LiH PES: R-Dependent PK (cap=1.5)', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plot_dir / 'task11_pes_rdependent.png', dpi=150)
        plt.close(fig)
        print(f"Plot 3 saved: debug/plots/task11_pes_rdependent.png")

    except ImportError:
        print("matplotlib not available — skipping plots.")

    print("\nDone.")


if __name__ == '__main__':
    main()
