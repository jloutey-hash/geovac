"""
PK Saturation Test: R_eq Convergence at Higher l_max.

Task 8: Extend l-dependent PK PES sweep to l_max = 2..8 and determine
whether R_eq converges (saturates) or keeps diverging.

Method: Adiabatic approximation (same as Task 7 pes_sweep.py).
  E_elec(R) = min_{R_e} U(R_e)  where  U(R_e) = [mu_0(R_e;R) + 15/8] / R_e^2
  E_total(R) = E_elec(R) + Z_A_eff * Z_B / R

Only l-dependent PK mode (delta_{l_i,0} per electron) is tested -- established
as best mode in Tasks 6-7.

References:
  Paper 17 Sec III.B (PK pseudopotential), Sec V.A (l_max divergence)
  Paper 15 Sec IV (adiabatic separation), Sec V (multichannel expansion)
  Tasks 1-7 diagnostic arc (asymmetric_bond_diagnostic_summary.md)
"""

import sys
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.level4_multichannel import (
    _channel_list,
    build_angular_hamiltonian,
)

# === Physical parameters ===
Z_A = 2.69            # Clementi-Raimondi Z_eff for Li
Z_B = 1.0             # H
PK_A = 6.93           # Ha*bohr^2 (PK amplitude)
PK_B_param = 7.00     # bohr^-2 (PK Gaussian exponent)
R_EXP = 3.015         # Experimental R_eq (bohr)

# l_max sweep range
L_MAX_VALUES = [2, 3, 4, 5, 6, 7, 8]

# PES scan grid: 15 points from 2.0 to 8.0 bohr
N_R = 15
R_VALUES = np.linspace(2.0, 8.0, N_R)

# Hyperradial grid: 25 points, log-spaced 0.3 to 5.0 bohr
N_RE = 25
R_E_GRID = np.logspace(np.log10(0.3), np.log10(5.0), N_RE)

# Angular grid
N_ALPHA = 50

# Diagnostic R values for decomposition (Plot D)
R_SMALL = 3.71   # R_eq at l_max=2
R_LARGE = 5.0    # far from equilibrium

# Fallback grid for cost management
FALLBACK_N_R = 10
FALLBACK_R_VALUES = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0])
FALLBACK_N_RE = 20
FALLBACK_R_E_GRID = np.logspace(np.log10(0.3), np.log10(5.0), FALLBACK_N_RE)
FALLBACK_N_ALPHA = 40

# Timing threshold to trigger fallback (seconds per PES curve)
FALLBACK_THRESHOLD = 120.0


def is_homonuclear() -> bool:
    return Z_A == Z_B


def get_channels(l_max: int) -> list:
    return _channel_list(l_max, homonuclear=is_homonuclear())


def solve_angular_mu0(R: float, R_e: float, l_max: int,
                      n_alpha: int) -> float:
    """
    Solve angular eigenvalue problem at (R, R_e) with l-dependent PK.
    Returns mu_0 (lowest angular eigenvalue).
    """
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_param,
        'atom': 'A',
        'channel_mode': 'l_dependent',
    }]

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=pk_pots,
    )
    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def compute_pes_curve(l_max: int, R_grid: np.ndarray, Re_grid: np.ndarray,
                      n_alpha: int) -> tuple:
    """
    Compute E_total(R) for a given l_max with l-dependent PK.
    Returns (R_values, E_total, E_elec, V_nn).
    """
    n_R = len(R_grid)
    E_elec = np.zeros(n_R)
    V_nn = np.zeros(n_R)

    for ir, R in enumerate(R_grid):
        U_vals = np.zeros(len(Re_grid))
        for ire, R_e in enumerate(Re_grid):
            mu0 = solve_angular_mu0(R, R_e, l_max, n_alpha)
            U_vals[ire] = (mu0 + 15.0 / 8.0) / (R_e ** 2)

        E_elec[ir] = np.min(U_vals)
        V_nn[ir] = Z_A * Z_B / R

    E_total = E_elec + V_nn
    return R_grid.copy(), E_total, E_elec, V_nn


def find_r_eq(R_vals: np.ndarray, E_vals: np.ndarray) -> tuple:
    """
    Find R_eq by cubic interpolation near the minimum.
    Returns (R_eq, E_min, success).
    """
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


def compute_e_elec_at_R(R: float, l_max: int, Re_grid: np.ndarray,
                        n_alpha: int) -> float:
    """Compute E_elec at a single R value (for decomposition analysis)."""
    U_vals = np.zeros(len(Re_grid))
    for ire, R_e in enumerate(Re_grid):
        mu0 = solve_angular_mu0(R, R_e, l_max, n_alpha)
        U_vals[ire] = (mu0 + 15.0 / 8.0) / (R_e ** 2)
    return float(np.min(U_vals))


def exp_saturation(x, R_inf, c, beta):
    """Saturating exponential: R_inf - c * exp(-beta * x)."""
    return R_inf - c * np.exp(-beta * np.asarray(x))


def linear_model(x, a, b):
    """Linear: a + b * x."""
    return a + b * np.asarray(x)


def exp_decay(x, A, gamma):
    """Exponential decay for drift increments: A * exp(-gamma * x)."""
    return A * np.exp(-gamma * np.asarray(x))


def main():
    print("=" * 76)
    print("PK Saturation Test: R_eq Convergence at Higher l_max")
    print("=" * 76)
    print(f"Z_A_eff = {Z_A}, Z_B = {Z_B}")
    print(f"PK: A = {PK_A} Ha*bohr^2, B = {PK_B_param} bohr^-2")
    print(f"PK mode: l-dependent (delta_{{l_i,0}} per electron)")
    print(f"Experimental R_eq = {R_EXP} bohr")
    print()

    # Current grid settings
    r_grid = R_VALUES.copy()
    re_grid = R_E_GRID.copy()
    n_alpha = N_ALPHA
    use_fallback = False

    # ===================================================================
    # Timing test: run l_max=2 to estimate total compute time
    # ===================================================================
    print("Timing calibration (l_max=2)...", flush=True)
    t0 = time.time()
    R_test, E_test, _, _ = compute_pes_curve(2, r_grid, re_grid, n_alpha)
    t_lmax2 = time.time() - t0
    print(f"  l_max=2: {t_lmax2:.1f}s for {len(r_grid)} R-points x "
          f"{len(re_grid)} R_e-points")

    # Estimate l_max=8 time: matrix size scales as (l_max+1)^4 * n_alpha^3
    # (eigenvalue decomposition of n_ch*n_alpha matrix)
    # More conservatively: time ~ n_ch^2 * n_alpha
    n_ch_2 = len(get_channels(2))
    n_ch_8 = len(get_channels(8))
    scale = (n_ch_8 / n_ch_2) ** 2
    est_lmax8 = t_lmax2 * scale
    total_est = sum(t_lmax2 * (len(get_channels(lm)) / n_ch_2) ** 2
                    for lm in L_MAX_VALUES)
    print(f"  Estimated l_max=8 time: {est_lmax8:.0f}s")
    print(f"  Estimated total time: {total_est:.0f}s")

    if est_lmax8 > FALLBACK_THRESHOLD:
        print(f"\n  ** l_max=8 estimated > {FALLBACK_THRESHOLD}s, "
              f"switching to fallback grid **")
        r_grid = FALLBACK_R_VALUES.copy()
        re_grid = FALLBACK_R_E_GRID.copy()
        n_alpha = FALLBACK_N_ALPHA
        use_fallback = True

        # Re-run timing with fallback grid
        t0 = time.time()
        R_test, E_test, _, _ = compute_pes_curve(2, r_grid, re_grid, n_alpha)
        t_lmax2 = time.time() - t0
        scale = (n_ch_8 / n_ch_2) ** 2
        est_lmax8 = t_lmax2 * scale
        total_est = sum(t_lmax2 * (len(get_channels(lm)) / n_ch_2) ** 2
                        for lm in L_MAX_VALUES)
        print(f"  Fallback l_max=2: {t_lmax2:.1f}s")
        print(f"  Fallback estimated total: {total_est:.0f}s")

    print(f"\nGrid: {len(r_grid)} R-points [{r_grid[0]:.1f}, {r_grid[-1]:.1f}], "
          f"{len(re_grid)} R_e-points [{re_grid[0]:.2f}, {re_grid[-1]:.2f}], "
          f"n_alpha={n_alpha}")
    print()

    # ===================================================================
    # Main computation: PES curves at each l_max
    # ===================================================================
    all_results = {}
    max_lmax_computed = 2  # Track highest l_max we actually compute

    for l_max in L_MAX_VALUES:
        n_ch = len(get_channels(l_max))
        mat_size = n_ch * n_alpha
        print(f"[l_max={l_max}] {n_ch} channels, matrix {mat_size}x{mat_size} ...",
              end=" ", flush=True)

        t0 = time.time()

        # Skip l_max=2 if already computed in timing test (and grid matches)
        if l_max == 2 and not use_fallback:
            R_vals, E_total = R_test, E_test
            _, _, E_elec_2, V_nn_2 = compute_pes_curve(2, r_grid, re_grid, n_alpha)
            E_elec = E_elec_2
            V_nn = V_nn_2
            # Recompute properly
            R_vals, E_total, E_elec, V_nn = compute_pes_curve(
                l_max, r_grid, re_grid, n_alpha)
        else:
            R_vals, E_total, E_elec, V_nn = compute_pes_curve(
                l_max, r_grid, re_grid, n_alpha)

        elapsed = time.time() - t0
        R_eq, E_min, success = find_r_eq(R_vals, E_total)

        all_results[l_max] = {
            'n_ch': n_ch,
            'R': R_vals.tolist(),
            'E_total': E_total.tolist(),
            'E_elec': E_elec.tolist(),
            'V_nn': V_nn.tolist(),
            'R_eq': R_eq,
            'E_min': E_min,
            'R_eq_found': success,
            'time_s': elapsed,
        }

        max_lmax_computed = l_max

        status = f"R_eq = {R_eq:.3f} bohr" if success else "NO MINIMUM"
        print(f"{status}, E_min = {E_min:.6f} Ha, {elapsed:.1f}s")

        # Safety: abort if single curve takes > 10 minutes
        if elapsed > 600:
            print(f"\n** ABORT: l_max={l_max} took {elapsed:.0f}s (>600s limit). "
                  f"Stopping here. **")
            L_MAX_VALUES_ACTUAL = [lm for lm in L_MAX_VALUES if lm <= l_max]
            break
    else:
        L_MAX_VALUES_ACTUAL = L_MAX_VALUES

    # ===================================================================
    # Decomposition: E_elec at R_small and R_large vs l_max
    # ===================================================================
    print(f"\nComputing E_elec decomposition at R={R_SMALL} and R={R_LARGE} ...",
          flush=True)

    e_elec_small = {}
    e_elec_large = {}
    for l_max in L_MAX_VALUES_ACTUAL:
        e_elec_small[l_max] = compute_e_elec_at_R(
            R_SMALL, l_max, re_grid, n_alpha)
        e_elec_large[l_max] = compute_e_elec_at_R(
            R_LARGE, l_max, re_grid, n_alpha)
        print(f"  l_max={l_max}: E_elec(R={R_SMALL}) = {e_elec_small[l_max]:.6f}, "
              f"E_elec(R={R_LARGE}) = {e_elec_large[l_max]:.6f}")

    # ===================================================================
    # Analysis
    # ===================================================================
    print("\n" + "=" * 76)
    print("RESULTS TABLE")
    print("=" * 76)
    print(f"{'l_max':>5}  {'N_ch':>5}  {'R_eq (bohr)':>12}  {'DeltaR_eq':>10}  "
          f"{'E_min (Ha)':>12}  {'DeltaE_min':>10}  {'Time(s)':>8}")
    print("-" * 76)

    r_eqs = []
    e_mins = []
    for l_max in L_MAX_VALUES_ACTUAL:
        res = all_results[l_max]
        r_eqs.append(res['R_eq'])
        e_mins.append(res['E_min'])

    for i, l_max in enumerate(L_MAX_VALUES_ACTUAL):
        res = all_results[l_max]
        dr = r_eqs[i] - r_eqs[i-1] if i > 0 else 0.0
        de = e_mins[i] - e_mins[i-1] if i > 0 else 0.0
        marker = "" if res['R_eq_found'] else " *"
        print(f"{l_max:>5}  {res['n_ch']:>5}  {res['R_eq']:>12.4f}  "
              f"{dr:>+10.4f}  {res['E_min']:>12.6f}  {de:>+10.6f}  "
              f"{res['time_s']:>8.1f}{marker}")

    # Drift increments
    print("\n" + "=" * 76)
    print("DRIFT INCREMENTS")
    print("=" * 76)

    drift_lmax = []
    drift_vals = []
    for i in range(1, len(L_MAX_VALUES_ACTUAL)):
        dl = r_eqs[i] - r_eqs[i-1]
        drift_lmax.append(L_MAX_VALUES_ACTUAL[i])
        drift_vals.append(dl)
        print(f"  l_max {L_MAX_VALUES_ACTUAL[i-1]}->{L_MAX_VALUES_ACTUAL[i]}: "
              f"DeltaR_eq = {dl:+.4f} bohr")

    # ===================================================================
    # Fit analysis
    # ===================================================================
    print("\n" + "=" * 76)
    print("FIT ANALYSIS")
    print("=" * 76)

    x_fit = np.array(L_MAX_VALUES_ACTUAL, dtype=float)
    y_fit = np.array(r_eqs)

    # Linear fit
    try:
        popt_lin, _ = curve_fit(linear_model, x_fit, y_fit)
        y_pred_lin = linear_model(x_fit, *popt_lin)
        ss_res_lin = np.sum((y_fit - y_pred_lin) ** 2)
        ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
        r2_lin = 1.0 - ss_res_lin / ss_tot if ss_tot > 0 else 0.0
        print(f"Linear: R_eq = {popt_lin[0]:.4f} + {popt_lin[1]:.4f} * l_max")
        print(f"  R^2 = {r2_lin:.6f}")
    except Exception as e:
        print(f"Linear fit failed: {e}")
        r2_lin = -1
        popt_lin = None

    # Saturating exponential fit
    try:
        # Initial guess: R_inf ~ last R_eq + 0.5, c ~ spread, beta ~ 0.3
        p0 = [r_eqs[-1] + 0.5, r_eqs[-1] - r_eqs[0] + 0.5, 0.3]
        popt_exp, pcov_exp = curve_fit(exp_saturation, x_fit, y_fit, p0=p0,
                                        maxfev=10000)
        y_pred_exp = exp_saturation(x_fit, *popt_exp)
        ss_res_exp = np.sum((y_fit - y_pred_exp) ** 2)
        r2_exp = 1.0 - ss_res_exp / ss_tot if ss_tot > 0 else 0.0
        R_inf = popt_exp[0]
        print(f"\nSaturating exponential: R_eq = {R_inf:.4f} - "
              f"{popt_exp[1]:.4f} * exp(-{popt_exp[2]:.4f} * l_max)")
        print(f"  R_inf = {R_inf:.4f} bohr")
        print(f"  R^2 = {r2_exp:.6f}")
        print(f"  R_inf - R_exp = {R_inf - R_EXP:+.4f} bohr "
              f"({100*abs(R_inf - R_EXP)/R_EXP:.1f}% error)")
    except Exception as e:
        print(f"Exponential fit failed: {e}")
        r2_exp = -1
        popt_exp = None
        R_inf = None

    # Drift increment fit
    if len(drift_vals) >= 3:
        abs_drifts = np.array([abs(d) for d in drift_vals])
        x_drift = np.array(drift_lmax, dtype=float)
        try:
            p0_d = [abs_drifts[0] * 3, 0.3]
            popt_drift, _ = curve_fit(exp_decay, x_drift, abs_drifts,
                                       p0=p0_d, maxfev=10000)
            print(f"\nDrift increments: |DeltaR_eq| = {popt_drift[0]:.4f} * "
                  f"exp(-{popt_drift[1]:.4f} * l_max)")
            print(f"  Decay rate gamma = {popt_drift[1]:.4f}")
            if popt_drift[1] > 0.1:
                print(f"  -> Geometric decrease confirmed (gamma > 0.1)")
            else:
                print(f"  -> Weak or no decrease (gamma < 0.1)")
        except Exception as e:
            print(f"Drift increment fit failed: {e}")
            popt_drift = None
    else:
        print("\nInsufficient data points for drift increment fit")
        popt_drift = None

    # ===================================================================
    # Bottom-line verdict
    # ===================================================================
    print("\n" + "=" * 76)
    print("BOTTOM LINE")
    print("=" * 76)

    if popt_exp is not None and r2_exp > r2_lin and r2_exp > 0.95:
        print(f"R_eq SATURATES to R_inf = {R_inf:.3f} bohr "
              f"({100*abs(R_inf - R_EXP)/R_EXP:.1f}% from experiment).")
        print(f"Exponential fit (R^2={r2_exp:.4f}) is better than linear "
              f"(R^2={r2_lin:.4f}).")
        print(f"The l-dependent PK model converges; remaining "
              f"{100*abs(R_inf - R_EXP)/R_EXP:.1f}% error is from the "
              f"composed geometry approximation.")
    elif popt_lin is not None and r2_lin > 0.95:
        print(f"R_eq DIVERGES linearly at rate {popt_lin[1]:.3f} bohr/l_max.")
        print(f"Linear fit (R^2={r2_lin:.4f}) dominates.")
        print(f"The Z_eff/PK model is fundamentally flawed for l_max > 2.")
    else:
        print("Inconclusive: neither linear nor exponential fit is dominant.")
        if len(drift_vals) >= 2:
            if abs(drift_vals[-1]) < abs(drift_vals[0]):
                print("Drift increments are decreasing, suggesting eventual "
                      "saturation.")
            else:
                print("Drift increments are NOT decreasing.")

    # E_min trend
    if len(e_mins) >= 3:
        de_last = e_mins[-1] - e_mins[-2]
        de_first = e_mins[1] - e_mins[0]
        if abs(de_last) < abs(de_first) * 0.5:
            print(f"E_min stabilizing: first increment {de_first:+.6f}, "
                  f"last increment {de_last:+.6f} Ha.")
        else:
            print(f"E_min NOT stabilizing: first increment {de_first:+.6f}, "
                  f"last increment {de_last:+.6f} Ha.")

    # Decomposition analysis
    print(f"\nDifferential correlation analysis:")
    de_small_total = e_elec_small[L_MAX_VALUES_ACTUAL[-1]] - e_elec_small[L_MAX_VALUES_ACTUAL[0]]
    de_large_total = e_elec_large[L_MAX_VALUES_ACTUAL[-1]] - e_elec_large[L_MAX_VALUES_ACTUAL[0]]
    print(f"  E_elec drop at R={R_SMALL}: {de_small_total:+.6f} Ha "
          f"(l_max={L_MAX_VALUES_ACTUAL[0]}->{L_MAX_VALUES_ACTUAL[-1]})")
    print(f"  E_elec drop at R={R_LARGE}: {de_large_total:+.6f} Ha "
          f"(l_max={L_MAX_VALUES_ACTUAL[0]}->{L_MAX_VALUES_ACTUAL[-1]})")
    if abs(de_large_total) > abs(de_small_total):
        ratio = abs(de_large_total) / abs(de_small_total) if abs(de_small_total) > 1e-10 else float('inf')
        print(f"  Large-R stabilizes {ratio:.2f}x faster -> confirms "
              f"differential correlation mechanism.")
    else:
        print(f"  Small-R stabilizes faster -> differential correlation "
              f"NOT confirmed.")

    # ===================================================================
    # Save data
    # ===================================================================
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)

    output = {
        'parameters': {
            'Z_A_eff': Z_A, 'Z_B': Z_B,
            'PK_A': PK_A, 'PK_B': PK_B_param,
            'PK_mode': 'l_dependent',
            'R_exp': R_EXP,
            'n_alpha': n_alpha, 'n_re': len(re_grid),
            'n_R': len(r_grid),
            'fallback_used': use_fallback,
        },
        'results': {},
        'decomposition': {
            'R_small': R_SMALL,
            'R_large': R_LARGE,
            'E_elec_small': {str(k): v for k, v in e_elec_small.items()},
            'E_elec_large': {str(k): v for k, v in e_elec_large.items()},
        },
        'fits': {
            'linear': {
                'params': popt_lin.tolist() if popt_lin is not None else None,
                'R2': r2_lin,
            },
            'exponential': {
                'params': popt_exp.tolist() if popt_exp is not None else None,
                'R2': r2_exp,
                'R_inf': R_inf,
            },
            'drift_decay': {
                'params': popt_drift.tolist() if popt_drift is not None else None,
            },
        },
        'summary': {
            'l_max_values': L_MAX_VALUES_ACTUAL,
            'R_eq_values': r_eqs,
            'E_min_values': e_mins,
            'drift_increments': drift_vals,
        },
    }
    for l_max in L_MAX_VALUES_ACTUAL:
        output['results'][str(l_max)] = all_results[l_max]

    json_path = data_dir / 'pk_saturation_test.json'
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

        # --- Plot A (money plot): R_eq vs l_max ---
        fig, ax = plt.subplots(1, 1, figsize=(9, 6))

        ax.plot(L_MAX_VALUES_ACTUAL, r_eqs, 'o-', color='#1f77b4',
                linewidth=2.5, markersize=10, label='l-dependent PK', zorder=5)

        ax.axhline(R_EXP, color='gray', linestyle='--', linewidth=1.5,
                    alpha=0.7, label=f'Experiment ({R_EXP} bohr)')

        # Overlay fits
        x_fine = np.linspace(min(L_MAX_VALUES_ACTUAL) - 0.5,
                             max(L_MAX_VALUES_ACTUAL) + 1.5, 100)
        if popt_lin is not None:
            ax.plot(x_fine, linear_model(x_fine, *popt_lin), '--',
                    color='#d62728', alpha=0.6,
                    label=f'Linear (R²={r2_lin:.4f})')
        if popt_exp is not None:
            ax.plot(x_fine, exp_saturation(x_fine, *popt_exp), '-.',
                    color='#2ca02c', alpha=0.8, linewidth=2,
                    label=f'Exp. saturation (R²={r2_exp:.4f})\n'
                          f'R$_\\infty$ = {R_inf:.3f} bohr')
            ax.axhline(R_inf, color='#2ca02c', linestyle=':',
                        alpha=0.4)

        ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
        ax.set_ylabel('$R_{\\rm eq}$ (bohr)', fontsize=14)
        ax.set_title('LiH $R_{\\rm eq}$ vs $l_{\\rm max}$ — '
                     'l-dependent PK Saturation Test', fontsize=13)
        ax.legend(fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(L_MAX_VALUES_ACTUAL)
        fig.tight_layout()
        fig.savefig(plot_dir / 'pk_saturation_req_vs_lmax.png', dpi=150)
        plt.close(fig)
        print(f"Plot A saved to debug/plots/pk_saturation_req_vs_lmax.png")

        # --- Plot B: Drift increments on log scale ---
        if len(drift_vals) >= 2:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))

            abs_drifts = [abs(d) for d in drift_vals]
            ax.semilogy(drift_lmax, abs_drifts, 'o-', color='#1f77b4',
                        linewidth=2, markersize=10)

            if popt_drift is not None:
                x_d_fine = np.linspace(min(drift_lmax), max(drift_lmax) + 1, 50)
                ax.semilogy(x_d_fine, exp_decay(x_d_fine, *popt_drift), '--',
                            color='#d62728', alpha=0.7,
                            label=f'Fit: A·exp(-{popt_drift[1]:.3f}·l_max)')
                ax.legend(fontsize=11)

            ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
            ax.set_ylabel('$|\\Delta R_{\\rm eq}|$ (bohr)', fontsize=14)
            ax.set_title('Drift Increments vs $l_{\\rm max}$ (log scale)',
                         fontsize=13)
            ax.grid(True, alpha=0.3, which='both')
            ax.set_xticks(drift_lmax)
            fig.tight_layout()
            fig.savefig(plot_dir / 'pk_saturation_drift_increments.png', dpi=150)
            plt.close(fig)
            print(f"Plot B saved to debug/plots/pk_saturation_drift_increments.png")

        # --- Plot C: PES curves for l_max = 2, 4, 6, 8 (or available) ---
        plot_lmax = [lm for lm in [2, 4, 6, 8] if lm in all_results]
        if len(plot_lmax) >= 2:
            fig, ax = plt.subplots(1, 1, figsize=(9, 6))
            colors_pes = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                          '#9467bd', '#8c564b']

            for i, l_max in enumerate(plot_lmax):
                res = all_results[l_max]
                ax.plot(res['R'], res['E_total'], '-o',
                        color=colors_pes[i % len(colors_pes)],
                        markersize=4, linewidth=1.5,
                        label=f'$l_{{\\rm max}}={l_max}$')
                if res['R_eq_found']:
                    ax.axvline(res['R_eq'], color=colors_pes[i % len(colors_pes)],
                               alpha=0.25, linestyle='--')

            ax.axvline(R_EXP, color='gray', linestyle='--', alpha=0.5,
                       label=f'Expt ({R_EXP})')
            ax.set_xlabel('$R$ (bohr)', fontsize=14)
            ax.set_ylabel('$E_{\\rm total}$ (Ha)', fontsize=14)
            ax.set_title('LiH PES Evolution with $l_{\\rm max}$ '
                         '(l-dependent PK)', fontsize=13)
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(plot_dir / 'pk_saturation_pes_curves.png', dpi=150)
            plt.close(fig)
            print(f"Plot C saved to debug/plots/pk_saturation_pes_curves.png")

        # --- Plot D: E_elec at R_small and R_large vs l_max ---
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        lm_list = L_MAX_VALUES_ACTUAL
        e_s = [e_elec_small[lm] for lm in lm_list]
        e_l = [e_elec_large[lm] for lm in lm_list]

        # Normalize to show relative change from l_max=2
        e_s_ref = e_s[0]
        e_l_ref = e_l[0]
        de_s = [e - e_s_ref for e in e_s]
        de_l = [e - e_l_ref for e in e_l]

        ax.plot(lm_list, de_s, 'o-', color='#1f77b4', linewidth=2,
                markersize=8, label=f'$R = {R_SMALL}$ bohr (near eq.)')
        ax.plot(lm_list, de_l, 's-', color='#d62728', linewidth=2,
                markersize=8, label=f'$R = {R_LARGE}$ bohr (extended)')

        ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
        ax.set_ylabel('$\\Delta E_{\\rm elec}$ from $l_{\\rm max}=2$ (Ha)',
                      fontsize=14)
        ax.set_title('Differential Correlation: E_elec Change vs $l_{\\rm max}$',
                     fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(lm_list)
        fig.tight_layout()
        fig.savefig(plot_dir / 'pk_saturation_differential_corr.png', dpi=150)
        plt.close(fig)
        print(f"Plot D saved to debug/plots/pk_saturation_differential_corr.png")

        # --- Plot E (bonus): E_min vs l_max ---
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        ax.plot(L_MAX_VALUES_ACTUAL, e_mins, 'o-', color='#1f77b4',
                linewidth=2, markersize=10)
        ax.set_xlabel('$l_{\\rm max}$', fontsize=14)
        ax.set_ylabel('$E_{\\rm min}$ (Ha)', fontsize=14)
        ax.set_title('Minimum Energy vs $l_{\\rm max}$ (l-dependent PK)',
                     fontsize=13)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(L_MAX_VALUES_ACTUAL)
        fig.tight_layout()
        fig.savefig(plot_dir / 'pk_saturation_emin_vs_lmax.png', dpi=150)
        plt.close(fig)
        print(f"Plot E saved to debug/plots/pk_saturation_emin_vs_lmax.png")

    except ImportError:
        print("matplotlib not available -- skipping plots.")

    print("\nDone.")


if __name__ == '__main__':
    main()
