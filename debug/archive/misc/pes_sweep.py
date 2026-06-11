"""
PES Sweep and R_eq Determination for Three PK Modes.

Task 7: Compute Born-Oppenheimer potential energy curves E(R) for LiH
using three PK pseudopotential modes (channel-blind, l-dependent, projected)
at l_max = 2, 3, 4.  Extract R_eq for each and compare against the
Paper 17 l_max divergence.

Method: Adiabatic approximation (Option B).
  At each R, sweep R_e, compute U(R_e) = [mu_0(R_e;R) + 15/8] / R_e^2,
  then E_elec(R) = min_{R_e} U(R_e).
  E_total(R) = E_elec(R) + Z_A_eff * Z_B / R.

References:
  Paper 15 Sec IV (Eq. 28, hyperradial equation)
  Paper 17 Sec III.B (PK parameters), Sec V.A (l_max divergence)
  Task 6 projected_pk.py (angular solver with three PK modes)
"""

import sys
import json
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.level4_multichannel import (
    _channel_list,
    build_angular_hamiltonian,
    compute_pk_pseudopotential,
)

# === Physical parameters ===
Z_A = 2.69            # Clementi-Raimondi Z_eff for Li
Z_B = 1.0             # H
PK_A = 6.93           # Ha*bohr^2 (PK amplitude)
PK_B_param = 7.00     # bohr^-2 (PK Gaussian exponent)
R_EXP = 3.015         # Experimental R_eq (bohr)

# PES scan grid (extended to 8.0 to capture channel-blind minima at high l_max)
R_VALUES = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0])

# Hyperradial grid: 30 points, logarithmic spacing from 0.3 to 6.0 bohr
N_RE = 30
R_E_GRID = np.logspace(np.log10(0.3), np.log10(6.0), N_RE)

# Angular grid (reduced for speed)
N_ALPHA = 50

# l_max values
L_MAX_VALUES = [2, 3, 4]

# PK modes
PK_MODES = ['channel_blind', 'l_dependent', 'projected']


def is_homonuclear() -> bool:
    return Z_A == Z_B


def get_channels(l_max: int) -> list:
    return _channel_list(l_max, homonuclear=is_homonuclear())


def solve_angular_mu0(R: float, R_e: float, l_max: int,
                      pk_mode: str) -> float:
    """
    Solve the angular eigenvalue problem at (R, R_e) for a given PK mode.
    Returns mu_0 (lowest angular eigenvalue).
    """
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))

    h = (np.pi / 2) / (N_ALPHA + 1)
    alpha = (np.arange(N_ALPHA) + 1) * h

    if pk_mode == 'projected':
        return _solve_projected(R, R_e, l_max, rho, z0, alpha)

    # Channel-blind or l-dependent: direct solve
    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_param,
        'atom': 'A',
        'channel_mode': pk_mode,
    }]

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=pk_pots,
    )
    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def _solve_projected(R: float, R_e: float, l_max: int,
                     rho: float, z0: float, alpha: np.ndarray) -> float:
    """Projected PK: solve without PK, get weights, re-solve with weighted PK."""
    n_alpha = len(alpha)
    channels = get_channels(l_max)
    n_ch = len(channels)

    # Step 1: Solve without PK
    H_base = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )
    evals_nopk, evecs_nopk = eigh(H_base)
    ground_vec = evecs_nopk[:, 0]

    # Step 2: Per-channel weight from ground eigenstate
    per_ch_weight = np.zeros(n_ch)
    for ic in range(n_ch):
        chunk = ground_vec[ic * n_alpha: (ic + 1) * n_alpha]
        per_ch_weight[ic] = np.sum(chunk ** 2)
    total = np.sum(per_ch_weight)
    if total > 1e-30:
        per_ch_weight /= total

    # Fractional l=0 content for each electron
    w1_total = sum(per_ch_weight[ic] for ic, ch in enumerate(channels)
                   if ch[0] == 0)
    w2_total = sum(per_ch_weight[ic] for ic, ch in enumerate(channels)
                   if ch[1] == 0)

    # Step 3: Compute PK potential profiles
    rho_A_re = (R / 2.0 + z0) / R_e
    rho_B_re = (R / 2.0 - z0) / R_e

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_param,
        'atom': 'A',
    }]
    V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
        alpha, rho_A_re, R_e, pk_pots, rho_B=rho_B_re,
    )

    # Step 4: Add projected PK
    H = H_base.copy()
    for ic, ch in enumerate(channels):
        l1, l2 = ch[0], ch[1]
        w1 = w1_total if l1 == 0 else 0.0
        w2 = w2_total if l2 == 0 else 0.0
        if w1 == 0.0 and w2 == 0.0:
            continue
        for i in range(n_alpha):
            ii = ic * n_alpha + i
            H[ii, ii] += R_e * w1 * V_pk_e1[i]
            H[ii, ii] += R_e * w2 * V_pk_e2[i]

    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def compute_pes_curve(l_max: int, pk_mode: str) -> tuple:
    """
    Compute E_total(R) for a given l_max and PK mode.

    Returns (R_values, E_total_values, E_elec_values, V_nn_values).
    """
    E_elec = np.zeros(len(R_VALUES))
    V_nn = np.zeros(len(R_VALUES))

    for ir, R in enumerate(R_VALUES):
        # Sweep R_e, compute adiabatic potential U(R_e)
        U_vals = np.zeros(N_RE)
        for ire, R_e in enumerate(R_E_GRID):
            mu0 = solve_angular_mu0(R, R_e, l_max, pk_mode)
            U_vals[ire] = (mu0 + 15.0 / 8.0) / (R_e ** 2)

        # Adiabatic approximation: E_elec = min U(R_e)
        E_elec[ir] = np.min(U_vals)

        # Nuclear repulsion (using Z_A_eff * Z_B for consistency with
        # electronic Hamiltonian)
        V_nn[ir] = Z_A * Z_B / R

    E_total = E_elec + V_nn
    return R_VALUES.copy(), E_total, E_elec, V_nn


def find_r_eq(R_vals: np.ndarray, E_vals: np.ndarray) -> tuple:
    """
    Find R_eq by cubic interpolation near the minimum.

    Returns (R_eq, E_min, success).
    """
    i_min = np.argmin(E_vals)

    # Check if minimum is at an edge (no clear well)
    if i_min == 0:
        return R_vals[0], E_vals[0], False
    if i_min == len(R_vals) - 1:
        return R_vals[-1], E_vals[-1], False

    # Use 5 points around minimum for cubic spline (or all if fewer)
    i_lo = max(0, i_min - 2)
    i_hi = min(len(R_vals), i_min + 3)
    R_sub = R_vals[i_lo:i_hi]
    E_sub = E_vals[i_lo:i_hi]

    if len(R_sub) < 3:
        return R_vals[i_min], E_vals[i_min], True

    cs = CubicSpline(R_sub, E_sub)
    R_fine = np.linspace(R_sub[0], R_sub[-1], 500)
    E_fine = cs(R_fine)
    i_fine_min = np.argmin(E_fine)

    return float(R_fine[i_fine_min]), float(E_fine[i_fine_min]), True


def main():
    print("=" * 76)
    print("PES Sweep: R_eq vs l_max for Three PK Modes (LiH)")
    print("=" * 76)
    print(f"Z_A_eff = {Z_A}, Z_B = {Z_B}")
    print(f"PK: A = {PK_A} Ha*bohr^2, B = {PK_B_param} bohr^-2")
    print(f"R grid: {R_VALUES} bohr")
    print(f"R_e grid: {N_RE} points, {R_E_GRID[0]:.2f} to {R_E_GRID[-1]:.2f} bohr (log-spaced)")
    print(f"n_alpha = {N_ALPHA}")
    print(f"V_NN = Z_A_eff * Z_B / R = {Z_A} * {Z_B} / R")
    print(f"Experimental R_eq = {R_EXP} bohr")
    print()

    # Diagnostic: print projected PK weights at each l_max
    print("--- Projected PK weight diagnostic (at R=3.015, R_e=1.5) ---")
    for l_max in L_MAX_VALUES:
        channels = get_channels(l_max)
        n_ch = len(channels)
        rho_diag = R_EXP / (2.0 * 1.5)
        z0_diag = R_EXP * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
        h = (np.pi / 2) / (N_ALPHA + 1)
        alpha_diag = (np.arange(N_ALPHA) + 1) * h
        H_diag = build_angular_hamiltonian(
            alpha_diag, rho_diag, 1.5, l_max=l_max, Z=1.0,
            m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0_diag,
            pk_potentials=None,
        )
        _, evecs_diag = eigh(H_diag)
        gv = evecs_diag[:, 0]
        per_ch = np.array([np.sum(gv[ic*N_ALPHA:(ic+1)*N_ALPHA]**2)
                           for ic in range(n_ch)])
        per_ch /= per_ch.sum()
        w1 = sum(per_ch[ic] for ic, ch in enumerate(channels) if ch[0] == 0)
        w2 = sum(per_ch[ic] for ic, ch in enumerate(channels) if ch[1] == 0)
        print(f"  l_max={l_max}: w1_l0={w1:.4f}, w2_l0={w2:.4f}, "
              f"channels={channels}, weights={[f'{w:.3f}' for w in per_ch]}")
    print()

    # Storage for all results
    all_results = {}

    total_curves = len(PK_MODES) * len(L_MAX_VALUES)
    curve_num = 0

    for pk_mode in PK_MODES:
        all_results[pk_mode] = {}
        for l_max in L_MAX_VALUES:
            curve_num += 1
            print(f"[{curve_num}/{total_curves}] Computing PES: "
                  f"{pk_mode}, l_max={l_max} ...", flush=True)

            R_vals, E_total, E_elec, V_nn = compute_pes_curve(l_max, pk_mode)
            R_eq, E_min, success = find_r_eq(R_vals, E_total)

            all_results[pk_mode][l_max] = {
                'R': R_vals.tolist(),
                'E_total': E_total.tolist(),
                'E_elec': E_elec.tolist(),
                'V_nn': V_nn.tolist(),
                'R_eq': R_eq,
                'E_min': E_min,
                'R_eq_found': success,
            }

            status = f"R_eq = {R_eq:.3f}" if success else "NO MINIMUM FOUND"
            print(f"       {status}  (E_min = {E_min:.6f} Ha)")

    # ===================================================================
    # Summary table
    # ===================================================================
    print("\n" + "=" * 76)
    print("SUMMARY: R_eq vs l_max for Each PK Mode")
    print("=" * 76)
    print(f"{'PK Mode':>15}  {'l_max':>5}  {'R_eq (bohr)':>12}  "
          f"{'Delta from expt':>16}  {'% error':>8}  {'E_min (Ha)':>12}")
    print("-" * 76)

    for pk_mode in PK_MODES:
        for l_max in L_MAX_VALUES:
            res = all_results[pk_mode][l_max]
            R_eq = res['R_eq']
            delta = R_eq - R_EXP
            pct = 100.0 * abs(delta) / R_EXP
            found = res['R_eq_found']
            marker = "" if found else " *"
            print(f"{pk_mode:>15}  {l_max:>5}  {R_eq:>12.3f}  "
                  f"{delta:>+16.3f}  {pct:>7.1f}%  {res['E_min']:>12.6f}{marker}")
        print()

    print("* = minimum at edge of R grid (not a true equilibrium)")
    print(f"\nExperimental R_eq = {R_EXP} bohr")

    # ===================================================================
    # Convergence assessment
    # ===================================================================
    print("\n" + "=" * 76)
    print("CONVERGENCE ASSESSMENT: R_eq trend with l_max")
    print("=" * 76)

    for pk_mode in PK_MODES:
        r_eqs = [all_results[pk_mode][lm]['R_eq'] for lm in L_MAX_VALUES]
        spread = max(r_eqs) - min(r_eqs)
        drift = r_eqs[-1] - r_eqs[0]

        if drift > 0.3:
            verdict = f"DIVERGENT (R_eq drifts outward by {drift:.2f} bohr)"
        elif drift < -0.3:
            verdict = f"CONVERGENT (R_eq drifts inward by {abs(drift):.2f} bohr)"
        elif spread < 0.15:
            verdict = f"STABLE (spread = {spread:.3f} bohr)"
        else:
            verdict = f"MODERATE DRIFT (spread = {spread:.3f} bohr)"

        print(f"  {pk_mode:>15}: {verdict}")
        print(f"    R_eq values: {', '.join(f'{r:.3f}' for r in r_eqs)}")

    # ===================================================================
    # Save data
    # ===================================================================
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)

    output = {
        'parameters': {
            'Z_A_eff': Z_A, 'Z_B': Z_B,
            'PK_A': PK_A, 'PK_B': PK_B_param,
            'R_exp': R_EXP,
            'n_alpha': N_ALPHA, 'n_re': N_RE,
            'V_NN': 'Z_A_eff * Z_B / R',
            'method': 'adiabatic approximation (min_Re U(Re))',
        },
        'results': {},
    }
    for pk_mode in PK_MODES:
        output['results'][pk_mode] = {}
        for l_max in L_MAX_VALUES:
            output['results'][pk_mode][str(l_max)] = all_results[pk_mode][l_max]

    with open(data_dir / 'pes_sweep_pk_modes.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nData saved to debug/data/pes_sweep_pk_modes.json")

    # ===================================================================
    # Plots
    # ===================================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plot_dir = Path(__file__).parent / 'plots'
        plot_dir.mkdir(exist_ok=True)

        colors = {'channel_blind': '#d62728', 'l_dependent': '#1f77b4',
                  'projected': '#2ca02c'}
        ls_map = {2: '-', 3: '--', 4: ':'}
        marker_map = {2: 'o', 3: 's', 4: '^'}

        # --- Plot A: PES curves, one panel per PK mode ---
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

        for ax, pk_mode in zip(axes, PK_MODES):
            for l_max in L_MAX_VALUES:
                res = all_results[pk_mode][l_max]
                R = res['R']
                E = res['E_total']
                ax.plot(R, E, ls_map[l_max], color=colors[pk_mode],
                        marker=marker_map[l_max], markersize=5,
                        label=f'$l_{{\\rm max}}={l_max}$', linewidth=1.5)

                # Mark R_eq
                if res['R_eq_found']:
                    ax.axvline(res['R_eq'], color=colors[pk_mode],
                               alpha=0.3, linestyle=ls_map[l_max])

            ax.axvline(R_EXP, color='gray', linestyle='--', alpha=0.5,
                       label=f'Expt ({R_EXP})')
            ax.set_xlabel('$R$ (bohr)', fontsize=12)
            ax.set_title(pk_mode.replace('_', '-'), fontsize=13)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

        axes[0].set_ylabel('$E_{\\rm total}$ (Ha)', fontsize=12)
        fig.suptitle('LiH PES: PK Mode Comparison at Different $l_{\\rm max}$',
                     fontsize=14)
        fig.tight_layout()
        fig.savefig(plot_dir / 'pes_sweep_pk_modes.png', dpi=150)
        plt.close(fig)
        print(f"Plot A saved to debug/plots/pes_sweep_pk_modes.png")

        # --- Plot B: R_eq vs l_max (the money plot) ---
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        for pk_mode in PK_MODES:
            r_eqs = [all_results[pk_mode][lm]['R_eq'] for lm in L_MAX_VALUES]
            ax.plot(L_MAX_VALUES, r_eqs, '-o', color=colors[pk_mode],
                    label=pk_mode.replace('_', '-'),
                    linewidth=2, markersize=8)

        ax.axhline(R_EXP, color='gray', linestyle='--', alpha=0.7,
                    label=f'Expt ({R_EXP} bohr)')
        ax.set_xlabel('$l_{\\rm max}$', fontsize=13)
        ax.set_ylabel('$R_{\\rm eq}$ (bohr)', fontsize=13)
        ax.set_title('LiH $R_{\\rm eq}$ vs $l_{\\rm max}$ for Three PK Modes',
                     fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(L_MAX_VALUES)
        fig.tight_layout()
        fig.savefig(plot_dir / 'pes_sweep_req_vs_lmax.png', dpi=150)
        plt.close(fig)
        print(f"Plot B saved to debug/plots/pes_sweep_req_vs_lmax.png")

    except ImportError:
        print("matplotlib not available -- skipping plots.")


if __name__ == '__main__':
    main()
