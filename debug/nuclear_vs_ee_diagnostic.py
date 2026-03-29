#!/usr/bin/env python3
"""
Task 5: Nuclear vs E-E Channel Spreading Diagnostic + Higher l_max Convergence.

Part A: Decomposes wavefunction spreading into nuclear vs electron-electron
contributions by solving the coupled-channel alpha-equation three ways:
  1. Nuclear only (W_ee = 0)
  2. E-E only (W_nuc diagonal only, no off-diagonal nuclear coupling)
  3. Full coupling (both)

Part B: Higher l_max convergence extrapolation for 6:1 system.

Reuses the coupled alpha-equation solver from angular_wavefunction_diagnostic.py.
"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import json
import sys
import time
from typing import List, Tuple, Dict, Optional

sys.path.insert(0, str(Path(__file__).parent))
from eigenchannel_diagnostic import gaunt, build_channel_list
from angular_wavefunction_diagnostic import (
    f_k_vec, ee_coupling_element, nuclear_coupling_vec,
    richardson_extrapolate,
)

PLOT_DIR = Path(__file__).parent / "plots"
DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Modified coupling builder with separate nuclear/ee control
# --------------------------------------------------------------------------- #

def build_alpha_coupling_decomposed(
    channels: List[Tuple[int, int]],
    alpha_grid: np.ndarray,
    rho: float,
    Z_A: float, Z_B: float,
    mode: str = 'full',
) -> np.ndarray:
    """
    Build W_ij(alpha) coupling array with decomposition control.

    mode:
      'full'         - full nuclear + e-e coupling (standard)
      'nuclear_only' - nuclear coupling only, W_ee = 0
      'ee_only'      - e-e coupling + diagonal-only nuclear (no off-diagonal nuclear)
    """
    Z_tot = Z_A + Z_B
    rho_A = rho * 2.0 * Z_B / Z_tot
    rho_B = rho * 2.0 * Z_A / Z_tot

    n_ch = len(channels)
    n_alpha = len(alpha_grid)
    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    W = np.zeros((n_ch, n_ch, n_alpha))

    for i, (l1p, l2p) in enumerate(channels):
        for j, (l1, l2) in enumerate(channels):
            # Nuclear coupling
            if mode in ('full', 'nuclear_only'):
                # Full nuclear: all matrix elements
                if l2p == l2:
                    W[i, j, :] += nuclear_coupling_vec(
                        l1p, l1, cos_a, Z_A, rho_A, Z_B, rho_B)
                if l1p == l1:
                    W[i, j, :] += nuclear_coupling_vec(
                        l2p, l2, sin_a, Z_A, rho_A, Z_B, rho_B)
            elif mode == 'ee_only':
                # Diagonal nuclear only (no off-diagonal nuclear coupling)
                if i == j:
                    if l2p == l2:
                        W[i, j, :] += nuclear_coupling_vec(
                            l1p, l1, cos_a, Z_A, rho_A, Z_B, rho_B)
                    if l1p == l1:
                        W[i, j, :] += nuclear_coupling_vec(
                            l2p, l2, sin_a, Z_A, rho_A, Z_B, rho_B)

            # E-e coupling
            if mode in ('full', 'ee_only'):
                W[i, j, :] += ee_coupling_element(l1p, l2p, l1, l2, cos_a, sin_a)

            # Diagonal: centrifugal + constant
            if i == j:
                centrifugal = np.zeros(n_alpha)
                if l1 > 0:
                    centrifugal += l1 * (l1 + 1) / (2.0 * cos_a**2)
                if l2 > 0:
                    centrifugal += l2 * (l2 + 1) / (2.0 * sin_a**2)
                W[i, i, :] += centrifugal - 2.0

    return W


def solve_alpha_decomposed(
    channels: List[Tuple[int, int]],
    rho: float,
    Z_A: float, Z_B: float,
    N_alpha: int = 200,
    mode: str = 'full',
    n_eig: int = 5,
) -> Dict:
    """
    Solve the coupled-channel alpha eigenvalue problem with decomposition control.
    """
    n_ch = len(channels)
    alpha_grid = np.linspace(0, np.pi / 2, N_alpha + 2)[1:-1]
    da = alpha_grid[1] - alpha_grid[0]

    W = build_alpha_coupling_decomposed(
        channels, alpha_grid, rho, Z_A, Z_B, mode)

    N_total = n_ch * N_alpha
    kinetic_diag = 1.0 / da**2
    kinetic_off = -0.5 / da**2

    use_sparse = N_total > 5000

    if use_sparse:
        H = lil_matrix((N_total, N_total))
    else:
        H = np.zeros((N_total, N_total))

    # Kinetic energy (tridiagonal per channel)
    for i_ch in range(n_ch):
        off = i_ch * N_alpha
        for k in range(N_alpha):
            H[off + k, off + k] += kinetic_diag
            if k > 0:
                H[off + k, off + k - 1] += kinetic_off
            if k < N_alpha - 1:
                H[off + k, off + k + 1] += kinetic_off

    # Potential coupling
    for i_ch in range(n_ch):
        for j_ch in range(n_ch):
            w_ij = W[i_ch, j_ch, :]
            if np.max(np.abs(w_ij)) < 1e-15:
                continue
            i_off = i_ch * N_alpha
            j_off = j_ch * N_alpha
            for k in range(N_alpha):
                H[i_off + k, j_off + k] += w_ij[k]

    # Solve
    actual_neig = min(n_eig, N_total - 2)
    if use_sparse:
        H_csr = csr_matrix(H)
        eigenvalues, eigenvectors = eigsh(H_csr, k=actual_neig, which='SA')
        idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
    else:
        eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, actual_neig - 1])

    mu0 = eigenvalues[0]
    psi0 = eigenvectors[:, 0]

    # Channel weights
    weights = np.zeros(n_ch)
    for i_ch in range(n_ch):
        offset = i_ch * N_alpha
        u_i = psi0[offset:offset + N_alpha]
        weights[i_ch] = np.sum(u_i**2) * da

    total_weight = np.sum(weights)
    if total_weight > 1e-15:
        weights /= total_weight

    return {
        'mu0': float(mu0),
        'eigenvalues': [float(e) for e in eigenvalues],
        'channels': channels,
        'weights': weights,
        'n_ch': n_ch,
    }


def extract_weight_stats(channels, weights):
    """Extract summary statistics from channel weights."""
    w00 = 0.0
    w_dipolar = 0.0
    for i, (l1, l2) in enumerate(channels):
        if l1 == 0 and l2 == 0:
            w00 = weights[i]
        if (l1 == 0 and l2 == 1) or (l1 == 1 and l2 == 0):
            w_dipolar += weights[i]

    sorted_w = np.sort(weights)[::-1]
    cum_w = np.cumsum(sorted_w)
    n90 = int(np.searchsorted(cum_w, 0.90) + 1)
    n95 = int(np.searchsorted(cum_w, 0.95) + 1)
    n99 = int(np.searchsorted(cum_w, 0.99) + 1)

    return {
        'w00': float(w00),
        'w_dipolar': float(w_dipolar),
        'n90': n90,
        'n95': n95,
        'n99': n99,
    }


# --------------------------------------------------------------------------- #
# Part A: Nuclear vs E-E decomposition
# --------------------------------------------------------------------------- #

def run_part_a():
    """Nuclear vs E-E channel spreading decomposition for 6:1 system."""
    print("=" * 80)
    print("PART A: NUCLEAR vs E-E CHANNEL SPREADING DECOMPOSITION")
    print("System: 6:1 (Z_A=6, Z_B=1), charge-center origin, rho=0.7, l_max=4")
    print("=" * 80)

    Z_A, Z_B = 6.0, 1.0
    rho = 0.7
    l_max = 4
    N_alpha = 300

    channels = build_channel_list(l_max, homonuclear=False)
    n_ch = len(channels)
    print(f"Channels: {n_ch} (l_max={l_max})")

    modes = [
        ('nuclear_only', 'Nuclear only (W_ee=0)'),
        ('ee_only', 'E-E only (diagonal nuclear)'),
        ('full', 'Full coupling'),
    ]

    results = {}
    for mode, label in modes:
        t0 = time.time()
        print(f"\n  Solving: {label} ... ", end='', flush=True)
        r = solve_alpha_decomposed(
            channels, rho, Z_A, Z_B, N_alpha=N_alpha, mode=mode)
        dt = time.time() - t0
        stats = extract_weight_stats(channels, r['weights'])
        results[mode] = {**r, **stats, 'label': label, 'time': dt}
        print(f"mu0={r['mu0']:.6f}  [{dt:.1f}s]")

    # Summary table
    print("\n" + "=" * 90)
    print("PART A SUMMARY TABLE")
    print("=" * 90)
    header = (f"{'Configuration':<30} {'mu0':<12} {'w(0,0)':<10} "
              f"{'w(dipolar)':<12} {'N(90%)':<8} {'N(95%)':<8}")
    print(header)
    print("-" * len(header))

    for mode, label in modes:
        r = results[mode]
        print(f"{label:<30} {r['mu0']:<12.6f} {r['w00']:<10.4f} "
              f"{r['w_dipolar']:<12.4f} {r['n90']:<8} {r['n95']:<8}")

    # Analysis
    nuc = results['nuclear_only']
    ee = results['ee_only']
    full = results['full']

    print("\n--- Analysis ---")
    print(f"  Nuclear spreading: w(0,0) = {nuc['w00']:.4f} "
          f"({nuc['n90']} ch for 90%)")
    print(f"  E-E spreading:     w(0,0) = {ee['w00']:.4f} "
          f"({ee['n90']} ch for 90%)")
    print(f"  Full spreading:    w(0,0) = {full['w00']:.4f} "
          f"({full['n90']} ch for 90%)")

    nuc_spreading = 1.0 - nuc['w00']
    ee_spreading = 1.0 - ee['w00']
    full_spreading = 1.0 - full['w00']

    if full_spreading > 0.01:
        nuc_frac = nuc_spreading / full_spreading * 100
        ee_frac = ee_spreading / full_spreading * 100
        print(f"\n  Nuclear contributes ~{nuc_frac:.0f}% of total spreading")
        print(f"  E-E contributes ~{ee_frac:.0f}% of total spreading")

        if nuc_frac > 70:
            print("  => NUCLEAR ASYMMETRY dominates channel spreading")
            print("  => Higher l_max SHOULD improve convergence")
        elif ee_frac > 70:
            print("  => E-E CORRELATION dominates channel spreading")
            print("  => Higher l_max has diminishing returns")
        else:
            print("  => Both nuclear and e-e contribute significantly")

    return results


# --------------------------------------------------------------------------- #
# Part B: Higher l_max convergence extrapolation
# --------------------------------------------------------------------------- #

def run_part_b():
    """Higher l_max convergence extrapolation for 6:1 system."""
    print("\n\n" + "=" * 80)
    print("PART B: HIGHER l_max CONVERGENCE EXTRAPOLATION")
    print("System: 6:1 (Z_A=6, Z_B=1), charge-center origin, rho=0.7")
    print("=" * 80)

    Z_A, Z_B = 6.0, 1.0
    rho = 0.7

    # Determine feasible l_max range
    # Matrix size = (l_max+1)^2 * N_alpha
    # l_max=10: 121 * N_alpha; at N_alpha=200, size = 24200 (feasible with sparse)
    # l_max=8: 81 * 200 = 16200 (feasible dense)
    # Try up to l_max=10 with N_alpha=200

    N_alpha = 200  # Reduced for computational feasibility
    l_max_values = list(range(0, 11))  # 0 through 10

    convergence = []
    all_weights = {}

    for l_max in l_max_values:
        channels = build_channel_list(l_max, homonuclear=False)
        n_ch = len(channels)
        N_total = n_ch * N_alpha

        t0 = time.time()
        print(f"  l_max={l_max:>2}: {n_ch:>4} ch, matrix {N_total:>6}x{N_total:<6} ",
              end='', flush=True)

        # Check feasibility
        if N_total > 60000:
            print("SKIPPED (too large)")
            continue

        r = solve_alpha_decomposed(
            channels, rho, Z_A, Z_B, N_alpha=N_alpha, mode='full')
        dt = time.time() - t0

        stats = extract_weight_stats(channels, r['weights'])
        mu0 = r['mu0']
        convergence.append({
            'l_max': l_max,
            'n_ch': n_ch,
            'mu0': mu0,
            'w00': stats['w00'],
            'w_dipolar': stats['w_dipolar'],
            'n90': stats['n90'],
            'n95': stats['n95'],
            'n99': stats['n99'],
            'time': dt,
        })

        # Store weights for comparison
        all_weights[l_max] = {
            str(channels[i]): float(r['weights'][i])
            for i in range(n_ch)
        }

        print(f"mu0={mu0:>12.6f}  w(0,0)={stats['w00']:.4f}  [{dt:.1f}s]")

    # Compute convergence metrics
    mu_vals = [c['mu0'] for c in convergence]
    l_vals = [c['l_max'] for c in convergence]

    print("\n" + "=" * 100)
    print("PART B SUMMARY TABLE")
    print("=" * 100)
    header = (f"{'l_max':<7} {'N_ch':<7} {'mu0':<14} {'|Delta_mu|':<14} "
              f"{'Est mu(inf)':<14} {'w(0,0)':<10} {'N(90%)':<8} {'N(95%)':<8}")
    print(header)
    print("-" * len(header))

    for idx, c in enumerate(convergence):
        delta_mu = abs(c['mu0'] - convergence[idx-1]['mu0']) if idx > 0 else float('nan')

        # Richardson extrapolation using last 3 points up to current
        if idx >= 2:
            l_arr = np.array([convergence[k]['l_max'] for k in range(idx+1)])
            mu_arr = np.array([convergence[k]['mu0'] for k in range(idx+1)])
            mu_inf = richardson_extrapolate(l_arr, mu_arr)
        else:
            mu_inf = float('nan')

        c['delta_mu'] = delta_mu
        c['mu_inf'] = mu_inf

        print(f"{c['l_max']:<7} {c['n_ch']:<7} {c['mu0']:<14.8f} "
              f"{delta_mu:<14.2e} {mu_inf:<14.8f} {c['w00']:<10.4f} "
              f"{c['n90']:<8} {c['n95']:<8}")

    # Convergence analysis
    if len(convergence) >= 4:
        deltas = [abs(convergence[i]['mu0'] - convergence[i-1]['mu0'])
                  for i in range(1, len(convergence))]
        l_deltas = [convergence[i]['l_max'] for i in range(1, len(convergence))]

        # Check if convergence is exponential or algebraic
        # Fit log(delta) vs l_max (exponential) and log(delta) vs log(l_max) (algebraic)
        valid = [(l, d) for l, d in zip(l_deltas, deltas) if d > 1e-15 and l > 0]
        if len(valid) >= 3:
            l_fit = np.array([v[0] for v in valid])
            d_fit = np.array([v[1] for v in valid])

            # Exponential: log(d) = a - b*l => d ~ exp(-b*l)
            try:
                exp_fit = np.polyfit(l_fit, np.log(d_fit), 1)
                exp_rate = -exp_fit[0]
            except Exception:
                exp_rate = 0.0

            # Algebraic: log(d) = a - p*log(l) => d ~ l^(-p)
            try:
                alg_fit = np.polyfit(np.log(l_fit), np.log(d_fit), 1)
                alg_power = -alg_fit[0]
            except Exception:
                alg_power = 0.0

            print(f"\n--- Convergence Analysis ---")
            print(f"  Exponential fit: |Delta_mu| ~ exp(-{exp_rate:.3f} * l_max)")
            print(f"  Algebraic fit:   |Delta_mu| ~ l_max^(-{alg_power:.2f})")

            # R-squared for each fit
            log_d = np.log(d_fit)
            exp_pred = np.polyval(exp_fit, l_fit)
            alg_pred = np.polyval(alg_fit, np.log(l_fit))
            ss_tot = np.sum((log_d - np.mean(log_d))**2)
            if ss_tot > 0:
                r2_exp = 1 - np.sum((log_d - exp_pred)**2) / ss_tot
                r2_alg = 1 - np.sum((log_d - alg_pred)**2) / ss_tot
                print(f"  R² (exponential): {r2_exp:.4f}")
                print(f"  R² (algebraic):   {r2_alg:.4f}")

                if r2_exp > r2_alg + 0.05:
                    print("  => Convergence is EXPONENTIAL — higher l_max is efficient")
                elif r2_alg > r2_exp + 0.05:
                    print("  => Convergence is ALGEBRAIC — higher l_max has diminishing returns")
                else:
                    print("  => Convergence type is ambiguous in this range")

    # Weight distribution evolution
    if len(convergence) >= 2:
        print(f"\n--- Weight Distribution Evolution ---")
        print(f"  {'l_max':<7} {'w(0,0)':<10} {'w(dipolar)':<12} "
              f"{'N(90%)/N_ch':<14} {'N(95%)/N_ch':<14}")
        for c in convergence:
            ratio90 = c['n90'] / c['n_ch'] if c['n_ch'] > 0 else 0
            ratio95 = c['n95'] / c['n_ch'] if c['n_ch'] > 0 else 0
            print(f"  {c['l_max']:<7} {c['w00']:<10.4f} {c['w_dipolar']:<12.4f} "
                  f"{ratio90:<14.4f} {ratio95:<14.4f}")

        # Does the proportional spread narrow?
        first_ratio = convergence[1]['n90'] / convergence[1]['n_ch'] if convergence[1]['n_ch'] > 0 else 1
        last_ratio = convergence[-1]['n90'] / convergence[-1]['n_ch'] if convergence[-1]['n_ch'] > 0 else 1

        if last_ratio < first_ratio * 0.8:
            print("  => Weight distribution NARROWS proportionally — convergence is efficient")
        elif last_ratio > first_ratio * 1.2:
            print("  => Weight distribution SPREADS proportionally — more channels keep filling")
        else:
            print("  => Weight distribution is roughly stable proportionally")

    return convergence, all_weights


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def make_plots(part_a_results, convergence, all_weights):
    """Generate all diagnostic plots."""

    # ---- Plot 1: Part A — channel weight comparison (3 panels) ---- #
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    modes_order = ['nuclear_only', 'ee_only', 'full']
    mode_labels = {
        'nuclear_only': 'Nuclear Only',
        'ee_only': 'E-E Only (diag nuc)',
        'full': 'Full Coupling',
    }

    for idx, mode in enumerate(modes_order):
        ax = axes[idx]
        r = part_a_results[mode]
        channels = r['channels']
        weights = r['weights']
        labels = [f'({l1},{l2})' for l1, l2 in channels]

        colors = ['steelblue' if w > 0.01 else 'lightcoral' for w in weights]
        ax.bar(range(len(channels)), weights,
               color=colors, alpha=0.8, edgecolor='black', linewidth=0.3)
        ax.set_xticks(range(0, len(channels), max(1, len(channels)//15)))
        ax.set_xticklabels(
            [labels[i] for i in range(0, len(channels), max(1, len(channels)//15))],
            rotation=90, fontsize=6)
        ax.set_ylabel('Channel weight')
        ax.set_title(f'{mode_labels[mode]}\nmu0={r["mu0"]:.4f}, w(0,0)={r["w00"]:.4f}')
        ax.set_ylim(0, 1.05)
        ax.axhline(0.01, color='red', ls=':', alpha=0.5)
        ax.grid(axis='y', alpha=0.3)

    plt.suptitle('Part A: Nuclear vs E-E Channel Weight Decomposition\n'
                 '(6:1 system, rho=0.7, l_max=4, charge-center origin)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'nuc_vs_ee_channel_weights.png', dpi=150)
    plt.close()
    print(f"\nSaved: {PLOT_DIR / 'nuc_vs_ee_channel_weights.png'}")

    # ---- Plot 2: mu0 vs l_max ---- #
    if len(convergence) >= 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        l_vals = [c['l_max'] for c in convergence]
        mu_vals = [c['mu0'] for c in convergence]

        ax1.plot(l_vals, mu_vals, 'o-', color='red', markersize=8, linewidth=2)
        ax1.set_xlabel('l_max', fontsize=12)
        ax1.set_ylabel('Ground-state eigenvalue mu_0', fontsize=12)
        ax1.set_title('mu_0 vs l_max (6:1 system)', fontsize=13)
        ax1.grid(True, alpha=0.3)

        # Richardson extrapolation estimate
        mu_inf_vals = [c.get('mu_inf', np.nan) for c in convergence]
        valid_inf = [(l, m) for l, m in zip(l_vals, mu_inf_vals) if not np.isnan(m)]
        if valid_inf:
            final_inf = valid_inf[-1][1]
            ax1.axhline(final_inf, color='green', ls='--', alpha=0.7,
                       label=f'Richardson est: {final_inf:.6f}')
            ax1.legend(fontsize=10)

        # Delta_mu vs l_max (log scale)
        if len(convergence) >= 2:
            deltas = [abs(convergence[i]['mu0'] - convergence[i-1]['mu0'])
                      for i in range(1, len(convergence))]
            l_deltas = [convergence[i]['l_max'] for i in range(1, len(convergence))]

            ax2.semilogy(l_deltas, deltas, 'o-', color='red',
                        markersize=8, linewidth=2, label='|Delta mu|')

            # Fit lines
            valid = [(l, d) for l, d in zip(l_deltas, deltas)
                     if d > 1e-15 and l > 0]
            if len(valid) >= 3:
                l_fit = np.array([v[0] for v in valid])
                d_fit = np.array([v[1] for v in valid])

                # Exponential fit
                try:
                    p = np.polyfit(l_fit, np.log(d_fit), 1)
                    l_plot = np.linspace(l_fit[0], l_fit[-1] + 2, 50)
                    ax2.semilogy(l_plot, np.exp(np.polyval(p, l_plot)),
                                '--', color='blue', alpha=0.7,
                                label=f'Exp fit: rate={-p[0]:.3f}')
                except Exception:
                    pass

                # Algebraic fit
                try:
                    p2 = np.polyfit(np.log(l_fit), np.log(d_fit), 1)
                    ax2.semilogy(l_plot, np.exp(np.polyval(p2, np.log(l_plot))),
                                ':', color='green', alpha=0.7,
                                label=f'Power fit: p={-p2[0]:.2f}')
                except Exception:
                    pass

            ax2.set_xlabel('l_max', fontsize=12)
            ax2.set_ylabel('|mu_0(l_max) - mu_0(l_max-1)|', fontsize=12)
            ax2.set_title('Convergence Rate (6:1 system)', fontsize=13)
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)

        plt.suptitle('Part B: Higher l_max Convergence (6:1 system, rho=0.7)',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(PLOT_DIR / 'lmax_convergence_extrapolation.png', dpi=150)
        plt.close()
        print(f"Saved: {PLOT_DIR / 'lmax_convergence_extrapolation.png'}")

    # ---- Plot 3: Weight distribution at multiple l_max values ---- #
    if len(all_weights) >= 2:
        # Pick representative l_max values
        available = sorted(all_weights.keys())
        show_lmax = [l for l in [0, 2, 4, 6, 8, 10] if l in available]
        if not show_lmax:
            show_lmax = available[:4]

        n_panels = len(show_lmax)
        fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))
        if n_panels == 1:
            axes = [axes]

        for idx, l_max in enumerate(show_lmax):
            ax = axes[idx]
            w_dict = all_weights[l_max]
            channels_sorted = sorted(w_dict.keys(),
                                     key=lambda x: w_dict[x], reverse=True)
            weights_sorted = [w_dict[c] for c in channels_sorted]

            n_show = min(20, len(channels_sorted))
            colors = ['steelblue' if w > 0.01 else 'lightcoral'
                      for w in weights_sorted[:n_show]]
            ax.bar(range(n_show), weights_sorted[:n_show],
                   color=colors, alpha=0.8, edgecolor='black', linewidth=0.3)
            ax.set_xticks(range(n_show))
            ax.set_xticklabels(channels_sorted[:n_show], rotation=90, fontsize=5)
            ax.set_ylabel('Weight')
            ax.set_title(f'l_max={l_max} ({len(w_dict)} ch)')
            ax.set_ylim(0, 1.05)
            ax.axhline(0.01, color='red', ls=':', alpha=0.5)
            ax.grid(axis='y', alpha=0.3)

        plt.suptitle('Part B: Weight Distribution Evolution (top 20 channels)\n'
                     '(6:1 system, rho=0.7, charge-center origin)',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(PLOT_DIR / 'lmax_weight_evolution.png', dpi=150)
        plt.close()
        print(f"Saved: {PLOT_DIR / 'lmax_weight_evolution.png'}")

    # ---- Plot 4: Richardson extrapolation visualization ---- #
    if len(convergence) >= 3:
        fig, ax = plt.subplots(figsize=(10, 6))
        l_vals = [c['l_max'] for c in convergence]
        mu_vals = [c['mu0'] for c in convergence]
        mu_inf_vals = [c.get('mu_inf', np.nan) for c in convergence]

        ax.plot(l_vals, mu_vals, 'o-', color='red', markersize=8,
                linewidth=2, label='mu_0(l_max)')
        valid_inf = [(l, m) for l, m in zip(l_vals, mu_inf_vals) if not np.isnan(m)]
        if valid_inf:
            ax.plot([v[0] for v in valid_inf], [v[1] for v in valid_inf],
                    's--', color='blue', markersize=6, alpha=0.7,
                    label='Richardson estimate mu(inf)')

        ax.set_xlabel('l_max', fontsize=12)
        ax.set_ylabel('mu_0', fontsize=12)
        ax.set_title('Richardson Extrapolation of mu_0(l_max -> inf)\n'
                     '(6:1 system, rho=0.7)', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / 'richardson_extrapolation.png', dpi=150)
        plt.close()
        print(f"Saved: {PLOT_DIR / 'richardson_extrapolation.png'}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    t_start = time.time()

    # Part A
    part_a_results = run_part_a()

    # Part B
    convergence, all_weights = run_part_b()

    # Plots
    print("\n\nGenerating plots...")
    make_plots(part_a_results, convergence, all_weights)

    # Save all data
    save_data = {
        'part_a': {},
        'part_b': {
            'convergence': convergence,
        },
        'metadata': {
            'Z_A': 6.0, 'Z_B': 1.0,
            'rho': 0.7,
            'origin': 'charge_center',
            'date': '2026-03-26',
        },
    }

    for mode in ['nuclear_only', 'ee_only', 'full']:
        r = part_a_results[mode]
        save_data['part_a'][mode] = {
            'label': r['label'],
            'mu0': r['mu0'],
            'w00': r['w00'],
            'w_dipolar': r['w_dipolar'],
            'n90': r['n90'],
            'n95': r['n95'],
            'n99': r['n99'],
            'weights': {
                str(r['channels'][i]): float(r['weights'][i])
                for i in range(r['n_ch'])
            },
        }

    # Store weight evolution from Part B
    save_data['part_b']['weights_by_lmax'] = {
        str(k): v for k, v in all_weights.items()
    }

    data_path = DATA_DIR / 'nuclear_vs_ee_decomposition.json'
    with open(data_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nData saved: {data_path}")

    dt_total = time.time() - t_start
    print(f"\nTotal runtime: {dt_total:.1f}s")

    # Final verdict
    print("\n" + "=" * 80)
    print("FINAL VERDICT")
    print("=" * 80)

    nuc = part_a_results['nuclear_only']
    ee = part_a_results['ee_only']
    full = part_a_results['full']

    nuc_spread = 1.0 - nuc['w00']
    ee_spread = 1.0 - ee['w00']
    full_spread = 1.0 - full['w00']

    print(f"\n  Channel spreading (1 - w(0,0)):")
    print(f"    Nuclear only:  {nuc_spread:.4f}")
    print(f"    E-E only:      {ee_spread:.4f}")
    print(f"    Full:          {full_spread:.4f}")

    if full_spread > 0.01:
        if nuc_spread > 0.7 * full_spread:
            print(f"\n  => NUCLEAR ASYMMETRY is the primary source of spreading")
            print(f"     ({nuc_spread/full_spread*100:.0f}% of full spreading)")
            print(f"  => Strategy: higher l_max SHOULD converge the eigenvalue")
        else:
            print(f"\n  => E-E CORRELATION contributes substantially to spreading")
            print(f"  => Higher l_max has limited benefit beyond capturing nuclear part")

    if len(convergence) >= 3:
        last = convergence[-1]
        mu_inf = last.get('mu_inf', np.nan)
        if not np.isnan(mu_inf):
            gap = abs(last['mu0'] - mu_inf) / abs(mu_inf) * 100
            print(f"\n  Highest computed l_max: {last['l_max']}")
            print(f"  mu_0 at l_max={last['l_max']}: {last['mu0']:.8f}")
            print(f"  Richardson estimate mu(inf): {mu_inf:.8f}")
            print(f"  Remaining gap: {gap:.2f}%")

            if gap < 1.0:
                print(f"  => Eigenvalue is WELL CONVERGED at l_max={last['l_max']}")
            elif gap < 5.0:
                print(f"  => Eigenvalue is REASONABLY converged, ~{gap:.1f}% gap")
            else:
                print(f"  => Eigenvalue is NOT yet converged ({gap:.1f}% gap)")
                print(f"     Higher l_max or different basis needed")


if __name__ == '__main__':
    main()
