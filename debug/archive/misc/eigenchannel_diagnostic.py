#!/usr/bin/env python3
"""
Eigenchannel Diagnostic: Nuclear coupling eigenspectrum analysis.

Computes the nuclear coupling matrix C_nuc in the (l1, l2) channel basis
and analyzes its eigenvalue spectrum as a function of charge asymmetry.
Determines whether eigenchannel rotation can improve convergence for
asymmetric bonds (Paper 15, Eq. 34).

All matrix elements are algebraic (Gaunt integrals via Wigner 3j symbols).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sympy.physics.wigner import wigner_3j
from typing import List, Tuple, Dict
import json
import os

# --------------------------------------------------------------------------- #
# Channel basis construction
# --------------------------------------------------------------------------- #

def build_channel_list(l_max: int, homonuclear: bool) -> List[Tuple[int, int]]:
    """Build (l1, l2) channel list for sigma (m=0) states."""
    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue  # gerade constraint
            channels.append((l1, l2))
    channels.sort(key=lambda c: (c[0] + c[1], c[0]))
    return channels


# --------------------------------------------------------------------------- #
# Gaunt integral via Wigner 3j (exact, algebraic)
# --------------------------------------------------------------------------- #

_gaunt_cache: Dict[Tuple[int, int, int], float] = {}

def gaunt(lp: int, k: int, l: int) -> float:
    """
    Gaunt integral G(l', k, l) = integral of P_{l'}(x) P_k(x) P_l(x) dx
    over [-1, 1].

    Equals 2 * (l' k l; 0 0 0)^2 via Wigner 3j symbol.
    Selection rules: l'+k+l even, |l'-l| <= k <= l'+l.
    """
    key = (lp, k, l)
    if key in _gaunt_cache:
        return _gaunt_cache[key]

    # Selection rules
    if (lp + k + l) % 2 != 0:
        _gaunt_cache[key] = 0.0
        return 0.0
    if k < abs(lp - l) or k > lp + l:
        _gaunt_cache[key] = 0.0
        return 0.0

    w3j = float(wigner_3j(lp, k, l, 0, 0, 0))
    val = 2.0 * w3j**2
    _gaunt_cache[key] = val
    return val


# --------------------------------------------------------------------------- #
# Radial factor f_k (split-region Legendre expansion)
# --------------------------------------------------------------------------- #

def f_k(s: float, rho: float, k: int) -> float:
    """
    f_k(s, rho) = (min(s,rho)/max(s,rho))^k / max(s,rho)

    Split-region: converges for both s > rho and s < rho.
    """
    if s <= 0.0 or rho <= 0.0:
        return 0.0
    mn = min(s, rho)
    mx = max(s, rho)
    return (mn / mx)**k / mx


# --------------------------------------------------------------------------- #
# Nuclear coupling matrix element
# --------------------------------------------------------------------------- #

def nuclear_coupling_element(
    lp: int, l: int,
    s: float,
    Z_A: float, rho_A: float,
    Z_B: float, rho_B: float,
) -> float:
    """
    Compute single-electron nuclear coupling: integral of
    P_{l'}(x) V_nuc(x; s, rho) P_l(x) dx

    V_nuc = -Z_A sum_k f_k(s, rho_A) P_k(x)
            -Z_B sum_k (-1)^k f_k(s, rho_B) P_k(x)

    Uses Gaunt integral G(l', k, l) for exact angular integration.
    """
    k_min = abs(lp - l)
    k_max = lp + l
    # Parity: l' + k + l must be even, so k has same parity as l'+l
    required_parity = (lp + l) % 2

    val = 0.0
    for k in range(k_min, k_max + 1):
        if k % 2 != required_parity:
            continue
        g = gaunt(lp, k, l)
        if g == 0.0:
            continue
        c_A = -Z_A * f_k(s, rho_A, k)
        c_B = -Z_B * ((-1)**k) * f_k(s, rho_B, k)
        val += g * (c_A + c_B)

    # Normalization: sqrt((2l'+1)(2l+1))
    norm = np.sqrt((2 * lp + 1) * (2 * l + 1))
    return norm * val


# --------------------------------------------------------------------------- #
# Full coupling matrix
# --------------------------------------------------------------------------- #

def build_coupling_matrix(
    channels: List[Tuple[int, int]],
    alpha: float,
    rho: float,
    Z_A: float, Z_B: float,
) -> np.ndarray:
    """
    Build the nuclear coupling matrix C_nuc in channel basis at fixed (alpha, rho).

    C_nuc[(l1',l2'), (l1,l2)] =
        delta_{l2',l2} * N1 * integral(P_{l1'} V_1 P_{l1} dx)
      + delta_{l1',l1} * N2 * integral(P_{l2'} V_2 P_{l2} dx)

    Charge-center origin: z0 = R(Z_A - Z_B) / (2(Z_A + Z_B))
    rho_A = rho * Z_B / (Z_A + Z_B) * 2  [distance to A from charge center]
    rho_B = rho * Z_A / (Z_A + Z_B) * 2  [distance to B from charge center]

    Actually: rho = R/(2*R_e), z0 = R*(Z_A-Z_B)/(2*(Z_A+Z_B))
    rho_A = R_A/R_e where R_A = R/2 - z0 = R*Z_B/(Z_A+Z_B)
    rho_B = R_B/R_e where R_B = R/2 + z0 = R*Z_A/(Z_A+Z_B)

    So rho_A/rho = 2*Z_B/(Z_A+Z_B) and rho_B/rho = 2*Z_A/(Z_A+Z_B).
    """
    Z_tot = Z_A + Z_B
    rho_A = rho * 2.0 * Z_B / Z_tot  # charge-center: A closer if Z_A > Z_B
    rho_B = rho * 2.0 * Z_A / Z_tot

    s1 = np.cos(alpha)
    s2 = np.sin(alpha)

    n_ch = len(channels)
    C = np.zeros((n_ch, n_ch))

    for i, (l1p, l2p) in enumerate(channels):
        for j, (l1, l2) in enumerate(channels):
            val = 0.0

            # Electron 1 coupling: delta_{l2', l2}
            if l2p == l2:
                val += nuclear_coupling_element(l1p, l1, s1, Z_A, rho_A, Z_B, rho_B)

            # Electron 2 coupling: delta_{l1', l1}
            if l1p == l1:
                val += nuclear_coupling_element(l2p, l2, s2, Z_A, rho_A, Z_B, rho_B)

            C[i, j] = val

    return C


# --------------------------------------------------------------------------- #
# Eigenvalue analysis
# --------------------------------------------------------------------------- #

def analyze_spectrum(C: np.ndarray) -> Dict:
    """Compute eigenvalues and cumulative Frobenius fraction."""
    eigenvalues = np.linalg.eigvalsh(C)
    sorted_eig = np.sort(np.abs(eigenvalues))[::-1]

    frob_norm = np.linalg.norm(C, 'fro')
    if frob_norm < 1e-15:
        return {
            'eigenvalues': eigenvalues,
            'sorted_abs': sorted_eig,
            'frob_norm': 0.0,
            'cumulative_fraction': np.ones(len(sorted_eig)),
            'n_90': 1, 'n_95': 1, 'n_99': 1,
        }

    # Cumulative fraction of Frobenius norm captured by top N eigenchannels
    # Frobenius norm squared = sum of eigenvalue^2 (for symmetric matrix)
    eig_sq = sorted_eig**2
    cum_sq = np.cumsum(eig_sq)
    total_sq = np.sum(eig_sq)
    cum_frac = np.sqrt(cum_sq / total_sq) if total_sq > 0 else np.ones(len(eig_sq))

    n_90 = int(np.searchsorted(cum_frac, 0.90) + 1)
    n_95 = int(np.searchsorted(cum_frac, 0.95) + 1)
    n_99 = int(np.searchsorted(cum_frac, 0.99) + 1)

    return {
        'eigenvalues': eigenvalues,
        'sorted_abs': sorted_eig,
        'frob_norm': frob_norm,
        'cumulative_fraction': cum_frac,
        'n_90': min(n_90, len(sorted_eig)),
        'n_95': min(n_95, len(sorted_eig)),
        'n_99': min(n_99, len(sorted_eig)),
    }


# --------------------------------------------------------------------------- #
# Main diagnostic
# --------------------------------------------------------------------------- #

def run_diagnostic():
    l_max = 4

    systems = [
        ('H2 (1:1)', 1.0, 1.0, True),
        ('HeH+ (2:1)', 2.0, 1.0, False),
        ('Hypo (6:1)', 6.0, 1.0, False),
    ]

    alpha_val = np.pi / 6  # asymmetric electron config
    rho_vals = [0.5, 1.0, 1.5]

    results = {}

    print("=" * 80)
    print("EIGENCHANNEL DIAGNOSTIC: Nuclear Coupling Eigenspectrum")
    print(f"l_max = {l_max}, m = 0 (sigma only), alpha = pi/6")
    print("=" * 80)

    for name, Z_A, Z_B, homonuclear in systems:
        channels = build_channel_list(l_max, homonuclear)
        n_ch = len(channels)
        ratio = Z_A / Z_B
        print(f"\n{'-' * 60}")
        print(f"System: {name}  (Z_A={Z_A}, Z_B={Z_B}, {n_ch} channels)")
        print(f"Channels: {channels}")
        print(f"{'-' * 60}")

        results[name] = {}

        for rho in rho_vals:
            C = build_coupling_matrix(channels, alpha_val, rho, Z_A, Z_B)
            analysis = analyze_spectrum(C)

            results[name][f'rho={rho}'] = {
                'n_channels': n_ch,
                'rho': rho,
                'frob_norm': float(analysis['frob_norm']),
                'n_90': analysis['n_90'],
                'n_95': analysis['n_95'],
                'n_99': analysis['n_99'],
                'eigenvalues': analysis['eigenvalues'].tolist(),
                'cumulative_fraction': analysis['cumulative_fraction'].tolist(),
            }

            print(f"\n  rho = {rho}:")
            print(f"    Frobenius norm:  {analysis['frob_norm']:.6f}")
            print(f"    Top 5 |eig|:     {analysis['sorted_abs'][:5]}")
            print(f"    N for 90%:       {analysis['n_90']}")
            print(f"    N for 95%:       {analysis['n_95']}")
            print(f"    N for 99%:       {analysis['n_99']}")

    # ---- Summary table ---- #
    print("\n\n" + "=" * 80)
    print("SUMMARY: Eigenchannels needed to capture Frobenius norm fraction")
    print("=" * 80)
    header = f"{'System':<20} {'rho':>5} {'N_ch':>5} {'||C||_F':>10} {'90%':>5} {'95%':>5} {'99%':>5}"
    print(header)
    print("-" * len(header))
    for name, _, _, _ in systems:
        for rho in rho_vals:
            r = results[name][f'rho={rho}']
            print(f"{name:<20} {rho:>5.1f} {r['n_channels']:>5} "
                  f"{r['frob_norm']:>10.4f} {r['n_90']:>5} {r['n_95']:>5} {r['n_99']:>5}")

    # ---- Plots ---- #
    plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Plot 1: Eigenvalue spectra for each system at rho=1.0
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for idx, (name, Z_A, Z_B, homonuclear) in enumerate(systems):
        ax = axes[idx]
        r = results[name]['rho=1.0']
        eigs = np.sort(np.abs(r['eigenvalues']))[::-1]
        ax.bar(range(1, len(eigs) + 1), eigs, color='steelblue', alpha=0.8)
        ax.set_xlabel('Eigenchannel index')
        ax.set_ylabel('|eigenvalue|')
        ax.set_title(f'{name}, rho=1.0')
        ax.set_yscale('log')
        # Set floor for log scale
        nonzero = eigs[eigs > 0]
        if len(nonzero) > 0:
            ax.set_ylim(bottom=nonzero[-1] * 0.1)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'eigenchannel_spectra_rho1.png'), dpi=150)
    plt.close()
    print(f"\nSaved: {os.path.join(plot_dir, 'eigenchannel_spectra_rho1.png')}")

    # Plot 2: Cumulative fraction for all (system, rho) combinations
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = {'H2 (1:1)': 'blue', 'HeH+ (2:1)': 'orange', 'Hypo (6:1)': 'red'}
    for idx, rho in enumerate(rho_vals):
        ax = axes[idx]
        for name, _, _, _ in systems:
            r = results[name][f'rho={rho}']
            cf = r['cumulative_fraction']
            ax.plot(range(1, len(cf) + 1), cf, 'o-',
                    color=colors[name], label=name, markersize=5)
        ax.axhline(0.90, color='gray', linestyle='--', alpha=0.5, label='90%')
        ax.axhline(0.95, color='gray', linestyle=':', alpha=0.5, label='95%')
        ax.axhline(0.99, color='gray', linestyle='-.', alpha=0.5, label='99%')
        ax.set_xlabel('Number of eigenchannels')
        ax.set_ylabel('Fraction of ||C||_F captured')
        ax.set_title(f'rho = {rho}')
        ax.legend(fontsize=8)
        ax.set_ylim(0.0, 1.05)
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'eigenchannel_convergence.png'), dpi=150)
    plt.close()
    print(f"Saved: {os.path.join(plot_dir, 'eigenchannel_convergence.png')}")

    # Plot 3: Cumulative fraction comparison at rho=1.0 (key diagnostic)
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, _, _, _ in systems:
        r = results[name]['rho=1.0']
        cf = r['cumulative_fraction']
        n = range(1, len(cf) + 1)
        ax.plot(n, cf, 'o-', color=colors[name], label=name,
                markersize=7, linewidth=2)
    ax.axhline(0.90, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(0.95, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(0.99, color='gray', linestyle='-.', alpha=0.5)
    ax.set_xlabel('Number of eigenchannels retained', fontsize=12)
    ax.set_ylabel('Fraction of ||C||_F captured', fontsize=12)
    ax.set_title('Nuclear Coupling Eigenchannel Convergence (rho=1.0, alpha=pi/6)',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'eigenchannel_key_diagnostic.png'), dpi=150)
    plt.close()
    print(f"Saved: {os.path.join(plot_dir, 'eigenchannel_key_diagnostic.png')}")

    # Plot 4: Heatmap of coupling matrix at rho=1.0 for each system
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for idx, (name, Z_A, Z_B, homonuclear) in enumerate(systems):
        ax = axes[idx]
        channels = build_channel_list(l_max, homonuclear)
        C = build_coupling_matrix(channels, alpha_val, 1.0, Z_A, Z_B)
        vmax = np.max(np.abs(C))
        im = ax.imshow(C, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title(f'{name}', fontsize=12)
        labels = [f'({l1},{l2})' for l1, l2 in channels]
        ax.set_xticks(range(len(channels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticks(range(len(channels)))
        ax.set_yticklabels(labels, fontsize=6)
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle('Nuclear Coupling Matrix C_nuc (rho=1.0, alpha=pi/6)', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'eigenchannel_coupling_matrices.png'), dpi=150)
    plt.close()
    print(f"Saved: {os.path.join(plot_dir, 'eigenchannel_coupling_matrices.png')}")

    # Plot 5: Eigenvalue decay rate comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, _, _, homonuclear in systems:
        r = results[name]['rho=1.0']
        eigs = np.sort(np.abs(r['eigenvalues']))[::-1]
        nonzero = eigs[eigs > 1e-15]
        if len(nonzero) > 1:
            normalized = nonzero / nonzero[0]
            ax.semilogy(range(1, len(normalized) + 1), normalized, 'o-',
                        color=colors[name], label=name, markersize=6, linewidth=2)
    ax.set_xlabel('Eigenchannel index', fontsize=12)
    ax.set_ylabel('|eigenvalue| / |eigenvalue_max|', fontsize=12)
    ax.set_title('Normalized Eigenvalue Decay (rho=1.0)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'eigenchannel_decay_rate.png'), dpi=150)
    plt.close()
    print(f"Saved: {os.path.join(plot_dir, 'eigenchannel_decay_rate.png')}")

    # Save data
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(data_dir, exist_ok=True)
    data_path = os.path.join(data_dir, 'eigenchannel_diagnostic.json')
    with open(data_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved data: {data_path}")


if __name__ == '__main__':
    run_diagnostic()
