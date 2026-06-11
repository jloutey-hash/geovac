#!/usr/bin/env python3
"""
Angular Wavefunction Channel Decomposition Diagnostic.

Directly solves the coupled-channel angular eigenvalue problem (Paper 15, Eq. 30)
at fixed rho to determine whether the partial-wave (l1, l2) expansion of the
angular *wavefunction* is the convergence bottleneck for asymmetric bonds.

Context: Task 1 showed eigenchannel rotation doesn't help (flat eigenvalue spectrum).
Task 2 showed Legendre expansion of the potential converges exponentially (l_max=2
represents the potential to <0.1% for 6:1 with charge-center origin). Therefore the
17% geometry error for 6:1 bonds must come from slow wavefunction convergence in
the channel basis, not from potential truncation.

Equations (Paper 15):
  Single-channel (l1=l2=0):
    -1/2 u''(alpha) + V_eff(alpha; rho) u(alpha) = mu u(alpha)
  with:
    V_eff(alpha) = -2 + l1(l1+1)/(2 cos^2 alpha) + l2(l2+1)/(2 sin^2 alpha)
                   + R_e * C_mol^{l1,l2}(alpha; rho)
  Multichannel:
    -1/2 u_i''(alpha) + sum_j W_ij(alpha; rho) u_j(alpha) = mu u_i(alpha)
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
from typing import List, Tuple, Dict

# Reuse algebraic coupling machinery from eigenchannel_diagnostic
sys.path.insert(0, str(Path(__file__).parent))
from eigenchannel_diagnostic import (
    gaunt, build_channel_list,
)

PLOT_DIR = Path(__file__).parent / "plots"
DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Vectorized f_k: operates on arrays
# --------------------------------------------------------------------------- #

def f_k_vec(s: np.ndarray, rho: float, k: int) -> np.ndarray:
    """Vectorized f_k(s, rho) = (min/max)^k / max."""
    s = np.asarray(s, dtype=float)
    result = np.zeros_like(s)
    mask = s > 0
    mn = np.minimum(s[mask], rho)
    mx = np.maximum(s[mask], rho)
    result[mask] = (mn / mx)**k / mx
    return result


# --------------------------------------------------------------------------- #
# E-e coupling element (sigma channels, m=0) — vectorized
# --------------------------------------------------------------------------- #

def ee_coupling_element(
    l1p: int, l2p: int, l1: int, l2: int,
    cos_a: np.ndarray, sin_a: np.ndarray,
) -> np.ndarray:
    """
    Electron-electron repulsion coupling between (l1', l2') and (l1, l2) channels.
    Vectorized over alpha grid.
    """
    min_sc = np.minimum(cos_a, sin_a)
    max_sc = np.maximum(cos_a, sin_a)

    W = np.zeros(len(cos_a))

    k_min = max(abs(l1p - l1), abs(l2p - l2))
    k_max_val = min(l1p + l1, l2p + l2)

    for k in range(k_min, k_max_val + 1):
        g1 = gaunt(l1p, k, l1)
        g2 = gaunt(l2p, k, l2)
        if abs(g1) < 1e-15 or abs(g2) < 1e-15:
            continue
        fk = (min_sc / max_sc) ** k
        W += g1 * g2 * fk / max_sc

    norm = np.sqrt((2*l1p+1) * (2*l2p+1) * (2*l1+1) * (2*l2+1)) / 4.0
    return norm * W


# --------------------------------------------------------------------------- #
# Nuclear coupling for one electron — vectorized
# --------------------------------------------------------------------------- #

def nuclear_coupling_vec(
    lp: int, l: int,
    s: np.ndarray,
    Z_A: float, rho_A: float,
    Z_B: float, rho_B: float,
) -> np.ndarray:
    """Vectorized nuclear coupling integral for one electron."""
    k_min = abs(lp - l)
    k_max = lp + l
    req_parity = (lp + l) % 2

    val = np.zeros_like(s)
    for k in range(k_min, k_max + 1):
        if k % 2 != req_parity:
            continue
        g = gaunt(lp, k, l)
        if g == 0.0:
            continue
        c_A = -Z_A * f_k_vec(s, rho_A, k)
        c_B = -Z_B * ((-1)**k) * f_k_vec(s, rho_B, k)
        val += g * (c_A + c_B)

    norm = np.sqrt((2 * lp + 1) * (2 * l + 1))
    return norm * val


# --------------------------------------------------------------------------- #
# Build full coupling potential W_ij(alpha) on the alpha grid
# --------------------------------------------------------------------------- #

def build_alpha_coupling(
    channels: List[Tuple[int, int]],
    alpha_grid: np.ndarray,
    rho: float,
    Z_A: float, Z_B: float,
    include_ee: bool = True,
) -> np.ndarray:
    """
    Build the full W_ij(alpha_k) coupling array.
    Returns array of shape (n_ch, n_ch, n_alpha).
    Uses charge-center origin for heteronuclear systems.
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
            if l2p == l2:
                W[i, j, :] += nuclear_coupling_vec(l1p, l1, cos_a,
                                                   Z_A, rho_A, Z_B, rho_B)
            if l1p == l1:
                W[i, j, :] += nuclear_coupling_vec(l2p, l2, sin_a,
                                                   Z_A, rho_A, Z_B, rho_B)

            # E-e coupling
            if include_ee:
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


# --------------------------------------------------------------------------- #
# Solve the multichannel alpha-equation
# --------------------------------------------------------------------------- #

def solve_alpha_equation(
    channels: List[Tuple[int, int]],
    rho: float,
    Z_A: float, Z_B: float,
    N_alpha: int = 200,
    include_ee: bool = True,
) -> Dict:
    """
    Solve the coupled-channel alpha eigenvalue problem at fixed rho.

    Discretizes: -1/2 u_i'' + sum_j W_ij(alpha) u_j = mu u_i
    with Dirichlet BCs u(0) = u(pi/2) = 0 on a finite-difference grid.
    """
    n_ch = len(channels)
    alpha_grid = np.linspace(0, np.pi / 2, N_alpha + 2)[1:-1]
    da = alpha_grid[1] - alpha_grid[0]

    # Build coupling potential
    W = build_alpha_coupling(channels, alpha_grid, rho, Z_A, Z_B, include_ee)

    # Assemble full matrix
    N_total = n_ch * N_alpha
    kinetic_diag = 1.0 / da**2
    kinetic_off = -0.5 / da**2

    # Use sparse for large systems
    use_sparse = N_total > 5000

    if use_sparse:
        H = lil_matrix((N_total, N_total))
    else:
        H = np.zeros((N_total, N_total))

    # Kinetic energy blocks (tridiagonal per channel)
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
    if use_sparse:
        H_csr = csr_matrix(H)
        n_eig = min(5, N_total - 2)
        eigenvalues, eigenvectors = eigsh(H_csr, k=n_eig, which='SA')
        idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
    else:
        n_eig = min(5, N_total)
        eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, n_eig - 1])

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
        'mu0': mu0,
        'eigenvalues': eigenvalues,
        'channels': channels,
        'weights': weights,
        'alpha_grid': alpha_grid,
        'n_ch': n_ch,
    }


# --------------------------------------------------------------------------- #
# V_eff for single (0,0) channel (for Plot D)
# --------------------------------------------------------------------------- #

def compute_Veff_00(
    alpha_grid: np.ndarray,
    rho: float, Z_A: float, Z_B: float,
    include_ee: bool = True,
) -> np.ndarray:
    """Compute V_eff(alpha) for the (0,0) channel."""
    channels = [(0, 0)]
    W = build_alpha_coupling(channels, alpha_grid, rho, Z_A, Z_B, include_ee)
    return W[0, 0, :]


# --------------------------------------------------------------------------- #
# Richardson extrapolation
# --------------------------------------------------------------------------- #

def richardson_extrapolate(l_vals: np.ndarray, mu_vals: np.ndarray) -> float:
    """Aitken delta-squared extrapolation using last 3 points."""
    if len(l_vals) < 3:
        return mu_vals[-1]
    m1, m2, m3 = mu_vals[-3], mu_vals[-2], mu_vals[-1]
    denom = (m3 - m2) - (m2 - m1)
    if abs(denom) < 1e-15:
        return m3
    return m3 - (m3 - m2)**2 / denom


# --------------------------------------------------------------------------- #
# Test systems
# --------------------------------------------------------------------------- #

def get_test_systems() -> List[Dict]:
    return [
        {
            'name': 'H2 (1:1)',
            'Z_A': 1.0, 'Z_B': 1.0,
            'R': 1.4, 'R_e': 1.0,
            'homonuclear': True,
            'color': 'blue',
        },
        {
            'name': 'HeH+ (2:1)',
            'Z_A': 2.0, 'Z_B': 1.0,
            'R': 1.46, 'R_e': 1.0,
            'homonuclear': False,
            'color': 'orange',
        },
        {
            'name': 'Hypo (6:1)',
            'Z_A': 6.0, 'Z_B': 1.0,
            'R': 1.4, 'R_e': 1.0,
            'homonuclear': False,
            'color': 'red',
        },
    ]


# --------------------------------------------------------------------------- #
# Main diagnostic
# --------------------------------------------------------------------------- #

def run_diagnostic():
    systems = get_test_systems()

    # Use smaller grid for convergence sweep (accuracy check at end)
    N_alpha_sweep = 200
    N_alpha_detail = 300

    # l_max values: go to 8 for homonuclear, 6 for heteronuclear (cost)
    l_max_homo = list(range(0, 9))
    l_max_hetero = list(range(0, 7))

    all_results = {}
    convergence_data = {}

    print("=" * 80)
    print("ANGULAR WAVEFUNCTION CHANNEL DECOMPOSITION DIAGNOSTIC")
    print("Solving coupled-channel alpha-equation (Paper 15, Eq. 30)")
    print(f"N_alpha = {N_alpha_sweep} (sweep), {N_alpha_detail} (detail)")
    print("Charge-center origin for heteronuclear systems")
    print("=" * 80)

    # ---------------------------------------------------------------------- #
    # Convergence sweep: mu0 vs l_max for each system
    # ---------------------------------------------------------------------- #
    for sys_info in systems:
        name = sys_info['name']
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        l_max_list = l_max_homo if sys_info['homonuclear'] else l_max_hetero

        print(f"\n{'='*60}")
        print(f"System: {name}  (Z_A={sys_info['Z_A']}, Z_B={sys_info['Z_B']}, rho={rho:.4f})")
        print(f"{'='*60}")

        mu_vs_lmax = []

        for l_max in l_max_list:
            channels = build_channel_list(l_max, sys_info['homonuclear'])
            n_ch = len(channels)
            if n_ch == 0:
                continue

            t0 = time.time()
            print(f"  l_max={l_max}: {n_ch:>3} ch ", end='', flush=True)
            result = solve_alpha_equation(
                channels, rho, sys_info['Z_A'], sys_info['Z_B'],
                N_alpha=N_alpha_sweep, include_ee=True,
            )
            dt = time.time() - t0
            mu0 = result['mu0']
            mu_vs_lmax.append((l_max, n_ch, mu0))

            w00 = result['weights'][0]
            print(f"mu0={mu0:>12.6f}  w(0,0)={w00:.4f}  [{dt:.1f}s]")

        convergence_data[name] = mu_vs_lmax
        all_results[name] = {'convergence': mu_vs_lmax}

    # ---------------------------------------------------------------------- #
    # Detailed weight analysis at l_max=4 for all systems
    # ---------------------------------------------------------------------- #
    print("\n\n" + "=" * 80)
    print("DETAILED CHANNEL WEIGHTS (l_max=4)")
    print("=" * 80)

    detail_results = {}
    for sys_info in systems:
        name = sys_info['name']
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        channels = build_channel_list(4, sys_info['homonuclear'])

        result = solve_alpha_equation(
            channels, rho, sys_info['Z_A'], sys_info['Z_B'],
            N_alpha=N_alpha_detail, include_ee=True,
        )
        detail_results[name] = result

        print(f"\n  {name} ({len(channels)} channels):")
        sorted_idx = np.argsort(result['weights'])[::-1]
        for rank, idx in enumerate(sorted_idx):
            w = result['weights'][idx]
            ch = result['channels'][idx]
            if w > 1e-6:
                print(f"    {rank+1:>2}. ({ch[0]},{ch[1]}): {w:.6f}")

        all_results[name]['weights_lmax4'] = {
            str(result['channels'][i]): float(result['weights'][i])
            for i in range(len(result['channels']))
        }

    # ---------------------------------------------------------------------- #
    # Summary table
    # ---------------------------------------------------------------------- #
    print("\n\n" + "=" * 100)
    print("SUMMARY TABLE")
    print("=" * 100)
    header = (f"{'System':<16} {'w(0,0)':<10} {'N(90%)':<8} {'N(95%)':<8} "
              f"{'N(99%)':<8} {'mu(l4)':<12} {'mu(l6)':<12} {'mu(best)':<12} "
              f"{'mu(inf)':<12}")
    print(header)
    print("-" * len(header))

    for sys_info in systems:
        name = sys_info['name']
        r4 = detail_results[name]
        w00 = r4['weights'][0]

        sorted_w = np.sort(r4['weights'])[::-1]
        cum_w = np.cumsum(sorted_w)
        n90 = int(np.searchsorted(cum_w, 0.90) + 1)
        n95 = int(np.searchsorted(cum_w, 0.95) + 1)
        n99 = int(np.searchsorted(cum_w, 0.99) + 1)

        conv = convergence_data[name]
        mu_dict = {lm: mu for lm, nc, mu in conv}
        mu4 = mu_dict.get(4, np.nan)
        mu6 = mu_dict.get(6, np.nan)
        mu_best = conv[-1][2] if conv else np.nan
        best_lmax = conv[-1][0] if conv else 0

        l_arr = np.array([lm for lm, _, _ in conv if lm >= 2])
        mu_arr = np.array([mu for lm, _, mu in conv if lm >= 2])
        mu_inf = richardson_extrapolate(l_arr, mu_arr) if len(l_arr) >= 3 else np.nan

        print(f"{name:<16} {w00:<10.4f} {n90:<8} {n95:<8} {n99:<8} "
              f"{mu4:<12.6f} {mu6:<12.6f} {mu_best:<12.6f} {mu_inf:<12.6f}")

        all_results[name]['summary'] = {
            'w00': float(w00),
            'n_90': int(n90), 'n_95': int(n95), 'n_99': int(n99),
            'mu_lmax4': float(mu4), 'mu_lmax6': float(mu6),
            'mu_best': float(mu_best), 'best_lmax': int(best_lmax),
            'mu_inf': float(mu_inf),
        }

    # ---------------------------------------------------------------------- #
    # Plot A: Channel weight distribution at l_max=4
    # ---------------------------------------------------------------------- #
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for idx, sys_info in enumerate(systems):
        ax = axes[idx]
        name = sys_info['name']
        result = detail_results[name]
        labels = [f'({l1},{l2})' for l1, l2 in result['channels']]
        colors_bar = ['steelblue' if w > 0.01 else 'lightcoral'
                      for w in result['weights']]
        ax.bar(range(len(result['channels'])), result['weights'],
               color=colors_bar, alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(result['channels'])))
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_ylabel('Channel weight')
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        ax.set_title(f'{name}\nrho={rho:.3f}, {len(result["channels"])} ch')
        ax.set_ylim(0, 1.05)
        ax.axhline(0.01, color='red', ls=':', alpha=0.5, label='1% threshold')
        ax.legend(fontsize=7)
        ax.grid(axis='y', alpha=0.3)

    plt.suptitle('Plot A: Channel Weight Distribution (l_max=4, with e-e coupling)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'angular_channel_weights.png', dpi=150)
    plt.close()
    print(f"\nSaved: {PLOT_DIR / 'angular_channel_weights.png'}")

    # ---------------------------------------------------------------------- #
    # Plot B: Eigenvalue convergence mu0 vs l_max
    # ---------------------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(10, 7))
    for sys_info in systems:
        name = sys_info['name']
        conv = convergence_data[name]
        l_vals = [lm for lm, nc, mu in conv]
        mu_vals = [mu for lm, nc, mu in conv]
        ax.plot(l_vals, mu_vals, 'o-', color=sys_info['color'], label=name,
                markersize=7, linewidth=2)

    ax.set_xlabel('l_max', fontsize=12)
    ax.set_ylabel('Ground-state eigenvalue mu_0', fontsize=12)
    ax.set_title('Plot B: Angular Eigenvalue Convergence vs Channel Truncation',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'angular_eigenvalue_convergence.png', dpi=150)
    plt.close()
    print(f"Saved: {PLOT_DIR / 'angular_eigenvalue_convergence.png'}")

    # ---------------------------------------------------------------------- #
    # Plot C: Cumulative channel weight
    # ---------------------------------------------------------------------- #
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for idx, sys_info in enumerate(systems):
        ax = axes[idx]
        name = sys_info['name']
        result = detail_results[name]
        sorted_w = np.sort(result['weights'])[::-1]
        cum_w = np.cumsum(sorted_w)

        ax.plot(range(1, len(cum_w) + 1), cum_w, 'o-', color=sys_info['color'],
                markersize=7, linewidth=2)
        ax.axhline(0.90, color='gray', ls='--', alpha=0.5, label='90%')
        ax.axhline(0.95, color='gray', ls=':', alpha=0.5, label='95%')
        ax.axhline(0.99, color='gray', ls='-.', alpha=0.5, label='99%')
        ax.set_xlabel('Channels (sorted by weight)', fontsize=10)
        ax.set_ylabel('Cumulative weight', fontsize=10)
        ax.set_title(f'{name}', fontsize=12)
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Plot C: Cumulative Channel Weight (l_max=4)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'angular_cumulative_weight.png', dpi=150)
    plt.close()
    print(f"Saved: {PLOT_DIR / 'angular_cumulative_weight.png'}")

    # ---------------------------------------------------------------------- #
    # Plot D: V_eff(alpha) for (0,0) channel
    # ---------------------------------------------------------------------- #
    alpha_plot = np.linspace(0.05, np.pi / 2 - 0.05, 300)
    fig, ax = plt.subplots(figsize=(10, 7))
    for sys_info in systems:
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        Veff = compute_Veff_00(alpha_plot, rho, sys_info['Z_A'], sys_info['Z_B'],
                               include_ee=True)
        ax.plot(np.degrees(alpha_plot), Veff, '-', color=sys_info['color'],
                label=sys_info['name'], linewidth=2)

    ax.set_xlabel('alpha (degrees)', fontsize=12)
    ax.set_ylabel('V_eff(alpha) for (0,0) channel', fontsize=12)
    ax.set_title('Plot D: Effective Potential in Alpha for (0,0) Channel', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'angular_Veff_00.png', dpi=150)
    plt.close()
    print(f"Saved: {PLOT_DIR / 'angular_Veff_00.png'}")

    # ---------------------------------------------------------------------- #
    # Save data
    # ---------------------------------------------------------------------- #
    json_results = {}
    for name, data in all_results.items():
        json_results[name] = {
            'convergence': [
                {'l_max': int(lm), 'n_ch': int(nc), 'mu0': float(mu)}
                for lm, nc, mu in data['convergence']
            ],
            'weights_lmax4': data.get('weights_lmax4', {}),
        }
        if 'summary' in data:
            json_results[name]['summary'] = data['summary']

    data_path = DATA_DIR / 'angular_wavefunction_decomposition.json'
    with open(data_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"\nData saved: {data_path}")

    # ---------------------------------------------------------------------- #
    # Key findings
    # ---------------------------------------------------------------------- #
    print("\n" + "=" * 80)
    print("KEY FINDINGS")
    print("=" * 80)

    for sys_info in systems:
        name = sys_info['name']
        s = all_results[name].get('summary', {})
        w00 = s.get('w00', 0)
        n90 = s.get('n_90', 0)
        n99 = s.get('n_99', 0)
        mu_best = s.get('mu_best', 0)
        mu_inf = s.get('mu_inf', 0)

        delta = abs(mu_best - mu_inf) / abs(mu_inf) * 100 if abs(mu_inf) > 1e-10 else 0
        print(f"\n  {name}:")
        print(f"    (0,0) channel carries {w00*100:.1f}% of wavefunction")
        print(f"    {n90} channels for 90%, {n99} channels for 99% of norm")
        print(f"    mu0 best: {mu_best:.6f}, estimated mu(inf): {mu_inf:.6f}")
        print(f"    Relative gap to extrapolated: {delta:.2f}%")

    # Comparison
    h2_s = all_results['H2 (1:1)'].get('summary', {})
    heh_s = all_results['HeH+ (2:1)'].get('summary', {})
    hyp_s = all_results['Hypo (6:1)'].get('summary', {})

    print(f"\n  Wavefunction localization comparison:")
    print(f"    {'System':<16} {'w(0,0)':<10} {'N(90%)':<8} {'N(99%)':<8}")
    for label, s in [('H2 (1:1)', h2_s), ('HeH+ (2:1)', heh_s), ('Hypo (6:1)', hyp_s)]:
        print(f"    {label:<16} {s.get('w00',0):<10.4f} "
              f"{s.get('n_90',0):<8} {s.get('n_99',0):<8}")

    if hyp_s.get('w00', 1) < h2_s.get('w00', 0):
        print("\n  => Asymmetric systems spread wavefunction across MORE channels")
        print("  => CONFIRMS: angular wavefunction convergence is the bottleneck")
        print("  => The potential converges fast (Task 2), but the wavefunction")
        print("     localizes near the heavy nucleus, requiring many partial waves")
    else:
        print("\n  => Asymmetric wavefunction is NOT more spread — revisit hypothesis")


if __name__ == '__main__':
    run_diagnostic()
