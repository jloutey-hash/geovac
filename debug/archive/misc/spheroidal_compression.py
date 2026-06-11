#!/usr/bin/env python3
"""
Prolate Spheroidal Basis Compression Test (Task 4).

Hypothesis: Prolate spheroidal harmonics — eigenfunctions of the one-electron
two-center angular equation — should concentrate the two-electron wavefunction
into fewer channels than the Legendre basis, because they are the natural
angular basis for a bond.

Method:
  Part 1: Compute prolate spheroidal harmonics S_n(eta; c, b) as Legendre
          expansions by diagonalizing the tridiagonal/pentadiagonal matrix.
  Part 2: Build the transformation matrix T (Legendre -> spheroidal).
  Part 3: Transform the Task 3 wavefunctions into spheroidal basis.
  Part 4: Measure compression (channel counts for 90/95/99% norm).
"""

import numpy as np
from scipy.linalg import eigh
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import sys
import time
from typing import List, Tuple, Dict

# Reuse the alpha-equation solver and coupling machinery from Task 3
sys.path.insert(0, str(Path(__file__).parent))
from eigenchannel_diagnostic import gaunt, build_channel_list
from angular_wavefunction_diagnostic import (
    solve_alpha_equation, build_alpha_coupling,
)

PLOT_DIR = Path(__file__).parent / "plots"
DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)


# =========================================================================== #
# Part 1: Prolate spheroidal harmonics via Legendre expansion
# =========================================================================== #

def eta2_matrix_element(l: int, lp: int) -> float:
    """
    Compute <P_l | eta^2 | P_{l'}> over [-1, 1].

    Uses the recurrence eta P_l = [(l+1)P_{l+1} + l P_{l-1}] / (2l+1) twice.
    Non-zero for |l - l'| = 0 or 2 (pentadiagonal in l).
    """
    # <P_l | eta^2 | P_{l'}> = sum over intermediates from applying eta twice
    # eta P_{l'} = a_{l'} P_{l'+1} + b_{l'} P_{l'-1}
    # where a_{l'} = (l'+1)/(2l'+1), b_{l'} = l'/(2l'+1)
    #
    # eta^2 P_{l'} = a_{l'} [a_{l'+1} P_{l'+2} + b_{l'+1} P_{l'}]
    #              + b_{l'} [a_{l'-1} P_{l'} + b_{l'-1} P_{l'-2}]
    #
    # <P_l | eta^2 | P_{l'}> picks out the coefficient of P_l in this expansion,
    # times the norm 2/(2l+1).

    if abs(l - lp) > 2:
        return 0.0

    norm = 2.0 / (2 * l + 1)  # <P_l | P_l>

    def a(n: int) -> float:
        return (n + 1) / (2 * n + 1) if n >= 0 else 0.0

    def b(n: int) -> float:
        return n / (2 * n + 1) if n >= 0 else 0.0

    val = 0.0

    # Coefficient of P_l in eta^2 P_{l'}:
    # From a_{l'} * eta P_{l'+1}: contributes a_{l'} * [a_{l'+1} delta_{l,l'+2} + b_{l'+1} delta_{l,l'}]
    # From b_{l'} * eta P_{l'-1}: contributes b_{l'} * [a_{l'-1} delta_{l,l'} + b_{l'-1} delta_{l,l'-2}]

    if l == lp + 2:
        val = a(lp) * a(lp + 1)
    elif l == lp:
        val = a(lp) * b(lp + 1) + b(lp) * a(lp - 1)
    elif l == lp - 2:
        val = b(lp) * b(lp - 1)

    return val * norm


def eta1_matrix_element(l: int, lp: int) -> float:
    """
    Compute <P_l | eta | P_{l'}> over [-1, 1].

    Non-zero for |l - l'| = 1 (tridiagonal).
    """
    norm = 2.0 / (2 * l + 1)  # <P_l | P_l>

    if lp == l + 1:
        # eta P_{l+1} has coefficient b(l+1) for P_l
        return (l + 1) / (2 * (l + 1) + 1) * norm
    elif lp == l - 1:
        # eta P_{l-1} has coefficient a(l-1) for P_l
        return l / (2 * l - 1) * norm if l >= 1 else 0.0
    else:
        return 0.0


def build_spheroidal_matrix(L_max: int, c: float, b: float = 0.0) -> np.ndarray:
    """
    Build the matrix M for the prolate spheroidal angular equation:
      d/deta [(1-eta^2) dS/deta] + [A - c^2 eta^2 - b eta - m^2/(1-eta^2)] S = 0

    For m=0, expanding S(eta) = sum_l d_l P_l(eta), the secular equation is:
      M d = A d
    where M_{ll'} = l(l+1) delta_{ll'} + c^2 <P_l|eta^2|P_{l'}> + b <P_l|eta|P_{l'}>

    Note: the angular equation has +A as eigenvalue and -c^2 eta^2 as potential,
    so M = l(l+1) + c^2 * <eta^2> + b * <eta>.

    Parameters
    ----------
    L_max : int
        Maximum angular momentum in expansion.
    c : float
        Spheroidal parameter (c = rho for our context).
    b : float
        Charge asymmetry parameter (b = (Z_A - Z_B) * R / 2 or similar).

    Returns
    -------
    M : ndarray, shape (L_max+1, L_max+1)
        Matrix whose eigenvectors give the Legendre expansion coefficients.
    """
    N = L_max + 1
    M = np.zeros((N, N))

    for l in range(N):
        for lp in range(N):
            if l == lp:
                M[l, lp] += l * (l + 1)
            M[l, lp] += c**2 * eta2_matrix_element(l, lp)
            if b != 0.0:
                M[l, lp] += b * eta1_matrix_element(l, lp)

    return M


def compute_spheroidal_harmonics(
    L_max: int, c: float, b: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute prolate spheroidal harmonics as Legendre expansions.

    Returns
    -------
    eigenvalues : ndarray, shape (L_max+1,)
        Separation constants A_n(c, b).
    T : ndarray, shape (L_max+1, L_max+1)
        Transformation matrix: T[n, l] = d_l^{(n)}(c, b).
        Columns of the eigenvector matrix, transposed so rows = spheroidal index.
    """
    M = build_spheroidal_matrix(L_max, c, b)
    eigenvalues, eigvecs = eigh(M)
    # eigvecs[:, n] = coefficients d_l for the n-th spheroidal harmonic
    # T[n, l] = eigvecs[l, n] => T = eigvecs.T
    T = eigvecs.T
    return eigenvalues, T


# =========================================================================== #
# Part 2: Solve alpha-equation and get full wavefunctions
# =========================================================================== #

def solve_alpha_with_wavefunctions(
    channels: List[Tuple[int, int]],
    rho: float, Z_A: float, Z_B: float,
    N_alpha: int = 300,
) -> Dict:
    """
    Solve the coupled-channel alpha-equation and return the full wavefunction
    u_i(alpha) for each channel, not just the integrated weights.
    """
    from scipy.sparse.linalg import eigsh
    from scipy.sparse import lil_matrix, csr_matrix

    n_ch = len(channels)
    alpha_grid = np.linspace(0, np.pi / 2, N_alpha + 2)[1:-1]
    da = alpha_grid[1] - alpha_grid[0]

    # Build coupling potential
    W = build_alpha_coupling(channels, alpha_grid, rho, Z_A, Z_B, include_ee=True)

    # Assemble full matrix
    N_total = n_ch * N_alpha
    kinetic_diag = 1.0 / da**2
    kinetic_off = -0.5 / da**2

    use_sparse = N_total > 5000
    if use_sparse:
        H = lil_matrix((N_total, N_total))
    else:
        H = np.zeros((N_total, N_total))

    for i_ch in range(n_ch):
        off = i_ch * N_alpha
        for k in range(N_alpha):
            H[off + k, off + k] += kinetic_diag
            if k > 0:
                H[off + k, off + k - 1] += kinetic_off
            if k < N_alpha - 1:
                H[off + k, off + k + 1] += kinetic_off

    for i_ch in range(n_ch):
        for j_ch in range(n_ch):
            w_ij = W[i_ch, j_ch, :]
            if np.max(np.abs(w_ij)) < 1e-15:
                continue
            i_off = i_ch * N_alpha
            j_off = j_ch * N_alpha
            for k in range(N_alpha):
                H[i_off + k, j_off + k] += w_ij[k]

    if use_sparse:
        H_csr = csr_matrix(H)
        eigenvalues, eigenvectors = eigsh(H_csr, k=min(5, N_total - 2), which='SA')
        idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
    else:
        n_eig = min(5, N_total)
        eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, n_eig - 1])

    mu0 = eigenvalues[0]
    psi0 = eigenvectors[:, 0]

    # Extract per-channel wavefunctions u_i(alpha)
    u_channels = np.zeros((n_ch, N_alpha))
    for i_ch in range(n_ch):
        offset = i_ch * N_alpha
        u_channels[i_ch, :] = psi0[offset:offset + N_alpha]

    # Compute weights
    weights = np.sum(u_channels**2, axis=1) * da
    total_weight = np.sum(weights)
    if total_weight > 1e-15:
        weights /= total_weight

    return {
        'mu0': mu0,
        'channels': channels,
        'weights': weights,
        'alpha_grid': alpha_grid,
        'da': da,
        'u_channels': u_channels,  # shape (n_ch, N_alpha)
        'total_norm': total_weight,
    }


# =========================================================================== #
# Part 3: Transform wavefunction to spheroidal basis
# =========================================================================== #

def transform_to_spheroidal(
    result: Dict,
    T: np.ndarray,
    L_max: int,
    homonuclear: bool,
) -> Dict:
    """
    Transform the two-electron channel wavefunction from Legendre (l1,l2)
    to spheroidal (n1,n2) basis.

    For each alpha point:
      u_{n1,n2}(alpha) = sum_{l1,l2} T[n1,l1] T[n2,l2] u_{l1,l2}(alpha)

    Then compute spheroidal channel weights.
    """
    channels = result['channels']
    u_channels = result['u_channels']  # (n_ch_leg, N_alpha)
    da = result['da']
    n_ch_leg = len(channels)
    N_alpha = u_channels.shape[1]

    # Build the spheroidal channel list (same structure as Legendre channels)
    sph_channels = []
    for n1 in range(L_max + 1):
        for n2 in range(L_max + 1):
            if homonuclear and (n1 + n2) % 2 != 0:
                continue  # preserve gerade constraint
            sph_channels.append((n1, n2))
    sph_channels.sort(key=lambda c: (c[0] + c[1], c[0]))

    n_ch_sph = len(sph_channels)

    # Build the two-electron transformation matrix
    # For each (n1, n2) spheroidal channel and (l1, l2) Legendre channel:
    # T2e[(n1,n2), (l1,l2)] = T[n1,l1] * T[n2,l2]
    T2e = np.zeros((n_ch_sph, n_ch_leg))
    for i_sph, (n1, n2) in enumerate(sph_channels):
        for j_leg, (l1, l2) in enumerate(channels):
            if n1 <= L_max and n2 <= L_max and l1 <= L_max and l2 <= L_max:
                T2e[i_sph, j_leg] = T[n1, l1] * T[n2, l2]

    # Transform: u_sph(alpha) = T2e @ u_leg(alpha) at each alpha point
    u_sph = T2e @ u_channels  # (n_ch_sph, N_alpha)

    # Compute spheroidal channel weights
    weights_sph = np.sum(u_sph**2, axis=1) * da
    total = np.sum(weights_sph)
    if total > 1e-15:
        weights_sph /= total

    return {
        'channels': sph_channels,
        'weights': weights_sph,
        'u_channels': u_sph,
        'T2e': T2e,
        'total_norm': total,
    }


# =========================================================================== #
# Test systems (consistent with Task 3)
# =========================================================================== #

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


# =========================================================================== #
# Main diagnostic
# =========================================================================== #

def run_diagnostic():
    systems = get_test_systems()
    L_max = 20  # Spheroidal harmonic expansion truncation
    l_max_alpha = 4  # Channel truncation for alpha-equation (matches Task 3)
    N_alpha = 300

    print("=" * 80)
    print("PROLATE SPHEROIDAL BASIS COMPRESSION TEST (Task 4)")
    print("=" * 80)
    print(f"Spheroidal L_max = {L_max}")
    print(f"Alpha-equation l_max = {l_max_alpha}, N_alpha = {N_alpha}")

    all_results = {}

    # ================================================================== #
    # Part 1: Compute spheroidal harmonics for each system
    # ================================================================== #
    print("\n" + "=" * 80)
    print("PART 1: PROLATE SPHEROIDAL HARMONICS")
    print("=" * 80)

    spheroidal_data = {}
    for sys_info in systems:
        name = sys_info['name']
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        Z_A, Z_B = sys_info['Z_A'], sys_info['Z_B']

        # Spheroidal parameter c = rho (natural scale)
        c = rho

        # Charge asymmetry: b proportional to (Z_A - Z_B)
        # In prolate spheroidal coordinates, b = (Z_B - Z_A) * R / 2
        # For the angular equation on eta in [-1,1], the asymmetry term is b*eta
        # Using charge-center scaled: b = (Z_A - Z_B) * rho (simplified)
        b = (Z_A - Z_B) * rho

        print(f"\n  {name}: c = {c:.4f}, b = {b:.4f}")

        eigenvalues, T = compute_spheroidal_harmonics(L_max, c, b)

        # Verify orthogonality: T^T T = I
        ortho_err = np.max(np.abs(T @ T.T - np.eye(L_max + 1)))
        print(f"    Orthogonality check: max|T T^T - I| = {ortho_err:.2e}")

        # Show first few eigenvalues
        print(f"    First 6 separation constants A_n: "
              f"{', '.join(f'{a:.4f}' for a in eigenvalues[:6])}")

        # Show how concentrated the first few harmonics are in Legendre basis
        for n in range(min(4, L_max + 1)):
            coeffs = T[n, :]
            dom_l = np.argmax(np.abs(coeffs))
            max_coeff = coeffs[dom_l]
            n_sig = np.sum(coeffs**2 > 0.01)  # Legendre components > 1%
            print(f"    S_{n}: dominant l={dom_l} (coeff={max_coeff:.4f}), "
                  f"{n_sig} significant Legendre components")

        spheroidal_data[name] = {
            'c': c, 'b': b,
            'eigenvalues': eigenvalues,
            'T': T,
        }

    # ================================================================== #
    # Part 2-3: Solve alpha-equation, transform, measure compression
    # ================================================================== #
    print("\n" + "=" * 80)
    print("PARTS 2-3: WAVEFUNCTION TRANSFORMATION & COMPRESSION")
    print("=" * 80)

    leg_results = {}
    sph_results = {}

    for sys_info in systems:
        name = sys_info['name']
        rho = sys_info['R'] / (2 * sys_info['R_e'])
        channels = build_channel_list(l_max_alpha, sys_info['homonuclear'])

        print(f"\n  {name}: {len(channels)} Legendre channels at l_max={l_max_alpha}")

        t0 = time.time()
        result = solve_alpha_with_wavefunctions(
            channels, rho, sys_info['Z_A'], sys_info['Z_B'], N_alpha=N_alpha,
        )
        dt = time.time() - t0
        print(f"    Alpha solver: mu0 = {result['mu0']:.6f}  [{dt:.1f}s]")

        leg_results[name] = result

        # Transform to spheroidal basis
        # Use only l_max_alpha for T (must match the channel space)
        T_full = spheroidal_data[name]['T']
        # Restrict T to l_max_alpha x l_max_alpha block
        T_sub = T_full[:l_max_alpha + 1, :l_max_alpha + 1]
        # Re-orthogonalize the sub-block (it's already orthogonal if L_max >= l_max_alpha)
        # Actually T_sub is a submatrix of an orthogonal matrix — need to re-check.
        # T_full is (L_max+1) x (L_max+1) orthogonal. T_sub is (l_max+1) x (l_max+1)
        # which is NOT necessarily orthogonal. We need to re-diagonalize at l_max_alpha.

        c = spheroidal_data[name]['c']
        b = spheroidal_data[name]['b']
        eigenvalues_sub, T_sub = compute_spheroidal_harmonics(l_max_alpha, c, b)

        ortho_err = np.max(np.abs(T_sub @ T_sub.T - np.eye(l_max_alpha + 1)))
        print(f"    T_sub orthogonality: max|T T^T - I| = {ortho_err:.2e}")

        sph_result = transform_to_spheroidal(
            result, T_sub, l_max_alpha, sys_info['homonuclear'],
        )
        sph_results[name] = sph_result

        # Print comparison
        print(f"\n    {'Basis':<12} {'w(0,0)':<10} {'N(90%)':<8} {'N(95%)':<8} {'N(99%)':<8}")
        print(f"    {'-'*46}")

        for label, res in [('Legendre', result), ('Spheroidal', sph_result)]:
            w = res['weights']
            w00 = w[0]
            sorted_w = np.sort(w)[::-1]
            cum_w = np.cumsum(sorted_w)
            n90 = int(np.searchsorted(cum_w, 0.90) + 1)
            n95 = int(np.searchsorted(cum_w, 0.95) + 1)
            n99 = int(np.searchsorted(cum_w, 0.99) + 1)
            print(f"    {label:<12} {w00:<10.4f} {n90:<8} {n95:<8} {n99:<8}")

    # ================================================================== #
    # Part 4: Summary table
    # ================================================================== #
    print("\n\n" + "=" * 100)
    print("SUMMARY TABLE: LEGENDRE vs SPHEROIDAL BASIS COMPRESSION")
    print("=" * 100)

    header = (f"{'System':<16} | {'--- Legendre ---':<36} | {'--- Spheroidal ---':<36} | {'Ratio':<8}")
    print(header)
    sub = (f"{'':16} | {'w(0,0)':<9} {'N(90%)':<8} {'N(95%)':<8} {'N(99%)':<8} "
           f"| {'w(0,0)':<9} {'N(90%)':<8} {'N(95%)':<8} {'N(99%)':<8} | {'N95L/S':<8}")
    print(sub)
    print("-" * 105)

    summary_data = {}
    for sys_info in systems:
        name = sys_info['name']
        leg = leg_results[name]
        sph = sph_results[name]

        def get_stats(w):
            sorted_w = np.sort(w)[::-1]
            cum = np.cumsum(sorted_w)
            return {
                'w00': float(w[0]),
                'n90': int(np.searchsorted(cum, 0.90) + 1),
                'n95': int(np.searchsorted(cum, 0.95) + 1),
                'n99': int(np.searchsorted(cum, 0.99) + 1),
            }

        ls = get_stats(leg['weights'])
        ss = get_stats(sph['weights'])
        ratio = ls['n95'] / ss['n95'] if ss['n95'] > 0 else float('inf')

        print(f"{name:<16} | {ls['w00']:<9.4f} {ls['n90']:<8} {ls['n95']:<8} {ls['n99']:<8} "
              f"| {ss['w00']:<9.4f} {ss['n90']:<8} {ss['n95']:<8} {ss['n99']:<8} | {ratio:<8.2f}")

        summary_data[name] = {
            'legendre': ls,
            'spheroidal': ss,
            'compression_ratio_95': ratio,
        }

    # ================================================================== #
    # Plot A: Side-by-side bar charts (Legendre vs Spheroidal)
    # ================================================================== #
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))

    for idx, sys_info in enumerate(systems):
        name = sys_info['name']

        # Legendre weights
        ax = axes[idx, 0]
        leg = leg_results[name]
        labels_leg = [f'({l1},{l2})' for l1, l2 in leg['channels']]
        colors_leg = ['steelblue' if w > 0.01 else 'lightgray' for w in leg['weights']]
        ax.bar(range(len(leg['channels'])), leg['weights'],
               color=colors_leg, alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(leg['channels'])))
        ax.set_xticklabels(labels_leg, rotation=90, fontsize=6)
        ax.set_ylabel('Channel weight')
        ax.set_title(f'{name} — Legendre basis', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.axhline(0.01, color='red', ls=':', alpha=0.5)
        ax.grid(axis='y', alpha=0.3)

        # Spheroidal weights
        ax = axes[idx, 1]
        sph = sph_results[name]
        labels_sph = [f'({n1},{n2})' for n1, n2 in sph['channels']]
        colors_sph = ['darkorange' if w > 0.01 else 'lightgray' for w in sph['weights']]
        ax.bar(range(len(sph['channels'])), sph['weights'],
               color=colors_sph, alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(sph['channels'])))
        ax.set_xticklabels(labels_sph, rotation=90, fontsize=6)
        ax.set_ylabel('Channel weight')
        ax.set_title(f'{name} — Spheroidal basis', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.axhline(0.01, color='red', ls=':', alpha=0.5)
        ax.grid(axis='y', alpha=0.3)

    plt.suptitle('Plot A: Channel Weights — Legendre (left) vs Spheroidal (right)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'spheroidal_compression_weights.png', dpi=150)
    plt.close()
    print(f"\nSaved: {PLOT_DIR / 'spheroidal_compression_weights.png'}")

    # ================================================================== #
    # Plot B: Cumulative weight curves — Legendre vs Spheroidal
    # ================================================================== #
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, sys_info in enumerate(systems):
        ax = axes[idx]
        name = sys_info['name']

        for label, res, color, ls in [
            ('Legendre', leg_results[name], sys_info['color'], '-'),
            ('Spheroidal', sph_results[name], 'darkgreen', '--'),
        ]:
            sorted_w = np.sort(res['weights'])[::-1]
            cum_w = np.cumsum(sorted_w)
            ax.plot(range(1, len(cum_w) + 1), cum_w, f'o{ls}',
                    color=color, markersize=5, linewidth=2, label=label)

        ax.axhline(0.90, color='gray', ls='--', alpha=0.4, label='90%')
        ax.axhline(0.95, color='gray', ls=':', alpha=0.4, label='95%')
        ax.axhline(0.99, color='gray', ls='-.', alpha=0.4, label='99%')
        ax.set_xlabel('Channels (sorted by weight)')
        ax.set_ylabel('Cumulative weight')
        ax.set_title(name)
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Plot B: Cumulative Channel Weight — Legendre vs Spheroidal',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'spheroidal_compression_cumulative.png', dpi=150)
    plt.close()
    print(f"Saved: {PLOT_DIR / 'spheroidal_compression_cumulative.png'}")

    # ================================================================== #
    # Plot C: Spheroidal harmonics S_n(eta) for 6:1 system
    # ================================================================== #
    eta_grid = np.linspace(-1, 1, 500)
    sys_61 = [s for s in systems if '6:1' in s['name']][0]
    sd = spheroidal_data[sys_61['name']]
    T_full = sd['T']

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # Top row: first 5 spheroidal harmonics S_n(eta)
    for n in range(5):
        ax = axes[0, n] if n < 3 else axes[1, n - 3]
        coeffs = T_full[n, :]  # d_l^{(n)} coefficients

        # Evaluate S_n(eta) = sum_l d_l P_l(eta)
        S_n = np.zeros_like(eta_grid)
        for l in range(L_max + 1):
            if abs(coeffs[l]) < 1e-15:
                continue
            # Evaluate P_l(eta) using numpy
            Pl = np.polynomial.legendre.Legendre.basis(l)(eta_grid)
            S_n += coeffs[l] * Pl

        ax.plot(eta_grid, S_n, 'r-', linewidth=2, label=f'$S_{n}(\\eta)$')

        # Overlay corresponding Legendre polynomial P_n(eta)
        Pn = np.polynomial.legendre.Legendre.basis(n)(eta_grid)
        ax.plot(eta_grid, Pn, 'b--', linewidth=1.5, alpha=0.6, label=f'$P_{n}(\\eta)$')

        ax.set_xlabel('$\\eta$')
        ax.set_title(f'$S_{n}$ (A={sd["eigenvalues"][n]:.2f})')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-1.5, 1.5)

    # Last panel: Legendre coefficient magnitudes for first 5 harmonics
    ax = axes[1, 2]
    for n in range(5):
        ax.semilogy(range(min(12, L_max + 1)), np.abs(T_full[n, :12]) + 1e-16,
                     'o-', markersize=4, label=f'$S_{n}$')
    ax.set_xlabel('Legendre index $l$')
    ax.set_ylabel('$|d_l^{(n)}|$')
    ax.set_title(f'Legendre coefficients (6:1, c={sd["c"]:.3f}, b={sd["b"]:.3f})')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1e-8, 2)

    plt.suptitle(f'Plot C: Spheroidal Harmonics for {sys_61["name"]} '
                 f'(c={sd["c"]:.3f}, b={sd["b"]:.3f})',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'spheroidal_harmonics_6to1.png', dpi=150)
    plt.close()
    print(f"Saved: {PLOT_DIR / 'spheroidal_harmonics_6to1.png'}")

    # ================================================================== #
    # Save JSON data
    # ================================================================== #
    json_out = {}
    for sys_info in systems:
        name = sys_info['name']
        sd = spheroidal_data[name]
        json_out[name] = {
            'spheroidal_params': {
                'c': sd['c'],
                'b': sd['b'],
                'L_max': L_max,
            },
            'separation_constants': sd['eigenvalues'][:10].tolist(),
            'summary': summary_data[name],
            'legendre_weights': {
                str(ch): float(w) for ch, w in
                zip(leg_results[name]['channels'], leg_results[name]['weights'])
            },
            'spheroidal_weights': {
                str(ch): float(w) for ch, w in
                zip(sph_results[name]['channels'], sph_results[name]['weights'])
            },
        }

    data_path = DATA_DIR / 'spheroidal_compression.json'
    with open(data_path, 'w') as f:
        json.dump(json_out, f, indent=2)
    print(f"\nData saved: {data_path}")

    # ================================================================== #
    # Bottom line
    # ================================================================== #
    print("\n" + "=" * 80)
    print("BOTTOM LINE")
    print("=" * 80)

    hyp = summary_data.get('Hypo (6:1)', {})
    h2 = summary_data.get('H2 (1:1)', {})

    hyp_leg = hyp.get('legendre', {})
    hyp_sph = hyp.get('spheroidal', {})
    ratio_61 = hyp.get('compression_ratio_95', 1.0)

    h2_leg = h2.get('legendre', {})
    h2_sph = h2.get('spheroidal', {})
    ratio_h2 = h2.get('compression_ratio_95', 1.0)

    print(f"\n  H2 (1:1):  Legendre N(95%)={h2_leg.get('n95','?')}, "
          f"Spheroidal N(95%)={h2_sph.get('n95','?')}, "
          f"compression ratio={ratio_h2:.2f}x")
    print(f"  Hypo (6:1): Legendre N(95%)={hyp_leg.get('n95','?')}, "
          f"Spheroidal N(95%)={hyp_sph.get('n95','?')}, "
          f"compression ratio={ratio_61:.2f}x")

    if ratio_61 > 1.5:
        print("\n  => CONFIRMED: Spheroidal basis provides significant compression")
        print("     for asymmetric systems. The spheroidal harmonics concentrate the")
        print("     wavefunction because they encode the two-center potential geometry.")
        print("     This suggests using spheroidal channel functions instead of Legendre")
        print("     polynomials in the Level 4 multichannel expansion for heteronuclear")
        print("     molecules would substantially improve convergence.")
    elif ratio_61 > 1.1:
        print("\n  => MODEST compression from spheroidal basis for asymmetric systems.")
        print("     The improvement exists but may not justify the implementation cost.")
        print("     The wavefunction spreading is partly due to electron-electron")
        print("     correlation (not just nuclear potential asymmetry).")
    else:
        print("\n  => NO significant compression from spheroidal basis.")
        print("     The angular wavefunction spreading is NOT primarily due to using")
        print("     the wrong one-electron basis — it is fundamentally a two-electron")
        print("     correlation effect. The Legendre basis is already near-optimal for")
        print("     the coupled two-electron problem, and a better one-electron angular")
        print("     basis cannot help much. This suggests the path forward is:")
        print("     (a) better treatment of e-e cusp/correlation, or")
        print("     (b) non-uniform grids in alpha, or")
        print("     (c) accepting that asymmetric bonds need more channels.")


if __name__ == '__main__':
    run_diagnostic()
