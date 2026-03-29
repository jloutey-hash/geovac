"""
Projected PK Pseudopotential Diagnostic for LiH l_max Convergence.

Tests three PK modes:
  Mode 1: Channel-blind -- same V_PK for all channels (current default).
  Mode 2: l-dependent (binary) -- V_PK weighted by delta_{l_i,0} per electron.
  Mode 3: Projected -- V_PK weighted by fractional l=0 content of each
           eigenchannel from nuclear coupling diagonalization.

Diagnostic: solve angular eigenvalue problem at fixed R, R_e for l_max = 2..6,
compare mu_0 convergence across modes.

References:
  Paper 17 Sec III.B (PK pseudopotential), Sec V.C (l_max divergence)
  Paper 15 Sec V (multichannel coupling)
"""

import sys
import json
import numpy as np
from scipy.linalg import eigh
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.level4_multichannel import (
    _channel_list,
    build_angular_hamiltonian,
    compute_nuclear_coupling,
    compute_pk_pseudopotential,
)


# === Physical parameters ===
# LiH at experimental R_eq
R_LiH = 3.015        # bohr (experimental)
Z_A = 2.69            # Clementi-Raimondi Z_eff for Li (after 1s^2 partial screening)
Z_B = 1.0             # H
R_e_test = 1.5        # bohr (reasonable valence hyperradius)

# PK parameters from Paper 17 Table I (Li core)
PK_A = 6.93           # Ha*bohr^2
PK_B_param = 7.00     # bohr^-2 (Gaussian exponent beta)

# Charge-center origin for asymmetric Z_A >> Z_B
# z0 = R * (Z_A - Z_B) / (2 * (Z_A + Z_B))
z0 = R_LiH * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))

# Nuclear positions in R_e units (rho = R/(2*R_e))
rho = R_LiH / (2.0 * R_e_test)

# For heteronuclear: A is at +rho_A, B at -rho_B from origin
# With z0 = 0 and equal Z_eff, rho_A = rho_B = rho
rho_A = rho / 2.0  # will be adjusted by build_angular_hamiltonian
rho_B = rho / 2.0

# l_max values to test
L_MAX_VALUES = [2, 3, 4, 5]

# Grid resolution (reduced for speed; 100 is adequate for eigenvalue convergence)
N_ALPHA = 100


def is_homonuclear() -> bool:
    """Check if system is treated as homonuclear by build_angular_hamiltonian."""
    return Z_A == Z_B


def get_channels(l_max: int) -> list:
    """Get channel list consistent with build_angular_hamiltonian."""
    return _channel_list(l_max, homonuclear=is_homonuclear())


def solve_angular_no_pk(l_max: int) -> tuple:
    """Solve angular problem with NO PK (baseline for nuclear coupling)."""
    h = (np.pi / 2) / (N_ALPHA + 1)
    alpha = (np.arange(N_ALPHA) + 1) * h

    H = build_angular_hamiltonian(
        alpha, rho, R_e_test, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )
    evals, evecs = eigh(H)
    channels = get_channels(l_max)
    return evals, evecs, alpha, channels


def solve_angular_with_pk(l_max: int, pk_channel_mode: str) -> tuple:
    """Solve angular problem with PK in specified channel mode."""
    h = (np.pi / 2) / (N_ALPHA + 1)
    alpha = (np.arange(N_ALPHA) + 1) * h

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_param,
        'atom': 'A',
        'channel_mode': pk_channel_mode,
    }]

    H = build_angular_hamiltonian(
        alpha, rho, R_e_test, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=pk_pots,
    )
    evals, evecs = eigh(H)
    channels = get_channels(l_max)
    return evals, evecs, alpha, channels


def compute_channel_l0_weights(evecs: np.ndarray, channels: list,
                               n_alpha: int) -> np.ndarray:
    """
    Compute the l=0 content of each eigenchannel.

    For each eigenvector, compute the fraction of the wavefunction norm
    that resides in channels where l1=0 or l2=0.

    Returns per-eigenchannel array of shape (n_eig,) with l=0 fractions.
    """
    n_ch = len(channels)
    n_eig = evecs.shape[0] if evecs.ndim == 2 else 1

    # Each eigenvector has shape (n_ch * n_alpha,).
    # Channel ic occupies indices [ic*n_alpha : (ic+1)*n_alpha].
    weights = np.zeros(n_eig)
    for ie in range(n_eig):
        vec = evecs[:, ie] if evecs.ndim == 2 else evecs
        total_norm = 0.0
        l0_norm = 0.0
        for ic, ch in enumerate(channels):
            l1, l2 = ch[0], ch[1]
            chunk = vec[ic * n_alpha: (ic + 1) * n_alpha]
            ch_norm = np.sum(chunk ** 2)
            total_norm += ch_norm
            # Count l=0 content: (delta_{l1,0} + delta_{l2,0}) / 2
            l0_frac = ((1.0 if l1 == 0 else 0.0) + (1.0 if l2 == 0 else 0.0)) / 2.0
            l0_norm += l0_frac * ch_norm
        if total_norm > 1e-30:
            weights[ie] = l0_norm / total_norm
    return weights


def solve_angular_projected_pk(l_max: int) -> tuple:
    """
    Solve angular problem with projected PK.

    Strategy:
    1. Solve WITHOUT PK to get nuclear-coupling eigenstates.
    2. Compute l=0 content of each channel in each eigenstate.
    3. Build a modified PK that weights each channel by its l=0 content
       from the ground-state eigenvector.
    4. Re-solve with these projected PK weights.

    This is a one-shot projection (not self-consistent).
    """
    # Step 1: Solve without PK to get nuclear eigenstates
    evals_nopk, evecs_nopk, alpha, channels = solve_angular_no_pk(l_max)

    n_ch = len(channels)
    n_alpha = len(alpha)

    # Step 2: Compute per-channel l=0 weight from ground-state eigenvector
    ground_vec = evecs_nopk[:, 0]
    per_channel_weight = np.zeros(n_ch)
    for ic, ch in enumerate(channels):
        chunk = ground_vec[ic * n_alpha: (ic + 1) * n_alpha]
        per_channel_weight[ic] = np.sum(chunk ** 2)

    # Normalize
    total = np.sum(per_channel_weight)
    if total > 1e-30:
        per_channel_weight /= total

    # Step 3: Build Hamiltonian with projected PK weights
    # Start with the base Hamiltonian (no PK)
    H_base = build_angular_hamiltonian(
        alpha, rho, R_e_test, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )

    # Compute PK potential profiles
    # Need rho_A, rho_B in R_e units as used by build_angular_hamiltonian
    # The build_angular_hamiltonian computes rho_A, rho_B internally from
    # rho, Z_A, Z_B, z0. Let's replicate that logic.
    R = 2.0 * rho * R_e_test
    rho_A_re = (R / 2.0 + z0) / R_e_test  # distance of A from origin in R_e units
    rho_B_re = (R / 2.0 - z0) / R_e_test  # distance of B from origin in R_e units

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_param,
        'atom': 'A',
    }]
    V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
        alpha, rho_A_re, R_e_test, pk_pots, rho_B=rho_B_re,
    )

    # Step 4: Add projected PK to base Hamiltonian
    H = H_base.copy()

    def idx(ch: int, i: int) -> int:
        return ch * n_alpha + i

    for ic, ch in enumerate(channels):
        l1, l2 = ch[0], ch[1]

        # Per-electron l=0 indicator (same as l-dependent mode)
        has_l0_e1 = 1.0 if l1 == 0 else 0.0
        has_l0_e2 = 1.0 if l2 == 0 else 0.0

        # Project: weight by the channel's contribution to the ground state
        # The idea: channels with more weight in the ground eigenstate
        # should get PK proportional to their l=0 character AND their
        # weight in the eigenstate.
        #
        # Actually, re-reading the task spec more carefully:
        # The projected PK weights each channel by the fractional l=0
        # content of the nuclear coupling eigenvectors.
        #
        # For the ground eigenvector, compute how much l=0 character it has
        # in this specific channel's contribution.
        #
        # The per-channel PK weight is:
        #   w_PK(l1, l2) = (delta_{l1,0} + delta_{l2,0}) / 2
        # But in projected mode, we use the ground-state eigenvector to
        # compute the effective overlap with l=0 character.
        #
        # Simpler interpretation: the PK weight for each channel is the
        # overlap of that channel with the l=0 subspace, weighted by how
        # much the eigenstate occupies that channel.

        # Per-electron weight: fractional l=0 content from eigenvector
        # For electron 1: sum over all channels with l1=0 of their weight
        # For electron 2: sum over all channels with l2=0 of their weight
        pass

    # Recompute: the projected weights for each electron
    # w1_proj = sum of per_channel_weight[ic] for all channels with l1=0
    # w2_proj = sum of per_channel_weight[ic] for all channels with l2=0
    w1_total = sum(per_channel_weight[ic] for ic, ch in enumerate(channels)
                   if ch[0] == 0)
    w2_total = sum(per_channel_weight[ic] for ic, ch in enumerate(channels)
                   if ch[1] == 0)

    # Now apply PK with per-channel projected weights
    for ic, ch in enumerate(channels):
        l1, l2 = ch[0], ch[1]

        # Per-electron projected weight:
        # Electron 1 feels PK proportional to: l1==0 AND the overall l=0
        # fraction of the eigenstate for electron 1's angular momentum.
        # If l1 != 0, this electron is orthogonal to the core regardless.
        # If l1 == 0, the weight is the fraction of the eigenstate in
        # l1=0 channels (= how much l=0 character the eigenstate has).
        w1 = w1_total if l1 == 0 else 0.0
        w2 = w2_total if l2 == 0 else 0.0

        if w1 == 0.0 and w2 == 0.0:
            continue

        for i in range(n_alpha):
            ii = idx(ic, i)
            H[ii, ii] += R_e_test * w1 * V_pk_e1[i]
            H[ii, ii] += R_e_test * w2 * V_pk_e2[i]

    evals, evecs = eigh(H)
    return evals, evecs, alpha, channels


def compute_per_channel_pk_weights(l_max: int) -> dict:
    """
    Compute per-channel PK weights for all three modes at given l_max.

    Returns dict with keys 'channel_blind', 'l_dependent', 'projected',
    each mapping to a list of (channel_label, w1, w2) tuples.
    """
    channels = get_channels(l_max)
    n_ch = len(channels)

    # Mode 1: Channel-blind
    blind = [(str(ch), 1.0, 1.0) for ch in channels]

    # Mode 2: l-dependent (binary)
    ldep = []
    for ch in channels:
        l1, l2 = ch[0], ch[1]
        w1 = 1.0 if l1 == 0 else 0.0
        w2 = 1.0 if l2 == 0 else 0.0
        ldep.append((str(ch), w1, w2))

    # Mode 3: Projected -- need eigenvector from no-PK solve
    evals, evecs, alpha, _ = solve_angular_no_pk(l_max)
    ground_vec = evecs[:, 0]
    n_alpha = len(alpha)

    per_ch_weight = np.zeros(n_ch)
    for ic in range(n_ch):
        chunk = ground_vec[ic * n_alpha: (ic + 1) * n_alpha]
        per_ch_weight[ic] = np.sum(chunk ** 2)
    total = np.sum(per_ch_weight)
    if total > 1e-30:
        per_ch_weight /= total

    w1_total = sum(per_ch_weight[ic] for ic, ch in enumerate(channels)
                   if ch[0] == 0)
    w2_total = sum(per_ch_weight[ic] for ic, ch in enumerate(channels)
                   if ch[1] == 0)

    proj = []
    for ic, ch in enumerate(channels):
        l1, l2 = ch[0], ch[1]
        w1 = w1_total if l1 == 0 else 0.0
        w2 = w2_total if l2 == 0 else 0.0
        proj.append((str(ch), w1, w2))

    return {
        'channel_blind': blind,
        'l_dependent': ldep,
        'projected': proj,
    }


def main():
    print("=" * 72)
    print("Projected PK Pseudopotential Diagnostic -- LiH l_max Convergence")
    print("=" * 72)
    print(f"R = {R_LiH} bohr (experimental LiH R_eq)")
    print(f"R_e = {R_e_test} bohr (fixed test hyperradius)")
    print(f"Z_A_eff = {Z_A}, Z_B = {Z_B}")
    print(f"PK: A = {PK_A} Ha*bohr^2, B = {PK_B_param} bohr^-2")
    print(f"rho = R/(2*R_e) = {rho:.4f}")
    print(f"n_alpha = {N_ALPHA}")
    print()

    # ===================================================================
    # Part 1: mu_0 vs l_max for all three PK modes
    # ===================================================================
    results = {mode: {} for mode in ['channel_blind', 'l_dependent', 'projected']}

    for l_max in L_MAX_VALUES:
        channels = get_channels(l_max)
        n_ch = len(channels)
        print(f"--- l_max = {l_max} ({n_ch} channels: {channels}) ---")

        # Mode 1: Channel-blind
        evals_cb, _, _, _ = solve_angular_with_pk(l_max, 'channel_blind')
        mu0_cb = evals_cb[0]
        results['channel_blind'][l_max] = float(mu0_cb)
        print(f"  Mode 1 (channel-blind):  mu_0 = {mu0_cb:.6f}")

        # Mode 2: l-dependent (binary)
        evals_ld, _, _, _ = solve_angular_with_pk(l_max, 'l_dependent')
        mu0_ld = evals_ld[0]
        results['l_dependent'][l_max] = float(mu0_ld)
        print(f"  Mode 2 (l-dependent):    mu_0 = {mu0_ld:.6f}")

        # Mode 3: Projected
        evals_pr, _, _, _ = solve_angular_projected_pk(l_max)
        mu0_pr = evals_pr[0]
        results['projected'][l_max] = float(mu0_pr)
        print(f"  Mode 3 (projected):      mu_0 = {mu0_pr:.6f}")

        print()

    # ===================================================================
    # Summary table
    # ===================================================================
    print("\n" + "=" * 72)
    print("SUMMARY TABLE: mu_0 vs l_max")
    print("=" * 72)
    print(f"{'l_max':>5}  {'Ch-blind':>12}  {'l-dependent':>12}  {'Projected':>12}"
          f"  {'Delta(cb)':>8}  {'Delta(ld)':>8}  {'Delta(pr)':>8}")
    print("-" * 72)

    ref_cb = results['channel_blind'][L_MAX_VALUES[0]]
    ref_ld = results['l_dependent'][L_MAX_VALUES[0]]
    ref_pr = results['projected'][L_MAX_VALUES[0]]

    for l_max in L_MAX_VALUES:
        cb = results['channel_blind'][l_max]
        ld = results['l_dependent'][l_max]
        pr = results['projected'][l_max]
        dcb = cb - ref_cb
        dld = ld - ref_ld
        dpr = pr - ref_pr
        print(f"{l_max:>5}  {cb:>12.6f}  {ld:>12.6f}  {pr:>12.6f}"
              f"  {dcb:>+8.4f}  {dld:>+8.4f}  {dpr:>+8.4f}")

    print("-" * 72)
    print("Delta = mu_0(l_max) - mu_0(l_max=2). Positive Delta -> PES shifts outward (divergent).")
    print("                               Negative Delta -> PES shifts inward (convergent).")

    # ===================================================================
    # Convergence assessment
    # ===================================================================
    print("\n" + "=" * 72)
    print("CONVERGENCE ASSESSMENT")
    print("=" * 72)

    for mode_name, mode_key in [("Channel-blind", "channel_blind"),
                                  ("l-dependent", "l_dependent"),
                                  ("Projected", "projected")]:
        vals = [results[mode_key][lm] for lm in L_MAX_VALUES]
        deltas = [v - vals[0] for v in vals]
        monotonic_up = all(deltas[i+1] >= deltas[i] for i in range(len(deltas)-1))
        monotonic_down = all(deltas[i+1] <= deltas[i] for i in range(len(deltas)-1))
        spread = max(vals) - min(vals)

        if monotonic_up and spread > 0.1:
            verdict = "DIVERGENT (mu_0 increasing monotonically)"
        elif monotonic_down and spread > 0.1:
            verdict = "CONVERGENT (mu_0 decreasing monotonically)"
        elif spread < 0.01:
            verdict = "STABLE (mu_0 variation < 0.01)"
        elif spread < 0.1:
            verdict = "WEAKLY CONVERGENT (mu_0 variation < 0.1)"
        else:
            verdict = f"OSCILLATORY (spread = {spread:.4f})"

        print(f"  {mode_name:15s}: {verdict}")

    # ===================================================================
    # Part 2: Per-channel PK weights at l_max=4
    # ===================================================================
    print("\n" + "=" * 72)
    print("PER-CHANNEL PK WEIGHTS at l_max = 4")
    print("=" * 72)

    weights = compute_per_channel_pk_weights(4)
    channels_4 = get_channels(4)

    print(f"{'Channel':>12}  {'Blind(w1,w2)':>14}  {'l-dep(w1,w2)':>14}  {'Proj(w1,w2)':>14}")
    print("-" * 60)
    for i, ch in enumerate(channels_4):
        cb_w = weights['channel_blind'][i]
        ld_w = weights['l_dependent'][i]
        pr_w = weights['projected'][i]
        print(f"{str(ch):>12}  ({cb_w[1]:.2f},{cb_w[2]:.2f})      "
              f"({ld_w[1]:.2f},{ld_w[2]:.2f})      "
              f"({pr_w[1]:.2f},{pr_w[2]:.2f})")

    # ===================================================================
    # Save data
    # ===================================================================
    output = {
        'parameters': {
            'R': R_LiH,
            'R_e': R_e_test,
            'Z_A_eff': Z_A,
            'Z_B': Z_B,
            'PK_A': PK_A,
            'PK_B': PK_B_param,
            'rho': rho,
            'n_alpha': N_ALPHA,
        },
        'mu0_vs_lmax': results,
        'per_channel_weights_lmax4': {
            mode: [(ch, w1, w2) for ch, w1, w2 in wlist]
            for mode, wlist in weights.items()
        },
    }

    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)
    with open(data_dir / 'projected_pk_diagnostic.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nData saved to debug/data/projected_pk_diagnostic.json")

    # ===================================================================
    # Plots
    # ===================================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plot_dir = Path(__file__).parent / 'plots'
        plot_dir.mkdir(exist_ok=True)

        # --- Plot A: mu_0 vs l_max ---
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        for mode_label, mode_key, color, marker in [
            ('Channel-blind', 'channel_blind', '#d62728', 'o'),
            ('l-dependent', 'l_dependent', '#1f77b4', 's'),
            ('Projected', 'projected', '#2ca02c', '^'),
        ]:
            vals = [results[mode_key][lm] for lm in L_MAX_VALUES]
            ax.plot(L_MAX_VALUES, vals, f'-{marker}', color=color,
                    label=mode_label, linewidth=2, markersize=8)

        ax.set_xlabel('$l_{\\rm max}$', fontsize=13)
        ax.set_ylabel('$\\mu_0$ (ground angular eigenvalue)', fontsize=13)
        ax.set_title(
            f'LiH Projected PK Diagnostic: $\\mu_0$ vs $l_{{\\rm max}}$\n'
            f'$R = {R_LiH}$ bohr, $R_e = {R_e_test}$ bohr, '
            f'$Z_{{A,eff}} = {Z_A}$, $Z_B = {Z_B}$',
            fontsize=11,
        )
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(L_MAX_VALUES)
        fig.tight_layout()
        fig.savefig(plot_dir / 'projected_pk_mu0_vs_lmax.png', dpi=150)
        plt.close(fig)
        print(f"Plot A saved to debug/plots/projected_pk_mu0_vs_lmax.png")

        # --- Plot B: Per-channel PK weights at l_max=4 ---
        fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
        ch_labels = [str(ch) for ch in channels_4]
        x = np.arange(len(ch_labels))
        bar_w = 0.35

        for ax, (mode_label, mode_key) in zip(axes, [
            ('Channel-blind', 'channel_blind'),
            ('l-dependent', 'l_dependent'),
            ('Projected', 'projected'),
        ]):
            w1s = [weights[mode_key][i][1] for i in range(len(channels_4))]
            w2s = [weights[mode_key][i][2] for i in range(len(channels_4))]
            ax.bar(x - bar_w/2, w1s, bar_w, label='$w_1$ (electron 1)',
                   color='#1f77b4', alpha=0.8)
            ax.bar(x + bar_w/2, w2s, bar_w, label='$w_2$ (electron 2)',
                   color='#ff7f0e', alpha=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(ch_labels, rotation=45, ha='right', fontsize=8)
            ax.set_title(mode_label, fontsize=12)
            ax.set_xlabel('Channel $(l_1, l_2)$', fontsize=10)
            if ax is axes[0]:
                ax.set_ylabel('PK weight', fontsize=11)
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, alpha=0.2, axis='y')

        fig.suptitle(f'Per-channel PK weights at $l_{{\\rm max}} = 4$', fontsize=13)
        fig.tight_layout()
        fig.savefig(plot_dir / 'projected_pk_channel_weights.png', dpi=150)
        plt.close(fig)
        print(f"Plot B saved to debug/plots/projected_pk_channel_weights.png")

    except ImportError:
        print("matplotlib not available -- skipping plots.")


if __name__ == '__main__':
    main()
