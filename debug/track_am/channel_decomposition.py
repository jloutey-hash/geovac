"""Channel decomposition analysis at l_max=2 for the 4-electron LiH solver.

Extracts the angular eigenvector at specific R values and decomposes it into
raw (l1,l2,l3,l4) channel contributions.
"""
import sys
sys.path.insert(0, 'c:/Users/jlout/OneDrive/Desktop/Project_Geometric')

import numpy as np
import time
from geovac.n_electron_solver import (
    solve_angular_4e_multichannel,
    _enumerate_channels_4e,
    _build_s4_22_projector,
    _s4_channel_basis,
)


def channel_weights_at_R(R, R_e, l_max=2, n_grid=6):
    """Compute channel weight decomposition at a given (R, R_e) point.

    Returns dict with channel labels and their weights.
    """
    Z_A, Z_B = 3.0, 1.0
    z0 = 0.0  # midpoint origin

    R_A = R / 2.0 - z0
    R_B = R / 2.0 + z0
    rho_A = R_A / R_e
    rho_B = R_B / R_e

    channels = _enumerate_channels_4e(l_max)
    n_ch = len(channels)

    # Get the S4 [2,2] channel basis
    ch_basis = _s4_channel_basis(channels)  # (n_ch, rank)
    rank = ch_basis.shape[1]

    # Solve for eigenvector
    evals, evecs = solve_angular_4e_multichannel(
        rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
        s4_projection=True, symmetry='s4', n_states=1,
    )

    if evecs is None:
        return None

    psi = evecs[0]  # shape: (rank * n_grid^3,)
    n_grid_total = n_grid ** 3

    # Reshape to (rank, n_grid^3)
    psi_block = psi.reshape(rank, n_grid_total)

    # Weight of S4 basis vector p: integral of |psi_p(grid)|^2
    s4_weights = np.sum(psi_block**2, axis=1)  # (rank,)

    # Map back to raw channels: weight of raw channel a
    # psi_raw_a = sum_p ch_basis[a, p] * psi_p
    # weight_a = sum_grid |psi_raw_a|^2 = sum_p,q ch_basis[a,p]*ch_basis[a,q] * (psi_p . psi_q)
    # For simplicity, compute the overlap matrix in S4 space
    overlap = psi_block @ psi_block.T  # (rank, rank)

    raw_weights = {}
    for a, ch in enumerate(channels):
        w = 0.0
        for p in range(rank):
            for q in range(rank):
                w += ch_basis[a, p] * ch_basis[a, q] * overlap[p, q]
        if abs(w) > 1e-10:
            raw_weights[ch] = w

    return {
        'mu': evals[0],
        'R': R,
        'R_e': R_e,
        'rank': rank,
        'n_ch': n_ch,
        's4_weights': s4_weights,
        'raw_weights': raw_weights,
        'channels': channels,
        'ch_basis': ch_basis,
    }


def main():
    l_max = 2
    n_grid = 6

    # R_e values near the minimum of the adiabatic curve (from l_max=2 results)
    # The adiabatic curve minimum is typically at R_e ~ 2-3 bohr
    R_e_target = 2.5  # approximate R_e at the adiabatic minimum

    R_values = [1.0, 1.5, 3.0]

    print("=" * 70)
    print(f"Channel decomposition at l_max={l_max}, n_grid={n_grid}")
    print("=" * 70)

    channels = _enumerate_channels_4e(l_max)
    ch_basis = _s4_channel_basis(channels)
    rank = ch_basis.shape[1]

    print(f"\nS4 [2,2] subspace: {rank} basis vectors from {len(channels)} raw channels")
    print(f"\nS4 basis vectors (channel decomposition):")
    for p in range(rank):
        nonzero = [(channels[a], ch_basis[a, p])
                   for a in range(len(channels)) if abs(ch_basis[a, p]) > 0.05]
        nonzero.sort(key=lambda x: -abs(x[1]))
        print(f"  S4 vector {p}: ", end="")
        for ch, w in nonzero[:6]:
            print(f"  {ch}:{w:.3f}", end="")
        print()

    for R in R_values:
        print(f"\n{'='*60}")
        print(f"R = {R:.3f} bohr")
        print(f"{'='*60}")

        t0 = time.time()
        result = channel_weights_at_R(R, R_e_target, l_max=l_max, n_grid=n_grid)
        t1 = time.time()

        if result is None:
            print("  No solution found (zero dim)")
            continue

        print(f"  mu = {result['mu']:.4f}")
        print(f"  Time: {t1-t0:.1f}s")

        # Sort raw weights
        sorted_weights = sorted(result['raw_weights'].items(), key=lambda x: -x[1])
        total_w = sum(w for _, w in sorted_weights)

        print(f"\n  Channel weights (total={total_w:.6f}):")
        print(f"  {'Channel':>20s}  {'Weight':>10s}  {'% total':>8s}")
        print(f"  {'-'*20}  {'-'*10}  {'-'*8}")
        cumulative = 0.0
        for ch, w in sorted_weights[:20]:
            pct = w / total_w * 100
            cumulative += pct
            print(f"  {str(ch):>20s}  {w:10.6f}  {pct:7.2f}%  (cum: {cumulative:.1f}%)")

        # Group by total L = l1+l2+l3+l4
        L_weights = {}
        for ch, w in result['raw_weights'].items():
            L = sum(ch)
            L_weights[L] = L_weights.get(L, 0.0) + w

        print(f"\n  Weight by total L = l1+l2+l3+l4:")
        for L in sorted(L_weights.keys()):
            pct = L_weights[L] / total_w * 100
            print(f"    L={L}: {L_weights[L]:.6f} ({pct:.2f}%)")

        # Group by max l
        lmax_weights = {}
        for ch, w in result['raw_weights'].items():
            ml = max(ch)
            lmax_weights[ml] = lmax_weights.get(ml, 0.0) + w

        print(f"\n  Weight by max(l_i):")
        for ml in sorted(lmax_weights.keys()):
            pct = lmax_weights[ml] / total_w * 100
            print(f"    max_l={ml}: {lmax_weights[ml]:.6f} ({pct:.2f}%)")


if __name__ == '__main__':
    main()
