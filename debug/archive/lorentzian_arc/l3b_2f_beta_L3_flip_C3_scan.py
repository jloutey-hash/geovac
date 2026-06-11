"""Scan: compute per-harmonic flip-Lichnerowicz constant
||[D_GV, M_spat_flip_N]||_op / ||W||_op  vs N over (2,3) and (3,5).

Goal: find a clean closed-form bound (analogous to Paper 38 L3's
C3^(N) = sqrt((N-1)/(N+1)) for natural generators).
"""

from __future__ import annotations

import numpy as np
import json
from collections import defaultdict
from pathlib import Path

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import CompactTemporalTruncatedOperatorSystem


def op_norm(A):
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def scan(n_max, N_t, T=2*np.pi):
    print(f"\n=== panel ({n_max}, {N_t}) ===")
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(K)
    J = K.J
    J_inv = np.linalg.inv(J)
    D_L_diag = 0.5 * (D_L + J @ D_L @ J_inv)
    D_L_off = 0.5 * (D_L - J @ D_L @ J_inv)

    # Recover D_GV ⊗ I_{N_t} via D_L_off / i. Actually D_L = i*(gamma0 ⊗ D_t + D_GV ⊗ I).
    # So D_L_off = i * D_GV ⊗ I_{N_t}.
    # We can extract D_GV from D_L_off by averaging over the temporal block.

    # D_GV acts on spatial part only.
    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    d_w = O.dim_spatial // 2
    d_s = O.dim_spatial

    # max-N for the operator system
    max_N = max(N for (N, L, M) in O.spat_labels)
    print(f"max spatial label N = {max_N}")
    print(f"||D_GV||_op = {op_norm(D_L_off):.4f}")

    # For each spatial generator, compute Paper 38-style ratio:
    # nat: ||[D_GV (spatial only), M_spat_nat]||_op / ||M_spat_nat||_op
    # flip: ||[D_GV (spatial only), M_spat_flip]||_op / ||M_spat_flip||_op

    # We need D_GV restricted to spatial. Build D_GV from spatial-side data.
    # D_L_off lives on dim_K = d_s * N_t. The off-block satisfies
    # D_L_off = i * D_GV ⊗ I_{N_t}. So restrict to a single temporal mode k=0:
    # take first d_s rows/cols of the (d_s*N_t, d_s*N_t) matrix, but actually
    # we need the block at temporal index 0 only.
    # The pure-tensor layout is kron(spat, temp), so for temporal index k = 0
    # (which is the central index for our symmetric construction):
    # D_GV_spat = - i * D_L_off[k*d_s:(k+1)*d_s, k*d_s:(k+1)*d_s].
    # Actually, since D_L_off = i * D_GV ⊗ I, the full matrix has block-diag
    # structure on the temporal index. We can read off D_GV as the (0,0) block.
    # Layout from build_chirality_flipping_generators is kron(M_spat, M_temp),
    # so the temporal index varies fastest. We need to extract the
    # spatial-only block. Easiest: D_GV_spat = - i * D_L_off[:d_s, :d_s] / ||M_temp_0[0,0]||
    # but M_temp_0 = diag(omega_k^0) = I.
    # Wait, layout: np.kron(M_spat, M_temp) has structure
    #   [M_spat[0,0] * M_temp,  M_spat[0,1] * M_temp, ...]
    # So for I_{N_t}, kron(D_GV, I_{N_t}) has block structure
    #   block[i,j] = D_GV[i,j] * I_{N_t}.
    # Hence D_GV[i,j] = D_L_off[i*N_t, j*N_t] / i (for each i, j in spatial).
    # Easier: average over diagonal to extract D_GV.
    D_GV_spat = np.zeros((d_s, d_s), dtype=np.complex128)
    for i in range(d_s):
        for j in range(d_s):
            # Average over k
            block_diag = np.zeros(N_t, dtype=np.complex128)
            for k in range(N_t):
                block_diag[k] = D_L_off[i*N_t + k, j*N_t + k] / 1j
            D_GV_spat[i, j] = block_diag.mean()
    print(f"||D_GV_spat||_op = {op_norm(D_GV_spat):.4f}")

    # Per-harmonic ratios
    by_N_nat = defaultdict(list)
    by_N_flip = defaultdict(list)

    for i, (N, L, M) in enumerate(O.spat_labels):
        M_spat_nat = O._spat_matrices[i]
        W = M_spat_nat[:d_w, :d_w].copy()
        M_spat_flip = np.zeros_like(M_spat_nat)
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        W_norm = op_norm(W)
        M_spat_flip_norm = op_norm(M_spat_flip)
        # By construction, M_spat_flip has same op_norm as W (it's block-diag).
        if W_norm == 0:
            continue

        # Spatial commutator (using extracted D_GV_spat):
        comm_nat = D_GV_spat @ M_spat_nat - M_spat_nat @ D_GV_spat
        comm_flip = D_GV_spat @ M_spat_flip - M_spat_flip @ D_GV_spat

        ratio_nat = op_norm(comm_nat) / op_norm(M_spat_nat) if op_norm(M_spat_nat) > 0 else 0
        ratio_flip = op_norm(comm_flip) / op_norm(M_spat_flip) if op_norm(M_spat_flip) > 0 else 0

        by_N_nat[N].append((L, M, ratio_nat))
        by_N_flip[N].append((L, M, ratio_flip))

    print(f"\n{'N':4s} {'#gens':6s} {'max(nat_ratio)':16s} {'C3_38(N)=sqrt((N-1)/(N+1))':30s} {'max(flip_ratio)':16s} {'flip/nat':10s}")
    for N in sorted(by_N_nat.keys()):
        max_nat = max(r for (_, _, r) in by_N_nat[N])
        max_flip = max(r for (_, _, r) in by_N_flip[N])
        C3_38 = np.sqrt(max((N-1)/(N+1), 0)) if N >= 1 else 0
        ratio_fn = max_flip / max_nat if max_nat > 1e-10 else float('inf')
        print(f"{N:4d} {len(by_N_nat[N]):6d} {max_nat:16.4f} {C3_38:30.4f} {max_flip:16.4f} {ratio_fn:10.4f}")

    return by_N_nat, by_N_flip


print("==" * 40)
scan(2, 3)
print("==" * 40)
scan(3, 5)
