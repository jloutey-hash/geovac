"""Diagnostic: compute Term A, Term B for the saturating flip generators,
and check what the correct Lichnerowicz constants are.
"""

from __future__ import annotations

import numpy as np
import json
from pathlib import Path

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import CompactTemporalTruncatedOperatorSystem


def op_norm(A):
    if A.size == 0:
        return 0.0
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def diagnose(n_max, N_t, T=2*np.pi):
    print(f"\n--- panel ({n_max}, {N_t}) ---")

    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L = lorentzian_dirac_compact_matrix(K)
    J = K.J
    J_inv = np.linalg.inv(J)
    D_L_diag = 0.5 * (D_L + J @ D_L @ J_inv)  # block-diagonal under J (gamma^0 ⊗ D_t piece)
    D_L_off = 0.5 * (D_L - J @ D_L @ J_inv)   # off-block (D_GV ⊗ I piece)

    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    d_w = O.dim_spatial // 2

    # Look at first few natural and flip generators.
    print(f"||D_L_diag|| = {op_norm(D_L_diag):.4f}, ||D_L_off|| = {op_norm(D_L_off):.4f}")

    # For each natural generator at p=0:
    print(f"{'label':12s} {'||W||':8s} {'||M_t_0||':10s} {'[D_GV,M_nat]':14s} {'[D_GV,M_flip]':14s} {'[gamma0,M_flip]_op':18s} {'TermA':10s} {'TermB_nat':10s} {'TermB_flip':10s} {'L_op_nat':10s} {'L_op_flip':10s}")
    for i, (N, L, M) in enumerate(O.spat_labels[:5]):
        M_spat_nat = O._spat_matrices[i]
        W = M_spat_nat[:d_w, :d_w].copy()
        M_spat_flip = np.zeros_like(M_spat_nat)
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        W_norm = op_norm(W)

        for p in range(min(2, N_t)):
            M_temp_p = O._temp_matrices[p]
            M_temp_norm = op_norm(M_temp_p)

            # Natural and flip multipliers:
            a_nat = np.kron(M_spat_nat, M_temp_p)
            a_flip = np.kron(M_spat_flip, M_temp_p)

            # Term A for flip: [D_L_diag, a_flip]
            term_A_flip = op_norm(D_L_diag @ a_flip - a_flip @ D_L_diag)
            term_B_flip = op_norm(D_L_off @ a_flip - a_flip @ D_L_off)
            L_op_flip = op_norm(D_L @ a_flip - a_flip @ D_L)

            # Term B for natural: [D_L_off, a_nat]
            term_A_nat = op_norm(D_L_diag @ a_nat - a_nat @ D_L_diag)
            term_B_nat = op_norm(D_L_off @ a_nat - a_nat @ D_L_off)
            L_op_nat = op_norm(D_L @ a_nat - a_nat @ D_L)

            # Bare gamma0, M_flip anti-commutator (factor in Term A):
            # gamma0 acts on the spatial 2-block: gamma0_spat * diag(W, -W) + diag(W, -W) * gamma0_spat
            # = [[0, I], [I, 0]] * diag(W, -W) = [[0, -W], [W, 0]]
            # + diag(W, -W) * [[0, I], [I, 0]] = [[0, W], [-W, 0]]
            # = 0 (the anti-commutator). But the commutator is [gamma0, diag(W, -W)] = [[0, I], [I, 0]] diag(W, -W) - diag(W, -W) [[0, I], [I, 0]]
            # = [[0, -W], [W, 0]] - [[0, W], [-W, 0]]
            # = [[0, -2W], [2W, 0]]
            # ||[gamma0, M_spat_flip]||_op = 2 ||W||_op
            gamma0 = np.kron(np.array([[0, 1], [1, 0]]), np.eye(d_w))
            comm_g0_Mflip = gamma0 @ M_spat_flip - M_spat_flip @ gamma0
            n_comm_g0_Mflip = op_norm(comm_g0_Mflip)

            print(f"({N},{L},{M},p={p}) {W_norm:.4f}   {M_temp_norm:.4f}    {op_norm(np.kron(M_spat_nat, M_temp_p)*0 + 0):.4f}        {op_norm(D_L_off @ a_flip - a_flip @ D_L_off):.4f}        {n_comm_g0_Mflip:.4f}            {term_A_flip:.4f}   {term_B_nat:.4f}    {term_B_flip:.4f}    {L_op_nat:.4f}   {L_op_flip:.4f}")

    # NEW: check the ratio L_op(a^flip) / [(2*||D_t|| + ||[D_GV,M^flip]||/||W||) * ||W|| * ||M_temp||]
    # This is a SHARPER bound: per-flip-generator C3.
    print("\n--- Sharper bound (per-flip-generator C3, includes ||[D_GV, M_spat_flip]||) ---")
    print(f"{'label':18s} {'L_op_flip':10s} {'2||D_t||*||W||*||Mt||':23s} {'||[D_GV,M_flip]||*||Mt||':25s} {'sum':10s} {'ratio':10s}")
    for i, (N, L, M) in enumerate(O.spat_labels[:10]):
        M_spat_nat = O._spat_matrices[i]
        W = M_spat_nat[:d_w, :d_w].copy()
        M_spat_flip = np.zeros_like(M_spat_nat)
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        W_norm = op_norm(W)

        # Use p=0 for canonical example.
        for p in [0, 1]:
            if p >= N_t:
                continue
            M_temp_p = O._temp_matrices[p]
            M_temp_norm = op_norm(M_temp_p)
            a_flip = np.kron(M_spat_flip, M_temp_p)

            # Direct LHS:
            L_op_flip = op_norm(D_L @ a_flip - a_flip @ D_L)

            # Detailed contributions:
            Dt_norm = op_norm(D_L_diag)  # = ||gamma0 ⊗ D_t|| = ||D_t||
            # ||[D_GV, M_spat_flip]||_op: extract D_GV from D_L_off / I
            # We don't have D_GV directly; use D_L_off ⊗ I^{-1} not well-defined.
            # Instead measure [D_L_off, a_flip] / ||M_temp_p||:
            comm_off = D_L_off @ a_flip - a_flip @ D_L_off
            comm_off_norm = op_norm(comm_off)
            # The factor we want: ||[D_GV, M_spat_flip]||_op (single tensor factor).
            # Since D_L_off = i * D_GV ⊗ I_{N_t}, [D_L_off, M_spat_flip o M_temp] = i [D_GV, M_spat_flip] o M_temp.
            # So ||comm_off|| = ||[D_GV, M_spat_flip]||_op * ||M_temp||_op.
            DGV_comm_flip_norm = comm_off_norm / M_temp_norm if M_temp_norm > 0 else 0

            # The actual Lichnerowicz constant we need:
            # L_op(a^flip) <= [2 ||D_t|| + DGV_comm_flip_norm/||W||] * ||W|| * ||M_temp||
            # i.e. C3 ~ 2 ||D_t||_op + (DGV_comm_flip_norm / ||W||)
            # Note (DGV_comm_flip_norm / ||W||) is a per-spatial-label per-harmonic constant
            # analogous to C3^(N) but for the chirality-flip multiplier.

            term_2Dt = 2.0 * Dt_norm * W_norm * M_temp_norm
            term_DGV = DGV_comm_flip_norm * M_temp_norm
            sum_bound = term_2Dt + term_DGV
            ratio = L_op_flip / sum_bound if sum_bound > 0 else float('inf')

            print(f"({N},{L},{M},p={p}) {L_op_flip:.4f}    {term_2Dt:.4f}                 {term_DGV:.4f}                  {sum_bound:.4f}  {ratio:.4f}")

    return D_L, D_L_diag, D_L_off, O


print("=" * 80)
diagnose(2, 3)
print("=" * 80)
diagnose(3, 5)
