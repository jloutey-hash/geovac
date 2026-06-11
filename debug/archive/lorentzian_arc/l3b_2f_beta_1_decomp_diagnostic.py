"""Decompose [D_L, a^flip] more finely to refine the L3b-2f-alpha §4.2
analytical prediction at p=0.

Per the alpha §4.2 reasoning:
  [D_L, a^flip_p] = i [gamma^0, M^spat_flip] (x) (D_t M^temp_p)
                  + i [D_GV, M^spat_flip] (x) M^temp_p
(both other terms vanish: [D_t, M^temp_p] = 0 since diagonal, and
[I, M^temp_p] = 0).

At p=0: M^temp_0 = I_{N_t}, so D_t M^temp_0 = D_t, which is non-zero
(diagonal of 2 pi i k / T). The alpha §4.2 claim that this term vanishes
was wrong; it confused "d/dt of the constant function" with the matrix
product D_t * I.

This driver verifies the corrected decomposition numerically.
"""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path
from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    lorentzian_dirac_compact_matrix,
    fourier_d_dt_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
    compact_temporal_multiplier_matrices,
)
from geovac.full_dirac_operator_system import camporesi_higuchi_full_dirac_matrix


def op_norm(A):
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def fro_norm(A):
    return float(np.linalg.norm(A))


n_max = 2
N_t = 3
T = 2.0 * np.pi

K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)

# Build the building blocks.
d_w = O.dim_spatial // 2
gamma0 = K.J_spatial  # = [[0, I], [I, 0]]
D_GV = camporesi_higuchi_full_dirac_matrix(K.basis_spatial)
D_t = fourier_d_dt_matrix(N_t, T)
I_t = np.eye(N_t, dtype=np.complex128)
I_spat = np.eye(O.dim_spatial, dtype=np.complex128)

print(f"D_t diagonal = {np.diag(D_t)}")
print(f"||D_t||_op = {op_norm(D_t):.4f}")

# Pick a spatial generator and build M^spat_flip = diag(W, -W).
for spat_idx in range(min(3, len(O._spat_matrices))):
    full_natural = O._spat_matrices[spat_idx]
    W = full_natural[:d_w, :d_w]
    M_spat_flip = np.zeros((O.dim_spatial, O.dim_spatial), dtype=np.complex128)
    M_spat_flip[:d_w, :d_w] = W
    M_spat_flip[d_w:, d_w:] = -W

    norm_W = op_norm(W)
    print(f"\nGenerator {O.spat_labels[spat_idx]} (||W||={norm_W:.4f}):")

    # Commutators.
    comm_gamma0_Mflip = gamma0 @ M_spat_flip - M_spat_flip @ gamma0
    comm_DGV_Mflip = D_GV @ M_spat_flip - M_spat_flip @ D_GV
    print(f"  ||[gamma^0, M^flip]||_op = {op_norm(comm_gamma0_Mflip):.4f}")
    print(f"  ||[D_GV, M^flip]||_op = {op_norm(comm_DGV_Mflip):.4f}")

    # For each p, compute the two-term decomposition.
    temp_mats = compact_temporal_multiplier_matrices(N_t, T)
    for p in range(N_t):
        M_temp_p = temp_mats[p]
        # Term A: i [gamma^0, M^spat_flip] (x) (D_t M^temp_p)
        # Term B: i [D_GV, M^spat_flip] (x) M^temp_p
        # Each kron is between spatial (dim_spat^2) and temporal (N_t^2).
        Dt_Mp = D_t @ M_temp_p
        term_A = 1j * np.kron(comm_gamma0_Mflip, Dt_Mp)
        term_B = 1j * np.kron(comm_DGV_Mflip, M_temp_p)
        full_comm = term_A + term_B

        # Cross-check: compute [D_L, a^flip_p] directly.
        a_flip_p = np.kron(M_spat_flip, M_temp_p)
        D_L = lorentzian_dirac_compact_matrix(K)
        actual = D_L @ a_flip_p - a_flip_p @ D_L
        residual = fro_norm(full_comm - actual)

        norm_A = op_norm(term_A)
        norm_B = op_norm(term_B)
        norm_full = op_norm(full_comm)
        norm_actual = op_norm(actual)
        print(
            f"  p={p}: ||term_A (time-piece)||={norm_A:.4f}, "
            f"||term_B (space-piece)||={norm_B:.4f}, "
            f"||sum||={norm_full:.4f}, ||actual||={norm_actual:.4f}, "
            f"reconstruction_residual={residual:.3e}"
        )
        # At p=0, M_temp_0 = I, D_t M_temp_0 = D_t, ||D_t||~?
        if p == 0:
            print(f"      D_t M_temp_0 = D_t, ||D_t||={op_norm(D_t):.4f}")

# Also dump key quantities to JSON for memo.
out = {
    "sprint": "L3b-2f-beta.1-decomp",
    "n_max": n_max, "N_t": N_t, "T": T,
    "D_t_diagonal": [complex(x).imag for x in np.diag(D_t)],
    "D_t_op_norm": op_norm(D_t),
    "structural_finding": (
        "Alpha §4.2 claim that [D_L_diag, a^flip_p=0] = 0 is WRONG. "
        "At p=0, M^temp_0 = I, but D_t * M^temp_0 = D_t (non-zero "
        "diagonal of 2 pi i k / T). The alpha §4.2 prose confused "
        "'d/dt applied to the constant function = 0' with the actual "
        "matrix product D_t * M^temp_p, which at p=0 equals D_t. "
        "The structural-identity break occurs at ALL p, not just p >= 1."
    ),
}
out_file = Path(__file__).parent / "data" / "l3b_2f_beta_1_decomp.json"
out_file.parent.mkdir(parents=True, exist_ok=True)
with open(out_file, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nWrote {out_file}")
