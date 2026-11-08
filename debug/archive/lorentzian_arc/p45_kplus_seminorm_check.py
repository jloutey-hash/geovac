"""P45 hardening check: does the K+ compression annihilate the Lipschitz seminorm?

Independent main-session verification of the Phase-1b adversarial finding
(debug/p45_adversarial_L5_memo.md, vector V1), run at the Paper 45 panel
cell (n_max, N_t) = (2, 3) and cross-checked at (3, 5).

Checks:
  C1. {J_spatial, D_GV} = 0   (forced by Krein-self-adjointness of
      D_L = i*(J_s (x) D_t + D_GV (x) I); the i-prefactor on the Hermitian
      spatial Dirac requires J-anticommutation).
  C2. P+ (i * D_GV (x) I_t) P+ = 0   (compression annihilates the spatial
      Dirac block).
  C3. K+-restricted Lipschitz seminorm s_restr(a) = ||[P+ D_L P+, P+ a P+]||
      for EVERY operator-system basis multiplier a.  Adversarial claim:
      identically zero on the entire system (kernel = O, not C*1).
  C4. Full-space seminorm s_full(a) = ||[D_L, a]||: count kernel elements.
      Adversarial claim: every pure-temporal multiplier (spatial label
      trivial) is already in the kernel BEFORE restriction.

Verdict criteria: claims confirmed iff max_a s_restr(a) == 0 to float
precision while max_a s_full(a) > 0 (i.e., the unrestricted system does
see the spatial Dirac, and the restriction kills everything).
"""

from __future__ import annotations

import numpy as np

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    lorentzian_dirac_compact_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)


def run_cell(n_max: int, N_t: int, T: float = 2.0 * np.pi) -> dict:
    op_sys = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    krein = op_sys.krein
    D_L = lorentzian_dirac_compact_matrix(krein)

    dim_K = krein.dim
    I_K = np.eye(dim_K, dtype=np.complex128)
    J = krein.J
    P_plus = 0.5 * (I_K + J)

    # --- C1: spatial anticommutation -----------------------------------
    # Reconstruct the spatial part of D_L: D_L = i*(time + space) with
    # time = J_s (x) D_t and space = D_GV (x) I_t.  Recover space part by
    # subtracting the time part built from krein attributes.
    from geovac.lorentzian_dirac_compact import fourier_d_dt_matrix

    D_t = fourier_d_dt_matrix(N_t, T)
    time_part = np.kron(krein.J_spatial, D_t)
    space_part = D_L / 1j - time_part  # = D_GV (x) I_t
    D_GV_kron = space_part

    anti = krein.J @ D_GV_kron + D_GV_kron @ krein.J
    c1_residual = float(np.linalg.norm(anti))

    # --- C2: compression kills spatial Dirac ---------------------------
    c2_residual = float(np.linalg.norm(P_plus @ (1j * D_GV_kron) @ P_plus))
    c2_scale = float(np.linalg.norm(1j * D_GV_kron))

    # --- C3 / C4: seminorms over the full multiplier basis -------------
    s_full = []
    s_restr = []
    D_restr = P_plus @ D_L @ P_plus
    for label, a in zip(op_sys.multiplier_labels, op_sys.multiplier_matrices):
        full = float(np.linalg.norm(D_L @ a - a @ D_L, 2))
        a_r = P_plus @ a @ P_plus
        restr = float(np.linalg.norm(D_restr @ a_r - a_r @ D_restr, 2))
        s_full.append((label, full))
        s_restr.append((label, restr))

    full_vals = np.array([v for _, v in s_full])
    restr_vals = np.array([v for _, v in s_restr])

    tol = 1e-12
    n_kernel_full = int(np.sum(full_vals < tol))
    n_kernel_restr = int(np.sum(restr_vals < tol))
    n_total = len(full_vals)

    return {
        "cell": (n_max, N_t),
        "dim_K": dim_K,
        "n_multipliers": n_total,
        "C1_anticommutator_residual": c1_residual,
        "C2_compressed_spatial_dirac_norm": c2_residual,
        "C2_uncompressed_spatial_dirac_norm": c2_scale,
        "C3_max_restricted_seminorm": float(restr_vals.max()),
        "C3_n_kernel_restricted": n_kernel_restr,
        "C4_max_full_seminorm": float(full_vals.max()),
        "C4_n_kernel_full": n_kernel_full,
    }


def main() -> None:
    print("=" * 72)
    print("P45 K+ seminorm annihilation check (Phase 1b independent verify)")
    print("=" * 72)
    confirmed_cells = 0
    cells = [(2, 3), (3, 5)]
    for n_max, N_t in cells:
        r = run_cell(n_max, N_t)
        print(f"\nCell (n_max, N_t) = {r['cell']}   dim_K = {r['dim_K']}, "
              f"multipliers = {r['n_multipliers']}")
        print(f"  C1 ||{{J, D_GV(x)I}}||          = {r['C1_anticommutator_residual']:.3e}"
              f"   (claim: 0)")
        print(f"  C2 ||P+ (i D_GV(x)I) P+||     = {r['C2_compressed_spatial_dirac_norm']:.3e}"
              f"   (uncompressed: {r['C2_uncompressed_spatial_dirac_norm']:.3e})")
        print(f"  C3 max_a ||[P+ D_L P+, a+]||  = {r['C3_max_restricted_seminorm']:.3e}"
              f"   kernel: {r['C3_n_kernel_restricted']}/{r['n_multipliers']}"
              f"   (claim: kernel = ALL)")
        print(f"  C4 max_a ||[D_L, a]||         = {r['C4_max_full_seminorm']:.3e}"
              f"   kernel: {r['C4_n_kernel_full']}/{r['n_multipliers']}")
        confirmed = (
            r["C1_anticommutator_residual"] < 1e-12
            and r["C2_compressed_spatial_dirac_norm"] < 1e-12
            and r["C3_max_restricted_seminorm"] < 1e-12
            and r["C4_max_full_seminorm"] > 1e-6
        )
        print(f"  --> adversarial claim {'CONFIRMED' if confirmed else 'NOT confirmed'}"
              f" at this cell")
        confirmed_cells += int(confirmed)

    print("\n" + "=" * 72)
    if confirmed_cells == len(cells):
        print("VERDICT: CONFIRMED at all cells -- the K+-restricted Lipschitz")
        print("seminorm is identically zero on the entire operator system.")
        print("thm:main's LHS is degenerate in the cited framework.")
    else:
        print(f"VERDICT: confirmed at {confirmed_cells}/{len(cells)} cells -- "
              "adversarial finding NOT fully reproduced; investigate.")
    print("=" * 72)


if __name__ == "__main__":
    main()
