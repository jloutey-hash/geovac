"""Diagnostic: probe a specific generator that shows L_op != L_block.

The chirality-doubling result (K+ block = K- block bit-exact) holds,
but L_block != L_op in general. Why? Because [D_L, a] has off-block-
diagonal content (chirality-flipping) even though a commutes with J.
The off-block-diagonal piece comes from the i*D_GV (x) I part of D_L
which anti-commutes with gamma^0 = J.

This diagnostic verifies the mathematical structure:
- D_L = D_L^diag + D_L^off, where D_L^diag = i*gamma^0 (x) d_t
  (block-diag in K+/K-) and D_L^off = i*D_GV (x) I (off-block-diag).
- For a commuting with J: [D_L, a] = [D_L^diag, a] + [D_L^off, a].
- The diag piece [D_L^diag, a] is block-diag; off piece is off-block-diag.
- L_op picks up both; L_block picks up only the diag piece (or its
  block-restricted norm = full diag-block-restricted norm).
"""

from __future__ import annotations

import numpy as np

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    fourier_d_dt_matrix,
    lorentzian_dirac_compact_matrix,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)
from geovac.full_dirac_operator_system import camporesi_higuchi_full_dirac_matrix


def main() -> None:
    n_max, N_t = 2, 3
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t)
    P_plus, P_minus = K.positive_negative_split()

    # Build the two pieces of D_L explicitly.
    D_GV = camporesi_higuchi_full_dirac_matrix(K.basis_spatial)
    D_t = fourier_d_dt_matrix(N_t, K.T)
    I_t = np.eye(N_t, dtype=np.complex128)

    D_L_diag = 1j * np.kron(K.J_spatial, D_t)
    D_L_off = 1j * np.kron(D_GV, I_t)
    D_L_total = D_L_diag + D_L_off

    # Sanity check: matches lorentzian_dirac_compact_matrix.
    D_L_ref = lorentzian_dirac_compact_matrix(K)
    decomp_residual = np.linalg.norm(D_L_total - D_L_ref)
    print(f"D_L decomposition residual: {decomp_residual:.3e}")

    # Verify the block structure of each piece.
    # D_L_diag should commute with J (block-diag in K+/K-).
    comm_diag_J = np.linalg.norm(K.J @ D_L_diag - D_L_diag @ K.J)
    print(f"||[J, D_L_diag]|| = {comm_diag_J:.3e} (should be 0)")

    # D_L_off should anti-commute with J (off-block-diag).
    anti_off_J = np.linalg.norm(K.J @ D_L_off + D_L_off @ K.J)
    print(f"||{{J, D_L_off}}|| = {anti_off_J:.3e} (should be 0)")

    # Pick a generator with L_op > L_block.
    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)

    # Find the first generator where L_op - L_block is large.
    for (label, M) in O.basis_matrices:
        comm = D_L_total @ M - M @ D_L_total
        L_op = float(np.linalg.svd(comm, compute_uv=False)[0])
        comm_pp = P_plus @ comm @ P_plus
        comm_mm = P_minus @ comm @ P_minus
        L_pp = float(np.linalg.svd(comm_pp, compute_uv=False)[0])
        L_mm = float(np.linalg.svd(comm_mm, compute_uv=False)[0])
        L_block = max(L_pp, L_mm)

        if L_op - L_block > 1e-3:
            print(f"\nGenerator {label}: L_op={L_op:.6f}, "
                  f"L_block={L_block:.6f}, L_pp={L_pp:.6f}, "
                  f"L_mm={L_mm:.6f}")

            # Decompose the commutator.
            comm_from_diag = D_L_diag @ M - M @ D_L_diag
            comm_from_off = D_L_off @ M - M @ D_L_off

            # Verify decomposition.
            comm_decomp_res = np.linalg.norm(
                (comm_from_diag + comm_from_off) - comm
            )
            print(f"  [D_L, M] = [D_L_diag,M] + [D_L_off,M] residual: "
                  f"{comm_decomp_res:.3e}")

            L_diag_op = float(np.linalg.svd(comm_from_diag,
                                            compute_uv=False)[0])
            L_off_op = float(np.linalg.svd(comm_from_off,
                                           compute_uv=False)[0])
            print(f"  ||[D_L_diag, M]||_op = {L_diag_op:.6f} "
                  f"(should equal L_block = {L_block:.6f})")
            print(f"  ||[D_L_off, M]||_op = {L_off_op:.6f}")

            # Block restrictions of the off-diagonal piece.
            comm_off_pm = P_plus @ comm_from_off @ P_minus
            comm_off_mp = P_minus @ comm_from_off @ P_plus
            L_off_cross_pm = float(np.linalg.svd(comm_off_pm,
                                                  compute_uv=False)[0])
            L_off_cross_mp = float(np.linalg.svd(comm_off_mp,
                                                  compute_uv=False)[0])
            print(f"  ||P_+ [D_L_off, M] P_-||_op = "
                  f"{L_off_cross_pm:.6f}")
            print(f"  ||P_- [D_L_off, M] P_+||_op = "
                  f"{L_off_cross_mp:.6f}")

            # Verify [D_L_diag, M] is block-diagonal (preserves K+, K-).
            cross_pm_diag = P_plus @ comm_from_diag @ P_minus
            cross_mp_diag = P_minus @ comm_from_diag @ P_plus
            print(f"  ||P_+ [D_L_diag, M] P_-||_op = "
                  f"{float(np.linalg.svd(cross_pm_diag, compute_uv=False)[0]):.3e} "
                  f"(should be 0)")
            print(f"  ||P_- [D_L_diag, M] P_+||_op = "
                  f"{float(np.linalg.svd(cross_mp_diag, compute_uv=False)[0]):.3e} "
                  f"(should be 0)")

            # Verify [D_L_off, M] is off-block-diagonal (no K+/K+ component).
            diag_pp_off = P_plus @ comm_from_off @ P_plus
            diag_mm_off = P_minus @ comm_from_off @ P_minus
            print(f"  ||P_+ [D_L_off, M] P_+||_op = "
                  f"{float(np.linalg.svd(diag_pp_off, compute_uv=False)[0]):.3e} "
                  f"(should be 0)")
            print(f"  ||P_- [D_L_off, M] P_-||_op = "
                  f"{float(np.linalg.svd(diag_mm_off, compute_uv=False)[0]):.3e} "
                  f"(should be 0)")

            break

    # And a generator where L_block = 0 (so L_op-L_block = L_op).
    print("\n\nSearching for L_block=0 case ...")
    for (label, M) in O.basis_matrices:
        comm = D_L_total @ M - M @ D_L_total
        L_op = float(np.linalg.svd(comm, compute_uv=False)[0])
        comm_pp = P_plus @ comm @ P_plus
        comm_mm = P_minus @ comm @ P_minus
        L_pp = float(np.linalg.svd(comm_pp, compute_uv=False)[0])
        L_block = max(L_pp, float(np.linalg.svd(comm_mm,
                                                 compute_uv=False)[0]))
        if L_block < 1e-12 and L_op > 1e-3:
            print(f"\nGenerator {label}: L_op={L_op:.6f}, "
                  f"L_block={L_block:.3e}")
            # This is a generator where the commutator is purely
            # off-block-diagonal: a commutes with both D_L_diag and
            # with itself, but [D_L_off, a] is purely off-block-diag.
            comm_from_diag = D_L_diag @ M - M @ D_L_diag
            comm_from_off = D_L_off @ M - M @ D_L_off
            L_diag_op = float(np.linalg.svd(comm_from_diag,
                                            compute_uv=False)[0])
            L_off_op = float(np.linalg.svd(comm_from_off,
                                           compute_uv=False)[0])
            print(f"  ||[D_L_diag, M]||_op = {L_diag_op:.3e}")
            print(f"  ||[D_L_off, M]||_op = {L_off_op:.6f}")
            print(f"  Total: L_op = {L_op:.6f}")
            print(f"  Block-diag piece of [D_L,M] = "
                  f"[D_L_diag, M] has norm {L_diag_op:.3e} (the "
                  f"L_block content).")
            break

    # Spatial multipliers like the identity at L=0,M=0: at L>0,M=0,
    # the spatial multiplier might commute with d_t (it's constant
    # in temporal slot) and might or might not commute with D_GV.
    # Let's understand the structure.

    # Concrete: every multiplier M = M^spat (x) M^temp (with M^temp
    # diagonal in momentum). Then
    # [D_L_diag, M] = i*[gamma^0 (x) d_t, M^spat (x) M^temp]
    #               = i*(gamma^0 M^spat (x) d_t M^temp - M^spat gamma^0
    #                    (x) M^temp d_t)
    # Since [gamma^0, M^spat] = 0 (Prop 5.1, M^spat is block-diagonal
    # in chirality), this simplifies to
    #               = i*gamma^0 M^spat (x) [d_t, M^temp]
    # which vanishes iff [d_t, M^temp] = 0, i.e., M^temp is constant.
    # M^temp = diag(omega_k^p): commutes with d_t (both diagonal in
    # momentum). So [D_L_diag, M] = 0 always!
    #
    # Hence the diag commutator is identically zero, and L_block(a) = 0
    # for every generator a where the spatial M^spat commutes with
    # the off-diagonal D_GV. But [D_L_off, M] = i*[D_GV, M^spat] (x) M^temp,
    # which is non-zero whenever M^spat does not commute with D_GV.

    print("\n\n==== STRUCTURAL ANALYSIS ====")
    # Verify [D_L_diag, M] = 0 for every generator.
    max_diag_comm = 0.0
    for (label, M) in O.basis_matrices:
        c = D_L_diag @ M - M @ D_L_diag
        r = float(np.linalg.norm(c))
        if r > max_diag_comm:
            max_diag_comm = r
    print(f"max ||[D_L_diag, M]||_F over all M in O^L: "
          f"{max_diag_comm:.3e}")
    print("This should be very small if [D_L_diag, M] = 0 holds.")

    # And verify [D_L_off, M] is purely off-block-diag for every M.
    max_off_diag_block = 0.0
    for (label, M) in O.basis_matrices:
        c_off = D_L_off @ M - M @ D_L_off
        # Project onto block-diagonal subspace: this should be ~0.
        c_off_pp = P_plus @ c_off @ P_plus
        c_off_mm = P_minus @ c_off @ P_minus
        r = max(float(np.linalg.norm(c_off_pp)),
                float(np.linalg.norm(c_off_mm)))
        if r > max_off_diag_block:
            max_off_diag_block = r
    print(f"max ||P_+ [D_L_off, M] P_+||_F (and P_- ... P_-) over all M: "
          f"{max_off_diag_block:.3e}")
    print("This should be very small: [D_L_off, M] is purely off-block-diag.")


if __name__ == "__main__":
    main()
