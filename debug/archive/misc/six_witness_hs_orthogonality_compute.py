"""Six-witness HS-orthogonality probe.

Sprint: Lorentzian-extension dictionary-completion (post-h_local_residual_pslq).
Date:   2026-05-17.

Goal
----

A prior sprint (debug/h_local_residual_pslq_memo.md) found that on the BW
canonical Lorentzian Krein wedge,

    <H_local, D_W^L>_HS  =  Tr(H_local^dagger D_W^L)  =  0      (bit-exact)

so the H_local-vs-D_W^L residual decomposes Pythagoreanly as

    r^2  =  ||H_local||_F^2 + ||D_W^L||_F^2.

The L2-E memo (Sprint L2-E synthesis) established that the SIX witnesses
(BW, HH at M=1, HH at M=2, Sewell at M=1, Unruh at a=1, Unruh at a=2) are
bit-identical at the algebra-action level because rho_W^L is
beta-independent under the BW choice H_local := K_L_alpha^W / beta.  But
the HS inner product <H_local, D_W^L> has been computed only ONCE (BW
canonical).

This sprint tests whether HS-orthogonality is UNIVERSAL across the six
witnesses.  If it is, the orthogonality promotes from "L2-E coincidence"
to "theorem of the truncated spectral triple under any wedge-KMS modular
construction."

Method
------

For each witness (BW, HH_M1, HH_M2, Sew_M1, Unruh_a1, Unruh_a2) at
n_max in {1, 2, 3} at N_t = 1:

  1.  Build the LorentzianModularHamiltonian via the existing factory.
  2.  Extract H_local := K_L_alpha^W / beta on the wedge.
  3.  Build the wedge-restricted Lorentzian Dirac D_W^L via
      restrict_operator_to_wedge_block(lorentzian_dirac_matrix(krein), ...).
  4.  Compute <H_local, D_W^L>_HS = Tr(H_local^dagger D_W^L) and the
      Frobenius norms ||H_local||_F, ||D_W^L||_F, ||H_local - D_W^L||_F.
  5.  Check Pythagoras: r^2 == ||H_local||^2 + ||D_W^L||^2 (bit-exact iff
      <H_local, D_W^L>_HS = 0).

Plus one spot-check at (n_max, N_t) = (3, 11) on BW to confirm
N_t > 1 persistence.

Bit-exactness discipline
------------------------

Machine precision in float64 inner products of D x D matrices scales as
~ eps * D^2 in the worst case.  We adopt the threshold

    bit_exact_zero(D)  :=  eps_machine * D^2  (approximately)

with eps_machine ~ 2.2e-16.  At D = dim_W_L for our panel:

  n_max = 1: dim_W_L = 2,  threshold ~ 2e-15
  n_max = 2: dim_W_L = 8,  threshold ~ 1.4e-14
  n_max = 3: dim_W_L = 20, threshold ~ 9e-14

Any residual at or below this scale is structural zero; anything above
is structural non-zero.

Output
------

  - debug/data/six_witness_hs_orthogonality.json
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from geovac.krein_space_construction import KreinSpace
from geovac.lorentzian_dirac import lorentzian_dirac_matrix
from geovac.modular_hamiltonian_lorentzian import (
    LorentzianModularHamiltonian,
    for_bisognano_wichmann_lorentzian,
    for_hartle_hawking_lorentzian,
    for_sewell_lorentzian,
    for_unruh_lorentzian,
    restrict_operator_to_wedge_block,
)


# ---------------------------------------------------------------------------
# Core HS-inner-product computation
# ---------------------------------------------------------------------------


def compute_hs_orthogonality_for_witness(
    lmh: LorentzianModularHamiltonian, witness_label: str,
) -> dict:
    """Compute <H_local, D_W^L>_HS and Pythagoras consistency for one witness.

    Parameters
    ----------
    lmh : LorentzianModularHamiltonian
        Already-built modular Hamiltonian for the witness (kappa_g already
        encoded in the factory).
    witness_label : str
        Human-readable witness identifier.

    Returns
    -------
    dict with HS inner product, Frobenius norms, Pythagoras consistency.
    """
    # H_local := K_L_alpha^W / beta on the wedge (dim_W_L x dim_W_L)
    H_local_wedge = lmh.K_L_alpha_W / lmh.beta

    # D_L^W: wedge-restricted Lorentzian Dirac
    D_L_full = lorentzian_dirac_matrix(lmh.krein)
    D_L_W = restrict_operator_to_wedge_block(
        D_L_full, lmh.krein, lmh.wedge,
    )

    dim_W_L = H_local_wedge.shape[0]
    assert D_L_W.shape == H_local_wedge.shape, (
        f"shape mismatch: H_local {H_local_wedge.shape} "
        f"vs D_L_W {D_L_W.shape}"
    )

    # HS inner product <H_local, D_L^W>_HS = Tr(H_local^dagger D_L^W)
    inner = np.trace(H_local_wedge.conj().T @ D_L_W)
    inner_re = float(np.real(inner))
    inner_im = float(np.imag(inner))
    inner_abs = float(np.abs(inner))

    # Frobenius norms (squared)
    norm_H_sq = float(np.linalg.norm(H_local_wedge, "fro") ** 2)
    norm_D_sq = float(np.linalg.norm(D_L_W, "fro") ** 2)

    # r = ||H_local - D_W^L||_F  (computed directly)
    r_direct = float(np.linalg.norm(H_local_wedge - D_L_W, "fro"))
    r_sq_direct = r_direct ** 2

    # Pythagoras predicted (assuming inner == 0):
    #   ||H - D||^2 = ||H||^2 + ||D||^2 - 2 Re<H, D>
    # So bit-exact Pythagoras r^2 = ||H||^2 + ||D||^2  iff  Re<H, D> = 0.
    r_sq_pythagoras = norm_H_sq + norm_D_sq
    pyth_residual = float(np.abs(r_sq_direct - r_sq_pythagoras))

    # Threshold for "bit-exact zero" of the HS inner product:
    # machine precision scaled by problem size
    eps_machine = np.finfo(np.float64).eps
    bit_zero_threshold = eps_machine * dim_W_L ** 2

    # Verdict
    is_orthogonal = inner_abs <= max(bit_zero_threshold, 1e-13)
    pyth_consistent = pyth_residual <= max(
        bit_zero_threshold * max(r_sq_direct, 1.0), 1e-12
    )

    return {
        "witness": witness_label,
        "n_max": lmh.krein.n_max,
        "N_t": lmh.krein.N_t,
        "kappa_g": float(lmh.kappa_g),
        "beta": float(lmh.beta),
        "dim_W_L": dim_W_L,
        "inner_re": inner_re,
        "inner_im": inner_im,
        "inner_abs": inner_abs,
        "norm_H_local_sq": norm_H_sq,
        "norm_D_L_W_sq": norm_D_sq,
        "r_sq_direct": r_sq_direct,
        "r_sq_pythagoras": r_sq_pythagoras,
        "pythagoras_residual": pyth_residual,
        "bit_zero_threshold": bit_zero_threshold,
        "is_orthogonal_HS": bool(is_orthogonal),
        "pythagoras_consistent": bool(pyth_consistent),
    }


def make_witness(label: str, n_max: int, N_t: int = 1) -> LorentzianModularHamiltonian:
    """Build a LorentzianModularHamiltonian for the named witness."""
    if label == "BW":
        return for_bisognano_wichmann_lorentzian(n_max, N_t=N_t, T_max=1.0)
    if label == "HH_M1":
        return for_hartle_hawking_lorentzian(n_max, N_t=N_t, T_max=1.0, M=1.0)
    if label == "HH_M2":
        return for_hartle_hawking_lorentzian(n_max, N_t=N_t, T_max=1.0, M=2.0)
    if label == "Sew_M1":
        return for_sewell_lorentzian(n_max, N_t=N_t, T_max=1.0, M=1.0)
    if label == "Unruh_a1":
        return for_unruh_lorentzian(n_max, N_t=N_t, T_max=1.0, a=1.0)
    if label == "Unruh_a2":
        return for_unruh_lorentzian(n_max, N_t=N_t, T_max=1.0, a=2.0)
    raise ValueError(f"unknown witness label: {label}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


WITNESSES = ["BW", "HH_M1", "HH_M2", "Sew_M1", "Unruh_a1", "Unruh_a2"]
N_MAX_PANEL = [1, 2, 3]


def run_panel(N_t: int = 1) -> list:
    """Compute six-witness HS-orthogonality at every n_max in the panel."""
    rows = []
    for n_max in N_MAX_PANEL:
        for w in WITNESSES:
            t0 = time.time()
            lmh = make_witness(w, n_max, N_t=N_t)
            row = compute_hs_orthogonality_for_witness(lmh, w)
            row["elapsed_sec"] = float(time.time() - t0)
            rows.append(row)
            print(
                f"  [{w:9s}] n_max={n_max} N_t={N_t}  "
                f"|<H,D>| = {row['inner_abs']:.4e}  "
                f"orthogonal = {row['is_orthogonal_HS']}  "
                f"pyth_res = {row['pythagoras_residual']:.4e}"
            )
    return rows


def main():
    print("=" * 78)
    print("SIX-WITNESS HS-ORTHOGONALITY PROBE")
    print("=" * 78)
    print()
    print("Main panel: 6 witnesses x 3 n_max at N_t = 1 (Riemannian limit)")
    print()
    main_panel = run_panel(N_t=1)

    print()
    print("Spot check: BW at (n_max, N_t) = (3, 11)  (N_t > 1 persistence)")
    print()
    lmh_spot = for_bisognano_wichmann_lorentzian(3, N_t=11, T_max=1.0)
    spot = compute_hs_orthogonality_for_witness(lmh_spot, "BW")
    spot["elapsed_sec"] = 0.0
    print(
        f"  [BW] n_max=3 N_t=11  |<H,D>| = {spot['inner_abs']:.4e}  "
        f"orthogonal = {spot['is_orthogonal_HS']}  "
        f"pyth_res = {spot['pythagoras_residual']:.4e}"
    )

    # Aggregate verdict
    all_orthogonal = all(r["is_orthogonal_HS"] for r in main_panel)
    all_pythagoras = all(r["pythagoras_consistent"] for r in main_panel)
    max_inner = max(r["inner_abs"] for r in main_panel)
    max_pyth = max(r["pythagoras_residual"] for r in main_panel)

    print()
    print("=" * 78)
    print("AGGREGATE VERDICT")
    print("=" * 78)
    print(f"  All 18 HS inner products bit-exact zero: {all_orthogonal}")
    print(f"  All 18 Pythagoras consistent:            {all_pythagoras}")
    print(f"  max |<H,D>|_HS over panel:               {max_inner:.4e}")
    print(f"  max Pythagoras residual over panel:      {max_pyth:.4e}")
    print(f"  Spot (3, 11) BW orthogonal:              "
          f"{spot['is_orthogonal_HS']}")

    out = {
        "sprint": "six_witness_hs_orthogonality",
        "date": "2026-05-17",
        "main_panel_N_t_1": main_panel,
        "spot_check_3_11_BW": spot,
        "summary": {
            "all_18_orthogonal": all_orthogonal,
            "all_18_pythagoras_consistent": all_pythagoras,
            "max_inner_abs": max_inner,
            "max_pythagoras_residual": max_pyth,
            "spot_check_orthogonal": spot["is_orthogonal_HS"],
            "spot_check_inner_abs": spot["inner_abs"],
        },
    }

    out_path = Path("debug/data/six_witness_hs_orthogonality.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
