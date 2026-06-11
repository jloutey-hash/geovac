"""Quick (2, 3) panel compute for Sprint L3b-2f-alpha — minimal scoping.

Drives the structural diagnostic + Lambda estimate at the (n_max=2, N_t=3)
cell only.  Computes flip vs natural commutator content at full precision
but skips the expensive propagation-number computation.

Output:
  debug/data/l3b_2f_alpha_enlarged_substrate.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

_LOG = Path(__file__).parent / "data" / "l3b_2f_alpha_progress.log"
_LOG.parent.mkdir(parents=True, exist_ok=True)
_LOG.write_text("", encoding="utf-8")


def L(msg):
    print(msg, flush=True)
    with open(_LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n"); f.flush()


L(f"[start] {time.strftime('%H:%M:%S')} importing geovac modules ...")

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)

L(f"[done] {time.strftime('%H:%M:%S')} imports done")


def op_norm(A):
    if A.size == 0: return 0.0
    return float(np.linalg.svd(A, compute_uv=False)[0])

def fro(A): return float(np.linalg.norm(A))


def run_panel(n_max, N_t, T, paper45_lambda):
    L(f"\n[panel] (n_max={n_max}, N_t={N_t}, T={T:.4f})  {time.strftime('%H:%M:%S')}")
    t0 = time.time()

    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    dim_K = K.dim
    P_plus, P_minus = K.positive_negative_split()
    d_plus, d_minus = K.split_dimensions()
    J = K.J
    L(f"  dim K = {dim_K}, K^+ = {d_plus}, K^- = {d_minus}, "
      f"elapsed = {time.time()-t0:.1f}s")

    D_L = lorentzian_dirac_compact_matrix(K)
    # Decompose: D_L_diag = (1/2)(D_L + J D_L J^-1), D_L_off = (1/2)(D_L - J D_L J^-1)
    Jinv = J  # J^2 = +I
    D_L_diag = 0.5 * (D_L + J @ D_L @ Jinv)
    D_L_off  = 0.5 * (D_L - J @ D_L @ Jinv)
    L(f"  D_L_op  = {op_norm(D_L):.4f}, "
      f"D_L_diag = {op_norm(D_L_diag):.4f}, "
      f"D_L_off  = {op_norm(D_L_off):.4f}, "
      f"sanity [J,D_L_diag]={fro(J@D_L_diag - D_L_diag@J):.1e}, "
      f"{{J,D_L_off}}={fro(J@D_L_off + D_L_off@J):.1e}, "
      f"elapsed = {time.time()-t0:.1f}s")

    O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)
    d_w = O.dim_spatial // 2
    L(f"  natural O^L: dim = {O.dim}, "
      f"dim_Weyl = {d_w}, n_gens = {len(O.basis_matrices)}, "
      f"elapsed = {time.time()-t0:.1f}s")

    # Construct flip generators in Paper 44's chiral basis:
    #   M^spat_flip = diag(W, -W),  {gamma^0, M^spat_flip} = 0.
    flip_gens = []
    for (N, Ls, M), full_nat in zip(O.spat_labels, O._spat_matrices):
        W = full_nat[:d_w, :d_w].copy()
        M_spat_flip = np.zeros_like(full_nat)
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        for p, gp in enumerate(O._temp_matrices):
            flip_gens.append(((int(N), int(Ls), int(M), int(p)),
                               np.kron(M_spat_flip, gp)))
    L(f"  built {len(flip_gens)} flip generators, "
      f"elapsed = {time.time()-t0:.1f}s")

    # Verify anti-commutation with J at the SPATIAL level.
    # (Tensoring with diagonal temporal matrices preserves anticommutation.)
    max_anti = 0.0; max_comm = 0.0
    for (lbl, Mf) in flip_gens:
        a = fro(J @ Mf + Mf @ J)
        c = fro(J @ Mf - Mf @ J)
        if a > max_anti: max_anti = a
        if c > max_comm: max_comm = c
    L(f"  max ||{{J, M_flip}}|| = {max_anti:.3e} "
      f"(expect ~0 if M_flip = diag(W,-W) anti-commutes), "
      f"max ||[J, M_flip]|| = {max_comm:.3e} "
      f"(expect >0, structurally non-zero)")

    # Compute structural-identity diagnostic for ALL flip generators.
    # Also for ALL natural generators (for baseline).
    flip_diag_norms = []
    flip_off_norms  = []
    flip_L_op_list  = []
    flip_L_block_P45_list = []
    flip_Cpp_list   = []
    flip_Cpm_list   = []
    nat_diag_norms = []
    nat_off_norms  = []
    nat_L_op_list  = []
    nat_L_block_P45_list = []
    nat_Cpp_list   = []

    L(f"  computing structural diagnostics for {len(flip_gens)} flip "
      f"+ {len(O.basis_matrices)} natural generators ...")

    # Cache D_L_proj_plus and D_L_proj_minus for L_block_P45.
    D_L_proj_plus  = P_plus @ D_L @ P_plus
    D_L_proj_minus = P_minus @ D_L @ P_minus

    t_loop = time.time()
    # Natural
    for k, (lbl, M) in enumerate(O.basis_matrices):
        c_DL    = D_L @ M - M @ D_L
        c_diag  = D_L_diag @ M - M @ D_L_diag
        c_off   = D_L_off  @ M - M @ D_L_off
        nat_diag_norms.append(op_norm(c_diag))
        nat_off_norms.append (op_norm(c_off))
        nat_L_op_list.append(op_norm(c_DL))
        # K^+ block restriction commutator.
        Cpp = P_plus @ c_DL @ P_plus
        nat_Cpp_list.append(op_norm(Cpp))
        # Paper 45 K^+-weak-form:
        ap = P_plus @ M @ P_plus
        L_P45 = op_norm(D_L_proj_plus @ ap - ap @ D_L_proj_plus)
        nat_L_block_P45_list.append(L_P45)

    L(f"    natural diagnostics done, elapsed = {time.time()-t_loop:.1f}s")

    t_loop = time.time()
    sample_flip = None
    for k, (lbl, M) in enumerate(flip_gens):
        c_DL    = D_L @ M - M @ D_L
        c_diag  = D_L_diag @ M - M @ D_L_diag
        c_off   = D_L_off  @ M - M @ D_L_off
        flip_diag_norms.append(op_norm(c_diag))
        flip_off_norms.append (op_norm(c_off))
        flip_L_op_list.append(op_norm(c_DL))
        Cpp = P_plus @ c_DL @ P_plus
        Cpm = P_plus @ c_DL @ P_minus
        flip_Cpp_list.append(op_norm(Cpp))
        flip_Cpm_list.append(op_norm(Cpm))
        ap = P_plus @ M @ P_plus
        # NB: for flip generators, P_+ M P_+ should vanish since M flips chirality.
        L_P45 = op_norm(D_L_proj_plus @ ap - ap @ D_L_proj_plus)
        flip_L_block_P45_list.append(L_P45)
        if sample_flip is None:
            sample_flip = {
                "label": list(lbl),
                "L_op": flip_L_op_list[-1],
                "norm_comm_D_L_diag": flip_diag_norms[-1],
                "norm_comm_D_L_off": flip_off_norms[-1],
                "norm_Cpp_block_plus": flip_Cpp_list[-1],
                "norm_Cpm_cross": flip_Cpm_list[-1],
                "L_block_P45_K+_norm": L_P45,
                "P_plus_M_P_plus_fro": fro(ap),
            }
    L(f"    flip diagnostics done, elapsed = {time.time()-t_loop:.1f}s")

    # Summary statistics
    max_nat_diag  = max(nat_diag_norms) if nat_diag_norms else 0.0
    max_nat_off   = max(nat_off_norms)  if nat_off_norms else 0.0
    max_nat_L_op  = max(nat_L_op_list) if nat_L_op_list else 0.0
    max_nat_L_P45 = max(nat_L_block_P45_list) if nat_L_block_P45_list else 0.0
    max_flip_diag = max(flip_diag_norms)
    max_flip_off  = max(flip_off_norms)
    max_flip_L_op = max(flip_L_op_list)
    max_flip_L_P45 = max(flip_L_block_P45_list)
    L(f"  STRUCTURAL SUMMARY at (n_max={n_max}, N_t={N_t}):")
    L(f"    max ||[D_L_diag, a_nat]||    = {max_nat_diag:.4e}  (expect ~0 by L3b-2a §3.3)")
    L(f"    max ||[D_L_diag, a_flip]||   = {max_flip_diag:.4e}  (expect >0 -- diagnostic key)")
    L(f"    max ||[D_L_off,  a_nat]||    = {max_nat_off:.4e}")
    L(f"    max ||[D_L_off,  a_flip]||   = {max_flip_off:.4e}")
    L(f"    max L_op(a_nat)              = {max_nat_L_op:.4e}")
    L(f"    max L_op(a_flip)             = {max_flip_L_op:.4e}")
    L(f"    max L_block_P45(a_nat)       = {max_nat_L_P45:.4e}")
    L(f"    max L_block_P45(a_flip)      = {max_flip_L_P45:.4e}  (expect ~0 if P_+ M_flip P_+ = 0)")

    # Verify chirality-projection on flip generators.
    sample_proj_fro = sample_flip["P_plus_M_P_plus_fro"] if sample_flip else 0.0
    L(f"    sample P_+ M_flip P_+ Frob   = {sample_proj_fro:.4e} "
      f"(expect ~0 since M_flip is off-block-diagonal in K^+/K^-)")

    # Lambda^enlarged estimate.
    # For natural: L_op = ||D_L_off, a|| (diag piece vanishes).
    # For flip: L_op has BOTH diag and off content; conservative
    # upper bound L_op ≤ sqrt(diag^2 + off^2).
    max_flip_L_op_bound = max(
        float(np.sqrt(d*d + o*o)) for d, o in zip(flip_diag_norms, flip_off_norms)
    )
    ratio = max_flip_L_op_bound / max(max_nat_L_op, 1e-30)
    lambda_enlarged = paper45_lambda * max(1.0, ratio)
    L(f"  Lambda ESTIMATE (Paper 45 = {paper45_lambda}):")
    L(f"    ratio max_flip_L_op_bound / max_nat_L_op = {ratio:.4f}")
    L(f"    Lambda^enlarged estimate                  = {lambda_enlarged:.4f}")

    elapsed = time.time() - t0
    L(f"  panel done in {elapsed:.1f}s")

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": dim_K,
        "dim_Kplus": d_plus,
        "dim_Kminus": d_minus,
        "dim_Weyl": d_w,
        "dim_natural": O.dim,
        "num_natural_generators": len(O.basis_matrices),
        "num_flip_generators": len(flip_gens),
        "achievable_env_natural": d_w * d_w * N_t,
        "achievable_env_enlarged": 2 * d_w * d_w * N_t,
        "full_envelope_dim": dim_K * dim_K,
        "D_L_decomposition": {
            "D_L_op_norm": op_norm(D_L),
            "D_L_diag_op_norm": op_norm(D_L_diag),
            "D_L_off_op_norm": op_norm(D_L_off),
            "comm_J_D_L_diag_fro": fro(J@D_L_diag - D_L_diag@J),
            "anti_J_D_L_off_fro": fro(J@D_L_off + D_L_off@J),
        },
        "max_anti_J_M_flip_fro": max_anti,
        "max_comm_J_M_flip_fro": max_comm,
        "structural_summary": {
            "max_norm_comm_D_L_diag_natural":  max_nat_diag,
            "max_norm_comm_D_L_diag_flip":     max_flip_diag,
            "max_norm_comm_D_L_off_natural":   max_nat_off,
            "max_norm_comm_D_L_off_flip":      max_flip_off,
            "max_L_op_natural":                max_nat_L_op,
            "max_L_op_flip":                   max_flip_L_op,
            "max_L_op_flip_upper_bound":       max_flip_L_op_bound,
            "max_L_block_P45_natural":         max_nat_L_P45,
            "max_L_block_P45_flip":            max_flip_L_P45,
        },
        "lambda_estimate": {
            "paper45_lambda": paper45_lambda,
            "lambda_enlarged_estimate": lambda_enlarged,
            "ratio_flip_over_natural": ratio,
            "max_natural_L_op": max_nat_L_op,
            "max_flip_L_op_upper_bound": max_flip_L_op_bound,
        },
        "sample_diagnostic": sample_flip,
        "propagation": {
            "prop_achievable": "ANALYTICAL_DEFERRED",
            "prop_full": "ANALYTICAL_DEFERRED",
            "notes": (
                "Propagation number on the enlarged substrate is "
                "analytically expected to remain prop=2 (achievable) "
                "because the enlarged generators are also chirality-doubled "
                "scalar multipliers (up to a sign), preserving the same "
                "spectral-block structure as the natural substrate.  Full "
                "envelope expected to remain prop=infinity (as on the "
                "natural substrate).  Computational verification deferred "
                "to a follow-on sub-sprint."
            ),
        },
        "elapsed_seconds": elapsed,
    }


def main():
    T = 2.0 * np.pi
    p45 = {(2, 3): 2.0746, (3, 5): 1.6101}

    panels = []
    # (3, 5) construction is observed to be very slow due to sympy-arithmetic
    # spinor 3-Y integrals; deferred to follow-on sub-sprint.  (2, 3) is the
    # representative cell for this scoping memo.
    test_cells = [(2, 3)]
    out_file = _LOG.parent / "l3b_2f_alpha_enlarged_substrate.json"

    for (n_max, N_t) in test_cells:
        panel = run_panel(n_max, N_t, T, p45[(n_max, N_t)])
        panels.append(panel)
        # incremental dump
        payload = {"sprint": "L3b-2f-alpha", "panels": panels}
        with open(out_file, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=str)
        L(f"[partial] wrote {out_file}")

    # Verdict
    ratios = [p["lambda_estimate"]["ratio_flip_over_natural"] for p in panels]
    lambdas = [p["lambda_estimate"]["lambda_enlarged_estimate"] for p in panels]
    min_r = min(ratios); max_r = max(ratios)
    if min_r > 1.1:
        verdict = "POSITIVE-GO"
    elif max_r <= 1.01:
        verdict = "NEGATIVE"
    elif min_r > 1.0:
        verdict = "POSITIVE-GO (narrow margin)"
    else:
        verdict = "PARTIAL"

    L("\n" + "=" * 70)
    L("VERDICT")
    L("=" * 70)
    L(f"  verdict           = {verdict}")
    L(f"  ratios            = {ratios}")
    L(f"  lambda_enlarged   = {lambdas}")
    L(f"  min ratio         = {min_r:.4f}")
    L(f"  max ratio         = {max_r:.4f}")

    payload = {
        "sprint": "L3b-2f-alpha",
        "panels": panels,
        "verdict": {
            "verdict": verdict,
            "ratios": ratios,
            "min_ratio": min_r,
            "max_ratio": max_r,
            "lambda_enlarged_estimates": lambdas,
        },
    }
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)
    L(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
