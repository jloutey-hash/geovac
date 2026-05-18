"""Numerical verification driver for L3b-2 Sub-sprint C (joint Berezin).

Verifies the four L4 properties of the joint Berezin reconstruction
map at panel cells (n_max, N_t) in {(2,3), (3,5), (4,7)} for ~5 joint
test functions per cell (mix of separable and non-separable).

Outputs:
  debug/data/l3b_2_sub_sprint_C.json
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

# Add project root to sys.path
PROJ = Path(__file__).resolve().parents[1]
if str(PROJ) not in sys.path:
    sys.path.insert(0, str(PROJ))

from geovac.joint_berezin_compact_temporal import (
    JointBerezinReconstruction,
    JointTestFunction,
    joint_lipschitz_inf_approx,
    joint_panel,
    make_joint_test_function,
    pure_tensor_function,
    temporal_lipschitz_inf,
)
from geovac.r25_l3_lipschitz_bound import (
    TestFunction,
    lipschitz_norm_inf_test_function,
    make_test_function,
)


def evaluate_cell(n_max: int, N_t: int, T: float = 2.0 * np.pi) -> dict:
    """Run the L4 four-property check at one (n_max, N_t) cell."""
    print(f"\n========== (n_max={n_max}, N_t={N_t}) ==========")
    t0 = time.time()

    jb = JointBerezinReconstruction(n_max=n_max, N_t=N_t, T=T)
    print(
        f"  dim_K={jb.dim_K}, dim_spatial={jb.dim_spatial}, "
        f"K_max_u1={jb.plancherel.K_max_u1}"
    )
    print(f"  ||hat_K^joint||_inf = {jb.plancherel.linfty_norm()} = "
          f"{float(jb.plancherel.linfty_norm()):.6f}")

    panel = joint_panel(n_max, N_t)
    print(f"  panel size: {len(panel)}")

    results: List[dict] = []
    for f in panel:
        is_pure_tensor = f.is_pure_tensor()
        try:
            B = jb.apply(f)
            op_norm = float(np.linalg.norm(B, ord=2)) if B.size > 0 else 0.0
        except Exception as e:
            results.append({
                "name": f.name,
                "error": f"apply failed: {e}",
            })
            continue

        # Lipschitz norm via triangle bound
        try:
            lip_inf = joint_lipschitz_inf_approx(f, T=T, metric="L2")
        except Exception as e:
            lip_inf = float("nan")

        # ||f||_infty rough upper bound via triangle on coefficients
        f_inf_upper = sum(abs(c) for _, c in f.coeff_dict.items())

        # (a) Positivity (only meaningful for known-positive functions)
        is_PSD_applicable = (
            "constant" in f.name or "axisymmetric_positive" in f.name
        )
        if is_PSD_applicable:
            psd_ok, min_eig = jb.verify_positivity(f)
        else:
            psd_ok, min_eig = (None, None)

        # (b) Contractivity
        contractive_ok, _, contract_ratio = jb.verify_contractivity(
            f, f_inf_upper
        )

        # (c) Approximate identity
        ai_residual, _ = jb.approximate_identity_residual(f)
        if lip_inf > 1e-12:
            ai_ratio = ai_residual / lip_inf
        else:
            ai_ratio = 0.0

        # (d) L3 compatibility
        try:
            l3_ok, l3_comm_norm, l3_ratio = jb.verify_l3_compatibility(
                f, lip_inf
            )
        except Exception as e:
            l3_ok, l3_comm_norm, l3_ratio = (None, float("nan"), float("nan"))

        # K^+ preservation (structural finding)
        kp_ok, kp_residual = jb.verify_krein_positive_preservation(f)

        # Factor check (for pure tensors)
        factor_ok = None
        factor_residual = None
        if is_pure_tensor and f.coeff_dict:
            # Construct f_s and f_t
            spatial_coeffs = {}
            temporal_coeffs = {}
            for ((N, L, M), q), c in f.coeff_dict.items():
                spatial_coeffs[(N, L, M)] = spatial_coeffs.get((N, L, M), 0)
                temporal_coeffs[q] = temporal_coeffs.get(q, 0)
            # For a rank-1 matrix c_{(N,L,M), q} = a_{NLM} b_q, pick
            # f_s as the first row (q with nonzero contribution)
            # and f_t as the column at that f_s key.
            # Simpler: use the marginal coefficients consistent with
            # rank-1 factorization.  If single term, trivial factorization.
            terms = list(f.coeff_dict.items())
            if len(terms) == 1:
                ((N, L, M), q), c = terms[0]
                f_s = make_test_function(
                    f"factor_s_{N}_{L}_{M}", {(N, L, M): c}
                )
                f_t = {q: 1.0}
                try:
                    factor_ok_b, factor_residual = jb.factor_check(f_s, f_t)
                    factor_ok = factor_ok_b
                except Exception:
                    factor_ok = None
                    factor_residual = None
            else:
                # Multi-term pure tensor: skip factor check (covered by
                # single-term cases)
                factor_ok = None
                factor_residual = None

        result = {
            "name": f.name,
            "is_pure_tensor": bool(is_pure_tensor),
            "is_PSD_applicable": bool(is_PSD_applicable),
            "op_norm": op_norm,
            "f_inf_upper": f_inf_upper,
            "lip_inf_L2": lip_inf,
            # (a) positivity
            "psd_ok": psd_ok,
            "min_eig": min_eig,
            # (b) contractivity
            "contractive_ok": contractive_ok,
            "contract_ratio": contract_ratio,
            # (c) approx-identity
            "approx_identity_residual": ai_residual,
            "approx_identity_ratio_to_lip": ai_ratio,
            # (d) L3 compatibility
            "l3_compat_ok": l3_ok,
            "l3_comm_norm": l3_comm_norm,
            "l3_ratio": l3_ratio,
            # K^+ preservation
            "krein_positive_ok": kp_ok,
            "krein_positive_residual": kp_residual,
            # factor check
            "factor_ok": factor_ok,
            "factor_residual": factor_residual,
        }
        results.append(result)

        psd_str = "n/a" if psd_ok is None else ("PSD" if psd_ok else "FAIL")
        l3_str = "n/a" if l3_ok is None else ("OK" if l3_ok else "FAIL")
        kp_str = "OK" if kp_ok else "FAIL"
        print(
            f"  {f.name[:38]:<38} "
            f"op={op_norm:6.3f}  ai={ai_residual:6.3f}  l3={l3_ratio:6.3f} "
            f"  PSD={psd_str}  L3={l3_str}  K+={kp_str}"
        )

    # Riemannian limit check (only at this n_max)
    if n_max <= 4:
        try:
            jb_1 = JointBerezinReconstruction(n_max=n_max, N_t=1, T=T)
            test_f_s = make_test_function(
                "test_Y3_(2,0,0)", {(2, 0, 0): 1.0}
            )
            rie_ok, rie_details = jb_1.reduce_to_paper38_at_N_t_1(test_f_s)
            rie_residual = rie_details["residual_F_norm"]
        except Exception as e:
            rie_ok = False
            rie_residual = float("nan")
            rie_details = {"error": str(e)}
    else:
        rie_ok = None
        rie_residual = None
        rie_details = None

    cell_summary = {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": jb.dim_K,
        "dim_spatial": jb.dim_spatial,
        "K_max_u1": jb.plancherel.K_max_u1,
        "linfty_norm_joint_K_hat": float(jb.plancherel.linfty_norm()),
        "cb_norm_su2_factor": float(jb.plancherel.weight_su2(n_max)),
        "cb_norm_u1_factor": float(jb.plancherel.weight_u1(0)),
        "panel_size": len(panel),
        "n_results": len(results),
        "results": results,
        "riemannian_limit_at_N_t_1": {
            "ok": rie_ok,
            "residual_F_norm": rie_residual,
            "details": rie_details,
        },
        "n_psd_applicable_pass": sum(
            1 for r in results if r.get("psd_ok") is True
        ),
        "n_contractivity_pass": sum(
            1 for r in results if r.get("contractive_ok") is True
        ),
        "n_l3_compat_pass": sum(
            1 for r in results if r.get("l3_compat_ok") is True
        ),
        "n_krein_positive_pass": sum(
            1 for r in results if r.get("krein_positive_ok") is True
        ),
        "max_approx_identity_residual": max(
            (r.get("approx_identity_residual", 0.0) for r in results),
            default=0.0,
        ),
        "max_l3_ratio": max(
            (r.get("l3_ratio", 0.0) for r in results
             if r.get("l3_ratio") is not None
             and not np.isnan(r.get("l3_ratio", float("nan")))),
            default=0.0,
        ),
        "wall_seconds": time.time() - t0,
    }
    return cell_summary


def main() -> int:
    panel_cells = [(2, 3), (3, 5), (4, 7)]
    T = 2.0 * np.pi

    print("=" * 60)
    print("L3b-2 Sub-sprint C: Joint Berezin reconstruction (L4)")
    print("=" * 60)

    cells: List[dict] = []
    for n_max, N_t in panel_cells:
        cell = evaluate_cell(n_max=n_max, N_t=N_t, T=T)
        cells.append(cell)

    out = {
        "sprint": "L3b-2 Sub-sprint C",
        "date": "2026-05-17",
        "T": T,
        "panel_cells": [{"n_max": n, "N_t": Nt} for n, Nt in panel_cells],
        "cells": cells,
        "summary": {
            "n_cells": len(cells),
            "all_riemannian_limits_pass": all(
                c["riemannian_limit_at_N_t_1"]["ok"] for c in cells
                if c["riemannian_limit_at_N_t_1"]["ok"] is not None
            ),
            "total_panel_entries": sum(c["panel_size"] for c in cells),
            "total_l3_compat_pass": sum(c["n_l3_compat_pass"] for c in cells),
            "total_contractivity_pass": sum(
                c["n_contractivity_pass"] for c in cells
            ),
            "total_krein_positive_pass": sum(
                c["n_krein_positive_pass"] for c in cells
            ),
            "total_psd_applicable_pass": sum(
                c["n_psd_applicable_pass"] for c in cells
            ),
        },
    }

    out_path = PROJ / "debug" / "data" / "l3b_2_sub_sprint_C.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # Final summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for c in cells:
        print(
            f"  (n_max={c['n_max']}, N_t={c['N_t']}): "
            f"panel={c['panel_size']}, "
            f"contractive={c['n_contractivity_pass']}/{c['n_results']}, "
            f"L3-compat={c['n_l3_compat_pass']}/{c['n_results']}, "
            f"K+={c['n_krein_positive_pass']}/{c['n_results']}, "
            f"PSD-applicable={c['n_psd_applicable_pass']}, "
            f"Riemannian limit OK={c['riemannian_limit_at_N_t_1']['ok']}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
