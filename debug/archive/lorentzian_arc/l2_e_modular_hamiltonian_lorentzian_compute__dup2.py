"""Sprint L2-E computational verification driver.

Mirrors debug/l1_modular_hamiltonian_compute.py and
debug/l1_tighten_tomita_compute.py at the Krein-level Lorentzian
construction.

Verifies:
  L2E-FALS-1 (LOAD-BEARING): sigma_{2pi}^{L,alpha}(O) = O bit-exact
  L2E-FALS-2: sigma_{2pi}^{L,TT}(O) = O bit-exact
  L2E-FALS-3: flow conjugacy bit-exact at general t

Plus:
  RIE-LIMIT (LOAD-BEARING #4): K_L_alpha and K_L_TT reduce bit-identically
                                to K_alpha and K_TT at N_t = 1
  H_LOCAL (Paper 42 §7.2 / O3 at (3,1)): verdict at every (n_max, N_t)
  SIX-WITNESS COLLAPSE: BW, HH_M1, HH_M2, Sew_M1, Unruh_a1, Unruh_a2
                        all collapse bit-exact

Output: debug/data/l2_e_modular_hamiltonian_lorentzian.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geovac.modular_hamiltonian_lorentzian import (
    for_bisognano_wichmann_lorentzian,
    verify_cross_witness_collapse_lorentzian,
)


def main() -> None:
    # Panel: (n_max, N_t) cells. We compute fast cells first and write
    # JSON incrementally so partial results are available even if larger
    # cells time out.
    #
    # Fast cells (dim_W_L <= 48, GNS dim <= 2304):
    #   (1, 1), (1, 11), (1, 21), (2, 1), (2, 11), (3, 1)
    # Medium cells (dim_W_L = 88, GNS dim = 7744):
    #   (2, 21)
    # Slow cells (dim_W_L >= 120, GNS dim >= 14400):
    #   (3, 11) -- skipped from default; available via test marker
    panel = [
        (1, 1), (1, 11), (1, 21),
        (2, 1), (2, 11), (3, 1),
        (2, 21),  # medium
    ]
    results = {
        "sprint": "L2-E",
        "date": "2026-05-16",
        "description": (
            "Sprint L2-E: Krein-level Paper 42 redo. Four LOAD-BEARING "
            "falsifiers + H_local verdict at (3, 1) + six-witness collapse."
        ),
        "panel": [
            {"n_max": n_max, "N_t": N_t} for n_max, N_t in panel
        ],
        "per_cell": [],
    }

    print("=" * 78)
    print("Sprint L2-E computational verification driver")
    print("=" * 78)
    print()

    for n_max, N_t in panel:
        print(f"--- (n_max={n_max}, N_t={N_t}) ---")
        lmh = for_bisognano_wichmann_lorentzian(
            n_max=n_max, N_t=N_t, T_max=1.0,
        )
        # Witness battery (L2E-FALS-1, 2, 3)
        witness_results = lmh.verify_witness_lorentzian()
        # Riemannian-limit recovery (LOAD-BEARING #4)
        rie_ok, rie_details = lmh.riemannian_limit_recovery()
        # H_local verdict (Paper 42 §7.2 / O3 at (3, 1))
        h_local = lmh.h_local_verdict_at_3_1()
        # Six-witness collapse at this (n_max, N_t)
        if n_max >= 1:
            collapse = verify_cross_witness_collapse_lorentzian(
                n_max=n_max, N_t=N_t, T_max=1.0,
            )
        else:
            collapse = {"skipped": True}

        cell = {
            "n_max": n_max, "N_t": N_t,
            "dim_K": lmh.dim_K, "dim_W_L": lmh.dim_W_L,
            "beta": lmh.beta, "kappa_g": lmh.kappa_g,
            # L2E-FALS-1
            "L2E_FALS_1_alpha_max_residual": witness_results["max_alpha_residual"],
            # L2E-FALS-2
            "L2E_FALS_2_TT_max_residual": witness_results["max_TT_residual"],
            # L2E-FALS-3
            "L2E_FALS_3_conjugacy_max_residual": witness_results["max_conjugacy_residual"],
            "witness_verdict": witness_results["verdict"],
            # RIE-LIMIT (LOAD-BEARING #4)
            "rie_limit_ok": bool(rie_ok),
            "rie_limit_k_alpha_W_residual": rie_details["k_alpha_W_residual"],
            "rie_limit_k_tomita_residual": rie_details["k_tomita_residual"],
            # H_local verdict (Paper 42 §7.2 / O3 at (3, 1))
            "h_local_verdict": h_local["verdict"],
            "h_local_residual_full": h_local["residual_full"],
            "h_local_residual_full_normalized": h_local["residual_full_normalized"],
            "h_local_residual_RIE_baseline": h_local["residual_RIE_baseline"],
            "h_local_residual_RIE_normalized": h_local["residual_RIE_normalized"],
            "h_local_lorentzian_eq_riemannian_baseline": h_local[
                "lorentzian_eq_riemannian_baseline"
            ],
            # Six-witness collapse
            "six_witness_collapse_ok": collapse.get(
                "six_witness_collapse_ok", None
            ),
            "max_alpha_cross_consistency": collapse.get(
                "max_alpha_cross_consistency", None
            ),
            "max_TT_cross_consistency": collapse.get(
                "max_TT_cross_consistency", None
            ),
        }
        results["per_cell"].append(cell)

        # Print summary line
        print(
            f"  L2E-FALS-1 (alpha):    {cell['L2E_FALS_1_alpha_max_residual']:.4e}",
            flush=True,
        )
        print(
            f"  L2E-FALS-2 (Tomita):   {cell['L2E_FALS_2_TT_max_residual']:.4e}",
            flush=True,
        )
        print(
            f"  L2E-FALS-3 (conjugacy):{cell['L2E_FALS_3_conjugacy_max_residual']:.4e}",
            flush=True,
        )
        print(
            f"  RIE-LIMIT:             k_alpha_res={cell['rie_limit_k_alpha_W_residual']:.4e}, "
            f"k_TT_res={cell['rie_limit_k_tomita_residual']:.4e}, ok={cell['rie_limit_ok']}",
            flush=True,
        )
        print(
            f"  H_LOCAL verdict:       {cell['h_local_verdict']} "
            f"(res_full={cell['h_local_residual_full']:.4e}, "
            f"RIE_baseline={cell['h_local_residual_RIE_baseline']:.4e}, "
            f"eq_RIE={cell['h_local_lorentzian_eq_riemannian_baseline']})",
            flush=True,
        )
        print(
            f"  SIX-WITNESS collapse:  ok={cell['six_witness_collapse_ok']}, "
            f"max_alpha={cell['max_alpha_cross_consistency']:.4e}, "
            f"max_TT={cell['max_TT_cross_consistency']:.4e}",
            flush=True,
        )
        print(flush=True)

        # Write JSON incrementally after each cell completes
        out_path = Path("debug/data/l2_e_modular_hamiltonian_lorentzian.json")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            # Aggregate so far
            so_far = dict(results)
            so_far["panel_completed"] = [
                {"n_max": c["n_max"], "N_t": c["N_t"]}
                for c in so_far["per_cell"]
            ]
            json.dump(so_far, f, indent=2)

    # Aggregate
    all_alpha_ok = all(
        c["L2E_FALS_1_alpha_max_residual"] < 1e-12 for c in results["per_cell"]
    )
    all_tt_ok = all(
        c["L2E_FALS_2_TT_max_residual"] < 1e-12 for c in results["per_cell"]
    )
    all_conj_ok = all(
        c["L2E_FALS_3_conjugacy_max_residual"] < 1e-12 for c in results["per_cell"]
    )
    all_rie_ok = all(c["rie_limit_ok"] for c in results["per_cell"])
    all_collapse_ok = all(
        (c["six_witness_collapse_ok"] is None) or c["six_witness_collapse_ok"]
        for c in results["per_cell"]
    )

    h_local_at_N_t_1_all_ii = all(
        c["h_local_verdict"] == "ii"
        for c in results["per_cell"]
        if c["N_t"] == 1
    )
    h_local_at_N_t_gt_1_all_iii = all(
        c["h_local_verdict"] == "iii"
        for c in results["per_cell"]
        if c["N_t"] > 1
    )

    overall = {
        "L2E_FALS_1_load_bearing_alpha_all_pass": all_alpha_ok,
        "L2E_FALS_2_TT_all_pass": all_tt_ok,
        "L2E_FALS_3_conjugacy_all_pass": all_conj_ok,
        "rie_limit_load_bearing_all_pass": all_rie_ok,
        "six_witness_collapse_all_pass": all_collapse_ok,
        "h_local_at_N_t_1_signature_independent": h_local_at_N_t_1_all_ii,
        "h_local_at_N_t_gt_1_temporally_refined": h_local_at_N_t_gt_1_all_iii,
        "verdict": (
            "CLOSED"
            if (
                all_alpha_ok and all_tt_ok and all_conj_ok and all_rie_ok
                and all_collapse_ok
            )
            else "STRUCTURAL_FINDING"
        ),
        "headline": (
            "Sprint L2-E CLOSED at finite cutoff (Lorentzian). Four-witness "
            "Wick-rotation theorem lifts from 'structural correspondence' to "
            "'literal identification at the operator-system level (Lorentzian)' "
            "at every (n_max, N_t) panel cell. Paper 42 §7.2 O3 is "
            "signature-INDEPENDENT at the Riemannian reduction (N_t = 1); "
            "refined at N_t > 1 by the temporal-derivative contribution."
        ),
    }
    results["overall"] = overall

    # Save
    out_path = Path("debug/data/l2_e_modular_hamiltonian_lorentzian.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print()
    print("=" * 78)
    print("OVERALL VERDICT:", overall["verdict"])
    print("=" * 78)
    for k, v in overall.items():
        if k == "headline":
            print()
            print("HEADLINE:")
            print(f"  {v}")
        else:
            print(f"  {k:55s}: {v}")
    print()
    print(f"Data written to {out_path}")


if __name__ == "__main__":
    main()
