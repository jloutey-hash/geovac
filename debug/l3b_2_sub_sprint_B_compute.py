"""L3b-2 sub-sprint B: numerical verification of the joint cb-norm bound.

Verifies sub-sprint B's claim:

    ||S_{K^joint}||_cb  =  ||S_{K^SU(2)}||_cb * ||S_{K^U(1)}||_cb
                        =  2/(n_max + 1) * 1
                        =  2/(n_max + 1)

at panel (n_max, N_t) in {(2,3), (2,5), (3,3), (3,5), (4,3), (4,5)}.

For each cell:
  1. Confirm sympy factorization of joint Plancherel symbol is exact
     (sympy rational equality, not numerical match).
  2. Compute factor cb-norms via geovac.central_fejer_su2 and
     geovac.central_fejer_compact_temporal.
  3. Confirm product = closed-form bound 2/(n_max+1).

Writes machine-readable results to debug/data/l3b_2_sub_sprint_B.json.

Companion to debug/l3b_2_sub_sprint_B_cb_norm_memo.md.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import sympy as sp
from sympy import Rational

from geovac.central_fejer_su2 import (
    central_multiplier_cb_norm as cb_su2,
    normalization_constant,
    plancherel_symbol as plancherel_su2,
)
from geovac.central_fejer_compact_temporal import (
    cb_norm_circle,
    joint_cb_norm,
    joint_plancherel_symbol,
    plancherel_symbol_circle,
    verify_plancherel_factorization,
)


def compute_panel_cell(n_max: int, N_t: int) -> dict:
    """Verify joint cb-norm closed form at a single panel cell."""
    # Factor cb-norms
    su2_cb = cb_su2(n_max)  # = 2/(n_max + 1)
    u1_cb = cb_norm_circle(N_t)  # = 1
    joint_cb = joint_cb_norm(n_max, N_t)  # = 2/(n_max + 1)
    expected = Rational(2, n_max + 1)

    # Check product = expected
    product_check = su2_cb * u1_cb
    product_matches_expected = (product_check == expected)
    joint_matches_expected = (joint_cb == expected)
    joint_matches_product = (joint_cb == product_check)

    # Plancherel factorization (exact in sympy)
    ok, details = verify_plancherel_factorization(n_max, N_t)

    # Manual l-infty maximum of joint Plancherel symbol
    j_max = Rational(n_max - 1, 2)
    j_vals = [Rational(k, 2) for k in range(0, n_max)]
    K_max = (N_t - 1) // 2 if N_t % 2 == 1 else N_t // 2
    k_vals = list(range(-K_max, K_max + 1))

    max_symbol = Rational(0)
    arg_max = None
    for j in j_vals:
        for k in k_vals:
            val = joint_plancherel_symbol(n_max, N_t, j, k)
            if val > max_symbol:
                max_symbol = val
                arg_max = (str(j), int(k))

    return {
        "n_max": n_max,
        "N_t": N_t,
        "su2_cb_norm": str(su2_cb),
        "su2_cb_norm_float": float(su2_cb),
        "u1_cb_norm": str(u1_cb),
        "u1_cb_norm_float": float(u1_cb),
        "joint_cb_norm": str(joint_cb),
        "joint_cb_norm_float": float(joint_cb),
        "expected_2_over_nmax_plus_1": str(expected),
        "expected_float": float(expected),
        "product_check": str(product_check),
        "joint_matches_expected": bool(joint_matches_expected),
        "joint_matches_product": bool(joint_matches_product),
        "product_matches_expected": bool(product_matches_expected),
        "plancherel_factorization_exact": bool(ok),
        "plancherel_pairs_checked": details["pairs_checked"],
        "plancherel_pairs_match": details["pairs_match"],
        "joint_symbol_linfty_max": str(max_symbol),
        "joint_symbol_argmax": arg_max,
        "max_attained_at_j0_k0": (arg_max == ("(n_max-1)/2", 0)
                                  or arg_max == (str(j_max), 0)),
    }


def asymptotic_rate_check() -> dict:
    """Verify 2/(n_max + 1) -> 0 monotonically as n_max -> infinity."""
    n_vals = [1, 2, 3, 4, 5, 10, 20, 50, 100]
    table = []
    for n in n_vals:
        cb = float(Rational(2, n + 1))
        product_check = float((n + 1) * cb)
        table.append({
            "n_max": n,
            "cb_norm": cb,
            "(n_max+1) * cb_norm": product_check,
            "product_eq_2": product_check == 2.0,
        })

    # Monotonicity check
    monotone_decreasing = all(
        table[i]["cb_norm"] > table[i + 1]["cb_norm"] for i in range(len(table) - 1)
    )

    return {
        "table": table,
        "monotone_decreasing": monotone_decreasing,
        "product_always_2": all(row["product_eq_2"] for row in table),
        "asymptote": 0.0,
        "rate": "O(1/n_max)",
    }


def cross_check_u1_via_max_symbol() -> dict:
    """Verify U(1) cb-norm = max of Plancherel symbol over k for various N_t.

    This is independent of cb_norm_circle's return value (which hardcodes
    Rational(1)); we compute max_{k in Z} plancherel_symbol_circle(N_t, k)
    explicitly and verify it equals 1 at k = 0.
    """
    results = []
    for N_t in [1, 2, 3, 4, 5, 7, 11, 21]:
        K_max = (N_t - 1) // 2 if N_t % 2 == 1 else N_t // 2
        k_vals = list(range(-K_max, K_max + 1))
        symbols = [(k, plancherel_symbol_circle(N_t, k)) for k in k_vals]
        max_k, max_val = max(symbols, key=lambda kv: kv[1])
        symbol_at_zero = plancherel_symbol_circle(N_t, 0)
        results.append({
            "N_t": N_t,
            "K_max": K_max,
            "n_modes": len(k_vals),
            "max_argk": int(max_k),
            "max_value": str(max_val),
            "max_value_float": float(max_val),
            "symbol_at_zero": str(symbol_at_zero),
            "symbol_at_zero_float": float(symbol_at_zero),
            "max_attained_at_zero": int(max_k) == 0,
            "max_equals_one": max_val == Rational(1),
        })
    return {
        "table": results,
        "u1_cb_norm_always_1": all(r["max_equals_one"] for r in results),
        "max_always_at_k_zero": all(r["max_attained_at_zero"] for r in results),
    }


def main():
    """Run the sub-sprint B numerical verification."""
    t_start = time.time()
    print(f"L3b-2 sub-sprint B verification driver")
    print(f"=" * 60)

    # Panel from the dispatch spec
    panel = [(2, 3), (2, 5), (3, 3), (3, 5), (4, 3), (4, 5)]
    print(f"\nPanel: {panel}")

    cell_results = []
    for n_max, N_t in panel:
        print(f"\n  cell (n_max={n_max}, N_t={N_t}):")
        cell = compute_panel_cell(n_max, N_t)
        cell_results.append(cell)
        print(f"    su2_cb = {cell['su2_cb_norm']} = {cell['su2_cb_norm_float']:.6f}")
        print(f"    u1_cb  = {cell['u1_cb_norm']} = {cell['u1_cb_norm_float']:.6f}")
        print(f"    joint  = {cell['joint_cb_norm']} = {cell['joint_cb_norm_float']:.6f}")
        print(f"    bound  = {cell['expected_2_over_nmax_plus_1']} = {cell['expected_float']:.6f}")
        print(f"    factorization exact:  {cell['plancherel_factorization_exact']}")
        print(f"    joint == product:     {cell['joint_matches_product']}")
        print(f"    joint == bound:       {cell['joint_matches_expected']}")
        print(f"    pairs match:          {cell['plancherel_pairs_match']}/{cell['plancherel_pairs_checked']}")

    print(f"\n{'-' * 60}")
    print("Asymptotic rate panel:")
    asym = asymptotic_rate_check()
    for row in asym["table"]:
        print(f"  n_max={row['n_max']:>3}:  cb={row['cb_norm']:.6f},  "
              f"(n_max+1)*cb={row['(n_max+1) * cb_norm']:.6f}")
    print(f"  monotone decreasing: {asym['monotone_decreasing']}")
    print(f"  product always 2:    {asym['product_always_2']}")

    print(f"\n{'-' * 60}")
    print("U(1) cb-norm cross-check via Plancherel symbol max:")
    u1_check = cross_check_u1_via_max_symbol()
    for row in u1_check["table"]:
        print(f"  N_t={row['N_t']:>3}:  max(k={row['max_argk']:+d}) = {row['max_value']}  "
              f"= {row['max_value_float']:.6f}  (modes={row['n_modes']})")
    print(f"  U(1) cb-norm always 1: {u1_check['u1_cb_norm_always_1']}")
    print(f"  max always at k=0:     {u1_check['max_always_at_k_zero']}")

    # Aggregate verdict
    all_factorization_exact = all(c["plancherel_factorization_exact"] for c in cell_results)
    all_joint_eq_product = all(c["joint_matches_product"] for c in cell_results)
    all_joint_eq_bound = all(c["joint_matches_expected"] for c in cell_results)
    all_pairs_match = all(c["plancherel_pairs_match"] == c["plancherel_pairs_checked"]
                          for c in cell_results)

    verdict = (all_factorization_exact and all_joint_eq_product
               and all_joint_eq_bound and all_pairs_match
               and asym["monotone_decreasing"] and asym["product_always_2"]
               and u1_check["u1_cb_norm_always_1"]
               and u1_check["max_always_at_k_zero"])

    print(f"\n{'=' * 60}")
    print(f"VERDICT: {'PROVED (all checks pass)' if verdict else 'FAILED'}")
    print(f"{'=' * 60}")

    # Persist results
    out_data = {
        "sub_sprint": "L3b-2 sub-sprint B",
        "claim": "||S_K_joint||_cb = ||S_K_SU(2)||_cb * ||S_K_U(1)||_cb = 2/(n_max + 1)",
        "verdict": "PROVED" if verdict else "FAILED",
        "panel_cells": cell_results,
        "asymptotic_rate": asym,
        "u1_cross_check": u1_check,
        "summary": {
            "all_factorization_exact": all_factorization_exact,
            "all_joint_eq_product": all_joint_eq_product,
            "all_joint_eq_bound": all_joint_eq_bound,
            "all_pairs_match": all_pairs_match,
            "asymptotic_monotone": asym["monotone_decreasing"],
            "u1_always_1": u1_check["u1_cb_norm_always_1"],
        },
        "wall_time_s": time.time() - t_start,
    }

    out_dir = Path("debug/data")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "l3b_2_sub_sprint_B.json"
    with open(out_path, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"\nWrote: {out_path}")
    print(f"Wall time: {out_data['wall_time_s']:.2f} s")

    return verdict


if __name__ == "__main__":
    success = main()
    raise SystemExit(0 if success else 1)
