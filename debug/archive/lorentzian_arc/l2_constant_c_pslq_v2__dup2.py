"""PSLQ v2: rerun PSLQ at the higher-precision c from the precision_push sprint.

Reads c_best from l2_constant_c_precision_push.json and runs the same PSLQ
campaigns at coefficient ceilings 10^4, 10^5, 10^6 against the same seven test
bases as l2_constant_c_identification.py.

This is the definitive PSLQ run; the original sprint had only ~14 digits of
effective precision in c.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, pi

# Reuse all basis builders from the original script
from l2_constant_c_identification import (  # noqa: E402
    build_basis_basic,
    build_basis_dirichlet_L,
    build_basis_multi_hurwitz,
    build_basis_stein_weiss,
    build_basis_wildcard,
    deduplicate_basis,
    run_pslq,
)


DPS_PSLQ = 250
COEFFICIENT_CEILINGS = [10**4, 10**5, 10**6]


def load_c_high_precision():
    """Load c at highest available precision.

    Uses MR-C 80-digit string value directly. The actual *digit precision*
    of c is bounded by the truncation order of the asymptotic fit; both
    MR-C K=5 vs K=6 disagree at digit 12, and the precision_push sprint
    K=6 disagrees with MR-C K=6 at digit 13. So c is reliable to ~12
    digits, with the 80-digit string carrying the panel precision but
    not the genuine analytical precision.

    PSLQ at ceiling 10^M against B-element basis needs ~M*B digits to
    avoid spurious relations; conversely a NULL result at higher digits
    of mismatch is a strong negative.
    """
    # MUST set mp.dps before mpf(<long string>) to preserve precision.
    mp.dps = 250
    p_v1 = PROJ / "debug" / "data" / "mr_c_l2_subleading.json"
    if p_v1.exists():
        with open(p_v1) as f:
            d1 = json.load(f)
        return mpf(d1["c_best"]), int(d1.get("best_K", 6)), len(d1.get("panel", []))
    p_v2 = PROJ / "debug" / "data" / "l2_constant_c_precision_push.json"
    if p_v2.exists():
        with open(p_v2) as f:
            d = json.load(f)
        return mpf(d["c_best_250dps"]), int(d["best_K"]), int(d["best_n_points"])
    raise FileNotFoundError("No c data found - run l2_constant_c_*.py scripts first")


def main():
    print("=" * 72)
    print("L2 constant c: PSLQ campaign v2 (high-precision c)")
    print("=" * 72)
    c, best_K, best_n_pts = load_c_high_precision()
    mp.dps = DPS_PSLQ
    print(f"c (50 dps) : {mpmath.nstr(c, 50)}")
    print(f"  best_K   : {best_K}")
    print(f"  panel pts: {best_n_pts}")
    print(f"  ceilings : {COEFFICIENT_CEILINGS}")
    print(f"  dps_PSLQ : {DPS_PSLQ}")
    print()

    # Build bases
    basis_basic = build_basis_basic(DPS_PSLQ)
    basis_dirichletL = build_basis_dirichlet_L(DPS_PSLQ)
    basis_multihurwitz = build_basis_multi_hurwitz(DPS_PSLQ)
    basis_steinweiss = build_basis_stein_weiss(DPS_PSLQ)
    basis_wildcard = build_basis_wildcard(DPS_PSLQ)

    test_bases = [
        ("MR-C_baseline", basis_basic),
        ("DirichletL_addon", deduplicate_basis(basis_basic + basis_dirichletL)),
        ("MultiHurwitz_addon", deduplicate_basis(basis_basic + basis_multihurwitz)),
        ("SteinWeissIBP_addon", deduplicate_basis(basis_basic + basis_steinweiss)),
        ("Wildcard_addon", deduplicate_basis(basis_basic + basis_wildcard)),
        ("DirichletL_StrictMin", deduplicate_basis(
            [("1", mpf(1)), ("pi", pi), ("4/pi", 4 / pi)] + basis_dirichletL)),
        ("Union_All", deduplicate_basis(
            basis_basic + basis_dirichletL + basis_multihurwitz +
            basis_steinweiss + basis_wildcard)),
    ]

    out = {
        "c_50dps": mpmath.nstr(c, 50),
        "best_K": best_K,
        "best_n_points": best_n_pts,
        "dps_PSLQ": DPS_PSLQ,
        "coefficient_ceilings": COEFFICIENT_CEILINGS,
        "results_by_basis": {},
    }

    for basis_name, basis in test_bases:
        print(f"--- Basis: {basis_name} ({len(basis)} elements) ---")
        results = []
        for M in COEFFICIENT_CEILINGS:
            r = run_pslq(c, basis, max_coeff=M, dps=DPS_PSLQ)
            results.append(r)
            if r.get("found"):
                cf = " + ".join(f"({q})*{lbl}" for q, lbl in r["closed_form_terms"])
                tag = "PASS_80" if r.get("verify_pass_1e_minus_80") else (
                      "PASS_50" if r.get("verify_pass_1e_minus_50") else (
                      "PASS_30" if r.get("verify_pass_1e_minus_30") else "FAIL"))
                print(f"  [M={M}, {tag}] c = {cf}")
                print(f"    abs err = {r.get('verify_abs_error')}")
            else:
                note = r.get("note", "no relation")
                print(f"  [M={M}] {note}")
        out["results_by_basis"][basis_name] = results
        print()

    out_path = PROJ / "debug" / "data" / "l2_constant_c_pslq_v2.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"-> {out_path}")
    print()

    # Summary
    print("=" * 72)
    print("SUMMARY (PSLQ v2)")
    print("=" * 72)
    any_pass = False
    for basis_name, results in out["results_by_basis"].items():
        passes_80 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_80")]
        passes_50 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_50")]
        passes_30 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_30")]
        if passes_80:
            best = passes_80[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"])
            print(f"  {basis_name:24s} PASS_80 at M={best['max_coeff']}: c = {cf}")
            any_pass = True
        elif passes_50:
            best = passes_50[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"])
            print(f"  {basis_name:24s} PASS_50 at M={best['max_coeff']}: c = {cf}")
            any_pass = True
        elif passes_30:
            best = passes_30[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"])
            print(f"  {basis_name:24s} PASS_30 at M={best['max_coeff']}: c = {cf}")
        else:
            largest_M = max(r["max_coeff"] for r in results)
            print(f"  {basis_name:24s} NULL at M up to {largest_M}")
    print()
    if any_pass:
        print("VERDICT: identification(s) found above.")
    else:
        print(f"VERDICT: NULL across all {len(test_bases)} test bases at all ceilings up to {max(COEFFICIENT_CEILINGS)}.")


if __name__ == "__main__":
    main()
