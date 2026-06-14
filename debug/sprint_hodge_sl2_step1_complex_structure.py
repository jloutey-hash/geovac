"""Sprint Hodge-SL2 — Step 1 diagnostic: a D-compatible complex structure.

The Hodge realization of GeoVac's SL_2 needs a weight-1 polarized Hodge
structure, whose complex structure is an operator with square -1. The
natural candidate is the Kramers real structure J (J^2 = -I, the
SU(2) = S^3 quaternionic structure). The codebase flags JD = -DJ for the
Kramers J as requiring CG-weighted adjacency ("future work"). This driver
resolves that open item: build the Kramers J on both the uniform (baseline,
known to fail) and CG-weighted (candidate) Dirac, and test JD = -DJ
exactly in sympy.

Decision gate (per scoping memo step 1):
  PASS     J^2 = -I AND JD + DJ = 0 bit-exact on CG substrate -> KO-dim-3
           real structure confirmed; a D-compatible complex structure
           exists; proceed to polarization (step 2).
  PARTIAL  J^2 = -I but JD + DJ != 0: characterize the residual
           (Frobenius, per-sector, constraint test); Hodge prerequisite
           compromised -> diagnose before proceeding.
  FAIL     structure absent.

Run: python debug/sprint_hodge_sl2_step1_complex_structure.py
"""

from __future__ import annotations

import json
from typing import Dict

import sympy as sp

from geovac.spectral_triple import FockSpectralTriple


def probe(n_max: int, adjacency: str) -> Dict:
    """Build Kramers J on the given adjacency and test J^2 and JD = -DJ."""
    st = FockSpectralTriple(
        n_max=n_max, j_type="kramers", adjacency_weights=adjacency
    )

    j2_ok, eps = st.check_J_squared()
    jd = st.check_kramers_D_relation()           # R = JD + DJ
    decomp = st.kramers_D_residual_analysis()    # diagonal / offdiag / constraint

    return {
        "n_max": n_max,
        "adjacency": adjacency,
        "dim_H": st.dim_H,
        "J2_is_minus_I": bool(j2_ok and eps == -1),
        "J2_epsilon": int(eps),
        "JD_anticommute_exact_zero": bool(jd["exact_zero"]),
        "JD_residual_n_nonzero": int(jd["n_nonzero"]),
        "JD_residual_frobenius": float(jd["frobenius_norm"]),
        "JD_residual_max_entry": str(jd["max_entry"]),
        "diag_residual": jd_safe(decomp["diagonal_residual"]),
        "offdiag_residual": jd_safe(decomp["offdiag_residual"]),
        "constraint_JA_plus_2LamJ_zero": bool(decomp["constraint_test"]["exact_zero"]),
        "constraint_frobenius": float(decomp["constraint_test"]["frobenius_norm"]),
    }


def jd_safe(d: Dict) -> Dict:
    """Coerce a residual sub-dict to JSON-safe primitives."""
    return {
        "exact_zero": bool(d["exact_zero"]),
        "n_nonzero": int(d["n_nonzero"]),
        "frobenius_norm": float(d["frobenius_norm"]),
    }


def verdict(cg2: Dict) -> str:
    """Decision gate on the CG-substrate result at n_max = 2."""
    if cg2["J2_is_minus_I"] and cg2["JD_anticommute_exact_zero"]:
        return "PASS"
    if cg2["J2_is_minus_I"]:
        return "PARTIAL"
    return "FAIL"


def main() -> None:
    results = {}
    # n_max = 2 (16-dim, fast) is the decision substrate; n_max = 3
    # (40-dim) is confirmation.
    for n_max in (2, 3):
        for adjacency in ("uniform", "cg"):
            key = f"nmax{n_max}_{adjacency}"
            print(f"=== {key} ===", flush=True)
            res = probe(n_max, adjacency)
            results[key] = res
            print(
                f"  J^2 = -I: {res['J2_is_minus_I']} (eps={res['J2_epsilon']});  "
                f"JD+DJ exact zero: {res['JD_anticommute_exact_zero']};  "
                f"residual nnz={res['JD_residual_n_nonzero']}, "
                f"Frob={res['JD_residual_frobenius']:.4g}",
                flush=True,
            )
            print(
                f"  diag JLam+LamJ zero: {res['diag_residual']['exact_zero']};  "
                f"offdiag JA+AJ zero: {res['offdiag_residual']['exact_zero']};  "
                f"constraint {{J,A}}+2LamJ zero: "
                f"{res['constraint_JA_plus_2LamJ_zero']}",
                flush=True,
            )

    gate = verdict(results["nmax2_cg"])
    cg3 = results.get("nmax3_cg", {})
    gate3 = (
        cg3.get("J2_is_minus_I") and cg3.get("JD_anticommute_exact_zero")
    )
    summary = {
        "step1_verdict_nmax2_cg": gate,
        "step1_confirm_nmax3_cg": bool(gate3),
        "baseline_uniform_fails_as_expected": (
            not results["nmax2_uniform"]["JD_anticommute_exact_zero"]
        ),
    }
    results["_summary"] = summary

    print("\n=== STEP 1 VERDICT ===", flush=True)
    print(f"  CG n_max=2: {gate}", flush=True)
    print(f"  CG n_max=3 confirms PASS: {bool(gate3)}", flush=True)
    print(
        f"  uniform baseline fails JD=-DJ as expected: "
        f"{summary['baseline_uniform_fails_as_expected']}",
        flush=True,
    )

    out = "debug/data/sprint_hodge_sl2_step1.json"
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
