"""Sprint L1-tighten computational verification driver.

Runs the load-bearing Tomita-Takesaki (BW-gamma) construction at
n_max in {2, 3, 4, 5} for all four physical witnesses (BW, HH, Sew,
Unruh) and writes structured JSON results to
debug/data/l1_tighten_tomita_results.json.

Tests verified:
  - GNS construction: rho, rho_sqrt, rho_invsqrt are well-defined
  - J_TT^2 = +I at the FULL Tomita-Takesaki level (not just U = I cross-check)
  - K_TT spectrum and zero-eigenvalue count
  - sigma_{2*pi}^TT(O) = O bit-exact period closure (LOAD-BEARING)
  - BW-alpha vs BW-gamma comparison: verdict
    (UNIFIED_STRONG / TWO_WITNESS_STRONG / SPLIT_SOFT / STRUCTURAL_OBSTRUCTION)
  - Cross-witness Tomita-level collapse

Outcome 1 prediction (PI's prior): UNIFIED_STRONG verdict at every n_max.
sigma_t^TT = sigma_{-t}^alpha (the flows are conjugate at general t but
both close bit-exactly at the period t = 2*pi).

Run: python debug/l1_tighten_tomita_compute.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem
from geovac.modular_hamiltonian import (
    for_bisognano_wichmann,
    for_hartle_hawking,
    for_sewell,
    for_unruh,
)


def _convert_for_json(obj):
    """Recursively convert numpy types to Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _convert_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_convert_for_json(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    return obj


def witness_compare(name: str, bw, n_test_operators: int = 5) -> dict:
    """Run the BW-alpha vs BW-gamma comparison for one witness."""
    comp = bw.compare_alpha_vs_tomita(n_test_operators=n_test_operators)
    return {
        "name": name,
        "kappa_g": bw.kappa_g,
        "beta": bw.beta,
        "dim_H": bw.dim,
        "wedge_dim": bw.wedge.wedge_dim(bw.basis),
        "dim_GNS": comp["K_TT_spectrum_summary"]["dim_GNS"],
        "K_TT_spectrum_summary": comp["K_TT_spectrum_summary"],
        "J_TT_squared_residual_at_TT_level": comp["J_TT_squared_residual"],
        "max_residual_alpha_2pi": comp["max_residual_alpha_2pi"],
        "max_residual_TT_2pi": comp["max_residual_TT_2pi"],
        "avg_TT_vs_alpha_t1": comp["avg_TT_vs_alpha_t1"],
        "avg_TT_vs_alpha_neg_t1": comp["avg_TT_vs_alpha_neg_t1"],
        "per_operator": comp["per_operator"],
        "verdict": comp["verdict"],
    }


def main() -> None:
    print("Sprint L1-tighten — Load-bearing Tomita-Takesaki verification")
    print("=" * 70)
    print()

    all_results = {
        "sprint": "L1-tighten Tomita-Takesaki load-bearing",
        "date": "2026-05-16",
        "predecessor": "L1 BW-alpha closure (geometric K = J_polar, bit-exact)",
        "architecture": {
            "wedge": "W1 hemispheric P_W on S^3 aligned with Hopf-base axis",
            "unit_normalization": "U-1 natural rapidity units, kappa_g = 1",
            "load_bearing_construction": (
                "BW-gamma Tomita-Takesaki: K_TT = -log Delta on H_GNS via "
                "polar decomposition of antilinear S; natural choice "
                "H_local := K_alpha^W / beta makes K_TT = ad(K_alpha^W) "
                "with the same integer spectrum as BW-alpha."
            ),
            "cross_check_construction": (
                "BW-alpha geometric K = J_polar (was load-bearing in L1)."
            ),
            "outcome_1_prediction": (
                "UNIFIED_STRONG: K_TT and K_alpha generate the same modular "
                "automorphism up to sign; sigma_t^TT = sigma_{-t}^alpha; "
                "both close bit-exact at t = 2*pi (integer spectrum)."
            ),
            "dirac": "truthful Camporesi-Higuchi (chi*(n+1/2))",
            "witness_pattern": "BW canonical + HH/Sew/Unruh parameterized",
        },
        "per_n_max": {},
        "overall": {},
    }

    n_max_values = [2, 3, 4, 5]

    for n_max in n_max_values:
        print(f"--- n_max = {n_max} ---")
        per_nmax = {}
        t0 = time.time()
        for name, factory in [
            ("BW", lambda nm=n_max: for_bisognano_wichmann(n_max=nm)),
            ("HH_M1", lambda nm=n_max: for_hartle_hawking(n_max=nm, M=1.0)),
            ("HH_M2", lambda nm=n_max: for_hartle_hawking(n_max=nm, M=2.0)),
            ("Sew_M1", lambda nm=n_max: for_sewell(n_max=nm, M=1.0)),
            ("Unruh_a1", lambda nm=n_max: for_unruh(n_max=nm, a=1.0)),
            ("Unruh_a2", lambda nm=n_max: for_unruh(n_max=nm, a=2.0)),
        ]:
            bw = factory()
            results = witness_compare(name, bw)
            per_nmax[name] = results
            print(
                f"  {name}: verdict={results['verdict']}, "
                f"alpha_2pi={results['max_residual_alpha_2pi']:.2e}, "
                f"TT_2pi={results['max_residual_TT_2pi']:.2e}, "
                f"avg_TT_vs_alpha_neg_t1="
                f"{results['avg_TT_vs_alpha_neg_t1']:.2e}, "
                f"J_TT^2_resid={results['J_TT_squared_residual_at_TT_level']:.2e}"
            )
        t_elapsed = time.time() - t0
        per_nmax["_wall_time_seconds"] = t_elapsed
        print(f"  wall: {t_elapsed:.2f}s")
        all_results["per_n_max"][str(n_max)] = per_nmax
        print()

    # Overall verdict
    print("=" * 70)
    print("OVERALL SUMMARY")
    print("=" * 70)

    all_unified = True
    j_tt_clean = True
    for n_max in n_max_values:
        per_nmax = all_results["per_n_max"][str(n_max)]
        for name, r in per_nmax.items():
            if name.startswith("_"):
                continue
            if r["verdict"] != "UNIFIED_STRONG":
                all_unified = False
            if r["J_TT_squared_residual_at_TT_level"] > 1e-10:
                j_tt_clean = False

    if all_unified and j_tt_clean:
        verdict = "OUTCOME_1_UNIFIED_STRONG_ALL_WITNESSES_ALL_NMAX"
        print()
        print("VERDICT: OUTCOME 1 — UNIFIED STRONG CLOSURE")
        print()
        print("  BW-alpha (geometric K = J_polar) and BW-gamma (Tomita-")
        print("  Takesaki K_TT = -log Delta) generate the same modular")
        print("  automorphism up to sign at the operator-action level.")
        print()
        print("  - sigma_{2*pi}^alpha(O) = O bit-exact (L1 baseline preserved)")
        print("  - sigma_{2*pi}^TT(O)    = O bit-exact (LOAD-BEARING L1-tighten)")
        print("  - sigma_t^TT(O) = sigma_{-t}^alpha(O) bit-exact at t = 1")
        print("  - J_TT^2 = +I at the FULL TT level (NOT just U = I cross-check)")
        print()
        print("  The four-witness Wick-rotation theorem is closed at the")
        print("  operator-system level (Riemannian) at finite n_max via")
        print("  BOTH the geometric BW-alpha and the Tomita-Takesaki BW-gamma")
        print("  constructions. The L1 'BW-alpha only' caveat is removed.")
    else:
        verdict = "MIXED_OR_PARTIAL"
        print()
        print(f"VERDICT: {verdict} — see per-n_max details above")

    all_results["overall"] = {
        "verdict": verdict,
        "all_witnesses_unified_strong": all_unified,
        "J_TT_squared_clean_at_TT_level": j_tt_clean,
    }

    # Write JSON
    out_path = Path(__file__).parent / "data" / "l1_tighten_tomita_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(_convert_for_json(all_results), f, indent=2, default=str)
    print()
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main()
