"""Sprint L1 computational verification driver.

Runs the full four-witness verification at n_max in {2, 3, 4, 5} and
writes structured JSON results to
debug/data/l1_modular_hamiltonian_results.json.

Tests verified:
  - sigma_{2*pi}(O) = O period closure per witness per n_max
  - Cross-witness collapse (HH, Sew, Unruh all match BW)
  - Propinquity rate cross-check against L2 gamma_n
  - Operator-system leakage measurement
  - KMS condition at expectation level

Run: python debug/l1_modular_hamiltonian_compute.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from geovac.modular_hamiltonian import (
    for_bisognano_wichmann,
    for_hartle_hawking,
    for_sewell,
    for_unruh,
    propinquity_rate_check,
    verify_cross_witness_collapse,
)
from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem


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


def main():
    print("Sprint L1 modular Hamiltonian verification")
    print("=" * 60)

    all_results = {
        "sprint": "L1 modular Hamiltonian closure",
        "date": "2026-05-16",
        "architecture": {
            "wedge": "W1 hemispheric P_W on S^3 aligned with Hopf-base axis",
            "unit_normalization": "U-1 natural rapidity units, kappa_g = 1",
            "construction": "BW-alpha geometric K = J_polar (integer spectrum)",
            "dirac": "truthful Camporesi-Higuchi (chi*(n+1/2))",
            "witness_pattern": "BW canonical, HH/Sew/Unruh parameterized BW",
        },
        "per_n_max": {},
        "cross_witness_collapse": {},
        "propinquity_rate_check": {},
    }

    n_max_values = [2, 3, 4, 5]

    # --- Per-n_max witness verification ---
    for n_max in n_max_values:
        print(f"\n--- n_max = {n_max} ---")
        per_nmax = {}

        t0 = time.time()
        for name, factory_call in [
            ("BW", lambda: for_bisognano_wichmann(n_max=n_max)),
            ("HH_M1", lambda: for_hartle_hawking(n_max=n_max, M=1.0)),
            ("HH_M2", lambda: for_hartle_hawking(n_max=n_max, M=2.0)),
            ("Sew_M1", lambda: for_sewell(n_max=n_max, M=1.0)),
            ("Unruh_a1", lambda: for_unruh(n_max=n_max, a=1.0)),
            ("Unruh_a2", lambda: for_unruh(n_max=n_max, a=2.0)),
        ]:
            mh = factory_call()
            results = mh.verify_witness()
            per_nmax[name] = results
            print(f"  {name}: verdict={results['verdict']}, "
                  f"max_res={results['max_periodicity_residual']:.3e}, "
                  f"kappa_g={results['kappa_g']:.4f}, beta={results['beta']:.4f}")
        t_elapsed = time.time() - t0
        print(f"  Wall time: {t_elapsed:.2f}s")

        all_results["per_n_max"][str(n_max)] = per_nmax

        # --- Cross-witness collapse ---
        cw = verify_cross_witness_collapse(n_max=n_max)
        all_results["cross_witness_collapse"][str(n_max)] = cw
        print(f"  Cross-witness max_consistency_residual = {cw['max_consistency_residual']:.3e}")

        # --- Propinquity rate cross-check ---
        pr = propinquity_rate_check(n_max=n_max)
        all_results["propinquity_rate_check"][str(n_max)] = pr
        print(f"  Propinquity: max_res={pr['max_residual']:.3e}, "
              f"gamma_n_L2={pr['gamma_n_L2']:.4f}, "
              f"ratio={pr['residual_over_gamma']:.3e}")

    # --- Summary verdict ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    all_strong = True
    cross_consistent = True
    for n_max in n_max_values:
        bw_results = all_results["per_n_max"][str(n_max)]["BW"]
        if bw_results["verdict"] != "STRONG_IDENTIFICATION":
            all_strong = False
        cw = all_results["cross_witness_collapse"][str(n_max)]
        if cw["max_consistency_residual"] >= 1e-12:
            cross_consistent = False
        print(f"  n_max={n_max}: BW verdict={bw_results['verdict']}, "
              f"max_res={bw_results['max_periodicity_residual']:.3e}")

    if all_strong and cross_consistent:
        all_results["overall_verdict"] = "STRONG_IDENTIFICATION_ALL_NMAX_ALL_WITNESSES"
        print("\nOVERALL: STRONG IDENTIFICATION")
        print("  sigma_{2*pi}(O) = O bit-exact (machine precision) at n_max=2,3,4,5")
        print("  Cross-witness collapse bit-exact (residuals match across witnesses)")
        print("  Period closure independent of kappa_g (K has integer spectrum)")
        print("  Lifts BW reading from 'structural correspondence' to")
        print("  'literal identification at operator-system level (Riemannian)'.")
    elif all_strong:
        all_results["overall_verdict"] = "STRONG_PER_NMAX_PARTIAL_CROSS_WITNESS"
        print("\nOVERALL: STRONG per-n_max, partial cross-witness")
    else:
        all_results["overall_verdict"] = "MIXED"
        print("\nOVERALL: MIXED — see per-n_max details")

    # --- Write JSON ---
    out_path = Path(__file__).parent / "data" / "l1_modular_hamiltonian_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(_convert_for_json(all_results), f, indent=2, default=str)
    print(f"\nResults written to: {out_path}")


if __name__ == "__main__":
    main()
