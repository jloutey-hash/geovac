"""G4a n_max=3 verification — extend Connes SM axiom checks to n_max=3.

Constructs the full Connes Standard Model almost-commutative triple
T = T_GV(n_max=3) x T_F^SM and verifies all load-bearing axioms.

Compares residuals across n_max in {1, 2, 3} to check whether violations
scale toward zero with increasing cutoff.

Output: debug/data/g4a_nmax3_verification.json
"""

import json
import time
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.standard_model_triple import (
    sm_one_generation,
    sm_gauge_only,
)


def run_axiom_verification(n_max: int) -> dict:
    """Run full axiom verification at given n_max, return results dict."""
    print(f"\n{'='*60}")
    print(f"  n_max = {n_max}")
    print(f"{'='*60}")

    t0 = time.time()
    T = sm_one_generation(n_max)
    t_build = time.time() - t0
    print(f"  Build time: {t_build:.2f}s")
    print(f"  dim_GV = {T.dim_GV}, dim_F = {T.dim_F}, dim_H = {T.dim_H}")

    # --- Axiom verification ---
    t0 = time.time()
    axioms = T.verify_axioms()
    t_axioms = time.time() - t0
    print(f"  Axiom verification time: {t_axioms:.2f}s")

    for key, val in axioms.items():
        status = "PASS" if val < 1e-10 else ("WARN" if val < 1.0 else "FAIL")
        print(f"    {key}: {val:.6e}  [{status}]")

    # --- Gauge group census ---
    t0 = time.time()
    if n_max >= 2:
        census = T.gauge_group_census(n_samples=30, seed=99)
        t_gauge = time.time() - t0
        print(f"  Gauge census time: {t_gauge:.2f}s")
        print(f"    Gauge group: {census['gauge_group']}")
        print(f"    SU(2) lepton: {census['SU2_lepton']}")
        print(f"    SU(2) quark:  {census['SU2_quark']}")
        print(f"    SU(3):        {census['SU3']}")
        print(f"    U(1):         {census['U1']}")
        print(f"    Lepton-quark mixing: {census['lepton_quark_mixing']}")
    else:
        census = {"gauge_group": "trivial (n_max=1)", "note": "skipped"}
        t_gauge = 0.0

    # --- Falsifier (Higgs check) ---
    t0 = time.time()
    higgs_zero_yukawa, reason_y, data_y = T.check_natural_negative(
        n_random_generators=30, seed=42
    )
    t_falsifier = time.time() - t0
    print(f"  Falsifier time: {t_falsifier:.2f}s")
    print(f"    With Yukawa: higgs_zero={higgs_zero_yukawa}")
    print(f"    Reason: {reason_y}")

    T0 = sm_gauge_only(n_max)
    higgs_zero_no_yukawa, reason_n, data_n = T0.check_natural_negative(
        n_random_generators=30, seed=42
    )
    print(f"    Without Yukawa: higgs_zero={higgs_zero_no_yukawa}")

    # --- Additional axiom checks not in verify_axioms ---
    # gamma^2, {gamma, D}, J gamma relation
    D = T.dirac_combined()
    gamma = T.chirality_combined()
    U_J = T.real_structure_combined()
    I_H = np.eye(T.dim_H, dtype=np.complex128)

    gamma_sq_residual = float(np.linalg.norm(gamma @ gamma - I_H))
    gamma_D_anticomm = float(np.linalg.norm(gamma @ D + D @ gamma))
    J_gamma = U_J @ np.conj(gamma)
    gamma_J = gamma @ U_J
    j_gamma_residual = float(np.linalg.norm(J_gamma - gamma_J))

    print(f"    gamma^2 = I residual: {gamma_sq_residual:.6e}")
    print(f"    {{gamma, D}} residual: {gamma_D_anticomm:.6e}")
    print(f"    J gamma - gamma J residual: {j_gamma_residual:.6e}")

    return {
        "n_max": n_max,
        "dim_GV": T.dim_GV,
        "dim_F": T.dim_F,
        "dim_H": T.dim_H,
        "build_time_s": round(t_build, 3),
        "axiom_time_s": round(t_axioms, 3),
        "gauge_time_s": round(t_gauge, 3),
        "falsifier_time_s": round(t_falsifier, 3),
        "axioms": {k: float(v) for k, v in axioms.items()},
        "extra_axioms": {
            "gamma_squared_residual": gamma_sq_residual,
            "gamma_D_anticommutator_residual": gamma_D_anticomm,
            "J_gamma_relation_residual": j_gamma_residual,
        },
        "gauge_census": {
            k: (bool(v) if isinstance(v, (bool, np.bool_)) else v)
            for k, v in census.items()
        },
        "falsifier_with_yukawa": {
            "higgs_zero": bool(higgs_zero_yukawa),
            "higgs_max": data_y["higgs_max"],
            "gauge_max": data_y["gauge_max"],
        },
        "falsifier_no_yukawa": {
            "higgs_zero": bool(higgs_zero_no_yukawa),
            "higgs_max": data_n["higgs_max"],
            "gauge_max": data_n["gauge_max"],
        },
    }


def main():
    print("G4a n_max=3 Connes SM verification")
    print("=" * 60)

    results = {}
    for n_max in [1, 2, 3]:
        results[f"n_max_{n_max}"] = run_axiom_verification(n_max)

    # --- Comparison table ---
    print("\n" + "=" * 60)
    print("  COMPARISON TABLE")
    print("=" * 60)
    print(f"{'Axiom':<40} {'n_max=1':>12} {'n_max=2':>12} {'n_max=3':>12}")
    print("-" * 76)

    axiom_keys = list(results["n_max_1"]["axioms"].keys())
    for key in axiom_keys:
        v1 = results["n_max_1"]["axioms"][key]
        v2 = results["n_max_2"]["axioms"][key]
        v3 = results["n_max_3"]["axioms"][key]
        print(f"{key:<40} {v1:>12.3e} {v2:>12.3e} {v3:>12.3e}")

    extra_keys = list(results["n_max_1"]["extra_axioms"].keys())
    for key in extra_keys:
        v1 = results["n_max_1"]["extra_axioms"][key]
        v2 = results["n_max_2"]["extra_axioms"][key]
        v3 = results["n_max_3"]["extra_axioms"][key]
        print(f"{key:<40} {v1:>12.3e} {v2:>12.3e} {v3:>12.3e}")

    print(f"\n{'dim_H':<40} {results['n_max_1']['dim_H']:>12} {results['n_max_2']['dim_H']:>12} {results['n_max_3']['dim_H']:>12}")

    # --- Save results ---
    outpath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data", "g4a_nmax3_verification.json"
    )
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {outpath}")


if __name__ == "__main__":
    main()
