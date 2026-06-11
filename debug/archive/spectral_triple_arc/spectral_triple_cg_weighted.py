"""
CG-weighted spectral triple comparison: uniform vs CG adjacency.

Compares the FockSpectralTriple at n_max=2 and n_max=3 with uniform
(kappa = -1/16) and CG-weighted (Wigner 3j) edge weights.  Main
diagnostic: does CG weighting improve the Kramers JD = -DJ relation?

Mathematical background:
  JD = -DJ requires BOTH:
    (a) JLambda = -LambdaJ  (diagonal part anticommutes with J)
    (b) J*(kappa*A) = -(kappa*A)*J  (adjacency part anticommutes)

  But JLambda = +LambdaJ always holds (Lambda is chi*|lam|, J preserves
  eigenvalue magnitude, so JLambdaJ^{-1} = +Lambda).  Therefore full
  anticommutation JD = -DJ is impossible unless Lambda = 0.

  The diagnostic therefore characterizes the residual JD + DJ and
  tests whether CG weights improve the off-diagonal part.

Output: JSON report at debug/data/spectral_triple_cg_weighted.json
"""

import json
import sys
import os
import time

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sympy import Rational
from geovac.spectral_triple import FockSpectralTriple


def analyze_spectral_triple(n_max: int, adjacency_weights: str, j_type: str) -> dict:
    """Build spectral triple and collect diagnostics."""
    t0 = time.time()
    st = FockSpectralTriple(
        n_max=n_max,
        j_type=j_type,
        adjacency_weights=adjacency_weights,
    )
    build_time = time.time() - t0

    result = {
        "n_max": n_max,
        "dim_H": st.dim_H,
        "n_sectors": st.n_sectors,
        "j_type": j_type,
        "adjacency_weights": adjacency_weights,
        "build_time_s": round(build_time, 3),
    }

    # Basic axiom checks
    result["D_selfadjoint"] = st.check_selfadjoint()
    result["gamma_squared_I"] = st.check_grading_square()

    j2_ok, j2_eps = st.check_J_squared()
    result["J_squared_ok"] = j2_ok
    result["J_squared_epsilon"] = j2_eps

    jd_ok, jd_msg = st.check_J_D_relation()
    result["J_D_relation_ok"] = jd_ok
    result["J_D_relation_msg"] = jd_msg

    result["pi_free"] = st.verify_pi_free()

    # Kramers residual analysis
    kr = st.check_kramers_D_relation()
    result["kramers_residual"] = {
        "exact_zero": kr["exact_zero"],
        "n_nonzero": kr["n_nonzero"],
        "frobenius_norm": kr["frobenius_norm"],
        "max_entry": str(kr["max_entry"]),
    }

    kra = st.kramers_D_residual_analysis()
    result["kramers_decomposition"] = {
        "diagonal_residual": kra["diagonal_residual"],
        "offdiag_residual": kra["offdiag_residual"],
        "constraint_test": kra["constraint_test"],
        "per_sector_max": [
            {"sector": str(s), "max_val": v}
            for s, v in kra["per_sector_max"]
        ],
    }

    # CG weight matrix stats (if applicable)
    if adjacency_weights == "cg" and st._cg_weights is not None:
        W = st._cg_weights
        N = st.dim_H
        n_nz = sum(1 for i in range(N) for j in range(N) if W[i, j] != 0)
        result["cg_weight_stats"] = {
            "n_nonzero_entries": n_nz,
            "n_edges": n_nz // 2,
        }

    # Spectrum summary
    spec = st.dirac_spectrum_sorted()
    result["spectrum"] = [
        {"eigenvalue": str(ev), "multiplicity": mult}
        for ev, mult in spec
    ]

    return result


def main():
    print("=" * 70)
    print("CG-weighted Spectral Triple Comparison")
    print("=" * 70)

    results = {}

    for n_max in [2, 3]:
        print(f"\n{'=' * 50}")
        print(f"  n_max = {n_max}")
        print(f"{'=' * 50}")

        for j_type in ["permutation", "kramers"]:
            for weights in ["uniform", "cg"]:
                key = f"n{n_max}_{j_type}_{weights}"
                print(f"\n  Building: n_max={n_max}, j_type={j_type}, "
                      f"adjacency_weights={weights} ...")

                r = analyze_spectral_triple(n_max, weights, j_type)
                results[key] = r

                print(f"    dim_H = {r['dim_H']}, build = {r['build_time_s']}s")
                print(f"    D self-adjoint: {r['D_selfadjoint']}")
                print(f"    J^2 epsilon: {r['J_squared_epsilon']}")
                print(f"    JD relation: {r['J_D_relation_msg']}")
                print(f"    pi-free: {r['pi_free']}")
                kr = r["kramers_residual"]
                print(f"    Kramers residual: exact_zero={kr['exact_zero']}, "
                      f"||R||_F={kr['frobenius_norm']:.6f}, "
                      f"n_nz={kr['n_nonzero']}")
                kd = r["kramers_decomposition"]
                print(f"      Diagonal part: ||JLam+LamJ||_F="
                      f"{kd['diagonal_residual']['frobenius_norm']:.6f}")
                print(f"      Off-diag part: ||JA+AJ||_F="
                      f"{kd['offdiag_residual']['frobenius_norm']:.6f}")
                print(f"      Constraint {'{'}J,A{'}'}+2LamJ: ||R||_F="
                      f"{kd['constraint_test']['frobenius_norm']:.6f}")

    # Comparison summary
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)

    for n_max in [2, 3]:
        print(f"\n  n_max = {n_max}:")
        for j_type in ["kramers"]:
            uni_key = f"n{n_max}_{j_type}_uniform"
            cg_key = f"n{n_max}_{j_type}_cg"
            uni = results[uni_key]["kramers_residual"]
            cg = results[cg_key]["kramers_residual"]
            print(f"    Kramers ||JD+DJ||_F: uniform={uni['frobenius_norm']:.6f}, "
                  f"cg={cg['frobenius_norm']:.6f}")
            if uni['frobenius_norm'] > 0:
                ratio = cg['frobenius_norm'] / uni['frobenius_norm']
                print(f"    CG/uniform ratio: {ratio:.4f}")

            uni_d = results[uni_key]["kramers_decomposition"]
            cg_d = results[cg_key]["kramers_decomposition"]
            print(f"    Off-diag ||JA+AJ||_F: uniform="
                  f"{uni_d['offdiag_residual']['frobenius_norm']:.6f}, "
                  f"cg={cg_d['offdiag_residual']['frobenius_norm']:.6f}")
            print(f"    Constraint ||{{J,A}}+2LamJ||_F: uniform="
                  f"{uni_d['constraint_test']['frobenius_norm']:.6f}, "
                  f"cg={cg_d['constraint_test']['frobenius_norm']:.6f}")

    # Save results
    out_path = os.path.join(os.path.dirname(__file__), "data",
                            "spectral_triple_cg_weighted.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
