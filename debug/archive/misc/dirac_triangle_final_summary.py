"""
Final summary verification: both (INT) and (SI) on all four panels.

Outputs counts and writes to JSON.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product, weights_of_irrep,
)


def check_pair(la, lam, lam_pr):
    """Check (INT) and (SI) for a pair."""
    rank = la.rank
    rho = la.rho_omega
    rho_sq = la.rho_sq
    lam_pr_dual = la.dual(lam_pr)
    decomp = tensor_product(la, lam, lam_pr_dual)
    weights_dual = weights_of_irrep(la, lam_pr_dual)

    lam_diff = tuple(lam[i] - lam_pr[i] for i in range(rank))
    lam_diff_sq = la._inner_omega(lam_diff, lam_diff)
    lam_diff_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    lam_diff_plus_rho_sq = la._inner_omega(lam_diff_plus_rho, lam_diff_plus_rho)

    sigmas_data = []
    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho = tuple(sigma[i] + rho[i] for i in range(rank))
        sigma_plus_rho_sq = la._inner_omega(sigma_plus_rho, sigma_plus_rho)
        c_sigma = sigma_plus_rho_sq - rho_sq

        int_holds = (lam_diff_sq <= c_sigma)
        si_holds = (lam_diff_plus_rho_sq <= sigma_plus_rho_sq)

        # Lemma 1 check
        mu_0 = tuple(sigma[i] - lam[i] for i in range(rank))
        lemma1_holds = (mu_0 in weights_dual)

        # Q in Q+
        Q = tuple(sigma[i] + lam_pr[i] - lam[i] for i in range(rank))
        # Convert to alpha-basis to check Q+ membership
        # In omega basis, simple roots are columns of Cartan; check Q can be written as nonneg combo
        # For simplicity, check |Q|^2 <= 2 <lambda', Q>
        Q_norm_sq = la._inner_omega(Q, Q)
        Q_dot_lam_pr = la._inner_omega(Q, lam_pr)
        lemma1_norm_holds = (Q_norm_sq <= 2 * Q_dot_lam_pr)

        sigmas_data.append({
            "sigma": list(sigma),
            "mult": mult,
            "int_holds": int_holds,
            "si_holds": si_holds,
            "lemma1_holds": lemma1_holds,
            "lemma1_norm_holds": lemma1_norm_holds,
        })

    return sigmas_data


def main():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "dirac_triangle_proof_verification.json")

    panels = [
        ("SU(3) p+q<=5", build_A(2), panel_dominant_weights(2, 5)),
        ("SU(4) a+b+c<=3", build_A(3), panel_dominant_weights(3, 3)),
        ("Sp(2)=C_2 a+b<=3", build_C2(), panel_dominant_weights(2, 3)),
        ("G_2 a+b<=2", build_G2(), panel_dominant_weights(2, 2)),
    ]

    results = {
        "description": "Final verification of (INT), (SI), Lemma 1, and Lemma 1 norm bound.",
        "summary": [],
    }

    print(f"{'Panel':<25} {'#sigmas':>9} {'INT':>6} {'SI':>6} {'L1':>6} {'L1norm':>8}")
    print("-" * 70)

    for label, la, panel in panels:
        total = 0
        int_pass = 0
        si_pass = 0
        l1_pass = 0
        l1norm_pass = 0

        for lam in panel:
            for lam_pr in panel:
                sigmas_data = check_pair(la, lam, lam_pr)
                for s in sigmas_data:
                    total += 1
                    if s["int_holds"]:
                        int_pass += 1
                    if s["si_holds"]:
                        si_pass += 1
                    if s["lemma1_holds"]:
                        l1_pass += 1
                    if s["lemma1_norm_holds"]:
                        l1norm_pass += 1

        print(f"{label:<25} {total:>9} {int_pass:>6} {si_pass:>6} {l1_pass:>6} {l1norm_pass:>8}")
        results["summary"].append({
            "label": label,
            "total_sigmas": total,
            "int_pass": int_pass,
            "si_pass": si_pass,
            "lemma1_pass": l1_pass,
            "lemma1_norm_pass": l1norm_pass,
        })

    total = sum(p["total_sigmas"] for p in results["summary"])
    int_pass = sum(p["int_pass"] for p in results["summary"])
    si_pass = sum(p["si_pass"] for p in results["summary"])
    l1_pass = sum(p["lemma1_pass"] for p in results["summary"])
    l1n_pass = sum(p["lemma1_norm_pass"] for p in results["summary"])
    print("-" * 70)
    print(f"{'TOTAL':<25} {total:>9} {int_pass:>6} {si_pass:>6} {l1_pass:>6} {l1n_pass:>8}")
    results["overall"] = {
        "total": total, "int_pass": int_pass, "si_pass": si_pass,
        "lemma1_pass": l1_pass, "lemma1_norm_pass": l1n_pass,
    }

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
