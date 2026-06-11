"""
Refined claim: For every mu in weights of V_{lambda'^*} such that lambda + mu + rho
is REGULAR (no zero coordinate / not on any Weyl wall):
    |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2

This is the right claim because BK contributions from wall-mu's are zero.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, weights_of_irrep,
)


def is_regular(weight):
    """A weight (in omega basis) is regular iff no coordinate is zero."""
    return all(c != 0 for c in weight)


def test_regular_witness(la, lam, lam_pr):
    """
    For all mu in weights of V_{lam_pr^*} with lambda+mu+rho regular:
    check |lam + mu + rho|^2 >= |lam - lam' + rho|^2.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    weights_dual = weights_of_irrep(la, lam_pr_dual)

    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    target_sq = la._inner_omega(lam_minus_lp_plus_rho, lam_minus_lp_plus_rho)

    # Also check if target itself is on wall
    target_regular = is_regular(lam_minus_lp_plus_rho)

    fails = []
    for mu in weights_dual:
        lam_mu_rho = tuple(lam[i] + mu[i] + rho[i] for i in range(rank))
        if not is_regular(lam_mu_rho):
            continue  # Skip wall weights
        norm_sq = la._inner_omega(lam_mu_rho, lam_mu_rho)
        if norm_sq < target_sq:
            fails.append({"mu": mu, "norm_sq": norm_sq, "target_sq": target_sq,
                          "lam_mu_rho": lam_mu_rho})
    return fails, target_regular


def main():
    print("Testing REFINED universal claim:")
    print("  For all mu in weights of V_{lambda'^*} with lambda+mu+rho REGULAR:")
    print("      |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2")
    print()

    panels = [
        ("SU(3) p+q<=4", build_A(2), panel_dominant_weights(2, 4)),
        ("SU(4) a+b+c<=2", build_A(3), panel_dominant_weights(3, 2)),
        ("Sp(2)=C_2 a+b<=2", build_C2(), panel_dominant_weights(2, 2)),
        ("G_2 a+b<=2", build_G2(), panel_dominant_weights(2, 2)),
    ]

    for label, la, panel in panels:
        total_pairs = 0
        total_failures = 0
        first_failures = []
        for lam in panel:
            for lam_pr in panel:
                total_pairs += 1
                fails, target_regular = test_regular_witness(la, lam, lam_pr)
                if fails:
                    total_failures += 1
                    if len(first_failures) < 8:
                        first_failures.append((lam, lam_pr, fails[:3], target_regular))
        print(f"{label}: {total_pairs} pairs, {total_failures} with regular-mu failures")
        for lam, lam_pr, fails, tr in first_failures:
            print(f"  lam={lam}, lam'={lam_pr}, target_regular={tr}:")
            for f in fails:
                print(f"    mu={f['mu']} -> lam+mu+rho={f['lam_mu_rho']}, "
                      f"|lam+mu+rho|^2={f['norm_sq']}, target={f['target_sq']}")


if __name__ == "__main__":
    main()
