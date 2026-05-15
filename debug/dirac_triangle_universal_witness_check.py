"""
Test the UNIVERSAL stronger claim:

  For every weight mu in V_{lambda'^*}, we have:
      |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2

This is independent of which sigma mu witnesses. If true, the proof is
elementary: just use that -lambda' is a weight of V_{lambda'^*}.

In fact, this is exactly the statement that
  |lambda + mu + rho|^2 is minimized over mu in V_{lambda'^*} at mu = -lambda'.

Equivalently:
  <lambda + rho, mu> is minimized at mu = -lambda' (and the value is -<lambda+rho, lambda'>).
And then by quadratic norm:
  |mu|^2 must compensate.

Actually let me just verify it numerically first.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, weights_of_irrep,
)


def test_universal_witness(la, lam, lam_pr):
    """
    For all mu in weights of V_{lam_pr^*}, check |lam + mu + rho|^2 >= |lam - lam' + rho|^2.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    weights_dual = weights_of_irrep(la, lam_pr_dual)

    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    target_sq = la._inner_omega(lam_minus_lp_plus_rho, lam_minus_lp_plus_rho)

    fail_examples = []
    for mu in weights_dual:
        lam_mu_rho = tuple(lam[i] + mu[i] + rho[i] for i in range(rank))
        norm_sq = la._inner_omega(lam_mu_rho, lam_mu_rho)
        if norm_sq < target_sq:
            fail_examples.append({
                "mu": mu, "norm_sq": norm_sq, "target_sq": target_sq,
            })
    return fail_examples


def main():
    print("Testing universal weight-witness claim:")
    print("  For all mu in weights of V_{lambda'^*}:")
    print("      |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2")
    print()
    print("If true, the proof is elementary.")
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
                fails = test_universal_witness(la, lam, lam_pr)
                if fails:
                    total_failures += 1
                    if len(first_failures) < 5:
                        first_failures.append((lam, lam_pr, fails[:3]))
        print(f"{label}: {total_pairs} pairs, {total_failures} pairs with at least one failing mu")
        if first_failures:
            for lam, lam_pr, fails in first_failures:
                print(f"  lam={lam}, lam'={lam_pr}:")
                for f in fails:
                    print(f"    mu={f['mu']}, |lam+mu+rho|^2={f['norm_sq']}, "
                          f"target={f['target_sq']}")


if __name__ == "__main__":
    main()
