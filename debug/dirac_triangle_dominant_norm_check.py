"""
Test the stronger inequality:
    |sigma + rho|^2 >= |lambda - lambda' + rho|^2 ??

If yes, this is the analytic proof: |sigma+rho|^2 - |rho|^2 = C(sigma) >=
|lambda - lambda' + rho|^2 - |rho|^2 = |lam - lam'|^2 + 2<lam-lam', rho>.

But we need C(sigma) >= |lam-lam'|^2, so we'd need 2<lam-lam', rho> <= 0,
which often fails.

The right direction: try
    |sigma + rho|^2 >= |lam - lam'|^2 + |rho|^2  (equivalent to INT)
which means
    |sigma + rho|^2 - |lam - lam' + rho|^2 >= -2<lam-lam', rho>
But this is just C(sigma) - |lam - lam'|^2 >= 0, the original inequality.

The cleanest reformulation, due to the dominance property of sigma+rho:

  Since sigma+rho is the DOMINANT representative of its Weyl orbit (the W-orbit
  of lam+mu+rho), and rho is regular dominant, we have

      <sigma+rho, lam+rho> >= <lam+mu+rho, lam+rho>
      (the dominant rep maximizes inner products with dominant weights)

  More useful: <sigma+rho, sigma+rho> >= <lam-lam', lam-lam'> ???

  This is the question. Test it.
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product,
)
import sympy as sp


def check_lam_minus_lp_in_decomp(la, lam, lam_pr):
    """
    For each sigma in V_lam (x) V_{lam_pr}^*, check various norm inequalities.
    """
    rank = la.rank
    rho = la.rho_omega
    rho_sq = la.rho_sq
    lam_pr_dual = la.dual(lam_pr)
    decomp = tensor_product(la, lam, lam_pr_dual)

    lam_minus_lp = tuple(lam[i] - lam_pr[i] for i in range(rank))
    lam_minus_lp_sq = la._inner_omega(lam_minus_lp, lam_minus_lp)

    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    lam_minus_lp_plus_rho_sq = la._inner_omega(lam_minus_lp_plus_rho, lam_minus_lp_plus_rho)

    results = []
    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho = tuple(sigma[i] + rho[i] for i in range(rank))
        sigma_plus_rho_sq = la._inner_omega(sigma_plus_rho, sigma_plus_rho)
        c_sigma = sigma_plus_rho_sq - rho_sq

        # Two questions:
        # 1. Is C(sigma) >= |lam - lam'|^2? (INT)
        # 2. Is |sigma+rho|^2 >= |lam - lam' + rho|^2? (stronger, would imply INT
        #    if <lam-lam', rho> <= 0)
        int_holds = lam_minus_lp_sq <= c_sigma
        stronger_holds = sigma_plus_rho_sq >= lam_minus_lp_plus_rho_sq

        results.append({
            "sigma": sigma,
            "mult": mult,
            "c_sigma": c_sigma,
            "lam_diff_sq": lam_minus_lp_sq,
            "int": int_holds,
            "lam_diff_rho_sq": lam_minus_lp_plus_rho_sq,
            "sigma_rho_sq": sigma_plus_rho_sq,
            "stronger": stronger_holds,
        })

    return results


def main():
    a2 = build_A(2)

    print("Testing whether |sigma+rho|^2 >= |lam-lam'+rho|^2 for all sigma in tensor product")
    print()

    panel = panel_dominant_weights(2, 4)
    print(f"SU(3) panel: {len(panel)} irreps")

    total_sigma = 0
    int_failures = 0
    stronger_failures = []
    for lam in panel:
        for lam_pr in panel:
            results = check_lam_minus_lp_in_decomp(a2, lam, lam_pr)
            for r in results:
                total_sigma += 1
                if not r["int"]:
                    int_failures += 1
                if not r["stronger"]:
                    stronger_failures.append({
                        "lam": lam, "lam_pr": lam_pr, "result": r,
                    })

    print(f"\nTotal sigmas: {total_sigma}")
    print(f"INT failures: {int_failures}")
    print(f"Stronger inequality failures: {len(stronger_failures)}")
    if stronger_failures:
        print("\nFirst few failures of |sigma+rho|^2 >= |lam-lam'+rho|^2:")
        for f in stronger_failures[:8]:
            r = f["result"]
            print(f"  lam={f['lam']}, lam'={f['lam_pr']}, sigma={r['sigma']}: "
                  f"|sigma+rho|^2={r['sigma_rho_sq']} < |lam-lam'+rho|^2={r['lam_diff_rho_sq']}")


if __name__ == "__main__":
    main()
