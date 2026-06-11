"""
Check: for all mu in V_{lambda'^*} that CONTRIBUTE to BK (i.e., lambda+mu+rho is
W-regular -- can be reflected to strict dominant interior without hitting walls),
does |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2?

This is the right "regular" criterion for the proof.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, weights_of_irrep,
)
from dirac_triangle_proof_verification import reflect


def is_w_regular_and_contributes(la, weight):
    """
    Check whether weight (in omega basis) is W-regular, i.e., its Weyl orbit
    contains a strictly dominant interior point.

    We reflect using the standard algorithm; if at any step a coordinate hits 0,
    the weight is on a wall.

    Returns (contributes_bool, dominant_image_if_contributes, sign).
    """
    rank = la.rank
    cur = list(weight)
    sign = 1
    max_iter = 200
    for _ in range(max_iter):
        if any(c == 0 for c in cur):
            return False, None, sign
        if all(c > 0 for c in cur):
            return True, tuple(cur), sign
        i_neg = -1
        for i in range(rank):
            if cur[i] < 0:
                i_neg = i
                break
        if i_neg < 0:
            break
        cur = list(reflect(la, tuple(cur), i_neg))
        sign *= -1
    return False, None, sign


def test_bk_contributing(la, lam, lam_pr):
    """
    For all BK-contributing mu in V_{lam_pr^*}: check |lam+mu+rho|^2 >= target.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    weights_dual = weights_of_irrep(la, lam_pr_dual)

    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    target_sq = la._inner_omega(lam_minus_lp_plus_rho, lam_minus_lp_plus_rho)

    fails = []
    for mu in weights_dual:
        lam_mu_rho = tuple(lam[i] + mu[i] + rho[i] for i in range(rank))
        contributes, dom_image, sign = is_w_regular_and_contributes(la, lam_mu_rho)
        if not contributes:
            continue
        norm_sq = la._inner_omega(lam_mu_rho, lam_mu_rho)
        if norm_sq < target_sq:
            fails.append({"mu": mu, "lam_mu_rho": lam_mu_rho,
                          "dom_image": dom_image, "sign": sign,
                          "norm_sq": norm_sq, "target_sq": target_sq})
    return fails


def main():
    print("Testing BK-CONTRIBUTING witness claim:")
    print("  For all mu in V_{lambda'^*} that contribute via BK:")
    print("      |lambda + mu + rho|^2 >= |lambda - lambda' + rho|^2")
    print()

    panels = [
        ("SU(3) p+q<=5", build_A(2), panel_dominant_weights(2, 5)),
        ("SU(4) a+b+c<=2", build_A(3), panel_dominant_weights(3, 2)),
        ("Sp(2)=C_2 a+b<=3", build_C2(), panel_dominant_weights(2, 3)),
        ("G_2 a+b<=2", build_G2(), panel_dominant_weights(2, 2)),
    ]

    for label, la, panel in panels:
        total_pairs = 0
        total_failures = 0
        first_failures = []
        for lam in panel:
            for lam_pr in panel:
                total_pairs += 1
                fails = test_bk_contributing(la, lam, lam_pr)
                if fails:
                    total_failures += 1
                    if len(first_failures) < 5:
                        first_failures.append((lam, lam_pr, fails[:3]))
        print(f"{label}: {total_pairs} pairs, {total_failures} with BK-contributing failures")
        for lam, lam_pr, fails in first_failures:
            print(f"  lam={lam}, lam'={lam_pr}:")
            for f in fails:
                print(f"    mu={f['mu']} -> lam+mu+rho={f['lam_mu_rho']} "
                      f"dom={f['dom_image']} sign={f['sign']}, "
                      f"norm^2={f['norm_sq']}, target^2={f['target_sq']}")


if __name__ == "__main__":
    main()
