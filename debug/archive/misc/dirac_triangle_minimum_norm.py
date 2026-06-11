"""
Test: is sigma_min := dominant_image(lambda - lambda' + rho) - rho
always in the decomposition with positive multiplicity, and does it have
the MINIMUM |sigma+rho|^2 norm among all sigmas in the decomposition?

This would close the proof: every sigma has |sigma+rho|^2 >= |sigma_min+rho|^2 =
dominant image norm = |lambda - lambda' + rho|^2.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product,
)
from dirac_triangle_proof_verification import reflect


def dominant_image(la, weight):
    """
    Return the dominant chamber image of `weight` and the sign of the Weyl element
    that brings it there.
    Returns (dom_weight_or_None_if_wall, sign).
    """
    rank = la.rank
    cur = list(weight)
    sign = 1
    max_iter = 200
    for _ in range(max_iter):
        if any(c == 0 for c in cur):
            return None, sign  # On wall
        if all(c > 0 for c in cur):
            return tuple(cur), sign
        # find first negative
        i_neg = -1
        for i in range(rank):
            if cur[i] < 0:
                i_neg = i
                break
        if i_neg < 0:
            return tuple(cur), sign
        cur = list(reflect(la, tuple(cur), i_neg))
        sign *= -1
    return None, sign


def test_minimum_norm(la, lam, lam_pr):
    """
    1. Compute sigma_min := dominant_image(lambda - lambda' + rho) - rho.
    2. Verify sigma_min is in V_lambda (x) V_{lambda'}^* with positive multiplicity.
    3. Verify all other sigmas in the decomp have |sigma+rho|^2 >= |sigma_min+rho|^2.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    decomp = tensor_product(la, lam, lam_pr_dual)

    # Compute sigma_min via dominant image
    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    dom_image, sign = dominant_image(la, lam_minus_lp_plus_rho)

    if dom_image is None:
        # On wall - sigma_min not well-defined; can still ask whether OTHER sigmas
        # all satisfy |sigma+rho|^2 >= |lambda - lambda' + rho|^2 (this is what we want)
        target_norm_sq = la._inner_omega(lam_minus_lp_plus_rho, lam_minus_lp_plus_rho)
        results = {"on_wall": True, "target_norm_sq": str(target_norm_sq),
                   "min_sigma_norm_sq": None, "min_sigma_present": False,
                   "all_sigma_meet_target": True, "decomp_size": 0}
        min_norm_sq_in_decomp = None
        for sigma, mult in decomp.items():
            if mult <= 0:
                continue
            results["decomp_size"] += 1
            sigma_plus_rho_sq = la._inner_omega(
                tuple(sigma[i] + rho[i] for i in range(rank)),
                tuple(sigma[i] + rho[i] for i in range(rank))
            )
            if min_norm_sq_in_decomp is None or sigma_plus_rho_sq < min_norm_sq_in_decomp:
                min_norm_sq_in_decomp = sigma_plus_rho_sq
            if sigma_plus_rho_sq < target_norm_sq:
                results["all_sigma_meet_target"] = False
        results["min_sigma_norm_sq"] = str(min_norm_sq_in_decomp)
        return results

    sigma_min = tuple(dom_image[i] - rho[i] for i in range(rank))
    sigma_min_norm_sq = la._inner_omega(dom_image, dom_image)

    # Check sigma_min in decomp?
    sigma_min_present = (sigma_min in decomp and decomp[sigma_min] > 0)

    # Check all sigmas have norm >= sigma_min_norm_sq
    all_meet_target = True
    min_norm_in_decomp = None
    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho_sq = la._inner_omega(
            tuple(sigma[i] + rho[i] for i in range(rank)),
            tuple(sigma[i] + rho[i] for i in range(rank))
        )
        if min_norm_in_decomp is None or sigma_plus_rho_sq < min_norm_in_decomp:
            min_norm_in_decomp = sigma_plus_rho_sq
        if sigma_plus_rho_sq < sigma_min_norm_sq:
            all_meet_target = False

    return {
        "on_wall": False,
        "sigma_min": list(sigma_min),
        "sigma_min_sign": sign,
        "sigma_min_present_in_decomp": sigma_min_present,
        "sigma_min_norm_sq": str(sigma_min_norm_sq),
        "all_sigma_meet_target": all_meet_target,
        "min_norm_in_decomp": str(min_norm_in_decomp),
        "decomp_size": sum(1 for s, m in decomp.items() if m > 0),
    }


def main():
    print("Testing minimum-norm hypothesis on all panels.")
    print()
    print("Hypothesis: for every (lam, lam'), there is a 'sigma_min' = dominant_image")
    print("(lam - lam' + rho) - rho such that:")
    print("  (a) sigma_min is in V_lam (x) V_{lam'}^* with positive multiplicity")
    print("  (b) every other sigma in the decomp has |sigma+rho|^2 >= |sigma_min+rho|^2")
    print()

    panels = [
        ("SU(3) p+q<=5", build_A(2), panel_dominant_weights(2, 5)),
        ("SU(4) a+b+c<=3", build_A(3), panel_dominant_weights(3, 3)),
        ("Sp(2)=C_2 a+b<=3", build_C2(), panel_dominant_weights(2, 3)),
        ("G_2 a+b<=2", build_G2(), panel_dominant_weights(2, 2)),
    ]

    for label, la, panel in panels:
        sigma_min_absent = 0
        not_meet_target_count = 0
        total_pairs = 0
        on_wall_count = 0

        for lam in panel:
            for lam_pr in panel:
                total_pairs += 1
                res = test_minimum_norm(la, lam, lam_pr)
                if res.get("on_wall"):
                    on_wall_count += 1
                    if not res["all_sigma_meet_target"]:
                        not_meet_target_count += 1
                else:
                    if not res["sigma_min_present_in_decomp"]:
                        sigma_min_absent += 1
                    if not res["all_sigma_meet_target"]:
                        not_meet_target_count += 1

        print(f"{label}: {total_pairs} pairs")
        print(f"  on_wall (sigma_min not well-defined): {on_wall_count}")
        print(f"  sigma_min present in decomp: {total_pairs - on_wall_count - sigma_min_absent}/{total_pairs - on_wall_count}")
        print(f"  All sigmas meet target |lam-lam'+rho|^2: {total_pairs - not_meet_target_count}/{total_pairs}")
        print()


if __name__ == "__main__":
    main()
