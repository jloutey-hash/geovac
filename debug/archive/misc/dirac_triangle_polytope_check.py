"""
Investigate the structure: the set of sigmas appearing in V_lambda (x) V_{lambda'^*}
form a convex polytope. We check:

1. The minimum norm |sigma + rho|^2 over the polytope = |lambda - lambda' + rho|^2
   exactly (= dominant image norm).

2. The function sigma -> |sigma + rho|^2 is convex on the polytope.

3. The polytope is contained in [lambda + conv(W * lambda'^*)] intersected with
   the dominant chamber.

If we can prove (3) and use convexity, then maybe the minimum is at a vertex
and we can identify the vertex.

But actually all we need is the inequality |sigma + rho|^2 >= |lambda - lambda' + rho|^2.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product, weights_of_irrep
)
from dirac_triangle_minimum_norm import dominant_image
from dirac_triangle_proof_verification import weyl_orbit


def analyze_polytope(la, lam, lam_pr):
    """
    Compute the support polytope and the minimum norm.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    decomp = tensor_product(la, lam, lam_pr_dual)

    sigmas = [s for s, m in decomp.items() if m > 0]
    norms = []
    for s in sigmas:
        s_plus_rho = tuple(s[i] + rho[i] for i in range(rank))
        norms.append(la._inner_omega(s_plus_rho, s_plus_rho))

    if not norms:
        return None

    min_norm = min(norms)
    min_sigma = sigmas[norms.index(min_norm)]

    target_norm = la._inner_omega(
        tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank)),
        tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    )

    # Check: does min_sigma = dominant_image(lambda - lambda' + rho) - rho?
    lam_minus_lp_plus_rho = tuple(lam[i] - lam_pr[i] + rho[i] for i in range(rank))
    dom_image, sign = dominant_image(la, lam_minus_lp_plus_rho)
    if dom_image is not None:
        expected_min_sigma = tuple(dom_image[i] - rho[i] for i in range(rank))
    else:
        expected_min_sigma = None

    return {
        "min_norm_in_decomp": min_norm,
        "min_sigma": min_sigma,
        "target_norm": target_norm,
        "match": (min_norm == target_norm),
        "expected_min_sigma": expected_min_sigma,
    }


def main():
    print("Polytope minimum-norm verification:")
    print()

    la = build_A(2)
    panel = panel_dominant_weights(2, 4)

    print(f"SU(3) panel: {len(panel)} irreps")

    matches = 0
    non_matches = 0
    for lam in panel:
        for lam_pr in panel:
            res = analyze_polytope(la, lam, lam_pr)
            if res is None:
                continue
            if res["match"]:
                matches += 1
            else:
                non_matches += 1
                if non_matches <= 5:
                    print(f"  Non-match: lam={lam}, lam'={lam_pr}, "
                          f"min_norm={res['min_norm_in_decomp']}, target={res['target_norm']}")

    print(f"\nMatches: {matches}, Non-matches: {non_matches}")
    print()

    # SU(4) panel
    print("SU(4) panel:")
    la = build_A(3)
    panel = panel_dominant_weights(3, 3)
    matches = 0
    non_matches = 0
    for lam in panel:
        for lam_pr in panel:
            res = analyze_polytope(la, lam, lam_pr)
            if res is None:
                continue
            if res["match"]:
                matches += 1
            else:
                non_matches += 1
                if non_matches <= 5:
                    print(f"  Non-match: lam={lam}, lam'={lam_pr}, "
                          f"min_norm={res['min_norm_in_decomp']}, target={res['target_norm']}")
    print(f"  Matches: {matches}, Non-matches: {non_matches}")


if __name__ == "__main__":
    main()
