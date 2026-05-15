"""
Test: for every Brauer-Klimyk witness (w, mu) of a contributing sigma,
does |lam + mu + rho|^2 >= |lam - lam'|^2 + |rho|^2 hold?

If yes for ALL witnesses, the proof is trivial. If yes only for SOME witnesses,
we need to identify the "good" witness systematically.

This script runs the test exhaustively on small panels.
"""

import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product,
    weights_of_irrep,
)
from dirac_triangle_proof_verification import (
    all_weyl_elements_as_reduced_words, apply_word,
)
from dirac_triangle_steinberg_check import steinberg_witnesses
import sympy as sp


def test_witness_norm_bound(la, lam, lam_pr, weyl_words):
    """
    For each sigma in the decomp, find all Brauer-Klimyk witnesses (w, mu).
    Check whether AT LEAST ONE witness has |lam+mu+rho|^2 >= |lam-lam'|^2 + |rho|^2.
    Equivalently: |sigma+rho|^2 >= |lam-lam'|^2 + |rho|^2, i.e., C(sigma) >= |lam-lam'|^2.

    But since |w(x)|^2 = |x|^2 (Weyl invariance), |lam+mu+rho|^2 = |sigma+rho|^2 for
    every witness. So this check is equivalent to (INT) holding, NOT to any
    per-witness condition.

    The real question is different: which witness contributes to sigma, and does
    that witness give a Cauchy-Schwarz-like bound on |lam-mu|?

    Let's also check: among the contributing mu's, is at least one of them
    "extremal" (i.e., lies in W*(lam_pr^*) = W*(-w_0 lam_pr)).
    """
    rank = la.rank
    rho = la.rho_omega
    rho_sq = la.rho_sq
    lam_pr_dual = la.dual(lam_pr)
    weights_dual = weights_of_irrep(la, lam_pr_dual)
    extremal_dual_weights = set()
    # Extremal weights of V_{lam_pr^*} are the Weyl orbit of lam_pr^* (the HW).
    from dirac_triangle_proof_verification import weyl_orbit
    for w in weyl_orbit(la, lam_pr_dual):
        if w in weights_dual:
            extremal_dual_weights.add(w)

    lam_pr_diff = tuple(lam[i] - lam_pr[i] for i in range(rank))
    lam_diff_sq = la._inner_omega(lam_pr_diff, lam_pr_diff)

    decomp = tensor_product(la, lam, lam_pr_dual)

    sigma_data = []
    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho = tuple(sigma[i] + rho[i] for i in range(rank))
        sigma_plus_rho_sq = la._inner_omega(sigma_plus_rho, sigma_plus_rho)
        c_sigma = sigma_plus_rho_sq - rho_sq

        # Get all witnesses
        witnesses = steinberg_witnesses(la, lam, lam_pr, sigma, weyl_words)

        # Check: among witnesses, is mu extremal (i.e., in W*lam_pr^*)?
        n_extremal_witnesses = sum(1 for w in witnesses if tuple(w["mu"]) in extremal_dual_weights)
        # Check: the maximum |lam+mu+rho|^2 over witnesses (all equal to sigma+rho)
        # OK so by Weyl invariance, all witnesses have |lam+mu+rho|^2 = sigma_plus_rho_sq.

        sigma_data.append({
            "sigma": list(sigma),
            "mult": mult,
            "c_sigma": str(c_sigma),
            "lam_diff_sq": str(lam_diff_sq),
            "int_holds": (lam_diff_sq <= c_sigma),
            "n_witnesses": len(witnesses),
            "n_extremal": n_extremal_witnesses,
        })

    return sigma_data, lam_diff_sq


def main():
    a2 = build_A(2)
    rho = a2.rho_omega

    print("Computing W(A_2)...")
    weyl_words = all_weyl_elements_as_reduced_words(a2)
    print(f"  |W| = {len(weyl_words)}")

    # Small panel test
    panel = panel_dominant_weights(2, 3)
    print(f"\nPanel: {len(panel)} irreps, {len(panel)**2} pairs")

    total_sigma = 0
    sigma_no_extremal = []  # sigma's with no extremal witness
    sigma_int_failures = 0  # should be 0

    for lam in panel:
        for lam_pr in panel:
            sigma_data, lam_diff_sq = test_witness_norm_bound(a2, lam, lam_pr, weyl_words)
            for s in sigma_data:
                total_sigma += 1
                if not s["int_holds"]:
                    sigma_int_failures += 1
                if s["n_extremal"] == 0:
                    sigma_no_extremal.append((lam, lam_pr, s))

    print(f"\nTotal sigmas in decompositions: {total_sigma}")
    print(f"INT failures: {sigma_int_failures}")
    print(f"Sigmas with NO extremal witness: {len(sigma_no_extremal)}")
    if sigma_no_extremal:
        print("\nFirst few sigmas with no extremal witness:")
        for lam, lam_pr, s in sigma_no_extremal[:15]:
            print(f"  lam={lam}, lam_pr={lam_pr}, sigma={s['sigma']}, mult={s['mult']}, "
                  f"INT holds: {s['int_holds']}")


if __name__ == "__main__":
    main()
