"""
Test the WEAKER Steinberg-formula claim:

  For every sigma appearing in V_lambda (x) V_{lambda'}^* with positive
  multiplicity, there exists at least one mu in the weights of V_{lambda'^*}
  (not just the Weyl orbit of the HW) and one w in W such that

      w(sigma + rho) = lambda + mu + rho.

This is GUARANTEED by the Brauer-Klimyk/Steinberg formula: if mult(sigma) > 0,
some uncancelled term must contribute.

We test:
  (a) That this weaker claim holds (it must, structurally).
  (b) Whether the *contributing* mu satisfies |lambda + mu + rho|^2 >=
      |lambda - lambda'|^2 + |rho|^2.

If (b) holds always, that's the analytic proof.

Then we test alternatively:
  (c) For the lambda+mu+rho that's NORM-INVARIANT to sigma+rho (one of them
      must be, since they differ by a Weyl reflection), does the NORM
      satisfy |lambda+mu+rho|^2 - |lambda+mu+rho|*0 >= |lambda-lambda'|^2 + |rho|^2?

Actually, this is equivalent. The point is just that |w(sigma+rho)|^2 =
|sigma+rho|^2 always.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product,
    weights_of_irrep,
)
from dirac_triangle_proof_verification import (
    weyl_orbit, all_weyl_elements_as_reduced_words,
    apply_word, reflect,
)
import sympy as sp


def steinberg_witnesses(la, lam, lam_pr, sigma, weyl_words):
    """
    For (lam, lam'), find all (w, mu) pairs such that
        w(sigma + rho) = lam + mu + rho
    where mu is in the weights of V_{lam_pr^*}.

    Returns list of dicts: {w_word, mu, mult_of_mu, sign, lam_plus_mu_plus_rho_sq}.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_pr_dual = la.dual(lam_pr)
    weights_dual = weights_of_irrep(la, lam_pr_dual)

    sigma_plus_rho = tuple(sigma[i] + rho[i] for i in range(rank))
    witnesses = []

    for word in weyl_words:
        w_sigma_plus_rho = apply_word(la, word, sigma_plus_rho)
        # We need: w_sigma_plus_rho = lam + mu + rho for some mu in weights_dual
        # i.e., mu = w_sigma_plus_rho - lam - rho
        mu = tuple(w_sigma_plus_rho[i] - lam[i] - rho[i] for i in range(rank))
        if mu in weights_dual:
            mult = weights_dual[mu]
            # Sign from Weyl word length
            sign = (-1) ** len(word)
            lam_plus_mu_plus_rho = tuple(lam[i] + mu[i] + rho[i] for i in range(rank))
            norm_sq = la._inner_omega(lam_plus_mu_plus_rho, lam_plus_mu_plus_rho)
            witnesses.append({
                "w_word": list(word),
                "mu": list(mu),
                "mult_mu": str(mult),
                "sign": sign,
                "lam_plus_mu_plus_rho": list(lam_plus_mu_plus_rho),
                "norm_sq": str(norm_sq),
            })

    return witnesses


def main():
    a2 = build_A(2)
    rho = a2.rho_omega

    print("Computing Weyl group of A_2...")
    weyl_words = all_weyl_elements_as_reduced_words(a2)
    print(f"  |W(A_2)| = {len(weyl_words)}")

    # Test on a non-trivial case where PRV claim fails
    # lam=(0,1)=anti-fundamental, lam_pr=(0,2)
    lam = (0, 1)
    lam_pr = (0, 2)
    lam_pr_dual = a2.dual(lam_pr)
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"\nlam={lam}, lam_pr={lam_pr}, lam_pr_dual={lam_pr_dual}")
    print(f"Decomp: {decomp}")

    lam_diff = tuple(lam[i] - lam_pr[i] for i in range(2))
    lam_diff_sq = a2._inner_omega(lam_diff, lam_diff)
    rho_sq = a2.rho_sq
    print(f"|lam-lam'|^2 = {lam_diff_sq}")
    print(f"|rho|^2 = {rho_sq}")
    print(f"|lam-lam'|^2 + |rho|^2 = {lam_diff_sq + rho_sq}")

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho_sq = a2._inner_omega(
            tuple(sigma[i] + rho[i] for i in range(2)),
            tuple(sigma[i] + rho[i] for i in range(2))
        )
        print(f"\n  sigma={sigma}, mult={mult}, |sigma+rho|^2={sigma_plus_rho_sq}")
        print(f"    C(sigma) = {sigma_plus_rho_sq - rho_sq}")
        print(f"    INT holds: {lam_diff_sq <= sigma_plus_rho_sq - rho_sq}")
        # Find all Steinberg witnesses
        witnesses = steinberg_witnesses(a2, lam, lam_pr, sigma, weyl_words)
        print(f"    Steinberg witnesses ({len(witnesses)}):")
        for wit in witnesses:
            print(f"      mu={wit['mu']}, mult_mu={wit['mult_mu']}, "
                  f"sign={wit['sign']}, |lam+mu+rho|^2={wit['norm_sq']}")

    print()
    # Another case: lam=(0,1), lam_pr=(3,0)
    lam = (0, 1)
    lam_pr = (3, 0)
    lam_pr_dual = a2.dual(lam_pr)
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"\nlam={lam}, lam_pr={lam_pr}, lam_pr_dual={lam_pr_dual}")
    print(f"Decomp: {decomp}")
    lam_diff = tuple(lam[i] - lam_pr[i] for i in range(2))
    lam_diff_sq = a2._inner_omega(lam_diff, lam_diff)
    print(f"|lam-lam'|^2 = {lam_diff_sq}")
    print(f"|lam-lam'|^2 + |rho|^2 = {lam_diff_sq + rho_sq}")

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho_sq = a2._inner_omega(
            tuple(sigma[i] + rho[i] for i in range(2)),
            tuple(sigma[i] + rho[i] for i in range(2))
        )
        print(f"\n  sigma={sigma}, mult={mult}, |sigma+rho|^2={sigma_plus_rho_sq}")
        print(f"    C(sigma) = {sigma_plus_rho_sq - rho_sq}")
        print(f"    INT holds: {lam_diff_sq <= sigma_plus_rho_sq - rho_sq}")
        witnesses = steinberg_witnesses(a2, lam, lam_pr, sigma, weyl_words)
        print(f"    Steinberg witnesses ({len(witnesses)}):")
        for wit in witnesses:
            print(f"      mu={wit['mu']}, mult_mu={wit['mult_mu']}, "
                  f"sign={wit['sign']}, |lam+mu+rho|^2={wit['norm_sq']}")


if __name__ == "__main__":
    main()
