"""Check whether the PRV-style kernel claim fails in non-trivial (lam != lam') cases."""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product
)
from dirac_triangle_proof_verification import (
    weyl_orbit, all_weyl_elements_as_reduced_words,
    apply_word, verify_kernel_claim,
)


def main():
    a2 = build_A(2)
    rho = a2.rho_omega

    print("Computing Weyl group of A_2...")
    weyl_words = all_weyl_elements_as_reduced_words(a2)
    print(f"  |W(A_2)| = {len(weyl_words)}")

    # All pairs with p+q <= 3 (small panel)
    panel = panel_dominant_weights(2, 3)
    print(f"Panel: {len(panel)} irreps")

    n_pairs_tested = 0
    n_failures_nontrivial = 0
    n_failures_trivial = 0  # lam == lam_pr
    n_failures_pair = 0
    failure_examples = []

    for lam in panel:
        for lam_pr in panel:
            n_pairs_tested += 1
            claim_ok, summary, diags = verify_kernel_claim(a2, lam, lam_pr, weyl_words)
            if not claim_ok:
                n_failures_pair += 1
                is_trivial = (lam == lam_pr)
                if is_trivial:
                    n_failures_trivial += 1
                else:
                    n_failures_nontrivial += 1
                    if len(failure_examples) < 20:
                        failed_sigmas = [d for d in diags if not d["claim_holds"]]
                        failure_examples.append({
                            "lam": lam, "lam_pr": lam_pr,
                            "failed_sigmas": failed_sigmas,
                        })

    print(f"\nTotal pairs: {n_pairs_tested}")
    print(f"Pairs with PRV-claim failure: {n_failures_pair}")
    print(f"  ... of which lam == lam_pr (trivial INT): {n_failures_trivial}")
    print(f"  ... of which lam != lam_pr (non-trivial!): {n_failures_nontrivial}")

    if failure_examples:
        print("\nNon-trivial failure examples:")
        for ex in failure_examples[:10]:
            print(f"  lam={ex['lam']}, lam_pr={ex['lam_pr']}:")
            for d in ex["failed_sigmas"]:
                # Compute C(sigma) and |lam-lam_pr|^2
                sigma = tuple(d["sigma"])
                C = a2.casimir(sigma)
                lam_diff = tuple(ex["lam"][i] - ex["lam_pr"][i] for i in range(2))
                lam_diff_sq = a2._inner_omega(lam_diff, lam_diff)
                holds = lam_diff_sq <= C
                print(f"    sigma={sigma}, mult={d['mult']}, "
                      f"|lam-lam'|^2={lam_diff_sq}, C={C}, INT holds: {holds}")


if __name__ == "__main__":
    main()
