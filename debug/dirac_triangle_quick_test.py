"""Quick test of the PRV-style kernel claim on SU(3) small panel."""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, panel_dominant_weights, tensor_product
)
from dirac_triangle_proof_verification import (
    weyl_orbit, all_weyl_elements_as_reduced_words,
    apply_word, verify_kernel_claim,
)


def main():
    a2 = build_A(2)
    rho = a2.rho_omega

    # SU(3) Weyl group has order 6
    print("Computing Weyl group of A_2...")
    t0 = time.time()
    weyl_words = all_weyl_elements_as_reduced_words(a2)
    print(f"  |W(A_2)| = {len(weyl_words)} (expected 6) in {time.time()-t0:.2f}s")
    for w in weyl_words:
        print(f"    word {w}: action on rho = {apply_word(a2, w, rho)}")
    print()

    # Test the worst case from previous sprint: lambda=(0,0), lambda'=(0,5)
    lam = (0, 0)
    lam_pr = (0, 5)
    print(f"Testing kernel claim at lam={lam}, lam'={lam_pr}")
    lam_pr_dual = a2.dual(lam_pr)
    print(f"  lam_pr_dual = {lam_pr_dual}")
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"  Decomposition: {decomp}")

    # -lambda' = (0, -5); compute its Weyl orbit
    minus_lp = tuple(-x for x in lam_pr)
    print(f"  -lam_pr = {minus_lp}")
    orbit = weyl_orbit(a2, minus_lp)
    print(f"  Orbit of -lam_pr: |orbit|={len(orbit)}; {orbit}")

    # Now verify the kernel claim
    claim_ok, summary, diags = verify_kernel_claim(a2, lam, lam_pr, weyl_words)
    print(f"  Claim holds: {claim_ok}")
    print(f"  Summary: {summary}")
    for d in diags:
        print(f"    sigma={tuple(d['sigma'])}, mult={d['mult']}, "
              f"claim={d['claim_holds']}, n_witnesses={d['n_witnesses']}")

    # Test a more typical case: lam=(2,0), lam'=(0,1)
    print()
    lam = (2, 0)
    lam_pr = (0, 1)
    print(f"Testing kernel claim at lam={lam}, lam'={lam_pr}")
    lam_pr_dual = a2.dual(lam_pr)
    print(f"  lam_pr_dual = {lam_pr_dual}")
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"  Decomposition: {decomp}")
    claim_ok, summary, diags = verify_kernel_claim(a2, lam, lam_pr, weyl_words)
    print(f"  Claim holds: {claim_ok}")
    for d in diags:
        print(f"    sigma={tuple(d['sigma'])}, mult={d['mult']}, "
              f"claim={d['claim_holds']}, n_witnesses={d['n_witnesses']}")

    # 8 x 8 case: lam=(1,1), lam'=(1,1)
    print()
    lam = (1, 1)
    lam_pr = (1, 1)
    print(f"Testing kernel claim at lam={lam}, lam'={lam_pr}")
    lam_pr_dual = a2.dual(lam_pr)
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"  Decomposition: {decomp}")
    claim_ok, summary, diags = verify_kernel_claim(a2, lam, lam_pr, weyl_words)
    print(f"  Claim holds: {claim_ok}")
    for d in diags:
        print(f"    sigma={tuple(d['sigma'])}, mult={d['mult']}, "
              f"claim={d['claim_holds']}, n_witnesses={d['n_witnesses']}")


if __name__ == "__main__":
    main()
