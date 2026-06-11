"""Check INT inequality for the 8x8 case where PRV claim fails."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, tensor_product
)
import sympy as sp


def main():
    a2 = build_A(2)
    rho = a2.rho_omega

    # 8 x 8 case
    lam = (1, 1)
    lam_pr = (1, 1)
    lam_pr_dual = a2.dual(lam_pr)
    decomp = tensor_product(a2, lam, lam_pr_dual)
    print(f"8 x 8 decomp = {decomp}")

    lam_diff = tuple(lam[i] - lam_pr[i] for i in range(2))
    lam_diff_sq = a2._inner_omega(lam_diff, lam_diff)
    print(f"|lam-lam_pr|^2 = {lam_diff_sq}")

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        C = a2.casimir(sigma)
        print(f"  sigma={sigma}, mult={mult}, C(sigma)={C}, "
              f"|lam-lam'|^2 <= C? {lam_diff_sq <= C}")

    # The 8 (sigma=(1,1)) with mult=2: but |lam-lam_pr|^2 = 0 since lam=lam_pr,
    # so trivially holds. INT holds trivially!
    print()
    print("Note: lam = lam_pr means |lam-lam_pr|^2 = 0, INT is trivial.")
    print()

    # Find a case where PRV claim might fail AND (INT) is non-trivial.
    # Try non-trivial mu-overlapping cases.
    cases = [
        ((1, 1), (2, 2)),
        ((2, 1), (1, 2)),
        ((1, 2), (2, 1)),
        ((3, 0), (1, 1)),
        ((1, 1), (3, 0)),
    ]
    for lam, lam_pr in cases:
        lam_pr_dual = a2.dual(lam_pr)
        decomp = tensor_product(a2, lam, lam_pr_dual)
        lam_diff = tuple(lam[i] - lam_pr[i] for i in range(2))
        lam_diff_sq = a2._inner_omega(lam_diff, lam_diff)
        print(f"\nlam={lam}, lam_pr={lam_pr}, |diff|^2={lam_diff_sq}, decomp={decomp}")


if __name__ == "__main__":
    main()
