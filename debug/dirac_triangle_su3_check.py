"""
Numerical verification of the Dirac-triangle inequality on SU(3).

After P-A established that the Casimir triangle inequality
    |C(pi) - C(pi')| <= C(sigma)         (CASIMIR-TRIANGLE; FALSE)
fails generically, we test the corrected L3 reformulation:

    |D(pi) - D(pi')| <= sqrt(C(sigma))    (DIRAC-TRIANGLE; CONJECTURED)

where D is the Kostant cubic Dirac eigenvalue and the Casimir convention is the
same as P-A's casimir_triangle_su3.py:

    C(p, q) = (p^2 + p*q + q^2 + 3p + 3q) / 3

In this normalization, C(adjoint) = C(1,1) = 3, and the Kostant cubic Dirac on
the bi-invariant compact Lie group SU(3) has

    |D(lambda)|^2 = <lambda + rho, lambda + rho>

which in agent's Casimir normalization works out to

    |D(lambda)|^2 = C(lambda) + <rho, rho>  with <rho, rho> = 1
                  = C(lambda) + 1

(Derivation: C(lambda) = <lambda+rho, lambda+rho> - <rho, rho>; adjoint has
lambda = theta = rho, so C(adjoint) = <2rho, 2rho> - <rho, rho> = 3<rho, rho>.
Since C(adjoint) = 3 in this normalization, <rho, rho> = 1.)

For SU(2) the analog is |D(j)|^2 = C(j) + 1/4 = (j + 1/2)^2, matching Paper 38.

We test the inequality on the same panel P-A used and report pass/fail with
ratios.
"""

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, nsimplify
from itertools import product
from typing import List, Tuple, Dict

# Import P-A's machinery
from casimir_triangle_su3 import (
    casimir_su3,
    dim_su3,
    tensor_product_su3,
    weights_of_irrep_su3,
)


# ---------------------------------------------------------------------------
# Dirac eigenvalue
# ---------------------------------------------------------------------------

def D_squared_su3(p: int, q: int) -> Rational:
    """
    Kostant cubic Dirac eigenvalue squared on SU(3) irrep (p, q).

    In agent's Casimir normalization (C(adjoint) = 3),
        |D(lambda)|^2 = C(lambda) + <rho, rho> = C(lambda) + 1.
    """
    return casimir_su3(p, q) + Rational(1)


def D_su3(p, q):
    """|D(lambda)| as sympy expression (may be irrational)."""
    return sqrt(D_squared_su3(p, q))


# ---------------------------------------------------------------------------
# Dirac-triangle inequality check
# ---------------------------------------------------------------------------

def verify_dirac_triangle(lam: Tuple[int, int], lam_prime: Tuple[int, int]
                          ) -> Tuple[bool, list, list]:
    """
    Test |D(pi) - D(pi')| <= sqrt(C(sigma)) for every sigma in pi (x) pi'^*.

    Returns (all_pass, violations, all_ratios) where each violation is
    (sigma, sqrt_C_sigma, |D - D'|, ratio).
    """
    p, q = lam
    p_pr, q_pr = lam_prime
    lam_prime_dual = (q_pr, p_pr)

    D = D_su3(p, q)
    D_pr = D_su3(p_pr, q_pr)
    dirac_diff = simplify(sp.Abs(D - D_pr))

    decomp = tensor_product_su3(lam, lam_prime_dual)

    violations = []
    all_ratios = []
    all_pass = True

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        C_sigma = casimir_su3(*sigma)
        if C_sigma == 0:
            # trivial rep; sqrt(0) = 0; only valid if D = D' (i.e., pi = pi')
            if dirac_diff != 0:
                all_pass = False
                violations.append((sigma, Rational(0), dirac_diff, sp.oo))
            else:
                all_ratios.append(Rational(0))
            continue
        sqrt_C_sigma = sqrt(C_sigma)
        ratio = dirac_diff / sqrt_C_sigma
        ratio_num = float(ratio)
        all_ratios.append(ratio_num)
        if ratio_num > 1.0 + 1e-12:
            all_pass = False
            violations.append((sigma, sqrt_C_sigma, dirac_diff, ratio_num))

    return all_pass, violations, all_ratios


# ---------------------------------------------------------------------------
# Run on P-A's panel
# ---------------------------------------------------------------------------

def run_panel():
    print("=" * 78)
    print("SU(3) Dirac-triangle inequality verification")
    print("=" * 78)

    # P-A's panel was all (p,q) with p+q <= 3, 100 ordered pairs total
    irreps = []
    for p in range(4):
        for q in range(4):
            if p + q <= 3:
                irreps.append((p, q))

    print(f"\nPanel: {len(irreps)} irreps with p+q <= 3:")
    print(f"  {irreps}")
    print(f"  Total ordered pairs: {len(irreps) ** 2} (incl. diagonal)")

    # Sanity check
    print("\n--- Sanity: Dirac eigenvalues ---")
    for lam in irreps:
        p, q = lam
        D2 = D_squared_su3(p, q)
        D = D_su3(p, q)
        print(f"  ({p},{q}): C={casimir_su3(p,q)}, D^2={D2}, |D|={float(D):.4f}")

    print("\n--- Comparing to P-A's known Casimir-triangle counterexamples ---")
    pa_counterexamples = [
        ((3, 0), (1, 0)),  # P-A headline; Casimir ratio 1.40
        ((0, 2), (0, 3)),  # Casimir ratio 2.0
        ((3, 0), (2, 0)),  # Casimir ratio 2.0
    ]
    for lam, lam_pr in pa_counterexamples:
        passes, viols, ratios = verify_dirac_triangle(lam, lam_pr)
        d = simplify(sp.Abs(D_su3(*lam) - D_su3(*lam_pr)))
        print(f"  lambda={lam}, lambda'={lam_pr}: |D-D'|={float(d):.4f}")
        decomp = tensor_product_su3(lam, (lam_pr[1], lam_pr[0]))
        for sigma, mult in decomp.items():
            if mult > 0:
                Cs = casimir_su3(*sigma)
                if Cs > 0:
                    ratio = float(d / sqrt(Cs))
                    flag = "PASS" if ratio <= 1.0 else "FAIL"
                    print(f"      sigma={sigma}, sqrt(C)={float(sqrt(Cs)):.4f}, "
                          f"ratio={ratio:.4f} [{flag}]")

    # Full panel scan
    print("\n--- Full panel scan ---")
    fail_count = 0
    pass_count = 0
    max_ratio = 0.0
    max_ratio_case = None
    sup_ratio_per_pair = []  # max ratio over sigma in this pair

    for lam in irreps:
        for lam_pr in irreps:
            passes, viols, ratios = verify_dirac_triangle(lam, lam_pr)
            if not passes:
                fail_count += 1
                print(f"  FAIL: lambda={lam}, lambda'={lam_pr}")
                for v in viols:
                    sigma, sqrtC, diff, ratio = v
                    print(f"    sigma={sigma}: ratio={ratio:.4f}")
            else:
                pass_count += 1
            if ratios:
                pair_sup = max(ratios) if ratios else 0
                if pair_sup > max_ratio:
                    max_ratio = pair_sup
                    max_ratio_case = (lam, lam_pr)
                sup_ratio_per_pair.append((lam, lam_pr, pair_sup))

    print("\n--- Summary ---")
    print(f"  Pass: {pass_count} / {pass_count + fail_count}")
    print(f"  Fail: {fail_count}")
    print(f"  Max ratio in panel: {max_ratio:.6f}")
    if max_ratio_case:
        lam, lam_pr = max_ratio_case
        print(f"  Achieved at: lambda={lam}, lambda'={lam_pr}")
        decomp = tensor_product_su3(lam, (lam_pr[1], lam_pr[0]))
        print(f"  Decomposition: {decomp}")
        D = D_su3(*lam)
        D_pr = D_su3(*lam_pr)
        diff = simplify(sp.Abs(D - D_pr))
        print(f"  |D - D'| = {diff} = {float(diff):.6f}")
        for sigma, mult in decomp.items():
            if mult > 0 and casimir_su3(*sigma) > 0:
                Cs = casimir_su3(*sigma)
                ratio = float(diff / sqrt(Cs))
                print(f"    sigma={sigma}: C={Cs}, sqrt(C)={float(sqrt(Cs)):.4f}, "
                      f"ratio={ratio:.4f}")

    # Asymptotic behavior: family (n, 0) (x) (1, 0) for n = 1, 2, ..., 8
    print("\n--- Asymptotic family: (n, 0) vs (1, 0) ---")
    print(f"  {'n':>3} {'|D - D_1|':>12} {'sqrt(C(sigma_min))':>20} {'ratio':>10}")
    for n in range(1, 9):
        lam = (n, 0)
        lam_pr = (1, 0)
        decomp = tensor_product_su3(lam, (0, 1))  # (n,0) (x) (0,1)
        # find sigma with smallest C(sigma) > 0
        smallest_C = None
        smallest_sigma = None
        for sigma, mult in decomp.items():
            if mult > 0:
                Cs = casimir_su3(*sigma)
                if Cs > 0 and (smallest_C is None or Cs < smallest_C):
                    smallest_C = Cs
                    smallest_sigma = sigma
        if smallest_C is None:
            continue
        D = D_su3(n, 0)
        D_pr = D_su3(1, 0)
        diff = float(simplify(sp.Abs(D - D_pr)))
        sqrt_Cs = float(sqrt(smallest_C))
        ratio = diff / sqrt_Cs
        print(f"  {n:>3} {diff:>12.4f} {sqrt_Cs:>20.4f} {ratio:>10.4f}    "
              f"(sigma_min={smallest_sigma})")

    return pass_count, fail_count, max_ratio


if __name__ == "__main__":
    run_panel()
