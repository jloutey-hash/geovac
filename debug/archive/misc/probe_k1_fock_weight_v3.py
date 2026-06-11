"""
Probe v3: The coupling^2 pattern and exact rational formula.

From v2, the coupling^2 values are:
  l=0: 1/16 (all n)
  l=1: 1/24, 1/19.2, 9/160, 1/17.14...
  l=2: 1/32, 1/22.86, ...
  etc.

Let me compute these as exact rationals and find the closed form.
Also: verify the Fock Jacobian derivation of kappa more rigorously.
"""

from __future__ import annotations

import json
import os

import sympy as sp
from sympy import (
    Rational, Symbol, sqrt, pi, simplify, factorial, gamma, gegenbauer,
    integrate, cancel, together, nsimplify
)

u = Symbol('u', real=True)


def gegenbauer_norm(k, lam):
    """Exact Gegenbauer norm h_k^lam."""
    return (pi * 2**(1 - 2*lam) * gamma(k + 2*lam)
            / (factorial(k) * (k + lam) * gamma(lam)**2))


def coupling_squared_exact(n, l):
    """
    Exact rational value of the Fock weight coupling^2.

    coupling = (1/2) * a_{k,k+1} * sqrt(h_{k+1}/h_k)
    coupling^2 = (1/4) * a_{k,k+1}^2 * h_{k+1}/h_k

    where k = n-l-1, lam = l+1,
    a_{k,k+1} = (k+1)/(2(k+lam)) = (n-l)/(2n),
    h_{k+1}/h_k = (k+2*lam)*(k+1) / ((k+1+lam)*(k+lam)) [from Gamma ratio formula]
                = (k+2l+2)(k+1) / ((k+l+2)(k+l+1))
                = (n+l+1)(n-l) / ((n+1)*n)

    So:
    coupling^2 = (1/4) * [(n-l)/(2n)]^2 * (n+l+1)(n-l) / ((n+1)*n)
               = (n-l)^3 * (n+l+1) / (16 * n^3 * (n+1))
    """
    k = n - l - 1
    lam = l + 1

    # Recursion coefficient squared
    a_sq = Rational(k + 1, 2 * (k + lam))**2

    # Norm ratio h_{k+1} / h_k
    # h_k = pi * 2^{1-2lam} * Gamma(k+2lam) / (k! * (k+lam) * Gamma(lam)^2)
    # h_{k+1}/h_k = Gamma(k+1+2lam)/Gamma(k+2lam) * k!/(k+1)! * (k+lam)/(k+1+lam)
    #             = (k+2lam) * 1/(k+1) * (k+lam)/(k+1+lam)
    # Wait, Gamma(k+1+2lam)/Gamma(k+2lam) = k+2lam (when incrementing by 1)
    h_ratio = Rational(k + 2*lam, 1) / Rational(k + 1, 1) * Rational(k + lam, k + 1 + lam)

    # Actually let me be more careful:
    # h_k = pi * 2^{1-2lam} * Gamma(k+2lam) / (k! * (k+lam) * Gamma(lam)^2)
    # h_{k+1} = pi * 2^{1-2lam} * Gamma(k+1+2lam) / ((k+1)! * (k+1+lam) * Gamma(lam)^2)
    #
    # h_{k+1}/h_k = Gamma(k+1+2lam)/Gamma(k+2lam) * k!/(k+1)! * (k+lam)/(k+1+lam)
    #             = (k+2*lam) * 1/(k+1) * (k+lam)/(k+1+lam)

    h_ratio_exact = Rational(k + 2*lam, k + 1) * Rational(k + lam, k + 1 + lam)

    csq = Rational(1, 4) * a_sq * h_ratio_exact

    # Substituting k = n-l-1, lam = l+1:
    # a = (k+1)/(2(k+lam)) = (n-l)/(2n)
    # h_ratio = (k+2lam)/(k+1) * (k+lam)/(k+1+lam)
    #         = (n-l-1+2l+2)/(n-l) * (n-l-1+l+1)/(n-l-1+l+1+1)
    #         = (n+l+1)/(n-l) * n/(n+1)

    csq_formula = Rational(1, 4) * (Rational(n-l, 2*n))**2 * Rational(n+l+1, n-l) * Rational(n, n+1)
    csq_simplified = simplify(csq_formula)

    assert simplify(csq - csq_simplified) == 0, f"Mismatch at n={n}, l={l}"

    return csq_simplified


def main():
    print("=" * 80)
    print("Exact rational coupling^2 values and closed-form formula")
    print("=" * 80)
    print()
    print("Formula: coupling^2(n,l) = (n-l)^2 * (n+l+1) / (16 * n^2 * (n+1))")
    print()
    print(f"{'(n,l)->(n+1,l)':<20} {'coupling^2':<20} {'= p/q':<15} {'16*n^2*(n+1)*c^2':<20} {'check (n-l)^2*(n+l+1)'}")
    print("-" * 95)

    results = {}

    for n in range(1, 10):
        for l in range(n):
            csq = coupling_squared_exact(n, l)
            factor = 16 * n**2 * (n + 1)
            numerator = csq * factor

            label = f"({n},{l})->({n+1},{l})"
            check = (n-l)**2 * (n+l+1)

            print(f"{label:<20} {str(csq):<20} = {float(csq):.8f}  "
                  f"{str(simplify(numerator)):<20} {check}")

            results[label] = {
                "coupling_squared": str(csq),
                "float": float(csq),
                "numerator": int(simplify(numerator)),
                "formula_check": int(check),
            }

    # Now check if the AVERAGE over all l at fixed n gives 1/16
    print("\n" + "=" * 80)
    print("Per-shell (2l+1)-weighted average of coupling^2")
    print("=" * 80)
    print()
    print("If the average equals 1/16, then kappa = -1/16 is the")
    print("degeneracy-weighted average of the Fock coupling^2.")
    print()
    print(f"{'n':<5} {'avg coupling^2':<20} {'ratio to 1/16':<15}")
    print("-" * 40)

    for n in range(1, 15):
        total = Rational(0)
        total_weight = 0
        for l in range(n):
            csq = coupling_squared_exact(n, l)
            weight = 2*l + 1
            total += weight * csq
            total_weight += weight

        avg = simplify(total / total_weight)
        ratio = simplify(avg / Rational(1, 16))

        print(f"{n:<5} {str(avg):<20} {str(ratio):<15}")
        results[f"n_{n}_avg"] = {"avg_csq": str(avg), "ratio_to_kappa": str(ratio)}

    # Check: is the formula for the average
    # <coupling^2>_n = (1/n^2) sum_{l=0}^{n-1} (2l+1) * (n-l)^2 * (n+l+1) / (16 n^2 (n+1))
    # = 1/(16 n^4 (n+1)) * sum_{l=0}^{n-1} (2l+1)(n-l)^2(n+l+1)

    print("\n" + "=" * 80)
    print("Closed-form for the sum S(n) = sum_{l=0}^{n-1} (2l+1)(n-l)^2(n+l+1)")
    print("=" * 80)

    n_sym = Symbol('n', positive=True, integer=True)
    l_sym = Symbol('l', positive=True, integer=True)

    for n in range(1, 12):
        S = sum((2*l+1) * (n-l)**2 * (n+l+1) for l in range(n))
        avg = Rational(S, 16 * n**4 * (n+1))
        print(f"  n={n}: S(n) = {S}, avg = S/(16*n^4*(n+1)) = {avg} = {float(avg):.8f}")

    # Check if S(n) = n^4 * (n+1) for all n (which would give avg = 1/16)
    print("\nCheck: is S(n) = n^4 * (n+1)?")
    for n in range(1, 12):
        S = sum((2*l+1) * (n-l)**2 * (n+l+1) for l in range(n))
        target = n**4 * (n+1)
        print(f"  n={n}: S={S}, n^4*(n+1)={target}, ratio = {S/target:.8f}")

    # The ratio is NOT 1 for n > 1, so the average is NOT 1/16.
    # Let me find the actual pattern.

    print("\nLet me find S(n) / n^2:")
    for n in range(1, 12):
        S = sum((2*l+1) * (n-l)**2 * (n+l+1) for l in range(n))
        print(f"  n={n}: S(n) = {S}, S/n^2 = {Rational(S, n**2)}, S/n^3 = {Rational(S, n**3)}, S/n^4 = {Rational(S, n**4)}")

    # Let me try polynomial fit for S(n)
    print("\nPolynomial structure of S(n):")
    sn_values = []
    for n in range(1, 12):
        S = sum((2*l+1) * (n-l)**2 * (n+l+1) for l in range(n))
        sn_values.append(S)
        # Check: S(n) = a*n^5 + b*n^4 + c*n^3 + ...
    print(f"  S(n) values: {sn_values}")

    # Use sympy to find the closed form
    n = Symbol('n')
    # S(n) = sum_{l=0}^{n-1} (2l+1)(n-l)^2(n+l+1)
    # Expand: (2l+1)(n-l)^2(n+l+1) = (2l+1)(n^2 - 2nl + l^2)(n+l+1)
    # This is a polynomial in l of degree 4. The sum of l^k from 0 to n-1
    # gives known Faulhaber formulas.

    # Let me compute symbolically
    l = Symbol('l')
    expr = (2*l + 1) * (n - l)**2 * (n + l + 1)
    expr_expanded = sp.expand(expr)
    print(f"\n  Expanded integrand: {expr_expanded}")

    # Compute the sum using known formulas
    S_closed = sp.summation(expr_expanded, (l, 0, n - 1))
    S_closed = sp.simplify(S_closed)
    S_closed = sp.factor(S_closed)
    print(f"  Closed form: S(n) = {S_closed}")

    # Verify
    for n_val in range(1, 8):
        computed = int(S_closed.subs(n, n_val))
        direct = sum((2*ll+1) * (n_val-ll)**2 * (n_val+ll+1) for ll in range(n_val))
        assert computed == direct, f"Mismatch at n={n_val}: {computed} vs {direct}"
    print(f"  (Verified for n=1..7)")

    # The average coupling^2
    avg_formula = sp.simplify(S_closed / (16 * n**2 * n * (n + 1)))  # weight is n^2 = sum (2l+1)
    avg_formula = sp.factor(avg_formula)
    print(f"\n  Average coupling^2 = S(n) / (16 * n^2 * (n+1)) = {avg_formula}")

    # Check limit as n -> infinity
    limit = sp.limit(avg_formula, n, sp.oo)
    print(f"  Limit as n -> infinity: {limit}")

    print()
    print("=" * 80)
    print("THE EXACT RESULT")
    print("=" * 80)
    print()
    print(f"  S(n) = sum_{{l=0}}^{{n-1}} (2l+1)(n-l)^2(n+l+1) = {S_closed}")
    print(f"  <coupling^2>_n = S(n) / (16 * n^2 * (n+1)) = {avg_formula}")
    print(f"  Limit as n->inf: {limit}")
    print()

    # Factor the relationship
    if limit == Rational(1, 16):
        print("  *** THE DEGENERACY-WEIGHTED AVERAGE OF COUPLING^2 APPROACHES 1/16 AS n->inf ***")
        print("  This means kappa = -1/16 is the ASYMPTOTIC average Fock weight coupling.")
    else:
        print(f"  The limit is {limit}, NOT 1/16 = {Rational(1,16)}")
        print(f"  Ratio limit/(1/16) = {simplify(limit * 16)}")

    # Save
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)
    json_path = os.path.join(data_dir, "probe_k1_fock_weight.json")

    results["closed_form_S"] = str(S_closed)
    results["avg_formula"] = str(avg_formula)
    results["limit_n_inf"] = str(limit)

    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {json_path}")

    return results


if __name__ == "__main__":
    main()
