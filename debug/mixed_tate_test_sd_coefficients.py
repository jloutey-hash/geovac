"""Mixed-Tate-period test for GeoVac discrete Seeley-DeWitt coefficients on S^3.

Sprint: Mixed-Tate Test (2026-06-03).

Tests whether the discrete Seeley-DeWitt (SD) coefficients of the squared
Dirac operator D^2 and the scalar Laplacian Delta on unit S^3 sit inside
the mixed-Tate-period ring of Fathizadeh-Marcolli arXiv:1611.01815, or
break out of it.

Two normalizations of the heat-trace expansion are tracked separately:

1. Raw heat-trace coefficients tilde_a_k defined by
       K(t) = sum_{k>=0} tilde_a_k * t^{k - d/2}
   For d=3 on S^3, the leading Dirac coefficients are
       tilde_a_0^D2 = sqrt(pi)/2,  tilde_a_1^D2 = -sqrt(pi)/4.
   These contain sqrt(pi), which is NOT a mixed-Tate period over Q.

2. Standard CC volume-normalized coefficients a_k = (4 pi)^{d/2} * tilde_a_k.
   For d=3 the (4 pi)^{3/2} prefactor exactly cancels the sqrt(pi), giving
       a_0^D2 = 4 pi^2,  a_1^D2 = -2 pi^2,  a_k^D2 = 0 for k>=2.
   For the scalar Laplacian, a_k^Delta = 2 pi^2 / k! for all k>=0.
   All of these lie in Q[pi^2] subset Q[(2 pi i)^{+/-1}] subset
   mixed-Tate-period ring over Q.

This is the convention used by Fathizadeh-Marcolli and by Connes-Chamseddine
when classifying the spectral action expansion coefficients.

Conclusion: the GeoVac discrete SD sector is VOLUME-NORMALIZATION DEPENDENT
on the mixed-Tate question, with the standard convention (matching F-M)
landing POSITIVE.
"""
from __future__ import annotations

import sympy as sp
from sympy import (
    Integer, Rational, factorial, oo, pi, simplify, sqrt, summation, symbols, zeta,
)


def dirac_zeta_D2_at_integer_s(s_int):
    """zeta_{D^2}(s) at integer s via direct degeneracy-weighted sum.

    Uses Camporesi-Higuchi spectrum |lambda_n| = n + 3/2 with degeneracy
    g_n = 2(n+1)(n+2) on unit S^3 (Paper 51 eq:CH_spec).
    """
    n = symbols("n", integer=True, nonnegative=True)
    return summation(2 * (n + 1) * (n + 2) * (n + Rational(3, 2))**(-2 * s_int), (n, 0, oo))


def dirac_raw_heat_trace_coefficients():
    """Raw heat-trace coefficients tilde_a_k for Dirac on unit S^3.

    Paper 51 Corollary 2.1 (two-term exactness):
        K_{D^2}(t) = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}
                     + O(e^{-pi^2/t}).
    """
    return {0: sqrt(pi) / 2, 1: -sqrt(pi) / 4}


def dirac_volume_normalized_SD():
    """Standard CC volume-normalized SD coefficients a_k = (4 pi)^{3/2} tilde_a_k."""
    pref = (4 * pi)**Rational(3, 2)
    return {k: simplify(v * pref) for k, v in dirac_raw_heat_trace_coefficients().items()}


def scalar_volume_normalized_SD(k_max=6):
    """Scalar Laplacian SD coefficients a_k^Delta = 2 pi^2 / k! on unit S^3.

    Paper 51 Theorem 3.1.
    """
    return {k: 2 * pi**2 / factorial(k) for k in range(k_max + 1)}


def classify_mixed_tate(expr):
    """Best-effort classifier for whether expr lies in mixed-Tate-period ring over Q.

    Returns one of: 'rational', 'Q[pi^2]', 'Q[pi]', 'Q[sqrt(pi)]', 'unknown'.
    """
    expr = sp.expand(expr)
    if expr.is_rational:
        return "rational (subset Q ⊂ mixed-Tate over Q)"
    # Substitute pi -> x and check whether x appears only with even integer exponents.
    x = symbols("x", positive=True)
    sub = expr.subs(pi, x)
    poly = sp.Poly(sub, x)
    powers = poly.monoms()
    flat = [p[0] for p in powers]
    if all(isinstance(p, int) and p % 2 == 0 for p in flat):
        return "Q[pi^2] (subset mixed-Tate over Q)"
    if all(isinstance(p, sp.Rational) and p.q == 1 for p in flat):
        return "Q[pi] (subset mixed-Tate over Q)"
    if any(isinstance(p, sp.Rational) and p.q == 2 for p in flat):
        return "Q[sqrt(pi)] (NOT mixed-Tate over Q; algebraic extension by sqrt(pi))"
    return "unknown"


if __name__ == "__main__":
    print("=" * 72)
    print("Mixed-Tate Test - Sprint canonical computation (2026-06-03)")
    print("=" * 72)
    print()

    print("(A) zeta_{D^2}(s) at integer s (Paper 28 T9 verification)")
    print("-" * 72)
    for s_int in range(2, 6):
        val = dirac_zeta_D2_at_integer_s(s_int)
        print(f"  zeta_D2({s_int}) = {simplify(val)}")
        print(f"    classification: {classify_mixed_tate(val)}")
    print()

    print("(B) Raw heat-trace coefficients K(t) = sum_k tilde_a_k t^{k-3/2}")
    print("-" * 72)
    raw = dirac_raw_heat_trace_coefficients()
    for k, v in raw.items():
        print(f"  tilde_a_{k}^D2 = {v}")
        print(f"    classification: {classify_mixed_tate(v**2)}  (squared, to drop sqrt)")
        # Note: classify works on polys in pi; sqrt(pi) shows as half-integer power
        if sqrt(pi) in v.atoms(sp.Pow) or v.has(sqrt(pi)):
            print("    contains sqrt(pi) literally  ==>  NOT mixed-Tate at raw level")
    print("  tilde_a_k = 0 for k >= 2 (two-term exact, Paper 51 Cor 2.1)")
    print()

    print("(C) Volume-normalized SD coefficients (standard CC / F-M convention)")
    print("-" * 72)
    for k, v in dirac_volume_normalized_SD().items():
        print(f"  a_{k}^Dirac = {v}")
        print(f"    classification: {classify_mixed_tate(v)}")
    print("  a_k^Dirac = 0 for k >= 2")
    print()
    print("Scalar Laplacian SD coefficients a_k^Delta = 2 pi^2 / k!:")
    for k, v in scalar_volume_normalized_SD(5).items():
        print(f"  a_{k}^Delta = {v}")
        print(f"    classification: {classify_mixed_tate(v)}")
    print()

    print("(D) Verdict")
    print("-" * 72)
    print("Raw heat-trace coefficients: contain sqrt(pi) literally.")
    print("  ==> NOT mixed-Tate over Q at the raw level.")
    print()
    print("Volume-normalized SD coefficients (standard CC + F-M convention):")
    print("  Dirac: a_0 = 4 pi^2, a_1 = -2 pi^2 (then zero).")
    print("  Scalar: a_k = 2 pi^2 / k! for all k.")
    print("  ==> All in Q[pi^2] ⊂ mixed-Tate over Q.")
    print()
    print("Conclusion: the GeoVac discrete SD sector matches F-M's mixed-Tate")
    print("classification under the standard convention, with the additional")
    print("structural property of two-term exactness (Dirac) or geometric-series")
    print("collapse (scalar) that the F-M continuum R-W expansion does not have.")
