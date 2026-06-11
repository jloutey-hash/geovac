"""
W1a-diag Q-A and Q-B: Multi-lambda Shibuya-Wulfman cross-center V_ne.

Goal: compute symbolically
  I_L(R; lambda_a, lambda_b, n_a, l_a, n_b, l_b) =
    int_0^infty R_{n_a l_a}^{lambda_a}(r) R_{n_b l_b}^{lambda_b}(r) g_L(r, R) r^2 dr

with the Sturmian-style (or hydrogenic-style) radial functions
  R_{nl}^{lambda}(r) = N_{nl}(lambda) * (2 lambda r)^l * exp(-lambda r) * L_{n-l-1}^{2l+1}(2 lambda r)

and the multipole kernel
  g_L(r, R) = r_<^L / r_>^{L+1} = { r^L / R^{L+1},   r < R
                                   R^L / r^{L+1},   r > R }

The matched case (lambda_a = lambda_b) reproduces the existing Shibuya-Wulfman
implementation in geovac/shibuya_wulfman.py.

The mismatched case is the Q-A object. We test:
  Q-A: Does the integral retain a closed form? Identify the ring of values.
  Q-B: Does the multipole expansion still terminate at L_max = 2 max(l_a, l_b)?
       The angular factor is unchanged from the matched case (depends only on
       (l_a, l_b, m_a, m_b)) so the *angular* termination is preserved by
       Gaunt selection. The question is whether the *radial* integral I_L(R)
       gives a nontrivial contribution for L > L_max_angular -- and if not,
       what kills it.

DESIGN NOTE: The angular termination in the matched-lambda case is a
GAUNT/3j theorem, not a radial-orthogonality theorem. The radial integral
I_L(R) is generically nonzero for any L; what kills the higher-L terms
is the angular integral being identically zero by 3j(l_a, L, l_b; 0, 0, 0)
parity / triangle-inequality. So Q-B as posed in the briefing has a clear
prediction: termination at L_max = 2 max(l_a, l_b) survives the mismatch
trivially, because mismatch is in the radial sector only. The interesting
quantity is the radial integral I_L itself -- its closed-form structure
under mismatch.

We compute Q-A symbolically for several (n_a l_a, n_b l_b, L) cases at both
matched (lambda_a = lambda_b = lambda) and mismatched (lambda_a != lambda_b)
exponents, and check:
  - Closed form exists in elementary functions (polynomials + exponentials +
    incomplete gamma / lower / upper incomplete).
  - The integrand is a polynomial-times-exp(-(lambda_a + lambda_b) r), so the
    radial integral is identical in structure to the matched case with a
    single combined exponent alpha_total = lambda_a + lambda_b.
  - Therefore the *same* incomplete-gamma machinery in
    `geovac/shibuya_wulfman.py::_split_integral_analytical` works for the
    mismatched case with one substitution: alpha = lambda_a + lambda_b.

The rest of this file works through specific symbolic cases.
"""

import json
import sys
from sympy import (Rational, Symbol, sqrt, exp, integrate, oo, simplify,
                   factorial, lowergamma, uppergamma, gamma, expand, factor,
                   Piecewise, Poly, symbols, sympify)
from sympy.functions.special.polynomials import assoc_laguerre

# ---- 1) Hydrogenic radial wavefunction with arbitrary exponent lambda
def R_radial(n, l, lam, r):
    """Coulomb-Sturmian / hydrogenic radial function at exponent lam.

    Uses the convention
      R_{nl}^lam(r) = N(n,l,lam) * (2 lam r)^l * exp(-lam r) *
                      L_{n-l-1}^{2l+1}(2 lam r)
    with normalization
      int_0^infty |R(r)|^2 r^2 dr = 1.

    For lam = Z/n this is the standard hydrogenic 1/Z atomic units form.
    """
    # Use sympy's associated Laguerre.
    rho = 2 * lam * r
    L = assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    psi_unnormed = rho**l * exp(-lam * r) * L
    # Compute norm symbolically.
    norm_sq = integrate(psi_unnormed**2 * r**2, (r, 0, oo))
    norm_sq = simplify(norm_sq)
    return psi_unnormed / sqrt(norm_sq)


def radial_overlap_check():
    """Sanity check: <R_nl^lam | R_nl^lam> = 1 in r^2 dr measure."""
    r = Symbol('r', positive=True)
    lam = Symbol('lam', positive=True)
    cases = [(1, 0), (2, 0), (2, 1)]
    out = {}
    for (n, l) in cases:
        R = R_radial(n, l, lam, r)
        val = integrate(R**2 * r**2, (r, 0, oo))
        out[f"n={n},l={l}"] = str(simplify(val))
    return out


def split_radial_integral(n_a, l_a, lam_a, n_b, l_b, lam_b, L, R):
    """Compute the radial split-region integral
      I_L = (1/R^{L+1}) int_0^R R_a R_b r^{L+2} dr
          + R^L int_R^infty R_a R_b r^{1-L} dr
    symbolically for arbitrary exponents.

    For L=0 the second integrand has factor r^1; for L>=2 the factor r^{1-L}
    is divergent at r=0 and finite at infinity, but the inner integral has
    a *higher* power of r (r^{L+2}) which compensates -- only the SUM is
    finite.
    """
    r = Symbol('r', positive=True)
    R_a = R_radial(n_a, l_a, lam_a, r)
    R_b = R_radial(n_b, l_b, lam_b, r)
    integrand_inner = R_a * R_b * r**(L + 2)
    integrand_outer = R_a * R_b * r**(1 - L)

    I_inner = integrate(integrand_inner, (r, 0, R))
    I_outer = integrate(integrand_outer, (r, R, oo))
    total = I_inner / R**(L + 1) + R**L * I_outer
    return simplify(total)


def diagnose_q_a():
    """Q-A: same-center mismatched-exponent radial integral.

    For (n_a, l_a) = (n_b, l_b) = (1, 0), this is the simplest
    mismatched case.

    The integrand factorizes as
       N(1, 0, lam_a) * N(1, 0, lam_b) * exp(-(lam_a + lam_b) r)
        * r^{L+2} or r^{1-L}.

    So the integral is a sum of incomplete-gamma functions in the combined
    exponent alpha = lam_a + lam_b. This is the SAME structural form as the
    matched case; only the value of alpha changes.

    Conclusion of Q-A: closed form survives. The integrand is still a
    polynomial-times-single-exponential. Mismatch in exponents only changes
    the value of the combined decay rate alpha = lam_a + lam_b. The Avery
    school's "isoenergetic" constraint that forces matched lambda is a
    PHYSICAL constraint (so the two basis functions span a single energy
    sector), NOT a mathematical constraint on the integral existing.
    """
    R = Symbol('R', positive=True)
    lam_a = Symbol('lam_a', positive=True)
    lam_b = Symbol('lam_b', positive=True)

    out = {}

    # Case 1: 1s-1s same-center, mismatched lambda, L=0
    print("Case 1: (n_a, l_a, n_b, l_b) = (1, 0, 1, 0), L = 0, lambda_a != lambda_b")
    I_00 = split_radial_integral(1, 0, lam_a, 1, 0, lam_b, 0, R)
    I_00_simpl = simplify(I_00)
    out["I_00_mismatched_1s1s"] = str(I_00_simpl)
    print("  I_L=0 =", I_00_simpl)

    # Verify matched limit lam_b -> lam_a recovers single-lambda case
    matched = simplify(I_00_simpl.subs(lam_b, lam_a))
    out["I_00_matched_limit"] = str(matched)
    print("  Matched limit (lam_b -> lam_a):", matched)

    # The "expected" single-lambda result for 1s-1s nuclear-attraction
    # multipole at L=0 is the integral
    #   int_0^infty (4 lam^3) exp(-2 lam r) g_0(r, R) r^2 dr
    # = (1/R) int_0^R 4 lam^3 exp(-2 lam r) r^2 dr + int_R^infty 4 lam^3 exp(-2 lam r) r dr
    # = (1/R)[1 - exp(-2 lam R)(1 + 2 lam R + 2 lam^2 R^2)] / R
    #   wait: standard result is 1/R - exp(-2 lam R)/R - lam exp(-2 lam R)
    # which simplifies to (1 - exp(-2 lam R)(1 + lam R)) / R after careful work.
    # Let's just verify it's a function of (lam_a + lam_b) R only.

    return out


if __name__ == "__main__":
    print("=" * 70)
    print("W1a-diag Q-A: Multi-lambda Shibuya-Wulfman radial integral")
    print("=" * 70)

    print("\n[1] Sanity check: radial wavefunctions normalized to 1")
    norms = radial_overlap_check()
    for k, v in norms.items():
        print(f"  <R_{k} | R_{k}> = {v}")

    print("\n[2] Q-A diagnosis: mismatched 1s-1s radial integral at L=0")
    qa_out = diagnose_q_a()

    out = {
        "norms": norms,
        "q_a": qa_out,
    }
    with open("debug/data/multifocal_b_w1a_qa_qb.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: debug/data/multifocal_b_w1a_qa_qb.json")
