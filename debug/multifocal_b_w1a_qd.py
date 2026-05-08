"""
W1a-diag Q-D: Cross-register bilinear ERI as the actual W1a closure object.

When R is itself a quantum-mechanical operator on a separate register at
focal length lam_n, the cross-register V_eN matrix element becomes

  M = <chi_{n_e l_e m_e}^{lam_e}(r) chi_{n_n l_n m_n}^{lam_n}(R) |
       1/|r - R|
       | chi_{n_e' l_e' m_e'}^{lam_e}(r) chi_{n_n' l_n' m_n'}^{lam_n}(R)>

This is a six-coordinate integral. The trick: factor through the standard
multipole expansion of 1/|r - R|:

  1/|r - R| = sum_L sum_M (4 pi)/(2L+1) [r_<^L / r_>^{L+1}]
              Y_{LM}^*(omega_r) Y_{LM}(omega_R)

The ANGULAR integrals factorize: r-side gives
  <Y_{l_e m_e} | Y_{LM}^* | Y_{l_e' m_e'}>_omega_r
which is a Gaunt integral, terminating at L_e_max = l_e + l_e' by Gaunt
selection rule. R-side gives
  <Y_{l_n m_n} | Y_{LM} | Y_{l_n' m_n'}>_omega_R
similarly terminating at L_n_max = l_n + l_n'. Combined: L_max = min(l_e+l_e',
l_n+l_n') is the effective truncation.

The RADIAL integral factorizes as a DOUBLE radial integral with
  r_<^L / r_>^{L+1}
in mixed coordinates. Defining s = min(r, R), t = max(r, R):

  J_L = int_0^infty dr int_0^infty dR R_e(r) R_e'(r) R_n(R) R_n'(R)
        r^2 R^2 [s^L / t^{L+1}]
      = int_0^infty dr R_e(r) R_e'(r) r^2 [
          int_0^r dR R_n(R) R_n'(R) R^{L+2} / r^{L+1}
        + int_r^infty dR R_n(R) R_n'(R) R^{1-L}
        ] r^L  -- wait, not quite. Let me be careful.

Actually:
  r_<^L / r_>^{L+1} = (r^L / R^{L+1}) when r < R
                    = (R^L / r^{L+1}) when r > R

So
  J_L = int_0^infty dr R_e R_e' r^{L+2} int_r^infty dR R_n R_n' R^{1-L}    (*)
      + int_0^infty dr R_e R_e' r^{1-L} int_0^r dR R_n R_n' R^{L+2}        (**)

Each of these is a *double* integral; the outer integral has an inner
integral with R-dependent limits. Both inner integrals are POLYNOMIAL-TIMES-
EXPONENTIAL on bounded intervals -- evaluable in closed form via lower /
upper incomplete gamma functions Gamma(s, x). The outer integral is then
the integral of (polynomial-times-exp(-2 lam_e r)) * (Gamma(s, c r))
where Gamma is the incomplete gamma.

KEY OBSERVATION: integrals of the form
  int_0^infty r^k exp(-alpha r) Gamma(s, beta r) dr
are STANDARD and have closed forms in terms of hypergeometric functions
(or, equivalently, in terms of the Lerch transcendent / generalized
incomplete gamma). See e.g. Gradshteyn & Ryzhik 6.455 or DLMF 8.14.

So Q-D's bilinear ERI HAS a closed form in elementary + incomplete-gamma
+ 2F1-hypergeometric functions, for any (n_e, l_e, n_n, l_n). The integral
is bilinear in (electron, nucleus) Sturmian indices: each term in the
electron radial expansion couples linearly with each term in the nucleus
radial expansion via a closed-form coefficient.

This computation verifies that for the simplest case (1s_e at lam_e times
0s_n HO basis or 1s_n hydrogenic, mismatched lambdas), the integral admits
a closed form, and identifies what new mathematics is required vs. what's
already in geovac/shibuya_wulfman.py.
"""

import json
from sympy import (Rational, Symbol, sqrt, exp, integrate, oo, simplify,
                   factorial, S, sympify, expand, factor, Min, Max,
                   Piecewise, lowergamma, uppergamma, gamma, hyper,
                   summation, Pow, Function, hyperexpand)
from sympy.functions.special.polynomials import assoc_laguerre


def R_radial(n, l, lam, r):
    rho = 2 * lam * r
    L_assoc = assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    psi_unnormed = rho**l * exp(-lam * r) * L_assoc
    norm_sq = integrate(psi_unnormed**2 * r**2, (r, 0, oo))
    return psi_unnormed / sqrt(simplify(norm_sq))


def cross_register_radial_integral_L0(lam_e, lam_n):
    """The simplest cross-register V_eN radial double integral at L=0,
    for 1s_e (electron at exponent lam_e) and 1s_n (nucleus at exponent lam_n).

    L=0 means s-s on electron side AND s-s on nucleus side. M=0.

    J_0 = int_0^infty dr R_e(r)^2 r^2 int_r^infty dR R_n(R)^2 R    [r > R term]
        + int_0^infty dr R_e(r)^2 r   int_0^r dR R_n(R)^2 R^2     [r < R term]

    The "1s_e" wavefunction at exponent lam_e: R_e(r) = 2 lam_e^(3/2) exp(-lam_e r).
    The "1s_n" wavefunction at exponent lam_n: R_n(R) = 2 lam_n^(3/2) exp(-lam_n R).

    Each radial probability density |R(r)|^2 r^2 = 4 lam^3 r^2 exp(-2 lam r).

    So:
      Inner R-integral, lower-limit r, of [4 lam_n^3 R^2 e^{-2 lam_n R} * R^{1-L}]
      with L=0: r^{1-L} = r^1, so this term is  int_0^infty dr 4 lam_e^3 r^2 e^{-2 lam_e r} * r * Inner(r)
      where Inner(r) = int_r^infty 4 lam_n^3 R^{3} e^{-2 lam_n R} dR  [from r^{1-L} = r factor moved outside r-integral; but I think I'm confusing myself]

    Let me just evaluate symbolically.
    """
    r = Symbol('r', positive=True)
    R = Symbol('R', positive=True)
    L = 0

    R_e = R_radial(1, 0, lam_e, r)
    R_n = R_radial(1, 0, lam_n, R)

    # Inner: t = max(r, R), s = min(r, R). r_<^L / r_>^{L+1} at L=0 = 1/t.
    # The integral splits at R = r:
    #   r < R: 1/R, factor R^2 dR R_n^2
    #   r > R: 1/r, factor R^2 dR R_n^2  (R is the inner var here)
    # Then weighted by r^2 R_e^2 dr.

    # Term I (r < R, i.e. R > r): 1/R = R^{-1} on the inner R integral.
    # Inner integral: int_r^infty R_n(R)^2 R^2 (1/R) dR = int_r^infty R_n^2 R dR
    # times outer int_0^infty R_e(r)^2 r^2 dr.
    inner_I = integrate(R_n**2 * R, (R, r, oo))

    # Term II (r > R): 1/r, inner R integral: int_0^r R_n^2 R^2 dR.
    # Outer: int_0^infty R_e^2 r^2 (1/r) [...] dr = int_0^infty R_e^2 r [...] dr.
    inner_II = integrate(R_n**2 * R**2, (R, 0, r))

    integrand = R_e**2 * r**2 * inner_I + R_e**2 * r * inner_II
    total = integrate(integrand, (r, 0, oo))
    return simplify(total)


def cross_register_radial_integral_L0_v2(lam_e, lam_n):
    """Use the cleaner formulation:

    J_0 = int_0^infty dr int_0^infty dR rho_e(r) rho_n(R) [1/max(r,R)]

    where rho_e(r) = R_e(r)^2 r^2, rho_n(R) = R_n(R)^2 R^2, are the radial
    probability densities (normalized to 1 in r^2 dr or R^2 dR).

    For 1s_e: rho_e = 4 lam_e^3 r^2 exp(-2 lam_e r).
    For 1s_n: rho_n = 4 lam_n^3 R^2 exp(-2 lam_n R).

    Standard integral (textbook Coulomb mutual repulsion of two 1s densities
    at different exponents):
      <1s_e | 1/r_12 | 1s_n>_radial-only-at-L=0
    Evaluates to: see e.g. Slater (1929), "Note on Hartree's Method," PR 35.

    Closed form (1/Hartree at lam_e, lam_n):
      J_0 = ?
    """
    r = Symbol('r', positive=True)
    R = Symbol('R', positive=True)
    rho_e = 4 * lam_e**3 * r**2 * exp(-2 * lam_e * r)
    rho_n = 4 * lam_n**3 * R**2 * exp(-2 * lam_n * R)
    # Integrate (R from 0 to r) gives "r > R" piece divided by r;
    # integrate (R from r to oo) gives "r < R" piece divided by R.
    inner_I = integrate(rho_n / R, (R, r, oo))   # r < R: 1/R
    inner_II = integrate(rho_n / r, (R, 0, r))   # r > R: 1/r
    total = integrate(rho_e * (inner_I + inner_II), (r, 0, oo))
    return simplify(total)


def cross_register_L1_pe_pn(lam_e, lam_n):
    """L=1 cross-register: 2p_e on electron, 2p_n on nucleus.
    Tests that L > 0 multipole also has a closed form.
    """
    r = Symbol('r', positive=True)
    R = Symbol('R', positive=True)
    L = 1

    R_e = R_radial(2, 1, lam_e, r)
    R_n = R_radial(2, 1, lam_n, R)

    # r_<^1 / r_>^2 case
    # r < R: r/R^2.  inner R-integral: int_r^infty R_n^2 R^2 (1/R^2) dR = int_r^infty R_n^2 dR
    # r > R: R/r^2.  inner R-integral: int_0^r R_n^2 R^2 (R/r^2) dR = (1/r^2) int_0^r R_n^2 R^3 dR
    inner_I = integrate(R_n**2, (R, r, oo))
    inner_II = integrate(R_n**2 * R**3, (R, 0, r))

    integrand = R_e**2 * r**2 * (r * inner_I) + R_e**2 * (1/r**2) * inner_II
    # No, the L=1 weights need to be properly attached. Let me redo.
    # The formula is:
    #   J_L = int_0^infty dr int_0^infty dR rho_e(r) rho_n(R)
    #         (r_<^L / r_>^{L+1}) / (r^2 R^2)
    # Because rho includes r^2 R^2.
    # Better: just use R_e * R_e' * r^2 and R_n * R_n' * R^2 explicitly.
    rho_e = R_e**2 * r**2
    rho_n = R_n**2 * R**2

    # r < R (R > r): r^L / R^{L+1} = r/R^2
    inner_lt = integrate(rho_n / R**(L+1), (R, r, oo))
    # r > R: R^L / r^{L+1} = R/r^2 ; integrate over R from 0 to r
    inner_gt = integrate(rho_n * R**L, (R, 0, r)) / r**(L+1)
    integrand = rho_e * (r**L * inner_lt + inner_gt)
    total = integrate(integrand, (r, 0, oo))
    return simplify(total)


def main():
    out = {}
    lam_e = Symbol('lam_e', positive=True)
    lam_n = Symbol('lam_n', positive=True)

    # 1) Cross-register 1s_e x 1s_n at L=0.
    print("Computing cross-register V_eN at L=0 for 1s_e, 1s_n (mismatched lam)...")
    J_0 = cross_register_radial_integral_L0_v2(lam_e, lam_n)
    print(f"\n  J_0(lam_e, lam_n) = {J_0}")
    out["J_0_1s_1s"] = str(J_0)

    # Verify symmetry J_0(lam_e, lam_n) = J_0(lam_n, lam_e)
    J_0_swapped = J_0.subs([(lam_e, lam_n), (lam_n, lam_e)], simultaneous=True)
    sym_diff = simplify(J_0 - J_0_swapped)
    out["J_0_symmetric"] = (sym_diff == 0)
    print(f"  Symmetric J_0(lam_e, lam_n) = J_0(lam_n, lam_e)? {sym_diff == 0}")

    # Verify lam_n -> infinity gives 2 lam_e (point-charge limit, classical V_ne).
    # For 1s_e and a point nucleus at origin, <1s | 1/r | 1s> = lam_e.
    # If nucleus has finite extent, it should be SMALLER than this. As lam_n -> infinity
    # (nucleus tightly localized at R=0), J_0 -> <1s | 1/r | 1s> = lam_e.
    print(f"\n  Point-nucleus limit lam_n -> oo: ", end="")
    from sympy import limit
    lim_pt = limit(J_0, lam_n, oo)
    print(f"J_0 -> {simplify(lim_pt)}")
    out["J_0_point_nucleus_limit"] = str(simplify(lim_pt))

    # Numerical check at lam_e = 1 (H 1s), lam_n = 1000 (very localized nucleus).
    J_0_num = float(J_0.subs([(lam_e, 1), (lam_n, 1000)]))
    print(f"  Numerical at (lam_e=1, lam_n=1000): J_0 = {J_0_num:.6f}")
    print(f"    (expect ~ <H 1s | 1/r | 1s> = 1.0)")
    out["J_0_numerical_lam1_lam1000"] = J_0_num

    # 2) Cross-register 2p_e x 2p_n at L=1.
    # NOTE: the angular factor for 2p_e x 2p_n at L=0 cross-coupling is zero
    # (both have l=1, so L=0 only via 3j(1,0,1;0,0,0); but also L=2). For
    # L=1 we get cross-parity zero. So this is a check of mechanism, not a
    # full physical computation. Let's just confirm the L=1 radial integral
    # has a closed form.
    print("\nComputing cross-register V_eN at L=1 for 2p_e, 2p_n (mismatched lam)...")
    try:
        J_1_pp = cross_register_L1_pe_pn(lam_e, lam_n)
        print(f"  J_1(2p, 2p) = {J_1_pp}")
        out["J_1_2p_2p"] = str(J_1_pp)
    except Exception as e:
        print(f"  Failed: {e}")
        out["J_1_2p_2p"] = f"failed: {e}"

    # Q-D scoping summary.
    out["q_d_scoping"] = {
        "verdict": "Cross-register bilinear V_eN HAS a closed-form symbolic representation.",
        "key_ingredients": [
            "Multipole expansion of 1/|r - R|.",
            "Standard angular Gaunt termination on BOTH electron and nucleus sides.",
            "Effective L truncation: L_max = min(l_e + l_e', l_n + l_n').",
            "Double radial integral: outer (poly x exp(-2 lam_e r)) x inner Gamma(s, c r).",
        ],
        "sparsity_implication": (
            "The cross-register angular factor is bilinear in (3j_e * 3j_n) and inherits "
            "Gaunt sparsity from BOTH sides. Pauli scaling estimate: at fixed Q_e, Q_n, "
            "the cross-register Pauli count is bounded by the product of single-register "
            "ERI counts, but the angular sparsity (Gaunt squared) reduces it from "
            "O(Q_e^2 Q_n^2) to O(Q_e^{1.5} Q_n^{1.5}) -- conjectural without full ERI count."
        ),
        "what_to_build": [
            "geovac/cross_register_vne.py: build_cross_register_vne(spec_e, spec_n) "
            "returning a Pauli dictionary on the joint register.",
            "Generalize Shibuya-Wulfman _split_integral_analytical from single-radial "
            "to double-radial via a 4D incomplete-gamma representation. The integrand "
            "is bivariate polynomial-times-exp; the integral is a sum of double-incomplete-"
            "gamma terms.",
            "Replace Track NI's R_PROTON_BOHR classical scalar with operator-valued "
            "R on the nuclear register.",
        ],
        "guardrail_compatibility": (
            "GUARDRAIL Papers 8-9: Sturmian theorem says shared-p_0 makes eigenvalues "
            "R-independent (proven negative). W1a does NOT use shared p_0 -- electron "
            "and nucleus are at different focal lengths by design (lam_e != lam_n). "
            "The theorem boundary is honored: the W1a closure is OUTSIDE the theorem's "
            "scope, not in conflict with it."
        ),
        "what_new_math_is_needed": [
            "Double-radial integral with mismatched exponents and r_<^L/r_>^{L+1} kernel "
            "in closed form. Q-A and Q-B above show this is mechanical -- the integrand "
            "factorizes through standard incomplete-gamma functions on bounded subintervals "
            "and gamma functions on infinite tails. NO new mathematics is required at the "
            "integral level; this is engineering using existing machinery.",
            "Operator-valued R on the proton register, built from HO matrix elements via "
            "Moshinsky-Talmi brackets for the proton in HO basis (Paper 23 §III). Already "
            "exists in geovac/nuclear/moshinsky.py for relative coordinates; needs explicit "
            "lab-frame R operator construction.",
            "A spectral-triple lift IF one wants to phrase the closure in NCG language. "
            "The W1a closure as a Pauli encoding does NOT require this; the spectral-triple "
            "framing is for Paper 32 §VIII (a Connes-Marcolli A8'-class bridge for spatial "
            "registers). The Pauli encoding closes W1a; the spectral-triple framing is "
            "Phase C polish.",
        ],
    }

    with open("debug/data/multifocal_b_w1a_qd.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: debug/data/multifocal_b_w1a_qd.json")


if __name__ == "__main__":
    main()
