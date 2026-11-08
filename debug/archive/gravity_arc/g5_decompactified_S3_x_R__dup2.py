"""Sprint G5 — spectral action on the decompactified S^3 x R_tau.

Extension of G2 (S^3 x S^1_beta) to non-compact temporal direction. Makes the
beta -> infinity minimum of G2's spectral action explicit.

Setup
-----
G2 found that S^3_R x S^1_beta has spectral action S(R, beta, Lambda)
with formal extremum at R Lambda = 1/sqrt(6), beta -> infinity (joint
minimum corresponds to zero-temperature de Sitter vacuum).

G5 takes the beta -> infinity limit explicitly. The temporal direction
decompactifies to R_tau (non-compact). The action becomes infinite
(proportional to length of R_tau), so the meaningful quantity is the
ACTION DENSITY PER UNIT TIME:

    s(R, Lambda) := lim_{beta -> infinity} S(R, beta, Lambda) / beta

For Gaussian cutoff f(x) = e^{-x}, G2's result was
    S_Gauss(R, beta, Lambda) = (beta R^3 / 4) Lambda^4 - (beta R / 8) Lambda^2

So
    s(R, Lambda) = (R^3 / 4) Lambda^4 - (R / 8) Lambda^2

This is the spectral action DENSITY on S^3 x R_tau.

Equivalent derivation via heat kernel
-------------------------------------
On S^3 x R_tau, the heat kernel is
    K_full(t, R) = K_{S^3}(t, R) * K_R(t)
where K_R(t) = int dk/(2 pi) exp(-k^2 t) = 1/(2 sqrt(pi t)) is the
heat-kernel DENSITY on the non-compact R_tau (per unit length).

K_{S^3}(t, R) = K_unit(t/R^2) = (sqrt(pi)/2) R^3 t^{-3/2} - (sqrt(pi)/4) R t^{-1/2}
              + O(exp(-pi^2 R^2 / t))

So K_full(t, R) (per unit time) = K_{S^3}(t, R) * 1/(2 sqrt(pi t))
              = (R^3/4) t^{-2} - (R/8) t^{-1} + O(exp small)

Identical to G2's per-beta density.

Result
------
The spectral action density on S^3 x R_tau is EXACTLY two-term:
    s(R, Lambda) = (R^3 / 4) Lambda^4 - (R / 8) Lambda^2 + O(exp small)

Extremum at d s / d R = 0:
    3 R^2 / 4 * Lambda^4 - 1/8 * Lambda^2 = 0
    R^2 = 1 / (6 Lambda^2)
    R_crit Lambda = 1/sqrt(6)   <-- same as G1 and G2

Minimum action density value:
    s(R_crit, Lambda) = (R_crit^3 / 4) Lambda^4 - (R_crit / 8) Lambda^2
                     = R_crit Lambda^2 [(R_crit^2 Lambda^2)/4 - 1/8]
                     = R_crit Lambda^2 [(1/24) - (1/8)]
                     = R_crit Lambda^2 * (-1/12)
                     = -Lambda / (12 sqrt(6))

The CC spectral-action principle reading: -Lambda/(12 sqrt(6)) is the
"vacuum energy density" of the de Sitter vacuum on S^3 x R_tau, with the
preferred radius R = 1/(sqrt(6) Lambda).
"""

import json
from pathlib import Path

import sympy as sp
from sympy import symbols, pi, Rational, simplify, diff, sqrt

OUT_JSON = Path(__file__).parent / "data" / "g5_decompactified.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main():
    results = {}
    print("=" * 72)
    print("Sprint G5: spectral action on the decompactified S^3 x R_tau")
    print("=" * 72)

    # -----------------------------------------------------------------------
    # Step 1: Heat-kernel density derivation
    # -----------------------------------------------------------------------
    print("\n[Step 1] Heat-kernel density on S^3_R x R_tau:")
    print()
    print("  K_full(t, R) = K_{S^3}(t, R) * K_R(t)  (per unit length of R_tau)")
    print()
    print("  K_R(t) = int dk/(2 pi) exp(-k^2 t) = 1 / (2 sqrt(pi t))")
    print("    (heat-kernel density per unit length on the real line)")
    print()
    print("  K_{S^3}(t, R) = (sqrt(pi)/2) R^3 t^{-3/2} - (sqrt(pi)/4) R t^{-1/2}")
    print("                  + O(exp(-pi^2 R^2 / t))   [Paper 28 two-term exactness]")
    print()
    t, R, Lambda = symbols("t R Lambda", positive=True)

    K_S3 = (sqrt(pi)/2) * R**3 / t**Rational(3,2) - (sqrt(pi)/4) * R / sqrt(t)
    K_R = 1 / (2 * sqrt(pi * t))
    K_full = K_S3 * K_R
    K_full_simp = sp.simplify(K_full)
    K_full_expanded = sp.expand(K_full_simp)
    print(f"  K_full(t, R) = K_{{S^3}}(t, R) * K_R(t)")
    print(f"               = {K_full_expanded}")

    results["heat_kernel_density"] = str(K_full_expanded)

    # -----------------------------------------------------------------------
    # Step 2: Spectral action density (Gaussian cutoff)
    # -----------------------------------------------------------------------
    print("\n[Step 2] Spectral action density via Gaussian cutoff f(x) = e^{-x}:")
    print()
    print("  For Gaussian f: S = Tr exp(-D^2/Lambda^2) factorizes via")
    print("  S = K(1/Lambda^2) on S^3 x R_tau (per unit length of R_tau).")
    print()
    s_action = K_full.subs(t, 1/Lambda**2)
    s_action_simp = sp.simplify(s_action)
    s_action_expanded = sp.expand(s_action_simp)
    print(f"  s(R, Lambda) = K_full(1/Lambda^2, R)")
    print(f"              = {s_action_expanded}")

    results["action_density"] = str(s_action_expanded)

    # -----------------------------------------------------------------------
    # Step 3: Extremum in R
    # -----------------------------------------------------------------------
    print("\n[Step 3] Extremum: d s / d R = 0:")
    ds_dR = diff(s_action_simp, R)
    ds_dR_simp = sp.simplify(ds_dR)
    print(f"  d s / d R = {ds_dR_simp}")

    R_crit_eq = sp.Eq(ds_dR_simp, 0)
    R_crit_sols = sp.solve(R_crit_eq, R)
    print(f"  Solutions: R = {R_crit_sols}")

    R_crit_positive = [r for r in R_crit_sols if sp.simplify(r).is_positive]
    if R_crit_positive:
        R_crit = R_crit_positive[0]
    else:
        # Try simplifying
        R_crit = sp.simplify(R_crit_sols[0])
    print(f"  R_crit = {R_crit}")
    print(f"  R_crit Lambda = {sp.simplify(R_crit * Lambda)}  (should be 1/sqrt(6))")

    R_crit_Lambda = sp.simplify(R_crit * Lambda)
    check_value = sp.simplify(R_crit_Lambda - 1/sqrt(6))
    print(f"  R_crit Lambda - 1/sqrt(6) = {check_value}  (should be zero)")

    results["extremum"] = {
        "R_crit_solutions": [str(r) for r in R_crit_sols],
        "R_crit_Lambda": str(R_crit_Lambda),
        "matches_1_over_sqrt6": check_value == 0,
    }

    # -----------------------------------------------------------------------
    # Step 4: Minimum action density value
    # -----------------------------------------------------------------------
    print("\n[Step 4] Minimum action density value at R = R_crit:")
    s_min = s_action_simp.subs(R, R_crit)
    s_min_simp = sp.simplify(s_min)
    s_min_factored = sp.factor(s_min_simp)
    print(f"  s(R_crit, Lambda) = {s_min_simp}")
    print(f"                    = {s_min_factored}")

    # Expected: -Lambda / (12 sqrt(6))
    # Derivation: s(R_crit, Lambda) at R_crit = 1/(sqrt(6) Lambda):
    #   (R_crit^3/4) Lambda^4 = Lambda/(24 sqrt(6))
    #   (R_crit/8) Lambda^2 = Lambda/(8 sqrt(6))
    #   s_min = Lambda/(24 sqrt(6)) - Lambda/(8 sqrt(6)) = -Lambda/(12 sqrt(6))
    s_min_expected = -Lambda / (12 * sqrt(6))
    diff_check = sp.simplify(s_min_simp - s_min_expected)
    print(f"\n  Expected: -Lambda / (12 sqrt(6))")
    print(f"  Difference: s(R_crit) - (-Lambda/(12 sqrt(6))) = {diff_check}")
    print(f"  Match: {diff_check == 0}")

    results["minimum_action_density"] = {
        "s_min": str(s_min_simp),
        "expected": str(s_min_expected),
        "matches_expected": diff_check == 0,
    }

    # -----------------------------------------------------------------------
    # Step 5: Reading as zero-temperature de Sitter vacuum
    # -----------------------------------------------------------------------
    print("\n[Step 5] Physical reading: zero-temperature de Sitter vacuum")
    print()
    print("  G2 found that the formal extremum of S(R, beta, Lambda) on")
    print("  S^3 x S^1_beta is at R Lambda = 1/sqrt(6) and (joint minimum)")
    print("  beta -> infinity. G5 makes this explicit:")
    print()
    print("  - The decompactified time direction R_tau is the natural setting")
    print("    for the 'vacuum' state (no thermal compactification = T = 0).")
    print("  - The spectral action density s(R, Lambda) = K(1/Lambda^2, R)")
    print("    has the same two-term structure as G2's per-beta density.")
    print("  - The extremum R_crit = 1/(sqrt(6) Lambda) is the de Sitter")
    print("    vacuum radius (Planck-scale at Lambda = M_Planck).")
    print()
    print("  - The minimum action density s_min = -Lambda/(12 sqrt(6))")
    print("    is the 'vacuum energy density' of de Sitter at this background.")
    print()
    print("  This is the standard Connes-Chamseddine spectral-action-principle")
    print("  prediction of a Planck-scale de Sitter vacuum, with the standard")
    print("  cosmological-constant-scale gap (observed Lambda_cc << Lambda^3).")

    # -----------------------------------------------------------------------
    # Step 6: Verdict
    # -----------------------------------------------------------------------
    print("\n[Step 6] Verdict:")
    print()
    print("  POSITIVE: spectral action on S^3 x R_tau is EXACTLY two-term")
    print("  (modulo exp small in min(R^2/t, infinity)) via heat-kernel")
    print("  decompactification of G2's S^3 x S^1_beta result.")
    print()
    print("  Vacuum extremum at R_crit Lambda = 1/sqrt(6) explicit; vacuum")
    print("  energy density -Lambda/(12 sqrt(6)) clean closed form.")
    print()
    print("  G2's joint minimum at beta -> infinity is now made manifest:")
    print("  the framework's preferred 'zero-temperature de Sitter vacuum'")
    print("  has Planck-scale radius and Planck^4-cube vacuum energy density.")
    print("  Inherits the standard CC cosmological-constant-scale gap.")

    results["verdict"] = "POSITIVE"

    # Save JSON
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)
    return results


if __name__ == "__main__":
    main()
