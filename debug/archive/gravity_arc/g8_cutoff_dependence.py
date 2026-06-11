"""Sprint G8 — Cutoff-function dependence of G_eff and Lambda_cc.

Extends G7 (which used Gaussian cutoff only) to arbitrary cutoff functions.
Identifies which framework predictions are cutoff-independent (true GeoVac
predictions) and which depend on the cutoff (calibration data).

Setup
-----
From G2: spectral action on S^3_R x S^1_beta has the form
    S(R, beta, Lambda) = c_4 phi_4 Lambda^4 V + c_2 phi_2 Lambda^2 int R sqrt(g)
where:
- c_4, c_2 are fixed by GeoVac's spectral structure
- phi_4, phi_2 are MOMENTS of the cutoff function f

For Mellin transform phi(s) = int_0^infty f(x) x^{s-1} dx:
- phi(2) = "Lambda^4 moment" (cosmological-constant-scale)
- phi(1) = "Lambda^2 moment" (Einstein-Hilbert-scale)

Standard CC: G_eff comes from c_2 phi(1) coefficient; Lambda_cc from c_4 phi(2).

Question
--------
Which combinations of G_eff and Lambda_cc are cutoff-independent?

Answer (derived below):
- G_eff depends on phi(1) only
- Lambda_cc depends on phi(2)/phi(1) ratio
- The DIMENSIONLESS combination Lambda_cc * G_eff is cutoff-DEPENDENT
- The combination G_eff * Lambda^2 = 6 pi / phi(1) (independent of which sphere)
- The combination Lambda_cc / Lambda^2 = 6 phi(2) / phi(1) (depends on phi-moments)

Specific cutoffs
----------------
1. Gaussian f(x) = e^{-x}: phi(s) = Gamma(s); phi(1) = 1, phi(2) = 1
2. Polynomial f(x) = e^{-x^2}: phi(s) = (1/2) Gamma(s/2); phi(1) = sqrt(pi)/2, phi(2) = 1/2
3. Sharp f(x) = Theta(1 - x): phi(s) = 1/s; phi(1) = 1, phi(2) = 1/2

Comparison
----------
Verify cutoff-dependent quantities and identify universals.
"""

import json
from pathlib import Path

import sympy as sp
from sympy import (
    symbols, pi, Rational, simplify, sqrt, gamma,
    exp, integrate, oo, Symbol,
)

OUT_JSON = Path(__file__).parent / "data" / "g8_cutoff_dependence.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main():
    results = {}
    print("=" * 72)
    print("Sprint G8: cutoff-function dependence of G_eff and Lambda_cc")
    print("=" * 72)

    # Symbols
    Lambda, R, beta, phi1, phi2 = symbols("Lambda R beta phi_1 phi_2", positive=True)
    s = symbols("s")

    # -----------------------------------------------------------------------
    # Step 1: Derive G_eff and Lambda_cc in terms of phi-moments
    # -----------------------------------------------------------------------
    print("\n[Step 1] G2 spectral action with general cutoff f:")
    print()
    print("  S(R, beta, Lambda) = phi(2) * (beta R^3 / 4) * Lambda^4")
    print("                     + phi(1) * (-beta R / 8) * Lambda^2")
    print()
    print("  where phi(s) = int_0^infty f(x) x^{s-1} dx is the Mellin transform of f")
    print("  evaluated at s = 2 (Lambda^4 coefficient) and s = 1 (Lambda^2 coefficient).")
    print()
    print("  For Gaussian f(x) = e^{-x}: phi(s) = Gamma(s), so phi(1) = phi(2) = 1.")
    print("  G7's Gaussian-cutoff result: G_eff = 6 pi/Lambda^2, Lambda_cc = 6 Lambda^2")
    print()
    print("  Generalization to arbitrary f:")
    print("    Lambda^4 coefficient matches Lambda_cc/(8 pi G_eff) * Vol")
    print("    Lambda^2 coefficient matches -1/(16 pi G_eff) * int R_scalar sqrt(g)")
    print()
    print("  Working out the algebra (Vol = 2 pi^2 R^3 beta, int R = 12 pi^2 R beta):")

    # Matching Lambda^4 term:
    # phi(2) * beta R^3 / 4 * Lambda^4 = Lambda_cc/(8 pi G_eff) * 2 pi^2 R^3 beta
    # => Lambda_cc / G_eff = phi(2) Lambda^4 / pi
    print()
    print("    Lambda^4: phi(2) * (beta R^3 / 4) * Lambda^4 = Lambda_cc/(8 pi G_eff) * 2 pi^2 R^3 beta")
    print("           => Lambda_cc / G_eff = phi(2) Lambda^4 / pi")

    # Matching Lambda^2 term:
    # phi(1) * (-beta R / 8) * Lambda^2 = -1/(16 pi G_eff) * 12 pi^2 R beta
    # => 1/G_eff = 1/(6 pi) * phi(1) Lambda^2
    # => G_eff = 6 pi / (phi(1) Lambda^2)
    print()
    print("    Lambda^2: phi(1) * (-beta R / 8) * Lambda^2 = -1/(16 pi G_eff) * 12 pi^2 R beta")
    print("           => G_eff = 6 pi / (phi(1) Lambda^2)")

    # Solve for Lambda_cc
    # Lambda_cc = G_eff * phi(2) Lambda^4 / pi
    #           = (6 pi / (phi(1) Lambda^2)) * phi(2) Lambda^4 / pi
    #           = 6 phi(2) / phi(1) * Lambda^2
    G_eff_general = 6 * pi / (phi1 * Lambda**2)
    Lambda_cc_general = 6 * phi2 / phi1 * Lambda**2

    print()
    print(f"  General results:")
    print(f"    G_eff(f, Lambda) = 6 pi / (phi(1) Lambda^2) = {G_eff_general}")
    print(f"    Lambda_cc(f, Lambda) = 6 phi(2) / phi(1) * Lambda^2 = {Lambda_cc_general}")

    results["general_formulas"] = {
        "G_eff": "6 pi / (phi(1) Lambda^2)",
        "Lambda_cc": "6 phi(2) / phi(1) * Lambda^2",
    }

    # -----------------------------------------------------------------------
    # Step 2: Specific cutoffs
    # -----------------------------------------------------------------------
    print("\n[Step 2] Three specific cutoffs:")

    cutoffs = []

    # Gaussian
    print("\n  (a) Gaussian f(x) = e^{-x}: phi(s) = Gamma(s)")
    phi1_gauss = sp.Integer(1)  # Gamma(1) = 1
    phi2_gauss = sp.Integer(1)  # Gamma(2) = 1
    G_eff_gauss = G_eff_general.subs([(phi1, phi1_gauss), (phi2, phi2_gauss)])
    Lambda_cc_gauss = Lambda_cc_general.subs([(phi1, phi1_gauss), (phi2, phi2_gauss)])
    print(f"    phi(1) = Gamma(1) = 1, phi(2) = Gamma(2) = 1")
    print(f"    G_eff = {G_eff_gauss}")
    print(f"    Lambda_cc = {Lambda_cc_gauss}")
    cutoffs.append({
        "name": "Gaussian",
        "f": "e^{-x}",
        "phi_1": str(phi1_gauss),
        "phi_2": str(phi2_gauss),
        "G_eff": str(G_eff_gauss),
        "Lambda_cc": str(Lambda_cc_gauss),
    })

    # Polynomial e^{-x^2}
    print("\n  (b) Polynomial f(x) = e^{-x^2}: phi(s) = (1/2) Gamma(s/2)")
    phi1_poly = Rational(1, 2) * gamma(Rational(1, 2))  # (1/2) Gamma(1/2)
    phi2_poly = Rational(1, 2) * gamma(1)  # (1/2) Gamma(1) = 1/2
    phi1_poly_simp = sp.simplify(phi1_poly)
    G_eff_poly = G_eff_general.subs([(phi1, phi1_poly), (phi2, phi2_poly)])
    G_eff_poly_simp = sp.simplify(G_eff_poly)
    Lambda_cc_poly = Lambda_cc_general.subs([(phi1, phi1_poly), (phi2, phi2_poly)])
    Lambda_cc_poly_simp = sp.simplify(Lambda_cc_poly)
    print(f"    phi(1) = (1/2) Gamma(1/2) = sqrt(pi)/2 = {phi1_poly_simp}")
    print(f"    phi(2) = (1/2) Gamma(1) = 1/2 = {phi2_poly}")
    print(f"    G_eff = {G_eff_poly_simp}")
    print(f"    Lambda_cc = {Lambda_cc_poly_simp}")
    cutoffs.append({
        "name": "Polynomial e^{-x^2}",
        "f": "e^{-x^2}",
        "phi_1": str(phi1_poly_simp),
        "phi_2": str(phi2_poly),
        "G_eff": str(G_eff_poly_simp),
        "Lambda_cc": str(Lambda_cc_poly_simp),
    })

    # Sharp Theta(1-x)
    print("\n  (c) Sharp f(x) = Theta(1 - x): phi(s) = 1/s")
    phi1_sharp = sp.Integer(1)  # 1/1
    phi2_sharp = Rational(1, 2)  # 1/2
    G_eff_sharp = G_eff_general.subs([(phi1, phi1_sharp), (phi2, phi2_sharp)])
    Lambda_cc_sharp = Lambda_cc_general.subs([(phi1, phi1_sharp), (phi2, phi2_sharp)])
    print(f"    phi(1) = 1, phi(2) = 1/2")
    print(f"    G_eff = {G_eff_sharp}")
    print(f"    Lambda_cc = {Lambda_cc_sharp}")
    cutoffs.append({
        "name": "Sharp",
        "f": "Theta(1 - x)",
        "phi_1": str(phi1_sharp),
        "phi_2": str(phi2_sharp),
        "G_eff": str(G_eff_sharp),
        "Lambda_cc": str(Lambda_cc_sharp),
    })

    results["specific_cutoffs"] = cutoffs

    # -----------------------------------------------------------------------
    # Step 3: Identify cutoff-independent invariants
    # -----------------------------------------------------------------------
    print("\n[Step 3] Cutoff-independent invariants:")
    print()
    print("  Consider the product:")
    print("    G_eff * Lambda^2 = 6 pi / phi(1)  -- depends on phi(1)")
    print("    Lambda_cc / Lambda^2 = 6 phi(2)/phi(1)  -- depends on phi(2)/phi(1)")
    print()
    print("  Therefore the product Lambda_cc * G_eff * Lambda^{-2} = ?")
    print()
    print("  Compute: Lambda_cc * G_eff / Lambda^2:")

    # Lambda_cc * G_eff / Lambda^2:
    # = (6 phi(2)/phi(1) Lambda^2) * (6 pi/(phi(1) Lambda^2)) / Lambda^2
    # = 36 pi phi(2) / (phi(1)^2 Lambda^2)
    invariant_candidate = Lambda_cc_general * G_eff_general / Lambda**2
    invariant_candidate_simp = sp.simplify(invariant_candidate)
    print(f"    Lambda_cc * G_eff / Lambda^2 = {invariant_candidate_simp}")
    print()
    print("  Hmm, depends on phi(1)^2 / phi(2). NOT cutoff-independent.")
    print()
    print("  Consider G_eff * Lambda^2:")
    G_eff_Lambda2 = sp.simplify(G_eff_general * Lambda**2)
    print(f"    G_eff * Lambda^2 = {G_eff_Lambda2}  (depends on phi(1))")
    print()
    print("  Consider Lambda_cc / Lambda^2:")
    Lambda_cc_over_Lambda2 = sp.simplify(Lambda_cc_general / Lambda**2)
    print(f"    Lambda_cc / Lambda^2 = {Lambda_cc_over_Lambda2}  (depends on phi(2)/phi(1))")

    # The RATIO G_eff * Lambda^2 / (Lambda_cc / Lambda^2):
    # = (6 pi / phi(1)) / (6 phi(2)/phi(1))
    # = pi / phi(2)
    # Still depends on phi(2)
    ratio_candidate = G_eff_Lambda2 / Lambda_cc_over_Lambda2
    ratio_candidate_simp = sp.simplify(ratio_candidate)
    print()
    print(f"  Ratio: (G_eff * Lambda^2) / (Lambda_cc / Lambda^2) = {ratio_candidate_simp}")
    print(f"    Still depends on phi(2). Not cutoff-independent.")
    print()

    print("  Conclusion: NO simple cutoff-independent combination of G_eff and Lambda_cc")
    print("  Both depend on phi-moments (calibration data).")
    print()
    print("  HOWEVER: the RATIO R_crit Lambda (from G2's extremum):")
    print("    R_crit * Lambda = sqrt(phi(1)/(6 phi(2)))")
    print()
    # Verify: R_crit Lambda from extremum of (R^3 phi(2)/4) Lambda^4 - (R phi(1)/8) Lambda^2
    # d/dR = 3R^2 phi(2)/4 * Lambda^4 - phi(1)/8 * Lambda^2 = 0
    # R^2 = phi(1)/(6 phi(2) Lambda^2)
    # R_crit Lambda = sqrt(phi(1)/(6 phi(2)))
    R_crit_Lambda = sp.sqrt(phi1 / (6 * phi2))
    print(f"    R_crit Lambda = sqrt(phi(1)/(6 phi(2))) = {R_crit_Lambda}")
    print()
    print("  For Gaussian: R_crit Lambda = 1/sqrt(6) ~ 0.408 (matches G1)")
    print("  For Polynomial: R_crit Lambda = sqrt((sqrt(pi)/2) / (6 * 1/2)) = sqrt(sqrt(pi)/6)")
    R_crit_poly = R_crit_Lambda.subs([(phi1, phi1_poly), (phi2, phi2_poly)])
    print(f"                                  ~ {float(sp.simplify(R_crit_poly)):.4f}")
    print(f"  For Sharp: R_crit Lambda = sqrt(1/(6 * 1/2)) = sqrt(1/3) ~ 0.577")

    results["cutoff_dependence_analysis"] = {
        "G_eff_Lambda2": str(G_eff_Lambda2),
        "Lambda_cc_over_Lambda2": str(Lambda_cc_over_Lambda2),
        "R_crit_Lambda_general": str(R_crit_Lambda),
        "all_cutoff_dependent": True,
        "no_simple_universal_combination": True,
    }

    # -----------------------------------------------------------------------
    # Step 4: Implication: framework predictions are NOT cutoff-independent
    # -----------------------------------------------------------------------
    print("\n[Step 4] Implication for framework predictions:")
    print()
    print("  The framework's predictions for Newton constant, cosmological")
    print("  constant, and preferred radius ALL depend on the choice of cutoff")
    print("  function f. This is standard in CC literature.")
    print()
    print("  The cutoff function moments phi(1), phi(2) are part of the framework's")
    print("  CALIBRATION DATA, not derivable from first principles. They are Class")
    print("  1 calibration data per memory/external_input_three_class_partition.")
    print()
    print("  This sharpens the structural-skeleton-scope reading: the framework's")
    print("  gravity predictions are PROPORTIONAL TO specific combinations of phi-")
    print("  moments. Without first-principles selection of the cutoff, the framework")
    print("  gives numerical predictions only up to overall cutoff calibration.")
    print()
    print("  CONSISTENT WITH STANDARD CC: this is the standard situation in CC's")
    print("  spectral-action framework. The cutoff function is external input.")

    results["framework_implication"] = {
        "all_gravity_predictions_cutoff_dependent": True,
        "cutoff_calibration_class": "Class 1 (external)",
        "consistent_with_standard_CC": True,
        "framework_predictions_proportional_to_phi_moments": True,
    }

    # -----------------------------------------------------------------------
    # Step 5: Verdict
    # -----------------------------------------------------------------------
    print("\n[Step 5] Verdict:")
    print()
    print("  POSITIVE-CALIBRATION-CLARIFICATION:")
    print("  - General formulas derived: G_eff = 6 pi/(phi(1) Lambda^2)")
    print("                              Lambda_cc = 6 phi(2)/phi(1) Lambda^2")
    print("                              R_crit Lambda = sqrt(phi(1)/(6 phi(2)))")
    print()
    print("  - All depend on cutoff function moments (Class 1 calibration data)")
    print()
    print("  - G7's specific numerical values are for Gaussian only")
    print()
    print("  Honest scope: framework's gravity predictions are PROPORTIONAL TO phi-")
    print("  moment combinations; the cutoff function is external input. This is")
    print("  consistent with standard CC.")
    print()
    print("  The framework's STRUCTURAL findings (two-term exactness, closed forms,")
    print("  extremality consistency) are cutoff-INDEPENDENT. The NUMERICAL predictions")
    print("  for G_eff, Lambda_cc are cutoff-DEPENDENT.")

    results["verdict"] = "POSITIVE-CALIBRATION-CLARIFICATION"

    # Save JSON
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)
    return results


if __name__ == "__main__":
    main()
