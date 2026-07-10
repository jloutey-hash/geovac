"""
g2_c2_verification.py — Verification of the c₂ identification for g-2 curvature expansion on S³

Documents the cross-invariant formula:
    c₂ = (2 - B*Δ - F*Δ - F/B) / 5 = 19/100 - 41*π²/25200

where B=42, F=π²/6, Δ=1/40 are the Paper 2 invariants.

This script is standalone — no geovac imports required.
"""

import json
from pathlib import Path

import mpmath
import sympy
from sympy import Rational, pi, sqrt, simplify, expand, nsimplify


def run_verification():
    mp = mpmath.mp
    mp.dps = 50

    results = {}

    # =====================================================================
    # Section 1: Paper 2 invariants
    # =====================================================================
    print("=" * 70)
    print("Section 1: Paper 2 invariants")
    print("=" * 70)

    B = 42
    F_exact = pi ** 2 / 6  # sympy
    Delta = Rational(1, 40)
    c1 = Rational(1, 2)

    B_mp = mpmath.mpf(42)
    F_mp = mp.pi ** 2 / 6
    Delta_mp = mpmath.mpf(1) / 40
    c1_mp = mpmath.mpf(1) / 2

    print(f"  B     = {B}")
    print(f"  F     = pi^2/6 = {mp.nstr(F_mp, 30)}")
    print(f"  Delta = 1/40   = {mp.nstr(Delta_mp, 30)}")
    print(f"  c_1   = 1/2    = {mp.nstr(c1_mp, 30)}")
    print()

    results["B"] = B
    results["F"] = str(F_mp)
    results["Delta"] = str(Delta_mp)
    results["c1"] = 0.5

    # =====================================================================
    # Section 2: Cross-invariant c₂ formula
    # =====================================================================
    print("=" * 70)
    print("Section 2: Cross-invariant c_2 formula")
    print("=" * 70)

    # mpmath computation
    c2_cross_mp = (2 - B_mp * Delta_mp - F_mp * Delta_mp - F_mp / B_mp) / 5
    print(f"  c_2 = (2 - B*Delta - F*Delta - F/B) / 5")
    print(f"       = {mp.nstr(c2_cross_mp, 40)}")

    # sympy exact computation
    B_sym = sympy.Integer(42)
    c2_cross_sym = (2 - B_sym * Delta - F_exact * Delta - F_exact / B_sym) / 5
    c2_cross_expanded = expand(c2_cross_sym)
    print(f"  Exact sympy: {c2_cross_expanded}")

    # Target form: 19/100 - 41*pi^2/25200
    c2_target = Rational(19, 100) - Rational(41, 25200) * pi ** 2
    diff_sym = simplify(c2_cross_expanded - c2_target)
    print(f"  Target:      19/100 - 41*pi^2/25200")
    print(f"  Difference:  {diff_sym}")
    assert diff_sym == 0, f"Symbolic verification FAILED: difference = {diff_sym}"
    print(f"  --> VERIFIED: c_2 = 19/100 - 41*pi^2/25200 (exact)")
    print()

    # Verify the rational decomposition step by step
    print("  Rational decomposition verification:")
    # c2 = (2 - 21/20 - pi^2/240 - pi^2/252) / 5
    inner = 2 - Rational(21, 20) - pi ** 2 / 240 - pi ** 2 / 252
    inner_simplified = expand(inner)
    print(f"    2 - B*Delta - F*Delta - F/B = {inner_simplified}")
    # Should be 19/20 - 41*pi^2/5040
    inner_target = Rational(19, 20) - Rational(41, 5040) * pi ** 2
    assert simplify(inner_simplified - inner_target) == 0
    print(f"    = 19/20 - 41*pi^2/5040  (verified)")
    print(f"    / 5 = 19/100 - 41*pi^2/25200  (verified)")
    print()

    results["c2_cross_mpmath"] = str(c2_cross_mp)
    results["c2_cross_exact"] = "19/100 - 41*pi^2/25200"
    results["c2_symbolic_match"] = True

    # =====================================================================
    # Section 3: Load computed data
    # =====================================================================
    print("=" * 70)
    print("Section 3: Load computed data (n_ext=1, n_int=0..25 + tail)")
    print("=" * 70)

    data_path = Path(__file__).parent / "data" / "g2_extended_nint_v2.json"
    with open(data_path, "r") as f:
        data = json.load(f)

    delta_computed = data["corrected_delta"]
    F2_over_S = data["corrected_F2_over_S"]
    tail_est = data["tail_estimate"]

    print(f"  Data file: {data_path}")
    print(f"  corrected_delta    = {delta_computed}")
    print(f"  corrected_F2/S     = {F2_over_S}")
    print(f"  tail_estimate      = {tail_est:.6e}")

    # At n_ext=1: lambda = n_ext + 3/2 = 5/2, so x = 1/lambda^2 = 4/25
    lam = mpmath.mpf(5) / 2
    x = 1 / lam ** 2  # 4/25 = 0.16
    print(f"  lambda (n_ext=1)   = {mp.nstr(lam, 10)}")
    print(f"  x = 1/lambda^2     = {mp.nstr(x, 10)}")

    delta_mp = mpmath.mpf(delta_computed)

    # Extract c2_apparent: delta = c1*x + c2*x^2 + c3*x^3 + ...
    # So c2_apparent = (delta - c1*x) / x^2
    c2_apparent = (delta_mp - c1_mp * x) / x ** 2
    print(f"  c_2 apparent       = (delta - c1*x) / x^2")
    print(f"                     = {mp.nstr(c2_apparent, 20)}")
    print()

    results["delta_computed"] = delta_computed
    results["lambda_n1"] = 2.5
    results["x_n1"] = float(x)
    results["c2_apparent"] = float(c2_apparent)

    # =====================================================================
    # Section 4: Candidate comparison table
    # =====================================================================
    print("=" * 70)
    print("Section 4: Candidate comparison table")
    print("=" * 70)

    candidates = {}

    # Cross-invariant
    candidates["cross-invariant"] = c2_cross_mp

    # PSLQ rational+pi^2: (6*pi^2 - 58)/7
    candidates["PSLQ (6pi^2-58)/7"] = (6 * mp.pi ** 2 - 58) / 7

    # 25/144 = c1^2 * 25/36
    candidates["25/144"] = mpmath.mpf(25) / 144

    # 7/40 = 7*Delta
    candidates["7/40 = 7*Delta"] = mpmath.mpf(7) / 40

    # 19/100 (rational part only)
    candidates["19/100 (rational)"] = mpmath.mpf(19) / 100

    print(f"  c_2 apparent = {mp.nstr(c2_apparent, 20)}")
    print()
    print(f"  {'Candidate':<25s} {'Value':>25s} {'|diff|':>15s} {'Match digits':>12s}")
    print(f"  {'-'*25} {'-'*25} {'-'*15} {'-'*12}")

    comparison_results = []
    for name, val in candidates.items():
        diff = abs(val - c2_apparent)
        if diff > 0:
            match_digits = -mpmath.log10(diff / abs(c2_apparent))
            md_str = f"{float(match_digits):.1f}"
        else:
            md_str = "exact"
        print(f"  {name:<25s} {mp.nstr(val, 20):>25s} {mp.nstr(diff, 6):>15s} {md_str:>12s}")
        comparison_results.append({
            "name": name,
            "value": str(val),
            "abs_diff": str(diff),
            "match_digits": md_str,
        })

    print()
    results["candidate_comparison"] = comparison_results

    # =====================================================================
    # Section 5: Implied c₃
    # =====================================================================
    print("=" * 70)
    print("Section 5: Implied c_3 from cross-invariant c_2")
    print("=" * 70)

    # delta_pred = c1*x + c2_cross*x^2
    delta_pred = c1_mp * x + c2_cross_mp * x ** 2
    residual = delta_mp - delta_pred
    c3_estimate = residual / x ** 3

    print(f"  delta (computed)   = {mp.nstr(delta_mp, 20)}")
    print(f"  delta (predicted)  = c1*x + c2_cross*x^2")
    print(f"                     = {mp.nstr(delta_pred, 20)}")
    print(f"  residual           = {mp.nstr(residual, 10)}")
    print(f"  c_3 estimate       = residual / x^3 = {mp.nstr(c3_estimate, 10)}")
    print(f"  |c_3|              = {mp.nstr(abs(c3_estimate), 6)}")
    c3_small = abs(c3_estimate) < 2e-7
    print(f"  |c_3| < 2e-7?     {'YES' if c3_small else 'NO'}")
    print()

    results["delta_predicted"] = str(delta_pred)
    results["residual"] = str(residual)
    results["c3_estimate"] = str(c3_estimate)
    results["c3_bound_satisfied"] = c3_small

    # =====================================================================
    # Section 6: Structural decomposition
    # =====================================================================
    print("=" * 70)
    print("Section 6: Structural decomposition")
    print("=" * 70)

    # B*Delta
    BD = B_sym * Delta
    print(f"  B * Delta  = 42 * 1/40 = {BD} = {float(BD)}")
    assert BD == Rational(21, 20)
    print(f"             = 21/20  (verified)")

    # F*Delta
    FD = F_exact * Delta
    FD_expanded = expand(FD)
    print(f"  F * Delta  = (pi^2/6) * (1/40) = {FD_expanded}")
    assert simplify(FD_expanded - pi ** 2 / 240) == 0
    print(f"             = pi^2/240  (verified)")

    # F/B
    FB = F_exact / B_sym
    FB_expanded = expand(FB)
    print(f"  F / B      = (pi^2/6) / 42 = {FB_expanded}")
    assert simplify(FB_expanded - pi ** 2 / 252) == 0
    print(f"             = pi^2/252  (verified)")

    # Numerator: 2 - BD - FD - FB
    num = 2 - BD - FD_expanded - FB_expanded
    num_expanded = expand(num)
    print(f"  2 - B*Delta - F*Delta - F/B = {num_expanded}")

    # Check rational part
    # Collect rational and pi^2 coefficients
    rat_part = 2 - Rational(21, 20)
    assert rat_part == Rational(19, 20)
    print(f"    Rational part: 2 - 21/20 = {rat_part} = 19/20  (verified)")

    # pi^2 coefficient: -1/240 - 1/252
    pi2_coeff = -Rational(1, 240) - Rational(1, 252)
    # LCM(240, 252) = ?
    # 240 = 2^4 * 3 * 5, 252 = 2^2 * 3^2 * 7, LCM = 2^4 * 3^2 * 5 * 7 = 5040
    pi2_coeff_simplified = pi2_coeff
    print(f"    pi^2 coeff: -1/240 - 1/252 = {pi2_coeff_simplified}")
    assert pi2_coeff_simplified == Rational(-41, 5040)
    print(f"              = -41/5040  (verified)")

    # Full numerator check
    num_target = Rational(19, 20) - Rational(41, 5040) * pi ** 2
    assert simplify(num_expanded - num_target) == 0
    print(f"    Full: 19/20 - 41*pi^2/5040  (verified)")

    # Divide by 5
    c2_final = num_target / 5
    c2_final_expanded = expand(c2_final)
    c2_final_target = Rational(19, 100) - Rational(41, 25200) * pi ** 2
    assert simplify(c2_final_expanded - c2_final_target) == 0
    print(f"  Dividing by 5:")
    print(f"    c_2 = 19/100 - 41*pi^2/25200  (verified)")
    print()

    results["BD"] = "21/20"
    results["FD"] = "pi^2/240"
    results["FB"] = "pi^2/252"
    results["pi2_coeff_numerator"] = "-41/5040"
    results["c2_rational_part"] = "19/100"
    results["c2_pi2_coeff"] = "-41/25200"

    # =====================================================================
    # Section 7: T9 consistency check
    # =====================================================================
    print("=" * 70)
    print("Section 7: T9 consistency check")
    print("=" * 70)

    a_rational = Rational(19, 100)
    b_rational = Rational(-41, 25200)

    print(f"  c_2 = a + b * pi^2")
    print(f"  a = {a_rational} = {float(a_rational)}")
    print(f"  b = {b_rational} = {float(b_rational):.10f}")
    print()
    print(f"  T9 theorem: one-loop QED on S^3 produces only pi^{{even}}")
    print(f"  transcendentals (rational multiples of pi^{{2k}}).")
    print()
    print(f"  c_2 content:")
    print(f"    - Rational part:   19/100          (T9 compatible: trivially)")
    print(f"    - pi^2 part:       -41/25200 * pi^2 (T9 compatible: pi^{{even}})")
    print(f"    - Odd zeta:        NONE")
    print(f"    - Catalan G:       NONE")
    print(f"    - Dirichlet beta:  NONE")
    print(f"  --> CONSISTENT with T9 theorem")
    print()

    results["T9_consistent"] = True
    results["transcendental_content"] = "pi^2 only (T9 compatible)"

    # =====================================================================
    # Section 8: Save results
    # =====================================================================
    print("=" * 70)
    print("Section 8: Save results")
    print("=" * 70)

    output_path = Path(__file__).parent / "data" / "g2_c2_verification.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"  Results saved to: {output_path}")
    print()

    # Final summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  c_2 = (2 - B*Delta - F*Delta - F/B) / 5")
    print(f"       = 19/100 - 41*pi^2/25200")
    print(f"       = {mp.nstr(c2_cross_mp, 30)}")
    print(f"  c_2 apparent (from data) = {mp.nstr(c2_apparent, 20)}")
    print(f"  |c_3| estimate           = {mp.nstr(abs(c3_estimate), 6)}")
    print(f"  T9 consistency:            YES (pi^{{even}} only)")
    print(f"  Symbolic verification:     PASSED")
    print()

    return results


if __name__ == "__main__":
    run_verification()
