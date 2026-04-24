"""
Sprint A v2 — focused verification of the v1 partial finding:
  R_cubic ≈ π³α³ to within 0.25%, with apparent next-order C ~ 1/3.

This script tests whether the v1 claim conflates:
  (A) R_cubic = K - 1/α  (which equals α² trivially from the cubic
                          α³ - Kα + 1 = 0, since 1/α = K - α²)
  (B) R_predict = K - 1/α - α²  (the residual AFTER subtracting the
                                  trivial cubic offset, i.e. the
                                  α³-and-higher part)

Then asks: is R_predict / (π³ · α³) close to 1?

Outputs a JSON with the verified high-precision residual values and a
short structural reading.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


def s(x: mp.mpf, n: int = 30) -> str:
    """Format mpf to n significant decimals."""
    return mp.nstr(x, n, strip_zeros=False)


def main() -> None:
    mp.mp.dps = 80  # 80 dps for headroom on α^3 ~ 10^-7

    # -----------------------------------------------------------------
    # Inputs (Paper 2 §III, IV.B/C)
    # -----------------------------------------------------------------
    B = mp.mpf(42)
    F = mp.pi**2 / 6                    # ζ_R(2) (Phase 4F α-J)
    Delta = mp.mpf(1) / 40              # 1/g_3^Dirac (Phase 4H SM-D)
    K = mp.pi * (B + F - Delta)         # Paper 2 Eq. (10) inner

    # CODATA 2022: α^{-1} = 137.035 999 084(21)
    alpha_inv_codata = mp.mpf("137.035999084")
    alpha_codata = 1 / alpha_inv_codata

    # -----------------------------------------------------------------
    # Cubic α³ - Kα + 1 = 0  (Paper 2 Eq. 9)
    # Smallest positive root = formula prediction for α
    # -----------------------------------------------------------------
    coeffs = [mp.mpf(1), mp.mpf(0), -K, mp.mpf(1)]
    roots = mp.polyroots(coeffs)
    real_roots = sorted(r.real for r in roots if abs(r.imag) < mp.mpf(10) ** (-50))
    alpha_formula = next(r for r in real_roots if r > 0)
    alpha_inv_formula = 1 / alpha_formula

    # -----------------------------------------------------------------
    # SUBTASK 1: clarify R_cubic vs R_predict
    # -----------------------------------------------------------------
    # The cubic α³ - Kα + 1 = 0 ⇒ K = α² + 1/α exactly.  So if you
    # plug *the formula's own* α back into K - 1/α you get exactly α²
    # (no information; tautology).
    #
    # The PHYSICAL residual is between the formula's α and CODATA α.
    # Two natural definitions of "residual":
    #
    # R_cubic_codata := K - 1/α_codata
    #   = α_codata² + (formula_correction)
    #   = α_codata² + (1/α_formula - 1/α_codata)
    #
    # R_predict := R_cubic_codata - α_codata²
    #   = 1/α_formula - 1/α_codata
    #   (= the "raw" first-order residual the v1 paper reports as 8.8e-8
    #      relative to K)
    #
    # The v1 claim "R_cubic ≈ π³α³ to 0.25%" is plausible only if it
    # means R_predict (i.e., AFTER stripping the trivial α² piece).
    # Verify directly.

    R_cubic_codata = K - alpha_inv_codata               # ≈ α² (large)
    alpha_sq = alpha_codata**2
    R_predict = R_cubic_codata - alpha_sq               # the "post-α² residual"

    # Equivalent: R_predict = 1/α_formula - 1/α_codata, by the cubic identity
    # for the formula α (K = α_f² + 1/α_f), so R_predict =
    # K - α_codata² - 1/α_codata
    #              =  α_f² + 1/α_f - α_codata² - 1/α_codata
    # (only equals 1/α_formula - 1/α_codata if α_f = α_codata; for the actual
    # offset the two terms both contribute)
    R_predict_check = (1 / alpha_formula) - alpha_inv_codata + (
        alpha_formula**2 - alpha_codata**2
    )
    # Sanity: should be exactly K - α_codata² - 1/α_codata up to 80 dps
    assert abs(R_predict_check - R_predict) < mp.mpf("1e-70")

    # -----------------------------------------------------------------
    # The π³α³ comparison
    # -----------------------------------------------------------------
    pi3_alpha3 = mp.pi**3 * alpha_codata**3
    ratio_R_predict_over_pi3alpha3 = R_predict / pi3_alpha3

    # Also test plain α^3 (no π³)
    alpha_cubed = alpha_codata**3
    ratio_R_predict_over_alpha3 = R_predict / alpha_cubed

    # -----------------------------------------------------------------
    # SUBTASK 3: extract the next-order coefficient C in
    #            R_predict = π³α³ (1 + Cα + O(α²))
    # If R_predict / (π³α³) = 1 + Cα + ..., then C ≈ (ratio - 1) / α
    # -----------------------------------------------------------------
    delta_from_unity = ratio_R_predict_over_pi3alpha3 - 1
    C_estimate = delta_from_unity / alpha_codata

    # Compare to candidate exact values for C:
    candidates_C = {
        "1/3": mp.mpf(1) / 3,
        "1/dim_S3=1/3": mp.mpf(1) / 3,    # alias
        "pi/9": mp.pi / 9,
        "1/2": mp.mpf(1) / 2,
        "1/4": mp.mpf(1) / 4,
        "1": mp.mpf(1),
        "pi/3": mp.pi / 3,
        "1/(2pi)": 1 / (2 * mp.pi),
        "ln_2": mp.log(2),
        "1/n_max=1/3": mp.mpf(1) / 3,     # alias
    }
    C_diffs = {
        name: float(abs(C_estimate - val) / abs(C_estimate))
        for name, val in candidates_C.items()
    }

    # -----------------------------------------------------------------
    # SUBTASK 2: structural interpretation of π³ prefactor
    # -----------------------------------------------------------------
    # Three natural candidates for π³:
    #  (i)   Vol(S^5) = π³  (spectral-action a_4 on a 5-dim manifold;
    #        matches Paper 24's Bargmann-Segal S^5 lattice).
    #  (ii)  Three external π factors stacked (one per geometric layer:
    #        Hopf base π, Fock-Dirichlet π, Dirac-degeneracy π — though
    #        Δ contributes 1/π, not π).  The composition rule already
    #        contains one outer π in K = π(...); π³α³ would then come
    #        from K³α³ ∼ (π × O(1))³ α³.
    #  (iii) Three-loop QED scaling: each loop carries one 1/π from the
    #        loop integral, so a 3-loop diagram has α³/π³.  An *inverse*
    #        suppression by π³ at the same α³ order is the signature of
    #        a SUM over three loops with a π³ enhancement that overcomes
    #        the loop-factor suppression.

    pi_cubed = mp.pi**3
    vol_S5 = pi_cubed                              # Vol(S^5) at unit radius
    K_cubed_times_alpha3 = K**3 * alpha_codata**3  # ≈ (π·B)³ α³ ~ π³B³α³
    K_cubed = K**3
    pi_to_3 = mp.pi**3
    # Order-of-magnitude check that π³ is the right prefactor (not, e.g., π² or π⁴):
    pi2_alpha3_ratio = R_predict / (mp.pi**2 * alpha_codata**3)
    pi4_alpha3_ratio = R_predict / (mp.pi**4 * alpha_codata**3)

    # -----------------------------------------------------------------
    # Pack & write
    # -----------------------------------------------------------------
    results = {
        "context": {
            "framing": "Sprint A v2 — verify v1 partial finding R_cubic ≈ π³α³",
            "B": int(B),
            "F_symbolic": "π²/6",
            "Delta_symbolic": "1/40",
            "K_symbolic": "π(42 + π²/6 - 1/40)",
            "K_value_50dps": s(K, 50),
            "alpha_inv_CODATA_2022": "137.035999084(21)",
            "precision_dps": 80,
            "v1_claim": "R_cubic ≈ π³α³ to 0.25%; next-order C ~ 1/3",
        },
        "subtask_1_residual_clarification": {
            "definition_R_cubic_codata": "K - 1/α_codata (≈ α² trivially via cubic identity)",
            "definition_R_predict": "(K - 1/α_codata) - α² (the post-cubic residual)",
            "alpha_inv_codata_50dps": s(alpha_inv_codata, 50),
            "alpha_inv_formula_50dps": s(alpha_inv_formula, 50),
            "alpha_formula_50dps": s(alpha_formula, 50),
            "K_minus_alphainv_codata": s(R_cubic_codata, 30),
            "alpha_sq_codata": s(alpha_sq, 30),
            "R_predict_post_alphasq": s(R_predict, 30),
            "v1_finding_assessment": (
                "v1 partial finding does NOT mean K - 1/α (which is "
                "trivially ≈ α² ≈ 5.33e-5).  The substantive claim is "
                "about R_predict = K - 1/α - α², the post-α² residual."
            ),
        },
        "subtask_2_pi3_alpha3_test": {
            "alpha_cubed_codata": s(alpha_cubed, 30),
            "pi3_alpha3_codata": s(pi3_alpha3, 30),
            "R_predict_30dps": s(R_predict, 30),
            "ratio_R_predict_over_pi3_alpha3": s(ratio_R_predict_over_pi3alpha3, 20),
            "ratio_R_predict_over_alpha3_no_pi": s(ratio_R_predict_over_alpha3, 20),
            "pi2_alpha3_ratio": s(pi2_alpha3_ratio, 12),
            "pi4_alpha3_ratio": s(pi4_alpha3_ratio, 12),
            "fractional_distance_from_pi3alpha3": s(
                abs(ratio_R_predict_over_pi3alpha3 - 1), 12
            ),
            "fractional_distance_pct": s(
                100 * abs(ratio_R_predict_over_pi3alpha3 - 1), 6
            ),
        },
        "subtask_3_next_order_coefficient_C": {
            "ansatz": "R_predict = π³α³ (1 + Cα + O(α²))",
            "C_estimate": s(C_estimate, 12),
            "C_estimate_in_dimensionless_units": s(C_estimate, 8),
            "v1_claim_C_one_third": "1/3 (suggested as 1/dim(S³))",
            "C_diffs_relative": C_diffs,
            "interpretation": (
                "C is the leading correction coefficient.  The v1 claim "
                "C ~ 1/3 = 1/dim(S³) is suggestive but the numerical fit "
                "cannot fix it to a unique closed form from a single data point."
            ),
        },
        "subtask_4_pi3_structural_reading": {
            "primary_candidate": "Vol(S^5) = π³",
            "primary_justification": (
                "S^5 is the natural ambient sphere for the Bargmann-Segal "
                "lattice (Paper 24).  Spectral-action a_4 coefficient on "
                "an odd-dim sphere of dim d carries Vol(S^d).  The Hopf "
                "fibration S^1 → S^3 → S^2 sits inside S^5 ⊃ S^3 by the "
                "complex Hopf S^1 → S^3 → CP^1 = S^2 viewed in C².  The "
                "next term in the K-expansion is naturally Vol(S^5)·α^3 "
                "from this enclosing geometry."
            ),
            "secondary_candidate": "K³ ≈ (π × 43.6)³ ≈ π³ × 8.3e4",
            "secondary_justification": (
                "Iterating the cubic α^3 = Kα - 1 once more:  α^3 = Kα - 1, "
                "α^6 = (Kα-1)² ≈ K²α² for small α.  At order α^3 the "
                "iteration carries a factor K which contains exactly one "
                "outer π.  Three iterations of the cubic at next order "
                "would give π³ from three K factors collapsing — though "
                "the dimensional accounting needs care."
            ),
            "tertiary_candidate": "Three-loop QED ~ α^3 with π enhancement",
            "tertiary_justification": (
                "Standard QED 3-loop magnetic moment coefficient is "
                "~1.18 × (α/π)^3, i.e. α^3/π^3.  The OBSERVED π³α³ enters "
                "at the OPPOSITE sign of the loop counting — it is an "
                "ENHANCEMENT, not a suppression — consistent with a "
                "geometric (volume-of-ambient-sphere) origin rather than "
                "a loop-integral origin."
            ),
            "preferred": (
                "Vol(S^5) = π³.  Justification: spectral-action expansion "
                "on S^3 ↪ R^4 ↪ S^4/S^5 enclosing geometry; matches "
                "Paper 24 Bargmann-Segal S^5 lattice; sign and magnitude "
                "consistent with enhancement from a higher-dim spectral "
                "trace, not loop suppression."
            ),
        },
        "subtask_5_paper2_integration_recommendation": {
            "where": "Paper 2 §IV.G (Open questions, item 'Understand the residual 8.8×10⁻⁸')",
            "what_to_add": (
                "(a) The residual K - 1/α has a trivial component α² from "
                "the cubic; the substantive residual is R_predict = "
                "K - 1/α - α² ≈ {0:s}.  "
                "(b) Numerical comparison: R_predict / (π³α³) = {1:s}, "
                "i.e. {2:s}% from unity at 80-dps precision.  "
                "(c) Conjectural reading: π³ = Vol(S^5) suggests the next "
                "K-expansion term originates in the ambient geometry of "
                "the Bargmann-Segal S^5 lattice (Paper 24).  "
                "(d) The next-order coefficient C in "
                "R_predict = π³α³(1 + Cα + ...) is C ≈ {3:s}; the "
                "tentative reading C ~ 1/3 = 1/dim(S³) cannot be fixed by "
                "a single CODATA data point.  Status: structural hint, "
                "NOT a derivation."
            ).format(
                s(R_predict, 8),
                s(ratio_R_predict_over_pi3alpha3, 8),
                s(100 * abs(ratio_R_predict_over_pi3alpha3 - 1), 4),
                s(C_estimate, 6),
            ),
            "scope_warnings": [
                "Phases 4B-4I closed nine mechanisms for K (CLAUDE.md §3).",
                "This is a STRUCTURAL HINT for the residual, not a derivation.",
                "Paper 2 stays conjectural; do NOT remove that label.",
                "Single CODATA data point cannot fix C uniquely — ratio test only.",
            ],
        },
    }

    out_dir = Path(__file__).resolve().parents[2] / "debug" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "alpha_error_budget_v2.json"
    with open(out_file, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"Wrote {out_file}")

    # Headline summary (ASCII-only for Windows console compatibility)
    print()
    print("=== Sprint A v2 headline ===")
    print(f"K (50 dps)                  = {s(K, 50)}")
    print(f"alpha^-1 CODATA             = {s(alpha_inv_codata, 30)}")
    print(f"K - 1/alpha (~ alpha^2)     = {s(R_cubic_codata, 16)}  (alpha^2 = {s(alpha_sq, 16)})")
    print(f"R_predict = above - alpha^2 = {s(R_predict, 16)}")
    print(f"pi^3 * alpha^3              = {s(pi3_alpha3, 16)}")
    print(f"R_predict / (pi^3 alpha^3)  = {s(ratio_R_predict_over_pi3alpha3, 12)}")
    print(f"|ratio - 1| (%)             = {s(100 * abs(ratio_R_predict_over_pi3alpha3 - 1), 6)} %")
    print(f"R_predict / alpha^3 (no pi) = {s(ratio_R_predict_over_alpha3, 12)}")
    print(f"C estimate                  = {s(C_estimate, 12)}  (vs claimed 1/3 = 0.3333)")


if __name__ == "__main__":
    main()
