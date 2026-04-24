"""
Sprint A — α-program reframe: high-precision error budget for K = π(B + F − Δ).

Computes K to 50+ dps, residual R = K − 1/α (with cubic correction), relative
error, and order-of-magnitude estimates for candidate next-order corrections.

Outputs: debug/data/alpha_error_budget.json
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


def main() -> None:
    mp.mp.dps = 60  # 60 decimal places throughout

    # -----------------------------------------------------------------
    # Fixed inputs
    # -----------------------------------------------------------------
    B = mp.mpf(42)
    F = mp.pi**2 / 6                                # ζ_R(2)
    Delta = mp.mpf(1) / 40                          # 1/g_3^Dirac

    # K = π(B + F − Δ)
    inner = B + F - Delta
    K = mp.pi * inner

    # CODATA 2022: α^{-1} = 137.035 999 084(21)
    alpha_inv_codata = mp.mpf("137.035999084")
    alpha_codata = 1 / alpha_inv_codata

    # -----------------------------------------------------------------
    # Subtask 1: residual R₀ = K − α^{-1}_CODATA   (raw, no cubic correction)
    # -----------------------------------------------------------------
    R0 = K - alpha_inv_codata

    # Paper 2 actually uses the cubic α³ − Kα + 1 = 0  (Eq. 9, 10 of paper_2_alpha.tex):
    # α^{-1} + α² = K  ⇒  α^{-1} = K − α²
    # so the formula's α^{-1} = K − α². Compute that explicitly.
    alpha_sq_codata = alpha_codata**2
    alpha_inv_formula = K - alpha_sq_codata
    R = alpha_inv_formula - alpha_inv_codata

    # Confirm by solving the cubic α³ − Kα + 1 = 0 at high precision and taking
    # the smallest positive root.
    # mpmath.polyroots on coefficients [1, 0, -K, 1]
    coeffs = [mp.mpf(1), mp.mpf(0), -K, mp.mpf(1)]
    roots = mp.polyroots(coeffs)
    real_roots = sorted(
        [r.real for r in roots if abs(r.imag) < mp.mpf(10) ** (-40)]
    )
    # Physical root: smallest positive real root (perturbative QED, α ≪ 1)
    alpha_formula = next(r for r in real_roots if r > 0)
    alpha_inv_formula_cubic = 1 / alpha_formula
    R_cubic = alpha_inv_formula_cubic - alpha_inv_codata

    # -----------------------------------------------------------------
    # Subtask 2: relative error
    # -----------------------------------------------------------------
    rel_err_raw = R0 / alpha_inv_codata
    rel_err_formula = R / alpha_inv_codata
    rel_err_cubic = R_cubic / alpha_inv_codata

    # -----------------------------------------------------------------
    # Subtask 3: candidate next-order corrections
    # -----------------------------------------------------------------
    # (a) Higher Seeley-DeWitt — a_3, a_4 enter at order (m_γ/Λ)^{2k}.
    #     On unit S³ the natural cutoff Λ corresponds to the highest mode
    #     n_max = 3 (Paper 2 selection).  Mass scale is set by the photon
    #     loop's IR cutoff = lowest non-zero Hodge-1 eigenvalue, μ₁ = 3.
    #     Take m_γ ~ 1 (unit S³ natural units), Λ ~ n_max = 3, so
    #     (m_γ/Λ)² ~ 1/9.  The Seeley-DeWitt expansion factor is
    #     a_2 = √π/8 (Paper 28 §III); the next term carries an additional
    #     coefficient ~ 1/(48π²) (vacuum-polarization prefactor).
    seeley_a2 = mp.sqrt(mp.pi) / 8
    vacpol_prefactor = 1 / (48 * mp.pi**2)
    sd_estimate = vacpol_prefactor * (mp.mpf(1) / 9) ** 1   # k=1 (next a_3 term)
    # Express as fraction of K
    cand_a = sd_estimate                                     # ~ 2.3e-4 (raw)
    cand_a_rel = cand_a / K

    # (b) Higher-m Casimir truncation — B(4) − B(3).
    #     B(m) = m(m-1)(m+1)(m+2)(2m+1)/20.  B(3)=42, B(4)=162.
    #     Increment ΔB = 120.  Suppression factor: in spectral-action
    #     framing the contribution of higher modes is Boltzmann-suppressed
    #     by exp(−|λ_n|/Λ).  At n=4, |λ_4| = 15, Λ ~ |λ_3| = 8,
    #     so suppression ~ exp(−15/8) ≈ 0.153.  Contribution:
    delta_B = mp.mpf(162 - 42)  # 120
    suppression_b = mp.exp(-mp.mpf(15) / 8)
    cand_b = mp.pi * delta_B * suppression_b      # if it entered K linearly
    cand_b_rel = cand_b / K
    # In a soft-cutoff spectral-action expansion the suppression can be much
    # stronger (Gaussian e^{−λ²/Λ²} not e^{−λ/Λ}); record both:
    suppression_b_gauss = mp.exp(-(mp.mpf(15) / 8) ** 2)
    cand_b_gauss = mp.pi * delta_B * suppression_b_gauss
    cand_b_gauss_rel = cand_b_gauss / K

    # (c) Two-loop QED correction — α/(2π) is the famous Schwinger anomaly
    #     factor.  Naively (1/α) · α/(2π) = 1/(2π) ≈ 0.159 — gigantic.
    #     For this to match observed R ~ 1e-5 (raw) or 1e-7 (with cubic),
    #     we need a suppression factor s such that
    #     (1/(2π)) · s ≈ |R|.  Solve for s:
    schwinger = 1 / (2 * mp.pi)
    s_suppression_raw = abs(R0) / schwinger
    s_suppression_formula = abs(R) / schwinger
    s_suppression_cubic = abs(R_cubic) / schwinger

    # (d) Non-abelian SU(2) Wilson contribution (Paper 30) — magnitude unknown;
    #     we record an "open" placeholder.  At weak coupling β → ∞ the
    #     plaquette expectation ⟨W⟩ → 1 and the action contribution scales
    #     like β·N_plaq.  In the maximal-torus limit (Paper 30 Prop 1) the
    #     SU(2) action reduces to U(1), so the *new* SU(2) content is
    #     suppressed by α² (next-to-leading non-abelian).  Estimate:
    cand_d = alpha_codata**2
    cand_d_rel = cand_d / K  # would-be relative weight if added to K linearly
    # Actually-additive estimate (since K already has dimensions of α^{-1}):
    cand_d_in_R = K * alpha_codata**2  # K times α² ≈ 137 · 5.3e-5 ≈ 7.3e-3

    # (e) Off-shell / non-Fock-shell contribution.
    #     The Fock projection picks up the n_max=3 cutoff exactly.  Off-shell
    #     contributions from n ≥ 4 modes scale as (n_max/n)² = (3/4)², (3/5)², …
    #     leading-correction sum:
    off_shell = sum(mp.mpf(3) / n for n in range(4, 50)) ** 2
    # That overestimates; the natural off-shell weight is the spectral
    # density tail beyond the truncation:
    off_shell_tail = sum(
        mp.mpf(n**2) * mp.exp(-mp.mpf(n) / 3) for n in range(4, 50)
    )
    cand_e = mp.pi * off_shell_tail
    cand_e_rel = cand_e / K

    # -----------------------------------------------------------------
    # Subtask 5: predicted form of next correction
    # -----------------------------------------------------------------
    # If the cubic-corrected residual R_cubic is at the 1e-7 level and
    # K ~ 137, then ΔK ~ R_cubic ~ 1e-5.  Compare to candidate
    # magnitudes:
    target_dK_cubic = R_cubic                       # the ΔK we need
    target_dK_raw = R0
    target_log_cubic = mp.log10(abs(target_dK_cubic))
    target_log_raw = mp.log10(abs(target_dK_raw))

    # -----------------------------------------------------------------
    # Pack results
    # -----------------------------------------------------------------
    def s(x: mp.mpf, n: int = 30) -> str:
        """Format mpf to n significant decimals."""
        return mp.nstr(x, n, strip_zeros=False)

    results = {
        "context": {
            "framing": "Sprint A — high-precision residual + form prediction",
            "B": 42,
            "F_symbolic": "pi^2 / 6",
            "F_value_30dps": s(F, 30),
            "Delta": "1/40",
            "K_symbolic": "pi * (42 + pi^2/6 - 1/40)",
            "alpha_inv_CODATA_2022": "137.035999084(21)",
            "precision_dps": 60,
        },
        "subtask_1_K_and_residual": {
            "K_50dps": s(K, 50),
            "alpha_inv_codata_50dps": s(alpha_inv_codata, 50),
            "R0_raw_residual_K_minus_alphainv": s(R0, 30),
            "R_cubic_corrected_residual_50dps": s(R_cubic, 50),
            "alpha_inv_formula_cubic_50dps": s(alpha_inv_formula_cubic, 50),
            "alpha_formula_cubic": s(alpha_formula, 50),
        },
        "subtask_2_relative_errors": {
            "rel_err_raw_K_only": s(rel_err_raw, 20),
            "rel_err_first_order_K_minus_alphasq": s(rel_err_formula, 20),
            "rel_err_cubic_root": s(rel_err_cubic, 20),
            "Paper2_quoted_value": "8.8e-8",
            "verified": abs(rel_err_cubic) < mp.mpf("1e-7"),
        },
        "subtask_3_candidate_corrections": {
            "candidate_a_seeley_dewitt_a3": {
                "form": "(1/(48 pi^2)) * (m_gamma / Lambda)^2 with m_gamma=1, Lambda=n_max=3",
                "magnitude": s(cand_a, 12),
                "rel_to_K": s(cand_a_rel, 12),
                "log10": s(mp.log10(abs(cand_a)), 6),
            },
            "candidate_b_higher_Casimir": {
                "form_linear": "pi * (B(4)-B(3)) * exp(-|lambda_4|/|lambda_3|)",
                "magnitude_linear_suppression": s(cand_b, 12),
                "log10_linear": s(mp.log10(abs(cand_b)), 6),
                "form_gauss": "pi * Delta_B * exp(-(|lambda_4|/|lambda_3|)^2)",
                "magnitude_gauss": s(cand_b_gauss, 12),
                "log10_gauss": s(mp.log10(abs(cand_b_gauss)), 6),
            },
            "candidate_c_two_loop_QED": {
                "form": "Schwinger 1/(2 pi) times suppression factor s",
                "schwinger_factor_inv2pi": s(schwinger, 12),
                "required_suppression_for_R_raw": s(s_suppression_raw, 8),
                "required_suppression_for_R_cubic": s(s_suppression_cubic, 8),
                "interpretation_raw": "log10(s_raw) ~= " + s(mp.log10(abs(s_suppression_raw)), 6),
                "interpretation_cubic": "log10(s_cubic) ~= " + s(mp.log10(abs(s_suppression_cubic)), 6),
            },
            "candidate_d_su2_wilson": {
                "form": "alpha^2 (next-to-leading non-abelian beyond Cartan torus)",
                "magnitude_alpha_squared": s(cand_d, 12),
                "magnitude_K_times_alphasq": s(cand_d_in_R, 12),
                "log10_alphasq": s(mp.log10(abs(cand_d)), 6),
                "log10_K_alphasq": s(mp.log10(abs(cand_d_in_R)), 6),
            },
            "candidate_e_off_shell": {
                "form": "pi * sum_{n>=4} n^2 exp(-n/n_max)",
                "magnitude": s(cand_e, 12),
                "rel_to_K": s(cand_e_rel, 12),
                "log10": s(mp.log10(abs(cand_e)), 6),
            },
        },
        "subtask_4_best_match_table": {
            "target_DeltaK_raw_log10": s(target_log_raw, 6),
            "target_DeltaK_cubic_log10": s(target_log_cubic, 6),
            "candidate_logs_DeltaK": {
                "a_seeley_dewitt": s(mp.log10(abs(cand_a)), 6),
                "b_casimir_linear":  s(mp.log10(abs(cand_b)), 6),
                "b_casimir_gauss":   s(mp.log10(abs(cand_b_gauss)), 6),
                "c_required_suppression_raw":  s(mp.log10(abs(s_suppression_raw)), 6),
                "c_required_suppression_cubic": s(mp.log10(abs(s_suppression_cubic)), 6),
                "d_alphasq":         s(mp.log10(abs(cand_d)), 6),
                "d_K_alphasq":       s(mp.log10(abs(cand_d_in_R)), 6),
                "e_off_shell":       s(mp.log10(abs(cand_e)), 6),
            },
        },
        "subtask_5_form_prediction": {
            "expected_factor_alpha_sq": "yes — α² ≈ 5.3e-5; matches RAW R0 well",
            "expected_factor_inv_Lambda_sq": "candidate (a) suggests 1/Λ² with Λ ~ n_max",
            "expected_factor_pi4": "F = π²/6 is leading; next zeta is ζ(4)=π⁴/90, suppressed by π²/15",
            "expected_factor_zeta3": "appears in Dirac sector at d_max=4 (Phase 4I D3); could enter K at next order",
        },
        "next_zeta_F_correction": {
            "form": "ζ(4) = π^4/90, would replace F=ζ(2) at next packing exponent d_max+2=6",
            "magnitude_zeta4": s(mp.zeta(4), 12),
            "rel_to_F": s(mp.zeta(4) / F, 12),
        },
    }

    out_dir = Path(__file__).resolve().parents[2] / "debug" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "alpha_error_budget.json"
    with open(out_file, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"Wrote {out_file}")

    # Also print the headline numbers to stdout for the report-back summary
    print("\n=== Sprint A headline ===")
    print(f"K (50 dps)            = {s(K, 50)}")
    print(f"alpha^-1 CODATA       = {s(alpha_inv_codata, 30)}")
    print(f"R0 = K - alpha^-1     = {s(R0, 20)}")
    print(f"R_cubic (cubic root)  = {s(R_cubic, 20)}")
    print(f"rel err raw           = {s(rel_err_raw, 12)}")
    print(f"rel err cubic         = {s(rel_err_cubic, 12)}")
    print(f"alpha^2               = {s(alpha_codata**2, 12)}")
    print(f"|R0| / alpha^2        = {s(abs(R0) / alpha_codata**2, 6)}  (ratio if R0 ~ alpha^2)")
    print(f"|R_cubic| / alpha^4   = {s(abs(R_cubic) / alpha_codata**4, 6)}  (ratio if R_cubic ~ alpha^4)")
    print(f"|R_cubic| / alpha^3   = {s(abs(R_cubic) / alpha_codata**3, 6)}  (ratio if R_cubic ~ alpha^3)")


if __name__ == "__main__":
    main()
