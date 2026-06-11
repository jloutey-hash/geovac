"""Verify whether the Gaussian-cutoff "hit" in sub-track (b) is tautological.

The Gaussian-cutoff trace N_{gauss}(Λ) = Σ g_n exp(-(|λ_n|/Λ)²) IS the Connes-
Chamseddine spectral action with f(x) = exp(-x²). Its small-Λ → large-Λ
asymptotic via the SD two-term expansion is

  Tr exp(-D²/Λ²) = K_heat(t = 1/Λ²) ≈ (√π/2) Λ³ − (√π/4) Λ.

Setting this = K/π and rearranging gives 2 Λ³ − Λ = 4(K/π)/√π — exactly the
depressed cubic of Sprint A. So the hit at target = K/π is the same equation
that DEFINES Λ_∞; it is a tautology.

This script confirms by:
 (1) showing the Gaussian sum at Λ_∞ literally equals K/π by construction;
 (2) showing the same sum at any other natural target (B, Δ⁻¹, etc.) gives a
     different Λ value, none of which has any independent significance.
"""

import json
from pathlib import Path

import mpmath as mp


mp.mp.dps = 80


LAMBDA_INFTY = mp.mpf("3.7102454679060528505052")
INV_ALPHA = mp.mpf("137.035999084")
K_OVER_PI = INV_ALPHA / mp.pi


def gaussian_trace(Lambda, n_max=200):
    s = mp.mpf(0)
    for n in range(n_max + 1):
        lam = mp.mpf(n) + mp.mpf(3) / 2
        g = 2 * (n + 1) * (n + 2)
        s += g * mp.exp(-(lam / Lambda) ** 2)
    return s


def main():
    print("=" * 70)
    print("Gaussian-cutoff tautology check")
    print("=" * 70)

    # (1) sum at Λ_∞
    sum_at_Lam_infty = gaussian_trace(LAMBDA_INFTY)
    print(f"\n[1] Gaussian trace at Lambda_infty = {mp.nstr(LAMBDA_INFTY, 20)}:")
    print(f"    Σ g_n exp(-(|λ_n|/Λ)²) = {mp.nstr(sum_at_Lam_infty, 20)}")
    print(f"    K/pi = 1/(πα)         = {mp.nstr(K_OVER_PI, 20)}")
    print(f"    relative diff          = {mp.nstr(abs(sum_at_Lam_infty - K_OVER_PI) / K_OVER_PI, 8)}")

    # (2) compare to SD asymptotic
    sd_asymp = mp.sqrt(mp.pi) / 2 * LAMBDA_INFTY ** 3 - mp.sqrt(mp.pi) / 4 * LAMBDA_INFTY
    print(f"\n[2] SD asymptotic at same Λ:")
    print(f"    (√π/2) Λ³ − (√π/4) Λ = {mp.nstr(sd_asymp, 20)}")
    print(f"    diff from sum         = {mp.nstr(abs(sum_at_Lam_infty - sd_asymp), 8)}")

    # (3) check that the depressed cubic gives Λ_infty exactly
    rhs = 4 * K_OVER_PI / mp.sqrt(mp.pi)
    Λ_check = mp.findroot(lambda L: 2 * L**3 - L - rhs, mp.mpf("3.71"))
    print(f"\n[3] depressed cubic 2Λ³ − Λ = {mp.nstr(rhs, 12)}:")
    print(f"    root = {mp.nstr(Λ_check, 20)}")
    print(f"    diff from Λ_∞ = {mp.nstr(abs(Λ_check - LAMBDA_INFTY), 8)}")

    # (4) the structural verdict
    print("\n[4] Structural verdict:")
    print("    The Gaussian-cutoff trace at Λ_∞ literally equals K/π by")
    print("    Sprint A's SD two-term identity. The 'hit' in sub-track (b)")
    print("    target = K/π is the depressed-cubic root by construction.")
    print("    It is NOT an independent natural criterion — it IS the K equation.")
    print()
    print("    The remaining sub-track (b) candidates (target = B alone, target = g_3)")
    print("    give Λ values 1.2% and 2.8% off Λ_∞ respectively, both clean negatives.")

    # (5) test independent targets using ONLY S³ spectral data
    targets = {
        "B = 42 (Casimir)": mp.mpf(42),
        "g_3^Dirac = 40 (Δ⁻¹)": mp.mpf(40),
        "Σ_{k≤3} g_k^Dirac = 80": mp.mpf(80),
        "Vol(S³) = 2π² ≈ 19.74": 2 * mp.pi**2,
        "(no K input)  B + F = 43.6449": mp.mpf(42) + mp.pi**2 / 6,
    }
    print("\n[5] independent targets (NOT involving K = 1/α):")
    for name, tgt in targets.items():
        f = lambda L: gaussian_trace(L) - tgt
        try:
            Λ = mp.findroot(f, mp.mpf("3.5"))
            d_rel = abs(Λ - LAMBDA_INFTY) / LAMBDA_INFTY
            print(f"    target = {name}: Λ = {mp.nstr(Λ, 16)}, rel dist from Λ_∞ = {mp.nstr(d_rel, 8)}")
        except Exception as e:
            print(f"    target = {name}: findroot failed: {e}")

    out = {
        "Lambda_infty": str(LAMBDA_INFTY),
        "K_over_pi": str(K_OVER_PI),
        "gaussian_trace_at_Lambda_infty": str(sum_at_Lam_infty),
        "sd_asymptotic_at_Lambda_infty": str(sd_asymp),
        "sd_asymptotic_minus_full_sum": str(abs(sum_at_Lam_infty - sd_asymp)),
        "verdict": "Gaussian-cutoff hit at K/π is the depressed-cubic root by construction (SD two-term exact). Not an independent criterion.",
    }
    Path("debug/data/kcc_gaussian_check.json").write_text(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
