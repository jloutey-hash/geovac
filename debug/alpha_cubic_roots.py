"""
Alpha Cubic Root Survey & QED Running
======================================

Part A: Solve α³ - Kα + 1 = 0 for all three roots and compare
        transformations against known physical constants.

Part B: One-loop QED running from the paper's α to find the scale
        where it matches CODATA α = 1/137.035999084.

Reference: Paper 2 — "The Fine Structure Constant from Spectral
Geometry of the Hopf Fibration" (papers/conjectures/paper_2_alpha.tex)
"""

import numpy as np
from numpy.polynomial import polynomial as P


# ============================================================
# Constants
# ============================================================

# Spectral invariants of the Hopf bundle
B = 42                          # Degeneracy-weighted Casimir trace (S² base, n_max=3)
F = np.pi**2 / 6                # ζ(2), spectral zeta of S¹ fiber
Delta = 1 / 40                  # Boundary correction from S³ total space

K = np.pi * (B + F - Delta)     # Coupling constant

# CODATA 2018
ALPHA_CODATA = 1 / 137.035999084
M_E_MEV = 0.51099895            # Electron mass in MeV

# Known physical constants for comparison
KNOWN_CONSTANTS = {
    "α_em (CODATA)":                    ALPHA_CODATA,
    "α_em⁻¹ (CODATA)":                 137.035999084,
    "sin²θ_W (MS-bar, M_Z)":           0.23122,
    "sin²θ_W (low energy)":            0.2387,
    "α_s(M_Z)":                         0.1179,
    "α_s(1 GeV)":                       0.50,
    "G_F (GeV⁻²)":                     1.1663788e-5,
    "y_e (electron Yukawa)":            2.94e-6,
    "α(M_Z) = 1/127.95":               1 / 127.95,
    "sin(θ_Cabibbo)":                   0.22500,
    "cos(θ_Cabibbo)":                   0.97436,
    "tan(θ_W)":                         0.54860,      # tan of Weinberg angle
    "α_em / π":                         ALPHA_CODATA / np.pi,
    "2π":                               2 * np.pi,
    "4π":                               4 * np.pi,
    "π²":                               np.pi**2,
    "1/(4π)":                           1 / (4 * np.pi),
    "e² (= 4πα)":                       4 * np.pi * ALPHA_CODATA,
    "m_e/m_μ":                          0.00483633,
    "m_e/m_τ":                          0.000287585,
    "m_u/m_t":                          1.29e-5,
    "Euler γ":                          0.5772156649,
}


# ============================================================
# Part A: Cubic root survey
# ============================================================

def solve_cubic() -> np.ndarray:
    """Solve α³ - Kα + 1 = 0 using numpy (coefficients in ascending order)."""
    # numpy.polynomial.polynomial uses ascending powers: c0 + c1*x + c2*x² + c3*x³
    coeffs = [1, -K, 0, 1]  # 1 - Kα + 0·α² + α³
    roots = P.polyroots(coeffs)
    # Return only real parts (all roots are real for this cubic)
    return np.sort(np.real(roots))


def compute_transformations(r: complex) -> dict:
    """Compute natural transformations of a root."""
    r = float(np.real(r))
    transforms = {}
    transforms["r"] = r
    transforms["1/r"] = 1 / r if r != 0 else np.inf
    transforms["|r|"] = abs(r)
    transforms["1/|r|"] = 1 / abs(r) if r != 0 else np.inf
    transforms["r²"] = r**2
    transforms["1/r²"] = 1 / r**2 if r != 0 else np.inf
    transforms["r³"] = r**3
    transforms["√|r|"] = np.sqrt(abs(r))
    transforms["sin²(arctan(r))"] = np.sin(np.arctan(r))**2
    transforms["cos²(arctan(r))"] = np.cos(np.arctan(r))**2
    transforms["r/(1+r²)"] = r / (1 + r**2)
    transforms["r²/(1+r²)"] = r**2 / (1 + r**2)
    transforms["1/(1+r²)"] = 1 / (1 + r**2)
    transforms["r/π"] = r / np.pi
    transforms["r·π"] = r * np.pi
    transforms["r/(2π)"] = r / (2 * np.pi)
    transforms["r/(4π)"] = r / (4 * np.pi)
    transforms["4π·r"] = 4 * np.pi * r
    transforms["ln|r|"] = np.log(abs(r)) if r != 0 else -np.inf
    transforms["exp(-1/|r|)"] = np.exp(-1 / abs(r)) if r != 0 else 0
    transforms["-r"] = -r
    transforms["1/(-r)"] = -1/r if r != 0 else np.inf
    transforms["sin(r)"] = np.sin(r)
    transforms["cos(r)"] = np.cos(r)
    return transforms


def find_matches(transforms: dict, threshold: float = 0.01) -> list:
    """Find matches between transformed values and known constants."""
    matches = []
    for tname, tval in transforms.items():
        if not np.isfinite(tval) or tval == 0:
            continue
        for cname, cval in KNOWN_CONSTANTS.items():
            if cval == 0:
                continue
            rel_err = abs(tval - cval) / abs(cval)
            if rel_err < threshold:
                matches.append((tname, tval, cname, cval, rel_err))
    return sorted(matches, key=lambda x: x[4])


def part_a():
    """Run Part A: cubic root survey."""
    print("=" * 78)
    print("PART A: CUBIC ROOT SURVEY")
    print("=" * 78)
    print()

    # Display K
    print(f"Spectral invariants:")
    print(f"  B (Casimir trace)     = {B}")
    print(f"  F (ζ(2))              = {F:.10f}")
    print(f"  Δ (boundary)          = {Delta:.10f}")
    print(f"  B + F - Δ             = {B + F - Delta:.10f}")
    print(f"  K = π(B + F - Δ)      = {K:.10f}")
    print()

    # Verify discriminant
    disc = 4 * K**3 - 27
    print(f"Discriminant 4K³ - 27   = {disc:.6f}  (> 0 → three real roots)")
    print()

    roots = solve_cubic()
    print(f"Roots of α³ - Kα + 1 = 0:")
    print("-" * 78)

    for i, r in enumerate(roots):
        rv = float(np.real(r))
        print(f"\n  Root {i+1}: α = {rv:+.12f}")
        print(f"          1/α = {1/rv:+.10f}")
        print(f"          α²  = {rv**2:.12f}")

        # Verify it satisfies the cubic
        residual = rv**3 - K * rv + 1
        print(f"          Residual (α³ - Kα + 1) = {residual:.2e}")

    print()

    # Self-consistency form: 1/α + α² = K
    print("Self-consistency 1/α + α² = K:")
    for i, r in enumerate(roots):
        rv = float(np.real(r))
        sc = 1/rv + rv**2
        print(f"  Root {i+1}: 1/α + α² = {sc:.10f}  (K = {K:.10f})")

    print()
    print("=" * 78)
    print("TRANSFORMATION SURVEY — matches to known constants (< 1% error)")
    print("=" * 78)

    for i, r in enumerate(roots):
        rv = float(np.real(r))
        print(f"\n{'─' * 78}")
        print(f"Root {i+1}: α = {rv:+.12f}  (1/α = {1/rv:+.10f})")
        print(f"{'─' * 78}")

        transforms = compute_transformations(rv)
        matches = find_matches(transforms, threshold=0.01)

        if not matches:
            print("  No matches below 1% threshold.")
            continue

        print(f"  {'Transform':<24s} {'Value':>16s}   {'Matches':>28s} {'Known':>16s} {'RelErr':>12s}")
        print(f"  {'─'*24} {'─'*16}   {'─'*28} {'─'*16} {'─'*12}")

        for tname, tval, cname, cval, rel_err in matches:
            flag = "***" if rel_err < 0.001 else " * " if rel_err < 0.01 else "   "
            print(f"  {tname:<24s} {tval:>16.10f}   {cname:>28s} {cval:>16.10f} {rel_err:>11.2e} {flag}")


# ============================================================
# Part B: QED Running
# ============================================================

def alpha_qed_running(alpha_0: float, mu_mev: float, m_e_mev: float = M_E_MEV) -> float:
    """
    One-loop QED running coupling.

    α(μ) = α₀ / [1 - (α₀/(3π)) ln(μ²/m_e²)]

    Parameters
    ----------
    alpha_0 : float
        Value of α at the electron mass scale.
    mu_mev : float
        Energy scale μ in MeV.
    m_e_mev : float
        Electron mass in MeV.

    Returns
    -------
    float
        α(μ) at scale μ.
    """
    log_ratio = np.log(mu_mev**2 / m_e_mev**2)
    denom = 1 - (alpha_0 / (3 * np.pi)) * log_ratio
    if denom <= 0:
        return np.inf  # Landau pole
    return alpha_0 / denom


def part_b():
    """Run Part B: QED running analysis."""
    print()
    print()
    print("=" * 78)
    print("PART B: QED RUNNING")
    print("=" * 78)
    print()

    roots = solve_cubic()

    # The paper's α is the smallest positive root
    positive_roots = [float(np.real(r)) for r in roots if np.real(r) > 0]
    alpha_paper = min(positive_roots)

    print(f"Paper's α (smallest positive root) = {alpha_paper:.12f}")
    print(f"Paper's α⁻¹                        = {1/alpha_paper:.10f}")
    print(f"CODATA α⁻¹                         = {1/ALPHA_CODATA:.10f}")
    print(f"Paper vs CODATA relative error      = {abs(alpha_paper - ALPHA_CODATA)/ALPHA_CODATA:.4e}")
    print()

    # Scan μ from 0.1 MeV to 1 TeV
    mu_values = np.logspace(np.log10(0.1), np.log10(1e6), 10000)  # MeV

    alpha_running = np.array([alpha_qed_running(alpha_paper, mu) for mu in mu_values])

    # Find where running α best matches CODATA
    valid = np.isfinite(alpha_running) & (alpha_running > 0)
    residuals = np.abs(alpha_running[valid] - ALPHA_CODATA)
    best_idx_valid = np.argmin(residuals)

    # Map back to full array
    valid_indices = np.where(valid)[0]
    best_idx = valid_indices[best_idx_valid]

    mu_star = mu_values[best_idx]
    alpha_star = alpha_running[best_idx]

    print("─" * 78)
    print("OPTIMAL MATCH TO CODATA")
    print("─" * 78)
    print(f"  μ* = {mu_star:.6f} MeV")
    print(f"  μ*/m_e = {mu_star / M_E_MEV:.6f}")
    print(f"  α(μ*) = {alpha_star:.12f}")
    print(f"  α(μ*)⁻¹ = {1/alpha_star:.10f}")
    print(f"  Relative error vs CODATA = {abs(alpha_star - ALPHA_CODATA)/ALPHA_CODATA:.4e}")
    print()

    # Check specific scales
    print("─" * 78)
    print("α AT RECOGNIZED ENERGY SCALES")
    print("─" * 78)

    scales = {
        "μ = m_e (0.511 MeV)":         M_E_MEV,
        "μ = 2m_e (threshold)":         2 * M_E_MEV,
        "μ = m_μ (105.66 MeV)":         105.6583755,
        "μ = m_τ (1776.9 MeV)":         1776.86,
        "μ = m_W (80379 MeV)":          80379.0,
        "μ = M_Z (91188 MeV)":          91187.6,
        "μ = 0 (Thomson, q²→0)":        None,      # special case
    }

    print(f"  {'Scale':<30s} {'μ (MeV)':>12s} {'α(μ)':>16s} {'α(μ)⁻¹':>14s} {'vs CODATA':>12s}")
    print(f"  {'─'*30} {'─'*12} {'─'*16} {'─'*14} {'─'*12}")

    for name, mu in scales.items():
        if mu is None:
            # Thomson limit: α at q²→0 is just α(m_e) itself
            alpha_val = alpha_paper
            print(f"  {name:<30s} {'→0':>12s} {alpha_val:>16.12f} {1/alpha_val:>14.10f} {abs(alpha_val - ALPHA_CODATA)/ALPHA_CODATA:>11.4e}")
        else:
            alpha_val = alpha_qed_running(alpha_paper, mu)
            if np.isfinite(alpha_val):
                print(f"  {name:<30s} {mu:>12.2f} {alpha_val:>16.12f} {1/alpha_val:>14.10f} {abs(alpha_val - ALPHA_CODATA)/ALPHA_CODATA:>11.4e}")
            else:
                print(f"  {name:<30s} {mu:>12.2f} {'Landau pole':>16s}")

    print()

    # Check α(M_Z) against 1/127.95
    alpha_mz = alpha_qed_running(alpha_paper, 91187.6)
    alpha_mz_expt = 1 / 127.952  # Experimental α(M_Z)
    print("─" * 78)
    print("CHECK: α(M_Z)")
    print("─" * 78)
    print(f"  Paper α run to M_Z (1-loop QED only): α⁻¹ = {1/alpha_mz:.6f}")
    print(f"  Experimental α(M_Z)⁻¹:                      127.952")
    print(f"  Relative error: {abs(alpha_mz - alpha_mz_expt)/alpha_mz_expt:.4e}")
    print(f"  (Note: 1-loop QED only; full SM running includes hadronic + weak)")
    print()

    # Landau pole
    # α₀/(3π) ln(μ²/m_e²) = 1  →  μ_Landau = m_e exp(3π/(2α₀))
    mu_landau = M_E_MEV * np.exp(3 * np.pi / (2 * alpha_paper))
    print("─" * 78)
    print("LANDAU POLE (1-loop)")
    print("─" * 78)
    print(f"  μ_Landau = m_e × exp(3π/(2α)) = {mu_landau:.4e} MeV")
    print(f"           = {mu_landau / 1e3:.4e} GeV")
    print(f"  (Well above Planck scale — QED is perturbative everywhere relevant)")
    print()

    # Summary: what scale does the paper's value correspond to?
    print("─" * 78)
    print("INTERPRETATION")
    print("─" * 78)

    # The paper's value is slightly LARGER than CODATA (α⁻¹ = 137.035987 < 137.035999)
    # Running α increases with μ, so the paper's α corresponds to a scale slightly
    # above m_e if it's meant to match CODATA at q²→0
    if alpha_paper > ALPHA_CODATA:
        print(f"  Paper α⁻¹ = {1/alpha_paper:.6f} < CODATA α⁻¹ = {1/ALPHA_CODATA:.6f}")
        print(f"  Paper α is LARGER than CODATA → corresponds to μ > m_e")
        print(f"  The paper's formula gives α at μ* = {mu_star:.4f} MeV = {mu_star/M_E_MEV:.4f} m_e")
    else:
        print(f"  Paper α⁻¹ = {1/alpha_paper:.6f} > CODATA α⁻¹ = {1/ALPHA_CODATA:.6f}")
        print(f"  Paper α is SMALLER than CODATA → running cannot reach CODATA at higher μ")
        print(f"  The paper's value is very close to the Thomson limit (q²→0)")
        # Check: what if running DOWN from some UV scale?
        # α_UV / [1 - (α_UV/3π) ln(μ²/m_e²)] = α_CODATA at μ = m_e
        # This is already handled by the scan


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    part_a()
    part_b()
