"""Sprint LS-6a: Hydrogen Lamb shift in the Eides §3.2 canonical convention.

Clean A/B re-derivation of the LS-1 one-loop Lamb shift in the standard
Eides-Grotch-Shelyuto (2001) Phys. Rep. 342, 63 §3.2 convention. This is
NOT a recomputation of new physics; it is a convention fix for the +38/45
self-energy constant LS-1 used (which inadvertently subtracts the Uehling
kernel constant 4/15 from the canonical Eides constant 10/9), replacing
it with the canonical Eides §3.2 form.

Hypothesis (from LS-5 scoping)
-------------------------------
LS-5 found that the +29.26 MHz LS-3 residual against experimental 1057.85
MHz decomposes into:
  - +24.7 MHz from the LS-1 SE convention choice (one-loop, NOT physics)
  - +5 MHz from genuine alpha^5 multi-loop ceiling

Re-deriving the same one-loop Lamb shift in the canonical Eides §3.2
convention should bring the predicted Lamb shift to ~1052 MHz, leaving a
residual of ~+5 to +6 MHz that is the genuine alpha^5 multi-loop test
target for the eventual two-loop sprint (LS-7).

What changes from LS-1
----------------------
LS-1 §2.2 used:
    DeltaE_SE(2S_{1/2}) = (alpha^3 Z^4 / pi n^3) [
        (4/3) ln(1/(Z*alpha)^2) - (4/3) ln(k_0/Ry) + 38/45
    ] Ha

LS-6a uses (Eides Eq. 3.32):
    DeltaE_SE(nS_{1/2}) = (alpha^3 Z^4 / pi n^3) [
        (4/3) ln(1/(Z*alpha)^2) - (4/3) ln(k_0/Ry) + 10/9
    ] Ha

The constant difference is +10/9 - 38/45 = +12/45 = +4/15. The 4/15
is precisely the Uehling kernel constant -- LS-1 inadvertently
subtracted it from the canonical Eides 10/9 constant. The correct
canonical form keeps Uehling as a SEPARATE vacuum-polarization
contribution (which it is) and uses 10/9 alone for the SE.

The 10/9 constant in Eides Eq. 3.32 already includes:
  - Bethe constant 5/6 (from the original Bethe 1947 derivation)
  - Karplus-Klein-Darwin contribution 1/3 - 5/30 = 1/6 (combined: 5/6 + 5/18 = ...)
  - The s-state Schwinger AMM contribution at j=1/2 (1/6)
The 10/9 = 5/6 + 5/18 = 15/18 + 5/18 = 20/18 = 10/9. So the 10/9
encodes ALL the one-loop SE finite constants for nS_{1/2} states.
There is NO separate AMM term to add for s-states in the canonical
Eides convention; the AMM is already inside 10/9.

For 2P_{1/2}, the canonical Eides form keeps the textbook -1/6
constant (BS Eq. 21.5, Itzykson-Zuber p. 345), unchanged from LS-1:
    DeltaE_SE(nP_{1/2}) = (alpha^3 Z^4 / pi n^3) [
        -(4/3) ln(k_0(n,1)/Ry) - 1/6
    ] Ha
The -1/6 here is the spin-orbit + AMM combination for j = l - 1/2 = 1/2,
which IS the j-dependent piece (l=1 has no Bethe-Karplus-Klein-Darwin
contact contribution; only AMM through spin-orbit).

Numerical impact (one-loop SE 2S, Z=1, n=2, alpha = 1/137.036):
    common_2S = alpha^3 Z^4 / (pi n^3) Ha
    common_2S * HA_TO_MHZ = 101.733 MHz / dimensionless unit
    Shift = (4/15) * 101.733 = +27.129 MHz

Direction: + raises Lamb shift. Predicted LS-6a Lamb = 1025.06 + 27.13 = 1052.19 MHz.
Residual vs experimental 1057.85 = +5.65 MHz (vs LS-5 prediction +5..+11 MHz).

References
----------
- M. I. Eides, H. Grotch, V. A. Shelyuto, "Theory of light hydrogenic
  bound states", Phys. Rep. 342 (2001) 63. §3.2 Eq. (3.32) for the
  canonical one-loop SE form.
- P. J. Mohr, G. Plunien, G. Soff, Phys. Rep. 293 (1998) 227. §V for the
  same canonical form.
- H. A. Bethe, Phys. Rev. 72 (1947) 339. The original Bethe derivation
  (with constant 5/6 in his Eq. 11; the 10/9 emerges after including
  Karplus-Klein-Darwin).
"""

from __future__ import annotations

import json
import math
from pathlib import Path


# ---------------------------------------------------------------------------
# Physical constants (atomic units = Hartree, a_0)
# ---------------------------------------------------------------------------

# Fine structure constant (CODATA 2018) -- same as LS-1
ALPHA = 1.0 / 137.035999084

# Hartree to MHz conversion
HA_TO_MHZ = 6_579_683_920.502  # MHz / Ha

# Nuclear charge (hydrogen)
Z = 1

# Experimental Lamb shift 2S_{1/2} - 2P_{1/2} in hydrogen (CODATA)
LAMB_EXP_MHZ = 1057.845

# Bethe logarithms (Drake & Swainson 1990) -- same as LS-1
BETHE_LOG_2S = 2.8117698931  # ln k_0(2, 0)
BETHE_LOG_2P = -0.0300167089  # ln k_0(2, 1)


# ---------------------------------------------------------------------------
# Vacuum polarization / Uehling shift  -- IDENTICAL to LS-1
# (the convention issue is in self-energy, not VP)
# ---------------------------------------------------------------------------

def vacuum_polarization_shift_hartree(n: int, l: int, Z: int = 1) -> float:
    """Uehling vacuum polarization shift (UNCHANGED from LS-1).

    DeltaE_VP(nl) = -(4 alpha^3 Z^4) / (15 pi n^3) * delta_{l,0}   [Hartree]

    Note: the 4/15 prefactor here is the same 4/15 = 10/9 - 38/45 that
    LS-1 inadvertently subtracted from the SE 10/9 constant. In the
    canonical Eides convention, Uehling enters here (as a separate VP
    contribution) and the SE constant remains 10/9.
    """
    if l != 0:
        return 0.0
    return -4.0 * ALPHA ** 3 * (Z ** 4) / (15.0 * math.pi * n ** 3)


# ---------------------------------------------------------------------------
# Self-energy in EIDES §3.2 CANONICAL CONVENTION (the LS-6a change)
# ---------------------------------------------------------------------------
#
# The canonical Eides §3.2 form for the one-loop nS_{1/2} self-energy is
# Eq. (3.32):
#
#   DeltaE_SE(nS_{1/2}) = (alpha (Z*alpha)^4 / pi n^3) m_e c^2
#                       * { (4/3) ln[(Z alpha)^{-2}]
#                         - (4/3) ln(k_0(n,0) / Ry)
#                         + 10/9
#                         }
#
# The 10/9 constant is the COMBINED Bethe + Karplus-Klein-Darwin + AMM
# constant for an s-state. There is NO separate AMM term to add.
#
# In atomic units (m_e c^2 = 1/alpha^2 Ha):
#   prefactor_SE  = alpha^3 Z^4 / (pi n^3)  Ha
#
# For nP_{1/2} (l=1, j=1/2), the canonical form is (Itzykson-Zuber p.345):
#
#   DeltaE_SE(nP_{1/2}) = (alpha^3 Z^4 / pi n^3) [
#                            -(4/3) ln(k_0(n,1) / Ry) - 1/6
#                         ] Ha
#
# The -1/6 here is the AMM/spin-orbit combination for j = l - 1/2.


def self_energy_shift_eides_hartree(n: int, l: int, j: float, Z: int = 1,
                                    ln_k0_overRy: float = None) -> dict:
    """One-loop self-energy in Eides §3.2 canonical convention.

    For nS_{1/2}:
      total = (alpha^3 Z^4 / pi n^3) [(4/3) ln[(Z alpha)^{-2}]
                                     - (4/3) ln(k_0(n,0)/Ry) + 10/9]

    For nP_{1/2}:
      total = (alpha^3 Z^4 / pi n^3) [-(4/3) ln(k_0(n,1)/Ry) - 1/6]
    """
    if ln_k0_overRy is None:
        raise ValueError("Must supply Bethe logarithm")

    alpha = ALPHA
    Za = Z * alpha
    n3 = n ** 3
    common = alpha ** 3 * (Z ** 4) / (math.pi * n3)

    if l == 0 and abs(j - 0.5) < 1e-9:
        # nS_{1/2}: canonical Eides Eq. 3.32 form
        # Bethe-log piece + 10/9 combined SE constant (NO separate AMM)
        bethe_log_piece = (4.0 / 3.0) * (
            math.log(1.0 / (Za ** 2)) - ln_k0_overRy
        )
        constant_piece = 10.0 / 9.0
        bracket = bethe_log_piece + constant_piece
        total = common * bracket
        return {
            "bethe_log_piece_Ha": common * bethe_log_piece,
            "constant_piece_Ha": common * constant_piece,
            "total_Ha": total,
            "convention": "Eides 3.32 (nS_{1/2}); 10/9 includes Bethe+KK+Darwin+AMM",
            "constant": "10/9",
        }

    elif l == 1 and abs(j - 0.5) < 1e-9 and n == 2:
        # 2P_{1/2}: Itzykson-Zuber form, unchanged from LS-1
        bethe_log_piece = -(4.0 / 3.0) * ln_k0_overRy
        constant_piece = -1.0 / 6.0
        bracket = bethe_log_piece + constant_piece
        total = common * bracket
        return {
            "bethe_log_piece_Ha": common * bethe_log_piece,
            "constant_piece_Ha": common * constant_piece,
            "total_Ha": total,
            "convention": "Itzykson-Zuber p.345 (2P_{1/2}); -1/6 is AMM/SO",
            "constant": "-1/6",
        }

    else:
        raise NotImplementedError(
            f"Eides convention only implemented for nS_{{1/2}} and 2P_{{1/2}}; "
            f"got n={n}, l={l}, j={j}"
        )


# ---------------------------------------------------------------------------
# LS-1 self-energy formula (for direct A/B comparison)
# ---------------------------------------------------------------------------

def self_energy_shift_ls1_hartree(n: int, l: int, j: float, Z: int = 1,
                                  ln_k0_overRy: float = None) -> float:
    """LS-1 self-energy (the +38/45 lumped form) for direct A/B comparison.

    For 2S_{1/2}:
      DeltaE_SE = (alpha^3 Z^4 / pi n^3) [(4/3) ln(1/(Z alpha)^2)
                                         - (4/3) ln(k_0/Ry) + 38/45]
    For 2P_{1/2}:
      DeltaE_SE = (alpha^3 Z^4 / pi n^3) [-(4/3) ln(k_0/Ry) - 1/6]
    """
    if ln_k0_overRy is None:
        raise ValueError("Must supply Bethe log")

    alpha = ALPHA
    Za = Z * alpha
    common = alpha ** 3 * (Z ** 4) / (math.pi * n ** 3)

    if l == 0:
        bracket = (4.0 / 3.0) * (math.log(1.0 / (Za ** 2)) - ln_k0_overRy) + 38.0 / 45.0
        return common * bracket
    elif l == 1 and n == 2:
        bracket = -(4.0 / 3.0) * ln_k0_overRy - 1.0 / 6.0
        return common * bracket
    else:
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Main computation: clean A/B between LS-1 and LS-6a (Eides)
# ---------------------------------------------------------------------------

def compute_lamb_shift_ab() -> dict:
    """Compute the Lamb shift in BOTH conventions for direct A/B."""
    # Vacuum polarization (identical in both)
    VP_2S = vacuum_polarization_shift_hartree(n=2, l=0, Z=Z)
    VP_2P = vacuum_polarization_shift_hartree(n=2, l=1, Z=Z)

    # ---- LS-1 self-energy (the +38/45 lumped form) ----
    SE_2S_ls1 = self_energy_shift_ls1_hartree(
        n=2, l=0, j=0.5, Z=Z, ln_k0_overRy=BETHE_LOG_2S)
    SE_2P_ls1 = self_energy_shift_ls1_hartree(
        n=2, l=1, j=0.5, Z=Z, ln_k0_overRy=BETHE_LOG_2P)

    total_2S_ls1 = SE_2S_ls1 + VP_2S
    total_2P_ls1 = SE_2P_ls1 + VP_2P
    lamb_ls1 = (total_2S_ls1 - total_2P_ls1) * HA_TO_MHZ

    # ---- LS-6a self-energy (Eides §3.2 canonical convention) ----
    SE_2S_eides = self_energy_shift_eides_hartree(
        n=2, l=0, j=0.5, Z=Z, ln_k0_overRy=BETHE_LOG_2S)
    SE_2P_eides = self_energy_shift_eides_hartree(
        n=2, l=1, j=0.5, Z=Z, ln_k0_overRy=BETHE_LOG_2P)

    total_2S_eides = SE_2S_eides["total_Ha"] + VP_2S
    total_2P_eides = SE_2P_eides["total_Ha"] + VP_2P
    lamb_eides = (total_2S_eides - total_2P_eides) * HA_TO_MHZ

    # ---- A/B comparison ----
    # The constant difference is purely 10/9 - 38/45 = 12/45 = 4/15
    se_2s_diff_const = 10.0 / 9.0 - 38.0 / 45.0  # = 4/15
    common_2S = ALPHA ** 3 * Z ** 4 / (math.pi * 8)
    expected_se_shift_2S = se_2s_diff_const * common_2S * HA_TO_MHZ
    actual_se_shift_2S = (SE_2S_eides["total_Ha"] - SE_2S_ls1) * HA_TO_MHZ

    # ---- Residuals ----
    err_ls1_mhz = lamb_ls1 - LAMB_EXP_MHZ
    err_eides_mhz = lamb_eides - LAMB_EXP_MHZ

    err_ls1_pct = 100.0 * err_ls1_mhz / LAMB_EXP_MHZ
    err_eides_pct = 100.0 * err_eides_mhz / LAMB_EXP_MHZ

    # Verdict: LS-5 reframing PASS if Lamb in [1051, 1057] AND residual in [3.5, 10.7]
    genuine_residual = -err_eides_mhz
    in_lamb_band = 1051 <= lamb_eides <= 1057
    in_residual_band = 3.5 <= genuine_residual <= 10.7
    if in_lamb_band and in_residual_band:
        verdict = "POSITIVE -- LS-5 reframing CONFIRMED"
    elif 1048 <= lamb_eides <= 1060:
        verdict = "POSITIVE PARTIAL -- LS-5 reframing approximately confirmed"
    else:
        verdict = "NEGATIVE -- LS-5 reframing not confirmed by computation"

    result = {
        "experimental_MHz": LAMB_EXP_MHZ,
        "alpha": ALPHA,
        "Z": Z,
        "n": 2,
        "ha_to_mhz": HA_TO_MHZ,

        # LS-1 (lumped +38/45 convention)
        "ls1": {
            "convention": "+38/45 lumped form (LS-1 §2.2)",
            "SE_2S_MHz": SE_2S_ls1 * HA_TO_MHZ,
            "SE_2P_MHz": SE_2P_ls1 * HA_TO_MHZ,
            "VP_2S_MHz": VP_2S * HA_TO_MHZ,
            "VP_2P_MHz": VP_2P * HA_TO_MHZ,
            "total_2S_MHz": total_2S_ls1 * HA_TO_MHZ,
            "total_2P_MHz": total_2P_ls1 * HA_TO_MHZ,
            "lamb_shift_MHz": lamb_ls1,
            "error_MHz": err_ls1_mhz,
            "error_pct": err_ls1_pct,
        },

        # LS-6a (Eides §3.2 canonical convention)
        "ls6a_eides": {
            "convention": "Eides 3.32 canonical (10/9 for nS_{1/2}, no separate AMM)",
            # 2S breakdown
            "SE_2S_bethe_log_MHz": SE_2S_eides["bethe_log_piece_Ha"] * HA_TO_MHZ,
            "SE_2S_constant_MHz":  SE_2S_eides["constant_piece_Ha"] * HA_TO_MHZ,
            "SE_2S_total_MHz":     SE_2S_eides["total_Ha"] * HA_TO_MHZ,
            "SE_2S_constant_used": SE_2S_eides["constant"],
            # 2P breakdown
            "SE_2P_bethe_log_MHz": SE_2P_eides["bethe_log_piece_Ha"] * HA_TO_MHZ,
            "SE_2P_constant_MHz":  SE_2P_eides["constant_piece_Ha"] * HA_TO_MHZ,
            "SE_2P_total_MHz":     SE_2P_eides["total_Ha"] * HA_TO_MHZ,
            "SE_2P_constant_used": SE_2P_eides["constant"],
            # Vacuum pol (same as LS-1)
            "VP_2S_MHz": VP_2S * HA_TO_MHZ,
            "VP_2P_MHz": VP_2P * HA_TO_MHZ,
            # Totals
            "total_2S_MHz": total_2S_eides * HA_TO_MHZ,
            "total_2P_MHz": total_2P_eides * HA_TO_MHZ,
            "lamb_shift_MHz": lamb_eides,
            "error_MHz": err_eides_mhz,
            "error_pct": err_eides_pct,
        },

        # Convention shift A/B
        "ab_comparison": {
            "shift_lamb_eides_minus_ls1_MHz": lamb_eides - lamb_ls1,
            "shift_se_2s_eides_minus_ls1_MHz": actual_se_shift_2S,
            "constant_diff_2S_eides_vs_ls1": "10/9 - 38/45 = 12/45 = 4/15",
            "constant_diff_2S_value": se_2s_diff_const,
            "expected_se_2s_shift_MHz_from_const": expected_se_shift_2S,
            "shift_matches_pure_const_diff":
                abs(actual_se_shift_2S - expected_se_shift_2S) < 1e-6,
            "interpretation": (
                "The +4/15 = 10/9 - 38/45 is exactly the Uehling kernel constant. "
                "LS-1 inadvertently subtracted Uehling from the canonical 10/9; "
                "Eides convention keeps Uehling separate (as VP) and uses 10/9 for SE."
            ),
        },

        # LS-5 reframing test
        "ls5_test": {
            "ls5_predicted_genuine_multi_loop_MHz_central": 7.10,
            "ls5_predicted_alpha5_band_MHz": "+3.5 to +10.7",
            "ls5_predicted_se_convention_shift_MHz": 24.7,
            "ls6a_actual_se_convention_shift_MHz": lamb_eides - lamb_ls1,
            "ls6a_actual_genuine_residual_MHz": genuine_residual,
            "ls5_reframing_confirmed": in_residual_band,
            "ls6a_in_target_lamb_band_1051_1057": in_lamb_band,
        },

        "verdict": verdict,
    }

    return result


def print_result(r: dict) -> None:
    """Print the A/B comparison."""
    print("=" * 78)
    print("Sprint LS-6a: Hydrogen Lamb shift in Eides §3.2 convention")
    print("=" * 78)
    print(f"\nFine structure alpha = {r['alpha']:.10e}")
    print(f"Z = {r['Z']}, n = {r['n']}")
    print(f"Experimental Lamb shift = {r['experimental_MHz']:.3f} MHz\n")

    print("-" * 78)
    print("LS-1 (lumped +38/45 convention)")
    print("-" * 78)
    ls1 = r["ls1"]
    print(f"  SE 2S      = {ls1['SE_2S_MHz']:+12.4f} MHz")
    print(f"  VP 2S      = {ls1['VP_2S_MHz']:+12.4f} MHz")
    print(f"  Total 2S   = {ls1['total_2S_MHz']:+12.4f} MHz")
    print(f"  SE 2P      = {ls1['SE_2P_MHz']:+12.4f} MHz")
    print(f"  VP 2P      = {ls1['VP_2P_MHz']:+12.4f} MHz")
    print(f"  Total 2P   = {ls1['total_2P_MHz']:+12.4f} MHz")
    print(f"  Lamb shift = {ls1['lamb_shift_MHz']:+12.4f} MHz")
    print(f"  Error      = {ls1['error_MHz']:+12.4f} MHz "
          f"({ls1['error_pct']:+.3f}%)")

    print()
    print("-" * 78)
    print("LS-6a (Eides §3.2 canonical: 10/9 for nS, -1/6 for 2P)")
    print("-" * 78)
    e = r["ls6a_eides"]
    print(f"  SE 2S Bethe log = {e['SE_2S_bethe_log_MHz']:+12.4f} MHz")
    print(f"  SE 2S const     = {e['SE_2S_constant_MHz']:+12.4f} MHz "
          f"(uses {e['SE_2S_constant_used']})")
    print(f"  SE 2S total     = {e['SE_2S_total_MHz']:+12.4f} MHz")
    print(f"  VP 2S           = {e['VP_2S_MHz']:+12.4f} MHz")
    print(f"  Total 2S        = {e['total_2S_MHz']:+12.4f} MHz")
    print()
    print(f"  SE 2P Bethe log = {e['SE_2P_bethe_log_MHz']:+12.4f} MHz")
    print(f"  SE 2P const     = {e['SE_2P_constant_MHz']:+12.4f} MHz "
          f"(uses {e['SE_2P_constant_used']})")
    print(f"  SE 2P total     = {e['SE_2P_total_MHz']:+12.4f} MHz")
    print(f"  VP 2P           = {e['VP_2P_MHz']:+12.4f} MHz")
    print(f"  Total 2P        = {e['total_2P_MHz']:+12.4f} MHz")
    print()
    print(f"  Lamb shift  = {e['lamb_shift_MHz']:+12.4f} MHz")
    print(f"  Error       = {e['error_MHz']:+12.4f} MHz "
          f"({e['error_pct']:+.3f}%)")

    print()
    print("-" * 78)
    print("A/B comparison: Eides - LS-1")
    print("-" * 78)
    ab = r["ab_comparison"]
    print(f"  Lamb shift change: {ab['shift_lamb_eides_minus_ls1_MHz']:+10.4f} MHz")
    print(f"  SE 2S change:      {ab['shift_se_2s_eides_minus_ls1_MHz']:+10.4f} MHz")
    print(f"  Constant difference: {ab['constant_diff_2S_eides_vs_ls1']}")
    print(f"  Predicted from constant only: "
          f"{ab['expected_se_2s_shift_MHz_from_const']:+10.4f} MHz")
    print(f"  Pure constant difference accounts for shift: "
          f"{ab['shift_matches_pure_const_diff']}")
    print(f"  Interpretation: {ab['interpretation']}")

    print()
    print("-" * 78)
    print("LS-5 reframing test")
    print("-" * 78)
    test = r["ls5_test"]
    print(f"  LS-5 predicted convention shift   : "
          f"{test['ls5_predicted_se_convention_shift_MHz']:+8.2f} MHz")
    print(f"  LS-6a actual convention shift     : "
          f"{test['ls6a_actual_se_convention_shift_MHz']:+8.2f} MHz")
    print(f"  LS-5 predicted alpha^5 band       : "
          f"{test['ls5_predicted_alpha5_band_MHz']} MHz")
    print(f"  LS-5 predicted central value      : "
          f"{test['ls5_predicted_genuine_multi_loop_MHz_central']:+8.2f} MHz")
    print(f"  LS-6a actual genuine residual     : "
          f"{test['ls6a_actual_genuine_residual_MHz']:+8.2f} MHz")
    print(f"  In target Lamb band [1051, 1057]  : "
          f"{test['ls6a_in_target_lamb_band_1051_1057']}")
    print(f"  LS-5 reframing CONFIRMED?         : "
          f"{test['ls5_reframing_confirmed']}")

    print()
    print(f"VERDICT: {r['verdict']}")
    print()


def save_data(result: dict, path: Path) -> None:
    """Save result to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    serializable = json.loads(json.dumps(result, default=str))
    with open(path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"Saved: {path}")


def main() -> dict:
    """Run the LS-6a A/B comparison."""
    print("Sprint LS-6a: Eides §3.2 canonical convention")
    print()

    result = compute_lamb_shift_ab()
    print_result(result)

    out_path = Path(__file__).parent / "data" / "ls6a_eides_convention.json"
    save_data(result, out_path)

    return result


if __name__ == "__main__":
    main()
