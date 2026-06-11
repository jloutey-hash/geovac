"""BR-C: He 2^3P fine-structure multiplet benchmark via Breit-Pauli perturbation theory.

Goal
----
Compute the He 2^3P_J (J = 0, 1, 2) multiplet energies using Breit-Pauli
perturbation theory with spin-orbit, spin-spin tensor, and spin-other-orbit
corrections. Compare the multiplet splittings against NIST (Drake 2006).

Reference values (Drake 2006 / CODATA, converted from MHz):
    E(2^3P_0) - E(2^3P_1) =  +29,616.951 MHz  =  +4.501e-9 Ha (INVERTED multiplet!)
    E(2^3P_1) - E(2^3P_2) =   +2,291.178 MHz  =  +3.482e-10 Ha
    full span                 31,908.129 MHz  =   4.850e-9 Ha

Note the J ordering of the 2^3P multiplet is INVERTED relative to the normal
Lande rule: E(J=0) > E(J=1) > E(J=2). This inversion is driven by the two-body
spin-spin tensor interaction, which dominates over the one-body spin-orbit for
the (1s)(np) configuration in helium.

T8 baseline
-----------
Pure single-particle H_SO gives the NORMAL Lande ordering E(J=0) < E(J=1) <
E(J=2) -- i.e. SIGN WRONG relative to NIST -- with magnitude alpha^2 Z^4/32
for the doublet splitting. For He 2^3P this is ~10,949 MHz (66% of NIST span
magnitude, wrong sign).

Approach
--------
1. Use the closed-form LS-coupled matrix elements for (1s)(np) ^3P_J
   from the Breit-Pauli formalism (Bethe-Salpeter §39; Drake 1971; Cowan 1981).
2. Evaluate the three parameters (zeta, M, N) using GeoVac T1/BR-B infrastructure.
3. Compare with NIST.

BR-B limitation (flagged):
The BR-B `compute_rk_breit_retarded_algebraic` returns 0 for many convergent
integrals -- the region-splitting in `_T_kernel_breit_retarded` fails when
m1 < 0 AND m2 < 0 even though the combined integral should converge. This is
a separate issue from the present benchmark. We fall back to direct sympy
integration of the Breit radial Slater integrals here, and document the
discrepancy.

Outputs
-------
    debug/data/br_c_2P_benchmark.json  -- computed splittings and breakdown
    debug/br_c_he_fine_structure_memo.md  -- human-readable writeup

Author
------
GeoVac Development Team, April 2026
Sprint 2 Track BR-C
"""

from __future__ import annotations

import json
import math
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Rational, sqrt, simplify, oo, integrate, exp, symbols, Float

# Ensure UTF-8 output on Windows
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

# Project imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.spin_orbit import so_diagonal_matrix_element
from geovac.dirac_matrix_elements import Z_sym, alpha_sym


# ===========================================================================
# Constants
# ===========================================================================

ALPHA_CODATA = 7.2973525693e-3
HA_TO_MHZ = 6.5796839204e9  # 1 Ha -> MHz (frequency)

# NIST / Drake 2006 reference values for He 2^3P multiplet
# Source: G.W.F. Drake, Phys. Scr. T120 (2005); CODATA; NIST ASD.
# The He 2^3P multiplet is INVERTED: E(J=0) > E(J=1) > E(J=2)
NIST_REF_MHZ = {
    "E(P0) - E(P1)": 29616.951,
    "E(P1) - E(P2)":  2291.178,
    "E(P0) - E(P2)": 29616.951 + 2291.178,
}

NIST_REF_HA = {k: v / HA_TO_MHZ for k, v in NIST_REF_MHZ.items()}


# ===========================================================================
# Hydrogenic orbital radial wavefunctions (sympy exact)
# ===========================================================================

def hydrogenic_R(n: int, l: int, Z_sym_eff, r_sym) -> sp.Expr:
    """Normalized hydrogenic radial wavefunction R_{n,l}(r) at charge Z.

    Returns R(r) such that integral_0^inf R^2 r^2 dr = 1.
    """
    from sympy import sqrt, factorial, assoc_laguerre, exp, Rational

    # Standard hydrogenic normalization
    # R_{n,l}(r) = sqrt[(2Z/n)^3 (n-l-1)!/(2n (n+l)!)]
    #              * exp(-Zr/n) (2Zr/n)^l L_{n-l-1}^{2l+1}(2Zr/n)
    rho = 2 * Z_sym_eff * r_sym / n
    norm = sqrt((2 * Z_sym_eff / n) ** 3 * factorial(n - l - 1)
                / (2 * n * factorial(n + l)))
    return norm * exp(-rho / 2) * rho ** l * assoc_laguerre(n - l - 1, 2 * l + 1, rho)


def breit_retarded_slater(
    n1: int, l1: int, n2: int, l2: int,
    n3: int, l3: int, n4: int, l4: int,
    l: int, Z: int,
) -> sp.Rational:
    """Compute the Breit-Pauli retarded Slater integral via direct sympy integration:

        R_BP^l = int_0^inf int_0^inf R_{n1 l1}(r1) R_{n3 l3}(r1)
                                      R_{n2 l2}(r2) R_{n4 l4}(r2)
                                      (r_<^l / r_>^(l+3)) r1^2 r2^2 dr1 dr2

    Expressed as exact sympy rational.

    Uses the region-splitting r1 < r2 and r1 > r2.
    """
    r1, r2 = symbols('r1 r2', positive=True, real=True)
    R13 = hydrogenic_R(n1, l1, Integer(Z), r1) * hydrogenic_R(n3, l3, Integer(Z), r1)
    R24 = hydrogenic_R(n2, l2, Integer(Z), r2) * hydrogenic_R(n4, l4, Integer(Z), r2)

    # Region I: r1 <= r2,  r_< = r1, r_> = r2
    # Integrand I: R13(r1) R24(r2) r1^l / r2^(l+3) r1^2 r2^2
    integrand_I = R13 * R24 * r1 ** l / r2 ** (l + 3) * r1 ** 2 * r2 ** 2
    # integrate r1 from 0 to r2, then r2 from 0 to oo
    inner_I = integrate(integrand_I, (r1, 0, r2))
    outer_I = integrate(inner_I, (r2, 0, oo))

    # Region II: r2 <= r1
    integrand_II = R13 * R24 * r2 ** l / r1 ** (l + 3) * r1 ** 2 * r2 ** 2
    inner_II = integrate(integrand_II, (r2, 0, r1))
    outer_II = integrate(inner_II, (r1, 0, oo))

    total = sp.simplify(outer_I + outer_II)
    return total


# ===========================================================================
# Closed-form He 1s 2p ^3P Breit-Pauli matrix elements (Bethe-Salpeter §39)
# ===========================================================================
#
# Reference: H.A. Bethe & E.E. Salpeter, "Quantum Mechanics of One- and Two-Electron
# Atoms" (1957), §39 "Fine Structure of Helium".
#
# For the 1s np ^3P_J states of He, the Breit-Pauli fine-structure operator in
# first-order perturbation theory has diagonal elements:
#
#   E(^3P_J) = E_SO(J) + E_SS(J) + E_SOO(J)
#
# where (Eqs. 39.12, 39.16, 39.17 in Bethe-Salpeter, or Drake 1971 Tables I,II):
#
#   E_SO(^3P_J)  = (zeta/2) * X(J)    with X(J) = J(J+1) - L(L+1) - S(S+1)
#                                         = J(J+1) - 4
#   E_SS(^3P_J)  = A_SS * f_SS(J)
#   E_SOO(^3P_J) = A_SOO * f_SOO(J)
#
# where for L=1, S=1:
#
#   X(J=0) = -4,  X(J=1) = -2,  X(J=2) = +2       (note X/2 = <L.S>_J)
#
# The two-body angular coefficients are DIFFERENT from Lande:
# Bethe-Salpeter §39 gives (for 1snp ^3P, standard convention):
#
#   f_SS(J=0)  = -2,    f_SS(J=1)  = +1,     f_SS(J=2)  = -1/5
#   f_SOO(J=0) = +2,    f_SOO(J=1) = +1,     f_SOO(J=2) = -1
#
# The radial amplitudes A_SS and A_SOO for (1s)(np) (Bethe-Salpeter eq. 39.14-15):
#
#   A_SS  = (3/2) alpha^2 * M_ret            where M_ret = <1s np | r_<^2/r_>^5 | 1s np>
#   A_SOO = (1/2) alpha^2 * [M_ret + Delta]  where Delta is a related integral
#
# The one-body SO parameter zeta_{np}, using hydrogenic orbitals with effective
# charge Z_eff:
#
#   zeta_{np} = alpha^2 * Z_nuc * <1/r^3>_{np}
#             = alpha^2 * Z_nuc * Z_eff^3 / (n^3 l (l+1/2) (l+1))
#             = alpha^2 * Z_nuc * Z_eff^3 / (3 n^3)              (for l=1)
#
# For He at Z_nuc=2 with Z_eff=1 for the unscreened 2p electron:
#   zeta_{2p} = alpha^2 * 2 * 1 / 24 = alpha^2 / 12
#
# In atomic units:
#   alpha^2 = 5.326e-5
#   zeta_{2p} (bare Z_eff=1) = 4.44e-6 Ha ~ 29.2 GHz
#
# ---------------------------------------------------------------------------
# ANGULAR-FACTOR NOTE
# ---------------------------------------------------------------------------
# The coefficients (f_SS, f_SOO) above come from Bethe-Salpeter §39 with a
# specific sign convention.  We VERIFY them numerically against the expected
# He 2^3P inversion pattern (E(P_0) > E(P_1) > E(P_2)) and document any sign
# flips needed for this convention.
# ---------------------------------------------------------------------------


# Bethe-Salpeter §39 J-dependent coefficients for (1s)(np) ^3P
# Format: f_SS[J], f_SOO[J] for J in {0, 1, 2}
F_SS_BS = {0: -2.0, 1: +1.0, 2: -Fraction(1, 5)}
F_SOO_BS = {0: +2.0, 1: +1.0, 2: -1.0}


def he_2P_matrix_elements(
    Z_nuc: int = 2,
    Z_eff_2p: float = 1.0,
    alpha: float = ALPHA_CODATA,
    use_direct_sympy: bool = True,
) -> Dict[str, float]:
    """Closed-form Breit-Pauli matrix elements for He 1s 2p ^3P_J.

    Parameters
    ----------
    Z_nuc : int
        Nuclear charge Z for the -Z/r potential (=2 for He).
    Z_eff_2p : float
        Effective charge for the 2p hydrogenic orbital (1s is at Z_nuc).
        Captures screening of the 2p by the 1s core.
    alpha : float
        Fine-structure constant.
    use_direct_sympy : bool
        If True, use direct sympy integration for the Breit-Pauli retarded
        Slater integrals (robust). If False, try the BR-B machinery (currently
        buggy for some cases — see BR-B divergence in region splitting).
    """
    # ---------- zeta_{2p}: one-body SO parameter ----------
    r3_inv_2p = (Z_eff_2p ** 3) / (2 ** 3 * 1 * 1.5 * 2)  # <1/r^3>_{2p}
    zeta_2p = alpha ** 2 * Z_nuc * r3_inv_2p

    # ---------- BP-retarded Slater integrals ----------
    # For the (1s)(2p) ^3P configuration, the SS tensor involves the l=2
    # partial wave of (1/r_12^3) (P_2(costheta) / r_12^3), which after
    # Breit-Pauli retardation has radial kernel r_<^2 / r_>^5.
    # The SOO mixes l=0 and l=2 partial waves.
    #
    # The "1s 2p" orbitals are hydrogenic at Z=Z_nuc for 1s, and
    # ... for the 2p, we need to account for screening via Z_eff_2p.
    # For simplicity here we use hydrogenic at Z_nuc for both and
    # separately note the Z_eff_2p correction applies as ~(Z_eff/Z_nuc)^3
    # multiplier to the 2p radial moments.

    # Direct sympy integration of the Breit-Pauli retarded integrals.
    # M_ret_l = <1s 2p | (r_<^l / r_>^(l+3)) | 1s 2p> with the "1s" on electron 1
    #                                                    and "2p" on electron 2
    if use_direct_sympy:
        # Use hydrogenic orbitals for both 1s (Z=Z_nuc) and 2p (Z=Z_eff_2p)
        # NOTE: mixing Z for different electrons requires matching the orbital
        # definitions; for the benchmark we use Z=Z_nuc for both and apply
        # Z_eff correction only to zeta_{2p} above.
        try:
            M2_ret = breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 2, Z_nuc)
            M2_ret_f = float(M2_ret)
            print(f"    M2_ret (l=2, direct sympy) = {M2_ret} = {M2_ret_f:.6e}")
        except Exception as e:
            print(f"    [warn] M2_ret direct integration failed: {e}")
            M2_ret_f = 0.0
        try:
            M0_ret = breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 0, Z_nuc)
            M0_ret_f = float(M0_ret)
            print(f"    M0_ret (l=0, direct sympy) = {M0_ret} = {M0_ret_f:.6e}")
        except Exception as e:
            print(f"    [warn] M0_ret direct integration failed: {e}")
            M0_ret_f = 0.0
    else:
        M2_ret_f = 0.0
        M0_ret_f = 0.0

    # ---------- SS and SOO amplitudes (Bethe-Salpeter §39.14-15) ----------
    # A_SS = (3/2) alpha^2 * M_ret(l=2)       (standard Bethe-Salpeter)
    # A_SOO contains both l=0 and l=2 contributions
    A_SS = 1.5 * alpha ** 2 * M2_ret_f
    A_SOO = 0.5 * alpha ** 2 * (M2_ret_f - 2 * M0_ret_f)   # approximate B-S form

    # ---------- Assemble E(^3P_J) ----------
    E_SO_J = {J: (zeta_2p / 2.0) * (J * (J + 1) - 4) / 1.0
              for J in [0, 1, 2]}
    # NOTE: factor 1/2 is included in zeta/2*X definition where X = J(J+1)-4;
    # X/2 = <L.S>, so zeta/2*X = zeta*<L.S>. We use E_SO = zeta/2 * X as per
    # Bethe-Salpeter eq. 39.12.

    E_SS_J = {J: A_SS * float(F_SS_BS[J]) for J in [0, 1, 2]}
    E_SOO_J = {J: A_SOO * float(F_SOO_BS[J]) for J in [0, 1, 2]}

    E_total_J = {J: E_SO_J[J] + E_SS_J[J] + E_SOO_J[J]
                  for J in [0, 1, 2]}

    # ---------- Splittings ----------
    split_ha = {
        "E(P0) - E(P1)": E_total_J[0] - E_total_J[1],
        "E(P1) - E(P2)": E_total_J[1] - E_total_J[2],
        "E(P0) - E(P2)": E_total_J[0] - E_total_J[2],
    }
    split_mhz = {k: v * HA_TO_MHZ for k, v in split_ha.items()}

    rel_err = {k: (split_mhz[k] - NIST_REF_MHZ[k]) / NIST_REF_MHZ[k]
                for k in NIST_REF_MHZ}

    # ---------- SO-only baseline ----------
    SO_only_ha = {
        "E(P0) - E(P1)": E_SO_J[0] - E_SO_J[1],
        "E(P1) - E(P2)": E_SO_J[1] - E_SO_J[2],
        "E(P0) - E(P2)": E_SO_J[0] - E_SO_J[2],
    }
    SO_only_mhz = {k: v * HA_TO_MHZ for k, v in SO_only_ha.items()}
    SO_only_rel_err = {k: (SO_only_mhz[k] - NIST_REF_MHZ[k]) / NIST_REF_MHZ[k]
                         for k in NIST_REF_MHZ}

    return {
        "Z_nuc": Z_nuc,
        "Z_eff_2p": Z_eff_2p,
        "alpha": alpha,
        "zeta_2p_Ha": zeta_2p,
        "zeta_2p_MHz": zeta_2p * HA_TO_MHZ,
        "M2_retarded": M2_ret_f,
        "M0_retarded": M0_ret_f,
        "A_SS": A_SS,
        "A_SOO": A_SOO,
        "E_SO_J_Ha": E_SO_J,
        "E_SS_J_Ha": E_SS_J,
        "E_SOO_J_Ha": E_SOO_J,
        "E_total_J_Ha": E_total_J,
        "splittings_computed_Ha": split_ha,
        "splittings_computed_MHz": split_mhz,
        "NIST_reference_MHz": NIST_REF_MHZ,
        "rel_errors": rel_err,
        "SO_only_baseline_MHz": SO_only_mhz,
        "SO_only_rel_errors": SO_only_rel_err,
    }


# ===========================================================================
# Sanity check: verify f_SS, f_SOO vs published He ^3P tables
# ===========================================================================

def verify_angular_coefficients() -> Dict[str, any]:
    """Verify the Breit-Pauli J-coefficients for He ^3P against published values.

    Bethe-Salpeter §39.12-17 for (1s)(np) ^3P:
      L=1, S=1, so X(J) = J(J+1) - 4:  -4, -2, +2
      f_SS(J):   -2,    +1,    -1/5
      f_SOO(J):  +2,    +1,    -1

    Cowan (1981) Table 9.10 "f-factors for LS-coupled tensors":
      Same values up to sign convention.

    We also check the constraint that the three splittings predicted from
    these coefficients, with *known* ratio 29617:2291 in He, consistent with
    A_SS, A_SOO being both positive and with A_SS >> zeta ~ A_SOO.
    """
    out = {}
    # Check: for zero zeta and zero SOO, pure SS should invert the multiplet
    # with ratio SS:(J=0)/(J=2) = -2 / (-1/5) = +10
    # Observed NIST: (P0-P2) = 31908 = sum of (P0-P1) + (P1-P2)
    #               (P0-P1)/(P1-P2) = 29617/2291 ~ 12.9  (close to 13)
    # This is consistent with the Breit-Pauli Lande-plus-inversion pattern
    # where SS dominates and inverts the multiplet.
    out["X_of_J"] = {J: J * (J + 1) - 4 for J in [0, 1, 2]}
    out["f_SS"] = {J: float(F_SS_BS[J]) for J in [0, 1, 2]}
    out["f_SOO"] = {J: float(F_SOO_BS[J]) for J in [0, 1, 2]}
    # Ratio test for SS only
    f_ratio_SS = (F_SS_BS[0] - F_SS_BS[1]) / (F_SS_BS[1] - F_SS_BS[2])
    out["SS_only_splitting_ratio_(P0-P1)/(P1-P2)"] = float(f_ratio_SS)
    f_ratio_SOO = (F_SOO_BS[0] - F_SOO_BS[1]) / (F_SOO_BS[1] - F_SOO_BS[2])
    out["SOO_only_splitting_ratio_(P0-P1)/(P1-P2)"] = float(f_ratio_SOO)
    out["NIST_splitting_ratio_(P0-P1)/(P1-P2)"] = NIST_REF_MHZ["E(P0) - E(P1)"] / NIST_REF_MHZ["E(P1) - E(P2)"]
    # ratio SS = (-2-1)/(1-(-1/5)) = -3 / 1.2 = -2.5
    # ratio SOO = (2-1)/(1-(-1)) = 1/2 = 0.5
    # NIST = 12.9.  Neither pure SS nor pure SOO reproduces 12.9, confirming
    # that the three-parameter fit (zeta, SS, SOO) is needed.
    return out


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 78)
    print("Track BR-C: He 2^3P Fine-Structure Multiplet Benchmark")
    print("=" * 78)

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Verify angular coefficients
    print("\n--- Breit-Pauli angular coefficients (Bethe-Salpeter §39) ---")
    angular = verify_angular_coefficients()
    print(f"  X(J):     {angular['X_of_J']}")
    print(f"  f_SS(J):  {angular['f_SS']}")
    print(f"  f_SOO(J): {angular['f_SOO']}")
    print(f"  SS-only ratio (P0-P1)/(P1-P2):   {angular['SS_only_splitting_ratio_(P0-P1)/(P1-P2)']:+.2f}")
    print(f"  SOO-only ratio (P0-P1)/(P1-P2):  {angular['SOO_only_splitting_ratio_(P0-P1)/(P1-P2)']:+.2f}")
    print(f"  NIST ratio (P0-P1)/(P1-P2):       {angular['NIST_splitting_ratio_(P0-P1)/(P1-P2)']:+.2f}")

    # Run computation at Z_eff=1.0 (standard approximation for excited 2p)
    print("\n--- Computation at Z_nuc=2, Z_eff_2p=1.0 (standard screening) ---")
    r1 = he_2P_matrix_elements(Z_nuc=2, Z_eff_2p=1.0)
    _print_result(r1)

    # Run at Z_eff=1.34 (Clementi-Raimondi for He 1s2p ^3P)
    print("\n--- Computation at Z_nuc=2, Z_eff_2p=1.34 (Clementi-Raimondi) ---")
    r2 = he_2P_matrix_elements(Z_nuc=2, Z_eff_2p=1.34)
    _print_result(r2)

    # Try with different Breit amplitude scaling
    # (to test whether the J-pattern is correct even if overall scale is off)
    print("\n--- Sensitivity check: amplitude scan for (zeta_SO, A_SS, A_SOO) ---")
    sensitivity_results = amplitude_sensitivity_scan(Z_nuc=2, Z_eff_2p=1.0)

    # Alternative convention test
    print("\n--- Alternative sign convention test (Drake 1971 / Johnson) ---")
    alt_results = try_alternative_sign_conventions(Z_nuc=2, Z_eff_2p=1.0)

    # Output
    data = {
        "track": "BR-C",
        "date": "2026-04-15",
        "description": "He 2^3P fine-structure multiplet benchmark via Breit-Pauli PT.",
        "method": "Closed-form LS-coupled matrix elements (Bethe-Salpeter §39) + direct sympy radial integration.",
        "reference_source": "Drake 2006 / NIST ASD",
        "reference_values_MHz": NIST_REF_MHZ,
        "reference_values_Ha": NIST_REF_HA,
        "angular_coefficients": angular,
        "result_Zeff_1.0": r1,
        "result_Zeff_1.34": r2,
        "sensitivity_scan": sensitivity_results,
        "alternative_conventions": alt_results,
        "known_issue_BR_B": (
            "The BR-B compute_rk_breit_retarded_algebraic function returns 0 for "
            "several convergent Breit-Pauli retarded integrals (e.g. (1s,1s;1s,1s) "
            "l=2 gives 0 vs published 83/640 from Bethe-Salpeter §38). The region "
            "splitting in _T_kernel_breit_retarded fails when m1<0 AND m2<0 "
            "even though the combined integral is convergent. This script "
            "bypasses BR-B by using direct sympy integration."
        ),
    }
    out_path = out_dir / "br_c_2P_benchmark.json"
    with out_path.open("w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


def _print_result(r: Dict):
    print(f"  zeta_2p = {r['zeta_2p_Ha']:.4e} Ha = {r['zeta_2p_MHz']:.2f} MHz")
    print(f"  M2_retarded  = {r['M2_retarded']:.6e}")
    print(f"  M0_retarded  = {r['M0_retarded']:.6e}")
    print(f"  A_SS  = {r['A_SS']:.4e}")
    print(f"  A_SOO = {r['A_SOO']:.4e}")
    print(f"\n  Energy contributions per J (Ha):")
    print(f"  {'J':>3} | {'E_SO':>12} | {'E_SS':>12} | {'E_SOO':>12} | {'E_total':>12}")
    for J in [0, 1, 2]:
        print(f"  {J:>3} | "
              f"{r['E_SO_J_Ha'][J]:>+12.4e} | "
              f"{r['E_SS_J_Ha'][J]:>+12.4e} | "
              f"{r['E_SOO_J_Ha'][J]:>+12.4e} | "
              f"{r['E_total_J_Ha'][J]:>+12.4e}")

    print(f"\n  Splittings (GeoVac SO+SS+SOO vs NIST):")
    for k in ["E(P0) - E(P1)", "E(P1) - E(P2)", "E(P0) - E(P2)"]:
        c = r['splittings_computed_MHz'][k]
        ref = NIST_REF_MHZ[k]
        err = r['rel_errors'][k] * 100
        sign = "OK" if (c > 0) == (ref > 0) else "FLIP"
        print(f"    {k}: {c:>+14.2f} MHz  (NIST {ref:>+10.2f}, rel err {err:+.1f}% [{sign}])")

    print(f"\n  SO-only baseline (T8 pattern):")
    for k in ["E(P0) - E(P1)", "E(P1) - E(P2)", "E(P0) - E(P2)"]:
        c = r['SO_only_baseline_MHz'][k]
        ref = NIST_REF_MHZ[k]
        err = r['SO_only_rel_errors'][k] * 100
        sign = "OK" if (c > 0) == (ref > 0) else "FLIP"
        print(f"    {k}: {c:>+14.2f} MHz  (NIST {ref:>+10.2f}, rel err {err:+.1f}% [{sign}])")


def amplitude_sensitivity_scan(Z_nuc: int = 2, Z_eff_2p: float = 1.0) -> Dict:
    """Sensitivity scan: what A_SS amplitude would exactly match NIST for a given zeta?

    This tells us whether the J-pattern (angular 9j structure) is correct,
    i.e. whether the NIST splittings CAN be reproduced with this choice of
    (f_SS, f_SOO) by adjusting A_SS and A_SOO.
    """
    from scipy.optimize import least_squares

    alpha = ALPHA_CODATA
    r3_inv_2p = (Z_eff_2p ** 3) / (2 ** 3 * 1 * 1.5 * 2)
    zeta_2p = alpha ** 2 * Z_nuc * r3_inv_2p

    # Free parameters: A_SS, A_SOO.  zeta is fixed from single-particle.
    def residuals(params):
        A_SS_fit, A_SOO_fit = params
        E_SO_J = {J: (zeta_2p / 2.0) * (J * (J + 1) - 4) for J in [0, 1, 2]}
        E_SS_J = {J: A_SS_fit * float(F_SS_BS[J]) for J in [0, 1, 2]}
        E_SOO_J = {J: A_SOO_fit * float(F_SOO_BS[J]) for J in [0, 1, 2]}
        E_tot = {J: E_SO_J[J] + E_SS_J[J] + E_SOO_J[J] for J in [0, 1, 2]}
        # Computed splittings in Ha
        dP0P1 = E_tot[0] - E_tot[1]
        dP1P2 = E_tot[1] - E_tot[2]
        # NIST in Ha
        ref_P0P1 = NIST_REF_HA["E(P0) - E(P1)"]
        ref_P1P2 = NIST_REF_HA["E(P1) - E(P2)"]
        return [dP0P1 - ref_P0P1, dP1P2 - ref_P1P2]

    # Initial guess
    sol = least_squares(residuals, [1e-5, 1e-6])
    A_SS_fit, A_SOO_fit = sol.x

    E_SO_J = {J: (zeta_2p / 2.0) * (J * (J + 1) - 4) for J in [0, 1, 2]}
    E_SS_J = {J: A_SS_fit * float(F_SS_BS[J]) for J in [0, 1, 2]}
    E_SOO_J = {J: A_SOO_fit * float(F_SOO_BS[J]) for J in [0, 1, 2]}
    E_tot = {J: E_SO_J[J] + E_SS_J[J] + E_SOO_J[J] for J in [0, 1, 2]}
    fit_split_ha = {
        "E(P0) - E(P1)": E_tot[0] - E_tot[1],
        "E(P1) - E(P2)": E_tot[1] - E_tot[2],
        "E(P0) - E(P2)": E_tot[0] - E_tot[2],
    }
    fit_split_mhz = {k: v * HA_TO_MHZ for k, v in fit_split_ha.items()}

    print(f"  Fit A_SS = {A_SS_fit:.4e} Ha, A_SOO = {A_SOO_fit:.4e} Ha")
    print(f"  [sanity: A_SS/zeta_2p = {A_SS_fit/zeta_2p:.2f}, "
          f"A_SOO/zeta_2p = {A_SOO_fit/zeta_2p:.2f}]")
    print(f"  Fit splittings (MHz):")
    for k in ["E(P0) - E(P1)", "E(P1) - E(P2)", "E(P0) - E(P2)"]:
        c = fit_split_mhz[k]
        ref = NIST_REF_MHZ[k]
        err = (c - ref) / ref * 100
        print(f"    {k}: {c:>+14.2f} MHz  (NIST {ref:>+10.2f}, rel err {err:+.3f}%)")

    # Express A_SS, A_SOO in units of alpha^2 to see what M and N they correspond to
    M_ret_implied = A_SS_fit / (1.5 * alpha ** 2)
    N_implied = 2.0 * A_SOO_fit / alpha ** 2

    # Compare with directly computed radial integrals (no fitting):
    # If the Bethe-Salpeter formula is literal, A_SS = 1.5 * alpha^2 * M2_ret.
    # With our direct sympy M2_ret = 6.93e-2, A_SS_direct = +5.54e-6 Ha.
    # The fit requires A_SS = -1.20e-6 Ha.
    # Ratio A_SS_fit / A_SS_direct = -0.22 -- the sign AND magnitude differ.
    # This tells us the convention is not literal B-S and/or there are
    # additional contributions not captured.
    alpha_num = ALPHA_CODATA
    try:
        M2_ret_direct = float(breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 2, Z_nuc))
        M0_ret_direct = float(breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 0, Z_nuc))
    except Exception:
        M2_ret_direct = 0.0
        M0_ret_direct = 0.0
    A_SS_direct_BS = 1.5 * alpha_num ** 2 * M2_ret_direct
    A_SOO_direct_BS = 0.5 * alpha_num ** 2 * (M2_ret_direct - 2 * M0_ret_direct)

    return {
        "zeta_2p_fixed_Ha": zeta_2p,
        "A_SS_fit_Ha": float(A_SS_fit),
        "A_SOO_fit_Ha": float(A_SOO_fit),
        "M_ret_implied_from_A_SS": float(M_ret_implied),
        "N_implied_from_A_SOO": float(N_implied),
        "fit_splittings_Ha": {k: float(v) for k, v in fit_split_ha.items()},
        "fit_splittings_MHz": {k: float(v) for k, v in fit_split_mhz.items()},
        "fit_residual_max": float(max(abs(sol.fun))),
        "fit_cost": float(sol.cost),
        "direct_sympy_M2_ret": M2_ret_direct,
        "direct_sympy_M0_ret": M0_ret_direct,
        "direct_BS_A_SS": A_SS_direct_BS,
        "direct_BS_A_SOO": A_SOO_direct_BS,
        "fit_vs_direct_A_SS_ratio": float(A_SS_fit / A_SS_direct_BS) if A_SS_direct_BS else None,
        "fit_vs_direct_A_SOO_ratio": float(A_SOO_fit / A_SOO_direct_BS) if A_SOO_direct_BS else None,
    }


def try_alternative_sign_conventions(Z_nuc: int = 2, Z_eff_2p: float = 1.0) -> Dict:
    """Test whether the Drake 1971 / Condon-Shortley sign convention flips
    match NIST better than our Bethe-Salpeter convention.

    Alternative convention trial:
      f_SS_alt(J)  = -f_SS_BS(J)  (sign-flipped)
      f_SOO_alt(J) = -f_SOO_BS(J)
      E_SO_alt(J)  = -E_SO_BS(J)   (sign-flipped; corresponds to using
                                    H_SO = +Z alpha^2 (kappa+1)/(4 n^3 ...)
                                    rather than the T2 minus sign)

    The NIST He 2^3P multiplet is INVERTED (E(P0) > E(P2)), and standard
    single-particle SO gives the NORMAL ordering (E(P0) < E(P2)). The
    inversion MUST come from two-body SS (or SOO), not from SO. With
    A_SS = +|something|, the SS alone pattern is:
      E_SS(J=0, f=-2) is largest negative -> E(P0) is pushed DOWN  [NORMAL]
      To invert, A_SS must be NEGATIVE (or f signs flipped).

    We test:
      (a) B-S convention with direct radial (wrong pattern)
      (b) Sign-flipped f_SS (equivalent to negative radial) (expected pattern)
      (c) All signs flipped
    """
    alpha = ALPHA_CODATA
    r3_inv_2p = (Z_eff_2p ** 3) / (2 ** 3 * 1 * 1.5 * 2)
    zeta_2p = alpha ** 2 * Z_nuc * r3_inv_2p

    try:
        M2_ret = float(breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 2, Z_nuc))
        M0_ret = float(breit_retarded_slater(1, 0, 1, 0, 2, 1, 2, 1, 0, Z_nuc))
    except Exception:
        M2_ret = 0.0
        M0_ret = 0.0

    A_SS_BS = 1.5 * alpha ** 2 * M2_ret
    A_SOO_BS = 0.5 * alpha ** 2 * (M2_ret - 2 * M0_ret)

    # Conventions to test
    conventions = [
        ("BS_direct",   +1, +1, +1, "Bethe-Salpeter §39.12-17 as cited, with positive A_SS, A_SOO"),
        ("SS_flipped",  +1, -1, +1, "Sign-flip on f_SS only"),
        ("SOO_flipped", +1, +1, -1, "Sign-flip on f_SOO only"),
        ("SS_SOO_flipped", +1, -1, -1, "Sign-flip on both f_SS and f_SOO"),
        ("all_flipped", -1, -1, -1, "Sign-flip on all three terms"),
    ]

    results = {}
    best_conv = None
    best_max_rel = float('inf')
    for name, s_so, s_ss, s_soo, desc in conventions:
        E_SO_J = {J: s_so * (zeta_2p / 2.0) * (J * (J + 1) - 4)
                  for J in [0, 1, 2]}
        E_SS_J = {J: s_ss * A_SS_BS * float(F_SS_BS[J])
                  for J in [0, 1, 2]}
        E_SOO_J = {J: s_soo * A_SOO_BS * float(F_SOO_BS[J])
                   for J in [0, 1, 2]}
        E_tot = {J: E_SO_J[J] + E_SS_J[J] + E_SOO_J[J]
                 for J in [0, 1, 2]}
        dP0P1 = (E_tot[0] - E_tot[1]) * HA_TO_MHZ
        dP1P2 = (E_tot[1] - E_tot[2]) * HA_TO_MHZ
        dP0P2 = (E_tot[0] - E_tot[2]) * HA_TO_MHZ
        err_P0P1 = (dP0P1 - NIST_REF_MHZ["E(P0) - E(P1)"]) / NIST_REF_MHZ["E(P0) - E(P1)"]
        err_P1P2 = (dP1P2 - NIST_REF_MHZ["E(P1) - E(P2)"]) / NIST_REF_MHZ["E(P1) - E(P2)"]
        err_P0P2 = (dP0P2 - NIST_REF_MHZ["E(P0) - E(P2)"]) / NIST_REF_MHZ["E(P0) - E(P2)"]
        max_rel = max(abs(err_P0P1), abs(err_P1P2), abs(err_P0P2))
        results[name] = {
            "desc": desc,
            "signs": (s_so, s_ss, s_soo),
            "dP0P1_MHz": dP0P1,
            "dP1P2_MHz": dP1P2,
            "dP0P2_MHz": dP0P2,
            "err_P0P1_rel": err_P0P1,
            "err_P1P2_rel": err_P1P2,
            "err_P0P2_rel": err_P0P2,
            "max_rel_err": max_rel,
        }
        print(f"  {name:>18}: dP0P1={dP0P1:+10.1f} MHz, dP1P2={dP1P2:+10.1f} MHz, "
              f"dP0P2={dP0P2:+10.1f} MHz, max_rel_err={max_rel*100:+.1f}%")
        if max_rel < best_max_rel:
            best_max_rel = max_rel
            best_conv = name

    results["best_convention"] = best_conv
    results["best_max_rel_err_pct"] = best_max_rel * 100
    print(f"\n  Best convention: {best_conv} with max rel err {best_max_rel*100:.1f}%")
    return results


if __name__ == "__main__":
    main()
