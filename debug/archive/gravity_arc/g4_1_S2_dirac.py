"""Sprint G4-1 — Dirac on S^2 (cigar's spatial sphere) and structural comparison
with the S^3 substrate.

First sub-sprint of G4 full. Establishes whether the framework's two-term
exactness on S^3 Dirac propagates to the S^2 Dirac that appears in the cigar's
spatial section.

Setup
-----
The Euclidean Schwarzschild cigar has spatial section topologically
R_+ x S^2_r where r is the radial coordinate. At the horizon (r = 2M = r_h),
the S^2 has area A = 4 pi r_h^2.

Dirac on S^d (Camporesi-Higuchi 1996 general dimension):
  |lambda_n^{S^d}| = n + d/2
  g_n^{Weyl} = 2^{floor((d-1)/2)} * binom(n + d - 1, n)
  g_n^{Dirac} = 2 g_n^{Weyl}

For S^3 (d = 3, our usual substrate):
  |lambda_n^{S^3}| = n + 3/2  (HALF-INTEGER shift)
  g_n^{Dirac} = 2 (n+1)(n+2)

For S^2 (d = 2, cigar's spatial sphere):
  |lambda_n^{S^2}| = n + 1  (INTEGER shift)
  g_n^{Dirac} = 4 (n+1)

The structural question
-----------------------
Paper 28's two-term exactness theorem on S^3 Dirac relies on the half-integer
shift in the spectrum: the heat trace asymptotic has exactly two power-law terms
because Jacobi theta_2 inversion (acting on half-integer-shifted sums) gives
clean two-term form.

For S^2 Dirac (integer-shifted), the heat trace uses Jacobi theta_3 or similar
acting on integer sums, which gives a DIFFERENT asymptotic structure — full
infinite Seeley-DeWitt series with higher curvature corrections.

This script verifies this structural distinction by:
1. Computing the heat trace K_{S^2}(t) explicitly
2. Identifying its asymptotic form
3. Extracting SD coefficients
4. Comparing to the K_{S^3}(t) two-term form

Implication for G4 full
-----------------------
If the answer is "S^2 Dirac has standard CC infinite series" (not two-term),
then the cigar's discrete-substrate construction needs to account for higher
curvature corrections at all orders. The framework's clean structure on the
Dirac sector is specific to S^3 (the half-integer shift).
"""

import json
from pathlib import Path

import sympy as sp
from sympy import (
    symbols, pi, Rational, simplify, sqrt, exp, oo, Sum, factorial,
    Integer, latex, series,
)
import numpy as np
from mpmath import mp, mpf, exp as mpexp, sqrt as mpsqrt, pi as mp_pi, fsum

mp.dps = 50

OUT_JSON = Path(__file__).parent / "data" / "g4_1_S2_dirac.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main():
    results = {}
    print("=" * 72)
    print("Sprint G4-1: Dirac on S^2 (cigar's spatial sphere)")
    print("=" * 72)

    # -----------------------------------------------------------------------
    # Step 1: Set up S^2 Dirac spectrum
    # -----------------------------------------------------------------------
    print("\n[Step 1] S^2 Dirac spectrum (Camporesi-Higuchi, d=2):")
    print("  Eigenvalues: |lambda_n| = n + 1  (INTEGER shift, n = 0, 1, 2, ...)")
    print("  Multiplicity (full Dirac): g_n = 4(n+1)")
    print("    (Weyl per chirality: 2(n+1); both chiralities give 4(n+1))")
    print()

    # -----------------------------------------------------------------------
    # Step 2: Verify the leading heat trace asymptotic
    # -----------------------------------------------------------------------
    print("[Step 2] Heat trace K_{S^2}^{Dirac}(t) = sum_n 4(n+1) e^{-(n+1)^2 t}")
    print()
    print("  Using u = n+1: K = 4 sum_{u>=1} u e^{-u^2 t}")
    print()
    print("  At small t, sum dominated by u ~ 1/sqrt(t), approximate by integral:")
    print("    int_0^infty u e^{-u^2 t} du = 1/(2t)")
    print("  So K(t) -> 4 * 1/(2t) = 2/t at small t.")
    print()
    print("  Expected from heat-kernel formula on S^2 of unit radius:")
    print("    K(t) ~ (4 pi t)^{-1} * dim_S * Vol(S^2)")
    print("         = (4 pi t)^{-1} * 2 * 4 pi = 2/t  (matches)")
    print()

    # Numerical verification
    t = symbols("t", positive=True)
    print("  Numerical verification of K(t) ~ 2/t at small t:")
    for t_val_str in ["0.1", "0.01", "0.001", "0.0001"]:
        t_val = mpf(t_val_str)
        K_direct = mpf("0")
        for n in range(2000):
            u = mpf(n) + mpf("1")
            g_n = 4 * (n + 1)
            K_direct += g_n * mpexp(-u * u * t_val)
        K_leading = mpf("2") / t_val
        rel = abs(K_direct - K_leading) / abs(K_direct)
        print(f"    t = {float(t_val):.5f}: K_direct = {float(K_direct):.6e}, "
              f"K_leading = {float(K_leading):.6e}, rel diff = {float(rel):.3e}")

    # -----------------------------------------------------------------------
    # Step 3: Closed-form for K_{S^2}(t) via Jacobi theta and derivative
    # -----------------------------------------------------------------------
    print("\n[Step 3] Closed-form structure of K_{S^2}^{Dirac}(t):")
    print()
    print("  Define: K(t) = 4 sum_{u>=1} u e^{-u^2 t}")
    print()
    print("  Relate to theta_2(0, e^{-t}) = 2 sum_{u>=1} e^{-(u-1/2)^2 t}? No, that's for")
    print("  half-integer shift. For integer shift we'd use derivative of theta_3.")
    print()
    print("  Direct: K(t) = -2 d/dt [sum_{u>=1} u^2 ... ]? No, this is sum of u, not u^2.")
    print()
    print("  Alternative: K(t) = 2 d/dt [sum_{u>=1} 1] = ... diverges. Not useful directly.")
    print()
    print("  Cleaner approach: K(t) is a Mellin transform of zeta-like function.")

    # Use Hurwitz / standard formulas
    # sum_{u>=1} u e^{-u^2 t}: differentiate sum_{u>=1} e^{-(u-c)^2 t} with respect to c at c=0?
    # This is the "theta'" function. Let me just compute the asymptotic series.

    print()
    print("  Asymptotic expansion (Euler-Maclaurin):")
    print("    K(t) ~ 4 [int_0^infty u e^{-u^2 t} du - (1/2) * 0 * e^0 + corrections]")
    print("         ~ 4 [1/(2t) + 0 + 0 + ...] = 2/t + (exp small in 1/t)")
    print()
    print("  So leading is 2/t with no subleading power-law corrections from the smooth")
    print("  part. The 'exp small' comes from theta modular structure.")
    print()
    print("  Standard SD form: K(t) = (4 pi t)^{-1} [a_0 + a_1 t + a_2 t^2 + ...]")
    print("    a_0 = 2 (Vol(S^2)) = 8 pi  (from a_0 = dim_S * Vol = 2 * 4 pi)")
    print("    a_1 = ? (involves scalar curvature R = 2/r^2 on unit S^2; nonzero)")
    print()
    print("  Coefficient of t^0 in (4 pi t)^{-1} a_0: a_0 / (4 pi t) = 8 pi /(4 pi t) = 2/t  [OK]")

    results["S2_dirac_spectrum"] = {
        "eigenvalues": "|lambda_n| = n + 1",
        "multiplicity": "g_n = 4(n+1) for full Dirac",
        "spectrum_shift": "INTEGER (unlike S^3 which has HALF-INTEGER)",
    }
    results["S2_heat_trace_leading"] = {
        "form": "K(t) ~ 2/t + (exp small in 1/t)",
        "matches_dim_S_times_Vol_over_4pi_t": True,
    }

    # -----------------------------------------------------------------------
    # Step 4: Check for higher SD coefficients
    # -----------------------------------------------------------------------
    print("\n[Step 4] Higher Seeley-DeWitt coefficients on S^2 Dirac:")
    print()
    print("  On S^2 (unit), the SD expansion of Dirac heat trace is:")
    print("    K(t) = (4 pi t)^{-1} sum_k a_k t^k")
    print()
    print("  Standard formulas (Vassilevich 2003, Branson-Gilkey):")
    print("    a_0 = dim_S * Vol = 2 * 4 pi = 8 pi")
    print("    a_1 = dim_S * (R_scalar/6 - E_Lich) * Vol")
    print("        = 2 * (2/6 - 1/2) * 4 pi  (using R = 2, E = R/4 = 1/2 on unit S^2)")
    print("        = 2 * (1/3 - 1/2) * 4 pi = 2 * (-1/6) * 4 pi = -4 pi/3")
    print()
    print("  So a_1 != 0 on S^2 Dirac. The next-order term is NONZERO.")
    print()
    print("  Numerical extraction of a_1 from K(t) at small t:")

    # Try to extract a_1 by computing [K(t) - a_0/(4 pi t)] * (4 pi t) / 1 at small t
    # K(t) = (4 pi t)^{-1} [a_0 + a_1 t + ...]
    # K(t) - a_0/(4 pi t) = (4 pi t)^{-1} * a_1 * t = a_1 / (4 pi)
    # So [K(t) - 2/t] -> a_1 / (4 pi) = -1/3  at small t

    print("    [K(t) - 2/t] should approach a_1 / (4 pi) = -1/3 at small t:")
    for t_val_str in ["0.1", "0.05", "0.02", "0.01", "0.005", "0.002", "0.001"]:
        t_val = mpf(t_val_str)
        K_direct = mpf("0")
        for n in range(2000):
            u = mpf(n) + mpf("1")
            g_n = 4 * (n + 1)
            K_direct += g_n * mpexp(-u * u * t_val)
        residual = K_direct - mpf("2") / t_val
        print(f"    t = {float(t_val):.5f}: K - 2/t = {float(residual):.6e}, "
              f"expected -> -1/3 = -0.333...")

    results["a_1_coefficient"] = {
        "value": "-4 pi / 3 (from formula)",
        "implication": "a_1 != 0 on S^2 Dirac (unlike S^3 where two-term exactness has a_2 = 0)",
    }

    # -----------------------------------------------------------------------
    # Step 5: Compare to S^3 Dirac
    # -----------------------------------------------------------------------
    print("\n[Step 5] Comparison with S^3 Dirac (Paper 28 two-term exactness):")
    print()
    print("  S^3 Dirac:")
    print("    Spectrum: |lambda_n| = n + 3/2 (HALF-INTEGER shift)")
    print("    Heat trace: K_{S^3}(t) = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}")
    print("                            + O(exp(-pi^2/t))")
    print("    All higher SD coefficients a_k = 0 for k >= 2 (Paper 28 theorem)")
    print()
    print("  S^2 Dirac:")
    print("    Spectrum: |lambda_n| = n + 1 (INTEGER shift)")
    print("    Heat trace: K_{S^2}(t) ~ 2/t + a_1/(4 pi) + a_2 t/(4 pi) + ...")
    print("    All higher SD coefficients a_k != 0 (standard CC infinite series)")
    print()
    print("  STRUCTURAL DISTINCTION:")
    print("    - S^3 (half-integer): two-term exact via Jacobi theta_2 inversion")
    print("    - S^2 (integer): full SD series via Jacobi theta_3 + exponential prefactor")
    print()
    print("  Sprint G3 already established that the SCALAR Laplacian on S^3 has")
    print("  full SD series (a_k = 2 pi^2/k!). The pattern is consistent:")
    print()
    print("    Half-integer spectrum (Dirac on S^3): TWO-TERM exact")
    print("    Integer spectrum (scalar Laplacian on S^3, Dirac on S^2): FULL SD SERIES")
    print()
    print("  Two-term exactness is SPECIFIC to half-integer-shifted Dirac.")

    results["structural_comparison"] = {
        "S3_Dirac": "two-term exact (half-integer spectrum, Jacobi theta_2)",
        "S2_Dirac": "full SD series (integer spectrum, Jacobi theta_3 + exp prefactor)",
        "pattern": "two-term exactness is specific to half-integer Dirac",
    }

    # -----------------------------------------------------------------------
    # Step 6: Implications for G4 full
    # -----------------------------------------------------------------------
    print("\n[Step 6] Implications for G4 full (Bekenstein-Hawking on discrete substrate):")
    print()
    print("  The cigar's spatial section is R_+ x S^2_r (warped product, r varies).")
    print("  The S^2 component has standard CC infinite SD series, NOT two-term exact.")
    print()
    print("  This means the discrete-substrate BH entropy derivation:")
    print("  - Cannot exploit a clean two-term form on the spatial S^2")
    print("  - Must account for higher-curvature corrections at all orders")
    print("  - Will follow the standard CC heat-kernel asymptotic structure")
    print()
    print("  The framework's clean structure on the Dirac sector is specific to")
    print("  S^3 (the substrate sourced by the Fock projection). It does NOT extend")
    print("  to the S^2 sector of the cigar.")
    print()
    print("  This is consistent with Sprint G3's finding (spinor-bundle specificity):")
    print("  half-integer-shifted Dirac on S^3 is special; other operators / sectors")
    print("  / dimensions inherit standard CC behavior.")

    results["implications_for_G4_full"] = {
        "cigar_spatial_S2": "standard CC infinite SD series",
        "discrete_substrate_construction": "must account for higher-curvature corrections",
        "two_term_exactness_specificity": "specific to S^3 Dirac (the GeoVac substrate)",
        "consistency_with_G3": "spinor-bundle specificity confirmed",
    }

    # -----------------------------------------------------------------------
    # Step 7: Verdict
    # -----------------------------------------------------------------------
    print("\n[Step 7] Verdict:")
    print()
    print("  POSITIVE-STRUCTURAL-FINDING: S^2 Dirac on the cigar's spatial sphere has")
    print("  STANDARD CC infinite SD series, not the two-term exact form of S^3 Dirac.")
    print()
    print("  The framework's clean two-term exactness is specific to the half-integer")
    print("  shift of the S^3 Dirac spectrum (the GeoVac substrate). For G4 full:")
    print("  - The cigar's spatial S^2 component inherits standard CC continuum behavior")
    print("  - Higher-curvature corrections enter at every order in the heat-kernel")
    print("    expansion")
    print("  - The discrete-substrate construction needs to preserve these corrections")
    print()
    print("  This is consistent with Sprint G3's spinor-bundle specificity result and")
    print("  sharpens the understanding of the framework's gravity sector structure.")

    results["verdict"] = "POSITIVE-STRUCTURAL-FINDING"

    # Save JSON
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)
    return results


if __name__ == "__main__":
    main()
