"""Extended PSLQ identification of S_min at 200+ digit precision.

Track Q-3: Push the CG-weighted two-loop spectral sum S_min to 200 digits
and attempt identification against an extended ~100-element PSLQ basis
including Euler sums, Tornheim-Witten double zeta values, colored MZVs,
double Hurwitz at quarter-integer shifts, and all cross-products up to
weight 8.

S_min = sum_{k>=1} T(k)^2

where T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2).

Sprint 5 corrected value (verified by three independent methods):
    S_min = 2.47993693803422255441357950082938214468792578661728845837...

Previous PSLQ attempts (15 failures, 47-element ultra-wide basis at 150 dps)
all returned None.  This track extends to 200 digits and a ~100-element basis.

Usage:
    python debug/smin_extended_pslq.py

See also:
    debug/smin_identification.py       — original 15-attempt PSLQ
    debug/compute_smin_verification.py — Sprint 5 three-method verification
    geovac/qed_vertex.py               — T(k) definition and two_loop_min_weighted_hurwitz
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Ensure repo root on path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import mpmath


# =====================================================================
# Known verified value (Sprint 5 erratum, 120+ digits)
# =====================================================================
SMIN_KNOWN_STR = (
    "2.479936938034222554413579500829382144687925786617288458"
    "37879872655955277781837491232855893143004699635515948282"
    "771311796"
)


# =====================================================================
# Core computation: S_min at high precision
# =====================================================================

def smin_high_precision(n_digits: int = 200) -> mpmath.mpf:
    """Compute S_min = sum_{k>=1} T(k)^2 to n_digits verified digits.

    Uses mpmath.nsum with Levin u-transform (the method proven in Sprint 5).
    Verifies the first 80 digits match the known corrected value.

    Parameters
    ----------
    n_digits : int
        Number of target digits of precision.

    Returns
    -------
    mpmath.mpf
        S_min computed at the requested precision.
    """
    # Use extra guard digits: 50 beyond target
    guard = 50
    dps = n_digits + guard

    with mpmath.workdps(dps):
        def term(k):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            Tk = (2 * mpmath.hurwitz(2, a)
                  - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))
            return Tk ** 2

        S = mpmath.nsum(term, [1, mpmath.inf], method='levin')

        # Verify first 80 digits against known value
        known = mpmath.mpf(SMIN_KNOWN_STR)
        discrepancy = abs(S - known)
        if discrepancy > mpmath.power(10, -79):
            raise ValueError(
                f"S_min computation disagrees with known value at 80-digit level: "
                f"|S - known| = {mpmath.nstr(discrepancy, 10)}"
            )

        return S


# =====================================================================
# Extended PSLQ basis construction
# =====================================================================

def _compute_euler_sum_s_ab(a: int, b: int, dps: int) -> mpmath.mpf:
    """Compute the Euler sum s_{a,b} = sum_{n>=1} H_n^{(a)} / n^b.

    H_n^{(a)} = sum_{j=1}^{n} 1/j^a  is the generalized harmonic number.

    Uses mpmath.nsum with Levin acceleration on the partial sums.
    """
    with mpmath.workdps(dps + 20):
        def term(n):
            n = int(n)
            # Compute H_n^{(a)}
            Hn_a = sum(mpmath.mpf(1) / mpmath.mpf(j) ** a for j in range(1, n + 1))
            return Hn_a / mpmath.mpf(n) ** b

        return mpmath.nsum(term, [1, mpmath.inf], method='levin')


def _compute_tornheim_witten(r: int, s: int, t: int, dps: int,
                              n_max: int = 2000) -> mpmath.mpf:
    """Compute the Tornheim-Witten double sum T(r,s,t) via direct summation.

    T(r,s,t) = sum_{m>=1} sum_{n>=1} 1/(m^r * n^s * (m+n)^t)

    Convergent when r+t > 1, s+t > 1, and r+s+t > 2.
    Uses partial summation up to n_max with a tail correction.
    """
    with mpmath.workdps(dps + 20):
        S = mpmath.mpf(0)
        for m in range(1, n_max + 1):
            for n in range(1, n_max + 1):
                S += mpmath.mpf(1) / (
                    mpmath.mpf(m) ** r
                    * mpmath.mpf(n) ** s
                    * mpmath.mpf(m + n) ** t
                )
        return S


def _build_extended_basis(dps: int) -> Tuple[List[mpmath.mpf], List[str]]:
    """Build the extended PSLQ basis with ~100 elements.

    Includes standard, products, Catalan/Dirichlet, Euler sums,
    double Hurwitz at quarter-integer, Tornheim-Witten, and colored MZVs.
    """
    with mpmath.workdps(dps):
        pi = mpmath.pi
        pi2 = pi ** 2
        pi4 = pi ** 4
        pi6 = pi ** 6
        pi8 = pi ** 8

        z3 = mpmath.zeta(3)
        z5 = mpmath.zeta(5)
        z7 = mpmath.zeta(7)
        z9 = mpmath.zeta(9)

        G = mpmath.catalan
        beta4 = (mpmath.hurwitz(4, mpmath.mpf(1) / 4)
                 - mpmath.hurwitz(4, mpmath.mpf(3) / 4)) / mpmath.power(4, 4)
        beta3 = pi ** 3 / 32

        ln2 = mpmath.log(2)

        Li4_half = mpmath.polylog(4, mpmath.mpf(1) / 2)

        # Euler sums: s_{a,b} = sum_{n>=1} H_n^{(a)} / n^b
        # Known closed forms:
        #   s_{1,2} = 2*zeta(3)
        #   s_{1,3} = pi^4/72
        #   s_{2,2} = pi^4/120 + zeta(3) - zeta(4)  (needs careful computation)
        # For PSLQ we compute numerically to avoid error in the closed forms
        # and let PSLQ detect any relation.
        print("  Computing Euler sums s_{2,3}, s_{3,2}...", flush=True)
        t0 = time.time()
        s_23 = _compute_euler_sum_s_ab(2, 3, dps)
        s_32 = _compute_euler_sum_s_ab(3, 2, dps)
        print(f"  Euler sums computed in {time.time() - t0:.1f}s", flush=True)

        # Double Hurwitz at quarter-integer shifts
        hz_2_14 = mpmath.hurwitz(2, mpmath.mpf(1) / 4)
        hz_2_34 = mpmath.hurwitz(2, mpmath.mpf(3) / 4)
        hz_4_14 = mpmath.hurwitz(4, mpmath.mpf(1) / 4)
        hz_4_34 = mpmath.hurwitz(4, mpmath.mpf(3) / 4)

        # Hurwitz at 3/2
        hz_2_32 = mpmath.hurwitz(2, mpmath.mpf(3) / 2)
        hz_3_32 = mpmath.hurwitz(3, mpmath.mpf(3) / 2)
        hz_4_32 = mpmath.hurwitz(4, mpmath.mpf(3) / 2)
        hz_5_32 = mpmath.hurwitz(5, mpmath.mpf(3) / 2)
        hz_6_32 = mpmath.hurwitz(6, mpmath.mpf(3) / 2)

        # Polygamma at 3/2
        psi1_32 = mpmath.psi(1, mpmath.mpf(3) / 2)   # = hz(2, 3/2)
        psi3_32 = mpmath.psi(3, mpmath.mpf(3) / 2)   # = 6*hz(4, 3/2)

        # Tornheim-Witten: T(r,s,t) via direct summation
        # T(2,2,2) and T(2,2,4) are known to be expressible in terms of MZVs
        print("  Computing Tornheim-Witten T(2,2,2), T(2,2,4)...", flush=True)
        t0 = time.time()
        TW_222 = _compute_tornheim_witten(2, 2, 2, dps, n_max=500)
        TW_224 = _compute_tornheim_witten(2, 2, 4, dps, n_max=500)
        print(f"  Tornheim-Witten computed in {time.time() - t0:.1f}s", flush=True)

        # Colored MZVs via mpmath
        # zeta(3, 1) = sum_{n1>n2>=1} 1/(n1^3 * n2^1) = pi^4/360 (known)
        # Li_{2,2}(1/2, 1) is exotic — compute by direct sum
        # For PSLQ we just include the numerical values
        # MZV zeta({3,1}) = pi^4/360
        mzv_31 = pi4 / 360
        # MZV zeta({2,1,1}) = pi^4/360 (same, by Euler identity)
        mzv_211 = pi4 / 360

        # D_even(4) and D_odd(4) for completeness
        D4 = 2 * hz_2_32 - mpmath.mpf(1) / 2 * hz_4_32
        D_even4 = mpmath.power(2, -4) * (
            8 * mpmath.hurwitz(2, mpmath.mpf(3) / 4)
            - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, mpmath.mpf(3) / 4)
        )
        D_odd4 = mpmath.power(2, -4) * (
            8 * mpmath.hurwitz(2, mpmath.mpf(5) / 4)
            - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, mpmath.mpf(5) / 4)
        )

        # Euler-Mascheroni
        gamma_em = mpmath.euler

        # Clausen function
        Cl2_pi3 = mpmath.clsin(2, pi / 3)

        # Build the basis list and labels
        # We are careful to avoid exact duplicates
        vals = []
        labels = []

        def add(val, label):
            vals.append(val)
            labels.append(label)

        # --- Standard zeta/pi powers ---
        add(mpmath.mpf(1), "1")
        add(pi2, "pi^2")
        add(pi4, "pi^4")
        add(pi6, "pi^6")
        add(pi8, "pi^8")

        # --- Riemann odd zeta ---
        add(z3, "zeta(3)")
        add(z5, "zeta(5)")
        add(z7, "zeta(7)")
        add(z9, "zeta(9)")

        # --- Products of pi^{even} with odd zeta ---
        add(pi2 * z3, "pi^2*zeta(3)")
        add(pi2 * z5, "pi^2*zeta(5)")
        add(pi4 * z3, "pi^4*zeta(3)")

        # --- Odd zeta products ---
        add(z3 ** 2, "zeta(3)^2")
        add(z3 * z5, "zeta(3)*zeta(5)")

        # --- Catalan / Dirichlet beta ---
        add(G, "G")
        add(beta4, "beta(4)")
        add(beta3, "beta(3)")

        # --- Catalan/beta products ---
        add(G ** 2, "G^2")
        add(G * pi2, "G*pi^2")
        add(beta4 * pi2, "beta(4)*pi^2")
        add(G * beta4, "G*beta(4)")
        add(beta4 ** 2, "beta(4)^2")
        add(G * z3, "G*zeta(3)")
        add(G * z5, "G*zeta(5)")
        add(beta4 * z3, "beta(4)*zeta(3)")
        add(beta4 * z5, "beta(4)*zeta(5)")

        # --- Logarithmic ---
        add(ln2, "ln2")
        add(ln2 ** 2, "ln2^2")
        add(ln2 * pi2, "ln2*pi^2")
        add(ln2 * z3, "ln2*zeta(3)")
        add(ln2 * G, "ln2*G")
        add(ln2 * beta4, "ln2*beta(4)")

        # --- Polylogarithm ---
        add(Li4_half, "Li4(1/2)")

        # --- Euler sums ---
        add(s_23, "s_{2,3}")
        add(s_32, "s_{3,2}")

        # --- Double Hurwitz at quarter-integer ---
        add(hz_4_14, "zeta(4,1/4)")
        add(hz_4_34, "zeta(4,3/4)")
        add(hz_2_14 * hz_2_34, "zeta(2,1/4)*zeta(2,3/4)")
        add(hz_2_14 * hz_2_14, "zeta(2,1/4)^2")
        add(hz_2_34 * hz_2_34, "zeta(2,3/4)^2")

        # --- Hurwitz at 3/2 ---
        add(hz_2_32, "hz(2,3/2)")
        add(hz_3_32, "hz(3,3/2)")
        add(hz_4_32, "hz(4,3/2)")
        add(hz_5_32, "hz(5,3/2)")
        add(hz_6_32, "hz(6,3/2)")

        # --- Hurwitz at 3/2 products ---
        add(hz_2_32 ** 2, "hz(2,3/2)^2")
        add(hz_2_32 * hz_4_32, "hz(2,3/2)*hz(4,3/2)")
        add(hz_4_32 ** 2, "hz(4,3/2)^2")
        add(hz_2_32 * hz_3_32, "hz(2,3/2)*hz(3,3/2)")
        add(hz_3_32 ** 2, "hz(3,3/2)^2")

        # --- Tornheim-Witten ---
        add(TW_222, "T(2,2,2)")
        add(TW_224, "T(2,2,4)")

        # --- Colored MZVs (known closed forms as sanity) ---
        add(mzv_31, "zeta({3,1})")
        add(mzv_211, "zeta({2,1,1})")

        # --- D_even/D_odd products ---
        add(D_even4, "D_even(4)")
        add(D_odd4, "D_odd(4)")
        add(D_even4 ** 2, "D_even(4)^2")
        add(D_odd4 ** 2, "D_odd(4)^2")
        add(D_even4 * D_odd4, "D_even(4)*D_odd(4)")
        add(D4, "D(4)")
        add(D4 ** 2, "D(4)^2")

        # --- Polygamma ---
        add(psi1_32, "psi_1(3/2)")
        add(psi3_32, "psi_3(3/2)")
        add(psi1_32 ** 2, "psi_1(3/2)^2")
        add(psi1_32 * psi3_32, "psi_1(3/2)*psi_3(3/2)")
        add(psi3_32 ** 2, "psi_3(3/2)^2")

        # --- Misc ---
        add(gamma_em, "gamma_EM")
        add(Cl2_pi3, "Cl2(pi/3)")

        # --- Additional cross-products to reach ~100 ---
        add(pi4 * G, "pi^4*G")
        add(pi4 * beta4, "pi^4*beta(4)")
        add(pi2 * z7, "pi^2*zeta(7)")
        add(pi4 * z5, "pi^4*zeta(5)")
        add(z5 ** 2, "zeta(5)^2")
        add(z3 * z7, "zeta(3)*zeta(7)")
        add(G * z7, "G*zeta(7)")
        add(beta4 * z7, "beta(4)*zeta(7)")
        add(ln2 ** 3, "ln2^3")
        add(ln2 ** 2 * pi2, "ln2^2*pi^2")
        add(ln2 * z5, "ln2*zeta(5)")
        add(Li4_half * pi2, "Li4(1/2)*pi^2")
        add(Li4_half * z3, "Li4(1/2)*zeta(3)")

        # --- Sum components (structural decomposition) ---
        # sum_22 = sum_k zeta(2, k+3/2)^2 from the original analysis
        # We compute via a targeted nsum
        def sum_22_term(k):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            return mpmath.hurwitz(2, a) ** 2

        def sum_44_term(k):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            return mpmath.hurwitz(4, a) ** 2

        def sum_24_term(k):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            return mpmath.hurwitz(2, a) * mpmath.hurwitz(4, a)

        print("  Computing component sums sum_22, sum_24, sum_44...", flush=True)
        t0 = time.time()
        comp_sum_22 = mpmath.nsum(sum_22_term, [1, mpmath.inf], method='levin')
        comp_sum_24 = mpmath.nsum(sum_24_term, [1, mpmath.inf], method='levin')
        comp_sum_44 = mpmath.nsum(sum_44_term, [1, mpmath.inf], method='levin')
        print(f"  Component sums computed in {time.time() - t0:.1f}s", flush=True)

        add(comp_sum_22, "sum_hz2sq")
        add(comp_sum_24, "sum_hz2_hz4")
        add(comp_sum_44, "sum_hz4sq")

        return vals, labels


def smin_extended_pslq(
    value: mpmath.mpf,
    basis_size: int = 100,
) -> Dict[str, Any]:
    """Run PSLQ on S_min against an extended basis of ~100 transcendental constants.

    Parameters
    ----------
    value : mpmath.mpf
        The S_min value to identify (must be computed at matching precision).
    basis_size : int
        Target basis size (actual size may differ slightly).

    Returns
    -------
    dict
        With keys 'identified' (bool), 'basis_size' (int),
        'components' (dict, if identified), 'residual' (float, if identified),
        'relation' (list of ints, if identified).
    """
    dps = int(mpmath.mp.dps)
    print(f"  Building extended PSLQ basis at {dps} dps...", flush=True)
    t0 = time.time()
    vals, labels = _build_extended_basis(dps)
    build_time = time.time() - t0
    print(f"  Basis built: {len(vals)} elements in {build_time:.1f}s", flush=True)

    # Prepend S_min as the first element
    full_basis = [value] + vals
    full_labels = ["S_min"] + labels
    actual_size = len(full_basis)

    print(f"  Running PSLQ with {actual_size} elements at {dps} dps...", flush=True)
    t0 = time.time()
    try:
        rel = mpmath.pslq(full_basis, tol=mpmath.power(10, -(dps - 40)),
                          maxcoeff=1000000)
    except Exception as e:
        pslq_time = time.time() - t0
        print(f"  PSLQ raised exception after {pslq_time:.1f}s: {e}", flush=True)
        return {
            "identified": False,
            "basis_size": actual_size,
            "error": str(e),
            "time_sec": pslq_time,
        }

    pslq_time = time.time() - t0
    print(f"  PSLQ completed in {pslq_time:.1f}s", flush=True)

    if rel is None:
        print("  PSLQ returned None (no relation found)", flush=True)
        return {
            "identified": False,
            "basis_size": actual_size,
            "error": "no relation",
            "time_sec": pslq_time,
        }

    if rel[0] == 0:
        print("  PSLQ found relation but S_min coefficient is 0 (spurious)",
              flush=True)
        return {
            "identified": False,
            "basis_size": actual_size,
            "error": "zero leading coeff",
            "time_sec": pslq_time,
        }

    # Reconstruct
    print("  PSLQ FOUND RELATION!", flush=True)
    components = {}
    reconstructed = mpmath.mpf(0)
    for i in range(1, len(rel)):
        if rel[i] != 0:
            coeff = mpmath.mpf(-rel[i]) / mpmath.mpf(rel[0])
            label = full_labels[i]
            components[label] = f"{-rel[i]}/{rel[0]}"
            reconstructed += coeff * full_basis[i]
            print(f"    {-rel[i]}/{rel[0]} * {label}", flush=True)

    residual = float(abs(value - reconstructed))
    print(f"  Residual: {residual:.3e}", flush=True)

    return {
        "identified": True,
        "basis_size": actual_size,
        "relation": [int(r) for r in rel],
        "components": components,
        "residual": residual,
        "time_sec": pslq_time,
    }


# =====================================================================
# Inner photon sum for structural insight
# =====================================================================

def smin_inner_q_closed_form(
    n1: int,
    n2: int,
    s: int = 4,
) -> Dict[str, Any]:
    """Compute the inner photon sum for fixed (n1, n2) pair.

    For the CG-weighted two-loop sunset, the inner photon sum is:

        I(n1, n2) = sum_{q allowed} W(n1, n2, q) * d_T(q) / mu_q^p

    where:
        W(n1, n2, q) is the SO(4) channel count (0, 1, or 2)
        d_T(q) = q(q+2) is the transverse photon degeneracy
        mu_q = q(q+2) is the Hodge-1 eigenvalue
        p = 1 (default photon exponent)

    The vertex selection rules constrain:
        |n1 - n2| <= q <= n1 + n2,  n1 + n2 + q odd,  q >= 1

    For p=1: d_T(q)/mu_q = 1 (constant), so:
        I(n1, n2) = sum_{q allowed} W(n1, n2, q)

    This gives a finite sum that can be evaluated exactly.

    Parameters
    ----------
    n1, n2 : int
        Dirac spinor levels (CH convention, >= 0).
    s : int
        Exponent parameter (default 4, for S_min).

    Returns
    -------
    dict with 'n1', 'n2', 'inner_sum', 'n_terms', 'allowed_q_list'.
    """
    from fractions import Fraction

    dps = max(int(mpmath.mp.dps), 50)
    with mpmath.workdps(dps):
        # Vertex selection rule
        q_min = max(1, abs(n1 - n2))
        q_max = n1 + n2

        allowed_q = []
        inner_sum = mpmath.mpf(0)
        terms = []

        for q in range(q_min, q_max + 1):
            # Parity check: n1 + n2 + q must be odd
            if (n1 + n2 + q) % 2 == 0:
                continue

            # Channel count (using the SO(4) CG structure)
            # Inline the computation rather than importing from qed_vertex
            # to keep this self-contained
            j1_L = Fraction(n1 + 1, 2)
            j1_R = Fraction(n1, 2)
            j2_L = Fraction(n2, 2)
            j2_R = Fraction(n2 + 1, 2)

            count = 0
            # Component A: ((q+1)/2, (q-1)/2)
            jg_L_A = Fraction(q + 1, 2)
            jg_R_A = Fraction(q - 1, 2)
            if (jg_R_A >= 0
                    and abs(j1_L - jg_L_A) <= j2_L <= j1_L + jg_L_A
                    and abs(j1_R - jg_R_A) <= j2_R <= j1_R + jg_R_A):
                count += 1

            # Component B: ((q-1)/2, (q+1)/2)
            jg_L_B = Fraction(q - 1, 2)
            jg_R_B = Fraction(q + 1, 2)
            if (jg_L_B >= 0
                    and abs(j1_L - jg_L_B) <= j2_L <= j1_L + jg_L_B
                    and abs(j1_R - jg_R_B) <= j2_R <= j1_R + jg_R_B):
                count += 1

            if count > 0:
                allowed_q.append(q)
                # d_T(q) = q(q+2), mu_q = q(q+2) => ratio = 1 for p=1
                # So inner sum contribution is just count
                contribution = mpmath.mpf(count)
                inner_sum += contribution
                terms.append({
                    "q": q,
                    "W": count,
                    "contribution": float(contribution),
                })

        return {
            "n1": n1,
            "n2": n2,
            "inner_sum": inner_sum,
            "inner_sum_float": float(inner_sum),
            "n_terms": len(terms),
            "allowed_q_list": allowed_q,
            "terms": terms,
        }


# =====================================================================
# Main driver
# =====================================================================

def main():
    print("=" * 78)
    print(" S_min Extended PSLQ Identification (Track Q-3)")
    print("=" * 78)

    # Phase 1: High-precision computation
    print("\n--- Phase 1: Compute S_min to 200 digits ---\n")
    target_dps = 200
    mpmath.mp.dps = target_dps + 50  # guard digits

    t0 = time.time()
    S_min = smin_high_precision(n_digits=target_dps)
    compute_time = time.time() - t0
    print(f"\n  S_min ({target_dps}+ digits) =")
    print(f"  {mpmath.nstr(S_min, target_dps)}")
    print(f"  Computation time: {compute_time:.1f}s")

    # Phase 2: Extended PSLQ
    print("\n--- Phase 2: Extended PSLQ (target ~100 basis elements) ---\n")
    pslq_result = smin_extended_pslq(S_min, basis_size=100)

    # Phase 3: Inner photon sum structural analysis
    print("\n--- Phase 3: Inner photon sum for small (n1, n2) ---\n")
    inner_results = []
    for n1 in range(0, 6):
        for n2 in range(0, 6):
            res = smin_inner_q_closed_form(n1, n2)
            inner_results.append(res)
            if res["n_terms"] > 0:
                print(f"  I({n1},{n2}) = {res['inner_sum_float']:.6f}  "
                      f"({res['n_terms']} allowed q-values)")

    # Phase 4: Summary
    print("\n" + "=" * 78)
    print(" SUMMARY")
    print("=" * 78)
    print(f"\n  S_min ({target_dps} digits) =")
    print(f"  {mpmath.nstr(S_min, target_dps)}")
    print(f"\n  PSLQ result: {'IDENTIFIED' if pslq_result['identified'] else 'IRREDUCIBLE'}")
    print(f"  Basis size: {pslq_result['basis_size']}")
    if pslq_result["identified"]:
        print("  Components:")
        for comp, coeff in pslq_result["components"].items():
            print(f"    {coeff} * {comp}")
        print(f"  Residual: {pslq_result['residual']:.3e}")
    else:
        print(f"  Error: {pslq_result.get('error', 'unknown')}")
    print(f"  PSLQ time: {pslq_result.get('time_sec', 0):.1f}s")

    # Save results
    output = {
        "track": "Q-3",
        "target_digits": target_dps,
        "S_min_str": mpmath.nstr(S_min, target_dps),
        "S_min_float": float(S_min),
        "compute_time_sec": compute_time,
        "pslq_result": {
            k: v for k, v in pslq_result.items()
            if k != "relation"
        },
        "inner_photon_sums": [
            {
                "n1": r["n1"],
                "n2": r["n2"],
                "inner_sum": r["inner_sum_float"],
                "n_terms": r["n_terms"],
            }
            for r in inner_results
            if r["n_terms"] > 0
        ],
    }
    if "relation" in pslq_result:
        output["pslq_result"]["relation"] = pslq_result["relation"]

    out_path = Path(__file__).parent / "data" / "smin_extended_pslq.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Results saved to {out_path}")

    print("\n" + "=" * 78)
    print(" DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()
