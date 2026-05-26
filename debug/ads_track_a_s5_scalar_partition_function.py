"""Track AdS-A extension to S^5 — Conformally coupled scalar partition function.

Goal
----
Test whether the master Mellin engine M2/M3 decomposition of Paper 50
extends from S^3 (CFT_3) to S^5 (CFT_4 on round 5-sphere).

Setup
-----
Round unit S^5:
  - Scalar Laplacian eigenvalues: lambda_n = n(n+4) for n = 0, 1, 2, ...
  - Multiplicity: m_n = (2n+4)(n+1)(n+2)(n+3)/6  [dim of degree-n harmonic
    polynomials on R^6, restricted to S^5]
  - Conformal mass shift: xi * R_scalar for d=5
    xi = (d-2)/(4(d-1)) = 3/16
    R_scalar (unit S^5) = d(d-1) = 20
    Shift = 3/16 * 20 = 15/4
  - Conformally coupled spectrum: (n+3/2)(n+5/2)
    [analogous to (n+1/2)(n+3/2) on S^3]

Re-index m = n + 1, m = 1, 2, ...:
  - lambda_m^2 = (m+1/2)(m+3/2) = m^2 + 2m + 3/4
                              = (m+1)^2 - 1/4
  - m_n = (2m+2)(m)(m+1)(m+2)/6 = (m+1)·m(m+1)(m+2)/3

Hmm let me recompute. Re-index n -> m via m = n+1:
  lambda_n = (n+3/2)(n+5/2) -> at m=n+1: (m+1/2)(m+3/2)
  m_n = (2n+4)(n+1)(n+2)(n+3)/6 -> at m=n+1: (2m+2)·m·(m+1)·(m+2)/6
                                            = 2(m+1)·m(m+1)(m+2)/6
                                            = m(m+1)^2(m+2)/3

So zeta_S5(s) = sum_{m=1}^infty m(m+1)^2(m+2)/3 / [(m+1/2)(m+3/2)]^s

This is a 4th-degree polynomial in m divided by (m+1/2)(m+3/2) to the s.

Cleaner: keep original n indexing.
  zeta_S5(s) = sum_{n=0}^infty m_n / [(n+3/2)(n+5/2)]^s
             = sum_{n=0}^infty (2n+4)(n+1)(n+2)(n+3)/6 / [(n+3/2)(n+5/2)]^s

For Hurwitz expansion, write u = n + 5/2 (centered at midpoint of (n+3/2)(n+5/2)):
  (n+3/2)(n+5/2) = (u-1)·u = u^2 - u
                = ... hmm doesn't factor nicely

Or use v = n + 2:
  (n+3/2)(n+5/2) = (v-1/2)(v+1/2) = v^2 - 1/4

So in v = n+2 (v = 2, 3, 4, ...):
  eigenvalue = v^2 - 1/4
  multiplicity = (2(v-2)+4)(v-1)(v)(v+1)/6 = 2v(v-1)v(v+1)/6 = v(v-1)(v+1)·2v/6 = v^2(v-1)(v+1)/3
              = v^2(v^2-1)/3

So zeta_S5_conf(s) = sum_{v=2}^infty v^2(v^2-1)/3 / (v^2 - 1/4)^s

Interesting structure! The multiplicity v^2(v^2-1)/3 starts at v=2.

For Hurwitz expansion:
  v^2(v^2-1)/3 = (v^4 - v^2)/3

At large v: (v^4 - v^2)/3 · (v^2 - 1/4)^{-s}
        = (v^4 - v^2)/3 · v^{-2s} · (1 - 1/(4v^2))^{-s}
        = (1/3) [v^{4-2s} - v^{2-2s}] · sum_k C(s+k-1, k)/(4^k v^{2k})

Summing over v >= 2 (need to handle v=0, 1 corrections):
  sum_{v=2}^infty v^{a} = zeta_R(-a) - 1 - 2^a  (for a < -1)
  ...but for a near 0 we need analytic continuation.

Easier: extend the sum to v=1 (where v^2-1 = 0, so the term vanishes) and
v=0 (where v^2(v^2-1) = 0 also vanishes). So:
  zeta_S5_conf(s) = sum_{v=0}^infty v^2(v^2-1)/3 / (v^2 - 1/4)^s
                  = sum_{v=2}^infty v^2(v^2-1)/3 / (v^2 - 1/4)^s

(v=0 and v=1 give 0 multiplicity, so no contribution. Good.)

Now we can use the same Hurwitz machinery as Paper 50:
  v^a · (v^2 - 1/4)^{-s} expansion gives Hurwitz/Riemann zeta at shifted args

For F_S^5 = -(1/2) zeta_S5_conf'(0).

Numerical strategy
------------------
1. Compute zeta_S5_conf(s) and its derivative at s=0 via Hurwitz expansion to k_max ~ 100, high dps
2. PSLQ against {log 2, zeta(3)/pi^2, zeta(5)/pi^4, pi^2, pi^4, 1} basis
3. Identify if it matches known KPS S^5 closed form
4. Check master Mellin engine M2/M3 decomposition extends or breaks
"""

import mpmath as mp
import sympy as sp
import json
from pathlib import Path


def zeta_S5_conf_prime_at_zero_hurwitz(k_max: int = 100, dps: int = 200) -> mp.mpf:
    """zeta'(0) for conformally coupled scalar on round S^5 via Hurwitz expansion.

    zeta(s) = sum_{v=2}^infty (v^4 - v^2)/3 · (v^2 - 1/4)^{-s}

    Expand (v^2 - 1/4)^{-s} via binomial series:
        (v^2 - 1/4)^{-s} = v^{-2s} sum_{k>=0} g_k(s) v^{-2k}

    where g_k(s) = C(s+k-1, k) / 4^k.

    Then zeta(s) = (1/3) sum_{k>=0} g_k(s) [zeta_R(2s+2k-4) - zeta_R(2s+2k-2)]

    (where the [zeta_R(2s+2k-4) - zeta_R(2s+2k-2)] comes from
     (v^4 - v^2)/v^{2s+2k} = v^{4-2s-2k} - v^{2-2s-2k} summed over v=1..infty.
     Sum starts at v=1 since v=0 gives v^0 type term not v^a; but v=0 doesn't
     contribute anyway because of v^2(v^2-1) = 0 at v=0 AND v=1.)

    Note: this Hurwitz form sums over v=1..infty (the standard Riemann form),
    not v=2..infty. But v=1 contributes 1·(1^2 - 1/4)^{-s} · (v^4 - v^2)/3 |_{v=1}
    = 0, so v=1 contributes zero anyway. v=0 same. So the Hurwitz expansion is
    valid.
    """
    mp.mp.dps = dps

    # ζ(s) = (1/3) Σ g_k(s) [ζ_R(2s + 2k - 4) - ζ_R(2s + 2k - 2)]
    #
    # At s = 0: g_0(0) = 1, g_k(0) = 0 for k >= 1.
    # So ζ(0) = (1/3) [ζ_R(-4) - ζ_R(-2)] = (1/3) [0 - 0] = 0 ✓
    #
    # ζ'(0):
    # Derivative of g_k(s) [ζ_R(2s+2k-4) - ζ_R(2s+2k-2)] at s=0
    #   = g_k'(0) [ζ_R(2k-4) - ζ_R(2k-2)]
    #   + g_k(0) · 2 [ζ_R'(2k-4) - ζ_R'(2k-2)]
    #
    # g_0(0) = 1, g_0'(0) = 0
    # g_k(0) = 0 for k >= 1, g_k'(0) = 1/(k * 4^k) for k >= 1

    # So ζ'(0) = (1/3) {
    #     g_0(0) · 2 [ζ_R'(-4) - ζ_R'(-2)]
    #     + Σ_{k>=1} g_k'(0) [ζ_R(2k-4) - ζ_R(2k-2)]
    # }
    #          = (1/3) {
    #     2 [ζ_R'(-4) - ζ_R'(-2)]
    #     + Σ_{k>=1} [ζ_R(2k-4) - ζ_R(2k-2)] / (k * 4^k)
    # }

    # First two terms: 2 [ζ_R'(-4) - ζ_R'(-2)]
    zeta_p_neg4 = mp.zeta(-4, derivative=1)
    zeta_p_neg2 = mp.zeta(-2, derivative=1)
    leading = 2 * (zeta_p_neg4 - zeta_p_neg2)

    # Identities (standard from Riemann functional equation):
    #   ζ_R'(-2) = -ζ(3)/(4π²)
    #   ζ_R'(-4) = +3 ζ(5)/(4π^4)
    # Verify:
    pred_neg2 = -mp.zeta(3) / (4 * mp.pi**2)
    pred_neg4 = 3 * mp.zeta(5) / (4 * mp.pi**4)

    err_neg2 = abs(zeta_p_neg2 - pred_neg2)
    err_neg4 = abs(zeta_p_neg4 - pred_neg4)
    print(f"  Verification: zeta_R'(-2) = -zeta(3)/(4 pi^2): err = {mp.nstr(err_neg2, 3)}")
    print(f"  Verification: zeta_R'(-4) = +3 zeta(5)/(4 pi^4): err = {mp.nstr(err_neg4, 3)}")

    # Sum: Σ_{k>=1} [ζ_R(2k-4) - ζ_R(2k-2)] / (k * 4^k)
    series_sum = mp.mpf(0)
    for k in range(1, k_max + 1):
        term = (mp.zeta(2*k - 4) - mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
        series_sum += term

    total = (leading + series_sum) / 3
    return total


def pslq_decomposition_S5(framework_value: mp.mpf, dps: int = 200) -> dict:
    """PSLQ-identify framework's zeta'(0) on S^5 against extended basis.

    Expected transcendentals on S^5: log 2, zeta(3)/pi^2, zeta(5)/pi^4
    (these are the M2 + M3 ring extended to include zeta(5)).
    """
    mp.mp.dps = dps
    target = framework_value

    basis = [
        mp.log(2),
        mp.zeta(3) / mp.pi**2,
        mp.zeta(5) / mp.pi**4,
    ]
    basis_names = ["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4"]

    vec = [target] + basis
    tol_exp = 40  # generous for the expansion
    relation = mp.pslq(vec, tol=mp.mpf('1e-' + str(tol_exp)), maxcoeff=10**6)

    if relation is None:
        return {"status": "no_relation_found", "tol": f"1e-{tol_exp}"}

    a0 = relation[0]
    if a0 == 0:
        return {"status": "trivial_relation", "raw": list(relation)}

    coeffs = {}
    for i, name in enumerate(basis_names):
        rat = sp.Rational(-relation[i+1], a0)
        coeffs[name] = str(rat)

    return {
        "status": "relation_found",
        "integer_relation": list(relation),
        "normalized_coefficients": coeffs,
    }


def main():
    print("=" * 70)
    print("Track AdS-A extension to S^5: Conformally coupled scalar")
    print("=" * 70)
    print()
    print("Setup:")
    print("  Spectrum: (n+3/2)(n+5/2) for n=0,1,2,...")
    print("  Multiplicity: (2n+4)(n+1)(n+2)(n+3)/6")
    print("  (Re-indexed: eigenvalue v^2-1/4, multiplicity v^2(v^2-1)/3 for v=2,3,...)")
    print()

    DPS = 200

    # Compute framework value at various k_max
    print("Step 1: Convergence of Hurwitz expansion")
    print(f"  {'k_max':>8}  {'zeta_S5_conf_prime(0)':>32}  {'F_S5_scalar':>32}")

    values = {}
    for k_max in (20, 40, 60, 80, 100):
        v = zeta_S5_conf_prime_at_zero_hurwitz(k_max=k_max, dps=DPS)
        F_val = -v / 2
        values[k_max] = (v, F_val)
        print(f"  {k_max:>8}  {mp.nstr(v, 25):>32}  {mp.nstr(F_val, 25):>32}")

    print()

    # PSLQ identification
    print("Step 2: PSLQ structural identification")
    best_zeta_prime = values[100][0]
    pslq_result = pslq_decomposition_S5(best_zeta_prime, dps=DPS)

    if pslq_result['status'] == "relation_found":
        print(f"  Integer relation found: {pslq_result['integer_relation']}")
        print(f"  Normalized coefficients (framework's zeta'(0) = ...):")
        for basis_name, coeff in pslq_result['normalized_coefficients'].items():
            print(f"    coeff of {basis_name:>20}: {coeff}")
    else:
        print(f"  {pslq_result['status']}")
        print(f"  Could not find clean integer relation in {{log(2), zeta(3)/pi^2, zeta(5)/pi^4}} basis.")
        print(f"  Maybe need more basis elements (zeta(2k+1)/pi^{{2k}}, log pi, etc.)")
    print()

    # Master Mellin engine reading
    print("Step 3: Master Mellin engine reading")
    print()
    print("  If the relation involves only {log(2), zeta(3)/pi^2, zeta(5)/pi^4}:")
    print("    Each is in a natural M-engine ring:")
    print("    - log(2): M2 (Seeley-DeWitt extension)")
    print("    - zeta(3)/pi^2, zeta(5)/pi^4: M3 (half-integer Hurwitz / odd-zeta)")
    print("  S^5 introduces zeta(5)/pi^4 as a new M3 transcendental (vs S^3 only zeta(3)/pi^2).")
    print("  Pattern: M3 ring on S^d has odd-zeta over even-pi powers up to zeta(d-2).")
    print()

    # Save results
    out_path = Path("debug/data/ads_track_a_s5_scalar_partition_function.json")
    out_path.parent.mkdir(exist_ok=True)

    results = {
        "track": "AdS-A extension",
        "system": "conformally coupled scalar on round unit S^5",
        "spectrum_eigenvalues": "(n+3/2)(n+5/2) for n=0,1,2,...",
        "spectrum_multiplicities": "(2n+4)(n+1)(n+2)(n+3)/6 for n=0,1,2,...",
        "dps": DPS,
        "convergence": [
            {
                "k_max": k,
                "zeta_prime_at_0": str(v[0]),
                "F_scalar_S5": str(v[1]),
            }
            for k, v in values.items()
        ],
        "pslq_decomposition": pslq_result,
    }

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()
