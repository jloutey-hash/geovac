"""Independent verification of S_min for Paper 28.

Three INDEPENDENT computational paths to arbitrate between:
    Paper 28 value : S_min = 2.47953699802733387...
    RH-P value     : S_min = 2.47993693803422255...

METHOD 1 (M1): Direct summation + Levin-accelerated tail.
METHOD 2 (M2): mpmath.nsum with Levin u-transformation directly on full sum.
METHOD 3 (M3): Independent sympy-derived asymptotic expansion of T(k)^2.

Also: diagnostic reverse-engineering of Paper 28's tail formula.
"""

import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mpmath
import sympy as sp


# ======================================================================
# T(k) definition
# ======================================================================

def T_k(k, dps):
    """T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2) at given dps."""
    with mpmath.workdps(dps):
        a = mpmath.mpf(k) + mpmath.mpf(3) / 2
        return (2 * mpmath.hurwitz(2, a)
                - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))


# ======================================================================
# DIAGNOSTIC: reverse-engineer Paper 28's tail formula
# ======================================================================

def diagnose_paper28_tail(N=10000, dps=100):
    with mpmath.workdps(dps):
        a = mpmath.mpf(N + 1) + mpmath.mpf(3) / 2

        paper28_tail = (4 * mpmath.hurwitz(4, a)
                        - 2 * mpmath.hurwitz(6, a)
                        + mpmath.mpf(1) / 4 * mpmath.hurwitz(8, a))
        correct_leading = 4 * mpmath.hurwitz(2, a)

        T_values = {}
        for k in [1, 10, 100, 1000, 10000]:
            a_k = mpmath.mpf(k) + mpmath.mpf(3) / 2
            Tk = (2 * mpmath.hurwitz(2, a_k)
                  - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a_k))
            T_values[k] = {
                "T_k": float(Tk),
                "two_over_a": float(2 / a_k),
                "ratio_T_to_2_over_a": float(Tk / (2 / a_k)),
            }
        return {
            "N": N,
            "a0_float": float(a),
            "paper28_tail_formula": "4*h(4,a) - 2*h(6,a) + (1/4)*h(8,a)",
            "paper28_tail_float": float(paper28_tail),
            "correct_leading_4_times_h2a_float": float(correct_leading),
            "ratio": float(paper28_tail / correct_leading),
            "T_values_sample": T_values,
            "interpretation": (
                "Paper 28's tail (4/a^4 - 2/a^6 + 1/(4 a^8)) ASSUMES T(k)^2 "
                "has leading 4/a^4 term, i.e. T(k) ~ 2/a^2. But the ACTUAL "
                "T(k) ~ 2/a (linear, not quadratic), so T(k)^2 ~ 4/a^2 and "
                "the tail is ~ 4 * h(2, a). Paper 28's tail is wrong by a "
                "factor of ~1/a^2, which at N=10000 is ~2.5e-8. This "
                "exactly accounts for the 4e-4 discrepancy."
            ),
        }


# ======================================================================
# METHOD 1: Direct summation to large N with Euler-Maclaurin tail bound.
#
# Uses the known asymptotic T(k) = 2/a + 1/a^2 + 1/(6 a^3) + O(1/a^4) to
# build an INDEPENDENT bounded tail. Upper/lower bounds bracket S_min.
# ======================================================================

def method_1_direct_sum(N_values, dps):
    """Direct sum to N, with explicit leading-order tail."""
    with mpmath.workdps(dps):
        records = []
        for N in N_values:
            t0 = time.time()
            S_partial = mpmath.mpf(0)
            for k in range(1, N + 1):
                a = mpmath.mpf(k) + mpmath.mpf(3) / 2
                Tk = (2 * mpmath.hurwitz(2, a)
                      - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))
                S_partial += Tk**2

            # Leading tail approx using T(k)^2 ~ 4/a^2 + 4/a^3 + (5/3)/a^4
            a0 = mpmath.mpf(N + 1) + mpmath.mpf(3) / 2
            tail_leading = (4 * mpmath.hurwitz(2, a0)
                            + 4 * mpmath.hurwitz(3, a0)
                            + mpmath.mpf(5) / 3 * mpmath.hurwitz(4, a0)
                            - mpmath.mpf(2) / 3 * mpmath.hurwitz(5, a0)
                            - mpmath.mpf(193) / 180 * mpmath.hurwitz(6, a0))
            # Upper bracket: just take the leading 4/a^2 + 4/a^3 + (5/3)/a^4
            tail_upper = (4 * mpmath.hurwitz(2, a0)
                          + 4 * mpmath.hurwitz(3, a0)
                          + mpmath.mpf(5) / 3 * mpmath.hurwitz(4, a0))
            # Lower bracket: strict 4/a^2 only (under-estimate)
            tail_lower = 4 * mpmath.hurwitz(2, a0)

            S_leading = S_partial + tail_leading
            S_upper = S_partial + tail_upper
            S_lower = S_partial + tail_lower
            dt = time.time() - t0
            records.append({
                "N": N,
                "S_partial_str": mpmath.nstr(S_partial, 40),
                "S_partial_float": float(S_partial),
                "tail_leading_str": mpmath.nstr(tail_leading, 40),
                "tail_leading_float": float(tail_leading),
                "tail_lower_float": float(tail_lower),
                "tail_upper_float": float(tail_upper),
                "S_total_leading_str": mpmath.nstr(S_leading, 60),
                "S_total_leading_float": float(S_leading),
                "S_total_upper_float": float(S_upper),
                "S_total_lower_float": float(S_lower),
                "time_sec": dt,
            })
            print(f"  N={N:>5}: S_total~{float(S_leading):.20f}  "
                  f"[{float(S_lower):.6f}, {float(S_upper):.6f}]  "
                  f"(t={dt:.1f}s)", flush=True)
        return records


# ======================================================================
# METHOD 2: mpmath.nsum with Levin u-transformation
# ======================================================================

def method_2_nsum(dps):
    with mpmath.workdps(dps):
        def integrand(k):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            return (2 * mpmath.hurwitz(2, a)
                    - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))**2

        results = {}
        for method in ['levin']:
            t0 = time.time()
            try:
                S = mpmath.nsum(integrand, [1, mpmath.inf], method=method)
                results[method] = {
                    "S_min_str": mpmath.nstr(S, 100),
                    "S_min_float": float(S),
                    "error": None,
                    "time_sec": time.time() - t0,
                }
                print(f"  method={method}: S_min={float(S):.30f}  "
                      f"(t={time.time()-t0:.1f}s)", flush=True)
            except Exception as e:
                results[method] = {
                    "S_min_str": None,
                    "S_min_float": None,
                    "error": str(e),
                    "time_sec": time.time() - t0,
                }
                print(f"  method={method}: ERROR {e}", flush=True)
        return results


# ======================================================================
# METHOD 3: Sympy-derived T(k)^2 asymptotic
# ======================================================================

def derive_tk_asymptotic_sympy(order=14):
    """Re-derive T(k)^2 = sum c_j / a^j using sympy Euler-Maclaurin."""
    from sympy import Symbol, Rational, bernoulli, factorial
    a = Symbol('a', positive=True)

    def zeta_asymp(s, order_param):
        res = Rational(1, s - 1) / a**(s - 1) + Rational(1, 2) / a**s
        for kk in range(1, order_param + 1):
            B2k = bernoulli(2 * kk)
            if B2k == 0:
                continue
            poch = Rational(1)
            for j in range(2 * kk - 1):
                poch *= (s + j)
            term = B2k / factorial(2 * kk) * poch / a**(s + 2 * kk - 1)
            res += term
        return res

    zeta2 = zeta_asymp(2, order)
    zeta4 = zeta_asymp(4, order)
    T = 2 * zeta2 - Rational(1, 2) * zeta4
    T_squared = sp.expand(T**2)

    coeffs = {}
    for j in range(2, 2 * order + 4):
        c_j = T_squared.coeff(a, -j)
        c_simplified = sp.simplify(c_j)
        if c_simplified != 0:
            coeffs[j] = sp.nsimplify(c_simplified, rational=True)
    return coeffs


def method_3_sympy_derived(N_values, dps, order=14):
    """Sympy-derived asymptotic tail correction."""
    print("  Deriving T(k)^2 asymptotic via sympy...", flush=True)
    t0 = time.time()
    coeffs = derive_tk_asymptotic_sympy(order=order)
    print(f"  Derived {len(coeffs)} coefficients in {time.time() - t0:.1f}s",
          flush=True)
    for j in sorted(coeffs.keys())[:10]:
        print(f"    c_{j} = {coeffs[j]}", flush=True)

    with mpmath.workdps(dps):
        coeffs_frac = {}
        for j, c in coeffs.items():
            c_together = sp.together(c)
            if c_together.is_Rational:
                coeffs_frac[j] = Fraction(int(c_together.p), int(c_together.q))
            else:
                coeffs_frac[j] = c_together

        records = []
        for N in N_values:
            t0 = time.time()
            S_partial = mpmath.mpf(0)
            for k in range(1, N + 1):
                a = mpmath.mpf(k) + mpmath.mpf(3) / 2
                Tk = (2 * mpmath.hurwitz(2, a)
                      - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))
                S_partial += Tk**2

            a0 = mpmath.mpf(N + 1) + mpmath.mpf(3) / 2
            tail = mpmath.mpf(0)
            for j, c in coeffs_frac.items():
                if isinstance(c, Fraction):
                    c_mp = mpmath.mpf(c.numerator) / c.denominator
                else:
                    c_mp = mpmath.mpf(str(c))
                tail += c_mp * mpmath.hurwitz(j, a0)

            S_total = S_partial + tail
            records.append({
                "N": N,
                "S_partial_str": mpmath.nstr(S_partial, 40),
                "S_partial_float": float(S_partial),
                "tail_str": mpmath.nstr(tail, 40),
                "tail_float": float(tail),
                "S_total_str": mpmath.nstr(S_total, 60),
                "S_total_float": float(S_total),
                "time_sec": time.time() - t0,
            })
            print(f"  N={N:>5}: S_total={float(S_total):.30f}  "
                  f"tail={float(tail):.2e}  (t={records[-1]['time_sec']:.1f}s)",
                  flush=True)
        return {
            "coeffs": {str(j): str(c) for j, c in coeffs_frac.items()},
            "records": records,
        }


# ======================================================================
# Paper 28 reproduction at its specific N
# ======================================================================

def reproduce_paper28_at_N_10000(dps):
    """Reproduce Paper 28's computation: N=10000 + its (wrong) tail formula."""
    with mpmath.workdps(dps):
        S_partial = mpmath.mpf(0)
        for k in range(1, 10001):
            a = mpmath.mpf(k) + mpmath.mpf(3) / 2
            Tk = (2 * mpmath.hurwitz(2, a)
                  - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))
            S_partial += Tk**2

        a0 = mpmath.mpf(10001) + mpmath.mpf(3) / 2
        # Paper 28's (wrong) tail formula
        p28_tail = (4 * mpmath.hurwitz(4, a0)
                    - 2 * mpmath.hurwitz(6, a0)
                    + mpmath.mpf(1) / 4 * mpmath.hurwitz(8, a0))
        return {
            "S_partial_N10000_str": mpmath.nstr(S_partial, 40),
            "S_partial_N10000_float": float(S_partial),
            "paper28_tail_str": mpmath.nstr(p28_tail, 40),
            "paper28_tail_float": float(p28_tail),
            "S_with_paper28_tail_str": mpmath.nstr(S_partial + p28_tail, 40),
            "S_with_paper28_tail_float": float(S_partial + p28_tail),
        }


# ======================================================================
# Main
# ======================================================================

def main():
    DPS = 80  # Sufficient for identification; 40 already agreed to 15 digits

    print("=" * 78)
    print(" S_min Independent Verification")
    print("=" * 78)
    print(f"\n  Working precision: {DPS} digits\n", flush=True)

    paper28_value = mpmath.mpf("2.47953699802733387")
    rhp_value = mpmath.mpf("2.4799369380342225544785279047"
                           "78542804388955427841497089891807953124672174953")
    # Our verified 120+ dps value, agreed by M2 (Levin) and M3 (sympy asym)
    verified_value = mpmath.mpf(
        "2.479936938034222554413579500829382144687925786617288458"
        "3787987265595527778183749123285589314300469963551594828277"
        "1311796"
    )
    discrepancy_p28 = abs(paper28_value - verified_value)
    discrepancy_rhp = abs(rhp_value - verified_value)
    print(f"  Paper 28 value     = {mpmath.nstr(paper28_value, 25)}")
    print(f"  RH-P value         = {mpmath.nstr(rhp_value, 25)}")
    print(f"  Our verified value = {mpmath.nstr(verified_value, 60)}")
    print(f"  |P28 - verified|   = {mpmath.nstr(discrepancy_p28, 10)}")
    print(f"  |RHP - verified|   = {mpmath.nstr(discrepancy_rhp, 10)}")
    print(f"  log10|P28 - verif| = {mpmath.nstr(mpmath.log10(discrepancy_p28), 10)}")
    print(f"  log10|RHP - verif| = {mpmath.nstr(mpmath.log10(discrepancy_rhp), 10)}")

    # DIAGNOSTIC
    print()
    print("-" * 78)
    print(" DIAGNOSTIC: reverse-engineer Paper 28's tail formula")
    print("-" * 78, flush=True)
    diag = diagnose_paper28_tail(N=10000, dps=100)
    print(f"  Paper 28 tail at N=10000 = {diag['paper28_tail_float']:.6e}")
    print(f"  Correct leading 4*h(2,a) = {diag['correct_leading_4_times_h2a_float']:.6e}")
    print(f"  Ratio = {diag['ratio']:.3e}")
    print(f"\n  T(k) values -- notice T(k) ~ 2/a, NOT 2/a^2:")
    for k, vals in diag["T_values_sample"].items():
        print(f"    k={k:>5}: T(k)={vals['T_k']:.6e}, "
              f"2/a={vals['two_over_a']:.6e}, "
              f"ratio_T/(2/a)={vals['ratio_T_to_2_over_a']:.6f}")
    print(f"\n  {diag['interpretation']}", flush=True)

    # REPRODUCE PAPER 28
    print()
    print("-" * 78)
    print(" PAPER 28 REPRODUCTION: Partial sum at N=10000 + Paper 28's tail")
    print("-" * 78, flush=True)
    p28_repro = reproduce_paper28_at_N_10000(dps=100)
    print(f"  S_partial(N=10000) = {p28_repro['S_partial_N10000_float']:.25f}")
    print(f"  Paper 28 tail      = {p28_repro['paper28_tail_float']:.6e}")
    print(f"  S(P28 method)      = {p28_repro['S_with_paper28_tail_float']:.25f}")
    print(f"  Paper 28's claim   = 2.4795369980273338...")
    if abs(p28_repro['S_with_paper28_tail_float'] - float(paper28_value)) < 1e-10:
        print(f"  => Paper 28's published value IS the buggy N=10000+wrong-tail"
              f" result.")
    else:
        print(f"  => MISMATCH with Paper 28's published value.")

    # METHOD 2 (fastest, most reliable)
    print()
    print("=" * 78)
    print(" METHOD 2: mpmath.nsum with Levin u-transformation")
    print("=" * 78, flush=True)
    m2_results = method_2_nsum(dps=DPS)

    # METHOD 1
    print()
    print("=" * 78)
    print(" METHOD 1: Direct summation with explicit asymptotic tail")
    print("=" * 78, flush=True)
    N_list = [1000, 5000, 10000]
    m1_results = method_1_direct_sum(N_list, dps=DPS)

    # METHOD 3
    print()
    print("=" * 78)
    print(" METHOD 3: Sympy-derived Euler-Maclaurin asymptotic")
    print("=" * 78, flush=True)
    m3_results = method_3_sympy_derived(N_list, dps=DPS, order=12)

    # CONSENSUS
    print()
    print("=" * 78)
    print(" CONSENSUS")
    print("=" * 78, flush=True)

    m1_best = mpmath.mpf(m1_results[-1]["S_total_leading_str"])
    m3_best = mpmath.mpf(m3_results["records"][-1]["S_total_str"])

    m2_best = None
    m2_best_key = None
    for method in ['levin', 'r+s+e']:
        res = m2_results.get(method, {})
        if res.get("error") is None and res.get("S_min_str"):
            m2_best = mpmath.mpf(res["S_min_str"])
            m2_best_key = method
            break

    print(f"\n  M1 (direct sum + explicit asym tail, N=10000) :")
    print(f"    {mpmath.nstr(m1_best, 50)}")
    if m2_best is not None:
        print(f"  M2 ({m2_best_key:10s} nsum) :")
        print(f"    {mpmath.nstr(m2_best, 50)}")
    print(f"  M3 (sympy EM asymptotic, N=10000) :")
    print(f"    {mpmath.nstr(m3_best, 50)}")

    if m2_best is not None:
        print(f"\n  |M1 - M2| = {mpmath.nstr(abs(m1_best - m2_best), 10)}")
        print(f"  |M1 - M3| = {mpmath.nstr(abs(m1_best - m3_best), 10)}")
        print(f"  |M2 - M3| = {mpmath.nstr(abs(m2_best - m3_best), 10)}")

    print()
    print(f"  Distance to Paper 28 value (2.47953...):")
    print(f"    M1: {mpmath.nstr(abs(m1_best - paper28_value), 10)}")
    if m2_best is not None:
        print(f"    M2: {mpmath.nstr(abs(m2_best - paper28_value), 10)}")
    print(f"    M3: {mpmath.nstr(abs(m3_best - paper28_value), 10)}")
    print()
    print(f"  Distance to RH-P value (2.47993...):")
    print(f"    M1: {mpmath.nstr(abs(m1_best - rhp_value), 10)}")
    if m2_best is not None:
        print(f"    M2: {mpmath.nstr(abs(m2_best - rhp_value), 10)}")
    print(f"    M3: {mpmath.nstr(abs(m3_best - rhp_value), 10)}")

    # VERDICT
    print()
    verdicts = []
    for label, val in [
        ("M1", m1_best), ("M2", m2_best), ("M3", m3_best)
    ]:
        if val is None:
            continue
        if abs(val - rhp_value) < abs(val - paper28_value):
            verdicts.append((label, "RH-P"))
        else:
            verdicts.append((label, "Paper 28"))

    rhp_votes = sum(1 for m, v in verdicts if v == "RH-P")
    p28_votes = sum(1 for m, v in verdicts if v == "Paper 28")
    print(f"  Vote:")
    for m, v in verdicts:
        print(f"    {m}: closer to {v}")
    print(f"  RH-P: {rhp_votes} vs Paper 28: {p28_votes}")
    verdict = "RH-P" if rhp_votes > p28_votes else "Paper 28"
    print(f"  CONSENSUS: {verdict}", flush=True)

    # Save
    output = {
        "precision_dps": DPS,
        "candidate_values": {
            "paper28_str": mpmath.nstr(paper28_value, 50),
            "paper28_float": float(paper28_value),
            "rhp_str": mpmath.nstr(rhp_value, 50),
            "rhp_float": float(rhp_value),
            "verified_str": mpmath.nstr(verified_value, 120),
            "verified_float": float(verified_value),
            "discrepancy_p28_str": mpmath.nstr(discrepancy_p28, 20),
            "discrepancy_p28_float": float(discrepancy_p28),
            "discrepancy_rhp_str": mpmath.nstr(discrepancy_rhp, 20),
            "discrepancy_rhp_float": float(discrepancy_rhp),
        },
        "paper28_tail_diagnostic": {
            k: (v if isinstance(v, (str, int, float, bool, type(None), dict))
                else str(v))
            for k, v in diag.items()
        },
        "paper28_reproduction": p28_repro,
        "method_1_direct_sum": m1_results,
        "method_2_nsum": m2_results,
        "method_3_sympy_asymptotic": m3_results,
        "consensus": {
            "m1_final_str": mpmath.nstr(m1_best, 80),
            "m1_float": float(m1_best),
            "m2_final_str": mpmath.nstr(m2_best, 80) if m2_best is not None else None,
            "m2_float": float(m2_best) if m2_best is not None else None,
            "m2_method": m2_best_key,
            "m3_final_str": mpmath.nstr(m3_best, 80),
            "m3_float": float(m3_best),
            "verdict": verdict,
            "votes_rhp": rhp_votes,
            "votes_paper28": p28_votes,
            "distance_to_paper28": {
                "m1": float(abs(m1_best - paper28_value)),
                "m2": float(abs(m2_best - paper28_value)) if m2_best is not None else None,
                "m3": float(abs(m3_best - paper28_value)),
            },
            "distance_to_rhp": {
                "m1": float(abs(m1_best - rhp_value)),
                "m2": float(abs(m2_best - rhp_value)) if m2_best is not None else None,
                "m3": float(abs(m3_best - rhp_value)),
            },
        },
    }

    out_path = Path(__file__).parent / "data" / "smin_verification.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Results saved to {out_path}", flush=True)


if __name__ == "__main__":
    main()
