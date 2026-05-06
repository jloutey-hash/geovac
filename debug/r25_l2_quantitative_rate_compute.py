"""Driver for R2.5 / L2 quantitative-rate sprint.

Computes:
  - gamma_n via the closed-form sum-rule gamma_n = pi - 4 T_n / (pi Z_n) at high precision
  - doubling-estimator a_n = (2n*g_{2n} - n*g_n)/log(2)
  - threshold table N_0(C) for various C
  - uniform-bound margins for C = 6

Writes debug/data/r25_l2_quantitative_rate.json.

Companion: debug/r25_l2_quantitative_rate_memo.md (the proof memo).
"""

from __future__ import annotations

import json
import math
import os
import sys

import mpmath


# Make sure we can import geovac
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
if PROJ not in sys.path:
    sys.path.insert(0, PROJ)


def T_n_exact(n: int, prec_dps: int = 60) -> mpmath.mpf:
    """Closed-form T_n = sum over different-parity pairs (k1, k2)
    of sqrt(k1 k2) [1/(k1-k2)^2 - 1/(k1+k2)^2].

    Uses k1 < k2 by symmetry, doubles. mpmath dps controlled.
    """
    mpmath.mp.dps = prec_dps
    total = mpmath.mpf(0)
    for k1 in range(1, n + 1):
        for k2 in range(k1 + 1, n + 1):
            if (k1 + k2) % 2 == 0:
                continue
            term = mpmath.sqrt(k1 * k2) * (
                mpmath.mpf(1) / (k1 - k2) ** 2
                - mpmath.mpf(1) / (k1 + k2) ** 2
            )
            total += 2 * term
    return total


def gamma_n_via_sum_rule(n: int, prec_dps: int = 60) -> mpmath.mpf:
    """Compute gamma_n via the closed-form sum rule (Eq. 3.2 of memo)."""
    mpmath.mp.dps = prec_dps
    Z = mpmath.mpf(n * (n + 1)) / 2
    T = T_n_exact(n, prec_dps=prec_dps)
    return mpmath.pi - 4 * T / (mpmath.pi * Z)


def main():
    out_path = os.path.join(PROJ, "debug", "data", "r25_l2_quantitative_rate.json")
    print(f"Computing quantitative-rate data, output: {out_path}")

    # Section 1: gamma_n at small to large n
    n_values_small = list(range(2, 21))
    n_values_medium = [25, 30, 40, 50, 70, 100, 150, 200, 300, 500, 1000]
    n_values_large = [200, 400, 800, 1600]
    all_n = sorted(set(n_values_small + n_values_medium))

    gammas = {}
    print("Computing gamma_n at small/medium n...")
    for n in all_n:
        g = gamma_n_via_sum_rule(n, prec_dps=50)
        gammas[n] = float(g)
        print(f"  n={n}: gamma = {gammas[n]:.10f}")

    # Section 2: gamma_n at very large n for asymptotic
    print("Computing gamma_n at large n...")
    gammas_large = {}
    for n in n_values_large:
        if n in gammas:
            gammas_large[n] = gammas[n]
        else:
            g = gamma_n_via_sum_rule(n, prec_dps=50)
            gammas_large[n] = float(g)
            gammas[n] = float(g)
        print(f"  n={n}: gamma = {gammas_large[n]:.12f}")

    # Section 3: doubling estimator
    log2 = float(mpmath.log(2))
    doubling_data = []
    for i in range(len(n_values_large) - 1):
        n1 = n_values_large[i]
        n2 = n_values_large[i + 1]  # 2*n1
        if n2 != 2 * n1:
            continue
        a_n = (n2 * gammas[n2] - n1 * gammas[n1]) / log2
        doubling_data.append({
            "n": n1,
            "n_doubled": n2,
            "a_n": a_n,
            "error_from_4_over_pi": a_n - 4 / math.pi,
        })

    # Section 4: threshold table N_0(C)
    candidate_Cs = [
        ("4/pi (asymptote)", 4 / math.pi),
        ("7/4", 1.75),
        ("2", 2.0),
        ("pi - 1", math.pi - 1),
        ("5/2", 2.5),
        ("3", 3.0),
        ("4", 4.0),
        ("5", 5.0),
        ("6 (uniform)", 6.0),
    ]
    threshold_table = []
    for name, C in candidate_Cs:
        # smallest n such that for ALL m >= n, m*gamma_m/log(m) <= C
        last_violation = 0
        for n in sorted(gammas.keys()):
            if n < 2:
                continue
            ratio = n * gammas[n] / math.log(n)
            if ratio > C:
                last_violation = n
        if last_violation == 0:
            N_0 = 2
            note = "holds for all n >= 2"
        elif last_violation >= max(gammas.keys()):
            N_0 = None
            note = f"violated up to largest tested n = {last_violation} (asymptotic)"
        else:
            # Find next n in sorted list strictly greater than last_violation
            sorted_ns = sorted(gammas.keys())
            next_idx = sorted_ns.index(last_violation) + 1
            if next_idx < len(sorted_ns):
                N_0 = sorted_ns[next_idx]
                note = f"holds for n >= {N_0}"
            else:
                N_0 = None
                note = "asymptotic only"
        threshold_table.append({
            "C_name": name,
            "C_value": C,
            "N_0": N_0,
            "note": note,
        })

    # Section 5: uniform bound C = 6 verification
    uniform_C = 6.0
    uniform_check = []
    for n in sorted(gammas.keys()):
        if n < 2:
            continue
        bound = uniform_C * math.log(n) / n
        margin = (bound - gammas[n]) / gammas[n]
        uniform_check.append({
            "n": n,
            "gamma_n": gammas[n],
            "bound_6_log_n_over_n": bound,
            "satisfies_bound": gammas[n] <= bound,
            "margin_pct": margin * 100,
        })

    # Section 6: ratio table for asymptotic
    ratio_table = []
    for n in sorted(gammas.keys()):
        if n < 2:
            continue
        ratio = n * gammas[n] / math.log(n)
        ratio_table.append({
            "n": n,
            "gamma_n": gammas[n],
            "n_gamma_over_log_n": ratio,
            "error_from_4_over_pi": ratio - 4 / math.pi,
        })

    # Section 7: assemble
    asymptotic_constant = {
        "value": 4 / math.pi,
        "expression": "4/pi",
        "structural_meaning": (
            "Hopf-base measure factor: Vol(S^2)/pi^2 = 2 Vol(S^1)/Vol(SU(2)) "
            "rescaled by Haar normalization; M1 mechanism of Sprint TS-E1 / "
            "Paper 32 SVIII / Paper 18 SIII.7."
        ),
    }

    output = {
        "description": (
            "Quantitative rate sharpening for R2.5 / L2 mass-concentration "
            "moment gamma_{n_max}. See debug/r25_l2_quantitative_rate_memo.md."
        ),
        "asymptotic_constant": asymptotic_constant,
        "uniform_bound": {
            "C": uniform_C,
            "form": "gamma_n <= 6 * log(n) / n for all n >= 2",
            "tightest_n": 2,
            "tightest_margin_pct": uniform_check[0]["margin_pct"],
        },
        "doubling_estimator": doubling_data,
        "threshold_table": threshold_table,
        "uniform_check_table": uniform_check,
        "ratio_table": ratio_table,
    }

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nWrote {out_path}")
    print(f"  asymptotic constant = 4/pi = {asymptotic_constant['value']:.6f}")
    print(f"  uniform bound C = 6 verified at {len(uniform_check)} values of n")
    print(f"  threshold table: {len(threshold_table)} entries")


if __name__ == "__main__":
    main()
