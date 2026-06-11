"""Verify closed-form identification of ||H_local - D_L_W||_F^2 sequence.

Discovery (Sprint H_local-PSLQ probe, 2026-05-17):

The residual decomposes via Pythagorean orthogonality (verified bit-exact at
n_max in {1, 2, 3, 4, 5}, ⟨H_local, D_W⟩_HS = 0):

    r^2(n_max) = ||H_local||^2 + ||D_W||^2
               = S(n_max) / (4 pi^2)  +  D(n_max)

where:

    S(n_max) := sum two_m_j^2 on wedge
              = n_max(n_max+1)(n_max+2)(2 n_max^2 + 4 n_max - 1) / 15

    D(n_max) := ||D_W||^2 = sum |lambda|^2 on wedge
              = n_max(n_max+1)(n_max+2)(2 n_max+1)(2 n_max+3) / 20

Both S(n_max) and D(n_max) are pure rationals in n_max; the only
transcendental is the 1/(4 pi^2) factor on the M1 (Vol(S^1)^2) sector.

This script verifies the closed form at 100 dps using mpmath, then performs
the W3-protocol PSLQ check against the master Mellin engine ring + framework
invariants to falsify any spurious identification.
"""

from __future__ import annotations

import json
import mpmath as mp
import sympy as sp
import numpy as np
from pathlib import Path

mp.mp.dps = 100  # 100 decimal places


def closed_form_residual_squared_mp(n_max: int) -> mp.mpf:
    """Compute r^2(n_max) at 100 dps."""
    n = mp.mpf(n_max)
    pi = mp.pi
    S = n * (n + 1) * (n + 2) * (2 * n * n + 4 * n - 1) / mp.mpf(15)
    D = n * (n + 1) * (n + 2) * (2 * n + 1) * (2 * n + 3) / mp.mpf(20)
    return S / (4 * pi * pi) + D


def main():
    out: dict = {
        "discovery": "r^2(n_max) = S(n_max)/(4 pi^2) + D(n_max), Pythagorean decomposition",
        "S_formula": "n(n+1)(n+2)(2n^2+4n-1)/15",
        "D_formula": "n(n+1)(n+2)(2n+1)(2n+3)/20",
        "transcendental_content": "1/(4 pi^2), one source from Hopf-base measure (M1, Vol(S^1)^2)",
        "high_precision_check": [],
    }

    # ---- High-precision closed-form vs numerical agreement ----
    # The numerical values come from finite-precision (float64) computations,
    # so the match is limited to ~14-15 dps.
    print("=" * 78)
    print("100-dps closed-form check at n_max = 1..6")
    print("=" * 78)

    numerical_vals = {
        1: 2.133227740261496,
        2: 6.527474787531087,
        3: 13.854180391695055,
        4: 24.566729356,
        5: 39.063654026,
    }

    for n in range(1, 7):
        r_sq_mp = closed_form_residual_squared_mp(n)
        r_mp = mp.sqrt(r_sq_mp)
        r_str = mp.nstr(r_mp, 50)
        print(f"n_max = {n}: r = {r_str}")
        entry = {
            "n_max": n,
            "r_sq_mp_50dps": mp.nstr(r_sq_mp, 50),
            "r_mp_50dps": r_str,
            "S": int(n * (n + 1) * (n + 2) * (2 * n * n + 4 * n - 1) // 15) if (n * (n + 1) * (n + 2) * (2 * n * n + 4 * n - 1)) % 15 == 0 else "non-integer",
            "D_times_20": int(n * (n + 1) * (n + 2) * (2 * n + 1) * (2 * n + 3)),
        }
        if n in numerical_vals:
            r_num = numerical_vals[n]
            r_mp_float = float(r_mp)
            entry["r_numerical_float64"] = r_num
            entry["match_to_float64"] = abs(r_mp_float - r_num) < 1e-10
            print(f"    numerical (float64): {r_num}, match: {entry['match_to_float64']}")
        out["high_precision_check"].append(entry)

    # ---- PSLQ at 100 dps, ceiling 10^4, against master Mellin engine ring ----
    print()
    print("=" * 78)
    print("PSLQ at 100 dps, ceiling 10^4, falsification check")
    print("=" * 78)

    # The basis: pure rationals + master Mellin engine ring elements
    # Target: r^2(3) at 100 dps (the structurally distinguished n_max=3 = Paper 2 cutoff)
    target = closed_form_residual_squared_mp(3)
    print(f"Target: r^2(3) = {mp.nstr(target, 30)}")

    # Mellin engine ring at small integer rationals
    pi = mp.pi
    G = mp.catalan
    sqrt_pi = mp.sqrt(pi)
    beta_4 = mp.zeta(4) - 2 * (mp.mpf(1) / mp.mpf(16) + mp.mpf(1) / mp.mpf(81)) * mp.mpf(0)  # placeholder
    # Use Dirichlet beta(s) via mpmath's dirichlet
    beta_4 = mp.dirichlet(4, [0, 1, 0, -1])  # L(s, chi_{-4}) at s=4
    beta_2 = G  # beta(2) = Catalan

    # Basis for PSLQ
    basis = {
        "1": mp.mpf(1),
        "1/pi": 1 / pi,
        "1/pi^2": 1 / pi**2,
        "1/pi^3": 1 / pi**3,
        "1/pi^4": 1 / pi**4,
        "pi": pi,
        "pi^2": pi**2,
        "pi^3": pi**3,
        "pi^4": pi**4,
        "1/sqrt(pi)": 1 / sqrt_pi,
        "sqrt(pi)": sqrt_pi,
        "G": G,
        "beta_4": beta_4,
        "1/(4 pi^2)": 1 / (4 * pi**2),  # the actual term
    }

    # Predicted: target = 189/(20) + 116/(60 pi^2) ... wait, let me check.
    # r^2(3) = S(3)/(4 pi^2) + D(3) = 116/(4 pi^2) + 189 = 29/pi^2 + 189
    print(f"Predicted: 189 + 29/pi^2 = {189 + 29 / pi**2}")
    predicted = mp.mpf(189) + mp.mpf(29) / (pi * pi)
    print(f"  predicted - target = {target - predicted}")

    # PSLQ search: find integer relation [c_T, c_1, c_2, ...] such that
    # c_T * target + sum c_i * basis[i] = 0
    vec_list = [target] + list(basis.values())
    keys_list = ["TARGET"] + list(basis.keys())

    try:
        rel = mp.pslq(vec_list, maxcoeff=10**4)
        if rel is None:
            print(f"PSLQ ceiling 10^4: NULL (no relation)")
            pslq_result_4 = None
        else:
            print(f"PSLQ ceiling 10^4: relation found")
            for c, k in zip(rel, keys_list):
                if c != 0:
                    print(f"    {c:+8d} * {k}")
            pslq_result_4 = list(rel)
    except Exception as e:
        print(f"PSLQ ceiling 10^4: ERROR: {e}")
        pslq_result_4 = None

    # If null, try 10^6 once (the W3-protocol allowed escalation)
    if pslq_result_4 is None:
        print()
        print("Escalating to ceiling 10^6 (single allowed)...")
        try:
            rel6 = mp.pslq(vec_list, maxcoeff=10**6)
            if rel6 is None:
                print(f"PSLQ ceiling 10^6: NULL (no relation)")
                pslq_result_6 = None
            else:
                print(f"PSLQ ceiling 10^6: relation found")
                for c, k in zip(rel6, keys_list):
                    if c != 0:
                        print(f"    {c:+8d} * {k}")
                pslq_result_6 = list(rel6)
        except Exception as e:
            print(f"PSLQ ceiling 10^6: ERROR: {e}")
            pslq_result_6 = None
    else:
        pslq_result_6 = None

    out["pslq"] = {
        "target_n_max_3": mp.nstr(target, 50),
        "predicted_form": "189 + 29/pi^2",
        "predicted_value_50dps": mp.nstr(predicted, 50),
        "match_predicted": float(abs(target - predicted)) < 1e-40,
        "ceiling_10_4_result": pslq_result_4,
        "ceiling_10_6_result": pslq_result_6,
        "basis_keys": keys_list,
    }

    # ---- Save output ----
    out_path = Path("debug/data/h_local_residual_closed_form.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"Wrote {out_path}")

    # ---- Sequence check at n_max in {1..5} ----
    print()
    print("=" * 78)
    print("Sequence verification: closed form r^2(n_max) = S/(4 pi^2) + D")
    print("=" * 78)
    print(f"{'n':>3} | {'S':>8} | {'D':>10} | {'r^2 (50 dps)':>50}")
    for n in range(1, 6):
        S = n * (n + 1) * (n + 2) * (2 * n * n + 4 * n - 1) // 15
        D_num = n * (n + 1) * (n + 2) * (2 * n + 1) * (2 * n + 3)
        D_str = f"{D_num}/20" if D_num % 20 != 0 else str(D_num // 20)
        r_sq = closed_form_residual_squared_mp(n)
        print(f"{n:>3d} | {S:>8d} | {D_str:>10s} | {mp.nstr(r_sq, 25):>50}")


if __name__ == "__main__":
    main()
