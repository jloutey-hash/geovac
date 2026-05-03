"""
g2_c3_symbolic_B.py -- Compute B(n_int) as exact rationals using sympy.

B(n_int) for the g-2 vertex correction on S³ is a sum of products of
CG coefficients (rational) and propagator factors (rational). Every
B(n_int) is therefore an exact rational number.

This script computes B(1), B(2), B(3) exactly to find the rational
structure, then attempts to identify c₃ algebraically.
"""

import json
import time
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational as R, pi, sqrt, nsimplify, simplify, sympify

DATA_DIR = Path(__file__).parent / "data"


@lru_cache(maxsize=2_000_000)
def cg_exact(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
    """CG coefficient as sympy Rational. Arguments are 2*j, 2*m integers."""
    if m1_2 + m2_2 != M_2:
        return R(0)
    if abs(j1_2 - j2_2) > J_2 or J_2 > j1_2 + j2_2:
        return R(0)
    if abs(m1_2) > j1_2 or abs(m2_2) > j2_2 or abs(M_2) > J_2:
        return R(0)
    val = clebsch_gordan(R(j1_2, 2), R(j2_2, 2), R(J_2, 2),
                         R(m1_2, 2), R(m2_2, 2), R(M_2, 2))
    return val


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, jsL2, jsR2, jtL2, jtR2):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL2, jgR2 in [(q + 1, q - 1), (q - 1, q + 1)]:
        if jgL2 < 0 or jgR2 < 0:
            continue
        if (abs(jsL2 - jgL2) <= jtL2 <= jsL2 + jgL2 and
                abs(jsR2 - jgR2) <= jtR2 <= jsR2 + jgR2):
            chs.append((jgL2, jgR2))
    return chs


def vertex_amp_exact(jsL2, jsR2, js2, mjs2,
                     jtL2, jtR2, jt2, mjt2,
                     jgL2, jgR2, mgL2, mgR2):
    """Vertex amplitude using exact sympy CG coefficients."""
    total = R(0)
    for mL1_2 in range(-jsL2, jsL2 + 1, 2):
        mR1_2 = mjs2 - mL1_2
        if abs(mR1_2) > jsR2:
            continue
        mL2_2 = mL1_2 + mgL2
        if abs(mL2_2) > jtL2:
            continue
        mR2_2 = mR1_2 + mgR2
        if abs(mR2_2) > jtR2:
            continue
        if mL2_2 + mR2_2 != mjt2:
            continue
        c1v = cg_exact(jsL2, jsR2, js2, mL1_2, mR1_2, mjs2)
        if c1v == 0:
            continue
        c2v = cg_exact(jtL2, jtR2, jt2, mL2_2, mR2_2, mjt2)
        if c2v == 0:
            continue
        c3v = cg_exact(jsL2, jgL2, jtL2, mL1_2, mgL2, mL2_2)
        c4v = cg_exact(jsR2, jgR2, jtR2, mR1_2, mgR2, mR2_2)
        total += c1v * c2v * c3v * c4v
    return total


def compute_B_exact(n_ext, n_int):
    """Compute B(n_int) = vertex(mj=+1/2) - vertex(mj=-1/2) as exact rational."""
    result = R(0)
    q_probe = 1

    for mj_ext_2 in [+1, -1]:
        sign = 1 if mj_ext_2 == +1 else -1

        jE_L_2 = n_ext + 1
        jE_R_2 = n_ext
        j_ext_2 = 1

        jI_L_2 = n_int + 1
        jI_R_2 = n_int

        # lambda = (2*n_int + 3)/2, lambda^4 = ((2*n_int+3)/2)^4
        lam_num = 2 * n_int + 3
        lam4 = R(lam_num**4, 16)  # (lam_num/2)^4 = lam_num^4/16

        if not vertex_allowed(n_int, n_int, q_probe):
            continue
        probe_chs = get_channels(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
        if not probe_chs:
            continue

        j_int_min_2 = 1
        j_int_max_2 = 2 * n_int + 1

        subtotal = R(0)

        for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
            for mj_int_2 in range(-j_int_2, j_int_2 + 1, 2):
                for mj_int_prime_2 in range(-j_int_2, j_int_2 + 1, 2):

                    probe_amp = R(0)
                    for jpL2, jpR2 in probe_chs:
                        for mpL2 in range(-jpL2, jpL2 + 1, 2):
                            for mpR2 in range(-jpR2, jpR2 + 1, 2):
                                pa = vertex_amp_exact(
                                    jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                    jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                    jpL2, jpR2, mpL2, mpR2)
                                probe_amp += pa

                    if probe_amp == 0:
                        continue

                    q_lo = max(1, abs(n_ext - n_int))
                    q_hi = n_ext + n_int

                    for q_loop in range(q_lo, q_hi + 1):
                        if not vertex_allowed(n_ext, n_int, q_loop):
                            continue
                        mu_q = R(q_loop * (q_loop + 2))

                        chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
                        chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

                        for jgL1_2, jgR1_2 in chs1:
                            for jgL2_2, jgR2_2 in chs2:
                                for mgL2_v in range(-jgL1_2, jgL1_2 + 1, 2):
                                    for mgR2_v in range(-jgR1_2, jgR1_2 + 1, 2):
                                        if abs(mgL2_v) > jgL2_2 or abs(mgR2_v) > jgR2_2:
                                            continue

                                        v1 = vertex_amp_exact(
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                            jgL1_2, jgR1_2, mgL2_v, mgR2_v)
                                        if v1 == 0:
                                            continue

                                        v3 = vertex_amp_exact(
                                            jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jgL2_2, jgR2_2, mgL2_v, mgR2_v)
                                        if v3 == 0:
                                            continue

                                        subtotal += v1 * probe_amp * v3 / (lam4 * mu_q)

        result += sign * subtotal

    return result


def main():
    print("=" * 70)
    print("  EXACT RATIONAL B(n_int) COMPUTATION")
    print("=" * 70)

    # Load float values for comparison
    inv_file = DATA_DIR / "g2_c3_investigation.json"
    with open(inv_file) as f:
        data = json.load(f)
    all_B_float = {int(k): float(v) for k, v in data["all_B"].items()}
    V_mag_float = float(data["V_mag"])

    results = {}

    for n_int in range(1, 4):  # start small: n_int=1,2,3
        t0 = time.time()
        B_exact = compute_B_exact(1, n_int)
        dt = time.time() - t0

        B_simplified = simplify(B_exact)
        B_float_val = float(B_simplified)
        B_stored = all_B_float.get(n_int, 0)
        diff = abs(B_float_val - B_stored)

        is_rational = B_simplified.is_Rational
        if is_rational:
            p = B_simplified.p
            q = B_simplified.q
            results[n_int] = {"p": int(p), "q": int(q), "float": B_float_val, "rational": True}
            print(f"\n  n_int={n_int} ({dt:.1f}s):")
            print(f"    B = {p}/{q}")
            print(f"    q factored = {_factor_small(q)}")
        else:
            results[n_int] = {"expr": str(B_simplified), "float": B_float_val, "rational": False}
            print(f"\n  n_int={n_int} ({dt:.1f}s):")
            print(f"    B = {B_simplified}")

        print(f"    float = {B_float_val:.15e}")
        print(f"    stored = {B_stored:.15e}")
        print(f"    diff   = {diff:.2e}")

    # Now compute the partial sum and extract c3 structure
    print("\n" + "=" * 70)
    print("  PARTIAL RATIONAL SUM AND c3 STRUCTURE")
    print("=" * 70)

    # Cumulative sum
    rational_results = {k: v for k, v in results.items() if v.get("rational", False)}
    if len(rational_results) == len(results):
        cum_B = R(0)
        for n_int in sorted(results.keys()):
            cum_B += R(results[n_int]["p"], results[n_int]["q"])
        n_last = max(results.keys())
        print(f"\n  cum_B(1..{n_last}) = {cum_B.p}/{cum_B.q}")
        print(f"  float       = {float(cum_B):.15e}")
    else:
        print("\n  Not all B values are rational; skipping cumulative sum")

    # Save
    out = {}
    for k, v in results.items():
        if v.get("rational", False):
            out[str(k)] = {"p": v["p"], "q": v["q"], "float": v["float"]}
        else:
            out[str(k)] = {"expr": v.get("expr", "?"), "float": v["float"]}
    outfile = DATA_DIR / "g2_c3_symbolic_B.json"
    with open(outfile, 'w') as f:
        json.dump({"B_exact": out}, f, indent=2)
    print(f"\n  Saved to {outfile}")


def _factor_small(n):
    """Simple factorization for display."""
    if n <= 1:
        return str(n)
    factors = {}
    d = 2
    temp = abs(n)
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factors.items()))


if __name__ == "__main__":
    main()
