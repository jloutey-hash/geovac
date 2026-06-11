"""
g2_c3_vmag_structure.py -- Check if B(n_int)/V_mag simplifies to rational.

V_mag(n_ext) = (2/3)(sqrt((n+3)/(n+1)) - sqrt(n/(n+2)))

For n_ext=1: V_mag = (2/3)(sqrt(2) - 1/sqrt(3))

B(n_int) lives in Q(sqrt(2), sqrt(3), sqrt(5), sqrt(7), ...)

Key question: does F2/Schwinger = cum_B / (V_mag * Schwinger) simplify
to Q[pi^2] as T9 theorem demands?
"""

import time
from sympy.physics.wigner import clebsch_gordan
from sympy import (Rational as R, sqrt, simplify, radsimp, nsimplify,
                   sqrtdenest, pi, factor, cancel, collect, Symbol)
from functools import lru_cache


def V_magnetic_closed(n):
    """Exact symbolic V_mag for n_ext = n."""
    return R(2, 3) * (sqrt(R(n + 3, n + 1)) - sqrt(R(n, n + 2)))


@lru_cache(maxsize=2_000_000)
def cg_exact(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
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
    result = R(0)
    q_probe = 1

    for mj_ext_2 in [+1, -1]:
        sign = 1 if mj_ext_2 == +1 else -1

        jE_L_2 = n_ext + 1
        jE_R_2 = n_ext
        j_ext_2 = 1

        jI_L_2 = n_int + 1
        jI_R_2 = n_int

        lam_num = 2 * n_int + 3
        lam4 = R(lam_num**4, 16)

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
    n_ext = 1
    V_mag = V_magnetic_closed(n_ext)
    V_mag_simp = radsimp(V_mag)

    print("=" * 70)
    print("  V_mag AND B/V_mag RATIONALITY TEST")
    print("=" * 70)

    print(f"\n  V_mag(n_ext=1) = {V_mag}")
    print(f"  V_mag simplified = {V_mag_simp}")
    print(f"  V_mag float = {float(V_mag_simp.evalf(30)):.15e}")

    # Check: V_mag^2 should simplify nicely
    V2 = radsimp(V_mag**2)
    print(f"\n  V_mag^2 = {V2}")
    print(f"  V_mag^2 float = {float(V2.evalf(30)):.15e}")

    # Compute B(n_int)/V_mag for n_int=1..3
    print("\n" + "=" * 70)
    print("  B(n_int) / V_mag RATIONALITY CHECK")
    print("=" * 70)

    SCHWINGER = R(1)  # We'll divide by Schwinger=alpha/(2pi) separately

    cum_B = R(0)
    cum_ratio = R(0)

    for n_int in range(1, 4):
        t0 = time.time()
        B = compute_B_exact(n_ext, n_int)
        dt = time.time() - t0
        B_simp = radsimp(B)

        cum_B += B_simp

        # Key test: B / V_mag
        ratio = radsimp(B_simp / V_mag_simp)
        cum_ratio += ratio
        cum_ratio_simp = radsimp(cum_ratio)

        is_rat = ratio.is_Rational
        cum_is_rat = cum_ratio_simp.is_Rational

        print(f"\n  n_int={n_int} ({dt:.1f}s):")
        print(f"    B = {B_simp}")
        print(f"    B float = {float(B_simp.evalf(30)):.15e}")
        print(f"    B/V_mag = {ratio}")
        print(f"    B/V_mag rational? {is_rat}")
        print(f"    cum(B/V_mag) = {cum_ratio_simp}")
        print(f"    cum rational? {cum_is_rat}")

    # Also try B * V_mag (might cancel differently)
    print("\n" + "=" * 70)
    print("  ALTERNATIVE: cum_B * V_mag")
    print("=" * 70)
    product = radsimp(cum_B * V_mag_simp)
    print(f"  cum_B * V_mag = {product}")
    print(f"  rational? {product.is_Rational}")


if __name__ == "__main__":
    main()
