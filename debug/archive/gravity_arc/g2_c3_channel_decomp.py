"""
g2_c3_channel_decomp.py -- Decompose B(n_int) by (j_int, channel_pair).

Key structural simplification for n_ext=1:
  - Loop photon has exactly q = n_int for each n_int >= 1
  - Only 2 channels: (jgL, jgR) = (n_int+1, n_int-1) and (n_int-1, n_int+1)
  - j_int ranges from 1/2 to n_int + 1/2

Goal: check if irrationals cancel at the channel-sum level (summing over
channels at fixed n_int and j_int), which would enable closed-form B(n_int).
"""

import time
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational as R, sqrt, simplify, radsimp, S
from functools import lru_cache


@lru_cache(maxsize=2_000_000)
def cg_exact(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
    if m1_2 + m2_2 != M_2:
        return R(0)
    if abs(j1_2 - j2_2) > J_2 or J_2 > j1_2 + j2_2:
        return R(0)
    if abs(m1_2) > j1_2 or abs(m2_2) > j2_2 or abs(M_2) > J_2:
        return R(0)
    return clebsch_gordan(R(j1_2, 2), R(j2_2, 2), R(J_2, 2),
                          R(m1_2, 2), R(m2_2, 2), R(M_2, 2))


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


def decompose_B(n_ext, n_int):
    """
    Decompose B(n_int) by (j_int, channel_pair_v1, channel_pair_v3).

    Returns dict: (j_int_2, ch1_idx, ch3_idx) -> contribution to B
    """
    jE_L_2 = n_ext + 1  # 2
    jE_R_2 = n_ext       # 1
    j_ext_2 = 1
    jI_L_2 = n_int + 1
    jI_R_2 = n_int

    lam_num = 2 * n_int + 3
    lam4 = R(lam_num**4, 16)

    q_probe = 1
    if not vertex_allowed(n_int, n_int, q_probe):
        return {}
    probe_chs = get_channels(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
    if not probe_chs:
        return {}

    q_loop = n_int  # single q value for n_ext=1
    if not vertex_allowed(n_ext, n_int, q_loop):
        return {}
    mu_q = R(q_loop * (q_loop + 2))

    chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
    chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

    j_int_min_2 = 1
    j_int_max_2 = 2 * n_int + 1

    contributions = {}

    for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
        for ch1_idx, (jgL1_2, jgR1_2) in enumerate(chs1):
            for ch3_idx, (jgL2_2, jgR2_2) in enumerate(chs2):
                key = (j_int_2, ch1_idx, ch3_idx)
                subtotal = R(0)

                for mj_ext_2 in [+1, -1]:
                    sign = 1 if mj_ext_2 == +1 else -1

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

                            v_amp = R(0)
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
                                    v_amp += v1 * v3

                            subtotal += sign * v_amp * probe_amp / (lam4 * mu_q)

                if subtotal != 0:
                    contributions[key] = radsimp(subtotal)

    return contributions


def main():
    print("=" * 70)
    print("  B(n_int) CHANNEL DECOMPOSITION")
    print("=" * 70)

    n_ext = 1

    for n_int in range(1, 6):
        t0 = time.time()
        contribs = decompose_B(n_ext, n_int)
        dt = time.time() - t0

        q_loop = n_int
        lam = R(2 * n_int + 3, 2)
        mu_q = R(q_loop * (q_loop + 2))

        jI_L_2 = n_int + 1
        jI_R_2 = n_int
        jE_L_2 = n_ext + 1
        jE_R_2 = n_ext

        chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
        chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

        print(f"\n  n_int={n_int} (lam={lam}, q={q_loop}, mu_q={mu_q}, {dt:.1f}s):")
        print(f"    v1 channels: {chs1}")
        print(f"    v3 channels: {chs2}")

        B_total = R(0)

        # Group by j_int
        j_int_vals = sorted(set(k[0] for k in contribs.keys()))
        for j_int_2 in j_int_vals:
            j_sub = {k: v for k, v in contribs.items() if k[0] == j_int_2}
            j_sum = sum(j_sub.values())
            j_sum_simp = radsimp(j_sum)
            is_rat = j_sum_simp.is_Rational
            B_total += j_sum_simp

            print(f"\n    j_int={j_int_2}/2:")

            for (_, ch1, ch3), val in sorted(j_sub.items()):
                ch1_label = f"({chs1[ch1][0]},{chs1[ch1][1]})" if ch1 < len(chs1) else "?"
                ch3_label = f"({chs2[ch3][0]},{chs2[ch3][1]})" if ch3 < len(chs2) else "?"
                print(f"      ch1={ch1_label} ch3={ch3_label}: {val}")
                print(f"        rational? {val.is_Rational}  float={float(val.evalf(15)):.10e}")

            print(f"      SUM over channels: {j_sum_simp}")
            print(f"        rational? {is_rat}  float={float(j_sum_simp.evalf(15)):.10e}")

        B_total_simp = radsimp(B_total)
        print(f"\n    B(n_int={n_int}) total = {B_total_simp}")
        print(f"      rational? {B_total_simp.is_Rational}")
        print(f"      float = {float(B_total_simp.evalf(15)):.10e}")

        # Check B² rationality
        B2 = radsimp(B_total_simp**2)
        print(f"      B² = {B2}")
        print(f"      B² rational? {B2.is_Rational}")


if __name__ == "__main__":
    main()
