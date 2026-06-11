"""
Extended convergence of the normalized anomalous magnetic moment on S^3.

Tracks F2 / [alpha/(2*pi)] vs n_max to determine whether the ratio
approaches 1 (flat-space Schwinger) or a curvature-corrected value.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S
import sys, json, time

sys.path.insert(0, '.')


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL:
            continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR:
            continue
        if mL2 + mR2 != mj_t:
            continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0:
            continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0:
            continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, j_sL, j_sR, j_tL, j_tR):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL, jgR in [(Rational(q+1,2), Rational(q-1,2)),
                      (Rational(q-1,2), Rational(q+1,2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j_sL - jgL) <= j_tL <= j_sL + jgL and
            abs(j_sR - jgR) <= j_tR <= j_sR + jgR):
            chs.append((jgL, jgR))
    return chs


def compute_vertex_3pt_single_nint(n_ext, j_ext, mj_ext, n_int, q_probe=1):
    """Contribution from a single n_int level."""
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2*n_int + 3, 2)
    lam4 = lam**4

    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_chs = get_channels(n_int, n_int, q_probe,
                              jI_L, jI_R, jI_L, jI_R)
    if not probe_chs:
        return S.Zero

    total = S.Zero

    for j_int in half_ints(j_int_max):
        if j_int < j_int_min:
            continue
        for mj_int in half_ints(j_int):
            for mj_int_prime in half_ints(j_int):
                probe_amp = S.Zero
                for jpL, jpR in probe_chs:
                    for mpL in half_ints(jpL):
                        for mpR in half_ints(jpR):
                            pa = vertex_amp_pol(
                                jI_L, jI_R, j_int, mj_int,
                                jI_L, jI_R, j_int, mj_int_prime,
                                jpL, jpR, mpL, mpR)
                            probe_amp += pa
                if probe_amp == 0:
                    continue

                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int
                for q_loop in range(q_lo, q_hi + 1):
                    if not vertex_allowed(n_ext, n_int, q_loop):
                        continue
                    mu_q = Rational(q_loop * (q_loop + 2))
                    chs1 = get_channels(n_ext, n_int, q_loop,
                                        jE_L, jE_R, jI_L, jI_R)
                    chs2 = get_channels(n_int, n_ext, q_loop,
                                        jI_L, jI_R, jE_L, jE_R)
                    for jgL1, jgR1 in chs1:
                        for jgL2, jgR2 in chs2:
                            for mgL in half_ints(jgL1):
                                for mgR in half_ints(jgR1):
                                    if abs(mgL) > jgL2 or abs(mgR) > jgR2:
                                        continue
                                    v1 = vertex_amp_pol(
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jI_L, jI_R, j_int, mj_int,
                                        jgL1, jgR1, mgL, mgR)
                                    if v1 == 0:
                                        continue
                                    v3 = vertex_amp_pol(
                                        jI_L, jI_R, j_int, mj_int_prime,
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jgL2, jgR2, mgL, mgR)
                                    if v3 == 0:
                                        continue
                                    total += v1 * probe_amp * v3 / (lam4 * mu_q)

    return total


def tree_level_probe_magnetic(n, j, q_probe=1):
    """V_tree(+1/2) - V_tree(-1/2) for the tree-level probe coupling."""
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    if not vertex_allowed(n, n, q_probe):
        return S.Zero
    chs = get_channels(n, n, q_probe, jL, jR, jL, jR)
    if not chs:
        return S.Zero
    total_up = S.Zero
    total_dn = S.Zero
    for jpL, jpR in chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                total_up += vertex_amp_pol(jL, jR, j, Rational(1, 2),
                                           jL, jR, j, Rational(1, 2),
                                           jpL, jpR, mpL, mpR)
                total_dn += vertex_amp_pol(jL, jR, j, Rational(-1, 2),
                                           jL, jR, j, Rational(-1, 2),
                                           jpL, jpR, mpL, mpR)
    return total_up - total_dn


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    V_mag = tree_level_probe_magnetic(n_ext, j_ext, q_probe=1)
    V_mag_float = float(V_mag)

    print("=" * 70)
    print("Extended convergence: F2 / [alpha/(2*pi)] on S^3")
    print(f"External: n_CH={n_ext}, j={j_ext}")
    print(f"V_tree_magnetic = {V_mag} = {V_mag_float:.8f}")
    print(f"Schwinger = {schwinger:.6e}")
    print("=" * 70)

    results = []
    B_cumul = S.Zero

    n_max_target = 8

    for n_int in range(n_max_target + 1):
        t0 = time.time()

        contrib_up = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
        contrib_dn = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)

        B_level = contrib_up - contrib_dn
        B_cumul += B_level

        B_level_f = float(B_level)
        B_cumul_f = float(B_cumul)
        F2 = B_cumul_f / V_mag_float if V_mag_float != 0 else 0
        ratio = F2 / schwinger if schwinger != 0 else 0

        elapsed = time.time() - t0

        if B_level_f != 0:
            print(f"  n={n_int}: dB={B_level_f:+.4e}, B={B_cumul_f:+.6e}, "
                  f"F2/Schwinger={ratio:.6f} ({elapsed:.1f}s)")
        else:
            print(f"  n={n_int}: zero ({elapsed:.1f}s)")

        results.append({
            'n_int': n_int,
            'B_level': B_level_f,
            'B_cumul': B_cumul_f,
            'F2': F2,
            'F2_over_schwinger': ratio,
            'time': elapsed,
        })

    # Check if converging to some rational multiple of schwinger
    final_ratio = results[-1]['F2_over_schwinger']
    print(f"\n  Final F2/Schwinger = {final_ratio:.6f}")
    print(f"  Possible values near {final_ratio:.3f}:")
    for num in range(1, 20):
        for den in range(1, 20):
            frac = num / den
            if abs(frac - final_ratio) < 0.02:
                print(f"    {num}/{den} = {frac:.6f} (diff = {abs(frac-final_ratio):.4e})")

    output = {
        'n_ext': n_ext,
        'V_magnetic': V_mag_float,
        'schwinger': schwinger,
        'per_level': results,
        'final_ratio': final_ratio,
    }

    with open('debug/data/alpha_g_minus_2_extended.json', 'w') as f:
        json.dump(output, f, indent=2)
    print("\nSaved to debug/data/alpha_g_minus_2_extended.json")


if __name__ == '__main__':
    main()
