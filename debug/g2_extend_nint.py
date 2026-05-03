"""
Extend the g-2 vertex correction computation from n_int=12 to higher levels.
Goal: improve tail estimate and get 6-7 digit precision on delta.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt as ssqrt
import json, time, sys

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


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    # Load existing data
    with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
        existing = json.load(f)

    V_mag = existing['V_magnetic_float']
    existing_levels = existing['per_level']

    # Current cumulative B
    B_cumul = sum(l['B_level'] for l in existing_levels)
    n_int_done = max(l['n_int'] for l in existing_levels)

    print(f"Existing computation: n_int=0..{n_int_done}")
    print(f"Cumulative B = {B_cumul:.10e}")
    print(f"F2/Schwinger = {B_cumul / V_mag / schwinger:.10f}")
    print()

    # Extend computation
    new_levels = []
    for n_int in range(n_int_done + 1, n_int_done + 4):  # 3 more levels
        print(f"Computing n_int={n_int}...", flush=True)
        t0 = time.time()

        contrib_up = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
        contrib_dn = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)

        B_level = float(contrib_up - contrib_dn)
        B_cumul += B_level
        F2_S = B_cumul / V_mag / schwinger

        elapsed = time.time() - t0

        print(f"  n_int={n_int}: B_level={B_level:.6e}  "
              f"cum F2/S={F2_S:.12f}  ({elapsed:.1f}s)")

        new_levels.append({
            'n_int': n_int,
            'B_level': B_level,
            'F2_over_schwinger': F2_S,
            'time': elapsed,
        })

    # Improved tail estimate
    import numpy as np

    all_B = [l['B_level'] for l in existing_levels if l['n_int'] >= 1] + \
            [l['B_level'] for l in new_levels]
    all_n = [l['n_int'] for l in existing_levels if l['n_int'] >= 1] + \
            [l['n_int'] for l in new_levels]

    # Separate even/odd and fit power laws
    even_n = [n for n in all_n if n % 2 == 0 and n >= 4]
    even_B = [all_B[all_n.index(n)] for n in even_n]
    odd_n = [n for n in all_n if n % 2 == 1 and n >= 3]
    odd_B = [all_B[all_n.index(n)] for n in odd_n]

    log_ne = np.log(np.array(even_n, dtype=float))
    log_Be = np.log(np.abs(np.array(even_B)))
    Ae = np.column_stack([np.ones_like(log_ne), -log_ne])
    fit_e, _, _, _ = np.linalg.lstsq(Ae, log_Be, rcond=None)
    Ce, pe = np.exp(fit_e[0]), fit_e[1]

    log_no = np.log(np.array(odd_n, dtype=float))
    log_Bo = np.log(np.abs(np.array(odd_B)))
    Ao = np.column_stack([np.ones_like(log_no), -log_no])
    fit_o, _, _, _ = np.linalg.lstsq(Ao, log_Bo, rcond=None)
    Co, po = np.exp(fit_o[0]), fit_o[1]

    print(f"\nPower law fits:")
    print(f"  Even: B ~ {Ce:.4e} * n^(-{pe:.2f})")
    print(f"  Odd:  B ~ {Co:.4e} * n^(-{po:.2f})")

    # Compute tail
    n_start = max(all_n) + 1
    tail = 0.0
    for n in range(n_start, 100000):
        if n % 2 == 0:
            tail += Ce * n**(-pe)
        else:
            tail += Co * n**(-po)

    tail_f2s = tail / V_mag / schwinger
    corrected = B_cumul / V_mag / schwinger + tail_f2s
    delta_corrected = corrected - 1.0

    print(f"\nTail estimate (n_int={n_start}..100000):")
    print(f"  Tail B = {tail:.4e}")
    print(f"  Tail F2/S = {tail_f2s:.4e}")
    print(f"  Corrected F2/S = {corrected:.12f}")
    print(f"  Corrected delta = {delta_corrected:.12f}")

    # Save extended results
    output = {
        'extended_levels': new_levels,
        'cumulative_B': B_cumul,
        'cumulative_F2_over_S': B_cumul / V_mag / schwinger,
        'tail_estimate': tail_f2s,
        'corrected_F2_over_S': corrected,
        'corrected_delta': delta_corrected,
        'power_law_even': {'C': float(Ce), 'p': float(pe)},
        'power_law_odd': {'C': float(Co), 'p': float(po)},
    }

    with open('debug/data/g2_extended_nint.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to debug/data/g2_extended_nint.json")


if __name__ == '__main__':
    main()
