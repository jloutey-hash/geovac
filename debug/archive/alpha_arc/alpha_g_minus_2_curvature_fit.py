"""Compute F2/Schwinger at multiple n_ext values to fit curvature expansion.

F2/Schwinger = 1 + c1/lam^2 + c2/lam^4 + ...
where lam = n_ext + 3/2 is the Dirac eigenvalue of the external state.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S
import sys, time, json
import numpy as np
from fractions import Fraction

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


def tree_level_probe_magnetic(n, j, q_probe=1):
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


def compute_ratio(n_ext, n_max_internal):
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    V_mag = tree_level_probe_magnetic(n_ext, j_ext, q_probe=1)
    V_mag_f = float(V_mag)
    if V_mag_f == 0:
        print(f"  V_magnetic = 0 for n_ext={n_ext}, skipping")
        return None

    print(f"  V_magnetic = {V_mag_f:.8f}")

    B_cumul = S.Zero
    for n_int in range(n_max_internal + 1):
        t0 = time.time()
        cup = compute_vertex_3pt_single_nint(n_ext, j_ext, Rational(1, 2), n_int)
        cdn = compute_vertex_3pt_single_nint(n_ext, j_ext, Rational(-1, 2), n_int)
        B_level = cup - cdn
        B_cumul += B_level
        elapsed = time.time() - t0
        bf = float(B_level)
        if bf != 0:
            print(f"  n_int={n_int}: dB={bf:+.4e} ({elapsed:.1f}s)")
        else:
            print(f"  n_int={n_int}: zero ({elapsed:.1f}s)")

    B_f = float(B_cumul)
    F2 = B_f / V_mag_f
    ratio = F2 / schwinger
    return {
        'n_ext': n_ext,
        'lambda_ext': float(Rational(2*n_ext + 3, 2)),
        'V_magnetic': V_mag_f,
        'B': B_f,
        'F2': F2,
        'F2_over_schwinger': ratio,
        'correction': ratio - 1,
        'n_max_internal': n_max_internal,
    }


def main():
    results = []

    # n_ext=1 from prior computation (n_max=12, well converged)
    results.append({
        'n_ext': 1,
        'lambda_ext': 2.5,
        'F2_over_schwinger': 1.08445,
        'correction': 0.08445,
        'n_max_internal': 12,
    })

    # Compute n_ext=2 and n_ext=3 (n_max=5 should suffice given n^-7 convergence)
    for n_ext in [2, 3]:
        n_max = 5
        print(f"\n{'='*60}")
        print(f"n_ext = {n_ext}, n_max = {n_max}")
        print(f"{'='*60}")
        t_total = time.time()
        r = compute_ratio(n_ext, n_max)
        if r:
            print(f"\n  Result: F2/Schwinger = {r['F2_over_schwinger']:.6f}")
            print(f"  correction = {r['correction']:.6f}")
            print(f"  Total time: {time.time() - t_total:.1f}s")
            results.append(r)

    # Fit curvature expansion
    print(f"\n{'='*60}")
    print("Curvature expansion fit")
    print(f"{'='*60}")

    lam_vals = [r['lambda_ext'] for r in results]
    corr_vals = [r['correction'] for r in results]

    print("\nData:")
    for r in results:
        lam = r['lambda_ext']
        print(f"  n_ext={r['n_ext']}, lam={lam}, correction={r['correction']:.6f}, "
              f"1/(2*lam^2)={1/(2*lam**2):.6f}")

    # Fit: correction = c1/lam^2 + c2/lam^4
    n_pts = len(results)
    if n_pts >= 2:
        A = np.array([[1/l**2, 1/l**4] for l in lam_vals])
        b = np.array(corr_vals)
        if n_pts == 2:
            coeffs = np.linalg.solve(A[:2], b[:2])
        else:
            coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

        c1, c2 = coeffs[0], coeffs[1]
        print(f"\n2-term fit: correction = {c1:.6f}/lam^2 + {c2:.6f}/lam^4")
        print(f"  c1 = {c1:.6f}")
        c1_frac = Fraction(c1).limit_denominator(100)
        print(f"  c1 ~ {c1_frac} = {float(c1_frac):.6f}")
        print(f"  Parker-Toms predicts c1 = 1/2 = 0.500000")
        print(f"  c2 = {c2:.6f}")
        c2_frac = Fraction(c2).limit_denominator(1000)
        print(f"  c2 ~ {c2_frac} = {float(c2_frac):.6f}")

        # Check residuals
        for r in results:
            lam = r['lambda_ext']
            pred = c1/lam**2 + c2/lam**4
            print(f"  n_ext={r['n_ext']}: data={r['correction']:.6f}, fit={pred:.6f}, "
                  f"residual={r['correction']-pred:.2e}")

    # If 3 points, also fit 3-term
    if n_pts >= 3:
        A3 = np.array([[1/l**2, 1/l**4, 1/l**6] for l in lam_vals])
        coeffs3 = np.linalg.solve(A3, np.array(corr_vals))
        print(f"\n3-term fit: correction = {coeffs3[0]:.6f}/lam^2 + {coeffs3[1]:.6f}/lam^4 + {coeffs3[2]:.6f}/lam^6")
        for i, c in enumerate(coeffs3):
            frac = Fraction(c).limit_denominator(1000)
            print(f"  c{i+1} = {c:.6f} ~ {frac} = {float(frac):.6f}")

    output = {
        'results': results,
    }
    if n_pts >= 2:
        output['fit_2term'] = {'c1': float(c1), 'c2': float(c2)}
    if n_pts >= 3:
        output['fit_3term'] = {
            'c1': float(coeffs3[0]),
            'c2': float(coeffs3[1]),
            'c3': float(coeffs3[2]),
        }

    with open('debug/data/alpha_g_minus_2_curvature_expansion.json', 'w') as f:
        json.dump(output, f, indent=2)
    print("\nSaved to debug/data/alpha_g_minus_2_curvature_expansion.json")


if __name__ == '__main__':
    main()
