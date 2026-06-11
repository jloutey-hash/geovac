"""
Compute the full vertex correction at n_ext=2,3,4 to separate
curvature expansion coefficients c1, c2, c3, c4.

At higher n_ext, the convergence in n_int is faster (larger lambda),
so fewer n_int levels are needed.

Strategy: compute B(n_ext, n_int) for n_ext=2,3,4, n_int=1..8,
then combine with the n_ext=1 converged data for a 4-point curvature fit.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S
import json, time, sys
import numpy as np

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


def compute_V_magnetic(n_ext):
    """V_magnetic for the external state (n_ext, j=1/2)."""
    j_ext = Rational(1, 2)
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    lam_ext = Rational(2*n_ext + 3, 2)

    probe_chs = get_channels(n_ext, n_ext, 1, jE_L, jE_R, jE_L, jE_R)

    up_amp = S.Zero
    for jpL, jpR in probe_chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                a = vertex_amp_pol(
                    jE_L, jE_R, j_ext, Rational(1, 2),
                    jE_L, jE_R, j_ext, Rational(1, 2),
                    jpL, jpR, mpL, mpR)
                up_amp += a

    dn_amp = S.Zero
    for jpL, jpR in probe_chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                a = vertex_amp_pol(
                    jE_L, jE_R, j_ext, Rational(-1, 2),
                    jE_L, jE_R, j_ext, Rational(-1, 2),
                    jpL, jpR, mpL, mpR)
                dn_amp += a

    V_mag = float(up_amp - dn_amp)
    return V_mag


def main():
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * np.pi)

    # Load n_ext=1 converged data
    with open('debug/data/g2_extended_nint.json') as f:
        ext1_data = json.load(f)

    with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
        ratio_data = json.load(f)

    # n_ext=1: already converged
    delta_1 = ext1_data['corrected_F2_over_S'] - 1.0
    lam_1 = 2.5

    print("=" * 70)
    print("  MULTI-n_ext g-2 CURVATURE EXTRACTION")
    print("=" * 70)
    print(f"\n  n_ext=1 (converged): delta = {delta_1:.10f}, lam = {lam_1}")

    # Compute at n_ext=2,3,4
    results = {1: {'delta': delta_1, 'lam': lam_1, 'lam2': lam_1**2}}

    for n_ext in [2, 3, 4]:
        j_ext = Rational(1, 2)
        lam = (2*n_ext + 3) / 2.0
        lam2 = lam**2

        print(f"\n  Computing n_ext={n_ext} (lam={lam}, lam^2={lam2})...")

        # First compute V_magnetic for this n_ext
        t0 = time.time()
        V_mag = compute_V_magnetic(n_ext)
        print(f"    V_magnetic = {V_mag:.10f} ({time.time()-t0:.1f}s)")

        # Compute B for each n_int
        cum_B = 0.0
        n_int_max = 8  # faster convergence at higher n_ext

        for n_int in range(0, n_int_max + 1):
            t0 = time.time()
            contrib_up = compute_vertex_3pt_single_nint(
                n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
            contrib_dn = compute_vertex_3pt_single_nint(
                n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)

            B_level = float(contrib_up - contrib_dn)
            cum_B += B_level
            F2_over_S = cum_B / V_mag / schwinger

            elapsed = time.time() - t0

            if abs(B_level) > 1e-15 or n_int <= 1:
                print(f"    n_int={n_int:2d}: B={B_level:12.4e}  "
                      f"cum F2/S={F2_over_S:.10f}  ({elapsed:.1f}s)")

            # Early termination if well converged
            if n_int >= 4 and abs(B_level) < abs(cum_B) * 1e-6:
                print(f"    Converged at n_int={n_int} (B_level/cum_B < 1e-6)")
                break

        delta = cum_B / V_mag / schwinger
        results[n_ext] = {
            'delta': delta,
            'lam': lam,
            'lam2': lam2,
            'V_mag': V_mag,
            'cum_B': cum_B,
            'n_int_max': n_int,
        }
        print(f"    RESULT: F2/S - 1 = {delta:.10f}")

    # === Curvature fit ===
    print("\n" + "=" * 70)
    print("  CURVATURE EXPANSION FIT")
    print("=" * 70)

    n_ext_list = sorted(results.keys())
    x_vals = np.array([1.0 / results[n]['lam2'] for n in n_ext_list])
    y_vals = np.array([results[n]['delta'] for n in n_ext_list])

    print(f"\n  Data points:")
    for n in n_ext_list:
        r = results[n]
        print(f"    n_ext={n}: lam={r['lam']:.1f}, 1/lam^2={1/r['lam2']:.6f}, delta={r['delta']:.10f}")

    # 2-parameter fit: delta = c1*x + c2*x^2
    X2 = np.column_stack([x_vals, x_vals**2])
    c_2p, _, _, _ = np.linalg.lstsq(X2, y_vals, rcond=None)
    print(f"\n  2-parameter fit: delta = c1/lam^2 + c2/lam^4")
    print(f"    c1 = {c_2p[0]:.10f}")
    print(f"    c2 = {c_2p[1]:.10f}")

    # 3-parameter fit: delta = c1*x + c2*x^2 + c3*x^3
    X3 = np.column_stack([x_vals, x_vals**2, x_vals**3])
    c_3p, _, _, _ = np.linalg.lstsq(X3, y_vals, rcond=None)
    print(f"\n  3-parameter fit: delta = c1/lam^2 + c2/lam^4 + c3/lam^6")
    print(f"    c1 = {c_3p[0]:.10f}")
    print(f"    c2 = {c_3p[1]:.10f}")
    print(f"    c3 = {c_3p[2]:.10f}")

    # 4-parameter fit
    if len(n_ext_list) >= 4:
        X4 = np.column_stack([x_vals, x_vals**2, x_vals**3, x_vals**4])
        c_4p, _, _, _ = np.linalg.lstsq(X4, y_vals, rcond=None)
        print(f"\n  4-parameter fit: delta = c1/lam^2 + c2/lam^4 + c3/lam^6 + c4/lam^8")
        print(f"    c1 = {c_4p[0]:.10f}")
        print(f"    c2 = {c_4p[1]:.10f}")
        print(f"    c3 = {c_4p[2]:.10f}")
        print(f"    c4 = {c_4p[3]:.10f}")

    # === PSLQ on coefficients ===
    print("\n" + "=" * 70)
    print("  PSLQ IDENTIFICATION OF COEFFICIENTS")
    print("=" * 70)

    try:
        import mpmath
        mpmath.mp.dps = 30

        # Parker-Toms: c1 should be R/12 = 1/2
        c1_best = c_3p[0] if len(n_ext_list) >= 3 else c_2p[0]
        print(f"\n  c1 = {c1_best:.10f}")
        print(f"  c1 - 1/2 = {c1_best - 0.5:.6e}")
        print(f"  Parker-Toms: c1 = R/12 = 1/2 {'CONFIRMED' if abs(c1_best - 0.5) < 0.01 else 'CHECK'}")

        c2_best = c_3p[1] if len(n_ext_list) >= 3 else c_2p[1]
        c2_mp = mpmath.mpf(str(c2_best))

        print(f"\n  c2 = {c2_best:.10f}")

        # PSLQ: try c2 = rational
        rel = mpmath.pslq([c2_mp, mpmath.mpf(1)], tol=1e-5, maxcoeff=1000)
        if rel is not None and rel[0] != 0:
            ratio = -mpmath.mpf(rel[1]) / rel[0]
            print(f"  PSLQ rational: c2 = {ratio} = {float(ratio):.10f}")

        # PSLQ: c2 against curvature basis
        # On S^3: R=6, RicciSq=12, RiemannSq=12
        # Commonly: R^2/180, R_ab^2/180, etc.
        basis = [c2_mp, mpmath.mpf(1), mpmath.pi**2]
        labels = ['c2', '1', 'pi^2']
        rel2 = mpmath.pslq(basis, tol=1e-4, maxcoeff=500)
        if rel2 is not None and rel2[0] != 0:
            terms = [f"{r}*{l}" for r, l in zip(rel2, labels) if r != 0]
            print(f"  PSLQ (rational + pi^2): {' + '.join(terms)} = 0")

        # Try neat fractions near c2
        for denom in range(1, 201):
            numer = round(c2_best * denom)
            if numer > 0:
                approx = numer / denom
                err = abs(approx - c2_best)
                if err < 1e-4:
                    from math import gcd
                    g = gcd(numer, denom)
                    print(f"  Near fraction: {numer//g}/{denom//g} = {approx:.10f} (err {err:.2e})")

    except ImportError:
        print("  mpmath not available")

    # === Save results ===
    output = {
        'multi_next_data': {str(k): v for k, v in results.items()},
        'fit_2p': {'c1': float(c_2p[0]), 'c2': float(c_2p[1])},
    }
    if len(n_ext_list) >= 3:
        output['fit_3p'] = {'c1': float(c_3p[0]), 'c2': float(c_3p[1]), 'c3': float(c_3p[2])}
    if len(n_ext_list) >= 4:
        output['fit_4p'] = {'c1': float(c_4p[0]), 'c2': float(c_4p[1]),
                            'c3': float(c_4p[2]), 'c4': float(c_4p[3])}

    with open('debug/data/g2_multi_next.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to debug/data/g2_multi_next.json")


if __name__ == '__main__':
    main()
