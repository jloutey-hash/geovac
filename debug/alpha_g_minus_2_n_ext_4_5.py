"""Compute F2/Schwinger at n_ext=4 and n_ext=5 on unit S^3.

Extends the 3-point dataset to 5 points for power-law exponent determination.
Uses exact sympy CG coefficients, same computation as alpha_g_minus_2_curvature_fit.py.

Existing data:
  n_ext=1 (lam=2.5): F2/Schwinger = 1.08445  (n_max=12)
  n_ext=2 (lam=3.5): F2/Schwinger = 0.6133   (n_max=5)
  n_ext=3 (lam=4.5): F2/Schwinger = 0.3702   (n_max=5)
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S
import sys, time, json
import numpy as np
from fractions import Fraction
from scipy.optimize import curve_fit

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
        B_f = float(B_cumul)
        F2_running = B_f / V_mag_f
        ratio_running = F2_running / schwinger
        if bf != 0:
            print(f"  n_int={n_int}: dB={bf:+.4e}, B_cumul={B_f:.6e}, "
                  f"F2/Schwinger={ratio_running:.6f} ({elapsed:.1f}s)")
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
    # Prior results
    prior = [
        {'n_ext': 1, 'lambda_ext': 2.5, 'F2_over_schwinger': 1.08445, 'n_max_internal': 12},
        {'n_ext': 2, 'lambda_ext': 3.5, 'F2_over_schwinger': 0.6133, 'n_max_internal': 5},
        {'n_ext': 3, 'lambda_ext': 4.5, 'F2_over_schwinger': 0.3702, 'n_max_internal': 5},
    ]

    new_results = []

    # Compute n_ext=4 and n_ext=5
    for n_ext in [4, 5]:
        n_max = 5
        print(f"\n{'='*60}")
        print(f"n_ext = {n_ext}, n_max = {n_max}")
        print(f"lambda = {n_ext + 1.5}")
        print(f"{'='*60}")
        t_total = time.time()
        r = compute_ratio(n_ext, n_max)
        if r:
            print(f"\n  Result: F2/Schwinger = {r['F2_over_schwinger']:.6f}")
            print(f"  Total time: {time.time() - t_total:.1f}s")
            new_results.append(r)

    # Combine all 5 data points
    all_results = []
    for p in prior:
        all_results.append(p)
    for r in new_results:
        all_results.append({
            'n_ext': r['n_ext'],
            'lambda_ext': r['lambda_ext'],
            'F2_over_schwinger': r['F2_over_schwinger'],
            'n_max_internal': r['n_max_internal'],
            'V_magnetic': r.get('V_magnetic'),
            'B': r.get('B'),
            'F2': r.get('F2'),
        })

    print(f"\n{'='*60}")
    print("5-point power-law fit")
    print(f"{'='*60}")

    lam = np.array([r['lambda_ext'] for r in all_results])
    ratio = np.array([r['F2_over_schwinger'] for r in all_results])

    print("\nData:")
    for r in all_results:
        print(f"  n_ext={r['n_ext']}, lam={r['lambda_ext']:.1f}, "
              f"F2/Schwinger={r['F2_over_schwinger']:.6f}")

    # --- Power-law fit: F2/Schwinger = A * lam^p ---
    log_lam = np.log(lam)
    log_ratio = np.log(ratio)

    # Linear fit in log-log space
    coeffs = np.polyfit(log_lam, log_ratio, 1)
    p_fit = coeffs[0]
    A_fit = np.exp(coeffs[1])
    print(f"\nLog-log linear fit: F2/Schwinger = {A_fit:.4f} * lam^({p_fit:.4f})")

    # Nonlinear fit with uncertainties
    def power_law(x, A, p):
        return A * x**p
    popt, pcov = curve_fit(power_law, lam, ratio, p0=[A_fit, p_fit])
    A_nl, p_nl = popt
    A_err, p_err = np.sqrt(np.diag(pcov))
    print(f"\nNonlinear fit: F2/Schwinger = {A_nl:.6f} * lam^({p_nl:.6f})")
    print(f"  A = {A_nl:.6f} +/- {A_err:.6f}")
    print(f"  p = {p_nl:.6f} +/- {p_err:.6f}")

    # Residuals for the free-exponent fit
    pred_free = A_nl * lam**p_nl
    resid_free = ratio - pred_free
    rms_free = np.sqrt(np.mean(resid_free**2))
    print(f"\n  RMS residual (free p): {rms_free:.2e}")
    for i, r in enumerate(all_results):
        print(f"    n_ext={r['n_ext']}: data={ratio[i]:.6f}, "
              f"pred={pred_free[i]:.6f}, resid={resid_free[i]:.2e}")

    # --- Forced p = -2 fit: F2/Schwinger = A / lam^2 ---
    # Least squares: min sum (ratio_i - A/lam_i^2)^2
    inv_lam2 = 1.0 / lam**2
    A_forced = np.sum(ratio * inv_lam2) / np.sum(inv_lam2**2)
    pred_forced = A_forced * inv_lam2
    resid_forced = ratio - pred_forced
    rms_forced = np.sqrt(np.mean(resid_forced**2))
    print(f"\nForced p=-2 fit: F2/Schwinger = {A_forced:.6f} / lam^2")
    print(f"  RMS residual (p=-2): {rms_forced:.2e}")
    for i, r in enumerate(all_results):
        print(f"    n_ext={r['n_ext']}: data={ratio[i]:.6f}, "
              f"pred={pred_forced[i]:.6f}, resid={resid_forced[i]:.2e}")

    # Check if A_forced is close to a clean number
    A_frac = Fraction(A_forced).limit_denominator(100)
    print(f"\n  A (p=-2) = {A_forced:.6f}")
    print(f"  Nearest simple fraction: {A_frac} = {float(A_frac):.6f}")
    # Check some specific candidates
    candidates = [
        ('6', 6.0),
        ('5', 5.0),
        ('4', 4.0),
        ('15/2', 7.5),
        ('11/2', 5.5),
        ('7', 7.0),
        ('9/2', 4.5),
        ('3*pi/2', 3*3.14159265/2),
        ('2*pi', 2*3.14159265),
        ('sqrt(35)', 35**0.5),
        ('sqrt(30)', 30**0.5),
        ('sqrt(40)', 40**0.5),
    ]
    print("\n  Candidate A values (for p=-2):")
    for name, val in candidates:
        print(f"    A = {name} = {val:.6f}, ratio to fit = {val/A_forced:.6f}")

    # --- Curvature expansion fit: F2/Schwinger = c0 + c1/lam^2 + c2/lam^4 ---
    # where c0 should be 0 in the flat-space limit
    print(f"\n{'='*60}")
    print("Curvature expansion (correction = F2/Schwinger - 1)")
    print(f"{'='*60}")

    corrections = ratio - 1.0
    # 2-term fit: correction = c1/lam^2 + c2/lam^4
    A2 = np.column_stack([inv_lam2, 1.0/lam**4])
    coeffs2, res2, _, _ = np.linalg.lstsq(A2, corrections, rcond=None)
    c1_2, c2_2 = coeffs2
    pred2 = A2 @ coeffs2
    rms2 = np.sqrt(np.mean((corrections - pred2)**2))
    print(f"\n2-term: correction = {c1_2:.6f}/lam^2 + {c2_2:.6f}/lam^4")
    print(f"  c1 = {c1_2:.6f} (Parker-Toms: 1/2 = 0.500000)")
    print(f"  c2 = {c2_2:.6f}")
    print(f"  RMS = {rms2:.2e}")

    # 3-term fit: correction = c1/lam^2 + c2/lam^4 + c3/lam^6
    A3 = np.column_stack([inv_lam2, 1.0/lam**4, 1.0/lam**6])
    coeffs3, res3, _, _ = np.linalg.lstsq(A3, corrections, rcond=None)
    c1_3, c2_3, c3_3 = coeffs3
    pred3 = A3 @ coeffs3
    rms3 = np.sqrt(np.mean((corrections - pred3)**2))
    print(f"\n3-term: correction = {c1_3:.6f}/lam^2 + {c2_3:.6f}/lam^4 + {c3_3:.6f}/lam^6")
    print(f"  c1 = {c1_3:.6f}")
    print(f"  c2 = {c2_3:.6f}")
    print(f"  c3 = {c3_3:.6f}")
    print(f"  RMS = {rms3:.2e}")
    for i, r in enumerate(all_results):
        print(f"    n_ext={r['n_ext']}: data={corrections[i]:.6f}, "
              f"pred={pred3[i]:.6f}, resid={corrections[i]-pred3[i]:.2e}")

    # Is p=-2 consistent? Compare sigma_p to |p+2|
    sigma_consistent = abs(p_nl + 2) / p_err if p_err > 0 else float('inf')
    print(f"\n{'='*60}")
    print("Consistency check: is p = -2 exact?")
    print(f"{'='*60}")
    print(f"  p = {p_nl:.6f} +/- {p_err:.6f}")
    print(f"  |p - (-2)| = {abs(p_nl + 2):.6f}")
    print(f"  |p - (-2)| / sigma_p = {sigma_consistent:.2f} sigma")
    if sigma_consistent < 2:
        print(f"  => p = -2 is CONSISTENT (within 2 sigma)")
    elif sigma_consistent < 3:
        print(f"  => p = -2 is MARGINAL (2-3 sigma)")
    else:
        print(f"  => p = -2 is REJECTED (> 3 sigma)")

    # Save results
    output = {
        'data_points': [],
        'power_law_fit': {
            'A': float(A_nl),
            'A_err': float(A_err),
            'p': float(p_nl),
            'p_err': float(p_err),
            'rms_residual': float(rms_free),
        },
        'forced_p_minus_2': {
            'A': float(A_forced),
            'rms_residual': float(rms_forced),
            'p_minus_2_sigma': float(sigma_consistent),
        },
        'curvature_expansion_2term': {
            'c1': float(c1_2),
            'c2': float(c2_2),
            'rms': float(rms2),
        },
        'curvature_expansion_3term': {
            'c1': float(c1_3),
            'c2': float(c2_3),
            'c3': float(c3_3),
            'rms': float(rms3),
        },
    }
    for r in all_results:
        dp = {
            'n_ext': r['n_ext'],
            'lambda_ext': r['lambda_ext'],
            'F2_over_schwinger': r['F2_over_schwinger'],
            'n_max_internal': r['n_max_internal'],
        }
        if 'V_magnetic' in r and r['V_magnetic'] is not None:
            dp['V_magnetic'] = r['V_magnetic']
            dp['B'] = r['B']
            dp['F2'] = r['F2']
        output['data_points'].append(dp)

    with open('debug/data/alpha_g_minus_2_5point_fit.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to debug/data/alpha_g_minus_2_5point_fit.json")


if __name__ == '__main__':
    main()
