"""
Multi-n_ext converged vertex correction for c2 curvature extraction.

Computes F2/Schwinger at n_ext = 1, 2, 3, 4, 5 with enough n_int levels
to converge each one, then fits the curvature expansion:

    delta(n_ext) = F2/S - 1 = c1/lam^2 + c2/lam^4 + c3/lam^6 + ...

where lam = n_ext + 3/2 (Dirac eigenvalue on unit S3).

Key insight: higher n_ext has larger lam, so:
  (a) the curvature expansion converges faster (less c3 contamination)
  (b) the n_int sum converges faster (power-law tail ~n_int^{-7})

n_ext=1 data is loaded from the existing converged run (n_int=0..25).
n_ext=2,3,4,5 are computed fresh with adaptive n_int stopping.

Expected runtime: ~6-10 hours total.

Target: c2 at 10+ digit precision (currently 7 digits from n_ext=1 alone).
This resolves whether c2 = 7/40 or c2 = 19/100 - 41*pi^2/25200.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S
import json
import time
import sys
import os
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
    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                      (Rational(q-1, 2), Rational(q+1, 2))]:
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


def power_law_tail(all_B, n_start):
    """Estimate tail sum from power-law fit to B(n) for n >= 3."""
    ns = sorted(n for n in all_B.keys() if n >= 3)
    even_n = [n for n in ns if n % 2 == 0 and n >= 4]
    even_B = [all_B[n] for n in even_n]
    odd_n = [n for n in ns if n % 2 == 1 and n >= 3]
    odd_B = [all_B[n] for n in odd_n]

    if len(even_n) < 2 or len(odd_n) < 2:
        return 0.0, {}, {}

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

    tail = 0.0
    for n in range(n_start, 100001):
        if n % 2 == 0:
            tail += Ce * n**(-pe)
        else:
            tail += Co * n**(-po)

    return tail, {'C': float(Ce), 'p': float(pe)}, {'C': float(Co), 'p': float(po)}


def save_checkpoint(results, filename='debug/data/g2_multi_next_converged.json'):
    """Save all results after each n_ext completes."""
    output = {}
    for n_ext, r in sorted(results.items()):
        output[str(n_ext)] = {
            'delta': r['delta'],
            'delta_corrected': r.get('delta_corrected', r['delta']),
            'lam': r['lam'],
            'lam2': r['lam2'],
            'V_mag': r.get('V_mag', None),
            'cum_B': r.get('cum_B', None),
            'n_int_max': r.get('n_int_max', None),
            'all_B': r.get('all_B', None),
            'tail_B': r.get('tail_B', None),
            'power_law_even': r.get('power_law_even', None),
            'power_law_odd': r.get('power_law_odd', None),
        }

    # Add polynomial fits if we have >= 2 points
    n_ext_list = sorted(results.keys())
    if len(n_ext_list) >= 2:
        x_vals = np.array([1.0 / results[n]['lam2'] for n in n_ext_list])
        y_vals = np.array([results[n].get('delta_corrected', results[n]['delta'])
                           for n in n_ext_list])

        # With c1 = 1/2 fixed (exact Parker-Toms), fit residual
        y_residual = y_vals - 0.5 * x_vals
        x2 = x_vals**2

        if len(n_ext_list) >= 2:
            # 1-param fit: residual = c2 * x^2
            X1 = x2.reshape(-1, 1)
            c_1p, _, _, _ = np.linalg.lstsq(X1, y_residual, rcond=None)
            output['fit_c1_fixed_1p'] = {
                'c1': 0.5, 'c2': float(c_1p[0]),
                'note': 'c1=1/2 fixed, 1-param fit for c2'
            }

        if len(n_ext_list) >= 3:
            # 2-param fit: residual = c2 * x^2 + c3 * x^3
            X2 = np.column_stack([x2, x_vals**3])
            c_2p, _, _, _ = np.linalg.lstsq(X2, y_residual, rcond=None)
            output['fit_c1_fixed_2p'] = {
                'c1': 0.5, 'c2': float(c_2p[0]), 'c3': float(c_2p[1]),
                'note': 'c1=1/2 fixed, 2-param fit for c2, c3'
            }

        if len(n_ext_list) >= 4:
            # 3-param fit: residual = c2*x^2 + c3*x^3 + c4*x^4
            X3 = np.column_stack([x2, x_vals**3, x_vals**4])
            c_3p, _, _, _ = np.linalg.lstsq(X3, y_residual, rcond=None)
            output['fit_c1_fixed_3p'] = {
                'c1': 0.5, 'c2': float(c_3p[0]), 'c3': float(c_3p[1]),
                'c4': float(c_3p[2]),
                'note': 'c1=1/2 fixed, 3-param fit for c2, c3, c4'
            }

        if len(n_ext_list) >= 5:
            # 4-param fit
            X4 = np.column_stack([x2, x_vals**3, x_vals**4, x_vals**5])
            c_4p, _, _, _ = np.linalg.lstsq(X4, y_residual, rcond=None)
            output['fit_c1_fixed_4p'] = {
                'c1': 0.5, 'c2': float(c_4p[0]), 'c3': float(c_4p[1]),
                'c4': float(c_4p[2]), 'c5': float(c_4p[3]),
                'note': 'c1=1/2 fixed, 4-param fit'
            }

        # Also free fits (c1 not fixed)
        if len(n_ext_list) >= 2:
            X_free = np.column_stack([x_vals, x2])
            c_free, _, _, _ = np.linalg.lstsq(X_free, y_vals, rcond=None)
            output['fit_free_2p'] = {'c1': float(c_free[0]), 'c2': float(c_free[1])}

        if len(n_ext_list) >= 3:
            X_free3 = np.column_stack([x_vals, x2, x_vals**3])
            c_free3, _, _, _ = np.linalg.lstsq(X_free3, y_vals, rcond=None)
            output['fit_free_3p'] = {
                'c1': float(c_free3[0]), 'c2': float(c_free3[1]),
                'c3': float(c_free3[2])
            }

        if len(n_ext_list) >= 4:
            X_free4 = np.column_stack([x_vals, x2, x_vals**3, x_vals**4])
            c_free4, _, _, _ = np.linalg.lstsq(X_free4, y_vals, rcond=None)
            output['fit_free_4p'] = {
                'c1': float(c_free4[0]), 'c2': float(c_free4[1]),
                'c3': float(c_free4[2]), 'c4': float(c_free4[3])
            }

    with open(filename, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"  [checkpoint saved to {filename}]", flush=True)


def main():
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * np.pi)

    print("=" * 70)
    print("  MULTI-n_ext CONVERGED g-2 CURVATURE EXTRACTION")
    print("  Target: c2 at 10+ digit precision")
    print("=" * 70)
    print(flush=True)

    results = {}

    # === n_ext=1: load from existing converged data (n_int=0..25) ===
    try:
        with open('debug/data/g2_extended_nint_v2.json') as f:
            v2_data = json.load(f)

        with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
            ratio_data = json.load(f)

        V_mag_1 = ratio_data['V_magnetic_float']
        all_B_1 = {int(k): v for k, v in v2_data['all_B'].items()}
        cum_B_1 = sum(all_B_1.values())

        # Use tail-corrected value
        delta_1 = v2_data['corrected_delta']
        delta_raw_1 = cum_B_1 / V_mag_1 / schwinger - 1.0

        results[1] = {
            'delta': delta_raw_1,
            'delta_corrected': delta_1,
            'lam': 2.5,
            'lam2': 6.25,
            'V_mag': V_mag_1,
            'cum_B': cum_B_1,
            'n_int_max': 25,
            'all_B': {str(k): v for k, v in all_B_1.items()},
            'tail_B': v2_data['tail_estimate'] * V_mag_1 * schwinger,
            'power_law_even': v2_data['power_law_even'],
            'power_law_odd': v2_data['power_law_odd'],
        }

        print(f"  n_ext=1 loaded (n_int=0..25, tail-corrected):")
        print(f"    V_mag = {V_mag_1:.12e}")
        print(f"    cum_B = {cum_B_1:.12e}")
        print(f"    delta_raw = {delta_raw_1:.12f}")
        print(f"    delta_corrected = {delta_1:.12f}")
        print(f"    lam = 2.5, lam^2 = 6.25")
        print(flush=True)

    except FileNotFoundError:
        print("  ERROR: n_ext=1 data not found. Run g2_extend_nint_v2.py first.")
        sys.exit(1)

    # === n_ext=2,3,4,5: compute fresh ===
    # Convergence strategy: at higher n_ext, the B(n_int) tail falls faster
    # because the vertex selection rule restricts q to [|n_ext-n_int|, n_ext+n_int]
    # and the Clebsch-Gordan coefficients suppress large angular momentum transfer.
    #
    # Adaptive stopping: run until |B_level| < |cum_B| * 1e-10 for 3 consecutive levels,
    # OR until n_int reaches a maximum.
    #
    # Conservative n_int_max per n_ext (based on convergence rate):
    #   n_ext=2: n_int_max=20 (lam=3.5, ~6x faster convergence than n_ext=1)
    #   n_ext=3: n_int_max=18 (lam=4.5)
    #   n_ext=4: n_int_max=15 (lam=5.5)
    #   n_ext=5: n_int_max=15 (lam=6.5)

    n_ext_configs = [
        (2, 20),
        (3, 18),
        (4, 15),
        (5, 15),
    ]

    for n_ext, n_int_max_target in n_ext_configs:
        j_ext = Rational(1, 2)
        lam = (2*n_ext + 3) / 2.0
        lam2 = lam**2

        print(f"\n{'='*70}")
        print(f"  n_ext={n_ext}, lam={lam}, lam^2={lam2}")
        print(f"  n_int_max target = {n_int_max_target}")
        print(f"{'='*70}", flush=True)

        # V_magnetic for this n_ext
        t0 = time.time()
        V_mag = compute_V_magnetic(n_ext)
        print(f"  V_magnetic = {V_mag:.12e}  ({time.time()-t0:.1f}s)", flush=True)

        # Compute B per n_int
        all_B = {}
        cum_B = 0.0
        consecutive_small = 0

        for n_int in range(0, n_int_max_target + 1):
            t0 = time.time()

            contrib_up = compute_vertex_3pt_single_nint(
                n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
            contrib_dn = compute_vertex_3pt_single_nint(
                n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)

            B_level = float(contrib_up - contrib_dn)
            cum_B += B_level
            all_B[n_int] = B_level

            elapsed = time.time() - t0

            F2_S = cum_B / V_mag / schwinger if V_mag != 0 else 0
            delta = F2_S - 1.0

            print(f"  n_int={n_int:2d}: B={B_level:12.4e}  "
                  f"cum_B={cum_B:12.6e}  "
                  f"delta={delta:.12f}  ({elapsed:.1f}s)", flush=True)

            # Adaptive stopping
            if n_int >= 6 and cum_B != 0:
                ratio = abs(B_level / cum_B)
                if ratio < 1e-10:
                    consecutive_small += 1
                else:
                    consecutive_small = 0

                if consecutive_small >= 3:
                    print(f"  CONVERGED at n_int={n_int} "
                          f"(3 consecutive |B/cum| < 1e-10)", flush=True)
                    break

        # Power-law tail correction
        tail_B, plaw_even, plaw_odd = power_law_tail(all_B, max(all_B.keys()) + 1)
        tail_f2s = tail_B / V_mag / schwinger if V_mag != 0 else 0

        delta_raw = cum_B / V_mag / schwinger - 1.0 if V_mag != 0 else 0
        delta_corrected = delta_raw + tail_f2s

        print(f"\n  RESULT n_ext={n_ext}:")
        print(f"    delta_raw = {delta_raw:.12f}")
        print(f"    tail_correction = {tail_f2s:.6e}")
        print(f"    delta_corrected = {delta_corrected:.12f}")
        if plaw_even:
            print(f"    power_law even: C={plaw_even['C']:.4e}, p={plaw_even['p']:.2f}")
            print(f"    power_law odd:  C={plaw_odd['C']:.4e}, p={plaw_odd['p']:.2f}")
        print(flush=True)

        results[n_ext] = {
            'delta': delta_raw,
            'delta_corrected': delta_corrected,
            'lam': lam,
            'lam2': lam2,
            'V_mag': V_mag,
            'cum_B': cum_B,
            'n_int_max': max(all_B.keys()),
            'all_B': {str(k): v for k, v in all_B.items()},
            'tail_B': tail_B,
            'power_law_even': plaw_even,
            'power_law_odd': plaw_odd,
        }

        # Save checkpoint after each n_ext
        save_checkpoint(results)

    # === Final analysis ===
    print("\n" + "=" * 70)
    print("  FINAL CURVATURE EXPANSION ANALYSIS")
    print("=" * 70)

    n_ext_list = sorted(results.keys())
    for n in n_ext_list:
        r = results[n]
        print(f"  n_ext={n}: lam={r['lam']:.1f}, 1/lam^2={1/r['lam2']:.8f}, "
              f"delta={r.get('delta_corrected', r['delta']):.12f}")

    x_vals = np.array([1.0 / results[n]['lam2'] for n in n_ext_list])
    y_vals = np.array([results[n].get('delta_corrected', results[n]['delta'])
                       for n in n_ext_list])

    # Fixed c1 = 1/2 fits
    y_residual = y_vals - 0.5 * x_vals

    for n_params in range(1, len(n_ext_list)):
        X = np.column_stack([x_vals**(k+2) for k in range(n_params)])
        c_fit, _, _, _ = np.linalg.lstsq(X, y_residual, rcond=None)
        labels = [f"c{k+2}" for k in range(n_params)]
        print(f"\n  {n_params}-param fit (c1=1/2 fixed):")
        for label, val in zip(labels, c_fit):
            print(f"    {label} = {val:.12f}")

    # PSLQ on the best c2
    try:
        import mpmath
        mpmath.mp.dps = 30

        # Use the highest-order fit available
        if len(n_ext_list) >= 4:
            X3 = np.column_stack([x_vals**2, x_vals**3, x_vals**4])
            c_best, _, _, _ = np.linalg.lstsq(X3, y_residual, rcond=None)
            c2_best = c_best[0]
            c3_best = c_best[1]
            print(f"\n  Best c2 (3-param, c1 fixed) = {c2_best:.12f}")
            print(f"  Best c3 = {c3_best:.12f}")
        elif len(n_ext_list) >= 3:
            X2 = np.column_stack([x_vals**2, x_vals**3])
            c_best, _, _, _ = np.linalg.lstsq(X2, y_residual, rcond=None)
            c2_best = c_best[0]
            c3_best = c_best[1]
            print(f"\n  Best c2 (2-param, c1 fixed) = {c2_best:.12f}")
            print(f"  Best c3 = {c3_best:.12f}")
        else:
            X1 = x_vals.reshape(-1, 1)**2
            c_best, _, _, _ = np.linalg.lstsq(X1, y_residual, rcond=None)
            c2_best = c_best[0]
            c3_best = None
            print(f"\n  Best c2 (1-param, c1 fixed) = {c2_best:.12f}")

        c2_mp = mpmath.mpf(str(c2_best))

        print("\n  === CANDIDATE TESTS ===")
        candidates = [
            ("7/40 = 7*Delta", mpmath.mpf(7)/40),
            ("19/100 - 41*pi^2/25200", mpmath.mpf(19)/100 - 41*mpmath.pi**2/25200),
            ("1/6", mpmath.mpf(1)/6),
            ("4/23", mpmath.mpf(4)/23),
        ]
        for label, val in candidates:
            gap = float(c2_mp - val)
            print(f"  {label} = {float(val):.12f}, gap = {gap:.6e}")

        # PSLQ
        print("\n  === PSLQ ===")
        basis_sets = [
            ("c2, 1, pi^2", [c2_mp, mpmath.mpf(1), mpmath.pi**2]),
            ("c2, 1, B*Delta, F*Delta, F/B",
             [c2_mp, mpmath.mpf(1), mpmath.mpf(42)/40,
              mpmath.pi**2/240, mpmath.pi**2/252]),
        ]
        for label, basis in basis_sets:
            rel = mpmath.pslq(basis, tol=1e-8, maxcoeff=500)
            if rel is not None:
                terms = []
                for r, l in zip(rel, label.split(', ')):
                    if r != 0:
                        terms.append(f"({r})*{l}")
                print(f"  {label}: {' + '.join(terms)} = 0")
            else:
                print(f"  {label}: no relation found")

    except Exception as e:
        print(f"  PSLQ error: {e}")

    save_checkpoint(results)
    print("\n  DONE.", flush=True)


if __name__ == '__main__':
    main()
