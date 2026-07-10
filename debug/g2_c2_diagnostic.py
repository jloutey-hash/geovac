"""
Diagnostic: understand the multi-n_ext normalization failure.

The j=1/2 vertex correction at n_ext >= 2 gives F2/Schwinger << 1,
approaching 0 as n_ext increases. This means either:
  (a) The j=1/2 state at higher levels has suppressed F2 (physical, but
      then the curvature expansion coefficients are j-dependent), or
  (b) The normalization V_magnetic is wrong for n_ext >= 2.

This diagnostic computes:
1. V_magnetic and the one-loop vertex for ALL j values at n_ext=2
2. A traced (j-summed) quantity that should have universal coefficients
3. Compares the j=1/2 and j=3/2 sectors

If the j=3/2 sector gives a DIFFERENT ratio, the coefficients are j-dependent
and we need the trace. If they're similar, the code has a bug.
"""

import sys, json, time
import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt

sys.path.insert(0, '.')


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_allowed(n1, n2, q):
    return (n1 + n2 + q) % 2 == 1 and abs(n1 - n2) <= q <= n1 + n2


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        cg1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if cg1 == 0:
            continue
        for mL2 in half_ints(j_tL):
            mR2 = mj_t - mL2
            if abs(mR2) > j_tR:
                continue
            cg2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
            if cg2 == 0:
                continue
            mL_g = mL1 - mL2
            mR_g = mR1 - mR2
            if abs(mL_g) > jgL or abs(mR_g) > jgR:
                continue
            if abs(mgL - mL_g) > Rational(1, 100):
                continue
            if abs(mgR - mR_g) > Rational(1, 100):
                continue
            cgL = clebsch_gordan(j_tL, jgL, j_sL, mL2, mgL, mL1)
            cgR = clebsch_gordan(j_tR, jgR, j_sR, mR2, mgR, mR1)
            total += cg1 * cg2 * cgL * cgR
    return total


def get_channels(n1, n2, q, j1L, j1R, j2L, j2R):
    if not vertex_allowed(n1, n2, q):
        return []
    chs = []
    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                      (Rational(q-1, 2), Rational(q+1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j1L - jgL) <= j2L <= j1L + jgL and
                abs(j1R - jgR) <= j2R <= j1R + jgR):
            chs.append((jgL, jgR))
    return chs


def compute_V_magnetic_j(n_ext, j_ext, mj_up, mj_dn):
    """Tree-level magnetic vertex for arbitrary j state."""
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)

    probe_chs = get_channels(n_ext, n_ext, 1, jE_L, jE_R, jE_L, jE_R)

    up_amp = S.Zero
    for jpL, jpR in probe_chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                a = vertex_amp_pol(
                    jE_L, jE_R, j_ext, mj_up,
                    jE_L, jE_R, j_ext, mj_up,
                    jpL, jpR, mpL, mpR)
                up_amp += a

    dn_amp = S.Zero
    for jpL, jpR in probe_chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                a = vertex_amp_pol(
                    jE_L, jE_R, j_ext, mj_dn,
                    jE_L, jE_R, j_ext, mj_dn,
                    jpL, jpR, mpL, mpR)
                dn_amp += a

    return float(up_amp - dn_amp)


def compute_vertex_nint(n_ext, j_ext, mj_ext, n_int, q_probe=1):
    """One-loop vertex contribution from a single n_int level."""
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
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * np.pi)

    print("=" * 70)
    print("  MULTI-j DIAGNOSTIC FOR g-2 NORMALIZATION")
    print("=" * 70)

    # === n_ext=1: only j=1/2, baseline check ===
    print("\n--- n_ext=1 (n_CH=0, only j=1/2) ---")
    j_half = Rational(1, 2)
    V1_half = compute_V_magnetic_j(1, j_half, Rational(1,2), Rational(-1,2))
    print(f"  V_mag(j=1/2) = {V1_half:.10f}")

    # Compute first few n_int levels
    B1_total = 0.0
    for ni in range(0, 6):
        t0 = time.time()
        up = compute_vertex_nint(1, j_half, Rational(1,2), ni)
        dn = compute_vertex_nint(1, j_half, Rational(-1,2), ni)
        B = float(up - dn)
        B1_total += B
        ratio = B1_total / (V1_half * schwinger) if V1_half != 0 else 0
        print(f"  n_int={ni}: B={B:12.6e}  cum={B1_total:12.6e}  "
              f"F2/S={ratio:.6f}  ({time.time()-t0:.1f}s)")

    # === n_ext=2: j=1/2 AND j=3/2 ===
    print("\n--- n_ext=2 (n_CH=1, j=1/2 and j=3/2) ---")

    # j=1/2 sector
    V2_half = compute_V_magnetic_j(2, j_half, Rational(1,2), Rational(-1,2))
    print(f"  V_mag(j=1/2, mj=+1/2 vs -1/2) = {V2_half:.10f}")

    # j=3/2 sector (use mj = +3/2 vs +1/2 for dipole)
    j_3half = Rational(3, 2)
    V2_3half_maxmin = compute_V_magnetic_j(2, j_3half, Rational(3,2), Rational(-3,2))
    V2_3half_mid = compute_V_magnetic_j(2, j_3half, Rational(3,2), Rational(1,2))
    print(f"  V_mag(j=3/2, mj=+3/2 vs -3/2) = {V2_3half_maxmin:.10f}")
    print(f"  V_mag(j=3/2, mj=+3/2 vs +1/2) = {V2_3half_mid:.10f}")

    # Traced V_mag: sum over all (j, m_j) weighted by m_j
    # For j=1/2: m_j * V(m_j) summed = (+1/2)*V(+1/2) + (-1/2)*V(-1/2)
    # This equals (1/2) * V_mag(j=1/2) if V is linear in m_j
    # For the trace, sum m_j * Λ(j, m_j) over all j, m_j
    print("\n  Computing traced V_mag (m_j-weighted sum)...")
    V_trace = 0.0
    for j_ext_val, j_ext_sym in [(0.5, Rational(1,2)), (1.5, Rational(3,2))]:
        for mj_val in half_ints(j_ext_sym):
            # Tree-level vertex for this (j, m_j) state
            jE_L = Rational(3, 2)  # n_ext=2
            jE_R = Rational(1, 1)
            probe_chs = get_channels(2, 2, 1, jE_L, jE_R, jE_L, jE_R)
            amp = S.Zero
            for jpL, jpR in probe_chs:
                for mpL in half_ints(jpL):
                    for mpR in half_ints(jpR):
                        a = vertex_amp_pol(
                            jE_L, jE_R, j_ext_sym, mj_val,
                            jE_L, jE_R, j_ext_sym, mj_val,
                            jpL, jpR, mpL, mpR)
                        amp += a
            V_trace += float(mj_val) * float(amp)
            if abs(float(mj_val)) > 0.01:
                print(f"    j={j_ext_val}, mj={float(mj_val):+.1f}: "
                      f"V_tree={float(amp):.10f}, mj*V={float(mj_val)*float(amp):.10f}")
    print(f"  V_mag_trace (Sigma mj * V_tree) = {V_trace:.10f}")

    # One-loop vertex for first few n_int, both j sectors
    print("\n  One-loop vertex, j=1/2 sector (n_ext=2):")
    B2_half_total = 0.0
    for ni in range(0, 6):
        t0 = time.time()
        up = compute_vertex_nint(2, j_half, Rational(1,2), ni)
        dn = compute_vertex_nint(2, j_half, Rational(-1,2), ni)
        B = float(up - dn)
        B2_half_total += B
        ratio = B2_half_total / (V2_half * schwinger) if V2_half != 0 else 0
        print(f"    n_int={ni}: B={B:12.6e}  cum={B2_half_total:12.6e}  "
              f"F2/S={ratio:.6f}  ({time.time()-t0:.1f}s)")

    print("\n  One-loop vertex, j=3/2 sector (n_ext=2, mj=+3/2 vs -3/2):")
    B2_3half_total = 0.0
    for ni in range(0, 6):
        t0 = time.time()
        up = compute_vertex_nint(2, j_3half, Rational(3,2), ni)
        dn = compute_vertex_nint(2, j_3half, Rational(-3,2), ni)
        B = float(up - dn)
        B2_3half_total += B
        ratio = B2_3half_total / (V2_3half_maxmin * schwinger) if V2_3half_maxmin != 0 else 0
        print(f"    n_int={ni}: B={B:12.6e}  cum={B2_3half_total:12.6e}  "
              f"F2/S={ratio:.6f}  ({time.time()-t0:.1f}s)")

    # Traced one-loop
    print("\n  One-loop vertex, TRACED (n_ext=2, first 4 n_int):")
    B_trace_total = 0.0
    for ni in range(0, 4):
        t0 = time.time()
        B_trace_level = 0.0
        for j_ext_val, j_ext_sym in [(0.5, Rational(1,2)), (1.5, Rational(3,2))]:
            for mj_val in half_ints(j_ext_sym):
                amp = compute_vertex_nint(2, j_ext_sym, mj_val, ni)
                B_trace_level += float(mj_val) * float(amp)
        B_trace_total += B_trace_level
        ratio = B_trace_total / (V_trace * schwinger) if V_trace != 0 else 0
        print(f"    n_int={ni}: B_trace={B_trace_level:12.6e}  "
              f"cum={B_trace_total:12.6e}  "
              f"F2/S_trace={ratio:.6f}  ({time.time()-t0:.1f}s)")

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  n_ext=1 j=1/2:   V_mag={V1_half:.6f}  cum_B(0..5)={B1_total:.6e}")
    print(f"  n_ext=2 j=1/2:   V_mag={V2_half:.6f}  cum_B(0..5)={B2_half_total:.6e}")
    print(f"  n_ext=2 j=3/2:   V_mag={V2_3half_maxmin:.6f}  cum_B(0..5)={B2_3half_total:.6e}")
    print(f"  n_ext=2 traced:  V_mag={V_trace:.6f}  cum_B(0..3)={B_trace_total:.6e}")

    r1 = B1_total / (V1_half * schwinger) if V1_half != 0 else 0
    r2h = B2_half_total / (V2_half * schwinger) if V2_half != 0 else 0
    r2_3h = B2_3half_total / (V2_3half_maxmin * schwinger) if V2_3half_maxmin != 0 else 0
    r2_tr = B_trace_total / (V_trace * schwinger) if V_trace != 0 else 0

    print(f"\n  F2/Schwinger ratios (partial sums, not converged):")
    print(f"    n_ext=1, j=1/2:          {r1:.6f}")
    print(f"    n_ext=2, j=1/2:          {r2h:.6f}")
    print(f"    n_ext=2, j=3/2:          {r2_3h:.6f}")
    print(f"    n_ext=2, traced:         {r2_tr:.6f}")

    print(f"\n  KEY QUESTION: does j=3/2 give a ratio closer to 1 than j=1/2?")
    print(f"  If yes -> curvature expansion is j-dependent, need trace")
    print(f"  If no  -> code bug in vertex computation for n_ext >= 2")

    # Save
    results = {
        'n_ext_1': {
            'V_mag_j_half': V1_half,
            'cum_B_j_half': B1_total,
            'ratio_j_half': r1,
        },
        'n_ext_2': {
            'V_mag_j_half': V2_half,
            'V_mag_j_3half': V2_3half_maxmin,
            'V_mag_trace': V_trace,
            'cum_B_j_half': B2_half_total,
            'cum_B_j_3half': B2_3half_total,
            'cum_B_trace': B_trace_total,
            'ratio_j_half': r2h,
            'ratio_j_3half': r2_3h,
            'ratio_trace': r2_tr,
        }
    }
    with open('debug/data/g2_c2_diagnostic.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved to debug/data/g2_c2_diagnostic.json")


if __name__ == '__main__':
    main()
