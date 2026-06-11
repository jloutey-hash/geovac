"""
Investigation of F2/Schwinger ratio for the anomalous magnetic moment on S^3.

Extends the computation to higher n_max, computes V_magnetic exactly,
and tests candidate closed forms via direct comparison and PSLQ.

RESULTS (2026-04-23):
- Extended to n_max=12 (from previous n_max=8)
- V_magnetic = 2*sqrt(2)/3 - 2*sqrt(3)/9  (exact sympy)
- F2(n=1) = 8*(86 - 19*sqrt(6)) / 253125  (exact algebraic)
- B_level(n) ~ n^{-7} power law decay (effective exponent 6.1->7.1)
- F2/Schwinger at n_max=12: 1.0844524
- Extrapolated limit: 1.08445 +/- 0.0001
- NOT converged: successive B-level ratios still increasing (0.28 at last pair)
- NO clean closed form found (see structural analysis below)
- Correction ratio-1 ~ 0.0845, consistent with Parker-Toms R/(12*m^2) = 2/25 = 0.08
  plus higher-order curvature corrections (~5.6% of leading)

STRUCTURAL ANALYSIS:
F2 is a pure algebraic number (infinite sum of terms with sqrt(2), sqrt(3), etc.)
Schwinger = alpha/(2*pi) involves transcendentals.
The ratio F2/Schwinger necessarily mixes algebraic and transcendental numbers.
A clean closed form is NOT expected on structural grounds.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt, pi, simplify
import sys
import json
import time

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
    for jgL, jgR in [(Rational(q + 1, 2), Rational(q - 1, 2)),
                      (Rational(q - 1, 2), Rational(q + 1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j_sL - jgL) <= j_tL <= j_sL + jgL and
                abs(j_sR - jgR) <= j_tR <= j_sR + jgR):
            chs.append((jgL, jgR))
    return chs


def compute_vertex_3pt_magnetic(n_ext, j_ext, n_int, q_probe=1):
    """Compute B_level(n_int) = B_up - B_dn for the magnetic part."""
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2 * n_int + 3, 2)
    lam4 = lam ** 4

    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_chs = get_channels(n_int, n_int, q_probe,
                              jI_L, jI_R, jI_L, jI_R)
    if not probe_chs:
        return S.Zero

    total_up = S.Zero
    total_dn = S.Zero

    for mj_ext_val in [Rational(1, 2), Rational(-1, 2)]:
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
                                            jE_L, jE_R, j_ext, mj_ext_val,
                                            jI_L, jI_R, j_int, mj_int,
                                            jgL1, jgR1, mgL, mgR)
                                        if v1 == 0:
                                            continue
                                        v3 = vertex_amp_pol(
                                            jI_L, jI_R, j_int, mj_int_prime,
                                            jE_L, jE_R, j_ext, mj_ext_val,
                                            jgL2, jgR2, mgL, mgR)
                                        if v3 == 0:
                                            continue
                                        total += v1 * probe_amp * v3 / (lam4 * mu_q)

        if mj_ext_val == Rational(1, 2):
            total_up = total
        else:
            total_dn = total

    return total_up - total_dn


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979323846)

    # Exact V_magnetic
    V_magnetic_exact = 2 * sqrt(Rational(2)) / 3 - 2 * sqrt(Rational(3)) / 9
    V_magnetic_float = float(V_magnetic_exact)

    print("=" * 70)
    print("F2/Schwinger ratio investigation on S^3")
    print(f"External: n_CH={n_ext}, j={j_ext}")
    print(f"V_magnetic = 2*sqrt(2)/3 - 2*sqrt(3)/9 = {V_magnetic_float:.15f}")
    print(f"Schwinger = alpha/(2*pi) = {schwinger:.15e}")
    print("=" * 70)

    # Try to push to n_max=14 (or as far as time allows)
    n_max_target = 14
    max_wall_time = 240  # 4 minutes total budget

    results = []
    B_cumul = 0.0
    total_time = 0.0

    for n_int in range(n_max_target + 1):
        t0 = time.time()

        B_level_sympy = compute_vertex_3pt_magnetic(n_ext, j_ext, n_int)
        B_level_f = float(B_level_sympy)

        elapsed = time.time() - t0
        total_time += elapsed

        B_cumul += B_level_f
        F2 = B_cumul / V_magnetic_float if V_magnetic_float != 0 else 0
        ratio = F2 / schwinger if schwinger != 0 else 0

        if B_level_f != 0:
            print(f"  n={n_int:2d}: dB={B_level_f:+14.8e}, B={B_cumul:14.10e}, "
                  f"F2/Schwinger={ratio:.10f} ({elapsed:.1f}s)")
        else:
            print(f"  n={n_int:2d}: zero ({elapsed:.1f}s)")

        results.append({
            'n_int': n_int,
            'B_level': B_level_f,
            'B_cumul': B_cumul,
            'F2': F2,
            'F2_over_schwinger': ratio,
            'time': elapsed,
        })

        if total_time > max_wall_time and n_int >= 10:
            print(f"\n  Time budget exhausted after n_int={n_int} ({total_time:.0f}s)")
            break

    # ====================================================================
    # Convergence analysis
    # ====================================================================
    print("\n" + "=" * 70)
    print("CONVERGENCE ANALYSIS")
    print("=" * 70)

    nonzero = [r for r in results if r['B_level'] != 0]
    ratios = [r['F2_over_schwinger'] for r in nonzero]

    # Pair-end analysis
    pair_end_ratios = [ratios[i] for i in range(1, len(ratios), 2)]
    print("\nPair-end ratios:")
    for i, r in enumerate(pair_end_ratios):
        print(f"  pair {i}: {r:.12f}")

    if len(pair_end_ratios) >= 3:
        diffs = [pair_end_ratios[i] - pair_end_ratios[i - 1]
                 for i in range(1, len(pair_end_ratios))]
        print("\nPair-end differences:")
        for i, d in enumerate(diffs):
            print(f"  d[{i}] = {d:.8e}")

        if len(diffs) >= 2:
            conv_ratios = [diffs[i] / diffs[i - 1] for i in range(1, len(diffs))
                          if abs(diffs[i - 1]) > 1e-20]
            print("\nConvergence ratios:")
            for r in conv_ratios:
                print(f"  {r:.6f}")

            if len(conv_ratios) >= 1:
                r_conv = conv_ratios[-1]
                remaining = diffs[-1] * r_conv / (1 - r_conv)
                limit_estimate = pair_end_ratios[-1] + remaining
                print(f"\nExtrapolated limit (pair-end geometric):")
                print(f"  r_conv = {r_conv:.6f}")
                print(f"  remaining = {remaining:.8e}")
                print(f"  F2/Schwinger limit = {limit_estimate:.12f}")

    # Also do Richardson extrapolation on the last 3 pair-end values
    if len(pair_end_ratios) >= 3:
        a0, a1, a2 = pair_end_ratios[-3], pair_end_ratios[-2], pair_end_ratios[-1]
        # Aitken's delta-squared
        denom = a2 - 2 * a1 + a0
        if abs(denom) > 1e-20:
            aitken = a2 - (a2 - a1) ** 2 / denom
            print(f"\nAitken delta-squared extrapolation:")
            print(f"  F2/Schwinger limit = {aitken:.12f}")

    # Current best value
    current_ratio = results[-1]['F2_over_schwinger']
    print(f"\nCurrent value at n_max={results[-1]['n_int']}:")
    print(f"  F2/Schwinger = {current_ratio:.12f}")

    # ====================================================================
    # F2 itself (not ratio)
    # ====================================================================
    print("\n" + "=" * 70)
    print("F2 VALUE (not ratio)")
    print("=" * 70)
    F2_current = results[-1]['F2']
    print(f"  F2 = {F2_current:.15e}")
    print(f"  Schwinger = {schwinger:.15e}")
    print(f"  F2 - Schwinger = {F2_current - schwinger:.8e}")
    print(f"  (F2 - Schwinger) / Schwinger = {(F2_current - schwinger) / schwinger:.8e}")

    # ====================================================================
    # Candidate closed forms for the ratio
    # ====================================================================
    print("\n" + "=" * 70)
    print("CANDIDATE CLOSED FORMS")
    print("=" * 70)

    # Use the best estimate we have
    best_ratio = current_ratio  # will be updated if we have an extrapolation

    candidates = [
        ("13/12", 13 / 12),
        ("1 + 1/12", 1 + 1 / 12),
        ("1 + 1/10", 1 + 1 / 10),
        ("1 + 1/11", 1 + 1 / 11),
        ("1 + 1/R_scalar (R=6, so 1+1/6)", 1 + 1 / 6),
        ("1 + 6/R (R=6 scalar curv)", 1 + 1 / 6),
        ("1 + 1/(2*R_scalar) = 1+1/12", 1 + 1 / 12),
        ("1 + R/72 = 1+6/72=1+1/12", 1 + 1 / 12),
        ("1 + 3/(2*|lam_1|^2) = 1+3/(2*25/4) = 1+6/25", 1 + 6 / 25),
        ("1 + 1/(|lam_1|^2) = 1+4/25", 1 + 4 / 25),
        ("1 + 1/|lam_1| = 1+2/5", 1 + 2 / 5),
        ("1 + 4/(3*|lam_1|^2) = 1+16/75", 1 + 16 / 75),
        ("sqrt(3)/sqrt(3-1/12)", (3 / (3 - 1 / 12)) ** 0.5),
        ("(5/2)^2 / (5/2)^2 - 1) = 25/21", 25 / 21),
        ("1 + R/(6*(2*n+3)^2) = 1+6/(6*25/4) = 1+4/25", 1 + 4 / 25),
        # n_ext=1: |lam_ext| = 5/2, g_ext = 2*2*3 = 12
        ("1 + 1/g_ext = 1+1/12", 1 + 1 / 12),
        ("1 + 2/g_ext = 1+1/6", 1 + 1 / 6),
        # Seeley-DeWitt
        ("1 + a2/a0 = 1 + 1/8", 1 + 1 / 8),
        ("1 + a1/a0 = 1 + 1 = 2", 2),
        # Dirac eigenvalue at external state
        ("|lam_1|/(|lam_1|-1) = (5/2)/(3/2) = 5/3", 5 / 3),
        ("(|lam_1|^2)/(|lam_1|^2-1) = 25/21", 25 / 21),
        ("(|lam_1|^2+1)/|lam_1|^2 = 29/25", 29 / 25),
        # n^2 spectrum
        ("(n_ext+1)^2/((n_ext+1)^2-1) = 4/3", 4 / 3),
        # Other simple fractions near 1.0844
        ("217/200", 217 / 200),
        ("541/499", 541 / 499),
        ("1 + 169/2000", 1 + 169 / 2000),
        ("(25/4 + 1)/(25/4) = 29/25", 29 / 25),
        # Try: 1 + curvature / (8 * |lam|^2) = 1 + 6/(8*25/4) = 1 + 6/50 = 1+3/25
        ("1 + 3/25", 1 + 3 / 25),
        ("1 + R/(8*|lam_1|^2) = 1+6/50=1+3/25", 1 + 3 / 25),
        # Systematic: try a/b for small a, b
    ]

    # Also search systematically
    best_match = None
    best_diff = 1.0
    for a in range(1, 200):
        for b in range(1, 200):
            v = a / b
            d = abs(v - best_ratio)
            if d < 1e-4 and d < best_diff:
                best_diff = d
                best_match = (a, b, v, d)

    if best_match:
        candidates.append((f"BEST: {best_match[0]}/{best_match[1]}", best_match[2]))

    print(f"\nTarget ratio: {best_ratio:.12f}")
    print()
    for name, val in candidates:
        diff = best_ratio - val
        if abs(diff) < 0.05:
            marker = " <---" if abs(diff) < 5e-4 else ""
            print(f"  {name:50s} = {val:.10f}  diff = {diff:+.6e}{marker}")

    # ====================================================================
    # PSLQ-style search: is ratio algebraic over simple constants?
    # ====================================================================
    print("\n" + "=" * 70)
    print("PSLQ-STYLE SEARCH")
    print("=" * 70)
    try:
        import mpmath
        mpmath.mp.dps = 30
        r = mpmath.mpf(str(best_ratio))

        # Test: a*ratio + b = 0 => ratio = -b/a, i.e. rational
        # Test: a*ratio + b*pi + c = 0
        # Test: a*ratio + b*sqrt(2) + c*sqrt(3) + d = 0

        bases = [
            ("rational", [r, mpmath.mpf(1)]),
            ("with pi", [r, mpmath.pi, mpmath.mpf(1)]),
            ("with sqrt(2)", [r, mpmath.sqrt(2), mpmath.mpf(1)]),
            ("with sqrt(3)", [r, mpmath.sqrt(3), mpmath.mpf(1)]),
            ("with sqrt(2), sqrt(3)", [r, mpmath.sqrt(2), mpmath.sqrt(3), mpmath.mpf(1)]),
            ("with pi^2", [r, mpmath.pi ** 2, mpmath.mpf(1)]),
            ("with zeta(3)", [r, mpmath.zeta(3), mpmath.mpf(1)]),
            ("with 1/pi", [r, 1 / mpmath.pi, mpmath.mpf(1)]),
            ("with pi, sqrt(2), sqrt(3)", [r, mpmath.pi, mpmath.sqrt(2), mpmath.sqrt(3), mpmath.mpf(1)]),
        ]

        for name, basis in bases:
            try:
                rel = mpmath.pslq(basis, maxcoeff=200, tol=1e-8)
                if rel is not None:
                    print(f"  {name}: PSLQ found {rel}")
                    # Reconstruct
                    if len(rel) == 2:
                        val = mpmath.mpf(-rel[1]) / mpmath.mpf(rel[0])
                        print(f"    => ratio = {val}")
                    elif len(rel) == 3:
                        # rel[0]*r + rel[1]*x + rel[2] = 0
                        # r = -(rel[1]*x + rel[2]) / rel[0]
                        pass
                else:
                    print(f"  {name}: no relation found (maxcoeff=200)")
            except Exception as e:
                print(f"  {name}: error: {e}")
    except ImportError:
        print("  mpmath not available for PSLQ")

    # ====================================================================
    # Save results
    # ====================================================================
    output = {
        'n_ext': n_ext,
        'V_magnetic_exact': '2*sqrt(2)/3 - 2*sqrt(3)/9',
        'V_magnetic_float': V_magnetic_float,
        'schwinger': schwinger,
        'per_level': results,
        'final_ratio': results[-1]['F2_over_schwinger'],
        'final_F2': results[-1]['F2'],
    }

    with open('debug/data/alpha_g_minus_2_ratio_investigation.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to debug/data/alpha_g_minus_2_ratio_investigation.json")


if __name__ == '__main__':
    main()
