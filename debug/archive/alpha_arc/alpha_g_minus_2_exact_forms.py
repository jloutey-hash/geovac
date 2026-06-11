"""
Compute exact sympy closed forms for V_magnetic(n_ext) and B(n_ext, n_int=1)
at each n_ext value (1 through 5).

The anomalous magnetic moment on S^3 depends on n_ext through the mode-dependent
vertex correction. We want to find whether F2_dominant / Schwinger * lambda^2
approaches a clean algebraic constant as n_ext -> infinity.

All computation uses sympy Rational arithmetic -- no floats until final display.
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt, simplify, nsimplify, pi, latex
import sys
import json
import time

sys.path.insert(0, '.')


def half_ints(j):
    """Generate all m values for angular momentum j."""
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    """Single-polarization vertex amplitude from SU(2)_L x SU(2)_R CG coefficients."""
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
    """Check SO(4) vertex selection rule: n1 + n2 + q must be odd."""
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, j_sL, j_sR, j_tL, j_tR):
    """Get allowed photon channels (jgL, jgR) for the transition."""
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


def tree_level_probe_magnetic_exact(n, j, q_probe=1):
    """Compute V_magnetic = V_tree(+1/2) - V_tree(-1/2) as exact sympy expression.

    Returns the exact symbolic result, not a float.
    """
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
    return simplify(total_up - total_dn)


def compute_B_magnetic_single_nint_exact(n_ext, j_ext, n_int, q_probe=1):
    """Compute B_magnetic(n_ext, n_int) = B(+1/2) - B(-1/2) as exact sympy expression.

    This is the magnetic part of the one-loop vertex correction with a single
    internal fermion mode n_int.
    """
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)
    lam_int = Rational(2 * n_int + 3, 2)
    lam4 = lam_int ** 4
    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_chs = get_channels(n_int, n_int, q_probe, jI_L, jI_R, jI_L, jI_R)
    if not probe_chs:
        return S.Zero

    # Compute B for m_j = +1/2 and -1/2
    results = {}
    for mj_ext in [Rational(1, 2), Rational(-1, 2)]:
        total = S.Zero
        for j_int in half_ints(j_int_max):
            if j_int < j_int_min:
                continue
            for mj_int in half_ints(j_int):
                for mj_int_prime in half_ints(j_int):
                    # Probe vertex (internal -> internal via probe)
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

                    # Loop photon: sum over allowed q
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

        results[mj_ext] = total

    return simplify(results[Rational(1, 2)] - results[Rational(-1, 2)])


def main():
    print("=" * 70)
    print("Exact algebraic forms for V_magnetic and B(n_int=1)")
    print("=" * 70)

    j_ext = Rational(1, 2)
    schwinger_exact_coeff = Rational(1, 1)  # F2/Schwinger is dimensionless

    all_results = []

    for n_ext in range(1, 8):
        lam = Rational(2 * n_ext + 3, 2)
        print(f"\n{'='*60}")
        print(f"n_ext = {n_ext},  lambda = {lam} = {float(lam)}")
        print(f"{'='*60}")

        # 1. Tree-level magnetic coupling
        t0 = time.time()
        V_mag = tree_level_probe_magnetic_exact(n_ext, j_ext, q_probe=1)
        t_vmag = time.time() - t0
        V_mag_simplified = simplify(V_mag)

        print(f"\n  V_magnetic (exact) = {V_mag_simplified}")
        print(f"  V_magnetic (float) = {float(V_mag_simplified):.12f}")
        print(f"  Time: {t_vmag:.1f}s")

        # 2. Vertex correction B at n_int=1 (the dominant contribution)
        #    n_int=0 has even parity so vertex_allowed(n_int=0, n_int=0, q=1) is false
        #    (0+0+1=1 is odd, so it IS allowed... let's check)
        n_int = 1
        t0 = time.time()
        B_mag = compute_B_magnetic_single_nint_exact(n_ext, j_ext, n_int, q_probe=1)
        t_bmag = time.time() - t0
        B_mag_simplified = simplify(B_mag)

        print(f"\n  B(n_int={n_int}) (exact) = {B_mag_simplified}")
        print(f"  B(n_int={n_int}) (float) = {float(B_mag_simplified):.12e}")
        print(f"  Time: {t_bmag:.1f}s")

        # Also compute n_int=0 for completeness
        t0 = time.time()
        B_mag_0 = compute_B_magnetic_single_nint_exact(n_ext, j_ext, 0, q_probe=1)
        t_b0 = time.time() - t0
        B_mag_0_simplified = simplify(B_mag_0)
        print(f"\n  B(n_int=0) (exact) = {B_mag_0_simplified}")
        print(f"  B(n_int=0) (float) = {float(B_mag_0_simplified):.12e}")
        print(f"  Time: {t_b0:.1f}s")

        # 3. F2_dominant = B(n_int=1) / V_magnetic
        if V_mag_simplified != 0 and B_mag_simplified != 0:
            F2_dom = simplify(B_mag_simplified / V_mag_simplified)
            print(f"\n  F2(n_int=1) = B/V_mag (exact) = {F2_dom}")
            print(f"  F2(n_int=1) (float) = {float(F2_dom):.12e}")

            # 4. Check lambda^2 scaling
            lam2_F2 = simplify(lam ** 2 * F2_dom)
            print(f"\n  lambda^2 * F2 (exact) = {lam2_F2}")
            print(f"  lambda^2 * F2 (float) = {float(lam2_F2):.12f}")

            # Also check lambda^(1.82) approximately
            lam_f = float(lam)
            F2_f = float(F2_dom)
            print(f"  lambda^1.82 * F2 (float) = {lam_f**1.82 * F2_f:.12f}")
            print(f"  lambda^2 * F2 (float) = {lam_f**2 * F2_f:.12f}")
        else:
            F2_dom = S.Zero
            lam2_F2 = S.Zero

        # Also compute B(n_int=0) contribution to F2
        if V_mag_simplified != 0 and B_mag_0_simplified != 0:
            F2_dom_0 = simplify(B_mag_0_simplified / V_mag_simplified)
            print(f"\n  F2(n_int=0) = B/V_mag (float) = {float(F2_dom_0):.12e}")

        entry = {
            'n_ext': n_ext,
            'lambda': str(lam),
            'lambda_float': float(lam),
            'V_magnetic_exact': str(V_mag_simplified),
            'V_magnetic_float': float(V_mag_simplified),
            'B_nint1_exact': str(B_mag_simplified),
            'B_nint1_float': float(B_mag_simplified),
            'B_nint0_exact': str(B_mag_0_simplified),
            'B_nint0_float': float(B_mag_0_simplified),
            'F2_nint1_exact': str(F2_dom),
            'F2_nint1_float': float(F2_dom),
            'lam2_F2_exact': str(lam2_F2),
            'lam2_F2_float': float(lam2_F2),
        }
        all_results.append(entry)

    # Summary table
    print(f"\n\n{'='*70}")
    print("SUMMARY TABLE")
    print(f"{'='*70}")
    print(f"{'n_ext':>5} {'lambda':>8} {'V_magnetic':>14} {'F2(n_int=1)':>14} {'lam^2*F2':>14}")
    print("-" * 60)
    for r in all_results:
        print(f"{r['n_ext']:5d} {r['lambda_float']:8.1f} {r['V_magnetic_float']:14.8f} "
              f"{r['F2_nint1_float']:14.6e} {r['lam2_F2_float']:14.8f}")

    # Check if lam^2 * F2 approaches a constant
    print(f"\n--- Checking lambda^2 * F2 trend ---")
    lam2_f2_values = [r['lam2_F2_float'] for r in all_results]
    for i, r in enumerate(all_results):
        print(f"  n_ext={r['n_ext']}: lam^2 * F2 = {r['lam2_F2_float']:.10f}")
    if len(lam2_f2_values) >= 2:
        ratio = lam2_f2_values[-1] / lam2_f2_values[-2] if lam2_f2_values[-2] != 0 else 0
        print(f"  Ratio of last two: {ratio:.6f}")

    # Print exact expressions for inspection
    print(f"\n--- Exact V_magnetic expressions ---")
    for r in all_results:
        print(f"  n_ext={r['n_ext']}: {r['V_magnetic_exact']}")

    print(f"\n--- Exact B(n_int=1) expressions ---")
    for r in all_results:
        print(f"  n_ext={r['n_ext']}: {r['B_nint1_exact']}")

    print(f"\n--- Exact lambda^2 * F2 expressions ---")
    for r in all_results:
        print(f"  n_ext={r['n_ext']}: {r['lam2_F2_exact']}")

    # Try to identify the limit using nsimplify
    print(f"\n--- Algebraic identification of lambda^2 * F2 ---")
    for r in all_results:
        val = r['lam2_F2_float']
        try:
            guess = nsimplify(val, rational=False, tolerance=1e-6)
            print(f"  n_ext={r['n_ext']}: {val:.10f} ~ {guess}")
        except Exception:
            print(f"  n_ext={r['n_ext']}: {val:.10f} (no identification)")

    # Save
    output = {
        'description': 'Exact algebraic forms for V_magnetic and B(n_int=1) at each n_ext',
        'results': all_results,
    }
    with open('debug/data/alpha_g_minus_2_exact_forms.json', 'w') as f:
        json.dump(output, f, indent=2)
    print("\nSaved to debug/data/alpha_g_minus_2_exact_forms.json")


if __name__ == '__main__':
    main()
