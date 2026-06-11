"""
Option C: Analytical approach to g-2 curvature coefficient c2.

Strategy:
1. Recompute B(n_int) keeping EXACT sympy rationals (not float)
2. Look for closed-form rational function pattern in n_int
3. If found, sum the series analytically to get exact delta

The B(n_int) values are exact rationals because:
- CG coefficients are algebraic (sqrt of rationals)
- Sum over all m_j of products of CGs gives 6j/9j symbols (rational)
- Division by lam^4 * mu_q is rational

Also tries:
- Richardson extrapolation on existing float data
- Wynn epsilon acceleration
- Rational function fitting
"""

from sympy.physics.wigner import clebsch_gordan, wigner_6j
from sympy import Rational, S, simplify, factor, nsimplify, Integer
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


def compute_vertex_exact(n_ext, j_ext, mj_ext, n_int, q_probe=1):
    """Same as compute_vertex_3pt_single_nint but returns EXACT sympy result."""
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


def wynn_epsilon(s_n):
    """Wynn's epsilon algorithm for sequence acceleration.
    Input: list of partial sums s_n.
    Returns accelerated estimates."""
    n = len(s_n)
    if n < 3:
        return s_n

    eps = np.zeros((n, n+1))
    eps[:, 0] = 0.0
    eps[:, 1] = np.array(s_n, dtype=float)

    for k in range(2, n+1):
        for i in range(n - k + 1):
            diff = eps[i+1, k-1] - eps[i, k-1]
            if abs(diff) < 1e-30:
                eps[i, k] = 1e30
            else:
                eps[i, k] = eps[i+1, k-2] + 1.0 / diff

    # Even columns are the accelerated values
    results = []
    for k in range(0, n, 2):
        if k+1 < n+1:
            results.append(eps[0, k+1])
    return results


def main():
    print("=" * 70)
    print("  ANALYTICAL APPROACH TO g-2 CURVATURE COEFFICIENT")
    print("=" * 70)

    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    # Load existing float data
    with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
        base_data = json.load(f)
    with open('debug/data/g2_extended_nint.json') as f:
        ext_data = json.load(f)

    V_mag = base_data['V_magnetic_float']

    all_B_float = {}
    for lev in base_data['per_level']:
        all_B_float[lev['n_int']] = lev['B_level']
    for lev in ext_data['extended_levels']:
        all_B_float[lev['n_int']] = lev['B_level']

    # ================================================================
    # Part 1: Exact rational B values for small n_int
    # ================================================================
    print("\n--- Part 1: Exact rational B(n_int) values ---")
    print("  Computing B(n_int) as exact sympy rationals...")

    exact_B = {}
    # n_int=0 is always zero (vertex_allowed fails)
    exact_B[0] = S.Zero

    for n_int in range(1, 7):  # 1..6, increasing cost
        t0 = time.time()
        print(f"\n  n_int={n_int}:", flush=True)

        c_up = compute_vertex_exact(n_ext, j_ext, Rational(1,2), n_int)
        c_dn = compute_vertex_exact(n_ext, j_ext, Rational(-1,2), n_int)

        B_exact = simplify(c_up - c_dn)
        elapsed = time.time() - t0

        exact_B[n_int] = B_exact

        # Check if rational
        is_rational = B_exact.is_rational
        if is_rational is None:
            # Try harder
            B_simplified = nsimplify(B_exact, rational=True)
            is_rational = B_simplified.is_rational
            if is_rational:
                B_exact = B_simplified
                exact_B[n_int] = B_exact

        B_float = float(B_exact)
        float_match = abs(B_float - all_B_float.get(n_int, 0))

        print(f"    B = {B_exact}")
        print(f"    float = {B_float:.10e}")
        print(f"    is_rational = {is_rational}")
        if n_int in all_B_float:
            print(f"    matches stored float: {float_match:.2e}")
        print(f"    ({elapsed:.1f}s)")

        if is_rational:
            from sympy import numer, denom
            p = B_exact.p if hasattr(B_exact, 'p') else numer(B_exact)
            q = B_exact.q if hasattr(B_exact, 'q') else denom(B_exact)
            print(f"    = {p}/{q}")
            print(f"    denominator = {q}")
            # Factor the denominator
            from sympy import factorint
            if isinstance(q, int) or q.is_integer:
                facs = factorint(int(q))
                print(f"    denom factored = {facs}")

    # ================================================================
    # Part 2: Pattern analysis on exact rationals
    # ================================================================
    print("\n\n--- Part 2: Pattern analysis ---")

    rational_B = {n: b for n, b in exact_B.items() if b != 0 and b.is_rational}

    if len(rational_B) >= 3:
        print("  Looking for rational function pattern B(n) = P(n)/Q(n)...")

        for n, b in sorted(rational_B.items()):
            lam = Rational(2*n + 3, 2)
            lam4 = lam**4
            # B * lam^4 removes the propagator factor
            b_reduced = simplify(b * lam4)
            print(f"    n={n}: B*lam^4 = {b_reduced}")

    # ================================================================
    # Part 3: Sequence acceleration on existing float data
    # ================================================================
    print("\n\n--- Part 3: Sequence acceleration on float data ---")

    # Build partial sum sequence
    ns = sorted(all_B_float.keys())
    partial_sums = []
    cumsum = 0.0
    for n in ns:
        cumsum += all_B_float[n]
        partial_sums.append(cumsum / V_mag / schwinger)  # F2/S

    print(f"  {len(partial_sums)} partial sums (n_int=0..{ns[-1]})")
    print(f"  Last partial sum F2/S = {partial_sums[-1]:.12f}")

    # Richardson extrapolation on partial sums
    print("\n  Richardson extrapolation on F2/S partial sums:")
    ps = np.array(partial_sums)
    # Use h = 1/n_int^2 as the expansion parameter (since B ~ n^{-6.5})
    # Actually, partial sum S(N) = S_inf + c/N^p + ... where p ~ 5.5
    # Use h_i = 1/N_i^5 as the parameter
    N_vals = np.array([float(n) for n in ns if n >= 1], dtype=float)
    S_vals = np.array([partial_sums[i] for i, n in enumerate(ns) if n >= 1])

    # Richardson with h = 1/n^5
    h_vals = 1.0 / N_vals**5
    m = min(len(h_vals), 8)  # use last 8 points
    h_use = h_vals[-m:]
    s_use = S_vals[-m:]

    table = np.zeros((m, m))
    table[:, 0] = s_use
    for j in range(1, m):
        for i in range(m - j):
            table[i, j] = (h_use[i+j] * table[i, j-1] - h_use[i] * table[i+1, j-1]) / (h_use[i+j] - h_use[i])

    print(f"  Using last {m} points with h = 1/n^5:")
    for j in range(m):
        print(f"    order {j}: F2/S = {table[0, j]:.12f}")

    rich_best = table[0, m-1]
    rich_delta = rich_best - 1.0
    print(f"\n  Best Richardson: delta = {rich_delta:.12f}")

    # Wynn epsilon acceleration
    print("\n  Wynn epsilon acceleration on F2/S partial sums:")
    wynn_results = wynn_epsilon(partial_sums[1:])  # skip n=0
    for i, val in enumerate(wynn_results[:6]):
        if abs(val) < 100:
            print(f"    epsilon_{2*i}: {val:.12f}")

    if len(wynn_results) >= 2:
        wynn_best = wynn_results[-1] if abs(wynn_results[-1]) < 100 else wynn_results[-2]
        wynn_delta = wynn_best - 1.0
        print(f"\n  Best Wynn: delta = {wynn_delta:.12f}")

    # Aitken delta-squared on last few partial sums
    print("\n  Aitken delta-squared on tail:")
    for i in range(max(0, len(partial_sums)-6), len(partial_sums)-2):
        s0, s1, s2 = partial_sums[i], partial_sums[i+1], partial_sums[i+2]
        denom = s2 - 2*s1 + s0
        if abs(denom) > 1e-20:
            aitken = s2 - (s2 - s1)**2 / denom
            print(f"    n={ns[i]},{ns[i+1]},{ns[i+2]}: {aitken:.12f}")

    # ================================================================
    # Part 4: Separate even/odd acceleration
    # ================================================================
    print("\n\n--- Part 4: Separate even/odd tail analysis ---")

    even_ns = [n for n in ns if n % 2 == 0 and n >= 2]
    odd_ns = [n for n in ns if n % 2 == 1 and n >= 1]

    even_B = [all_B_float[n] for n in even_ns]
    odd_B = [all_B_float[n] for n in odd_ns]

    # For each subsequence, compute n^7 * B(n) to see if it approaches a constant
    print("  n^7 * B(n) for even n:")
    for n, b in zip(even_ns, even_B):
        print(f"    n={n:3d}: n^7*B = {n**7 * b:.8f}")

    print("  n^7 * B(n) for odd n:")
    for n, b in zip(odd_ns, odd_B):
        print(f"    n={n:3d}: n^7*B = {n**7 * b:.8f}")

    # Try n^6 * B and n^7 * B to determine exact power
    print("\n  Even: checking decay power...")
    for p in [6, 6.5, 7, 7.5]:
        vals = [n**p * abs(b) for n, b in zip(even_ns[-5:], even_B[-5:])]
        ratio = vals[-1] / vals[0] if vals[0] > 0 else 0
        print(f"    n^{p:.1f} * B: {vals[0]:.6f} -> {vals[-1]:.6f} (ratio {ratio:.4f})")

    print("  Odd: checking decay power...")
    for p in [6, 6.5, 7, 7.5]:
        vals = [n**p * abs(b) for n, b in zip(odd_ns[-5:], odd_B[-5:])]
        ratio = vals[-1] / vals[0] if vals[0] > 0 else 0
        print(f"    n^{p:.1f} * B: {vals[0]:.6f} -> {vals[-1]:.6f} (ratio {ratio:.4f})")

    # ================================================================
    # Part 5: PSLQ with improved precision
    # ================================================================
    print("\n\n--- Part 5: PSLQ with accelerated values ---")

    try:
        import mpmath
        mpmath.mp.dps = 50

        # Use Richardson-accelerated delta
        delta_rich = mpmath.mpf(str(rich_delta))
        lam2 = mpmath.mpf('25/4')
        c1 = mpmath.mpf('1/2')
        c2_rich = (delta_rich - c1/lam2) * lam2**2

        print(f"  Richardson-accelerated values:")
        print(f"    delta = {float(delta_rich):.12f}")
        print(f"    c2 = {float(c2_rich):.10f}")

        # Compare with raw c2
        c2_raw = mpmath.mpf(str(ext_data['corrected_delta']))
        c2_raw = (c2_raw - float(c1/lam2)) * float(lam2**2)
        print(f"  Tail-corrected c2 = {float(c2_raw):.10f}")
        print(f"  Difference = {float(c2_rich - c2_raw):.6e}")

        # PSLQ on delta with tighter tolerance
        print(f"\n  PSLQ on Richardson delta = {float(delta_rich):.12f}:")

        bases = [
            (['delta', '1', 'pi^2'],
             [delta_rich, mpmath.mpf(1), mpmath.pi**2]),
            (['delta', '1', 'pi^2', 'pi^4'],
             [delta_rich, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4]),
            (['delta', '1', 'pi^2', 'log2'],
             [delta_rich, mpmath.mpf(1), mpmath.pi**2, mpmath.log(2)]),
            (['delta', '1', 'pi^2', 'zeta3'],
             [delta_rich, mpmath.mpf(1), mpmath.pi**2, mpmath.zeta(3)]),
            (['delta', '1', 'pi^2', 'G'],
             [delta_rich, mpmath.mpf(1), mpmath.pi**2, mpmath.catalan]),
            (['delta', '1', 'pi', 'pi^2'],
             [delta_rich, mpmath.mpf(1), mpmath.pi, mpmath.pi**2]),
        ]

        for labels, vals in bases:
            try:
                rel = mpmath.pslq(vals, tol=1e-6, maxcoeff=500)
                if rel is not None and rel[0] != 0:
                    terms = [f"{r}*{l}" for r, l in zip(rel, labels) if r != 0]
                    ratio = sum(-r*float(v) for r, v in zip(rel[1:], vals[1:])) / rel[0]
                    err = abs(ratio - float(delta_rich))
                    print(f"    {labels[1:]}: {' + '.join(terms)} = 0")
                    print(f"      => delta = {ratio:.12f} (err {err:.2e})")
                else:
                    print(f"    {labels[1:]}: no relation")
            except Exception as e:
                print(f"    {labels[1:]}: error ({e})")

    except ImportError:
        print("  (mpmath not available)")

    # ================================================================
    # Summary
    # ================================================================
    print("\n" + "=" * 70)
    print("  ANALYTICAL APPROACH SUMMARY")
    print("=" * 70)
    print(f"""
  Exact rational B values: computed for n_int=1..6
  Sequence acceleration: Richardson + Wynn epsilon applied
  PSLQ: attempted on accelerated values

  Key question: is B(n_int) a rational function of n_int?
  If yes, the series can be summed in closed form.
  If no, the best path is continued brute-force (Option A).
""")


if __name__ == '__main__':
    main()
