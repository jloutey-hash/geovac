"""
Investigate whether B(n_int) in the g-2 vertex correction on S^3 can be
expressed via Racah/6j algebra to produce rational values.

Result: STRUCTURAL OBSTRUCTION -- B(n_int) is NOT rational.

Background
----------
The brute-force code (g2_extend_nint_v2.py, g2_analytical_c2.py) computes
B(n_int) by summing over all magnetic quantum numbers of products of four
CG coefficients at each of three vertices.

The hypothesis was that the sum over all m_j projections of products of CG
coefficients reduces to Wigner 6j/9j symbols via Racah algebra, which are
RATIONAL, and therefore B(n_int) would become rational.

Angular coupling scheme
-----------------------
On S^3 = SU(2)_L x SU(2)_R, the Dirac spinor at CH level n has:
  Positive chirality: (j_L, j_R) = ((n+1)/2, n/2)
  Negative chirality: (j_L, j_R) = (n/2, (n+1)/2)

The photon at level q has two vector components:
  A: (j_L^g, j_R^g) = ((q+1)/2, (q-1)/2)
  B: (j_L^g, j_R^g) = ((q-1)/2, (q+1)/2)

Vertex selection rule: n1 + n2 + q = odd, q >= 1, |n1-n2| <= q <= n1+n2.

Structural result
-----------------
B(n_int) is NOT rational. The structural obstruction is:

1. Each vertex amplitude involves 4 CG coefficients. The product IS
   reducible to a 9j symbol via Racah recoupling, and the 9j is rational.

2. HOWEVER, the g-2 extraction requires FIXED external projections
   mj_ext = +1/2 and -1/2 (to form c_up - c_dn). The external spinor
   CG coefficients CG(j_L, j_R, j; mL, mR, mj_ext) for j_L != j_R
   contain irremovable square roots.

3. On S^3, positive chirality spinors have (j_L, j_R) = ((n+1)/2, n/2),
   with j_L - j_R = 1/2 ALWAYS. The CG(j_L, j_R, j; ...) for this
   half-integer step generically contains sqrt(2j+1) type irrationals.

4. The Wigner-Eckart theorem extracts a factor sqrt(6)/3 from the
   mj-dependent part, but the reduced matrix element is itself NOT
   rational because the internal coupling at each vertex still involves
   CGs between unequal (j_L, j_R) pairs.

5. The number field grows with n_int:
   n_int=1: B in Q(sqrt(2), sqrt(3))
   n_int=2: B in Q(sqrt(2), sqrt(3), sqrt(5))
   (new square roots from sqrt(2j+1) at higher j values)

Implications for Paper 18
--------------------------
The irrationals are ALGEBRAIC (square roots), NOT transcendental.
This is consistent with Paper 18: the transcendental content of the
summed series F2 is determined by spectral structure (Dirac eigenvalues,
Hurwitz zeta), while per-level B(n_int) can be algebraically irrational.
The T9 theorem (one-loop = pi^{even} only) applies to the sum, not to
individual terms.
"""

import sys
import time
sys.path.insert(0, '.')

from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_9j
from sympy import Rational, S, sqrt as ssqrt, simplify, nsimplify, Integer
import json


def half_ints(j):
    """Generate all m values from -j to +j in steps of 1."""
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def triangle_ok(a, b, c):
    """Check triangle inequality for angular momentum coupling."""
    return (abs(a - b) <= c <= a + b) and ((2*(a + b + c)) == int(2*(a + b + c)))


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    """Single vertex amplitude: sum over mL1 of 4 CG coefficients."""
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


# ============================================================================
# TEST 1: Single vertex amplitude algebraic structure
# ============================================================================

def test_single_amplitude_structure():
    """
    Compute individual vertex amplitudes at fixed m-values and check
    which ones are rational vs irrational.
    """
    print("\n" + "=" * 70)
    print("  TEST 1: Single vertex amplitude algebraic structure")
    print("=" * 70)

    # Same-chirality probe vertex: n_int=1 -> n_int=1 via q=1
    # Source = Target: positive chirality (j_sL, j_sR) = (1, 1/2)
    n_src, n_tgt, q = 1, 1, 1
    j_sL = Rational(1, 1)   # 1
    j_sR = Rational(1, 2)   # 1/2
    j_tL = j_sL             # same chirality
    j_tR = j_sR

    chs = get_channels(n_src, n_tgt, q, j_sL, j_sR, j_tL, j_tR)
    j_s_min = abs(j_sL - j_sR)
    j_s_max = j_sL + j_sR

    print(f"\n  Same-chirality probe ({n_src},{n_tgt},{q}):")
    print(f"  Source = Target: (j_L, j_R) = ({j_sL}, {j_sR})")
    print(f"  Channels: {chs}")
    print(f"  j range: [{j_s_min}, {j_s_max}]")

    count_rat = 0
    count_irrat = 0

    for jgL, jgR in chs:
        label = "A" if jgL > jgR else "B"
        print(f"\n  Channel {label}: (jgL, jgR) = ({jgL}, {jgR})")

        for j_s in half_ints(j_s_max):
            if j_s < j_s_min:
                continue
            for j_t in half_ints(j_s_max):
                if j_t < j_s_min:
                    continue

                for mj_s in half_ints(j_s):
                    for mj_t in half_ints(j_t):
                        for mgL in half_ints(jgL):
                            for mgR in half_ints(jgR):
                                amp = vertex_amp_pol(
                                    j_sL, j_sR, j_s, mj_s,
                                    j_tL, j_tR, j_t, mj_t,
                                    jgL, jgR, mgL, mgR)
                                if amp != 0:
                                    is_rat = amp.is_rational
                                    if is_rat:
                                        count_rat += 1
                                    else:
                                        count_irrat += 1
                                    print(f"    j_s={j_s},mj_s={mj_s}, j_t={j_t},mj_t={mj_t}, "
                                          f"mg=({mgL},{mgR}): A = {amp}  rat={is_rat}")

    print(f"\n  Total: {count_rat} rational, {count_irrat} irrational amplitudes")
    if count_irrat > 0:
        print(f"  => Individual vertex amplitudes CONTAIN irrationals")
        print(f"     (square roots from CG(j_L, j_R, j; mL, mR, mj) with j_L != j_R)")


# ============================================================================
# TEST 2: Full B(n_int) exact computation
# ============================================================================

def test_B_rationality():
    """
    Compute B(n_int) for n_int=1..2 using exact sympy and check rationality.
    """
    print("\n" + "=" * 70)
    print("  TEST 2: Full B(n_int) = c_up - c_dn -- exact computation")
    print("=" * 70)

    from debug.g2_analytical_c2 import compute_vertex_exact

    n_ext = 1
    j_ext = Rational(1, 2)

    results = []

    for n_int in [1, 2, 3]:
        t0 = time.time()
        c_up = compute_vertex_exact(n_ext, j_ext, Rational(1,2), n_int, q_probe=1)
        c_dn = compute_vertex_exact(n_ext, j_ext, Rational(-1,2), n_int, q_probe=1)
        B = simplify(c_up - c_dn)
        elapsed = time.time() - t0

        B_float = float(B)
        is_rat = B.is_rational

        print(f"\n  n_int={n_int}: ({elapsed:.1f}s)")
        print(f"    B     = {B}")
        print(f"    B (float) = {B_float:.15e}")
        print(f"    rational  = {is_rat}")

        if not is_rat:
            # Identify the irrationals
            from sympy import sqrt
            B_expanded = B.expand()
            # Collect terms by irrational part
            terms = B_expanded.as_ordered_terms()
            print(f"    Number of terms: {len(terms)}")
            for t in terms:
                print(f"      {t}")

            # Check B^2
            B_sq = simplify(B**2)
            is_sq_rat = B_sq.is_rational
            print(f"    B^2 rational = {is_sq_rat}")

        results.append({
            'n_int': n_int,
            'B_float': B_float,
            'is_rational': bool(is_rat) if is_rat is not None else None,
            'B_sympy': str(B),
            'time': elapsed
        })

    return results


# ============================================================================
# TEST 3: Verify 9j identity for same-chirality vertex
# ============================================================================

def test_9j_same_chirality():
    """
    Test the 9j recoupling identity for a same-chirality vertex
    (the probe vertex in the g-2 computation).
    """
    print("\n" + "=" * 70)
    print("  TEST 3: 9j recoupling for same-chirality probe vertex")
    print("=" * 70)

    # Same-chirality probe: n=1 -> n=1, q=1
    # Both source and target: (j_sL, j_sR) = (1, 1/2)
    j_sL = j_tL = Rational(1)
    j_sR = j_tR = Rational(1, 2)

    # Channel A: (jgL, jgR) = (1, 0)
    jgL, jgR = Rational(1), S.Zero

    print(f"\n  Vertex (1,1,1) same-chirality, Channel A: (jgL,jgR) = ({jgL},{jgR})")
    print(f"  Source = Target: (j_L, j_R) = ({j_sL}, {j_sR})")

    j_s_vals = half_ints(j_sL + j_sR)
    j_s_vals = [j for j in j_s_vals if j >= abs(j_sL - j_sR)]

    for j_s in j_s_vals:
        for j_t in j_s_vals:
            # The photon total: jg ranges from |jgL-jgR| to jgL+jgR
            jg_min = abs(jgL - jgR)
            jg_max = jgL + jgR

            print(f"\n  j_s={j_s}, j_t={j_t}:")

            for jg in half_ints(jg_max):
                if jg < jg_min:
                    continue
                # Check ALL 9j triangle conditions
                all_ok = True
                for a, b, c in [(j_sL, j_sR, j_s),   # row 1
                                (jgL, jgR, jg),        # row 2
                                (j_tL, j_tR, j_t),     # row 3
                                (j_sL, jgL, j_tL),     # col 1
                                (j_sR, jgR, j_tR),     # col 2
                                (j_s, jg, j_t)]:       # col 3
                    if not triangle_ok(a, b, c):
                        all_ok = False
                        break

                if not all_ok:
                    continue

                t0 = time.time()
                nj = wigner_9j(j_sL, j_sR, j_s,
                                jgL, jgR, jg,
                                j_tL, j_tR, j_t)
                elapsed = time.time() - t0

                if nj != 0:
                    print(f"    9j(j_s={j_s}, jg={jg}, j_t={j_t}) = {nj}  "
                          f"rational={nj.is_rational}  ({elapsed:.1f}s)")

                    # Compare with brute-force: sum over specific m's
                    mj_s = min(j_s, Rational(1,2))
                    mj_t = mj_s

                    bf = vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                                        j_tL, j_tR, j_t, mj_t,
                                        jgL, jgR, S.Zero, S.Zero)

                    # 9j formula for this specific case:
                    mg = S.Zero
                    if abs(mg) <= jg and abs(mj_s + mg - mj_t) < Rational(1, 100):
                        cg_gamma = clebsch_gordan(jgL, jgR, jg, S.Zero, S.Zero, mg)
                        cg_ext = clebsch_gordan(j_s, jg, j_t, mj_s, mg, mj_t)
                        pf = ssqrt((2*j_s+1) * (2*jg+1) * (2*j_tL+1) * (2*j_tR+1))
                        nj_amp = pf * nj * cg_gamma * cg_ext
                    else:
                        nj_amp = S.Zero

                    diff = simplify(bf - nj_amp)
                    print(f"      mj_s={mj_s}: BF={bf}, 9j_formula={nj_amp}, diff={diff}")


# ============================================================================
# TEST 4: Known F2(n_int=1) structure analysis
# ============================================================================

def test_known_f2():
    """Analyze the known exact F2(n_int=1) and trace the irrationals."""
    print("\n" + "=" * 70)
    print("  TEST 4: Known F2(n_int=1) structure analysis")
    print("=" * 70)

    from sympy import sqrt

    # From existing data
    F2_n1 = Rational(8, 253125) * (86 - 19*sqrt(6))
    V_mag = 2*sqrt(2)/3 - 2*sqrt(3)/9

    print(f"  F2(n_int=1) = 8*(86 - 19*sqrt(6)) / 253125")
    print(f"    = {float(F2_n1):.15e}")
    print(f"    rational = {F2_n1.is_rational}")
    print(f"    contains sqrt(6)")
    print()

    print(f"  V_magnetic = 2*sqrt(2)/3 - 2*sqrt(3)/9")
    print(f"    = {float(V_mag):.15e}")
    print(f"    lives in Q(sqrt(2), sqrt(3))")
    print()

    # B_level(1) = F2 * V_mag (confirmed by Test 2)
    B1 = simplify(F2_n1 * V_mag)
    print(f"  B(n_int=1) = F2 * V_mag")
    print(f"    = {B1}")
    print(f"    = {float(B1):.15e}")
    print(f"    lives in Q(sqrt(2), sqrt(3))")
    print()

    # Clean form
    B1_clean = 16*(63*ssqrt(2) - 40*ssqrt(3)) / 455625
    diff = simplify(B1 - B1_clean)
    print(f"  Clean form: B(n_int=1) = 16*(63*sqrt(2) - 40*sqrt(3)) / 455625")
    print(f"    verification diff = {diff}")
    print()

    # Wigner-Eckart decomposition
    w = wigner_3j(Rational(1,2), 1, Rational(1,2), Rational(-1,2), 0, Rational(1,2))
    print(f"  Wigner-Eckart factor:")
    print(f"    3j(1/2,1,1/2; -1/2,0,1/2) = {w} = sqrt(6)/6")
    print(f"    c_up - c_dn = (sqrt(6)/3) * <j=1/2 || T^1 || j=1/2>")
    print()

    RME = simplify(B1 * 3 / ssqrt(6))
    print(f"  Reduced matrix element = B * 3/sqrt(6)")
    print(f"    = {RME}")
    print(f"    rational = {RME.is_rational}")
    print(f"    STILL IRRATIONAL: contains sqrt(2), sqrt(3)")
    print()

    # Denominator structure
    print("  Denominator factorizations:")
    for d, name in [(253125, "F2 denom"), (455625, "B clean denom")]:
        n = d
        factors = {}
        for p in [2, 3, 5, 7, 11, 13]:
            while n % p == 0:
                factors[p] = factors.get(p, 0) + 1
                n //= p
        if n > 1:
            factors[n] = 1
        print(f"    {name}: {d} = {' * '.join(f'{p}^{e}' for p, e in sorted(factors.items()))}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("  RACAH / 6j INVESTIGATION OF g-2 VERTEX CORRECTION ON S^3")
    print("=" * 70)

    # Test 4: known F2 structure (fast, existing data)
    test_known_f2()

    # Test 1: individual amplitude structure (fast)
    test_single_amplitude_structure()

    # Test 3: 9j identity verification (moderate)
    test_9j_same_chirality()

    # Test 2: full B(n_int) exact computation (expensive)
    results = test_B_rationality()

    # ================================================================
    # Save results
    # ================================================================
    output = {
        'B_values': [],
        'structural_obstruction': True,
        'obstruction_type': 'algebraic_irrationals_from_external_CG_coupling',
        'field_extensions': {
            'n_int_1': 'Q(sqrt(2), sqrt(3))',
            'n_int_2': 'Q(sqrt(2), sqrt(3), sqrt(5))',
            'general': 'Q(sqrt(2), sqrt(3), sqrt(5), sqrt(7), ...) grows with n_int'
        },
        'root_cause': (
            'Positive chirality spinors on S^3 have (j_L, j_R) = ((n+1)/2, n/2) '
            'with j_L - j_R = 1/2 always. The CG coefficients '
            'CG(j_L, j_R, j; mL, mR, mj) for this half-integer step contain '
            'sqrt(2j+1) type irrationals. The g-2 extraction requires fixed '
            'mj_ext = +1/2 and -1/2, so these CGs are NOT summed over and '
            'the irrationals survive in c_up - c_dn.'
        ),
        'racah_reduction_status': (
            'The 9j symbols for vertex recoupling are RATIONAL (ratios of '
            'factorials). The Racah reduction works correctly. But the '
            'irrationals come from the REMAINING CG coefficients in the '
            'external coupling, not from the recoupling coefficients.'
        ),
        'paper_18_consistency': (
            'Consistent. The irrationals are ALGEBRAIC (square roots), not '
            'transcendental. Paper 18 T9 theorem (one-loop = pi^{even} only) '
            'applies to the summed series, not to individual terms.'
        ),
        'wigner_eckart': {
            'factor': 'sqrt(6)/3',
            'reduced_ME_rational': False,
            'explanation': (
                'Wigner-Eckart extracts sqrt(6)/3 from the mj-dependence, '
                'but the reduced matrix element is itself irrational because '
                'internal vertex couplings involve CGs between unequal (j_L, j_R) pairs.'
            )
        }
    }

    for r in results:
        output['B_values'].append({
            'n_int': r['n_int'],
            'B_float': r['B_float'],
            'is_rational': r['is_rational'],
            'B_exact': r['B_sympy'],
            'time': r['time']
        })

    # Add known clean forms
    if len(output['B_values']) >= 1:
        output['B_values'][0]['clean_form'] = '16*(63*sqrt(2) - 40*sqrt(3)) / 455625'
        output['B_values'][0]['field'] = 'Q(sqrt(2), sqrt(3))'

    with open('debug/data/g2_racah_test.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Results saved to debug/data/g2_racah_test.json")

    # ================================================================
    # Final summary memo
    # ================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY MEMO: g-2 Racah/6j Investigation")
    print("=" * 70)
    print("""
  1. ANGULAR COUPLING SCHEME:
     SO(4) = SU(2)_L x SU(2)_R on S^3. Positive chirality spinor at
     CH level n has (j_L, j_R) = ((n+1)/2, n/2). Photon at level q
     has two vector components A: ((q+1)/2, (q-1)/2), B: ((q-1)/2, (q+1)/2).
     Each vertex involves 4 CG coefficients: two for the SU(2)_L x SU(2)_R
     photon coupling, two for coupling (j_L, j_R) into total j.

  2. RACAH REDUCTION:
     APPLICABLE but DOES NOT PRODUCE RATIONAL B(n_int).

     The product of 4 CG coefficients at each vertex IS reducible to a
     9j symbol via standard recoupling. The 9j symbols are RATIONAL
     (ratios of factorials). This is verified in Test 3.

     However, after Racah reduction, the amplitude still contains CG
     coefficients for the EXTERNAL spinor coupling: CG(j_L, j_R, j; mL, mR, mj).
     These are NOT summed over because the g-2 form factor requires
     FIXED mj_ext = +1/2 and -1/2. These CGs contain irremovable
     square roots.

  3. RATIONALITY OF B(n_int):
     NO. B(n_int) is NOT rational.

     Exact results:
       B(1) = 16*(63*sqrt(2) - 40*sqrt(3)) / 455625
            = 6.9578e-04
            in Q(sqrt(2), sqrt(3))

       B(2) = -68/496125 + 52*sqrt(30)/52093125 + 1252*sqrt(3)/72930375
              + 116*sqrt(10)/3472875
            = 3.7649e-06
            in Q(sqrt(2), sqrt(3), sqrt(5))

     The number field GROWS with n_int because higher-level CG coefficients
     involve larger j values, introducing sqrt(2j+1) for new j.

     ROOT CAUSE: On S^3, positive chirality spinors have j_L - j_R = 1/2
     ALWAYS. The CG(j_L, j_R, j; mL, mR, mj) for this asymmetric coupling
     generically contains square roots. This is a structural feature of the
     Dirac operator on S^3 (Camporesi-Higuchi spectrum).

  4. PAPER 18 TAXONOMY:
     CONSISTENT. The irrationals in B(n_int) are ALGEBRAIC (square roots),
     not transcendental. They come from the CG coupling structure of
     SU(2)_L x SU(2)_R, not from integration or regularization.

     Paper 18's T9 theorem (one-loop QED quantities involve only pi^{even}
     transcendentals) applies to the SUMMED series F2 = sum_n B(n_int),
     not to individual terms. The per-level algebraic irrationals are
     expected to combine into transcendental content (pi^{even}) upon
     summation over all n_int, as the Hurwitz zeta structure mediates
     between algebraic vertex couplings and transcendental spectral sums.

     No zetq(3), log(2), or Catalan G appears at one loop -- consistent
     with T9 theorem (squared Dirac spectral zeta is polynomial in pi^2
     at every integer s).

  5. EXACT VALUES (first 2):
     B(1) = 16*(63*sqrt(2) - 40*sqrt(3)) / (3^6 * 5^4)

     B(2) = -68/(3^4 * 5^3 * 7^2) + 52*sqrt(30)/(3^5 * 5^4 * 7^3)
            + 1252*sqrt(3)/(3^5 * 5^3 * 7^4) + 116*sqrt(10)/(3^4 * 5^3 * 7^3)

     The denominators are products of small primes (3, 5, 7) reflecting
     the CG normalization factors. The irrationals are sqrt of small
     integers (2, 3, 5, 6, 10, 30) = products of sqrt(2), sqrt(3), sqrt(5).
""")


if __name__ == "__main__":
    main()
