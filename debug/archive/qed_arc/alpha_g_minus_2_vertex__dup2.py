"""
Three-point vertex correction on S^3 for anomalous magnetic moment.

The one-loop vertex correction has the topology:

    ext_in --[V1]-- int --[V2]-- ext_out
                |            |
                +--- gamma --+

where gamma is the LOOP photon and there is an EXTERNAL probe photon
at one vertex. The external photon breaks SO(4) -> SO(3), enabling
m_j-dependent amplitudes.

Physical structure:
  - External electron: |n_ext, j_ext, m_j> (positive chirality)
  - Internal electron: |n_int, j_int, m_j_int> (same chirality, gamma^mu preserves it)
  - Loop photon: mode q_loop with SO(4) vector-harmonic channels
  - Probe photon: mode q_probe = 1 (lowest mode, uniform B-field analog)
    with FIXED polarization (mg_L, mg_R) to break SO(4) symmetry

The vertex correction amplitude for specific external m_j is:

  Lambda(m_j) = sum_{n_int, j_int, m_j_int, q_loop, q_probe}
      V1(ext_in -> int; loop_photon) * V2(int -> ext_out; loop_photon)
      * V_probe(ext -> ext; probe_photon)
      / (lambda_int^4 * mu_loop)

The anomalous magnetic moment piece is B = Lambda(+1/2) - Lambda(-1/2).
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt, pi
import sys, json, time

sys.path.insert(0, '.')


def half_ints(j):
    """All m values from -j to +j."""
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    """CG amplitude for a SPECIFIC photon polarization (mgL, mgR).

    Returns the 4-CG product summed over internal (mL1) decomposition.
    """
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
    """Return list of allowed (jgL, jgR) photon channels."""
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


def compute_vertex_correction_3pt(n_ext, j_ext, mj_ext,
                                   n_max_int, q_probe=1):
    """Three-point vertex correction with external probe photon.

    The one-loop vertex has:
      - Left vertex: external electron absorbs probe photon AND loop photon
      - Right vertex: internal electron emits loop photon back to external
      - OR: factored as probe at one vertex, loop at another

    Physical topology (vertex correction of the current vertex):

        ext(m_j) --[probe+loop]--> int(m_j_int) --[loop]--> ext(m_j)

    But actually the correct one-loop vertex correction inserts the
    loop ON TOP of the external probe vertex. The simplest approach:

    Gamma^mu_{corrected}(p', p) = sum_int [
        V(ext -> int; loop) * G_int * V(int -> ext; loop)
    ] * V_probe(ext; mu)

    where V_probe(ext; mu) is the tree-level probe coupling that
    gives the m_j dependence.

    In practice on S^3, the probe photon at q=1 has (jg_L, jg_R)
    channels and definite polarization. We compute:

    Lambda(m_j) = sum over {int states, loop photon, probe polarization}
        of the vertex-corrected matrix element.

    For the anomalous magnetic moment, we need the piece of the vertex
    function that is proportional to sigma^{mu nu} q_nu, which on S^3
    corresponds to the m_j-dependent part of the probe coupling
    modified by the loop.

    Approach: compute the one-loop corrected vertex for the probe photon.
    The probe photon at q_probe=1 with specific (mgL, mgR) defines a
    "magnetic field direction". The vertex correction modifies this coupling.

    Full one-loop corrected probe vertex element:
        <ext, m_j'| Gamma^{1-loop}_{probe} |ext, m_j>
        = sum_{int, q_loop} [
            <ext, m_j'| V_{probe} |int, m_j_int>
            * 1/(lambda_int^4 * mu_loop)
            * sum over loop-photon polarizations of
              <int, m_j_int| V_{loop} |ext, m_j> * <ext, m_j| V_{loop} |int, m_j_int>
          ]

    But this double-counts. The standard vertex correction is:

    Gamma^mu(p', p) = integral [
        gamma^nu S(k) gamma^mu S(k-q) gamma_nu D(p-k)
    ]

    On S^3, this becomes a sum over internal modes and loop photon modes:

    Lambda(m_j_out, m_j_in; probe) = sum_{n_int, j_int, m_j_int, q_loop}
        V(ext_out -> int; loop) * V(int -> ext_in; probe) * V_loop_propagator
        / (lambda_int^4 * mu_loop)

    Wait -- this is a three-vertex diagram with one loop:
        ext_in(m_j) -> [V_probe] -> int1(m_j_1)
                                          |
                        [loop photon propagator]
                                          |
                     int2(m_j_2) -> [V_loop_out] -> ext_out(m_j')
                                          |
                         int1 = int2 (same internal line)

    No, the standard one-loop vertex correction is:

    ext_in(p) ---gamma^mu--- int(k) ---gamma^nu--- ext_out(p')
                   |                      |
                   +--- photon (p-k) -----+

    where gamma^mu is the PROBE vertex (external current) and gamma^nu are
    the TWO QED vertices where the loop photon attaches. The internal electron
    propagates between the two loop-photon vertices.

    On S^3 this is:

    ext_in(m_j) --[V1: absorb loop photon]--> int(m_j_int) --[V2: emit loop photon, absorb probe]--> ext_out(m_j')

    Actually the probe vertex is between the two loop-photon vertices.
    Let me use the standard ordering:

    ext_in --[loop vertex 1]--> int --[probe vertex]--> int --[loop vertex 2]--> ext_out

    No, in the one-loop vertex correction there are exactly three vertices
    and two propagators (one electron, one photon in the loop):

    ext_in(p) --[VERT: emit loop photon q]--> int(p-q) --[VERT: probe mu]--> int(p'-q) --[VERT: absorb loop photon q]--> ext_out(p')

    But wait, the internal momenta are different on the two sides of the
    probe vertex, so the internal propagator BEFORE the probe has momentum (p-q)
    and AFTER the probe has momentum (p'-q). For forward scattering p=p',
    the internal line is the same on both sides.

    For the ANOMALOUS MAGNETIC MOMENT, we need the piece that is linear in
    the probe photon momentum transfer. At zero momentum transfer (p=p'),
    we get the charge form factor F1, not F2. The anomalous moment comes from
    the first-order response to the probe photon's polarization/momentum.

    On S^3, the "momentum transfer" of the probe photon is its angular mode
    number q_probe. For q_probe=1 (dipole), the probe photon carries angular
    momentum 1, which CAN distinguish m_j=+1/2 from m_j=-1/2.

    Simplest correct structure: the vertex correction with a q=1 probe
    photon inserted at one vertex. The loop correction modifies the
    coupling of the external state to this probe. The m_j-dependent
    piece of this modified coupling IS the anomalous magnetic moment.

    So the corrected probe coupling for external state m_j is:

    Gamma(m_j) = sum_{n_int, j_int, m_j_int, q_loop}
        sum_{loop photon polarizations (mg)}
        sum_{probe photon polarization (mp)}
        V(ext,m_j -> int,m_j_int; loop,mg) * V(int,m_j_int -> ext,m_j; probe,mp AND loop,-mg)
        / (lambda_int^4 * mu_loop)

    But this isn't right either -- the probe and loop photon are at
    DIFFERENT vertices.

    Let me just implement the simplest physically correct structure:

    FORWARD vertex correction with probe:
    For m_j_out = m_j_in = m_j (diagonal, forward scattering):

    Lambda(m_j; mp_L, mp_R) = sum_{n_int, j_int, m_j_int, q_loop} [
        V1(ext,m_j | loop,mgL,mgR -> int,m_j_int)
      * V2(int,m_j_int | probe,mpL,mpR -> int,m_j_int')  [this changes m_j_int!]
      * V3(int,m_j_int' | loop,-mgL,-mgR -> ext,m_j)     [loop photon reabsorbed]
      / (lambda_int^4 * mu_loop)
    ]

    Hmm, this has THREE vertices and the probe in the middle. But actually
    the one-loop vertex correction diagram has only TWO internal vertices
    (where the loop photon attaches) and ONE external current insertion
    (the probe). The external current insertion doesn't create a new vertex
    in the loop -- it's on the EXTERNAL line.

    Let me reconsider. The standard QED vertex correction modifies the
    electromagnetic form factor. The TREE-LEVEL coupling of the external
    electron to a probe photon of mode q_probe is:

    V_tree(m_j, m_j'; probe_pol) = CG amplitude

    The ONE-LOOP CORRECTED coupling is:

    V_1loop(m_j, m_j'; probe_pol) = sum_{int, loop}
        V(ext,m_j -> int,m_j_int; loop)
        * G_int(n_int)  [= 1/lambda_int^2]
        * V_tree(int,m_j_int, m_j_int'; probe_pol)
        * G_int(n_int)  [= 1/lambda_int^2]
        * V(int,m_j_int' -> ext,m_j'; loop)
        * D_loop(q)     [= 1/mu_q]

    where V_tree in the middle is the tree-level coupling of the INTERNAL
    electron to the probe photon.

    For the form factor we take m_j' = m_j (forward diagonal):

    Lambda(m_j; probe) = sum_{n_int, j_int, m_j_int, m_j_int', q_loop}
        V(ext,m_j -> int,m_j_int; loop)
        * V_tree(int,m_j_int -> int,m_j_int'; probe)
        * V(int,m_j_int' -> ext,m_j; loop)
        / (lambda_int^4 * mu_loop)

    The internal state has TWO magnetic quantum numbers: m_j_int (before probe)
    and m_j_int' (after probe). They're connected by the probe photon's
    angular momentum.

    THIS is the structure that can give m_j-dependent results.
    The probe photon's polarization (mpL, mpR) determines which
    m_j_int -> m_j_int' transitions are allowed, and the CG coefficients
    at the loop-photon vertices depend on m_j, creating the asymmetry.
    """
    # External: positive chirality
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)

    # Probe photon at q_probe = 1
    # For q=1: channels are (jgL, jgR) = (1, 0) or (0, 1)
    # (1,0) means jgL=1, jgR=0 — a "left-handed" probe
    # (0,1) means jgL=0, jgR=1 — a "right-handed" probe
    # For a magnetic field (spin-sensitive probe), use one definite polarization

    total = S.Zero
    nonzero_count = 0

    for n_int in range(n_max_int + 1):
        # Internal: same chirality (gamma^mu preserves)
        jI_L = Rational(n_int + 1, 2)
        jI_R = Rational(n_int, 2)

        lam = Rational(2*n_int + 3, 2)
        lam4 = lam**4

        j_int_min = abs(jI_L - jI_R)
        j_int_max = jI_L + jI_R

        for j_int in half_ints(j_int_max):
            if j_int < j_int_min:
                continue

            for mj_int in half_ints(j_int):
                for mj_int_prime in half_ints(j_int):

                    # Probe photon coupling: int(m_j_int) -> int(m_j_int')
                    # This is a SAME-STATE coupling (n_int -> n_int, same chirality)
                    # with probe photon q_probe
                    # Need: vertex_allowed(n_int, n_int, q_probe) must be true
                    if not vertex_allowed(n_int, n_int, q_probe):
                        continue

                    probe_chs = get_channels(n_int, n_int, q_probe,
                                              jI_L, jI_R, jI_L, jI_R)
                    if not probe_chs:
                        continue

                    # Sum over probe photon polarizations
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

                    # Loop photon: ext -> int and int' -> ext
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
                                # Joint loop-photon polarization sum
                                for mgL in half_ints(jgL1):
                                    for mgR in half_ints(jgR1):
                                        if abs(mgL) > jgL2 or abs(mgR) > jgR2:
                                            continue

                                        # V1: ext(m_j) -> int(m_j_int) via loop
                                        v1 = vertex_amp_pol(
                                            jE_L, jE_R, j_ext, mj_ext,
                                            jI_L, jI_R, j_int, mj_int,
                                            jgL1, jgR1, mgL, mgR)
                                        if v1 == 0:
                                            continue

                                        # V3: int(m_j_int') -> ext(m_j) via loop
                                        v3 = vertex_amp_pol(
                                            jI_L, jI_R, j_int, mj_int_prime,
                                            jE_L, jE_R, j_ext, mj_ext,
                                            jgL2, jgR2, mgL, mgR)
                                        if v3 == 0:
                                            continue

                                        contrib = v1 * probe_amp * v3 / (lam4 * mu_q)
                                        total += contrib
                                        if contrib != 0:
                                            nonzero_count += 1

    return total, nonzero_count


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    n_max = 3  # smaller cutoff first for speed

    print("=" * 60)
    print("Three-point vertex correction on S^3")
    print("(one-loop correction to probe-photon coupling)")
    print(f"External: n_CH={n_ext}, j={j_ext}")
    print(f"Internal levels: n_CH = 0..{n_max}")
    print(f"Probe photon: q=1")
    print("=" * 60)

    results = {}

    for q_probe in [1]:
        print(f"\n--- Probe q={q_probe} ---")

        # First check: does the probe coupling itself exist?
        # n_int -> n_int with q_probe
        for n_test in range(n_max + 1):
            allowed = vertex_allowed(n_test, n_test, q_probe)
            print(f"  n_int={n_test}: probe coupling allowed = {allowed} "
                  f"(parity: {n_test}+{n_test}+{q_probe}={2*n_test+q_probe}, "
                  f"odd={bool((2*n_test+q_probe)%2)})")

        t0 = time.time()

        val_up, count_up = compute_vertex_correction_3pt(
            n_ext, j_ext, Rational(1, 2), n_max, q_probe=q_probe)
        val_dn, count_dn = compute_vertex_correction_3pt(
            n_ext, j_ext, Rational(-1, 2), n_max, q_probe=q_probe)

        elapsed = time.time() - t0

        B = val_up - val_dn

        print(f"\n  L(+1/2) = {float(val_up):.10e} ({count_up} nonzero terms)")
        print(f"  L(-1/2) = {float(val_dn):.10e} ({count_dn} nonzero terms)")
        print(f"  B = L(+1/2) - L(-1/2) = {float(B):.6e}")
        if B != 0:
            print(f"  B exact = {B}")
        print(f"  Time: {elapsed:.1f}s")

        results[f'q_probe_{q_probe}'] = {
            'val_up': float(val_up),
            'val_dn': float(val_dn),
            'B': float(B),
            'B_exact': str(B),
            'count_up': count_up,
            'count_dn': count_dn,
        }

    # Also try q_probe = 2 for comparison
    print("\n--- Probe q=2 ---")
    for n_test in range(n_max + 1):
        allowed = vertex_allowed(n_test, n_test, 2)
        print(f"  n_int={n_test}: probe coupling allowed = {allowed}")

    t0 = time.time()
    val_up_2, count_up_2 = compute_vertex_correction_3pt(
        n_ext, j_ext, Rational(1, 2), n_max, q_probe=2)
    val_dn_2, count_dn_2 = compute_vertex_correction_3pt(
        n_ext, j_ext, Rational(-1, 2), n_max, q_probe=2)
    elapsed = time.time() - t0
    B2 = val_up_2 - val_dn_2

    print(f"\n  L(+1/2) = {float(val_up_2):.10e} ({count_up_2} nonzero terms)")
    print(f"  L(-1/2) = {float(val_dn_2):.10e} ({count_dn_2} nonzero terms)")
    print(f"  B = L(+1/2) - L(-1/2) = {float(B2):.6e}")
    print(f"  Time: {elapsed:.1f}s")

    results['q_probe_2'] = {
        'val_up': float(val_up_2),
        'val_dn': float(val_dn_2),
        'B': float(B2),
        'count_up': count_up_2,
        'count_dn': count_dn_2,
    }

    results['n_ext'] = n_ext
    results['n_max'] = n_max

    with open('debug/data/alpha_g_minus_2_vertex.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved to debug/data/alpha_g_minus_2_vertex.json")


if __name__ == '__main__':
    main()
