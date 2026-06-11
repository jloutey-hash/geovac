"""
m_j-resolved vertex correction on S³ — anomalous magnetic moment probe.

Resolves the SO(4) CG coefficients at each vertex to compute
Λ(m_j=+1/2) and Λ(m_j=-1/2) separately. If they differ, the
difference B = L(+1/2) - L(-1/2) is the spin-dependent piece ∝ F₂.

Two modes:
  - "joint": same photon at both vertices (one-loop vertex correction)
  - "factored": independent photons (existing dirac_vertex.py structure)
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt
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
    The photon polarization is FIXED — no sum over (mgL, mgR).
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


def vertex_amp_summed(j_sL, j_sR, j_s, mj_s,
                      j_tL, j_tR, j_t, mj_t,
                      jgL, jgR):
    """CG amplitude summed over ALL photon polarizations."""
    total = S.Zero
    for mgL in half_ints(jgL):
        for mgR in half_ints(jgR):
            total += vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                                    j_tL, j_tR, j_t, mj_t,
                                    jgL, jgR, mgL, mgR)
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


def compute_single_contribution(n_ext, n_int, j_ext, mj_ext, mode='joint'):
    """Compute vertex correction from a SINGLE (n_ext, n_int) pair.

    mode='joint': one photon, joint polarization sum
    mode='factored': independent photons (existing code structure)
    """
    # Ext: positive chirality
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    # Int: positive chirality (same as external — gamma^mu preserves chirality)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2*n_int + 3, 2)
    lam4 = lam**4

    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    total = S.Zero

    for j_int in half_ints(j_int_max):
        if j_int < j_int_min:
            continue
        for mj_int in half_ints(j_int):

            if mode == 'joint':
                # One photon: same q, same polarization at both vertices
                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int
                for q in range(q_lo, q_hi + 1):
                    chs1 = get_channels(n_ext, n_int, q,
                                        jE_L, jE_R, jI_L, jI_R)
                    chs2 = get_channels(n_int, n_ext, q,
                                        jI_L, jI_R, jE_L, jE_R)
                    mu_q = Rational(q * (q + 2))

                    for jgL1, jgR1 in chs1:
                        for jgL2, jgR2 in chs2:
                            # Joint polarization sum
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
                                    v2 = vertex_amp_pol(
                                        jI_L, jI_R, j_int, mj_int,
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jgL2, jgR2, mgL, mgR)
                                    total += v1 * v2 / (lam4 * mu_q)

            elif mode == 'factored':
                # Independent photons (existing code structure)
                # V1 summed, V2 summed, product / lam^4
                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int

                v1_sum = S.Zero
                for q1 in range(q_lo, q_hi + 1):
                    for jgL, jgR in get_channels(n_ext, n_int, q1,
                                                  jE_L, jE_R, jI_L, jI_R):
                        v1_sum += vertex_amp_summed(
                            jE_L, jE_R, j_ext, mj_ext,
                            jI_L, jI_R, j_int, mj_int,
                            jgL, jgR)

                v2_sum = S.Zero
                for q2 in range(q_lo, q_hi + 1):
                    for jgL, jgR in get_channels(n_int, n_ext, q2,
                                                  jI_L, jI_R, jE_L, jE_R):
                        v2_sum += vertex_amp_summed(
                            jI_L, jI_R, j_int, mj_int,
                            jE_L, jE_R, j_ext, mj_ext,
                            jgL, jgR)

                total += v1_sum * v2_sum / lam4

    return total


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    n_max = 4

    print("=" * 60)
    print("m_j-resolved vertex correction on S^3")
    print(f"External: n_CH={n_ext}, j={j_ext}")
    print(f"Internal levels: n_CH = 0..{n_max}")
    print("=" * 60)

    results = {}

    for mode in ['joint', 'factored']:
        print(f"\n--- Mode: {mode} ---")
        t0 = time.time()

        val_up = S.Zero
        val_dn = S.Zero

        for n_int in range(n_max + 1):
            t1 = time.time()
            contrib_up = compute_single_contribution(
                n_ext, n_int, j_ext, Rational(1, 2), mode=mode)
            contrib_dn = compute_single_contribution(
                n_ext, n_int, j_ext, Rational(-1, 2), mode=mode)

            val_up += contrib_up
            val_dn += contrib_dn

            elapsed = time.time() - t1
            if contrib_up != 0 or contrib_dn != 0:
                print(f"  n_int={n_int}: up={float(contrib_up):.8f}, "
                      f"dn={float(contrib_dn):.8f}, "
                      f"diff={float(contrib_up - contrib_dn):.2e} "
                      f"({elapsed:.1f}s)")
            else:
                print(f"  n_int={n_int}: zero (forbidden by selection rules) "
                      f"({elapsed:.1f}s)")

        B = val_up - val_dn
        total_time = time.time() - t0

        print(f"\n  L(+1/2) = {float(val_up):.10f}")
        print(f"  L(-1/2) = {float(val_dn):.10f}")
        print(f"  B = L(+1/2) - L(-1/2) = {float(B):.2e}")
        if B != 0:
            print(f"  B exact = {B}")
        print(f"  Time: {total_time:.1f}s")

        results[mode] = {
            'val_up': str(val_up),
            'val_dn': str(val_dn),
            'val_up_float': float(val_up),
            'val_dn_float': float(val_dn),
            'B': float(B),
            'B_exact': str(B),
        }

    # Cross-check with shell-summed
    from geovac.dirac_vertex import vertex_correction_shell_summed
    shell_val = vertex_correction_shell_summed(n_ext, n_max)
    print(f"\nShell-summed L = {float(shell_val):.10f}")

    for mode in ['joint', 'factored']:
        r = results[mode]
        avg = (r['val_up_float'] + r['val_dn_float']) / 2
        if float(shell_val) != 0:
            ratio = avg / float(shell_val)
            print(f"  {mode} avg / shell = {ratio:.6f}")

    results['shell_summed'] = float(shell_val)
    results['n_ext'] = n_ext
    results['n_max'] = n_max

    with open('debug/data/alpha_g_minus_2_mj.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved to debug/data/alpha_g_minus_2_mj.json")


if __name__ == '__main__':
    main()
