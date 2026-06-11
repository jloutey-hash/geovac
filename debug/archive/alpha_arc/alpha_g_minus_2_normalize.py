"""
Normalize the three-point vertex correction B by the tree-level
probe coupling to extract the anomalous magnetic moment.

The physical quantity is F2 = B / (2 * V_tree_magnetic), where
V_tree_magnetic is the tree-level magnetic coupling (the piece of
the tree-level vertex that is proportional to sigma^{mu nu}).

On S^3, the tree-level coupling of the external electron to a
q=1 probe photon is:

V_tree(m_j; probe_pol) = sum over CG coefficients

The magnetic part is the m_j-dependent piece:
V_tree_magnetic = V_tree(+1/2) - V_tree(-1/2)

Then the anomalous magnetic moment is:
a_e = F2(0) = B / V_tree_magnetic

(up to normalization factors from the probe photon wavefunction).
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt
import sys, json, time

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


def tree_level_probe_coupling(n, j, mj, q_probe=1):
    """Tree-level coupling of |n, j, m_j> to a probe photon of mode q.

    The coupling is n -> n (same state) via the probe.
    Sum over ALL probe polarizations and ALL photon channels.
    """
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)

    if not vertex_allowed(n, n, q_probe):
        return S.Zero

    chs = get_channels(n, n, q_probe, jL, jR, jL, jR)
    if not chs:
        return S.Zero

    total = S.Zero
    for jpL, jpR in chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                amp = vertex_amp_pol(jL, jR, j, mj,
                                     jL, jR, j, mj,
                                     jpL, jpR, mpL, mpR)
                total += amp
    return total


def tree_level_probe_by_polarization(n, j, mj, q_probe=1):
    """Tree-level probe coupling broken down by polarization."""
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)

    if not vertex_allowed(n, n, q_probe):
        return {}

    chs = get_channels(n, n, q_probe, jL, jR, jL, jR)
    results = {}

    for jpL, jpR in chs:
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                amp = vertex_amp_pol(jL, jR, j, mj,
                                     jL, jR, j, mj,
                                     jpL, jpR, mpL, mpR)
                key = f"({float(jpL)},{float(jpR)})_({float(mpL)},{float(mpR)})"
                results[key] = float(amp)
    return results


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    print("=" * 70)
    print("Tree-level probe coupling and normalization analysis")
    print(f"External state: n_CH={n_ext}, j={j_ext}")
    print("=" * 70)

    # Tree-level coupling for m_j = +1/2 and -1/2
    print("\n--- Tree-level probe coupling (q_probe=1) ---")
    for mj in [Rational(1, 2), Rational(-1, 2)]:
        v_tree = tree_level_probe_coupling(n_ext, j_ext, mj, q_probe=1)
        print(f"  V_tree(m_j={float(mj):+.1f}) = {float(v_tree):.8f}")
        print(f"    exact: {v_tree}")

        # By polarization
        by_pol = tree_level_probe_by_polarization(n_ext, j_ext, mj, q_probe=1)
        for k, v in by_pol.items():
            if v != 0:
                print(f"    pol {k}: {v:.6f}")

    v_up = tree_level_probe_coupling(n_ext, j_ext, Rational(1, 2))
    v_dn = tree_level_probe_coupling(n_ext, j_ext, Rational(-1, 2))
    v_magnetic = float(v_up - v_dn)
    v_charge = float(v_up + v_dn) / 2

    print(f"\n  V_tree(+1/2) - V_tree(-1/2) = {v_magnetic:.8f}  (magnetic)")
    print(f"  [V_tree(+1/2) + V_tree(-1/2)]/2 = {v_charge:.8f}  (charge)")

    # Now the vertex correction B (from convergence run)
    B_value = 7.026510e-04  # from convergence run at n_max=5

    print(f"\n--- Normalization ---")
    print(f"  B (vertex correction) = {B_value:.6e}")
    print(f"  Schwinger target = {schwinger:.6e}")
    print(f"  B / Schwinger = {B_value/schwinger:.4f}")

    if v_magnetic != 0:
        F2 = B_value / v_magnetic
        print(f"\n  F2 = B / V_magnetic = {F2:.6e}")
        print(f"  F2 / (alpha/(2*pi)) = {F2/schwinger:.4f}")
        print(f"  F2 * 2*pi / alpha = {F2 * 2 * 3.14159265358979 / alpha:.4f}")

    if v_charge != 0:
        ratio_BC = B_value / v_charge
        print(f"\n  B / V_charge = {ratio_BC:.6e}")
        print(f"  (B/V_charge) / (alpha/(2*pi)) = {ratio_BC/schwinger:.4f}")

    # Also check: what if the probe should NOT be summed over polarizations?
    # In QED the magnetic moment comes from a SPECIFIC component of the current.
    # Let's check each polarization separately.
    print("\n--- Per-polarization analysis ---")
    jL = Rational(n_ext + 1, 2)
    jR = Rational(n_ext, 2)
    chs = get_channels(n_ext, n_ext, 1, jL, jR, jL, jR)
    print(f"  Probe channels at q=1: {[(float(a),float(b)) for a,b in chs]}")

    for jpL, jpR in chs:
        print(f"\n  Channel ({float(jpL)}, {float(jpR)}):")
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                amp_up = vertex_amp_pol(jL, jR, j_ext, Rational(1, 2),
                                        jL, jR, j_ext, Rational(1, 2),
                                        jpL, jpR, mpL, mpR)
                amp_dn = vertex_amp_pol(jL, jR, j_ext, Rational(-1, 2),
                                        jL, jR, j_ext, Rational(-1, 2),
                                        jpL, jpR, mpL, mpR)
                diff = amp_up - amp_dn
                if amp_up != 0 or amp_dn != 0:
                    print(f"    pol ({float(mpL):+.1f}, {float(mpR):+.1f}): "
                          f"up={float(amp_up):.6f}, dn={float(amp_dn):.6f}, "
                          f"diff={float(diff):.6f}")

    # Check different external states
    print("\n--- Different external states ---")
    for n_test in range(5):
        j_test = Rational(1, 2)
        if not vertex_allowed(n_test, n_test, 1):
            print(f"  n_CH={n_test}: probe forbidden (parity)")
            continue
        jL_t = Rational(n_test + 1, 2)
        jR_t = Rational(n_test, 2)
        j_min_t = abs(jL_t - jR_t)
        j_max_t = jL_t + jR_t
        # j_ext = 1/2 only exists if j_min <= 1/2 <= j_max
        if j_test < j_min_t or j_test > j_max_t:
            print(f"  n_CH={n_test}: j=1/2 not in range [{float(j_min_t)}, {float(j_max_t)}]")
            continue
        v_up_t = tree_level_probe_coupling(n_test, j_test, Rational(1, 2))
        v_dn_t = tree_level_probe_coupling(n_test, j_test, Rational(-1, 2))
        print(f"  n_CH={n_test}: V_up={float(v_up_t):.6f}, V_dn={float(v_dn_t):.6f}, "
              f"diff={float(v_up_t-v_dn_t):.6f}")


if __name__ == '__main__':
    main()
