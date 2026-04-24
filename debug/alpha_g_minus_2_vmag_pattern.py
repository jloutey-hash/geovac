"""
Analyze the V_magnetic pattern in detail.

The exact expressions show a beautiful structure:
  n_ext=1: V_mag = 2*sqrt(2)/3 - 2*sqrt(3)/9
  n_ext=2: V_mag = 2*sqrt(15)/9 - sqrt(2)/3
  n_ext=3: V_mag = sqrt(6)/3 - 2*sqrt(15)/15
  n_ext=4: V_mag = 2*sqrt(35)/15 - 2*sqrt(6)/9
  n_ext=5: V_mag = 4*sqrt(3)/9 - 2*sqrt(35)/21
  n_ext=6: V_mag = 2*sqrt(7)/7 - sqrt(3)/3
  n_ext=7: V_mag = sqrt(5)/3 - 2*sqrt(7)/9

There's a two-term structure: V_mag = A_n * sqrt(p_n) - B_n * sqrt(q_n).
Find the pattern in the coefficients and sqrt arguments.
"""

from sympy import Rational, sqrt, simplify, S, nsimplify, cancel
from sympy.physics.wigner import clebsch_gordan
import json


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


def tree_level_probe_per_channel(n, j, q_probe=1):
    """Return V_tree(+1/2) and V_tree(-1/2) per channel."""
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)

    if not vertex_allowed(n, n, q_probe):
        return {}

    chs = get_channels(n, n, q_probe, jL, jR, jL, jR)
    if not chs:
        return {}

    results = {}
    for jpL, jpR in chs:
        ch_up = S.Zero
        ch_dn = S.Zero
        for mpL in half_ints(jpL):
            for mpR in half_ints(jpR):
                ch_up += vertex_amp_pol(jL, jR, j, Rational(1, 2),
                                        jL, jR, j, Rational(1, 2),
                                        jpL, jpR, mpL, mpR)
                ch_dn += vertex_amp_pol(jL, jR, j, Rational(-1, 2),
                                        jL, jR, j, Rational(-1, 2),
                                        jpL, jpR, mpL, mpR)
        results[(jpL, jpR)] = {
            'up': simplify(ch_up),
            'dn': simplify(ch_dn),
            'magnetic': simplify(ch_up - ch_dn),
            'charge': simplify((ch_up + ch_dn) / 2),
        }
    return results


print("=" * 70)
print("V_MAGNETIC STRUCTURE ANALYSIS")
print("=" * 70)

j_ext = Rational(1, 2)

# Compute V per channel for n_ext=1..10
for n in range(1, 11):
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    lam = Rational(2 * n + 3, 2)

    print(f"\n--- n_ext={n}, (jL,jR)=({jL},{jR}), lambda={lam} ---")

    if not vertex_allowed(n, n, 1):
        print("  q_probe=1 not allowed (parity)")
        continue

    ch_data = tree_level_probe_per_channel(n, j_ext, q_probe=1)

    total_magnetic = S.Zero
    for (jpL, jpR), d in ch_data.items():
        print(f"  Channel ({jpL}, {jpR}):")
        print(f"    V_up   = {d['up']}")
        print(f"    V_dn   = {d['dn']}")
        print(f"    V_mag  = {d['magnetic']} = {float(d['magnetic']):.10f}")
        print(f"    V_chg  = {d['charge']} = {float(d['charge']):.10f}")
        total_magnetic += d['magnetic']

    total_magnetic = simplify(total_magnetic)
    print(f"  Total V_magnetic = {total_magnetic}")
    print(f"                   = {float(total_magnetic):.12f}")

    # Analyze the structure: V_mag * lambda
    vlam = simplify(total_magnetic * lam)
    print(f"  V_magnetic * lambda = {vlam}")
    print(f"                      = {float(vlam):.12f}")

    # V_mag^2
    v2 = simplify(total_magnetic**2)
    print(f"  V_magnetic^2 = {v2} = {float(v2):.12f}")

    # Try to express V_mag in terms of jL*(jL+1) and jR*(jR+1)
    # CG sum formulas for j=1/2, q=1 coupling should have a
    # closed form in terms of 6j or Racah coefficients
    jL_cas = jL * (jL + 1)
    jR_cas = jR * (jR + 1)
    dim_L = 2 * jL + 1
    dim_R = 2 * jR + 1
    dim_prod = dim_L * dim_R
    print(f"  jL(jL+1)={jL_cas}, jR(jR+1)={jR_cas}, dimL={dim_L}, dimR={dim_R}")

    # Try: V_mag ~ 1 / sqrt(jL*(jL+1))
    ratio = simplify(total_magnetic * sqrt(jL_cas))
    print(f"  V_mag * sqrt(jL*(jL+1)) = {float(ratio):.10f}")

    # Try: V_mag ~ 1 / (dim_L * dim_R) ?
    ratio2 = simplify(total_magnetic * dim_prod)
    print(f"  V_mag * (2jL+1)(2jR+1) = {float(ratio2):.10f}")

    # Try: V_mag ~ 1 / lambda
    ratio3 = simplify(total_magnetic * lam)
    print(f"  V_mag * lambda = {float(ratio3):.10f}")


print("\n\n--- Summary: V_magnetic * lambda ---")
vlam_vals = []
for n in range(1, 11):
    if not vertex_allowed(n, n, 1):
        continue
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    lam = Rational(2 * n + 3, 2)
    ch_data = tree_level_probe_per_channel(n, j_ext, q_probe=1)
    total_mag = S.Zero
    for d in ch_data.values():
        total_mag += d['magnetic']
    total_mag = simplify(total_mag)
    vlam = simplify(total_mag * lam)
    vlam_f = float(vlam)
    vlam_vals.append((n, float(lam), vlam_f))
    print(f"  n={n}: V_mag*lam = {vlam_f:.10f}")

# Check if V_mag * lambda -> constant
print(f"\n  Ratios of consecutive V_mag*lambda:")
for i in range(1, len(vlam_vals)):
    r = vlam_vals[i][2] / vlam_vals[i - 1][2]
    print(f"  n=({vlam_vals[i-1][0]},{vlam_vals[i][0]}): ratio = {r:.8f}")

# Check V_mag * lambda^(3/2) or similar
print(f"\n  V_mag * lambda^(3/2):")
for n, lam_f, vlam_f in vlam_vals:
    print(f"  n={n}: {vlam_f * lam_f**0.5:.10f}")

# Limit: V_mag * lambda -> 4/3?
from sympy import Rational as R
print(f"\n  Comparison: 4/3 = {4/3:.10f}")
print(f"  Last V_mag*lambda = {vlam_vals[-1][2]:.10f}")
print(f"  4/3 - V_mag*lambda[last] = {4/3 - vlam_vals[-1][2]:.6e}")
