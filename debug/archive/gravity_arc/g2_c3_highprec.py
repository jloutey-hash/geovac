"""
g2_c3_highprec.py -- High-precision c3 extraction using mpmath arithmetic.

Strategy: reuse sympy CG (exact rational) but accumulate in mpmath at 80 dps
instead of float64. This avoids the 7-digit cancellation loss in c3 extraction.

For n_int=0..50 (existing data), recompute the vertex sum at high precision.
Then extract c3 with the full 80-dps precision of the partial sum.
"""

import json
import sys
import time
from pathlib import Path
from functools import lru_cache

import mpmath

mpmath.mp.dps = 80

sys.path.insert(0, '.')

DATA_DIR = Path(__file__).parent / "data"

# Paper 2 invariants at high precision
B_HOPF = mpmath.mpf(42)
F_HOPF = mpmath.pi**2 / 6
DELTA_HOPF = mpmath.mpf(1) / 40
ALPHA = mpmath.mpf('7.2973525693e-3')
SCHWINGER = ALPHA / (2 * mpmath.pi)

C1 = mpmath.mpf(1) / 2
C2 = (2 - B_HOPF * DELTA_HOPF - F_HOPF * DELTA_HOPF - F_HOPF / B_HOPF) / 5

# n_ext=1 kinematics
LAM = mpmath.mpf(5) / 2
LAM2 = LAM**2
X = 1 / LAM2  # 4/25
X3 = X**3

print(f"c1 = {C1}")
print(f"c2 = {mpmath.nstr(C2, 20)}")
print(f"x  = {mpmath.nstr(X, 20)}  (= 4/25)")
print(f"x^3= {mpmath.nstr(X3, 20)}")

# CG from sympy (exact rational, cached)
from sympy.physics.wigner import clebsch_gordan as _cg_sympy
from sympy import Rational as R


@lru_cache(maxsize=2_000_000)
def cg_mp(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
    """CG coefficient as mpmath.mpf (exact rational converted at 80 dps)."""
    if m1_2 + m2_2 != M_2:
        return mpmath.mpf(0)
    if abs(j1_2 - j2_2) > J_2 or J_2 > j1_2 + j2_2:
        return mpmath.mpf(0)
    if abs(m1_2) > j1_2 or abs(m2_2) > j2_2 or abs(M_2) > J_2:
        return mpmath.mpf(0)
    val = _cg_sympy(R(j1_2, 2), R(j2_2, 2), R(J_2, 2),
                     R(m1_2, 2), R(m2_2, 2), R(M_2, 2))
    return mpmath.mpf(val.evalf(80))


def half_ints_2(two_j):
    return list(range(-two_j, two_j + 1, 2))


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, jsL2, jsR2, jtL2, jtR2):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL2, jgR2 in [(q + 1, q - 1), (q - 1, q + 1)]:
        if jgL2 < 0 or jgR2 < 0:
            continue
        if (abs(jsL2 - jgL2) <= jtL2 <= jsL2 + jgL2 and
                abs(jsR2 - jgR2) <= jtR2 <= jsR2 + jgR2):
            chs.append((jgL2, jgR2))
    return chs


def vertex_amp_mp(jsL2, jsR2, js2, mjs2,
                  jtL2, jtR2, jt2, mjt2,
                  jgL2, jgR2, mgL2, mgR2):
    """Vertex amplitude using mpmath arithmetic."""
    total = mpmath.mpf(0)
    for mL1_2 in half_ints_2(jsL2):
        mR1_2 = mjs2 - mL1_2
        if abs(mR1_2) > jsR2:
            continue
        mL2_2 = mL1_2 + mgL2
        if abs(mL2_2) > jtL2:
            continue
        mR2_2 = mR1_2 + mgR2
        if abs(mR2_2) > jtR2:
            continue
        if mL2_2 + mR2_2 != mjt2:
            continue
        c1v = cg_mp(jsL2, jsR2, js2, mL1_2, mR1_2, mjs2)
        if c1v == 0:
            continue
        c2v = cg_mp(jtL2, jtR2, jt2, mL2_2, mR2_2, mjt2)
        if c2v == 0:
            continue
        c3v = cg_mp(jsL2, jgL2, jtL2, mL1_2, mgL2, mL2_2)
        c4v = cg_mp(jsR2, jgR2, jtR2, mR1_2, mgR2, mR2_2)
        total += c1v * c2v * c3v * c4v
    return total


def compute_B_level_mp(n_ext, n_int):
    """Compute B(n_int) = vertex(mj=+1/2) - vertex(mj=-1/2) at mpmath precision."""
    result = mpmath.mpf(0)
    for mj_ext_2 in [+1, -1]:
        sign = 1 if mj_ext_2 == +1 else -1

        jE_L_2 = n_ext + 1
        jE_R_2 = n_ext
        j_ext_2 = 1

        jI_L_2 = n_int + 1
        jI_R_2 = n_int

        lam = mpmath.mpf(2 * n_int + 3) / 2
        lam4 = lam**4

        q_probe = 1
        if not vertex_allowed(n_int, n_int, q_probe):
            continue

        probe_chs = get_channels(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
        if not probe_chs:
            continue

        j_int_min_2 = 1
        j_int_max_2 = 2 * n_int + 1

        subtotal = mpmath.mpf(0)

        for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
            for mj_int_2 in range(-j_int_2, j_int_2 + 1, 2):
                for mj_int_prime_2 in range(-j_int_2, j_int_2 + 1, 2):

                    probe_amp = mpmath.mpf(0)
                    for jpL2, jpR2 in probe_chs:
                        for mpL2 in range(-jpL2, jpL2 + 1, 2):
                            for mpR2 in range(-jpR2, jpR2 + 1, 2):
                                pa = vertex_amp_mp(
                                    jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                    jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                    jpL2, jpR2, mpL2, mpR2)
                                probe_amp += pa

                    if probe_amp == 0:
                        continue

                    q_lo = max(1, abs(n_ext - n_int))
                    q_hi = n_ext + n_int

                    for q_loop in range(q_lo, q_hi + 1):
                        if not vertex_allowed(n_ext, n_int, q_loop):
                            continue
                        mu_q = mpmath.mpf(q_loop * (q_loop + 2))

                        chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
                        chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

                        for jgL1_2, jgR1_2 in chs1:
                            for jgL2_2, jgR2_2 in chs2:
                                for mgL2_v in range(-jgL1_2, jgL1_2 + 1, 2):
                                    for mgR2_v in range(-jgR1_2, jgR1_2 + 1, 2):
                                        if abs(mgL2_v) > jgL2_2 or abs(mgR2_v) > jgR2_2:
                                            continue
                                        v1 = vertex_amp_mp(
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                            jgL1_2, jgR1_2, mgL2_v, mgR2_v)
                                        if v1 == 0:
                                            continue
                                        v3 = vertex_amp_mp(
                                            jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jgL2_2, jgR2_2, mgL2_v, mgR2_v)
                                        if v3 == 0:
                                            continue
                                        subtotal += v1 * probe_amp * v3 / (lam4 * mu_q)

        result += sign * subtotal

    return result


def main():
    # Load existing float B values for comparison
    inv_file = DATA_DIR / "g2_c3_investigation.json"
    with open(inv_file) as f:
        data = json.load(f)
    all_B_float = {int(k): float(v) for k, v in data["all_B"].items()}
    V_mag = mpmath.mpf(data["V_mag"])

    print(f"\nV_mag = {mpmath.nstr(V_mag, 20)}")
    print(f"Schwinger = {mpmath.nstr(SCHWINGER, 20)}")

    # Verify against float values for n_int=1..3
    print("\n--- Verification: mpmath vs float for n_int=1..3 ---")
    for n_int in [1, 2, 3]:
        t0 = time.time()
        B_mp = compute_B_level_mp(1, n_int)
        dt = time.time() - t0
        B_fl = all_B_float.get(n_int, 0)
        diff = float(abs(B_mp - mpmath.mpf(B_fl)))
        print(f"  n_int={n_int}: B_mp={mpmath.nstr(B_mp, 25)}, B_float={B_fl:.15e}, diff={diff:.2e} ({dt:.1f}s)")

    # Full computation for n_int=0..50 at mpmath precision
    print("\n--- Full mpmath computation n_int=0..50 ---")
    all_B_mp = {}
    cum_B_mp = mpmath.mpf(0)

    for n_int in range(0, 51):
        t0 = time.time()
        if n_int == 0:
            B_mp = mpmath.mpf(0)
        else:
            B_mp = compute_B_level_mp(1, n_int)
        dt = time.time() - t0
        all_B_mp[n_int] = B_mp
        cum_B_mp += B_mp

        # Compare with float
        B_fl = all_B_float.get(n_int, 0)
        diff = float(abs(B_mp - mpmath.mpf(B_fl)))

        if n_int <= 5 or n_int % 5 == 0:
            print(f"  n_int={n_int:3d}: B={mpmath.nstr(B_mp, 20):>28s}  "
                  f"diff_vs_float={diff:.2e}  cum={mpmath.nstr(cum_B_mp, 25)}  ({dt:.1f}s)",
                  flush=True)

    # Extract c3 at high precision (no tail correction)
    print("\n--- c3 extraction (partial sum only, no tail) ---")
    F2_S_partial = cum_B_mp / (V_mag * SCHWINGER)
    delta_partial = F2_S_partial - 1
    c1_term = C1 * X
    c2_term = C2 * X**2
    residual_partial = delta_partial - c1_term - c2_term
    c3_partial = residual_partial / X3

    print(f"  cum_B       = {mpmath.nstr(cum_B_mp, 30)}")
    print(f"  F2/S        = {mpmath.nstr(F2_S_partial, 25)}")
    print(f"  delta       = {mpmath.nstr(delta_partial, 25)}")
    print(f"  c1*x        = {mpmath.nstr(c1_term, 25)}")
    print(f"  c2*x^2      = {mpmath.nstr(c2_term, 25)}")
    print(f"  residual    = {mpmath.nstr(residual_partial, 20)}")
    print(f"  c3 (partial)= {mpmath.nstr(c3_partial, 20)}")

    # Save high-precision B values
    out = {
        "dps": 80,
        "n_max": 50,
        "all_B": {str(k): mpmath.nstr(v, 40) for k, v in all_B_mp.items()},
        "cum_B": mpmath.nstr(cum_B_mp, 40),
        "c3_partial": mpmath.nstr(c3_partial, 20),
        "c2": mpmath.nstr(C2, 30),
    }
    outfile = DATA_DIR / "g2_c3_highprec.json"
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {outfile}")


if __name__ == "__main__":
    main()
