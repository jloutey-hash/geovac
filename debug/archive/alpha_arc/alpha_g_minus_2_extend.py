"""
Extend the B(n_ext, n_int=1) computation to n_ext=11..15 and
refine the asymptotic analysis.

Also: look at the product CG_sum * (n+1)*(n+2) more carefully.
The slope is converging (Richardson -> 0.2397 at n=10).
If q_inf = 12/50 = 6/25 = 0.24, then c0 = 12*q_inf/625 = 12*0.24/625 = 0.004608.
But the Shanks^3 gives ~0.00489. Let's see.

Actually the relation c0 = 12*q/625 was derived assuming V_mag ~ 4/(3n)
and mu_q ~ n^2. Let me be more precise.

V_mag(n) = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]
         = (2/3) * [2/n - 3/n^2 + 10/(3n^3) - ...] (from sympy series)
         = 4/(3n) - 2/n^2 + 20/(9n^3) - ...

mu_q = n*(n+2) = n^2 + 2n

prop = (5/2)^4 * n*(n+2) = 625/16 * (n^2+2n)

lam = n + 3/2

F2 = B / V_mag = CG_sum / (prop * V_mag)

lam^2 * F2 = (n+3/2)^2 * CG_sum / (625/16 * (n^2+2n) * V_mag)

For large n:
lam^2 ~ n^2
prop ~ 625/16 * n^2
V_mag ~ 4/(3n)
So lam^2*F2 ~ n^2 * CG_sum / (625/16 * n^2 * 4/(3n))
= n^2 * CG_sum * 3n / (625/16 * n^2 * 4)
= CG_sum * 3n * 16 / (625 * 4)
= CG_sum * 12n / 625

So c0 = lim_{n->inf} CG_sum * 12n / 625
     = 12/625 * lim_{n->inf} n * CG_sum(n)

And CG_sum*(n+1)*(n+2) ~ q*n => CG_sum ~ q*n/n^2 = q/n
=> n*CG_sum -> q
=> c0 = 12*q/625

But we need to be more careful about subleading terms.
"""

from sympy import (Rational, sqrt, simplify, S, sympify)
from sympy.physics.wigner import clebsch_gordan
import time


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def V_magnetic_closed(n):
    return Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR: continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL: continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR: continue
        if mL2 + mR2 != mj_t: continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0: continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0: continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


def compute_B_nint1(n_ext):
    """Compute B(n_ext, n_int=1) exactly."""
    n_int = 1
    j_ext = Rational(1, 2)
    jEL = Rational(n_ext + 1, 2)
    jER = Rational(n_ext, 2)
    jIL = Rational(1)
    jIR = Rational(1, 2)

    lam_int = Rational(5, 2)
    q = n_ext
    mu_q = Rational(q * (q + 2))
    prop = lam_int**4 * mu_q

    loop_channels = []
    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                      (Rational(q-1, 2), Rational(q+1, 2))]:
        if jgL < 0 or jgR < 0: continue
        if not (abs(jEL - jgL) <= jIL <= jEL + jgL): continue
        if not (abs(jER - jgR) <= jIR <= jER + jgR): continue
        loop_channels.append((jgL, jgR))

    probe_channels = []
    for jpL, jpR in [(Rational(1), Rational(0)), (Rational(0), Rational(1))]:
        if not (abs(jIL - jpL) <= jIL <= jIL + jpL): continue
        if not (abs(jIR - jpR) <= jIR <= jIR + jpR): continue
        probe_channels.append((jpL, jpR))

    j_int_vals = [Rational(1, 2), Rational(3, 2)]

    B_up = S.Zero
    B_dn = S.Zero

    for jgL, jgR in loop_channels:
        for mgL in half_ints(jgL):
            for mgR in half_ints(jgR):
                for jpL, jpR in probe_channels:
                    for mpL in half_ints(jpL):
                        for mpR in half_ints(jpR):
                            for j_int in j_int_vals:
                                for m_int in half_ints(j_int):
                                    for m_int_p in half_ints(j_int):
                                        v1u = vertex_amp_pol(
                                            jEL, jER, j_ext, Rational(1, 2),
                                            jIL, jIR, j_int, m_int,
                                            jgL, jgR, mgL, mgR)
                                        if v1u == 0: continue
                                        vp = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int,
                                            jIL, jIR, j_int, m_int_p,
                                            jpL, jpR, mpL, mpR)
                                        if vp == 0: continue
                                        v2u = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int_p,
                                            jEL, jER, j_ext, Rational(1, 2),
                                            jgL, jgR, mgL, mgR)
                                        B_up += v1u * vp * v2u

                                        v1d = vertex_amp_pol(
                                            jEL, jER, j_ext, Rational(-1, 2),
                                            jIL, jIR, j_int, m_int,
                                            jgL, jgR, mgL, mgR)
                                        if v1d == 0: continue
                                        v2d = vertex_amp_pol(
                                            jIL, jIR, j_int, m_int_p,
                                            jEL, jER, j_ext, Rational(-1, 2),
                                            jgL, jgR, mgL, mgR)
                                        B_dn += v1d * vp * v2d

    return simplify((B_up - B_dn) / prop)


# First reproduce n=1..10, then extend to 11..15
all_lam2_f2 = []

# Load cached results from the previous run
cached = {
    1: 0.023588460396, 2: 0.013757413244, 3: 0.010569894823,
    4: 0.009010317729, 5: 0.008089889121, 6: 0.007484099716,
    7: 0.007055800827, 8: 0.006737245018, 9: 0.006491188403,
    10: 0.006295488045,
}

for n in range(1, 11):
    lam2_f2 = cached[n]
    all_lam2_f2.append(lam2_f2)
    print(f"n={n:2d}: lam^2*F2 = {lam2_f2:.12f} (cached)")

# Compute n=11..15
for n in range(11, 16):
    t0 = time.time()
    B = compute_B_nint1(n)
    dt = time.time() - t0

    V_mag = V_magnetic_closed(n)
    F2 = simplify(B / V_mag)
    lam = Rational(2*n + 3, 2)
    lam2_F2 = float(lam**2 * F2)

    all_lam2_f2.append(lam2_F2)
    print(f"n={n:2d}: lam^2*F2 = {lam2_F2:.12f} ({dt:.1f}s)")


# Shanks analysis with 15 points
def shanks(seq):
    n = len(seq)
    if n < 3: return seq[-1:]
    results = []
    for i in range(n - 2):
        denom = (seq[i+2] - seq[i+1]) - (seq[i+1] - seq[i])
        if abs(denom) < 1e-30:
            results.append(seq[i+2])
        else:
            results.append(seq[i+2] - (seq[i+2] - seq[i+1])**2 / denom)
    return results


print("\n\n=== Shanks transforms (15 points) ===")
s1 = shanks(all_lam2_f2)
print("S1:", [f"{v:.10f}" for v in s1[-5:]])

s2 = shanks(s1) if len(s1) >= 3 else []
if s2:
    print("S2:", [f"{v:.10f}" for v in s2[-5:]])

s3 = shanks(s2) if len(s2) >= 3 else []
if s3:
    print("S3:", [f"{v:.10f}" for v in s3[-5:]])

s4 = shanks(s3) if len(s3) >= 3 else []
if s4:
    print("S4:", [f"{v:.10f}" for v in s4[-5:]])

s5 = shanks(s4) if len(s4) >= 3 else []
if s5:
    print("S5:", [f"{v:.10f}" for v in s5[-5:]])

# Best estimate
for label, seq in [("S5", s5), ("S4", s4), ("S3", s3), ("S2", s2), ("S1", s1)]:
    if seq:
        best = seq[-1]
        print(f"\nBest from {label}: c0 = {best:.15f}")
        print(f"1/c0 = {1/best:.10f}")
        print(f"c0 * 625/12 = {best*625/12:.10f}")
        break


# Also try Richardson on the last few points assuming error ~ a/n + b/n^2
print("\n\n=== Three-point Richardson (error ~ a/n + b/n^2) ===")
for i in range(2, len(all_lam2_f2)):
    n1 = i - 1
    n2 = i
    n3 = i + 1
    T1 = all_lam2_f2[i-2]
    T2 = all_lam2_f2[i-1]
    T3 = all_lam2_f2[i]
    # T = c0 + a/n + b/n^2
    # Three equations, three unknowns: c0, a, b
    # n1*n2*n3 * [T1*(n2-n3) + T2*(n3-n1) + T3*(n1-n2)] = some polynomial
    # Actually let's solve the 3x3 system
    # T1 = c0 + a/n1 + b/n1^2
    # T2 = c0 + a/n2 + b/n2^2
    # T3 = c0 + a/n3 + b/n3^2
    import numpy as np
    A = np.array([[1, 1./n1, 1./n1**2],
                  [1, 1./n2, 1./n2**2],
                  [1, 1./n3, 1./n3**2]])
    rhs = np.array([T1, T2, T3])
    sol = np.linalg.solve(A, rhs)
    c0_est = sol[0]
    if i >= len(all_lam2_f2) - 5:
        print(f"  n=({n1},{n2},{n3}): c0 = {c0_est:.12f}, a = {sol[1]:.6f}, b = {sol[2]:.6f}")


# Try candidates for c0
print("\n\n=== c0 candidate identification ===")
# Get best estimate
best_c0 = None
for seq in [s5, s4, s3, s2, s1]:
    if seq:
        best_c0 = seq[-1]
        break

import math
tests = [
    ("1/200", 1/200),
    ("1/204", 1/204),
    ("1/210", 1/210),
    ("1/216", 1/216),
    ("1/220", 1/220),
    ("1/225", 1/225),
    ("1/(6*pi^2)", 1/(6*math.pi**2)),
    ("1/(18*pi)", 1/(18*math.pi)),
    ("pi/645", math.pi/645),
    ("pi/648", math.pi/648),
    ("pi/650", math.pi/650),
    ("2/(9*pi^2)", 2/(9*math.pi**2)),
    ("4/(81*pi)", 4/(81*math.pi)),
    ("1/(3*pi*sqrt(7))", 1/(3*math.pi*math.sqrt(7))),
    ("48/(625*9)", 48/5625),
    ("48*6/(625*25)", 48*6/15625),
    ("12*6/(625*25)", 12*6/15625),
    ("1/(4*pi*4)", 1/(16*math.pi)),
    ("2*sqrt(2)/(9*pi*5)", 2*math.sqrt(2)/(45*math.pi)),
]
for name, val in tests:
    rel = (val - best_c0) / best_c0 if best_c0 else 0
    print(f"  {name:25s} = {val:.12f}, rel err = {rel:+.4e}")
