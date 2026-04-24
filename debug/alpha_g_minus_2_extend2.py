"""
Extend B(n_ext, n_int=1) to n_ext=16..25 and perform high-precision
Richardson/Neville extrapolation for c0 = lim lambda^2 * F2.

Also: attempt to identify c0 algebraically using PSLQ or nsimplify.

From 15-point data:
  8-point Richardson gives c0 ~ 0.004590077
  1/c0 ~ 217.86

Strategy: compute more points to sharpen the extrapolation, then
test a wider set of algebraic candidates.
"""

from sympy import (Rational, sqrt, simplify, S, nsimplify)
from sympy.physics.wigner import clebsch_gordan
import time
import json
import os


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


# Cached lam^2*F2 values from n=1..15 (verified correct, both j_int sectors)
cached_lam2_f2 = {
    1: 0.023588460396, 2: 0.013757413244, 3: 0.010569894823,
    4: 0.009010317729, 5: 0.008089889121, 6: 0.007484099716,
    7: 0.007055800827, 8: 0.006737245018, 9: 0.006491188403,
    10: 0.006295488045,
    11: 0.006136064889, 12: 0.006004176817, 13: 0.005893369437,
    14: 0.005799060929, 15: 0.005717830319,
}

all_lam2_f2 = []
all_cg_sum = []

# Load cached
for n in range(1, 16):
    all_lam2_f2.append(cached_lam2_f2[n])

# Compute n=16..25
print("=== Computing B(n_ext, n_int=1) for n_ext=16..25 ===")
new_data = {}

for n in range(16, 26):
    t0 = time.time()
    B = compute_B_nint1(n)
    dt = time.time() - t0

    V_mag = V_magnetic_closed(n)
    F2 = simplify(B / V_mag)
    lam = Rational(2*n + 3, 2)
    lam2_F2 = float(lam**2 * F2)

    prop = Rational(625, 16) * n * (n + 2)
    cg_sum = float(simplify(B * prop))

    all_lam2_f2.append(lam2_F2)
    new_data[n] = {'lam2_f2': lam2_F2, 'cg_sum': cg_sum, 'B': float(B)}

    print(f"n={n:2d}: lam^2*F2 = {lam2_F2:.12f}, CG_sum = {cg_sum:.12f} ({dt:.1f}s)")


# Save all data
output = {
    'lam2_f2': {str(i+1): v for i, v in enumerate(all_lam2_f2)},
    'new_data': {str(k): v for k, v in new_data.items()},
}
out_path = os.path.join(os.path.dirname(__file__), 'data', 'alpha_g_minus_2_extend2.json')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nData saved to {out_path}")


# Richardson/Neville extrapolation with all 25 points
print("\n\n=== Multi-point Richardson extrapolation (25 points) ===")
import numpy as np

N = len(all_lam2_f2)

# Three-point Richardson assuming T(n) = c0 + a/n + b/n^2
print("\n--- Three-point Richardson (error ~ a/n + b/n^2) ---")
for i in range(N - 2):
    n1 = i + 1
    n2 = i + 2
    n3 = i + 3
    T1 = all_lam2_f2[i]
    T2 = all_lam2_f2[i+1]
    T3 = all_lam2_f2[i+2]
    A = np.array([[1, 1./n1, 1./n1**2],
                  [1, 1./n2, 1./n2**2],
                  [1, 1./n3, 1./n3**2]])
    rhs = np.array([T1, T2, T3])
    sol = np.linalg.solve(A, rhs)
    if i >= N - 8:
        print(f"  n=({n1},{n2},{n3}): c0 = {sol[0]:.15f}")


# Four-point Richardson assuming T(n) = c0 + a/n + b/n^2 + c/n^3
print("\n--- Four-point Richardson (error ~ a/n + b/n^2 + c/n^3) ---")
for i in range(N - 3):
    n1 = i + 1
    n2 = i + 2
    n3 = i + 3
    n4 = i + 4
    T = [all_lam2_f2[i+k] for k in range(4)]
    ns = [n1, n2, n3, n4]
    A = np.array([[1, 1./n, 1./n**2, 1./n**3] for n in ns])
    sol = np.linalg.solve(A, T)
    if i >= N - 8:
        print(f"  n=({n1},{n2},{n3},{n4}): c0 = {sol[0]:.15f}")


# Five-point Richardson
print("\n--- Five-point Richardson (error ~ sum a_k/n^k, k=1..4) ---")
for i in range(N - 4):
    ns = [i + 1 + k for k in range(5)]
    T = [all_lam2_f2[i + k] for k in range(5)]
    A = np.array([[1] + [1./n**k for k in range(1, 5)] for n in ns])
    sol = np.linalg.solve(A, T)
    if i >= N - 6:
        print(f"  n=({ns[0]}..{ns[-1]}): c0 = {sol[0]:.15f}")


# Six-point Richardson
print("\n--- Six-point Richardson (error ~ sum a_k/n^k, k=1..5) ---")
for i in range(N - 5):
    ns = [i + 1 + k for k in range(6)]
    T = [all_lam2_f2[i + k] for k in range(6)]
    A = np.array([[1] + [1./n**k for k in range(1, 6)] for n in ns])
    sol = np.linalg.solve(A, T)
    if i >= N - 5:
        print(f"  n=({ns[0]}..{ns[-1]}): c0 = {sol[0]:.15f}")


# Seven-point Richardson
print("\n--- Seven-point Richardson (error ~ sum a_k/n^k, k=1..6) ---")
for i in range(N - 6):
    ns = [i + 1 + k for k in range(7)]
    T = [all_lam2_f2[i + k] for k in range(7)]
    A = np.array([[1] + [1./n**k for k in range(1, 7)] for n in ns])
    sol = np.linalg.solve(A, T)
    if i >= N - 4:
        print(f"  n=({ns[0]}..{ns[-1]}): c0 = {sol[0]:.15f}")


# Eight-point Richardson
print("\n--- Eight-point Richardson (error ~ sum a_k/n^k, k=1..7) ---")
for i in range(N - 7):
    ns = [i + 1 + k for k in range(8)]
    T = [all_lam2_f2[i + k] for k in range(8)]
    A = np.array([[1] + [1./n**k for k in range(1, 8)] for n in ns])
    sol = np.linalg.solve(A, T)
    if i >= N - 3:
        print(f"  n=({ns[0]}..{ns[-1]}): c0 = {sol[0]:.15f}")


# Use the BEST estimate from the highest-order Richardson at latest data
# and try to identify c0 algebraically
print("\n\n=== Algebraic identification of c0 ===")
# Compute best estimate from 8-point Richardson at latest window
ns = list(range(N - 7, N + 1))
T = all_lam2_f2[ns[0]-1:]
A_mat = np.array([[1] + [1./n**k for k in range(1, 8)] for n in ns])
best_sol = np.linalg.solve(A_mat, T)
c0_best = best_sol[0]
print(f"Best 8-point Richardson c0 = {c0_best:.15f}")
print(f"1/c0 = {1/c0_best:.12f}")
print(f"c0 * 625 = {c0_best * 625:.12f}")
print(f"c0 * 625/12 = {c0_best * 625/12:.12f}")
print(f"c0 * 2500 = {c0_best * 2500:.12f}")
print(f"c0 * 2500/48 = {c0_best * 2500/48:.12f}")

import math
print(f"\nc0 / (alpha/(2*pi)) = {c0_best / (1/137.036 / (2*math.pi)):.10f}")
print(f"c0 * 2*pi/alpha = {c0_best * 2 * math.pi * 137.036:.10f}")

# Extended candidate testing
print("\n--- Extended candidate testing ---")
alpha_em = 1/137.035999084
candidates = [
    ("1/218", 1/218),
    ("1/217", 1/217),
    ("1/219", 1/219),
    ("1/216", 1/216),
    ("1/220", 1/220),
    ("2/435", 2/435),
    ("2/437", 2/437),
    ("3/653", 3/653),
    ("3/654", 3/654),
    ("4/872", 4/872),
    ("4/871", 4/871),
    ("1/(6*pi^2)", 1/(6*math.pi**2)),
    ("1/(2*pi*sqrt(12))", 1/(2*math.pi*math.sqrt(12))),
    ("1/(4*pi*sqrt(3))", 1/(4*math.pi*math.sqrt(3))),
    ("pi/(4*625)", math.pi/2500),
    ("12/(625*pi)", 12/(625*math.pi)),
    ("48/(625*pi^2)", 48/(625*math.pi**2)),
    ("12/(625*pi^2)*pi", 12*math.pi/(625*math.pi**2)),
    ("alpha/(2*pi)", alpha_em/(2*math.pi)),
    ("2*alpha/pi", 2*alpha_em/math.pi),
    ("4*alpha/(2*pi)", 4*alpha_em/(2*math.pi)),
    ("alpha^2", alpha_em**2),
    ("alpha/(2*pi)*8/3", alpha_em/(2*math.pi)*8/3),
    ("12*0.24/625", 12*0.24/625),
    ("12*6/(625*25)", 12*6/15625),
    ("48/(625*pi*sqrt(3))", 48/(625*math.pi*math.sqrt(3))),
    ("1/(25*sqrt(pi)*sqrt(3))", 1/(25*math.sqrt(math.pi)*math.sqrt(3))),
    ("12/(5^5)", 12/3125),
    ("12/(5^4*5.0155)", 12/(625*5.0155)),
    ("2/(3*sqrt(5)*5^3)", 2/(3*math.sqrt(5)*125)),
    ("4/(3*5^4)", 4/(3*625)),
    ("8/(3*5^4)", 8/(3*625)),
    ("16/(3*5^5)", 16/(3*3125)),
    ("4/(5^4*sqrt(3))", 4/(625*math.sqrt(3))),
    ("8/(5^4*sqrt(7))", 8/(625*math.sqrt(7))),
    ("4/(3*5^3*sqrt(5))", 4/(3*125*math.sqrt(5))),
    ("48/(5^4*21)", 48/(625*21)),
    ("48/(5^4*20)", 48/12500),
    ("48/(625*20.944)", 48/(625*20.944)),
    ("48/(625*21)", 48/13125),
    ("48/(625*pi*sqrt(3)*4/3)", 48*3/(625*math.pi*math.sqrt(3)*4)),
]

# Sort by absolute relative error
results = []
for name, val in candidates:
    rel_err = abs(val - c0_best) / c0_best
    results.append((rel_err, name, val))
results.sort()

print(f"c0_best = {c0_best:.15f}")
print(f"\nTop 15 closest candidates:")
for rel_err, name, val in results[:15]:
    print(f"  {name:30s} = {val:.15f}, rel err = {rel_err:.6e}")


# nsimplify attempt
print("\n\n=== nsimplify on c0, 1/c0, c0*625, c0*2500/48 ===")
for label, val in [("c0", c0_best),
                   ("1/c0", 1/c0_best),
                   ("c0*625", c0_best*625),
                   ("c0*2500/48", c0_best*2500/48),
                   ("c0*625/12", c0_best*625/12),
                   ("c0*625*12", c0_best*625*12)]:
    try:
        result = nsimplify(val, tolerance=1e-4, rational=False)
        print(f"  {label:20s} = {val:.12f} -> nsimplify = {result} ({float(result):.12f})")
    except Exception as e:
        print(f"  {label:20s} = {val:.12f} -> nsimplify failed: {e}")


# Shanks transform on the full 25-point sequence
print("\n\n=== Shanks transforms (25 points) ===")
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

s1 = shanks(all_lam2_f2)
print(f"S1 (last 5): {[f'{v:.12f}' for v in s1[-5:]]}")
s2 = shanks(s1) if len(s1) >= 3 else []
if s2: print(f"S2 (last 5): {[f'{v:.12f}' for v in s2[-5:]]}")
s3 = shanks(s2) if len(s2) >= 3 else []
if s3: print(f"S3 (last 5): {[f'{v:.12f}' for v in s3[-5:]]}")
s4 = shanks(s3) if len(s3) >= 3 else []
if s4: print(f"S4 (last 5): {[f'{v:.12f}' for v in s4[-5:]]}")
s5 = shanks(s4) if len(s4) >= 3 else []
if s5: print(f"S5 (last 5): {[f'{v:.12f}' for v in s5[-5:]]}")
s6 = shanks(s5) if len(s5) >= 3 else []
if s6: print(f"S6 (last 5): {[f'{v:.12f}' for v in s6[-5:]]}")
s7 = shanks(s6) if len(s6) >= 3 else []
if s7: print(f"S7 (last 5): {[f'{v:.12f}' for v in s7[-5:]]}")

for label, seq in [("S7", s7), ("S6", s6), ("S5", s5), ("S4", s4), ("S3", s3)]:
    if seq:
        best = seq[-1]
        print(f"\nBest from {label}: c0 = {best:.15f}")
        print(f"  1/c0 = {1/best:.12f}")
        print(f"  c0*625/12 = {best*625/12:.12f}")
        break


# CG_sum analysis with full data
print("\n\n=== CG_sum * (n+1)*(n+2) analysis ===")
# Recompute CG_sum for all n
for n in range(1, N+1):
    lam = n + 1.5
    lam2_f2 = all_lam2_f2[n-1]
    V_mag = float(V_magnetic_closed(n))
    F2 = lam2_f2 / lam**2
    B = F2 * V_mag
    prop = 625.0/16.0 * n * (n+2)
    cg = B * prop
    cg_prod = cg * (n+1) * (n+2)
    if n >= N - 10:
        print(f"  n={n:2d}: CG*(n+1)(n+2) = {cg_prod:.12f}")
