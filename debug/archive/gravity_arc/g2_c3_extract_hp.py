"""
g2_c3_extract_hp.py -- High-precision c3 extraction from existing float64 B values.

Key insight: each B(n_int) is computed from exact sympy CG and stored as float64
(15-16 digits). The precision loss happens in the EXTRACTION step where we subtract
c1*x + c2*x^2 from the accumulated delta. By doing extraction at 80 dps, we
preserve all available digits.

Additionally: compute exact B(n_int) for n_int=1..5 using full sympy Rational
to verify rationality and measure float64 error.
"""

import json
import sys
import time
from pathlib import Path
from fractions import Fraction

import mpmath
mpmath.mp.dps = 80

sys.path.insert(0, '.')
DATA_DIR = Path(__file__).parent / "data"

# ------------------------------------------------------------------
# Paper 2 invariants (exact)
# ------------------------------------------------------------------
B = mpmath.mpf(42)
F = mpmath.pi**2 / 6
Delta = mpmath.mpf(1) / 40
ALPHA = mpmath.mpf('7.2973525693e-3')
SCHWINGER = ALPHA / (2 * mpmath.pi)

BD = B * Delta        # 21/20
FD = F * Delta        # pi^2/240
FB = F / B            # pi^2/252
D  = Delta            # 1/40

C1 = mpmath.mpf(1) / 2
C2 = (2 - BD - FD - FB) / 5

LAM = mpmath.mpf(5) / 2
LAM2 = LAM**2
X = 1 / LAM2
X3 = X**3
LAM6 = LAM**6

print("=" * 70)
print("  HIGH-PRECISION c3 EXTRACTION FROM EXISTING DATA")
print("=" * 70)
print(f"  c2 = {mpmath.nstr(C2, 25)}")
print(f"  x  = 4/25 = {mpmath.nstr(X, 20)}")
print(f"  x^3 = {mpmath.nstr(X3, 20)}")

# ------------------------------------------------------------------
# Load existing float64 B values
# ------------------------------------------------------------------
inv_file = DATA_DIR / "g2_c3_investigation.json"
with open(inv_file) as f:
    data = json.load(f)
V_mag = mpmath.mpf(str(data["V_mag"]))
all_B_str = data["all_B"]

print(f"\n  V_mag = {mpmath.nstr(V_mag, 20)}")
print(f"  Loaded {len(all_B_str)} B values")

# Convert to mpmath (preserving all float64 digits)
all_B_mp = {}
for k, v in all_B_str.items():
    n = int(k)
    all_B_mp[n] = mpmath.mpf(str(v))

# ------------------------------------------------------------------
# Step 1: High-precision extraction from cumulative float data
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 1: HP EXTRACTION (cumulative float B, mpmath arithmetic)")
print("=" * 70)

cum_B = mpmath.mpf(0)
for n in range(0, 51):
    cum_B += all_B_mp[n]

F2_S = cum_B / (V_mag * SCHWINGER)
delta = F2_S - 1
c1_term = C1 * X
c2_term = C2 * X**2
residual = delta - c1_term - c2_term
c3_partial = residual / X3

print(f"  cum_B (n=0..50)  = {mpmath.nstr(cum_B, 30)}")
print(f"  F2/Schwinger     = {mpmath.nstr(F2_S, 25)}")
print(f"  delta            = {mpmath.nstr(delta, 25)}")
print(f"  c1*x             = {mpmath.nstr(c1_term, 25)}")
print(f"  c2*x^2           = {mpmath.nstr(c2_term, 25)}")
print(f"  residual (no tail) = {mpmath.nstr(residual, 20)}")
print(f"  c3 (partial, no tail) = {mpmath.nstr(c3_partial, 15)}")

# Convergence study: c3 at each truncation
print(f"\n  Convergence of c3_partial (no tail):")
for n_cut in [10, 15, 20, 25, 30, 35, 40, 45, 50]:
    cb = sum(all_B_mp[n] for n in range(0, n_cut + 1))
    f2s = cb / (V_mag * SCHWINGER)
    d = f2s - 1
    r = d - c1_term - c2_term
    c3p = r / X3
    print(f"    n={n_cut:3d}: c3_partial = {mpmath.nstr(c3p, 15)}")

# ------------------------------------------------------------------
# Step 2: Tail correction with power-law fit
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 2: TAIL CORRECTION")
print("=" * 70)

import numpy as np

n_fit_min = 25
even_n = [n for n in range(n_fit_min, 51) if n % 2 == 0 and n > 0]
odd_n = [n for n in range(n_fit_min, 51) if n % 2 == 1]

def logfit(ns):
    logns = np.log(np.array(ns, dtype=float))
    logBs = np.log(np.abs(np.array([float(all_B_mp[n]) for n in ns])))
    A = np.column_stack([np.ones_like(logns), -logns])
    x, _, _, _ = np.linalg.lstsq(A, logBs, rcond=None)
    return np.exp(x[0]), x[1]

C_e, p_e = logfit(even_n)
C_o, p_o = logfit(odd_n)
print(f"  Even: C={C_e:.6e}, p={p_e:.4f}")
print(f"  Odd:  C={C_o:.6e}, p={p_o:.4f}")

def tail_mp(n_start, Ce, pe, Co, po):
    n_even = n_start if n_start % 2 == 0 else n_start + 1
    n_odd = n_start if n_start % 2 == 1 else n_start + 1
    te = mpmath.mpf(Ce) * mpmath.power(2, -mpmath.mpf(pe)) * mpmath.zeta(mpmath.mpf(pe), mpmath.mpf(n_even) / 2)
    to = mpmath.mpf(Co) * mpmath.power(2, -mpmath.mpf(po)) * mpmath.zeta(mpmath.mpf(po), mpmath.mpf(n_odd) / 2)
    return te + to

tail = tail_mp(51, C_e, p_e, C_o, p_o)
total_B = cum_B + tail
F2_S_full = total_B / (V_mag * SCHWINGER)
delta_full = F2_S_full - 1
residual_full = delta_full - c1_term - c2_term
c3_full = residual_full / X3

# Uncertainty from power-law fit
dp = 0.05
tail_hi = tail_mp(51, C_e, p_e - dp, C_o, p_o - dp)
tail_lo = tail_mp(51, C_e, p_e + dp, C_o, p_o + dp)
tail_unc = max(abs(tail_hi - tail), abs(tail_lo - tail))
c3_unc = tail_unc / (V_mag * SCHWINGER * X3)

print(f"\n  Tail (n>50)      = {mpmath.nstr(tail, 15)}")
print(f"  Tail unc         = {mpmath.nstr(tail_unc, 5)}")
print(f"  c3 (with tail)   = {mpmath.nstr(c3_full, 15)}")
print(f"  c3 uncertainty   = {mpmath.nstr(c3_unc, 5)}")
print(f"  c3 significance  = {float(abs(c3_full)/c3_unc):.1f} sigma")

# ------------------------------------------------------------------
# Step 3: Test specific algebraic hypotheses for c3
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 3: ALGEBRAIC HYPOTHESIS TESTING")
print("=" * 70)

c3_target = c3_full
c3_partial_target = c3_partial

# The PSLQ candidate from previous run: c3*LAM6 = 2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2
def formula_pslq(coeff_dict):
    """Evaluate a formula c3*LAM6 = sum(a_i * mono_i)."""
    monomials = {
        '1': mpmath.mpf(1),
        'BD': BD, 'FD': FD, 'FB': FB, 'D': D,
        'BD^2': BD**2, 'FD^2': FD**2, 'FB^2': FB**2, 'D^2': D**2,
        'BD*D': BD*D, 'FD*D': FD*D, 'FB*D': FB*D,
        'BD*FD': BD*FD, 'BD*FB': BD*FB, 'FD*FB': FD*FB,
        'c2': C2, 'c2^2': C2**2,
    }
    val = mpmath.mpf(0)
    for name, coeff in coeff_dict.items():
        val += mpmath.mpf(coeff) * monomials[name]
    return val / LAM6

# Candidate from PSLQ run
cand_pslq = formula_pslq({'1': 2, 'BD': -7, 'FD': -4, 'BD^2': 5, 'FD^2': 2, 'FB^2': -1})
err_pslq = abs((cand_pslq - c3_target) / c3_target)
print(f"\n  PSLQ: c3*L6 = 2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2")
print(f"    value = {mpmath.nstr(cand_pslq, 15)}")
print(f"    target = {mpmath.nstr(c3_target, 15)}")
print(f"    rel_err = {float(err_pslq):.3e}")
print(f"    abs_err = {mpmath.nstr(abs(cand_pslq - c3_target), 5)}")

# Test: can the error be absorbed into the tail?
pslq_abs_err = abs(cand_pslq - c3_target)
print(f"    error/unc = {float(pslq_abs_err / c3_unc):.2f}")

# Systematic search: c3 * N = polynomial(BD, FD, FB, D) for small N
print(f"\n  Systematic search: c3 * N = poly(BD, FD, FB, D)")
print(f"  Testing N = 1..15 (matching c2 pattern where N=5)")

# Separate T9: c3 = r0 + r2*pi^2 + r4*pi^4
# BD = 21/20 (rational)
# FD = pi^2/240 (pi^2*rational)
# FB = pi^2/252 (pi^2*rational)
# D = 1/40 (rational)
# BD^2 = 441/400 (rational)
# FD^2 = pi^4/57600 (pi^4*rational)
# FB^2 = pi^4/63504 (pi^4*rational)
# BD*FD = 21*pi^2/4800 (pi^2*rational)
# BD*FB = 21*pi^2/5040 = pi^2/240 (pi^2*rational) -- same as FD!
# FD*FB = pi^4/(240*252) = pi^4/60480 (pi^4*rational)
# etc.

# Note: BD*FB = (21/20)(pi^2/252) = 21*pi^2/5040 = pi^2/240 = FD !!
# This means BD*FB = FD, which constrains the search.

print(f"\n  Structural identity: BD*FB = {mpmath.nstr(BD*FB, 20)}")
print(f"                      FD    = {mpmath.nstr(FD, 20)}")
print(f"  BD*FB == FD: {abs(BD*FB - FD) < mpmath.mpf('1e-60')}")

# Define the monomials with their T9 decomposition
# Rational monomials:
rat_monos = {
    '1': mpmath.mpf(1),
    'BD': mpmath.mpf(21) / 20,
    'D': mpmath.mpf(1) / 40,
    'BD^2': (mpmath.mpf(21) / 20)**2,
    'D^2': (mpmath.mpf(1) / 40)**2,
    'BD*D': mpmath.mpf(21) / 800,
}
# pi^2 monomials (coefficient of pi^2):
pi2_monos = {
    'FD': mpmath.mpf(1) / 240,
    'FB': mpmath.mpf(1) / 252,
    'BD*FD': mpmath.mpf(21) / (20 * 240),
    'BD*FB': mpmath.mpf(21) / (20 * 252),
    'FD*D': mpmath.mpf(1) / (240 * 40),
    'FB*D': mpmath.mpf(1) / (252 * 40),
}
# pi^4 monomials (coefficient of pi^4):
pi4_monos = {
    'FD^2': mpmath.mpf(1) / (240**2),
    'FB^2': mpmath.mpf(1) / (252**2),
    'FD*FB': mpmath.mpf(1) / (240 * 252),
}

print(f"\n  Rational monomials: {list(rat_monos.keys())}")
print(f"  pi^2 monomials: {list(pi2_monos.keys())}")
print(f"  pi^4 monomials: {list(pi4_monos.keys())}")

# For a formula c3 = (sum a_i*mono_i) / N, decompose into T9:
# c3 = r0 + r2*pi^2 + r4*pi^4
# where r0, r2, r4 are rational

# At this precision, c3 partial (no tail) should give 7-8 digits
# c3 full (with tail) gives 2-3 digits due to tail uncertainty

# Use c3_partial for structural analysis (more precise, known offset from true c3)
print(f"\n  Using c3_partial (no tail) for structural search:")
print(f"    c3_partial = {mpmath.nstr(c3_partial, 15)}")
print(f"    c3_full    = {mpmath.nstr(c3_full, 15)}")
print(f"    difference = {mpmath.nstr(c3_full - c3_partial, 5)}")

# ------------------------------------------------------------------
# Step 4: Enhanced PSLQ on c3_partial (higher relative precision)
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 4: PSLQ ON c3_partial (no tail)")
print("=" * 70)

# The value c3_partial has ~7-8 significant digits of precision
# This is much better than c3_full (~2-3 digits)
# But c3_partial != c3_full; it's missing the tail

# PSLQ approach: try to identify c3_partial in the T9 ring
# c3_partial should also be in the T9 ring (it's a finite sum of rational
# quantities from CG products, so it's actually rational!)

# Wait - the PARTIAL sum (no tail) IS an exact quantity. It's the
# sum of B(n)/V_mag/Schwinger for n=0..50, minus 1 - c1*x - c2*x^2.
# Each B(n) is a specific algebraic number. If they're all rational,
# then c3_partial is of the form q + r*pi^2 + s*pi^4 where the pi
# dependence comes only from c2 (which has pi^2).

# Actually: cum_B is a sum of numbers from CG products (probably algebraic).
# V_mag involves Schwinger = alpha/(2*pi), so there's a 1/pi factor.
# delta = cum_B/(V_mag * alpha/(2pi)) - 1
# So delta involves cum_B * 2*pi / (V_mag * alpha)
# This introduces a pi factor...

# Let's check what V_mag is
print(f"\n  V_mag = {mpmath.nstr(V_mag, 20)}")
print(f"  1/(V_mag*Schwinger) = {mpmath.nstr(1/(V_mag*SCHWINGER), 20)}")

# V_mag = Schwinger * sum_probe  (from the normalization)
# Actually V_mag is the probe vertex magnitude, which is also a CG sum

# The issue is that alpha and pi enter through V_mag*Schwinger = V_mag*alpha/(2pi)
# So c3 = (cum_B_contribution) * 2*pi/(V_mag*alpha) / x^3 - c1/x^2 - c2/x - 1/x^3

# This is getting complex. Let me just try PSLQ directly on c3_partial.

# For PSLQ, I'll multiply c3_partial by various small integers N and check
# if N*c3_partial = rational + rational*pi^2 + rational*pi^4

# Scale c3_partial to O(1)
target = c3_partial
pi2 = mpmath.pi**2
pi4 = mpmath.pi**4

# Try PSLQ at increasing norms
print(f"\n  PSLQ tests on c3_partial:")
for N in [1, 5, 7, 10, 12, 15, 20, 25, 42, 100, 252, 840, 5040, 25200, 50400]:
    scaled = N * target
    vec = [scaled, mpmath.mpf(1), pi2, pi4]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**6, maxsteps=5000)
        if rel is not None and rel[0] != 0:
            # c3_partial = (-rel[1] - rel[2]*pi^2 - rel[3]*pi^4) / (rel[0] * N)
            check = sum(r * v for r, v in zip(rel, vec))
            if abs(check) < mpmath.mpf('1e-20'):
                r0 = mpmath.mpf(-rel[1]) / (rel[0] * N)
                r2 = mpmath.mpf(-rel[2]) / (rel[0] * N)
                r4 = mpmath.mpf(-rel[3]) / (rel[0] * N)
                val = r0 + r2 * pi2 + r4 * pi4
                err = abs(val - target)
                print(f"  N={N:6d}: rel={rel}, check={mpmath.nstr(check, 5)}")
                print(f"    c3 = {mpmath.nstr(r0,8)} + {mpmath.nstr(r2,8)}*pi^2 + {mpmath.nstr(r4,8)}*pi^4")
                print(f"    err = {mpmath.nstr(err, 5)}")
    except Exception:
        pass

# ------------------------------------------------------------------
# Step 5: Test the c2-pattern generalization
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 5: c2-PATTERN GENERALIZATION")
print("=" * 70)

# c2 = (2 - BD - FD - FB)/5 works.
# What about c3 = (a0 + a1*BD + a2*FD + a3*FB + a4*D +
#                  a5*BD^2 + a6*FD^2 + a7*FB^2 + a8*D^2 +
#                  a9*BD*FD + a10*BD*FB + ...) / N  ?

# Since c3 ~ 6e-7 and the monomials are O(1) to O(1e-4),
# we need massive cancellation. This means |a_i| must be large
# or N must be very large.

# For c2, the coefficients were all ±1 or 2, and N=5.
# For c3, if coefficients are O(10) and N=O(10^6), that's plausible.

# But a brute-force search over 15 integer coefficients is infeasible.
# Let me instead check if c3 * LAM6 or c3 * (some framework number)
# has a clean PSLQ decomposition.

# Test: c3 * various scalings
scalings = {
    'LAM6': LAM6,
    'LAM6*5': LAM6 * 5,
    'LAM6*7': LAM6 * 7,
    'LAM6*42': LAM6 * 42,
    '25200': mpmath.mpf(25200),  # denominator from c2
    '5040': mpmath.mpf(5040),    # 7!
    '2520': mpmath.mpf(2520),    # LCM(1..10)
    '840': mpmath.mpf(840),      # LCM(1..8)
    '40*LAM6': 40 * LAM6,
    '42*LAM6': 42 * LAM6,
}

print(f"\n  Testing scaled PSLQ: c3*S in basis [1, pi^2, pi^4, BD, FD, FB, D]")
basis_names = ['1', 'pi^2', 'pi^4', 'BD', 'FD', 'FB', 'D', 'BD^2', 'FD^2', 'FB^2', 'D^2']
basis_vals = [
    mpmath.mpf(1), pi2, pi4, BD, FD, FB, D,
    BD**2, FD**2, FB**2, D**2
]

for name, S in scalings.items():
    scaled = c3_partial * S
    vec = [scaled] + basis_vals
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**4, maxsteps=5000)
        if rel is not None and rel[0] != 0:
            check = sum(r * v for r, v in zip(rel, vec))
            if abs(check) < mpmath.mpf('1e-15'):
                print(f"\n  S={name}: PSLQ hit!")
                print(f"    relation: {rel}")
                for i, (r, n) in enumerate(zip(rel, ['c3*S'] + basis_names)):
                    if r != 0:
                        print(f"      {r:+d} * {n}")
                print(f"    check: {mpmath.nstr(check, 5)}")
    except Exception:
        pass

# ------------------------------------------------------------------
# Step 6: Exact computation for n_int=1 using sympy Rational
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  STEP 6: EXACT B(n_int=1) VIA SYMPY RATIONAL")
print("=" * 70)

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, sqrt, nsimplify, S

def vertex_amp_exact(jsL2, jsR2, js2, mjs2,
                     jtL2, jtR2, jt2, mjt2,
                     jgL2, jgR2, mgL2, mgR2):
    """Exact sympy vertex amplitude."""
    total = S.Zero
    for mL1_2 in range(-jsL2, jsL2 + 1, 2):
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
        c1v = clebsch_gordan(Rational(jsL2,2), Rational(jsR2,2), Rational(js2,2),
                              Rational(mL1_2,2), Rational(mR1_2,2), Rational(mjs2,2))
        if c1v == 0:
            continue
        c2v = clebsch_gordan(Rational(jtL2,2), Rational(jtR2,2), Rational(jt2,2),
                              Rational(mL2_2,2), Rational(mR2_2,2), Rational(mjt2,2))
        if c2v == 0:
            continue
        c3v = clebsch_gordan(Rational(jsL2,2), Rational(jgL2,2), Rational(jtL2,2),
                              Rational(mL1_2,2), Rational(mgL2,2), Rational(mL2_2,2))
        c4v = clebsch_gordan(Rational(jsR2,2), Rational(jgR2,2), Rational(jtR2,2),
                              Rational(mR1_2,2), Rational(mgR2,2), Rational(mR2_2,2))
        total += c1v * c2v * c3v * c4v
    return total

def vertex_allowed_s(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1

def get_channels_s(n_src, n_tgt, q, jsL2, jsR2, jtL2, jtR2):
    if not vertex_allowed_s(n_src, n_tgt, q):
        return []
    chs = []
    for jgL2, jgR2 in [(q + 1, q - 1), (q - 1, q + 1)]:
        if jgL2 < 0 or jgR2 < 0:
            continue
        if (abs(jsL2 - jgL2) <= jtL2 <= jsL2 + jgL2 and
                abs(jsR2 - jgR2) <= jtR2 <= jsR2 + jgR2):
            chs.append((jgL2, jgR2))
    return chs

def compute_B_exact(n_ext, n_int):
    """Compute B(n_int) using exact sympy arithmetic."""
    from sympy import S as Sym
    result = Sym.Zero

    jE_L_2 = n_ext + 1
    jE_R_2 = n_ext
    j_ext_2 = 1

    jI_L_2 = n_int + 1
    jI_R_2 = n_int

    lam = Rational(2 * n_int + 3, 2)
    lam4 = lam**4

    q_probe = 1
    if not vertex_allowed_s(n_int, n_int, q_probe):
        return Sym.Zero

    probe_chs = get_channels_s(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
    if not probe_chs:
        return Sym.Zero

    j_int_min_2 = 1
    j_int_max_2 = 2 * n_int + 1

    for mj_ext_2 in [+1, -1]:
        sign = Sym.One if mj_ext_2 == +1 else -Sym.One

        for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
            for mj_int_2 in range(-j_int_2, j_int_2 + 1, 2):
                for mj_int_prime_2 in range(-j_int_2, j_int_2 + 1, 2):

                    probe_amp = Sym.Zero
                    for jpL2, jpR2 in probe_chs:
                        for mpL2 in range(-jpL2, jpL2 + 1, 2):
                            for mpR2 in range(-jpR2, jpR2 + 1, 2):
                                pa = vertex_amp_exact(
                                    jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                    jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                    jpL2, jpR2, mpL2, mpR2)
                                probe_amp += pa

                    if probe_amp == 0:
                        continue

                    q_lo = max(1, abs(n_ext - n_int))
                    q_hi = n_ext + n_int

                    for q_loop in range(q_lo, q_hi + 1):
                        if not vertex_allowed_s(n_ext, n_int, q_loop):
                            continue
                        mu_q = q_loop * (q_loop + 2)

                        chs1 = get_channels_s(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
                        chs2 = get_channels_s(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

                        for jgL1_2, jgR1_2 in chs1:
                            for jgL2_2, jgR2_2 in chs2:
                                for mgL2_v in range(-jgL1_2, jgL1_2 + 1, 2):
                                    for mgR2_v in range(-jgR1_2, jgR1_2 + 1, 2):
                                        if abs(mgL2_v) > jgL2_2 or abs(mgR2_v) > jgR2_2:
                                            continue
                                        v1 = vertex_amp_exact(
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                            jgL1_2, jgR1_2, mgL2_v, mgR2_v)
                                        if v1 == 0:
                                            continue
                                        v3 = vertex_amp_exact(
                                            jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jgL2_2, jgR2_2, mgL2_v, mgR2_v)
                                        if v3 == 0:
                                            continue
                                        result += sign * v1 * probe_amp * v3 / (lam4 * mu_q)

    return result

# Compute B(1) exactly
t0 = time.time()
B1_exact = compute_B_exact(1, 1)
dt = time.time() - t0
B1_simplified = B1_exact.simplify()
print(f"  B(1) exact = {B1_simplified}")
print(f"  B(1) float = {float(B1_simplified):.15e}")
print(f"  B(1) stored = {float(all_B_mp[1]):.15e}")
print(f"  Difference = {abs(float(B1_simplified) - float(all_B_mp[1])):.2e}")
print(f"  Is rational? {B1_simplified.is_rational}")
print(f"  Time: {dt:.1f}s")

# Compute B(2) exactly
t0 = time.time()
B2_exact = compute_B_exact(1, 2)
dt = time.time() - t0
B2_simplified = B2_exact.simplify()
print(f"\n  B(2) exact = {B2_simplified}")
print(f"  B(2) float = {float(B2_simplified):.15e}")
print(f"  B(2) stored = {float(all_B_mp[2]):.15e}")
print(f"  Difference = {abs(float(B2_simplified) - float(all_B_mp[2])):.2e}")
print(f"  Is rational? {B2_simplified.is_rational}")
print(f"  Time: {dt:.1f}s")


if __name__ == "__main__":
    main()
