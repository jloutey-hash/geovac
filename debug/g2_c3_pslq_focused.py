"""
g2_c3_pslq_focused.py -- Focused algebraic identification of c3.

Strategy: since B(n_int) are algebraic irrationals in Q(sqrt(2), sqrt(3), sqrt(5)),
and PSLQ needs more precision than we have (~7-8 digits of c3_partial),
we use a multi-pronged approach:

1. Verify the structural candidates from previous runs against c3_partial
2. Test whether c3 * (small integer) lies in the T9 ring Q + Q*pi^2 + Q*pi^4
3. Compute exact B(1)..B(5) via sympy Rational CG, improving c3 precision
4. Run PSLQ on the improved c3 value
"""

import json
import sys
import time
from pathlib import Path
from fractions import Fraction
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

# Bilinear combinations
BD = B_HOPF * DELTA_HOPF   # 42/40 = 21/20
FD = F_HOPF * DELTA_HOPF   # pi^2/240
FB = F_HOPF / B_HOPF       # pi^2/252
D = DELTA_HOPF              # 1/40

# n_ext=1 kinematics
LAM = mpmath.mpf(5) / 2
LAM2 = LAM**2
X = 1 / LAM2  # 4/25
X3 = X**3
LAM6 = LAM**6  # 15625/64

pi2 = mpmath.pi**2
pi4 = mpmath.pi**4

# ------------------------------------------------------------------
# Load data and extract c3_partial at maximum float64 precision
# ------------------------------------------------------------------
print("=" * 70)
print("  FOCUSED c3 IDENTIFICATION")
print("=" * 70)

inv_file = DATA_DIR / "g2_c3_investigation.json"
with open(inv_file) as f:
    data = json.load(f)

all_B = {int(k): float(v) for k, v in data["all_B"].items()}
V_mag = mpmath.mpf(data["V_mag"])

# Accumulate cum_B at mpmath precision from float values
cum_B = mpmath.mpf(0)
for n in range(51):
    cum_B += mpmath.mpf(all_B.get(n, 0))

delta = cum_B / (V_mag * SCHWINGER) - 1
c1_term = C1 * X
c2_term = C2 * X**2
residual = delta - c1_term - c2_term
c3_partial = residual / X3

print(f"  c3_partial (n=0..50, no tail) = {mpmath.nstr(c3_partial, 20)}")

# Tail correction (from previous fits)
C_even = mpmath.mpf('2.907587e-02')
p_even = mpmath.mpf('7.6205')
C_odd = mpmath.mpf('7.951405e-02')
p_odd = mpmath.mpf('7.7052')

tail = mpmath.mpf(0)
for n in range(51, 10001):
    if n % 2 == 0:
        tail += C_even * mpmath.power(n, -p_even)
    else:
        tail += C_odd * mpmath.power(n, -p_odd)

tail_unc = abs(tail) * mpmath.mpf('0.226')  # 22.6% from fit uncertainties
tail_c3 = tail / (V_mag * SCHWINGER * X3)
c3_full = c3_partial + tail_c3
c3_unc = tail_unc / (V_mag * SCHWINGER * X3)

print(f"  c3_full    (with tail)        = {mpmath.nstr(c3_full, 20)}")
print(f"  c3 uncertainty                = {mpmath.nstr(c3_unc, 5)}")
print(f"  tail contribution to c3       = {mpmath.nstr(tail_c3, 5)}")

# ------------------------------------------------------------------
# Test structural candidates from previous sessions
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  CANDIDATE VERIFICATION")
print("=" * 70)

candidates = {
    "PSLQ_LAM6": (2 - 7*BD - 4*FD + 5*BD**2 + 2*FD**2 - FB**2) / LAM6,
    "(-8fd+2fb+3fd^2)/413361": (-8*FD + 2*FB + 3*FD**2) / 413361,
    "(-27fd*fb-12d^2)/85819": (-27*FD*FB - 12*D**2) / 85819,
    "(-25fb^2+12bd*d)/(-465656)": (-25*FB**2 + 12*BD*D) / (-465656),
    "(-4d-13fb^2)/201882": (-4*D - 13*FB**2) / 201882,
    "(-4fd-8fb+9d)/425532": (-4*FD - 8*FB + 9*D) / 425532,
    "-c2^2/(2*c1)": -C2**2 / (2*C1),
    "c2*(c2-c1)/(c1^2)": C2*(C2-C1)/(C1**2),
    "c2^3/c1^2": C2**3 / C1**2,
    "(c2^2 - c2/5)/c1": (C2**2 - C2/5)/C1,
}

# Also test simple c2-based combinations
for a in range(-5, 6):
    for b in range(-5, 6):
        for N in [1, 2, 3, 5, 7, 10, 12, 15, 20, 25, 30, 42, 50, 100, 120, 168, 200, 252, 840]:
            if a == 0 and b == 0:
                continue
            val = (a * C2**2 + b * C2) / N
            if abs(val) > 0 and abs(val) < 1e-4:
                err = abs(val - c3_full)
                if err < 3 * c3_unc:
                    name = f"({a}*c2^2 + {b}*c2)/{N}"
                    candidates[name] = val

# Test c1, c2, BD, FD, FB combinations for c3
for a in range(-3, 4):
    for b in range(-3, 4):
        for c in range(-3, 4):
            for N in [1, 5, 7, 10, 25, 42, 100, 252, 840, 5040, 25200]:
                val = (a*BD + b*FD + c*FB) / (N * LAM6)
                if abs(val) > 1e-10 and abs(val) < 1e-4:
                    err = abs(val - c3_full)
                    if err < 2 * c3_unc:
                        name = f"({a}BD+{b}FD+{c}FB)/({N}*L6)"
                        candidates[name] = val

print(f"\n  Testing {len(candidates)} candidates against c3_full:")
print(f"  c3_full = {mpmath.nstr(c3_full, 15)}")
print(f"  c3_unc  = {mpmath.nstr(c3_unc, 5)}")

results = []
for name, val in candidates.items():
    err = abs(val - c3_full)
    rel = float(err / abs(c3_full))
    sigma = float(err / c3_unc)
    results.append((sigma, rel, name, val))

results.sort()
print(f"\n  Top 20 candidates (sorted by sigma):")
for sigma, rel, name, val in results[:20]:
    print(f"    {name:50s}  val={mpmath.nstr(val,12):>18s}  "
          f"rel_err={rel:.3e}  sigma={sigma:.2f}")

# ------------------------------------------------------------------
# Check if c3 is simply related to c2
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  c3 vs c2 RELATIONSHIP")
print("=" * 70)

ratio = c3_full / C2
print(f"  c3/c2 = {mpmath.nstr(ratio, 15)}")
print(f"  c3/c2^2 = {mpmath.nstr(c3_full / C2**2, 15)}")
print(f"  c3*5/c2 = {mpmath.nstr(c3_full * 5 / C2, 15)}")
print(f"  c3/c1 = {mpmath.nstr(c3_full / C1, 15)}")
print(f"  c3/(c1*c2) = {mpmath.nstr(c3_full / (C1*C2), 15)}")
print(f"  c3*LAM6 = {mpmath.nstr(c3_full * LAM6, 15)}")

# Check if c3 = -c2^2 * (something rational/pi)
r = c3_full / (-C2**2)
print(f"\n  c3 / (-c2^2) = {mpmath.nstr(r, 15)}")
print(f"  This * 5 = {mpmath.nstr(r * 5, 15)}")
print(f"  This * 10 = {mpmath.nstr(r * 10, 15)}")
print(f"  This * 25 = {mpmath.nstr(r * 25, 15)}")
print(f"  This * 50 = {mpmath.nstr(r * 50, 15)}")

# Check the pattern: c1 = 1/2, c2 = (2-BD-FD-FB)/5
# Does c3 = (polynomial in BD, FD, FB, D) / N ?
# With LAM6 = (5/2)^6
# c3 * LAM6 = scaled_c3 ~ -0.0001454... this should be a poly in BD, FD, FB
scaled = c3_full * LAM6
print(f"\n  c3 * LAM6 = {mpmath.nstr(scaled, 15)}")
print(f"  c3 * LAM6 / BD = {mpmath.nstr(scaled / BD, 15)}")
print(f"  c3 * LAM6 / FD = {mpmath.nstr(scaled / FD, 15)}")
print(f"  c3 * LAM6 / D^2 = {mpmath.nstr(scaled / D**2, 15)}")

# ------------------------------------------------------------------
# c2 derivation pattern analysis
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  c2 PATTERN ANALYSIS (for extending to c3)")
print("=" * 70)

# c2 = (2 - BD - FD - FB)/5
# Let's verify this is really the formula
c2_check = (2 - BD - FD - FB) / 5
print(f"  c2 = (2 - BD - FD - FB)/5 = {mpmath.nstr(c2_check, 20)}")
print(f"  c2 (from formula)          = {mpmath.nstr(C2, 20)}")
print(f"  Match: {abs(c2_check - C2) < mpmath.mpf('1e-60')}")

# If c3 follows a similar pattern:
# c3 = (a + b*BD + c*FD + d*FB + e*BD^2 + f*FD^2 + g*FB^2 + h*BD*FD + ...) / N

# With c3 ~ -5.946e-7 and the monomials:
print(f"\n  Monomial values:")
monomials = {
    '1': mpmath.mpf(1),
    'BD': BD,
    'FD': FD,
    'FB': FB,
    'D': D,
    'D^2': D**2,
    'BD^2': BD**2,
    'FD^2': FD**2,
    'FB^2': FB**2,
    'BD*FD': BD*FD,
    'BD*FB': BD*FB,
    'FD*FB': FD*FB,
    'BD*D': BD*D,
    'FD*D': FD*D,
    'FB*D': FB*D,
}
for name, val in monomials.items():
    print(f"    {name:10s} = {mpmath.nstr(val, 15)}")

# ------------------------------------------------------------------
# PSLQ with relaxed tolerance on c3_partial
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  PSLQ WITH RELAXED TOLERANCE")
print("=" * 70)

# c3_partial has ~7-8 significant digits. Set dps to 25 for PSLQ.
mpmath.mp.dps = 25

target_partial = c3_partial
target_full = c3_full

# Test 1: c3 in T9 basis [1, pi^2, pi^4]
print("\n  PSLQ: c3 * N = a + b*pi^2 + c*pi^4")
for N in [1, 5, 7, 10, 25, 42, 100, 120, 168, 252, 840, 2520, 5040, 25200]:
    for use_full in [True, False]:
        target = target_full if use_full else target_partial
        label = "full" if use_full else "part"
        scaled = N * target
        vec = [scaled, mpmath.mpf(1), pi2, pi4]
        try:
            rel = mpmath.pslq(vec, maxcoeff=10**5, maxsteps=3000)
            if rel is not None and rel[0] != 0:
                check = sum(r * v for r, v in zip(rel, vec))
                if abs(check) < mpmath.mpf('1e-8'):
                    c3_recon = sum(-r*v for r, v in zip(rel[1:], vec[1:])) / (rel[0] * N)
                    err = abs(c3_recon - target)
                    print(f"  N={N:6d} ({label}): rel={rel}, check={mpmath.nstr(check, 3)}, "
                          f"c3={mpmath.nstr(c3_recon, 10)}, err={mpmath.nstr(err, 3)}")
        except Exception:
            pass

# Test 2: c3 * LAM6 in extended basis
print("\n  PSLQ: c3*LAM6 in [1, pi^2, pi^4, BD, FD, FB]")
mpmath.mp.dps = 25
for use_full in [True, False]:
    target = target_full if use_full else target_partial
    label = "full" if use_full else "part"
    scaled = target * LAM6
    vec = [scaled, mpmath.mpf(1), pi2, pi4, BD, FD, FB]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**4, maxsteps=5000)
        if rel is not None:
            check = sum(r * v for r, v in zip(rel, vec))
            labels = ['c3*L6', '1', 'pi^2', 'pi^4', 'BD', 'FD', 'FB']
            print(f"  ({label}) rel={rel}")
            print(f"    check={mpmath.nstr(check, 5)}")
            for r, lab in zip(rel, labels):
                if r != 0:
                    print(f"    {r:+d}*{lab}")
    except Exception as e:
        print(f"  ({label}) PSLQ error: {e}")

# Test 3: c3*LAM6 in degree-2 basis
print("\n  PSLQ: c3*LAM6 in degree-2 basis")
mpmath.mp.dps = 25
for use_full in [True, False]:
    target = target_full if use_full else target_partial
    label = "full" if use_full else "part"
    scaled = target * LAM6
    vec = [scaled, mpmath.mpf(1), pi2, pi4, BD, FD, FB, BD**2, FD**2, FB**2]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**3, maxsteps=5000)
        if rel is not None:
            check = sum(r * v for r, v in zip(rel, vec))
            labels = ['c3*L6', '1', 'pi^2', 'pi^4', 'BD', 'FD^2', 'FB', 'BD^2', 'FD^2', 'FB^2']
            print(f"  ({label}) rel={rel}")
            print(f"    check={mpmath.nstr(check, 5)}")
    except Exception as e:
        print(f"  ({label}) PSLQ error: {e}")

# Test 4: c3 in terms of c1, c2
print("\n  PSLQ: c3 in [c1, c2, c1^2, c2^2, c1*c2]")
mpmath.mp.dps = 25
for N in [1, 5, 10, 25, 50, 100]:
    for use_full in [True]:
        target = target_full if use_full else target_partial
        label = "full" if use_full else "part"
        vec = [N*target, C1, C2, C1**2, C2**2, C1*C2]
        try:
            rel = mpmath.pslq(vec, maxcoeff=10**4, maxsteps=5000)
            if rel is not None and rel[0] != 0:
                check = sum(r * v for r, v in zip(rel, vec))
                if abs(check) < mpmath.mpf('1e-5'):
                    labels = ['N*c3', 'c1', 'c2', 'c1^2', 'c2^2', 'c1*c2']
                    print(f"  N={N} ({label}): rel={rel}, check={mpmath.nstr(check, 3)}")
        except Exception:
            pass

# ------------------------------------------------------------------
# Direct test of the c2-pattern hypothesis
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  DIRECT HYPOTHESIS TESTS")
print("=" * 70)

mpmath.mp.dps = 80

# Hypothesis: c3 follows c2's pattern but with degree-2 terms
# c2 = (2 - BD - FD - FB)/5
# Parker-Toms: c1 = R/12 where R = 6 (scalar curvature of unit S^3)
# So c1 = 1/2 = R/12

# The general heat-kernel expansion on S^3 gives:
# c_n = sum of products of curvature invariants of total weight 2n
# For S^3: R = d(d-1) = 6, Ric = (d-1)g = 2g, Riem = 2(g^g - g^g) (max sym)
# So all curvature invariants collapse to powers of R = 6

# c1 = R/12 = 1/2  ✓
# c2 = (5R^2 + ...)/720 for spin-1/2 on round S^d (standard result)

# But our c2 involves B, F, Δ — not just R. So the "curvature expansion"
# is not the standard Seeley-DeWitt one; it's specific to the vertex
# correction on the Fock-projected S^3.

# Hypothesis A: c3 = polynomial(BD, FD, FB) / N for small integer N
# We know c3 ~ -5.946e-7, and the monomials are:
# Degree 1: BD~1.05, FD~0.0411, FB~0.0392, D~0.025
# Degree 2: BD^2~1.1025, FD^2~0.00169, FB^2~0.00154, etc.
# Degree 3 would give ~1e-3 for D^3

# For c3 ~ -6e-7, we need N ~ 1e6 if using degree-1 terms
# or N ~ 3000 if using degree-2 terms like FD^2 ~ 1.7e-3

# The PSLQ candidate: c3*LAM6 = 2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2
# LAM6 = 15625/64 ~ 244.14

cand_A = (2 - 7*BD - 4*FD + 5*BD**2 + 2*FD**2 - FB**2) / LAM6
err_A = abs(cand_A - c3_full)
sig_A = float(err_A / c3_unc)
print(f"\n  Candidate A (PSLQ LAM6):")
print(f"    c3 = (2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2) / LAM6")
print(f"    = {mpmath.nstr(cand_A, 15)}")
print(f"    c3_full = {mpmath.nstr(c3_full, 15)}")
print(f"    error = {mpmath.nstr(err_A, 5)}")
print(f"    sigma = {sig_A:.2f}")

# Expand in terms of B, F, Δ
# 2 - 7*(42/40) - 4*(π²/240) + 5*(42/40)² + 2*(π²/240)² - (π²/252)²
# = 2 - 7*21/20 - 4*π²/240 + 5*(21/20)² + 2*(π²/240)² - (π²/252)²
# = 2 - 147/20 - π²/60 + 5*441/400 + 2*π⁴/57600 - π⁴/63504
# = 2 - 7.35 - 0.1645... + 5.5125 + 0.000342... - 0.000155...
# = 0.0 + small terms

# Rational part: 2 - 147/20 + 2205/400 = 2 - 7.35 + 5.5125 = 0.1625
# Pi^2 part: -1/60 = -0.01667...
# Pi^4 part: 2/57600 - 1/63504 = (2*63504 - 57600)/(57600*63504) = (127008-57600)/... = 69408/3657830400

print(f"\n  Rational part: {float(2 - mpmath.mpf(147)/20 + mpmath.mpf(2205)/400)}")
print(f"  = 2 - 147/20 + 2205/400")
rat_part = Fraction(2) - Fraction(147, 20) + Fraction(2205, 400)
print(f"  = {rat_part} = {float(rat_part)}")

pi2_part = Fraction(-1, 60)
print(f"  Pi^2 coefficient: {pi2_part} = {float(pi2_part)}")

pi4_part = Fraction(2, 57600) - Fraction(1, 63504)
print(f"  Pi^4 coefficient: {pi4_part} = {float(pi4_part)}")

# So candidate A = (13/80 - π²/60 + 69408/(57600*63504)*π⁴) / LAM6
# LAM6 = 15625/64
# c3 = (13/80 - π²/60 + ...*π⁴) * 64/15625

c3_A_exact = (rat_part + pi2_part * float(pi2) + pi4_part * float(pi4)) * 64/15625
print(f"\n  c3_A (float check) = {c3_A_exact:.15e}")
print(f"  c3_full (mpmath)   = {float(c3_full):.15e}")

# Hypothesis B: try c2-analog pattern
# c2 = (2 - BD - FD - FB)/5
# Maybe c3 = (a0 + a1*BD + a2*FD + a3*FB) * (b0 + b1*BD + b2*FD + b3*FB) / N
# = c2-like * c2-like / N

# c2 * c1 = (2 - BD - FD - FB)/5 * 1/2 = (2 - BD - FD - FB)/10
c2c1 = C2 * C1
print(f"\n  c1*c2 = {mpmath.nstr(c2c1, 15)}")
print(f"  c3/(c1*c2) = {mpmath.nstr(c3_full/(c2c1), 15)}")

# c2^2 = ((2-BD-FD-FB)/5)^2
c2sq = C2**2
print(f"  c2^2 = {mpmath.nstr(c2sq, 15)}")
print(f"  c3/c2^2 = {mpmath.nstr(c3_full/c2sq, 15)}")

# Is c3 = c2 * (something simple)?
# c3/c2 ~ -3.42e-6. What's that in terms of invariants?
ratio_c3_c2 = c3_full / C2
print(f"\n  c3/c2 = {mpmath.nstr(ratio_c3_c2, 15)}")
print(f"  c3/c2 * 5 = {mpmath.nstr(ratio_c3_c2 * 5, 15)}")
print(f"  c3/c2 * 25 = {mpmath.nstr(ratio_c3_c2 * 25, 15)}")

# c3/c2 * 25 should be ~ -8.55e-5
# D^2 = 1/1600 = 6.25e-4
# FD = pi^2/240 ~ 0.0411
# c3/c2 * 25 / D^2 ~ -0.137
# c3/c2 * 25 / FD^2 ~ -0.0506

# Hypothesis C: Parker-Toms / DeWitt pattern
# On round S^d, the heat-kernel coefficients for spin-1/2 have a known structure
# c_n involves Bernoulli numbers and curvature traces
# For unit S^3 (R=6, K=1):
# c1 = R/12 = 1/2
# c2 = (5R^2 - 2R_{ab}R^{ab} + 2R_{abcd}R^{abcd} - 60*R*xi)/360 for scalar
# For S^3: R_{ab}R^{ab} = R^2/3 = 12, R_{abcd}R^{abcd} = 2R^2/3(3-1) = 12
# So c2_scalar = (180 - 24 + 24)/360 = 1/2
# But our c2 = 0.1739... ≠ 1/2, so this is NOT the standard heat-kernel expansion
# Our c2 involves the VERTEX CORRECTION, not just free propagation

# This is important: the curvature expansion for the vertex correction
# is NOT the standard DeWitt expansion. It's specific to the diagram topology.

print("\n" + "=" * 70)
print("  HIGHER-PRECISION c3 VIA EXACT B(1)")
print("=" * 70)

# Since B(1) = -128*sqrt(3)/91125 + 112*sqrt(2)/50625 is known exactly,
# we can compute its contribution to delta at arbitrary precision.
# B(1) dominates cum_B by >99%, so this dramatically improves precision.

from sympy import sqrt as sym_sqrt, Rational as R, nsimplify, N as sym_N

B1_exact_sympy = -128*sym_sqrt(3)/91125 + 112*sym_sqrt(2)/50625
B1_mp = mpmath.mpf(str(sym_N(B1_exact_sympy, 60)))

# Float vs exact difference
B1_float = mpmath.mpf(all_B[1])
B1_diff = B1_mp - B1_float
print(f"  B(1) exact (60 digits) = {mpmath.nstr(B1_mp, 30)}")
print(f"  B(1) float             = {mpmath.nstr(B1_float, 20)}")
print(f"  Difference             = {mpmath.nstr(B1_diff, 5)}")

# The difference tells us how much precision we gain
# B(1) ~ 6.96e-4, float error ~ B1_diff
# This propagates to delta as: delta_correction = B1_diff / (V_mag * SCHWINGER)
delta_correction = B1_diff / (V_mag * SCHWINGER)
print(f"  Delta correction       = {mpmath.nstr(delta_correction, 5)}")
print(f"  c3 correction          = {mpmath.nstr(delta_correction / X3, 5)}")

# Use exact B(1) and float B(2..50)
cum_B_improved = B1_mp  # exact B(1)
for n in range(2, 51):
    cum_B_improved += mpmath.mpf(all_B.get(n, 0))

delta_improved = cum_B_improved / (V_mag * SCHWINGER) - 1
residual_improved = delta_improved - c1_term - c2_term
c3_partial_improved = residual_improved / X3

print(f"\n  c3_partial (improved with exact B(1)) = {mpmath.nstr(c3_partial_improved, 20)}")
print(f"  c3_partial (original float)           = {mpmath.nstr(c3_partial, 20)}")
print(f"  Difference                            = {mpmath.nstr(c3_partial_improved - c3_partial, 5)}")

# Now add exact B(2) if available
B2_exact_text = "-68/496125 + 52*sqrt(30)/52093125 + 1252*sqrt(3)/72930375 + 116*sqrt(10)/3472875"
from sympy import sqrt, Rational
B2_exact_sympy = -Rational(68, 496125) + 52*sqrt(30)/52093125 + 1252*sqrt(3)/72930375 + 116*sqrt(10)/3472875
B2_mp = mpmath.mpf(str(sym_N(B2_exact_sympy, 60)))
B2_float = mpmath.mpf(all_B[2])
B2_diff = B2_mp - B2_float
print(f"\n  B(2) exact (60 digits) = {mpmath.nstr(B2_mp, 30)}")
print(f"  B(2) float             = {mpmath.nstr(B2_float, 20)}")
print(f"  B(2) difference        = {mpmath.nstr(B2_diff, 5)}")

# Use exact B(1) and B(2), float B(3..50)
cum_B_improved2 = B1_mp + B2_mp
for n in range(3, 51):
    cum_B_improved2 += mpmath.mpf(all_B.get(n, 0))

delta_improved2 = cum_B_improved2 / (V_mag * SCHWINGER) - 1
residual_improved2 = delta_improved2 - c1_term - c2_term
c3_partial_improved2 = residual_improved2 / X3

print(f"\n  c3_partial (exact B(1,2))  = {mpmath.nstr(c3_partial_improved2, 20)}")
print(f"  c3_partial (exact B(1))    = {mpmath.nstr(c3_partial_improved, 20)}")
print(f"  Difference                 = {mpmath.nstr(c3_partial_improved2 - c3_partial_improved, 5)}")

# ------------------------------------------------------------------
# V_mag precision check
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  V_mag PRECISION ANALYSIS")
print("=" * 70)

# V_mag enters as: delta = cum_B / (V_mag * Schwinger) - 1
# Relative error in V_mag → relative error in (1+delta) → absolute error in delta
# V_mag is computed as the probe vertex magnitude (a CG sum at n_int=n_ext=1)
# It should also be algebraic

# What is V_mag exactly?
# V_mag = |sum of vertex amplitudes for probe at n_int=n_ext=1|
# This is a sum of products of CG coefficients, so algebraic

# The stored V_mag has 15 digits. If we could compute it exactly,
# that would improve c3 by removing the V_mag precision limitation.

print(f"  V_mag = {mpmath.nstr(V_mag, 20)}")
print(f"  V_mag^2 = {mpmath.nstr(V_mag**2, 20)}")
print(f"  Schwinger = {mpmath.nstr(SCHWINGER, 20)}")
print(f"  V_mag * Schwinger = {mpmath.nstr(V_mag * SCHWINGER, 20)}")

# The V_mag precision limits c3 to about:
# delta ~ 0.0845, V_mag relative error ~ 1e-15
# error in delta ~ 0.0845 * 1e-15 ~ 8.5e-17
# error in residual ~ 8.5e-17
# error in c3 ~ 8.5e-17 / 4.1e-3 ~ 2.1e-14
# c3 ~ 6e-7, relative error ~ 3.5e-8 → about 7.5 significant digits

# Even with exact B values, V_mag at float64 limits c3 to ~7-8 digits.
# To break past this, we need V_mag at higher precision.
# V_mag is the probe normalization, computable from CG coefficients at n_int=1.

# But actually — ALPHA also has limited precision (10 digits)!
# Error in alpha ~ alpha * 1.5e-10
# Error in Schwinger ~ Schwinger * 1.5e-10
# Error in delta ~ delta * 1.5e-10 ~ 1.3e-11
# Error in residual ~ 1.3e-11
# Error in c3 ~ 1.3e-11 / 4.1e-3 ~ 3.2e-9
# That's about 2.2 significant digits of c3!

# This is the DOMINANT error source. Alpha's 10-digit precision
# limits c3 to ~2 significant digits.

# WAIT — that can't be right. The PREVIOUS extraction got c3 to 4-5 digits.
# Let me reconsider...

# Actually, the issue is more subtle. The vertex correction is computed as
# F2 = sum(B(n)) / V_mag, and then divided by Schwinger = alpha/(2*pi).
# The alpha factor cancels between the numerator (which includes alpha in the
# coupling) and the denominator (Schwinger = alpha/(2*pi)).

# In the computation, each B(n) is proportional to alpha (from the coupling),
# and we divide by Schwinger = alpha/(2*pi). So the alpha dependence cancels:
# F2/Schwinger = (alpha * vertex_sum) / (V_mag * alpha/(2*pi)) = 2*pi * vertex_sum / V_mag

# So actually F2/Schwinger does NOT depend on alpha! It's purely geometric.
# Let me verify...

# From g2_c3_fast.py:
# V_mag is computed with probe B(n_int=n_ext=1, q=1)
# F2 = (sum_n B(n)) / V_mag
# Schwinger = alpha / (2*pi)
# ratio = F2 / Schwinger = sum_n B(n) / (V_mag * alpha/(2*pi))

# But B(n) as computed in the code does NOT include an explicit alpha factor.
# It's purely a CG-weighted sum with propagator factors 1/(lam^4 * mu_q).
# And V_mag is also purely CG.

# So the alpha enters ONLY through Schwinger in the denominator.
# F2/Schwinger = sum_n B(n) / (V_mag * alpha/(2*pi))
# = 2*pi * sum_n B(n) / (V_mag * alpha)

# This DOES depend on alpha. The alpha doesn't cancel because B(n) doesn't
# include alpha — the code computes the BARE vertex sum at O(alpha^0),
# and the physical vertex correction is alpha times this.

# So the precision of alpha matters. At 10 digits, we lose:
# delta ~ 0.085, error ~ 0.085 * 1.5e-10 ~ 1.3e-11
# residual ~ -2.5e-9, error ~ 1.3e-11
# c3 = residual/x^3, error ~ 1.3e-11/4.1e-3 ~ 3.2e-9
# c3 ~ 6e-7, relative error ~ 5.3e-3 → about 2.3 digits

# This means ALPHA IS THE BOTTLENECK. We can't get more than ~2-3 digits
# of c3 without using a higher-precision alpha value.

# CODATA 2018: alpha = 7.2973525693(11) × 10^-3
# CODATA 2022: alpha = 7.29735256937(28) × 10^-3 (more precise but not yet in code)

# Let's use the CODATA 2022 value with full precision
ALPHA_HP = mpmath.mpf('7.2973525693e-3')  # only 10 digits available
# Actually, the CODATA recommended value is:
# 1/alpha = 137.035999084(21)
# alpha = 1/137.035999084 = 0.0072973525693...

# We can get more digits from 1/alpha:
alpha_inv = mpmath.mpf('137.035999084')
ALPHA_FROM_INV = 1 / alpha_inv
print(f"\n  alpha from direct    = {mpmath.nstr(ALPHA, 20)}")
print(f"  alpha from 1/alpha   = {mpmath.nstr(ALPHA_FROM_INV, 20)}")
print(f"  difference           = {mpmath.nstr(ALPHA - ALPHA_FROM_INV, 5)}")

# Actually 1/alpha = 137.035999084(21) has 12 significant digits
# vs alpha = 7.2973525693(11)e-3 has 10 significant digits

# Use 1/alpha for higher precision
SCHWINGER_HP = ALPHA_FROM_INV / (2 * mpmath.pi)

# Recompute c3 with higher-precision Schwinger
delta_hp = cum_B_improved2 / (V_mag * SCHWINGER_HP) - 1
residual_hp = delta_hp - c1_term - c2_term
c3_hp = residual_hp / X3

print(f"\n  c3 (hp alpha, exact B1+B2) = {mpmath.nstr(c3_hp, 20)}")
print(f"  c3 (original alpha)        = {mpmath.nstr(c3_partial_improved2, 20)}")
print(f"  Difference                 = {mpmath.nstr(c3_hp - c3_partial_improved2, 5)}")

# Now the precision is limited by V_mag (15 digits) and B(3..50) (15 digits)
# Effective c3 precision should be ~7-8 significant digits

# Add tail correction
c3_hp_full = c3_hp + tail_c3
print(f"\n  c3_hp_full (with tail) = {mpmath.nstr(c3_hp_full, 15)}")

# ------------------------------------------------------------------
# Final PSLQ with best precision
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  FINAL PSLQ WITH BEST AVAILABLE PRECISION")
print("=" * 70)

mpmath.mp.dps = 30
target = c3_hp_full

print(f"  Target c3 = {mpmath.nstr(target, 15)}")

# Test in pure T9 basis [1, pi^2, pi^4] at various scalings
print(f"\n  c3 * N in [1, pi^2, pi^4]:")
for N in [1, 5, 7, 10, 12, 25, 42, 120, 252, 840, 2520, 5040, 25200, 50400, 85819]:
    scaled = N * target
    vec = [scaled, mpmath.mpf(1), pi2, pi4]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**6, maxsteps=5000)
        if rel is not None and rel[0] != 0:
            check = sum(r * v for r, v in zip(rel, vec))
            if abs(check) < abs(target) * mpmath.mpf('0.01'):  # within 1% relative
                c3_r = -sum(r*v for r, v in zip(rel[1:], vec[1:])) / (rel[0] * N)
                err = abs(c3_r - target)
                rel_e = float(err / abs(target))
                print(f"  N={N:6d}: rel={rel}, check={mpmath.nstr(check, 3)}, "
                      f"rel_err={rel_e:.2e}")
    except Exception:
        pass

# Test c3*LAM6 in extended basis
print(f"\n  c3*LAM6 in [1, pi^2, pi^4, BD, FD, FB, BD^2, FD^2, FB^2, D^2]:")
scaled = target * LAM6
basis = [scaled, mpmath.mpf(1), pi2, pi4, BD, FD, FB, BD**2, FD**2, FB**2, D**2]
labels = ['c3*L6', '1', 'pi2', 'pi4', 'BD', 'FD', 'FB', 'BD^2', 'FD^2', 'FB^2', 'D^2']
try:
    rel = mpmath.pslq(basis, maxcoeff=500, maxsteps=10000)
    if rel is not None:
        check = sum(r * v for r, v in zip(rel, basis))
        print(f"  rel={rel}")
        print(f"  check={mpmath.nstr(check, 5)}")
        for r, lab in zip(rel, labels):
            if r != 0:
                print(f"    {r:+d}*{lab}")
except Exception as e:
    print(f"  Error: {e}")

# Test c3 in [c1, c2, c1*c2, c2^2, c1^2*c2]
print(f"\n  c3 in [c2, c2^2, c1*c2, c1^2]:")
for N in [1, 5, 10, 25, 50, 100, 252]:
    vec_c = [N*target, C2, C2**2, C1*C2, C1**2]
    try:
        rel = mpmath.pslq(vec_c, maxcoeff=10**4, maxsteps=5000)
        if rel is not None and rel[0] != 0:
            check = sum(r * v for r, v in zip(rel, vec_c))
            if abs(check) < abs(target) * mpmath.mpf('0.1'):
                print(f"  N={N}: rel={rel}, check={mpmath.nstr(check, 3)}")
    except Exception:
        pass

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
print(f"  c3_partial (no tail, float B)     = {mpmath.nstr(c3_partial, 12)}")
print(f"  c3_partial (exact B1+B2)          = {mpmath.nstr(c3_partial_improved2, 12)}")
print(f"  c3_hp (hp alpha, exact B1+B2)     = {mpmath.nstr(c3_hp, 12)}")
print(f"  c3_hp_full (with tail)            = {mpmath.nstr(c3_hp_full, 12)}")
print(f"  c3 uncertainty                    = {mpmath.nstr(c3_unc, 5)}")
print(f"\n  Candidate A: (2-7BD-4FD+5BD^2+2FD^2-FB^2)/LAM6")
print(f"    = {mpmath.nstr(cand_A, 12)}")
print(f"    sigma from c3_hp_full = {float(abs(cand_A - c3_hp_full)/c3_unc):.2f}")
print(f"\n  B(n_int) are algebraic irrationals in Q(sqrt(2), sqrt(3), sqrt(5))")
print(f"  Alpha precision (10 digits) is the dominant bottleneck for c3 extraction")
