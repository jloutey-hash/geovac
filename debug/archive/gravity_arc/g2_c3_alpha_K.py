"""
g2_c3_alpha_K.py -- Extract c3 using alpha_K = 1/K at 80 dps.

Key insight: the curvature coefficients c_k are alpha-independent geometric
quantities of S³. Alpha enters only through the overall Schwinger normalization.
Using alpha_K = 1/(pi(B+F-Delta)) at 80 dps eliminates the 10-digit CODATA
bottleneck, giving ~15 digits of c3 (limited by the float64 B values).
"""

import json
import sys
from pathlib import Path

import mpmath
mpmath.mp.dps = 80

sys.path.insert(0, '.')

DATA_DIR = Path(__file__).parent / "data"

# Paper 2 invariants at high precision
B_HOPF = mpmath.mpf(42)
F_HOPF = mpmath.pi**2 / 6
DELTA_HOPF = mpmath.mpf(1) / 40

# Alpha from Paper 2: K = pi(B + F - Delta)
K = mpmath.pi * (B_HOPF + F_HOPF - DELTA_HOPF)
ALPHA_K = 1 / K
SCHWINGER_K = ALPHA_K / (2 * mpmath.pi)

# Known curvature coefficients
C1 = mpmath.mpf(1) / 2
C2 = (2 - B_HOPF * DELTA_HOPF - F_HOPF * DELTA_HOPF - F_HOPF / B_HOPF) / 5

# n_ext=1 kinematics
LAM = mpmath.mpf(5) / 2
LAM2 = LAM**2
X = 1 / LAM2  # 4/25
X3 = X**3
LAM6 = LAM**6  # 15625/64

# Bilinear combinations
BD = B_HOPF * DELTA_HOPF        # 21/20
FD = F_HOPF * DELTA_HOPF        # pi^2/240
FB = F_HOPF / B_HOPF            # pi^2/252
D = DELTA_HOPF                  # 1/40
D2 = D**2                       # 1/1600

print("=" * 70)
print("  c3 EXTRACTION WITH alpha_K (Paper 2 self-consistent)")
print("=" * 70)
print(f"  K = pi(B + F - Delta) = {mpmath.nstr(K, 30)}")
print(f"  alpha_K = 1/K          = {mpmath.nstr(ALPHA_K, 25)}")
print(f"  Schwinger_K            = {mpmath.nstr(SCHWINGER_K, 25)}")
print(f"  c1 = {C1}")
print(f"  c2 = {mpmath.nstr(C2, 25)}")
print(f"  x  = 4/25 = {mpmath.nstr(X, 20)}")

# CODATA alpha for comparison
ALPHA_CODATA = mpmath.mpf('7.2973525693e-3')
print(f"\n  alpha_K      = {mpmath.nstr(ALPHA_K, 20)}")
print(f"  alpha_CODATA = {mpmath.nstr(ALPHA_CODATA, 20)}")
print(f"  difference   = {mpmath.nstr(ALPHA_K - ALPHA_CODATA, 5)}")
print(f"  rel diff     = {mpmath.nstr((ALPHA_K - ALPHA_CODATA)/ALPHA_K, 5)}")

# Load B values
inv_file = DATA_DIR / "g2_c3_investigation.json"
with open(inv_file) as f:
    data = json.load(f)
all_B_float = {int(k): float(v) for k, v in data["all_B"].items()}
V_mag_float = float(data["V_mag"])
V_mag = mpmath.mpf(data["V_mag"])

print(f"\n  V_mag = {mpmath.nstr(V_mag, 20)}")

# Cumulative B sum
cum_B = mpmath.mpf(0)
for n in range(51):
    cum_B += mpmath.mpf(all_B_float.get(n, 0))

print(f"  cum_B (n=0..50) = {mpmath.nstr(cum_B, 25)}")

# Tail correction (power-law fit from previous analysis)
p_even = mpmath.mpf('7.5229')
p_odd = mpmath.mpf('7.6252')
B_50 = mpmath.mpf(all_B_float[50])
B_49 = mpmath.mpf(all_B_float[49])
tail_even = B_50 * mpmath.mpf(50)**(p_even) / (p_even - 1) / mpmath.mpf(50)**(p_even - 1)
tail_odd = B_49 * mpmath.mpf(49)**(p_odd) / (p_odd - 1) / mpmath.mpf(49)**(p_odd - 1)
tail = tail_even + tail_odd
# Uncertainty estimate
tail_unc = abs(tail) * mpmath.mpf('0.2')

print(f"  tail correction       = {mpmath.nstr(tail, 10)}")
print(f"  tail uncertainty      = {mpmath.nstr(tail_unc, 5)}")

# ===================================================================
# EXTRACTION with alpha_K
# ===================================================================
print("\n" + "=" * 70)
print("  STEP 1: c3 extraction (no tail)")
print("=" * 70)

F2_S = cum_B / (V_mag * SCHWINGER_K)
delta = F2_S - 1
c1_term = C1 * X
c2_term = C2 * X**2
residual = delta - c1_term - c2_term
c3_partial = residual / X3

print(f"  F2/S        = {mpmath.nstr(F2_S, 25)}")
print(f"  delta       = {mpmath.nstr(delta, 20)}")
print(f"  c1*x        = {mpmath.nstr(c1_term, 20)}")
print(f"  c2*x^2      = {mpmath.nstr(c2_term, 20)}")
print(f"  residual    = {mpmath.nstr(residual, 15)}")
print(f"  c3_partial  = {mpmath.nstr(c3_partial, 15)}")

print("\n" + "=" * 70)
print("  STEP 2: c3 extraction (with tail)")
print("=" * 70)

cum_B_full = cum_B + tail
F2_S_full = cum_B_full / (V_mag * SCHWINGER_K)
delta_full = F2_S_full - 1
residual_full = delta_full - c1_term - c2_term
c3_full = residual_full / X3
c3_unc = (tail_unc / (V_mag * SCHWINGER_K)) / X3

print(f"  cum_B_full  = {mpmath.nstr(cum_B_full, 25)}")
print(f"  delta_full  = {mpmath.nstr(delta_full, 20)}")
print(f"  residual    = {mpmath.nstr(residual_full, 15)}")
print(f"  c3_full     = {mpmath.nstr(c3_full, 15)}")
print(f"  c3_unc      = {mpmath.nstr(c3_unc, 5)}")

# Alpha uncertainty contribution
dalpha = mpmath.mpf('2.1e-11')  # CODATA uncertainty on 1/alpha
c3_alpha_unc = abs(cum_B / (V_mag * SCHWINGER_K**2) * dalpha / (2 * mpmath.pi * K**2)) / X3
print(f"  c3_alpha_unc (CODATA) = {mpmath.nstr(c3_alpha_unc, 5)}")
print(f"  c3_alpha_unc (K)      = 0 (exact)")

# ===================================================================
# COMPARISON: alpha_K vs alpha_CODATA
# ===================================================================
print("\n" + "=" * 70)
print("  STEP 3: Compare alpha_K vs alpha_CODATA extraction")
print("=" * 70)

SCHWINGER_CODATA = ALPHA_CODATA / (2 * mpmath.pi)
F2_S_codata = cum_B_full / (V_mag * SCHWINGER_CODATA)
delta_codata = F2_S_codata - 1
residual_codata = delta_codata - c1_term - c2_term
c3_codata = residual_codata / X3

print(f"  c3 (alpha_K)      = {mpmath.nstr(c3_full, 15)}")
print(f"  c3 (alpha_CODATA) = {mpmath.nstr(c3_codata, 15)}")
print(f"  difference        = {mpmath.nstr(c3_full - c3_codata, 5)}")
print(f"  The c_k are alpha-independent, but the extraction precision differs.")

# Also with the 12-digit 1/alpha
alpha_hp = 1 / mpmath.mpf('137.035999084')
schwinger_hp = alpha_hp / (2 * mpmath.pi)
F2_S_hp = cum_B_full / (V_mag * schwinger_hp)
delta_hp = F2_S_hp - 1
residual_hp = delta_hp - c1_term - c2_term
c3_hp = residual_hp / X3

print(f"  c3 (1/alpha=137.035999084) = {mpmath.nstr(c3_hp, 15)}")
print(f"  diff K vs hp               = {mpmath.nstr(c3_full - c3_hp, 5)}")

# ===================================================================
# PSLQ IDENTIFICATION
# ===================================================================
print("\n" + "=" * 70)
print("  STEP 4: PSLQ identification of c3")
print("=" * 70)

target = c3_full
target_partial = c3_partial
print(f"  c3_full    = {mpmath.nstr(target, 20)}")
print(f"  c3_partial = {mpmath.nstr(target_partial, 20)}")

# Test Candidate A
cand_A = (2 - 7*BD - 4*FD + 5*BD**2 + 2*FD**2 - FB**2) / LAM6
print(f"\n  Candidate A = {mpmath.nstr(cand_A, 20)}")
print(f"  c3_full     = {mpmath.nstr(target, 20)}")
print(f"  diff        = {mpmath.nstr(target - cand_A, 10)}")
print(f"  sigma       = {float(abs(target - cand_A) / c3_unc):.2f}")

# PSLQ: c3*LAM6 in [1, pi^2, pi^4]
print("\n--- PSLQ: c3*LAM6 in [1, pi^2, pi^4] ---")
for label, c3_val in [("full", target), ("partial", target_partial)]:
    vec = [c3_val * LAM6, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**6, tol=mpmath.mpf(10)**(-20))
        if rel:
            check = sum(r * v for r, v in zip(rel, vec))
            if rel[0] != 0:
                print(f"  ({label}) rel={rel}")
                print(f"    check={mpmath.nstr(check, 5)}")
                r0, r1, r2, r3 = rel
                print(f"    c3*LAM6 = ({-r1}/{r0}) + ({-r2}/{r0})*pi^2 + ({-r3}/{r0})*pi^4")
            else:
                print(f"  ({label}) TRIVIAL rel={rel} (c3 coefficient = 0)")
        else:
            print(f"  ({label}) No relation found")
    except Exception as e:
        print(f"  ({label}) PSLQ error: {e}")

# PSLQ: c3*LAM6 in extended basis [1, pi^2, pi^4, BD, FD, FB]
print("\n--- PSLQ: c3*LAM6 in [1, pi^2, pi^4, BD, FD, FB] ---")
labels_ext = ['c3*L6', '1', 'pi2', 'pi4', 'BD', 'FD', 'FB']
for label, c3_val in [("full", target), ("partial", target_partial)]:
    vec = [c3_val * LAM6, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4,
           BD, FD, FB]
    try:
        rel = mpmath.pslq(vec, maxcoeff=10**4, tol=mpmath.mpf(10)**(-20))
        if rel:
            check = sum(r * v for r, v in zip(rel, vec))
            if rel[0] != 0:
                print(f"  ({label}) rel={rel}")
                print(f"    labels: {labels_ext}")
                print(f"    check={mpmath.nstr(check, 5)}")
                terms = []
                for i in range(1, len(rel)):
                    if rel[i] != 0:
                        terms.append(f"({-rel[i]}/{rel[0]})*{labels_ext[i]}")
                print(f"    c3*LAM6 = {' + '.join(terms)}")
            else:
                print(f"  ({label}) TRIVIAL rel={rel}")
                # Show what it found
                terms = [f"{r}*{l}" for r, l in zip(rel[1:], labels_ext[1:]) if r != 0]
                print(f"    Identity: {' + '.join(terms)} = 0")
        else:
            print(f"  ({label}) No relation found")
    except Exception as e:
        print(f"  ({label}) PSLQ error: {e}")

# PSLQ: c3*LAM6 in degree-2 Paper 2 basis
print("\n--- PSLQ: c3*LAM6 in degree-2 basis ---")
labels_d2 = ['c3*L6', '1', 'pi2', 'pi4', 'BD', 'FD', 'FB',
             'BD^2', 'FD^2', 'FB^2', 'BD*FD', 'BD*FB', 'FD*FB', 'D^2']
vals_d2 = [target * LAM6, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4,
           BD, FD, FB,
           BD**2, FD**2, FB**2, BD*FD, BD*FB, FD*FB, D2]
try:
    rel = mpmath.pslq(vals_d2, maxcoeff=10**4, tol=mpmath.mpf(10)**(-20))
    if rel:
        check = sum(r * v for r, v in zip(rel, vals_d2))
        print(f"  rel={rel}")
        print(f"  labels: {labels_d2}")
        print(f"  check={mpmath.nstr(check, 5)}")
        if rel[0] != 0:
            terms = []
            for i in range(1, len(rel)):
                if rel[i] != 0:
                    terms.append(f"({-rel[i]}/{rel[0]})*{labels_d2[i]}")
            print(f"  c3*LAM6 = {' + '.join(terms)}")
        else:
            terms = [f"{r}*{l}" for r, l in zip(rel[1:], labels_d2[1:]) if r != 0]
            print(f"  TRIVIAL: {' + '.join(terms)} = 0")
    else:
        print(f"  No relation found")
except Exception as e:
    print(f"  PSLQ error: {e}")

# PSLQ: c3 in [c2, c2^2, c1*c2]
print("\n--- PSLQ: c3 in [c2, c2^2, c1*c2] ---")
labels_c = ['c3', 'c2', 'c2^2', 'c1*c2']
vals_c = [target, C2, C2**2, C1*C2]
try:
    rel = mpmath.pslq(vals_c, maxcoeff=10**6, tol=mpmath.mpf(10)**(-15))
    if rel:
        check = sum(r * v for r, v in zip(rel, vals_c))
        print(f"  rel={rel}")
        print(f"  check={mpmath.nstr(check, 5)}")
        if rel[0] != 0:
            terms = []
            for i in range(1, len(rel)):
                if rel[i] != 0:
                    terms.append(f"({-rel[i]}/{rel[0]})*{labels_c[i]}")
            print(f"  c3 = {' + '.join(terms)}")
    else:
        print(f"  No relation found")
except Exception as e:
    print(f"  PSLQ error: {e}")

# PSLQ with c2 included in main basis
print("\n--- PSLQ: c3*LAM6 in [1, pi^2, pi^4, BD, FD, FB, c2, c2^2, c1*c2] ---")
labels_c2 = ['c3*L6', '1', 'pi2', 'pi4', 'BD', 'FD', 'FB', 'c2', 'c2^2', 'c1*c2']
vals_c2 = [target * LAM6, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4,
           BD, FD, FB, C2, C2**2, C1*C2]
try:
    rel = mpmath.pslq(vals_c2, maxcoeff=10**4, tol=mpmath.mpf(10)**(-20))
    if rel:
        check = sum(r * v for r, v in zip(rel, vals_c2))
        print(f"  rel={rel}")
        print(f"  labels: {labels_c2}")
        print(f"  check={mpmath.nstr(check, 5)}")
        if rel[0] != 0:
            terms = []
            for i in range(1, len(rel)):
                if rel[i] != 0:
                    terms.append(f"({-rel[i]}/{rel[0]})*{labels_c2[i]}")
            print(f"  c3*LAM6 = {' + '.join(terms)}")
        else:
            terms = [f"{r}*{l}" for r, l in zip(rel[1:], labels_c2[1:]) if r != 0]
            print(f"  TRIVIAL: {' + '.join(terms)} = 0")
    else:
        print(f"  No relation found")
except Exception as e:
    print(f"  PSLQ error: {e}")

# ===================================================================
# SYSTEMATIC CANDIDATE SCAN
# ===================================================================
print("\n" + "=" * 70)
print("  STEP 5: Systematic candidate scan")
print("=" * 70)

# Generate degree-2 candidates in Paper 2 invariants / LAM6
monomials = {
    '1': mpmath.mpf(1),
    'BD': BD, 'FD': FD, 'FB': FB, 'D': D,
    'BD^2': BD**2, 'FD^2': FD**2, 'FB^2': FB**2,
    'BD*FD': BD*FD, 'BD*FB': BD*FB, 'FD*FB': FD*FB,
    'D^2': D2, 'BD*D': BD*D, 'FD*D': FD*D, 'FB*D': FB*D,
}

# Full scan: c3*LAM6 = sum of a_i * m_i for small integer a_i
target_scaled = target * LAM6
print(f"  Target: c3*LAM6 = {mpmath.nstr(target_scaled, 20)}")

hits = []
keys = list(monomials.keys())
vals = [monomials[k] for k in keys]

# 2-term combinations
for i in range(len(keys)):
    for j in range(i, len(keys)):
        for a in range(-10, 11):
            if a == 0: continue
            for b in range(-10, 11):
                if b == 0 and i == j: continue
                cand = a * vals[i] + b * vals[j]
                if abs(cand) < 1e-15: continue
                rel_err = abs(float((cand - target_scaled) / target_scaled))
                if rel_err < 1e-4:
                    sigma = float(abs(cand - target_scaled) / (c3_unc * LAM6))
                    hits.append((rel_err, sigma, a, keys[i], b, keys[j], float(cand)))

# 3-term combinations with small coefficients
for i in range(len(keys)):
    for j in range(i, len(keys)):
        for k in range(j, len(keys)):
            for a in range(-5, 6):
                if a == 0: continue
                for b in range(-5, 6):
                    for c in range(-5, 6):
                        if b == 0 and c == 0: continue
                        cand = a * vals[i] + b * vals[j] + c * vals[k]
                        if abs(cand) < 1e-15: continue
                        rel_err = abs(float((cand - target_scaled) / target_scaled))
                        if rel_err < 1e-5:
                            sigma = float(abs(cand - target_scaled) / (c3_unc * LAM6))
                            hits.append((rel_err, sigma, a, keys[i], b, keys[j], c, keys[k], float(cand)))

hits.sort()
print(f"\n  Found {len(hits)} hits within 0.01% relative error")
print(f"  Top 30 (sorted by rel_err):")
for h in hits[:30]:
    if len(h) == 7:
        rel_err, sigma, a, ki, b, kj, val = h
        print(f"    {a}*{ki} + {b}*{kj} = {val:.15e}  rel_err={rel_err:.3e}  sigma={sigma:.2f}")
    else:
        rel_err, sigma, a, ki, b, kj, c, kk, val = h
        print(f"    {a}*{ki} + {b}*{kj} + {c}*{kk} = {val:.15e}  rel_err={rel_err:.3e}  sigma={sigma:.2f}")

# ===================================================================
# V_mag EXACT COMPUTATION
# ===================================================================
print("\n" + "=" * 70)
print("  STEP 6: V_mag precision analysis")
print("=" * 70)

# V_mag = probe vertex normalization. Can we compute it algebraically?
# For n_int=0 probe on n_ext=1: it's a CG-weighted sum
# V_mag should be computable from CG coefficients (exact rationals)

# First, quantify V_mag's contribution to c3 uncertainty
# c3 = (cum_B/(V_mag * S_K) - 1 - c1*x - c2*x^2) / x^3
# dc3/dV_mag = -cum_B/(V_mag^2 * S_K * x^3)
dc3_dVmag = -cum_B / (V_mag**2 * SCHWINGER_K * X3)
# V_mag uncertainty: last digit of float
V_mag_unc = mpmath.mpf('1e-15')
c3_Vmag_unc = abs(dc3_dVmag) * V_mag_unc
print(f"  V_mag = {mpmath.nstr(V_mag, 20)}")
print(f"  dc3/dV_mag = {mpmath.nstr(dc3_dVmag, 10)}")
print(f"  V_mag float uncertainty ~ 1e-15")
print(f"  c3 uncertainty from V_mag = {mpmath.nstr(c3_Vmag_unc, 5)}")
print(f"  c3 uncertainty from tail  = {mpmath.nstr(c3_unc, 5)}")

# Precision budget
print(f"\n  PRECISION BUDGET:")
print(f"    alpha_K: exact (80 dps)")
print(f"    B values: float64 (~15 digits)")
print(f"    V_mag: float64 (~15 digits)")
print(f"    c3 V_mag contrib: {mpmath.nstr(c3_Vmag_unc, 3)}")
print(f"    c3 tail contrib:  {mpmath.nstr(c3_unc, 3)}")
print(f"    c3_partial reliable digits: ~{max(0, int(-mpmath.log10(abs(c3_Vmag_unc / c3_partial))))}")
print(f"    c3_full reliable digits:    ~{max(0, int(-mpmath.log10(abs(c3_unc / c3_full))))}")

# ===================================================================
# SUMMARY
# ===================================================================
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
print(f"  c3_partial (alpha_K, no tail) = {mpmath.nstr(c3_partial, 15)}")
print(f"  c3_full    (alpha_K, w/ tail) = {mpmath.nstr(c3_full, 15)}")
print(f"  c3_unc (tail)                 = {mpmath.nstr(c3_unc, 5)}")
print(f"  c3_unc (V_mag)                = {mpmath.nstr(c3_Vmag_unc, 5)}")
print(f"")
print(f"  Candidate A: (2-7BD-4FD+5BD^2+2FD^2-FB^2)/LAM6")
print(f"    = {mpmath.nstr(cand_A, 20)}")
print(f"    diff from c3_full = {mpmath.nstr(target - cand_A, 10)}")
print(f"    sigma = {float(abs(target - cand_A) / c3_unc):.3f}")
