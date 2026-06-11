"""
g2_c3_candidate_test.py -- Test Candidate A and extract c3 with maximum precision.

Strategy:
1. Compute Candidate A = (2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2)/LAM^6 exactly
2. Use Richardson extrapolation on the c3(N) convergence sequence to tighten the value
3. Run PSLQ at the best available precision

The c3(N) convergence data from n_cutoff=15..50 shows c3 approaching a limit.
The tail correction dominates the uncertainty. Richardson extrapolation on the
sequence c3(N) can extract the limit more precisely than any single N.
"""

import json
import sys
import mpmath
import numpy as np
from pathlib import Path

mpmath.mp.dps = 80

DATA_DIR = Path(__file__).parent / "data"

# Paper 2 invariants at full mpmath precision
B = mpmath.mpf(42)
F = mpmath.pi**2 / 6
D = mpmath.mpf(1) / 40  # Delta
ALPHA = mpmath.mpf('7.2973525693e-3')
SCHWINGER = ALPHA / (2 * mpmath.pi)

C1 = mpmath.mpf(1) / 2
C2 = (2 - B*D - F*D - F/B) / 5

# Bilinears
BD = B * D   # = 21/20
FD = F * D   # = pi^2/240
FB = F / B   # = pi^2/252

# n_ext=1 kinematics
LAM = mpmath.mpf(5) / 2
LAM2 = LAM**2
X = 1 / LAM2  # 4/25
X3 = X**3
LAM6 = LAM**6  # 15625/64

print("=" * 70)
print("  c3 CANDIDATE TEST AND PRECISION EXTRACTION")
print("=" * 70)

print(f"\nPaper 2 invariants:")
print(f"  B  = {B}")
print(f"  F  = {mpmath.nstr(F, 25)}")
print(f"  D  = {D}")
print(f"  BD = {mpmath.nstr(BD, 25)} = {mpmath.nstr(mpmath.mpf(21)/20, 25)}")
print(f"  FD = {mpmath.nstr(FD, 25)}")
print(f"  FB = {mpmath.nstr(FB, 25)}")
print(f"  c2 = {mpmath.nstr(C2, 25)}")
print(f"  x  = {mpmath.nstr(X, 15)} = 4/25")
print(f"  x^3= {mpmath.nstr(X3, 15)}")
print(f"  LAM^6 = {mpmath.nstr(LAM6, 15)} = 15625/64")

# =====================================================================
# STEP 1: Candidate A evaluation
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 1: CANDIDATE A EVALUATION")
print("=" * 70)

# Candidate A: c3 * LAM^6 = 2 - 7*BD - 4*FD + 5*BD^2 + 2*FD^2 - FB^2
cand_A_scaled = 2 - 7*BD - 4*FD + 5*BD**2 + 2*FD**2 - FB**2
cand_A_c3 = cand_A_scaled / LAM6

print(f"\nCandidate A: c3*LAM^6 = 2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2")
print(f"  2        = {mpmath.nstr(mpmath.mpf(2), 25)}")
print(f"  -7BD     = {mpmath.nstr(-7*BD, 25)}")
print(f"  -4FD     = {mpmath.nstr(-4*FD, 25)}")
print(f"  +5BD^2   = {mpmath.nstr(5*BD**2, 25)}")
print(f"  +2FD^2   = {mpmath.nstr(2*FD**2, 25)}")
print(f"  -FB^2    = {mpmath.nstr(-FB**2, 25)}")
print(f"  SUM      = {mpmath.nstr(cand_A_scaled, 25)}")
print(f"  c3_A     = {mpmath.nstr(cand_A_c3, 20)}")

# Rational + pi^2 + pi^4 decomposition
# 2 - 7*(21/20) - 4*(pi^2/240) + 5*(21/20)^2 + 2*(pi^2/240)^2 - (pi^2/252)^2
# = 2 - 147/20 + 5*441/400 - pi^2/60 + 2*pi^4/57600 - pi^4/63504
# = 2 - 7.35 + 5.5125 - pi^2/60 + pi^4*(2/57600 - 1/63504)
# rational part: 2 - 147/20 + 2205/400 = (800 - 2940 + 2205)/400 = 65/400 = 13/80
r0 = mpmath.mpf(13) / 80  # rational part
# pi^2 part: -4*FD = -4*pi^2/240 = -pi^2/60
r2_coeff = mpmath.mpf(-1) / 60  # coefficient of pi^2
# pi^4 part: 2*(pi^2/240)^2 - (pi^2/252)^2 = pi^4*(2/57600 - 1/63504)
# 2/57600 = 1/28800; 1/63504
# LCM of 28800 and 63504:
# 28800 = 2^6 * 3^2 * 5^2
# 63504 = 2^4 * 3^2 * 441 = 2^4 * 3^2 * 21^2 = 2^4 * 3^2 * 3^2 * 7^2 = 2^4 * 3^4 * 7^2
# Actually: 63504 = 63504. Let me factor: 63504/2=31752/2=15876/2=7938/2=3969
# 3969 = 63^2 = (7*9)^2 = 9*441 = 9*21^2
# So 63504 = 16*3969 = 16*63^2 = 2^4 * 63^2 = 2^4 * (7*9)^2 = 2^4 * 7^2 * 3^4
# LCM = 2^6 * 3^4 * 5^2 * 7^2 = 64 * 81 * 25 * 49 = 64*81*1225 = 6350400
# 1/28800 = 220.5/6350400 ... hmm, let me just compute numerically
r4_coeff = 2 / mpmath.mpf(57600) - 1 / mpmath.mpf(63504)

# Exact fractions:
# 2/57600 = 1/28800
# 1/63504
# 1/28800 - 1/63504 = (63504 - 28800)/(28800*63504) = 34704/1828915200
# Simplify: gcd(34704, 1828915200)
# 34704 = 2^4 * 3 * 723 = 2^4 * 3 * 3 * 241 = 2^4 * 9 * 241 = 2^4 * 3^2 * 241
# 1828915200 = 28800 * 63504 = (2^6*3^2*5^2) * (2^4*3^4*7^2)
#            = 2^10 * 3^6 * 5^2 * 7^2
# gcd = 2^4 * 3^2 = 144
# 34704/144 = 241
# 1828915200/144 = 12700800
# So r4_coeff = 241/12700800

r4_exact = mpmath.mpf(241) / 12700800
print(f"\nDecomposition: c3*LAM^6 = r0 + r2*pi^2 + r4*pi^4")
print(f"  r0 = 13/80 = {mpmath.nstr(r0, 25)}")
print(f"  r2 = -1/60 = {mpmath.nstr(r2_coeff, 25)}")
print(f"  r4 = 241/12700800 = {mpmath.nstr(r4_exact, 25)}")

# Cross-check
cand_A_check = r0 + r2_coeff * mpmath.pi**2 + r4_exact * mpmath.pi**4
print(f"\n  Cross-check: r0+r2*pi^2+r4*pi^4 = {mpmath.nstr(cand_A_check, 25)}")
print(f"  Direct computation           = {mpmath.nstr(cand_A_scaled, 25)}")
print(f"  Difference                   = {mpmath.nstr(cand_A_check - cand_A_scaled, 5)}")

c3_candA = cand_A_check / LAM6
print(f"\n  c3 (Candidate A) = {mpmath.nstr(c3_candA, 25)}")

# =====================================================================
# STEP 2: Load measured c3 and compare
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 2: COMPARISON WITH MEASURED c3")
print("=" * 70)

with open(DATA_DIR / "g2_c3_investigation.json") as f:
    data = json.load(f)

all_B = {int(k): float(v) for k, v in data["all_B"].items()}
V_mag = mpmath.mpf(str(data["V_mag"]))

# Best measured c3 from convergence study
conv = data["c3_convergence"]
c3_meas = conv[-1]["c3"]  # n_cutoff=50
c3_unc = conv[-1]["c3_unc"]
c3_sigma_meas = conv[-1]["c3_sigma"]

print(f"\n  c3 (measured, n_cutoff=50) = {c3_meas:.10e}")
print(f"  c3 uncertainty            = {c3_unc:.3e}")
print(f"  sigma                     = {c3_sigma_meas:.1f}")

diff = float(c3_candA) - c3_meas
sigma_diff = abs(diff) / c3_unc
print(f"\n  c3 (Candidate A)          = {float(c3_candA):.10e}")
print(f"  Difference                = {diff:.3e}")
print(f"  Sigma                     = {sigma_diff:.2f}")

# =====================================================================
# STEP 3: Richardson extrapolation on c3(N) sequence
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 3: RICHARDSON EXTRAPOLATION")
print("=" * 70)

# The c3(N) values converge as c3_inf + A/N^p for some p.
# Use the convergence study data
ns = [c["n_int_cutoff"] for c in conv]
c3s = [c["c3"] for c in conv]

print(f"\nConvergence data:")
for i, (n, c3v) in enumerate(zip(ns, c3s)):
    print(f"  N={n:3d}: c3 = {c3v:18.10e}")

# Aitken's delta-squared for acceleration
# Works on sequences converging as x_n -> L + C*r^n
print("\nAitken delta-squared acceleration (last 3 points):")
for i in range(len(c3s) - 2):
    x0, x1, x2 = c3s[i], c3s[i+1], c3s[i+2]
    denom = x2 - 2*x1 + x0
    if abs(denom) > 1e-30:
        x_ait = x0 - (x1 - x0)**2 / denom
        print(f"  Aitken({ns[i]},{ns[i+1]},{ns[i+2]}): c3 = {x_ait:18.10e}")

# Richardson extrapolation assuming c3(N) = c3_inf + A*N^(-p)
# Use log-log fit on consecutive differences
print("\nRichardson analysis on consecutive differences:")
diffs = []
for i in range(1, len(c3s)):
    d = c3s[i] - c3s[i-1]
    diffs.append((ns[i], d))
    if i >= 3:
        # Estimate convergence rate from ratio of consecutive differences
        d_prev = c3s[i-1] - c3s[i-2]
        if abs(d_prev) > 1e-30 and abs(d) > 1e-30:
            ratio = d / d_prev
            n_ratio = ns[i-1] / ns[i-2]
            if abs(ratio) > 0 and abs(ratio) < 1:
                p_est = -np.log(abs(ratio)) / np.log(ns[i] / ns[i-1])
                print(f"  N={ns[i]}: delta={d:12.4e}  ratio={ratio:.4f}  p_est={p_est:.2f}")

# Use last 4 points for quadratic Richardson
print("\nQuadratic extrapolation (last 4 points):")
n_last = np.array(ns[-4:], dtype=float)
c_last = np.array(c3s[-4:], dtype=float)
# Fit c3(N) = c3_inf + a/N^p + b/N^(2p) using least squares
# First estimate p from the approach rate
# c3(50) - c3(45) ≈ A*(45^-p - 50^-p)
# c3(45) - c3(40) ≈ A*(40^-p - 45^-p)
d1 = c3s[-1] - c3s[-2]  # 50 vs 45
d2 = c3s[-2] - c3s[-3]  # 45 vs 40
d3 = c3s[-3] - c3s[-4]  # 40 vs 35

# For power-law tail B(n) ~ C*n^(-p), the partial sum tail goes as ~n^(1-p)
# So c3(N) ~ c3_inf + const * N^(1-p)
# Ratio of differences: d1/d2 ≈ (N2^(1-p) - N1^(1-p)) / (N3^(1-p) - N2^(1-p))
# This is a known analysis pattern

# Simple three-point extrapolation: assume c3(N) = L + A/N^q
# (c3(N2) - c3(N1)) / (c3(N3) - c3(N2)) ≈ (N1^-q - N2^-q) / (N2^-q - N3^-q)
# Try q=5,6,7 and see which gives consistent L
print("\nParametric extrapolation c3(N) = L + A/N^q:")
for q in [4, 5, 6, 7, 8]:
    # Two equations, two unknowns (L, A) from last 2 points
    N1, N2 = float(ns[-2]), float(ns[-1])
    c1, c2 = c3s[-2], c3s[-1]
    A_est = (c2 - c1) / (N2**(-q) - N1**(-q))
    L_est = c1 - A_est * N1**(-q)

    # Cross-check with third point
    N3 = float(ns[-3])
    c3_pred = L_est + A_est * N3**(-q)
    resid = c3_pred - c3s[-3]

    sigma_L = abs(L_est - float(c3_candA)) / c3_unc
    print(f"  q={q}: L={L_est:18.10e}  A={A_est:12.4e}  "
          f"resid(N={ns[-3]})={resid:10.2e}  "
          f"sigma_vs_A={sigma_L:.2f}")

# =====================================================================
# STEP 4: Best-estimate c3 and PSLQ
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 4: PSLQ WITH BEST c3 ESTIMATE")
print("=" * 70)

# Use Richardson-extrapolated c3 (q=6 or 7 tends to work for power-law convergence)
# But first, let's also try a direct Shanks transformation
print("\nShanks transformation (epsilon algorithm):")
# Shanks transform of the last few points
e_table = [c3s[-5:]]  # row 0
for k in range(1, 5):
    row = []
    for j in range(len(e_table[k-1]) - 1):
        denom = e_table[k-1][j+1] - e_table[k-1][j]
        if abs(denom) < 1e-30:
            break
        if k == 1:
            row.append(1.0 / denom)
        else:
            row.append(e_table[k-2][j+1] + 1.0 / denom)
    e_table.append(row)
    if len(row) == 0:
        break

print("  Epsilon table (even rows are estimates):")
for k, row in enumerate(e_table):
    if k % 2 == 0 and len(row) > 0:
        for j, val in enumerate(row):
            sigma = abs(val - float(c3_candA)) / c3_unc
            print(f"    e[{k}][{j}] = {val:18.10e}  (sigma_vs_A={sigma:.2f})")

# Try Wynn's epsilon algorithm more carefully
print("\nWynn epsilon algorithm on c3(N=15..50):")
s_list = c3s[:]
n_pts = len(s_list)
eps = [[mpmath.mpf(0)] * (n_pts + 1) for _ in range(n_pts + 1)]
for i in range(n_pts):
    eps[i][1] = mpmath.mpf(str(s_list[i]))

for k in range(2, n_pts + 1):
    for i in range(n_pts - k + 1):
        diff = eps[i+1][k-1] - eps[i][k-1]
        if abs(diff) < mpmath.mpf('1e-50'):
            eps[i][k] = mpmath.mpf('1e50')
        else:
            eps[i][k] = eps[i+1][k-2] + 1/diff

print("  Even-column Wynn estimates (series limit accelerants):")
for k in range(2, n_pts + 1, 2):
    for i in range(min(3, n_pts - k + 1)):
        val = float(eps[i][k])
        if abs(val) < 1e-3:
            sigma = abs(val - float(c3_candA)) / c3_unc
            print(f"    eps[{i}][{k}] = {val:18.10e}  (sigma_vs_A={sigma:.2f})")

# Best estimate: use the most converged Richardson/Wynn value
# For now, use the n_cutoff=50 measured value with its uncertainty
c3_best = mpmath.mpf(str(c3_meas))
c3_best_unc = c3_unc

# Try PSLQ with Candidate A's structure
# c3 = (13/80 - pi^2/60 + 241*pi^4/12700800) / (15625/64)
# = 64*(13/80 - pi^2/60 + 241*pi^4/12700800) / 15625
# = (64*13/80 - 64*pi^2/60 + 64*241*pi^4/12700800) / 15625

print("\n--- PSLQ identification ---")
pi2 = mpmath.pi**2
pi4 = mpmath.pi**4

# T9 basis: c3 = (a0 + a2*pi^2 + a4*pi^4) / D
# where a0, a2, a4 are rationals and D is integer
# We search for c3*LAM^6 = a0 + a2*pi^2 + a4*pi^4

c3L6 = c3_best * LAM6
print(f"\n  c3*LAM^6 = {mpmath.nstr(c3L6, 15)}")
print(f"  Candidate A: {mpmath.nstr(cand_A_scaled, 15)}")
print(f"  Diff: {mpmath.nstr(c3L6 - cand_A_scaled, 5)}")

# PSLQ on [c3*LAM^6, 1, pi^2, pi^4]
mpmath.mp.dps = 50
vec = mpmath.matrix([c3L6, mpmath.mpf(1), pi2, pi4])
rel = mpmath.pslq(vec, maxcoeff=10000, tol=mpmath.mpf('1e-8'))
print(f"\n  PSLQ [c3L6, 1, pi^2, pi^4] maxcoeff=10000:")
if rel:
    check = sum(float(rel[j]) * float(vec[j]) for j in range(len(vec)))
    print(f"    relation: {list(rel)}")
    print(f"    check: {check:.3e}")
    if rel[0] != 0:
        print(f"    c3L6 = ({-rel[1]}/{rel[0]})*1 + ({-rel[2]}/{rel[0]})*pi^2 + ({-rel[3]}/{rel[0]})*pi^4")
else:
    print(f"    No relation found")

# Try Paper 2 bilinear basis
# c3*LAM^6 might be a polynomial in BD, FD, FB, D
bases = {
    "pure_T9": {
        "labels": ["c3L6", "1", "pi^2", "pi^4"],
        "vals": [c3L6, 1, pi2, pi4],
    },
    "paper2_deg1": {
        "labels": ["c3L6", "1", "BD", "FD", "FB", "D"],
        "vals": [c3L6, 1, BD, FD, FB, D],
    },
    "paper2_deg2": {
        "labels": ["c3L6", "1", "BD", "FD", "FB", "BD^2", "FD^2", "FB^2", "BD*FD", "BD*FB", "FD*FB"],
        "vals": [c3L6, 1, BD, FD, FB, BD**2, FD**2, FB**2, BD*FD, BD*FB, FD*FB],
    },
    "c2_basis": {
        "labels": ["c3L6", "c2", "c2^2", "c1*c2", "BD*c2", "FD*c2", "FB*c2", "D*c2"],
        "vals": [c3L6, C2, C2**2, C1*C2, BD*C2, FD*C2, FB*C2, D*C2],
    },
    "with_D_and_c2": {
        "labels": ["c3L6", "1", "BD", "FD", "FB", "D", "c2", "BD^2", "D^2"],
        "vals": [c3L6, 1, BD, FD, FB, D, C2, BD**2, D**2],
    },
}

print(f"\n--- Systematic PSLQ search ---")
for bname, basis in bases.items():
    labels = basis["labels"]
    vals = basis["vals"]
    vec = mpmath.matrix(vals)

    for mc in [100, 1000, 10000]:
        rel = mpmath.pslq(vec, maxcoeff=mc, tol=mpmath.mpf('1e-6'))
        if rel:
            check = sum(float(rel[j]) * float(vec[j]) for j in range(len(vec)))
            # Check if c3L6 coefficient is nonzero
            is_trivial = (rel[0] == 0)
            tag = "TRIVIAL" if is_trivial else "HIT"
            print(f"  {bname} (mc={mc}): {tag}  rel={list(rel)}  check={check:.3e}")
            if not is_trivial:
                # Decode the relation
                terms = []
                for j in range(1, len(rel)):
                    if rel[j] != 0:
                        terms.append(f"{-rel[j]}/{rel[0]}*{labels[j]}")
                print(f"    c3L6 = {' + '.join(terms)}")

                # Compute the predicted c3 from the relation
                c3L6_pred = sum(-mpmath.mpf(rel[j]) * vals[j] for j in range(1, len(rel))) / rel[0]
                c3_pred = c3L6_pred / LAM6
                sigma = abs(float(c3_pred) - c3_meas) / c3_unc
                print(f"    c3_pred = {mpmath.nstr(c3_pred, 15)}")
                print(f"    sigma vs measured = {sigma:.2f}")
            break  # found at this maxcoeff, skip larger
        else:
            if mc == 10000:
                print(f"  {bname}: no relation at maxcoeff={mc}")

# =====================================================================
# STEP 5: Direct formula candidates from c2 structure
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 5: STRUCTURAL CANDIDATES FROM c2 PATTERN")
print("=" * 70)

# c2 = (2 - BD - FD - FB) / 5
# c2 can also be written as: c2 = (2 - BD - FD - FB) / 5
#
# Natural generalization patterns for c3:
# Pattern 1: polynomial in (BD, FD, FB) of degree 3
# Pattern 2: c3 involves D^2 terms (next order)
# Pattern 3: c3 = f(c1, c2, B, F, D)

# Let's check: does c3 = c2^2/k for some integer k?
ratio_c2sq = float(c3_best) / float(C2**2)
print(f"\n  c3/c2^2 = {ratio_c2sq:.6e}")
print(f"  c3/(c2*c1) = {float(c3_best) / float(C2*C1):.6e}")

# Does c3*LAM^6 = A*c2^2*LAM^4 + B_term?
# c2*x^2 = c2/LAM^4, so c2*LAM^4 is a natural scale
c2L4 = C2 * LAM**4
print(f"\n  c2*LAM^4 = {mpmath.nstr(c2L4, 15)}")
print(f"  (c2*LAM^4)^2 = {mpmath.nstr(c2L4**2, 15)}")
ratio = c3L6 / c2L4**2
print(f"  c3*LAM^6 / (c2*LAM^4)^2 = {mpmath.nstr(ratio, 15)}")

# Check specific structural formulas from heat kernel expansion
# On S3, R=6, R_ij R^ij = 6, R_ijkl R^ijkl = 12
# The Parker-Toms expansion has known structure at each order
# c1 = R/12 = 1/2 (confirmed)
# c2 involves R^2 terms and the electromagnetic field

# For c3, the heat kernel gives terms ~ R^3, R*Rij*Rij, Rijkl*Rijkl*R, etc.
# On maximally symmetric S3: everything reduces to powers of R

# Try: c3 = polynomial in c1, c2, and R-dependent terms
# R = 6 on unit S3, so R^k / (some normalization)
print(f"\n  Checking c3 = a*c1*c2*X + b*c2^2*X + ...  (higher-order heat kernel)")

# Actually, the curvature expansion should be:
# F2/S = 1 + c1*X + c2*X^2 + c3*X^3 + ...
# where X = 1/lambda^2 and the c_k are FIXED constants (curvature-dependent)
# On unit S3 they are specific numbers.

# Let's check: is c3 = -(c1*c2^2 + something simple)/denominator ?
test_vals = {
    "c1*c2": C1 * C2,
    "c2^2": C2**2,
    "c1^2*c2": C1**2 * C2,
    "c1*c2^2": C1 * C2**2,
    "c2/B": C2 / B,
    "c2*D": C2 * D,
    "c2*BD": C2 * BD,
    "c2*FD": C2 * FD,
    "c2*FB": C2 * FB,
}

print(f"\n  c3 = {mpmath.nstr(c3_best, 12)}")
print(f"\n  Ratios c3/quantity:")
for name, val in test_vals.items():
    r = float(c3_best) / float(val)
    print(f"    c3/({name}) = {r:.8e}")

# =====================================================================
# STEP 6: Alternative candidates from dimensional analysis
# =====================================================================
print("\n" + "=" * 70)
print("  STEP 6: SYSTEMATIC CANDIDATE GENERATION")
print("=" * 70)

# c3*LAM^6 ≈ -1.45e-4
# Look for combinations of {1, BD, FD, FB, D, BD^2, FD^2, FB^2, BD*FD, BD*FB, FD*FB, D^2}
# with SMALL integer coefficients that match

target = float(c3L6)
print(f"\n  Target: c3*LAM^6 = {target:.10e}")

# The rational part of Candidate A is 13/80 = 0.1625
# The pi^2 part is -pi^2/60 = -0.16449...
# The pi^4 part is 241*pi^4/12700800 = 0.001849...
# Total: 0.1625 - 0.16449 + 0.00185 = -0.000145...

# Let me check: is there a SIMPLER formula?
# c3*LAM^6 = (a + b*pi^2 + c*pi^4) / D for small integers a,b,c,D

print("\n  Searching c3*LAM^6 = (a + b*pi^2)/D for small a,b,D...")
best_hits = []
for D_val in range(1, 200):
    for a in range(-20, 21):
        for b in range(-20, 21):
            if a == 0 and b == 0:
                continue
            val = (a + b * float(pi2)) / D_val
            rel_err = abs(val - target) / abs(target) if target != 0 else abs(val)
            if rel_err < 0.01:  # within 1%
                best_hits.append((rel_err, a, b, D_val, val))

best_hits.sort()
print(f"  Found {len(best_hits)} hits within 1%:")
for rel_err, a, b, d, val in best_hits[:15]:
    sigma = abs(val - target) / c3_unc
    print(f"    ({a} + {b}*pi^2)/{d} = {val:.10e}  rel_err={rel_err:.4e}  sigma={sigma:.1f}")

# Also search with pi^4
print("\n  Searching c3*LAM^6 = (a + b*pi^2 + c*pi^4)/D for small a,b,c,D...")
best_hits_3 = []
for D_val in range(1, 100):
    for a in range(-10, 11):
        for b in range(-10, 11):
            for c in range(-5, 6):
                if a == 0 and b == 0 and c == 0:
                    continue
                val = (a + b * float(pi2) + c * float(pi4)) / D_val
                rel_err = abs(val - target) / abs(target) if target != 0 else abs(val)
                if rel_err < 1e-4:  # within 0.01%
                    best_hits_3.append((rel_err, a, b, c, D_val, val))

best_hits_3.sort()
print(f"  Found {len(best_hits_3)} hits within 0.01%:")
for rel_err, a, b, c, d, val in best_hits_3[:20]:
    sigma = abs(val - target) / c3_unc
    print(f"    ({a} + {b}*pi^2 + {c}*pi^4)/{d} = {val:.10e}  rel_err={rel_err:.6e}  sigma={sigma:.1f}")

# =====================================================================
# Save results
# =====================================================================
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)

print(f"\n  c3 (measured, N=50)  = {c3_meas:.10e} ± {c3_unc:.3e}")
print(f"  c3 (Candidate A)     = {float(c3_candA):.10e}")
print(f"  sigma (A vs meas)    = {sigma_diff:.2f}")
print(f"\n  Candidate A: c3*LAM^6 = 13/80 - pi^2/60 + 241*pi^4/12700800")
print(f"             = (2 - 7BD - 4FD + 5BD^2 + 2FD^2 - FB^2)")

results = {
    "c3_measured": c3_meas,
    "c3_unc": c3_unc,
    "c3_candidate_A": float(c3_candA),
    "sigma_A_vs_meas": sigma_diff,
    "candidate_A_formula": {
        "r0": "13/80",
        "r2": "-1/60",
        "r4": "241/12700800",
        "bilinear": "2 - 7*BD - 4*FD + 5*BD^2 + 2*FD^2 - FB^2",
    },
}

outfile = DATA_DIR / "g2_c3_candidate_test.json"
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {outfile}")


if __name__ == "__main__":
    main() if 'main' in dir() else None
