"""
g-2 curvature coefficient c2 vs Paper 2 spectral invariants: PSLQ identification.

The QED g-2 computation on S^3 produces a curvature expansion:
  F2/Schwinger = 1 + c1/lam^2 + c2/lam^4 + ...
where lam is the Dirac eigenvalue cutoff.

We know c1 = 1/2 exactly (Parker-Toms, R/12 = 6/12 for unit S^3).
We have c2 ~ 0.17394 at ~7 reliable digits from n_ext=1, n_int=0..25.

Question: is c2 built from the same spectral invariants that appear in
Paper 2's alpha formula K = pi(B + F - Delta)?

Paper 2 invariants:
  B = 42  (finite Casimir trace)
  F = pi^2/6 ~ 1.6449  (Fock degeneracy Dirichlet at d_max=4)
  Delta = 1/40 = 0.025  (Dirac mode degeneracy boundary term)
  K = pi(B + F - Delta) ~ 137.036  (the alpha prediction)
  alpha ~ 1/137.036
  c1 = 1/2

Also relevant spectral invariants:
  |lam_3| = 9/2  (Dirac eigenvalue at cutoff)
  g_3^Dirac = 40 = Delta^{-1}
  N(2) = 5  (cumulative state count)
  dim(S^3) = 3
"""

import json
import sys
import os
from math import gcd

sys.path.insert(0, '.')

import mpmath
mpmath.mp.dps = 60  # 60 digits working precision

# =====================================================================
# Load data and compute c2
# =====================================================================

with open('debug/data/g2_extended_nint_v2.json') as f:
    v2_data = json.load(f)

with open('debug/data/g2_extended_nint.json') as f:
    v1_data = json.load(f)

# Best available delta: from v2 (n_int=0..25, tail-corrected)
delta_v2 = mpmath.mpf(str(v2_data['corrected_delta']))
# Cross-check: from v1 (n_int=0..15, tail-corrected)
delta_v1 = mpmath.mpf(str(v1_data['corrected_delta']))

# lam = (2*n_ext + 3)/2 at n_ext=1 -> 5/2
lam = mpmath.mpf('5') / 2
lam2 = lam**2  # 25/4
lam4 = lam2**2  # 625/16

# c1 = R/12 = 1/2 (exact, Parker-Toms)
c1 = mpmath.mpf('1') / 2

# c2_apparent = (delta - c1/lam^2) * lam^4
# But this is contaminated by c3/lam^2 + c4/lam^4 + ...
# At n_ext=1, lam^2 = 25/4, so c3 contamination is c3 * 4/25 ~ 16% if c3 ~ c2
c2_v2 = (delta_v2 - c1 / lam2) * lam4
c2_v1 = (delta_v1 - c1 / lam2) * lam4

# Estimate c3 contamination from the two datasets
# c2_apparent = c2_true + c3/lam^2 + c4/lam^4 + ...
# The difference between v1 and v2 comes from tail corrections, not c3

print("=" * 70)
print("  c2 vs Paper 2 spectral invariants: PSLQ identification")
print("=" * 70)

print(f"\n  delta (v2, n_int=0..25) = {mpmath.nstr(delta_v2, 15)}")
print(f"  delta (v1, n_int=0..15) = {mpmath.nstr(delta_v1, 15)}")
print(f"  Difference: {float(abs(delta_v2 - delta_v1)):.2e}")
print(f"\n  c1 = 1/2 (exact, Parker-Toms)")
print(f"  c1/lam^2 = 2/25 = {float(c1/lam2):.10f}")
print(f"\n  c2 (from v2 delta) = {mpmath.nstr(c2_v2, 15)}")
print(f"  c2 (from v1 delta) = {mpmath.nstr(c2_v1, 15)}")
print(f"  c2 difference: {float(abs(c2_v2 - c2_v1)):.2e}")
print(f"\n  Using c2 = {mpmath.nstr(c2_v2, 15)} for PSLQ searches")
print(f"  Estimated precision: ~7 digits (tail + c3 contamination)")

c2 = c2_v2

# =====================================================================
# Paper 2 invariants at high precision
# =====================================================================

B = mpmath.mpf(42)
F = mpmath.pi**2 / 6
Delta = mpmath.mpf(1) / 40
K = mpmath.pi * (B + F - Delta)
alpha = 1 / K
c1_val = mpmath.mpf('1') / 2

# Additional spectral invariants
lam3 = mpmath.mpf('9') / 2  # |lambda_3| Dirac eigenvalue at cutoff
g3_dirac = mpmath.mpf(40)   # g_3^Dirac = Delta^{-1}
N2 = mpmath.mpf(5)          # cumulative state count
dim_S3 = mpmath.mpf(3)
R_scalar = mpmath.mpf(6)    # Ricci scalar of unit S^3

print(f"\n  Paper 2 invariants:")
print(f"    B = {B}")
print(f"    F = pi^2/6 = {mpmath.nstr(F, 15)}")
print(f"    Delta = 1/40 = {Delta}")
print(f"    K = pi(B+F-Delta) = {mpmath.nstr(K, 15)}")
print(f"    alpha = 1/K = {mpmath.nstr(alpha, 15)}")
print(f"    c1 = 1/2 = {c1_val}")
print(f"    |lam_3| = 9/2 = {lam3}")
print(f"    g_3^Dirac = 40 = {g3_dirac}")
print(f"    N(2) = 5 = {N2}")
print(f"    dim(S^3) = 3 = {dim_S3}")
print(f"    R = 6 = {R_scalar}")

# =====================================================================
# Specific candidates
# =====================================================================

print("\n" + "=" * 70)
print("  SPECIFIC CANDIDATES")
print("=" * 70)

candidates = [
    ("7*Delta = 7/40", 7 * Delta),
    ("c1*Delta*B = (1/2)(1/40)(42) = 21/40", c1_val * Delta * B),
    ("Delta*(B+F-Delta)/pi = Delta*K/pi^2", Delta * K / mpmath.pi**2),
    ("c1^2*Delta*B = (1/4)(1/40)(42) = 42/160", c1_val**2 * Delta * B),
    ("F / ((B-2)*Delta)", F / ((B - 2) * Delta)),
    ("c1*F/B", c1_val * F / B),
    ("Delta*B/pi", Delta * B / mpmath.pi),
    ("1/(2*B*Delta) = 1/2.1", 1 / (2 * B * Delta)),
    ("c1*Delta*(B-2)", c1_val * Delta * (B - 2)),
    ("Delta^2*B", Delta**2 * B),
    ("c1/dim_S3", c1_val / dim_S3),
    ("1/R", 1 / R_scalar),
    ("R*Delta/2", R_scalar * Delta / 2),
    ("c1*Delta*R", c1_val * Delta * R_scalar),
    ("1/(4*N2) = 1/20", mpmath.mpf(1) / (4 * N2)),
    ("Delta*R = 6/40 = 3/20", Delta * R_scalar),
    ("c1^2*R/12 = 1/8", c1_val**2 * R_scalar / 12),
    ("7/40", mpmath.mpf(7) / 40),
    ("c1^2 - Delta = 1/4 - 1/40 = 9/40", c1_val**2 - Delta),
    ("(c1 - Delta)/dim_S3", (c1_val - Delta) / dim_S3),
    ("c1*(c1-Delta) = (1/2)(19/40)", c1_val * (c1_val - Delta)),
]

print(f"\n  {'Candidate':>50s}  {'Value':>14s}  {'Error':>10s}  {'Rel%':>8s}")
print(f"  {'-'*50}  {'-'*14}  {'-'*10}  {'-'*8}")
for label, val in sorted(candidates, key=lambda x: abs(float(x[1]) - float(c2))):
    err = abs(float(val) - float(c2))
    rel = err / abs(float(c2)) * 100
    marker = " <--" if rel < 1 else ""
    print(f"  {label:>50s}  {float(val):14.10f}  {err:10.2e}  {rel:7.3f}%{marker}")

# =====================================================================
# PSLQ BASIS SET SEARCHES
# =====================================================================

print("\n" + "=" * 70)
print("  PSLQ SEARCHES (10 basis sets)")
print("=" * 70)

results = {}
pslq_tol = mpmath.mpf('1e-5')  # tolerance matched to ~7 digit precision
max_coeff = 1000

def run_pslq(basis_name, labels, vals, tol=pslq_tol, maxc=max_coeff):
    """Run PSLQ and report results."""
    print(f"\n  --- Basis {basis_name}: {labels} ---")
    try:
        rel = mpmath.pslq(vals, tol=float(tol), maxcoeff=maxc)
        if rel is not None and rel[0] != 0:
            # Build the formula string
            terms = []
            for r, l in zip(rel, labels):
                if r != 0:
                    terms.append(f"({r})*{l}")
            formula = " + ".join(terms) + " = 0"

            # Solve for c2
            c2_pred = sum(-mpmath.mpf(r) * v for r, v, l in
                         zip(rel[1:], vals[1:], labels[1:])) / rel[0]
            resid = abs(float(c2_pred) - float(c2))

            print(f"    HIT: {formula}")
            print(f"    => c2 = {mpmath.nstr(c2_pred, 15)}")
            print(f"    Residual: {resid:.2e}")
            print(f"    Coefficients: {rel}")

            result = {
                'hit': True,
                'relation': rel,
                'formula': formula,
                'c2_predicted': float(c2_pred),
                'residual': resid,
            }
        else:
            print(f"    No relation found (maxcoeff={maxc})")
            result = {'hit': False}
    except Exception as e:
        print(f"    Error: {e}")
        result = {'hit': False, 'error': str(e)}

    results[basis_name] = result
    return result


# Basis 1: {c2, 1, 1/B, Delta, 1/K} -- simple reciprocals
run_pslq("1", ["c2", "1", "1/B", "Delta", "1/K"],
         [c2, mpmath.mpf(1), 1/B, Delta, 1/K])

# Basis 2: {c2, 1, Delta, c1*Delta, Delta^2} -- powers of Delta
run_pslq("2", ["c2", "1", "Delta", "c1*Delta", "Delta^2"],
         [c2, mpmath.mpf(1), Delta, c1_val*Delta, Delta**2])

# Basis 3: {c2, 1, pi^2, Delta, 1/B} -- mixed
run_pslq("3", ["c2", "1", "pi^2", "Delta", "1/B"],
         [c2, mpmath.mpf(1), mpmath.pi**2, Delta, 1/B])

# Basis 4: {c2, 1, F, Delta, F*Delta, c1} -- F and Delta combinations
run_pslq("4", ["c2", "1", "F", "Delta", "F*Delta", "c1"],
         [c2, mpmath.mpf(1), F, Delta, F*Delta, c1_val])

# Basis 5: {c2, 1, 1/(B-2), Delta, F/B} -- compound ratios
run_pslq("5", ["c2", "1", "1/(B-2)", "Delta", "F/B"],
         [c2, mpmath.mpf(1), 1/(B-2), Delta, F/B])

# Basis 6: {c2, 1, c1^2, c1*Delta, Delta*F, c1/B} -- products
run_pslq("6", ["c2", "1", "c1^2", "c1*Delta", "Delta*F", "c1/B"],
         [c2, mpmath.mpf(1), c1_val**2, c1_val*Delta, Delta*F, c1_val/B])

# Basis 7: {c2, 1, alpha, alpha^2, pi^2*alpha} -- alpha itself
run_pslq("7", ["c2", "1", "alpha", "alpha^2", "pi^2*alpha"],
         [c2, mpmath.mpf(1), alpha, alpha**2, mpmath.pi**2*alpha])

# Basis 8: {c2, 1, R/6, (R/6)^2, Delta*R/6} where R=6
run_pslq("8", ["c2", "1", "R/6", "(R/6)^2", "Delta*R/6"],
         [c2, mpmath.mpf(1), R_scalar/6, (R_scalar/6)**2, Delta*R_scalar/6])

# =====================================================================
# Basis 9: Brute force rational p/q scan
# =====================================================================

print(f"\n  --- Basis 9: Brute force rational scan p/q, q <= 200 ---")
c2_float = float(c2)
best_rationals = []
for q in range(1, 201):
    p = round(c2_float * q)
    if p <= 0:
        continue
    g = gcd(p, q)
    if g > 1:
        continue
    approx = p / q
    err = abs(approx - c2_float)
    if err < 1e-3:
        best_rationals.append((p, q, approx, err))

best_rationals.sort(key=lambda x: x[3])
print(f"    Best rational approximants (err < 1e-3):")
for p, q, approx, err in best_rationals[:15]:
    rel_pct = err / c2_float * 100
    marker = " <--" if rel_pct < 0.1 else ""
    print(f"    {p}/{q} = {approx:.10f}  err = {err:.2e}  ({rel_pct:.4f}%){marker}")

results["9_rational_scan"] = {
    'best_rationals': [(p, q, approx, err) for p, q, approx, err in best_rationals[:15]]
}

# =====================================================================
# Basis 10: Exhaustive c2 = a + b*Delta + c*c1 + d*F for |coeffs| <= 20
# =====================================================================

print(f"\n  --- Basis 10: Exhaust c2 = a + b*Delta + c*c1 + d*F, |coeffs| <= 20 ---")
print(f"    (searching ~68,921 combinations...)")

best_linear = []
c2_target = float(c2)
Delta_f = float(Delta)
c1_f = float(c1_val)
F_f = float(F)

for a_num in range(-20, 21):
    for a_den in range(1, 21):
        a = a_num / a_den
        for b in range(-20, 21):
            for c_coeff in range(-20, 21):
                # Solve for d: c2 = a + b*Delta + c*c1 + d*F
                # d = (c2 - a - b*Delta - c*c1) / F
                remainder = c2_target - a - b * Delta_f - c_coeff * c1_f
                d_exact = remainder / F_f
                d_round = round(d_exact)
                if abs(d_round) > 20:
                    continue
                pred = a + b * Delta_f + c_coeff * c1_f + d_round * F_f
                err = abs(pred - c2_target)
                if err < 5e-4:
                    best_linear.append((a_num, a_den, b, c_coeff, d_round, pred, err))

best_linear.sort(key=lambda x: x[6])
print(f"    Found {len(best_linear)} combinations with err < 5e-4")
if best_linear:
    print(f"    Top 10:")
    for a_num, a_den, b, c_coeff, d, pred, err in best_linear[:10]:
        a_str = f"{a_num}/{a_den}" if a_den > 1 else f"{a_num}"
        rel_pct = err / c2_target * 100
        print(f"      c2 = {a_str} + {b}*Delta + {c_coeff}*c1 + {d}*F"
              f"  = {pred:.10f}  err = {err:.2e}  ({rel_pct:.4f}%)")

results["10_linear_scan"] = {
    'count': len(best_linear),
    'best': [(a_num, a_den, b, c_coeff, d, pred, err)
             for a_num, a_den, b, c_coeff, d, pred, err in best_linear[:10]]
}

# =====================================================================
# ADDITIONAL PSLQ SEARCHES -- extended invariants
# =====================================================================

print("\n" + "=" * 70)
print("  ADDITIONAL PSLQ SEARCHES")
print("=" * 70)

# 11: Against heat-kernel coefficients (known on S^3)
a0_SD = mpmath.sqrt(mpmath.pi)  # Seeley-DeWitt a_0 on unit S^3
a1_SD = mpmath.sqrt(mpmath.pi)  # a_1
a2_SD = mpmath.sqrt(mpmath.pi) / 8  # a_2
run_pslq("11_heat_kernel", ["c2", "1", "a0_SD^2", "a2_SD", "1/a0_SD^2"],
         [c2, mpmath.mpf(1), a0_SD**2, a2_SD, 1/a0_SD**2])

# 12: Against zeta values (Paper 18 taxonomy)
run_pslq("12_zeta", ["c2", "1", "pi^2", "zeta(3)", "G_catalan"],
         [c2, mpmath.mpf(1), mpmath.pi**2, mpmath.zeta(3), mpmath.catalan])

# 13: Against Dirac-specific invariants
run_pslq("13_dirac", ["c2", "1", "1/g3", "1/|lam3|", "1/N2", "1/(g3*|lam3|)"],
         [c2, mpmath.mpf(1), 1/g3_dirac, 1/lam3, 1/N2, 1/(g3_dirac*lam3)])

# 14: Against c1 and its powers with R
run_pslq("14_c1_R", ["c2", "1", "c1^2", "(R/12)^2", "c1*R/12"],
         [c2, mpmath.mpf(1), c1_val**2, (R_scalar/12)**2, c1_val*R_scalar/12])

# 15: Cross-products of Paper 2 invariants
run_pslq("15_cross", ["c2", "1", "B*Delta", "F*Delta", "B*Delta^2", "F/B"],
         [c2, mpmath.mpf(1), B*Delta, F*Delta, B*Delta**2, F/B])

# 16: c2 against K and its parts
run_pslq("16_K_parts", ["c2", "1", "1/K", "pi/K", "Delta/K", "B/K"],
         [c2, mpmath.mpf(1), 1/K, mpmath.pi/K, Delta/K, B/K])

# 17: Wider rational + transcendental basis
run_pslq("17_wide", ["c2", "1", "pi", "pi^2", "log(2)", "1/pi"],
         [c2, mpmath.mpf(1), mpmath.pi, mpmath.pi**2, mpmath.log(2), 1/mpmath.pi])

# 18: Try B-related fractions specifically
run_pslq("18_B_fracs", ["c2", "1", "1/B", "1/B^2", "Delta/B", "F/B^2"],
         [c2, mpmath.mpf(1), 1/B, 1/B**2, Delta/B, F/B**2])

# 19: Parker-Toms progression: c1 = R/12, c2 = ? * (R/12)^2
pt_ratio = c2 / (c1_val**2)
print(f"\n  --- Ratio c2/c1^2 = {mpmath.nstr(pt_ratio, 15)} ---")
print(f"  If c2 = f * c1^2 = f/4, then f = {mpmath.nstr(pt_ratio, 15)}")
run_pslq("19_c1_sq", ["c2/c1^2", "1", "pi^2", "1/B", "Delta"],
         [pt_ratio, mpmath.mpf(1), mpmath.pi**2, 1/B, Delta])

# 20: Delta and its neighbors
run_pslq("20_Delta_ext", ["c2", "1", "Delta", "Delta*c1", "Delta*dim_S3", "Delta*R/12"],
         [c2, mpmath.mpf(1), Delta, Delta*c1_val, Delta*dim_S3, Delta*R_scalar/12])

# =====================================================================
# DEEPER ANALYSIS: factor structure of c2
# =====================================================================

print("\n" + "=" * 70)
print("  DEEPER ANALYSIS")
print("=" * 70)

# How does c2 decompose relative to known quantities?
print(f"\n  Key ratios:")
print(f"    c2 / c1          = {mpmath.nstr(c2/c1_val, 12)}")
print(f"    c2 / c1^2        = {mpmath.nstr(c2/c1_val**2, 12)}")
print(f"    c2 / Delta       = {mpmath.nstr(c2/Delta, 12)}")
print(f"    c2 / (c1*Delta)  = {mpmath.nstr(c2/(c1_val*Delta), 12)}")
print(f"    c2 / F           = {mpmath.nstr(c2/F, 12)}")
print(f"    c2 / (F*Delta)   = {mpmath.nstr(c2/(F*Delta), 12)}")
print(f"    c2 / alpha       = {mpmath.nstr(c2/alpha, 12)}")
print(f"    c2 * B           = {mpmath.nstr(c2*B, 12)}")
print(f"    c2 * K           = {mpmath.nstr(c2*K, 12)}")
print(f"    c2 * 40          = {mpmath.nstr(c2*40, 12)}")
print(f"    c2 * 120         = {mpmath.nstr(c2*120, 12)}")
print(f"    c2 * 240         = {mpmath.nstr(c2*240, 12)}")
print(f"    c2 * 360         = {mpmath.nstr(c2*360, 12)}")
print(f"    c2 * 720         = {mpmath.nstr(c2*720, 12)}")
print(f"    c2 * 5040        = {mpmath.nstr(c2*5040, 12)}")
print(f"    c2 / (1/(4*pi^2))= {mpmath.nstr(c2/(1/(4*mpmath.pi**2)), 12)}")

# Test c2 * common denominators for near-integer results
print(f"\n  c2 * n for small n (looking for near-integers):")
for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 18, 20,
          24, 25, 30, 32, 36, 40, 42, 48, 50, 60, 72, 80, 84,
          90, 100, 120, 140, 144, 150, 160, 168, 180, 200, 210,
          240, 252, 280, 300, 336, 360, 420, 480, 504, 560, 600,
          720, 840, 1000, 1008, 1260, 1680, 2520, 5040]:
    val = float(c2) * n
    frac_part = val - round(val)
    if abs(frac_part) < 0.05:
        print(f"    c2 * {n:5d} = {val:12.6f}  (frac = {frac_part:+.6f}) <-- NEAR INTEGER")

# =====================================================================
# PSLQ on c2 * 40 (= c2 / Delta) against simple basis
# =====================================================================

c2_over_delta = c2 / Delta
print(f"\n  c2 / Delta = c2 * 40 = {mpmath.nstr(c2_over_delta, 12)}")
run_pslq("21_c2_over_delta", ["c2/Delta", "1", "pi^2", "pi", "log(2)"],
         [c2_over_delta, mpmath.mpf(1), mpmath.pi**2, mpmath.pi, mpmath.log(2)])

# c2 * B
c2_times_B = c2 * B
print(f"\n  c2 * B = c2 * 42 = {mpmath.nstr(c2_times_B, 12)}")
run_pslq("22_c2_times_B", ["c2*B", "1", "pi^2", "pi", "F"],
         [c2_times_B, mpmath.mpf(1), mpmath.pi**2, mpmath.pi, F])

# c2 * K
c2_times_K = c2 * K
print(f"\n  c2 * K = c2 * 137.036 = {mpmath.nstr(c2_times_K, 12)}")
run_pslq("23_c2_times_K", ["c2*K", "1", "pi^2", "pi", "B"],
         [c2_times_K, mpmath.mpf(1), mpmath.pi**2, mpmath.pi, B])

# =====================================================================
# PSLQ: does c2 live in the ring Q(pi^2)?
# =====================================================================

print(f"\n  --- Is c2 in Q(pi^2)? ---")
run_pslq("24_Qpi2", ["c2", "1", "pi^2", "pi^4"],
         [c2, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4],
         maxc=5000)

# PSLQ: does c2 live in Q(pi)?
run_pslq("25_Qpi", ["c2", "1", "pi", "pi^2", "pi^3"],
         [c2, mpmath.mpf(1), mpmath.pi, mpmath.pi**2, mpmath.pi**3],
         maxc=5000)

# =====================================================================
# Summary of all PSLQ results
# =====================================================================

print("\n" + "=" * 70)
print("  SUMMARY OF ALL RESULTS")
print("=" * 70)

c2_float = float(c2)
print(f"\n  c2 = {c2_float:.10f} (7 reliable digits)")
print(f"  c2 = {mpmath.nstr(c2, 15)}")

hits = {k: v for k, v in results.items() if isinstance(v, dict) and v.get('hit')}
misses = {k: v for k, v in results.items() if isinstance(v, dict) and not v.get('hit') and 'hit' in v}

if hits:
    print(f"\n  PSLQ HITS ({len(hits)}):")
    for name, r in sorted(hits.items()):
        print(f"    Basis {name}: {r['formula']}")
        print(f"      c2_pred = {r['c2_predicted']:.10f}  resid = {r['residual']:.2e}")
else:
    print(f"\n  NO PSLQ HITS across all {len(results)} basis sets")
    print(f"  (This is consistent with c2 being a new transcendental or")
    print(f"   the precision being insufficient for PSLQ identification)")

if misses:
    print(f"\n  PSLQ MISSES ({len(misses)}):")
    for name in sorted(misses.keys()):
        print(f"    Basis {name}")

# Report on specific candidates
print(f"\n  CLOSEST SPECIFIC CANDIDATES:")
close_candidates = []
for label, val in candidates:
    err = abs(float(val) - c2_float)
    rel = err / c2_float * 100
    close_candidates.append((label, float(val), err, rel))
close_candidates.sort(key=lambda x: x[2])
for label, val, err, rel in close_candidates[:5]:
    print(f"    {label:>50s}: {val:.10f}  err={err:.2e}  ({rel:.3f}%)")

# Rational scan results
if results.get("9_rational_scan", {}).get('best_rationals'):
    rats = results["9_rational_scan"]['best_rationals']
    print(f"\n  BEST RATIONAL APPROXIMANTS:")
    for p, q, approx, err in rats[:5]:
        rel = err / c2_float * 100
        print(f"    {p}/{q} = {approx:.10f}  err={err:.2e}  ({rel:.4f}%)")

# Honest assessment
print(f"""
  ASSESSMENT:
  c2 ~ {c2_float:.6f} at ~7 digit precision.

  The apparent value c2 = (delta - c1/lam^2) * lam^4 is CONTAMINATED
  by c3/lam^2 + c4/lam^4 + ... at the ~16% level (if c3 ~ c2).

  The TRUE c2 could differ from the apparent value by O(10^-2).

  At 7-digit precision, PSLQ can identify rational numbers with
  denominator up to ~1000, but NOT transcendental combinations
  (which need 30+ digits).

  NEXT STEPS to resolve c2:
  1. Compute F2/Schwinger at n_ext=2,3,4 with sufficient n_int
     to get multi-point data for a 3-parameter fit (c1, c2, c3)
  2. The n_int=1 exact forms at n_ext=1..7 give a PARTIAL view
     (99% of total at n_ext=1, but less dominant at higher n_ext)
  3. Richardson extrapolation on multi-n_ext data would yield
     c2 at 10+ digit precision, enabling definitive PSLQ
""")

# =====================================================================
# Save results
# =====================================================================

output = {
    'c2_apparent': float(c2),
    'c2_v1': float(c2_v1),
    'c2_v2': float(c2_v2),
    'delta_v2': float(delta_v2),
    'delta_v1': float(delta_v1),
    'c1': 0.5,
    'precision_digits': 7,
    'paper2_invariants': {
        'B': 42,
        'F': float(F),
        'Delta': 0.025,
        'K': float(K),
        'alpha': float(alpha),
    },
    'specific_candidates': [
        {'label': label, 'value': float(val),
         'error': abs(float(val) - float(c2)),
         'rel_pct': abs(float(val) - float(c2)) / abs(float(c2)) * 100}
        for label, val in candidates
    ],
    'pslq_results': {},
    'rational_scan': results.get("9_rational_scan", {}),
    'linear_scan': results.get("10_linear_scan", {}),
}

# Convert PSLQ results (skip scan results already handled)
for name, r in results.items():
    if name.startswith("9_") or name.startswith("10_"):
        continue
    if isinstance(r, dict):
        pslq_entry = {'hit': r.get('hit', False)}
        if r.get('hit'):
            pslq_entry['relation'] = r['relation']
            pslq_entry['formula'] = r['formula']
            pslq_entry['c2_predicted'] = r['c2_predicted']
            pslq_entry['residual'] = r['residual']
        output['pslq_results'][f'basis_{name}'] = pslq_entry

os.makedirs('debug/data', exist_ok=True)
with open('debug/data/g2_c2_paper2_pslq.json', 'w') as f:
    json.dump(output, f, indent=2, default=float)

print(f"\n  Results saved to debug/data/g2_c2_paper2_pslq.json")
print(f"  Script: debug/g2_c2_paper2_pslq.py")
