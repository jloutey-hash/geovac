"""
Phase 5: Systematic residual analysis of the Hopf bundle alpha formula.

Formula: alpha^3 - K*alpha + 1 = 0, K = pi*(42 + pi^2/6 - 1/40)
Residual: alpha_CODATA - alpha_formula ~ 8.8e-8 relative

Tests whether the residual matches:
  1. QED radiative corrections (alpha/pi, alpha^2, etc.)
  2. Spectral invariant corrections from the Hopf bundle
  3. Mathematical constants (brute force search)
  4. Self-consistency (what correction to K closes the gap?)
  5. CODATA vintage comparison (2014 vs 2018 vs 2022)
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from mpmath import mp, mpf, pi as mpi, zeta, log, sqrt, fabs, nstr, polyroots, power

# Set precision to 80 decimal digits
mp.dps = 80

# ==============================================================================
# CONSTANTS (high precision)
# ==============================================================================
B = mpf(42)
F = mpi**2 / 6  # zeta(2)
Delta = mpf(1) / 40
K = mpi * (B + F - Delta)

# CODATA values (alpha^{-1})
CODATA_2014_inv = mpf('137.035999139')
CODATA_2014_unc = mpf('0.000000031')
CODATA_2018_inv = mpf('137.035999084')
CODATA_2018_unc = mpf('0.000000021')
# CODATA 2022 (TAMU Cs measurement): published 2024
CODATA_2022_inv = mpf('137.035999177')
CODATA_2022_unc = mpf('0.000000013')

print("=" * 78)
print("PHASE 5: SYSTEMATIC RESIDUAL ANALYSIS")
print("=" * 78)

# ==============================================================================
# TASK 1: EXACT RESIDUAL
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 1: EXACT RESIDUAL (50+ digit precision)")
print("=" * 78)

print(f"\nK = pi*(42 + pi^2/6 - 1/40)")
print(f"K = {nstr(K, 50)}")

# Solve cubic alpha^3 - K*alpha + 1 = 0
# Using mpmath polyroots for high precision
# Polynomial: x^3 + 0*x^2 - K*x + 1
roots = polyroots([1, 0, -K, 1])

# Sort real roots
real_roots = sorted([r.real for r in roots if fabs(r.imag) < mpf(10)**(-70)])
assert len(real_roots) == 3, f"Expected 3 real roots, got {len(real_roots)}"

alpha_formula = real_roots[1]  # smallest positive root
alpha_inv_formula = 1 / alpha_formula

print(f"\nThree roots of alpha^3 - K*alpha + 1 = 0:")
print(f"  r1 = {nstr(real_roots[0], 50)}")
print(f"  r2 = {nstr(real_roots[1], 50)}  <-- alpha_formula")
print(f"  r3 = {nstr(real_roots[2], 50)}")

# Verify
residual_check = alpha_formula**3 - K * alpha_formula + 1
print(f"\nVerification: alpha^3 - K*alpha + 1 = {nstr(residual_check, 10)}")

print(f"\nalpha_formula     = {nstr(alpha_formula, 50)}")
print(f"alpha_formula^-1  = {nstr(alpha_inv_formula, 50)}")

# Compute residuals against all CODATA values
for label, codata_inv, unc in [
    ("CODATA 2014", CODATA_2014_inv, CODATA_2014_unc),
    ("CODATA 2018", CODATA_2018_inv, CODATA_2018_unc),
    ("CODATA 2022", CODATA_2022_inv, CODATA_2022_unc),
]:
    alpha_codata = 1 / codata_inv
    delta_alpha = alpha_codata - alpha_formula       # shift in alpha
    delta_inv = codata_inv - alpha_inv_formula        # shift in alpha^{-1}
    rel_err = fabs(delta_alpha) / alpha_codata
    sigma = fabs(delta_inv) / unc                     # how many sigma away

    print(f"\n--- {label}: alpha^{{-1}} = {nstr(codata_inv, 15)} +/- {nstr(unc, 2)} ---")
    print(f"  delta(alpha)     = {nstr(delta_alpha, 20)}  (signed)")
    print(f"  delta(alpha^-1)  = {nstr(delta_inv, 20)}  (signed)")
    print(f"  |relative error| = {nstr(rel_err, 10)}")
    print(f"  sigma distance   = {nstr(sigma, 4)}")

# Store key values for later use
alpha_codata_2018 = 1 / CODATA_2018_inv
delta_alpha = alpha_codata_2018 - alpha_formula
delta_inv = CODATA_2018_inv - alpha_inv_formula
rel_err = fabs(delta_alpha) / alpha_codata_2018
alpha_f = alpha_formula  # shorthand

print(f"\n*** PRIMARY RESIDUAL (CODATA 2018) ***")
print(f"  delta(alpha^-1) = {nstr(delta_inv, 30)}")
print(f"  This is the number we are trying to match.")

# ==============================================================================
# TASK 2: QED RADIATIVE CORRECTIONS
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 2: QED RADIATIVE CORRECTIONS")
print("=" * 78)

alpha_phys = alpha_codata_2018  # use physical alpha for QED terms

qed_candidates = {
    "alpha/(2*pi)  [Schwinger]":        alpha_phys / (2 * mpi),
    "alpha^2/(2*pi)":                   alpha_phys**2 / (2 * mpi),
    "alpha/pi":                         alpha_phys / mpi,
    "alpha^2/pi":                       alpha_phys**2 / mpi,
    "alpha/(3*pi)  [VP coeff]":         alpha_phys / (3 * mpi),
    "alpha^2/(3*pi)":                   alpha_phys**2 / (3 * mpi),
    "(alpha/pi)^2":                     (alpha_phys / mpi)**2,
    "alpha^2*ln(alpha^-1)":             alpha_phys**2 * log(1 / alpha_phys),
    "alpha^3":                          alpha_phys**3,
    "2*alpha^3":                        2 * alpha_phys**3,
}

print(f"\n{'Candidate':<35} {'Value':>20} {'delta_alpha':>20} {'delta_inv':>20} {'Match(alpha)':>14} {'Match(inv)':>14}")
print("-" * 125)

for name, val in qed_candidates.items():
    # Compare against delta_alpha (shift in alpha)
    if fabs(delta_alpha) > 0:
        match_alpha = val / fabs(delta_alpha)
    else:
        match_alpha = mpf(0)

    # Compare against delta_inv (shift in alpha^{-1})
    if fabs(delta_inv) > 0:
        match_inv = val / fabs(delta_inv)
    else:
        match_inv = mpf(0)

    print(f"  {name:<33} {nstr(val, 10):>20} {nstr(fabs(delta_alpha), 10):>20} "
          f"{nstr(fabs(delta_inv), 10):>20} {nstr(match_alpha, 8):>14} {nstr(match_inv, 8):>14}")

# Also compute delta_inv in units of these QED quantities
print(f"\n--- Ratios: delta_inv / QED_quantity ---")
for name, val in qed_candidates.items():
    ratio = delta_inv / val if val != 0 else mpf(0)
    print(f"  delta_inv / ({name}) = {nstr(ratio, 15)}")


# ==============================================================================
# TASK 3: SPECTRAL INVARIANT CORRECTIONS
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 3: SPECTRAL INVARIANT CORRECTIONS")
print("=" * 78)

B4 = mpf(162)  # B(4) = next Casimir trace after selection cutoff
zeta3 = zeta(3)

spectral_candidates = {
    "1/(K*B) = 1/(137*42)":            1 / (K * B),
    "Delta/K = (1/40)/137":             Delta / K,
    "zeta(3)/(2*pi^2) / K":            zeta3 / (2 * mpi**2) / K,
    "Delta^2 = (1/40)^2":              Delta**2,
    "F*Delta = (pi^2/6)*(1/40)":        F * Delta,
    "1/B^2 = 1/1764":                  1 / B**2,
    "pi/(K*B)":                        mpi / (K * B),
    "Delta/B = (1/40)/42":             Delta / B,
    "1/K^2":                           1 / K**2,
    "alpha_formula*Delta":             alpha_f * Delta,
    "B4/K = 162/137":                  B4 / K,
}

print(f"\n{'Candidate':<35} {'Value':>22} {'Ratio to delta_inv':>22} {'% match':>12}")
print("-" * 95)

for name, val in spectral_candidates.items():
    ratio = delta_inv / val if val != 0 else mpf(0)
    # "match" = how close ratio is to a simple integer or simple fraction
    print(f"  {name:<33} {nstr(val, 14):>22} {nstr(ratio, 14):>22}")

# Check direct matches
print(f"\n--- Direct match quality (|candidate - delta_inv| / |delta_inv|) ---")
for name, val in spectral_candidates.items():
    rel_match = fabs(val - fabs(delta_inv)) / fabs(delta_inv)
    print(f"  {name:<33} relative deviation = {nstr(rel_match, 8)}")


# ==============================================================================
# TASK 4: MATHEMATICAL CONSTANTS BRUTE FORCE
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 4: MATHEMATICAL CONSTANTS SEARCH")
print("=" * 78)

target = fabs(delta_inv)
print(f"\nTarget: |delta_inv| = {nstr(target, 30)}")

best_matches = []  # (match_quality, description, value)

# 4a: Simple fractions 1/n
print("\n--- 4a: Simple fractions 1/n ---")
for n in range(1, 100001):
    val = mpf(1) / n
    rel = fabs(val - target) / target
    if rel < 0.001:  # 0.1% match
        best_matches.append((float(rel), f"1/{n}", val))

# 4b: pi^k / n
print("--- 4b: pi^k / n ---")
for k in range(-3, 5):
    if k == 0:
        continue
    pi_k = mpi**k
    # n ~ pi^k / target
    n_approx = float(pi_k / target)
    for n in range(max(1, int(n_approx) - 2), int(n_approx) + 3):
        if n <= 0:
            continue
        val = pi_k / n
        rel = fabs(val - target) / target
        if rel < 0.001:
            best_matches.append((float(rel), f"pi^{k}/{n}", val))

# 4c: zeta(k) / n
print("--- 4c: zeta(k) / n ---")
for k in [2, 3, 4, 5]:
    zk = zeta(k)
    n_approx = float(zk / target)
    for n in range(max(1, int(n_approx) - 2), int(n_approx) + 3):
        if n <= 0:
            continue
        val = zk / n
        rel = fabs(val - target) / target
        if rel < 0.001:
            best_matches.append((float(rel), f"zeta({k})/{n}", val))

# 4d: ln(x) / n
print("--- 4d: ln(x) / n ---")
for base_name, base_val in [("2", log(2)), ("3", log(3)), ("pi", log(mpi))]:
    n_approx = float(base_val / target)
    for n in range(max(1, int(n_approx) - 2), int(n_approx) + 3):
        if n <= 0:
            continue
        val = base_val / n
        rel = fabs(val - target) / target
        if rel < 0.001:
            best_matches.append((float(rel), f"ln({base_name})/{n}", val))

# 4e: Combinations: pi^j * zeta(k) / n, etc.
print("--- 4e: Combinations ---")
for k in [2, 3]:
    zk = zeta(k)
    for j in [-2, -1, 1, 2]:
        combo = mpi**j * zk
        n_approx = float(combo / target)
        for n in range(max(1, int(n_approx) - 2), int(n_approx) + 3):
            if n <= 0:
                continue
            val = combo / n
            rel = fabs(val - target) / target
            if rel < 0.001:
                best_matches.append((float(rel), f"pi^{j}*zeta({k})/{n}", val))

# Also try: alpha_phys^k / n for small k
for k in [1, 2, 3]:
    ak = alpha_phys**k
    n_approx = float(ak / target)
    for n in range(max(1, int(n_approx) - 5), int(n_approx) + 6):
        if n <= 0:
            continue
        val = ak / n
        rel = fabs(val - target) / target
        if rel < 0.001:
            best_matches.append((float(rel), f"alpha^{k}/{n}", val))

# Sort and report
best_matches.sort(key=lambda x: x[0])
print(f"\n*** TOP 30 MATCHES (within 0.1% of |delta_inv|) ***")
print(f"{'Match %':>10} {'Expression':<30} {'Value':>25}")
print("-" * 68)
for rel, desc, val in best_matches[:30]:
    print(f"  {rel*100:.6f}% {desc:<30} {nstr(val, 18):>25}")


# ==============================================================================
# TASK 5: SELF-CONSISTENCY (epsilon correction to K)
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 5: SELF-CONSISTENCY — EPSILON CORRECTION TO K")
print("=" * 78)

# Find epsilon such that K_corrected = K + epsilon gives alpha = alpha_CODATA
# From the cubic: alpha^3 - (K+eps)*alpha + 1 = 0
# => eps = (alpha^3 - K*alpha + 1) / alpha  evaluated at alpha = alpha_CODATA
# But alpha_CODATA^3 - K*alpha_CODATA + 1 != 0 (it has a residual)

cubic_resid_at_codata = alpha_codata_2018**3 - K * alpha_codata_2018 + 1
epsilon = -cubic_resid_at_codata / alpha_codata_2018
# Check: (K+eps)*alpha = alpha^3 + 1, so eps = (alpha^3 + 1)/alpha - K

K_corrected = K + epsilon
print(f"\nK_original  = {nstr(K, 30)}")
print(f"K_corrected = {nstr(K_corrected, 30)}")
print(f"epsilon     = {nstr(epsilon, 30)}")
print(f"epsilon/K   = {nstr(epsilon / K, 15)}")

# Verify
roots_corrected = polyroots([1, 0, -K_corrected, 1])
real_corrected = sorted([r.real for r in roots_corrected if fabs(r.imag) < mpf(10)**(-70)])
print(f"\nVerification: smallest positive root of corrected cubic:")
print(f"  alpha_corrected^-1 = {nstr(1/real_corrected[1], 20)}")
print(f"  alpha_CODATA^-1    = {nstr(CODATA_2018_inv, 20)}")

# Test epsilon against spectral quantities
print(f"\n--- Does epsilon match any spectral quantity? ---")
epsilon_candidates = {
    "Delta^2 = 1/1600":                   Delta**2,
    "F*Delta = pi^2/240":                  F * Delta,
    "1/B = 1/42":                          1 / B,
    "zeta(3)/K":                           zeta3 / K,
    "spectral det gap (0.04276)":          mpf('0.04276'),
    "1/K":                                 1 / K,
    "Delta/pi":                            Delta / mpi,
    "1/(B*pi)":                            1 / (B * mpi),
    "alpha^2/pi":                          alpha_phys**2 / mpi,
    "alpha/(2*pi) [Schwinger]":            alpha_phys / (2 * mpi),
    "pi/K^2":                              mpi / K**2,
    "Delta/K":                             Delta / K,
    "1/B^2":                               1 / B**2,
}

print(f"\n{'Candidate':<35} {'Value':>22} {'epsilon/candidate':>22} {'Match quality':>15}")
print("-" * 97)
for name, val in epsilon_candidates.items():
    ratio = epsilon / val if val != 0 else mpf(0)
    rel_match = fabs(fabs(ratio) - 1) if fabs(ratio) > 0.5 and fabs(ratio) < 2 else fabs(ratio)
    print(f"  {name:<33} {nstr(val, 14):>22} {nstr(ratio, 14):>22} {nstr(rel_match, 8):>15}")

# Also: is epsilon close to a simple fraction?
print(f"\n--- epsilon as simple fraction ---")
inv_eps = 1 / epsilon
print(f"  1/epsilon = {nstr(inv_eps, 20)}")
print(f"  Nearest integers: {int(float(inv_eps))-1}, {int(float(inv_eps))}, {int(float(inv_eps))+1}")


# ==============================================================================
# TASK 6: CODATA VINTAGE COMPARISON
# ==============================================================================
print("\n" + "=" * 78)
print("TASK 6: CODATA VINTAGE COMPARISON")
print("=" * 78)

print(f"\nalpha_formula^-1 = {nstr(alpha_inv_formula, 20)}")
print(f"\n{'Vintage':<15} {'alpha^-1':>20} {'unc':>15} {'delta_inv':>20} {'sigma':>10} {'rel_err':>15}")
print("-" * 97)

for label, codata_inv, unc in [
    ("CODATA 2014", CODATA_2014_inv, CODATA_2014_unc),
    ("CODATA 2018", CODATA_2018_inv, CODATA_2018_unc),
    ("CODATA 2022", CODATA_2022_inv, CODATA_2022_unc),
]:
    dinv = codata_inv - alpha_inv_formula
    sigma = fabs(dinv) / unc
    rel = fabs(dinv) / codata_inv
    print(f"  {label:<13} {nstr(codata_inv, 15):>20} {nstr(unc, 2):>15} "
          f"{nstr(dinv, 15):>20} {nstr(sigma, 5):>10} {nstr(rel, 8):>15}")

print(f"\n--- Trend analysis ---")
d14 = CODATA_2014_inv - alpha_inv_formula
d18 = CODATA_2018_inv - alpha_inv_formula
d22 = CODATA_2022_inv - alpha_inv_formula
print(f"  2014 delta_inv = {nstr(d14, 15)}")
print(f"  2018 delta_inv = {nstr(d18, 15)}")
print(f"  2022 delta_inv = {nstr(d22, 15)}")
if d14 > d18 > 0 or d14 < d18 < 0:
    if fabs(d22) < fabs(d18):
        print("  TREND: residual DECREASING (measurements approaching formula) -- 2014->2018")
        if fabs(d22) > fabs(d18):
            print("  BUT: 2022 reverses trend (residual increasing)")
    else:
        print("  NOTE: 2022 value moves AWAY from formula")
else:
    print("  NOTE: non-monotonic movement")

if d22 > d18:
    print(f"  2022 moves AWAY from formula by {nstr(d22 - d18, 10)}")
else:
    print(f"  2022 moves TOWARD formula by {nstr(d18 - d22, 10)}")


# ==============================================================================
# SUMMARY
# ==============================================================================
print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)

print(f"""
Key findings:

1. EXACT RESIDUAL (CODATA 2018):
   delta(alpha^-1) = {nstr(delta_inv, 20)}
   relative error  = {nstr(rel_err, 10)}

2. QED CORRECTIONS: The residual is compared against standard QED expansion
   terms. Look at the ratios above to see if any are close to unity or
   simple integers.

3. SPECTRAL INVARIANTS: The residual is compared against combinations of
   B, F, Delta, K from the Hopf bundle formula.

4. MATHEMATICAL CONSTANTS: Top matches from brute-force search are listed.

5. EPSILON CORRECTION: K needs a correction of epsilon = {nstr(epsilon, 15)}
   to match CODATA 2018 exactly.

6. CODATA TREND: The three vintages show the direction of measurement drift
   relative to the formula prediction.
""")
