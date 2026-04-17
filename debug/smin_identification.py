"""Identify the irreducible CG-weighted sunset sum S_min.

S_min = sum_{k=1}^inf T(k)^2

where T(k) = 2*hurwitz(2, k+3/2) - (1/2)*hurwitz(4, k+3/2)

is the Dirac Dirichlet series tail starting at level k.  The original
18-element PSLQ basis failed.  This script tries extended bases including
polylogarithms, Euler sums, multiple zeta values, and cross-products.

Usage:
    python debug/smin_identification.py
"""

import sys
import os
import json
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mpmath

# Set very high precision for PSLQ
mpmath.mp.dps = 150

# ---------------------------------------------------------------------------
# Compute S_min at high precision
# ---------------------------------------------------------------------------

def compute_T_k(k, s=4):
    """T(k) = 2*hurwitz(s-2, k+3/2) - (1/2)*hurwitz(s, k+3/2)."""
    a = mpmath.mpf(k) + mpmath.mpf(3) / 2
    return 2 * mpmath.hurwitz(s - 2, a) - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, a)


def compute_S_min_with_tail(n_terms=5000, s=4):
    """Compute S_min = sum_{k=1}^n_terms T(k)^2, with tail correction."""
    S = mpmath.mpf(0)
    for k in range(1, n_terms + 1):
        S += compute_T_k(k, s)**2

    # Tail estimate: T(k) ~ 2/(k+3/2)^2 for large k
    # T(k)^2 ~ 4/(k+3/2)^4 - 2/(k+3/2)^6 + 1/4*(k+3/2)^8
    a_tail = mpmath.mpf(n_terms + 1) + mpmath.mpf(3) / 2
    tail = (4 * mpmath.hurwitz(4, a_tail)
            - 2 * mpmath.hurwitz(6, a_tail)
            + mpmath.mpf(1) / 4 * mpmath.hurwitz(8, a_tail))
    return S, S + tail, tail


print("=" * 70)
print("S_min IDENTIFICATION ATTEMPT")
print("=" * 70)
print("Working precision: %d digits" % mpmath.mp.dps)
print()

# Step 1: Compute S_min
print("Step 1: Computing S_min at high precision...")
t0 = time.time()
S_raw, S_corrected, tail = compute_S_min_with_tail(n_terms=5000, s=4)
t1 = time.time()
print("  Raw sum (5000 terms): %s" % mpmath.nstr(S_raw, 60))
print("  Tail correction:      %s" % mpmath.nstr(tail, 30))
print("  Corrected S_min:      %s" % mpmath.nstr(S_corrected, 60))
print("  Time: %.1fs" % (t1 - t0))

print("\n  Convergence check with 10000 terms...")
t0 = time.time()
S_raw2, S_corrected2, tail2 = compute_S_min_with_tail(n_terms=10000, s=4)
t1 = time.time()
print("  Raw sum (10000 terms): %s" % mpmath.nstr(S_raw2, 60))
print("  Corrected S_min:       %s" % mpmath.nstr(S_corrected2, 60))
print("  Diff(corrected):       %s" % mpmath.nstr(abs(S_corrected2 - S_corrected), 20))
print("  Time: %.1fs" % (t1 - t0))

S_min = S_corrected2
print("\n  Using S_min = %s" % mpmath.nstr(S_min, 80))

# ---------------------------------------------------------------------------
# Step 2: Compute basis constants
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Step 2: Computing basis constants...")

pi = mpmath.pi
pi2 = pi**2
pi4 = pi**4
pi6 = pi**6
pi8 = pi**8

z3 = mpmath.zeta(3)
z5 = mpmath.zeta(5)
z7 = mpmath.zeta(7)

G = mpmath.catalan
beta4 = (mpmath.hurwitz(4, mpmath.mpf(1)/4) - mpmath.hurwitz(4, mpmath.mpf(3)/4)) / 256
beta3 = pi**3 / 32

ln2 = mpmath.log(2)

Li2_half = mpmath.polylog(2, mpmath.mpf(1)/2)
Li3_half = mpmath.polylog(3, mpmath.mpf(1)/2)
Li4_half = mpmath.polylog(4, mpmath.mpf(1)/2)

gamma_em = mpmath.euler

# Hurwitz values at 3/2
hz_vals = {}
for s_val in [2, 3, 4, 5, 6]:
    hz_vals[s_val] = mpmath.hurwitz(s_val, mpmath.mpf(3)/2)

# D(4) from Hurwitz
D4 = 2 * hz_vals[2] - mpmath.mpf(1)/2 * hz_vals[4]

print("  D(4)    = %s" % mpmath.nstr(D4, 30))
print("  D(4)^2  = %s" % mpmath.nstr(D4**2, 30))
print("  S_min   = %s" % mpmath.nstr(S_min, 30))
print("  ratio   = %s" % mpmath.nstr(S_min / D4**2, 15))

# Clausen
Cl2_pi3 = mpmath.clsin(2, mpmath.pi/3)

# ---------------------------------------------------------------------------
# Step 3: PSLQ attempts
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Step 3: PSLQ identification attempts")
print("=" * 70)

results = {}

def try_pslq(name, basis_vals, basis_labels, tol=1e-50, maxcoeff=200000):
    """Try PSLQ with given basis, return result dict."""
    print("\n  Attempt: %s" % name)
    print("  Basis size: %d elements" % len(basis_labels))
    full_basis = [S_min] + list(basis_vals)
    full_labels = ["S_min"] + list(basis_labels)
    try:
        rel = mpmath.pslq(full_basis, tol=tol, maxcoeff=maxcoeff)
    except Exception as e:
        print("  PSLQ FAILED: %s" % e)
        results[name] = {"identified": False, "error": str(e)}
        return None

    if rel is None:
        print("  PSLQ returned None (no relation found)")
        results[name] = {"identified": False, "error": "no relation"}
        return None

    if rel[0] == 0:
        print("  PSLQ found relation but S_min coefficient is 0 (spurious)")
        results[name] = {"identified": False, "error": "zero leading coeff"}
        return None

    print("  PSLQ FOUND RELATION!")
    components = {}
    reconstructed = mpmath.mpf(0)
    for i in range(1, len(rel)):
        if rel[i] != 0:
            coeff = mpmath.mpf(-rel[i]) / mpmath.mpf(rel[0])
            label = full_labels[i]
            components[label] = "%d/%d" % (-rel[i], rel[0])
            reconstructed += coeff * full_basis[i]
            print("    %d/%d * %s" % (-rel[i], rel[0], label))

    residual = float(abs(S_min - reconstructed))
    print("  Residual: %.3e" % residual)

    result = {
        "identified": True,
        "relation": [int(r) for r in rel],
        "components": components,
        "residual": residual,
        "denominator": int(rel[0]),
    }
    results[name] = result
    return result


# --- Attempt 1: Original 18-element basis (reproduce failure) ---
try_pslq(
    "original_18",
    [mpmath.mpf(1), pi2, pi4, pi6, pi8, z3, z5, z7,
     pi2*z3, pi2*z5, G, beta4, pi2*G, pi2*beta4,
     ln2, pi2*ln2, z3**2, z3*z5],
    ["1", "pi2", "pi4", "pi6", "pi8", "z3", "z5", "z7",
     "pi2*z3", "pi2*z5", "G", "beta4", "pi2*G", "pi2*beta4",
     "ln2", "pi2*ln2", "z3^2", "z3*z5"]
)

# --- Attempt 2: Add polylogarithms ---
try_pslq(
    "with_polylog",
    [mpmath.mpf(1), pi2, pi4, pi6, z3, z5, G, beta4,
     pi2*G, pi2*beta4, Li2_half, Li3_half, Li4_half, ln2, ln2**2],
    ["1", "pi2", "pi4", "pi6", "z3", "z5", "G", "beta4",
     "pi2*G", "pi2*beta4", "Li2(1/2)", "Li3(1/2)", "Li4(1/2)", "ln2", "ln2^2"]
)

# --- Attempt 3: Add Euler-Mascheroni and cross products ---
try_pslq(
    "with_gamma_cross",
    [mpmath.mpf(1), pi2, pi4, z3, z5, G, beta4,
     G*z3, G**2, beta4*z3, G*beta4, gamma_em,
     pi2*z3, pi2*G, pi2*beta4, ln2],
    ["1", "pi2", "pi4", "z3", "z5", "G", "beta4",
     "G*z3", "G^2", "beta4*z3", "G*beta4", "gamma",
     "pi2*z3", "pi2*G", "pi2*beta4", "ln2"]
)

# --- Attempt 4: Add beta(3) and multiple zeta ---
try_pslq(
    "with_beta3_mzv",
    [mpmath.mpf(1), pi2, pi4, pi6, z3, z5, G, beta4, beta3,
     pi2*G, pi2*beta4, G**2, G*beta4, beta4**2,
     z3**2, z3*G, z3*beta4],
    ["1", "pi2", "pi4", "pi6", "z3", "z5", "G", "beta4", "beta3",
     "pi2*G", "pi2*beta4", "G^2", "G*beta4", "beta4^2",
     "z3^2", "z3*G", "z3*beta4"]
)

# --- Attempt 5: Hurwitz at 3/2 squares ---
try_pslq(
    "hurwitz_products",
    [mpmath.mpf(1), pi2, pi4, pi6, z3, z5,
     G, beta4,
     hz_vals[2], hz_vals[3], hz_vals[4], hz_vals[5], hz_vals[6],
     hz_vals[2]**2, hz_vals[2]*hz_vals[3], hz_vals[2]*hz_vals[4],
     hz_vals[3]**2, hz_vals[3]*hz_vals[4], hz_vals[4]**2],
    ["1", "pi2", "pi4", "pi6", "z3", "z5",
     "G", "beta4",
     "hz2", "hz3", "hz4", "hz5", "hz6",
     "hz2^2", "hz2*hz3", "hz2*hz4",
     "hz3^2", "hz3*hz4", "hz4^2"]
)

# --- Attempt 6: Decompose into three component sums ---
print("\n  Computing the three component sums A, B, C...")
sum_22 = mpmath.mpf(0)  # sum_k zeta(2,k+3/2)^2
sum_24 = mpmath.mpf(0)  # sum_k zeta(2,k+3/2)*zeta(4,k+3/2)
sum_44 = mpmath.mpf(0)  # sum_k zeta(4,k+3/2)^2

n_terms_comp = 10000
for k in range(1, n_terms_comp + 1):
    a = mpmath.mpf(k) + mpmath.mpf(3) / 2
    h2 = mpmath.hurwitz(2, a)
    h4 = mpmath.hurwitz(4, a)
    sum_22 += h2**2
    sum_24 += h2 * h4
    sum_44 += h4**2

print("  A = sum zeta(2,k+3/2)^2            = %s" % mpmath.nstr(sum_22, 40))
print("  B = sum zeta(2,k+3/2)*zeta(4,k+3/2) = %s" % mpmath.nstr(sum_24, 40))
print("  C = sum zeta(4,k+3/2)^2            = %s" % mpmath.nstr(sum_44, 40))
print("  S_min = 4A - 2B + C/4              = %s" % mpmath.nstr(4*sum_22 - 2*sum_24 + sum_44/4, 40))
print("  S_min (direct)                      = %s" % mpmath.nstr(S_min, 40))

# Try PSLQ on each component
for label, val in [("sum_hz2sq", sum_22), ("sum_hz2_hz4", sum_24), ("sum_hz4sq", sum_44)]:
    try_pslq(
        "component_%s" % label,
        [mpmath.mpf(1), pi2, pi4, pi6, pi8, z3, z5, z7,
         G, beta4, pi2*G, pi2*beta4,
         G**2, G*beta4, beta4**2,
         z3**2, z3*G, z3*beta4,
         ln2, pi2*z3, pi2*z5],
        ["1", "pi2", "pi4", "pi6", "pi8", "z3", "z5", "z7",
         "G", "beta4", "pi2*G", "pi2*beta4",
         "G^2", "G*beta4", "beta4^2",
         "z3^2", "z3*G", "z3*beta4",
         "ln2", "pi2*z3", "pi2*z5"],
    )

# --- Attempt 7: Polygamma / trigamma products ---
psi1_32 = mpmath.psi(1, mpmath.mpf(3)/2)
psi3_32 = mpmath.psi(3, mpmath.mpf(3)/2)

try_pslq(
    "with_polygamma",
    [mpmath.mpf(1), pi2, pi4, pi6, z3, z5,
     G, beta4, psi1_32, psi3_32,
     psi1_32**2, psi1_32*psi3_32, psi3_32**2,
     G**2, G*beta4, beta4**2],
    ["1", "pi2", "pi4", "pi6", "z3", "z5",
     "G", "beta4", "psi1(3/2)", "psi3(3/2)",
     "psi1(3/2)^2", "psi1(3/2)*psi3(3/2)", "psi3(3/2)^2",
     "G^2", "G*beta4", "beta4^2"]
)

# --- Attempt 8: Ultra-wide PSLQ (47 elements) ---
print("\n  Attempting ultra-wide PSLQ...")

ultra_basis = [
    mpmath.mpf(1),
    pi2, pi4, pi6, pi8,
    z3, z5, z7,
    G, beta4, beta3,
    pi2*z3, pi2*z5, pi2*z7,
    pi4*z3, pi4*z5,
    z3**2, z3*z5, z5**2,
    z3*G, z3*beta4, z5*G, z5*beta4,
    pi2*G, pi2*beta4,
    pi4*G, pi4*beta4,
    G**2, G*beta4, beta4**2,
    ln2, ln2**2,
    ln2*z3, ln2*G, ln2*beta4, ln2*pi2,
    Li3_half, Li4_half,
    hz_vals[2], hz_vals[3], hz_vals[4], hz_vals[5],
    hz_vals[2]**2, hz_vals[2]*hz_vals[4], hz_vals[4]**2,
    gamma_em,
    Cl2_pi3,
]

ultra_labels = [
    "1",
    "pi2", "pi4", "pi6", "pi8",
    "z3", "z5", "z7",
    "G", "beta4", "beta3",
    "pi2*z3", "pi2*z5", "pi2*z7",
    "pi4*z3", "pi4*z5",
    "z3^2", "z3*z5", "z5^2",
    "z3*G", "z3*beta4", "z5*G", "z5*beta4",
    "pi2*G", "pi2*beta4",
    "pi4*G", "pi4*beta4",
    "G^2", "G*beta4", "beta4^2",
    "ln2", "ln2^2",
    "ln2*z3", "ln2*G", "ln2*beta4", "ln2*pi2",
    "Li3(1/2)", "Li4(1/2)",
    "hz(2,3/2)", "hz(3,3/2)", "hz(4,3/2)", "hz(5,3/2)",
    "hz(2,3/2)^2", "hz(2,3/2)*hz(4,3/2)", "hz(4,3/2)^2",
    "gamma",
    "Cl2(pi/3)",
]

try_pslq(
    "ultra_wide_47",
    ultra_basis,
    ultra_labels,
    maxcoeff=1000000,
    tol=1e-40,
)

# --- Attempt 9: Ratio S_min/D(4)^2 ---
print("\n  Trying PSLQ on S_min / D(4)^2...")
ratio = S_min / D4**2
try_pslq(
    "ratio_Smin_D4sq",
    [mpmath.mpf(1), pi2, pi4, z3, z5, G, beta4,
     pi2*G, pi2*beta4, G**2, G*beta4, beta4**2,
     z3**2, z3*G, z3*beta4, ln2, Li3_half],
    ["1", "pi2", "pi4", "z3", "z5", "G", "beta4",
     "pi2*G", "pi2*beta4", "G^2", "G*beta4", "beta4^2",
     "z3^2", "z3*G", "z3*beta4", "ln2", "Li3(1/2)"],
)

# We need to replace S_min in the PSLQ call for the ratio
# The try_pslq function uses S_min global. Override:
print("\n  Trying PSLQ on the ratio directly...")
saved_smin = S_min
try:
    S_min = ratio
    try_pslq(
        "ratio_direct",
        [mpmath.mpf(1), pi2, pi4, z3, z5, G, beta4,
         pi2*G, pi2*beta4, G**2, G*beta4, beta4**2,
         z3**2, z3*G, z3*beta4, ln2, Li3_half],
        ["1", "pi2", "pi4", "z3", "z5", "G", "beta4",
         "pi2*G", "pi2*beta4", "G^2", "G*beta4", "beta4^2",
         "z3^2", "z3*G", "z3*beta4", "ln2", "Li3(1/2)"],
    )
finally:
    S_min = saved_smin

# --- Attempt 10: Try sum_22 / zeta(2,3/2)^2 ---
ratio_22 = sum_22 / hz_vals[2]**2
print("\n  sum_22 / zeta(2,3/2)^2 = %s" % mpmath.nstr(ratio_22, 40))
saved_smin = S_min
try:
    S_min = ratio_22
    try_pslq(
        "ratio_sum22_hz2sq",
        [mpmath.mpf(1), pi2, pi4, z3, z5, G, beta4,
         ln2, gamma_em, hz_vals[3]/hz_vals[2],
         hz_vals[4]/hz_vals[2]**2],
        ["1", "pi2", "pi4", "z3", "z5", "G", "beta4",
         "ln2", "gamma", "hz3/hz2", "hz4/hz2^2"],
    )
finally:
    S_min = saved_smin

# --- Attempt 11: Multiple Hurwitz zeta approach ---
# zeta_2(s1,s2;a) = sum_{n1>n2>=0} 1/((n1+a)^s1 * (n2+a)^s2)
# = (1/2)*[zeta(s1,a)*zeta(s2,a) - zeta(s1+s2,a)]
# Our sum_22 relates to a shifted version of this.
# Compute the double Hurwitz explicitly
dh_22_32 = (hz_vals[2]**2 - hz_vals[4]) / 2  # zeta_2(2,2; 3/2)
dh_24_32 = (hz_vals[2]*hz_vals[4] - hz_vals[6]) / 2  # zeta_2(2,4; 3/2)
dh_42_32 = dh_24_32  # symmetric
dh_44_32 = (hz_vals[4]**2 - mpmath.hurwitz(8, mpmath.mpf(3)/2)) / 2

print("\n  Double Hurwitz zeta values at a=3/2:")
print("  zeta_2(2,2;3/2) = %s" % mpmath.nstr(dh_22_32, 30))
print("  zeta_2(2,4;3/2) = %s" % mpmath.nstr(dh_24_32, 30))

# Now try: S_min in terms of double Hurwitz
saved_smin = S_min
S_min_val = saved_smin
try:
    S_min = S_min_val
    try_pslq(
        "double_hurwitz",
        [mpmath.mpf(1), pi2, pi4, pi6, pi8,
         hz_vals[2], hz_vals[4],
         dh_22_32, dh_24_32, dh_44_32,
         hz_vals[2]**2, hz_vals[4]**2, hz_vals[2]*hz_vals[4],
         G, beta4, z3, z5, G**2, beta4**2, G*beta4],
        ["1", "pi2", "pi4", "pi6", "pi8",
         "hz(2)", "hz(4)",
         "dh(2,2)", "dh(2,4)", "dh(4,4)",
         "hz(2)^2", "hz(4)^2", "hz(2)*hz(4)",
         "G", "beta4", "z3", "z5", "G^2", "beta4^2", "G*beta4"],
    )
finally:
    S_min = saved_smin

# --- Attempt 12: Hurwitz at quarter-integer shifts (even/odd decomposition) ---
# Since D_even/D_odd use zeta(s, 3/4) and zeta(s, 5/4),
# the min-weighted sum might involve products of these.
D_even4 = mpmath.power(2, -4) * (8 * mpmath.hurwitz(2, mpmath.mpf(3)/4)
                                  - mpmath.mpf(1)/2 * mpmath.hurwitz(4, mpmath.mpf(3)/4))
D_odd4 = mpmath.power(2, -4) * (8 * mpmath.hurwitz(2, mpmath.mpf(5)/4)
                                 - mpmath.mpf(1)/2 * mpmath.hurwitz(4, mpmath.mpf(5)/4))

print("\n  D_even(4) = %s" % mpmath.nstr(D_even4, 30))
print("  D_odd(4)  = %s" % mpmath.nstr(D_odd4, 30))
print("  D(4)      = %s" % mpmath.nstr(D_even4 + D_odd4, 30))
print("  D(4) ref  = %s" % mpmath.nstr(D4, 30))

saved_smin = S_min
try:
    S_min = saved_smin
    try_pslq(
        "even_odd_products",
        [mpmath.mpf(1), pi2, pi4, pi6, pi8,
         D_even4, D_odd4,
         D_even4**2, D_odd4**2, D_even4*D_odd4,
         G, beta4, z3, z5,
         G**2, beta4**2, G*beta4,
         D_even4*G, D_even4*beta4, D_odd4*G, D_odd4*beta4],
        ["1", "pi2", "pi4", "pi6", "pi8",
         "De4", "Do4",
         "De4^2", "Do4^2", "De4*Do4",
         "G", "beta4", "z3", "z5",
         "G^2", "beta4^2", "G*beta4",
         "De4*G", "De4*beta4", "Do4*G", "Do4*beta4"],
    )
finally:
    S_min = saved_smin

# ---------------------------------------------------------------------------
# Step 4: Numerical exploration
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Step 4: Numerical exploration")
print("=" * 70)

print("\n  S_min (100 digits) =")
print("  %s" % mpmath.nstr(saved_smin, 100))

print("\n  Ratios with known constants:")
for label, val in [
    ("pi^2", pi2), ("pi^4", pi4), ("zeta(3)", z3), ("zeta(5)", z5),
    ("G (Catalan)", G), ("beta(4)", beta4),
    ("G*zeta(3)", G*z3), ("G^2", G**2), ("beta(4)^2", beta4**2),
    ("D(4)^2", D4**2), ("Li3(1/2)", Li3_half),
    ("D_even(4)^2", D_even4**2), ("D_odd(4)^2", D_odd4**2),
    ("D_even(4)*D_odd(4)", D_even4*D_odd4),
]:
    r = saved_smin / val
    print("  S_min / %-20s = %s" % (label, mpmath.nstr(r, 25)))

# Check if S_min - some_combination is simpler
print("\n  Testing S_min - known combinations:")
for label, val in [
    ("D(4)^2", D4**2),
    ("2*D_even(4)*D_odd(4)", 2*D_even4*D_odd4),
    ("D_even(4)^2 + D_odd(4)^2", D_even4**2 + D_odd4**2),
    ("(De+Do)^2/2 + (De-Do)^2/2", (D_even4+D_odd4)**2/2 + (D_even4-D_odd4)**2/2),
]:
    diff = saved_smin - val
    print("  S_min - %-35s = %s" % (label, mpmath.nstr(diff, 30)))

# ---------------------------------------------------------------------------
# Step 5: Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Step 5: Summary")
print("=" * 70)

n_attempts = len(results)
n_identified = sum(1 for r in results.values() if r.get("identified", False))
print("\n  Total attempts: %d" % n_attempts)
print("  Identified: %d" % n_identified)
print("  Failed: %d" % (n_attempts - n_identified))

for name, res in results.items():
    status = "IDENTIFIED" if res.get("identified") else "FAILED"
    print("\n  [%s] %s" % (status, name))
    if res.get("identified"):
        for comp, coeff in res["components"].items():
            print("    %s * %s" % (coeff, comp))
        print("    Residual: %.3e" % res["residual"])
    elif "error" in res:
        print("    Error: %s" % res["error"])

# Save results
output = {
    "S_min_str": mpmath.nstr(saved_smin, 100),
    "S_min_float": float(saved_smin),
    "D4_squared_float": float(D4**2),
    "ratio_Smin_D4sq": float(saved_smin / D4**2),
    "precision_dps": 150,
    "n_terms": 10000,
    "component_sums": {
        "sum_22": float(sum_22),
        "sum_24": float(sum_24),
        "sum_44": float(sum_44),
    },
    "attempts": {},
}
for k, v in results.items():
    output["attempts"][k] = {kk: vv for kk, vv in v.items() if kk != "relation"}

outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "data", "smin_identification.json")
os.makedirs(os.path.dirname(outpath), exist_ok=True)
with open(outpath, "w") as f:
    json.dump(output, f, indent=2, default=str)
print("\n  Results saved to %s" % outpath)

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
