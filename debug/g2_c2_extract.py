"""
Extract c2 coefficient in the g-2 curvature expansion on S^3.

Strategy:
  1. Use the EXISTING converged data at n_ext=1 (from g2_extended_nint_v2.json).
  2. Compute F2/Schwinger at n_ext=2..7 using compute_anomalous_magnetic_moment
     with moderate n_int_max (convergence is much faster at higher n_ext since
     the series decays as ~ n_int^{-7}).
  3. Fit the curvature expansion delta = c1*x + c2*x^2 + c3*x^3 + ...
     with c1 = 1/2 fixed (Parker-Toms exact result).
  4. Identify c2 via PSLQ against rational and heat-kernel bases.

The expansion is:
    F2(n_ext) / Schwinger = 1 + c1/lambda^2 + c2/lambda^4 + ...
    where lambda = n_ext + 3/2 (Dirac eigenvalue on unit S^3)
    and x = 1/lambda^2.

Known: c1 = R/12 = 1/2 exactly (Parker-Toms 1985).
Target: c2 ~ 0.174 (approximate from single-point extraction).

Output: debug/data/g2_c2_analytical.json
"""

from __future__ import annotations
import json
import sys
import os
import time
from fractions import Fraction

sys.path.insert(0, '.')

import numpy as np
import mpmath
mpmath.mp.dps = 50

print("=" * 70)
print("  EXTRACT c2 FROM g-2 CURVATURE EXPANSION ON S^3")
print("=" * 70)

# ---- Step 1: Load the existing converged data at n_ext=1 ----
print("\nStep 1: Load existing converged data")

with open('debug/data/g2_extended_nint_v2.json') as f:
    v2_data = json.load(f)

delta_nex1 = float(v2_data['corrected_delta'])
f2_over_s_nex1 = float(v2_data['corrected_F2_over_S'])
print(f"  n_ext=1: F2/S = {f2_over_s_nex1:.15f} (n_int=0..25, tail-corrected)")
print(f"  n_ext=1: delta = {delta_nex1:.15f}")

# ---- Step 2: Compute at higher n_ext ----
print("\nStep 2: Compute vertex correction at n_ext = 2..7")
print("  (Higher n_ext converge faster - fewer n_int needed)")

from geovac.qed_anomalous_moment import compute_anomalous_magnetic_moment

# n_ext=2..7 with decreasing n_int_max (convergence ~ n_int^{-7})
# At n_ext=1, 99.3% comes from n_int <= 5.
# At n_ext >= 2, convergence is even faster.
targets = [
    (2, 8),
    (3, 6),
    (4, 5),
    (5, 4),
    (6, 4),
    (7, 3),
]

results = {
    1: {'n_ext': 1, 'lambda': 2.5, 'F2_over_S': f2_over_s_nex1,
        'delta': delta_nex1, 'n_int_max': 25, 'source': 'cached'},
}

for n_ext, n_int_max in targets:
    print(f"  Computing n_ext={n_ext}, n_int_max={n_int_max}...", end='', flush=True)
    t0 = time.time()
    try:
        r = compute_anomalous_magnetic_moment(n_ext, n_int_max)
        elapsed = time.time() - t0
        ratio = r['F2_over_schwinger']
        delta = ratio - 1.0
        lam = (2*n_ext + 3) / 2.0
        results[n_ext] = {
            'n_ext': n_ext,
            'lambda': lam,
            'F2_over_S': ratio,
            'delta': delta,
            'n_int_max': n_int_max,
            'elapsed_s': elapsed,
            'source': 'computed',
        }
        print(f"  {elapsed:.1f}s  F2/S = {ratio:.12f}  delta = {delta:.12f}")
    except Exception as e:
        elapsed = time.time() - t0
        print(f"  FAILED ({elapsed:.1f}s): {e}")

# ---- Step 3: Fit curvature expansion ----
print("\nStep 3: Fit curvature expansion (c1 = 1/2 fixed)")

# Collect good data
data_pts = []
for n_ext in sorted(results.keys()):
    r = results[n_ext]
    lam = r['lambda']
    x = 1.0 / lam**2  # x = 1/lambda^2
    delta = r['delta']
    data_pts.append((n_ext, lam, x, delta))
    print(f"  n_ext={n_ext}: lambda={lam:.1f}  x=1/lam^2={x:.8f}  delta={delta:.12f}")

n_pts = len(data_pts)
x_arr = np.array([p[2] for p in data_pts])
delta_arr = np.array([p[3] for p in data_pts])

# Subtract c1 = 1/2 contribution
c1_exact = 0.5
delta_adj = delta_arr - c1_exact * x_arr  # residual = c2*x^2 + c3*x^3 + ...

print(f"\n  After subtracting c1*x (c1 = 1/2):")
for i, (n_ext, lam, x, delta) in enumerate(data_pts):
    print(f"  n_ext={n_ext}: residual = {delta_adj[i]:.12e}  (= c2*x^2 + c3*x^3 + ...)")

# Polynomial fits at increasing degree
print("\n  Constrained fits (c1 = 1/2 fixed):")
best_c2_estimates = []

for degree in range(1, min(n_pts, 5)):
    # Fit residual = c2*x^2 + c3*x^3 + ... + c_{d+1}*x^{d+1}
    V = np.zeros((n_pts, degree))
    for i in range(n_pts):
        for j in range(degree):
            V[i, j] = x_arr[i]**(j+2)

    coeffs, residuals, rank, sv = np.linalg.lstsq(V, delta_adj, rcond=None)
    pred = V @ coeffs
    max_err = np.max(np.abs(delta_adj - pred))

    labels = ['c2', 'c3', 'c4', 'c5']
    print(f"\n  Degree-{degree} fit (terms x^2 ... x^{degree+1}):")
    for j in range(degree):
        print(f"    {labels[j]} = {coeffs[j]:.12f}")
    print(f"    Max residual: {max_err:.2e}")

    best_c2_estimates.append((degree, coeffs[0], max_err))

# Use Richardson extrapolation: compare c2 estimates at different degrees
print("\n  Richardson extrapolation of c2:")
for deg, c2_val, err in best_c2_estimates:
    print(f"    degree-{deg}: c2 = {c2_val:.12f}  (max residual {err:.2e})")

# Best estimate: highest degree with good conditioning
if len(best_c2_estimates) >= 2:
    c2_best = best_c2_estimates[-1][1]  # highest degree
    c2_low = best_c2_estimates[-2][1]   # second highest
    # Richardson: c2_true ~ c2_high + (c2_high - c2_low) * factor
    # This is approximate; the key is the convergence pattern.
    print(f"\n  Convergence: c2 estimates = {[f'{e[1]:.10f}' for e in best_c2_estimates]}")
else:
    c2_best = best_c2_estimates[0][1]

print(f"\n  BEST c2 estimate: {c2_best:.12f}")

# ---- Step 4: Analytical identification ----
print("\nStep 4: Analytical identification of c2")

c2_mp = mpmath.mpf(str(c2_best))

# PSLQ against rational basis
print(f"  c2 = {mpmath.nstr(c2_mp, 15)}")

# Simple rational approximants
print("\n  Best rational approximants (err < 0.5%):")
best_rats = []
for d in range(1, 1000):
    n = round(float(c2_mp) * d)
    if n > 0:
        from math import gcd
        g = gcd(n, d)
        nn, dd = n // g, d // g
        rat_val = nn / dd
        err = abs(rat_val - float(c2_mp))
        rel = err / float(c2_mp)
        if rel < 0.005:
            best_rats.append((nn, dd, rat_val, err, rel))

best_rats.sort(key=lambda x: x[3])
for nn, dd, val, err, rel in best_rats[:15]:
    print(f"    {nn}/{dd} = {val:.12f}  err = {err:.2e} ({rel*100:.4f}%)")

# PSLQ against various bases
print("\n  PSLQ identifications:")

# Pure rational
try:
    rel = mpmath.pslq([c2_mp, mpmath.mpf(1)], tol=1e-8, maxcoeff=10000)
    if rel and rel[0] != 0:
        p, q = -rel[1], rel[0]
        from math import gcd
        g = gcd(abs(p), abs(q))
        print(f"    Rational: c2 = {p//g}/{q//g} = {float(p)/float(q):.12f}")
except:
    print(f"    Rational PSLQ: no identification")

# With pi^2
try:
    rel = mpmath.pslq([c2_mp, mpmath.mpf(1), mpmath.pi**2], tol=1e-8, maxcoeff=500)
    if rel and rel[0] != 0:
        print(f"    PSLQ(1, pi^2): {rel[0]}*c2 + {rel[1]} + {rel[2]}*pi^2 = 0")
        val = (-rel[1] - rel[2]*mpmath.pi**2) / rel[0]
        print(f"    => c2 = {mpmath.nstr(val, 15)}")
except:
    pass

# With pi^2 and pi^4
try:
    rel = mpmath.pslq([c2_mp, mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4],
                       tol=1e-8, maxcoeff=200)
    if rel and rel[0] != 0:
        print(f"    PSLQ(1, pi^2, pi^4): {rel}")
except:
    pass

# Candidates from heat kernel analysis
print("\n  Heat kernel candidates for c2:")

# On unit S^3:
# a1_per_comp = R/6 = 1
# a2 (propagator) per comp from Vassilevich formula:
#   a2_density = (1/360)(5R^2 - 2Ric^2 + 2Riem^2) = 180/360 = 1/2 (geometric)
#   + E^2 - R*E/6 = 9/4 - 3/2 = 3/4 (Lichnerowicz E=R/4=3/2)
#   + (1/12)*Omega^2_per_comp (spin connection)
#
# IMPORTANT: The code in qed_vacuum_polarization uses a DIFFERENT convention
# (the propagator SD formula with -30RE + 60E^2, giving integrand=45).
# The Vassilevich Eq 4.3 formula gives:
#   a2_geom = 1/2, a2_E = 3/4, a2_Omega = to be determined
#
# Per component: a2_geom/dim_S = 1/8 (the code gives 1/8 for a2_density_per_comp)
#
# The full a2 per component (with spin connection) = 9/8 (computed in Section 2)

# With proper vertex weights:
# c1 = (1/2) * a1_per_comp = (1/2) * 1 = 1/2  [matches exact]
#
# For c2 we need the vertex proper-time structure.
# The Schwinger kernel for F2 is proportional to z(1-z).
# At order s^2, three types of curvature insertion:
#   (A) a2 from one propagator:  weight ~ integral of z^3(1-z) or z(1-z)^3
#   (B) a1*a1 from both:        weight ~ integral of z^2(1-z)^2
#
# The Parker-Toms factor 1/2 for c1 comes from:
#   c1 = a1 * [int z^2(1-z)dz] / [int z(1-z)dz] = a1 * (1/12)/(1/6) = a1/2
# (Only one propagator contributes to F2 magnetic part.)

# For c2 (magnetic part from one propagator):
# Type A: a2 * int z^3(1-z)dz / int z(1-z)dz = a2 * (1/20)/(1/6) = 3*a2/10
# Type B: (a1)^2 * int z^2(1-z)^2 dz / int z(1-z)dz
#        = (a1)^2 * (1/30)/(1/6) = (a1)^2/5
# But only ONE propagator carries the magnetic form factor:
# Type A: (3/10)*a2  (a2 from magnetic propagator only)
# Type B: need vertex projection -- half goes to F2, half to F1?

# Various candidate formulas for c2:
R = 6
a1 = 1  # a1_per_comp = R/6
a2_prop = 1/8  # a2 per comp from propagator SD (code: integrand=45, density=1/8)
a2_full = 9/8  # a2 per comp including spin connection

candidates = [
    # Pure vertex structure with propagator a2
    ("(3/10)*a2_prop + (1/5)*a1^2", 0.3*a2_prop + 0.2*a1**2),
    ("(3/10)*a2_prop", 0.3*a2_prop),
    ("(1/2)*a2_prop", 0.5*a2_prop),
    ("(1/2)*a2_prop + (1/10)*a1^2", 0.5*a2_prop + 0.1*a1**2),
    ("c1^2 = 1/4", 0.25),
    ("c1^2/2 = 1/8", 0.125),
    ("7/40 = 7*Delta", 7/40),
    ("9/40 (Parker-Toms vertex)", 9/40),
    ("25/144", 25/144),
    ("3/16", 3/16),
    ("5/28", 5/28),
    ("5/29", 5/29),
    # With full a2
    ("(3/10)*a2_full + (1/5)*a1^2", 0.3*a2_full + 0.2*a1**2),
    ("(3/10)*a2_full", 0.3*a2_full),
    ("(1/2)*a2_full", 0.5*a2_full),
    # Simple rationals near 0.174
    ("7/40", 7/40),
    ("1/6", 1/6),
    ("7/36", 7/36),
    ("13/75", 13/75),
    ("31/180", 31/180),
    ("29/168", 29/168),
    # Combinations with c1
    ("c1*(c1+1/5) = 0.5*0.7", 0.5*0.7),
    ("c1*(c1-1/6) = 1/6", 1/6),
    ("c1^2*(1-1/5) = 1/5", 0.25*0.8),
    ("c1^2*(2/3+1/30) = 7/30*c1^2 ... wrong", 7/30*0.25),
    # 25/144 is close, and = (5/12)^2
    ("(5/12)^2 = 25/144", 25/144),
    # Try rationals involving pi
    ("1/(2*pi) + 1/20", 1/(2*np.pi) + 0.05),
]

candidates.sort(key=lambda x: abs(x[1] - float(c2_mp)))
print(f"  c2 = {float(c2_mp):.12f}")
for label, val in candidates[:15]:
    err = abs(val - float(c2_mp))
    pct = err / float(c2_mp) * 100
    mark = " <-- MATCH" if pct < 0.05 else (" <-- CLOSE" if pct < 0.5 else "")
    print(f"    {label:50s} = {val:.12f}  ({pct:.3f}%){mark}")

# ---- Step 5: Save results ----
print("\nStep 5: Save results")

output = {
    'description': 'c2 extraction from g-2 curvature expansion on S^3',
    'method': 'multi-n_ext polynomial fit with c1=1/2 fixed',
    'c1_exact': '1/2',
    'c1_source': 'Parker-Toms: R/12 on unit S^3',
    'data_points': {},
    'fit_results': [],
    'c2_best': float(c2_best),
}

for n_ext in sorted(results.keys()):
    r = results[n_ext]
    output['data_points'][str(n_ext)] = {
        'n_ext': r['n_ext'],
        'lambda': r['lambda'],
        'F2_over_S': r['F2_over_S'],
        'delta': r['delta'],
        'n_int_max': r['n_int_max'],
        'source': r['source'],
    }

for deg, c2_val, err in best_c2_estimates:
    output['fit_results'].append({
        'degree': deg,
        'c2': float(c2_val),
        'max_residual': float(err),
    })

# Add candidate analysis
output['candidates'] = {}
for label, val in candidates[:10]:
    output['candidates'][label] = {
        'value': float(val),
        'error': abs(float(val) - float(c2_best)),
        'error_pct': abs(float(val) - float(c2_best)) / abs(float(c2_best)) * 100,
    }

with open('debug/data/g2_c2_analytical.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"  Saved to debug/data/g2_c2_analytical.json")

print("\n" + "=" * 70)
print("  DONE")
print("=" * 70)
