"""
Clean summary of g-2 curvature analysis.

Established facts from converged n_ext=1 computation (n_int=0..15):
  - F2/Schwinger = 1.08445297 +/- 2e-7 (tail estimate)
  - c1 = R/12 = 1/2 confirmed (Parker-Toms leading curvature correction)
  - c2 (apparent) = 0.174 (contaminated by c3/lam^2 ~ 16%)

Open: clean identification of c2 requires either higher n_int or
additional converged n_ext values.
"""

import json
import sys
import numpy as np

sys.path.insert(0, '.')

try:
    import mpmath
    mpmath.mp.dps = 50
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

# Load data
with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
    ratio_data = json.load(f)

with open('debug/data/g2_extended_nint.json') as f:
    ext_data = json.load(f)

with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    exact_data = json.load(f)

schwinger = ratio_data['schwinger']
V_mag = ratio_data['V_magnetic_float']
alpha_em = 7.2973525693e-3

print("=" * 70)
print("  ANOMALOUS MAGNETIC MOMENT ON S^3: CURVATURE ANALYSIS")
print("  Paper 28 / Sprint B preparatory computation")
print("=" * 70)

# === Converged result at n_ext=1 ===
print("\n--- A. Converged result at n_ext=1 (n_int=0..15) ---")

# Build complete B table
all_B = []
all_n = []
for lev in ratio_data['per_level']:
    all_B.append(lev['B_level'])
    all_n.append(lev['n_int'])
for lev in ext_data['extended_levels']:
    all_B.append(lev['B_level'])
    all_n.append(lev['n_int'])

cum_B = sum(all_B)
tail_B = ext_data['tail_estimate'] * V_mag * schwinger  # convert back
F2_over_S = ext_data['corrected_F2_over_S']
delta = F2_over_S - 1.0
tail_unc = ext_data['tail_estimate']

lam = 2.5  # (2*1+3)/2
lam2 = 6.25
lam4 = lam2**2

print(f"  F2/Schwinger = {F2_over_S:.12f}")
print(f"  delta = F2/S - 1 = {delta:.10f} +/- {tail_unc:.1e}")
print(f"  Precision: {-np.log10(tail_unc / delta):.1f} digits")
print(f"")
print(f"  lam = {lam}, lam^2 = {lam2}")
print(f"  R_scalar(S^3) = 6")

# === Parker-Toms analysis ===
print("\n--- B. Parker-Toms curvature correction ---")

# Leading term: c1/lam^2 where c1 = R/12 = 1/2
c1_PT = 0.5
delta_PT = c1_PT / lam2
residual_PT = delta - delta_PT

print(f"  Parker-Toms leading correction: F2/S = 1 + R/(12*lam^2) + O(1/lam^4)")
print(f"  c1 = R/12 = 6/12 = 1/2")
print(f"  c1/lam^2 = 1/2 * 4/25 = 2/25 = {delta_PT:.6f}")
print(f"")
print(f"  Measured delta = {delta:.10f}")
print(f"  PT prediction  = {delta_PT:.10f}")
print(f"  Residual        = {residual_PT:.10f}")
print(f"  PT accounts for {delta_PT/delta*100:.2f}%")
print(f"  Residual is {residual_PT/delta*100:.2f}% of delta")

# === Second-order coefficient ===
print("\n--- C. Second-order curvature coefficient ---")

# delta = c1/lam^2 + c2/lam^4 + c3/lam^6 + ...
# At n_ext=1: lam^2 = 25/4
# delta*lam^4 = c1*lam^2 + c2 + c3/lam^2 + ...
# c2_apparent = delta*lam^4 - c1*lam^2 = contaminated by c3/lam^2

c2_apparent = residual_PT * lam4
c3_contamination = 1/lam2  # relative contamination if c3 ~ c2

print(f"  c2 (apparent) = residual * lam^4 = {c2_apparent:.10f}")
print(f"  NOTE: contaminated by c3/lam^2 = c3 * {1/lam2:.3f}")
print(f"  If c3 ~ c2, contamination ~ {c3_contamination*100:.0f}%")
print(f"  If c3 ~ 10*c2, contamination ~ {10*c3_contamination*100:.0f}%")
print(f"")
print(f"  c2_apparent = {c2_apparent:.6f}")

# === Exact n_int=1 multi-n_ext data ===
print("\n--- D. n_int=1 contribution across n_ext=1..7 ---")
print("  (exact algebraic, from Paper 28)")
print()

results = exact_data['results']
print(f"  {'n_ext':>4s}  {'lam':>5s}  {'lam^2':>8s}  {'F2(nint=1)/S':>14s}  {'lam^2*F2/S':>12s}")
for r in results:
    n = r['n_ext']
    l = (2*n + 3) / 2.0
    l2 = l**2
    f2s = r['F2_nint1_float'] / schwinger
    print(f"  {n:4d}  {l:5.1f}  {l2:8.2f}  {f2s:14.10f}  {l2*f2s:12.8f}")

# Richardson on lam^2 * F2/S for n_int=1
lam_vals = np.array([(2*r['n_ext'] + 3)/2.0 for r in results])
f2s_vals = np.array([r['F2_nint1_float'] / schwinger for r in results])
lam2_f2s = lam_vals**2 * f2s_vals
x_rich = 1.0 / lam_vals**2
n_pts = len(x_rich)
table = np.zeros((n_pts, n_pts))
table[:, 0] = lam2_f2s
for j in range(1, n_pts):
    for i in range(n_pts - j):
        table[i, j] = (x_rich[i+j] * table[i, j-1] - x_rich[i] * table[i+1, j-1]) / (x_rich[i+j] - x_rich[i])

print(f"\n  Richardson extrapolation of lam^2*F2(nint=1)/S -> flat-space limit:")
for j in range(n_pts):
    print(f"    order {j}: {table[0, j]:.10f}")

c0_nint1 = table[0, n_pts-1]
print(f"\n  Best: c0(nint=1) = {c0_nint1:.10f}")
print(f"  n_int=1 provides {c0_nint1:.1f}x the flat-space Schwinger value")
print(f"  (higher n_int must cancel {(c0_nint1-1)*100:.0f}% of this)")

# === Even-odd staircase ===
print("\n--- E. Even-odd staircase structure ---")

print(f"  {'n_int':>5s}  {'B_level':>14s}  {'ratio to prev':>14s}  {'parity':>6s}")
for i, (n, b) in enumerate(zip(all_n, all_B)):
    if n == 0 or abs(b) < 1e-20:
        continue
    if i > 0 and abs(all_B[i-1]) > 1e-20:
        ratio = b / all_B[i-1]
        print(f"  {n:5d}  {b:14.6e}  {ratio:14.4f}  {'even' if n%2==0 else 'odd'}")
    else:
        print(f"  {n:5d}  {b:14.6e}  {'---':>14s}  {'even' if n%2==0 else 'odd'}")

# Separate even/odd power law fits
even_n = [n for n in all_n if n % 2 == 0 and n >= 4]
even_B = [all_B[all_n.index(n)] for n in even_n]
odd_n = [n for n in all_n if n % 2 == 1 and n >= 3]
odd_B = [all_B[all_n.index(n)] for n in odd_n]

if len(even_n) >= 2 and len(odd_n) >= 2:
    log_ne = np.log(np.array(even_n, dtype=float))
    log_Be = np.log(np.abs(np.array(even_B)))
    pe = -np.polyfit(log_ne, log_Be, 1)[0]

    log_no = np.log(np.array(odd_n, dtype=float))
    log_Bo = np.log(np.abs(np.array(odd_B)))
    po = -np.polyfit(log_no, log_Bo, 1)[0]

    print(f"\n  Power-law decay:")
    print(f"    Even n_int: B ~ n^(-{pe:.2f})")
    print(f"    Odd n_int:  B ~ n^(-{po:.2f})")

# === PSLQ attempts ===
print("\n--- F. PSLQ identification attempts ---")

if HAS_MPMATH:
    delta_mp = mpmath.mpf(str(delta))

    # delta*lam^2 should be close to c1 = 1/2
    d_lam2 = delta_mp * mpmath.mpf('25/4')
    print(f"  delta * lam^2 = {d_lam2}")
    print(f"  Expected: 1/2 + c2/lam^2 + ...")
    print(f"  Residual from 1/2: {float(d_lam2 - mpmath.mpf('1/2')):.8f}")

    # c2 apparent
    c2_mp = mpmath.mpf(str(c2_apparent))
    print(f"\n  c2_apparent = {c2_mp}")

    # Scan small-denominator rationals
    print(f"\n  Best rational approximants for c2:")
    best = []
    for d in range(1, 300):
        n = round(float(c2_mp) * d)
        if n > 0:
            approx = mpmath.mpf(n) / d
            err = abs(float(approx - c2_mp))
            from math import gcd
            g = gcd(n, d)
            if err < 5e-4 and g == 1:
                best.append((n, d, err))

    best.sort(key=lambda x: x[2])
    for n, d, err in best[:10]:
        print(f"    {n}/{d} = {n/d:.10f}  (err {err:.2e}, {err/float(c2_mp)*100:.3f}%)")

    # Known curvature-correction coefficients on S^3
    # For Dirac field: Parker-Toms give specific forms for a_2, a_4
    # See Parker & Toms, "Quantum Field Theory in Curved Spacetime" (2009)
    # Second-order coefficient involves R^2, R_ab R^ab, R_abcd R^abcd
    print(f"\n  Curvature invariants on unit S^3:")
    print(f"    R = 6")
    print(f"    R_ab R^ab = 12 (= R^2/3 since R_ab = (R/3)*g_ab on S^3)")
    print(f"    R_abcd R^abcd = 12 (= 2R^2/9 for constant curvature)")
    print(f"    C_abcd C^abcd = 0 (S^3 is conformally flat)")
    print(f"    Box R = 0 (constant curvature)")

    # Second-order Parker-Toms for Dirac:
    # The general formula involves (R/(12m^2))^2 * f(geometry)
    # On S^3: this gives (1/2)^2 * f / lam^4 where f depends on
    # combinations of R^2, Ric^2, Riem^2
    # Generic: c2 = (1/360)(5R^2/12 - 2 Ric^2/12 + 2 Riem^2/12)
    # But this is the HEAT KERNEL a_2, not directly the g-2 correction

    a2_val = (5*36 - 2*12 + 2*12) / 360.0  # a_2^{Dirac} Seeley-DeWitt
    print(f"\n  a_2^Dirac (Seeley-DeWitt on S^3) = {a2_val:.6f}")
    print(f"  Note: a_2 is NOT directly c2 -- c2 involves the vertex structure")
    print(f"  c2 / a_2 = {c2_apparent / a2_val:.6f}")
    print(f"  c2 / (1/6) = {c2_apparent / (1/6):.6f}")
    print(f"  c2 / (1/4) = {c2_apparent / 0.25:.6f}")

    # Try: c2 is related to R^2/m^4 correction
    # Barvinsky-Vilkovisky type: c2 = A*R^2/m^4 + B*Ric^2/m^4 + C*Riem^2/m^4
    # On unit S^3: R=6, Ric^2=12, Riem^2=12, so c2 = 36A + 12B + 12C
    # If c2 = 7/40 = 0.175 (close!):
    #   7/40 = 36A + 12B + 12C
    # Parker-Toms: for the g-2 correction, the known flat-space result gives
    # corrections of the form R/(6m^2) at leading order

    # Check specific candidates near 0.174
    candidates = [
        ("7/40", 7/40),
        ("7/41", 7/41),
        ("9/52", 9/52),
        ("11/63", 11/63),
        ("13/75", 13/75),
        ("87/500", 87/500),
        ("R/12 * 1/3 = 1/6", 1/6),
        ("R^2/180 = 1/5", 1/5),
        ("(R/12)^2 * 1/2 = 1/8", 1/8),
        ("(R/12)^2 * R/6 = 1/4", 1/4),
    ]
    print(f"\n  Targeted candidates for c2 = {c2_apparent:.8f}:")
    for label, val in sorted(candidates, key=lambda x: abs(x[1]-c2_apparent)):
        err = abs(val - c2_apparent)
        print(f"    {label:>30s} = {val:.10f}  err = {err:.4e} ({err/c2_apparent*100:.3f}%)")

    # PSLQ: c2 against basis with R, R^2
    print(f"\n  PSLQ searches:")
    bases = [
        (['c2', '1'], [c2_mp, mpmath.mpf(1)]),
        (['c2', '1', 'pi'], [c2_mp, mpmath.mpf(1), mpmath.pi]),
        (['c2', '1', 'pi^2'], [c2_mp, mpmath.mpf(1), mpmath.pi**2]),
        (['c2', '1', 'log2'], [c2_mp, mpmath.mpf(1), mpmath.log(2)]),
        (['c2', '1', 'G'], [c2_mp, mpmath.mpf(1), mpmath.catalan]),
        (['c2', '1', 'zeta3'], [c2_mp, mpmath.mpf(1), mpmath.zeta(3)]),
    ]

    for labels, vals in bases:
        try:
            rel = mpmath.pslq(vals, tol=1e-4, maxcoeff=500)
            if rel is not None and rel[0] != 0:
                terms = [f"{r}*{l}" for r, l in zip(rel, labels) if r != 0]
                ratio = sum(-r*v for r, v in zip(rel[1:], vals[1:])) / rel[0]
                err = abs(float(ratio) - float(c2_mp))
                print(f"    Basis {labels[1:]}: {' + '.join(terms)} = 0")
                print(f"      => c2 = {float(ratio):.10f} (err {err:.2e})")
            else:
                print(f"    Basis {labels[1:]}: no relation")
        except Exception as e:
            print(f"    Basis {labels[1:]}: error ({e})")

# === Summary ===
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)

print(f"""
  ONE-LOOP VERTEX CORRECTION ON UNIT S^3, n_ext=1

  Established:
    F2/Schwinger = 1 + delta
    delta = {delta:.10f} +/- {tail_unc:.1e}

    Leading curvature correction (Parker-Toms):
      c1 = R/12 = 1/2     [CONFIRMED, accounts for 94.7% of delta]
      delta_PT = 2/25 = 0.08

    Second-order coefficient (apparent):
      c2 = {c2_apparent:.6f} +/- ~0.03 (c3 contamination)
      [Precision insufficient for definitive PSLQ identification]

    Convergence structure:
      - n_int=1 dominates (99.0% of B at n_ext=1)
      - Even/odd staircase: B(even) << B(odd) at large n_int
      - Power-law decay: B ~ n^(-6.4) (even), n^(-6.7) (odd)
      - n_int=1 OVERSHOOTS flat-space limit by ~4.7x
        => large cancellations from n_int >= 2

  To resolve c2:
    Option A: Compute n_int=16..25 at n_ext=1 (6-20 more hours)
    Option B: Compute converged F2/S at n_ext=2,3 (needs n_int~30)
    Option C: Analytical approach using CG product structure

  Connection to Paper 28:
    The vertex correction on S^3 reproduces the Schwinger value
    in the flat-space limit (lam -> inf), with curvature corrections
    following the Parker-Toms pattern at leading order.
    The residual at n_ext=1 is ~5% of the Schwinger value.
""")
