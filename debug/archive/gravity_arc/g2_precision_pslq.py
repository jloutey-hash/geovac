"""
g-2 curvature analysis with extended data (n_int=0..15).
Extract curvature expansion coefficients with 6-7 digit precision.
Target: F2/Schwinger = 1 + c1/lam^2 + c2/lam^4 + ...
where lam = (2*n_ext+3)/2 and we study n_ext=1 (lam=5/2, lam^2=25/4).
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

# Load all data
with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
    ratio_data = json.load(f)

with open('debug/data/g2_extended_nint.json') as f:
    ext_data = json.load(f)

schwinger = ratio_data['schwinger']
V_mag = ratio_data['V_magnetic_float']

print("=" * 70)
print("  g-2 CURVATURE ANALYSIS — EXTENDED PRECISION")
print("=" * 70)

# === Part 1: Improved delta at n_ext=1 ===
print("\n--- 1. Improved delta at n_ext=1 ---")
delta_raw = ext_data['corrected_F2_over_S'] - 1.0
delta_tail_unc = ext_data['tail_estimate']
print(f"  F2/Schwinger = {ext_data['corrected_F2_over_S']:.12f}")
print(f"  delta = F2/S - 1 = {delta_raw:.10f}")
print(f"  Tail uncertainty: {delta_tail_unc:.2e}")
print(f"  Reliable digits: ~{-np.log10(delta_tail_unc / delta_raw):.1f}")

# At n_ext=1: lam = 5/2, lam^2 = 25/4 = 6.25
lam = 2.5
lam2 = 6.25
R_scalar = 6.0  # Ricci scalar of unit S^3

# Parker-Toms leading curvature correction for Dirac field
# c1 = R/(12*m^2) but on S^3 with lam = (2n+3)/2, the leading term is
# delta ~ c1 / lam^2 where c1 is related to R_scalar
c1_PT = R_scalar / 12.0  # = 0.5
delta_PT = c1_PT / lam2   # = 0.08
residual = delta_raw - delta_PT

print(f"\n  Parker-Toms prediction: c1 = R/12 = {c1_PT}")
print(f"  Leading term: c1/lam^2 = {delta_PT:.6f}")
print(f"  Measured delta = {delta_raw:.10f}")
print(f"  Residual = {residual:.10f}")
print(f"  Parker-Toms accounts for {delta_PT/delta_raw*100:.2f}%")

# === Part 2: Extract c2 from residual ===
print("\n--- 2. Second-order curvature coefficient ---")

# residual ~ c2 / lam^4
c2_measured = residual * lam2**2
print(f"  c2 (from residual * lam^4) = {c2_measured:.10f}")

# Try to identify c2 with PSLQ
print(f"\n  c2 = {c2_measured:.10f}")

# Curvature invariants on unit S^3
# R = 6, R_ij R^ij = 12, R_ijkl R^ijkl = 12, C_ijkl C^ijkl = 0 (conformally flat)
# Commonly appearing rational combinations:
print(f"  c2 / (R/12) = {c2_measured / 0.5:.8f}")  # relative to leading term
print(f"  c2 / (R^2/144) = {c2_measured / (36/144):.8f}")  # R^2/144 = 1/4
print(f"  c2 / (R_ij R^ij / 144) = {c2_measured / (12/144):.8f}")

# Check simple rationals
candidates_c2 = [
    ("1/4", 0.25),
    ("1/5", 0.2),
    ("1/6", 1/6),
    ("2/9", 2/9),
    ("5/18", 5/18),
    ("11/40", 11/40),
    ("4/15", 4/15),
    ("17/60", 17/60),
    ("3/10", 0.3),
    ("7/25", 7/25),
    ("41/150", 41/150),
    ("11/40", 11/40),
    ("53/192", 53/192),
    ("43/156", 43/156),
]

print(f"\n  Simple rational candidates for c2:")
for label, val in sorted(candidates_c2, key=lambda x: abs(x[1] - c2_measured)):
    err = abs(val - c2_measured)
    print(f"    {label:>8s} = {val:.10f}  err = {err:.2e}  ({err/c2_measured*100:.4f}%)")

# === Part 3: Direct PSLQ on delta ===
print("\n--- 3. PSLQ on delta ---")

if HAS_MPMATH:
    delta_mp = mpmath.mpf(str(delta_raw))

    # delta should be a rational function of lam^2 = 25/4
    # Try: delta = a + b/lam^2 + c/lam^4 + d/lam^6
    # But we only have one n_ext value, so we identify delta*lam^2 and delta*lam^4

    d_lam2 = delta_mp * mpmath.mpf('25/4')  # = c1 + c2/lam^2 + ...
    d_lam4 = delta_mp * mpmath.mpf('25/4')**2  # = c1*lam^2 + c2 + c3/lam^2 + ...

    print(f"  delta * lam^2 = {d_lam2}")
    print(f"  delta * lam^4 = {d_lam4}")

    # PSLQ: is delta*lam^2 a simple combination?
    # Expected: delta*lam^2 = 1/2 + c2/lam^2 + ... ~ 0.5278
    d_lam2_res = d_lam2 - mpmath.mpf('1/2')
    print(f"  delta*lam^2 - 1/2 = {d_lam2_res}")
    print(f"  This should be c2/lam^2 + higher = c2*(4/25) + ...")
    c2_from_lam2 = d_lam2_res * mpmath.mpf('25/4')
    print(f"  c2 estimate from this = {c2_from_lam2}")

    # PSLQ on c2 against rational basis
    print(f"\n  PSLQ on c2 = {c2_measured}:")

    # Test: c2 is a rational number
    c2_mp = mpmath.mpf(str(c2_measured))

    # Basis: {c2, 1}
    try:
        rel = mpmath.pslq([c2_mp, mpmath.mpf(1)], tol=1e-6, maxcoeff=10000)
        if rel is not None and rel[0] != 0:
            ratio = -mpmath.mpf(rel[1]) / rel[0]
            print(f"    PSLQ {rel}: c2 = {ratio} = {float(ratio):.10f}")
    except:
        pass

    # Basis: {c2, 1, pi^2}
    try:
        rel = mpmath.pslq([c2_mp, mpmath.mpf(1), mpmath.pi**2], tol=1e-6, maxcoeff=1000)
        if rel is not None and rel[0] != 0:
            terms = []
            labels = ['c2', '1', 'pi^2']
            for r, l in zip(rel, labels):
                if r != 0:
                    terms.append(f"{r}*{l}")
            print(f"    PSLQ (rational+pi^2): {' + '.join(terms)} = 0")
    except:
        pass

    # Basis: {c2, 1, 1/pi}
    try:
        rel = mpmath.pslq([c2_mp, mpmath.mpf(1), 1/mpmath.pi], tol=1e-6, maxcoeff=1000)
        if rel is not None and rel[0] != 0:
            terms = []
            labels = ['c2', '1', '1/pi']
            for r, l in zip(rel, labels):
                if r != 0:
                    terms.append(f"{r}*{l}")
            print(f"    PSLQ (rational+1/pi): {' + '.join(terms)} = 0")
    except:
        pass

    # More targeted: curvature invariant combinations
    # On S^3: R=6, R_ab R^ab=12, Weyl=0
    # Standard Dirac a_2 heat kernel on S^3:
    # a_2 = (1/360)*(5R^2 - 2 R_ab R^ab - 7 R_abcd R^abcd - 12 Delta R)
    #      = (1/360)*(5*36 - 2*12 - 7*12 - 0) = (180-24-84)/360 = 72/360 = 1/5
    # So a_2^Dirac = 1/5 on unit S^3
    a2_dirac = 1.0/5

    # DeWitt: a_2^scalar = (1/180)(R^2 - R_ab^2 + R_abcd^2 - 6 Delta R)
    #        = (1/180)(36 - 12 + 12 - 0) = 36/180 = 1/5
    a2_scalar = 1.0/5

    print(f"\n  Heat kernel a_2 (Dirac, unit S^3) = {a2_dirac}")
    print(f"  Heat kernel a_2 (scalar, unit S^3) = {a2_scalar}")
    print(f"  c2 / a_2 = {c2_measured / a2_dirac:.8f}")

    # The Schwinger term is a one-loop effect, so c2 might be
    # related to the next-order heat kernel coefficient
    # Check c2 = R^2/180 = 36/180 = 1/5
    print(f"\n  c2 vs 1/5 = 0.2: diff = {c2_measured - 0.2:.6e}")
    print(f"  c2 vs 1/4 = 0.25: diff = {c2_measured - 0.25:.6e}")
    print(f"  c2 vs 11/40 = 0.275: diff = {c2_measured - 0.275:.6e}")
    print(f"  c2 vs 7/25 = 0.28: diff = {c2_measured - 0.28:.6e}")

    # PSLQ on delta itself with extended basis
    print(f"\n  Extended PSLQ on delta = {delta_raw}:")

    pi_mp = mpmath.pi
    basis_labels = [
        'delta', '1', 'pi', 'pi^2', '1/pi', 'log(2)', 'sqrt(2)',
        'sqrt(3)', 'sqrt(6)', 'G', 'zeta(3)'
    ]
    G = mpmath.catalan
    basis_vals = [
        delta_mp, mpmath.mpf(1), pi_mp, pi_mp**2, 1/pi_mp,
        mpmath.log(2), mpmath.sqrt(2), mpmath.sqrt(3), mpmath.sqrt(6),
        G, mpmath.zeta(3)
    ]

    try:
        rel = mpmath.pslq(basis_vals, tol=1e-5, maxcoeff=10000)
        if rel is not None and rel[0] != 0:
            terms = [f"{r}*{l}" for r, l in zip(rel, basis_labels) if r != 0]
            print(f"    {' + '.join(terms)} = 0")
            # Solve for delta
            delta_formula = sum(-r*v for r, v, l in zip(rel, basis_vals, basis_labels)
                              if l != 'delta' and r != 0) / rel[0]
            print(f"    delta = {float(delta_formula):.12f}")
            print(f"    actual = {float(delta_mp):.12f}")
            print(f"    diff   = {float(abs(delta_formula - delta_mp)):.2e}")
        else:
            print(f"    No relation found")
    except Exception as e:
        print(f"    Error: {e}")

    # === Part 4: Multi-n_ext curvature fit ===
    print("\n--- 4. Multi-n_ext curvature fit (using n_int=1 data) ---")

    with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
        exact_data = json.load(f)

    # At each n_ext, we have F2(n_int=1) / Schwinger
    # This gives the dominant term of the curvature expansion
    # F2/S = c0 + c1/lam^2 + c2/lam^4 + ...
    # where c0 should be 1 in the flat-space limit

    results = exact_data['results']
    n_ext_vals = [r['n_ext'] for r in results]
    lam_vals = [(2*n + 3)/2.0 for n in n_ext_vals]
    lam2_vals = [l**2 for l in lam_vals]
    F2_nint1 = [r['F2_nint1_float'] for r in results]

    # F2/Schwinger at each n_ext (only n_int=1 contribution)
    F2_over_S_nint1 = [f2 / schwinger for f2 in F2_nint1]

    print(f"  n_ext  lam    lam^2   F2(nint=1)/S")
    for n, l, l2, fs in zip(n_ext_vals, lam_vals, lam2_vals, F2_over_S_nint1):
        print(f"  {n:3d}  {l:5.1f}  {l2:7.2f}  {fs:.10f}")

    # Fit: F2/S = sum_k c_k * (1/lam^2)^k
    # Use x = 1/lam^2 as variable
    x_vals = np.array([1.0/l2 for l2 in lam2_vals])
    y_vals = np.array(F2_over_S_nint1)

    # 4-parameter fit: c0 + c1*x + c2*x^2 + c3*x^3
    X = np.column_stack([np.ones_like(x_vals), x_vals, x_vals**2, x_vals**3])
    coeffs, _, _, _ = np.linalg.lstsq(X, y_vals, rcond=None)
    c0_fit, c1_fit, c2_fit, c3_fit = coeffs

    print(f"\n  4-parameter fit (n_int=1 only):")
    print(f"    c0 = {c0_fit:.10f}  (should be ~4.68 = n_int=1 fraction of flat-space)")
    print(f"    c1 = {c1_fit:.10f}")
    print(f"    c2 = {c2_fit:.10f}")
    print(f"    c3 = {c3_fit:.10f}")

    # Richardson extrapolation on lam^2 * F2/S
    print(f"\n  Richardson extrapolation of lam^2 * (F2/S)(nint=1):")
    lam2_F2_S = [l2 * fs for l2, fs in zip(lam2_vals, F2_over_S_nint1)]
    for i, (n, val) in enumerate(zip(n_ext_vals, lam2_F2_S)):
        print(f"    n_ext={n}: lam^2 * F2/S = {val:.10f}")

    # Neville table
    x_rich = np.array([1.0/l2 for l2 in lam2_vals])
    vals = np.array(lam2_F2_S)
    n_pts = len(x_rich)
    table = np.zeros((n_pts, n_pts))
    table[:, 0] = vals
    for j in range(1, n_pts):
        for i in range(n_pts - j):
            table[i, j] = (x_rich[i+j] * table[i, j-1] - x_rich[i] * table[i+1, j-1]) / (x_rich[i+j] - x_rich[i])

    print(f"\n  Richardson diagonal (c0*lam^2 limit):")
    for j in range(n_pts):
        print(f"    order {j}: {table[0, j]:.12f}")

    # === Part 5: Per-level structure for delta ===
    print("\n--- 5. Per-level B structure (extended) ---")

    all_levels = ratio_data['per_level']

    print(f"  n_int  B_level         cum_B           F2/S")
    cum_B = 0
    for lev in all_levels:
        n = lev['n_int']
        b = lev['B_level']
        cum_B += b
        fs = cum_B / V_mag / schwinger + 1.0  # approximate
        print(f"  {n:3d}  {b:15.6e}  {cum_B:15.6e}  {lev['F2_over_schwinger']:.12f}")

    for lev in ext_data['extended_levels']:
        n = lev['n_int']
        b = lev['B_level']
        cum_B += b
        print(f"  {n:3d}  {b:15.6e}  {cum_B:15.6e}  {lev['F2_over_schwinger']:.12f}")

    # === Part 6: delta * lam^4 PSLQ (the key test) ===
    print("\n--- 6. PSLQ on delta * lam^4 ---")

    d_lam4_val = delta_raw * lam2**2  # delta * (25/4)^2 = delta * 625/16
    print(f"  delta * lam^4 = {d_lam4_val:.10f}")
    print(f"  = c1*lam^2 + c2 + c3/lam^2 + ...")
    print(f"  = 0.5*6.25 + c2 + ...")
    print(f"  = 3.125 + c2 + ...")
    c2_direct = d_lam4_val - 3.125
    print(f"  c2 (direct) = {c2_direct:.10f}")

    # Same as before but now with extended precision
    c2_mp = mpmath.mpf(str(c2_direct))

    print(f"\n  PSLQ on c2 against rational basis:")
    for maxc in [100, 1000, 10000]:
        try:
            rel = mpmath.pslq([c2_mp, mpmath.mpf(1)], tol=1e-5, maxcoeff=maxc)
            if rel is not None and rel[0] != 0:
                ratio = -mpmath.mpf(rel[1]) / rel[0]
                err = abs(float(ratio) - float(c2_mp))
                print(f"    maxcoeff={maxc}: c2 = {ratio} = {float(ratio):.10f} (err {err:.2e})")
        except:
            pass

    # PSLQ against curvature invariants
    print(f"\n  PSLQ against curvature basis:")
    # On S^3: relevant invariants per unit volume
    # R = 6, R^2 = 36, R_ab R^ab = 12, R_abcd R^abcd = 12
    # Known one-loop Dirac g-2 corrections involve:
    # R/(12m^2) at O(1/m^2) (Parker-Toms)
    # At O(1/m^4): (R^2 - R_ab R^ab)/360 + ...
    # Let's try PSLQ with {c2, 1, R, R^2, R_ab^2, R_abcd^2}
    # = {c2, 1, 6, 36, 12, 12}

    basis_curv = [c2_mp, mpmath.mpf(1)]
    labels_curv = ['c2', '1']

    try:
        rel = mpmath.pslq(basis_curv, tol=1e-5, maxcoeff=200)
        if rel is not None and rel[0] != 0:
            ratio = -mpmath.mpf(rel[1]) / rel[0]
            print(f"    Rational: c2 = {ratio}")
    except:
        pass

    # === Part 7: Summary table ===
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)

    print(f"""
  At n_ext=1 (ground state, unit S^3):
    lam = 5/2, lam^2 = 25/4

    F2/Schwinger = 1 + delta
    delta = {delta_raw:.10f} +/- {delta_tail_unc:.1e}

  Curvature expansion: delta = c1/lam^2 + c2/lam^4 + ...
    c1 = R/12 = 1/2 (Parker-Toms, confirmed at {delta_PT/delta_raw*100:.1f}%)
    c2 = {c2_direct:.10f} (from residual)

  Best rational approximants for c2:
    c2 ~ 0.278 (between 5/18 = 0.2778 and 11/40 = 0.275)

  Precision: ~{-np.log10(delta_tail_unc / delta_raw):.1f} digits on delta
             ~{-np.log10(delta_tail_unc * lam2):.1f} digits on c2
""")

else:
    print("  mpmath not available — skipping PSLQ")


if __name__ == '__main__':
    pass
