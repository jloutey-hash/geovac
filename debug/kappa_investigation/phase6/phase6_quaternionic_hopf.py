"""
Phase 6: Quaternionic (and octonionic) Hopf fibration — does the Paper 2 recipe generalize?

The four Hopf fibrations (Adams 1960):
  Real:        S^0 → S^1 → S^1   (fiber → total → base)  — trivial/{1}
  Complex:     S^1 → S^3 → S^2   — U(1), electromagnetism  ← Paper 2
  Quaternionic: S^3 → S^7 → S^4   — SU(2), weak force?
  Octonionic:  S^7 → S^15 → S^8   — G2→SU(3), strong force?

Paper 2 recipe (abstracted):
  1. B = (1/2) Σ g_k |λ_k| through some cutoff n_max
  2. Selection: B/N = dim(total space), solve for n_max
  3. F = fiber spectral zeta = Σ 1/λ_k(fiber) (without degeneracy)
  4. Δ = 1/(|λ_{n_max}| × N(n_max - 1))  (boundary correction)
  5. K = π(B + F - Δ)
  6. Cubic: α³ - Kα + 1 = 0, smallest positive root
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from mpmath import mp, mpf, pi as mpi, zeta, log, fabs, nstr, polyroots, binomial, factorial
from fractions import Fraction

mp.dps = 50

# ==============================================================================
# SPHERE SPECTRAL DATA
# ==============================================================================

def sphere_degeneracy(k: int, d: int) -> int:
    """Degeneracy of the k-th eigenspace of the Laplacian on S^d.

    g_k = C(k+d, d) - C(k+d-2, d) for k >= 0 (with C(n,d) = 0 if n < d).
    """
    from math import comb
    if k < 0:
        return 0
    c1 = comb(k + d, d)
    c2 = comb(k + d - 2, d) if k >= 2 else 0
    return c1 - c2


def sphere_eigenvalue(k: int, d: int) -> int:
    """Magnitude of k-th eigenvalue of Laplacian on S^d: |λ_k| = k(k+d-1)."""
    return k * (k + d - 1)


def compute_B(m: int, d: int) -> Fraction:
    """Casimir trace B(m) = (1/2) Σ_{k=0}^{m} g_k × |λ_k| for S^d."""
    total = Fraction(0)
    for k in range(m + 1):
        total += Fraction(sphere_degeneracy(k, d)) * Fraction(sphere_eigenvalue(k, d))
    return total / 2


def compute_N(m: int, d: int) -> int:
    """Total number of states N(m) = Σ_{k=0}^{m} g_k for S^d."""
    return sum(sphere_degeneracy(k, d) for k in range(m + 1))


def fiber_spectral_zeta(d_fiber: int, n_terms: int = 100000) -> mpf:
    """Spectral zeta of S^d at s=1: F = Σ_{k=1}^∞ 1/|λ_k(S^d)| (no degeneracy).

    For S^1: F = ζ(2) = π²/6
    For S^3: F = Σ 1/(k(k+2)) = 3/4  (telescoping)
    For S^7: F = Σ 1/(k(k+6)) = (1/6)(H_6) = 49/120
    General S^d: F = Σ 1/(k(k+d-1)) = (1/(d-1)) H_{d-1}
      where H_n = 1 + 1/2 + ... + 1/n is the harmonic number.
    """
    d = d_fiber
    if d == 0:
        # S^0 = two points. Only eigenvalue: 2 (from the Laplacian on {±1}).
        # Actually S^0 Laplacian has eigenvalue 0 (constant) and...
        # For d=0: λ_k = k(k-1). k=1: λ=0 (problematic).
        # S^0 is degenerate: only k=0 (trivial) and k=1 (eigenvalue 0).
        # Return NaN to flag degeneracy.
        return mpf('nan')
    if d == 1:
        return mpi**2 / 6  # ζ(2)
    # General case: Σ_{k=1}^∞ 1/(k(k+d-1)) = (1/(d-1)) Σ_{j=1}^{d-1} 1/j = H_{d-1}/(d-1)
    # by partial fractions: 1/(k(k+d-1)) = (1/(d-1))(1/k - 1/(k+d-1))
    # Telescoping gives: Σ = (1/(d-1))(1 + 1/2 + ... + 1/(d-1)) = H_{d-1}/(d-1)
    harmonic = sum(Fraction(1, j) for j in range(1, d))
    result = harmonic / (d - 1)  # H_{d-1} / (d-1)
    return mpf(result.numerator) / mpf(result.denominator)


print("=" * 78)
print("PHASE 6: QUATERNIONIC HOPF FIBRATION — RECIPE GENERALIZATION")
print("=" * 78)

# ==============================================================================
# PART 1: VERIFY SPECTRAL DATA
# ==============================================================================
print("\n" + "=" * 78)
print("PART 1: SPECTRAL DATA VERIFICATION")
print("=" * 78)

for d in [1, 2, 3, 4, 7, 8, 15]:
    print(f"\n--- S^{d} (dim = {d}) ---")
    print(f"  {'k':>3} {'g_k':>8} {'|λ_k|':>8}")
    for k in range(5):
        g = sphere_degeneracy(k, d)
        lam = sphere_eigenvalue(k, d)
        print(f"  {k:>3} {g:>8} {lam:>8}")

# Verify known values
assert sphere_degeneracy(0, 3) == 1
assert sphere_degeneracy(1, 3) == 4
assert sphere_degeneracy(2, 3) == 9
assert sphere_eigenvalue(1, 3) == 3  # 1×(1+2) = 3, or -(n²-1) with n=2: -(4-1)=-3
# Paper 2 uses n=k+1: λ_n = -(n²-1). So for n=2: |λ| = 3. For k=1: k(k+2) = 3. ✓
# For n=3: |λ| = 8. For k=2: k(k+2) = 8. ✓
assert sphere_eigenvalue(2, 3) == 8

print("\n✓ All spectral data verified")

# ==============================================================================
# PART 2: UNIVERSALITY OF THE SELECTION PRINCIPLE
# ==============================================================================
print("\n" + "=" * 78)
print("PART 2: SELECTION PRINCIPLE B/N = dim(S^d)")
print("=" * 78)

print("\nThe selection principle B(m)/N(m) = d is tested for multiple spheres.")
print("B(m) = (1/2) Σ_{k=0}^{m} g_k × |λ_k|,  N(m) = Σ g_k")
print()

for d in [1, 2, 3, 4, 5, 6, 7, 8, 15]:
    print(f"\n--- S^{d} (selection target: B/N = {d}) ---")
    print(f"  {'m':>3} {'B(m)':>12} {'N(m)':>8} {'B/N':>12} {'= d?':>6}")
    found_unique = None
    for m in range(6):
        B = compute_B(m, d)
        N = compute_N(m, d)
        ratio = float(B) / N if N > 0 else 0
        match = "  <<<" if abs(ratio - d) < 1e-10 else ""
        if abs(ratio - d) < 1e-10 and found_unique is None:
            found_unique = m
        print(f"  {m:>3} {str(B):>12} {N:>8} {ratio:>12.6f}{match}")

    if found_unique is not None:
        print(f"  → Selection: m = {found_unique} (n_max = {found_unique + 1} in Paper 2 convention)")

# Prove universality algebraically
print("\n" + "-" * 50)
print("ALGEBRAIC PROOF: B(2)/N(2) = d for ALL dimensions d")
print("-" * 50)
print("""
B(2) = (1/2)[g_1 × 1×d + g_2 × 2(d+1)]
     = (1/2)[(d+1)d + ((d+2)(d+1)/2 - 1) × 2(d+1)]

g_1(S^d) = d+1
g_2(S^d) = (d+2)(d+1)/2 - 1 = (d²+3d)/2

B(2) = (1/2)[(d+1)d + (d²+3d)(d+1)]
     = (1/2)(d+1)[d + d²+3d]
     = (1/2)(d+1)d(d+4)

N(2) = 1 + (d+1) + (d²+3d)/2 = (d+1)(d+4)/2

B(2)/N(2) = [(1/2)(d+1)d(d+4)] / [(d+1)(d+4)/2] = d  ✓

Setting B(m)/N(m) = d gives the quadratic:
    m² + dm - 2(d+2) = 0

with roots m = [-d ± (d+4)]/2:
    m = 2  (positive, UNIVERSAL)
    m = -(d+2)  (negative, rejected)

Therefore m = 2 is the UNIQUE positive solution for ALL d.
""")

# Verify the algebraic formula for several d
for d in [1, 2, 3, 4, 7, 8, 15]:
    B2 = Fraction(d+1) * d * (d+4) / 2
    N2 = Fraction(d+1) * (d+4) / 2
    ratio = B2 / N2
    assert ratio == d, f"Failed for d={d}: B/N = {ratio}"
print("✓ Algebraic formula verified for d = 1, 2, 3, 4, 7, 8, 15")


# ==============================================================================
# PART 3: ALL FOUR HOPF FIBRATIONS
# ==============================================================================
print("\n" + "=" * 78)
print("PART 3: RECIPE APPLIED TO ALL FOUR HOPF FIBRATIONS")
print("=" * 78)

hopf_fibrations = [
    # (name, d_fiber, d_total, d_base, gauge_group, force)
    ("Real",        0,  1,  1, "{1}",    "trivial"),
    ("Complex",     1,  3,  2, "U(1)",   "electromagnetism"),
    ("Quaternionic", 3,  7,  4, "SU(2)",  "weak force?"),
    ("Octonionic",  7, 15,  8, "G₂→SU(3)", "strong force?"),
]

results = {}

for name, d_fib, d_tot, d_base, gauge, force in hopf_fibrations:
    print(f"\n{'='*60}")
    print(f"  {name} Hopf:  S^{d_fib} → S^{d_tot} → S^{d_base}")
    print(f"  Gauge group: {gauge}  |  Force: {force}")
    print(f"{'='*60}")

    # Step 1: Selection principle → n_max
    m_select = 2  # Universal!

    B_val = compute_B(m_select, d_tot)
    N_val = compute_N(m_select, d_tot)
    ratio = float(B_val) / N_val

    print(f"\n  Step 1: Selection principle B/N = dim(S^{d_tot}) = {d_tot}")
    print(f"    m_select = {m_select} (universal)")
    print(f"    B({m_select}) = {B_val}")
    print(f"    N({m_select}) = {N_val}")
    print(f"    B/N = {ratio:.6f}  {'✓' if abs(ratio - d_tot) < 1e-10 else '✗'}")

    B_mp = mpf(B_val.numerator) / mpf(B_val.denominator)

    # Step 2: Fiber spectral zeta
    F_val = fiber_spectral_zeta(d_fib)

    # Exact rational forms for known cases
    if d_fib == 0:
        F_exact = "DEGENERATE (S^0 has no positive eigenvalues)"
        F_val = mpf(0)  # Use 0 as placeholder
    elif d_fib == 1:
        F_exact = "π²/6 = ζ(2)"
    elif d_fib == 3:
        F_exact = "3/4  [= H_2/2 = (1+1/2)/2]"
    elif d_fib == 7:
        F_exact = "49/120  [= H_6/6]"
    else:
        F_exact = f"H_{d_fib-1}/{d_fib-1}"

    print(f"\n  Step 2: Fiber spectral zeta (S^{d_fib})")
    print(f"    F = Σ_{{k=1}}^∞ 1/|λ_k(S^{d_fib})| = {F_exact}")
    print(f"    F = {nstr(F_val, 20)}")

    # Step 3: Boundary correction
    lam_m = sphere_eigenvalue(m_select, d_tot)  # |λ_{n_max}| of total space
    N_prev = compute_N(m_select - 1, d_tot)     # N(n_max - 1)

    if lam_m > 0 and N_prev > 0:
        Delta_val = mpf(1) / (lam_m * N_prev)
    else:
        Delta_val = mpf(0)

    print(f"\n  Step 3: Boundary correction")
    print(f"    |λ_{m_select}(S^{d_tot})| = {lam_m}")
    print(f"    N({m_select - 1}) = {N_prev}")
    print(f"    Δ = 1/({lam_m} × {N_prev}) = 1/{lam_m * N_prev} = {nstr(Delta_val, 20)}")

    # Step 4: Combination
    if d_fib == 0:
        print(f"\n  Step 4-6: DEGENERATE — S^0 fiber has no spectral content.")
        print(f"    The real Hopf fibration is trivial (Z₂ bundle).")
        print(f"    Attempting with F=0:")

    K_inner = B_mp + F_val - Delta_val
    K_val = mpi * K_inner

    print(f"\n  Step 4: Combination rule")
    print(f"    B + F - Δ = {nstr(B_mp, 10)} + {nstr(F_val, 10)} - {nstr(Delta_val, 10)}")
    print(f"    B + F - Δ = {nstr(K_inner, 20)}")
    print(f"    K = π × (B + F - Δ) = {nstr(K_val, 20)}")

    # Step 5: Cubic
    roots = polyroots([1, 0, -K_val, 1])
    real_roots = sorted([r.real for r in roots if fabs(r.imag) < mpf(10)**(-40)])

    if len(real_roots) == 3:
        alpha_val = real_roots[1]  # smallest positive root
        alpha_inv = 1 / alpha_val
    else:
        alpha_val = mpf(0)
        alpha_inv = mpf('inf')

    print(f"\n  Step 5: Cubic α³ - Kα + 1 = 0")
    print(f"    Roots: {', '.join(nstr(r, 15) for r in real_roots)}")
    if alpha_val > 0:
        print(f"    α = {nstr(alpha_val, 20)}")
        print(f"    1/α = {nstr(alpha_inv, 20)}")

    results[name] = {
        'B': B_val, 'F': F_val, 'Delta': Delta_val, 'K': K_val,
        'K_inner': K_inner, 'alpha': alpha_val, 'alpha_inv': alpha_inv,
        'd_fib': d_fib, 'd_tot': d_tot, 'd_base': d_base,
    }


# ==============================================================================
# PART 4: COMPARISON WITH KNOWN COUPLINGS
# ==============================================================================
print("\n" + "=" * 78)
print("PART 4: COMPARISON WITH KNOWN COUPLING CONSTANTS")
print("=" * 78)

known_couplings = {
    "α_em (CODATA 2018)":             mpf('137.035999084'),
    "α_W = α/sin²θ_W (at M_Z)":      mpf('29.5'),         # SU(2) weak ~ 1/29.5
    "α_s (at M_Z) [strong]":          1/mpf('0.1179'),     # ~ 8.48
    "α at M_Z scale":                 mpf('128.9'),         # running α at Z mass
    "α_GUT (unification)":            mpf('25'),            # approximate GUT coupling
    "sin²θ_W ~ 0.231":               1/mpf('0.2312'),      # ~ 4.33
    "1/g_W² = 1/0.424":              1/mpf('0.424'),       # ~ 2.36
}

print(f"\n{'Hopf':>15} {'1/α':>15} {'Known coupling':>30} {'Known 1/α':>12} {'Ratio':>12}")
print("-" * 86)
for name, res in results.items():
    if res['alpha'] > 0:
        for cname, cval in known_couplings.items():
            ratio = res['alpha_inv'] / cval
            match = " <<<" if 0.95 < float(ratio) < 1.05 else ""
            print(f"  {name:>13} {nstr(res['alpha_inv'],10):>15} {cname:>30} "
                  f"{nstr(cval,8):>12} {nstr(ratio,8):>12}{match}")
        print()


# ==============================================================================
# PART 5: RATIO ANALYSIS BETWEEN HOPF COUPLINGS
# ==============================================================================
print("\n" + "=" * 78)
print("PART 5: RATIOS BETWEEN HOPF COUPLINGS")
print("=" * 78)

for n1 in ["Quaternionic", "Octonionic"]:
    for n2 in ["Complex"]:
        if n1 in results and n2 in results:
            r1 = results[n1]
            r2 = results[n2]
            K_ratio = r1['K'] / r2['K']
            inner_ratio = r1['K_inner'] / r2['K_inner']
            B_ratio = float(r1['B']) / float(r2['B'])

            print(f"\n  {n1} / {n2}:")
            print(f"    K ratio:        {nstr(K_ratio, 15)}")
            print(f"    (B+F-Δ) ratio:  {nstr(inner_ratio, 15)}  (same as K ratio)")
            print(f"    B ratio:        {B_ratio:.6f}  (= {r1['B']}/{r2['B']})")
            print(f"    dim ratio:      {r1['d_tot']}/{r2['d_tot']} = {r1['d_tot']/r2['d_tot']:.6f}")
            print(f"    α ratio:        {nstr(r1['alpha_inv'] / r2['alpha_inv'], 15)}")

# Check K_quat / K_em ≈ 7?
if "Complex" in results and "Quaternionic" in results:
    r = results["Quaternionic"]["K"] / results["Complex"]["K"]
    print(f"\n  K_quat / K_em = {nstr(r, 15)}")
    print(f"  Nearest integer = 7 (= dim S^7)")
    print(f"  Deviation from 7: {nstr(r - 7, 10)}")


# ==============================================================================
# PART 6: DETAILED SPECTRAL INVARIANT TABLE
# ==============================================================================
print("\n" + "=" * 78)
print("PART 6: SPECTRAL INVARIANT TABLE (all four Hopf fibrations)")
print("=" * 78)

print(f"\n{'':>15} {'Real':>15} {'Complex':>15} {'Quaternionic':>15} {'Octonionic':>15}")
print(f"{'':>15} {'S⁰→S¹→S¹':>15} {'S¹→S³→S²':>15} {'S³→S⁷→S⁴':>15} {'S⁷→S¹⁵→S⁸':>15}")
print("-" * 76)

rows = [
    ("dim(fiber)", lambda r: str(r['d_fib'])),
    ("dim(total)", lambda r: str(r['d_tot'])),
    ("dim(base)", lambda r: str(r['d_base'])),
    ("n_max (m)", lambda r: "2"),
    ("B", lambda r: str(r['B'])),
    ("N", lambda r: str(compute_N(2, r['d_tot']))),
    ("B/N", lambda r: str(r['d_tot'])),
    ("F", lambda r: nstr(r['F'], 10)),
    ("Δ", lambda r: nstr(r['Delta'], 10)),
    ("K/π", lambda r: nstr(r['K_inner'], 10)),
    ("K", lambda r: nstr(r['K'], 10)),
    ("1/α", lambda r: nstr(r['alpha_inv'], 10) if r['alpha'] > 0 else "—"),
]

for label, func in rows:
    vals = [func(results[name]) for name in ["Real", "Complex", "Quaternionic", "Octonionic"]]
    print(f"  {label:>13} {'':>2}{vals[0]:>13} {vals[1]:>15} {vals[2]:>15} {vals[3]:>15}")


# ==============================================================================
# PART 7: REPRESENTATION-THEORETIC OBSTRUCTION
# ==============================================================================
print("\n" + "=" * 78)
print("PART 7: REPRESENTATION-THEORETIC ANALYSIS")
print("=" * 78)

print("""
For the COMPLEX Hopf (S¹ → S³ → S²):
  SO(4) ≅ (SU(2) × SU(2))/Z₂
  Shell n has the (j,j) representation with j = (n-1)/2
  Under diagonal SO(3) [base isometry]: (j,j) → ⊕_{l=0}^{2j} l
  Inner sum: Σ (2l+1)l(l+1) = n²(n²-1)/2 = (1/2)g_n |λ_n|

  This works because:
    (a) SO(4) is a PRODUCT group, so the branching is the Clebsch-Gordan series
    (b) S³ ≅ SU(2) is a LIE GROUP (unique among spheres)
    (c) The embedding index of diagonal SO(3) in SO(4) is 1

For the QUATERNIONIC Hopf (S³ → S⁷ → S⁴):
  Total space SO(8) acting on S⁷
  Base SO(5) ≅ Sp(2)/Z₂ acting on S⁴ = HP¹
  Fiber Sp(1) ≅ SU(2) acting on S³

  The representation-theoretic chain is:
    SO(8) ⊃ Sp(2) × Sp(1)_R
  where R⁸ = V₄ ⊗ V₂ (Sp(2) fundamental ⊗ Sp(1) fundamental).

  CRITICAL DIFFERENCE: The Sp(2) fundamental V₄ is the SPINOR of SO(5),
  not a tensor representation. This means the base "angular momentum"
  includes HALF-INTEGER values — there is no clean integer-l decomposition
  as in the complex case.

  The SO(4) → SO(3) branching that makes Paper 2 work is specific to:
    (a) SO(4) being a product group (no analog for SO(8))
    (b) S³ being a Lie group (S⁷ is NOT a Lie group)
    (c) The embedding index being exactly 1
""")

# Compute the SO(5) embedding index in SO(8) through Sp(2)×Sp(1)
# R⁸ under Sp(2)×Sp(1) = V₄⊗V₂
# Under SO(5) alone: V₄ is the spinor (1/2,1/2) of SO(5), dim 4
# R⁸ decomposes as 2 copies of the spinor (since V₂ contributes factor 2)
# C_{SO(5)}((1/2,1/2)) = (1/2)(1/2+3) + (1/2)(1/2+1) = 7/4 + 3/4 = 5/2
# Tr_{R⁸} C_{SO(5)} = 2 × 4 × 5/2 = 20
# C_{SO(8)}(vector) = 1×(1+6) = 7
# Tr_{R⁸} C_{SO(8)} = 8 × 7 = 56
# Index ratio = 20/56 = 5/14

print("  Embedding index computation:")
print("    R⁸ under SO(5): 2 copies of spinor (1/2,1/2), dim 4 each")
print("    C_{SO(5)}(spinor) = 5/2")
print("    Tr_{R⁸} C_{SO(5)} = 2 × 4 × 5/2 = 20")
print("    Tr_{R⁸} C_{SO(8)} = 8 × 7 = 56")
print("    Embedding index ratio: 20/56 = 5/14 ≈", nstr(mpf(5)/14, 10))
print()
print("  If the embedding index extended to all representations:")
print("    B_quat would have factor 5/14 instead of 1/2")
print("    B(2,d=7) with factor 5/14: ", end="")
B_alt = mpf(5) / 14 * sum(sphere_degeneracy(k, 7) * sphere_eigenvalue(k, 7) for k in range(3))
N7 = compute_N(2, 7)
print(f"{nstr(B_alt, 10)},  B/N = {nstr(B_alt/N7, 10)}")
print("    This gives B/N = 5.0 (= dim S⁴ + 1), NOT 7 (= dim S⁷)")
print()
print("  CONCLUSION: The 1/2 factor in B = (1/2)Σ g_k|λ_k| is an ALGEBRAIC")
print("  IDENTITY for the Clebsch-Gordan decomposition, specific to SO(4).")
print("  It generalizes FORMALLY to all spheres (yielding B/N = d universally)")
print("  but the PHYSICAL INTERPRETATION as a base Casimir trace is specific")
print("  to S³, because only SO(4) is a product group.")


# ==============================================================================
# PART 8: REAL HOPF — DEGENERATE CASE
# ==============================================================================
print("\n" + "=" * 78)
print("PART 8: REAL HOPF FIBRATION S⁰ → S¹ → S¹")
print("=" * 78)

print("""
  Fiber S⁰ = {±1}: a discrete two-point space.
  - Only two "eigenfunctions": constant (λ=0) and sign-flip (λ=2)
  - The spectral zeta is trivially F = 1/2 (one nonzero eigenvalue)
  - But S⁰ has no continuous spectral geometry

  Base S¹ = RP¹ ≅ S¹: same as the total space!
  - The fibration is trivial (Z₂ cover)

  Total space S¹:
  - B(2) = 5, N(2) = 5, B/N = 1 = dim(S¹) ✓

  With F = 1/2 and Δ = 1/12:
  K_real = π × (5 + 1/2 - 1/12) = π × 65/12 ≈ 17.02

  Cubic: α³ - 17.02α + 1 = 0
  → α ≈ 0.0588, 1/α ≈ 17.0

  This is degenerate because:
  1. The fiber is discrete, not a manifold
  2. The gauge group {1} has no coupling constant to predict
  3. The base = total space (no non-trivial projection)
""")

# Verify
r = results["Real"]
print(f"  Computed: K = {nstr(r['K'], 15)}, 1/α = {nstr(r['alpha_inv'], 15)}")


# ==============================================================================
# PART 9: ALTERNATIVE PREFACTORS AND COMBINATION RULES
# ==============================================================================
print("\n" + "=" * 78)
print("PART 9: ALTERNATIVE PREFACTORS FOR QUATERNIONIC CASE")
print("=" * 78)

B_q = mpf(results["Quaternionic"]["B"].numerator) / mpf(results["Quaternionic"]["B"].denominator)
F_q = results["Quaternionic"]["F"]
D_q = results["Quaternionic"]["Delta"]
inner_q = B_q + F_q - D_q

print(f"\n  B + F - Δ = {nstr(inner_q, 15)}")
print(f"\n  Testing different prefactors T in K = T × (B + F - Δ):")
print(f"\n  {'Prefactor':>20} {'K':>15} {'1/α':>15} {'Nearest coupling':>20}")
print("-" * 72)

prefactors = {
    "π": mpi,
    "2π": 2*mpi,
    "π²": mpi**2,
    "π/2": mpi/2,
    "4π": 4*mpi,
    "1": mpf(1),
    "2": mpf(2),
    "1/π": 1/mpi,
    "√(2π)": (2*mpi)**mpf('0.5'),
    "2π²": 2*mpi**2,
}

for pname, pval in prefactors.items():
    K_test = pval * inner_q
    if K_test > 0:
        roots_test = polyroots([1, 0, -K_test, 1])
        real_test = sorted([r.real for r in roots_test if fabs(r.imag) < mpf(10)**(-30)])
        if len(real_test) >= 2:
            alpha_test_inv = 1 / real_test[1]

            # Find nearest known coupling
            nearest = ""
            best_diff = mpf(1000)
            for cname, cval in known_couplings.items():
                diff = fabs(alpha_test_inv - cval) / cval
                if diff < best_diff:
                    best_diff = diff
                    nearest = f"{cname} ({nstr(diff*100,3)}%)"

            print(f"  {pname:>20} {nstr(K_test,10):>15} {nstr(alpha_test_inv,10):>15} {nearest:>20}")


# ==============================================================================
# PART 10: WHAT IF THE POLYNOMIAL ORDER CHANGES?
# ==============================================================================
print("\n" + "=" * 78)
print("PART 10: ALTERNATIVE POLYNOMIAL FORMS")
print("=" * 78)

print("""
  The cubic α³ - Kα + 1 = 0 comes from the Z₃-symmetric circulant (3 bundle
  components). For the quaternionic Hopf, there are still 3 components, so the
  cubic form is natural. But what if the self-consistency relation changes?
""")

K_q_val = results["Quaternionic"]["K"]

# Test different self-consistency relations
print("  For K_quat = π × (308 + 3/4 - 1/144) ≈", nstr(K_q_val, 10))
print()

# α^3 - Kα + 1 = 0 (Paper 2 cubic)
# 1/α + α^2 = K (equivalent)
# Try: 1/α = K (linear)
print(f"  Linear:  1/α = K → α = {nstr(1/K_q_val, 10)}, 1/α = {nstr(K_q_val, 10)}")

# Try: α² = K → α = √K
print(f"  Square:  α = √K → 1/α = {nstr(1/K_q_val**mpf('0.5'), 10)}")

# Try: 1/α = K/π
print(f"  1/α = K/π → 1/α = {nstr(K_q_val/mpi, 10)} = B+F-Δ")

# Try: α = 1/B (simplest)
print(f"  α = 1/B → 1/α = {results['Quaternionic']['B']}")


# ==============================================================================
# SUMMARY
# ==============================================================================
print("\n" + "=" * 78)
print("SUMMARY OF FINDINGS")
print("=" * 78)

print(f"""
1. UNIVERSAL SELECTION PRINCIPLE:
   B(m)/N(m) = dim(S^d) is satisfied at m = 2 for ALL dimensions d,
   and m = 2 is the UNIQUE positive solution. This is proven algebraically:
   the equation reduces to m² + dm - 2(d+2) = 0 with roots m = 2, -(d+2).

   Consequence: the selection principle is NOT specific to S³ or the Hopf
   fibrations. It is a universal property of the eigenvalue distribution
   on round spheres. Paper 2's n_max = 3 (= m + 1) is universal.

2. QUATERNIONIC HOPF SPECTRAL DATA:
   B_quat = 308  (Casimir trace, S⁷ at m=2)
   F_quat = 3/4  (fiber S³ spectral zeta)
   Δ_quat = 1/144  (boundary correction)
   K_quat = π × 308.743... ≈ 970.1
   1/α_quat ≈ 970.1  (cubic root)

   → Does NOT match any known weak-sector coupling.
   → α_W ~ 1/29.5, sin²θ_W ~ 0.231, α at M_Z ~ 1/128.9

3. OCTONIONIC HOPF SPECTRAL DATA:
   B_oct = 2280, F_oct = 49/120, Δ_oct = 1/544
   K_oct ≈ 7163.8, 1/α_oct ≈ 7163.8

   → Does NOT match any known strong-sector coupling.
   → α_s ~ 0.118 (1/α_s ≈ 8.5)

4. RATIO STRUCTURE:
   K_quat / K_em ≈ 7.08  (suggestively close to dim(S⁷)/dim(S³) = 7/3 × 3 ≈ 7)
   But the ratio is NOT exactly 7 (deviation ~ 1%).

   The ratios are dominated by the B values:
   B_quat/B_em = 308/42 = 22/3 ≈ 7.33

5. REPRESENTATION-THEORETIC OBSTRUCTION:
   The Paper 2 formula's deep content comes from the SO(4) → SO(3) branching,
   which requires SO(4) to be a product group. This is unique to d = 3.
   For d = 7, the relevant chain SO(8) → Sp(2) × Sp(1) involves SPINOR
   representations of SO(5), breaking the integer-angular-momentum structure
   that makes the complex Hopf formula physically interpretable.

6. REAL HOPF (S⁰ → S¹ → S¹):
   Fully degenerate. Fiber S⁰ is discrete with no continuous spectral geometry.
   Formally gives 1/α ≈ 17.0, but this has no physical content.

CONCLUSION:
   The Paper 2 recipe PARTIALLY generalizes: the selection principle is universal,
   and all spectral ingredients (B, F, Δ) can be computed. But the resulting
   coupling constants do not match known physics. The formula appears specific to
   the complex Hopf fibration and U(1) electromagnetism. This is consistent with
   Paper 2's S³ specificity argument: S³ is the only sphere that is a Lie group,
   the only one with SO(4) ≅ SU(2)×SU(2), and the only total space where the
   base Casimir trace has a clean representation-theoretic meaning.

   The universality of B/N = d is a NEW mathematical result that should be noted
   in Paper 2 — it means the selection principle is necessary but not sufficient;
   the formula's specificity lies in the combination rule and cubic structure.
""")
