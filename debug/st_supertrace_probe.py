"""
Supertrace spectral action probe — ST-A/B combined.

Tests whether K/π = B + F − Δ arises from a spectral action supertrace
(scalar − spinor) on S³ at finite cutoff n_max=3.

Key hypothesis: the minus sign on Δ is the standard (-1)^F boson-fermion
sign from the Connes-Chamseddine spectral action.
"""

import numpy as np
import json
from fractions import Fraction
try:
    import mpmath
    mpmath.mp.dps = 50
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

# =============================================================
# Constants
# =============================================================
B_target = 42
F_target = float(mpmath.pi**2 / 6) if HAS_MPMATH else np.pi**2 / 6
Delta_target = Fraction(1, 40)
K_over_pi = B_target + F_target - float(Delta_target)  # ~43.6199
K_target = K_over_pi * np.pi  # ~137.036

print("=" * 72)
print("SUPERTRACE SPECTRAL ACTION PROBE")
print("=" * 72)
print(f"Targets: B = {B_target}, F = π²/6 = {F_target:.15f}")
print(f"         Δ = 1/40 = {float(Delta_target):.15f}")
print(f"         K/π = B + F - Δ = {K_over_pi:.15f}")
print(f"         K = π(B+F-Δ) = {K_target:.15f}")
print(f"         1/α (CODATA 2018) = 137.035999084")
print()

# =============================================================
# S³ spectra
# =============================================================

def scalar_eigenvalue(n):
    """Scalar Laplacian on S³: λ_n = n² - 1, n = 1, 2, 3, ..."""
    return n * n - 1

def scalar_degeneracy(n):
    """Degeneracy of scalar shell n: g_n = n²."""
    return n * n

def dirac_eigenvalue_sq(n_ch):
    """D² eigenvalue at CH index n: (n + 3/2)²."""
    return (n_ch + 1.5) ** 2

def dirac_degeneracy(n_ch):
    """Single-chirality Dirac degeneracy: g_n = 2(n+1)(n+2)."""
    return 2 * (n_ch + 1) * (n_ch + 2)

def dirac_degeneracy_total(n_ch):
    """Total D² degeneracy (both chiralities): 4(n+1)(n+2)."""
    return 4 * (n_ch + 1) * (n_ch + 2)

def casimir_weight_per_shell(n):
    """Per-shell Casimir sum for shell n: Σ_{l=0}^{n-1} (2l+1)l(l+1)."""
    return sum((2*l + 1) * l * (l + 1) for l in range(n))

def casimir_weight(m):
    """Cumulative Casimir B(m) = Σ_{n=1}^{m} Σ_{l=0}^{n-1} (2l+1)l(l+1)."""
    return sum(casimir_weight_per_shell(n) for n in range(1, m + 1))

def casimir_weight_closed(m):
    """Closed form: B(m) = m(m-1)(m+1)(m+2)(2m+1)/20."""
    return m * (m - 1) * (m + 1) * (m + 2) * (2 * m + 1) // 20

# Verify closed form
for n in range(1, 10):
    assert casimir_weight(n) == casimir_weight_closed(n), f"Mismatch at n={n}"
print("✓ Casimir closed form verified for n=1..9")

# =============================================================
# TEST 1: Euler-Maclaurin on Dirac mode count
# =============================================================
print("\n" + "=" * 72)
print("TEST 1: Euler-Maclaurin boundary correction for Dirac mode count")
print("=" * 72)

# g(n) = 4(n+1)(n+2) for CH index n = 0, 1, 2, ...
# Exact sum: N_D(N) = Σ_{n=0}^{N} 4(n+1)(n+2) = 4(N+1)(N+2)(N+3)/3

N_cut = 3  # CH cutoff

# Exact sum
exact_sum = sum(dirac_degeneracy_total(n) for n in range(N_cut + 1))
exact_formula = 4 * (N_cut + 1) * (N_cut + 2) * (N_cut + 3) // 3
print(f"Exact sum Σ_{{n=0}}^{{{N_cut}}} 4(n+1)(n+2) = {exact_sum}")
print(f"Closed form 4(N+1)(N+2)(N+3)/3 = {exact_formula}")
assert exact_sum == exact_formula

# Euler-Maclaurin: Σ_{n=0}^{N} g(n) = ∫₀^N g(x)dx + [g(0)+g(N)]/2 + EM corrections
# g(x) = 4(x+1)(x+2) = 4(x² + 3x + 2)
# ∫₀^N g(x)dx = 4[x³/3 + 3x²/2 + 2x]₀^N = 4(N³/3 + 3N²/2 + 2N)
integral = 4 * (N_cut**3 / 3 + 3 * N_cut**2 / 2 + 2 * N_cut)
g_0 = dirac_degeneracy_total(0)  # 4·1·2 = 8
g_N = dirac_degeneracy_total(N_cut)  # 4·4·5 = 80
boundary_correction = (g_0 + g_N) / 2  # (8 + 80)/2 = 44

# g'(x) = 4(2x + 3), g''(x) = 8, g'''(x) = 0 (polynomial of degree 2)
# EM corrections: B₂/2! [g'(N) - g'(0)] = (1/6)/2 · [4(2N+3) - 4·3]
#               = 1/12 · 4 · 2N = 2N/3
g_prime_N = 4 * (2 * N_cut + 3)  # 4·9 = 36
g_prime_0 = 4 * 3  # 12
em_B2 = Fraction(1, 12) * (g_prime_N - g_prime_0)  # (1/12)(36-12) = 2

# Higher Bernoulli terms: g''(x) = 8, g'''(x) = 0
# B₄/4! [g'''(N) - g'''(0)] = 0 (g is degree 2, all derivatives ≥ 3 vanish)

print(f"\nEM decomposition of Σ_{{n=0}}^{{{N_cut}}} g(n) = {exact_sum}:")
print(f"  Integral ∫₀^{N_cut} g(x)dx = {integral:.6f}")
print(f"  Boundary [g(0)+g({N_cut})]/2 = ({g_0}+{g_N})/2 = {boundary_correction:.1f}")
print(f"  B₂ correction = {float(em_B2):.6f}")
print(f"  Higher terms = 0 (g is degree 2)")
print(f"  Total EM = {integral + boundary_correction + float(em_B2):.6f}")
print(f"  Exact    = {exact_sum}")

print(f"\n  g({N_cut})/2 = {g_N/2:.1f}")
print(f"  g₃^Dirac = {dirac_degeneracy(N_cut)} = Δ⁻¹ = {1/float(Delta_target):.0f}")
print(f"  MATCH: g({N_cut})/2 = g₃^Dirac? {g_N/2 == dirac_degeneracy(N_cut)}")

# The single-chirality boundary correction
print(f"\n  Single-chirality g({N_cut}) = {dirac_degeneracy(N_cut)}")
print(f"  Half of that = {dirac_degeneracy(N_cut)/2}")
print(f"  Total (both chiralities) g({N_cut}) = {g_N}")
print(f"  Half of total = {g_N/2} = {int(g_N/2)}")

# KEY: the EM boundary correction at the upper endpoint is g(N)/2
# For the TOTAL (both chiralities): g(3)/2 = 80/2 = 40 = g₃^Dirac = Δ⁻¹
# For SINGLE chirality: 2(4)(5)/2 = 20 ≠ 40
print(f"\n  ★ EM upper boundary term g_total(3)/2 = {g_N//2} = Δ⁻¹ = {int(1/float(Delta_target))}")
print(f"    VERDICT: The Euler-Maclaurin boundary correction at n=3 for the")
print(f"    D² mode count IS exactly Δ⁻¹ = g₃^Dirac = 40.")

# =============================================================
# TEST 2: Spectral action scan — supertrace
# =============================================================
print("\n" + "=" * 72)
print("TEST 2: Supertrace spectral action scan")
print("=" * 72)

# Cutoff functions
def f_sharp(x):
    return 1.0 if x <= 1.0 else 0.0

def f_cc(x):
    """Connes-Chamseddine (1-x)²₊"""
    return max(0.0, 1.0 - x) ** 2

def f_exp(x):
    return np.exp(-x)

def f_rational(x):
    return 1.0 / (1.0 + x) ** 2

cutoff_funcs = {
    'sharp': f_sharp,
    'CC (1-x)²': f_cc,
    'exp(-x)': f_exp,
    '1/(1+x)²': f_rational,
}

N_sum = 5000  # terms for convergence

def spectral_action_scalar(f, Lambda_sq, N=N_sum):
    """S_scalar = Σ_{n=1}^{N} n² · f((n²-1)/Λ²)"""
    total = 0.0
    for n in range(1, N + 1):
        x = (n * n - 1) / Lambda_sq
        total += n * n * f(x)
    return total

def spectral_action_dirac(f, Lambda_sq, N=N_sum):
    """S_dirac = Σ_{n=0}^{N} 4(n+1)(n+2) · f((n+3/2)²/Λ²)"""
    total = 0.0
    for n in range(N + 1):
        x = (n + 1.5) ** 2 / Lambda_sq
        total += 4 * (n + 1) * (n + 2) * f(x)
    return total

# Scan over Λ² values
Lambda_sq_candidates = np.concatenate([
    np.array([8.0, 9.0, 81/4, 20.25, 3.75]),  # specific candidates
    np.linspace(1.0, 30.0, 100),
])
Lambda_sq_candidates = np.unique(Lambda_sq_candidates)

# Normalization constants to try
norm_constants = {
    '1': 1.0,
    '1/2': 0.5,
    '1/4': 0.25,
    '1/8': 0.125,
    '1/3': 1/3,
}

best_hits = []

for fname, f in cutoff_funcs.items():
    for Lsq in Lambda_sq_candidates:
        S_s = spectral_action_scalar(f, Lsq, N=2000)
        S_d = spectral_action_dirac(f, Lsq, N=2000)

        for cname, c in norm_constants.items():
            supertrace = S_s - c * S_d

            # Check against targets
            for target_name, target_val in [
                ('K/π', K_over_pi),
                ('B+F', B_target + F_target),
                ('B', float(B_target)),
                ('K', K_target),
            ]:
                if abs(target_val) > 1e-10:
                    rel_err = abs(supertrace - target_val) / abs(target_val)
                    if rel_err < 0.05:  # within 5%
                        best_hits.append({
                            'f': fname,
                            'Λ²': float(Lsq),
                            'c': cname,
                            'S_scalar': S_s,
                            'S_dirac': S_d,
                            'supertrace': supertrace,
                            'target': target_name,
                            'target_val': target_val,
                            'rel_err': rel_err,
                        })

# Also check ratios
ratio_hits = []
for fname, f in cutoff_funcs.items():
    for Lsq in Lambda_sq_candidates:
        S_s = spectral_action_scalar(f, Lsq, N=2000)
        S_d = spectral_action_dirac(f, Lsq, N=2000)
        if abs(S_d) > 1e-10:
            ratio = S_s / S_d
            for target_name, target_val in [('K/π', K_over_pi), ('B', 42.0), ('1/4', 0.25)]:
                if abs(target_val) > 1e-10:
                    rel_err = abs(ratio - target_val) / abs(target_val)
                    if rel_err < 0.02:
                        ratio_hits.append({
                            'f': fname, 'Λ²': float(Lsq),
                            'ratio': ratio, 'target': target_name,
                            'rel_err': rel_err,
                        })

# Sort by relative error
best_hits.sort(key=lambda h: h['rel_err'])
print(f"\nTop 10 supertrace hits (within 5% of targets):")
print(f"{'f':<15} {'Λ²':<8} {'c':<6} {'super':<15} {'target':<8} {'val':<15} {'rel_err':<10}")
for h in best_hits[:10]:
    print(f"{h['f']:<15} {h['Λ²']:<8.3f} {h['c']:<6} {h['supertrace']:<15.6f} "
          f"{h['target']:<8} {h['target_val']:<15.6f} {h['rel_err']:<10.6f}")

if ratio_hits:
    ratio_hits.sort(key=lambda h: h['rel_err'])
    print(f"\nRatio hits (S_scalar/S_dirac within 2% of targets):")
    for h in ratio_hits[:5]:
        print(f"  f={h['f']}, Λ²={h['Λ²']:.3f}: ratio={h['ratio']:.6f} "
              f"vs {h['target']}={h.get('target', '')}, err={h['rel_err']:.4%}")

# =============================================================
# TEST 3: Per-shell supertrace at n_max=3
# =============================================================
print("\n" + "=" * 72)
print("TEST 3: Per-shell analysis at n_max=3")
print("=" * 72)

print("\nScalar shells (Fock n = 1, 2, 3):")
print(f"  {'n':<4} {'λ_n':<8} {'g_n=n²':<8} {'B_n':<8} {'B_n/n²':<10}")
B_cumul = 0
N_cumul = 0
for n in range(1, 4):
    lam = scalar_eigenvalue(n)
    g = scalar_degeneracy(n)
    Bn_shell = casimir_weight_per_shell(n)
    B_cumul += Bn_shell
    N_cumul += g
    print(f"  {n:<4} {lam:<8} {g:<8} {Bn_shell:<8} {Bn_shell/g if g > 0 else 0:<10.4f}")
print(f"  Cumulative: N = {N_cumul}, B = {B_cumul}, B/N = {B_cumul/N_cumul:.6f}")

print(f"\nDirac shells (CH n = 0, 1, 2, 3):")
print(f"  {'n_CH':<6} {'|λ|':<8} {'g_D':<8} {'g_D_total':<10}")
D_cumul = 0
for n in range(4):
    lam = n + 1.5
    g = dirac_degeneracy(n)
    gt = dirac_degeneracy_total(n)
    D_cumul += gt
    print(f"  {n:<6} {lam:<8.1f} {g:<8} {gt:<10}")
print(f"  Cumulative D² modes: {D_cumul}")

# The per-shell supertrace with various normalizations
print(f"\nPer-shell supertrace δ_n = B_n - c · g_n^Dirac:")
for c_val, c_name in [(1, '1'), (Fraction(1,2), '1/2'), (Fraction(1,4), '1/4'),
                       (Fraction(21,40), '21/40'), (Fraction(1,20), '1/20')]:
    total = 0
    for n_fock in range(1, 4):
        n_ch = n_fock - 1  # CH index
        Bn = casimir_weight_per_shell(n_fock)
        gD = dirac_degeneracy(n_ch)
        delta = Bn - float(c_val) * gD
        total += delta
    # Add tail estimate for F
    tail_F = F_target - sum(1/n**2 for n in range(1, 4))  # F - H_3
    result = total + tail_F
    print(f"  c={c_name}: Σδ = {total:.4f}, + F_tail = {result:.6f} "
          f"(vs K/π = {K_over_pi:.6f}, diff = {result - K_over_pi:.6f})")

# =============================================================
# TEST 4: Casimir-weighted spectral action tail
# =============================================================
print("\n" + "=" * 72)
print("TEST 4: Casimir-weighted spectral zeta tail")
print("=" * 72)

# Z_C(s) = Σ_{n=1}^∞ B_n / (n²-1)^s
# Z_C(s, 3) = B_1/(0)^s + B_2/3^s + B_3/8^s
# Note: B_1 = 0 (only l=0 in n=1 shell, Casimir = 0)
# So Z_C(s, 3) = B_2/3^s + B_3/8^s = 6/3^s + 36/8^s

if HAS_MPMATH:
    print("\nCasimir spectral zeta Z_C(s) = Σ B_n / (n²-1)^s:")
    print(f"  B_1={casimir_weight(1)}, B_2={casimir_weight(2)}, B_3={casimir_weight(3)}")

    for s in [1, 2, 3, 4]:
        # Truncated (n=2,3 only; n=1 has λ=0)
        Z_trunc = sum(casimir_weight(n) / (n*n - 1)**s for n in range(2, 4))
        # Tail (n=4 to large N)
        Z_tail = sum(casimir_weight_closed(n) / (n*n - 1)**s for n in range(4, 10000))
        Z_total = Z_trunc + Z_tail
        print(f"  s={s}: Z_trunc(3) = {Z_trunc:.10f}, Z_tail = {Z_tail:.10f}, "
              f"Z_total = {Z_total:.10f}")

        # Check against F and B targets
        if abs(Z_tail - F_target) / F_target < 0.5:
            print(f"        ★ Z_tail / F = {Z_tail / F_target:.6f}")
        if abs(Z_total - (B_target + F_target)) < 10:
            print(f"        ★ Z_total - (B+F) = {Z_total - (B_target + F_target):.6f}")

# =============================================================
# TEST 5: Fock degeneracy Dirichlet — EM decomposition
# =============================================================
print("\n" + "=" * 72)
print("TEST 5: Euler-Maclaurin on ζ(2) = Σ 1/n² with N=3 cutoff")
print("=" * 72)

# f(n) = 1/n²
# Σ_{n=1}^{3} 1/n² = 1 + 1/4 + 1/9 = 49/36
H3 = Fraction(1, 1) + Fraction(1, 4) + Fraction(1, 9)
print(f"Partial sum H_3 = {H3} = {float(H3):.15f}")
print(f"ζ(2) = π²/6 = {F_target:.15f}")
print(f"Tail = ζ(2) - H_3 = {F_target - float(H3):.15f}")

# EM decomposition of Σ_{n=1}^{N} 1/n²
# ∫₁^N 1/x² dx = 1 - 1/N
# [f(1) + f(N)]/2 = (1 + 1/N²)/2
# B₂/2! [f'(N) - f'(1)] = (1/12)[-2/N³ + 2] = (1/6)(1 - 1/N³)
# B₄/4! [f'''(N) - f'''(1)] = (-1/720)[-24/N⁵ + 24] = (1/30)(1 - 1/N⁵)

N = 3
em_integral = 1.0 - 1.0/N  # = 2/3
em_boundary = (1.0 + 1.0/N**2) / 2  # = (1 + 1/9)/2 = 5/9
em_B2 = (1.0/6) * (1.0 - 1.0/N**3)  # = (1/6)(1 - 1/27) = 26/162
em_B4 = (1.0/30) * (1.0 - 1.0/N**5)  # = (1/30)(1 - 1/243) = 242/7290

em_sum = em_integral + em_boundary + em_B2 + em_B4
print(f"\nEM decomposition of Σ_{{n=1}}^{{{N}}} 1/n²:")
print(f"  Integral ∫₁^{N} 1/x² dx = {em_integral:.10f}")
print(f"  Boundary [f(1)+f({N})]/2 = {em_boundary:.10f}")
print(f"  B₂ correction = {em_B2:.10f}")
print(f"  B₄ correction = {em_B4:.10f}")
print(f"  EM total (4 terms) = {em_sum:.10f}")
print(f"  Exact H_3 = {float(H3):.10f}")
print(f"  EM error = {abs(em_sum - float(H3)):.2e}")

# Check EM boundary term against Δ
print(f"\n  EM upper boundary f(N)/2 = 1/(2·{N}²) = {1/(2*N**2):.10f}")
print(f"  Δ = 1/40 = {float(Delta_target):.10f}")
print(f"  MATCH? {abs(1/(2*N**2) - float(Delta_target)) < 1e-10}")
print(f"  f(N)/2 = 1/18 vs Δ = 1/40 — NO MATCH")

# =============================================================
# TEST 6: Direct spectral action at specific Λ values
# =============================================================
print("\n" + "=" * 72)
print("TEST 6: Spectral action values at specific Λ")
print("=" * 72)

# For exp(-x) cutoff, the spectral action has a closed-form heat kernel expansion
# S(Λ) = Σ g_n exp(-λ_n/Λ²)

# Scalar: S_s(Λ) = Σ_{n=1}^∞ n² exp(-(n²-1)/Λ²) = e^{1/Λ²} Σ n² e^{-n²/Λ²}
# This is related to the Jacobi theta function θ₃

# Just compute numerically
print(f"\nexp(-x) cutoff — scalar and Dirac spectral actions:")
print(f"{'Λ²':<10} {'S_scalar':<18} {'S_dirac':<18} {'S_s/S_d':<12} {'S_s - S_d/4':<15}")
for Lsq in [4, 8, 9, 12, 16, 20.25, 25, 36, 81/4]:
    S_s = spectral_action_scalar(f_exp, Lsq, N=500)
    S_d = spectral_action_dirac(f_exp, Lsq, N=500)
    print(f"{Lsq:<10.2f} {S_s:<18.8f} {S_d:<18.8f} {S_s/S_d:<12.8f} {S_s - S_d/4:<15.8f}")

# =============================================================
# TEST 7: The key insight test — B/N = dim(S³)
# =============================================================
print("\n" + "=" * 72)
print("TEST 7: Selection principle B(m)/N(m) = dim(S³) = 3")
print("=" * 72)

print(f"\n{'m':<4} {'N(m)':<10} {'B(m)':<10} {'B/N':<12} {'3(m-1)(m+2)/10':<18}")
for m in range(1, 9):
    Nm = m * (m + 1) * (2 * m + 1) // 6
    Bm = casimir_weight_closed(m)
    ratio = Bm / Nm if Nm > 0 else 0
    formula = 3 * (m - 1) * (m + 2) / 10
    match = "★ = dim(S³)" if abs(ratio - 3.0) < 1e-10 else ""
    print(f"{m:<4} {Nm:<10} {Bm:<10} {ratio:<12.6f} {formula:<18.6f} {match}")

# =============================================================
# TEST 8: "What spectral quantity gives B=42 as truncated + F as tail?"
# =============================================================
print("\n" + "=" * 72)
print("TEST 8: Looking for spectral sum that decomposes as B + F")
print("=" * 72)

# We need a sequence a_n such that:
# Σ_{n=1}^{3} a_n = 42 = B
# Σ_{n=4}^{∞} a_n = π²/6 = F (after regularization)
#
# One candidate: a_n = B_n for n ≤ 3, a_n = 1/n² for n ≥ 4?
# But that's two different sequences stitched together.
#
# Better question: is there a single function g(n) such that
# g(n) ~ B_n for small n and g(n) ~ 1/n² for large n?
#
# B_n = n(n-1)(n+1)(n+2)(2n+1)/20 ~ n⁵/10 for large n
# 1/n² ~ n⁻²
# These have totally different asymptotics. No single g(n) works.

# BUT: what if the spectral action WEIGHT f((n²-1)/Λ²) converts B_n to 1/n²?
# Then: B_n · f((n²-1)/Λ²) ~ 1/n² for large n
# f(x) ~ 1/(n⁷) for x ~ n² (to kill the n⁵ growth of B_n and get n⁻²)
# f(x) ~ x⁻⁷/² — a specific power-law cutoff

# Test: Casimir-weighted spectral action with power-law cutoff f(x) = x^{-s}
# S_C(s) = Σ_{n=2}^∞ B_n · (n²-1)^{-s}
# For convergence: B_n ~ n⁵, (n²-1)^{-s} ~ n^{-2s}, need 5-2s < -1 → s > 3
print("\nCasimir-weighted spectral zeta Σ B_n · (n²-1)^{-s}:")
print(f"{'s':<6} {'Z_C(s,3)':<18} {'Z_C(s,∞)':<18} {'tail':<18} {'tail/F':<12}")

for s_val in [3.5, 4, 4.5, 5, 6]:
    Z3 = sum(casimir_weight_closed(n) / (n*n - 1)**s_val for n in range(2, 4))
    Z_inf = sum(casimir_weight_closed(n) / (n*n - 1)**s_val for n in range(2, 20000))
    tail = Z_inf - Z3
    print(f"{s_val:<6.1f} {Z3:<18.10f} {Z_inf:<18.10f} {tail:<18.10f} {tail/F_target:<12.6f}")

# =============================================================
# TEST 9: Seeley-DeWitt coefficient supertrace
# =============================================================
print("\n" + "=" * 72)
print("TEST 9: Seeley-DeWitt coefficient analysis")
print("=" * 72)

# Scalar Laplacian on unit S³ (d=3, R=6):
# a₀^s = vol(S³) / (4π)^{3/2} = 2π² / (4π)^{3/2} = 2π² / (8π√π) = π^{1/2}/4 = √π/4
# a₁^s = (R/6) · a₀^s = 1 · √π/4 = √π/4

# Dirac D² on unit S³:
# a₀^d = √π (from qed_vacuum_polarization.py)
# a₁^d = √π
# a₂^d = √π/8

a0_s = np.sqrt(np.pi) / 4
a1_s = np.sqrt(np.pi) / 4  # R/6 = 1 on unit S³
a0_d = np.sqrt(np.pi)
a1_d = np.sqrt(np.pi)
a2_d = np.sqrt(np.pi) / 8

print(f"Scalar: a₀ = √π/4 = {a0_s:.10f}, a₁ = √π/4 = {a1_s:.10f}")
print(f"Dirac:  a₀ = √π   = {a0_d:.10f}, a₁ = √π   = {a1_d:.10f}, a₂ = √π/8 = {a2_d:.10f}")
print(f"Ratio a₀^d/a₀^s = {a0_d/a0_s:.1f} (factor of 4 — Dirac has 4× modes per eigenvalue)")

# Supertrace of SD coefficients: Δa_k = a_k^s - c · a_k^d
# For the ratio to be meaningful, need correct normalization
# The standard CC spectral action: Tr[f(D²/Λ²)] = Σ f_k a_k Λ^{d-2k}
# On S³ (d=3): f₁ Λ³ a₀ + f₀ Λ a₁ + f_{-1} Λ⁻¹ a₂ + ...

print(f"\nSD expansion S(Λ) = Σ f_k Λ^{{3-2k}} a_k for exp(-x) cutoff:")
print(f"  f_k = ∫₀^∞ e^{{-u}} u^{{(1-2k)/2}} du = Γ((3-2k)/2)")
print(f"  f₁ = Γ(1/2) = √π = {np.sqrt(np.pi):.10f}")
print(f"  f₀ = Γ(-1/2) = -2√π = {-2*np.sqrt(np.pi):.10f}")

f1_moment = np.sqrt(np.pi)  # Γ(1/2)
f0_moment = -2 * np.sqrt(np.pi)  # Γ(-1/2)

print(f"\nSD spectral action for exp(-x) cutoff:")
print(f"{'Λ':<8} {'S_s(SD)':<15} {'S_d(SD)':<15} {'S_s-S_d/4':<15}")
for Lambda in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]:
    S_s_sd = f1_moment * Lambda**3 * a0_s + f0_moment * Lambda * a1_s
    S_d_sd = f1_moment * Lambda**3 * a0_d + f0_moment * Lambda * a1_d
    print(f"{Lambda:<8.1f} {S_s_sd:<15.6f} {S_d_sd:<15.6f} {S_s_sd - S_d_sd/4:<15.6f}")

# =============================================================
# TEST 10: The "inverse solve" — what Λ gives K/π from each spectral action?
# =============================================================
print("\n" + "=" * 72)
print("TEST 10: Inverse solve — what Λ makes S(Λ) = K/π?")
print("=" * 72)

from scipy.optimize import brentq

for fname, f in [('exp(-x)', f_exp), ('CC (1-x)²', f_cc), ('1/(1+x)²', f_rational)]:
    def objective_scalar(log_Lsq):
        Lsq = np.exp(log_Lsq)
        return spectral_action_scalar(f, Lsq, N=1000) - K_over_pi

    def objective_dirac(log_Lsq):
        Lsq = np.exp(log_Lsq)
        return spectral_action_dirac(f, Lsq, N=1000) - K_over_pi

    # Find Λ² where scalar spectral action = K/π
    try:
        log_Lsq_s = brentq(objective_scalar, -2, 10)
        Lsq_s = np.exp(log_Lsq_s)
        print(f"  {fname}: S_scalar(Λ²={Lsq_s:.6f}) = K/π = {K_over_pi:.6f}")
        print(f"    Λ = {np.sqrt(Lsq_s):.6f}, Λ² = {Lsq_s:.6f}")
        # Check if Λ² is a recognizable number
        for name, val in [('8=n²-1@3', 8), ('9=n²@3', 9), ('81/4=(9/2)²', 20.25),
                          ('3/2', 1.5), ('15/4', 3.75), ('4', 4), ('16', 16)]:
            if abs(Lsq_s - val) / val < 0.01:
                print(f"    ★ Close to {name} = {val}")
    except (ValueError, Exception) as e:
        print(f"  {fname}: scalar solve failed ({e})")

    try:
        log_Lsq_d = brentq(objective_dirac, -2, 10)
        Lsq_d = np.exp(log_Lsq_d)
        print(f"  {fname}: S_dirac(Λ²={Lsq_d:.6f}) = K/π = {K_over_pi:.6f}")
    except (ValueError, Exception) as e:
        pass

# =============================================================
# SUMMARY
# =============================================================
print("\n" + "=" * 72)
print("SUMMARY OF PROBE RESULTS")
print("=" * 72)

results = {
    'test1_em_boundary': {
        'verdict': 'POSITIVE',
        'finding': 'EM boundary correction at n=3 for D² mode count = g_3^Dirac/2 per chirality; '
                   'total (both chiralities) g(3)/2 = 40 = Δ⁻¹ exactly',
    },
    'test2_supertrace_scan': {
        'verdict': 'PENDING',
        'finding': f'Best hits from scan (see table above), {len(best_hits)} within 5%',
    },
    'test5_em_on_zeta2': {
        'verdict': 'NEGATIVE',
        'finding': 'EM boundary correction 1/(2N²) = 1/18 ≠ Δ = 1/40',
    },
    'test7_selection': {
        'verdict': 'CONFIRMED',
        'finding': 'B(3)/N(3) = 3 = dim(S³); unique integer hit',
    },
    'test8_decomposition': {
        'verdict': 'NEGATIVE',
        'finding': 'No single spectral quantity decomposes as B (truncated) + F (tail) — '
                   'B_n ~ n⁵ and 1/n² have incompatible asymptotics',
    },
}

for test, r in results.items():
    print(f"  {test}: {r['verdict']} — {r['finding']}")

# Save results
output = {
    'targets': {'B': B_target, 'F': F_target, 'Delta': float(Delta_target),
                'K_over_pi': K_over_pi, 'K': K_target},
    'test1_em': {
        'g_N_total': int(g_N), 'boundary_half': int(g_N // 2),
        'delta_inverse': int(1 / float(Delta_target)),
        'match': g_N // 2 == int(1 / float(Delta_target)),
    },
    'best_supertrace_hits': best_hits[:20],
    'results': results,
}

with open('debug/data/st_supertrace_probe.json', 'w') as fout:
    json.dump(output, fout, indent=2, default=str)
print(f"\nData saved to debug/data/st_supertrace_probe.json")
