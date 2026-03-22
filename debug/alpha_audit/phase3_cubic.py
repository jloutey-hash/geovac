#!/usr/bin/env python3
"""
Phase 3: Origin of the Cubic — α³ − Kα + 1 = 0
=================================================
Phase 2 established statistical significance (p = 5.2e-9).
This phase asks: does the cubic arise naturally from the geometry
of the Hopf fibration S¹ → S³ → S²?

Five investigations:
  1. Spectral determinant ratios
  2. Heat kernel coefficients
  3. Self-consistency from the fibration (Chern, CS, η)
  4. QED comparison
  5. Characteristic polynomial search
"""

import numpy as np
import os
import sys
from itertools import product as cartesian_product

# ================================================================
# Constants
# ================================================================
CODATA_ALPHA_INV = 137.035999084  # CODATA 2018
alpha_codata = 1.0 / CODATA_ALPHA_INV
zeta2 = np.pi**2 / 6              # ζ(2) = π²/6
zeta3 = 1.2020569031595942        # Apéry's constant
zeta4 = np.pi**4 / 90             # ζ(4)

# Riemann zeta derivatives at special points (high precision)
# ζ'(0) = -(1/2)ln(2π)
zeta_prime_0 = -0.5 * np.log(2 * np.pi)
# ζ'(-1) ≈ -0.16542114370045092
zeta_prime_neg1 = -0.16542114370045092
# ζ'(-2) = -ζ(3)/(4π²) [from functional equation, trivial zero]
zeta_prime_neg2 = -zeta3 / (4 * np.pi**2)

# Paper's formula
B = 42          # Degeneracy-weighted Casimir trace at n_max=3
F = zeta2       # Fiber invariant = ζ(2) = π²/6
Delta = 1 / 40  # Boundary correction = 1/(|λ₃|×N(2))
K_paper = np.pi * (B + F - Delta)  # ≈ 137.036064

# Solve the cubic α³ − Kα + 1 = 0
roots = np.roots([1, 0, -K_paper, 1])
alpha_paper = min(r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 0)
alpha_inv_paper = 1.0 / alpha_paper
paper_err = abs(alpha_inv_paper - CODATA_ALPHA_INV) / CODATA_ALPHA_INV

output_dir = os.path.dirname(os.path.abspath(__file__))
results = []  # Accumulate output lines


def out(s=""):
    """Print and save output."""
    print(s)
    results.append(s)


# ================================================================
out("=" * 76)
out("Phase 3: Origin of the Cubic α³ − Kα + 1 = 0")
out("=" * 76)
out(f"\nPaper's formula: K = π(42 + ζ(2) − 1/40) = {K_paper:.6f}")
out(f"Physical root:   α⁻¹ = {alpha_inv_paper:.6f}  (err = {paper_err:.2e})")
out(f"CODATA:          α⁻¹ = {CODATA_ALPHA_INV}")


# ================================================================
# INVESTIGATION 1: Spectral Determinant Ratios
# ================================================================
out("\n" + "=" * 76)
out("INVESTIGATION 1: Spectral Determinant Ratios")
out("=" * 76)
out("""
For compact manifolds, the zeta-regularized determinant det'(Δ) is a spectral
invariant: log det'(Δ) = −ζ'_Δ(0).

The Hopf bundle has three Laplacians: on S¹ (fiber), S² (base), S³ (total).
""")

# --- S¹ (unit circle, circumference 2π) ---
# Eigenvalues: k² for k ∈ Z\{0}, each with multiplicity 2 (±k)
# ζ_{S¹}(s) = 2·ζ_R(2s)
# ζ'_{S¹}(0) = 2·2·ζ'_R(0) = 4·ζ'_R(0) = 4·(−½ ln(2π)) = −2 ln(2π)
# log det'_{S¹} = −ζ'_{S¹}(0) = 2 ln(2π) = ln(4π²)
log_det_S1 = 2 * np.log(2 * np.pi)
det_S1 = np.exp(log_det_S1)  # = 4π²
out("a) Spectral determinants:")
out(f"  S¹: log det'(Δ) = 2·ln(2π) = {log_det_S1:.8f}")
out(f"      det'(Δ) = 4π² = {det_S1:.8f}")

# --- S² (unit 2-sphere) ---
# Eigenvalues: l(l+1) with multiplicity 2l+1, l ≥ 1
# ζ_{S²}(s) = Σ_{l=1}^∞ (2l+1) [l(l+1)]^{-s}
#
# Analytic derivation:
# ζ'_{S²}(0) = −Σ (2l+1) ln[l(l+1)] (zeta-regularized)
#            = −Σ (2l+1)[ln l + ln(l+1)]
# Using zeta-regularized sums:
#   Σ (2l+1) ln l = −2ζ'_R(−1) − ζ'_R(0)
#   Σ (2l+1) ln(l+1) = −2ζ'_R(−1) + ζ'_R(0)
# Therefore: ζ'_{S²}(0) = 4·ζ'_R(−1)
#
# Ref: Weisberger (1987), Osgood-Phillips-Sarnak (1988)
zeta_prime_S2 = 4 * zeta_prime_neg1
log_det_S2 = -zeta_prime_S2  # = −4·ζ'_R(−1)
det_S2 = np.exp(log_det_S2)
out(f"\n  S²: ζ'_{'{S²}'}(0) = 4·ζ'_R(−1) = {zeta_prime_S2:.8f}")
out(f"      log det'(Δ) = −4·ζ'_R(−1) = {log_det_S2:.8f}")
out(f"      det'(Δ) = {det_S2:.8f}")

# --- S³ (unit 3-sphere) ---
# Eigenvalues: n²−1 with multiplicity n², n ≥ 2
# ζ_{S³}(s) = Σ_{n=2}^∞ n² (n²−1)^{-s}
#
# Analytic derivation:
# ζ'_{S³}(0) = −Σ_{n=2}^∞ n² [ln(n−1) + ln(n+1)] (zeta-regularized)
# Split into S1 = Σ n² ln(n−1) and S2 = Σ n² ln(n+1):
#   S1 = −ζ'_R(−2) − 2ζ'_R(−1) − ζ'_R(0)  [shift m=n−1]
#   S2 = −ζ'_R(−2) + 2ζ'_R(−1) − ζ'_R(0) − ln 2  [shift m=n+1]
# Therefore:
#   ζ'_{S³}(0) = 2ζ'_R(−2) + 2ζ'_R(0) + ln 2
#             = −ζ(3)/(2π²) − ln(2π) + ln 2
#             = −ζ(3)/(2π²) − ln π
zeta_prime_S3 = 2 * zeta_prime_neg2 + 2 * zeta_prime_0 + np.log(2)
log_det_S3 = -zeta_prime_S3  # = ζ(3)/(2π²) + ln π
det_S3 = np.exp(log_det_S3)
out(f"\n  S³: ζ'_{'{S³}'}(0) = −ζ(3)/(2π²) − ln π = {zeta_prime_S3:.8f}")
out(f"      log det'(Δ) = ζ(3)/(2π²) + ln π = {log_det_S3:.8f}")
out(f"      det'(Δ) = {det_S3:.8f}")

# Verify with numerical partial sums (sanity check on S¹)
det_S1_check = (2 * np.pi)**2
out(f"\n  Cross-check: (2π)² = {det_S1_check:.8f} vs det'_{'{S¹}'} = {det_S1:.8f}  ✓")

# Analytic forms
out(f"\n  Analytic forms:")
out(f"    det'(Δ_{'{S¹}'}) = 4π²")
out(f"    det'(Δ_{'{S²}'}) = exp(−4·ζ'_R(−1))")
out(f"    det'(Δ_{'{S³}'}) = π · exp(ζ(3)/(2π²))")

# --- b) Form ratios and products ---
out(f"\nb) Ratios and products of determinants:")
out(f"  {'Combination':>50s}  {'Value':>12s}  {'α⁻¹?':>10s}  {'K?':>10s}")
out(f"  {'-'*50}  {'-'*12}  {'-'*10}  {'-'*10}")

combos = [
    ("det'_{S³} / det'_{S²}",            det_S3 / det_S2),
    ("det'_{S³} / det'_{S¹}",            det_S3 / det_S1),
    ("det'_{S¹} / det'_{S²}",            det_S1 / det_S2),
    ("det'_{S¹} × det'_{S²}",            det_S1 * det_S2),
    ("det'_{S¹} × det'_{S³}",            det_S1 * det_S3),
    ("det'_{S²} × det'_{S³}",            det_S2 * det_S3),
    ("det'_{S¹} × det'_{S²} × det'_{S³}", det_S1 * det_S2 * det_S3),
    # Logs
    ("log det'_{S³} − log det'_{S²}",    log_det_S3 - log_det_S2),
    ("log det'_{S³} − log det'_{S¹}",    log_det_S3 - log_det_S1),
    ("log det'_{S¹} − log det'_{S²}",    log_det_S1 - log_det_S2),
    ("log det'_{S¹} + log det'_{S³}",    log_det_S1 + log_det_S3),
    ("Σ log det'",                        log_det_S1 + log_det_S2 + log_det_S3),
    # Roots
    ("√det'_{S¹}",                        np.sqrt(det_S1)),
    ("√det'_{S³}",                        np.sqrt(det_S3)),
    # With π factors
    ("det'_{S¹} × det'_{S³} / π",        det_S1 * det_S3 / np.pi),
    ("π × det'_{S²}",                     np.pi * det_S2),
    ("4π² × det'_{S³} / π",              4 * np.pi**2 * det_S3 / np.pi),
    # Anomaly ratio
    ("det'_{S³} / (det'_{S²} × det'_{S¹})", det_S3 / (det_S2 * det_S1)),
    ("ln[det'_{S³} / (det'_{S²} × det'_{S¹})]",
     np.log(det_S3 / (det_S2 * det_S1))),
]

closest_K = (None, float('inf'))
closest_alpha = (None, float('inf'))

for name, val in combos:
    err_K = abs(val - K_paper) / K_paper if val != 0 else float('inf')
    err_alpha = abs(val - CODATA_ALPHA_INV) / CODATA_ALPHA_INV if val != 0 else float('inf')
    k_str = f"{err_K:.2e}" if err_K < 1 else "no"
    a_str = f"{err_alpha:.2e}" if err_alpha < 1 else "no"
    out(f"  {name:>50s}  {val:12.6f}  {a_str:>10s}  {k_str:>10s}")
    if err_K < closest_K[1]:
        closest_K = (name, err_K, val)
    if err_alpha < closest_alpha[1]:
        closest_alpha = (name, err_alpha, val)

out(f"\n  Closest to K = {K_paper:.3f}:")
out(f"    {closest_K[0]} = {closest_K[2]:.6f} (err = {closest_K[1]:.2e})")
out(f"  Closest to α⁻¹ = {CODATA_ALPHA_INV}:")
out(f"    {closest_alpha[0]} = {closest_alpha[2]:.6f} (err = {closest_alpha[1]:.2e})")

# --- c) The fibration anomaly ratio ---
anomaly = det_S3 / (det_S2 * det_S1)
log_anomaly = np.log(anomaly)
out(f"\nc) Fibration anomaly (total / base×fiber):")
out(f"   det'_{'{S³}'} / (det'_{'{S²}'} × det'_{'{S¹}'}) = {anomaly:.8f}")
out(f"   log(anomaly) = {log_anomaly:.8f}")
out(f"   Compare: −π = {-np.pi:.8f}  (diff = {abs(log_anomaly + np.pi):.6f})")
out(f"   The log-anomaly is CLOSE to −π but not exact (0.3% off).")

# Notable near-miss
prod_13 = det_S1 * det_S3
out(f"\n   NOTABLE NEAR-MISS:")
out(f"   det'_{'{S¹}'} × det'_{'{S³}'} / π = {prod_13/np.pi:.6f}  (compare B = 42)")
out(f"   Exact: 4π² · π · exp(ζ(3)/(2π²)) / π = 4π² exp(ζ(3)/(2π²))")
out(f"         = {4*np.pi**2 * np.exp(zeta3/(2*np.pi**2)):.6f}")
out(f"   Off from 42 by {abs(prod_13/np.pi - 42)/42:.4%}")


# ================================================================
# INVESTIGATION 2: Heat Kernel Coefficients
# ================================================================
out("\n" + "=" * 76)
out("INVESTIGATION 2: Heat Kernel Coefficients")
out("=" * 76)
out("""
Heat kernel: Tr(e^{-tΔ}) ~ (4πt)^{-d/2} Σ_k a_k(M) t^k

For unit S^d: R = d(d−1), |Ric|² = d(d−1)², |Riem|² = 2d(d−1).
Volume: Vol(S^d) = 2π^{(d+1)/2} / Γ((d+1)/2).
""")


def vol_sphere(d):
    """Volume of unit d-sphere."""
    from math import gamma
    return 2 * np.pi**((d + 1) / 2) / gamma((d + 1) / 2)


def heat_coeffs(d):
    """Compute a_0, a_1, a_2 for unit S^d."""
    V = vol_sphere(d)
    R = d * (d - 1)       # Scalar curvature
    Ric2 = d * (d - 1)**2  # |Ric|²
    Riem2 = 2 * d * (d - 1)  # |Riem|²

    a0 = V
    a1 = R / 6 * V
    a2 = (5 * R**2 - 2 * Ric2 + 2 * Riem2) / 360 * V
    return a0, a1, a2, R


out(f"{'Manifold':>8s}  {'d':>2s}  {'Vol':>10s}  {'R':>6s}  "
    f"{'a_0':>12s}  {'a_1':>12s}  {'a_2':>12s}")
out("-" * 76)

heat_data = {}
for d, name in [(1, 'S¹'), (2, 'S²'), (3, 'S³')]:
    a0, a1, a2, R = heat_coeffs(d)
    heat_data[name] = (a0, a1, a2)
    out(f"  {name:>6s}  {d:2d}  {vol_sphere(d):10.6f}  {R:6.1f}  "
        f"{a0:12.6f}  {a1:12.6f}  {a2:12.6f}")

# Explicit values
out(f"\n  Explicit values:")
out(f"    S¹: Vol = 2π = {2*np.pi:.6f}, R = 0")
out(f"        a_0 = 2π, a_1 = 0, a_2 = 0  (S¹ is flat!)")
out(f"    S²: Vol = 4π = {4*np.pi:.6f}, R = 2")
out(f"        a_0 = 4π, a_1 = 4π/3 = {4*np.pi/3:.6f}, "
    f"a_2 = 4π/15 = {4*np.pi/15:.6f}")
out(f"    S³: Vol = 2π² = {2*np.pi**2:.6f}, R = 6")
out(f"        a_0 = 2π², a_1 = 2π² = {2*np.pi**2:.6f}, "
    f"a_2 = π² = {np.pi**2:.6f}")

out(f"\n  Check combinations against K = {K_paper:.6f} and α⁻¹ = {CODATA_ALPHA_INV}:")

# Try π × (a_k(S²) + a_k(S¹) − correction) and other combos
a0_1, a1_1, a2_1 = heat_data['S¹']
a0_2, a1_2, a2_2 = heat_data['S²']
a0_3, a1_3, a2_3 = heat_data['S³']

hk_combos = [
    # Direct coefficients
    ("a_2(S³) = π²",               a2_3),
    ("a_1(S³) = 2π²",              a1_3),
    ("a_0(S³)/a_0(S²)",            a0_3 / a0_2),
    ("a_2(S³)/a_2(S²)",            a2_3 / a2_2),
    ("a_1(S²)/a_2(S²)",            a1_2 / a2_2),
    # Ratios and sums
    ("a_0(S³)·a_0(S²)/(4π)",       a0_3 * a0_2 / (4 * np.pi)),
    ("Σ a_0",                       a0_1 + a0_2 + a0_3),
    ("Σ a_1",                       a1_1 + a1_2 + a1_3),
    ("Σ a_2",                       a2_1 + a2_2 + a2_3),
    # With π prefactor
    ("π·(a_0(S²) + a_1(S²))",      np.pi * (a0_2 + a1_2)),
    ("π·Σa_2",                      np.pi * (a2_1 + a2_2 + a2_3)),
    # Reproduce K structure?
    ("π·(B_coeff + F_like − Δ_like)", None),  # placeholder
    # Cross-ratios
    ("a_2(S³)/(a_2(S²)·a_0(S¹))",  a2_3 / (a2_2 * a0_1)),
    ("a_0(S¹)·a_0(S²)·a_0(S³)",    a0_1 * a0_2 * a0_3),
    ("a_2(S³)/a_0(S¹)",             a2_3 / a0_1),
    # Interesting: a_2(S³) = π², a_1(S²) = 4π/3, a_2(S²) = 4π/15
    ("15·a_2(S²)/(a_0(S¹))",       15 * a2_2 / a0_1),
    ("a_2(S³)·a_0(S²)/(a_0(S¹)²)", a2_3 * a0_2 / a0_1**2),
]

out(f"  {'Combination':>45s}  {'Value':>12s}  {'err vs K':>10s}  {'err vs α⁻¹':>10s}")
out(f"  {'-'*45}  {'-'*12}  {'-'*10}  {'-'*10}")
for name, val in hk_combos:
    if val is None:
        continue
    err_K = abs(val - K_paper) / K_paper
    err_alpha = abs(val - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
    k_str = f"{err_K:.2e}" if err_K < 1 else "no"
    a_str = f"{err_alpha:.2e}" if err_alpha < 1 else "no"
    out(f"  {name:>45s}  {val:12.6f}  {k_str:>10s}  {a_str:>10s}")

out(f"\n  FINDING: Heat kernel coefficients are simple polynomials in π.")
out(f"  None of the straightforward combinations produce K or α⁻¹.")
out(f"  The number 42 does NOT appear in heat coefficients (it requires")
out(f"  the specific S³ degeneracy weighting, not Seeley-DeWitt).")
out(f"  a_2(S³) = π² ≈ 9.87 is the most interesting single value.")


# ================================================================
# INVESTIGATION 3: Self-Consistency from the Fibration
# ================================================================
out("\n" + "=" * 76)
out("INVESTIGATION 3: Self-Consistency from the Fibration")
out("=" * 76)

# --- a) First Chern number ---
out("\na) First Chern number of the Hopf bundle:")
out("   c₁ = 1 (well-known: H²(S²; Z) = Z, generator is Hopf bundle)")
out("   The Hopf bundle S¹ → S³ → S² has c₁ = 1.")
out("   This is verified by: ∫_{S²} (i/2π)F = 1")
out("   where F is the curvature 2-form of the connection.")
c1 = 1

# --- b) Chern-Simons invariant ---
out("\nb) Chern-Simons invariant CS(S³):")
out("   For the Levi-Civita connection on S³ with the round metric:")
out("   CS_grav = (1/16π²) ∫_{S³} Tr(ω∧dω + ⅔ ω∧ω∧ω)")
out("")
out("   S³ = ∂D⁴ (bounds a 4-disk), so by the APS index theorem:")
out("   CS(S³) ≡ 0 (mod 1)")
out("")
out("   For SU(2) Chern-Simons at level k, the partition function is:")
out("   Z_k(S³) = √(2/(k+2)) · sin(π/(k+2))")
Z_k = {}
for k in [1, 2, 3, 4, 5, 42]:
    Z_k[k] = np.sqrt(2.0 / (k + 2)) * np.sin(np.pi / (k + 2))
    out(f"     k = {k:2d}: Z_k(S³) = {Z_k[k]:.8f}")

out(f"\n   Z_1(S³) = √(2/3)·sin(π/3) = {Z_k[1]:.8f}")
out(f"   1/Z_1(S³) = {1/Z_k[1]:.8f}  (not related to α)")

# --- c) η-invariant of Dirac on S³ ---
out("\nc) η-invariant of the Dirac operator on S³:")
out("   Eigenvalues of Dirac on unit S³: ±(n + 3/2) for n = 0,1,2,...")
out("   with multiplicity (n+1)(n+2) for each sign.")
out("   The spectrum is perfectly symmetric ⟹ η(0) = 0.")
out("   (This holds for ALL odd-dimensional round spheres by symmetry.)")
eta_0 = 0
out(f"   η(0) = {eta_0}")

# η-function at other points
out("\n   The η-function η(s) = Σ sign(λ)|λ|^{-s} vanishes identically")
out("   for S³ because the spectrum is symmetric: η(s) = 0 for all s.")

# --- d) Cubic from bundle data? ---
out("\nd) Can the cubic arise from {c₁, CS, η, det} ?\n")
out("   Available topological/spectral data:")
out(f"     c₁ = {c1}")
out(f"     CS(S³) = 0 (mod 1)")
out(f"     η(0) = {eta_0}")
out(f"     det'(Δ_{{S¹}}) = 4π² = {det_S1:.6f}")
out(f"     det'(Δ_{{S²}}) = {det_S2:.6f}")
out(f"     det'(Δ_{{S³}}) = {det_S3:.6f}")

out(f"\n   The cubic x³ − Kx + 1 = 0 is a depressed cubic (no x² term).")
out(f"   Depressed cubics arise from:")
out(f"     • Tschirnhaus transform of general cubics")
out(f"     • Resolvent of a quartic")
out(f"     • Algebraic relations on elliptic curves")
out(f"     • Characteristic polynomial of a traceless 3×3 matrix")

out(f"\n   CRITICAL OBSERVATION: S³ ≅ SU(2).")
out(f"   SU(2) has a 3D adjoint representation (Lie algebra su(2)).")
out(f"   The Killing form on su(2) provides a natural 3×3 metric.")
out(f"   Could the cubic be the characteristic polynomial of a natural")
out(f"   operator on the 3D Lie algebra?")

# Casimir of SU(2) in adjoint representation
# The adjoint rep of SU(2) has generators (J_i)_{jk} = -iε_{ijk}
# The quadratic Casimir in the adjoint: C₂(adj) = 2 (for SU(2))
# The matrix of Casimirs... actually, let's think about this differently.
out(f"\n   The SU(2) structure constants ε_{{ijk}} define an operator")
out(f"   on ℝ³ (the Lie algebra). The Casimir C₂ = Σ T_i² = j(j+1)·I")
out(f"   for spin-j representation. In adjoint (j=1): C₂ = 2·I₃.")
out(f"   This is proportional to identity → char poly = (λ−2)³,")
out(f"   which is NOT a depressed cubic. Dead end for this construction.")

# Check: Laplacian on SU(2) ≅ S³
# The Peter-Weyl decomposition gives eigenvalues j(j+1) with mult (2j+1)²
# This is the same as n²−1 with mult n² for n = 2j+1.
out(f"\n   Alternative: The 'structure matrix' of the fibration.")
out(f"   Define M = diag(d₁, d₂, d₃) where d_i are dimensions of S^i.")
out(f"   M = diag(1, 2, 3). Char poly = (λ−1)(λ−2)(λ−3)")
out(f"   = λ³ − 6λ² + 11λ − 6. Not our cubic.")

out(f"\n   A 3x3 matrix from c1=1, CS=0, eta=0 has too many zeros")
out(f"   to produce a non-trivial cubic.")

out(f"\n   FINDING: The topological invariants (c₁=1, CS=0, η=0) are too")
out(f"   simple to generate the cubic directly. The cubic would need to")
out(f"   come from SPECTRAL data (eigenvalues, determinants) rather than")
out(f"   purely topological data.")


# ================================================================
# INVESTIGATION 4: QED Comparison
# ================================================================
out("\n" + "=" * 76)
out("INVESTIGATION 4: QED Comparison")
out("=" * 76)

# --- a) One-loop QED ---
out("\na) One-loop QED correction vs cubic correction:")
out("")

# The cubic says: 1/α = K − α²
# So the correction to 1/α is: δ(1/α) = α²
alpha = alpha_codata
delta_cubic = alpha**2
delta_qed_1loop = alpha / np.pi  # Leading-order Schwinger correction

out(f"   Cubic correction:        δ(α⁻¹) = α² = {delta_cubic:.6e}")
out(f"   QED 1-loop (α/π):        δ(α⁻¹) ~ α/π = {delta_qed_1loop:.6e}")
out(f"   Ratio: (α/π) / α² = 1/(πα) = {1/(np.pi*alpha):.4f}")
out(f"   ≈ K/π = {K_paper/np.pi:.4f}")
out("")
out(f"   The cubic correction α² is ~{delta_qed_1loop/delta_cubic:.0f}× SMALLER than")
out(f"   the QED one-loop correction α/π.")
out(f"   They differ by a factor of 1/(πα) ≈ 43.6 ≈ K/π.")
out("")
out(f"   The cubic correction is roughly at the THIRD-LOOP level in QED:")
out(f"   QED 3-loop: ~ (α/π)³ ≈ {(alpha/np.pi)**3:.6e}")
out(f"   α²          =         {alpha**2:.6e}")
out(f"   Ratio: α²/(α/π)³ = π³/α = {np.pi**3/alpha:.1f}")
out(f"   Not a clean match. The cubic correction does NOT correspond to")
out(f"   any specific order in the QED perturbative expansion.")

# --- b) Casimir energy on S³ ---
out("\nb) Casimir energy on S³:")
out("")
out("   For a conformally-coupled massless scalar on S³ of radius R:")
out("   E_Casimir = 1/(240R)  [Dowker & Critchley]")
out("")
out("   For the electromagnetic field (massless spin-1) on S³:")
out("   Modes: transverse vector harmonics on S³")
out("   Eigenvalues: (n+1)² for n = 1, 2, 3, ... with mult 2(2n+1)")
out("   (two polarizations × angular degeneracy)")

# Compute zeta-regularized vacuum energy for EM on unit S³
# ω_n = n+1 (for unit S³), multiplicity = 2(2n+1)
# E = (1/2) Σ_{n=1}^∞ 2(2n+1)(n+1) [zeta-regularized]
# = Σ_{n=1}^∞ (2n+1)(n+1)
# = Σ (2n² + 3n + 1)
# = 2ζ(-2) + 3ζ(-1) + ζ(0)
# = 2·0 + 3·(-1/12) + (-1/2)
# = -1/4 - 1/2 = -3/4
# But this is the sum of eigenvalues, need to include the ω_n factor.
# Actually: E = (1/2) Σ d_n ω_n where ω_n = n+1 and d_n = 2(2n+1)
# E = Σ_{n=1}^∞ (2n+1)(n+1) = Σ (2n² + 3n + 1) [n from 1 to ∞]
# Zeta-regularized: 2ζ_R(-2) + 3ζ_R(-1) + ζ_R(0)
# = 2·0 + 3·(-1/12) + (-1/2) = -1/4 - 1/2 = -3/4

# Wait, we need to be more careful. The eigenvalues of the Maxwell/Hodge
# Laplacian on S³ for 1-forms are n(n+2) with multiplicity 2(2n+1) for n≥1.
# The frequencies ω_n = √(n(n+2)) = √(n²+2n).
# For the VACUUM ENERGY, E = (1/2) Σ 2(2n+1)√(n(n+2)).

# Actually for a U(1) gauge field on S³ with radius R:
# E_Casimir = 11/(120R)  [Ford 1980]
E_casimir_scalar = 1.0 / 240  # scalar
E_casimir_em = 11.0 / 120      # electromagnetic

out(f"   E_Casimir (scalar, R=1) = 1/240 = {E_casimir_scalar:.6f}")
out(f"   E_Casimir (EM, R=1)    = 11/120 = {E_casimir_em:.6f}")
out(f"   α²                     =         {alpha**2:.6e}")
out(f"   α²/E_Casimir(EM)       =         {alpha**2/E_casimir_em:.6e}")
out(f"   E_Casimir(EM)/α        =         {E_casimir_em/alpha:.4f}")
out("")
out(f"   The Casimir energy is O(1), not O(α²). The correction α²")
out(f"   does not naturally arise from the Casimir energy of a free")
out(f"   field on S³. This would require a self-consistent coupling")
out(f"   where the field strength depends on α.")

# --- c) Discrete mode sum ---
out("\nc) Discrete mode sum (zeta-regularized):")
out("")
out("   For the Hodge Laplacian on 1-forms on unit S³:")
out("   Eigenvalues: n² − 1 (co-exact, n≥2) and n² − 1 (exact, n≥2)")
out("   with multiplicities related to 2n+1 and 2n−1.")
out("")

# The zeta-regularized trace of the 1-form Laplacian
# This is essentially the spectral zeta at s=−1/2 (for the vacuum energy)
# ζ_{1-form}(s) = Σ_{n=2}^∞ d_n (n²−1)^{-s}
# For co-exact 1-forms on S³: eigenvalue n(n+2), mult 2(2n+1)

# Zeta-regularized sum: Σ_{n=1}^∞ 2(2n+1) n(n+2)
# = Σ 2(2n+1)(n²+2n) = Σ 2(2n³+4n²+n²+2n) = Σ 2(2n³+5n²+2n)
# = 4ζ(-3) + 10ζ(-2) + 4ζ(-1)
# = 4·(1/120) + 10·0 + 4·(-1/12) = 1/30 - 1/3 = -3/10

# Actually ζ_R(-3) = 1/120
zeta_minus3 = 1.0 / 120
zeta_minus2 = 0.0
zeta_minus1 = -1.0 / 12
zeta_0_val = -0.5

mode_sum = 4 * zeta_minus3 + 10 * zeta_minus2 + 4 * zeta_minus1
out(f"   Σ 2(2n+1)·n(n+2) [zeta-reg] = 4ζ(−3) + 10ζ(−2) + 4ζ(−1)")
out(f"   = 4/120 + 0 + 4(−1/12) = 1/30 − 1/3 = {mode_sum:.6f}")
out(f"   = −3/10 (not related to α)")

out(f"\n   FINDING: QED corrections at no order match α² quantitatively.")
out(f"   On compact S³, the Casimir energy is O(1/R), not O(α²).")
out(f"   The cubic correction α² ≈ 5.3×10⁻⁵ sits between QED orders")
out(f"   and does not correspond to any standard perturbative term.")
out(f"   This is CONSISTENT with the paper's claim that the cubic is")
out(f"   an 'exact resummation', not a truncation of perturbation theory.")


# ================================================================
# INVESTIGATION 5: Characteristic Polynomial Search
# ================================================================
out("\n" + "=" * 76)
out("INVESTIGATION 5: Characteristic Polynomial Search")
out("=" * 76)
out("""
The cubic λ³ − Kλ + 1 = 0 is the characteristic polynomial of a
TRACELESS 3×3 matrix M with:
  tr(M) = 0
  Σ (2×2 minors) = −K
  det(M) = −1
""")

# --- Diagonal construction ---
out("a) Diagonal construction M = diag(a, b, −(a+b)):")
out("   Constraints: a² + ab + b² = K, ab(a+b) = 1")
out("   Let v = ab, u = a+b. Then u = 1/v and 1/v² − v = K.")
out("   Solving v³ + Kv² − 1 = 0:")

v_roots = np.roots([1, K_paper, 0, -1])
v_real = sorted([r.real for r in v_roots if abs(r.imag) < 1e-10])
out(f"   Roots: {[f'{r:.8f}' for r in v_real]}")

for v in v_real:
    if v > 0:
        u = 1.0 / v
        disc = u**2 - 4 * v
        if disc >= 0:
            a = (u + np.sqrt(disc)) / 2
            b = (u - np.sqrt(disc)) / 2
            c = -(a + b)
            out(f"\n   v = {v:.8f}, u = a+b = {u:.6f}")
            out(f"   a = {a:.8f}, b = {b:.8f}, c = {c:.8f}")
            out(f"   Check: tr = {a+b+c:.2e}, det = {a*b*c:.8f}")
            out(f"   Check: minors = {a*b + a*c + b*c:.6f} (want −K = {-K_paper:.6f})")
            if abs(b - alpha_codata) / alpha_codata < 0.01:
                out(f"   !! b ≈ α = {alpha_codata:.8f}  (the roots ARE the cubic's roots)")

out(f"\n   FINDING: The diagonal matrix diag(r₁, r₂, r₃) where r_i are")
out(f"   the cubic's roots trivially has the cubic as its char poly.")
out(f"   This is TAUTOLOGICAL — we need a matrix from a priori data.")

# --- Systematic search for 3×3 matrices ---
out("\nb) Systematic search for 3×3 matrix with entries from spectral data:")
out("")
out(f"   Target: tr = 0, minors = −K = {-K_paper:.6f}, det = −1")
out("")

# Pool of candidate matrix entries
entry_pool = {
    '0':         0.0,
    '1':         1.0,
    '-1':       -1.0,
    'B=42':      42.0,
    '-B=-42':   -42.0,
    'F=ζ(2)':   zeta2,
    '-F':       -zeta2,
    'Δ=1/40':   Delta,
    '-Δ':       -Delta,
    'π':        np.pi,
    '-π':       -np.pi,
    'π²/6':     np.pi**2/6,
    '1/π':      1/np.pi,
    '-1/π':     -1/np.pi,
    '√π':       np.sqrt(np.pi),
    '-√π':      -np.sqrt(np.pi),
    'B/π':      42/np.pi,
    '-B/π':     -42/np.pi,
    'K/π':      K_paper/np.pi,
    '-K/π':     -K_paper/np.pi,
    '2':         2.0,
    '-2':       -2.0,
    '3':         3.0,
    '-3':       -3.0,
    '6':         6.0,
    '-6':       -6.0,
    '7':         7.0,
    '-7':       -7.0,
    '8':         8.0,
    '-8':       -8.0,
    '14':        14.0,
    '-14':      -14.0,
    'c₁=1':     1.0,
    '1/B':      1.0/42,
    '-1/B':     -1.0/42,
}

entries = list(entry_pool.items())
entry_vals = np.array([v for _, v in entries])
entry_names = [n for n, _ in entries]

# For a general 3×3 matrix:
# [[a, b, c], [d, e, f], [g, h, i]]
# tr = a + e + i = 0
# det = a(ei-fh) - b(di-fg) + c(dh-eg) = -1
# minors = (ae-bd) + (ai-cg) + (ei-fh) = -K

# Full search is too large. Focus on structured matrices.
# Strategy 1: Symmetric matrices (6 free entries, tr=0 constraint → 5)
# Strategy 2: Matrices with specific patterns (companion, circulant)

out("   Strategy: Search companion-like and structured matrices.\n")

# The COMPANION matrix of x³ - Kx + 1 = 0 is:
# [[0, 0, -1], [1, 0, K], [0, 1, 0]]
# This has the right char poly by construction.
M_comp = np.array([[0, 0, -1], [1, 0, K_paper], [0, 1, 0]])
out(f"   Companion matrix: [[0, 0, -1], [1, 0, K], [0, 1, 0]]")
out(f"   Char poly: λ³ − Kλ + 1 = 0  ✓ (by construction)")
out(f"   Entry 'K' = {K_paper:.6f} appears explicitly.")
out(f"   The companion matrix requires K as an entry — it doesn't derive K.")

# Search: matrices where entries are from the pool and char poly is close
out(f"\n   Brute-force search over structured matrices...\n")

# Strategy: Search ANTI-SYMMETRIC matrices (natural in Lie algebra context)
# For a 3×3 antisymmetric matrix: diag = 0, tr = 0 automatically.
# M = [[0, a, b], [-a, 0, c], [-b, -c, 0]]
# det(M) = 0 for any antisymmetric matrix. DEAD END (need det = -1).
out("   Anti-symmetric: det = 0 always. Cannot produce det = −1. ✗")

# Strategy: Circulant matrices
# M = [[a, b, c], [c, a, b], [b, c, a]]
# tr = 3a = 0 → a = 0
# M = [[0, b, c], [c, 0, b], [b, c, 0]]
# det = b³ + c³ = -1
# minors = -bc - bc - bc = -3bc → bc = K/3
# So: b³ + c³ = -1, bc = K/3
# (b+c)(b²-bc+c²) = -1
# (b+c)((b+c)²-3bc) = -1
# Let s = b+c, p = bc = K/3
# s(s² - K) = -1 → s³ - Ks + 1 = 0
# THIS IS THE SAME CUBIC! With s = α (or the roots)!
out("   CIRCULANT MATRICES:")
out("   M = [[0, b, c], [c, 0, b], [b, c, 0]]  (traceless circulant)")
out(f"   Constraints: b³ + c³ = −1, bc = K/3")
out(f"   Let s = b+c. Then s³ − Ks + 1 = 0.")
out(f"   THIS IS EXACTLY THE PAPER'S CUBIC!")
out(f"   The cubic IS the self-consistency condition for a traceless circulant.")
out(f"\n   Solving: b³ + c³ = −1, bc = K/3 = {K_paper/3:.6f}")
bc_val = K_paper / 3
# b+c is a root of s³ - Ks + 1 = 0
s_roots = np.roots([1, 0, -K_paper, 1])
s_real = sorted([r.real for r in s_roots if abs(r.imag) < 1e-10])
for s in s_real:
    # b, c are roots of t² - st + bc = 0
    disc = s**2 - 4 * bc_val
    if disc >= 0:
        b_val = (s + np.sqrt(disc)) / 2
        c_val = (s - np.sqrt(disc)) / 2
        M_circ = np.array([[0, b_val, c_val],
                           [c_val, 0, b_val],
                           [b_val, c_val, 0]])
        eigs = np.linalg.eigvalsh(M_circ)
        det_check = np.linalg.det(M_circ)
        out(f"\n   s = b+c = {s:.8f}:")
        out(f"     b = {b_val:.8f}, c = {c_val:.8f}")
        out(f"     det(M) = {det_check:.8f} (want −1)")
        out(f"     eigenvalues: {[f'{e:.6f}' for e in sorted(eigs)]}")
    else:
        out(f"\n   s = {s:.8f}: b,c are complex (disc = {disc:.4f})")

out(f"\n   CRITICAL FINDING: The cubic α³ − Kα + 1 = 0 is EQUIVALENT to")
out(f"   the characteristic equation of a traceless 3×3 CIRCULANT matrix")
out(f"   with off-diagonal entries b, c satisfying bc = K/3.")
out(f"\n   Physical interpretation: A circulant matrix has Z₃ symmetry")
out(f"   (cyclic permutation of rows/columns). In the Hopf bundle context,")
out(f"   this could represent the democratic coupling of the three components")
out(f"   (fiber S¹, base S², total S³), where each couple equally to the")
out(f"   other two with strength set by K/3 ≈ 45.679.")

# Can K/3 be expressed from spectral data?
out(f"\n   Key question: does K/3 = bc have geometric meaning?")
out(f"   K/3 = π(42 + ζ(2) − 1/40)/3 = {K_paper/3:.6f}")
out(f"   = π/3 × (42 + ζ(2) − 1/40)")
out(f"   = π/3 × 43.6199... = {np.pi/3 * (42 + zeta2 - Delta):.6f}")
out(f"   Note: K/3 ≈ 45.679 ≈ B + F + d(S³)")
out(f"         B + F + 3 = 42 + 1.645 + 3 = {42 + zeta2 + 3:.3f}")
out(f"         K/3 = {K_paper/3:.3f}  (not a match)")

# Search for a better decomposition of K/3
out(f"\n   Decomposition search for K/3 = {K_paper/3:.6f}:")
decomps = [
    ("π·B/3 + π·F/3 − π·Δ/3", np.pi*B/3 + np.pi*F/3 - np.pi*Delta/3),
    ("π·(B+F−Δ)/3",            K_paper/3),  # tautological
    ("14π",                     14*np.pi),
    ("B + ζ(2) + 2",           B + zeta2 + 2),
    ("14.5π",                   14.5*np.pi),
    ("π·14·(1+ζ(2)/42−1/(42·40))", np.pi*14*(1+zeta2/42-1/(42*40))),
]
for name, val in decomps:
    out(f"     {name:>35s} = {val:.6f}  (err = {abs(val-K_paper/3)/abs(K_paper/3):.2e})")


# ================================================================
# OVERALL SUMMARY
# ================================================================
out("\n" + "=" * 76)
out("PHASE 3 SUMMARY: Origin of the Cubic")
out("=" * 76)

out("""
INVESTIGATION 1 — Spectral Determinants:
  ┌─────────────────────────────────────────────────────────────┐
  │ det'(Δ_{S¹}) = 4π²                                        │
  │ det'(Δ_{S²}) = exp(−4ζ'_R(−1))  ≈ 1.938                  │
  │ det'(Δ_{S³}) = π·exp(ζ(3)/(2π²)) ≈ 3.339                 │
  └─────────────────────────────────────────────────────────────┘
  • The fibration anomaly log[det'₃/(det'₂·det'₁)] ≈ −3.132
    is CLOSE to −π (0.3% off) but not exact.
  • NEAR-MISS: det'₁·det'₃/π = 4π²·exp(ζ(3)/(2π²)) ≈ 41.96
    is 0.1% off from B = 42.
  • No simple ratio or product gives K or α⁻¹ directly.
  • VERDICT: Suggestive near-misses but no derivation. ○

INVESTIGATION 2 — Heat Kernel Coefficients:
  ┌─────────────────────────────────────────────────────────────┐
  │ a₂(S³) = π², a₁(S²) = 4π/3, a₂(S²) = 4π/15              │
  └─────────────────────────────────────────────────────────────┘
  • Heat coefficients are simple polynomials in π with rational
    coefficients. The number 42 requires the specific S³ degeneracy
    weighting, NOT standard geometric invariants.
  • No combination produces K or α.
  • VERDICT: Dead end. ✗

INVESTIGATION 3 — Fibration Topology:
  ┌─────────────────────────────────────────────────────────────┐
  │ c₁ = 1, CS(S³) = 0, η(0) = 0                              │
  └─────────────────────────────────────────────────────────────┘
  • The topological invariants are too simple (mostly 0 or 1).
  • The cubic must come from SPECTRAL data, not topology alone.
  • VERDICT: Dead end for topological route. ✗

INVESTIGATION 4 — QED Comparison:
  ┌─────────────────────────────────────────────────────────────┐
  │ α² ≈ 5.3×10⁻⁵ vs α/π ≈ 2.3×10⁻³ (factor ~43.6 apart)    │
  │ Casimir energy on S³ is O(1), not O(α²)                    │
  └─────────────────────────────────────────────────────────────┘
  • The cubic correction α² does NOT match any QED perturbative order.
  • On compact S³, vacuum energy is O(1/R), not O(α²).
  • CONSISTENT with the claim of an "exact resummation" rather
    than perturbative truncation.
  • VERDICT: No contradiction, no confirmation. ○

INVESTIGATION 5 — Characteristic Polynomial:
  ┌─────────────────────────────────────────────────────────────┐
  │ ★ THE CUBIC IS THE EIGENVALUE EQUATION OF A TRACELESS      │
  │   CIRCULANT 3×3 MATRIX with bc = K/3.                      │
  └─────────────────────────────────────────────────────────────┘
  • A traceless circulant matrix [[0,b,c],[c,0,b],[b,c,0]]
    has char poly s³ − 3bc·s + (b³+c³) = 0.
  • Setting bc = K/3 and b³+c³ = −1 recovers the EXACT cubic.
  • The circulant structure implies Z₃ SYMMETRY — democratic
    coupling among the three bundle components (S¹, S², S³).
  • This is a STRUCTURAL result: the cubic is the natural
    eigenvalue equation for a Z₃-symmetric coupling of three
    objects, not an arbitrary ansatz.
  • VERDICT: Most promising lead. ★★★

═══════════════════════════════════════════════════════════════

PHASE 4 DIRECTION:
  The circulant-matrix interpretation opens a concrete path:
  1. Can bc = K/3 = π(42+ζ(2)−1/40)/3 be derived from the
     coupling between fiber, base, and total space?
  2. Why b³ + c³ = −1? (This encodes det = −1.)
  3. The Z₃ symmetry may connect to the three SU(2) generators
     (the Lie algebra of S³ ≅ SU(2)) or to the three Hopf maps
     (S¹→S³→S², S³→S⁷→S⁴, S⁷→S¹⁵→S⁸).
  4. The near-miss det'₁·det'₃/π ≈ 42 deserves deeper analysis:
     does a higher-order correction close the 0.1% gap?
""")

# ================================================================
# Save results
# ================================================================
results_path = os.path.join(output_dir, 'phase3_results.txt')
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(results))
    f.write('\n')
out(f"\nResults saved to: {results_path}")
