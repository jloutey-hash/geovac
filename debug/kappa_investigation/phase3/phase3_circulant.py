"""
Phase 3: Top-down circulant analysis of the cubic alpha^3 - K*alpha + 1 = 0.

The cubic is the characteristic polynomial of the traceless Z3-symmetric circulant:
    M = [[0, b, c],
         [c, 0, b],
         [b, c, 0]]
with bc = K/3 and b^3 + c^3 = -1.

Key question: can b and c be identified with spectral invariants of individual
Hopf bundle components, giving K = 3bc without the additive rule?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import io

# Force UTF-8 output
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

OUTPUT_DIR = Path(__file__).parent

# ==============================================================================
# Constants
# ==============================================================================
pi = np.pi
ALPHA_PHYSICAL = 1 / 137.035999084  # CODATA 2018
B = 42
F = pi**2 / 6
Delta = 1 / 40
K = pi * (B + F - Delta)  # 137.03606...

print("=" * 70)
print("PHASE 3: CIRCULANT STRUCTURE OF THE CUBIC")
print("=" * 70)
print(f"\nK = pi*(B + F - Delta) = {K:.10f}")
print(f"1/alpha_physical       = {1/ALPHA_PHYSICAL:.10f}")
print(f"K vs 1/alpha: diff     = {K - 1/ALPHA_PHYSICAL:.6e}")

# ==============================================================================
# TASK 1: Circulant Algebra -- Solve for b and c
# ==============================================================================
print("\n" + "=" * 70)
print("TASK 1: CIRCULANT ALGEBRA")
print("=" * 70)

p = K / 3  # bc = K/3
print(f"\nbc = K/3 = {p:.10f}")
print(f"b^3 + c^3 = -1")

# The cubic for s = b + c:
# s^3 - 3p*s + 1 = 0  =>  s^3 - K*s + 1 = 0
# This is the SAME cubic as for alpha! (with s in place of alpha)
print(f"\ns = b + c satisfies: s^3 - K*s + 1 = 0")
print("This is the SAME cubic as alpha^3 - K*alpha + 1 = 0!")
print("=> b + c must be one of the three roots of the cubic.")

# Solve the cubic s^3 - K*s + 1 = 0
# np.roots takes coefficients in DESCENDING power order:
# s^3 + 0*s^2 + (-K)*s + 1 = 0
cubic_coeffs = [1, 0, -K, 1]  # s^3 + 0*s^2 - K*s + 1
roots = np.roots(cubic_coeffs)

# All roots should be real (discriminant > 0 for K >> 1)
print(f"\nThree roots of s^3 - K*s + 1 = 0:")
real_roots = []
for i, r in enumerate(roots):
    if np.abs(r.imag) < 1e-10:
        real_roots.append(r.real)
    else:
        print(f"  WARNING: complex root {r}")
        real_roots.append(r.real)

real_roots.sort()
s1, s2, s3 = real_roots
print(f"  s_1 (large negative) = {s1:.15f}")
print(f"  s_2 (small positive) = {s2:.15f}")
print(f"  s_3 (large positive) = {s3:.15f}")

# Verify roots
for i, s in enumerate([s1, s2, s3]):
    resid = s**3 - K*s + 1
    print(f"  Check s_{i+1}: s^3 - K*s + 1 = {resid:.6e}")

print(f"\n  1/s_2 = {1/s2:.10f}  (cf. 1/alpha_phys = {1/ALPHA_PHYSICAL:.10f})")
print(f"  s_2   = {s2:.15e}")
print(f"  alpha = {ALPHA_PHYSICAL:.15e}")
print(f"  s_2 - alpha = {s2 - ALPHA_PHYSICAL:.6e}")

# Vieta relations check
print(f"\nVieta relations:")
print(f"  s1 + s2 + s3 = {s1+s2+s3:.15e}  (should be 0)")
print(f"  s1*s2 + s1*s3 + s2*s3 = {s1*s2+s1*s3+s2*s3:.10f}  (should be -K = {-K:.10f})")
print(f"  s1*s2*s3 = {s1*s2*s3:.10f}  (should be -1)")

# For each root s_i, find b, c from:
#   b + c = s_i
#   b * c = K/3
# => b, c are roots of t^2 - s_i*t + K/3 = 0
# Discriminant = s_i^2 - 4*K/3

print(f"\n--- Checking which root gives real b, c ---")
print(f"Need discriminant s_i^2 - 4K/3 >= 0")
print(f"4K/3 = {4*K/3:.10f}")

real_solutions = []
complex_solutions = []

for i, s in enumerate([s1, s2, s3]):
    disc = s**2 - 4*K/3
    print(f"\n  s_{i+1} = {s:.10f}")
    print(f"    discriminant = {disc:.10f}")
    if disc >= 0:
        b = (s + np.sqrt(disc)) / 2
        c = (s - np.sqrt(disc)) / 2
        print(f"    b = {b:.15f}")
        print(f"    c = {c:.15f}")
        print(f"    VERIFICATION:")
        print(f"      b*c = {b*c:.10f}  (should be {K/3:.10f})")
        print(f"      b^3 + c^3 = {b**3 + c**3:.10f}  (should be -1)")
        print(f"      b + c = {b+c:.10f}  (should be {s:.10f})")
        print(f"      b/c = {b/c:.10f}")
        real_solutions.append((i+1, s, b, c))
    else:
        # Complex b, c (conjugate pair since s is real)
        b = (s + 1j*np.sqrt(-disc)) / 2
        c = (s - 1j*np.sqrt(-disc)) / 2
        print(f"    b = {b.real:.10f} + {b.imag:.10f}i")
        print(f"    c = {c.real:.10f} + {c.imag:.10f}i")
        print(f"    |b| = |c| = {abs(b):.10f}")
        print(f"    arg(b) = {np.angle(b):.10f} rad")
        print(f"    VERIFICATION:")
        print(f"      b*c = {(b*c).real:.10f}  (should be {K/3:.10f})")
        print(f"      b^3 + c^3 = {(b**3 + c**3).real:.10f}  (should be -1)")
        complex_solutions.append((i+1, s, b, c))

# ==============================================================================
# Detailed analysis of real solutions
# ==============================================================================
print("\n" + "=" * 70)
print("DETAILED ANALYSIS OF REAL SOLUTIONS")
print("=" * 70)

for idx, s, b, c in real_solutions:
    print(f"\nSolution from s_{idx} = {s:.10f}:")
    print(f"  b = {b:.15f}")
    print(f"  c = {c:.15f}")
    print(f"  b - c = {b-c:.15f}")
    print(f"  b/c = {b/c:.15f}")
    print(f"  Derived quantities:")
    print(f"    b*c     = {b*c:.10f}  = K/3")
    print(f"    b+c     = {b+c:.10f}  = s_{idx}")
    print(f"    b-c     = {b-c:.10f}")
    print(f"    b^2     = {b**2:.10f}")
    print(f"    c^2     = {c**2:.10f}")
    print(f"    b^2+c^2 = {b**2+c**2:.10f}")
    print(f"    b^3     = {b**3:.10f}")
    print(f"    c^3     = {c**3:.10f}")
    print(f"    b^3+c^3 = {b**3+c**3:.10f}  (should be -1)")

# ==============================================================================
# TASK 1c: Test b and c against spectral invariants
# ==============================================================================
print("\n" + "=" * 70)
print("TASK 1c: SPECTRAL INVARIANT IDENTIFICATION")
print("=" * 70)

# Spectral determinant values from Phase 2
det_prime_S1_times_det_prime_S3_over_pi = 41.957  # near-miss vs B=42

candidates = {
    'sqrt(K/3)': np.sqrt(K/3),
    '-sqrt(K/3)': -np.sqrt(K/3),
    'pi*sqrt(B/3)': pi*np.sqrt(B/3),
    'sqrt(pi*B)': np.sqrt(pi*B),
    'pi*F': pi*F,
    'pi*Delta': pi*Delta,
    'B/pi': B/pi,
    'F*pi': F*pi,
    'pi*B': pi*B,
    'B': float(B),
    'F': F,
    'Delta': Delta,
    'K/pi': K/pi,
    'pi': pi,
    'pi^2': pi**2,
    'sqrt(K)': np.sqrt(K),
    '-1/cbrt(2)': -1/2**(1/3),
    'cbrt(K/3)': (K/3)**(1/3),
    '-cbrt(K/3)': -(K/3)**(1/3),
    'B + F': B + F,
    'B - Delta': B - Delta,
    'F - Delta': F - Delta,
    'B + F - Delta': B + F - Delta,
    'K/3': K/3,
    'sqrt(B)': np.sqrt(B),
    'sqrt(F)': np.sqrt(F),
    'pi*sqrt(B)': pi*np.sqrt(B),
    'B/3': B/3,
    '(B+F)/3': (B+F)/3,
    'pi*B/3': pi*B/3,
    'pi*14': pi*14,
    '14': 14.0,
    'sqrt(pi*B/3)': np.sqrt(pi*B/3),
    'K^(1/3)': K**(1/3),
    'pi*sqrt(14)': pi*np.sqrt(14),
    'sqrt(3*K)': np.sqrt(3*K),
    '2*sqrt(K/3)': 2*np.sqrt(K/3),
    'K/pi^2': K/pi**2,
    'B*pi/K': B*pi/K,
    '3/K': 3/K,
    '1/K': 1/K,
    'sqrt(K)/pi': np.sqrt(K)/pi,
    'pi/sqrt(K)': pi/np.sqrt(K),
    '2*pi': 2*pi,
    '4*pi': 4*pi,
    'exp(1)': np.e,
    'exp(pi)': np.exp(pi),
    'K - B*pi': K - B*pi,
    '(F-Delta)*pi': (F-Delta)*pi,
    'F*pi - Delta*pi': F*pi - Delta*pi,
}

for idx, s, b, c in real_solutions:
    print(f"\nSolution s_{idx}: b = {b:.10f}, c = {c:.10f}")

    for label, val_to_test in [("b", b), ("c", c)]:
        print(f"\n  Testing {label} = {val_to_test:.10f} against candidates:")
        matches = []
        for name, val in candidates.items():
            if abs(val_to_test) > 1e-15:
                rel_err = abs(val_to_test - val) / abs(val_to_test)
            else:
                rel_err = abs(val_to_test - val)
            matches.append((rel_err, name, val))
        matches.sort()
        print(f"    {'Candidate':35s} {'Value':>15s}  {'Rel Error':>12s}")
        print(f"    {'-'*35:35s} {'-'*15:>15s}  {'-'*12:>12s}")
        for err, name, val in matches[:10]:
            marker = " <-- CLOSE" if err < 0.01 else ""
            print(f"    {name:35s} {val:15.10f}  {err:12.6e}{marker}")

# Also test for complex solutions
for idx, s, b, c in complex_solutions:
    print(f"\nComplex solution s_{idx}: |b| = |c| = {abs(b):.10f}")
    print(f"  Testing |b| against candidates:")
    abs_b = abs(b)
    matches = []
    for name, val in candidates.items():
        rel_err = abs(abs_b - abs(val)) / abs_b
        matches.append((rel_err, name, val))
    matches.sort()
    for err, name, val in matches[:5]:
        marker = " <-- MATCH" if err < 1e-6 else ""
        print(f"    {name:35s} |val|={abs(val):15.10f}  rel_err = {err:.6e}{marker}")

# ==============================================================================
# TASK 2: THE SELF-REFERENTIAL STRUCTURE
# ==============================================================================
print("\n" + "=" * 70)
print("TASK 2: SELF-REFERENTIAL STRUCTURE")
print("=" * 70)

print(f"\nThe cubic s^3 - K*s + 1 = 0 has three roots:")
print(f"  s_1 = {s1:.15f}  (large negative)")
print(f"  s_2 = {s2:.15f}  (small positive ~= 1/K ~= alpha)")
print(f"  s_3 = {s3:.15f}  (large positive)")

print(f"\nFor REAL b, c, we need s^2 >= 4K/3 = {4*K/3:.6f}")
print(f"  s_1^2 = {s1**2:.6f}  {'>=':4s} {4*K/3:.6f}  -> {'YES' if s1**2 >= 4*K/3 else 'NO'}")
print(f"  s_2^2 = {s2**2:.10f}  {'>=':4s} {4*K/3:.6f}  -> {'YES' if s2**2 >= 4*K/3 else 'NO'}")
print(f"  s_3^2 = {s3**2:.6f}  {'>=':4s} {4*K/3:.6f}  -> {'YES' if s3**2 >= 4*K/3 else 'NO'}")

# 2a: Identify which root
print(f"\n2a: b + c = alpha (s_2) requires COMPLEX b,c")
print(f"    b + c = s_1 (large negative, ~-11.7) gives REAL b,c")
print(f"    b + c = s_3 (large positive, ~+11.7) gives REAL b,c")

# 2b: Algebraic property — k=0 Fourier mode
print(f"\n2b: The k=0 Fourier mode of circulant [[0,b,c],[c,0,b],[b,c,0]] is")
print(f"    lambda_0 = 0 + b*1 + c*1 = b + c.")
print(f"    So b+c is ALWAYS an eigenvalue of M. This is structural.")
print(f"    The 'self-referential' property (b+c satisfies same cubic as alpha)")
print(f"    is simply the statement that b+c IS one of the eigenvalues of M,")
print(f"    and ALL eigenvalues of M satisfy the same characteristic polynomial.")

# 2c: Physical meaning
print(f"\n2c: The two real solutions:")
for idx, s, b, c in real_solutions:
    print(f"    s_{idx} = {s:.10f} => b = {b:.10f}, c = {c:.10f}")

print(f"\n    The large roots are ~ +/- sqrt(K) = +/- {np.sqrt(K):.6f}")
print(f"    Making b and c individually O(1) to O(10)")

# What are the eigenvalues when b+c = s1 (large negative)?
print(f"\n    CRITICAL: NO real solutions exist!")
print(f"    s_1^2 = {s1**2:.6f} < 4K/3 = {4*K/3:.6f}")
print(f"    s_3^2 = {s3**2:.6f} < 4K/3 = {4*K/3:.6f}")
print(f"    All three roots give COMPLEX CONJUGATE b, c pairs.")
print(f"    |b| = |c| = sqrt(K/3) = {np.sqrt(K/3):.10f} for ALL three solutions.")

# ==============================================================================
# TASK 3: EIGENVALUES AND FOURIER DIAGONALIZATION
# ==============================================================================
print("\n" + "=" * 70)
print("TASK 3: EIGENVALUES AND FOURIER DIAGONALIZATION")
print("=" * 70)

omega = np.exp(2j*pi/3)

print("\n  Since b,c are always complex conjugates, the circulant M is complex.")
print("  But the characteristic polynomial still has the three real roots s_1, s_2, s_3.")
print("")
print("  For each solution (b+c = s_i, b*c = K/3):")

for idx, s, b, c in complex_solutions:
    print(f"\n--- Solution s_{idx}: b+c = {s:.10f} ---")
    print(f"  b = {b.real:.10f} + {b.imag:.10f}i")
    print(f"  c = {c.real:.10f} + {c.imag:.10f}i")

    # Build complex M
    M = np.array([[0, b, c],
                   [c, 0, b],
                   [b, c, 0]])

    # Eigenvalues
    eigs_M = np.linalg.eigvals(M)
    eigs_sorted = sorted(eigs_M, key=lambda x: x.real)
    print(f"  Eigenvalues of M:")
    for j, e in enumerate(eigs_sorted):
        print(f"    eig_{j+1} = {e.real:+.10f} + {e.imag:+.6e}i")

    # Check against cubic roots
    print(f"  Cubic roots: {s1:.10f}, {s2:.10f}, {s3:.10f}")
    print(f"  Match check:")
    for j, e in enumerate(eigs_sorted):
        best_j = min(range(3), key=lambda k: abs(e - [s1, s2, s3][k]))
        best_s = [s1, s2, s3][best_j]
        print(f"    eig = {e.real:+.10f}{e.imag:+.6e}i  closest to s_{best_j+1} = {best_s:.10f}  |diff| = {abs(e - best_s):.6e}")

    # Circulant formula eigenvalues
    print(f"\n  Circulant formula: lambda_k = b*omega^k + c*omega^(2k)")
    for k in range(3):
        lam = b * omega**k + c * omega**(2*k)
        print(f"    k={k}: {lam.real:+.10f} + {lam.imag:+.6e}i")

print("\n3c: KEY STRUCTURAL RESULT")
print("  b and c are ALWAYS complex conjugates: b = r*exp(i*theta), c = r*exp(-i*theta)")
print(f"  with r = sqrt(K/3) = {np.sqrt(K/3):.10f}")
print(f"")
print(f"  The three solutions correspond to three phases theta_1, theta_2, theta_3")
print(f"  related by the constraint cos(3*theta) = -1/(2*r^3):")
r_val = np.sqrt(K/3)
cos3theta_val = -1 / (2 * r_val**3)
print(f"  cos(3*theta) = -1/(2*(K/3)^(3/2)) = {cos3theta_val:.15f}")
base_angle = np.arccos(cos3theta_val)
print(f"  3*theta_0 = arccos({cos3theta_val:.10f}) = {base_angle:.10f}")
print(f"  Three solutions: theta = theta_0/3, (theta_0 + 2*pi)/3, (theta_0 + 4*pi)/3")
for k in range(3):
    theta_k = (base_angle + 2*pi*k) / 3
    print(f"    theta_{k} = {theta_k:.10f} rad = {np.degrees(theta_k):.4f} deg")
    b_k = r_val * np.exp(1j*theta_k)
    c_k = r_val * np.exp(-1j*theta_k)
    s_k = (b_k + c_k).real
    print(f"      b+c = 2*r*cos(theta) = {s_k:.10f}")

# Also the negative arccos branch gives 3 more, but cos(3theta) has period 2pi/3
# so we get the same 3 values
print(f"\n  Matching to cubic roots:")
for k in range(3):
    theta_k = (base_angle + 2*pi*k) / 3
    s_k = 2 * r_val * np.cos(theta_k)
    best_j = min(range(3), key=lambda j: abs(s_k - [s1, s2, s3][j]))
    best_s = [s1, s2, s3][best_j]
    print(f"    theta_{k}: b+c = {s_k:.10f} matches s_{best_j+1} = {best_s:.10f} (diff {abs(s_k-best_s):.2e})")

print(f"\n  INTERPRETATION:")
print(f"  All three roots of the cubic are b+c = 2*sqrt(K/3)*cos(theta)")
print(f"  where theta cycles through three values separated by ~2*pi/3.")
print(f"  The 'alpha' root s_2 corresponds to theta ~ pi/2 (nearly pure imaginary b).")
print(f"  The large roots s_1, s_3 correspond to theta ~ pi + delta, delta")
print(f"  (nearly real b, with opposite signs).")

# ==============================================================================
# TASK 4: REPRESENTATION-THEORETIC CONTENT
# ==============================================================================
print("\n" + "=" * 70)
print("TASK 4: REPRESENTATION THEORY")
print("=" * 70)

# 4a: Z3 action
print("\n4a: Z3 symmetry")
print("  The Z3 is algebraic, not geometric (S1, S2, S3 have different dims).")
print("  It permutes ROLES (fiber, base, total) cyclically.")
print("  Possible origins:")
print("    - Three SU(2) generators (Pauli matrices)")
print("    - Three Hopf maps: S3->S2, S7->S4, S15->S8")
print("    - Triality of SO(8) (which governs the S7 Hopf fibration)")
print("    - The three embeddings S1 x S1 x S1 -> S3 (Clifford tori)")

# 4b: Deformation from SU(2) isotropic case
print("\n4b: Deformation from SU(2) structure constants")
print("  SU(2): b = c = 1")
print("  Characteristic polynomial: s^3 - 3s - 2 = 0 = (s-2)(s+1)^2")
su2_roots = np.roots([1, 0, -3, -2])
su2_roots_real = np.sort(np.real(su2_roots))
print(f"  Roots: {su2_roots_real}  (degenerate: -1, -1, 2)")

print(f"\n  Physical: b != c, K = {K:.6f}")
print(f"  SU(2): b = c = 1, K_SU2 = 3")
print(f"  'Coupling enhancement': K/K_SU2 = {K/3:.6f}")

# 4c: Isotropy analysis
print("\n4c: Isotropy analysis")
print(f"  Isotropic case: b = c")
print(f"    bc = b^2 = K/3  =>  b = sqrt(K/3) = {np.sqrt(K/3):.10f}")
print(f"    b^3 + c^3 = 2b^3 = -1  =>  b = -1/cbrt(2) = {-1/2**(1/3):.10f}")
print(f"  These are INCOMPATIBLE. The isotropic point does NOT exist.")
print(f"  The constraint surface bc = K/3 and b^3 + c^3 = -1 forces anisotropy.")

print(f"\n  Since b = r*exp(i*theta), c = r*exp(-i*theta):")
print(f"    r = sqrt(K/3) = {np.sqrt(K/3):.10f}")
print(f"    For the alpha solution (s_2): theta ~ pi/2")
print(f"    b/c = exp(2i*theta) ~ exp(i*pi) = -1  (nearly anti-phase)")
print(f"    For SU(2): b = c = 1 (real), so theta = 0, r = 1")
print(f"    Deformation: r: 1 -> {np.sqrt(K/3):.4f},  theta: 0 -> pi/2")

# ==============================================================================
# ADDITIONAL: Explore the complex solutions more carefully
# ==============================================================================
print("\n" + "=" * 70)
print("ADDITIONAL: COMPLEX SOLUTIONS (b+c = alpha)")
print("=" * 70)

for idx, s, b, c in complex_solutions:
    if abs(s - s2) < 1e-6:
        print(f"\n  s_{idx} = s_2 = alpha ~ {s:.15f}")
        print(f"  b = {b.real:.10f} + {b.imag:.10f}i")
        print(f"  c = {c.real:.10f} - {c.imag:.10f}i  (conjugate)")
        print(f"  |b| = |c| = {abs(b):.10f}")
        print(f"  arg(b) = {np.angle(b):.10f} rad = {np.degrees(np.angle(b)):.6f} deg")
        print(f"")
        print(f"  Note: |b| = |c| = sqrt(K/3) = {np.sqrt(K/3):.10f}")
        print(f"  Because |b|^2 = b*conj(b) = b*c = K/3 (since c = conj(b))")
        print(f"  So the complex solution has |b| = |c| = sqrt(K/3) EXACTLY.")
        print(f"  The modulus encodes the coupling strength; the phase encodes")
        print(f"  the asymmetry that makes b^3+c^3 = -1 (determinant condition).")
        print(f"")
        # Phase analysis
        theta = np.angle(b)
        print(f"  Phase theta = {theta:.10f}")
        print(f"  pi/2 = {pi/2:.10f}")
        print(f"  theta - pi/2 = {theta - pi/2:.10f}")
        # b = r*exp(i*theta), c = r*exp(-i*theta), r = sqrt(K/3)
        # b^3 + c^3 = 2*r^3*cos(3*theta) = -1
        # cos(3*theta) = -1/(2*r^3)
        r_complex = abs(b)
        cos3theta = -1 / (2 * r_complex**3)
        # arccos gives the principal value in [0, pi], but 3*theta may be
        # in a different branch. Check all three: arccos, 2*pi - arccos, arccos + 2*pi
        base_angle = np.arccos(cos3theta)
        print(f"  From det condition: cos(3*theta) = -1/(2*r^3) = {cos3theta:.15f}")
        print(f"  arccos({cos3theta:.10f}) = {base_angle:.10f}")
        # Three solutions for theta in [0, 2*pi):
        # 3*theta = base_angle, base_angle + 2*pi, base_angle + 4*pi
        # OR 3*theta = 2*pi - base_angle (negative branch)
        # More correctly: cos(3*theta) = cos3theta has solutions
        #   3*theta = base_angle + 2*k*pi  or  3*theta = -base_angle + 2*k*pi
        print(f"  Possible theta values:")
        found_match = False
        for k in range(-3, 4):
            for sign in [1, -1]:
                t = (sign * base_angle + 2*pi*k) / 3
                if abs(np.cos(3*t) - cos3theta) < 1e-10:
                    match = abs(t - theta) < 1e-6
                    if match:
                        found_match = True
                    marker = " <-- MATCH" if match else ""
                    if abs(t - theta) < 0.5:  # only show nearby solutions
                        print(f"    theta = {t:.10f}  (k={k}, sign={sign:+d})  cos(3t) = {np.cos(3*t):.10f}{marker}")
        if found_match:
            print(f"  Theta verification: PASSED")
        print(f"")
        print(f"  INTERPRETATION: In the complex (b+c = alpha) solution,")
        print(f"  the circulant parameters are b = sqrt(K/3) * exp(i*theta)")
        print(f"  where theta is fixed by the determinant condition.")
        print(f"  The 'coupling strength' K/3 sets the modulus.")
        print(f"  The 'chirality' (det M = -1) sets the phase.")

# Also check: what is s_2 exactly?
print(f"\n  Precision check on s_2:")
print(f"  s_2 = {s2:.15e}")
print(f"  1/K = {1/K:.15e}")
print(f"  s_2 - 1/K = {s2 - 1/K:.6e}")
print(f"  s_2 ~ 1/K + 1/K^3 = {1/K + 1/K**3:.15e}")
print(f"  (From series: if s is small, s^3 ~ 0, so s ~ 1/K)")

# ==============================================================================
# PLOTS
# ==============================================================================
print("\n" + "=" * 70)
print("GENERATING PLOTS")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: The cubic and its roots
s_range = np.linspace(-15, 15, 1000)
y = s_range**3 - K*s_range + 1

ax = axes[0, 0]
ax.plot(s_range, y, 'b-', linewidth=2)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)
for sr in [s1, s2, s3]:
    ax.plot(sr, 0, 'ro', markersize=8)
    if abs(sr) > 1:
        ax.annotate(f'{sr:.4f}', (sr, 0), textcoords="offset points",
                    xytext=(0, 15), ha='center', fontsize=9)
    else:
        ax.annotate(f'{sr:.6f}\n(~alpha)', (sr, 0), textcoords="offset points",
                    xytext=(40, 25), ha='center', fontsize=9,
                    arrowprops=dict(arrowstyle='->', color='red'))
ax.set_xlabel('s')
ax.set_ylabel('f(s) = s^3 - Ks + 1')
ax.set_title('Cubic equation: three real roots')
ax.set_ylim(-500, 500)
ax.grid(True, alpha=0.3)

# Plot 2: Three complex solutions in polar form
ax = axes[0, 1]
r_val = np.sqrt(K/3)
# Circle of radius r
theta_circle = np.linspace(0, 2*pi, 100)
ax.plot(r_val * np.cos(theta_circle), r_val * np.sin(theta_circle), 'k--', alpha=0.2)

colors_sol = ['red', 'green', 'blue']
for (idx_c, s, b, c), col in zip(complex_solutions, colors_sol):
    ax.plot(b.real, b.imag, 'o', color=col, markersize=10, label=f'b (s_{idx_c}={s:.4f})')
    ax.plot(c.real, c.imag, 's', color=col, markersize=8, alpha=0.5)
    ax.plot([0, b.real], [0, b.imag], '-', color=col, alpha=0.3)
    ax.plot([0, c.real], [0, c.imag], '--', color=col, alpha=0.3)

ax.set_xlabel('Re(b)')
ax.set_ylabel('Im(b)')
ax.set_title(f'All three (b,c) pairs\n|b|=|c|=sqrt(K/3)={r_val:.3f}')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

# Plot 3: Root evolution with K
ax = axes[1, 0]
K_range = np.linspace(3.5, 150, 500)
roots_vs_K = []
for Ki in K_range:
    ri = np.roots([1, 0, -Ki, 1])
    ri_real = np.sort(np.real(ri))
    roots_vs_K.append(ri_real)
roots_vs_K = np.array(roots_vs_K)

ax.plot(K_range, roots_vs_K[:, 0], 'r-', label='s_1 (large negative)', linewidth=2)
ax.plot(K_range, roots_vs_K[:, 1], 'g-', label='s_2 ~ 1/K (alpha)', linewidth=2)
ax.plot(K_range, roots_vs_K[:, 2], 'b-', label='s_3 (large positive)', linewidth=2)
ax.axvline(x=K, color='gray', linestyle='--', alpha=0.5, label=f'K = {K:.2f}')
ax.axvline(x=3, color='orange', linestyle='--', alpha=0.5, label='K = 3 (SU(2))')
ax.set_xlabel('K (coupling strength)')
ax.set_ylabel('Roots of s^3 - Ks + 1 = 0')
ax.set_title('Root evolution: SU(2) (K=3) to Physical (K~137)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 4: Complex solution phase diagram
ax = axes[1, 1]
# Show b, c in complex plane for complex solutions
theta_range = np.linspace(0, 2*pi, 100)
r_val = np.sqrt(K/3)
circle_x = r_val * np.cos(theta_range)
circle_y = r_val * np.sin(theta_range)
ax.plot(circle_x, circle_y, 'b--', alpha=0.3, label=f'|z| = sqrt(K/3) = {r_val:.3f}')

for idx_c, s, b, c in complex_solutions:
    if abs(s - s2) < 1e-6:
        ax.plot(b.real, b.imag, 'ro', markersize=10, label=f'b (s={s:.6f})')
        ax.plot(c.real, c.imag, 'bs', markersize=10, label=f'c = conj(b)')
        # Draw angle
        ax.plot([0, b.real], [0, b.imag], 'r-', alpha=0.5)
        ax.plot([0, c.real], [0, c.imag], 'b-', alpha=0.5)
        theta_b = np.angle(b)
        ax.annotate(f'theta = {theta_b:.4f} rad\n= {np.degrees(theta_b):.2f} deg',
                   (b.real, b.imag), textcoords="offset points",
                   xytext=(10, 10), fontsize=9)

ax.set_xlabel('Re')
ax.set_ylabel('Im')
ax.set_title('Complex solution: b, c in complex plane\n(b+c = alpha, |b|=|c|=sqrt(K/3))')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'circulant_eigenvalues.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'circulant_eigenvalues.png'}")

# Plot 5: Eigenvalue structure
fig2, ax = plt.subplots(figsize=(10, 4))
# Number line showing all three roots
for i, (sr, col, lab) in enumerate(zip([s1, s2, s3],
                                        ['red', 'green', 'blue'],
                                        [f's_1={s1:.4f}', f's_2={s2:.6f} (alpha)', f's_3={s3:.4f}'])):
    ax.plot(sr, 0, 'o', color=col, markersize=14, zorder=5)
    ax.annotate(lab, (sr, 0), textcoords="offset points",
               xytext=(0, 20), ha='center', fontsize=10,
               arrowprops=dict(arrowstyle='->', color=col))

# Mark sqrt(K) for reference
ax.axvline(x=np.sqrt(K), color='gray', linestyle=':', alpha=0.5)
ax.axvline(x=-np.sqrt(K), color='gray', linestyle=':', alpha=0.5)
ax.text(np.sqrt(K), 0.3, f'sqrt(K)={np.sqrt(K):.2f}', ha='center', fontsize=8, color='gray')
ax.text(-np.sqrt(K), 0.3, f'-sqrt(K)', ha='center', fontsize=8, color='gray')

ax.axhline(y=0, color='k', linewidth=0.5)
ax.set_xlabel('Value')
ax.set_title('Three eigenvalues of the circulant M = Three roots of the cubic')
ax.set_yticks([])
ax.set_ylim(-0.5, 0.8)
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'eigenvalue_number_line.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'eigenvalue_number_line.png'}")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
print("\n" + "=" * 70)
print("SUMMARY OF KEY FINDINGS")
print("=" * 70)

print(f"""
1. CIRCULANT ALGEBRA:
   K = {K:.10f}, bc = K/3 = {K/3:.10f}
   b + c satisfies the SAME cubic as alpha (self-referential).

   Three roots: s_1 = {s1:.10f}, s_2 = {s2:.10f} (~alpha), s_3 = {s3:.10f}

   ALL THREE roots give COMPLEX b,c (s^2 < 4K/3 = {4*K/3:.2f} for all):
     s_1^2 = {s1**2:.2f}, s_2^2 ~ 0, s_3^2 = {s3**2:.2f}  (all < {4*K/3:.2f})

   b and c are ALWAYS complex conjugates: b = r*exp(i*theta), c = r*exp(-i*theta)
   with r = sqrt(K/3) = {np.sqrt(K/3):.10f} for ALL solutions.
   b + c = 2*r*cos(theta) gives the three roots via three theta values.""")

r_val = np.sqrt(K/3)
cos3t = -1 / (2 * r_val**3)
base_a = np.arccos(cos3t)
for k in range(3):
    tk = (base_a + 2*pi*k) / 3
    sk = 2*r_val*np.cos(tk)
    best = min(range(3), key=lambda j: abs(sk - [s1,s2,s3][j]))
    print(f"     theta_{k} = {tk:.6f} rad -> b+c = {sk:.10f} = s_{best+1}")

print(f"""
2. SELF-REFERENTIAL STRUCTURE:
   The k=0 Fourier mode of any circulant [[0,b,c],[c,0,b],[b,c,0]] is b+c.
   So b+c is ALWAYS an eigenvalue. The self-reference is structural, not
   a coincidence: ALL eigenvalues satisfy the characteristic polynomial.

   Vieta: s_1 + s_2 + s_3 = 0
          s_1 * s_2 * s_3 = -1
          s_1*s_2 + s_1*s_3 + s_2*s_3 = -K

   So s_1 * s_3 ~ -1/alpha (the large roots carry the 'inverse alpha' scale).

3. SPECTRAL INVARIANT IDENTIFICATION:
   No exact match found for b or c individually against spectral invariants.
   The off-diagonal parameters are ALGEBRAICALLY determined by K through
   the cubic, not independently identifiable with B, F, or Delta.

   However, the COMPLEX solution (b+c = alpha) is elegant:
     |b| = |c| = sqrt(K/3)  [coupling strength]
     arg(b) fixed by det(M) = -1  [chirality]
   This separates K into modulus (physics) and phase (topology).

4. ISOTROPY:
   The isotropic point b = c does NOT exist on the constraint surface
   (sqrt(K/3) != -1/cbrt(2)). b and c are always complex conjugates.
   For the alpha solution, theta ~ pi/2 means b ~ +i*r, c ~ -i*r (anti-phase).
   The 'anisotropy' is encoded as a phase rotation, not a magnitude difference.""")

print(f"""
5. KEY STRUCTURAL INSIGHT:
   The circulant explains WHY the self-consistency equation is a depressed
   cubic with unit constant term (tracelessness + Z3 + det=-1).
   It does NOT independently determine K.

   The path forward is not finding b,c from spectral invariants (they're
   algebraically determined by K), but understanding WHY the Z3 circulant
   structure arises from the Hopf bundle. The det(M) = -1 condition
   (fermionic sign / Hopf linking number) is the most physically motivated
   constraint; the bc = K/3 condition encodes the coupling strength.

   The COMPLEX solution b = sqrt(K/3)*exp(i*theta) with cos(3*theta) = -1/(2*(K/3)^(3/2))
   may be the most natural: it says the coupling lives on a circle of
   radius sqrt(K/3) in the complex plane, with the phase fixed by chirality.
""")

plt.close('all')
print("Phase 3 computation complete.")
