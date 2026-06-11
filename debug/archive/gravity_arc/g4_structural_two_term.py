"""Structural investigation: WHY is the S^3 Dirac spectral action two-term exact?

The mechanism: zeta_unit(-k) = 0 for all k >= 0. This happens because the
Camporesi-Higuchi spectrum on S^3 has:
  eigenvalues |lambda_n| = n + 3/2 (half-integer shifted linear)
  degeneracies g_n = 2(n+1)(n+2) (quadratic polynomial in n)

Paper 51 G3 established: scalar Laplacian on S^3 does NOT have two-term
exactness (a_k^Delta = 2*pi^2/k! for scalar, non-zero). So it's spinor-specific.

Questions:
1. Is the identity zeta(-k) = 0 sensitive to the SHIFT (3/2)?
2. Is it sensitive to the DEGENERACY pattern g_n = 2(n+1)(n+2)?
3. What minimal conditions produce the cancellation?
4. Are there OTHER spectra (not S^3 Dirac) that also have this property?
"""

import numpy as np
from fractions import Fraction
import sympy as sp


def spectral_zeta_neg_k(eigenvalues, degeneracies, k: int, N_max: int = 500):
    """Compute zeta(-k) = sum g_n * |lambda_n|^{2k} for a spectral triple.

    At s = -k, |lambda|^{-2s} = |lambda|^{2k}, so this is a polynomial sum.
    """
    total = 0.0
    for n in range(N_max):
        if n < len(eigenvalues):
            lam = eigenvalues[n]
            g = degeneracies[n]
        else:
            break
        total += g * abs(lam) ** (2 * k)
    return total


def zeta_s3_dirac(k: int, N_max: int = 200):
    """zeta_unit(-k) for CH Dirac on S^3: |lam_n| = n+3/2, g_n = 2(n+1)(n+2)."""
    # Use exact rational arithmetic
    total = Fraction(0)
    for n in range(N_max):
        lam = Fraction(2*n + 3, 2)  # n + 3/2
        g = 2 * (n + 1) * (n + 2)
        total += g * lam ** (2 * k)
    return total


def zeta_s3_scalar(k: int, N_max: int = 200):
    """zeta(-k) for scalar Laplacian on S^3: |lam_n|^2 = n(n+2), g_n = (n+1)^2.

    Note: scalar uses SQUARED eigenvalues n(n+2) = (n+1)^2 - 1,
    so |lam_n|^{2k} = (n(n+2))^k.
    """
    total = Fraction(0)
    for n in range(1, N_max):  # start at n=1 (n=0 is zero mode)
        lam_sq = n * (n + 2)
        g = (n + 1) ** 2
        total += g * lam_sq ** k
    return total


def zeta_general_linear(k: int, shift, deg_coeffs, N_max: int = 200):
    """zeta(-k) for a general spectrum |lam_n| = n + shift, g_n = poly(n).

    deg_coeffs: coefficients of degeneracy polynomial g_n = c0 + c1*n + c2*n^2 + ...
    """
    total = Fraction(0)
    shift = Fraction(shift) if not isinstance(shift, Fraction) else shift
    for n in range(N_max):
        lam = n + shift
        g = sum(Fraction(c) * n**i for i, c in enumerate(deg_coeffs))
        total += g * lam ** (2 * k)
    return total


def main():
    print("=" * 72)
    print("STRUCTURAL: Why is the S^3 Dirac spectral action two-term exact?")
    print("=" * 72)

    # Test 1: Verify zeta_unit(-k) = 0 for S^3 Dirac
    print("\n  === Test 1: zeta_unit(-k) for S^3 Dirac ===")
    print(f"  Spectrum: |lam_n| = n + 3/2, g_n = 2(n+1)(n+2)")
    print(f"  {'k':>3}  {'zeta(-k)':>20}  {'zero?':>6}")
    for k in range(8):
        z = zeta_s3_dirac(k)
        print(f"  {k:3d}  {float(z):20.6f}  {'YES' if z == 0 else 'NO'}")

    # Test 2: S^3 scalar (should NOT be zero)
    print("\n  === Test 2: zeta(-k) for S^3 scalar Laplacian ===")
    print(f"  Spectrum: lam_sq = n(n+2), g_n = (n+1)^2")
    print(f"  {'k':>3}  {'zeta(-k)':>20}")
    for k in range(6):
        z = zeta_s3_scalar(k)
        print(f"  {k:3d}  {float(z):20.6f}")

    # Test 3: What if we change the SHIFT from 3/2?
    print("\n  === Test 3: Sensitivity to the spectral shift ===")
    print(f"  Keep g_n = 2(n+1)(n+2), vary shift in |lam_n| = n + shift")
    print(f"  {'shift':>8}  {'z(-0)':>12}  {'z(-1)':>12}  {'z(-2)':>12}  "
          f"{'z(-3)':>12}  {'all_zero':>10}")
    for shift_num, shift_den in [(1,2), (1,1), (3,2), (2,1), (5,2), (3,1)]:
        shift = Fraction(shift_num, shift_den)
        deg_coeffs = [Fraction(4), Fraction(6), Fraction(2)]  # 2(n+1)(n+2) = 2n^2+6n+4
        zs = [zeta_general_linear(k, shift, deg_coeffs) for k in range(4)]
        all_zero = all(z == 0 for z in zs)
        print(f"  {float(shift):8.3f}  {float(zs[0]):12.4f}  {float(zs[1]):12.4f}  "
              f"{float(zs[2]):12.4f}  {float(zs[3]):12.4f}  "
              f"{'YES' if all_zero else 'NO'}")

    # Test 4: What if we change the DEGENERACY?
    print("\n  === Test 4: Sensitivity to the degeneracy pattern ===")
    print(f"  Keep shift=3/2, vary g_n")
    print(f"  {'g_n pattern':>25}  {'z(-0)':>12}  {'z(-1)':>12}  {'z(-2)':>12}  "
          f"{'all_zero':>10}")

    degeneracy_tests = [
        ("2(n+1)(n+2) [CH S^3]", [Fraction(4), Fraction(6), Fraction(2)]),
        ("(n+1)^2", [Fraction(1), Fraction(2), Fraction(1)]),
        ("n+1 [linear]", [Fraction(1), Fraction(1)]),
        ("(n+1)(n+2)(n+3)/3", [Fraction(2), Fraction(4), Fraction(8,3), Fraction(2,3)]),
        ("4(n+1) [S^2 Dirac]", [Fraction(4), Fraction(4)]),
        ("1 [constant]", [Fraction(1)]),
    ]
    for label, coeffs in degeneracy_tests:
        shift = Fraction(3, 2)
        zs = [zeta_general_linear(k, shift, coeffs) for k in range(4)]
        all_zero = all(z == 0 for z in zs)
        print(f"  {label:>25}  {float(zs[0]):12.4f}  {float(zs[1]):12.4f}  "
              f"{float(zs[2]):12.4f}  {'YES' if all_zero else 'NO'}")

    # Test 5: The algebraic identity
    # g_n = 2(n+1)(n+2) = 2[(n+3/2)^2 - 1/4] when shift = 3/2
    # So zeta(s) = 2*sum (n+3/2)^{-2s} * [(n+3/2)^2 - 1/4]
    #            = 2*[zeta_H(2s-2, 3/2) - (1/4)*zeta_H(2s, 3/2)]
    # At s=-k: 2*[zeta_H(-2k-2, 3/2) - (1/4)*zeta_H(-2k, 3/2)]
    # = 2*[-B_{2k+3}(3/2)/(2k+3) + (1/4)*B_{2k+1}(3/2)/(2k+1)]
    print("\n  === Test 5: The algebraic mechanism ===")
    print(f"  g_n = 2(n+1)(n+2) = 2[(n+3/2)^2 - 1/4]")
    print(f"  zeta(s) = 2[zeta_H(2s-2, 3/2) - (1/4)zeta_H(2s, 3/2)]")
    print(f"  At s=-k: need -2B_{{2k+3}}(3/2)/(2k+3) + (1/2)B_{{2k+1}}(3/2)/(2k+1) = 0")
    print()
    x = sp.Rational(3, 2)
    print(f"  {'k':>3}  {'B_{2k+1}(3/2)':>16}  {'B_{2k+3}(3/2)':>16}  "
          f"{'term1':>12}  {'term2':>12}  {'sum':>12}")
    for k in range(6):
        B_odd_1 = sp.bernoulli(2*k+1, x)
        B_odd_3 = sp.bernoulli(2*k+3, x)
        term1 = -2 * B_odd_3 / (2*k+3)
        term2 = sp.Rational(1, 2) * B_odd_1 / (2*k+1)
        total = term1 + term2
        print(f"  {k:3d}  {str(B_odd_1):>16}  {str(B_odd_3):>16}  "
              f"{str(term1):>12}  {str(term2):>12}  {str(total):>12}")

    # Test 6: Can we characterize ALL (shift, degeneracy) pairs that give
    # two-term exactness?
    print("\n  === Test 6: Search for other two-term-exact spectra ===")
    print(f"  Testing |lam_n| = n + shift, g_n = a + b*n + c*n^2")
    print(f"  Searching shifts in [1/4, 7/2] and degree-2 degeneracies")
    print(f"  {'shift':>8}  {'(a,b,c)':>20}  {'z(-0)':>10}  {'z(-1)':>10}  "
          f"{'z(-2)':>10}  {'exact?':>6}")

    found = []
    for shift_4 in range(1, 15):  # shift = shift_4/4
        shift = Fraction(shift_4, 4)
        # For each shift, find if there's a (a,b,c) that makes all zeta(-k)=0
        # Set up linear system: zeta(-k) = sum_n (a+b*n+c*n^2)*(n+shift)^{2k} = 0
        # This is linear in (a,b,c). Three equations (k=0,1,2) determine (a,b,c).
        # Check if k=3,4,... also vanish.

        # Build coefficient matrix for k=0,1,2
        N = 100
        A_mat = np.zeros((3, 3))
        for k in range(3):
            for n in range(N):
                lam_2k = float(n + shift) ** (2*k)
                A_mat[k, 0] += lam_2k  # coefficient of a
                A_mat[k, 1] += n * lam_2k  # coefficient of b
                A_mat[k, 2] += n**2 * lam_2k  # coefficient of c

        # The system A_mat @ [a,b,c]^T = 0 has nontrivial solutions
        # when rank < 3. Check singular values.
        U, S, Vt = np.linalg.svd(A_mat)
        if S[-1] / S[0] < 1e-10:
            # Rank deficient — nontrivial null space
            null_vec = Vt[-1]
            # Normalize so c=1 (quadratic leading coeff)
            if abs(null_vec[2]) > 1e-10:
                null_vec = null_vec / null_vec[2]
                a_val, b_val, c_val = null_vec
                # Verify at k=3,4,5
                zs = []
                for k in range(6):
                    z = sum((a_val + b_val*n + c_val*n**2) *
                            float(n + shift)**(2*k) for n in range(N))
                    zs.append(z)
                all_small = all(abs(z) < 1e-6 for z in zs)
                if all_small and c_val > 0:  # positive leading = physical
                    found.append((float(shift), a_val, b_val, c_val))
                    print(f"  {float(shift):8.4f}  ({a_val:.3f},{b_val:.3f},"
                          f"{c_val:.3f})  {zs[0]:10.2e}  {zs[1]:10.2e}  "
                          f"{zs[2]:10.2e}  {'YES' if all_small else 'NO'}")

    if not found:
        print("  No additional two-term-exact spectra found in the search range.")
    else:
        print(f"\n  Found {len(found)} two-term-exact spectra:")
        for shift, a, b, c in found:
            print(f"    shift={shift:.4f}, g_n = {a:.4f} + {b:.4f}*n + {c:.4f}*n^2")
            # Compare to CH: shift=1.5, g = 4+6n+2n^2 = 2(n+1)(n+2)
            if abs(shift - 1.5) < 0.01:
                print(f"    ^ THIS IS S^3 CH Dirac (up to normalization)")

    print("\n" + "=" * 72)
    print("SYNTHESIS")
    print("=" * 72)


if __name__ == "__main__":
    main()
