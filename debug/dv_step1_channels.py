"""Sprint 5 DV Step 1: Enumerate all (k_1, k_2) channels contributing to
the He (1s)(2p) ^3P_J direct and exchange matrix elements of C^(K)/r_12^{K+1}.

For each channel, compute the angular coefficient (9j × reduced ME) and
print it symbolically. The angular coefficients combine with the (channel-
specific) radial integrals to give the full matrix element.
"""
from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan


def red_C(l, lp, k):
    """Reduced matrix element <l || C^(k) || l'> (Racah normalization)."""
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * wigner_3j(l, k, lp, 0, 0, 0)


def enumerate_channels(K, la, lb, lc, ld, L=1, Lp=1, k_max=4):
    """List all (k1, k2) with non-zero contribution to
    <(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(K) || (l_c l_d) L'>.
    """
    result = []
    for k1 in range(0, k_max + 1):
        for k2 in range(0, k_max + 1):
            # Triangle (k1, k2, K)
            if abs(k1 - k2) > K or k1 + k2 < K:
                continue
            # Parity: <l || C^(k) || l'> requires k + l + l' even
            if (k1 + la + lc) % 2 != 0:
                continue
            if (k2 + lb + ld) % 2 != 0:
                continue
            # Triangle (l_a, k1, l_c) and (l_b, k2, l_d)
            if abs(la - lc) > k1 or la + lc < k1:
                continue
            if abs(lb - ld) > k2 or lb + ld < k2:
                continue
            red1 = red_C(la, lc, k1)
            red2 = red_C(lb, ld, k2)
            n9 = wigner_9j(la, lb, L, lc, ld, Lp, k1, k2, K)
            # Coupling of k1, k2 to K along 0: <k1 0 k2 0 | K 0>
            # (projective selection rule; zero unless k1+k2+K is even, already enforced)
            cg = clebsch_gordan(k1, k2, K, 0, 0, 0)

            if simplify(red1) == 0 or simplify(red2) == 0:
                continue
            if simplify(n9) == 0:
                continue
            # Angular coefficient (without radial):
            # Spatial reduced ME  (l_a l_b) L || [C^(k1)(1) x C^(k2)(2)]^(K) || (l_c l_d) L' contains:
            #   sqrt((2L+1)(2L'+1)(2K+1)) * 9j * red1 * red2
            # The bipolar expansion coefficient of Y^K/r_12^{K+1} in [C^(k1)(1) x C^(k2)(2)]^(K) contains:
            #   sqrt((2k1+1)(2k2+1)) <k1 0 k2 0 | K 0> / sqrt(2K+1)   (from Sack/Brink-Satchler)
            # Combining:
            spatial_me = (sqrt(Integer((2 * L + 1) * (2 * Lp + 1) * (2 * K + 1))) *
                          n9 * red1 * red2)
            # Bipolar-decomposition angular coefficient (without radial kernel):
            bipolar_ang = sqrt(Integer((2 * k1 + 1) * (2 * k2 + 1))) * cg / sqrt(Integer(2 * K + 1))

            ang_total = simplify(spatial_me * bipolar_ang)
            result.append((k1, k2, simplify(cg), simplify(red1), simplify(red2),
                           simplify(n9), ang_total))
    return result


print("=" * 76)
print("Sprint 5 DV Step 1: Channel enumeration")
print("=" * 76)

for K, label in [(2, "SS"), (1, "SOO"), (0, "monopole")]:
    print(f"\n\n=== K = {K} ({label}) ===")

    print(f"\n--- DIRECT channel: (l_a, l_b, l_c, l_d) = (0, 1, 0, 1), L=L'=1 ---")
    for ch in enumerate_channels(K, 0, 1, 0, 1):
        k1, k2, cg, r1, r2, n9, ang = ch
        print(f"  (k1={k1}, k2={k2}): <.|C^{k1}|.>={r1}, <.|C^{k2}|.>={r2}, "
              f"9j={n9}, <k1 0 k2 0|K 0>={cg}")
        print(f"    angular coefficient (with bipolar prefactor / sqrt(2K+1)): {ang}")

    print(f"\n--- EXCHANGE channel: (l_a, l_b, l_c, l_d) = (0, 1, 1, 0), L=L'=1 ---")
    for ch in enumerate_channels(K, 0, 1, 1, 0):
        k1, k2, cg, r1, r2, n9, ang = ch
        print(f"  (k1={k1}, k2={k2}): <.|C^{k1}|.>={r1}, <.|C^{k2}|.>={r2}, "
              f"9j={n9}, <k1 0 k2 0|K 0>={cg}")
        print(f"    angular coefficient: {ang}")
