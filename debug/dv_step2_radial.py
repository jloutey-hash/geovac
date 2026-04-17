"""Sprint 5 DV Step 2: Numerically evaluate the bipolar kernel radial
integrals for the (k_1, k_2) channels contributing to He (1s)(2p) ^3P.

Goal: verify whether the channel-specific bipolar kernel integrals
     I(k1, k2, K; a, b; c, d)
evaluate to the SAME radial integral as Drake's M^K_{dir/exch} (up to
known angular factors), or to DIFFERENT integrals.

The bipolar expansion of Y^K_Q(hat r_12) / r_12^{K+1} (from
Varshalovich §5.17 or Brink-Satchler App. 5 / Sack 1964) is a
sum of terms with channel-specific (k_1, k_2) radial kernels. For
the "leading" bipolar expansion where k_1 + k_2 = K:

    Y^K_Q(hat r_12) / r_12^{K+1}
        = sum_{k1+k2=K}  beta(k1, k2, K) * [Y^{k1}(1) (x) Y^{k2}(2)]^K_Q
          * r_1^{k1} r_2^{k2} / r_>^{2K+1}  [leading term]

with beta(k1, k2, K) a specific rational coefficient. But this is only
the leading part; the full expansion has HIGHER-(k_1, k_2) channels
with different radial kernels (see Brink-Satchler 4.6.13).

For the Breit-Pauli radial kernel 1/r_12^{K+3} (which is what appears in
Drake's M^K), we need to account for the r_12^2 in the denominator
beyond Y^K/r_12^{K+1}.

Actually, the Breit-Pauli SS operator is:
    H_SS ∝ [(s_1 · hat r_12)(s_2 · hat r_12) - (1/3) s_1 · s_2] / r_12^3
          ∝ [Y^(2)(hat r_12) · ...] / r_12^3

So it involves Y^(K=2)(hat r_12) / r_12^{K+1=3}.  Match!

The bipolar expansion of Y^(K)(hat r_12) / r_12^{K+1} is the object we
need. Reference formulae:

--- Brink-Satchler "Angular Momentum" App. 5 / Eq. 4.6.13 ---

Y^K_Q(hat r_12) / r_12^{K+1}
  = (1/4pi) sum_{k1,k2} (2K+1)! / [(2k1+1)!! (2k2+1)!! (2K)!]^{1/2}
    * ... complex expression

This is unwieldy. Let me instead use a DIFFERENT strategy: directly
COMPUTE the needed radial integrals from scratch for each (k1, k2)
channel and see how they relate to M^K_dir and M^K_exch.

For the channel (k1, k2, K) contributing to the DIRECT path:
    I_dir(k1, k2, K) = int_0^inf int_0^inf R_{1s}(r_1)^2 R_{2p}(r_2)^2
                         * g_{k1 k2 K}(r_1, r_2) * r_1^2 r_2^2 dr_1 dr_2
where g_{k1 k2 K} is the bipolar radial kernel.

STANDARD FORM (Rose "Elementary Theory of Angular Momentum" 1957 §7):

The bipolar expansion of a general f(r_12) is
    f(r_12) = sum_L F_L(r_1, r_2) P_L(cos theta_12)

For f(r_12) = 1/r_12^n, the F_L are given explicitly in Sack 1964.
For f(r_12) = Y^K(hat r_12) / r_12^{K+1}, the bipolar form is:

   Y^K_M(hat r_12) / r_12^{K+1}
     = (4pi / (2K+1))^{1/2} * sum_{k1,k2: triangle(k1,k2,K)}
         <k1, k2; K || K> * r_<^K / r_>^{K+1}
         * sum_{q1, q2} <k1 q1 k2 q2 | K M> Y^{k1}_{q1}(1) Y^{k2}_{q2}(2)

where <k1, k2; K || K> is a specific angular coefficient and r_<^K/r_>^{K+1}
is the UNIVERSAL radial kernel (!). But this form is only valid for
CERTAIN (k1, k2) — not all.

The specific identity we need is:

BRINK-SATCHLER 4.6.13 (for the Coulomb kernel Y^0/r_12 = 1/(4pi r_12)):
   1/r_12 = 4pi sum_{l,m} (1/(2l+1)) r_<^l / r_>^{l+1} Y^l_m(1)* Y^l_m(2)

Generalizing to Y^K(hat r_12) / r_12^{K+1}:

   Y^K_M(hat r_12) / r_12^{K+1}
     = (4pi / (2K+1)) sum_{l1,l2: triangle(l1, l2, K)} (-1)^{?}
         * sqrt((2l1+1)(2l2+1)/((2K+1)(4pi)^2))
         * <l1 0 l2 0 | K 0>
         * r_<^? / r_>^?
         * [Y^(l1)(hat r_1) (x) Y^(l2)(hat r_2)]^(K)_M

(I'm losing myself in convention differences. Let me approach this
differently — brute force numeric check.)

BRUTE FORCE CHECK:

I'll evaluate the MATRIX ELEMENT of C^(2)_0(hat r_12) / r_12^3 between
coupled (1s, 2p) ^3P states directly via numerical angular integration
AND via the bipolar-channel decomposition, and see which matches.
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
# Force UTF-8 for stdout
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

from geovac.breit_integrals import breit_ss_radial

# ============================================================
# Drake's M^K integrals (at Z=1, symbolic closed form)
# ============================================================
M2_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1)
M2_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1)
M1_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 1, Z=1)
M1_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 1, Z=1)

print("Drake M^K integrals at Z=1:")
print(f"  M^2_dir  = {simplify(M2_dir)}")
print(f"  M^2_dir  ≈ {float(M2_dir):.6e}")
print(f"  M^2_exch = {simplify(M2_exch)}")
print(f"  M^2_exch ≈ {float(M2_exch):.6e}")
print(f"  M^1_dir  = {simplify(M1_dir)}")
print(f"  M^1_dir  ≈ {float(M1_dir):.6e}")
print(f"  M^1_exch = {simplify(M1_exch)}")
print(f"  M^1_exch ≈ {float(M1_exch):.6e}")

# ============================================================
# Drake's A_SS at Z=1
# ============================================================
A_SS_drake = Rational(3, 50) * M2_dir - Rational(2, 5) * M2_exch
A_SOO_drake = Rational(3, 2) * M1_dir - M1_exch

print(f"\nA_SS (Drake, Z=1, no alpha^2): {simplify(A_SS_drake)}")
print(f"  Float: {float(A_SS_drake):.6e}")
print(f"A_SOO (Drake, Z=1): {simplify(A_SOO_drake)}")
print(f"  Float: {float(A_SOO_drake):.6e}")

# ============================================================
# Now: compute <^3P_1, M_J=1 | H_SS | ^3P_1, M_J=1> from scratch
# using the correct Racah form with 9j recoupling.
#
# The matrix element at J=1 corresponds to A_SS (f_SS(1)=1).
#
# Wigner-Eckart on scalar tensor (Edmonds 7.1.6):
#
# <(L S) J M | [T^(K)(space) . U^(K)(spin)]^(0) | (L S) J M>
#   = (-1)^{L+S+J} 6j{L S J; S L K} / sqrt(2J+1)
#       * <L||T||L> * <S||U||S>
#
# where the ".^(0)" scalar product divides by sqrt(2K+1), and the
# (-1)^K factor might also appear depending on convention.
#
# Drake's A_SS = <M_J=1| H_SS |M_J=1> at L=S=1, J=1.
# With 6j{1,1,1;1,1,2} = -1/6 and (-1)^{1+1+1} = -1:
#
#   (-1) * (-1/6) / sqrt(3) = 1/(6 sqrt 3) ~ 0.0962
#
# So A_SS = (1/(6 sqrt 3)) * <L=1||T^(2)||L=1> * <S=1||U^(2)||S=1>
#
# We know <S=1||U^(2)||S=1> = sqrt(5)/2.
#
# H_SS = -alpha^2 sqrt(24 pi / 5) [s1 x s2]^(2) . Y^(2)(hat r_12)/r_12^3
#      = -alpha^2 sqrt(6) [s1 x s2]^(2) . C^(2)(hat r_12)/r_12^3
# (using Y^K = sqrt((2K+1)/(4pi)) C^K, so sqrt(24 pi/5) * sqrt(5/(4pi)) = sqrt(6))
#
# So T^(2) = -sqrt(6) C^(2)(hat r_12)/r_12^3 * alpha^2 (where alpha^2 is pulled out)
#
# <L=1||T^(2)||L=1> (spatial) = -sqrt(6) alpha^2 * <L=1 || C^(2)(hat r_12)/r_12^3 || L=1>
#
# <L=1 || C^(2)(hat r_12)/r_12^3 || L=1>  needs bipolar decomposition.
#
# PART A: Direct channel (L=L'=1, (l_a, l_b)=(0,1)=(l_c, l_d))
# PART B: Exchange (L=L'=1, (l_a, l_b)=(0,1), (l_c, l_d)=(1,0))
#
# Total = (1/sqrt 2)^2 * (Direct - Exchange)  [for triplet spatial antisymm.]
# But the factor is actually  (1/2) * (A_dir - A_exch - A_exch + A_dir) if
# we unpack the antisymmetrized product — it's easier to just handle the
# (direct - exchange) form naturally.

# ============================================================
# STEP: Compute spatial reduced matrix element channel-by-channel
# ============================================================
#
# <(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(K) || (l_c l_d) L'>
#   = sqrt((2L+1)(2L'+1)(2K+1)) 9j{l_a l_b L; l_c l_d L'; k1 k2 K}
#       * <l_a||C^(k1)||l_c> <l_b||C^(k2)||l_d>
#
# and the bipolar expansion of C^(K)(hat r_12)/r_12^{K+1} into
# [C^(k1)(1) x C^(k2)(2)]^(K) form is:
#
# C^(K)_Q(hat r_12) / r_12^{K+1}
#   = sum_{k1,k2: triangle(k1,k2,K), k1+k2+K even}
#       B(k1, k2, K; r_1, r_2) [C^(k1)(1) x C^(k2)(2)]^(K)_Q
#
# with B(k1, k2, K; r_1, r_2) the channel radial kernel.
#
# THE SACK THEOREM (Sack 1964 J. Math. Phys. 5, 245):
# Expansion of r_12^{-nu} (for nu > 0) in bipolar form:
#
# r_12^{-nu} = sum_{n, l} A(nu, n, l; r_1, r_2) P_l(cos theta_12)
#
# with A(nu, n, l) specific radial functions (Gegenbauer expansion).
#
# For r_12^{-nu} P_K(cos theta_12) = r_12^{-nu} * (4pi/(2K+1)) Y^K_0(hat r_12)* Y^K_0(hat ?) -- wait
#
# The KEY FORMULA WE NEED (Rose 1957, or Edmonds):
#
# Y^K_M(hat r_12) / r_12^{K+1} factorizes as:
#
# Y^K_M(hat r_12) / r_12^{K+1}
#   = sqrt(4pi/(2K+1)) * sum_{l1, l2} <l1 0 l2 0 | K 0> * sqrt((2l1+1)(2l2+1)/(4pi))
#       * g_{l1 l2 K}(r_1, r_2) * [Y^(l1)(hat r_1) (x) Y^(l2)(hat r_2)]^(K)_M
#
# where g_{l1 l2 K}(r_1, r_2) depends on (l1, l2, K) AND on the relative
# size of r_1 and r_2. For the LEADING Sack term l1+l2 = K:
#
# g_{l1 l2 K}(r_1, r_2) = r_1^{l1} r_2^{l2} / r_>^{2K+1}  (if r_1 < r_2 means r_> = r_2)
#                        ...symmetrically if r_2 < r_1.
#
# But this is NOT the full story — there are additional terms at l1+l2 = K+2, K+4, ...
# HOWEVER, these higher terms involve k1+k2+K > K*2, which for our (l_a, l_b, l_c, l_d)
# selection rules cannot contribute to <(0,1) L || | | (0,1) L> or <(0,1) L || | | (1,0) L>.
#
# So for OUR specific problem, ONLY the leading Sack term (k1+k2=K) contributes:
#
# DIRECT (K=2): (k1=0, k2=2). Kernel: r_1^0 r_2^2 / r_>^5 = r_2^2 / r_>^5.
#   When r_1 < r_2: r_> = r_2. Kernel = r_2^2 / r_2^5 = 1/r_2^3.
#   When r_1 > r_2: r_> = r_1. Kernel = r_2^2 / r_1^5.
#
# EXCHANGE (K=2): (k1=1, k2=1). Kernel: r_1 r_2 / r_>^5.
#   r_1 < r_2: r_1 r_2 / r_2^5 = r_1 / r_2^4.
#   r_1 > r_2: r_1 r_2 / r_1^5 = r_2 / r_1^4.
#
# Compare to Drake's M^K = integral r_<^K / r_>^{K+3}:
#   M^K direct/exchange = integral of r_<^2 / r_>^5.
#   r_1 < r_2: r_<^2/r_>^5 = r_1^2 / r_2^5.
#   r_1 > r_2: r_<^2/r_>^5 = r_2^2 / r_1^5.

# These are DIFFERENT radial kernels. The bipolar (k1=0, k2=2) kernel is
# r_2^2/r_>^5, while Drake's M^K is r_<^2/r_>^5.
#
# Let's compute each channel's radial integral numerically and express
# it in terms of M^K.

# ============================================================
# STEP: Numerically compute channel-specific radial integrals
# ============================================================
import numpy as np
from scipy import integrate

# Hydrogenic radial functions at Z=1:
# R_{n,l}(r) = N_{n,l} r^l e^{-r/n} * L^{(2l+1)}_{n-l-1}(2r/n)
#
# Normalized 1s: R_{1,0}(r) = 2 e^{-r}
# Normalized 2p: R_{2,1}(r) = (1/(2 sqrt(6))) r e^{-r/2}


def R_1s(r):
    return 2 * np.exp(-r)


def R_2p(r):
    return (1 / (2 * np.sqrt(6))) * r * np.exp(-r / 2)


def radial_integral(f_r1, f_r2, g_kernel, r_max=30.0):
    """Compute ∫∫ f_r1(r1) f_r2(r2) g_kernel(r1, r2) r1^2 r2^2 dr1 dr2."""
    def integrand(r1, r2):
        return f_r1(r1) * f_r2(r2) * g_kernel(r1, r2) * r1**2 * r2**2
    val, err = integrate.dblquad(integrand, 0, r_max, 0, r_max, epsabs=1e-10, epsrel=1e-8)
    return val, err


# Drake M^2_dir: direct = <1s, 2p | ... | 1s, 2p>, kernel = r_<^2 / r_>^5
# Here "direct" means orbital on electron 1 is 1s (unchanged) and on e2 is 2p (unchanged).
# So f_r1(r1) = R_1s(r1)^2, f_r2(r2) = R_2p(r2)^2, kernel = r_<^2 / r_>^5
def kernel_rlt_rgt(K_lt, K_gt):
    """Return kernel function r_<^K_lt / r_>^K_gt."""
    def g(r1, r2):
        if r1 < r2:
            return r1**K_lt / r2**K_gt
        else:
            return r2**K_lt / r1**K_gt
    return g


def kernel_explicit(p1, p2, q1, q2):
    """Kernel r_1^p1 r_2^p2 / r_>^{q} for r_1 < r_2 (= r_2^{q1}),
    and r_1^p1 r_2^p2 / r_>^q for r_1 > r_2 (= r_1^{q2}).

    Actually: kernel = r_1^p1 r_2^p2 / r_>^Q where Q = 2K+1."""
    def g(r1, r2):
        if r1 < r2:
            return r1**p1 * r2**p2 / r2**(q1)
        else:
            return r1**p1 * r2**p2 / r1**(q2)
    return g


# Drake's M^2_dir: kernel r_<^2 / r_>^5
val_M2d_num, _ = radial_integral(
    lambda r: R_1s(r) ** 2,
    lambda r: R_2p(r) ** 2,
    kernel_rlt_rgt(2, 5),
)
print(f"\nM^2_dir (numerical): {val_M2d_num:.6e}")
print(f"M^2_dir (symbolic):  {float(M2_dir):.6e}")

# Drake's M^2_exch: direct = <1s, 2p | ... | 2p, 1s> exchange (orbital swap)
# kernel r_<^2 / r_>^5 but with <1s|.|2p> on r_1 and <2p|.|1s> on r_2
# = R_1s(r_1) R_2p(r_1) R_2p(r_2) R_1s(r_2)
val_M2e_num, _ = radial_integral(
    lambda r: R_1s(r) * R_2p(r),
    lambda r: R_2p(r) * R_1s(r),
    kernel_rlt_rgt(2, 5),
)
print(f"M^2_exch (numerical): {val_M2e_num:.6e}")
print(f"M^2_exch (symbolic):  {float(M2_exch):.6e}")

# Bipolar (k1=0, k2=2, K=2) kernel for DIRECT:
# kernel = r_1^0 r_2^2 / r_>^(2K+1) = r_2^2 / r_>^5
# r_1 < r_2: r_> = r_2, kernel = r_2^2 / r_2^5 = 1/r_2^3
# r_1 > r_2: r_> = r_1, kernel = r_2^2 / r_1^5
def kernel_bipolar_direct_k1k2_K(k1, k2, K):
    def g(r1, r2):
        r_gt = max(r1, r2)
        return r1**k1 * r2**k2 / r_gt**(2 * K + 1)
    return g


val_Idir_02_2, _ = radial_integral(
    lambda r: R_1s(r) ** 2,
    lambda r: R_2p(r) ** 2,
    kernel_bipolar_direct_k1k2_K(0, 2, 2),
)
print(f"\nBipolar (k1=0, k2=2, K=2) DIRECT radial I: {val_Idir_02_2:.6e}")
print(f"  Kernel: r_2^2 / r_>^5")
print(f"  Ratio to M^2_dir: {val_Idir_02_2 / val_M2d_num:.6f}")

# Bipolar (k1=1, k2=1, K=2) kernel for EXCHANGE:
# kernel = r_1 r_2 / r_>^5
# We also need to be careful which density this multiplies. For EXCHANGE:
# < 1s(1) 2p(2) | kernel * [angular] | 2p(1) 1s(2) >
# radial part = R_1s(r_1) R_2p(r_1) R_2p(r_2) R_1s(r_2) * kernel
val_Iexch_11_2, _ = radial_integral(
    lambda r: R_1s(r) * R_2p(r),
    lambda r: R_2p(r) * R_1s(r),
    kernel_bipolar_direct_k1k2_K(1, 1, 2),
)
print(f"\nBipolar (k1=1, k2=1, K=2) EXCHANGE radial I: {val_Iexch_11_2:.6e}")
print(f"  Kernel: r_1 r_2 / r_>^5")
print(f"  Ratio to M^2_exch: {val_Iexch_11_2 / val_M2e_num:.6f}")

# For K=1 SOO: checked empty spatial bipolar; skip.
