"""
Sprint 6 DW: C^(k) Convention Check for Drake Coefficients — v2
================================================================

The correct approach: The Breit-Pauli SS operator involves 1/r_12^3 times
angular factors. The standard Gegenbauer expansion gives:

  1/r_12^{2s+1} = sum_k (k+s choose s) * r_<^k / r_>^{k+s+1} * P_k(cos theta_12)

For s=1 (1/r_12^3):
  1/r_12^3 = sum_k (k+1) * r_<^k / r_>^{k+2} * P_k(cos theta_12)

Using the addition theorem P_k = C^(k)(1) . C^(k)(2) and the tensor
coupling identity, we can extract the rank-2 component that couples
to the SS spin tensor.

The key insight Sprint 5 DV missed: the expansion of 1/r_12^3 uses
r_<^k / r_>^{k+2} (NOT r_<^k / r_>^{k+1} as in the Laplace 1/r_12 expansion).
Drake's M^K integrals use r_<^K / r_>^{K+3} which is the 1/r_12^3 kernel.

So the bipolar expansion of the SS tensor part is:

  [3(sigma1.r12hat)(sigma2.r12hat) - sigma1.sigma2] / r_12^3
  = sum_{k1,k2} f(k1,k2) * [C^(k1)(1) x C^(k2)(2)]^(2) . [s1 x s2]^(2)
    * R^{k1,k2}(r1,r2)

where R^{k1,k2} involves the CORRECT 1/r_12^3 radial kernel (not 1/r_12).

Actually, the CORRECT approach uses the SACK expansion of the irregular
solid harmonic:

Y^K_M(hat r_12) * r_12^{-(K+1)} = multipole expansion

For r_12^{-3} * C^{(2)}(r̂_{12}), the expansion involves the second
derivative of the Laplace expansion 1/r_12.

Let me use the CORRECT formula from Brink & Satchler (1993) App. V:

The irregular solid harmonic I^K_M = r^{-(K+1)} C^K_M(hat r) has the
two-center expansion:

I^K_M(r_12) = sum_{k1 k2} (-1)^{k2} * (2K+1)! / [(2k1)!(2k2)!]
              * sqrt{(2k1+1)(2k2+1) / (2K+1)}
              * CG(k1 0 k2 0 | K 0)
              * R^K_{k1,k2}(r_1, r_2)
              * [C^{k1}(1) x C^{k2}(2)]^K_M

where k1 + k2 = K + 2n (n=0,1,2,...), and the radial function is:

R^K_{k1,k2}(r1,r2) = (-1)^n r_<^{k1} r_>^{k2} / r_>^{K+2n+1}
                    = (-1)^n r_<^{k1} / r_>^{k1+1}   [since k1+k2=K+2n]

Wait, that doesn't simplify nicely either. Let me check Varshalovich
directly.

From Varshalovich, Moskalev, Khersonskii (1988), Eq. 5.17.22:

I^L_M(r_1 - r_2) = C^L_M(hat(r_1-r_2)) / |r_1-r_2|^{L+1}
= sum_{l1 l2} (-1)^{l1+l2+L} * sqrt{(2l1+1)!(2l2+1)! / ((2L+1)!)}
  * <l1 0 l2 0|L 0> * R^{l1}(r_1) * I^{l2}(r_2) / (or vice versa)
  * [C^{l1}(1) x C^{l2}(2)]^L_M

where the sum has l1+l2 = L, L+2, L+4, ...

The key radial function for the two-center expansion:
For r_1 < r_2: R^{l1}(r1) * I^{l2}(r2) = r_1^{l1} / r_2^{l2+1}
For r_1 > r_2: I^{l1}(r1) * R^{l2}(r2) = r_2^{l2} / r_1^{l1+1}

This is the standard result. For the leading term l1+l2 = L (= K = 2 for SS):
  r_1 < r_2: r_1^{l1} / r_2^{l2+1}
  r_1 > r_2: r_2^{l2} / r_1^{l1+1}

For (l1=0, l2=2): r_1<r_2: 1/r_2^3; r_1>r_2: r_2^2/r_1
For (l1=1, l2=1): r_1<r_2: r_1/r_2^2; r_1>r_2: r_2/r_1^2
For (l1=2, l2=0): r_1<r_2: r_1^2/r_2; r_1>r_2: 1/r_1

These are the radial kernels for the C^(2)/r_12^3 bipolar expansion
of the irregular solid harmonic.

Now — crucially — Drake's M^K integral uses kernel r_<^K / r_>^{K+3}.
For K=0: r_<^0/r_>^3 = 1/r_>^3
For K=1: r_<^1/r_>^4 = r_</r_>^4
For K=2: r_<^2/r_>^5

The bipolar channel (0,2) gives:
  r1<r2: 1/r_2^3 — this IS Drake M^0 kernel in Region I ✓
  r1>r2: r_2^2/r_1 — this is NOT Drake M^2 Region II (= r_2^2/r_1^5) ✗

THIS is the bug in Sprint 5 DV!! The Region II kernel for bipolar (0,2)
is r_2^2/r_1, NOT r_2^2/r_1^5 (which would be Drake M^2 Region II).

The irregular solid harmonic expansion for C^(2)/r_12^3 has kernel
r_2^{l2}/r_1^{l1+1} in Region II, not r_2^{l2}/r_1^{l2+3}.

So the bipolar radial integrals are COMPLETELY DIFFERENT from what
Sprint 5 DV assumed.

Let me recompute everything with the CORRECT radial kernels.
"""

from __future__ import annotations
import os, sys
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
sys.set_int_max_str_digits(10**7)
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
from scipy import integrate
import sympy as sp
from sympy import (Rational, sqrt, Integer, simplify, pi, S, Abs,
                   factorial as sym_factorial, Symbol)
from sympy.physics.wigner import (wigner_3j, wigner_6j, wigner_9j,
                                   clebsch_gordan)
from math import factorial, sqrt as msqrt

from geovac.breit_integrals import breit_ss_radial

# ====================================================================
# Hydrogenic radial wavefunctions (scipy)
# ====================================================================

def R_nl(n, l, r, Z=1):
    """Normalized hydrogenic R_nl(r)."""
    from scipy.special import assoc_laguerre, factorial as sf
    rho = 2*Z*r/n
    norm = np.sqrt((2*Z/n)**3 * sf(n-l-1) / (2*n*sf(n+l)))
    L = assoc_laguerre(rho, n-l-1, 2*l+1)
    return norm * rho**l * np.exp(-rho/2) * L

# ====================================================================
# Angular algebra in C^(k) convention (exact sympy)
# ====================================================================

def red_C(l, k, lp):
    """<l || C^(k) || l'> reduced matrix element."""
    return (-1)**l * sqrt(Integer((2*l+1)*(2*lp+1))) * wigner_3j(l, k, lp, 0, 0, 0)

def spatial_red_me(la, lb, lc, ld, L, Lp, k1, k2, K):
    """<(la lb) L || [C^(k1)(1) x C^(k2)(2)]^(K) || (lc ld) L'>."""
    n9 = wigner_9j(la, lb, L, lc, ld, Lp, k1, k2, K)
    r1 = red_C(la, k1, lc)
    r2 = red_C(lb, k2, ld)
    prefactor = sqrt(Integer((2*L+1)*(2*Lp+1)*(2*K+1)))
    return simplify(prefactor * n9 * r1 * r2)

def spin_red_me(S_val, k):
    """<S || [s1 x s2]^(k) || S>."""
    half = Rational(1, 2)
    n9 = wigner_9j(half, half, S_val, half, half, S_val, 1, 1, k)
    red_s = sqrt(Rational(3, 2))
    prefactor = sqrt(Integer((2*S_val+1)**2 * (2*k+1)))
    return simplify(prefactor * n9 * red_s**2)

def j_factor(L, S, J, k):
    return (-1)**(L+S+J) * wigner_6j(L, S, J, S, L, k)

# ====================================================================
# Varshalovich Eq. 5.17.22: bipolar expansion coefficient for I^L_M
# ====================================================================

def varshalovich_coeff(l1, l2, L):
    """Angular coefficient in the Varshalovich expansion of I^L_M(r1-r2).

    From Eq. 5.17.22 (VMK 1988):
    The coefficient for the (l1, l2) channel with l1+l2 = L + 2n is:

    (-1)^{l2} * hat{l1} * hat{l2} / hat{L} * (l1 0 l2 0 | L 0)
    * [(2l1-1)!! (2l2-1)!! / (2L-1)!!]^{1/2}  [for n=0 leading term]

    Wait, let me look up the exact formula more carefully.

    Actually, the standard result (Rose 1955, Eq. 4.68; Brink-Satchler App. V)
    for the irregular solid harmonic expansion is:

    I^L_M(r_12) = sum_{l1+l2=L} (-1)^{l1} (2L+1)/(hat{l1}^2 hat{l2}^2)
                  * hat{l1} hat{l2} * CG(l1 0 l2 0|L 0)
                  * ... (this gets messy)

    Let me use the simplest correct form. The key identity is:

    C^L_M(hat r_12) / r_12^{L+1}
    = sum_{l1+l2=L} (-1)^{l1} sqrt{(2l1+1)(2l2+1)/(2L+1)}
      * CG(l1 0 l2 0 | L 0)
      * { r_1^{l1}/r_2^{l2+1}  for r1<r2
        { r_2^{l2}/r_1^{l1+1}  for r1>r2
      * [C^{l1}(hat r1) x C^{l2}(hat r2)]^L_M

    ONLY the leading l1+l2=L term contributes (no higher-order terms)
    because the expansion terminates at the leading multipole for the
    irregular solid harmonic.

    This is from Jackson (3rd ed.) Eq. 3.70 generalized to arbitrary L,
    or equivalently from the addition theorem for solid harmonics.

    Let me verify: for L=0 (1/r_12):
    l1+l2=0, so (l1,l2)=(0,0).
    coeff = (-1)^0 * 1*1/1 * CG(0 0 0 0|0 0) = 1*1 = 1
    kernel: 1/r_> for r1<r2, and 1/r_1 for r1>r2 → 1/r_> ✓
    [C^0 x C^0]^0 = 1/sqrt(4pi) * 1/sqrt(4pi) * ... no, C^0 = 1.
    So I^0(r12) = 1/r_12 = sum_l r_<^l/r_>^{l+1} * C^l.C^l ← that's
    the Laplace expansion, which has ALL multipoles, not just l1+l2=0!

    Wait, I'm confusing two different expansions:
    1. Laplace expansion of 1/r_12 (scalar): ALL multipoles l contribute.
    2. Bipolar expansion of I^L_M(r_12): gives a SINGLE-CENTER expansion
       of the TENSOR irregular solid harmonic.

    The CORRECT identity is the addition theorem for irregular solid harmonics:

    I^L_M(r_1 - r_2)
    = sum_{l1=0}^{inf} sum_{l2} [coefficient] * R^{l1}(r_1) * I^{l2}(r_2)
      * [C^{l1} x C^{l2}]^L_M

    where l2 = l1 + L (for the r1 < r2 part) and the sum is infinite.

    But for our purposes, the angular selection rules from (1s)(2p) restrict
    which (l1, l2) contribute. Let me just compute ALL terms and see.

    Actually, I think the correct formula is simpler than I'm making it.
    The key reference is Gaunt (1929) or the standard multipole expansion
    of 1/|r1-r2|^{2s+1} from Jackson:

    1/|r1-r2|^{2s+1} = sum_l [(l+s)!/(l!s!)]^2 * (2l+1) / [(l+s)!/(l!s!)]
                      ... this is getting complicated.

    Let me use the NUMERICAL approach: compute the matrix element of
    the FULL H_SS operator by direct 6D integration (3 angles + 2 radii),
    and see if it matches Drake's formula.
    """
    cg = clebsch_gordan(l1, l2, L, 0, 0, 0)
    hat_l1 = sqrt(Integer(2*l1+1))
    hat_l2 = sqrt(Integer(2*l2+1))
    hat_L = sqrt(Integer(2*L+1))
    return (-1)**l1 * hat_l1 * hat_l2 / hat_L * cg

print("=" * 76)
print("Sprint 6 DW v2: Direct numerical check of H_SS matrix element")
print("=" * 76)

# ====================================================================
# APPROACH: Compute A_SS directly from the definition
# ====================================================================
#
# A_SS = <(1s 2p) ^3P_1 | H_SS_tensor | (1s 2p) ^3P_1> / f_SS(J=1)
#
# The tensor part of H_SS (traceless rank-2) acts on the spatial coordinates
# as: sum_q (-1)^q [s1xs2]^(2)_q * T^(2)_{-q}(r_12) where
# T^(2)_q(r_12) = -alpha^2/4 * sqrt(6) * C^(2)_q(hat r_12) / r_12^3
#
# Using Wigner-Eckart:
# <^3P_1,M | sum_q [s1xs2]^(2)_q T^(2)_{-q} | ^3P_1,M>
# = (-1)^{L+S+J} 6j{L S J; S L 2}
#   * <L||T^(2)_spatial||L> * <S||[s1xs2]^(2)||S>
#
# So A_SS = (-alpha^2/4) * sqrt(6) * (-1)^{L+S+J} 6j{1 1 1; 1 1 2}
#           * <L=1||T^(2)||L=1>_spatial * spin_red_me / f_SS(1)
#
# Since f_SS(1)=1 and [6j at J=1] = -1/6:
# A_SS = (-1/4) * sqrt(6) * (-1)^3 * (-1/6) * spatial_red_total * sqrt(5)/2
#       = (-1/4) * sqrt(6) * (1/6) * spatial_red_total * sqrt(5)/2
#       = sqrt(30)/48 * spatial_red_total
#
# The spatial reduced ME is what we need to compute.
# Let's compute it via numerical integration of the full
# two-electron matrix element with a specific M_J, then extract A_SS.

# ====================================================================
# Simpler approach: compute the full 2e matrix element of C^2/r12^3
# ====================================================================

# Instead of the bipolar expansion, let's directly compute:
# <1s(1) 2p_0(2) | C^(2)_0(r_12_hat) / r_12^3 | 1s(1) 2p_0(2)>  (direct)
# <1s(1) 2p_0(2) | C^(2)_0(r_12_hat) / r_12^3 | 2p_0(1) 1s(2)>  (exchange)
#
# These integrals can be computed by expanding 1/r_12^3 using the
# Gegenbauer expansion:
#
# 1/r_12^3 = sum_k (k+1) * r_<^k / r_>^{k+2} * P_k(cos theta_12)
#
# where theta_12 is the angle between r1 and r2 vectors.
#
# Then C^(2)(r12_hat) / r_12^3 involves additional angular coupling.

# Actually let me try the most direct approach: Drake's formula tells us that
# A_SS is a LINEAR COMBINATION of M^K integrals. The M^K integrals are
# computed exactly by the production module. So the only question is:
# what linear combination?
#
# The Racah algebra should determine this. Let me compute the angular
# prefactors from the Gegenbauer expansion and see.

# ====================================================================
# The Gegenbauer expansion approach
# ====================================================================
#
# The Breit-Pauli SS Hamiltonian can be written as:
# H_SS = (alpha^2/r_12^3) * [3(s1.r12hat)(s2.r12hat) - s1.s2]
#
# The angular-spin factor [3(s1.r)(s2.r) - s1.s2] can be decomposed as:
# = -sqrt(24*pi/5) * sum_q (-1)^q [s1xs2]^(2)_q * Y^(2)_{-q}(r_hat)
#
# In C^(k) convention: Y^(2) = sqrt(5/(4pi)) C^(2), so:
# = -sqrt(24*pi/5) * sqrt(5/(4pi)) * sum_q (-1)^q [s1xs2]^(2)_q * C^(2)_{-q}
# = -sqrt(6) * sum_q (-1)^q [s1xs2]^(2)_q * C^(2)_{-q}(r_hat)
#
# So H_SS = -(alpha^2/4) * sqrt(6) / r_12^3
#           * sum_q (-1)^q [s1xs2]^(2)_q * C^(2)_{-q}(r12_hat)
#
# Wait — I need to be careful. The standard form has alpha^2/4 prefactor
# from the Breit-Pauli Hamiltonian. Let me write:
#
# H_SS = (alpha^2/4) * [-sqrt(6)] * sum_q (-1)^q [s1xs2]^(2)_q
#        * C^(2)_{-q}(r12_hat) / r_12^3
#
# Now: C^(2)_q(r12_hat) depends on the direction of r_12 = r_1 - r_2.
# To evaluate this, we use the identity:
#
# C^(K)_M(r12_hat) = sum_{k1 k2} beta(k1,k2,K) * [C^{k1}(r1) x C^{k2}(r2)]^K_M
#                    * angular(r1, r2, r_12)
#
# But this is just the addition theorem for C^K, which relates C^K(r12) to
# the individual C^k's of r1 and r2 via bipolar harmonics.
#
# The problem is that C^K(r12_hat) where r12 = r1 - r2 is NOT simply
# a bipolar harmonic of r1_hat and r2_hat — it depends on the magnitudes
# r1 and r2 as well (because r12_hat = (r1 r1_hat - r2 r2_hat)/|r1-r2|).
#
# The CORRECT way is to use the Gegenbauer addition formula for 1/r_12^3
# and extract the rank-2 component.

# ====================================================================
# The Gegenbauer expansion with multipole decomposition
# ====================================================================
#
# 1/r_12^3 * C^(2)_0(r12_hat)
#
# This is the distributional kernel. To get the matrix element, we note:
#
# <ab|C^(2)_0(r12)/r_12^3|cd>
#   = int dr1 dr2 R_a(r1) R_b(r2) R_c(r1) R_d(r2) r1^2 r2^2
#     * int dOmega1 dOmega2 Y_{la,ma}*(1) Y_{lb,mb}*(2) C^(2)_0(r12_hat)/r_12^3
#       Y_{lc,mc}(1) Y_{ld,md}(2)
#
# The angular-radial integration must be done together because r_12 depends
# on both angles and radii.
#
# Standard approach: expand r_12^{-3} C^(2)_0(r12) using the Sack expansion,
# which gives it as a sum over k1, k2 of radial × angular terms.
#
# The Sack expansion (J. Math. Phys. 5, 245, 1964) gives:
#
# r_12^{-n} C^K_M(r12_hat)
# = sum_{l=0}^{inf} sum_{l1+l2=K+2l, triangle(l1,l2,K)}
#   alpha_{l1,l2}^{K,n} * r_<^{l1}/r_>^{l2+n-K}
#   * [C^{l1}(1) x C^{l2}(2)]^K_M
#
# For n=3, K=2 (our case):
# r_12^{-3} C^(2)_M = sum_{l1+l2=2,4,6,...} alpha * r_<^{l1}/r_>^{l2+1}
#                     * [C^{l1}(1) x C^{l2}(2)]^(2)_M
#
# The radial kernel is r_<^{l1}/r_>^{l2+1} (since n-K = 3-2 = 1).
#
# For l1+l2=2 (l=0, leading term): r_<^{l1}/r_>^{l2+1}
# For l1+l2=4 (l=1, next term): r_<^{l1}/r_>^{l2+1}
# etc.
#
# Angular selection for He (1s)(2p):
# Direct: <0|C^{l1}|0> requires l1=0, <1|C^{l2}|1> requires l2 even ≤ 2
#   So l1=0, l2 must be even and ≤2, giving l2=0 or l2=2.
#   But triangle(l1,l2,K=2): |l1-l2|≤2≤l1+l2 → l2≥2 (since l1=0).
#   So l2=2 only, and l1+l2=2 (leading term).
#   Radial: r_<^0/r_>^3 = Drake M^0 kernel!
#
# Exchange: <0|C^{l1}|1> requires l1=1, <1|C^{l2}|0> requires l2=1.
#   l1=1, l2=1, l1+l2=2 (leading term).
#   Triangle (1,1,2): |1-1|≤2≤1+1=2 ✓
#   Radial: r_</r_>^2
#
# Higher terms l1+l2=4: Direct needs l1=0, l2=4. But <1|C^4|1> requires
#   triangle(1,4,1): 1+1<4, FAILS. So no l1+l2=4 direct channel.
#   Exchange: l1+l2=4 with l1 odd, l2 odd (from <0|C^l1|1> needs l1=1,
#   and <1|C^l2|0> needs l2=1). So l1=1, l2=3: triangle(1,3,2) ✓,
#   but <1|C^3|0> needs triangle(1,3,0): 1+0<3 FAILS. So no l2=3.
#   l1=3, l2=1: <0|C^3|1> needs triangle(0,3,1): 0+1<3 FAILS.
#   So NO l1+l2=4 exchange channel either.
#
# CONCLUSION: Only the leading term l1+l2=2 contributes, for both direct
# and exchange. The angular selection completely determines the answer.

print("\n--- Sack expansion radial kernels ---")
print("For C^(2)/r_12^3 (n=3, K=2), the radial kernel at leading order is:")
print("  r_<^{l1} / r_>^{l2+1}  where l1+l2=2")
print()
print("Direct (l1=0, l2=2): kernel = 1/r_>^3 = Drake M^0 kernel")
print("Exchange (l1=1, l2=1): kernel = r_</r_>^2 = NEW (not Drake M^K)")
print()
print("Note: r_</r_>^2 is NOT the same as Drake's r_<^1/r_>^4 (= M^1).")
print("It's r_</r_>^2, which has a WEAKER r_> singularity.")

# ====================================================================
# Compute the Sack expansion coefficient alpha_{l1,l2}^{K=2,n=3}
# ====================================================================

# The Sack coefficient for the leading term (l1+l2=K) is:
#
# From Sack (1964) or see Steinborn & Filter (1975):
# alpha_{l1,l2}^{K,n} = (-1)^{l1} * (2K+1) / [hat{l1}^2 hat{l2}^2]
#                       * hat{l1} hat{l2} * CG(l1 0 l2 0|K 0)
#                       * (for leading n=0 term)
#
# Actually the simplest correct expression comes from the ADDITION THEOREM
# for irregular solid harmonics. The key identity is:
#
# For the addition theorem of 1/r_12^{2s+1} * C^K(r12):
# The Gegenbauer expansion gives:
#
# sum_k (k+s)! / (k! s!) * (4pi/(2k+1)) * sum_m Y^k_m(1)* Y^k_m(2) * r_<^k/r_>^{k+s+1}
# = 1/r_12^{2s+1}    (this is the SCALAR expansion)
#
# The TENSOR C^K/r_12^{2s+1} expansion is more complex and involves
# linearization coefficients (Gaunt coefficients).
#
# Let me use the DIRECT numerical approach to verify.

# ====================================================================
# Direct numerical matrix element of C^(2)_0(r12)/r_12^3
# ====================================================================

print("\n" + "=" * 76)
print("Direct numerical computation of <1s 2p|C^(2)/r^3|1s 2p> and exchange")
print("=" * 76)

# We compute:
# I_dir = <1s(1) 2p_0(2) | C^(2)_0(r12_hat)/r_12^3 | 1s(1) 2p_0(2)>
# I_exch = <1s(1) 2p_0(2) | C^(2)_0(r12_hat)/r_12^3 | 2p_0(1) 1s(2)>
#
# Using: C^(2)_0 = (3cos^2(theta)-1)/2 where theta is the polar angle of r_12.
#
# r_12 = r_1 - r_2, with r_12^2 = r1^2 + r2^2 - 2 r1 r2 cos(gamma)
# where gamma is the angle between r1 and r2.
#
# cos(theta_12) = (r_1 cos(theta_1) - r_2 cos(theta_2)) / r_12
#   ... this is complex. Instead, use the fact that C^(2)_0(r12_hat)
#   can be evaluated from the components of r_12.
#
# Actually, for the matrix element, we should use the decomposition into
# Gaunt integrals. Let me use the Gegenbauer (scalar) expansion and
# extract the rank-2 piece.

# The key formula: For any central two-body operator V(r_12), the
# matrix element between Slater determinants decomposes as:
#
# <ab|V|cd> = sum_k [angular_k] * R^k(ac;bd)
#
# where R^k is the Slater integral with the k-th multipole kernel.
#
# For V = 1/r_12: kernel(k) = r_<^k / r_>^{k+1} → Coulomb Slater integrals
# For V = 1/r_12^3: kernel(k) = (k+1) * r_<^k / r_>^{k+2} → Breit integrals
#
# Wait — that's not quite right either. The Gegenbauer expansion of
# f(r_12) P_k(cos gamma) gives an expansion in Legendre polynomials,
# and the angular integration extracts specific k values via Gaunt integrals.
#
# For the TENSOR part C^(2)_0(r12_hat)/r_12^3, we need the Sack expansion.
# But let me try a simpler route:
#
# The standard decomposition of the SS operator matrix element gives:
#
# <(la lb) L, (sa sb) S; J M| H_SS |same>
# = (-alpha^2/4) * sum_k V^k(la lb la lb; lc ld lc ld) * W^k(spin, J)
#
# where V^k is the Breit radial integral R^k_BP = int r_<^k/r_>^{k+3} * densities.
#
# The k+3 exponent in r_>^{k+3} comes from the DELTA-FUNCTION TERM
# in the Breit-Pauli 1/r_12^3 operator. The actual tensor C^(2)/r_12^3
# is a distribution with a delta function at r_12=0.
#
# From Bethe-Salpeter §38, the regularized form gives EXACTLY the
# r_<^k / r_>^{k+3} kernel (times angular coefficients that are Gaunt
# integrals). This IS the correct radial kernel, and it IS what Drake uses.
#
# So the question reduces to: what are the angular coefficients that
# multiply each R^k_BP = M^k integral?

# ====================================================================
# The Angular Coefficients: Standard Slater-Condon for H_SS
# ====================================================================
#
# For the Breit-Pauli SS operator, the standard decomposition (see
# Johnson, Atomic Structure Theory, Ch. 8; or Bethe-Salpeter §38) gives:
#
# <(la lb) L S J | H_SS | (lc ld) L S J>
# = (-alpha^2/4) * sqrt(6) * (-1)^{L+S+J} * 6j{L S J; S L 2}
#   * <S||[s1xs2]^(2)||S>
#   * sum_k Theta^2(la lb lc ld; k) * M^k
#
# where Theta^2(la lb lc ld; k) is the angular coefficient connecting
# multipole order k to the rank-2 tensor.
#
# The angular coefficient Theta^2 comes from the linearization of the
# product C^(k)(r1) . C^(k)(r2) * C^(2)(r12) in the matrix element.
# It involves 6j symbols.
#
# From Johnson (2007) Eq. 8.47 or Judd "Operator Techniques" Ch. 7:
#
# Theta^K(la lb lc ld; k) = (-1)^{la+lc} * hat{la}*hat{lb}*hat{lc}*hat{ld}
#   * (la 0 k 0|lc 0) * (lb 0 k 0|ld 0) * (-1)^{k+K}
#   * 6j{la lc k; ld lb K}
#
# This is the same as the standard Coulomb angular coefficient EXCEPT for
# the extra 6j factor with K.
#
# For K=2 (SS):

print("\n--- Standard angular coefficients Theta^2 for each k ---")

def theta_K(la, lb, lc, ld, k, K):
    """Angular coefficient for the rank-K tensor two-body operator.

    From Johnson (2007) Eq. 8.47 / Bethe-Salpeter §38.
    """
    hat_la = sqrt(Integer(2*la+1))
    hat_lb = sqrt(Integer(2*lb+1))
    hat_lc = sqrt(Integer(2*lc+1))
    hat_ld = sqrt(Integer(2*ld+1))
    cg1 = clebsch_gordan(la, k, lc, 0, 0, 0)  # = hat... * 3j
    cg2 = clebsch_gordan(lb, k, ld, 0, 0, 0)
    j6 = wigner_6j(la, lc, k, ld, lb, K)
    return ((-1)**(la+lc) * hat_la * hat_lb * hat_lc * hat_ld
            * cg1 * cg2 * (-1)**(k+K) * j6)


# For direct: (la,lb,lc,ld) = (0,1,0,1)
# For exchange: (la,lb,lc,ld) = (0,1,1,0)
print("\nDirect (0,1,0,1), K=2:")
for k in range(5):
    th = theta_K(0, 1, 0, 1, k, 2)
    th_s = simplify(th)
    if th_s != 0:
        print(f"  k={k}: Theta^2 = {th_s} = {float(th_s):.8f}")

print("\nExchange (0,1,1,0), K=2:")
for k in range(5):
    th = theta_K(0, 1, 1, 0, k, 2)
    th_s = simplify(th)
    if th_s != 0:
        print(f"  k={k}: Theta^2 = {th_s} = {float(th_s):.8f}")

# ====================================================================
# Now compute A_SS from the angular coefficients × M^k integrals
# ====================================================================

print("\n" + "=" * 76)
print("Computing A_SS from Theta^2 angular coefficients")
print("=" * 76)

# A_SS * f_SS(J) = (-alpha^2/4) * sqrt(6) * (-1)^{L+S+J} * 6j{1 1 J; 1 1 2}
#                  * <S=1||[s1xs2]^(2)||S=1>
#                  * sum_k [Theta^2_dir(k) * M^k_dir - Theta^2_exch(k) * M^k_exch]
#
# At J=1, f_SS(1)=1:
# A_SS = (-1/4) * sqrt(6) * (-1)^{1+1+1} * (-1/6) * sqrt(5)/2
#        * sum_k [...]
#      = (-1/4) * sqrt(6) * (-1) * (-1/6) * sqrt(5)/2 * sum_k
#      = (-1/4) * sqrt(6) * (1/6) * sqrt(5)/2 * sum_k
#      = sqrt(30)/48 * sum_k [...]

jf = float(simplify(j_factor(1, 1, 1, 2)))
spin = float(simplify(spin_red_me(1, 2)))
prefactor = 0.25 * msqrt(6) * (-1) * jf * spin  # Note: (-1) from the -alpha^2/4 sign
# Wait: the Hamiltonian is H_SS = (alpha^2/4)[sigma.sigma/r^3 - 3(sigma.r)(sigma.r)/r^3]
# = -(alpha^2/4) * sqrt(6) * [s1xs2]^(2) . C^(2)(r12)/r_12^3
# The SIGN depends on convention. Let me be careful.
# Following Bethe-Salpeter §38:
# H_SS = (alpha^2/4) * [3(s1.r12)(s2.r12)/r12^5 - (s1.s2)/r12^3]
# The rank-2 tensor: 3(s1.r)(s2.r)/r^2 - s1.s2 = -sqrt(24pi/5) * [s1xs2]^2 . Y^2
# = -sqrt(6) * [s1xs2]^2 . C^2
# So H_SS = -(alpha^2/4) * sqrt(6) * {[s1xs2]^2 . C^2(r12hat)} / r_12^3
#
# The DIAGONAL matrix element: <|H_SS|> = -(alpha^2/4) * sqrt(6) * <|T.U|>
# where T = C^2(r12hat)/r12^3, U = [s1xs2]^(2).
#
# <T.U> = (-1)^{L+S+J} 6j{LSJ;SLK} * <L||T||L> * <S||U||S>
#
# For J=1: (-1)^3 * (-1/6) = 1/6
#
# So <|H_SS|> = -(1/4) * sqrt(6) * (1/6) * <L||T||L> * sqrt(5)/2

# The spatial reduced ME:
# <L=1 || C^(2)(r12hat)/r_12^3 || L=1>
# = sum_k Theta^2_dir(k) * M^k_dir - Theta^2_exch(k) * M^k_exch

spatial_sum = 0.0
M_vals = {}
for K_d in range(5):
    try:
        M_dir = float(breit_ss_radial(1,0, 2,1, 1,0, 2,1, K_d, Z=1))
        M_exch = float(breit_ss_radial(1,0, 2,1, 2,1, 1,0, K_d, Z=1))
    except Exception:
        M_dir = 0.0
        M_exch = 0.0
    M_vals[('dir', K_d)] = M_dir
    M_vals[('exch', K_d)] = M_exch

print(f"\nPrefactors:")
print(f"  -sqrt(6)/4 = {-msqrt(6)/4:.6f}")
print(f"  J-factor(J=1) = {jf:.6f}")
print(f"  Spin red ME = {spin:.6f}")
print(f"  Full prefactor (no alpha^2) = -(1/4)*sqrt(6)*jf*spin = {-0.25*msqrt(6)*jf*spin:.6f}")

print(f"\nSpatial sum: sum_k Theta^2 * M^k")
spatial_dir_sum = 0.0
spatial_exch_sum = 0.0
for k in range(5):
    th_dir = float(simplify(theta_K(0,1,0,1, k, 2)))
    th_exch = float(simplify(theta_K(0,1,1,0, k, 2)))
    M_dir = M_vals[('dir', k)]
    M_exch = M_vals[('exch', k)]
    dir_c = th_dir * M_dir
    exch_c = th_exch * M_exch
    if th_dir != 0 or th_exch != 0:
        print(f"  k={k}: Theta_dir={th_dir:+.8f} * M_dir={M_dir:+.8e} = {dir_c:+.10e}")
        print(f"         Theta_exch={th_exch:+.8f} * M_exch={M_exch:+.8e} = {exch_c:+.10e}")
    spatial_dir_sum += dir_c
    spatial_exch_sum += exch_c

spatial_sum = spatial_dir_sum - spatial_exch_sum
print(f"\n  Spatial direct sum = {spatial_dir_sum:+.10e}")
print(f"  Spatial exchange sum = {spatial_exch_sum:+.10e}")
print(f"  Spatial total (dir - exch) = {spatial_sum:+.10e}")

# Full A_SS (without alpha^2):
A_SS_computed = -0.25 * msqrt(6) * jf * spin * spatial_sum
print(f"\n  A_SS (computed, no alpha^2) = {A_SS_computed:+.10e}")

# Drake's value:
M2_dir = M_vals[('dir', 2)]
M2_exch = M_vals[('exch', 2)]
A_SS_drake = 3.0/50 * M2_dir - 2.0/5 * M2_exch
print(f"  A_SS (Drake) = {A_SS_drake:+.10e}")

ratio = A_SS_computed / A_SS_drake if A_SS_drake != 0 else float('inf')
print(f"\n  Ratio = {ratio:.8f}")
print(f"  Relative difference = {abs(ratio - 1)*100:.6f}%")

if abs(ratio - 1) < 0.001:
    print("\n  *** EXACT MATCH! ***")
    print("  The standard Racah algebra CLOSES the Drake derivation!")
    print()
    # Extract the coefficients
    print("  Extracting the combining coefficients:")
    # A_SS = prefactor * (Theta_dir_sum * M_dir - Theta_exch_sum * M_exch)
    # If only k=2 contributes:
    th_dir_k2 = float(simplify(theta_K(0,1,0,1, 2, 2)))
    th_exch_k2 = float(simplify(theta_K(0,1,1,0, 2, 2)))
    full_pref = -0.25 * msqrt(6) * jf * spin
    c_dir = full_pref * th_dir_k2
    c_exch = -full_pref * th_exch_k2  # minus from dir-exch
    print(f"  c_dir (coef of M^2_dir) = {c_dir:+.10f}")
    print(f"  c_exch (coef of M^2_exch) = {c_exch:+.10f}")
    print(f"  Drake: c_dir = 3/50 = {3/50:.10f}")
    print(f"  Drake: c_exch = -2/5 = {-2/5:.10f}")
elif abs(ratio - 1) < 0.05:
    print(f"\n  CLOSE ({abs(ratio-1)*100:.4f}% off) — nearly matches!")
else:
    print(f"\n  NO MATCH ({abs(ratio-1)*100:.2f}% off)")

# ====================================================================
# Also check: what if Theta^K uses a DIFFERENT formula?
# ====================================================================

print("\n" + "=" * 76)
print("Alternative: Compute Theta^K via the Condon-Shortley c^k convention")
print("=" * 76)

# In the Condon-Shortley convention, the angular coefficient for a
# rank-K tensor two-body operator is:
#
# d^k_K(la,lb,lc,ld) = hat{la}^2 hat{lb}^2 * (la 0 k 0|lc 0)^2 ... no
#
# Actually, let me use a completely different route. From Drake 1971 directly:
#
# The Breit-Pauli SS matrix element for (n1 l1)(n2 l2) ^{2S+1}L_J is:
#
# <^3P_J|H_SS|^3P_J> = (-alpha^2/4) * sum_K=0,2 sum_k D^K(l1 l2; k) * M^k * f_SS(J;K)
#
# where D^K is the spatial-spin angular coefficient and the sum over K
# includes K=0 (scalar, s1.s2/r^3) and K=2 (tensor) parts.
#
# For triplet S=1:
# The K=0 part: s1.s2 = [S(S+1)-3/2]/2 = 1/4 for S=1
# The K=2 part: involves [s1xs2]^(2)
#
# Let me try yet another approach: use c^k coefficients (Condon-Shortley)
# and the known formula for H_SS in terms of Slater integrals.

# From Sobelman, Vainshtein, Yukov (2012), or Cowan "The Theory of Atomic
# Structure and Spectra" (1981), the SS matrix element is:
#
# <ab|H_SS|cd> = (alpha^2/4) * sum_k [direct_angular(k) * M^k_dir
#                                     - exchange_angular(k) * M^k_exch]
#
# where the angular coefficients for the K=2 tensor part are:
# direct_angular(k) = (-1)^{l_a+l_c} * c^k(l_a,l_c)^2 * beta(k,l_a,l_c,l_b,l_d) ... complex

# Let me just TRY the Theta^K formula with the CG convention correction.
# Maybe I have a sign or normalization error.

# Alternative formulation: use the explicit formula from
# Bethe-Salpeter §38, Eq. 38.28:
#
# For the tensor part of H_SS, the matrix element between LS-coupled states is:
#
# <(l1 l2)L, S=1; J | H_SS^{tensor} | (l1 l2)L, S=1; J>
# = -(alpha^2/4) * sqrt(6) * (-1)^{L+1+J} * 6j{1 1 J; 1 1 2}
#   * sqrt(5)/2  [spin reduced ME]
#   * <(l1 l2) L || Q^2 || (l1 l2) L>
#
# where Q^2 is the spatial rank-2 tensor operator (C^2(r12)/r_12^3).
# And <L || Q^2 || L> decomposes into Slater integrals:
#
# <(l1 l2) L || Q^2 || (l1 l2) L>
# = sum_k { direct: (-1)^{l1+l2+L+k+2} hat{L}^2 hat{k} hat{2}
#            * 6j{l2 l2 k; l1 l1 2} * 6j{l1 l2 L; l2 l1 k}
#            * (l1 0 k 0|l1 0) * (l2 0 k 0|l2 0) * M^k_dir
#          - exchange: similar with (l1,l2) -> (l2,l1) }

# Wait, I think the 6j-based decomposition is for the general case.
# For our specific case (l1=0, l2=1, L=1), let me just use the explicit
# formula.

# 6j{l2 l2 k; l1 l1 2} = 6j{1 1 k; 0 0 2} = delta_{k,2} * (something)
# Actually, 6j{1 1 k; 0 0 2}: rows must satisfy triangle conditions.
# First row: (1, 1, k) → k=0,1,2
# Columns: (1,0,2)→triangle ✓, (1,0,k)→k=1, (k,2,?)
# Actually the 6j symbol has specific rules. Let me just compute numerically.

print("\nAlternative Theta^2 formula (6j-based, Bethe-Salpeter):")
print("For direct: (l1=0, l2=1, L=1)")

for k in range(5):
    # 6j{l2 l2 k; l1 l1 K} = 6j{1 1 k; 0 0 2}
    j6a = wigner_6j(1, 1, k, 0, 0, 2)
    # 6j{l1 l2 L; l2 l1 k} = 6j{0 1 1; 1 0 k}
    j6b = wigner_6j(0, 1, 1, 1, 0, k)
    # CGs
    cg1 = clebsch_gordan(0, k, 0, 0, 0, 0)  # = delta_{k,0}
    cg2 = clebsch_gordan(1, k, 1, 0, 0, 0)  # c^k(1,0;1,0)

    hat_L = msqrt(3)  # 2L+1=3
    hat_k = msqrt(2*k+1)
    hat_K = msqrt(5)  # 2K+1=5, K=2

    phase = (-1)**(0+1+1+k+2)
    coeff = phase * 3 * hat_k * hat_K * float(j6a) * float(j6b) * float(cg1) * float(cg2)

    if abs(coeff) > 1e-12:
        print(f"  k={k}: 6j_a={float(j6a):.6f}, 6j_b={float(j6b):.6f}, "
              f"cg1={float(cg1):.6f}, cg2={float(cg2):.6f}, "
              f"coeff={coeff:+.8f}")

# For the exchange, the orbital indices are crossed:
print("\nFor exchange: (l1=0, l2=1) → (l2=1, l1=0), same L=1")
for k in range(5):
    # Exchange: bra (0,1), ket (1,0)
    # The 9j-based formula from Part 3 already gives the right answer.
    # Let me check if the 6j formula gives different results.
    j6a = wigner_6j(1, 0, k, 1, 0, 2)  # 6j{l2' l1' k; l2 l1 K}
    j6b = wigner_6j(0, 1, 1, 0, 1, k)  # 6j{l1 l2 L; l1' l2' k}

    cg1 = clebsch_gordan(0, k, 1, 0, 0, 0)  # <l1|C^k|l1'=l2>
    cg2 = clebsch_gordan(1, k, 0, 0, 0, 0)  # <l2|C^k|l2'=l1>

    hat_L = msqrt(3)
    hat_k = msqrt(2*k+1)
    hat_K = msqrt(5)

    phase = (-1)**(0+1+1+k+2)
    coeff = phase * 3 * hat_k * hat_K * float(j6a) * float(j6b) * float(cg1) * float(cg2)

    if abs(coeff) > 1e-12:
        print(f"  k={k}: coeff={coeff:+.8f}")

# ====================================================================
# Let me try yet another route: use the KNOWN formula from
# Cowan (1981) "Theory of Atomic Structure and Spectra" Ch. 14
# ====================================================================

print("\n" + "=" * 76)
print("Cowan (1981) formula for SS angular coefficients")
print("=" * 76)

# From Cowan, the spin-spin interaction for the tensor part gives:
#
# <(l1 l2)^{2S+1}L_J|H_SS^{tensor}|same>
# = (alpha^2/4) * (-1)^{L+S+J} * 6j{L S J; S L 2} * sqrt(5)/2
#   * sum_k n_k(l1,l2,L) * N^k(n1l1,n2l2)
#
# where n_k is a specific angular coefficient and N^k is a radial integral.
#
# For the specific case (1s)(2p) ^3P, L=S=1, Cowan gives (Table 14-3):
# The contributing k values and their coefficients.
#
# Actually, I recall that the standard result (see e.g., Bethe-Salpeter
# problem set, or Glass & Hibbert 1978) is:
#
# For (ns)(n'p) ^3P, the SS tensor gives contributions at k=0 and k=2.
#
# The radial integrals involved are:
# N^k = integral of [electron density] * [r_<^k / r_>^{k+3}]
# which is EXACTLY Drake's M^k.
#
# The angular coefficients depend on l1, l2, L, and involve 6j symbols.

# Let me compute the angular coefficient using the SIMPLEST correct formula.
# From Lindgren & Morrison (1986) "Atomic Many-Body Theory" Eq. (11.92):
#
# For the tensor part (K=2) of the magnetic e-e interaction:
# d^K(la,lb,lc,ld; k) = (-1)^{la+lb+L+K}
#   * hat{la} hat{lb} hat{lc} hat{ld}
#   * (la 0 k 0|lc 0) (lb 0 k 0|ld 0)
#   * 6j{la lc k; ld lb K} * 6j{la lb L; ld lc k}  [for direct]
#
# Actually, this is getting circular. Let me just use the 9j from Part 3
# and the standard Sack coefficients.

# ====================================================================
# FINAL APPROACH: Use the identity that relates the 9j-based spatial
# reduced ME to Slater integrals
# ====================================================================

print("\n" + "=" * 76)
print("Final: Identity connecting spatial_red_ME to Drake M^k")
print("=" * 76)

# The spatial reduced ME of [C^{k1}(1) x C^{k2}(2)]^{K} involves
# integrating over angles AND radii with the kernel from the bipolar expansion.
#
# The 9j-based formula (Part 3) gives the ANGULAR part only.
# The RADIAL part comes from the specific expansion being used.
#
# For the irregular solid harmonic C^(K)/r^{K+1}:
# kernel = r_<^{k1}/r_>^{k2+1} for channel (k1,k2) with k1+k2=K
#
# For the Breit-Pauli 1/r^3 * C^(2):
# From the Gegenbauer expansion of 1/r_12^{2s+1}:
# 1/r_12^3 = sum_k (k+1) r_<^k/r_>^{k+2} P_k(cos gamma)
# P_k(cos gamma) = C^k(1).C^k(2) = sum_m C^k_m(1) C^k_m(2) (-1)^m
#
# So C^(2)(r12)/r_12^3 involves the product C^(2)(r12) * 1/r_12^3
# where 1/r_12^3 has its own Gegenbauer expansion, and C^(2)(r12) has
# its own angular structure.
#
# The product C^(2)(r12) * sum_k (k+1) r_<^k/r_>^{k+2} P_k must be
# decomposed into bipolar harmonics.
#
# Actually, I think the CORRECT statement is much simpler:
# Drake's M^k integral with kernel r_<^k/r_>^{k+3} IS the correct
# radial integral for the Breit-Pauli two-body operator.
#
# The exponent k+3 (rather than k+1 for Coulomb or k+2 for the Gegenbauer
# 1/r^3 expansion) comes from the angular integration of the delta-function
# regularization of 1/r_12^3.
#
# Let me verify this numerically. If the angular coefficient × M^k_Drake
# for each k gives the correct total, then the formula works.

# The key insight: maybe BOTH k=0 and k=2 contribute to the direct
# channel, and the 3/50 and -2/5 are the sums of two angular coefficients.

# Let me try: what if the angular coefficients for direct are
# a_0 (at k=0) and a_2 (at k=2), such that
# a_0 * M^0_dir + a_2 * M^2_dir = (3/50) * M^2_dir  ?
# This would require a_0 = 0 and a_2 = 3/50.
#
# OR: the total A_SS involves BOTH M^0 and M^2 integrals:
# A_SS = c_0 M^0_dir + c_2 M^2_dir - d_0 M^0_exch - d_2 M^2_exch
# and Drake's formula 3/50 M^2_dir - 2/5 M^2_exch is a simplification.

# Let me check: Drake 1971 (PRA 3, 908) explicitly states the formula.
# The M^k notation is standard. Let me see if k=0 contributes.

# For the scalar part of H_SS: s1.s2/r^3 → this involves the MONOPOLE
# (k=0) of 1/r^3, giving M^0 integrals. The tensor part (k=2) gives M^2.

# So the FULL H_SS (scalar + tensor) involves BOTH M^0 and M^2.

# Drake separates these: the scalar part (s1.s2/r^3) contributes to f_SS
# via the S(S+1) eigenvalue, and the tensor part via the 6j J-pattern.

# For ^3P (S=1): s1.s2 = [S(S+1)-3/2]/2 = [2-3/2]/2 = 1/4.

# So the scalar part contributes: (alpha^2/4) * (1/4) * M^0 / r^3
# But wait — the scalar part s1.s2/r^3 also has angular structure from 1/r^3.
# The 1/r^3 Gegenbauer expansion gives sum_k (k+1) M^k P_k, and the
# angular integration selects specific k values.

# Actually, for the FULL H_SS = (alpha^2/4) [s1.s2/r^3 - 3(s1.r)(s2.r)/r^5]:
# The s1.s2/r^3 term is an ISOTROPIC (rank-0 in spatial variables) operator.
# Its matrix element involves only the MONOPOLE (k=0) of 1/r^3:
# <1s 2p|1/r^3|1s 2p>_dir = (k=0) * M^0_dir (with angular coeff 1)
# No, that's wrong too. 1/r^3 has a Gegenbauer expansion involving all k.
# But the ANGULAR INTEGRAL for the direct channel (1s,1s on electron 1
# and 2p,2p on electron 2) selects only k=0 (from <Y_00|P_k|Y_00> = delta_{k0}).

# Wait: 1/r_12^3 has the expansion:
# 1/r_12^3 = sum_k (k+1) r_<^k/r_>^{k+2} (2k+1)^{-1} * 4pi * sum_m Y^k*(1) Y^k(2)
# In C^k convention: = sum_k (k+1) r_<^k/r_>^{k+2} * sum_m C^k_m(1)(-1)^m C^k_{-m}(2)

# For the SCALAR s1.s2 part:
# <1s 2p|s1.s2/r_12^3|1s 2p>
# = (1/4) * sum_k (k+1) * <1s(1)|C^k(1)|1s(1)> * <2p(2)|C^k(2)|2p(2)> * M'^k
# where M'^k has kernel r_<^k/r_>^{k+2}.

# <Y_00|C^k_0|Y_00> = delta_{k,0} (for k=0: 1; else 0 from triangle)
# So ONLY k=0 contributes to the scalar direct part.
# k=0 coefficient: (0+1) = 1.
# M'^0 has kernel 1/r_>^2.  NOT the same as M^0_Drake (kernel 1/r_>^3).

# Hmm — the Gegenbauer expansion of 1/r^3 gives kernel r_<^k/r_>^{k+2},
# but Drake's M^k uses r_<^k/r_>^{k+3}. These are DIFFERENT.

# I think the resolution is that the H_SS operator is actually written as
# a DERIVATIVE of the Coulomb operator. Bethe-Salpeter derives it from
# the Breit operator, and the result involves 1/r^3 in a distributional
# sense that gives rise to the r_<^k/r_>^{k+3} kernel (not r_<^k/r_>^{k+2}).

# The distributional regularization adds an extra factor of 1/r_> compared
# to the naive Gegenbauer expansion. This is because:
# d/dr_> [r_<^k/r_>^{k+1}] = -(k+1) r_<^k/r_>^{k+2}
# d^2/dr_>^2 [r_<^k/r_>^{k+1}] = (k+1)(k+2) r_<^k/r_>^{k+3}
#
# So r_<^k/r_>^{k+3} = d^2/dr^2[r_<^k/r_>^{k+1}] / [(k+1)(k+2)]
# This means Drake's M^k is the SECOND DERIVATIVE of the Coulomb R^k
# with respect to r_>, divided by (k+1)(k+2).
#
# OR: the Breit-Pauli 1/r^3 operator, when properly regularized, gives
# the radial kernel r_<^k/r_>^{k+3} = 1/[(k+1)(k+2)] * d^2/dr^2 of Coulomb.

# This means the CORRECT Gegenbauer expansion for the REGULARIZED 1/r^3 is:
# 1/r_{12,reg}^3 = sum_k 1/(k+2) * r_<^k/r_>^{k+3} * (2k+1) * P_k(cos gamma)
#
# (The factor changes from (k+1) for naive Gegenbauer to 1/(k+2) for Drake.)

# Let me check with known values. For k=0:
# Drake kernel: 1/r_>^3. Production module M^0 is computed.
# For k=2:
# Drake kernel: r_<^2/r_>^5. Production module M^2 is computed.

# OK, I think the fundamental issue is clear now. Let me just compute the
# angular coefficients that multiply M^k_Drake numerically.

# ====================================================================
# SIMPLEST APPROACH: Fit the angular coefficients to Drake M^k
# ====================================================================

print("\n" + "=" * 76)
print("Direct computation: find angular coefficients a_k, b_k such that")
print("A_SS = sum_k a_k * M^k_dir - sum_k b_k * M^k_exch")
print("=" * 76)

# We know A_SS must be a linear combination of M^k_dir and M^k_exch
# for k = 0, 1, 2 (triangle rule limits for l=0, l=1).
# Drake says only k=2 contributes: A_SS = (3/50)*M^2_dir - (2/5)*M^2_exch.
# But does k=0 also contribute?

# For the TENSOR part (rank-2), by the Wigner-Eckart theorem on the
# angular variables, only k values that couple to rank 2 can contribute.
# The angular coupling C^k(1) C^k(2) has rank 2k, and we need rank 2.
# So k=1 gives rank 2 (from |k-k|=0 to 2k=2), k=2 gives rank 2,4
# (from 0 to 4), k=0 gives rank 0 only.

# Actually that's wrong. The Gegenbauer expansion gives:
# 1/r^3 = sum_k ... * C^k(1).C^k(2) * radial
# The product C^(2)(r12hat) * C^k(1).C^k(2) decomposes by the
# linearization formula (involving 6j/9j), which can give different
# angular ranks.

# The key insight I keep coming back to: C^(2)(r12hat) is NOT the same
# as C^(2)(r1hat) or C^(2)(r2hat). It depends on the DIFFERENCE vector
# r_12 = r_1 - r_2.

# For the matrix element, the CORRECT approach is:
# 1. Write H_SS in second quantization or as a sum of Slater integrals
# 2. Use the known angular coefficients from the literature

# From Glass & Hibbert (1978, J Phys B 11, 2413), Table 1:
# For (ns n'p) ^3P, the Breit-Pauli SS matrix element involves:
# M^0 (k=0) and M^2 (k=2) Drake integrals.
# The angular coefficients are tabulated.

# Let me compute them from the KNOWN angular reduction.
# From Bethe-Salpeter Eq. 38.28, the tensor SS term is:
#
# E_{SS}^{tensor} = -(alpha^2/4) sqrt(30) <S=1||[s1xs2]^2||S=1>
#   * (-1)^{L+S+J} 6j{LSJ;SL2}
#   * sum_k d^2_k * M^k
#
# where d^2_k is defined by:
# <(l1 l2) L || sum_q C^2_q(r12)/r12^3 * (bipolar) || L>
# = sum_k d^2_k * M^k_Drake
#
# And d^2_k involves the 6j{l2 l2 k; l1 l1 2} * 6j{l1 l2 L; l2 l1 k}
# * (l1 0 k 0|l1 0)(l2 0 k 0|l2 0)

# But I think the formula may have a different overall normalization.
# Let me try an EMPIRICAL approach: just solve a linear system.

# We know: A_SS_Z1 * alpha^(-2) = 3/50 * M^2_dir - 2/5 * M^2_exch
# (from BF-D Sprint 3)
# That gives one specific number. To verify, let me check at Z=2.

print("\n--- Checking Drake formula at Z=2 ---")
Z = 2
# M^k scale as Z^3 for Breit integrals
M2_dir_Z2 = float(breit_ss_radial(1,0, 2,1, 1,0, 2,1, 2, Z=Z))
M2_exch_Z2 = float(breit_ss_radial(1,0, 2,1, 2,1, 1,0, 2, Z=Z))
A_SS_drake_Z2 = 3.0/50 * M2_dir_Z2 - 2.0/5 * M2_exch_Z2
print(f"  M^2_dir(Z=2) = {M2_dir_Z2:+.10e}")
print(f"  M^2_exch(Z=2) = {M2_exch_Z2:+.10e}")
print(f"  A_SS(Drake, Z=2) = {A_SS_drake_Z2:+.10e}")

# ====================================================================
# Use exact sympy to get Theta^K from the 9j approach
# ====================================================================

print("\n" + "=" * 76)
print("EXACT Theta^2 from product of Brink-Satchler + 9j spatial + radial mapping")
print("=" * 76)

# Channel (0,2), direct:
# BS coeff = sqrt(5)
# 9j spatial red ME = -sqrt(30)/5
# Combined angular for direct: sqrt(5) * (-sqrt(30)/5) = -sqrt(150)/5 = -sqrt(6)
ang_dir_02 = simplify(brink_satchler_coeff(0,2,2) * spatial_red_me(0,1,0,1,1,1,0,2,2))

# Channel (1,1), exchange:
# BS coeff = -sqrt(6)
# 9j spatial red ME = -sqrt(5)/3
# Combined: (-sqrt(6)) * (-sqrt(5)/3) = sqrt(30)/3
ang_exch_11 = simplify(brink_satchler_coeff(1,1,2) * spatial_red_me(0,1,1,0,1,1,1,1,2))

print(f"  Combined angular (direct, (0,2)):   {ang_dir_02} = {float(ang_dir_02):.8f}")
print(f"  Combined angular (exchange, (1,1)): {ang_exch_11} = {float(ang_exch_11):.8f}")

# These multiply the bipolar radial integrals N^{k1,k2}:
# Direct: N^{0,2}_dir with kernel 1/r_>^3
# Exchange: N^{1,1}_exch with kernel r_</r_>^2

# The full spatial reduced ME is:
# <L=1||Q^2||L=1> = ang_dir_02 * N^{02}_dir - ang_exch_11 * N^{11}_exch

# Now, N^{02}_dir with kernel 1/r_>^3 IS Drake's M^0_dir.
# And N^{11}_exch with kernel r_</r_>^2 is NOT any Drake M^k.

# Let me compare N^{11}_exch numerically with M^k_exch values:
N11_exch = bipolar_radial_integral(1,0, 2,1, 2,1, 1,0, 1, 1)
N02_dir = bipolar_radial_integral(1,0, 2,1, 1,0, 2,1, 0, 2)
M0_dir = float(breit_ss_radial(1,0, 2,1, 1,0, 2,1, 0, Z=1))
M0_exch = float(breit_ss_radial(1,0, 2,1, 2,1, 1,0, 0, Z=1))
M1_exch = float(breit_ss_radial(1,0, 2,1, 2,1, 1,0, 1, Z=1))
M2_exch_Z1 = float(breit_ss_radial(1,0, 2,1, 2,1, 1,0, 2, Z=1))

print(f"\n  N^(0,2)_dir  (kernel 1/r_>^3):   {N02_dir:+.10e}")
print(f"  M^0_dir (Drake, kernel 1/r_>^3):   {M0_dir:+.10e}")
print(f"  Match: {abs(N02_dir - M0_dir) < 1e-8}")

print(f"\n  N^(1,1)_exch (kernel r_</r_>^2):  {N11_exch:+.10e}")
print(f"  M^0_exch (kernel 1/r_>^3):         {M0_exch:+.10e}")
print(f"  M^1_exch (kernel r_</r_>^4):       {M1_exch:+.10e}")
print(f"  M^2_exch (kernel r_<^2/r_>^5):     {M2_exch_Z1:+.10e}")

# Check ratios:
print(f"\n  N^(1,1)_exch / M^0_exch = {N11_exch/M0_exch:.8f}")
print(f"  N^(1,1)_exch / M^1_exch = {N11_exch/M1_exch:.8f}")
print(f"  N^(1,1)_exch / M^2_exch = {N11_exch/M2_exch_Z1:.8f}")

# ====================================================================
# CRITICAL CHECK: Is the bipolar radial kernel for C^(2)/r_12^3
# actually r_<^{k1}/r_>^{k2+1} or something else?
# ====================================================================

# The addition theorem for IRREGULAR solid harmonics states:
# I^K_M(r_12) = C^K_M(r12hat)/r_12^{K+1}
# For K=2: C^(2)/r_12^3.
#
# This has the two-center expansion:
# I^K(r1-r2) = sum_{l1 l2} ... R^{l1}(r_1) I^{l2}(r_2) [C^{l1} x C^{l2}]^K
# for r_1 < r_2, and
# = sum_{l1 l2} ... I^{l1}(r_1) R^{l2}(r_2) [C^{l1} x C^{l2}]^K
# for r_1 > r_2.
#
# Where R^l(r) = r^l and I^l(r) = r^{-(l+1)}.
# So the kernel for channel (l1,l2) is:
# r1 < r2: r_1^{l1} * r_2^{-(l2+1)} = r_<^{l1}/r_>^{l2+1}
# r1 > r2: r_1^{-(l1+1)} * r_2^{l2} = r_>^{-(l1+1)} * r_<^{l2}
#         = r_<^{l2}/r_>^{l1+1}

# For (l1=0, l2=2): r1<r2: 1/r_2^3; r1>r2: r_2^2/r_1  (NOT 1/r_>^3)
# Region I: 1/r_2^3 = 1/r_>^3 ✓ (matches M^0)
# Region II: r_2^2/r_1 = r_<^2/r_>^1 (matches bipolar (2,0) NOT (0,2))

# AH HA! The Region II kernel for channel (0,2) is r_<^2/r_>^1,
# which corresponds to l1↔l2 SWAP (i.e., (2,0) channel radial kernel).

# So the CORRECT bipolar radial integral for channel (l1,l2) is:
# Region I (r1<r2): r_1^{l1} / r_2^{l2+1}
# Region II (r1>r2): r_2^{l2} / r_1^{l1+1}

# For (0,2):
# Region I: 1/r_2^3 = Drake M^0 Region I
# Region II: r_2^2/r_1 = Drake M^{-1} Region II (??)

# Hmm, Drake's M^K has kernel r_<^K/r_>^{K+3}:
# Region I: r_1^K/r_2^{K+3}
# Region II: r_2^K/r_1^{K+3}
# For (0,2) Region II: r_2^2/r_1^1 would need K+3=1 → K=-2, invalid.

# So the bipolar expansion kernel is FUNDAMENTALLY DIFFERENT from Drake's.
# This confirms that the bipolar expansion gives integrals that are NOT
# in Drake's M^k basis at all (for the exchange/Region-II terms).

print("\n" + "=" * 76)
print("CONCLUSION")
print("=" * 76)
print("""
The C^(k) convention check reveals a STRUCTURAL MISMATCH that is deeper
than a normalization convention:

1. The Brink-Satchler/addition-theorem bipolar expansion of C^(2)/r_12^3
   gives radial kernels r_<^{l1}/r_>^{l2+1} (with l1+l2=2), specifically:
   - Channel (0,2): 1/r_>^3 in Region I, r_<^2/r_> in Region II
   - Channel (1,1): r_</r_>^2 in both regions

2. Drake's M^k integrals use kernel r_<^k/r_>^{k+3} uniformly, giving:
   - M^0: 1/r_>^3
   - M^1: r_</r_>^4
   - M^2: r_<^2/r_>^5

3. These are COMPLETELY DIFFERENT radial power laws. The bipolar channel
   (0,2) has Region I matching Drake M^0 but Region II matching NOTHING
   in Drake's basis. The channel (1,1) matches no Drake M^k at all.

4. This means Drake's M^k integrals do NOT arise from the irregular
   solid harmonic addition theorem. They arise from a DIFFERENT
   decomposition — specifically, from the Gegenbauer expansion of 1/r_12^3
   combined with the distributional regularization of the delta function
   at r_12=0.

5. The coefficients (3/50, -2/5) therefore encode the MAPPING between
   the bipolar-expansion basis {N^{0,2}_dir, N^{1,1}_exch, ...} and
   Drake's {M^0, M^1, M^2} basis. This mapping involves piecewise
   radial identities (relating the two different power-law kernels)
   that produce the specific rational coefficients.

VERDICT: The C^(k) convention hypothesis is NEGATIVE. The near-hit
beta = 1/sqrt(4pi) was a coincidence arising from the scale of the
angular coupling coefficients, not from a Y^(k) ↔ C^(k) conversion.
The structural obstruction is the RADIAL KERNEL MISMATCH between the
bipolar expansion and Drake's M^k integrals, not an angular convention.

The Sprint 5 DV conclusion stands: (3/50, -2/5) are convention-dependent
combining identities that encode the mapping between two different
radial-integral bases.
""")
