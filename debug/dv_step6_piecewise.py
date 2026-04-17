"""Sprint 5 DV Step 6: Use the production module's piecewise region
integrals to decompose bipolar channels into Drake M^K basis.

The production module `geovac.breit_integrals` exposes `_t_kernel_region_I`
which returns the Region I (r_1 < r_2) piece of the Drake M^K integral.
Region II is Region I with electron labels swapped.

For the bipolar channel (k_1, k_2, K=2), the kernel is
r_1^{k_1} r_2^{k_2} / r_>^5. In:
  Region I (r_1 < r_2, r_> = r_2): kernel = r_1^{k_1} r_2^{k_2 - 5}
  Region II (r_1 > r_2, r_> = r_1): kernel = r_1^{k_1 - 5} r_2^{k_2}

For Drake's M^K (standard Breit-Pauli r_<^K / r_>^{K+3}):
  Region I: r_1^K / r_2^{K+3}
  Region II: r_1^{-(K+3)} r_2^K

Matching bipolar (k_1, k_2) K=2 Region I to Drake K_d Region I requires:
  k_1 = K_d, k_2 - 5 = -(K_d + 3)
  => K_d = k_1 and K_d = k_2 - 5 + 3 = k_2 - 2
  so k_1 + 2 = k_2, or k_2 = k_1 + 2. (k_1, k_2) = (0,2), (1,3)... BUT
  we ALSO need triangle (k_1, k_2, K=2) which restricts k_1 ∈ {0,1,2} with
  k_2 = 2 - k_1 for the LEADING term.

OK so (k_1 + k_2 = K = 2, k_1 = k_2 - 2) with k_1 + k_2 = 2 gives
2 k_1 = 0, so k_1 = 0, k_2 = 2.
  Bipolar(0, 2) Region I = Drake M^{K_d=0} Region I ✓

For bipolar(1, 1) Region I: kernel = r_1 / r_2^4 = Drake K_d=1 Region I.
For bipolar(2, 0) Region I: kernel = r_1^2 / r_2^5 = Drake K_d=2 Region I.

Summary of the channel-to-Drake-piece map:
  Bipolar(0, 2) dir = Drake M^0 dir_I  + Drake M^2 dir_II
  Bipolar(1, 1) dir = Drake M^1 dir_I  + Drake M^1 dir_II = Drake M^1 dir
  Bipolar(2, 0) dir = Drake M^2 dir_I  + Drake M^0 dir_II
  Bipolar(0, 2) exch = Drake M^0 exch_I + Drake M^2 exch_II
  Bipolar(1, 1) exch = Drake M^1 exch_I + Drake M^1 exch_II = Drake M^1 exch
  Bipolar(2, 0) exch = Drake M^2 exch_I + Drake M^0 exch_II

For exchange density (r_1 <-> r_2 symmetric):
  M^K_exch_II = M^K_exch_I (by relabeling symmetry)
  So M^K_exch = M^K_exch_I + M^K_exch_II = 2 * M^K_exch_I
  Bipolar(0, 2) exch = M^0_exch/2 + M^2_exch/2
  Bipolar(1, 1) exch = M^1_exch
  Bipolar(2, 0) exch = M^2_exch/2 + M^0_exch/2
"""
from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
from pathlib import Path
from fractions import Fraction

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, symbols, expand, collect, together, nsimplify

from geovac.breit_integrals import (
    breit_ss_radial, _t_kernel_region_I, _expand_product, _fraction_sqrt,
    _simplify_log_form,
)


def production_M_K_piecewise(n1, l1, n3, l3, n2, l2, n4, l4, k):
    """Return (M^k Region I, M^k Region II) for the specified radial integral.

    Uses the production module's _t_kernel_region_I directly.
    """
    terms13, nsq13 = _expand_product(n1, l1, n3, l3)
    terms24, nsq24 = _expand_product(n2, l2, n4, l4)
    norm = _fraction_sqrt(nsq13 * nsq24)
    norm_s = Rational(norm.numerator, norm.denominator)

    # Region I
    integral_I = Rational(0)
    for c13, p13, a13 in terms13:
        for c24, p24, a24 in terms24:
            if c13 == 0 or c24 == 0:
                continue
            integral_I = integral_I + c13 * c24 * _t_kernel_region_I(
                p13, p24, a13, a24, "breit", k
            )
    I = simplify(norm_s * integral_I)

    # Region II: swap (p13, a13) <-> (p24, a24)
    integral_II = Rational(0)
    for c13, p13, a13 in terms13:
        for c24, p24, a24 in terms24:
            if c13 == 0 or c24 == 0:
                continue
            integral_II = integral_II + c13 * c24 * _t_kernel_region_I(
                p24, p13, a24, a13, "breit", k
            )
    II = simplify(norm_s * integral_II)
    return I, II


print("=" * 76)
print("Sprint 5 DV Step 6: Bipolar channel integrals via piecewise Drake M^K")
print("=" * 76)

# ============================================================
# Compute Drake M^K pieces (K = 0, 1, 2)
# ============================================================
pieces = {}
for K_d in (0, 1, 2):
    I_d, II_d = production_M_K_piecewise(1, 0, 2, 1, 1, 0, 2, 1, K_d)
    I_e, II_e = production_M_K_piecewise(1, 0, 2, 1, 2, 1, 1, 0, K_d)
    pieces[(K_d, 'dir', 'I')] = _simplify_log_form(I_d)
    pieces[(K_d, 'dir', 'II')] = _simplify_log_form(II_d)
    pieces[(K_d, 'exch', 'I')] = _simplify_log_form(I_e)
    pieces[(K_d, 'exch', 'II')] = _simplify_log_form(II_e)

    # Total and cross-check
    total_d = simplify(I_d + II_d)
    total_e = simplify(I_e + II_e)
    prod_d = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, K_d, Z=1)
    prod_e = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, K_d, Z=1)
    print(f"\n  M^{K_d}_dir: I = {float(I_d):+.6e}, II = {float(II_d):+.6e}, "
          f"total = {float(total_d):+.6e} (vs production: {float(prod_d):+.6e}) "
          f"diff={float(simplify(total_d - prod_d)):+.2e}")
    print(f"  M^{K_d}_exch: I = {float(I_e):+.6e}, II = {float(II_e):+.6e}, "
          f"total = {float(total_e):+.6e} (vs production: {float(prod_e):+.6e}) "
          f"diff={float(simplify(total_e - prod_e)):+.2e}")
    # Check exchange symmetry: I_e should equal II_e
    diff_exch_sym = simplify(I_e - II_e)
    if diff_exch_sym == 0:
        print(f"    [exchange r_1<->r_2 symmetry confirmed: I = II]")
    else:
        print(f"    [asymmetry: I - II = {diff_exch_sym}]")

# ============================================================
# Now: Bipolar(k_1, k_2, K=2) channel integrals
# ============================================================
#
# For channel (k_1, k_2, K=2), the kernel is r_1^{k_1} r_2^{k_2} / r_>^5.
#
# Region I (r_1 < r_2, r_>=r_2): kernel = r_1^{k_1} r_2^{k_2 - 5}
# Region II (r_1 > r_2, r_>=r_1): kernel = r_1^{k_1 - 5} r_2^{k_2}
#
# We compute this directly using _t_kernel_region_I with modified (a, b, k_pow, kernel_type).
# Actually, the production module's kernel is r_<^k / r_>^{k+3}, so we need to
# fake it. The simplest: compute as a linear combination of known Drake pieces.
#
# Matching:
#   Bipolar Region I kernel r_1^{k_1} r_2^{k_2 - 5} = r_1^K / r_2^{K+3}
#     ==> K = k_1 and K+3 = 5 - k_2, i.e., K = 2 - k_2 and K = k_1 so k_1 + k_2 = 2. ✓
#
# So for each (k_1, k_2) with k_1 + k_2 = 2 (which is K):
#   Region I of bipolar(k_1, k_2, K=2) = Region I of Drake M^{k_1}
#   Region II of bipolar(k_1, k_2, K=2) = Region II of Drake M^{k_2}

print("\n\n=== Bipolar channel integrals (K=2, k_1 + k_2 = 2) ===")
bipolar_in_drake = {}
for k1, k2 in [(0, 2), (1, 1), (2, 0)]:
    # Direct channel
    dir_bipolar = pieces[(k1, 'dir', 'I')] + pieces[(k2, 'dir', 'II')]
    exch_bipolar = pieces[(k1, 'exch', 'I')] + pieces[(k2, 'exch', 'II')]
    bipolar_in_drake[(k1, k2, 'dir')] = simplify(dir_bipolar)
    bipolar_in_drake[(k1, k2, 'exch')] = simplify(exch_bipolar)
    print(f"\n  Bipolar(k_1={k1}, k_2={k2}, K=2):")
    print(f"    Direct  = M^{k1}_dir_I + M^{k2}_dir_II")
    print(f"            = {float(dir_bipolar):+.6e}")
    print(f"    Exchange = M^{k1}_exch_I + M^{k2}_exch_II")
    print(f"            = {float(exch_bipolar):+.6e}")

# ============================================================
# Step 7: Assemble the full <^3P_1 | H_SS | ^3P_1> = A_SS matrix element
# and extract (3/50, -2/5) from the bipolar decomposition.
# ============================================================
#
# <^3P_1 | H_SS | ^3P_1> = A_SS
#
# From Racah decomposition:
#   A_SS = sum over (k_1, k_2) channels of (direct - exchange)
#          * (angular 9j factor)
#          * (spin reduced ME factor)
#          * (bipolar radial integral for that (k_1, k_2))
#
# Using scalar product decomposition:
#   <[T^(K)(space) . U^(K)(spin)]^(0)>_{J=L=S=1} = (-1)^{L+S+J} 6j{LSJ;SLK}/sqrt(2J+1)
#     * <L||T_space||L> <S||U_spin||S>
#
# With spatial reduced ME expanded bipolarly:
#   <L=1 || T_space^(K) || L=1>_{channel (k1,k2)}_dir or exch
#     = sqrt((2L+1)(2L+1)(2K+1)) * 9j{l_a l_b L; l_c l_d L; k1 k2 K}
#       * <l_a||C^(k1)||l_c> <l_b||C^(k2)||l_d>
#       * [bipolar radial integral]
#
# Triplet spatial wavefunction:
#   |triplet> = (1/sqrt 2) [|1s 2p> - |2p 1s>] coupled to L=1
# Matrix element of scalar operator V:
#   <triplet | V | triplet> = (1/2) [<1s 2p|V|1s 2p> - <1s 2p|V|2p 1s>
#                                    - <2p 1s|V|1s 2p> + <2p 1s|V|2p 1s>]
# For scalar V that's Hermitian and spin-independent:
#   <1s 2p|V|1s 2p> = <2p 1s|V|2p 1s>  (by relabeling)
#   <1s 2p|V|2p 1s> = <2p 1s|V|1s 2p>
# So <triplet | V | triplet> = <1s 2p|V|1s 2p> - <1s 2p|V|2p 1s> = direct - exchange
# which matches the standard rule.

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, symbols, nsimplify
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan


def red_C(l, lp, k):
    """Reduced matrix element <l || C^(k) || l'> (Racah normalization)."""
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * wigner_3j(l, k, lp, 0, 0, 0)


# Spin reduced ME:
# <S=1 || [s_1 (x) s_2]^(2) || S=1> = sqrt(5)/2
S = 1
k_spin = 2
n9_spin = wigner_9j(Rational(1, 2), Rational(1, 2), S, 1, 1, k_spin, Rational(1, 2), Rational(1, 2), S)
red_spin = sqrt(Integer((2 * S + 1) * (2 * S + 1) * (2 * k_spin + 1))) * n9_spin * Rational(3, 2)
red_spin = simplify(red_spin)
print(f"\n<S=1 || [s_1 (x) s_2]^(2) || S=1> = {red_spin}")

# ============================================================
# Now compute direct and exchange contributions to A_SS
# ============================================================
K = 2
la, lb = 0, 1

# 6j for J-pattern:
L_val, S_val, J_val = 1, 1, 1
phase = (-1) ** (L_val + S_val + J_val)
j6 = wigner_6j(L_val, S_val, J_val, S_val, L_val, K)
J_factor = phase * j6 / sqrt(Integer(2 * J_val + 1))
J_factor = simplify(J_factor)
print(f"\nJ-factor (phase * 6j / sqrt(2J+1)): {J_factor}")

# Hamiltonian prefactor: H_SS = -alpha^2 sqrt(24 pi / 5) * [s1 x s2]^(2) . Y^(2)(r_hat_12)/r_12^3
# With Y^K = sqrt((2K+1)/(4pi)) C^K:
#   Y^(2)_Q = sqrt(5/(4pi)) C^(2)_Q
# So H_SS = -alpha^2 sqrt(24 pi / 5) * sqrt(5/(4pi)) * [s1 x s2]^(2) . C^(2)(r_hat_12)/r_12^3
#         = -alpha^2 * sqrt(6) * [s1 x s2]^(2) . C^(2)(r_hat_12)/r_12^3
H_SS_prefactor = -sqrt(6)  # alpha^2 pulled out

# Bipolar expansion of C^(2)(hat r_12) / r_12^3 in the "leading" form
# (k_1 + k_2 = K = 2):
#   C^(2)_Q(hat r_12) / r_12^3 = sum_{k_1+k_2=2} beta(k_1, k_2, 2)
#                                  [C^(k_1)(1) (x) C^(k_2)(2)]^(2)_Q
#                                  * r_1^{k_1} r_2^{k_2} / r_>^5
# where beta(k_1, k_2, K) is the bipolar amplitude.
#
# The beta coefficient (from Sack 1964 / Brink-Satchler 4.6.13, adapted to C^k):
#   beta(k_1, k_2, K) = sqrt[ (2K+1)! / ((2k_1)! (2k_2)!) ] * (-1)^? * (1 / sqrt(2K+1))
#
# Let's DERIVE it by matching to the Laplace expansion K=0:
#   1/r_12 = sum_l (r_<^l / r_>^{l+1}) P_l(cos theta_12)
#          = sum_l (r_<^l / r_>^{l+1}) C^l(1) . C^l(2)
#          = sum_l (r_<^l / r_>^{l+1}) sum_m (-1)^m C^l_m(1) C^l_{-m}(2)
#
# With scalar product of rank-0 tensors:
#   [C^l(1) (x) C^l(2)]^(0)_0 = -1/sqrt(2l+1) * C^l(1) . C^l(2)
# So
#   1/r_12 = -sum_l (r_<^l / r_>^{l+1}) sqrt(2l+1) [C^l(1) (x) C^l(2)]^(0)_0
#
# But 1/r_12 = Y^0(hat r_12) / r_12^1 * sqrt(4pi) = C^0(hat r_12) / r_12 (trivial).
# For the target form C^K(hat r_12) / r_12^{K+1}, we would have to decompose into
# [C^{k_1}(1) x C^{k_2}(2)]^K rather than ^0.

# Let me use the "concrete" approach: just compute the ANGULAR reduced ME
# for each (k_1, k_2) channel and TRUST that the radial integral is the
# one I derived above. The OVERALL NORMALIZATION (beta prefactor) will
# come out from matching to a known K=0 identity.

# ============================================================
# Spatial reduced ME (direct channel):
#   <(0,1) L=1 || T^(K=2)_space |(direct)| (0,1) L=1>
# Expanded bipolarly as sum over (k_1, k_2) with ONE contribution:
#   (k_1=0, k_2=2) since <0||C^0||0>=1 only for k_1=0, and <1||C^2||1> for k_2 even, triangle with K=2
# ============================================================

# For DIRECT: (l_a, l_b, l_c, l_d) = (0, 1, 0, 1)
# Only contributing (k_1, k_2): (0, 2)
k1, k2 = 0, 2
red1_dir = red_C(0, 0, k1)   # 1
red2_dir = red_C(1, 1, k2)   # <1||C^2||1> = -sqrt(30)/5
n9_dir = wigner_9j(0, 1, 1, 0, 1, 1, k1, k2, K)
norm_L = sqrt(Integer((2 * 1 + 1) * (2 * 1 + 1) * (2 * K + 1)))  # sqrt(3*3*5) = sqrt(45) = 3 sqrt(5)
spatial_red_dir = norm_L * n9_dir * red1_dir * red2_dir
spatial_red_dir = simplify(spatial_red_dir)
print(f"\nDirect spatial reduced ME (k_1=0, k_2=2) [no radial factor]: {spatial_red_dir}")
# Expected: sqrt((2L+1)^2(2K+1)) = sqrt(45) = 3 sqrt(5)
#   * 9j (0,1,1; 0,1,1; 0,2,2) = sqrt(5)/15 (from Step 1)
#   * red_1 = 1
#   * red_2 = -sqrt(30)/5
# = 3 sqrt(5) * sqrt(5)/15 * (-sqrt(30)/5) = 3*5/15 * (-sqrt(30)/5) = 1 * (-sqrt(30)/5) = -sqrt(30)/5

# For EXCHANGE: (l_a, l_b, l_c, l_d) = (0, 1, 1, 0)
# Only contributing (k_1, k_2): (1, 1)
k1, k2 = 1, 1
red1_exch = red_C(0, 1, k1)  # <0||C^1||1> = -1
red2_exch = red_C(1, 0, k2)  # <1||C^1||0> = 1
n9_exch = wigner_9j(0, 1, 1, 1, 0, 1, k1, k2, K)
spatial_red_exch = norm_L * n9_exch * red1_exch * red2_exch
spatial_red_exch = simplify(spatial_red_exch)
print(f"Exchange spatial reduced ME (k_1=1, k_2=1) [no radial factor]: {spatial_red_exch}")

# ============================================================
# Assemble: A_SS (amplitude coefficient, no alpha^2)
# ============================================================
# A_SS = <^3P_1 | H_SS / alpha^2 | ^3P_1>
#      = (-sqrt 6) * [sum over channels ] * J_factor * red_spin
#      = (-sqrt 6) * (direct - exchange) * J_factor * red_spin
#
# where "direct - exchange" includes:
#   direct: spatial_red_dir(k_1=0, k_2=2) * beta(0, 2, 2) * I_bipolar(0,2)_dir
#   exchange: spatial_red_exch(k_1=1, k_2=1) * beta(1, 1, 2) * I_bipolar(1,1)_exch
#
# beta(k_1, k_2, K) = UNKNOWN; we'll calibrate it by matching to the
# Coulomb expansion identity.

# For now: express A_SS as a LINEAR combination of Drake pieces with
# unknown bipolar prefactors. The linearization:
#
# A_SS = CONST * (
#   spatial_red_dir * beta(0,2,2) * (M^0_dir_I + M^2_dir_II)
#   - spatial_red_exch * beta(1,1,2) * (M^1_exch_I + M^1_exch_II)
# )
#
# where CONST = (-sqrt 6) * J_factor * red_spin.
# Drake's form: A_SS = 3/50 M^2_dir - 2/5 M^2_exch
#                    = 3/50 (M^2_dir_I + M^2_dir_II) - 2/5 (M^2_exch_I + M^2_exch_II)
#
# These can match if and only if:
#   Bipolar-derived coefficient of M^2_dir_I = 0 (not present in Drake's M^2_dir_I)
#   Wait -- Drake's M^2_dir = M^2_dir_I + M^2_dir_II has coef 3/50 on TOTAL.
#   But bipolar derived expression has coef on M^0_dir_I (from (0,2) channel) and
#   M^2_dir_II (from (0,2) channel) - these are DIFFERENT pieces.
#
# So the Drake formula (in terms of M^K totals) CAN ONLY be correct if other
# terms in the bipolar expansion exist that reconstruct the Drake combination.
#
# In particular, we also need the (2, 0) channel:
#   Bipolar(2, 0, K=2) dir = M^2_dir_I + M^0_dir_II
#   Bipolar(2, 0, K=2) exch = M^2_exch_I + M^0_exch_II
#
# But (k_1=2, k_2=0) contributes to which 9j? For DIRECT (0,1,0,1),
# we need <0||C^2||0> which is ZERO (parity mismatch). So Bipolar(2, 0)
# is FORBIDDEN on the direct channel.
#
# The only contributing direct (k_1, k_2) is (0, 2).
# The only contributing exchange (k_1, k_2) is (1, 1).
#
# So the answer to A_SS from bipolar is:
#   A_SS = CONST * [C_d * (M^0_dir_I + M^2_dir_II) - C_e * M^1_exch]
# where C_d = spatial_red_dir * beta(0,2,2), C_e = spatial_red_exch * beta(1,1,2).
#
# For this to equal 3/50 M^2_dir - 2/5 M^2_exch, we need the pieces to rearrange.
# But M^0_dir_I + M^2_dir_II is NOT a simple rational multiple of M^2_dir.
#
# THIS SUGGESTS OUR BIPOLAR EXPANSION IS MISSING TERMS.

# ============================================================
# The resolution: the bipolar expansion of Y^K(r_hat_12) / r_12^{K+1}
# has MORE terms than just k_1 + k_2 = K. In particular, the HIGHER
# rank contributions from r_<^m r_>^n for (m, n) != (K, -K-1) form
# additional channels.
#
# Let me consult the CORRECT Sack/Wen formula.
# ============================================================

# === Sack 1964 formula (Appendix) ===
# f(r_12) has bipolar expansion:
# f(r_12) = sum_{L} f_L(r_1, r_2) P_L(cos theta_12)
# where f_L(r_1, r_2) = (2L+1)/2 integrate[ f(r_12) P_L(cos theta) sin theta, {theta, 0, pi}]
#
# This gives the RADIAL projections of f onto P_L. These are Legendre-expansion
# coefficients.
#
# For f(r_12) = 1/r_12:
#   f_L(r_1, r_2) = r_<^L / r_>^{L+1}  (Laplace expansion)
#
# For f(r_12) = 1/r_12^3 (our Breit-Pauli r_12 denominator):
#   f_L(r_1, r_2) = ?
#
# Compute numerically: f(r_12) = 1/r_12^3 where r_12 = sqrt(r_1^2 + r_2^2 - 2 r_1 r_2 cos theta).
# Legendre project.

# === BRINK-SATCHLER 4.6.13 EXACT FORM (for Y^K(r_hat_12)/r_12^{K+1}) ===
#
# I have to be more careful. Let me consult a verified formula.
#
# The Fontana-Hansen expansion (see Kay & Silverstone 1969, JCP 50, 3137):
#
# Y^K_M(hat r_12) / r_12^{K+1} = (4pi)^{-1/2} * sum_{l1,l2}
#   [ (2l1+1)(2l2+1) / (2K+1) ]^{1/2} * <l1 0 l2 0 | K 0> * Q(l1, l2, K; r_1, r_2)
#   * [Y^{l1}(hat r_1) (x) Y^{l2}(hat r_2)]^K_M
#
# where Q(l1, l2, K; r_1, r_2) contains the (r_1, r_2) dependence.
# The form of Q depends on the triangle (l1, l2, K).
#
# For the case K = l_1 + l_2 (leading term):
#   Q(l_1, l_2, K=l_1+l_2; r_1, r_2) = r_1^{l_1} r_2^{l_2} / r_>^{2K+1}
#     * (2K)! / ((2l_1)! (2l_2)!)^{1/2}    [Brink-Satchler 4.6.13]
#
# For lower-rank (l_1 + l_2 = K - 2, K - 4, ...), there are additional
# terms from the Sack expansion. But these ADDITIONAL terms involve higher
# derivatives and are typically zero for the specific (l_1, l_2, K) where
# triangle is satisfied by BOTH (l_1, K, l_1') and (l_2, K, l_2').
#
# For our He problem with (l_a, l_b, l_c, l_d) = (0, 1, 0, 1) direct or
# (0, 1, 1, 0) exchange, the allowed (l_1, l_2) on the direct side are
# ONLY (0, 2) [NOT (0, 0) because K=2 triangle forbidden, NOT (2, 2) because
# <l=0||C^2||0>=0]. On exchange side, only (1, 1).
#
# In both cases, l_1 + l_2 = K = 2, so we are in the LEADING term regime
# and Q = r_<^{l_1} r_>^{l_2} / r_>^{2K+1}... wait not exactly, Q is r_1^{l_1} r_2^{l_2} / r_>^{2K+1}.

# Given this convention, the beta coefficient is:
#   beta(l_1, l_2, K=l_1+l_2) = [(2K)! / ((2l_1)! (2l_2)!)]^{1/2} / (4pi)^{1/2}
#                               * [(2l_1+1)(2l_2+1)/(2K+1)]^{1/2}
#
# For (l_1=0, l_2=2, K=2):
#   beta(0, 2, 2) = sqrt(4! / (0! * 4!)) / sqrt(4pi) * sqrt(1*5/5) = 1/sqrt(4pi)
#                 * 1 = 1/sqrt(4pi)
#
# For (l_1=1, l_2=1, K=2):
#   beta(1, 1, 2) = sqrt(4! / (2! * 2!)) / sqrt(4pi) * sqrt(3*3/5)
#                 = sqrt(6) / sqrt(4pi) * 3/sqrt(5)
#                 = 3 sqrt(6) / sqrt(20 pi)
#                 = 3 sqrt(6/5) / (2 sqrt(pi))
#   numerical: 3 * sqrt(1.2) / (2 * 1.7725) = 3.286 / 3.545 = 0.927

# But these are amplitudes for Y^(L1) x Y^(L2). In Racah C^k convention:
#   Y^L = sqrt((2L+1)/4pi) C^L
#   C^L = sqrt(4pi/(2L+1)) Y^L
# So the bipolar [Y^{l1}(1) x Y^{l2}(2)]^K is related to [C^{l1}(1) x C^{l2}(2)]^K by
#   [Y^{l1} x Y^{l2}]^K = sqrt((2l_1+1)(2l_2+1)/(4pi)^2) [C^{l1} x C^{l2}]^K

# Converting beta from the Y^L form to C^L form:
#   beta_C = beta_Y * sqrt((2l_1+1)(2l_2+1)/(4pi)^2)   [wait, this is in the OTHER direction]
#
# Let me just write out the full equation in Y^L form and then re-express.

# ============================================================
# INSTEAD: instead of chasing convention, let me ASSUME the  form
#   C^(K)(hat r_12) / r_12^{K+1} = sum_{k_1 + k_2 = K} alpha(k_1, k_2, K)
#       [C^{k_1}(1) (x) C^{k_2}(2)]^K_Q * r_1^{k_1} r_2^{k_2} / r_>^{2K+1}
# with alpha(k_1, k_2, K) to be determined, and calibrate alpha by
# matching to the Coulomb expansion K = 0.
# ============================================================
print("\n\n=== Calibrate bipolar coefficient via K=0 Coulomb identity ===")
# For K=0, the Laplace expansion gives
#   1/r_12 = sum_l (r_<^l / r_>^{l+1}) C^l(1) . C^l(2)
#          = sum_l (r_<^l / r_>^{l+1}) sum_m (-1)^m C^l_m(1) C^l_{-m}(2)
#
# In [C^l(1) x C^l(2)]^(0) form: since [C^l(1) x C^l(2)]^(0)_0 = 1/sqrt(2l+1) *
#   sum_m C^l_m(1) C^l_{-m}(2) (-1)^(l-m)... careful about conventions
#
# The scalar product: T^(k).U^(k) = sum_q (-1)^q T^k_q U^k_{-q} = -(2k+1)^{1/2} [T x U]^(0)_0
# [T x U]^(0)_0 = 1/(sqrt(2k+1)) sum_q (-1)^(k-q) T^k_q U^k_{-q} ... not standard anyway.
#
# Standard: [T^k(1) x U^k(2)]^(0) = sum_q <k q k -q | 0 0> T^k_q(1) U^k_{-q}(2)
#   <k q k -q | 0 0> = (-1)^(k-q) / sqrt(2k+1)
# So [T x U]^(0)_0 = (-1)^{-k} / sqrt(2k+1) sum_q (-1)^q T^k_q U^k_{-q}
#                  = (-1)^k / sqrt(2k+1) T.U   (using sum_q (-1)^q T^k_q U^k_{-q} = T.U)
# Equivalently, T.U = (-1)^k sqrt(2k+1) [T x U]^(0)_0

# Thus
#   1/r_12 = sum_l (r_<^l / r_>^{l+1}) C^l(1) . C^l(2)
#          = sum_l (r_<^l / r_>^{l+1}) (-1)^l sqrt(2l+1) [C^l(1) x C^l(2)]^(0)_0

# Is this of the form sum_{k_1+k_2=K=0} alpha [C^{k_1}(1) x C^{k_2}(2)]^0 * r_1^{k_1} r_2^{k_2} / r_>^{2*0+1=1}?
#
# The "k_1 + k_2 = K = 0" restriction forces k_1 = k_2 = 0.
# But our Laplace expansion has ALL l. These are different expansions!
#
# CONCLUSION: the bipolar expansion of Y^K(r_hat_12)/r_12^{K+1} is NOT
# restricted to k_1 + k_2 = K. It includes all (k_1, k_2) with triangle
# (k_1, k_2, K) AND parity (k_1 + k_2 + K) even.
#
# The leading "k_1 + k_2 = K" is only ONE term; for K=0 the full expansion
# has ALL l (with k_1 = l = k_2, since (l, l, 0) triangle).
#
# Let me write the CORRECT expansion for general K.

# THE CORRECT FORMULA (I'll derive it by numerical fitting):
# Y^K_M(hat r_12) / r_12^{K+1} = sum_{k_1, k_2 triangle (k_1, k_2, K),
#                                    k_1+k_2+K even}
#                                  alpha(k_1, k_2, K) [Y^{k_1}(1) x Y^{k_2}(2)]^K_M
#                                  * g(k_1, k_2, K; r_1, r_2)
#
# with alpha and g to be found. For K=0 we know:
# 1/r_12 = sum_l (r_<^l / r_>^{l+1}) (-1)^l sqrt(2l+1) [C^l(1) x C^l(2)]^(0)_0

print("\nBipolar expansion for general K is more complex than 'k_1+k_2=K'")
print("It requires a Sack-Wen-Steinborn-Filter style derivation, not just leading term.")
