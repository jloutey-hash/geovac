"""Sprint 5 DV Step 11: Pure-float computation, no symbolic printing of M^K.

Uses production integrals (floats via float(breit_ss_radial(..., Z=1))),
piecewise extracted via _t_kernel_region_I, and scipy for bipolar integral
verification. Avoids sympy's huge-rational printing problem.
"""
from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
sys.set_int_max_str_digits(10**7)
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fractions import Fraction
import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, pi, Symbol
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

from geovac.breit_integrals import breit_ss_radial, _t_kernel_region_I, _expand_product, _fraction_sqrt


def piecewise_BP_float(n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, k):
    """Region-split Drake M^K integral. Returns (I, II) FLOAT tuple."""
    terms_e1, nsq_e1 = _expand_product(n_a, l_a, n_c, l_c)
    terms_e2, nsq_e2 = _expand_product(n_b, l_b, n_d, l_d)
    norm = _fraction_sqrt(nsq_e1 * nsq_e2)
    norm_s = Rational(norm.numerator, norm.denominator)

    I_I = Rational(0)
    I_II = Rational(0)
    for c13, p13, a13 in terms_e1:
        for c24, p24, a24 in terms_e2:
            if c13 == 0 or c24 == 0:
                continue
            I_I = I_I + c13 * c24 * _t_kernel_region_I(p13, p24, a13, a24, "breit", k)
            I_II = I_II + c13 * c24 * _t_kernel_region_I(p24, p13, a24, a13, "breit", k)
    # AVOID expensive simplify; use direct float eval.
    return float(norm_s * I_I), float(norm_s * I_II)


# Compute all Drake M^K pieces as floats at Z=1
pieces = {}
for K_d in (0, 1, 2):
    pieces[(K_d, 'dir', 'I')], pieces[(K_d, 'dir', 'II')] = piecewise_BP_float(
        1, 0, 2, 1, 1, 0, 2, 1, K_d
    )
    pieces[(K_d, 'exch', 'I')], pieces[(K_d, 'exch', 'II')] = piecewise_BP_float(
        1, 0, 2, 1, 2, 1, 1, 0, K_d
    )

# Print piecewise (float) values
print("Drake M^K piecewise (Z=1, floats):")
for K_d in (0, 1, 2):
    print(f"  M^{K_d}_dir:  I = {pieces[(K_d, 'dir', 'I')]:+.6e}, II = {pieces[(K_d, 'dir', 'II')]:+.6e}, "
          f"total = {pieces[(K_d, 'dir', 'I')] + pieces[(K_d, 'dir', 'II')]:+.6e}")
    print(f"  M^{K_d}_exch: I = {pieces[(K_d, 'exch', 'I')]:+.6e}, II = {pieces[(K_d, 'exch', 'II')]:+.6e}, "
          f"total = {pieces[(K_d, 'exch', 'I')] + pieces[(K_d, 'exch', 'II')]:+.6e}")

# Bipolar channel integrals
# Bipolar(k_1, k_2, K=2) dir = Drake M^{k_1} Region I + Drake M^{k_2} Region II
I_bipolar_02_dir = pieces[(0, 'dir', 'I')] + pieces[(2, 'dir', 'II')]
I_bipolar_11_dir = pieces[(1, 'dir', 'I')] + pieces[(1, 'dir', 'II')]
I_bipolar_20_dir = pieces[(2, 'dir', 'I')] + pieces[(0, 'dir', 'II')]
I_bipolar_02_exch = pieces[(0, 'exch', 'I')] + pieces[(2, 'exch', 'II')]
I_bipolar_11_exch = pieces[(1, 'exch', 'I')] + pieces[(1, 'exch', 'II')]
I_bipolar_20_exch = pieces[(2, 'exch', 'I')] + pieces[(0, 'exch', 'II')]

print(f"\nBipolar (K=2) channel integrals (floats):")
print(f"  Bipolar(0,2) dir:  {I_bipolar_02_dir:+.6e}")
print(f"  Bipolar(0,2) exch: {I_bipolar_02_exch:+.6e}")
print(f"  Bipolar(1,1) dir:  {I_bipolar_11_dir:+.6e}")
print(f"  Bipolar(1,1) exch: {I_bipolar_11_exch:+.6e}")
print(f"  Bipolar(2,0) dir:  {I_bipolar_20_dir:+.6e}")
print(f"  Bipolar(2,0) exch: {I_bipolar_20_exch:+.6e}")

# Drake's A_SS target:
M2_dir = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1))
M2_exch = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1))
A_SS_target = 3 / 50 * M2_dir - 2 / 5 * M2_exch
print(f"\nDrake A_SS target (Z=1, no alpha^2):")
print(f"  3/50 * M^2_dir = {3/50 * M2_dir:+.6e}")
print(f"  -2/5 * M^2_exch = {-2/5 * M2_exch:+.6e}")
print(f"  Sum = {A_SS_target:+.6e}")

# ============================================================
# Now apply the Racah formula with unknown (beta_02, beta_11)
# ============================================================
# <^3P_1, 1 | H_SS | ^3P_1, 1> = A_SS (with f_SS(1) = 1)
# Racah: = (-sqrt 6) * (-1/6) * (sqrt 5 / 2)
#          * [ sp_red_dir * beta_02 * I_bipolar_02_dir - sp_red_exch * beta_11 * I_bipolar_11_exch ]
#        = scaling * [ a * beta_02 + b * beta_11 ]
# where
sp_red_dir_f = -float(sqrt(30) / 5)
sp_red_exch_f = -float(sqrt(5) / 3)
spin_red_f = float(sqrt(5) / 2)
H_prefactor_f = -float(sqrt(6))
J_factor_f = -1 / 6

scaling = H_prefactor_f * J_factor_f * spin_red_f
a_coef = sp_red_dir_f * I_bipolar_02_dir
b_coef = -sp_red_exch_f * I_bipolar_11_exch

print(f"\nA_SS = scaling * (beta_02 * a + beta_11 * b)")
print(f"  scaling = {scaling:+.6e}")
print(f"  a = sp_red_dir * I_bipolar_02_dir = {a_coef:+.6e}")
print(f"  b = -sp_red_exch * I_bipolar_11_exch = {b_coef:+.6e}")

bracket_target = A_SS_target / scaling
print(f"\nbracket target = {bracket_target:+.6e}")
print(f"Need: beta_02 * a + beta_11 * b = bracket_target")
print(f"      beta_02 * {a_coef:+.6e} + beta_11 * {b_coef:+.6e} = {bracket_target:+.6e}")

# ============================================================
# Test various beta values
# ============================================================
print("\n=== Test specific (beta_02, beta_11) combinations ===")
import math
beta_candidates = {
    "1 (unit)": 1.0,
    "sqrt(2)": math.sqrt(2),
    "sqrt(3)": math.sqrt(3),
    "sqrt(5)": math.sqrt(5),
    "sqrt(6)": math.sqrt(6),
    "sqrt(3/5)": math.sqrt(3/5),
    "sqrt(5/3)": math.sqrt(5/3),
    "sqrt(6/5)": math.sqrt(6/5),
    "sqrt(1/5)": math.sqrt(1/5),
    "sqrt(1/20)": math.sqrt(1/20),
    "1/sqrt(30)": 1/math.sqrt(30),
    "sqrt(3/50)": math.sqrt(3/50),
    "sqrt(1/10)": math.sqrt(1/10),
    "1/2": 0.5,
    "1/sqrt(4pi)": 1/math.sqrt(4*math.pi),
    "sqrt(5)/2": math.sqrt(5)/2,
}

# Search for (beta_02, beta_11) such that A_SS_test = A_SS_target
print(f"\nSearching for beta values with ratio close to 1...")
best = None
for name02, v02 in beta_candidates.items():
    for name11, v11 in beta_candidates.items():
        A_test = scaling * (v02 * a_coef + v11 * b_coef)
        ratio = A_test / A_SS_target
        if abs(ratio - 1.0) < 0.02:
            print(f"  beta_02={name02} (={v02:.4f}), beta_11={name11} (={v11:.4f}): ratio={ratio:.6f}")

# ============================================================
# Alternative: solve the linear system with M^2_dir, M^2_exch as constraints
# ============================================================
print("\n=== ALTERNATIVE: express A_SS_derived as linear combination of M^K ===")
# We have:
#   A_SS_derived = scaling * beta_02 * sp_red_dir * I_bipolar_02_dir
#                + scaling * (-beta_11) * sp_red_exch * I_bipolar_11_exch
#   = scaling * beta_02 * sp_red_dir * [M^0_dir_I + M^2_dir_II]
#    + scaling * (-beta_11) * sp_red_exch * [M^1_exch]   (since I = II in exch)

# In the Drake basis:
#   A_SS = c_0d * M^0_dir_Region_I + c_2d_II * M^2_dir_Region_II + c_1e * M^1_exch_total
# where c_0d (coef of M^0_dir_I) = scaling * beta_02 * sp_red_dir
#       c_2d_II (coef of M^2_dir_II) = same as c_0d (both in direct channel)
#       c_1e (coef of M^1_exch) = scaling * (-beta_11) * (-sp_red_exch)

# Now we want this to equal 3/50 M^2_dir - 2/5 M^2_exch = 3/50 (M^2_dir_I + M^2_dir_II) - 2/5 (M^2_exch_I + M^2_exch_II)
# That requires:
#   coefficient of M^0_dir_I: we have c_0d, Drake has 0. So c_0d = 0 UNLESS there's
#     another channel that cancels it.
#   coefficient of M^2_dir_I: we have 0 (no channel contributes), Drake has 3/50.
# These are INCONSISTENT. Hence the bipolar (0, 2) leading form CANNOT reproduce Drake.

# This CONFIRMS that the bipolar expansion has more terms than just the leading k_1+k_2=K.

# Let's check: is there a bipolar (2, 0) channel that COULD contribute? For direct
# path, (k_1=2, k_2=0) requires <0||C^2||0>, which is ZERO. So no, (2, 0) direct is ZERO.
# What about (k_1=0, k_2=0)? <0||C^0||0> = 1, <1||C^0||1> = sqrt(3), but triangle (0, 0, 2)
# fails. So (0, 0) excluded. Only (0, 2) for direct.

# For exchange (1, 0, 1 <-> 0, 1, 0), only (1, 1) contributes.

# So the bipolar angular reduced matrix element is unique. And the INTEGRAL it
# gives is not proportional to M^2_{dir/exch}.

# CONCLUSION: the Drake (3/50, -2/5) cannot arise from the "leading bipolar"
# expansion with just k_1+k_2=K.
#
# The ACTUAL bipolar expansion must have additional terms. Let me check what
# Brink-Satchler or Varshalovich actually say about the (k_1+k_2 > K) terms.

print("\n=== OBSERVATION ===")
print("The leading bipolar expansion k_1+k_2=K cannot reproduce Drake's formula.")
print("The bipolar expansion of Y^K/r_12^{K+1} has additional terms at k_1+k_2 > K.")
print("These additional terms contribute to the DIRECT channel (forbidden at leading)")
print("and modify the effective radial kernel.")
print()
print("The FULL Sack 1964 expansion has:")
print("  Y^K_M(hat r_12) / r_12^{K+1}")
print("     = sum_{L>=K, (L-K) even} [f_L(r_1, r_2)] * [angular term at multipole L]")
print()
print("Where [angular term] = sum_{k_1, k_2 | triangle (k_1, k_2, K), k_1+k_2=L}")
print("                     [C^{k_1}(1) x C^{k_2}(2)]^(K)")
print()
print("So the leading L=K term has k_1+k_2=K, but next-to-leading L=K+2 term has")
print("k_1+k_2=K+2 with MORE (k_1, k_2) pairs (including, e.g., (K+2, 0), (0, K+2))")
print("AND DIFFERENT radial kernels f_L.")
print()
print("For our target — <^3P_J|H_SS|^3P_J> of He (1s)(2p) — the direct ANGULAR")
print("selection forbids k_1 != 0 or 2, and the exchange forbids k_1, k_2 != 1.")
print("So only (0, 2) direct and (1, 1) exchange can contribute AT ANY L.")
print()
print("But at L=K+2=4, we need k_1+k_2=4 with (k_1, k_2) in {(0, 4), (2, 2), (4, 0)}.")
print("For direct (l_a=0, l_b=1): <0|C^2|0> = 0 but <0|C^0|0> = 1 — but (0, 4) has k_1=0,")
print("and <1|C^4|1> triangle (1, 4, 1) requires 1+1 >= 4, which fails. So (0, 4) excluded.")
print("(2, 2) has <0|C^2|0> = 0, excluded. (4, 0) has <0|C^4|0> = 0, excluded.")
print()
print("So at higher L=K+2 there are NO new angular channels for our problem.")
print()
print("Conclusion: the bipolar expansion at L=K=2 IS the full answer. The discrepancy")
print("between bipolar and Drake's M^K basis comes from the DEFINITION OF DRAKE'S M^K.")
print("Drake uses a different radial integral decomposition than r_<^K/r_>^{K+3}.")

# ============================================================
# Let me check NUMERICALLY whether some specific beta gives Drake's answer
# ============================================================

print("\n=== Solving the bipolar = Drake identification ===")
# We have:
#   A_SS_derived = scaling * (beta_02 * sp_red_dir * I_bipolar_02_dir
#                            - beta_11 * sp_red_exch * I_bipolar_11_exch)
# A_SS_target = 3/50 * M^2_dir - 2/5 * M^2_exch
#
# I_bipolar_02_dir = M^0_dir_I + M^2_dir_II
# I_bipolar_11_exch = M^1_exch_I + M^1_exch_II = M^1_exch
#
# Substituting:
# A_SS_derived = scaling * beta_02 * sp_red_dir * (M^0_dir_I + M^2_dir_II)
#              + scaling * (-beta_11) * (-sp_red_exch) * M^1_exch

# Express in the Drake basis as a vector (c_{M^0_dir_I}, c_{M^2_dir_II}, c_{M^1_exch}):
# Derived: scaling * beta_02 * sp_red_dir * 1  (on M^0_dir_I)
#          scaling * beta_02 * sp_red_dir * 1  (on M^2_dir_II)
#          scaling * (-beta_11) * (-sp_red_exch) * 1  (on M^1_exch)

# Drake formula in the same basis:
# 3/50 M^2_dir - 2/5 M^2_exch = 3/50 (M^2_dir_I + M^2_dir_II) - 2/5 (M^2_exch_I + M^2_exch_II)
# In terms of our basis (M^0_dir_I, M^2_dir_II, M^1_exch), this is:
#  3/50 on M^2_dir_I  (NOT in our basis)
#  3/50 on M^2_dir_II
#  -2/5 on M^2_exch_I  (NOT in our basis)
#  -2/5 on M^2_exch_II  (NOT in our basis)
# Drake's basis is DIFFERENT from our bipolar basis!
#
# So they can't be made equal by any choice of (beta_02, beta_11), UNLESS
# there are ADDITIONAL bipolar channels (at higher L) that contribute.

# We've shown all angular-allowed channels contribute ONLY at L=K. So the bipolar
# leading-form expansion does not match Drake.

# There must be an error in my matching of bipolar (k_1, k_2, K) to Drake pieces.

# WAIT: let me recheck. Bipolar(0, 2, K=2) direct has kernel r_1^0 r_2^2 / r_>^5.
# In Region I (r_1 < r_2, r_> = r_2): kernel = r_2^2 / r_2^5 = 1/r_2^3
# In Region II (r_1 > r_2, r_> = r_1): kernel = r_2^2 / r_1^5
#
# Drake M^K in Region I has kernel r_<^K / r_>^{K+3} = r_1^K / r_2^{K+3}:
#   K=0: 1 / r_2^3  ✓ matches Region I of bipolar(0, 2) !
#   K=2: r_1^2 / r_2^5
# Drake M^K in Region II has kernel r_<^K / r_>^{K+3} = r_2^K / r_1^{K+3}:
#   K=0: 1 / r_1^3
#   K=2: r_2^2 / r_1^5  ✓ matches Region II of bipolar(0, 2) !
#
# So Bipolar(0, 2) direct Region I = M^0_dir_I, Bipolar(0, 2) direct Region II = M^2_dir_II.
# This is what I had. So the bipolar(0, 2) kernel is "hybrid":
#   direct I = same radial kernel as M^0 direct I
#   direct II = same radial kernel as M^2 direct II

# These match piecewise, but the TOTAL is NOT M^0_dir nor M^2_dir alone.

# So our derivation gives:
#   A_SS_derived = scaling * beta_02 * sp_red_dir * (M^0_dir_I + M^2_dir_II)
#                + scaling * (-beta_11) * (-sp_red_exch) * M^1_exch
#
# But Drake gives:
#   A_SS = 3/50 M^2_dir - 2/5 M^2_exch
#        = 3/50 (M^2_dir_I + M^2_dir_II) - 2/5 (M^2_exch_I + M^2_exch_II)
#
# These have DIFFERENT CONSTITUENT INTEGRALS. They can both be CORRECT but
# phrased differently: if the total A_SS (a physical observable) is the same,
# there exists some integration-by-parts or change-of-variables that relates
# M^0_dir_I + M^2_dir_II to 3/50 * (M^2_dir_I + M^2_dir_II) or similar.
#
# ALTERNATIVELY: the derivation approach I used for A_SS is incorrect — maybe
# missing a term or using the wrong Wigner-Eckart formula.

# Actually wait — I've been using (-1)^{L+S+J} * 6j as the J-factor, but Edmonds
# 7.1.6 says the RIGHT formula for scalar tensor product is
#   <(L S) J M | [T^K(space) · U^K(spin)]_{sc} | (L'S') J' M'>
#     = delta_{L L'} delta_{S S'} delta_{J J'} delta_{M M'} (-1)^{L'+S+J}
#        * {L S J; S' L' K} / ... check the explicit formula

# Let me re-check. Messiah (2nd ed.) Eq. (C.64) or Varshalovich 13.1.4:
# For scalar tensor product of rank-K tensors T^K(1) and U^K(2) acting on electrons 1 and 2:
# <J M | [T^K(1) · U^K(2)]_0 | J M> = (-1)^{L+S+J} W(L L J S S J; K) <L||T||L><S||U||S>
# where W is the Racah W-coefficient, W = (-1)^{l+m+n+p} 6j{l m n; p q r}...
# Formally W(abcd; ef) = (-1)^{...} 6j{...}.
# The STANDARD form is (-1)^{L+S+J} 6j{L S J; S L K}.

# Wait: Edmonds 7.1.6 gives (scalar invariant):
# <(j_1 j_2) j m | T^(K)(1) . U^(K)(2) | (j_1' j_2') j m>
#   = delta_{j j'} delta_{m m'} (-1)^{j_1 + j_2' + j}
#     * 6j{j_1 j_2 j; j_2' j_1' K}
#     * <j_1 || T || j_1'> <j_2 || U || j_2'>

# For diagonal (j_1 = j_1', j_2 = j_2'): this becomes (-1)^{j_1 + j_2 + j} 6j{j_1 j_2 j; j_2 j_1 K}.
# For our (L = j_1, S = j_2, J = j): (-1)^{L+S+J} 6j{L S J; S L K}.
# That's what I used. So J-factor is correct.

# But wait — the scalar product T.U convention: T.U = sum_q (-1)^q T_q U_{-q}.
# And [T (x) U]^(0)_0 = -sum_q <K q K -q | 0 0> T_q U_{-q} / (??). Let me re-check.

# <K q K -q | 0 0> = (-1)^{K-q} / sqrt(2K+1).
# So [T (x) U]^(0)_0 = sum_q <K q K -q | 0 0> T_q U_{-q}
#                   = sum_q (-1)^{K-q} / sqrt(2K+1) T_q U_{-q}
#                   = (-1)^K / sqrt(2K+1) sum_q (-1)^(-q) T_q U_{-q}  [note q=-q trick]
#                   = (-1)^K / sqrt(2K+1) sum_q (-1)^q T_q U_{-q}     [sum_q, so (-1)^(-q) = (-1)^q by re-indexing]
#                   = (-1)^K / sqrt(2K+1) T.U
# Thus T.U = (-1)^K * sqrt(2K+1) * [T x U]^(0)_0

# For K=2: T.U = sqrt(5) * [T x U]^(0)_0. (No sign change.)

# Let's verify: are we using the SCALAR PRODUCT T.U in the Hamiltonian or [T x U]^(0)?
# H_SS = -alpha^2 [ (s1 . s2) - 3 (s1 . hat r)(s2 . hat r) ] / r^3
# The angular factor is 3 (s1 . hat r)(s2 . hat r) - (s1 . s2) = ... this is the
# rank-2 spherical tensor (s1 . s2)_{rank-2} in the coupled form.

# Standard: [s1 (x) s2]^(2)_0 = (3 s_{1z} s_{2z} - s1 . s2) / sqrt(6)  (up to convention)
# Let me re-derive:
#   [s1 (x) s2]^(K)_Q = sum_{q1 q2} <1 q1 1 q2 | K Q> s_{1, q1} s_{2, q2}
#   For K=2, Q=0: [s1 x s2]^(2)_0 = sum_{q1+q2=0} <1 q1 1 -q1 | 2 0> s_{1, q1} s_{2, -q1}
#   q1=0: <1 0 1 0 | 2 0> s_{1,0} s_{2,0} = sqrt(2/3) s_{1z} s_{2z}
#   q1=+1: <1 +1 1 -1 | 2 0> s_{1, +1} s_{2, -1} = sqrt(1/6) * s_{1, +1} s_{2, -1}
#   q1=-1: similar, sqrt(1/6) * s_{1, -1} s_{2, +1}
# This is a specific combination. The Breit-Pauli H_SS is built from this.

# H_SS = -alpha^2 * sqrt(24 pi / 5) * sum_q (-1)^q [s1 x s2]^(2)_q Y^(2)_{-q}(hat r)/r^3
# This is WRITTEN as a scalar product [A x B]^(0) ∝ sum_q (-1)^q A_q B_{-q} = A.B.
# So H_SS ∝ [s1 x s2]^(2) . Y^(2)(hat r)/r^3.

# Using T.U = (-1)^K sqrt(2K+1) [T x U]^(0):
# H_SS = -alpha^2 sqrt(24 pi / 5) * sum_q (-1)^q [s1 x s2]^(2)_q Y^(2)_{-q}
#      = -alpha^2 sqrt(24 pi / 5) * {[s1 x s2]^(2) . Y^(2)}
# For the scalar product convention, <J|T.U|J> = (-1)^{L+S+J} 6j{L S J; S L K} <L||T||L><S||U||S>.

# OK so my formula was CORRECT. And the J-factor for SS (K=2, J=1) is
#   (-1)^(1+1+1) 6j{1 1 1; 1 1 2} = (-1) * (-1/6) = 1/6
# Not -1/6 as I had. Let me recompute.

j6_val = wigner_6j(1, 1, 1, 1, 1, 2)
phase = (-1)**(1+1+1)
J_factor_corrected = phase * j6_val
print(f"\nJ-factor re-check:")
print(f"  6j(1,1,1; 1,1,2) = {j6_val}")
print(f"  (-1)^(L+S+J) = {phase}")
print(f"  J-factor = {phase * j6_val} = {float(phase * j6_val):.6f}")

# With corrected J-factor (let me use whichever sign):
scaling2 = H_prefactor_f * float(phase * j6_val) * spin_red_f
print(f"\nCorrected scaling: {scaling2:+.6e}")
bracket_target2 = A_SS_target / scaling2
print(f"Corrected bracket target: {bracket_target2:+.6e}")
print(f"Need: beta_02 * {a_coef:+.6e} + beta_11 * {b_coef:+.6e} = {bracket_target2:+.6e}")

# ============================================================
# NUMERICAL SEARCH: find (beta_02, beta_11) that reproduce Drake
# ============================================================
# We need: beta_02 * a + beta_11 * b = bracket_target2
# This is ONE equation in TWO unknowns — infinitely many solutions.
# But if we assume beta_{k_1 k_2 K} depends only on the (k_1, k_2, K) triplet (not on
# the specific matrix element), we can pick a NATURAL FORM and check.

# Expected forms from Varshalovich §5.17:
#   beta(k_1, k_2, K) = sqrt((2K+1)!/((2k_1)!(2k_2)!)) * (CG or sqrt factor)
# For (0, 2, 2):  sqrt(4!/(0! 4!)) = 1
# For (1, 1, 2):  sqrt(4!/(2! 2!)) = sqrt(6)

print("\n=== Test natural forms for (beta_02, beta_11) ===")
for beta_02, name02 in [(1.0, "1"), (1/math.sqrt(4*math.pi), "1/sqrt(4pi)"),
                         (math.sqrt(5), "sqrt(5)"), (math.sqrt(5)/math.sqrt(4*math.pi), "sqrt(5/(4pi))"),
                         (5/math.sqrt(4*math.pi), "5/sqrt(4pi)")]:
    for beta_11, name11 in [(math.sqrt(6), "sqrt(6)"), (math.sqrt(6/5), "sqrt(6/5)"),
                             (math.sqrt(6/(4*math.pi)), "sqrt(6/(4pi))"),
                             (3*math.sqrt(6), "3 sqrt(6)"),
                             (math.sqrt(30), "sqrt(30)")]:
        A_test = scaling2 * (beta_02 * a_coef + beta_11 * b_coef)
        ratio = A_test / A_SS_target
        marker = " <-- " if abs(ratio - 1) < 0.05 else ""
        if abs(ratio - 1) < 1.0:
            print(f"  beta_02={name02} ({beta_02:.4f}), beta_11={name11} ({beta_11:.4f}): ratio={ratio:+.6f}{marker}")
