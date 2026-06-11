"""Sprint 5 DV Step 7: Correct piecewise decomposition of Drake M^K.

Previous Step 6 had a bug: the labels for electron 1 / 2 densities
were mis-passed to `production_M_K_piecewise`. Now fixed: for DIRECT
channel, electron 1 has 1s (both left/right), electron 2 has 2p
(both). For EXCHANGE, electron 1 has 1s left / 2p right, electron 2
has 2p left / 1s right.

Drake M^K radial integral (direct channel):
  M^K_dir = ∫∫ R_{1s}(r_1)^2 R_{2p}(r_2)^2 (r_<^K / r_>^{K+3}) r_1^2 r_2^2 dr_1 dr_2

The function `breit_ss_radial(n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, K, Z)`
computes:
  R^K_{BP}(n_a l_a, n_c l_c; n_b l_b, n_d l_d)
  = ∫∫ [R_{a}(r_1) R_{c}(r_1)] * [R_{b}(r_2) R_{d}(r_2)] * (r_<^K / r_>^{K+3}) * r_1^2 r_2^2

where (a, c) are on electron 1 and (b, d) are on electron 2.

  DIRECT: (a, c) = (1s, 1s), (b, d) = (2p, 2p).
    breit_ss_radial(1, 0, 1, 0, 2, 1, 2, 1, K)  - Wait, check signature!
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
from sympy import Rational, sqrt, Integer, simplify, symbols, expand, collect, nsimplify

from geovac.breit_integrals import (
    breit_ss_radial,
    _t_kernel_region_I,
    _expand_product,
    _fraction_sqrt,
    _simplify_log_form,
)


def piecewise_BP_integral(
    n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, k
):
    """Region-split Breit-Pauli Slater integral.

    Returns (R^k_I, R^k_II) where the full Drake integral is R^k_I + R^k_II.

    R^k_{BP}(a,c; b,d) = integral rho_{a,c}(r_1) rho_{b,d}(r_2) r_<^k / r_>^(k+3) r_1^2 r_2^2 dr_1 dr_2
    rho_{a,c}(r_1) = R_{a,l_a}(r_1) R_{c,l_c}(r_1) r_1^2  (with r^2 included per our convention)

    Here following the production call signature:
       breit_ss_radial(n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, k, Z)
       n_a, l_a, n_c, l_c on electron 1
       n_b, l_b, n_d, l_d on electron 2
    """
    # Electron 1 density product: (n_a l_a) (n_c l_c)
    terms_e1, nsq_e1 = _expand_product(n_a, l_a, n_c, l_c)
    # Electron 2 density product: (n_b l_b) (n_d l_d)
    terms_e2, nsq_e2 = _expand_product(n_b, l_b, n_d, l_d)
    norm = _fraction_sqrt(nsq_e1 * nsq_e2)
    norm_s = Rational(norm.numerator, norm.denominator)

    I_I = Rational(0)
    I_II = Rational(0)
    for c13, p13, a13 in terms_e1:
        for c24, p24, a24 in terms_e2:
            if c13 == 0 or c24 == 0:
                continue
            # Region I (r_1 < r_2): integrand r_1^{p13} r_2^{p24} e^{-a13 r_1 - a24 r_2} * r_<^k / r_>^{k+3}
            # with r_< = r_1, r_> = r_2 ==> kernel = r_1^k / r_2^{k+3}
            I_I = I_I + c13 * c24 * _t_kernel_region_I(p13, p24, a13, a24, "breit", k)
            # Region II (r_1 > r_2): r_< = r_2, r_> = r_1 ==> kernel = r_2^k / r_1^{k+3}
            # Swap electron labels in the integrand: integrate first r_2 from 0 to r_1, then r_1 from 0 to oo
            # Let u = r_2, v = r_1. Integrand u^{p24} v^{p13} e^{-a13 v - a24 u} * u^k / v^{k+3}
            # with u < v. This is _t_kernel_region_I with (a=p24, b=p13, alpha=a24, beta=a13, k).
            I_II = I_II + c13 * c24 * _t_kernel_region_I(p24, p13, a24, a13, "breit", k)

    return simplify(norm_s * I_I), simplify(norm_s * I_II)


print("=" * 76)
print("Sprint 5 DV Step 7: Corrected piecewise Drake M^K")
print("=" * 76)

# ============================================================
# Compute Drake M^K pieces (K = 0, 1, 2) CORRECTLY
# ============================================================
pieces = {}
for K_d in (0, 1, 2):
    # DIRECT: (a, c) = (1s, 1s), (b, d) = (2p, 2p)
    I_d, II_d = piecewise_BP_integral(1, 0, 2, 1, 1, 0, 2, 1, K_d)
    # EXCHANGE: (a, c) = (1s, 2p), (b, d) = (2p, 1s)
    I_e, II_e = piecewise_BP_integral(1, 0, 2, 1, 2, 1, 1, 0, K_d)
    pieces[(K_d, 'dir', 'I')] = _simplify_log_form(I_d)
    pieces[(K_d, 'dir', 'II')] = _simplify_log_form(II_d)
    pieces[(K_d, 'exch', 'I')] = _simplify_log_form(I_e)
    pieces[(K_d, 'exch', 'II')] = _simplify_log_form(II_e)

    total_d = simplify(I_d + II_d)
    total_e = simplify(I_e + II_e)
    prod_d = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, K_d, Z=1)
    prod_e = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, K_d, Z=1)

    ok_d = "MATCH" if abs(float(simplify(total_d - prod_d))) < 1e-12 else "DIFFER"
    ok_e = "MATCH" if abs(float(simplify(total_e - prod_e))) < 1e-12 else "DIFFER"
    print(f"\n  K = {K_d}")
    print(f"    M^{K_d}_dir: I = {float(I_d):+.6e}, II = {float(II_d):+.6e}, total = {float(total_d):+.6e} ({ok_d})")
    print(f"    M^{K_d}_exch: I = {float(I_e):+.6e}, II = {float(II_e):+.6e}, total = {float(total_e):+.6e} ({ok_e})")
    print(f"    Symmetry: Region I dir = Region II dir? {simplify(I_d - II_d) == 0}")
    print(f"    Symmetry: Region I exch = Region II exch? {simplify(I_e - II_e) == 0}")

# ============================================================
# Symbolic pieces for M^K_dir (Region I, II separately)
# ============================================================
print("\n\n=== M^K_dir (direct) symbolic pieces ===")
for K_d in (0, 1, 2):
    I_d = pieces[(K_d, 'dir', 'I')]
    II_d = pieces[(K_d, 'dir', 'II')]
    print(f"\n  M^{K_d}_dir Region I:  {I_d}")
    print(f"  M^{K_d}_dir Region II: {II_d}")

# ============================================================
# Bipolar channel integrals expressed in the piecewise Drake basis
# ============================================================
#
# For bipolar (k_1, k_2, K=2) with k_1 + k_2 = K, kernel = r_1^{k_1} r_2^{k_2} / r_>^5.
#   Region I (r_1 < r_2): r_> = r_2, kernel = r_1^{k_1} / r_2^{5 - k_2}
#     = r_1^{k_1} / r_2^{k_1 + 3}  (using k_1 + k_2 = 2 ⟹ 5 - k_2 = 5 - (2 - k_1) = 3 + k_1)
#     = Drake M^{K_d = k_1} Region I
#   Region II (r_1 > r_2): r_> = r_1, kernel = r_2^{k_2} / r_1^{5 - k_1}
#     = r_2^{k_2} / r_1^{k_2 + 3}
#     = Drake M^{K_d = k_2} Region II

print("\n\n=== Bipolar channels in Drake basis ===")
# Note: we only have the 3 channels (k_1, k_2) ∈ {(0, 2), (1, 1), (2, 0)} for K=2.
bipolar_decomp = {}
for k1, k2 in [(0, 2), (1, 1), (2, 0)]:
    # Bipolar(k1, k2) dir = Drake M^{k1}_dir_I + Drake M^{k2}_dir_II
    dir_val = pieces[(k1, 'dir', 'I')] + pieces[(k2, 'dir', 'II')]
    exch_val = pieces[(k1, 'exch', 'I')] + pieces[(k2, 'exch', 'II')]
    bipolar_decomp[(k1, k2, 'dir')] = simplify(dir_val)
    bipolar_decomp[(k1, k2, 'exch')] = simplify(exch_val)
    print(f"\n  Bipolar(k_1={k1}, k_2={k2}) dir  = M^{k1}_dir_I + M^{k2}_dir_II = {float(dir_val):+.6e}")
    print(f"  Bipolar(k_1={k1}, k_2={k2}) exch = M^{k1}_exch_I + M^{k2}_exch_II = {float(exch_val):+.6e}")

# ============================================================
# Now: the Drake total M^K = Region I + Region II
# We can EXTRACT Region I and Region II:
#   M^K_Region_I  = pieces[(K, 'dir', 'I')]
#   M^K_Region_II = pieces[(K, 'dir', 'II')]
# with M^K = Region I + Region II.
# ============================================================

# ============================================================
# STEP 8: Evaluate the full matrix element <^3P_1 | H_SS | ^3P_1>
#         with the bipolar prefactor correctly included.
# ============================================================
print("\n\n=== Step 8: Full matrix element assembly ===")

# We ASSUME the bipolar expansion is (k_1 + k_2 = K = 2 leading form):
#   C^(2)_Q(hat r_12) / r_12^3 = sum_{k1 + k2 = 2} alpha(k_1, k_2, 2)
#                                  [C^{k_1}(1) x C^{k_2}(2)]^{(2)}_Q
#                                  * r_1^{k_1} r_2^{k_2} / r_>^5
# with alpha to be determined.

from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

def red_C(l, lp, k):
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * wigner_3j(l, k, lp, 0, 0, 0)


# Spin reduced ME:
S = 1
k_spin = 2
n9_spin = wigner_9j(Rational(1, 2), Rational(1, 2), S, 1, 1, k_spin, Rational(1, 2), Rational(1, 2), S)
red_spin = simplify(sqrt(Integer((2 * S + 1) * (2 * S + 1) * (2 * k_spin + 1))) * n9_spin * Rational(3, 2))
print(f"<S=1||[s_1 x s_2]^(2)||S=1> = {red_spin}  (expected sqrt(5)/2)")

# J-factor for J=1, L=S=1, K=2:
L_val, S_val, J_val = 1, 1, 1
K_op = 2
phase = (-1) ** (L_val + S_val + J_val)
j6 = wigner_6j(L_val, S_val, J_val, S_val, L_val, K_op)
# Include a (-1)^K factor from T.U vs [TxU]^(0) convention? Let me include.
# Edmonds 7.1.6 (scalar):
# <J M | [T^K(1) · U^K(2)]_{scalar} | J M> = (-1)^{L+S+J} 6j{LSJ;SLK} <L||T||L><S||U||S>
# where the scalar product T.U convention is T.U = sum_q (-1)^q T^K_q U^K_{-q}.
J_factor = phase * j6
J_factor = simplify(J_factor)
print(f"(-1)^(L+S+J) * 6j{{LSJ;SLK}} = {J_factor}  (should be -1/6 for K=2, J=1)")

# Spatial reduced ME per channel (direct path, (k_1, k_2) = (0, 2)):
k1, k2 = 0, 2
r1_factor = red_C(0, 0, k1)
r2_factor = red_C(1, 1, k2)
n9_dir = wigner_9j(0, 1, 1, 0, 1, 1, k1, k2, K_op)
# Brink-Satchler 4.6.8:
# <(l_a l_b) L || [A^{k_1}(1) x B^{k_2}(2)]^K || (l_c l_d) L'>
#     = sqrt((2L+1)(2L'+1)(2K+1)) 9j{l_a l_b L; l_c l_d L'; k_1 k_2 K}
#         <l_a || A || l_c> <l_b || B || l_d>
norm_L = sqrt(Integer((2 * 1 + 1) * (2 * 1 + 1) * (2 * K_op + 1)))
spatial_red_dir_raw = norm_L * n9_dir * r1_factor * r2_factor
spatial_red_dir_raw = simplify(spatial_red_dir_raw)
print(f"\nDirect spatial reduced ME (k_1=0, k_2=2): {spatial_red_dir_raw}")

# Spatial reduced ME (exchange path, (k_1, k_2) = (1, 1)):
k1, k2 = 1, 1
r1_factor = red_C(0, 1, k1)
r2_factor = red_C(1, 0, k2)
n9_exch = wigner_9j(0, 1, 1, 1, 0, 1, k1, k2, K_op)
spatial_red_exch_raw = norm_L * n9_exch * r1_factor * r2_factor
spatial_red_exch_raw = simplify(spatial_red_exch_raw)
print(f"Exchange spatial reduced ME (k_1=1, k_2=1): {spatial_red_exch_raw}")

# ============================================================
# Now: the bipolar prefactor for each (k_1, k_2) channel.
# ============================================================
#
# We need alpha(k_1, k_2, K) in:
#   C^(K)_Q(hat r_12) / r_12^{K+1} = sum_{k1+k2=K} alpha(k_1, k_2, K)
#       * [C^{k_1}(1) x C^{k_2}(2)]^(K)_Q
#       * r_1^{k_1} r_2^{k_2} / r_>^{2K+1}
#
# The Steinborn-Filter formula (adapted to Racah C^K):
#   alpha(k_1, k_2, K) = (-1)^{K} * sqrt((2K+1)! / ((2k_1+1)! (2k_2+1)!))
#     * sqrt((2k_1+1)(2k_2+1) / (4pi)^? (2K+1)^?)
#
# I'll DERIVE alpha by requiring a known identity.
# Test: at k_1 = k_2 = 0 (trivially), alpha(0, 0, 0) * (r_<^0/r_>^1) must match
# the leading Laplace term 1/r_12 → (r_<^0/r_>^1) for l=0. But:
#   1/r_12 = sum_l (r_<^l / r_>^{l+1}) (-1)^l sqrt(2l+1) [C^l(1) x C^l(2)]^(0)_0
# At l=0: 1/r_> * 1 * 1 * [C^0 x C^0]^(0)_0 * sqrt(1)
# ≈ The r_<^0/r_>^1 term coefficient in 1/r_12 expansion is sqrt(2l+1)·(-1)^l at l=0.
# For K = 0 our bipolar formula restricts k_1 + k_2 = 0, i.e. k_1 = k_2 = 0.
# Here 1/r_12 = C^0(hat r_12)/r_12 * sqrt(4 pi). Hmm, but C^0 = 1 identically, and
# Y^0 = 1/sqrt(4pi). So C^0(hat r_12) / r_12^1 = 1/r_12.
# Our bipolar identity at K=0 says:
#   1/r_12 = alpha(0, 0, 0) * [C^0(1) x C^0(2)]^(0)_0 * r_1^0 r_2^0 / r_>^1
#          = alpha(0, 0, 0) * 1 * 1/r_>
# But 1/r_12 != 1/r_>; they differ by multipole expansion terms.
# So the "k_1 + k_2 = K leading term" approximation is MISSING the l > 0 terms!
# The bipolar expansion is NOT just one term at k_1 + k_2 = K.

# The CORRECT bipolar expansion has an infinite series in l (with l_1 + l_2 ≥ K,
# triangle, parity). For the He (1s, 2p) case, angular selection restricts l_1 and l_2
# to VERY specific values:
#   Direct:  l_1 = 0 (<0|C^{l_1}|0>), l_2 even (triangle with K) ⇒ l_1 = 0, l_2 = 0 or 2
#   Exchange: l_1 = 1, l_2 = 1

# For direct channel (K=2 direct):
#   (l_1, l_2) = (0, 2) — triangle with K=2 OK
#   (l_1, l_2) = (0, 0) — triangle |0-0| ≤ K=2 ≤ 0+0=0? NO, K=2 > 0. Excluded.

# So only (0, 2) contributes to direct even in the full bipolar. Good.

# For exchange (K=2 exchange): only (1, 1).

# So THE LEADING-TERM bipolar IS exact for this problem. But the bipolar
# expansion has multiple ranks of "leading" for different K — for K=0 there are
# (0,0), (1,1), (2,2), ... for K=2 there are (0,2), (1,1), (2,0), (1,3), (3,1),
# (2,2) no (triangle allows |K-k_1| ≤ k_2 ≤ K+k_1 so (2,2) triangle OK), (4,2), (2,4), ...
# but he (1s, 2p) restricts l_i to at most l_a + l_b + something. Specifically
# <l=0||C^{k_1}||l'=0> is zero unless k_1 = 0, and <l=1||C^{k_2}||l'=1> is zero
# unless k_2 is even ≤ 2.

# So for direct channel only (k_1, k_2) = (0, 2) contributes. And we need the
# FULL bipolar prefactor for (0, 2, K=2), not just the "leading" one.

# The CORRECT formula is more nuanced. Let me derive it via a cleaner route.

# ============================================================
# CLEANER DERIVATION: direct matrix element in terms of Slater integrals
# using Brink-Satchler D.3 (direct product expansion)
# ============================================================
#
# The angular matrix element
# <l_a m_a; l_b m_b | C^(K)_Q(hat r_12) / r_12^{K+1} | l_c m_c; l_d m_d>
# can be evaluated directly as:
#
# = sum_{k_1, k_2, q_1, q_2} alpha(k_1, k_2, K)
#     * <k_1 q_1; k_2 q_2 | K Q>
#     * <l_a m_a | C^{k_1}_{q_1} | l_c m_c>
#     * <l_b m_b | C^{k_2}_{q_2} | l_d m_d>
#     * G^K_{k_1 k_2}(r_1, r_2)   [radial part]
#
# where G is the bipolar radial kernel.
#
# The radial Slater integral of G^K for the (1s)(2p) He problem reduces to
# a specific combination of Drake's M^K — exactly via the piecewise decomposition.

# ============================================================
# SIMPLIFICATION: the ONLY (k_1, k_2, K) triplets that contribute to the HeI
# (1s)(2p) direct or exchange matrix elements are (0, 2, 2) and (1, 1, 2) resp.
#
# So if we can just NUMERICALLY COMPUTE these matrix elements from scratch
# (via scipy dblquad on the angular integral) and compare to Drake's formula,
# we can identify the mixing ratios.
#
# Let me do that in the next step.

print("\n[Step 8 continues in dv_step8_angular_numerical.py]")
