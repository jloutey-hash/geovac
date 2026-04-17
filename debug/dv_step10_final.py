"""Sprint 5 DV Step 10: Compute A_SS numerically via the Racah-coupled
channel decomposition + PSLQ-identify in Drake's M^K basis.

KEY EMPIRICAL FINDINGS:
  - Bipolar angular-allowed channels for He (1s, 2p) ^3P:
      Direct:  only (k_1=0, k_2=2)   [from Step 1 channel enumeration]
      Exchange: only (k_1=1, k_2=1)
  - Piecewise decomposition (Step 7 symbolic):
      M^K_dir Region I + Region II (with direct density) has I ≠ II.
      M^K_exch Region I = Region II by r_1 <-> r_2 exchange density symmetry.

THEOREM (Drake mixing ratios from leading-order bipolar + spin-coupling):
  A_SS = -alpha^2 sqrt(6) * J_factor * spin_red *
         [spatial_red_dir * beta(0,2,2) * I_bipolar(0,2,2)_dir
          - spatial_red_exch * beta(1,1,2) * I_bipolar(1,1,2)_exch]

SPECIFIC IDENTIFICATION:

Spatial reduced ME:
  spatial_red_dir  (k_1=0, k_2=2) = -sqrt(30)/5
  spatial_red_exch (k_1=1, k_2=1) = -sqrt(5)/3

J-factor = (-1)^(L+S+J) 6j{1,1,1;1,1,2} = (-1)^3 * (-1/6) = 1/6

Spin reduced ME: <S=1||[s_1 x s_2]^(2)||S=1> = sqrt(5)/2

Bipolar prefactor (unknown — calibrated from consistency):
  beta(k_1, k_2, K) for the expansion
    C^(K)_Q(hat r_12) / r_12^{K+1} = sum_{k_1, k_2} beta(k_1, k_2, K)
       [C^{k_1}(1) (x) C^{k_2}(2)]^(K)_Q * R(k_1, k_2, K; r_1, r_2)

Leading-term (Sack-type, k_1 + k_2 = K) radial kernel:
  R(k_1, k_2, K=k_1+k_2; r_1, r_2) = r_1^{k_1} r_2^{k_2} / r_>^{2K+1}

From Sack 1964, Varshalovich §5.17:
  beta(k_1, k_2, K) = sqrt((2K+1)!/((2k_1)!(2k_2)!))
    * sqrt((2k_1+1)(2k_2+1)/(4pi(2K+1)))
    * < k_1 0 k_2 0 | K 0 >

(This is in Y^K notation; for Racah C^K we adjust by sqrt((2K+1)/4pi) factors.)

For k_1 + k_2 = K:
  beta(k_1, k_2, K) (Y^K form) = sqrt((2K)!/((2k_1)!(2k_2)!))  — the "binomial sqrt"
    * sqrt((2k_1+1)(2k_2+1)/((4pi)^2 (2K+1)))

Let me verify this at K = 0, k_1 = k_2 = 0:
  beta(0, 0, 0) = sqrt(0!/(0! 0!)) * sqrt(1*1/((4pi)^2 * 1)) = 1 * 1/(4pi)

For Coulomb (K=0, k_1=k_2=0): expansion is 1/r_12 = sum_l terms with l=0, 1, 2, ...
  The l=0 leading term should be 1/r_> = 1/r_max_of(r_1, r_2). That's not quite right.

Actually the Laplace expansion:
  1/r_12 = sum_l (r_<^l / r_>^{l+1}) P_l(cos theta_12)
  At l=0: 1/r_> (not the full expansion at K=0 "leading" form).

So the "leading term" in the bipolar expansion is l=l(leading) only for the specific
K we're expanding. For K=0, the "leading" is l_1 = l_2 = 0 with kernel 1/r_>, which
is the LEADING Laplace term. Higher l_1 + l_2 = 2, 4, ... add the l>0 contributions.

In our Breit-Pauli K=2 case, the "leading" is k_1 + k_2 = 2 with kernel r_1^{k_1}r_2^{k_2}/r_>^5.
There are ALSO higher (k_1 + k_2 = 4, 6, ...) terms. However, for He (1s, 2p) matrix elements,
angular selection rules restrict ALL channels to just (k_1, k_2) = (0, 2) direct and (1, 1)
exchange. So higher-order (k_1 + k_2 = 4+) does NOT contribute.

Thus the bipolar expansion IS effectively single-term for our problem, with kernel
r_1^{k_1} r_2^{k_2} / r_>^5.

So the expression
  A_SS = (stuff) * [spatial_red_dir * beta(0,2,2) * I_bipolar(0,2)_dir
                   - spatial_red_exch * beta(1,1,2) * I_bipolar(1,1)_exch]

should DIRECTLY equal Drake's A_SS = 3/50 M^2_dir - 2/5 M^2_exch.

Let me verify by NUMERICALLY computing both and checking.
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
sys.set_int_max_str_digits(100000)
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, pi, sympify, Symbol, nsimplify
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

from geovac.breit_integrals import breit_ss_radial, _t_kernel_region_I, _expand_product, _fraction_sqrt


def piecewise_BP_integral(n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, k):
    """Region-split Drake M^K integral. Returns (I, II) symbolic tuple."""
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
    return simplify(norm_s * I_I), simplify(norm_s * I_II)


# ============================================================
# Pre-compute all needed integrals
# ============================================================
print("Computing symbolic Drake M^K (Z=1, piecewise)...")
pieces = {}
for K_d in (0, 1, 2):
    pieces[(K_d, 'dir', 'I')], pieces[(K_d, 'dir', 'II')] = piecewise_BP_integral(
        1, 0, 2, 1, 1, 0, 2, 1, K_d
    )
    pieces[(K_d, 'exch', 'I')], pieces[(K_d, 'exch', 'II')] = piecewise_BP_integral(
        1, 0, 2, 1, 2, 1, 1, 0, K_d
    )

# Bipolar channel integrals
# Bipolar(k_1, k_2, K=2) dir = Drake M^{k_1} Region I + Drake M^{k_2} Region II
I_bipolar_02_dir = pieces[(0, 'dir', 'I')] + pieces[(2, 'dir', 'II')]
I_bipolar_02_exch = pieces[(0, 'exch', 'I')] + pieces[(2, 'exch', 'II')]
I_bipolar_11_dir = pieces[(1, 'dir', 'I')] + pieces[(1, 'dir', 'II')]
I_bipolar_11_exch = pieces[(1, 'exch', 'I')] + pieces[(1, 'exch', 'II')]
I_bipolar_20_dir = pieces[(2, 'dir', 'I')] + pieces[(0, 'dir', 'II')]
I_bipolar_20_exch = pieces[(2, 'exch', 'I')] + pieces[(0, 'exch', 'II')]

print(f"\nBipolar (0, 2, K=2) direct:  {float(I_bipolar_02_dir):+.6e}")
print(f"Bipolar (0, 2, K=2) exchange: {float(I_bipolar_02_exch):+.6e}")
print(f"Bipolar (1, 1, K=2) direct:   {float(I_bipolar_11_dir):+.6e}")
print(f"Bipolar (1, 1, K=2) exchange: {float(I_bipolar_11_exch):+.6e}")

# ============================================================
# Assemble A_SS using the bipolar formula with unknown beta
# ============================================================
beta_02 = Symbol('beta_02')  # beta(k_1=0, k_2=2, K=2)
beta_11 = Symbol('beta_11')  # beta(k_1=1, k_2=1, K=2)

# Spatial reduced ME (already computed):
sp_red_dir = -sqrt(30) / 5   # (k_1=0, k_2=2) direct
sp_red_exch = -sqrt(5) / 3   # (k_1=1, k_2=1) exchange

# Spin reduced ME:
spin_red = sqrt(5) / 2

# J-factor (at J=1, L=S=1, K=2):
J_factor = Rational(-1, 6)  # (-1)^(L+S+J) * 6j = -1/6

# H_SS prefactor -alpha^2 sqrt(6):
H_prefactor = -sqrt(6)  # (alpha^2 pulled out)

# Apply the Wigner-Eckart scalar formula:
# <^3P_1, 1 | H_SS | ^3P_1, 1> = H_prefactor * J_factor * spin_red
#                                 * [sp_red_dir * beta_02 * I_bipolar_02_dir
#                                    - sp_red_exch * beta_11 * I_bipolar_11_exch]

A_SS_derived = H_prefactor * J_factor * spin_red * (
    sp_red_dir * beta_02 * I_bipolar_02_dir
    - sp_red_exch * beta_11 * I_bipolar_11_exch
)
A_SS_derived = simplify(A_SS_derived)
# (skip symbolic print — contains huge integers in log args)
print("\nA_SS derived (with beta_02, beta_11 symbolic): [skipped symbolic print; extremely long]")

# Numerical value of bracket [...] with unit beta:
unit_bracket_val = (
    float(sp_red_dir) * 1 * float(I_bipolar_02_dir)
    - float(sp_red_exch) * 1 * float(I_bipolar_11_exch)
)
print(f"Unit bracket (beta = 1): {unit_bracket_val:+.6e}")

# Target A_SS (Drake):
M2_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1)
M2_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1)
A_SS_target = Rational(3, 50) * M2_dir - Rational(2, 5) * M2_exch
print(f"\nA_SS target (Drake, Z=1): {float(A_SS_target):+.6e}")

# ============================================================
# Factor out the common prefactor from A_SS_derived
# ============================================================
# A_SS_derived = H_prefactor * J_factor * spin_red * (a * beta_02 + b * beta_11)
#              = scaling * (a * beta_02 + b * beta_11)
# where scaling = H_prefactor * J_factor * spin_red

scaling = H_prefactor * J_factor * spin_red
a_coef = simplify(sp_red_dir * I_bipolar_02_dir)
b_coef = simplify(-sp_red_exch * I_bipolar_11_exch)

print(f"\nA_SS = scaling * (beta_02 * a_coef + beta_11 * b_coef)")
print(f"  scaling = {float(scaling):.6e}  (symbolic: -sqrt(6)*(-1/6)*sqrt(5)/2 = sqrt(30)/12)")
print(f"  a_coef = sp_red_dir * I_bipolar_02_dir = {float(a_coef):+.6e}")
print(f"  b_coef = -sp_red_exch * I_bipolar_11_exch = {float(b_coef):+.6e}")

# Target: A_SS_target should equal scaling * (beta_02 * a_coef + beta_11 * b_coef)
# i.e. beta_02 * a_coef + beta_11 * b_coef = A_SS_target / scaling
bracket_target = simplify(A_SS_target / scaling)
print(f"\nbracket target = A_SS_target / scaling")
print(f"              = {float(bracket_target):+.6e}")
print(f"Solve: beta_02 * {float(a_coef):+.6e} + beta_11 * {float(b_coef):+.6e} = {float(bracket_target):+.6e}")

# ============================================================
# Sanity: what would (beta_02, beta_11) be if we used the EXPECTED Steinborn-Filter-type
# formula beta(k_1, k_2, K) = sqrt((2K+1)!/((2k_1)!(2k_2)!))... ?
# ============================================================
# For (k_1=0, k_2=2, K=2): (2K)! = 4! = 24, (2k_1)! = 1, (2k_2)! = 24
#   ratio = 24/24 = 1, sqrt = 1. So beta_02 (core form) = 1.
# For (k_1=1, k_2=1, K=2): (2k_1)! = 2, (2k_2)! = 2, ratio = 24/4 = 6, sqrt(6).
#   So beta_11 (core form) = sqrt(6).

# Then sol: with these beta values, what A_SS do we get?
beta_02_test = 1  # Guess: sqrt((2K+1)!/((2k_1)!(2k_2)!)) at (0,2,2) = 1
beta_11_test = sqrt(6)  # At (1,1,2) = sqrt(6)

A_SS_test = simplify(scaling * (a_coef * beta_02_test + b_coef * beta_11_test))
ratio_f = float(A_SS_test) / float(A_SS_target)
print(f"\nWith beta_02=1, beta_11=sqrt(6) (Steinborn-Filter form):")
print(f"  A_SS_test = {float(A_SS_test):+.6e}")
print(f"  A_SS_target = {float(A_SS_target):+.6e}")
print(f"  Ratio (test / target) = {ratio_f:.6f}")

# Try other natural forms
print("\nSearching for (beta_02, beta_11) that gives A_SS_test / A_SS_target = 1:")
for beta_02_try in [Rational(1, 2), Rational(1), sqrt(2), Rational(1, 2), Rational(2)]:
    for beta_11_try in [Rational(1), sqrt(6), sqrt(Rational(6, 5)), Rational(3), sqrt(Rational(3, 2)), sqrt(Rational(5, 2))]:
        A_test = simplify(scaling * (a_coef * beta_02_try + b_coef * beta_11_try))
        r_f = float(A_test) / float(A_SS_target)
        if abs(r_f - 1.0) < 0.05:
            print(f"  beta_02={beta_02_try}, beta_11={beta_11_try}: ratio = {r_f:.6f}")
