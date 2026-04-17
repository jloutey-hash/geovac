"""Sprint 5 DV Step 5: Analytic derivation of Drake mixing ratios via
bipolar channel decomposition.

KEY ALGEBRAIC IDENTITY (derived in Step 4):

For the bipolar kernel r_1^{k_1} r_2^{k_2} / r_>^{2K+1} at K=2 (so denom = r_>^5):
  Bipolar (k_1 = 0, k_2 = 2) Region I: r_1^0 r_2^{-3} = 1/r_2^3 = Drake M^0 Region I
  Bipolar (0, 2) Region II: r_1^{-5} r_2^2 = Drake M^2 Region II
  Bipolar (1, 1) Region I: r_1 r_2^{-4} = Drake M^1 Region I
  Bipolar (1, 1) Region II: r_1^{-4} r_2 = Drake M^1 Region II
  Bipolar (2, 0) Region I: r_1^2 r_2^{-5} = Drake M^2 Region I
  Bipolar (2, 0) Region II: r_1^{-3} = Drake M^0 Region II

So each bipolar channel is a HYBRID of two Drake channels (one each in
regions I and II). For the EXCHANGE density (which is r_1 <-> r_2
symmetric), Region I = Region II, so:
  Bipolar(0,2) exch = M^0_exch/2 + M^2_exch/2
  Bipolar(1,1) exch = M^1_exch
  Bipolar(2,0) exch = M^2_exch/2 + M^0_exch/2

For the DIRECT density (asymmetric: rho_1s(r_1) * rho_2p(r_2)), Region
I and II contribute different fractions. Let's compute them.

THE BIPOLAR EXPANSION
=====================

Sack 1964 / Brink-Satchler App 5: the bipolar expansion of Y^K_Q(r_hat_12) / r_12^{K+1}
has a specific form. In the "leading" version (k_1 + k_2 = K):

  Y^K_Q(hat r_12) / r_12^{K+1}
    = sum_{k_1=0}^{K} beta(k_1, k_2 = K-k_1, K)
      * [Y^{k_1}(hat r_1) (x) Y^{k_2}(hat r_2)]^K_Q
      * R(k_1, k_2 = K-k_1, K; r_1, r_2)

with R = r_1^{k_1} r_2^{k_2} / r_>^{2K+1} (!)
and beta(k_1, k_2, K) a specific rational coefficient.

The BETA coefficient comes from the Gegenbauer expansion:

  beta(k_1, k_2, K) = ((2K+1)!/(2k_1+1)! (2k_2+1)!)^{1/2}   [Brink-Satchler 4.6.13]

or equivalently

  beta(k_1, k_2, K) = 1/(4pi)^{1/2} * ((2K+1)!/(2k_1)! (2k_2)!)^{1/2}
                        * (1/(2K+1)^{1/2})   [Varshalovich 5.17]

There are MULTIPLE conventions floating around. Let me DERIVE the correct
one by matching to the Coulomb expansion K=0:

For K=0, the standard identity is:
  1/r_12 = sum_l (r_<^l / r_>^{l+1}) * (4pi/(2l+1)) sum_m Y^l_m(1)* Y^l_m(2)
         = sum_l (r_<^l / r_>^{l+1}) * (4pi/(2l+1)) sqrt((2l+1)/(4pi))^2
            * C^l(1) . C^l(2)   [in Racah norm]
         = sum_l (r_<^l / r_>^{l+1}) * (2l+1)/(4pi)  ??

Hmm, let me do it from scratch with Condon-Shortley phase and Racah's C^k tensor.

C^k_q(hat r) = sqrt(4pi/(2k+1)) Y^k_q(hat r)
Y^k_q(hat r) = sqrt((2k+1)/(4pi)) C^k_q(hat r)

Laplace expansion:
  1/r_12 = sum_l (r_<^l / r_>^{l+1}) P_l(cos theta_12)
         = sum_l (r_<^l / r_>^{l+1}) C^l(1) . C^l(2)

where C^l(1) . C^l(2) = sum_m (-1)^m C^l_m(1) C^l_{-m}(2) is the scalar
product of rank-l tensors.

In terms of [C^l(1) (x) C^l(2)]^(0)_0 = -1/sqrt(2l+1) * C^l(1) . C^l(2),
we have
  1/r_12 = -sum_l (r_<^l / r_>^{l+1}) sqrt(2l+1) [C^l(1) (x) C^l(2)]^(0)_0

Note 1/r_12 = 1/|r_12| is a scalar (rank 0). This is the Coulomb
expansion in bipolar form for Y^0(hat r_12)/r_12^1.

Now: Y^K_Q(hat r_12) / r_12^{K+1} is a rank-K spherical tensor in hat r_12.
It DOES NOT have a simple bipolar expansion in one-variable Y^K on each
electron, because Y^K(hat r_12) is not a product of functions of hat r_1
and hat r_2.

The correct identity (Wen 1973, J. Math. Phys. 14, 1320; Steinborn-Filter 1975):

  f(r_12) Y^K_M(hat r_12) = sum_{k_1, k_2} F_{k_1, k_2}^{f, K}(r_1, r_2)
                              * [Y^{k_1}(hat r_1) (x) Y^{k_2}(hat r_2)]^(K)_M

with triangle (k_1, k_2, K) and parity (k_1 + k_2 + K) even.

FOR OUR SPECIFIC CASE f(r_12) = 1/r_12^{K+1}, K=2 (breit-pauli SS):

The coefficients F_{k_1,k_2}^{1/r^{K+1}, K} are derived in Wen 1973 Eq. (27)
or equivalently in Sack 1964. For K=2 and 1/r_12^3:

  F_{0, 2} = sqrt(5/(4pi)) * (3/pi)^{1/2} * (r_<^2 / r_>^5)??

(I don't have the reference at hand; let me derive it via an alternative route.)

RECURSION APPROACH
==================

We know:
  1/r_12 = sum_l (r_<^l / r_>^{l+1}) C^l(1) . C^l(2)   [Laplace expansion]

Apply ∇_{r_12} to both sides:
  ∇(1/r_12) = -r_12 / r_12^3

Squared norm: (∇(1/r_12))^2 = r_12^2/r_12^6 = 1/r_12^4.

Not directly useful.

Alternative: differentiate
  1/r_12^n = -1/(n-2) ∇^2(1/r_12^{n-2})   [for n != 2]

in bipolar form. But applying Laplacians is complex.

PRAGMATIC APPROACH: directly numerically compute each bipolar channel
integral and express as linear combination of M^K. We did this in
Step 2 / Step 3. The key finding from Step 2 was:
  Bipolar(0, 2, K=2) direct = 2.591 * M^2_dir
  Bipolar(1, 1, K=2) exch   = 1.320 * M^2_exch

Let's carefully compute these ratios symbolically by splitting region I
and region II integrals of the hydrogenic densities directly.
"""
from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, symbols, exp, integrate, oo, log, Matrix

from geovac.breit_integrals import breit_ss_radial

r1, r2 = symbols('r_1 r_2', positive=True)

# Hydrogenic densities (with r^2 volume element included)
# At Z=1:
#   R_1s(r) = 2 e^{-r}, so (R_1s)^2 r^2 = 4 r^2 e^{-2r}
#   R_2p(r) = (1/(2 sqrt 6)) r e^{-r/2}, so (R_2p)^2 r^2 = (1/24) r^4 e^{-r}
#   R_1s R_2p r^2 = (1/sqrt 6) r^3 e^{-3r/2}
d_1s_sq = 4 * r1**2 * exp(-2 * r1)    # dens for electron 1 in "direct" channel
d_2p_sq = Rational(1, 24) * r2**4 * exp(-r2)  # dens for electron 2 in direct
d_1s2p_r1 = r1**3 * exp(-3 * r1 / 2) / sqrt(6)  # dens for electron 1 in exchange
d_2p1s_r2 = r2**3 * exp(-3 * r2 / 2) / sqrt(6)  # dens for electron 2 in exchange


def int_piecewise(dens_r1_expr, dens_r2_expr, p1, p2, Q):
    """∫∫ dens_r1_expr(r_1) dens_r2_expr(r_2) * r_1^p1 r_2^p2 / r_>^Q
    split into Region I (r_1 < r_2) and Region II (r_1 > r_2).

    Returns (I_region_I, I_region_II) each symbolically.
    """
    # Region I: r_1 < r_2, so r_> = r_2
    kernel_I = r1**p1 * r2**p2 / r2**Q
    inner = integrate(dens_r1_expr * dens_r2_expr * kernel_I, (r1, 0, r2))
    I_I = integrate(inner, (r2, 0, oo))
    # Region II: r_1 > r_2, so r_> = r_1
    kernel_II = r1**p1 * r2**p2 / r1**Q
    inner = integrate(dens_r1_expr * dens_r2_expr * kernel_II, (r1, r2, oo))
    I_II = integrate(inner, (r2, 0, oo))
    return simplify(I_I), simplify(I_II)


# ============================================================
# Compute each Drake M^K piecewise
# ============================================================
print("Computing Drake M^K piecewise integrals...")
pieces = {}
for K_d in (0, 1, 2):
    dir_I, dir_II = int_piecewise(d_1s_sq, d_2p_sq, K_d, 0, K_d + 3)
    dir_I_II, dir_II_II = int_piecewise(d_1s_sq, d_2p_sq, 0, K_d, K_d + 3)
    # Wait -- for M^K we need r_<^K in each region. Region I: r_< = r_1, so r_1^K.
    # Region II: r_< = r_2, so r_2^K. So for direct, we need two different integrands.
    # Let's just re-do properly.

    # DIRECT Region I: kernel = r_1^K / r_2^{K+3}
    dir_I = integrate(integrate(d_1s_sq * d_2p_sq * r1**K_d / r2**(K_d + 3), (r1, 0, r2)), (r2, 0, oo))
    # DIRECT Region II: kernel = r_2^K / r_1^{K+3}
    dir_II = integrate(integrate(d_1s_sq * d_2p_sq * r2**K_d / r1**(K_d + 3), (r1, r2, oo)), (r2, 0, oo))
    # EXCHANGE Region I: same form
    exch_I = integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r1**K_d / r2**(K_d + 3), (r1, 0, r2)), (r2, 0, oo))
    exch_II = integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r2**K_d / r1**(K_d + 3), (r1, r2, oo)), (r2, 0, oo))
    pieces[(K_d, 'dir', 'I')] = simplify(dir_I)
    pieces[(K_d, 'dir', 'II')] = simplify(dir_II)
    pieces[(K_d, 'exch', 'I')] = simplify(exch_I)
    pieces[(K_d, 'exch', 'II')] = simplify(exch_II)
    total_d = simplify(dir_I + dir_II)
    total_e = simplify(exch_I + exch_II)
    prod_d = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, K_d, Z=1)
    prod_e = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, K_d, Z=1)
    print(f"  M^{K_d}_dir: I = {float(dir_I):+.6e}, II = {float(dir_II):+.6e}, "
          f"total = {float(total_d):+.6e} (production: {float(prod_d):+.6e})")
    print(f"  M^{K_d}_exch: I = {float(exch_I):+.6e}, II = {float(exch_II):+.6e}, "
          f"total = {float(total_e):+.6e} (production: {float(prod_e):+.6e})")
    if K_d == 2:
        print("  M^2_dir Region I symbolic:  ", pieces[(2, 'dir', 'I')])
        print("  M^2_dir Region II symbolic: ", pieces[(2, 'dir', 'II')])

# ============================================================
# Verify the bipolar channel identities
# ============================================================
print("\n\n=== Bipolar = piecewise Drake identities ===")

# Bipolar (k_1, k_2, K=2): kernel r_1^{k_1} r_2^{k_2} / r_>^5

# Bipolar (0, 2) direct: kernel 1/r_2^3 in region I (= M^0 region I) + r_2^2/r_1^5 in region II (= M^2 region II)
# Predicted:
pred_02_dir = pieces[(0, 'dir', 'I')] + pieces[(2, 'dir', 'II')]
# Let's compute the actual bipolar integral
actual_02_dir = integrate(integrate(d_1s_sq * d_2p_sq * r1**0 * r2**2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                integrate(integrate(d_1s_sq * d_2p_sq * r1**0 * r2**2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
actual_02_dir = simplify(actual_02_dir)
pred_02_dir = simplify(pred_02_dir)
print(f"\nBipolar(0,2)_dir:")
print(f"  Predicted (M^0_dir_I + M^2_dir_II): {pred_02_dir} = {float(pred_02_dir):+.6e}")
print(f"  Actual                             : {actual_02_dir} = {float(actual_02_dir):+.6e}")
print(f"  Difference                         : {simplify(pred_02_dir - actual_02_dir)}")

pred_02_exch = pieces[(0, 'exch', 'I')] + pieces[(2, 'exch', 'II')]
actual_02_exch = integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r2**2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                 integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r2**2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
pred_02_exch = simplify(pred_02_exch)
actual_02_exch = simplify(actual_02_exch)
print(f"\nBipolar(0,2)_exch:")
print(f"  Predicted (M^0_exch_I + M^2_exch_II): {float(pred_02_exch):+.6e}")
print(f"  Actual                              : {float(actual_02_exch):+.6e}")
print(f"  Difference                          : {simplify(pred_02_exch - actual_02_exch)}")

# Bipolar (1, 1) direct: kernel r_1 r_2 / r_>^5 = M^1 region I + M^1 region II = M^1_total!
# (symmetric: r_1 r_2 is same in both regions)
pred_11_dir = pieces[(1, 'dir', 'I')] + pieces[(1, 'dir', 'II')]
actual_11_dir = integrate(integrate(d_1s_sq * d_2p_sq * r1 * r2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                integrate(integrate(d_1s_sq * d_2p_sq * r1 * r2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
pred_11_dir = simplify(pred_11_dir)
actual_11_dir = simplify(actual_11_dir)
print(f"\nBipolar(1,1)_dir:")
print(f"  Predicted (M^1_dir total): {float(pred_11_dir):+.6e}")
print(f"  Actual                   : {float(actual_11_dir):+.6e}")
print(f"  Difference               : {simplify(pred_11_dir - actual_11_dir)}")

pred_11_exch = pieces[(1, 'exch', 'I')] + pieces[(1, 'exch', 'II')]
actual_11_exch = integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r1 * r2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                 integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r1 * r2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
pred_11_exch = simplify(pred_11_exch)
actual_11_exch = simplify(actual_11_exch)
print(f"\nBipolar(1,1)_exch:")
print(f"  Predicted (M^1_exch total): {float(pred_11_exch):+.6e}")
print(f"  Actual                    : {float(actual_11_exch):+.6e}")
print(f"  Difference                : {simplify(pred_11_exch - actual_11_exch)}")

# Bipolar (2, 0) direct: kernel r_1^2 / r_>^5 = M^2 region I + M^0 region II
pred_20_dir = pieces[(2, 'dir', 'I')] + pieces[(0, 'dir', 'II')]
actual_20_dir = integrate(integrate(d_1s_sq * d_2p_sq * r1**2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                integrate(integrate(d_1s_sq * d_2p_sq * r1**2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
pred_20_dir = simplify(pred_20_dir)
actual_20_dir = simplify(actual_20_dir)
print(f"\nBipolar(2,0)_dir:")
print(f"  Predicted (M^2_dir_I + M^0_dir_II): {float(pred_20_dir):+.6e}")
print(f"  Actual                            : {float(actual_20_dir):+.6e}")
print(f"  Difference                        : {simplify(pred_20_dir - actual_20_dir)}")

pred_20_exch = pieces[(2, 'exch', 'I')] + pieces[(0, 'exch', 'II')]
actual_20_exch = integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r1**2 / r2**5, (r1, 0, r2)), (r2, 0, oo)) + \
                 integrate(integrate(d_1s2p_r1 * d_2p1s_r2 * r1**2 / r1**5, (r1, r2, oo)), (r2, 0, oo))
pred_20_exch = simplify(pred_20_exch)
actual_20_exch = simplify(actual_20_exch)
print(f"\nBipolar(2,0)_exch:")
print(f"  Predicted (M^2_exch_I + M^0_exch_II): {float(pred_20_exch):+.6e}")
print(f"  Actual                              : {float(actual_20_exch):+.6e}")
print(f"  Difference                          : {simplify(pred_20_exch - actual_20_exch)}")
