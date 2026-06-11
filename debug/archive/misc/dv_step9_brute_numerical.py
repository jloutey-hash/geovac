"""Sprint 5 DV Step 9: Brute-force numerical evaluation of <^3P_1 | H_SS | ^3P_1>

Compute A_SS = <^3P_1, M_J = 1 | H_SS | ^3P_1, M_J = 1> directly by 6-dim quadrature
(3 radial + 3 angular dims after using M_J and radial integration to reduce).

Strategy:
  - Use the ROTATIONAL symmetry to simplify: only M_J-diagonal matrix elements
    of a scalar operator are nonzero.
  - Choose M_J = 1 for cleanest CG coefficients.
  - Use the simplest ALIGNED-SPIN state: M_S = 1 (both spins up), M_L = 0 (spatial symmetry
    breaks at 2p_0).

|^3P_1, M_J = 1> has two components:
  a * |L=1 M_L=1, S=1 M_S=0> + b * |L=1 M_L=0, S=1 M_S=1>
with a = 1/sqrt(2), b = -1/sqrt(2).

We work with the (M_L=0, M_S=1) piece (b component), which has spins up on both electrons
(|up up>) and spatial |L=1 M_L=0> = (1/sqrt 2)[|1s(1) 2p_0(2)> - |2p_0(1) 1s(2)>].

H_SS has zero spin matrix element on |up up> ONLY for rank-2 tensor components Q=0:
  <up up | [s_1 (x) s_2]^(2)_0 | up up> = (nonzero by explicit computation)

Actually for triplet M_S = +1:
  [s_1 x s_2]^(2)_0 on |up up> = Y^2_0-like component, which gives
  <up up | [s_1 x s_2]^(2)_0 | up up> = ?
  Using s_1 · s_2 = (S^2 - s_1^2 - s_2^2)/2 = (2 - 3/4 - 3/4)/2 = 1/4
  For triplet, s_1 · s_2 = 1/4. That's scalar, but
  [s_1 (x) s_2]^(2) = coupled rank-2.
  By Wigner-Eckart: <S=1 M_S=1 | [s_1 x s_2]^(2)_Q | S=1 M_S'=1> nonzero only for Q=0.
  The value is <S=1 1; 2 0 | S=1 1> * <S=1 || [s1 x s2]^(2) || S=1> / sqrt(2S+1)
     = <1 1; 2 0 | 1 1> * (sqrt(5)/2) / sqrt(3)
  <1 1; 2 0 | 1 1> = sqrt(1/10)  (from CG tables)
  => ME = sqrt(1/10) * sqrt(5)/2 / sqrt(3) = sqrt(1/30) * 1/2

Let me just directly enumerate ALL contributions to A_SS numerically and avoid
convention issues.
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import mpmath as mp
mp.mp.dps = 30

from geovac.breit_integrals import breit_ss_radial

# ============================================================
# Reference values: Drake's A_SS at Z=1, no alpha^2
# ============================================================
import sympy as sp
from sympy import Rational, sqrt, Integer, simplify

M0_dir = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 0, Z=1))
M0_exc = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 0, Z=1))
M1_dir = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 1, Z=1))
M1_exc = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 1, Z=1))
M2_dir = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1))
M2_exc = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1))

A_SS_drake_f = 3/50 * M2_dir - 2/5 * M2_exc
print(f"Drake A_SS (no alpha^2): {A_SS_drake_f:+.6e}")
print(f"  = (3/50) M^2_dir - (2/5) M^2_exch")
print()
print(f"M^K values (Z=1):")
print(f"  M^0_dir  = {M0_dir:+.6e}")
print(f"  M^1_dir  = {M1_dir:+.6e}")
print(f"  M^2_dir  = {M2_dir:+.6e}")
print(f"  M^0_exch = {M0_exc:+.6e}")
print(f"  M^1_exch = {M1_exc:+.6e}")
print(f"  M^2_exch = {M2_exc:+.6e}")

# ============================================================
# BRUTE FORCE: compute <^3P_1 | H_SS | ^3P_1> at Z = 1 using numerical 6D quadrature.
# ============================================================
#
# The operator H_SS = -alpha^2 * { (s_1 . s_2) - 3 (s_1 . hat r_12)(s_2 . hat r_12) } / r_12^3
# We factor out -alpha^2 and compute <H_SS/(-alpha^2)> = A_SS / (-alpha^2).
#
# Using M_J = 1 with the M_L=0, M_S=1 component ONLY (the b=-1/sqrt(2) part):
# The FULL ^3P_1 state has cross-coupling between the two components. But since
# H_SS is spin-dependent, it can mix M_S=0 <-> M_S=1 via rank-2 tensor.
#
# To ISOLATE the direct matrix element involving only <M_L=0, M_S=1 | O | M_L=0, M_S=1>,
# take the J=1 M_J=1 expectation value with BOTH components.
#
# Actually, simpler: use M_J = 2 where |^3P_1, 2> does not exist (J=1 has |M|≤1).
# So M_J=2 gives no |^3P_1>. Let me stick with M_J=1.
#
# Alternative: compute <^3P_2, M_J=2 | H_SS | ^3P_2, M_J=2>.
#   f_SS(2) = -1/5. So A_SS = -5 * <^3P_2, 2 | H_SS | ^3P_2, 2>.
# And |^3P_2, 2> = |L=1 M_L=1, S=1 M_S=1> with CG = 1 (single component).
#
# <^3P_2, 2 | H_SS | ^3P_2, 2> = <|1s(1) 2p_+1(2) up up> - |2p_+1(1) 1s(2) up up>|  / sqrt 2
#                              * H_SS * same
# This has ONLY ONE SPATIAL STATE (M_L=1) and ONE SPIN STATE (M_S=1).
# The spin matrix element of H_SS on |up up> is with M_S unchanged — which requires Q=0.
#
# This simplifies enormously! Let's compute <^3P_2, M=2 | H_SS | ^3P_2, M=2>.

# |^3P_2, M=2> = (1/sqrt 2) [|1s(1) 2p_+1(2) up up> - |2p_+1(1) 1s(2) up up>]
#   (antisymmetric spatial × symmetric spin: triplet)
#
# The matrix element:
#   <^3P_2, 2 | H_SS | ^3P_2, 2> = direct - exchange
#   direct = <1s(1) 2p_+1(2) up up | H_SS | 1s(1) 2p_+1(2) up up>
#   exchange = <1s(1) 2p_+1(2) up up | H_SS | 2p_+1(1) 1s(2) up up>
#
# For |up up> both bra and ket, only Q=0 contribution of [s_1 x s_2]^(2)_Q survives.
#   <up up | [s_1 x s_2]^(2)_0 | up up> = ?
#
# By direct computation: [s_1 x s_2]^(2)_0 = sum_{q1, q2} <1 q1; 1 q2 | 2 0> s_{q1}(1) s_{q2}(2)
# On |up up> (m1 = m2 = +1/2):
#   Only q1 + q2 = 0 contributes (conservation), so q1 = q2 = 0 or q1=+1, q2=-1 or q1=-1, q2=+1
#   q1=q2=0: s_0(1)|up> = (1/2)|up>, s_0(2)|up>=(1/2)|up>; product = 1/4
#   q1=+1: s_+(1)|up>=0 (|up> is max), so 0.
#   q1=-1: s_-(1)|up> = something on |up> down in m. But s_-1 = (s_x-i s_y)/sqrt(2) lowers m,
#          so s_-1 |up> = (1/sqrt 2)(|down>). So we'd map |up up> -> (1/sqrt 2)|down up>.
#          This ISN'T |up up> so <up up | q1=-1, q2=+1 | up up> = 0.
#   So only q1=q2=0 contributes, with CG <1 0; 1 0 | 2 0> = sqrt(2/3).
#   => <up up | [s_1 x s_2]^(2)_0 | up up> = sqrt(2/3) * 1/4 = sqrt(2/3)/4 = 1/(2 sqrt 6)

# Y^(2)_0(hat r_12) = sqrt(5/(4 pi)) P_2(cos theta_12)
# = sqrt(5/(4 pi)) * (3 cos^2 theta_12 - 1)/2

# So the matrix element:
# <^3P_2, 2 | H_SS | ^3P_2, 2>
#   = -alpha^2 * sqrt(24 pi / 5) * (-1)^0 * <up up | [s_1 x s_2]^(2)_0 | up up>
#      * <spatial ket | Y^(2)_0(hat r_12)/r_12^3 | spatial ket>
# The spatial ket is (1/sqrt 2) [|1s 2p_+1> - |2p_+1 1s>].
# <spatial|O|spatial> = (direct - exchange) with O = Y^(2)_0(hat r_12)/r_12^3
# (note the minus sign for triplet spatial antisymmetry).

# Radial + Angular matrix elements of Y^(2)_0(hat r_12)/r_12^3:

# Direct: <1s(1) 2p_+1(2) | Y^(2)_0(hat r_12)/r_12^3 | 1s(1) 2p_+1(2)>
# Exchange: <1s(1) 2p_+1(2) | Y^(2)_0(hat r_12)/r_12^3 | 2p_+1(1) 1s(2)>

# Numerically: compute these via 6D integral.

# Simplification: Y^(2)_0 = sqrt(5/(4 pi)) P_2(cos theta_12).
# 1/r_12^3 * P_2(cos theta_12) has bipolar expansion that we can compute.

# Convergent formulation: introduce regularizer by computing 1/|r_1 - r_2|^3 * P_2(cos theta_12).
# This integral may have singularity at r_1 = r_2; check integrability.

# The integral
#  <1s(1) 2p_+1(2) | Y^(2)_0(hat r_12)/r_12^3 | 1s(1) 2p_+1(2)>
# = int int int int int int R_1s(r_1)^2 R_2p(r_2)^2 |Y_1,+1|^2 |Y_0,0|^2 ...
# = int int R_1s(r_1)^2 R_2p(r_2)^2 * r_1^2 r_2^2 * A(r_1, r_2) dr_1 dr_2
# where A(r_1, r_2) is the angular-integrated part:
#  A(r_1, r_2) = int int dΩ_1 dΩ_2 |Y_0,0(Ω_1)|^2 |Y_1,+1(Ω_2)|^2 Y^2_0(hat r_12) / |r_1 - r_2|^3

# Let's compute A for specific (r_1, r_2) numerically, then integrate radially.

import scipy.special

def compute_A(r1_val, r2_val, ngrid=40):
    """Compute A(r_1, r_2) for the DIRECT channel."""
    # A(r_1, r_2) = int_{0}^{pi} sin theta_1 dtheta_1 int_0^{2pi} dphi_1
    #               int_0^{pi} sin theta_2 dtheta_2 int_0^{2pi} dphi_2
    #               |Y_{0,0}(theta_1, phi_1)|^2 * |Y_{1,+1}(theta_2, phi_2)|^2
    #               * Y^2_0(theta_12, phi_12)
    #               / |r_1 - r_2|^3
    #
    # With |Y_{0,0}|^2 = 1/(4pi) (const), we can pull out factor 1/(4pi) and integrate over Ω_1.
    # Then integral over Ω_1 is trivial: 4pi * (const).
    # But wait, the Y^2_0(hat r_12)/r_12^3 DEPENDS on Ω_1 and Ω_2 via hat r_12 = (r_1 - r_2)/|r_1 - r_2|.
    # So we can't pull out.
    #
    # Parameterize: r_1 = (r1_val, theta_1, phi_1), r_2 = (r2_val, theta_2, phi_2).
    # r_12 = r_1 - r_2.  |r_12|^2 = r_1^2 + r_2^2 - 2 r_1 r_2 cos theta_{12} with
    #   cos theta_{12} = sin theta_1 sin theta_2 cos(phi_1 - phi_2) + cos theta_1 cos theta_2.
    # hat r_12 = r_{12} / |r_12|.

    # This is a 4D integral (theta_1, phi_1, theta_2, phi_2) that we evaluate via Gauss-Legendre.
    pass


# This gets complex fast. Let me use mpmath's quad for the reduced integral.
# Actually for a DIRECT channel, we can reduce to a 3D angular integral by fixing the
# orientation of e1 (say theta_1 = 0, phi_1 = 0) and rotating everything.

# Tip: for a SCALAR operator integrated over a rotationally-invariant state, we can
# simplify. But here |Y_1,+1|^2 is NOT rotationally invariant (preferred axis).

# Simpler approach: work entirely with the partial-wave expansion of 1/r_12^3 * Y^2_0.
# Expand in Legendre polynomials of cos theta_12:
#   Y^2_0(hat r_12)/r_12^3 = sqrt(5/(4 pi)) P_2(cos theta_12) / r_12^3
# and use
#   P_L(cos theta_12) = (4 pi / (2 L + 1)) sum_{M} Y^L_M(Omega_1)* Y^L_M(Omega_2)

# What's 1/r_12^3 expanded in Legendre? It's given by the Sack expansion:
#   1/r_12^n = sum_L A_L^{(n)}(r_1, r_2) P_L(cos theta_12)
# with A_L^{(n)}(r_1, r_2) specific functions.

# For n = 3 (our case):
#   A_L^{(3)}(r_1, r_2) = ???

# Actually Sack 1964 gives for 1/r_12^{2 nu}:
#   1/r_12^{2 nu} = sum_{l=0}^{inf} f_l^(nu)(r_1, r_2) P_l(cos theta_12)
# with f_l^(nu) a specific radial function involving Gegenbauer polynomials.

# The formula for f_l^(n) (where n = 1 is Coulomb, n = 3 is our case):
# f_l^(n)(r_1, r_2) = [gamma(n-l)!/(2l+1)!] * ... (non-trivial)

# This is getting too technical. Let me abandon the "derive it from scratch" approach
# and instead do a CROSS-CHECK: use the Drake formula to predict the physical result
# numerically and verify consistency.

# INSTEAD OF DERIVING: compute the matrix element via 6D scipy integration (which
# is slow but direct) and PSLQ-identify against Drake's M^K basis.

print("\nDirect 6D numerical integration implementation:")
print("This requires careful handling of the 1/r_12^3 singularity; deferring.")
print()

# ============================================================
# ALTERNATIVE: solve the linear system to derive (3/50, -2/5).
#
# If Drake's formula is
#    A_SS = c_dir * M^2_dir + c_exch * M^2_exch
# with the unknown (c_dir, c_exch), and we have the experimental A_SS-value
# reproduced to 0.20% by (c_dir, c_exch) = (3/50, -2/5):
#
# The rational structure can be reverse-engineered from WHICH bipolar channels
# contribute. Let me list the ONLY angular-allowed channels for He (1s, 2p) ^3P:
#   Direct  (K=2): (k_1=0, k_2=2)
#   Exchange (K=2): (k_1=1, k_2=1)
# No other channels contribute due to angular selection (parity + triangle).
#
# The direct bipolar integral (k_1=0, k_2=2) in terms of Drake's M^K:
#   I_bipolar(0,2,K=2) direct = M^0_dir_I + M^2_dir_II
# This is NOT a multiple of M^2_dir (nor any M^K total alone).
# So the Drake formula's (3/50 M^2_dir - 2/5 M^2_exch) is not the "natural" output
# of the bipolar expansion.
#
# CONCLUSION: the Drake formula uses M^K integrals that are REDEFINED from the naive
# r_<^K/r_>^{K+3} form. Specifically, Drake likely uses modified radial integrals that
# are region-specific, or uses a different angular integral definition.
# ============================================================

# Let me look at what Drake 1971 actually defines as M^K:
# In Drake 1971 (Phys Rev A 3, 908), the formula is in §III. For the Breit operators,
# the integrals are defined by specific radial amplitudes. From memory, Drake's M^K
# convention is:
#   M^K = (radial integral) with specific normalization conventions including angular factors.
#
# Rather than reverse-engineering from partial information, let me COMPUTE A_SS from
# the Racah-coupled formula using the bipolar (0, 2) direct and (1, 1) exchange
# radial integrals, and THEN extract what the formula looks like in Drake's basis.

# ============================================================
# Compute A_SS symbolically using the bipolar channel integrals
# ============================================================

# <^3P_1 | H_SS | ^3P_1> = A_SS * f_SS(1) = A_SS (at J=1, f_SS(1) = 1)
#
# With:
#   A_SS = (-alpha^2) sqrt(6) * (-1)^(L+S+J) 6j{LSJ;SLK} * red_spin *
#             [spatial_red_dir * beta(0,2,2) * I_bipolar(0,2,2)_dir
#              - spatial_red_exch * beta(1,1,2) * I_bipolar(1,1,2)_exch]
#
# Known quantities (at Z=1, no alpha^2):
#   Bipolar(0, 2, 2) direct integral I_d = M^0_dir_I + M^2_dir_II
#     Region I of M^0_dir has value: (from Step 7) I_d_I^0 = ... (symbolic value)
#     Region II of M^2_dir has value: (from Step 7) I_d_II^2 = ...
#   Bipolar(1, 1, 2) exchange integral I_e = M^1_exch_total = M^1_exch
#
# J-factor J_factor = (-1)^(L+S+J) 6j = -1/6 at J=1, L=S=1, K=2
# spin_red = sqrt 5 / 2
# spatial_red_dir = -sqrt(30)/5  (k_1=0, k_2=2 channel)
# spatial_red_exch = -sqrt(5)/3  (k_1=1, k_2=1 channel)
#
# beta(0, 2, 2), beta(1, 1, 2) are the UNKNOWN bipolar prefactors.
#
# The identity A_SS(computed) == A_SS(Drake) gives one equation in TWO unknowns
# beta(0,2,2) and beta(1,1,2). We need ANOTHER equation. Let me use a different
# state or channel to get another equation.
#
# Option: compute <^3P_2 | H_SS | ^3P_2> which gives f_SS(2) A_SS = -A_SS/5.
# This uses the SAME bipolar channels (0,2) direct and (1,1) exchange, just with
# the J-factor changed from -1/6 to (-1)^(1+1+2) 6j{1,1,2; 1,1,2} = ...
# This is a COPY of the same equation (up to the J-factor). Doesn't help.

# Alternative: use a different MULTIPOLE K (K=0 direct-Coulomb, K=1 exchange-Breit, etc.)
# to get another equation. But these involve DIFFERENT (k_1, k_2) channels.
#
# For K = 0 (Coulomb-like rank-0 tensor operator Y^(0)(hat r_12)/r_12^1 = 1/r_12):
#   Direct: (k_1=0, k_2=0) single channel
#   Exchange: (k_1=1, k_2=1) single channel (required by angular selection <0|C^1|1>, <1|C^1|0>)
#
# The K=0 matrix element is the FAMILIAR Coulomb matrix element:
#   <^3P_J | 1/r_12 | ^3P_J> = F^0(1s, 1s; 2p, 2p) - (1/3) G^1(1s, 2p; 2p, 1s)
# This is a STANDARD result:
#   Coulomb direct = F^0 = R^0(1s, 2p; 1s, 2p) = M^0_dir (in Coulomb, not Breit, version)
#   Coulomb exchange = G^1 = R^1(1s, 2p; 2p, 1s) = M^1_exch (coulomb)
#   Total: <^3P_J | 1/r_12 | ^3P_J> = F^0 - (1/3) G^1 (rank-0 scalar, J-independent)

# The coefficient is 1 (direct) and -1/3 (exchange).
#
# If we DERIVE THE SAME COEFFICIENTS via the same bipolar formula we used for SS,
# we can calibrate beta(0, 0, 0) and beta(1, 1, 0).

# Formula for K = 0:
#   <^3P_J | O^(0)(space) | ^3P_J> = (-1)^{L+S+J} 6j{LSJ;SL 0} * <L=1||O||L=1> (scalar spin trivial)
#   6j{L,S,J; S,L,0} = 1 / sqrt((2L+1)(2S+1)) * delta... — wait this has special form.
#   For k=0: 6j{L S J; S L 0} = (-1)^{L+S+J} / sqrt((2L+1)(2S+1))  (from Edmonds 6.3.2)
#   So (-1)^{L+S+J} * 6j = 1/sqrt((2L+1)(2S+1))

# For the COULOMB (rank-0 operator, so spin is identity rank-0),
#   <L=1 || O^(0) || L=1> = sum_{(k1, k2)} ...

# This is getting complicated. Let me ONE more time try a simpler approach:
# just solve the equation A_SS_numerical = c_d * M^2_dir + c_e * M^2_exch
# and see if (3/50, -2/5) come out.

# We need ONE equation with TWO unknowns (c_d, c_e). So we need another piece.
# Is there a way to separate direct from exchange contributions numerically?

# Yes: the DIRECT contribution alone is (from the Racah formula):
#   A_SS_direct = -alpha^2 sqrt(6) * J_factor * red_spin * spatial_red_dir *
#                   beta(0, 2, 2) * I_bipolar(0,2,2)_dir
# and the EXCHANGE contribution is similar with minus sign.

# Drake's formula splits A_SS = c_d * M^2_dir + c_e * M^2_exch, where c_d is from direct
# and c_e is from exchange. These should correspond to the direct and exchange
# contributions above, respectively.

# So if we can WRITE: A_SS_direct = (3/50) M^2_dir, what does it say about bipolar(0,2)_dir?
# I_bipolar(0,2)_dir = M^0_dir_Region_I + M^2_dir_Region_II (not just M^2_dir total!)

# This discrepancy says that our CHANNEL-AMPLITUDE IDENTITY is not correct — the
# bipolar expansion has terms we're missing.

# ============================================================
# REAL RESOLUTION: Drake 1971 uses a DIFFERENT definition of M^K than what
# `breit_ss_radial` computes!
# ============================================================

# Drake 1971, Eq. (16) defines the radial integrals as:
#   M^k(a, b; c, d) = int_0^inf dr_1 int_0^{r_1} dr_2 [DENS_1 DENS_2] (r_2^k / r_1^{k+3})
#                  + int_0^inf dr_1 int_{r_1}^inf dr_2 [DENS_1 DENS_2] (r_1^k / r_2^{k+3})
# The same as Breit: r_<^k / r_>^{k+3}. So no difference here.
#
# BUT, the DECOMPOSITION formula (A_SS = 3/50 M^2_dir - 2/5 M^2_exch) is NOT a
# "bipolar expansion" decomposition — it's a SELECTION of specific M^K values
# that together give the correct physical A_SS. There may be "implicit" cancellations
# in the bipolar expansion that aren't visible.
#
# Most likely: the bipolar expansion has MORE than just k_1 + k_2 = K terms.
# The higher-order bipolar terms (k_1 + k_2 = K + 2, K + 4, ...) contribute
# too, and their sum matches exactly (3/50 M^2_dir - 2/5 M^2_exch).

# Let me VERIFY this hypothesis: compute A_SS numerically using:
#   A_SS_bipolar_leading = (bipolar(0,2,2) direct minus bipolar(1,1,2) exchange)
#                           with appropriate prefactors.

# And cross-check against A_SS_drake = 3/50 M^2_dir - 2/5 M^2_exch.

# To do this, we need beta(0, 2, 2) and beta(1, 1, 2).

# For the calibration, I'll use the CASE K = 0 where Drake's formula is known to be:
#   A_00 = F^0 - (1/3) G^1  (coulomb direct minus 1/3 times exchange)
#   (hmm this is for a different operator)

# Cleaner CALIBRATION: use the K = 0 IDENTITY for the BREIT kernel.
# For K = 0 operator Y^0(hat r_12)/r_12^1 (trivial angular, just 1/r_12):
#   This ISN'T a Breit-Pauli operator. Different operator.

# OK let me try yet another angle. COMPUTE the full 6D integral numerically via mpmath
# and just see the number.

print("End of analytic step 9. Moving to direct 6D numerical integration in Step 10.")
