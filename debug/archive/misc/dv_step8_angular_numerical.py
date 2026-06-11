"""Sprint 5 DV Step 8: Compute <^3P_1, M=1 | H_SS | ^3P_1, M=1> = A_SS
directly via (angular_me * radial_integral) in the position-space view,
identifying the specific coefficients of M^K_dir and M^K_exch.

The SS operator (Breit-Pauli retarded form):
  H_SS = -alpha^2 * [ (s_1 · s_2) / r_12^3
                      - 3 (s_1 · hat r_12)(s_2 · hat r_12) / r_12^3 ]

In tensor (rank-2) form:
  H_SS = -alpha^2 * sqrt(24 pi / 5) sum_q (-1)^q [s_1 (x) s_2]^(2)_q
           * Y^(2)_{-q}(hat r_12) / r_12^3

The He ^3P_J, M_J = J matrix element at J=1 equals A_SS * f_SS(1) = A_SS
(since f_SS(1) = 1). So we can compute A_SS directly from M_J = 1.

STEP: Build the ^3P_1, M=1 wavefunction and apply H_SS.
  |^3P_1, M_J=1> = sum_{M_L + M_S = 1} <1 M_L; 1 M_S | 1 1> |1s 2p; M_L> |S=1, M_S>
  where |1s 2p; M_L> = |1s m_a=0> |2p m_b=M_L> for M_L ∈ {-1, 0, 1} (since
  the ANTISYMMETRIZED spatial state with 1s ⊗ 2p + 2p ⊗ 1s, but for triplet
  we need spatial ANTISYMMETRIC: (1/sqrt 2)(|1s 2p> - |2p 1s>))

CG coefficients:
  <1 1; 1 0 | 1 1> = 1/sqrt(2)
  <1 0; 1 1 | 1 1> = -1/sqrt(2)
  <1 -1; 1 0 | 1 1> = 0 (since 1 0; 1 0 -> M_J=0 only)

Actually, for S=1, M_S=1 (triplet with spins aligned) all spatial M_L's allowed. M_L+M_S=1 with M_S in {-1,0,1}:
  (M_L=2, M_S=-1): excluded (L=1 so M_L ≤ 1)
  (M_L=1, M_S=0): CG = <1 1; 1 0 | 1 1>
  (M_L=0, M_S=1): CG = <1 0; 1 1 | 1 1>

From Varshalovich tables:
  <1 1; 1 0 | 1 1> = 1/sqrt(2)
  <1 0; 1 1 | 1 1> = -1/sqrt(2)

STRATEGY:
=========
For each (M_L, M_S) pair with M_L + M_S = 1, build the spatial triplet L=1, M_L
state and the S=1, M_S spin state. Apply H_SS using the Y^2_{-q} operator; sum over q.

Spatial part: (1/sqrt(2))[1s(1) 2p_{M_L}(2) - 2p_{M_L}(1) 1s(2)].
Compute <spatial | Y^2_{-q}(hat r_12) / r_12^3 | spatial> directly.
(Note: the spatial matrix element may be complex in terms of Slater integrals.)

PRACTICAL IMPLEMENTATION:
=========================
Evaluate matrix elements as nested integrals over r_1, r_2, theta_1, theta_2,
phi_1, phi_2 (6-dimensional). Use numerical integration with mpmath for
high precision, then PSLQ-identify the result in terms of Drake M^K.

For efficiency, use the multipole expansion of 1/r_12^3:
  1/r_12^3 = (partial derivative of 1/r_12 with respect to something)^2/?
  Not clean. Instead use:
     1/r_12^n P_L(cos theta_12) = sum_l f_L^{n, l}(r_1, r_2) P_l(cos theta_1) P_l(cos theta_2) ?
     (these are Sack's radial projections). Not separable.

INSTEAD: use the ANGULAR factorization only on the ANGULAR part. Integrate
radially via sympy.

The angular integral of Y^K_Q(hat r_12) between hydrogenic states (l=0, 1 only
in our case):
  <0, 0; 1, m | Y^K_Q(hat r_12) | 0, 0; 1, m'>_{angular}  (direct)
  <0, 0; 1, m | Y^K_Q(hat r_12) | 1, m'; 0, 0>_{angular}  (exchange)

USE THE BIPOLAR EXPANSION WITH SPECIFIC (but verified) COEFFICIENT:

Rose's "Elementary Theory of Angular Momentum" (1957), Eq. 4.30:
  Y^K_M(hat r_12) = sum_{l_1, l_2} (4pi / (2K+1))^{-1/2}
                     * (2l_1+1) (2l_2+1) * <l_1 0 l_2 0 | K 0>
                     * (r_1^{l_1} r_2^{l_2} / r_12^{l_1 + l_2 + 1})
                     * [Y^{l_1}(1) x Y^{l_2}(2)]^K_M  / (2K+1)

I am confusing myself. Let me just DIRECTLY COMPUTE the matrix element
numerically via 6-dimensional integration and see what comes out.
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
from scipy import integrate, special
import mpmath as mp
mp.mp.dps = 30

# ============================================================
# Step A: Compute numerically <^3P_1, M=1 | H_SS | ^3P_1, M=1> at Z=1
# ============================================================
#
# H_SS = -alpha^2 * [(s_1 . s_2) - 3 (s_1 . hat r_12)(s_2 . hat r_12)] / r_12^3
#
# In the TENSOR form, the operator is
#
#  O(space) = Y^(2)_Q(hat r_12) / r_12^3
#  O(spin)  = [s_1 (x) s_2]^(2)_Q
#
# Combined as scalar product:
#  H_SS = -alpha^2 sqrt(24pi/5) sum_q (-1)^q [s_1 x s_2]^(2)_q * Y^(2)_{-q}(hat r_12) / r_12^3
#
# At M_J=J=1 (of ^3P_1), the matrix element A_SS = <...|H_SS|...> is a scalar.
#
# Let me compute <1s m_a=0; 2p m_b=1 | Y^(2)_Q(hat r_12)/r_12^3 | 1s m_a=0; 2p m_b=1 >
# with specific Q values.
#
# By spherical symmetry, only Q values where M_L = M_L' are nonzero.
# Here M_L = M_L' = 0+1 = 1, so Q = 0.
#
# Similarly for exchange:
#   <1s 0; 2p 1 | Y^(2)_0(hat r_12)/r_12^3 | 2p 1; 1s 0>
# with M_L = 1 on both sides.
#
# Since these are scalar contributions (Q=0), we can work with real Y^(2)_0 or just
# the P_2(cos theta_12) piece. Since Y^2_0 = sqrt(5/(4pi)) P_2(cos theta_12), we
# have
#
#  Y^(2)_0(hat r_12) / r_12^3  =  sqrt(5/(4pi)) * P_2(cos theta_12) / r_12^3

# So the RADIAL-ANGULAR COMBINED integral we need is:
#
#  I_dir(P2) = <1s(1) 2p_1(2) | P_2(cos theta_12) / r_12^3 | 1s(1) 2p_1(2)>
#  I_exc(P2) = <1s(1) 2p_1(2) | P_2(cos theta_12) / r_12^3 | 2p_1(1) 1s(2)>
#
# And then A_SS (direct part) = -alpha^2 * sqrt(24 pi/5) * sqrt(5/(4pi)) * (prefactor q=0 spin)
#                              * (I_dir - I_exc)
#                              = -alpha^2 * sqrt(6) * (prefactor q=0 spin) * (I_dir - I_exc)
#
# with "prefactor q=0 spin" the matrix element of [s_1 x s_2]^(2)_0 between
# the |S=1, M_S=0> state (since M_S was 0 to give M_L=M_J=1).

# Actually wait: A_SS is the coefficient such that <M_J=1|H|M_J=1> = A_SS * f_SS(1) = A_SS.
# Let me be more careful about the CG expansion of |^3P_1 M_J=1>.

# |^3P_1, M_J=1> = (1/sqrt 2) |L=1, M_L=1; S=1, M_S=0> + (-1/sqrt 2) |L=1, M_L=0; S=1, M_S=1>
#
# Using the ANTISYMMETRIZED spatial triplet |L=1, M_L>:
# |L=1 M_L> = (1/sqrt 2) [|1s m_a=0; 2p m_b=M_L> - |2p m_b=M_L; 1s m_a=0>]
# (where the ordering is (electron 1)(electron 2))
#
# But for M_L = 1: this is (1/sqrt 2)[|1s, 2p_1> - |2p_1, 1s>]
# For M_L = 0: (1/sqrt 2)[|1s, 2p_0> - |2p_0, 1s>]

# The spin part: |S=1, M_S>:
#   |S=1, M_S=0> = (1/sqrt 2)(|up down> + |down up>)
#   |S=1, M_S=1> = |up up>

# Note: In Condon-Shortley convention, spin-coupled states are:
#   |S=1, M_S=0> = (1/sqrt(2)) [|1/2, 1/2> |1/2, -1/2> + |1/2, -1/2> |1/2, 1/2>]

# Since H_SS is spin-DEPENDENT (via the [s_1 x s_2]^(2)), it couples different M_S
# sectors. But diagonal ME (M_S unchanged): <M_S=0| [s_1 x s_2]^(2)_0 |M_S=0> is nonzero.
# Cross-ME <M_S=0| [s_1 x s_2]^(2)_{+1} |M_S=1> is also nonzero.

# For M_J=1 evaluation, both (M_L=1, M_S=0) and (M_L=0, M_S=1) can mix via
# Y^2_Q with Q = +/-1.
#
# Given the CG coefficients (a=-1/sqrt(2), b=+1/sqrt(2)) for (M_L=1, M_S=0) and
# (M_L=0, M_S=1), the matrix element at M_J=1 is:
#  <1 1| H |1 1> = |a|^2 H_{MM0,MM0} + |b|^2 H_{ML0MS1, ML0MS1}
#                   + a^* b H_{ML1MS0, ML0MS1} + a b^* H_{ML0MS1, ML1MS0}
# (with ML1MS0 meaning M_L=1, M_S=0 etc.)
#
# With a = -1/sqrt(2), b = 1/sqrt(2):
#   |a|^2 = |b|^2 = 1/2
#   a^* b = -1/2
#   a b^* = -1/2
# So = (1/2) [H_{M_L1 M_S0} + H_{M_L0 M_S1} - H_{off} - H_{off^*}]

# By rotational invariance (H is a scalar), all M_J-diagonal elements of H for the same J
# are equal: H|^3P_1, M_J=0> = H|^3P_1, M_J=1> = same eigenvalue.
# So I can choose M_J=0 where CG is simpler... but let's stick with M_J=1.

# Actually easier: evaluate <M_J = 0 | H | M_J = 0> for J=1.
# |^3P_1, M_J=0> = sum_{M_L+M_S=0} CG |L=1 M_L; S=1 M_S>
# CG coefficients: <1 M_L; 1 M_S | 1 0>
#   M_L = -1, M_S = 1: <1 -1; 1 1 | 1 0> = 1/sqrt 2
#   M_L = 0, M_S = 0: <1 0; 1 0 | 1 0> = 0  (this CG is zero!)
#   M_L = 1, M_S = -1: <1 1; 1 -1 | 1 0> = -1/sqrt 2

# So |^3P_1, 0> = (1/sqrt 2) (|L=1 -1; S=1 1> - |L=1 1; S=1 -1>)

# For this state, M_L+M_S=0, and we can compute each contribution to <H|>.

# OK let me just CODE IT NUMERICALLY and see what A_SS at Z=1 comes out to.

# First, some numerical references:
from geovac.breit_integrals import breit_ss_radial
from fractions import Fraction
import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, pi

M2_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1)
M2_exc = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1)
M1_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 1, Z=1)
M1_exc = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 1, Z=1)
M0_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 0, Z=1)
M0_exc = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 0, Z=1)

A_SS_drake = Rational(3, 50) * M2_dir - Rational(2, 5) * M2_exc
print(f"A_SS (Drake, Z=1, no alpha^2):")
print(f"  Numerical: {float(A_SS_drake):+.6e}")
print(f"  = 3/50 * M^2_dir - 2/5 * M^2_exch")
print(f"  M^2_dir  = {float(M2_dir):+.6e}")
print(f"  M^2_exch = {float(M2_exc):+.6e}")

# ============================================================
# Step B: numerical 6D integration of the angular matrix element
# ============================================================

import numpy as np
# At M_J = 1 (chosen to be simpler with M_L=0, M_S=1 term dominant):
# |^3P_1, M_J=1> = (1/sqrt 2) (|L=1 M_L=1, S=1 M_S=0> - |L=1 M_L=0, S=1 M_S=1>)
# Wait, let me redo the CG. <1 1; 1 0 | 1 1> from Varshalovich:
#   for CG <j_1 m_1; j_2 m_2 | j m>: with j_1=1, j_2=1, m_1=1, m_2=0, j=1, m=1
#   This is <1, 1; 1, 0 | 1, 1> = 1/sqrt(2).
#   Similarly <1, 0; 1, 1 | 1, 1> = -1/sqrt(2).
# So |^3P_1, 1> = (1/sqrt 2) |1 1; 1 0> - (1/sqrt 2) |1 0; 1 1>

# Using <M_J=0 | H | M_J=0> for simplicity:
# |^3P_1, 0> = (1/sqrt 2) |1 -1; 1 1> - (1/sqrt 2) |1 1; 1 -1>
# <0 | H | 0> = 1/2 H(-1,1; -1,1) + 1/2 H(1,-1; 1,-1) - 1/2 H(-1,1; 1,-1) - 1/2 H(1,-1; -1,1)
# By symmetry: H(-1,1;-1,1) = H(1,-1;1,-1) and H(-1,1;1,-1) = H(1,-1;-1,1)^*

# Let me just code the simplest case: M_L = M_L' = 0, M_S = M_S' = 0.
# But M_J = M_L+M_S=0, and for ^3P_1 at M_J=0, CG(M_L=0, M_S=0 | J=1, M=0) = 0.
# So |^3P_1, 0> has no (0, 0) component. OK.

# Let me compute <^3P_0, M_J=0 | H_SS | ^3P_0, M_J=0>:
# |^3P_0, 0> = sum_{M_L+M_S=0} <1 M_L; 1 M_S | 0 0> |L=1, M_L; S=1, M_S>
# CG <1 M_L; 1 M_S | 0 0>:
#   (M_L=-1, M_S=1): 1/sqrt(3)
#   (M_L=0, M_S=0): -1/sqrt(3)
#   (M_L=1, M_S=-1): 1/sqrt(3)
# So |^3P_0, 0> = (1/sqrt 3) [|L=1 -1; S=1 1> - |L=1 0; S=1 0> + |L=1 1; S=1 -1>]

# At J=0, f_SS(0) = -2 in the Drake formula. So <^3P_0 | H | ^3P_0> / A_SS = -2
# and A_SS = -(1/2) <^3P_0 | H_SS | ^3P_0>.

# This is still complex. Let me use the EASIEST M_S = 1 case to isolate a single term.

# Actually, instead of numerical 6D integration, let me use the KNOWN VALUE of A_SS
# and the IDENTIFIED bipolar channel coefficients to ALGEBRAICALLY INVERT.

# Given:
#   A_SS = c_dir * M^K_dir + c_exch * M^K_exch + (lower-K contributions?)
#
# If we can identify the LINEAR MAP from the bipolar channels to the Drake basis,
# then Drake's c_dir = 3/50 and c_exch = -2/5 will just fall out as the projection
# coefficients.
#
# The bipolar decomposition (from Step 7) gave:
#   Bipolar(0, 2) dir = M^0_dir_I + M^2_dir_II   (asymmetric direct density)
#   Bipolar(1, 1) exch = M^1_exch_I + M^1_exch_II = M^1_exch  (symmetric exchange density)
#
# But Drake's c_dir and c_exch are FULL M^K totals. So the bipolar result is NOT in
# Drake's basis directly.
#
# UNLESS: the 9j coupling + bipolar prefactor generate multiple (k_1, k_2) channels
# that sum to give the M^K_total form. But we showed earlier that only (0, 2) direct
# and (1, 1) exchange are angular-allowed.
#
# Therefore, Drake's formula (in terms of M^K_total) does NOT emerge cleanly from
# the bipolar expansion. Drake's formula must be using a DIFFERENT decomposition.

# CONCLUSION (for Sprint 5 DV):
# The "Drake M^K integrals" are NOT just the naive r_<^K/r_>^{K+3}-kernel integrals.
# They're likely LINEAR COMBINATIONS of bipolar channel integrals, chosen so that
# the final formula is clean.
#
# Let me examine Drake 1971 (Phys Rev A 3, 908) directly for the specific M^K definition.

# NUMERICAL ANCHOR: compute A_SS from the raw 6-d integral and compare to Drake's.

# Define hydrogenic radial functions at Z=1:
def R_1s(r): return 2 * np.exp(-r)
def R_2p(r): return (1 / (2 * np.sqrt(6))) * r * np.exp(-r / 2)

# And angular Y^1_m for l=1:
def Y_10(theta, phi):
    return np.sqrt(3 / (4 * np.pi)) * np.cos(theta)
def Y_1p1(theta, phi):
    return -np.sqrt(3 / (8 * np.pi)) * np.sin(theta) * np.exp(1j * phi)
def Y_1m1(theta, phi):
    return np.sqrt(3 / (8 * np.pi)) * np.sin(theta) * np.exp(-1j * phi)

def Y_00(theta, phi):
    return np.sqrt(1 / (4 * np.pi))

# Matrix element <1s(1) 2p_1(2) | 1/r_12^3 * P_2(cos theta_12) | 1s(1) 2p_1(2)>
# = INT INT INT INT INT INT R_1s(r_1) R_2p(r_2) * Y_00^*(1) Y_1p1^*(2)
#   * 1/|r_1 - r_2|^3 * P_2(cos theta_12)
#   * R_1s(r_1) R_2p(r_2) * Y_00(1) Y_1p1(2)
#   * r_1^2 sin theta_1 dr_1 dtheta_1 dphi_1 * r_2^2 sin theta_2 dr_2 dtheta_2 dphi_2
#
# This is a 6D integral. Not easy.
#
# Cleaner approach: use the PARTIAL-WAVE EXPANSION of 1/r_12^n * P_K(cos theta_12)
# and evaluate each l-partial wave separately. Or just use the Drake paper directly.

print("\nNumerical 6D integration skipped — switching to Drake 1971 direct reference.")
