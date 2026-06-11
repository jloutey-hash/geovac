"""Sprint 4 Track DD: First-principles derivation of Drake 1971 combining
coefficients for the (1s)(2p) He 2^3P multiplet.

Target
------
Show that

  <^3P_J | H_SS  | ^3P_J> = f_SS(J)  * A_SS
  <^3P_J | H_SOO | ^3P_J> = f_SOO(J) * A_SOO

with f_SS(J) = (-2, +1, -1/5), f_SOO(J) = (+2, +1, -1) for J = 0, 1, 2 and

  A_SS  = alpha^2 * ( +3/50 * M^2_dir  - 2/5 * M^2_exch )
  A_SOO = alpha^2 * ( +3/2  * M^1_dir  - 1   * M^1_exch )

BF-D identified these coefficients by brute rational search. Our job is to
derive them symbolically via Wigner 9j algebra.

Approach
--------
We use the DIRECT Slater-determinant route with sympy-exact angular algebra.
For each J, build the |^3P_J, M_J = J> state as a coherent superposition of
two-electron spin-orbital Slater determinants. Then evaluate

  <^3P_J | V | ^3P_J>

where V is the Breit-Pauli SS or SOO tensor operator, using the Slater rules
for two-body matrix elements in the SD basis. Each SD two-body matrix element
reduces via multipole expansion to a sum of radial integrals M^k times
angular factors computed from Wigner 3j/6j/9j algebra.

Key BP operator definitions (Bethe-Salpeter §§38-39, with standard
Condon-Shortley phase and Racah tensor conventions):

  H_SS = -alpha^2 * (r_12)^{-3} * [ 3 (s_1 . r-hat_12)(s_2 . r-hat_12) - s_1 . s_2 ]

  H_SOO = alpha^2 * (r_12)^{-3} *
          { (s_1 + 2 s_2) . [r_12 x p_1]  +  (s_2 + 2 s_1) . [r_12 x p_2] }   # retarded version

The (rank-2 spatial) . (rank-2 spin) decomposition of H_SS is
  H_SS = -alpha^2 * sqrt(24 pi / 5) * sum_q (-1)^q *
            [s_1 (x) s_2]^(2)_q * Y^(2)_{-q}(r-hat_12) / r_12^3

and the corresponding rank-1 form of H_SOO has spatial tensor [r x p].

Reference constants
-------------------
The key reduced matrix elements needed:
  <1/2||s||1/2> = sqrt(3)/2   (Edmonds / Condon-Shortley)
  <l=1||L||l=1> = sqrt(6)     (rank-1 orbital angular momentum)

"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
from fractions import Fraction
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, expand, symbols, pi, I as sp_I
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

# ============================================================
# Symbolic radial integrals (no numerical values used)
# ============================================================
M0d, M1d, M2d = symbols("M0d M1d M2d", real=True)
M0e, M1e, M2e = symbols("M0e M1e M2e", real=True)
alpha_sym = symbols("alpha", positive=True, real=True)

# Map (k, channel) -> symbolic integral
M_sym = {
    (0, 'd'): M0d, (1, 'd'): M1d, (2, 'd'): M2d,
    (0, 'e'): M0e, (1, 'e'): M1e, (2, 'e'): M2e,
}


# ============================================================
# Elementary tensor matrix elements
# ============================================================

def c_tensor_me(l, m, k, q, lp, mp):
    """<l m | C^(k)_q | l' m'>, Racah-normalized tensor.

    Edmonds (5.4.1):
      <l m | C^(k)_q | l' m'> = (-1)^{l-m} sqrt((2l+1)(2l'+1)) *
                                 (l k l'; -m q m') * (l k l'; 0 0 0)
    """
    if m - q != mp:
        return Integer(0)
    ph = (-1) ** (l - m)
    return ph * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * \
           wigner_3j(l, k, lp, -m, q, mp) * \
           wigner_3j(l, k, lp, 0, 0, 0)


def s_tensor_me(m, q, mp):
    """<1/2 m | s_q | 1/2 m'> where s is the spin-1/2 rank-1 tensor.

    Use: s_{+1} = -(s_x + i s_y)/sqrt(2), s_0 = s_z, s_{-1} = (s_x - i s_y)/sqrt(2)
    On spin-1/2 states.
    """
    half = Rational(1, 2)
    if not (m in [half, -half] and mp in [half, -half]):
        return Integer(0)
    if m - q != mp:
        return Integer(0)
    if q == 0:
        return m if m == mp else Integer(0)
    elif q == 1:
        # <m=+1/2 | s_{+1} | m'=-1/2> = -1/sqrt(2)
        return -Integer(1) / sqrt(2) if (m == half and mp == -half) else Integer(0)
    elif q == -1:
        # <m=-1/2 | s_{-1} | m'=+1/2> = +1/sqrt(2)
        return Integer(1) / sqrt(2) if (m == -half and mp == half) else Integer(0)
    return Integer(0)


def L_orbital_tensor_me(l, m, q, lp, mp):
    """<l m | L^(1)_q | l' m'>, where L is the orbital angular momentum.

    For same l (= l'):
      <l m | L_0 | l m'> = m delta_{m m'}
      <l m | L_{+1} | l m'> = -sqrt((l-m+1)(l+m)/2) delta_{m', m-1}  (raises m by 1 in CS convention)

    Wait, conventions are tricky. Use spherical-tensor convention:
      L^(1)_0 = L_z
      L^(1)_{+1} = -(L_x + i L_y)/sqrt(2) = -L_+/sqrt(2)
      L^(1)_{-1} = (L_x - i L_y)/sqrt(2)  = L_-/sqrt(2)
    with L_+ |l m> = sqrt(l(l+1) - m(m+1)) |l m+1>
    and  L_- |l m> = sqrt(l(l+1) - m(m-1)) |l m-1>.
    """
    if l != lp:
        return Integer(0)
    if m - q != mp:
        return Integer(0)
    if q == 0:
        return Integer(m) if m == mp else Integer(0)
    elif q == 1:
        # L^(1)_{+1} |l m'> = -sqrt(l(l+1)-m'(m'+1))/sqrt(2) |l m'+1>
        # so <l m| L_{+1} |l m'> = -sqrt((l(l+1)-mp*(mp+1))/2) for m = mp+1
        if m == mp + 1:
            return -sqrt(Integer(l * (l + 1) - mp * (mp + 1))) / sqrt(2)
        return Integer(0)
    elif q == -1:
        # L^(1)_{-1} |l m'> = sqrt(l(l+1)-m'(m'-1))/sqrt(2) |l m'-1>
        # so <l m| L_{-1} |l m'> = sqrt((l(l+1)-mp*(mp-1))/2) for m = mp-1
        if m == mp - 1:
            return sqrt(Integer(l * (l + 1) - mp * (mp - 1))) / sqrt(2)
        return Integer(0)
    return Integer(0)


# ============================================================
# Step 1: Compute the 6-dim f_SS(J) and f_SOO(J) from the Racah 6j
# ============================================================
#
# Standard result (Edmonds 7.1.8 or Sobelman 6.1): for the diagonal matrix
# element of the scalar operator [T^(k)(space) . U^(k)(spin)]^(0) in the
# LS-coupled basis,
#
#   <(LS) J M | [T . U]^(0) | (LS) J M>
#     = (-1)^{L+S+J} * 6j{L S J; S L k} * <L||T||L> <S||U||S>
#
# NO reduced / non-reduced factor, NO extra sqrt(2J+1) for scalar operators.
# (A scalar operator has <J M | O | J M> = <J M' | O | J M'> for all M, M';
# the 6j pattern encodes the J-dependence.)
#
# We can thus predict f_SS(J), f_SOO(J) up to an overall J-independent
# normalization:


def f_J_6j(k):
    """Return { J : (-1)^{L+S+J} * 6j{L S J; S L k} } for L = S = 1, J = 0, 1, 2.

    This is (up to the J-independent factor <L||T||L><S||U||S>) the J-dependent
    coefficient of the scalar rank-(k,k) operator diagonal matrix element.
    """
    out = {}
    for J in (0, 1, 2):
        phase = (-1) ** (1 + 1 + J)   # (L + S + J), L = S = 1
        six_j = wigner_6j(1, 1, J, 1, 1, k)
        out[J] = sp.simplify(phase * six_j)
    return out


def verify_f_J_patterns():
    print("=" * 76)
    print("STEP 1: J-pattern from the rank-(k, k) scalar-operator 6j")
    print("=" * 76)

    target_SS = {0: Rational(-2), 1: Rational(1), 2: Rational(-1, 5)}
    target_SOO = {0: Rational(2), 1: Rational(1), 2: Rational(-1)}

    for k, target, label in [(2, target_SS, "SS"), (1, target_SOO, "SOO")]:
        g = f_J_6j(k)
        print(f"\n  {label} (k = {k}):")
        # Normalize so g(J=1) == target(J=1)  (the clean rational)
        norm = target[1] / g[1]
        print(f"    bare 6j with phase: {g}")
        print(f"    normalization factor (target(J=1)/g(J=1)) = {sp.simplify(norm)}")
        for J in (0, 1, 2):
            val = sp.simplify(norm * g[J])
            match = "OK" if sp.simplify(val - target[J]) == 0 else "MISMATCH"
            print(f"    J = {J}: norm * g = {val}, target = {target[J]}  [{match}]")

    return


# ============================================================
# Step 2: Build |^3P_J, M_J = J> as a sum over SDs
# ============================================================
# Canonical SD ordering: we treat each SD as a sorted pair of spin-orbitals.
# Spin-orbital key: (nlm, spin). For 1s, nl = (1, 0), m = 0. For 2p, nl = (2, 1).
# The spin-orbital tuple: ((n, l, m_l), m_s).


def spinorbital_key(n, l, m_l, m_s):
    """Canonical key for a spin-orbital."""
    return ((n, l, m_l), m_s)


def sd_key(sob1, sob2):
    """Canonical SD key: sorted tuple plus sign (+1 if sob1 < sob2, -1 else)."""
    if sob1 == sob2:
        return None, 0   # Pauli violation
    if sob1 < sob2:
        return (sob1, sob2), 1
    return (sob2, sob1), -1


def build_3P_J_state(J):
    """Build |^3P_J, M_J = J> in the SD basis.

    Returns a dict { (sob1, sob2) : coefficient }  where (sob1, sob2) is the
    canonically-ordered spin-orbital tuple.
    """
    half = Rational(1, 2)
    state = {}

    # |^3P_J, M_J = J> = sum_{M_L, M_S} <L=1 M_L S=1 M_S | J J> |L M_L; S M_S>
    for M_L in (-1, 0, 1):
        for M_S in (-1, 0, 1):
            if M_L + M_S != J:
                continue
            cg_LSJ = clebsch_gordan(1, 1, J, M_L, M_S, J)
            if sp.simplify(cg_LSJ) == 0:
                continue

            # |L=1 M_L; S=1 M_S> for (1s, 2p) config:
            # Spatial: l_a = 0, l_b = 1, coupled to L = 1, M_L
            # m_a = 0 (since l_a = 0), m_b = M_L
            cg_space = clebsch_gordan(0, 1, 1, 0, M_L, M_L)
            if sp.simplify(cg_space) == 0:
                continue

            # Spin: s_a = s_b = 1/2, coupled to S = 1, M_S
            for sa_2 in (-1, 1):
                for sb_2 in (-1, 1):
                    sa = Rational(sa_2, 2)
                    sb = Rational(sb_2, 2)
                    if sa + sb != M_S:
                        continue
                    cg_spin = clebsch_gordan(half, half, 1, sa, sb, M_S)
                    if sp.simplify(cg_spin) == 0:
                        continue

                    total = sp.simplify(cg_LSJ * cg_space * cg_spin)
                    if total == 0:
                        continue

                    # |phi_a(1) alpha(1); phi_b(2) beta(2)> with phi_a = 1s, phi_b = 2p_M_L
                    # This is a SINGLE orbital product, NOT antisymmetrized yet.
                    # Antisymmetric combination:
                    #   |ab>_AS = (1/sqrt(2)) [|a(1)b(2)> - |b(1)a(2)>]

                    sob_a = spinorbital_key(1, 0, 0, sa)
                    sob_b = spinorbital_key(2, 1, M_L, sb)

                    key, sign = sd_key(sob_a, sob_b)
                    if key is None:
                        continue   # Pauli
                    # Coefficient: total / sqrt(2) in the canonical SD basis
                    # (the 1/sqrt(2) comes from antisymmetrization normalization,
                    # and the sign from canonical ordering).
                    coef = sp.simplify(sign * total / sqrt(2))
                    state[key] = sp.simplify(state.get(key, Integer(0)) + coef)

    # Drop zeros
    state = {k: v for k, v in state.items() if sp.simplify(v) != 0}
    return state


def check_normalization(state):
    """Compute sum of |c|^2 (should be 1 for a properly normalized state)."""
    total = sum(sp.simplify(c * sp.conjugate(c)) for c in state.values())
    return sp.simplify(total)


# ============================================================
# Step 3: Two-body matrix element machinery
# ============================================================

def two_body_me_general(sd_bra, sd_ket, spatial_op, spin_op):
    """Compute <sd_bra | V | sd_ket> where V = spatial_op (x) spin_op is a
    two-body operator that factorizes into independent spatial and spin tensor
    parts, each acting on electrons 1 and 2 as separate tensor products.

    Here spatial_op and spin_op are functions of the form
        f((orbital_1_bra, orbital_2_bra), (orbital_1_ket, orbital_2_ket))  ->  number
    for the UNANTISYMMETRIZED ordering.

    sd_bra, sd_ket: canonical SD keys of the form (sob1, sob2) with sob1 < sob2.

    Slater rules for two-body: if sd_bra and sd_ket differ in 0, 1, or 2
    spin-orbitals, the matrix element is structured accordingly.
    For identical SDs:
        <ab|V|ab> - <ab|V|ba>

    Since we're computing diagonal ^3P_J matrix elements, we can compute the
    FULL <bra|V|ket> including off-diagonal SD contributions because the
    state superposition handles them.  So we compute the matrix element
    between two specific SDs.

    For simplicity, we break the SDs into their two spin-orbitals and compute
    all possible "direct" and "exchange" pairings.
    """
    (sob_a, sob_b) = sd_bra   # bra spin-orbitals (canonically ordered)
    (sob_c, sob_d) = sd_ket   # ket spin-orbitals

    # An antisymmetrized SD:
    #   |sob_a sob_b> = (1/sqrt(2)) [|sob_a(1) sob_b(2)> - |sob_b(1) sob_a(2)>]
    #
    # So <sob_a sob_b|V|sob_c sob_d> =
    #    (1/2) [ <sob_a(1) sob_b(2)|V|sob_c(1) sob_d(2)>
    #          - <sob_a(1) sob_b(2)|V|sob_d(1) sob_c(2)>
    #          - <sob_b(1) sob_a(2)|V|sob_c(1) sob_d(2)>
    #          + <sob_b(1) sob_a(2)|V|sob_d(1) sob_c(2)> ]
    #   = <sob_a(1) sob_b(2)|V|sob_c(1) sob_d(2)> - <sob_a(1) sob_b(2)|V|sob_d(1) sob_c(2)>
    # (using symmetry V(1,2) = V(2,1) for a symmetric two-body operator)
    # direct - exchange.

    direct = spatial_op(sob_a, sob_b, sob_c, sob_d) * spin_op(sob_a, sob_b, sob_c, sob_d)
    exchange = spatial_op(sob_a, sob_b, sob_d, sob_c) * spin_op(sob_a, sob_b, sob_d, sob_c)
    return direct - exchange


# ============================================================
# Step 4: Spatial and spin tensor operators
# ============================================================
#
# The SS operator:
#   H_SS = alpha^2 * sum_q (-1)^q * O^SS_q_(space)(r1,r2) * O^SS_{-q}_(spin)(s1,s2)
# where
#   O^SS_q_(space) = - sqrt(24 pi / 5) * Y^(2)_q(r-hat_12) / r_12^3
#   O^SS_q_(spin)  = [s_1 (x) s_2]^(2)_q
#
# The SOO operator:
#   H_SOO = alpha^2 * sum_q (-1)^q * O^SOO_q_(space)(r1,r2) * O^SOO_{-q}_(spin+orb)(s1,s2,p1,p2)
# which is more complex (spin and orbital mixed).
#
# For the SS rank-2 spatial operator, we use the multipole expansion:
#   -sqrt(24 pi / 5) * Y^(2)_q(r-hat_12) / r_12^3
#     = sum_{k1, k2} A(k1, k2) * [C^(k1)(1) (x) C^(k2)(2)]^(2)_q * R_{k1,k2}(r_1, r_2)
#
# where C^(k) = sqrt(4 pi / (2k+1)) Y^(k) is the Racah tensor, and R_{k1,k2}
# is a specific radial kernel.  The explicit Drake 1971 result is that for
# the (1s)(2p) ^3P state, only (k1, k2) = (0, 2), (2, 0), (1, 1) contribute.
#
# This is getting into the weeds. Let me use a DIFFERENT, simpler approach.


# ============================================================
# Simpler approach: use Drake's "N_k" form directly
# ============================================================
#
# Drake 1971 Eqs. (17), Bethe-Salpeter 39.13 give the closed form:
#
#   <^3P_J | H_SS | ^3P_J> = f_SS(J) * alpha^2 * N^(2)
#
# where
#   N^(2) = -(3/10) * M^2_exch  +  (9/100) * M^2_dir ???
#
# The specific combination depends on convention.  But the KEY POINT is that
# the J-dependent part f_SS(J) equals the 6j pattern (-2, +1, -1/5) UP TO an
# overall sign and normalization.
#
# In BF-D's convention, BF-D defines f_SS via
#   E_SS(J) = A_SS * f_SS(J)    with f_SS(0,1,2) = (-2, +1, -1/5)
# This means A_SS = <M_J = 1 | H_SS | M_J = 1> since f_SS(J=1) = 1.
#
# So A_SS is literally the J=1 matrix element of H_SS.  And we just need to
# compute <^3P_1, M_J=1 | H_SS | ^3P_1, M_J=1> symbolically and extract
# the coefficients of M^2_dir, M^2_exch.

# This is the clean plan. Let me implement it by explicitly building the
# two-electron SD state for J=1, M_J=1, and evaluating the H_SS matrix element
# via tensor decomposition.

# For the DIAGONAL matrix element of a rank-(k_spatial, k_spin) tensor operator
# coupled to scalar, the key identity is Racah (7.1.7):
#
# <(l_a s_a)(l_b s_b) L S J | [T^(k)(1) (x) U^(k)(2)]^(0) | (l_a s_a)(l_b s_b) L S J>
#   = (-1)^{L+S+J} * 6j{L L k; S S k; J J 0} * ...
#
# No wait, that's for a scalar operator built from (rank-k space) tensored with
# (rank-k spin).  For our case we have:
#
#   H_SS = T^(2)(space, electrons 1&2) . U^(2)(spin, electrons 1&2)
#
# where T^(2) is a SPACE rank-2 tensor FORMED FROM individual electron 1 and
# electron 2 spatial tensors, and U^(2) is analogous for spin.
#
# The DIAGONAL matrix element in the LS-coupled (l_1 l_2) L, (s_1 s_2) S basis:
#
# <(l_a l_b) L (s_a s_b) S J | T^(2)(space) . U^(2)(spin) | (l_a l_b) L (s_a s_b) S J>
#   = (-1)^{L+S+J} * 6j{L S J; S L 2} *
#        <(l_a l_b) L || T^(2)(space) || (l_a l_b) L>
#        <(s_a s_b) S || U^(2)(spin) || (s_a s_b) S>
#
# using Edmonds (7.1.1) -- reducing the L-S-J coupling.

# The spatial reduced matrix element is itself a 9j:
# <(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(2) || (l_a l_b) L>
#   = sqrt((2L+1)(2*2+1)(2L+1)) * 9j{l_a l_b L; k1 k2 2; l_a l_b L}
#     * <l_a || C^(k1) || l_a> * <l_b || C^(k2) || l_b>
#
# For (l_a, l_b) = (0, 1):
# - <0 || C^(k1) || 0> = delta_{k1, 0}
# - <1 || C^(k2) || 1> = sqrt(3) * delta_{k2,0}  OR  -sqrt(30)/5 * delta_{k2,2}
#
# So the "direct" spatial reduced m.e. picks up:
#   (k1, k2) = (0, 0)  -- rank 0, not 2; excluded
#   (k1, k2) = (0, 2)  -- rank 0+2=2, OK
# Therefore only (k1, k2) = (0, 2) contributes to the spatial "direct" channel.
#
# The spatial "exchange" reduced m.e.:
# <(0 1) L=1 || [C^(k1)(1) (x) C^(k2)(2)]^(2) || (1 0) L'=1>
#   = sqrt(9 * 5) * 9j{0 1 1; k1 k2 2; 1 0 1}
#     * <0||C^(k1)||1> <1||C^(k2)||0>
# With <0||C^(k1)||1> = -delta_{k1,1} (from reduced_C(0,1,1) = -1)
#      <1||C^(k2)||0> = +delta_{k2,1}
# So only (k1, k2) = (1, 1) contributes to exchange.
#
# PHYSICAL INTERPRETATION (this is THE key structural observation):
# * DIRECT SS channel: spatial rank-(0, 2) -- the "retarded Slater integral"
#   M^2_dir is the radial piece with k1=0 multipole on electron 1 and k2=2
#   multipole on electron 2 (effectively: orb density rho_{1s}(r_1) times
#   rank-2 spherical harmonic of orb density rho_{2p}(r_2)).
# * EXCHANGE SS channel: spatial rank-(1, 1) -- this is a rank-1 coupling
#   channel on each side with 1s-2p mixing.  BUT the 9j with rank 2 ensures
#   the TOTAL rank is 2.
#
# For SOO (rank 1):
# * DIRECT channel: (k1, k2) triangular to 1, with both same-orbital.
#   (k1, k2) = (0, 1): <0||C^0||0><1||C^1||1> triangular to 1 -- ALLOWED
#   (k1, k2) = (1, 0): <0||C^1||0><1||C^0||1> -- <0||C^1||0> = 0 because (l k l' = 0 1 0; 0 0 0) = 0
#   So only (k1, k2) = (0, 1) contributes to "direct" in SOO.
# * EXCHANGE channel: (k1, k2) triangular to 1 with l<->l swap.
#   (k1, k2) = (1, 0): <0||C^1||1><1||C^0||0> -- ALLOWED
#   (k1, k2) = (0, 1): <0||C^0||1><1||C^1||0> -- <0||C^0||1> = 0.
#   So only (k1, k2) = (1, 0) contributes to "exchange" in SOO.
#
# Now the spin operators:
#
# SS operator has spin part [s_1 x s_2]^(2)_q.  The reduced matrix element for
# S=1 coupled is <S=1 || [s_1 (x) s_2]^(2) || S=1>.
# SOO operator has spin part (s_1 + 2 s_2) which is rank-1 but a SUM of two
# rank-1 operators.  The reduced m.e. for S=1 -> S=1 of s_1 or s_2 alone is
# given by Edmonds 7.1.8:
#   <S || s_i || S> (for S = coupled s_1 + s_2 = 1)
#    = (-1)^{1/2 + 1/2 + S + 1} sqrt((2S+1)(2*1+1)) * 6j{1/2 1/2 S; 1 S 1/2}
#      * <1/2||s||1/2>
# This gives non-zero for i=1 AND i=2 (both contribute).

# Let's now evaluate these symbolically.


def spin_reduced_rankk_sigma(k):
    """<S=1 || [sigma_1 (x) sigma_2]^(k) || S=1> using Edmonds 7.1.7.

    Note: sigma = 2 s, so <1/2 || sigma || 1/2> = 2 <1/2 || s || 1/2> = 2 sqrt(3)/2 = sqrt(3).
    Wait -- <1/2 || s || 1/2> = sqrt(3)/2 in some conventions, sqrt(3/2) in others.

    Edmonds (uses S^2 = S(S+1)): <1/2 || s || 1/2> = sqrt(3/2)  (Edmonds 5.4)
    Wigner convention (uses J^2 = j(j+1)): same.

    So <1/2||s||1/2> = sqrt(3/2).

    Edmonds 7.1.7:
      <(j1 j2) J || [A^(k1)(1) (x) B^(k2)(2)]^(K) || (j'_1 j'_2) J'>
      = sqrt((2K+1)(2J+1)(2J'+1)) *
        9j{j1  j2  J; k1 k2 K; j'_1  j'_2  J'} * <j1||A||j'_1> <j2||B||j'_2>

    For j1 = j2 = j'_1 = j'_2 = 1/2, J = J' = 1, k1 = k2 = 1, K = k:
      = sqrt((2k+1) * 3 * 3) * 9j{1/2 1/2 1; 1 1 k; 1/2 1/2 1} * (sqrt(3/2))^2
      = 3 sqrt(2k+1) * 9j{...} * (3/2)
    """
    half = Rational(1, 2)
    red_s = sqrt(Rational(3, 2))  # <1/2||s||1/2> in Edmonds convention
    out = sqrt(Integer(3 * 3 * (2 * k + 1))) * \
          wigner_9j(half, half, 1, 1, 1, k, half, half, 1) * red_s * red_s
    return sp.simplify(out)


def spin_s1_reduced_S1():
    """<S=1 || s_1 || S=1> via Edmonds 7.1.8.

    Formula: <(j1 j2) J || O_1 || (j1 j2) J'> = delta_{J, J'} * (-1)^{j1 + j2 + J + 1}
      * sqrt((2J+1)(2J+1)) * 6j{J j2 j1; 1 j1 J} * <j1||O||j1>

    Wait, the general 7.1.8 for acting on electron 1 only:
      <(j1 j2) J || A_1^(k) || (j'_1 j_2) J'>
      = (-1)^{j1 + j2 + J' + k} sqrt((2J+1)(2J'+1)) 6j{j1 J j2; J' j'_1 k} * <j1||A||j'_1>
    """
    half = Rational(1, 2)
    S = 1
    j1, j2 = half, half
    Sp = 1   # J' = J
    k = 1
    red_s = sqrt(Rational(3, 2))
    out = (-1) ** (j1 + j2 + Sp + k) * sqrt(Integer((2 * S + 1) * (2 * Sp + 1))) * \
          wigner_6j(j1, S, j2, Sp, j1, k) * red_s
    return sp.simplify(out)


def spin_s2_reduced_S1():
    """<S=1 || s_2 || S=1> via Edmonds 7.1.8.  (By symmetry equals s_1 for S=1.)"""
    half = Rational(1, 2)
    S = 1
    j1, j2 = half, half
    Sp = 1
    k = 1
    red_s = sqrt(Rational(3, 2))
    # Acting on electron 2:
    out = (-1) ** (j1 + j2 + S + k) * sqrt(Integer((2 * S + 1) * (2 * Sp + 1))) * \
          wigner_6j(j2, S, j1, Sp, j2, k) * red_s
    return sp.simplify(out)


def show_spin_reduced_matrix_elements():
    print("\n" + "=" * 76)
    print("STEP 2: Spin reduced matrix elements for S = 1")
    print("=" * 76)
    print(f"  <1||s||1/2>  [single electron]   = {sqrt(Rational(3, 2))}")
    print(f"  <S=1||s_1||S=1>                 = {spin_s1_reduced_S1()}")
    print(f"  <S=1||s_2||S=1>                 = {spin_s2_reduced_S1()}")
    for k in (1, 2):
        val = spin_reduced_rankk_sigma(k)
        # Convert sigma-x-sigma to s-x-s: sigma = 2 s so [sigma x sigma]^k = 4 [s x s]^k
        val_s = val / Integer(4)
        print(f"  <S=1 || [sigma_1 (x) sigma_2]^({k}) || S=1> = {val}")
        print(f"  <S=1 || [s_1 (x) s_2]^({k}) || S=1>         = {val_s}")


# ============================================================
# Step 5: Full <^3P_J | H_SS | ^3P_J> via Edmonds reduction
# ============================================================
#
# The Racah reduction for the scalar operator [T^(k)(space) . U^(k)(spin)]^(0)
# in the (LS)J coupled basis is (Edmonds 7.1.8, with k_1 = k_2 = k and K = 0):
#
# <(LS) J || [T^(k) . U^(k)]^(0) || (LS) J> / sqrt(2J+1)
#   = (-1)^{L+S+J} * 6j{L S J; S L k} * <L||T||L> <S||U||S>
#
# The physical diagonal matrix element is:
# <(LS) J M | [T . U]^(0) | (LS) J M> = <(LS) J || [T . U]^(0) || (LS) J> / sqrt(2J+1)
#   = (-1)^{L+S+J} * 6j{L S J; S L k} * <L||T||L> <S||U||S> / sqrt(2J+1)
#
# Wait, that contradicts my earlier formula by a sqrt(2J+1).  Let me check.
#
# For a rank-0 tensor O^(0), Wigner-Eckart: <J M | O^(0) | J M'> = delta_{M M'} * <J||O||J> / sqrt(2J+1).
# So the diagonal (M = M') matrix element is <J||O||J> / sqrt(2J+1).
#
# The reduced m.e. of the coupled [T^(k) . U^(k)]^(0) is (from 7.1.7 with K = 0):
#   <(LS)J || [T^(k) . U^(k)]^(0) || (LS)J>
#   = sqrt(2J+1) * (-1)^{L+S+J} * 6j{L L k; S S k; J J 0}  -- actually no, the 9j{...;J J 0} reduces
#
# The 9j with a zero argument:
#   {a b c; d e f; g h 0} = delta_{g,h} delta_{c,f} * (-1)^{b+c+d+g} / sqrt((2c+1)(2g+1)) * 6j{a b c; e d g}
#
# Applied: {L L k; S S k; J J 0} = delta_{J,J} delta_{k,k} (-1)^{L+k+S+J} / sqrt((2k+1)(2J+1)) * 6j{L L k; S S J}
#  = (-1)^{L+S+J+k} / sqrt((2k+1)(2J+1)) * 6j{L L k; S S J}
#
# And 6j{L L k; S S J} = 6j{L S J; S L k} (symmetry of 6j).
#
# So <(LS)J || [T . U]^(0) || (LS)J>
#   = sqrt(2J+1) * (-1)^{L+S+J} * (-1)^{L+S+J+k} / sqrt((2k+1)(2J+1)) * 6j{L S J; S L k}
#     * <L||T||L> <S||U||S>
#   = (-1)^{2(L+S+J)+k} / sqrt(2k+1) * 6j{L S J; S L k} * <L||T||L> <S||U||S>
#   = (-1)^k / sqrt(2k+1) * 6j{L S J; S L k} * <L||T||L> <S||U||S>
#
# And the physical diagonal m.e. is:
# <JM|O|JM> = <J||O||J> / sqrt(2J+1)
#           = (-1)^k / sqrt((2k+1)(2J+1)) * 6j{L S J; S L k} * <L||T||L> <S||U||S>
#
# So the PHYSICAL J-dependence is
#   h(J, k) = (-1)^k / sqrt((2J+1)(2k+1)) * 6j{L S J; S L k}
#
# Let me verify this is the "f_SS(J)" pattern:

def physical_J_pattern(k):
    """Physical J-dependence h(J, k) = (-1)^k / sqrt((2J+1)(2k+1)) * 6j{L S J; S L k}"""
    out = {}
    L, S = 1, 1
    for J in (0, 1, 2):
        six_j = wigner_6j(L, S, J, S, L, k)
        h = (-1) ** k / sqrt(Integer((2 * J + 1) * (2 * k + 1))) * six_j
        out[J] = sp.simplify(h)
    return out


# ============================================================
# Step 6: Driver
# ============================================================


def main():
    print("=" * 76)
    print("Sprint 4 Track DD: Drake 1971 coefficient derivation")
    print("=" * 76)

    # Step 1: J-pattern from 6j
    verify_f_J_patterns()

    # Also print the physical h(J, k):
    print("\n" + "=" * 76)
    print("Physical J-pattern h(J, k) = (-1)^k / sqrt((2J+1)(2k+1)) * 6j{L S J; S L k}")
    print("=" * 76)
    for k in (1, 2):
        h = physical_J_pattern(k)
        print(f"\n  k = {k}:")
        for J in (0, 1, 2):
            print(f"    J = {J}: h = {h[J]}")
        # Ratios: h(J)/h(1)
        print(f"    ratios h(J)/h(1):")
        for J in (0, 1, 2):
            r = sp.simplify(h[J] / h[1])
            print(f"      J = {J}: {r}")

    # Step 2: Spin reduced m.e.
    show_spin_reduced_matrix_elements()

    # Step 3: |^3P_J> normalization
    print("\n" + "=" * 76)
    print("STEP 3: |^3P_J, M_J=J> state normalization check")
    print("=" * 76)
    for J in (0, 1, 2):
        state = build_3P_J_state(J)
        n = check_normalization(state)
        print(f"  |^3P_{J}, M_J={J}>: ||.||^2 = {n}  ({len(state)} SDs)")


if __name__ == "__main__":
    main()
