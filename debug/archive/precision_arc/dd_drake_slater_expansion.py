"""Sprint 4 Track DD: Slater-determinant / brute-force derivation of Drake 1971
combining coefficients A_SS, A_SOO for He 2^3P.

Approach
--------
We compute the diagonal matrix element <^3P_J, M_J=J | H_BP | ^3P_J, M_J=J>
of the Breit-Pauli operators H_SS (rank-2) and H_SOO (rank-1) by:

  1. Constructing the ^3P_J state as a linear combination of
     |L=1, M_L; S=1, M_S> states via Clebsch-Gordan.
  2. Decomposing each |L M_L; S M_S> into Slater determinants in the
     (1s, 2p_0, 2p_+1, 2p_-1) * (up, dn) basis.
  3. Applying the two-body operator matrix element formulas for SS and SOO.
  4. Collecting the coefficients of M^k_dir and M^k_exch.

The radial integrals are used ONLY symbolically (M^k_dir, M^k_exch as free
symbols). The angular integration is done exactly in sympy rationals.

The SS and SOO operators are treated in their standard rank-2 and rank-1
tensor forms respectively. The operator structure used:

  H_SS = -alpha^2 * (sqrt(24 pi / 5)) * sum_q (-1)^q
              [sigma_1 (x) sigma_2]^(2)_q * Y^(2)_{-q}(r-hat_12) / r_12^3

  (Bethe-Salpeter §38.14, or equivalently
  H_SS = alpha^2 * [ (s_1 . s_2 - 3 (s_1 . r-hat_12)(s_2 . r-hat_12))/r_12^3 ]
  up to an overall factor, with (s_1 . s_2) having no J-dependence beyond
  diagonal on S=1 eigenvalue and the rank-2 part carrying all J-dependence.)

  H_SOO = alpha^2 * sum_{i != j} (1/r_ij^3) [r_ij x p_i] . (s_i + 2 s_j)

  (Bethe-Salpeter §38.15, symmetrized.)

Convention check
----------------
We take the Drake 1971 / Johnson 2007 convention and identify the "direct"
channel with bra and ket Slater-determinant configurations both of the form
|1s(1), 2p(2); m_l_a, m_l_b, spin_a, spin_b>, and "exchange" as the
configuration swap |1s(2), 2p(1); ...>.

Radial conventions: the Slater integrals M^k follow the `geovac.breit_integrals`
convention
    M^k_dir  = R^k_BP(1s,1s; 2p,2p)   (direct: electron 1 holds 1s-1s, electron 2 holds 2p-2p)
    M^k_exch = R^k_BP(1s,2p; 2p,1s)   (exchange: each electron holds a 1s-2p bra-ket)

References
----------
- Bethe-Salpeter §§38-39 (spin-spin and spin-other-orbit tensor forms)
- Drake 1971 (explicit matrix elements for He excited states)
- Edmonds 1957 (tensor algebra conventions)
- Cowan 1981 §12 (SS in LS-coupled basis)
- sympy.physics.wigner (3j, 6j, 9j via Racah formula)
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, expand, symbols, Matrix, I as sp_I, pi, S as sp_S
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan, gaunt


# ============================================================
# Symbolic radial integrals
# ============================================================
M0d, M1d, M2d = symbols("M0d M1d M2d", real=True)
M0e, M1e, M2e = symbols("M0e M1e M2e", real=True)
alpha = symbols("alpha", positive=True, real=True)

# Dictionary: {(k, channel): symbol}  channel = 'd' or 'e'
M = {
    (0, 'd'): M0d, (1, 'd'): M1d, (2, 'd'): M2d,
    (0, 'e'): M0e, (1, 'e'): M1e, (2, 'e'): M2e,
}


# ============================================================
# LS -> (m_L, M_S) decomposition for ^3P_J (L=S=1)
# ============================================================

def expand_LSJ_state(L, S, J, M_J):
    """|L S J M_J> = sum_{M_L, M_S} <L M_L S M_S | J M_J> |L M_L; S M_S>.

    Returns dict {(M_L, M_S): coefficient}.
    """
    out = {}
    for M_L in range(-L, L + 1):
        for M_S_sympy in [Rational(-S), Rational(0), Rational(S)] if S == 1 else [
                Rational(-S), Rational(S)
        ]:
            # integer S = 1 here
            pass

    # Simpler: iterate m_L over integers and M_S over integers for S=1
    for M_L in range(-int(L), int(L) + 1):
        for M_S in range(-int(S), int(S) + 1):
            if M_L + M_S != M_J:
                continue
            cg = clebsch_gordan(L, S, J, M_L, M_S, M_J)
            if sp.simplify(cg) != 0:
                out[(M_L, M_S)] = cg
    return out


def ml_expand(l_a, l_b, L, M_L):
    """|l_a l_b L M_L> = sum_{m_a, m_b} <l_a m_a l_b m_b | L M_L> |l_a m_a; l_b m_b>.

    Returns dict {(m_a, m_b): coefficient}.
    """
    out = {}
    for m_a in range(-int(l_a), int(l_a) + 1):
        for m_b in range(-int(l_b), int(l_b) + 1):
            if m_a + m_b != M_L:
                continue
            cg = clebsch_gordan(l_a, l_b, L, m_a, m_b, M_L)
            if sp.simplify(cg) != 0:
                out[(m_a, m_b)] = cg
    return out


def ms_expand(s_a, s_b, S, M_S):
    """Similar for spins s_a = s_b = 1/2 coupled to S = 0 or 1."""
    out = {}
    half = Rational(1, 2)
    for sa_2 in (-1, 1):
        for sb_2 in (-1, 1):
            sa = Rational(sa_2, 2)
            sb = Rational(sb_2, 2)
            if sa + sb != M_S:
                continue
            cg = clebsch_gordan(half, half, S, sa, sb, M_S)
            if sp.simplify(cg) != 0:
                out[(sa, sb)] = cg
    return out


# ============================================================
# Rank-2 SS operator: matrix element in (m_a, m_b, sa, sb) basis
# ============================================================
#
# H_SS = C_SS * sum_q (-1)^q [sigma_1 (x) sigma_2]^(2)_q * Y^2_{-q}(r-hat_12) / r_12^3
#
# The operator factorizes as spin tensor x spatial tensor:
#   <sa_f sb_f|[sigma_1 x sigma_2]^(2)_q|sa_i sb_i> (rank 2 spin tensor)
#   <spatial_f|Y^2_{-q}(r-hat_12)/r_12^3|spatial_i> (rank 2 spatial tensor)
#
# For the spatial part on (1s, 2p) orbitals:
#   The multipole expansion of Y^2(r-hat_12)/r_12^3 = ...
#   (after expansion in terms of individual orbital tensors)
#   selects specific (k_1, k_2) multipole contributions.
#
# We use the Racah tensor form: the operator on electrons (1, 2) of the
# coupled tensor [T^(k1)(1) (x) T^(k2)(2)]^(K) has matrix elements
#   <f1 f2|[T^(k1)(1) x T^(k2)(2)]^(K)_Q|i1 i2>
#   = sum_{q1, q2} <k1 q1 k2 q2|K Q> <f1|T^(k1)_q1|i1> <f2|T^(k2)_q2|i2>
#
# For the Breit-Pauli SS operator, we use Bethe-Salpeter eq 38.14:
#   H_SS = -alpha^2 sqrt(24 pi / 5) * sum_q (-1)^q
#           [s_1 x s_2]^(2)_q * Y^(2)_{-q}(r-hat_12) / r_12^3
#
# where I use s_i (not sigma_i) to match Edmonds and BS.
# The factor in front matters for the overall A_SS magnitude.
#
# NOTE: Y^2(r-hat_12) / r_12^3 does NOT factorize into tensor operators on
# separate electrons directly.  It factorizes after the multipole expansion
# of 1/r_12^3 weighted by Y^2(r-hat_12).  This is the key complication.
#
# The standard result is (Edmonds 6.4):
#   Y^2_{-q}(r-hat_12) / r_12^3 = sum_{k1, k2} (coupling coefficient) *
#       [Y^(k1)(1-hat) (x) Y^(k2)(2-hat)]^(2)_{-q} * f(r_<,r_>)
#
# where f(r_<,r_>) is a rank-specific radial kernel and (k1, k2) are in a
# triangular relation with rank 2.  This allows the operator to be written
# as a direct sum of individual-electron tensor operators.
#
# For the rank-2 SS case specifically, the clean form is (Jucys-Levinson,
# Biedenharn 1981 §3; Judd 1963):
#   H_SS = -(5 * sqrt(6 pi / 5)) * alpha^2 * sum_q (-1)^q [s_1 x s_2]^(2)_q *
#              sum_{k} (-1)^k C(k, 2) * [Y^(k)(1) x Y^(k+1)(2)]^(2)_{-q} * f_k(r_1, r_2)
#
# but this is getting hairy.  We take a different route below.


# ============================================================
# Direct computational route: use Bethe-Salpeter §39 closed formulas
# ============================================================
#
# The BS §39 closed formulas for <^3P_J | H_SS | ^3P_J> and <^3P_J | H_SOO | ^3P_J>
# in the (1s)(np) configuration are (schematically):
#
#   <(1s)(np)^3P_J | H_SS | ^3P_J> = C_SS(J) * Q(1s, np)
#   <(1s)(np)^3P_J | H_SOO| ^3P_J> = C_SOO(J) * P(1s, np)
#
# where Q and P are specific radial Slater integrals and C_SS, C_SOO are
# J-dependent angular coefficients.
#
# For Bethe-Salpeter §39.14 (spin-spin):
#   C_SS(J) = ... // (J-dependent rank-2 factor)
#   Q = (3/50) alpha^2 * (M^2_direct - ... M^2_exchange)
#
# For §39.16 (spin-other-orbit):
#   C_SOO(J) = ... // (J-dependent rank-1 factor)
#   P = (1/2) alpha^2 * (M^1_direct - M^1_exchange)
#
# These formulas are stated in many references BUT the factor (3/50) and (1/2)
# come from specific 9j / Racah evaluations.  To derive them from first
# principles:
#
# STEP 1: Verify that the J-pattern f_SS(J) = (-2, 1, -1/5) IS the result
#         of the rank-2 angular coupling in (L=S=1) -> J.
#
# STEP 2: Verify that the J-pattern f_SOO(J) = (2, 1, -1) IS the result
#         of the rank-1 angular coupling.
#
# STEP 3: Compute the spatial + spin reduced matrix elements for rank-2 and
#         rank-1 operators in the (1s)(2p) spatial configuration.
#
# STEP 4: Combine to give A_SS, A_SOO.
#

def f_pattern_from_6j(L, S, k):
    """Return the f(J) pattern from the rank-k scalar tensor operator.

    The diagonal matrix element of a (k,k)-coupled-to-scalar operator
    in |(L S) J M> is, up to reduced factors,
      <(LS)J M | [T^(k)(space) . U^(k)(spin)]^(0) | (LS)J M>
      = (-1)^{L+S+J} {L L k; S S k; J J 0}(9j) <L||T||L> <S||U||S>
      = (-1)^{k} / sqrt((2J+1)(2k+1)) * 6j{L S J; S L k} * <L||T||L> <S||U||S>

    The J-dependence is entirely in the factor
      g(J, k) = (-1)^{k} * 6j{L S J; S L k} / sqrt(2J+1)

    (the sqrt(2k+1) factors into the overall normalization).
    """
    out = {}
    for J in range(int(abs(L - S)), int(L + S) + 1):
        six_j = wigner_6j(L, S, J, S, L, k)
        g = (-1) ** k / sqrt(2 * J + 1) * six_j
        out[J] = sp.simplify(g)
    return out


# ============================================================
# Step 1: J-pattern, rescaled to match (f_SS, f_SOO) normalizations
# ============================================================

def verify_J_pattern():
    print("=" * 76)
    print("STEP 1: J-pattern verification (f_SS / f_SOO)")
    print("=" * 76)
    # f_SS (k=2)
    g_SS = f_pattern_from_6j(1, 1, 2)
    # f_SOO (k=1)
    g_SOO = f_pattern_from_6j(1, 1, 1)

    # Expected targets:
    f_SS_target = {0: Rational(-2), 1: Rational(1), 2: Rational(-1, 5)}
    f_SOO_target = {0: Rational(2), 1: Rational(1), 2: Rational(-1)}

    print("\n  SS (k=2):")
    for J in (0, 1, 2):
        print(f"    J={J}: g = {g_SS[J]}, target f_SS(J) = {f_SS_target[J]}")
    # Normalization factor: divide target by g at J=0:
    norm_SS = f_SS_target[0] / g_SS[0]
    print(f"    Normalization (target J=0 / g(J=0)) = {sp.simplify(norm_SS)}")
    # Check that this normalization gives the right ratios for J=1, J=2
    for J in (0, 1, 2):
        norm_val = sp.simplify(norm_SS * g_SS[J])
        match = "OK" if sp.simplify(norm_val - f_SS_target[J]) == 0 else "MISMATCH"
        print(f"    J={J}: norm * g = {norm_val}, target = {f_SS_target[J]}  [{match}]")

    print("\n  SOO (k=1):")
    for J in (0, 1, 2):
        print(f"    J={J}: g = {g_SOO[J]}, target f_SOO(J) = {f_SOO_target[J]}")
    norm_SOO = f_SOO_target[0] / g_SOO[0]
    print(f"    Normalization (target J=0 / g(J=0)) = {sp.simplify(norm_SOO)}")
    for J in (0, 1, 2):
        norm_val = sp.simplify(norm_SOO * g_SOO[J])
        match = "OK" if sp.simplify(norm_val - f_SOO_target[J]) == 0 else "MISMATCH"
        print(f"    J={J}: norm * g = {norm_val}, target = {f_SOO_target[J]}  [{match}]")

    return norm_SS, norm_SOO


# ============================================================
# Step 2: Spatial reduced matrix element of the rank-k operator on (1s)(2p)
# ============================================================
#
# The key formula, for the DIAGONAL path (bra = ket = (l_a, l_b) = (0, 1), L=1)
# of a rank-k operator built from <l_a' | C^(k1) | l_a> <l_b' | C^(k2) | l_b>
# factorable tensor:
#
#   <(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(k) || (l_a' l_b') L'>
#   = sqrt((2L+1)(2k+1)(2L'+1)) *
#         9j{ l_a  l_b  L ; k1  k2  k ; l_a' l_b'  L' } *
#         <l_a || C^(k1) || l_a'> <l_b || C^(k2) || l_b'>
#
# For the SS/SOO radial-integral convention, we use k1 = k (with k2 = 0) on
# the DIRECT path (both 1s's stay on e1) and k1 = 1, k2 = k-1 cross on the
# EXCHANGE path.  WAIT -- this requires careful bookkeeping:
#
# For SS (rank-2):
#   Multipole expansion of Y^2(r-hat_12)/r_12^3 in terms of |(r_1, theta_1)|
#   and |(r_2, theta_2)| separately gives kernel contributions of the form
#       f_kk'(r_1, r_2) * [Y^(k)(1-hat) x Y^(k')(2-hat)]^(2)_{-q}
#   with (k, k') triangular to 2, i.e., (0,2), (1,1), (2,0), (2,2), (1,3), ...
#
# For SOO (rank-1):
#   Similar with (k, k') triangular to 1: (0, 1), (1, 0), (1, 1), ...
#
# CRITICAL INSIGHT: the M^k_dir and M^k_exch in the `breit_integrals` module
# are defined at a SPECIFIC k = multipole of the 1/r_12^{k+3} radial factor.
# They are NOT the (k, k') coupling indices of the spatial tensor operators!
#
# Drake 1971 uses the notation:
#   M^k = integral of r_<^k / r_>^(k+3) weighted by the orbital density product.
#
# This is the radial-only integral.  The angular part on each electron is
# already projected out.  So the 9j expansion above needs to be summed over
# the two independent (k1, k2) multipoles.
#
# Let me think again...


# ============================================================
# Actually, the cleanest approach is different.  We compute the <^3P_J|H|J>
# matrix element in the coupled Slater-determinant formalism.
#
# For two electrons in orbitals (phi_a, phi_b) = (1s, 2p), with L, S, J coupled,
# the antisymmetrized state is (schematically)
#
#   |^3P_J M_J> = (1/sqrt(2)) [ |phi_a(1) phi_b(2); spins singlet/triplet, LSJ>
#                             - |phi_a(2) phi_b(1); spins singlet/triplet, LSJ> ]
#
# The diagonal matrix element of any symmetric two-body operator V(1,2) is
#
#   <^3P_J|V|^3P_J> = <ab|V|ab> - <ab|V|ba>    (for triplet spin-antisymmetric
#                                                spatial combination
#                                                -- wait this is wrong)
#
# For spatial-antisymmetric + spin-symmetric (triplet):
#   <^3P_J|V|^3P_J> (spatial) = (1/2) [<ab|V|ab> + <ba|V|ba>] - (1/2) [<ab|V|ba> + <ba|V|ab>]
#                              = <ab|V|ab> - <ab|V|ba>
# (direct minus exchange)
#
# For spatial-symmetric + spin-antisymmetric (singlet): <ab|V|ab> + <ab|V|ba>
# (direct plus exchange)
#
# This direct-minus-exchange pattern applies UNIFORMLY to all spin-independent
# pieces of V.  But for SPIN-DEPENDENT V, the spin matrix element enters as a
# multiplicative factor and also affects the direct/exchange recoupling.
#
# For SS and SOO, V is spin-dependent: V = T(space) . U(spin).  The direct and
# exchange matrix elements of T(space) get multiplied by different spin matrix
# elements of U(spin).
#
# The GENERAL formula for the two-body matrix element with spin-dependent
# rank-(k,k) tensor operator is (Slater rules generalized for spin-tensor):
#
#   <phi_a alpha_1 phi_b alpha_2|V|phi_a alpha_1 phi_b alpha_2>
#   = direct_spatial_me * spin_direct_me
#   - exchange_spatial_me * spin_exchange_me
#
# where spin_direct_me = <alpha_1|U_1|alpha_1> * <alpha_2|U_2|alpha_2>
# and   spin_exchange_me = <alpha_1|U_1|alpha_2> * <alpha_2|U_2|alpha_1>
# (swapped spin slots).
#
# But we need to be CAREFUL: V = [T(space)(1,2) x U(spin)(1,2)]^(0) couples
# space and spin AS A WHOLE between the two electrons, so the direct and
# exchange are only well-defined after the state is written in a
# product-of-coordinates (not LSJ-coupled) basis.
#
# Let's just DO it:


# ============================================================
# Full Slater-determinant computation
# ============================================================

def build_3P_J_state(J):
    """Build the |^3P_J, M_J = J> state as a dictionary over SD amplitudes.

    SD basis: for 2 electrons, we use |orbital_a, spin_a; orbital_b, spin_b>_AS
    antisymmetrized under exchange.  The orbital index runs over (1s, 2p_{-1},
    2p_0, 2p_{+1}).

    Returns
    -------
    state : dict {(orb1, s1, orb2, s2): coefficient}
        orb_i in {('1s',), ('2p', m_l)} where m_l in {-1, 0, 1}, and s_i in
        {-1/2, +1/2}.  The dict is antisymmetric under (orb1,s1) <-> (orb2,s2)
        but we only store the "ordered" SD (e.g., by some canonical ordering)
        with the normalization sqrt(2) for non-equivalent pairs.

    In practice, we represent each basis state as a FROZEN SET
    {(orb, s), (orb', s')}.
    """
    state = {}
    L, S, M_J = 1, 1, J

    # Expand |(L=1, S=1) J M_J> in |LS M_L M_S>
    lsj_dict = {}
    for M_L in range(-L, L + 1):
        for M_S in range(-S, S + 1):
            if M_L + M_S != M_J:
                continue
            cg = clebsch_gordan(L, S, J, M_L, M_S, M_J)
            if sp.simplify(cg) != 0:
                lsj_dict[(M_L, M_S)] = cg

    # For each (M_L, M_S), expand into (m_a, m_b; sa, sb)
    # Note: for (1s)(2p), the "a" orbital is always 1s (l_a=0, m_a=0) and
    # "b" orbital is always 2p (l_b=1, m_b varying).
    for (M_L, M_S), cg_lsj in lsj_dict.items():
        # Spatial CG for coupling l_a=0 x l_b=1 -> L=1, M_L
        # m_a must be 0 (since l_a = 0), m_b = M_L
        if M_L < -1 or M_L > 1:
            continue
        m_a, m_b = 0, M_L
        cg_space = clebsch_gordan(0, 1, L, m_a, m_b, M_L)

        # Spin CG for coupling s_a = s_b = 1/2 -> S = 1, M_S
        half = Rational(1, 2)
        for sa_2 in (-1, 1):
            for sb_2 in (-1, 1):
                sa = Rational(sa_2, 2)
                sb = Rational(sb_2, 2)
                if sa + sb != M_S:
                    continue
                cg_spin = clebsch_gordan(half, half, S, sa, sb, M_S)
                total_cg = sp.simplify(cg_lsj * cg_space * cg_spin)
                if total_cg == 0:
                    continue
                # Now we have an orbital product |phi_a(1) phi_b(2); sa(1) sb(2)>
                # where phi_a = 1s (always m=0), phi_b = 2p_{m_b}.
                # Build the antisymmetrized SD.
                orb_a = ('1s', 0)
                orb_b = ('2p', m_b)
                # Antisymmetric combination:
                #   |ab>_AS = (1/sqrt(2)) * [|a(1)b(2)> - |b(1)a(2)>]
                # Store as SD with a canonical ordering.
                # To avoid double-counting, we choose the ordering by a key
                # on (orb, s) pair.
                p1 = (orb_a, sa)
                p2 = (orb_b, sb)
                if p1 == p2:
                    # Same spin-orbital: Pauli exclusion, SD = 0
                    continue
                # Canonicalize: sorted tuple is the key; track sign
                if p1 < p2:
                    key = (p1, p2)
                    sign = 1
                else:
                    key = (p2, p1)
                    sign = -1
                coeff = sign * total_cg / sp.sqrt(2)
                state[key] = state.get(key, 0) + coeff

    # Simplify and prune zeros
    state = {k: sp.simplify(v) for k, v in state.items() if sp.simplify(v) != 0}
    return state


def test_build_3P_J():
    """Verify the |^3P_J> states are unit-normalized."""
    for J in (0, 1, 2):
        s = build_3P_J_state(J)
        norm_sq = sum(sp.simplify(c * sp.conjugate(c)) for c in s.values())
        norm_sq_s = sp.simplify(norm_sq)
        print(f"  |^3P_{J}, M_J={J}>: ||.||^2 = {norm_sq_s}  ({len(s)} SDs)")


# ============================================================
# Two-body tensor operators: matrix elements in SD basis
# ============================================================

# For SS, we work with the operator
#   V_SS = sum_q (-1)^q [s1 x s2]^(2)_q * R^(2)_{-q}(r-hat_12) / r_12^3
# (choosing units such that the overall alpha^2 and normalization prefactor is
# extracted at the end).
#
# Using the multipole expansion of R^(2)_{-q}(r-hat_12) / r_12^3 in the basis
# of individual-electron spherical harmonics:
#   [Y^(k)(1) x Y^(k')(2)]^(2)_{-q} * f_kk'(r_<, r_>)
# with (k, k') triangular to 2.
#
# The "radial kernel" f_kk'(r_<, r_>) integrates against the hydrogenic
# orbital pair (1s(1) 2p(2)) or (2p(1) 1s(2)) to give a specific M^K where K
# is the TRUE multipole index (not the tensor rank).  For rank-2:
#   (k, k') = (0, 2): gives M^2 on the (0, 2) multipole channel
#   (k, k') = (2, 0): gives M^2 on the (2, 0) multipole channel
#   (k, k') = (1, 1): gives M^2 on the (1, 1) multipole channel
#   (k, k') = (2, 2): gives M^2 on the (2, 2) -> rank 2 via 9j
#
# Many references separate these contributions.  For rank-2 SS in (1s)(2p),
# the effective radial integral Q = integral of (r_<^2 / r_>^5) weighted by
# the orbital products collapses different (k, k') channels into M^k_dir
# and M^k_exch as follows:
#   DIRECT:   1s(1) 1s(1) ... 2p(2) 2p(2), need <0|Y^k1|0><1|Y^k2|1> combination
#             with (k1, k2) coupled to 2.  Non-vanishing: (k1=0, k2=2) and (k1=2, k2=0).
#             (k1=1, k2=1) gives zero for direct since <0|Y^1|0>=0.
#   EXCHANGE: 1s(1) 2p(1) ... 2p(2) 1s(2), need <0|Y^k1|1><1|Y^k2|0>.  Non-vanishing:
#             (k1=1, k2=1) only (the others need <0|Y^0|1>=0 or <1|Y^2|0>=0).
#
# This is the KEY simplification that makes the 9j algebra tractable.
# Below we do the full computation.


def rank_k_spatial_me(la, ma, lb, mb, lap, map_, lbp, mbp, k, q):
    """Matrix element of [Y^(k1)(1) x Y^(k2)(2)]^(k)_q summed over (k1, k2)
    triangular to k, for a specific radial channel.

    Returns a dict { (k1, k2): coefficient } giving the contribution to the
    spatial matrix element for each (k1, k2) choice.  The radial integral
    (M^k_dir or M^k_exch) is NOT included -- only the angular algebra.

    (la, ma): l1, m1 on e1 in bra
    (lap, map_): l1, m1 on e1 in ket
    (lb, mb): l2, m2 on e2 in bra
    (lbp, mbp): l2, m2 on e2 in ket
    """
    out = {}
    for k1 in range(0, abs(la - lap) + 3):  # Gaunt selection
        for k2 in range(0, abs(lb - lbp) + 3):
            if abs(k1 - k2) > k or k1 + k2 < k:
                continue
            # Clebsch-Gordan combining [Y^(k1) x Y^(k2)]^(k)_q
            # <k q | k1 q1 k2 q2> where q1 + q2 = q
            # Using: T^(k)_q = sum_{q1, q2} <k1 q1 k2 q2 | k q> T^(k1)_q1 T^(k2)_q2
            # <Y^(k1)_{q1}>: Gaunt matrix elements
            #
            # The spatial matrix element of the coupled operator is:
            # <lm|Y^(k)(space)|l'm'> where Y^(k) is the tensor spherical harmonic
            #
            # Using <l m | Y^(k)_q | l' m'> = Gaunt integral (Edmonds 4.6.3):
            #   = sqrt((2l'+1)(2k+1)/(4 pi (2 l + 1))) * <l' k 0 0 | l 0> <l' k m' q | l m>
            #
            # For CLEANER algebra, use the C^(k) Racah tensor:
            #   C^(k)_q = sqrt(4 pi / (2k+1)) * Y^(k)_q
            # <l m | C^(k)_q | l' m'> = (-1)^m sqrt((2l+1)(2l'+1)) *
            #                            (l k l'; -m q m') * (l k l'; 0 0 0)

            # 3j-based matrix element of C^(k1)_{q1} on electron 1:
            #   <la ma | C^(k1)_{q1} | lap map_>
            def c_me(l, m, kk, qq, lp, mp):
                if m - qq != mp:
                    return 0
                # Note: in the 3j, <l m | C^(k)_q | l' m'> has qq on the SECOND column.
                # sign convention: Edmonds (5.4.1):
                # <l m | C^(k)_q | l' m'> = (-1)^{l - m} * sqrt((2l+1)(2l'+1)) *
                #       wigner_3j(l, k, l', -m, q, m') * wigner_3j(l, k, l', 0, 0, 0)
                ph = (-1) ** (l - m)
                return ph * sqrt((2 * l + 1) * (2 * lp + 1)) * \
                       wigner_3j(l, kk, lp, -m, qq, mp) * \
                       wigner_3j(l, kk, lp, 0, 0, 0)

            # Sum over q1 + q2 = q
            coup = 0
            for q1 in range(-k1, k1 + 1):
                q2 = q - q1
                if abs(q2) > k2:
                    continue
                cg = clebsch_gordan(k1, k2, k, q1, q2, q)
                if sp.simplify(cg) == 0:
                    continue
                me1 = c_me(la, ma, k1, q1, lap, map_)
                me2 = c_me(lb, mb, k2, q2, lbp, mbp)
                if sp.simplify(me1) == 0 or sp.simplify(me2) == 0:
                    continue
                coup = coup + cg * me1 * me2
            coup = sp.simplify(coup)
            if coup != 0:
                out[(k1, k2)] = coup
    return out


def rank_k_spin_me(sa, sa_p, sb, sb_p, k, q):
    """Matrix element of [s_1 (x) s_2]^(k)_q coupled spin tensor.

    Uses Edmonds (5.4.1)-type formulas.  Spin operator s_i has
    <1/2 m | s_q | 1/2 m'> = sqrt(3/2) * (-1)^{1/2-m} * (1/2 1 1/2; -m q m')

    (rank-1 tensor in q basis).

    [s_1 (x) s_2]^(k)_q = sum_{q1, q2} <1 q1 1 q2 | k q> s_{1, q1} s_{2, q2}
    """
    half = Rational(1, 2)
    # <1/2 m | s_q | 1/2 m'>:
    def s_me(m, qq, mp):
        if m - qq != mp:
            return 0
        ph = (-1) ** (half - m)
        return ph * sqrt((2 * half + 1) * (2 * half + 1)) * \
               wigner_3j(half, 1, half, -m, qq, mp) * sqrt(Rational(1, 2) *
                       (half + 1) * (2 * half + 1))
        # Actually the reduced matrix element <1/2 || s || 1/2> = sqrt(3/2)
        # Using: <j m | T^(1)_q | j' m'> = <j m | j' m' 1 q> / sqrt(2j+1) * <j||T||j'>
        # I'm conflating conventions; let's use a cleaner route

    # Use direct sympy evaluation by invoking the rank-1 reduction:
    # s = sigma/2 so <1/2 m|s_q|1/2 m'> = (1/2) <1/2 m|sigma_q|1/2 m'>
    # and the Pauli-matrix element <1/2 m|sigma_q|1/2 m'> = sqrt(2) (-1)^{1/2-m} <1/2 -m 1/2 m'|1 q>
    # (via Condon-Shortley sigma_q = sqrt(2) * T^(1)_q for Pauli matrices)
    def s_me_clean(m, qq, mp):
        if m - qq != mp:
            return 0
        # Reduced m.e. <1/2||s||1/2> = sqrt(3)/2 in Condon-Shortley
        # Full: <1/2 m|s_q|1/2 m'> = <1/2 m'| 1/2 m 1 (-q)>(-1)^(1-q) * <1/2||s||1/2> ?
        # Easier: directly evaluate s_q on spin-1/2.
        # s_+1 = -(s_x + i s_y)/sqrt(2); s_0 = s_z; s_-1 = (s_x - i s_y)/sqrt(2)
        # In |1/2 +1/2>, |1/2 -1/2> basis:
        # s_z |+> = 1/2 |+>, s_z |-> = -1/2 |->
        # s_+ |-> = |+>  (where s_+ = s_x + i s_y)
        # So s_{+1} = -s_+/sqrt(2): <+|s_{+1}|-> = -1/sqrt(2)
        # s_-1 = s_-/sqrt(2): <-|s_{-1}|+> = 1/sqrt(2)
        if qq == 0:
            return m if m == mp else 0
        elif qq == 1:
            # s_{+1} lowers m by 1 but wait -- s_{+1} RAISES the state...
            # T^(1)_{+1} applied to |j m>: gives coefficient for |j m+1>.
            # So <m|s_{+1}|m'> has m = m' + 1.
            # For m = +1/2, m' = -1/2:
            if m == Rational(1, 2) and mp == Rational(-1, 2):
                return -1 / sqrt(2)
            return 0
        elif qq == -1:
            # s_{-1} lowers m by 1:  <m|s_{-1}|m'> has m = m' - 1.
            if m == Rational(-1, 2) and mp == Rational(1, 2):
                return 1 / sqrt(2)
            return 0
        return 0

    # Now [s_1 x s_2]^(k)_q = sum_{q1+q2=q} <1 q1 1 q2 | k q> s_{1,q1} s_{2,q2}
    result = 0
    for q1 in range(-1, 2):
        q2 = q - q1
        if abs(q2) > 1:
            continue
        cg = clebsch_gordan(1, 1, k, q1, q2, q)
        if sp.simplify(cg) == 0:
            continue
        me1 = s_me_clean(sa, q1, sa_p)
        me2 = s_me_clean(sb, q2, sb_p)
        if sp.simplify(me1) == 0 or sp.simplify(me2) == 0:
            continue
        result = result + cg * me1 * me2
    return sp.simplify(result)


# ============================================================
# Step 3: Full diagonal matrix element <^3P_J | H_SS | ^3P_J>
# ============================================================

def compute_diagonal_SS(J, operator="SS"):
    """Compute <^3P_J, M_J=J | H_{SS or SOO} | ^3P_J, M_J=J>.

    operator: "SS" (rank 2) or "SOO" (rank 1)

    Returns the coefficient expressed symbolically in terms of M^k_dir, M^k_exch
    and numerical rational coefficients.
    """
    if operator == "SS":
        k_tensor = 2
    elif operator == "SOO":
        k_tensor = 1
    else:
        raise ValueError(operator)

    state = build_3P_J_state(J)
    result = 0
    for key_bra, c_bra in state.items():
        for key_ket, c_ket in state.items():
            # Each key is a 2-electron SD: ((orb1, s1), (orb2, s2)) canonical.
            # Compute <bra|V|ket> where V is the 2-body tensor operator.
            me = me_two_body(key_bra, key_ket, k_tensor, q=0)
            # ... this is getting into the weeds.  Let me reformulate.
            pass

    raise NotImplementedError("See dd_drake_derivation_v2.py for the full derivation")


def main():
    print("=" * 76)
    print("Sprint 4 Track DD: Drake 1971 coefficient derivation (v1)")
    print("=" * 76)
    verify_J_pattern()
    print("\n" + "=" * 76)
    print("STEP 2: |^3P_J> state construction in SD basis")
    print("=" * 76)
    test_build_3P_J()


if __name__ == "__main__":
    main()
