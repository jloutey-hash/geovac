"""Sprint 4 Track DD: Drake 1971 coefficients via direct Slater-determinant
matrix element expansion.

Approach
--------
We build |^3P_J, M_J = J> as an antisymmetrized linear combination of
two-electron SDs in the |1s(1) 2p_m(2)> basis. For each J, we then compute
<^3P_J, M_J = J | H_{SS or SOO} | ^3P_J, M_J = J> by summing over pairs of
SDs and applying the two-body matrix element formula.

For each matrix element <ab | V | cd>, we use the multipole expansion of
V(1, 2) in terms of single-electron spherical tensor operators, yielding a
sum over multipole indices K weighted by the radial Drake integral M^K and
the angular (3j/9j) coefficient.

The target is to derive
  A_SS  = alpha^2 * ( 3/50 * M^2_dir - 2/5 * M^2_exch )
  A_SOO = alpha^2 * ( 3/2  * M^1_dir -       M^1_exch )
symbolically.

Strategy: identify the |^3P_1, M_J=1> state, compute <^3P_1|H|^3P_1>, and
extract the coefficients of M^k_dir and M^k_exch.  These ARE A_SS and A_SOO
since f_SS(J=1) = f_SOO(J=1) = 1.
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
from sympy import Rational, sqrt, Integer, simplify, symbols, pi, I as sp_I
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan


# Symbolic radial integrals
M0d, M1d, M2d = symbols("M0d M1d M2d", real=True)
M0e, M1e, M2e = symbols("M0e M1e M2e", real=True)
alpha_sym = symbols("alpha", positive=True, real=True)

# Map (K, channel) -> sympy symbol
M_sym = {
    (0, 'd'): M0d, (1, 'd'): M1d, (2, 'd'): M2d,
    (0, 'e'): M0e, (1, 'e'): M1e, (2, 'e'): M2e,
}


def reduced_C(l, lp, k):
    """<l || C^(k) || l'> (Racah normalized)."""
    three_j = wigner_3j(l, k, lp, 0, 0, 0)
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * three_j


def c_me(l, m, k, q, lp, mp):
    """<l m | C^(k)_q | l' m'>."""
    if m - q != mp:
        return Integer(0)
    return (-1) ** (l - m) * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * \
           wigner_3j(l, k, lp, -m, q, mp) * \
           wigner_3j(l, k, lp, 0, 0, 0)


def s_me(m, q, mp):
    """<1/2 m | s_q | 1/2 m'> where s = sigma/2 (spin operator).

    s_{+1} = -(s_x + i s_y)/sqrt(2) on spin-1/2:
      <1/2 +1/2 | s_{+1} | 1/2 -1/2> = -1/(2 sqrt(2))
      (rank-1 spherical tensor, spin-1/2 RME <1/2||s||1/2> = sqrt(3/2))
    """
    half = Rational(1, 2)
    m, mp = sp.Rational(m), sp.Rational(mp)
    if m - q != mp:
        return Integer(0)
    if q == 0:
        return m if m == mp else Integer(0)
    elif q == 1:
        # s_{+1} raises m by 1 via spherical tensor convention
        # <m = +1/2 | s_{+1} | mp = -1/2> = -1/(2 sqrt(2))
        if m == half and mp == -half:
            return -Integer(1) / (Integer(2) * sqrt(2))
        return Integer(0)
    elif q == -1:
        # s_{-1} lowers m by 1
        if m == -half and mp == half:
            return Integer(1) / (Integer(2) * sqrt(2))
        return Integer(0)
    return Integer(0)


# ============================================================
# Two-body matrix element via multipole expansion
# ============================================================
#
# For the SS operator:
#   H_SS = -alpha^2 sqrt(24 pi / 5) * sum_Q (-1)^Q * [s_1 (x) s_2]^(2)_Q
#           * Y^(2)_{-Q}(r-hat_12) / r_12^3
#
# We factorize Y^(2)_Q(r-hat_12) / r_12^3 in single-electron spherical tensor
# components.  The key identity (Brink & Satchler, "Angular Momentum" App. III):
#
# Y^(K)_Q(r-hat_12) / r_12^{K+1}
#   = sum_{k1 k2, triangle(k1,k2,K)} T_{k1 k2}^K(r_1, r_2)
#     * [Y^(k1)(1-hat) x Y^(k2)(2-hat)]^(K)_Q
#
# where T_{k1 k2}^K has a clean structure.  For the specific case K = 2 and
# the Breit-Pauli denominator r_12^3 (which equals r_12^{K+1} for K = 2):
#   T_{k1 k2}^2(r_1, r_2) = (4 pi / 5) * specific radial kernel
#
# And the specific radial kernel evaluates between hydrogenic orbitals to
# combinations of Drake's M^k integrals.  The exact formulas are:
#
# For K = 2 (SS rank-2) with Breit-Pauli kernel 1/r_12^3:
#   T_{0, 2}^2(r_1, r_2) = 1 / r_2^3  (for r_1 < r_2) or r_1^2 / r_2^5 (r_1 > r_2)?
# Wait, this requires careful work.  Let me instead take a KEY SHORTCUT:
#
# The angular projection of 1/r_12^3 * Y^(K)(r-hat_12) onto the two-electron
# orbital basis gives a specific linear combination of the Drake integrals.
# For the (1s)(2p) configuration at He, the explicit Drake 1971 formulas
# (Eqs. 17) ARE:
#
#   <(1s)(2p) ^3P_J | H_SS | ...> = (-alpha^2) * [
#      (specific coefficient) * sum_k f_SS(J) * M^2_{dir, or exch}
#   ]
#
# This is the direct-minus-exchange pattern with specific 9j coefficients.


# ============================================================
# NEW STRATEGY: do the matrix element in finite-dimensional m-basis
# ============================================================
#
# The observation: f_SS(J=1) = 1 and f_SOO(J=1) = 1.  So
#   A_SS  = <^3P_1, M_J=1 | H_SS | ^3P_1, M_J=1>  / alpha^2
#   A_SOO = <^3P_1, M_J=1 | H_SOO | ^3P_1, M_J=1> / alpha^2
#
# We compute these matrix elements by:
#   1. Writing |^3P_1, M_J=1> as a weighted sum of 2-electron product states
#      |1s m_a, s_a; 2p m_b, s_b> (NOT antisymmetrized yet).
#   2. Applying the 2-body Breit-Pauli operator term-by-term in this basis.
#   3. Collecting the direct-minus-exchange contributions.
#   4. Identifying which M^k_dir and M^k_exch each contribution corresponds to.
#
# For the SS operator, the tensor form is
#   H_SS = -alpha^2 * sum_q (-1)^q * sqrt(24 pi / 5) * [s_1 x s_2]^(2)_q *
#           Y^(2)_{-q}(r-hat_12) / r_12^3
#
# The angular matrix element of Y^(2)_{-q}(r-hat_12) / r_12^3 between
# |1s m_a ; 2p m_b> and |1s m_a' ; 2p m_b'> (DIRECT) or between
# |1s m_a ; 2p m_b> and |2p m_a' ; 1s m_b'> (EXCHANGE) is
#   DIRECT:   delta_{m_a, m_a'}(l_a=0) * <2p m_b | Y^(2)_{-q}(hat r) | 2p m_b'> / r_2^3
#             ... wait this doesn't separate like that.
#
# The RIGHT way: use the multipole expansion
#   Y^(K)_Q(r-hat_12) / r_12^3
#     = sum_{k1 k2} A(k1, k2; K) * C^(k1)_q1(1) * C^(k2)_q2(2) * <k1 q1 k2 q2 | K Q>
#       * g_{k1 k2}^K(r_1, r_2)
# where g_{k1 k2}^K is a specific radial kernel.  The Drake integrals M^k
# are just the RADIAL integrals after factoring out angular parts:
#   M^K_dir  = <1s 1s | g^K | 2p 2p> radial
#   M^K_exch = <1s 2p | g^K | 2p 1s> radial
# where the same K appears because the tensor has total rank K.
#
# The multipole decomposition (from standard angular algebra):
#
# Y^(K)_Q(r-hat_12) / r_12^{K+1}
#   = sqrt(4 pi / (2K+1)) * sum_{k1 k2 triangle(k1,k2,K)}
#        sqrt((2k1+1)(2k2+1)/(4pi)) * <k1 0 k2 0 | K 0> *
#        [C^(k1)(1) x C^(k2)(2)]^(K)_Q * r_<^K / r_>^{K+1}
#
# (Varshalovich Quantum Theory of Angular Momentum, Eq 5.17.2, adapted).
#
# For the Breit-Pauli factor 1/r_12^{K+2} at K=2 (so kernel 1/r_12^{K+3}=1/r_12^5),
# the radial piece is r_<^K / r_>^{K+3} = r_<^2 / r_>^5.  Good.
#
# Putting it together for the SS tensor with prefactor -sqrt(24 pi/5):
#
#   -sqrt(24 pi / 5) * Y^(2)_{-Q}(r-hat_12) / r_12^3
#     = -sqrt(24 pi/5) * sqrt(4 pi / 5) * sum_{k1 k2 tri(k1,k2,2)}
#        sqrt((2k1+1)(2k2+1)/(4pi)) * <k1 0 k2 0 | 2 0> *
#        [C^(k1)(1) x C^(k2)(2)]^(2)_{-Q} * r_<^2 / r_>^5
#     = -sqrt(24 pi^2 * 4 / 25) * sum ...
#     = -(4 pi sqrt(6) / 5) * sum ...
#     = -4 pi sqrt(6) / 5 * sum_{k1 k2}
#        sqrt((2k1+1)(2k2+1)/(4pi)) <k1 0 k2 0 | 2 0> [C(1) x C(2)]^(2)_{-Q}
#        * r_<^2 / r_>^5
#
# The 4pi factor will cancel against the 1/(4pi) normalization in the Drake
# radial integral convention.  Key: the angular coefficient is
#   A(k1, k2, K=2) = sqrt((2k1+1)(2k2+1)) * <k1 0 k2 0 | 2 0>
# (with appropriate 4pi factors absorbed).
#
# We take a simpler route: just use the Y^(K)(r-hat_12) multipole expansion
# form directly, plug into the angular matrix element in |l_a m_a> basis,
# and evaluate.
#

# ============================================================
# Direct computation: build matrix elements term by term
# ============================================================

def multipole_coef(k1, k2, K, Q):
    """CG combining [C^(k1)(1) x C^(k2)(2)]^(K)_Q.

    Y^(K)_Q(r-hat_12) / r_12^{K+1} expands as:
      sum_{k1 k2 tri(k1,k2,K)}  B(k1, k2, K) [C^(k1)(1) x C^(k2)(2)]^(K)_Q *
        r_<^K / r_>^{K+1}   [Coulomb kernel]
    or
      sum ... * r_<^K / r_>^{K+3}   [Breit-Pauli kernel]

    The angular coefficient B(k1, k2, K) = <k1 0 k2 0 | K 0> up to normalization.
    """
    # <k1 0 k2 0 | K 0> is a specific CG coefficient
    return clebsch_gordan(k1, k2, K, 0, 0, 0)


def drake_direct_matrix_element(J_ket, K):
    """Compute <^3P_J, M_J=J_ket | [T^(K)(space) . U^(K)(spin)]^(0) | same>
    where the spin tensor is [s_1 x s_2]^(K).

    Returns a sympy expression in M_sym symbols.
    """
    half = Rational(1, 2)

    # |^3P_J> expansion: sum over (M_L, M_S) with M_L + M_S = J_ket (=M_J)
    # and for each (M_L, M_S), expand spatial (l_a=0, l_b=1 -> L=1) and spin.
    # Each term is a spin-orbital product |1s m_a ; 2p m_b ; s_a, s_b>

    # Build up ORBITAL PRODUCT (not antisymmetrized) state:
    product_state = {}   # key: (m_a, m_b, s_a, s_b), value: coefficient
    for M_L in (-1, 0, 1):
        for M_S in (-1, 0, 1):
            if M_L + M_S != J_ket:
                continue
            cg_LSJ = clebsch_gordan(1, 1, J_ket, M_L, M_S, J_ket)
            if sp.simplify(cg_LSJ) == 0:
                continue
            # Spatial
            m_a, m_b = 0, M_L
            cg_space = clebsch_gordan(0, 1, 1, 0, M_L, M_L)
            if sp.simplify(cg_space) == 0:
                continue
            # Spin
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
                    key = (m_a, m_b, sa, sb)
                    product_state[key] = product_state.get(key, Integer(0)) + total

    product_state = {k: sp.simplify(v) for k, v in product_state.items() if sp.simplify(v) != 0}

    # For triplet (S=1, spatial antisymmetric), the ACTUAL wavefunction is
    #   |psi>_AS = (1/sqrt(2)) [|1s m_a ; 2p m_b> - |2p m_b ; 1s m_a>] * spin_triplet
    # So <psi_AS | V | psi_AS> = direct - exchange = <a b | V | a b> - <a b | V | b a>

    # Now compute <product_state | V | product_state>_direct and _exchange.
    # V = [T^(K)(space) . U^(K)(spin)]^(0) = sum_q (-1)^q / sqrt(2K+1) T^(K)_q U^(K)_{-q}
    # (the scalar product normalization for rank-K tensors)

    # T^(K)(space) = sum_{k1 k2, Q} A(k1 k2, K, Q) [C^(k1)(1) x C^(k2)(2)]^(K)_Q * kernel
    # U^(K)(spin) = [s_1 (x) s_2]^(K)

    me_total = Integer(0)

    # Direct contribution: <m_a m_b sa sb | V | m_a' m_b' sa' sb'> where the
    # orbitals on each electron are unchanged (1s on e1, 2p on e2).
    for (ma, mb, sa, sb), c_bra in product_state.items():
        for (map_, mbp, sap, sbp), c_ket in product_state.items():
            # DIRECT: both electrons keep their orbitals.
            # <1s m_a sa; 2p m_b sb | V | 1s m_a' sa'; 2p m_b' sb'>
            # V = sum_q (-1)^q / sqrt(2K+1) * T^(K)_q * U^(K)_{-q}
            # T^(K)_q = sum_{k1 k2} A(k1 k2 K) [C^(k1)(1) x C^(k2)(2)]^(K)_q * M^K_dir
            # U^(K)_{-q} = [s_1 x s_2]^(K)_{-q}

            # Sum over Q = q:
            me_dir = Integer(0)
            for Q in range(-K, K + 1):
                # Angular matrix element of spatial tensor:
                ang = Integer(0)
                for k1 in range(0, K + 1 + 1):
                    for k2 in range(0, K + 1 + 1):
                        if abs(k1 - k2) > K or k1 + k2 < K:
                            continue
                        # CG combining k1 and k2 into K
                        # [C^(k1)(1) x C^(k2)(2)]^(K)_Q = sum_{q1 q2} <k1 q1 k2 q2 | K Q> C^(k1)_q1(1) C^(k2)_q2(2)
                        for q1 in range(-k1, k1 + 1):
                            q2 = Q - q1
                            if abs(q2) > k2:
                                continue
                            cg = clebsch_gordan(k1, k2, K, q1, q2, Q)
                            if sp.simplify(cg) == 0:
                                continue
                            # <1s m_a | C^(k1)_q1 | 1s m_a'> for l=0 -> nonzero only for k1=0, q1=0
                            me1 = c_me(0, ma, k1, q1, 0, map_)
                            # <2p m_b | C^(k2)_q2 | 2p m_b'>
                            me2 = c_me(1, mb, k2, q2, 1, mbp)
                            if sp.simplify(me1) == 0 or sp.simplify(me2) == 0:
                                continue
                            # Multipole coefficient from Y^(K) expansion
                            # A(k1, k2, K) = sqrt((2k1+1)(2k2+1)) * <k1 0 k2 0 | K 0>
                            # (Varshalovich Eq. 5.17.2)
                            A_k1k2K = sqrt(Integer((2 * k1 + 1) * (2 * k2 + 1))) * \
                                      clebsch_gordan(k1, k2, K, 0, 0, 0)
                            ang = ang + cg * me1 * me2 * A_k1k2K
                ang = sp.simplify(ang)
                if ang == 0:
                    continue
                # Spin tensor m.e.: <sa sb | [s_1 x s_2]^(K)_{-Q} | sa' sb'>
                spin = Integer(0)
                for q1s in range(-1, 2):
                    q2s = -Q - q1s
                    if abs(q2s) > 1:
                        continue
                    cg_s = clebsch_gordan(1, 1, K, q1s, q2s, -Q)
                    if sp.simplify(cg_s) == 0:
                        continue
                    me_s1 = s_me(sa, q1s, sap)
                    me_s2 = s_me(sb, q2s, sbp)
                    if sp.simplify(me_s1) == 0 or sp.simplify(me_s2) == 0:
                        continue
                    spin = spin + cg_s * me_s1 * me_s2
                spin = sp.simplify(spin)
                if spin == 0:
                    continue
                # Phase from scalar product normalization: (-1)^Q
                me_dir = me_dir + (-1) ** Q * ang * spin
            me_dir = sp.simplify(me_dir)
            # Multiply by the DIRECT radial integral M^K_dir (symbolic) and the
            # state coefficients:
            if me_dir != 0:
                me_total = me_total + sp.conjugate(c_bra) * c_ket * me_dir * M_sym[(K, 'd')]

            # EXCHANGE: orbitals swapped on electrons 1 and 2.
            # <1s m_a sa; 2p m_b sb | V | 2p m_b' sb'; 1s m_a' sa'>  (with exchange of e1<->e2)
            # Effectively: <1s m_a| C^(k1) |2p m_b'> on electron 1
            #              <2p m_b| C^(k2) |1s m_a'> on electron 2
            # Similarly for spin.

            me_exc = Integer(0)
            for Q in range(-K, K + 1):
                ang_exc = Integer(0)
                for k1 in range(0, K + 1 + 1):
                    for k2 in range(0, K + 1 + 1):
                        if abs(k1 - k2) > K or k1 + k2 < K:
                            continue
                        for q1 in range(-k1, k1 + 1):
                            q2 = Q - q1
                            if abs(q2) > k2:
                                continue
                            cg = clebsch_gordan(k1, k2, K, q1, q2, Q)
                            if sp.simplify(cg) == 0:
                                continue
                            # Exchange: electron 1 is 1s-to-2p, electron 2 is 2p-to-1s
                            me1 = c_me(0, ma, k1, q1, 1, mbp)   # <1s m_a|C^(k1)|2p m_b'>
                            me2 = c_me(1, mb, k2, q2, 0, map_)  # <2p m_b|C^(k2)|1s m_a'>
                            if sp.simplify(me1) == 0 or sp.simplify(me2) == 0:
                                continue
                            A_k1k2K = sqrt(Integer((2 * k1 + 1) * (2 * k2 + 1))) * \
                                      clebsch_gordan(k1, k2, K, 0, 0, 0)
                            ang_exc = ang_exc + cg * me1 * me2 * A_k1k2K
                ang_exc = sp.simplify(ang_exc)
                if ang_exc == 0:
                    continue
                # Spin exchange: note spin labels are ALSO swapped in exchange
                # (we exchange both orbitals AND spins on electrons 1 and 2).
                spin_exc = Integer(0)
                for q1s in range(-1, 2):
                    q2s = -Q - q1s
                    if abs(q2s) > 1:
                        continue
                    cg_s = clebsch_gordan(1, 1, K, q1s, q2s, -Q)
                    if sp.simplify(cg_s) == 0:
                        continue
                    # Exchange spin: <sa| s_1^{q1}  |sbp>, <sb| s_2^{q2} |sap>
                    # Wait -- in the exchange, ELECTRONS swap, so spin slots
                    # become: electron 1 now has spin (sb), electron 2 has (sa)
                    # So <sa sb | V(exchange) | sa' sb'> =
                    #   sum_q T_q(space, 1 with ket b' vs bra a) T_-q(space, 2 with ket a' vs bra b) *
                    #         U_q1(spin, electron 1: <sa|..|sb'>) U_q2(spin, electron 2: <sb|..|sa'>)
                    me_s1 = s_me(sa, q1s, sbp)    # <sa | s_1^{q1}_{(electron 1)} | sb'>
                    me_s2 = s_me(sb, q2s, sap)    # <sb | s_2^{q2}_{(electron 2)} | sa'>
                    if sp.simplify(me_s1) == 0 or sp.simplify(me_s2) == 0:
                        continue
                    spin_exc = spin_exc + cg_s * me_s1 * me_s2
                spin_exc = sp.simplify(spin_exc)
                if spin_exc == 0:
                    continue
                me_exc = me_exc + (-1) ** Q * ang_exc * spin_exc
            me_exc = sp.simplify(me_exc)
            if me_exc != 0:
                # Minus sign for exchange (spatial antisymmetrization for triplet):
                me_total = me_total - sp.conjugate(c_bra) * c_ket * me_exc * M_sym[(K, 'e')]

    # The scalar tensor product normalization [T x U]^(0) = sum_q (-1)^q T_q U_{-q} / sqrt(2K+1)
    me_total = me_total / sqrt(Integer(2 * K + 1))

    return sp.simplify(me_total)


def main():
    print("=" * 76)
    print("Sprint 4 Track DD: Direct Slater-determinant matrix element")
    print("=" * 76)

    # Compute <^3P_1, M_J=1 | H_SS | ^3P_1, M_J=1>  (K=2)
    for K, label in [(2, "SS"), (1, "SOO")]:
        print(f"\n  {label} (K = {K}):")
        for J in (0, 1, 2):
            me = drake_direct_matrix_element(J, K)
            print(f"    <^3P_{J}, M_J = {J} | H_{label} | ^3P_{J}, M_J = {J}> = {me}")
            # Collect coefficients of M^K_dir and M^K_exch:
            # me = c_dir * M_sym[(K, 'd')] + c_exc * M_sym[(K, 'e')]
            c_d = me.coeff(M_sym[(K, 'd')])
            c_e = me.coeff(M_sym[(K, 'e')])
            # Residual after removing M^K_dir, M^K_exch contributions
            res = sp.simplify(me - c_d * M_sym[(K, 'd')] - c_e * M_sym[(K, 'e')])
            print(f"      coef of M^{K}_dir: {c_d}")
            print(f"      coef of M^{K}_exch: {c_e}")
            print(f"      residual: {res}")


if __name__ == "__main__":
    main()
