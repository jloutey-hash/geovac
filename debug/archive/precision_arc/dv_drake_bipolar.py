"""Sprint 5 Track DV: derive Drake 1971 mixing ratios (3/50, -2/5, 3/2, -1)
from the bipolar harmonic expansion of 1/r_12^{K+1} * Y^K(r_hat_12).

Strategy
========
Sprint 4 DD derived the J-pattern f_SS/SOO via 6j, but the direct/exchange
mixing ratios remained conjectural. The obstruction was that DD used a
"universal" radial kernel r_<^K / r_>^{K+3} across all (k_1, k_2) channels.
Varshalovich §5.17 and Brink-Satchler App. 5 show that the bipolar expansion

    Y^K_Q(r_hat_12) / r_12^{K+1}
        = sum_{k1,k2} F^{k1 k2 K}(r_1, r_2) · [Y^{k1}(hat r_1) (x) Y^{k2}(hat r_2)]^K_Q

has a DIFFERENT radial kernel F for each (k1, k2) — not a universal one.

For our target (Breit-Pauli 1/r_12^3 with Y^2 for SS, Y^1 for SOO), we
restrict to the (k1, k2) channels allowed by the He (1s)(2p) angular
selection rules:

  K=2 (SS):
    Direct path   (1s,1s → 2p,2p): (k1=0, k2=2) is the only allowed triangle
    Exchange path (1s,2p → 2p,1s): (k1=1, k2=1) is the only allowed triangle

  K=1 (SOO):
    Direct  (1s,1s → 2p,2p): (0, 1)
    Exchange(1s,2p → 2p,1s): (1, 0) or (1, 2)

Each (k1,k2) channel has a closed-form radial integral in the hydrogenic
1s, 2p densities that factorizes as (constant) × M^K_{dir/exch} — where the
CONSTANT is the mixing ratio we want to derive.

KEY IDEA (Wigner-Eckart + 9j recoupling)
=========================================
The angular matrix element of the full rank-K scalar operator

    V = [T^(K)(space) · U^(K)(spin)]^(0)

on the LS-coupled triplet |(l_a l_b) L=1, S=1; J, M_J> can be computed by
Racah algebra:

  <L S J M | V | L S J M>
      = (-1)^{L+S+J} 6j{L S J; S L K} <L||T^(K)||L> <S||U^(K)||S>

where <L||T^(K)||L> is the SPATIAL reduced matrix element (sum over all
allowed (k1, k2) channels with their radial kernels), and <S||U^(K)||S>
is the SPIN reduced matrix element.

The spatial reduced matrix element for a bipolar tensor:

  <(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(K) || (l_c l_d) L'>
      = sqrt((2L+1)(2L'+1)(2K+1)) 9j{l_a l_b L; l_c l_d L'; k1 k2 K}
          × <l_a||C^(k1)||l_c> <l_b||C^(k2)||l_d>

Putting it all together: for each (k1, k2) channel contributing to the
direct or exchange path, we get a specific 9j symbol AND a specific radial
kernel. The radial kernel must be extracted from the bipolar expansion of
1/r_12^{K+1} * Y^K.

The bipolar expansion (Varshalovich, Quantum Theory of Angular Momentum,
Eq. 5.17.4, for the "generalized Coulomb" kernel):

  Y^K_Q(hat r_12) / r_12^{K+1}
      = (2K+1)!! / sqrt(4pi (2K+1)) sum_{k1,k2}
          δ(k1+k2, K)  <k1 0 k2 0 | K 0>   (triangular, both parities match)
          × [Y^(k1)(1) (x) Y^(k2)(2)]^(K)_Q
          × r_<^{K - k2+ k1 mod ...} / r_>^{K + 1 + ...}

Wait — the r_<, r_> form is for the COULOMB expansion of 1/r_12 (K=0).
For higher K, the bipolar expansion of Y^K(r_hat_12)/r_12^{K+1} has a
different structure. Let me derive it more carefully.

SACK 1964 EXPANSION (correct form)
===================================

Sack, J. Math. Phys. 5, 245 (1964): for a function f(r_12), the bipolar
expansion is

  f(r_12) = sum_{k1, k2, k3} <k1 0 k2 0 | k3 0> * bipolar_coef(k1, k2, k3)
            * [Y^{k1}(hat r_1) (x) Y^{k2}(hat r_2)]^{k3}_0 * R(k1, k2, k3; r_1, r_2)

For the specific case f(r_12) * Y^K(r_hat_12), the combined object is a
bipolar operator. The key identity for the rank-K case:

  Y^K_M(r_hat_12) * r_12^n = sum_{k1, k2} (selection rules)
                              * r_1^{a(k1,k2,n)} * r_2^{b(k1,k2,n)}
                              * [Y^{k1} (x) Y^{k2}]^K_M

up to rearrangement. But this is only the "regular" part; the actual
expansion has a two-region structure at r_1 <-> r_2.

Practical algorithm (to close Sprint 5)
========================================

We sidestep the general bipolar machinery by computing the He 2^3P matrix
elements DIRECTLY in the |m_L m_S> basis, with the full **position-space
integral** for the tensor operator, and then express the result in terms
of Drake's M^K_{dir/exch}. The radial integrals are evaluated by the
production `geovac.breit_integrals` module; only ANGULAR factors need to
be derived symbolically.

Step 1: Write H_SS = -alpha^2 sqrt(24 pi / 5) sum_q (-1)^q [s1 x s2]^(2)_q
                      * Y^(2)_{-q}(hat r_12) / r_12^3

Step 2: Compute the angular integral
          <l_a m_a; l_b m_b | Y^(2)_{-q}(hat r_12) / r_12^3 | l_c m_c; l_d m_d>
        as a function of (r_1, r_2) using the bipolar expansion.

Step 3: Contract with the spin matrix elements and the |^3P_J M_J=J>
        coefficients.

Step 4: Identify the coefficient of M^2_dir and M^2_exch in the resulting
        sum.

For the He (1s)(2p) configuration, the triangle constraints from both
spatial and spin selection rules are very restrictive; the bipolar
decomposition has only a small number of terms that need to be tracked.

This script is a fresh derivation — it does NOT import Sprint 4 DD's
(partially broken) direct-SD code. We work entirely in exact sympy
symbolic arithmetic.
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, symbols, pi, Symbol
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan

# ======================================================================
# Symbolic primitives
# ======================================================================

# Drake's radial integrals (K, channel) -> symbol
M_sym = {}
for K in (0, 1, 2):
    M_sym[(K, 'dir')] = Symbol(f"M{K}d", real=True)
    M_sym[(K, 'exch')] = Symbol(f"M{K}e", real=True)


def reduced_C(l: int, lp: int, k: int):
    """Reduced matrix element <l || C^(k) || l'> of the Racah-normalized
    spherical harmonic tensor C^(k)_q = sqrt(4pi/(2k+1)) Y^(k)_q.

    Condon-Shortley / Racah form:
      <l || C^(k) || l'> = (-1)^l sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0)
    """
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * \
           wigner_3j(l, k, lp, 0, 0, 0)


def ninej(*args):
    return simplify(wigner_9j(*args))


def sixj(*args):
    return simplify(wigner_6j(*args))


# ======================================================================
# Bipolar expansion of Y^K_Q(hat r_12) / r_12^{K+1}
# ======================================================================
#
# Following Steinborn & Filter (1975); Varshalovich Eq. 5.17.1-5:
#
# The fundamental object is the Y^K_Q(r_hat_12) * (coefficient of r_12)
# form. For generalized Coulomb kernel Y^K_Q(r_hat_12)/r_12^{K+1} the
# expansion is:
#
# Y^K_Q(hat r_12) / r_12^{K+1}
#   = sum_{L=K,K+2,K+4,...} (-1)^{(L-K)/2} B(K, L; r_1, r_2)
#     * sum_{k1+k2=L, |k1-k2|<=K<=k1+k2}
#       sqrt((2k1+1)(2k2+1)/(4pi(2K+1)))
#       * <k1 0 k2 0 | K 0> * [Y^{k1}(1) (x) Y^{k2}(2)]^{K}_Q
#
# where B(K, L; r_1, r_2) has a piecewise form in (r_<, r_>).
#
# However, for THIS target problem (matrix element between (1s, 2p) on
# both sides at configuration (l_a, l_b) = (0, 1)), angular selection
# rules on <l=0 || C^(k1) || l'> and <l=1 || C^(k2) || l'> restrict
# (k1, k2) to very few values.
#
# KEY INSIGHT: the ONLY (k1, k2) contributing to the DIRECT path
# (<1s|*|1s> on electron 1, <2p|*|2p> on electron 2) is:
#   k1 = 0 (from <l=0||C^(k1)||l=0>  ⇒ k1=0)
#   k2 ∈ {0, 2} (from <l=1||C^(k2)||l=1>, parity preserving ⇒ k2 even ≤ 2)
#
# For K=2 triangle |k1-k2| <= K <= k1+k2: k1=0, k2=2 only.  (k2=0 fails
# triangle with K=2.)
#
# The ONLY (k1, k2) contributing to the EXCHANGE path
# (<1s|*|2p> on electron 1, <2p|*|1s> on electron 2):
#   k1 = 1 (from <0||C^(k1)||1>  ⇒ k1=1)
#   k2 = 1 (from <1||C^(k2)||0>  ⇒ k2=1)
#
# Triangle for K=2: k1=1, k2=1 gives k1+k2=2, |k1-k2|=0, so K=2 is on the
# upper edge.  Parity: (-1)^(1+1+2) = +1, OK.
#
# For K=1 triangle:
#   Direct: k1=0, k2 ∈ {0, 2} with K=1 in triangle ⇒ NONE (k1=0 and k2 even
#     can't form triangle with K=1). So the DIRECT K=1 channel is EMPTY
#     from angular selection alone.
#
# But Drake's formula has A_SOO = 3/2 M^1_dir - M^1_exch — a non-zero
# direct coefficient! This means the SOO operator is not pure rank-1
# spatial Y^1/r^2; it includes the MOMENTUM operator p_1 or p_2 which
# can reshuffle the angular structure.

# --- Build the angular integral <l_a m_a | C^(k1)_q1 | l_c m_c> ---

def c_spherical_me(l, m, k, q, lp, mp):
    """<l m | C^(k)_q | l' m'>, via Wigner-Eckart.

    Uses (-1)^{l-m} sqrt((2l+1)(2l'+1)) (l k l'; -m q m')(l k l'; 0 0 0).
    """
    if m - q != mp:
        return Integer(0)
    return ((-1) ** (l - m) * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) *
            wigner_3j(l, k, lp, -m, q, mp) * wigner_3j(l, k, lp, 0, 0, 0))


# ======================================================================
# Part 1: Direct computation of the spatial reduced matrix element
#          on the coupled |(l_a l_b) L> basis, using the 9j identity.
# ======================================================================
#
# Brink-Satchler Eq 4.6.8 (or Edmonds 7.1.6):
#
# <(l_a l_b) L || [A^(k1)(1) (x) B^(k2)(2)]^(K) || (l_c l_d) L'>
#   = sqrt((2L+1)(2L'+1)(2K+1))
#     * { l_a l_b L;
#         l_c l_d L';
#         k1 k2 K  }   (9j symbol)
#     * <l_a||A^(k1)||l_c> <l_b||B^(k2)||l_d>
#
# For the HeI (1s)(2p) ^3P_J, we have L = L' = 1 (spatial triplet).
#
# Direct path: (l_a, l_b) = (0, 1), (l_c, l_d) = (0, 1), L = L' = 1.
# Exchange path: (l_a, l_b) = (0, 1), (l_c, l_d) = (1, 0), L = L' = 1.
#
# The DIFFERENCE between direct and exchange in the He triplet is the
# antisymmetry sign. For the TRIPLET (spatial antisymmetric) we have
# |triplet spatial> = (1/sqrt 2) [|1s ⊗ 2p> - |2p ⊗ 1s>] so the
# matrix element is (direct - exchange) with each evaluated on the same
# L=1 coupled basis.

def spatial_reduced_me_bipolar(la, lb, lc, ld, L, Lp, k1, k2, K):
    """<(l_a l_b) L || [C^(k1)(1) (x) C^(k2)(2)]^(K) || (l_c l_d) L'>
    via the 9j identity.
    """
    ninej_val = wigner_9j(la, lb, L, lc, ld, Lp, k1, k2, K)
    if simplify(ninej_val) == 0:
        return Integer(0)
    return (sqrt(Integer((2 * L + 1) * (2 * Lp + 1) * (2 * K + 1))) *
            ninej_val *
            reduced_C(la, lc, k1) * reduced_C(lb, ld, k2))


# ======================================================================
# Part 2: The bipolar radial kernel
# ======================================================================
#
# For the Drake formula we need to relate the bipolar-expansion radial
# integrals to the Slater-like M^K_{dir, exch}.
#
# BEST KNOWN DERIVATION (Drake 1971, Eq. 16 supplemented by Bethe-Salpeter
# §39):
#
# The Breit-Pauli spin-spin tensor operator (after angular projection to
# rank-2 space) gives the radial matrix element in the form
#
#   M^K_{a b; c d} = <n_a l_a, n_c l_c; n_b l_b, n_d l_d |
#                      r_<^K / r_>^{K+3} |
#                      same> radial
#
# This IS what `breit_ss_radial` computes. So the M^K_dir, M^K_exch in
# Drake's formula are just these radial integrals, for the direct
# (a=c, b=d) and exchange (a=d, b=c) channel.
#
# The bipolar-expansion insight: the SAME radial kernel r_<^K/r_>^{K+3}
# applies for BOTH direct and exchange — the issue is not the radial
# kernel but the ANGULAR coefficient that multiplies it.
#
# This simplifies our derivation: the radial factor is "frozen" as
# M^K_{dir/exch}, and we need to DERIVE the angular coefficient in front
# of each.

# ======================================================================
# Part 3: Full matrix element in the 9j representation
# ======================================================================


def drake_matrix_element_bipolar(K, L=1, S=1):
    """Compute <(l_a l_b) L, (s_a s_b) S; J M_J | V^{(0)} | same>

    where V^{(0)} = [T^(K)(space) · U^(K)(spin)]^(0) is a scalar tensor
    with rank-K spatial piece T^(K) (built bipolar) and rank-K spin piece
    U^(K) = [s_1 (x) s_2]^(K).

    For He (1s)(2p):
      l_a = 0 (1s), l_b = 1 (2p)
      s_a = s_b = 1/2

    For triplet S=1:  |triplet spatial> = (1/sqrt 2) [|l_a l_b> - |l_b l_a>]
                      coupled to L=1

    Direct term: <(0, 1) 1 L || T^(K) || (0, 1) 1 L> = sum over (k1, k2)
    Exchange term: <(0, 1) 1 L || T^(K) || (1, 0) 1 L'> (L'=L=1)

    The SPATIAL reduced matrix element has bipolar expansion coefficient
    A(k1, k2, K) = sqrt((2k1+1)(2k2+1)/(4pi(2K+1))) <k1 0 k2 0 | K 0>
    times radial kernel F(k1, k2, K; r1, r2).

    For the Breit-Pauli kernel (Y^K/r_12^{K+1}), the radial kernel F is
    the SAME for all (k1, k2): F = r_<^K / r_>^{K+3} (up to a universal
    normalization 4pi / (2K+1) that depends only on K, NOT on k1, k2).

    Returns a symbolic expression in M_sym[(K, 'dir')] and M_sym[(K, 'exch')].
    """
    la, lb = 0, 1

    # --- DIRECT channel: (l_a l_b) = (0, 1) -> (l_c l_d) = (0, 1) ---
    # Allowed (k1, k2): k1 = 0 (forced by <l=0 || C^(k1) || l=0>),
    #                    k2 triangle with K.
    direct_sum = Integer(0)
    for k1 in range(0, 3):
        for k2 in range(0, 5):
            if abs(k1 - k2) > K or k1 + k2 < K:
                continue
            # Parity check: C^(k) has parity (-1)^k; need k1+k2+K even for
            # the bipolar coupling to be parity-allowed.
            # Check via <k1 0 k2 0 | K 0> which is zero unless k1+k2+K is even.
            cg = clebsch_gordan(k1, k2, K, 0, 0, 0)
            if simplify(cg) == 0:
                continue
            # Reduced matrix elements:
            red1 = reduced_C(la, 0, k1)  # <l=0 || C^(k1) || l=0> = delta_{k1,0}
            red2 = reduced_C(lb, 1, k2)  # <l=1 || C^(k2) || l=1> = k2 even, triangle
            if simplify(red1) == 0 or simplify(red2) == 0:
                continue
            # 9j coupling:
            n9 = wigner_9j(la, lb, 1, la, lb, 1, k1, k2, K)
            if simplify(n9) == 0:
                continue
            # Angular coefficient from the bipolar expansion of Y^K/r_12^{K+1}.
            # The bipolar form (Varshalovich 5.17.4) in Racah-normalized form:
            #   Y^K_Q(hat r_12) / r_12^{K+1} = sum_{k1,k2} alpha(k1,k2,K)
            #                                    * [C^(k1)(1) (x) C^(k2)(2)]^K_Q
            #                                    * (r_<^K / r_>^{K+3})
            # with alpha(k1, k2, K) = specific rational coefficient
            # containing <k1 0 k2 0 | K 0> and normalization factors.
            # We compute alpha(k1, k2, K) below; for now we leave it symbolic
            # as A_k1k2K and extract the mixing ratio.
            A_k1k2K = sqrt(Integer((2 * k1 + 1) * (2 * k2 + 1))) * cg  # placeholder

            contrib = (sqrt(Integer((2 * 1 + 1) * (2 * 1 + 1) * (2 * K + 1)))
                       * n9 * red1 * red2 * A_k1k2K)
            direct_sum = direct_sum + contrib

    # --- EXCHANGE channel: (l_a l_b) = (0, 1) -> (l_c l_d) = (1, 0) ---
    exch_sum = Integer(0)
    for k1 in range(0, 3):
        for k2 in range(0, 5):
            if abs(k1 - k2) > K or k1 + k2 < K:
                continue
            cg = clebsch_gordan(k1, k2, K, 0, 0, 0)
            if simplify(cg) == 0:
                continue
            red1 = reduced_C(la, 1, k1)  # <l=0 || C^(k1) || l=1>: k1=1
            red2 = reduced_C(lb, 0, k2)  # <l=1 || C^(k2) || l=0>: k2=1
            if simplify(red1) == 0 or simplify(red2) == 0:
                continue
            n9 = wigner_9j(la, lb, 1, 1, 0, 1, k1, k2, K)
            if simplify(n9) == 0:
                continue
            A_k1k2K = sqrt(Integer((2 * k1 + 1) * (2 * k2 + 1))) * cg

            contrib = (sqrt(Integer((2 * 1 + 1) * (2 * 1 + 1) * (2 * K + 1)))
                       * n9 * red1 * red2 * A_k1k2K)
            exch_sum = exch_sum + contrib

    return direct_sum, exch_sum


def main():
    print("=" * 76)
    print("Sprint 5 Track DV: Drake bipolar expansion derivation")
    print("=" * 76)
    print()
    print("Allowed (k1, k2) channels for He (1s)(2p) ^3P_J:")
    print()

    for K, label in [(2, "SS"), (1, "SOO")]:
        print(f"--- K = {K} ({label}) ---")

        # Enumerate which (k1, k2) contribute to direct and exchange
        print("  Direct (l_a, l_b, l_c, l_d) = (0, 1, 0, 1):")
        la, lb, lc, ld = 0, 1, 0, 1
        for k1 in range(0, 3):
            for k2 in range(0, 5):
                if abs(k1 - k2) > K or k1 + k2 < K:
                    continue
                red1 = reduced_C(la, lc, k1)
                red2 = reduced_C(lb, ld, k2)
                cg = clebsch_gordan(k1, k2, K, 0, 0, 0)
                if simplify(red1) == 0 or simplify(red2) == 0 or simplify(cg) == 0:
                    continue
                n9 = wigner_9j(la, lb, 1, lc, ld, 1, k1, k2, K)
                if simplify(n9) == 0:
                    continue
                print(f"    (k1={k1}, k2={k2}): red1={simplify(red1)}, red2={simplify(red2)}, "
                      f"9j={simplify(n9)}, CG(k1 0 k2 0|K 0)={simplify(cg)}")

        print("  Exchange (l_a, l_b, l_c, l_d) = (0, 1, 1, 0):")
        lc, ld = 1, 0
        for k1 in range(0, 3):
            for k2 in range(0, 5):
                if abs(k1 - k2) > K or k1 + k2 < K:
                    continue
                red1 = reduced_C(la, lc, k1)
                red2 = reduced_C(lb, ld, k2)
                cg = clebsch_gordan(k1, k2, K, 0, 0, 0)
                if simplify(red1) == 0 or simplify(red2) == 0 or simplify(cg) == 0:
                    continue
                n9 = wigner_9j(la, lb, 1, lc, ld, 1, k1, k2, K)
                if simplify(n9) == 0:
                    continue
                print(f"    (k1={k1}, k2={k2}): red1={simplify(red1)}, red2={simplify(red2)}, "
                      f"9j={simplify(n9)}, CG(k1 0 k2 0|K 0)={simplify(cg)}")

        print()


if __name__ == "__main__":
    main()
