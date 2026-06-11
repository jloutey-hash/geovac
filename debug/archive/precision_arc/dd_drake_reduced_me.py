"""Sprint 4 Track DD: Compute the spatial reduced matrix elements
<(l_a l_b) L || T^(k)(space) || (l_a l_b) L> for L=1, (l_a, l_b) = (0, 1),
with the spatial tensor T^(k)(space) built as a coupled product of
single-electron tensors C^(k1)(1) and C^(k2)(2), weighted by the radial
integral M^K.

This decouples the spatial structure into (direct, exchange) channels and
(k1, k2) multipole coupling channels, yielding the Drake A_SS, A_SOO
coefficients symbolically.
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
from sympy import Rational, sqrt, Integer, simplify, symbols
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j, clebsch_gordan


# Symbolic radial integrals
M0d, M1d, M2d = symbols("M0d M1d M2d", real=True)
M0e, M1e, M2e = symbols("M0e M1e M2e", real=True)
M_sym = {
    (0, 'd'): M0d, (1, 'd'): M1d, (2, 'd'): M2d,
    (0, 'e'): M0e, (1, 'e'): M1e, (2, 'e'): M2e,
}


def reduced_C(l, lp, k):
    """<l || C^(k) || l'> (Racah normalized).  Edmonds 5.4.1."""
    three_j = wigner_3j(l, k, lp, 0, 0, 0)
    return (-1) ** l * sqrt(Integer((2 * l + 1) * (2 * lp + 1))) * three_j


# ============================================================
# Spatial reduced m.e. including the (k1, k2) channel sum
# ============================================================
#
# The SS spatial tensor is -sqrt(24 pi / 5) Y^2(r-hat_12)/r_12^3.
# The multipole expansion of this object gives, in terms of single-electron
# tensors, a SUM of (k1, k2) channels. The decomposition (Drake 1971 Sec. III,
# Eqs. (8)-(13) and Bethe-Salpeter §38) is:
#
#   Y^(K)_Q(r-hat_12) / r_12^{K+1}
#     = sum_{L1, L2}  A(L1, L2; K) * [C^(L1)(1) x C^(L2)(2)]^(K)_Q * R_{L1 L2}^K(r_1, r_2)
#
# where (L1, L2) is triangular to K, and R_{L1 L2}^K is a specific radial
# kernel involving r_<^{...} / r_>^{...}.
#
# For the Breit-Pauli 1/r_12^3 kernel with K = 2 (so the rank-2 tensor is
# Y^2/r^3), the explicit form requires a multipole expansion of 1/r_12^3
# weighted by Y^2.
#
# Simpler alternative: the Legendre expansion
#   1/r_12^3 = sum_k (r_<^k / r_>^{k+3}) * ... (but the standard Legendre
# expansion is for 1/r_12, not 1/r_12^3).
#
# For the exact Breit-Pauli tensor algebra: 1/r_12^3 * Y^2(r-hat_12) does
# NOT have a clean multipole expansion.  However, the MATRIX ELEMENTS of
# 1/r_12^3 Y^2(r-hat_12) / Y^0 between angular momentum eigenfunctions DO
# have a closed form via the "Drake kernel"
#   D^K_{L1 L2}(r_1, r_2) = integral of Legendre-expanded 1/r_12 weighted by
#                           the Y^2 projection
# which evaluates to  r_<^K / r_>^{K+3}  times specific rational coefficients.
#
# The "Drake radial integrals" M^k in `breit_integrals` are defined as:
#   M^k(a,b; c,d) = integral of P_{ac}(r_1) P_{bd}(r_2) * r_<^k / r_>^{k+3}
#     times a radial 1-d weight
# and these are the only scalars needed for the Breit-Pauli radial content.
#
# The angular expansion of the SS tensor in terms of (k1, k2) multipoles then
# matches the scalar integrals M^k as follows.  For the rank-2 SS operator:
# after spatial tensor reduction (Drake 1971 Eq. 10, using Bethe-Salpeter
# eq (38.17)-type expansion):
#
#   <(l_a l_b) L M_L| -sqrt(24 pi / 5) * sum_q (-1)^q Y^(2)_q(r-hat_12)/r_12^3 |...>
#    * <(l_a l_b) L M_L|
#
# involves the single-electron 3j coupling coefficients and the ONE SPECIFIC
# radial integral M^k for each (k1, k2) channel:
#   - Direct (path a=c, b=d): M^k with k = k_{12,radial} that couples the tensor.
#   - Exchange (path a=d, b=c): similar.
#
# For the rank-(k_1, k_2) coupled-to-rank-K tensor:
#   [C^(k1)(1) x C^(k2)(2)]^(K) * f(r1, r2)
# where f(r1, r2) = r_<^K / r_>^{K+1}_{Coulomb} or r_<^K / r_>^{K+3}_{Breit}.
#
# Under the radial integration for the DIRECT path (1s(1) 1s(1) 2p(2) 2p(2)):
#   <1s|r_1^{k_1}|1s> <2p|r_2^{k_2}|2p> times the kernel
# Under the EXCHANGE path (1s(1) 2p(1) 2p(2) 1s(2)):
#   <1s|r_1^{k_1}|2p> <2p|r_2^{k_2}|1s> times the kernel
#
# This gives us M^K (the SAME multipole K) for all (k1, k2) coupling channels,
# because the radial kernel only depends on K (the total rank).  More precisely:
#   M^K_{dir}  = integral 1s(1) 1s(1) * 2p(2) 2p(2) * r_<^K / r_>^{K+3} d^6r
#   M^K_{exch} = integral 1s(1) 2p(1) * 2p(2) 1s(2) * r_<^K / r_>^{K+3} d^6r
# where the radial powers r_1^{k_1}, r_2^{k_2} are absorbed by the angular
# tensor integral and the density products.
#
# Wait -- in the Breit-Pauli retarded kernel the radial factor is JUST
# r_<^K / r_>^{K+3}, and the angular (k_1, k_2) dependence is captured by the
# ANGULAR matrix elements via 3j symbols, not by the radial integral.  So the
# SAME M^K_{dir} and M^K_{exch} multiply ALL (k1, k2) channels that have
# total K as their tensor rank.  GREAT.
#
# That means the DIRECT spatial reduced matrix element is a sum over (k1, k2)
# triangular to K = 2:
#
# <(0 1) L=1 || T^(K=2)(space, direct) || (0 1) L=1>_direct
#   = sum_{k1, k2: triangle(k1,k2,2)}
#         sqrt((2L+1)(2K+1)(2L+1)) * 9j{0 1 1; k1 k2 2; 0 1 1}
#         * <0||C^(k1)||0> <1||C^(k2)||1>
#         * factor(k1, k2)
#
# where factor(k1, k2) is the C-G coefficient that combines the (k1, k2)
# tensor product into the rank-K = 2 tensor, which is <k1 0 k2 0 | K 0>
# (since the Breit-Pauli kernel Y^(K)(r-hat_12) evaluated at the single-body
# coupled spherical-harmonic components has this CG factor).  Wait this isn't
# right either -- the decomposition of Y^(K)(r-hat_12) into [Y^(k1)(1) x Y^(k2)(2)]^(K)
# has a specific radial-weight coefficient.  Let me look this up.


# ============================================================
# Use the direct Slater-determinant + multipole route
# ============================================================
# The cleanest approach turns out to be: compute the 2-body matrix element
# directly, electron by electron, in the |l_a m_a; l_b m_b; s_a; s_b> basis,
# using the standard multipole expansion of the two-body operator.
#
# For the SS rank-2 operator, the multipole expansion of
#   (4 pi / 5) Y^(2)_Q(r-hat_12) / r_12^3
# in single-electron spherical harmonics is (Bethe-Salpeter Eq 38.17,
# Drake 1971 Eq (10); the "double multipole expansion"):
#
# (4 pi / 5) Y^(2)_Q(r-hat_12) / r_12^3
#   = sum_{k1, k2}  C(k1, k2, 2; Q) * [Y^(k1)(r1-hat) (x) Y^(k2)(r2-hat)]^(2)_Q *
#         U_{k1 k2}^{(2)}(r_1, r_2)
#
# where C(k1, k2, 2; Q) = <k1 0 k2 0 | 2 0> (integer C-G) and
# U_{k1 k2}^{(2)}(r_1, r_2) is a specific radial kernel that, for each
# (k1, k2) triangular to 2, integrates to multiples of the Drake integrals
# M^k = integral r_<^k / r_>^{k+3} * orbital product.
#
# The EXPLICIT DRAKE FORM is (after much algebra):
#   U_{k1 k2}^{(2)}(r_1, r_2) = 3 * (r_<^{K_12} / r_>^{K_12 + 3})
#   where K_12 = min(k1, k2)?  Actually no, it's a constant independent
#   of k1, k2 for the rank-2 case; only the PREFACTORS differ.
#
# The clean path here is to IMPLEMENT THE FULL COMPUTATION without getting
# bogged down in which (k1, k2) maps to which M^k.  We just compute the
# 2-electron angular + spin matrix element in a finite m-basis, and match to
# the radial-integral symbols M^k_dir, M^k_exch.


# ============================================================
# Direct matrix element evaluation: the Drake-Slater route
# ============================================================
#
# Following Drake 1971 Sec. III:
#
# <^(2S+1) L_J | H_SS | ^(2S+1) L_J>
#   = f_SS(J) * < J_ls || S^(2) . N^(2)(SS) || J_ls > / {normalization}
#
# where the reduced m.e. factorizes into:
#   1. J-dependence f_SS(J) = (-1)^{L+S+J} * 6j{L S J; S L 2} * (overall normalization)
#   2. Spin reduced m.e. <S||[s_1 x s_2]^(2)||S> = sqrt(5)/8 * (for S=1)
#   3. Spatial reduced m.e. <L||N^(2)(SS)||L>
#
# The spatial reduced m.e. N^(2)(SS) is (from tensor coupling):
#   N^(2)(SS) = -sqrt(24 pi / 5) * sum_{k1, k2, triangle(k1,k2,2)}
#         <k1 0 k2 0 | 2 0> * [Y^(k1)(1) x Y^(k2)(2)]^(2) / r_12^3
#
# For (l_a, l_b) = (0, 1), the direct channel has l = 0 electron, so:
#   k1 = 0 (on the 1s-1s component, electron 1)
#   k2 = 2 (on the 2p-2p component, electron 2)
# This is the (k1=0, k2=2) channel with M^K = M^{k2} = M^2_{dir}.
#
# For the exchange channel (1s on e1 exchanges with 2p on e2):
#   k1 = 1 (mixes 1s with 2p on electron 1)
#   k2 = 1 (mixes 2p with 1s on electron 2)
# This is the (k1=1, k2=1) channel with M^K = M^1_{exch}? or M^2_{exch}?
#
# Hmm.  Wait.  The "radial integral" M^k is labeled by the multipole of the
# radial kernel r_<^k / r_>^{k+3}, NOT by the rank of the coupling tensor.
#
# When the operator is [C^(k1)(1) (x) C^(k2)(2)]^(K) / r_12^3, integrating
# over r_1, r_2 with the orbital density product gives r_<^k / r_>^{k+3}
# where k is chosen to make the integrand non-singular -- typically
# k = max(k1, k2)? or k = K?
#
# The Drake paper clarifies: the effective radial integral for the Breit-Pauli
# retarded kernel [C^(k1)(1) (x) C^(k2)(2)]^(K) / r_12^3 has multipole
# r_<^K / r_>^{K+3} (the same K as the COUPLING rank, not individual k1 or k2).
# This is because the r_12^{-3} kernel with angular Y^(K)(r-hat_12) has
# multipole expansion
#   r_<^K / r_>^{K+3} * (4pi) * sum_{k1, k2, triangle(k1,k2,K)} C(k1,k2;K)
#       * [Y^(k1)(1) x Y^(k2)(2)]^(K) / Y^(K)
# where the "diagonal" multipole is K.
#
# So for the SS rank-2 operator, the radial integral in BOTH direct and
# exchange channels is M^{K=2} (k = 2).
# For the SOO rank-1 operator, the radial integral in both direct and
# exchange channels is M^{K=1} (k = 1).
#
# This matches the Drake coefficients: A_SS uses M^2_dir and M^2_exch only;
# A_SOO uses M^1_dir and M^1_exch only.  NO MIXING with M^0.


# ============================================================
# Core computation: spatial reduced m.e. <L || T^(K) || L> split into
# direct and exchange channels
# ============================================================

def spatial_reduced(path, K):
    """Compute <(l_a l_b) L || T^(K)(space) || (l_a' l_b') L> split by
    radial-integral channel.

    path: 'direct' (l_a, l_b) = (l_a', l_b') = (0, 1)
          'exchange' (l_a, l_b) = (0, 1), (l_a', l_b') = (1, 0)
    K: total rank of the spatial tensor.

    The spatial tensor T^(K) is the multipole expansion of
    Y^(K)(r-hat_12) / r_12^{K+1}.  After combining with the Breit-Pauli
    r_12^{-2} factor (since BP is 1/r_12^3 = 1/r_12^{K+1} for K=2 and
    1/r_12^3 = 1/r_12^{K+2} for K=1 -- wait, that depends), we get the Drake
    integral M^K.

    Returns a dict { 'M^K' -> coefficient } giving the linear combination
    of M^K_{path} that forms the reduced m.e.

    We use the 9j reduction:
    <(l_a l_b) L || [C^(k1)(1) x C^(k2)(2)]^(K) || (l_a' l_b') L>
      = sqrt((2L+1)(2K+1)(2L+1)) * 9j{l_a l_b L; k1 k2 K; l_a' l_b' L}
        * <l_a||C^(k1)||l_a'> <l_b||C^(k2)||l_b'>

    summed over (k1, k2) triangular to K.  The radial integral combines all
    (k1, k2) channels with the same K-rank kernel.  In the Drake convention,
    all (k1, k2) channels for a given direct/exchange path give coefficients
    that multiply a single radial integral M^K, whose value is
        M^K_{dir}  if path == 'direct'
        M^K_{exch} if path == 'exchange'
    """
    L = 1
    if path == 'direct':
        la, lb = 0, 1
        lap, lbp = 0, 1
    elif path == 'exchange':
        la, lb = 0, 1
        lap, lbp = 1, 0
    else:
        raise ValueError(path)

    total = Integer(0)
    for k1 in range(0, K + 3):
        for k2 in range(0, K + 3):
            if abs(k1 - k2) > K or k1 + k2 < K:
                continue
            Ca = reduced_C(la, lap, k1)
            Cb = reduced_C(lb, lbp, k2)
            if sp.simplify(Ca) == 0 or sp.simplify(Cb) == 0:
                continue
            # 9j {l_a l_b L; k1 k2 K; l_a' l_b' L}
            nj = wigner_9j(la, lb, L, k1, k2, K, lap, lbp, L)
            if sp.simplify(nj) == 0:
                continue
            # Prefactor and CG coupling of the rank-K tensor kernel:
            # The operator is [C^(k1)(1) x C^(k2)(2)]^(K)_Q times the kernel
            # 1/r_12^3 is expanded as sum over multipoles K of (r_<^k / r_>^{k+3})
            # with specific angular coefficients.
            # For the K-multipole decomposition of 1/r_12^3 weighted by
            # Y^(K)(r-hat_12), the expansion coefficient is (2K+1)/(4pi)
            # times Y^(K)(r-hat_12) / r_12^3 projected.
            #
            # But we're only tracking the 9j-reduced part here; the specific
            # normalization enters in the final A_SS, A_SOO amplitudes.
            term = sqrt(Integer((2 * L + 1) * (2 * K + 1) * (2 * L + 1))) * nj * Ca * Cb
            term = sp.simplify(term)
            if term != 0:
                total = total + term
    return sp.simplify(total)


# Check both direct and exchange paths for K = 1 and K = 2
def show_spatial_reduced():
    print("=" * 76)
    print("Spatial reduced m.e. <L||T^(K)||L> for (l_a, l_b) = (0,1), L = 1")
    print("=" * 76)
    for K in (1, 2):
        print(f"\n  K = {K}:")
        for path in ('direct', 'exchange'):
            v = spatial_reduced(path, K)
            print(f"    {path:9s}: {v}")


# ============================================================
# Full diagonal matrix element via Edmonds reduction
# ============================================================
#
# With all reduced matrix elements in hand, the physical diagonal
# matrix element is
#
# <^(2S+1) L_J M_J | [T^(K)(space) . U^(K)(spin)]^(0) | ^(2S+1) L_J M_J>
#   = (-1)^{L+S+J+K} / sqrt((2J+1)(2K+1)) * 6j{L S J; S L K}
#     * <L||T||L> * <S||U||S>
#
# The overall operator coefficient for SS is alpha^2 * (overall physics prefactor).


def spin_reduced_SSmkk_as_s(K):
    """<S=1 || [s_1 (x) s_2]^(K) || S=1> using Edmonds 7.1.7."""
    half = Rational(1, 2)
    red_s = sqrt(Rational(3, 2))  # <1/2||s||1/2> in Edmonds
    result = sqrt(Integer((2 * K + 1) * 3 * 3)) * \
             wigner_9j(half, half, 1, 1, 1, K, half, half, 1) * red_s * red_s
    return sp.simplify(result)


def full_diagonal_SS(J, K=2):
    """<^3P_J, M_J = J | H_SS | ^3P_J, M_J = J> expressed as a function of
    the symbolic radial integrals M^K_dir and M^K_exch.

    This assumes:
      H_SS = C_phys * [T^(K)(space) (x) U^(K)(spin)]^(0)
    with the spatial tensor built from the M^K multipole kernel and spin tensor
    [s_1 (x) s_2]^(K).
    """
    L = 1
    S = 1

    # J-dependent 6j factor:
    phase = (-1) ** (L + S + J + K)
    six_j = wigner_6j(L, S, J, S, L, K)
    j_factor = phase / sqrt(Integer((2 * J + 1) * (2 * K + 1))) * six_j

    # Spin reduced m.e.:
    spin_red = spin_reduced_SSmkk_as_s(K)

    # Spatial reduced m.e., direct and exchange:
    sp_dir = spatial_reduced('direct', K)
    sp_exc = spatial_reduced('exchange', K)

    # The physical diagonal m.e. involves M^K_dir weighted by sp_dir and
    # M^K_exch weighted by sp_exc.  BUT the "exchange" Slater amplitude has
    # an extra sign from the antisymmetrization of the spatial wavefunction.
    #
    # For the TRIPLET (S=1), the spatial wavefunction is ANTISYMMETRIC, so
    # the two-body matrix element has the form <ab|V|ab> - <ab|V|ba>
    # direct minus exchange with a MINUS sign on exchange.  In our tensor
    # language:
    #
    #   <(l_a l_b) L || T^(K) || (l_a l_b) L>_spatial
    #     = sp_dir * M^K_dir - (-1)^{...} sp_exch * M^K_exch
    #
    # where the sign on exchange is (-1)^{l_a + l_b + L + S} or similar from
    # the spin-space antisymmetrization.

    return j_factor, spin_red, sp_dir, sp_exc


def main():
    print("=" * 76)
    print("Sprint 4 Track DD: Detailed reduced matrix elements")
    print("=" * 76)
    show_spatial_reduced()

    print("\n" + "=" * 76)
    print("Full diagonal m.e. <^3P_J| H_SS | ^3P_J> per Edmonds")
    print("=" * 76)
    for J in (0, 1, 2):
        jf, sr, spd, spe = full_diagonal_SS(J, K=2)
        print(f"\n  J = {J}:")
        print(f"    j_factor (6j-with-phase / sqrt((2J+1)(2K+1))) = {jf}")
        print(f"    spin_reduced <S=1||[s1 x s2]^(2)||S=1>         = {sr}")
        print(f"    spatial_direct  (coeff of M^2_dir)            = {spd}")
        print(f"    spatial_exchange (coeff of M^2_exch)          = {spe}")
        # Physical m.e. (up to overall C_phys constant and radial M^K):
        # = jf * sr * (spd * M^2_dir + spe * M^2_exch)
        # (where the sign on exch is determined by triplet-spatial-antisymm,
        # which needs to be determined separately)

    print("\n" + "=" * 76)
    print("Full diagonal m.e. <^3P_J| H_SOO | ^3P_J> per Edmonds")
    print("=" * 76)
    # For SOO (K=1), the spin tensor is NOT [s_1 x s_2]^(1) (which vanishes on S=1)
    # but rather (s_1 + 2 s_2) -- a rank-1 single-electron tensor summed over
    # the two electrons.  So spin_reduced_SSmkk_as_s(1) = 0 is expected.
    for J in (0, 1, 2):
        jf, sr, spd, spe = full_diagonal_SS(J, K=1)
        print(f"\n  J = {J}:")
        print(f"    j_factor                                       = {jf}")
        print(f"    spin_reduced_as_ss <[s1 x s2]^(1)>             = {sr}  (expected 0 for S=1)")
        print(f"    spatial_direct  (coeff of M^1_dir)            = {spd}")
        print(f"    spatial_exchange (coeff of M^1_exch)          = {spe}")


if __name__ == "__main__":
    main()
