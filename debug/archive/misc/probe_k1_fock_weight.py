"""
Probe whether kappa = -1/16 emerges from the Fock momentum-space weight function.

The Fock projection maps hydrogen wavefunctions from R^3 to S^3 via
stereographic projection. In momentum space, the hydrogenic wavefunction
phi_{nlm}(p) is related to hyperspherical harmonics Y_{nlm}(Omega) on S^3 by:

    phi_{nlm}(p) = (2 p0)^{5/2} / (p^2 + p0^2)^2  *  Y_{nlm}(Omega(p))

where Omega(p) is the stereographic image of p on S^3. The factor
w(p) = (p^2 + p0^2)^{-2} is the Fock weight.

The GeoVac graph Hamiltonian has:
  - Diagonal: -Z^2/(2n^2)  (exact hydrogenic eigenvalue)
  - Off-diagonal: kappa * (-A) = (1/16) * A  where A connects (n,l,m)↔(n±1,l,m)

QUESTION: Do the off-diagonal Fock weight matrix elements between adjacent
shells equal 1/16 or a rational multiple thereof?

APPROACH:
  (A) Compute <n+1,l,m| w(p) |n,l,m> in momentum space using the known
      hydrogenic momentum wavefunctions.
  (B) Compute the Laplace-Beltrami operator matrix elements in the
      Fock-weighted inner product.
  (C) Compute matrix elements of (p^2+p0^2)^{-2} using S^3 hyperspherical
      harmonic recursion relations.

All computations use sympy for exact rational/symbolic arithmetic.
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

import sympy as sp
from sympy import (
    Rational, Symbol, symbols, sqrt, pi, oo, integrate,
    simplify, factorial, assoc_laguerre, exp, cos, sin,
    gegenbauer, binomial, cancel, together, apart, factor
)

# ============================================================================
# Momentum-space hydrogenic wavefunctions
# ============================================================================
#
# The normalized hydrogenic momentum wavefunction is (Bethe-Salpeter, Podolsky-Pauling):
#
#   phi_{nlm}(p) = Y_{lm}(theta_p, phi_p) * f_{nl}(p)
#
# where the radial part in momentum space is:
#
#   f_{nl}(p) = N_{nl} * (2 p0)^{5/2} * n^2 * (2 n p / p0)^l / (n^2 p^2/p0^2 + 1)^{l+2}
#             * C_{n-l-1}^{l+1}((n^2 p^2/p0^2 - 1)/(n^2 p^2/p0^2 + 1))
#
# with p0 = Z/n for shell n, and C_k^lambda is the Gegenbauer polynomial.
#
# We work at Z=1 (hydrogen). The momentum-space integration measure is
# p^2 dp * dOmega, and the angular part gives orthogonality in (l,m).
#
# For the radial overlap, we use the substitution:
#   t = p / p0, so p = p0 * t, dp = p0 dt
#   argument of Gegenbauer: (n^2 t^2 - 1)/(n^2 t^2 + 1) = x
#
# At p0 = 1 (ground state shell) or at general p0:
#   phi_{nl}(p) = N_{nl} * (2p0)^{5/2} * n^2 * (2np/p0)^l / (n^2 p^2/p0^2 + 1)^{l+2}
#                * C_{n-l-1}^{l+1}((n^2 p^2/p0^2 - 1)/(n^2 p^2/p0^2 + 1))
#
# The normalization is int |phi_{nlm}|^2 d^3p = 1.

p_var = Symbol('p', positive=True, real=True)
p0_var = Symbol('p0', positive=True, real=True)
t_var = Symbol('t', positive=True, real=True)  # t = p / p0
x_var = Symbol('x', real=True)  # Gegenbauer argument

def hydrogenic_momentum_radial_unnorm(n: int, l: int, p0_val):
    """
    Unnormalized radial momentum wavefunction f_{nl}(p) at general p0.

    Uses the Fock representation:
      f_{nl}(p) ∝ p^l * p0^{l+5/2} / (p^2 + p0^2/n^2)^{l+2}
                 * C_{n-l-1}^{l+1}(x)

    where x = (p^2 - p0^2/n^2) / (p^2 + p0^2/n^2) (Fock's change of variable).

    At the energy shell, E_n = -Z^2/(2n^2), p0 = Z = 1.
    The characteristic momentum for shell n is p_n = Z/n = p0/n.

    The conventional form uses (p0/n) as the momentum scale parameter.
    """
    # p_n = p0/n is the characteristic momentum for shell n
    p_n = p0_val / n

    # Fock's substitution variable
    # x = (p^2 - p_n^2) / (p^2 + p_n^2)
    denom_x = p_var**2 + p_n**2
    x_arg = (p_var**2 - p_n**2) / denom_x

    # Gegenbauer polynomial
    C = gegenbauer(n - l - 1, l + 1, x_arg)

    # Full unnormalized form (Podolsky-Pauling / Fock):
    # f_nl(p) ∝ p^l / (p^2 + p_n^2)^(l+2) * C_{n-l-1}^{l+1}(x)
    # The overall p0-dependent prefactor is (p_n)^{l + 5/2} * n^2
    # but for computing ratios and checking kappa, the overall constant
    # cancels. We keep it for correct normalization.

    prefactor = n**2 * (Rational(2) * p_n)**(l + Rational(5, 2))
    f_nl = prefactor * p_var**l / denom_x**(l + 2) * C

    return f_nl


def compute_radial_norm_sq(f_nl):
    """Compute 4*pi * int_0^inf |f_nl|^2 p^2 dp (radial norm squared)."""
    integrand = f_nl**2 * p_var**2
    result = 4 * pi * integrate(integrand, (p_var, 0, oo))
    return simplify(result)


def compute_radial_overlap_with_weight(f_n1l, f_n2l, weight_power=2):
    """
    Compute 4*pi * int_0^inf f_{n1,l}(p) * w(p) * f_{n2,l}(p) * p^2 dp

    where w(p) = (p^2 + p0^2)^{-weight_power}  (Fock weight at weight_power=2).

    Note: the Fock weight uses p0 (the overall energy scale), NOT p_n = p0/n.
    """
    w = (p_var**2 + p0_var**2)**(-weight_power)
    integrand = f_n1l * w * f_n2l * p_var**2
    result = 4 * pi * integrate(integrand, (p_var, 0, oo))
    return simplify(result)


def compute_radial_overlap(f_n1l, f_n2l):
    """
    Compute 4*pi * int_0^inf f_{n1,l}(p) * f_{n2,l}(p) * p^2 dp
    (bare overlap without Fock weight).
    """
    integrand = f_n1l * f_n2l * p_var**2
    result = 4 * pi * integrate(integrand, (p_var, 0, oo))
    return simplify(result)


# ============================================================================
# APPROACH A: Direct off-diagonal Fock weight matrix elements
# ============================================================================

def approach_A_direct_fock_weight():
    """
    Compute <n+1, l, m | w(p) | n, l, m> for the Fock weight w(p) = (p^2 + p0^2)^{-2}.

    Since the angular parts Y_{lm} are orthonormal, the matrix element
    reduces to the radial integral:

    <n',l| w |n,l> = 4*pi * int f_{n'l}(p) * w(p) * f_{nl}(p) * p^2 dp

    We compute for (n,l) -> (n+1,l) transitions (same l, same m — the
    inter-shell transitions that define the GeoVac adjacency).
    """
    print("=" * 70)
    print("APPROACH A: Direct off-diagonal Fock weight matrix elements")
    print("=" * 70)

    results = {}

    # Transitions to compute: (n, l) -> (n+1, l)
    transitions = [
        # Shell 1 -> 2
        (1, 0, 2, 0),  # (1,0,0) -> (2,0,0)
        # Shell 2 -> 3
        (2, 0, 3, 0),  # (2,0,0) -> (3,0,0)
        (2, 1, 3, 1),  # (2,1,m) -> (3,1,m)
        # Shell 3 -> 4
        (3, 0, 4, 0),  # (3,0,0) -> (4,0,0)
        (3, 1, 4, 1),  # (3,1,m) -> (4,1,m)
        (3, 2, 4, 2),  # (3,2,m) -> (4,2,m)
    ]

    # Also compute diagonal elements for normalization reference
    diag_cases = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]

    print("\n--- Diagonal elements <n,l|w|n,l> ---")
    print(f"{'(n,l)':<10} {'<n,l|w|n,l> (p0 free)':<40} {'at p0=1':<15}")
    print("-" * 65)

    for n, l in diag_cases:
        f_nl = hydrogenic_momentum_radial_unnorm(n, l, p0_var)
        norm_sq = compute_radial_norm_sq(f_nl)

        diag_weighted = compute_radial_overlap_with_weight(f_nl, f_nl)
        diag_ratio = simplify(diag_weighted / norm_sq)
        diag_at_1 = simplify(diag_ratio.subs(p0_var, 1))

        label = f"({n},{l})"
        print(f"{label:<10} {str(diag_ratio):<40} {str(diag_at_1):<15}")

        results[f"diag_({n},{l})_ratio_p0free"] = str(diag_ratio)
        results[f"diag_({n},{l})_at_p0=1"] = str(diag_at_1)

    print("\n--- Off-diagonal elements <n',l|w|n,l> / sqrt(<n',l|w|n',l> * <n,l|w|n,l>) ---")
    print(f"{'(n,l)->(n\',l)':<18} {'raw overlap':<35} {'normalized':<20} {'at p0=1':<15}")
    print("-" * 88)

    for n1, l1, n2, l2 in transitions:
        assert l1 == l2, "Must have same l for inter-shell transition"
        l = l1

        f_n1l = hydrogenic_momentum_radial_unnorm(n1, l, p0_var)
        f_n2l = hydrogenic_momentum_radial_unnorm(n2, l, p0_var)

        # Normalization squares
        norm1_sq = compute_radial_norm_sq(f_n1l)
        norm2_sq = compute_radial_norm_sq(f_n2l)

        # Weighted overlap
        offdiag_weighted = compute_radial_overlap_with_weight(f_n1l, f_n2l)

        # Normalized: <n'|w|n> / sqrt(<n'|w|n'> * <n|w|n>) — but this uses the
        # un-weighted norms. Let me compute a cleaner quantity:
        #
        # <n'|w|n> / sqrt(<n'|n'> * <n|n>) = Fock-weighted coupling between
        # normalized states
        offdiag_norm = simplify(offdiag_weighted / sqrt(norm1_sq * norm2_sq))
        offdiag_norm_p0_1 = simplify(offdiag_norm.subs(p0_var, 1))

        label = f"({n1},{l})->({n2},{l})"
        raw_str = str(simplify(offdiag_weighted))
        if len(raw_str) > 32:
            raw_str = raw_str[:32] + "..."
        print(f"{label:<18} {raw_str:<35} {str(offdiag_norm):<20} {str(offdiag_norm_p0_1):<15}")

        results[f"offdiag_({n1},{l})->({n2},{l})_norm_p0free"] = str(offdiag_norm)
        results[f"offdiag_({n1},{l})->({n2},{l})_norm_at_p0=1"] = str(offdiag_norm_p0_1)

    return results


# ============================================================================
# APPROACH B: Laplace-Beltrami in Fock-weighted inner product
# ============================================================================

def approach_B_laplacian_fock_weighted():
    """
    The Laplace-Beltrami operator on S^3 has eigenvalues -(n^2-1) with
    eigenfunctions Y_{nlm}. In the standard L^2(S^3) inner product,
    the operator is diagonal.

    In the FOCK-WEIGHTED inner product, with weight w(chi) depending on
    the hyperspherical angle chi, the matrix elements are:

    <Y_{n'lm}| Delta_{S3} |Y_{nlm}>_w = ∫ Y_{n'lm}* Delta_S3 Y_{nlm} * w(chi) d^3Omega

    On S^3, Fock's projection gives chi_n = 2 arctan(1/n) for shell n.
    The weight function in the S^3 picture is w(chi) = sin^2(chi/2) / [cos^2(chi/2) + sin^2(chi/2)]^2
    which is actually a smooth function of chi.

    Let's work with Gegenbauer polynomials instead. On S^3 (parametrized by chi in [0,pi]),
    the hyperspherical harmonics factorize as:

        Y_{nlm}(chi, theta, phi) = R_n^l(chi) * Y_{lm}(theta, phi)

    where R_n^l(chi) = N_n^l * sin^l(chi) * C_{n-l-1}^{l+1}(cos chi)

    The standard inner product on S^3 is d^3Omega = sin^2(chi) * sin(theta) d chi d theta d phi.

    The Fock projection corresponds to a specific weighting of this inner product
    via the Jacobian of the stereographic map.

    Rather than computing this directly on S^3, let's use momentum space and
    compute the Laplacian in that representation.
    """
    print("\n" + "=" * 70)
    print("APPROACH B: Laplace-Beltrami in Fock-weighted inner product")
    print("=" * 70)

    # On S^3, using chi as the "radial" hyperspherical angle:
    # The hyperspherical harmonic is proportional to sin^l(chi) * C_{n-l-1}^{l+1}(cos chi)
    # The S^3 measure is sin^2(chi) dchi d^2Omega (d^2Omega gives 4pi from spherical harmonics)

    # Fock's stereographic projection from momentum p to chi:
    #   p/p0 = tan(chi/2)   so   chi = 2 arctan(p/p0)
    #   cos(chi) = (p0^2 - p^2)/(p0^2 + p^2)   [NOTE: this is -x_Fock]
    #   sin(chi) = 2 p p0 / (p0^2 + p^2)
    #
    # The Jacobian of the Fock map is:
    #   dOmega_{S^3} = (2p0/(p^2+p0^2))^3 * (2p0) * d^3p
    #                = (2p0)^4 / (p^2+p0^2)^3 * p^2 dp dOmega_p
    #
    # So the S^3 inner product in momentum space is:
    #   <f|g>_{S^3} = int f*(Omega(p)) g(Omega(p)) (2p0)^4/(p^2+p0^2)^3 * p^2 dp dOmega_p
    #
    # The hydrogenic wavefunction phi_{nlm}(p) is related to Y_{nlm}(Omega) by:
    #   phi_{nlm}(p) = [(2p0)^{5/2}/(p^2+p0^2)^2] * Y_{nlm}(Omega)
    #   => Y_{nlm}(Omega) = [(p^2+p0^2)^2/(2p0)^{5/2}] * phi_{nlm}(p)
    #
    # So:
    #   <Y_{n'}|Y_n>_{S^3} = int Y* Y (2p0)^4/(p^2+p0^2)^3 p^2 dp dOmega_p
    #                       = int [(p^2+p0^2)^2/(2p0)^{5/2}]^2 * phi* phi * (2p0)^4/(p^2+p0^2)^3 * p^2 dp dOmega_p
    #                       = int phi* phi * (p^2+p0^2)/(2p0) * p^2 dp dOmega_p
    #
    # This means the standard S^3 inner product in momentum space representation is:
    #   <Y_{n'lm}|Y_{nlm}>_{S^3} = delta_{nn'} = int phi_{n'l}* phi_{nl} * (p^2+p0^2)/(2p0) * 4pi p^2 dp
    #
    # Now, the Laplace-Beltrami eigenvalue equation is Delta Y_n = -(n^2-1) Y_n.
    # In the standard inner product:
    #   <Y_{n'}| Delta |Y_n>_{S^3} = -(n^2-1) delta_{nn'}
    #
    # We want the matrix elements in a MODIFIED inner product that uses the
    # Fock weight w = (p^2+p0^2)^{-2} additionally:
    #   <Y_{n'}| Delta |Y_n>_w = int Y_{n'}* (Delta Y_n) * w * dOmega_{S^3}
    #                          = -(n^2-1) * int Y_{n'}* Y_n * w * dOmega_{S^3}
    #
    # Since Delta Y_n = -(n^2-1) Y_n, the operator is still diagonal times the
    # modified overlap! So the question reduces to computing:
    #   <Y_{n'}| Y_n>_w = int phi_{n'l}* phi_{nl} * (p^2+p0^2)/(2p0) * w(p) * 4pi p^2 dp
    #
    # where w(p) = (p^2+p0^2)^{-2}.
    #
    # This gives:
    #   <Y_{n'}| Y_n>_w = (1/(2p0)) * 4pi * int phi_{n'l} phi_{nl} * (p^2+p0^2)^{-1} p^2 dp
    #
    # This is different from Approach A! In Approach A we computed:
    #   int phi_{n'l} phi_{nl} * (p^2+p0^2)^{-2} p^2 dp
    #
    # Here we have:
    #   int phi_{n'l} phi_{nl} * (p^2+p0^2)^{-1} p^2 dp
    # (one less power of the weight).

    results = {}

    # Compute the modified overlap <Y_{n'}|Y_n>_w for adjacent shells
    # This is (1/(2p0)) * 4pi * int f_{n'l} * f_{nl} * (p^2+p0^2)^{-1} * p^2 dp

    transitions = [
        (1, 0, 2, 0),
        (2, 0, 3, 0),
        (2, 1, 3, 1),
        (3, 0, 4, 0),
        (3, 1, 4, 1),
        (3, 2, 4, 2),
    ]

    print("\nComputing Fock-weighted S^3 overlaps <Y_{n'}|Y_n>_w")
    print("where w = (p^2+p0^2)^{-2} on S^3 translates to (p^2+p0^2)^{-1} in mom. space")
    print()

    for n1, l1, n2, l2 in transitions:
        l = l1
        f_n1l = hydrogenic_momentum_radial_unnorm(n1, l, p0_var)
        f_n2l = hydrogenic_momentum_radial_unnorm(n2, l, p0_var)

        # Standard overlap (for normalization)
        # <Y|Y>_{S^3} = int phi phi * (p^2+p0^2)/(2p0) * 4pi p^2 dp
        # For normalized phi, this should give delta_{nn'}.
        # Since our phi's are unnormalized, compute the diagonal overlaps too.

        # Modified weight: (p^2+p0^2)^{-1} instead of (p^2+p0^2)^{-2}
        offdiag_w1 = compute_radial_overlap_with_weight(f_n1l, f_n2l, weight_power=1)

        # Diagonal for normalization
        norm1_sq = compute_radial_norm_sq(f_n1l)
        norm2_sq = compute_radial_norm_sq(f_n2l)

        # Normalized coupling
        coupling = simplify(offdiag_w1 / sqrt(norm1_sq * norm2_sq))
        coupling_p0_1 = simplify(coupling.subs(p0_var, 1))

        label = f"({n1},{l})->({n2},{l})"
        print(f"  {label}: coupling = {coupling}")
        print(f"  {label}: at p0=1: {coupling_p0_1}")

        results[f"B_({n1},{l})->({n2},{l})_coupling_p0free"] = str(coupling)
        results[f"B_({n1},{l})->({n2},{l})_coupling_at_p0=1"] = str(coupling_p0_1)

    return results


# ============================================================================
# APPROACH C: Direct S^3 Gegenbauer computation
# ============================================================================

def approach_C_gegenbauer_on_S3():
    """
    Work directly on S^3 with hyperspherical harmonics.

    The "radial" part of Y_{nlm} on S^3 (the chi-dependent part) is:

        R_n^l(chi) = N_n^l * sin^l(chi) * C_{n-l-1}^{l+1}(cos chi)

    The S^3 measure is sin^2(chi) dchi.

    The STANDARD inner product gives:
        int_0^pi R_{n'}^l(chi) R_n^l(chi) sin^2(chi) dchi = delta_{nn'} * (normalization)

    The Fock projection maps chi_n = 2 arctan(p/p0) at the energy shell p = p0/n.
    In the chi variable, the Fock weight w(p) = (p^2+p0^2)^{-2} transforms to:

        p = p0 tan(chi/2)
        p^2 + p0^2 = p0^2 / cos^2(chi/2)
        w(p) = cos^4(chi/2) / p0^4

    So the Fock weight on S^3 is simply cos^4(chi/2) (up to p0 factors).

    The off-diagonal coupling is then:
        <n',l| w |n,l>_{S^3} = N_n'^l * N_n^l * int_0^pi sin^{2l}(chi) C_{n'-l-1}^{l+1}(cos chi)
                                * C_{n-l-1}^{l+1}(cos chi) * cos^4(chi/2) * sin^2(chi) dchi

    Let u = cos(chi), then sin(chi) dchi = -du, sin^2(chi) = 1-u^2,
    cos^2(chi/2) = (1+u)/2, cos^4(chi/2) = (1+u)^2/4.

    So the integral becomes:
        (1/4) * int_{-1}^{1} (1-u^2)^l * C_{n'-l-1}^{l+1}(u) * C_{n-l-1}^{l+1}(u)
                              * (1+u)^2 * (1-u^2) du
      = (1/4) * int_{-1}^{1} (1-u)^{l+1} (1+u)^{l+3} C_{n'-l-1}^{l+1}(u) C_{n-l-1}^{l+1}(u) du
    """
    print("\n" + "=" * 70)
    print("APPROACH C: Direct S^3 Gegenbauer integration")
    print("=" * 70)

    u = Symbol('u', real=True)
    results = {}

    transitions = [
        (1, 0, 2, 0),
        (2, 0, 3, 0),
        (2, 1, 3, 1),
        (3, 0, 4, 0),
        (3, 1, 4, 1),
        (3, 2, 4, 2),
    ]

    # Gegenbauer normalization on S^3:
    # int_{-1}^1 (1-u^2)^{lambda-1/2} [C_k^lambda(u)]^2 du = pi 2^{1-2lambda} Gamma(k+2lambda) / (k! (k+lambda) [Gamma(lambda)]^2)

    print("\nComputing off-diagonal Fock-weight matrix elements on S^3")
    print("Weight w(chi) = cos^4(chi/2) = (1+cos(chi))^2 / 4\n")

    # Standard Gegenbauer normalization (for computing diagonal, needed for normalization)
    def gegenbauer_norm_sq(k, lam):
        """int_{-1}^1 (1-u^2)^{lam-1/2} [C_k^lam(u)]^2 du"""
        return pi * 2**(1 - 2*lam) * sp.gamma(k + 2*lam) / (
            factorial(k) * (k + lam) * sp.gamma(lam)**2
        )

    for n1, l1, n2, l2 in transitions:
        l = l1
        k1 = n1 - l - 1  # Gegenbauer index for state n1
        k2 = n2 - l - 1  # Gegenbauer index for state n2
        lam = l + 1       # Gegenbauer lambda parameter

        # Standard normalization integral for the S^3 hyperspherical harmonics:
        # int_{-1}^1 (1-u^2)^l * C_{k}^{l+1}(u) * C_{k}^{l+1}(u) * (1-u^2) du
        # = int_{-1}^1 (1-u^2)^{l+1} [C_k^{l+1}(u)]^2 du
        # This equals the Gegenbauer norm with lambda = l+1 and integration weight (1-u^2)^{l+1}
        # The standard Gegenbauer weight is (1-u^2)^{lambda-1/2} = (1-u^2)^{l+1/2}
        # So our weight has an extra (1-u^2)^{1/2} = sin(chi).
        #
        # Actually, let's just integrate directly with sympy for correctness.

        C_k1 = gegenbauer(k1, lam, u)
        C_k2 = gegenbauer(k2, lam, u)

        # Standard overlap (for normalization)
        # The S^3 chi-integral is: int_0^pi sin^{2l}(chi) C * C * sin^2(chi) dchi
        # = int_{-1}^1 (1-u^2)^l * C * C * (1-u^2) * |du/dchi| dchi
        # But u = cos(chi), du = -sin(chi) dchi, so dchi = -du/sin(chi), and
        # sin^2(chi) dchi = sin(chi) (-du) = sqrt(1-u^2) du
        #
        # Actually more carefully: d^3Omega_{S^3} = sin^2(chi) dchi d^2Omega
        # with chi in [0,pi], and d^2Omega integrates to 4pi * (angular Y_lm normalization).
        #
        # The chi integral alone: int_0^pi |R_n^l|^2 sin^2(chi) dchi
        # R_n^l = N sin^l(chi) C_{n-l-1}^{l+1}(cos chi)
        # So the integral is N^2 int_0^pi sin^{2l}(chi) [C(cos chi)]^2 sin^2(chi) dchi
        # Substituting u = cos chi: dchi = -du/sin(chi) = -du/sqrt(1-u^2)
        # sin^{2l+2}(chi) = (1-u^2)^{l+1}
        # So: N^2 int_{-1}^1 (1-u^2)^{l+1} [C(u)]^2 / sqrt(1-u^2) du  -- NO!
        #
        # Let me be more careful:
        # chi integral: int_0^pi f(chi) sin^2(chi) dchi
        # u = cos(chi), du = -sin(chi) dchi
        # sin^2(chi) dchi = sin(chi) * sin(chi) dchi / sin(chi) = sin(chi) * (-du) / sin(chi) WRONG
        # No: sin^2(chi) dchi, and dchi = -du/sin(chi)
        # So sin^2(chi) dchi = sin^2(chi) * (-du) / sin(chi) = sin(chi) * (-du) = sqrt(1-u^2) (-du)
        # = sqrt(1-u^2) du (flipping limits)
        #
        # So:
        # int_0^pi sin^{2l}(chi) C(cos chi)^2 sin^2(chi) dchi
        # = int_{-1}^1 (1-u^2)^l C(u)^2 sqrt(1-u^2) du
        # = int_{-1}^1 (1-u^2)^{l + 1/2} C(u)^2 du
        #
        # This IS the standard Gegenbauer orthogonality with lambda = l+1:
        # int_{-1}^1 (1-u^2)^{lambda-1/2} C_k^lambda(u)^2 du  with lambda = l+1

        # Standard diagonal norm: Gegenbauer orthogonality
        diag1_integrand = (1 - u**2)**(l + Rational(1, 2)) * C_k1**2
        diag2_integrand = (1 - u**2)**(l + Rational(1, 2)) * C_k2**2

        print(f"  Computing ({n1},{l})->({n2},{l}): k1={k1}, k2={k2}, lambda={lam}")

        diag1 = integrate(diag1_integrand, (u, -1, 1))
        diag1 = simplify(diag1)
        diag2 = integrate(diag2_integrand, (u, -1, 1))
        diag2 = simplify(diag2)

        # Off-diagonal with Fock weight cos^4(chi/2) = (1+u)^2 / 4:
        # int_0^pi sin^{2l}(chi) C_{k1}(cos chi) C_{k2}(cos chi) cos^4(chi/2) sin^2(chi) dchi
        # = int_{-1}^1 (1-u^2)^{l+1/2} C_{k1}(u) C_{k2}(u) * (1+u)^2/4 du

        offdiag_integrand = (1 - u**2)**(l + Rational(1, 2)) * C_k1 * C_k2 * (1 + u)**2 / 4
        offdiag = integrate(offdiag_integrand, (u, -1, 1))
        offdiag = simplify(offdiag)

        # Normalized coupling
        coupling = simplify(offdiag / sqrt(diag1 * diag2))

        label = f"({n1},{l})->({n2},{l})"
        print(f"    diag({n1},{l}) = {diag1}")
        print(f"    diag({n2},{l}) = {diag2}")
        print(f"    offdiag = {offdiag}")
        print(f"    coupling = {coupling}")

        # Check relationship to 1/16
        ratio_to_1_16 = simplify(coupling * 16)
        print(f"    coupling * 16 = {ratio_to_1_16}")
        print()

        results[f"C_({n1},{l})->({n2},{l})_diag1"] = str(diag1)
        results[f"C_({n1},{l})->({n2},{l})_diag2"] = str(diag2)
        results[f"C_({n1},{l})->({n2},{l})_offdiag"] = str(offdiag)
        results[f"C_({n1},{l})->({n2},{l})_coupling"] = str(coupling)
        results[f"C_({n1},{l})->({n2},{l})_coupling_x16"] = str(ratio_to_1_16)

    return results


# ============================================================================
# APPROACH D: Operator matrix elements — what gives 1/16?
# ============================================================================

def approach_D_operator_derivation():
    """
    The graph Hamiltonian is H = kappa * (D - A) where kappa = -1/16.
    For hydrogen (Z=1), H has eigenvalues -1/(2n^2) with degeneracy n^2.
    The graph Laplacian L = D - A has eigenvalues:
        lambda_n = -(n^2 - 1)  (Laplace-Beltrami on unit S^3)

    So H = kappa * L gives E_n = kappa * (-(n^2-1)) = (n^2-1)/16.

    But the actual eigenvalue is E_n = -1/(2n^2).

    The diagonal of H is -1/(2n^2) (exact hydrogenic).
    The off-diagonal is kappa * (-A) = (1/16) * A.

    So the question is: WHERE does 1/16 come from in the Fock projection?

    The Laplace-Beltrami on S^3 has eigenvalue -(n^2-1) for the n-th harmonic.
    In the graph discretization, the Laplacian L = D - A splits into:
        L_ii = d_i (degree of node i) — this is the diagonal
        L_ij = -A_ij for i ≠ j — this is the off-diagonal

    For the GeoVac graph, A_ij = 1 for adjacent states.
    The question then is: what fixes the relationship between the graph Laplacian
    eigenvalues -(n^2-1) and the physical eigenvalues -1/(2n^2)?

    E_n = kappa * lambda_n  =>  -1/(2n^2) = kappa * (-(n^2-1))

    This gives kappa = 1/(2n^2(n^2-1)), which is n-DEPENDENT, NOT -1/16.

    The resolution: the diagonal of H is NOT kappa * d_i. Instead, the
    diagonal is set independently to -Z^2/(2n^2) while the off-diagonal
    is kappa * (-A). The graph Laplacian is NOT producing the eigenvalues
    on its own — the diagonal is "by hand."

    BUT: the graph Laplacian L = D - A does have eigenvalues -(n^2-1).
    And kappa * L would give kappa * (-(n^2-1)). For this to equal -1/(2n^2)
    we would need kappa(n) = 1/(2n^2(n^2-1)), which is n-dependent.

    The ACTUAL construction uses:
        H_ii = -1/(2n^2)     [exact]
        H_ij = 1/16  for adjacent i,j

    kappa = -1/16 is the off-diagonal coupling scale. It's fixed by demanding
    that the FULL matrix (diagonal + off-diagonal) has the correct spectrum.
    For hydrogen at Z=1, the n=1 eigenvalue is -1/2, and the shift from the
    off-diagonal is what makes the spectrum exactly -1/(2n^2).

    Let me verify this by computing the full matrix spectrum for small n_max.
    And then check whether the Fock weight gives the same off-diagonal coupling.
    """
    print("\n" + "=" * 70)
    print("APPROACH D: Checking the graph spectrum and Fock weight connection")
    print("=" * 70)

    results = {}

    # For the SIMPLEST case: n_max=2, l=0 only (just s-orbitals).
    # States: (1,0,0) and (2,0,0)
    # H = [[-1/2, kappa], [kappa, -1/8]]
    # where kappa > 0 (since kappa_phys = -1/16, and H_ij = kappa_phys * (-A_ij) = 1/16 * A_ij = 1/16)

    # Eigenvalues of this 2x2 matrix:
    kappa_val = Rational(1, 16)  # off-diagonal coupling (positive, from kappa=-1/16 times -A=+1)

    H_11 = Rational(-1, 2)   # -Z^2/(2*1^2)
    H_22 = Rational(-1, 8)   # -Z^2/(2*2^2)
    H_12 = kappa_val

    # det(H - lambda I) = 0
    lam = Symbol('lambda')
    char_poly = (H_11 - lam) * (H_22 - lam) - H_12**2
    eigenvalues = sp.solve(char_poly, lam)

    print(f"\n2x2 matrix H for (1,0,0)-(2,0,0) with kappa = 1/16:")
    print(f"  H = [[{H_11}, {H_12}], [{H_12}, {H_22}]]")
    print(f"  Characteristic polynomial: {simplify(char_poly)}")
    print(f"  Eigenvalues: {eigenvalues}")
    print(f"  Exact H eigenvalues should be: -1/2, -1/8")
    print(f"  The off-diagonal 1/16 SHIFTS the eigenvalues from -1/2 and -1/8.")

    results["2x2_eigenvalues_with_kappa_1_16"] = [str(e) for e in eigenvalues]
    results["exact_eigenvalues"] = ["-1/2", "-1/8"]

    # So the eigenvalues are NOT exactly -1/2 and -1/8 with the off-diagonal!
    # The graph Hamiltonian with diagonal = -Z^2/(2n^2) and off-diagonal = 1/16
    # does NOT give exact hydrogenic eigenvalues when n_max is finite.
    # The eigenvalues converge to -1/(2n^2) only in the limit n_max -> infinity.

    # This means kappa = -1/16 is NOT determined by demanding exact eigenvalues
    # at finite n_max. It's determined by the CONTINUUM LIMIT where the graph
    # Laplacian converges to the Laplace-Beltrami operator.

    # In the continuum limit, the S^3 Laplace-Beltrami has:
    #   eigenvalues: -(n^2-1) for the n-th harmonic
    #   These are the eigenvalues of the adjacency-degree difference D - A
    #
    # The physical Hamiltonian is:
    #   H = -nabla^2/(2) - Z/r
    #   In momentum space with Fock projection:
    #   H = p^2/2 - (2p0/pi^2) * (delta_{S^3})/(p^2+p0^2)^2
    #     where delta_{S^3} is the S^3 angle-dependent term

    # Actually, the Fock relation is much simpler. Fock showed:
    #   The Schrodinger equation becomes, after projection:
    #   (p^2 + p0^2)^2 / (8pi^2 p0) * [Delta_{S^3} + (n^2-1)] psi = 0
    #
    #   on each energy shell with p0^2 = -2E = Z^2/n^2.
    #
    # The graph discretization replaces Delta_{S^3} with L = D - A.
    # The spectrum of L on the Fock graph IS -(n^2-1), EXACTLY for all n <= n_max.
    # This is proven in Paper 7.
    #
    # So the graph Laplacian already gives the EXACT eigenvalues -(n^2-1).
    # The factor kappa = -1/16 maps these to ENERGIES via:
    #   E_n = kappa * (-(n^2 - 1)) + E_shift
    #
    # But wait — that's NOT what the code does. The code uses:
    #   diagonal: -Z^2/(2n^2)   [NOT from the graph]
    #   off-diagonal: (1/16) * A   [FROM the graph]
    #
    # In the PURE graph approach (AtomicSolver), H = kappa * Z^2 * (D - A),
    # and the eigenvalues of (D-A) are -(n^2-1), so:
    #   E_n = (-Z^2/16) * (-(n^2-1)) = Z^2 (n^2-1)/16
    #
    # For n=1: E_1 = 0 (not -1/2!). The spectrum is shifted.
    # The hydrogenic spectrum is E_n = -Z^2/(2n^2) = -Z^2/2 + Z^2(n^2-1)/(2n^2).
    #
    # Actually E_n = -Z^2/(2n^2). And (n^2-1)/16 = (n^2-1)/16.
    # For n=1: (1-1)/16 = 0. For n=2: 3/16. For n=3: 8/16 = 1/2.
    #
    # While -1/(2n^2): n=1: -1/2. n=2: -1/8. n=3: -1/18.
    #
    # These are DIFFERENT sequences! The graph Laplacian gives (n^2-1)/16
    # while the exact spectrum is -1/(2n^2). So where does 1/16 come from?

    # The key is Paper 7: the graph Laplacian eigenvalues -(n^2-1) correspond
    # to the Laplace-Beltrami on S^3. The RELATIONSHIP between the S^3
    # eigenvalues and the physical energies involves the Fock projection factor:
    #
    # E_n = -Z^2/(2n^2)
    # lambda_n = -(n^2-1) = -(n^2) + 1 = -n^2 + 1
    #
    # E_n = Z^2/(2n^2) * (lambda_n / something)? No.
    #
    # The correct relation: on the energy shell, p0 = Z/n. The S^3 radius is
    # set by p0. The eigenvalue equation on S^3 at radius 1 gives n^2-1.
    # The physical energy is E_n = -p0^2/2 = -Z^2/(2n^2).
    #
    # So the mapping from S^3 eigenvalue to energy is:
    #   n^2 - 1  =>  E_n = -Z^2/(2n^2)
    #   (n^2-1) = n^2 - 1  =>  n^2 = (n^2-1) + 1 = |lambda_n| + 1
    #   E_n = -Z^2 / (2(|lambda_n|+1))
    #
    # This is a NONLINEAR mapping — not a simple scale factor!
    # So kappa is NOT the ratio E_n / lambda_n (which would be n-dependent).

    # Let me compute: for the pure graph Laplacian approach where H = kappa * (D-A),
    # what value of kappa gives the best approximation to -1/(2n^2)?

    print("\n--- Pure graph spectrum: kappa * (-(n^2-1)) vs -1/(2n^2) ---")
    for n in range(1, 6):
        lambda_n = -(n**2 - 1)
        E_exact = Rational(-1, 2*n**2)
        if n > 1:
            kappa_needed = E_exact / lambda_n
            print(f"  n={n}: lambda={lambda_n}, E_exact={E_exact}, kappa_needed={kappa_needed} = {float(kappa_needed):.6f}")
        else:
            print(f"  n={n}: lambda={lambda_n}, E_exact={E_exact}, kappa_needed=undefined (0/0)")

    # Let me check the n=2 case: kappa_needed = (-1/8) / (-3) = 1/24
    # n=3: (-1/18) / (-8) = 1/144
    # These are 1/(2n^2(n^2-1)/n^2) ... no: kappa = 1/(2n^2) * 1/(n^2-1) * n^2 = ...
    # kappa = E/lambda = 1/(2(n^2-1)n^2/n^2)... let me just be explicit:
    # kappa(n) = (-1/(2n^2)) / (-(n^2-1)) = 1/(2n^2(n^2-1))
    # n=2: 1/(2*4*3) = 1/24
    # n=3: 1/(2*9*8) = 1/144
    # n=4: 1/(2*16*15) = 1/480
    # These are NOT 1/16.

    # So the "pure" graph H = kappa*(D-A) with constant kappa does NOT
    # give the exact hydrogen spectrum.

    # What kappa = -1/16 achieves is: with DIAGONAL = -Z^2/(2n^2) set by hand,
    # the off-diagonal coupling of 1/16 gives eigenvalues that CONVERGE to
    # -1/(2n^2) as n_max -> infinity. But at any finite n_max, the eigenvalues
    # are shifted.

    # The GRAPH-NATIVE approach uses this hybrid: exact diagonal + graph off-diagonal.
    # kappa = -1/16 is then fixed by demanding that in the n_max -> infinity limit,
    # the full matrix reproduces the correct spectrum.

    # In the continuum limit, the off-diagonal coupling between adjacent shells
    # from the S^3 Laplace-Beltrami operator (in the Fock-projected formalism)
    # is what determines kappa. Let's compute this!

    # The S^3 Laplacian in the Gegenbauer basis has matrix elements:
    # <Y_{n'lm}| Delta_{S^3} |Y_{nlm}> = -(n^2-1) delta_{nn'}
    # It IS diagonal — the Gegenbauer polynomials are eigenfunctions.

    # But the GRAPH Laplacian D-A in the standard (hydrogenic) basis IS also
    # diagonal with eigenvalue -(n^2-1) for n=1..n_max.

    # So the off-diagonal coupling in the graph Hamiltonian h1 comes from
    # the MISMATCH between the hydrogenic (n,l,m) labeling and the graph
    # adjacency structure.

    # Specifically: the graph connects (n,l,m) to (n+1,l,m) with weight 1.
    # The Laplacian D-A has eigenvalue -(n^2-1) but the eigenvectors are
    # NOT the standard basis vectors |n,l,m>. The |n,l,m> states are
    # eigenvectors of both the S^3 Laplacian AND the degree matrix D,
    # so D-A having them as eigenvectors means A must also be simultaneously
    # diagonalizable in this basis.

    # Wait — that CAN'T be right if A connects different n. Let me think again.

    # The S^3 Laplacian eigenvalue -(n^2-1) has degeneracy n^2 (all l,m with
    # l=0..n-1, m=-l..l). WITHIN a given eigenspace, the states |n,l,m>
    # may be connected by the adjacency (m↔m±1 within same l), but the
    # inter-shell connections (n↔n±1) go BETWEEN different eigenspaces.

    # The graph Laplacian L = D - A is block-diagonal in the Fock eigenspaces
    # ONLY if A is. But A connects different n, so L is NOT block-diagonal.

    # RESOLUTION: Paper 7 proves that the spectrum of (D-A) on the GeoVac graph
    # IS {-(n^2-1) with degeneracy n^2} for all n<=n_max. This means the
    # eigenvectors of D-A mix different (l,m) within the same n-shell but
    # NOT different n-shells. The inter-shell connections from A are exactly
    # cancelled by the degree D to produce within-shell eigenvalues.

    # This is possible because the degree d_i of node (n,l,m) is:
    #   d_{n,l,m} = (number of m±1 neighbors) + (number of n±1 neighbors)
    #             = (min(1,l-m) + min(1,l+m)) + (1 if n>1) + (1 if n<n_max)
    # For interior nodes (1 < n < n_max, -l < m < l): d = 2 + 2 = 4
    # For edge cases it varies.

    # The key insight is that for the INFINITE graph (n_max -> infinity),
    # D-A has the spectrum -(n^2-1) with the degeneracy structure.

    # SO: kappa = -1/16 maps the S^3 Laplacian eigenvalues to the
    # hydrogenic energies via:
    #   H = kappa * Z^2 * (D - A) + shift
    #   E_n = kappa * Z^2 * (-(n^2-1))
    #
    # For n=1: E_1 = 0 (the shift is needed to get -1/2).
    # For n=2: E_2 = -Z^2 * 3/16 (vs exact -Z^2/8)

    # Hmm, but Z^2*3/16 = 3/16 while exact is 1/8 = 2/16. These are different!

    # Let me re-examine. The AtomicSolver code says:
    #   kinetic_scale = -1/16 * Z^2
    #   H = kinetic_scale * laplacian
    # where laplacian = D - A.
    #
    # eigenvalues of laplacian = -(n^2-1) for the GeoVac graph.
    # So E_n = (-Z^2/16) * (-(n^2-1)) = Z^2(n^2-1)/16.
    #
    # This gives: E_1 = 0, E_2 = 3/16, E_3 = 1/2, E_4 = 15/16, E_5 = 3/2
    # while exact: E_1 = -1/2, E_2 = -1/8, E_3 = -1/18, ...
    #
    # The relationship is: E_graph = Z^2(n^2-1)/16 = -1/(2n^2) + 1/2
    # Wait: -1/(2n^2) + 1/2 = (n^2-1)/(2n^2). That's not (n^2-1)/16.
    #
    # Unless... the Graph Laplacian eigenvalues are NOT -(n^2-1)?

    # Let me check with the actual code.

    print("\n--- Checking actual graph Laplacian eigenvalues ---")

    try:
        from geovac.lattice import GeometricLattice
        import numpy as np
        from scipy import sparse

        for n_max in [3, 5, 10]:
            lattice = GeometricLattice(max_n=n_max)
            A_mat = lattice.adjacency
            if sparse.issparse(A_mat):
                A_dense = A_mat.toarray()
            else:
                A_dense = np.array(A_mat)

            D_mat = np.diag(A_dense.sum(axis=1))
            L = D_mat - A_dense
            eigenvalues = np.sort(np.linalg.eigvalsh(L))

            # Group by approximate value
            unique_evals = []
            for ev in eigenvalues:
                if not unique_evals or abs(ev - unique_evals[-1][0]) > 1e-8:
                    unique_evals.append((ev, 1))
                else:
                    unique_evals[-1] = (unique_evals[-1][0], unique_evals[-1][1] + 1)

            print(f"\n  n_max={n_max}: graph Laplacian eigenvalues (value, degeneracy):")
            for val, deg in unique_evals[:8]:
                # Check if val = n^2-1 for some n
                n_sq_m1 = val
                n_check = sp.sqrt(n_sq_m1 + 1)
                print(f"    {val:.4f} (deg {deg}), n^2-1={val:.4f}, n={float(n_check):.4f}")

            # Compute H = (-1/16) * L
            H = (-1/16) * L
            H_eigenvalues = np.sort(np.linalg.eigvalsh(H))

            print(f"  H = (-1/16)*L eigenvalues:")
            unique_H = []
            for ev in H_eigenvalues:
                if not unique_H or abs(ev - unique_H[-1][0]) > 1e-8:
                    unique_H.append((ev, 1))
                else:
                    unique_H[-1] = (unique_H[-1][0], unique_H[-1][1] + 1)

            for val, deg in unique_H[:8]:
                # Compare to -1/(2n^2)
                print(f"    {val:.6f} (deg {deg})")

            print(f"  Exact -1/(2n^2) for n=1..{min(n_max,5)}:")
            for n in range(1, min(n_max+1, 6)):
                print(f"    -1/(2*{n}^2) = {-1/(2*n**2):.6f}")

    except ImportError:
        print("  (Could not import geovac.lattice — skipping numerical check)")

    results["approach_D_note"] = (
        "The graph Laplacian eigenvalues -(n^2-1) do NOT simply scale by kappa=-1/16 "
        "to give the hydrogen spectrum -1/(2n^2). The relationship is more subtle: "
        "the hybrid h1 matrix (exact diagonal + graph off-diagonal at kappa) gives "
        "a matrix whose spectrum converges to the correct values."
    )

    return results


# ============================================================================
# APPROACH E: Off-diagonal from Gegenbauer recursion
# ============================================================================

def approach_E_gegenbauer_recursion():
    """
    The Gegenbauer polynomials satisfy the three-term recursion:

        (k+1) C_{k+1}^lambda(x) = 2(k+lambda) x C_k^lambda(x) - (k+2*lambda-1) C_{k-1}^lambda(x)

    On S^3, the hyperspherical harmonic R_n^l(chi) = sin^l(chi) C_{n-l-1}^{l+1}(cos chi).

    The adjacency connects (n,l,m) to (n+1,l,m), which means we need the
    matrix element of the identity operator (or some weight) between
    R_n^l and R_{n+1}^l on S^3.

    But more directly: the Gegenbauer recursion with x = cos(chi) connects
    C_{k+1}^lambda with C_k^lambda via multiplication by cos(chi).
    So cos(chi) is the "adjacency" operator that couples adjacent Gegenbauer indices.

    The matrix element <R_{n+1}^l | cos(chi) | R_n^l>_{S^3} should be a clean rational.

    Note: cos(chi) in Fock's projection is:
        cos(chi) = (p0^2 - p^2) / (p0^2 + p^2)  [at p = p0 tan(chi/2)]

    So cos(chi) is related to the Fock energy variable!

    Also: cos^2(chi/2) = (1 + cos(chi))/2 = p0^2 / (p^2 + p0^2)

    The Fock weight (p^2+p0^2)^{-2} = 1/(p0^4) * cos^4(chi/2)
    = (1/p0^4) * ((1+cos chi)/2)^2.

    So the weight w = cos^4(chi/2) = (1+cos chi)^2/4, which involves
    1 + 2*cos(chi) + cos^2(chi). Using the Gegenbauer recursion,
    cos(chi) couples k to k±1, and cos^2(chi) couples to k±2 (and diagonal).

    The off-diagonal matrix elements of the weight are DETERMINED by the
    Gegenbauer recursion coefficients and the weight structure.
    """
    print("\n" + "=" * 70)
    print("APPROACH E: Gegenbauer recursion for cos(chi) couplings")
    print("=" * 70)

    results = {}
    u = Symbol('u', real=True)

    # Gegenbauer recursion: (k+1) C_{k+1}^lam(u) = 2(k+lam) u C_k^lam(u) - (k+2lam-1) C_{k-1}^lam(u)
    # => u C_k^lam(u) = [(k+1)/(2(k+lam))] C_{k+1}^lam(u) + [(k+2lam-1)/(2(k+lam))] C_{k-1}^lam(u)

    # Matrix element: <C_{k+1}^lam | u | C_k^lam>_w where w = (1-u^2)^{lam-1/2}
    # Using the recursion, this equals:
    # <C_{k+1} | [(k+1)/(2(k+lam))] C_{k+1} + [(k+2lam-1)/(2(k+lam))] C_{k-1}>_w
    # = [(k+1)/(2(k+lam))] * <C_{k+1}|C_{k+1}>_w + 0  (orthogonality kills second term)
    # = [(k+1)/(2(k+lam))] * h_{k+1}^lam
    #
    # where h_k^lam = int_{-1}^1 (1-u^2)^{lam-1/2} [C_k^lam(u)]^2 du = pi 2^{1-2lam} Gamma(k+2lam) / [k! (k+lam) Gamma(lam)^2]

    # For the S^3 inner product (with the extra sin^{2l}(chi) = (1-u^2)^l factor),
    # the weight is (1-u^2)^{l+1/2} which means lambda = l+1.

    print("\nGegenbauer recursion coefficients for cos(chi) coupling:")
    print(f"{'(n,l) -> (n+1,l)':<20} {'k':<5} {'lam':<5} {'a_{k,k+1} = (k+1)/(2(k+lam))':<40}")
    print("-" * 70)

    transitions = [
        (1, 0), (2, 0), (3, 0),  # l=0 transitions
        (2, 1), (3, 1),          # l=1 transitions
        (3, 2),                  # l=2 transition
    ]

    for n, l in transitions:
        k = n - l - 1     # Gegenbauer index
        lam = l + 1        # Gegenbauer lambda

        # Recursion coefficient for cos(chi) coupling k -> k+1
        a_kp1 = Rational(k + 1, 2 * (k + lam))

        # Also the k-1 coefficient (for completeness)
        if k > 0:
            a_km1 = Rational(k + 2*lam - 1, 2*(k + lam))
        else:
            a_km1 = None

        label = f"({n},{l})->({n+1},{l})"
        print(f"{label:<20} {k:<5} {lam:<5} {str(a_kp1):<40}")

        results[f"E_cos_chi_({n},{l})->({n+1},{l})"] = str(a_kp1)

    # Now compute the FULL coupling from the Fock weight w = (1+u)^2/4 = 1/4 + u/2 + u^2/4.
    #
    # <C_{k+1}^lam | w(u) | C_k^lam>_w = (1/4)<C_{k+1}|C_k>_w + (1/2)<C_{k+1}|u C_k>_w + (1/4)<C_{k+1}|u^2 C_k>_w
    #
    # Term 1: <C_{k+1}|C_k>_w = 0 by orthogonality (different indices)
    # Term 2: (1/2)<C_{k+1}|u C_k>_w = (1/2) * a_{k,k+1} * h_{k+1}  [from recursion above]
    # Term 3: (1/4)<C_{k+1}|u^2 C_k>_w — need u^2 expansion

    # For u^2: apply recursion twice, or use:
    # u^2 C_k^lam = u * [a_{k+1} C_{k+1} + a_{k-1} C_{k-1}]
    # = a_{k+1} [a'_{k+2} C_{k+2} + a'_{k} C_k] + a_{k-1} [a''_{k} C_k + a''_{k-2} C_{k-2}]
    # where a' coefficients are from the recursion applied to C_{k+1}, etc.

    # <C_{k+1} | u^2 | C_k>_w: the u^2 term couples C_k to C_{k-2}, C_k, C_{k+2}
    # The only term that survives in <C_{k+1}|...> is:
    # from u * a_{k,k+1} C_{k+1}: u C_{k+1} = a_{k+1,k+2} C_{k+2} + a_{k+1,k} C_k
    # => the C_{k+1} component of u^2 C_k is a_{k,k+1} * a_{k+1,k} * h_{k+1}? No.
    # Actually <C_{k+1}| u C_{k+1}> = a_{k+1,k+2} h_{k+2} / h_{k+1} * h_{k+1}
    # Wait, I need to be more careful with the inner product.

    # Let me just use the recursion directly:
    # u C_k = a_{k,k+1} C_{k+1} + a_{k,k-1} C_{k-1}
    # u * u C_k = u * a_{k,k+1} C_{k+1} + u * a_{k,k-1} C_{k-1}
    # = a_{k,k+1} [a_{k+1,k+2} C_{k+2} + a_{k+1,k} C_k] + a_{k,k-1} [a_{k-1,k} C_k + a_{k-1,k-2} C_{k-2}]
    #
    # The coefficient of C_{k+1} in u^2 C_k is ZERO! (The recursion u C_j = a C_{j+1} + b C_{j-1}
    # only connects j to j±1, so u^2 C_k connects to k, k±2, never k+1.)

    # THEREFORE: <C_{k+1} | u^2 C_k>_w = 0 !!!

    # So the Fock weight coupling <C_{k+1} | w | C_k>_w with w = (1+u)^2/4 is:
    # = 0 + (1/2) a_{k,k+1} h_{k+1} + 0
    # = (1/2) a_{k,k+1} h_{k+1}
    #
    # The NORMALIZED coupling (dividing by sqrt(h_k * h_{k+1})) is:
    # = (1/2) a_{k,k+1} sqrt(h_{k+1} / h_k)

    print("\n--- IMPORTANT RESULT ---")
    print("The u^2 term in the Fock weight w = (1+u)^2/4 = 1/4 + u/2 + u^2/4")
    print("gives ZERO off-diagonal coupling (u^2 couples k to k±2, not k to k+1).")
    print("The constant term 1/4 also gives zero (orthogonality).")
    print("Only the u/2 = cos(chi)/2 term contributes!")
    print()
    print("The Fock weight off-diagonal coupling between k and k+1 is:")
    print("  <k+1| w |k>_w = (1/2) * (k+1)/(2(k+lam)) * h_{k+1}^lam")
    print()

    # Gegenbauer norm: h_k^lam = pi * 2^{1-2lam} * Gamma(k+2lam) / (k! * (k+lam) * Gamma(lam)^2)
    def h_k_lam(k_val, lam_val):
        return (pi * 2**(1 - 2*lam_val) * sp.gamma(k_val + 2*lam_val)
                / (factorial(k_val) * (k_val + lam_val) * sp.gamma(lam_val)**2))

    print(f"{'(n,l)->(n+1,l)':<20} {'normalized coupling':<25} {'coupling * 16':<15}")
    print("-" * 60)

    for n, l in transitions:
        k = n - l - 1
        lam = l + 1

        a_kp1 = Rational(k + 1, 2 * (k + lam))

        h_k = simplify(h_k_lam(k, lam))
        h_kp1 = simplify(h_k_lam(k + 1, lam))

        # Unnormalized coupling: (1/2) * a_{k,k+1} * h_{k+1}
        unnorm_coupling = Rational(1, 2) * a_kp1 * h_kp1

        # Normalized coupling: divide by sqrt(h_k * h_{k+1})
        norm_coupling = simplify(unnorm_coupling / sqrt(h_k * h_kp1))
        norm_coupling = simplify(norm_coupling)

        # Simplify further
        norm_coupling2 = simplify(Rational(1, 2) * a_kp1 * sqrt(h_kp1 / h_k))
        norm_coupling2 = simplify(norm_coupling2)

        # x16
        coupling_x16 = simplify(norm_coupling * 16)

        label = f"({n},{l})->({n+1},{l})"
        print(f"{label:<20} {str(norm_coupling):<25} {str(coupling_x16):<15}")

        results[f"E_fock_coupling_({n},{l})->({n+1},{l})_normalized"] = str(norm_coupling)
        results[f"E_fock_coupling_({n},{l})->({n+1},{l})_x16"] = str(coupling_x16)

    # Let me also compute this in a fully explicit form
    print("\n--- Explicit Gegenbauer norm ratios ---")
    for n, l in transitions:
        k = n - l - 1
        lam = l + 1

        h_k = simplify(h_k_lam(k, lam))
        h_kp1 = simplify(h_k_lam(k + 1, lam))
        ratio = simplify(h_kp1 / h_k)

        print(f"  ({n},{l}): h_{k}^{lam} = {h_k}")
        print(f"  ({n},{l}): h_{k+1}^{lam} = {h_kp1}")
        print(f"  ({n},{l}): ratio h_{k+1}/h_{k} = {ratio}")
        print()

    return results


# ============================================================================
# Main
# ============================================================================

def main():
    all_results = {}

    # Approach C is the most direct and clean
    print("Starting with Approach C (direct S^3 Gegenbauer integration)")
    print("This is the most reliable approach.\n")

    results_C = approach_C_gegenbauer_on_S3()
    all_results["approach_C"] = results_C

    # Approach E uses Gegenbauer recursion (algebraic, no integration needed)
    results_E = approach_E_gegenbauer_recursion()
    all_results["approach_E"] = results_E

    # Approach D checks the actual graph spectrum
    results_D = approach_D_operator_derivation()
    all_results["approach_D"] = results_D

    # Save results
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)
    json_path = os.path.join(data_dir, "probe_k1_fock_weight.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {json_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print("The key finding is from Approach E: the Gegenbauer recursion shows that")
    print("the Fock weight cos^4(chi/2) = (1+cos chi)^2/4 has off-diagonal coupling")
    print("between adjacent Gegenbauer indices (k, k+1) given ENTIRELY by the")
    print("cos(chi)/2 term (the u^2 term gives zero coupling for adjacent indices).")
    print()
    print("The normalized coupling is: (1/2) * (k+1)/(2(k+lam)) * sqrt(h_{k+1}/h_k)")
    print("where k = n-l-1 and lam = l+1.")
    print()
    print("For the specific (n,l) combinations:")
    for n, l in [(1,0), (2,0), (3,0), (2,1), (3,1), (3,2)]:
        k = n - l - 1
        lam = l + 1
        a = Rational(k+1, 2*(k+lam))
        key = f"E_fock_coupling_({n},{l})->({n+1},{l})_normalized"
        if key in results_E:
            print(f"  ({n},{l})->({n+1},{l}): a = {a} = {float(a):.6f}, full coupling = {results_E[key]}")

    print()
    print("The cos(chi) recursion coefficient a_{k,k+1} = (k+1)/(2(k+lam)) is:")
    print("  (1,0)->2: 1/2")
    print("  (2,0)->3: 1/2")
    print("  (3,0)->4: 1/2")
    print("  (2,1)->3: 1/4")
    print("  (3,1)->4: 1/3")
    print("  (3,2)->4: 1/6")
    print()
    print("These are (n-l)/(2n) = k+1/(2(k+l+1)).")
    print("For l=0: always 1/2 (independent of n).")
    print("The coupling is n-dependent and l-dependent, NOT constant 1/16.")

    return all_results


if __name__ == "__main__":
    main()
