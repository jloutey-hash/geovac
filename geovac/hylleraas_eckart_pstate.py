"""
Hylleraas-Eckart P-state extension (Schwartz 1961 form)
=======================================================

P-state (L=1) two-electron trial wavefunctions for the Eckart double-alpha
basis. The Schwartz 1961 ansatz for L=1, M=0:

    Psi_{2^1P}^{M=0} = (z_1 + z_2) * chi_+(s, t, u)
                       + (z_1 - z_2) * chi_-(s, t, u)

where
  chi_+(s, t, u) = exp(-alpha s) * cosh(beta t) * s^l * t^{2m} * u^n
                   (singlet, r_1 <-> r_2 SYMMETRIC; "+" channel)
  chi_-(s, t, u) = exp(-alpha s) * sinh(beta t) * s^l * t^{2m} * u^n
                   (singlet, r_1 <-> r_2 ANTISYMMETRIC; "-" channel)

The (z_1 + z_2) factor is symmetric under r_1 <-> r_2 and pairs with the
symmetric radial piece (cosh, even t-power); the (z_1 - z_2) factor is
antisymmetric and pairs with sinh.

Track 5 of the Hylleraas-Eckart sprint (CLAUDE.md §2 backlog entry).
Headline application: He 2^1P -> 1^1S oscillator strength (Drake handbook
f = 0.276; extended-CI baseline at f = 0.286, +3.4%; Hylleraas r_12 explicit
correlation is the structural closure).

Session 1 scope (this module):
* P-state symmetric ("+") channel basis representation
* Overlap, V_ne, V_ee via index-shifted masters from
  hylleraas_eckart_recurrence.py
* Validate against a single-basis-function known reference

Future-session scope:
* Antisymmetric ("-") channel; full 2D singlet 2^1P basis
* Kinetic energy (Hartree-form expansion picks up extra terms from
  gradient acting on the (z_1 + z_2) factor)
* Cross-sector matrix elements <1^1S | z_1 + z_2 | 2^1P> -> dipole
* Variational driver for 2^1P
* Oscillator strength assembly

Angular reduction
-----------------

For Schwartz form Phi = (z_1 + z_2) chi(r_1, r_2, r_12), after SO(3)
averaging over the overall orientation (3 of the 6 angles), the
P-state matrix elements reduce to S-state-style master integrals
with an extra polynomial-in-(s,t,u) factor:

    < (z_1+z_2)^2 >_SO(3)  =  (1/3) (2 r_1^2 + 2 r_2^2 - r_12^2)
                            =  (1/3) (s^2 + t^2 - u^2)

(Using r_1^2 + r_2^2 = (s^2+t^2)/2 and r_12 = u.)

So for two basis functions Phi_p = (z_1+z_2) chi_p and Phi_q = (z_1+z_2) chi_q
in the SYMMETRIC channel (both cosh-factored, even t-power):

    <Phi_p | Phi_q> = (pi^2 / 3) * int (s^2 + t^2 - u^2) chi_p chi_q
                                          * u (s^2 - t^2) ds dt du
                    = (pi^2 / 3) * [I(L+2, M, N) + I(L, M+1, N) - I(L, M, N+2)]
                                   evaluated with the singlet cosh combination
                                   (1/2)[at B_+ + at B_-]

with L = l_p + l_q, M = m_p + m_q, N = n_p + n_q, and:
    I(L, M, N; alpha, B)
      = int exp(-2 alpha s) s^L u^{N+1} t^{2M} (s^2 - t^2) cosh(B t) dt du ds

V_ne and V_ee follow the same pattern with the J and K masters respectively
(see hylleraas_eckart_closed_forms.py).

References
----------
* Schwartz, C. (1961). "Lamb shift in the helium atom." Phys. Rev. 123, 1700.
  AND Schwartz, C. (1961). J. Math. Phys. 2, 568. (Hylleraas-Schwartz form
  for excited P-states.)
* Drake, G. W. F. (1996). Atomic, Molecular, and Optical Physics Handbook.
  AIP. Sec. 11. (Pekeris-Drake reductions for L > 0 states.)
* Pekeris, C. L. (1958). Phys. Rev. 112, 1649. (Perimetric coordinates;
  P-state extension in subsequent papers.)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from geovac.hylleraas_r12 import (
    HylleraasBasisFn,
    hylleraas_basis_total_degree,
    hylleraas_volume_element_factor,
    overlap_element_eckart,
    _ke_term_contributions,
)
from geovac.hylleraas_eckart_closed_forms import (
    I_HE_cosh_polynomial,
    J_HE_cosh_polynomial,
    K_HE_cosh_polynomial,
)
from geovac.hylleraas_eckart_recurrence import master_C_gen, master_S_gen


# ---------------------------------------------------------------------------
# P-state basis function
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HylleraasPStateBasisFn:
    """A single P-state basis function in the Schwartz 1961 form.

    phi_{l, m, n}^{P, sym}(s, t, u) = (z_1 + z_2) * exp(-alpha s)
                                       * cosh(beta t) * s^l * t^{2m} * u^n
    phi_{l, m, n}^{P, antisym}(s, t, u) = (z_1 - z_2) * exp(-alpha s)
                                            * sinh(beta t) * s^l * t^{2m} * u^n

    The (l, m, n) indices are the radial polynomial powers; the angular L=1
    M=0 content is in the (z_1 +/- z_2) prefactor.

    channel:
      'sym':    (z_1 + z_2) * cosh(beta t) * t^{2m}  (singlet symmetric)
      'antisym': (z_1 - z_2) * sinh(beta t) * t^{2m}  (singlet antisymmetric)

    Both channels can coexist in a single 2^1P trial wavefunction
    (Schwartz 1961).
    """

    l: int
    m: int
    n: int
    channel: str = 'sym'  # 'sym' or 'antisym'

    def __post_init__(self) -> None:
        if self.l < 0 or self.m < 0 or self.n < 0:
            raise ValueError(
                f"l, m, n must be >= 0; got ({self.l}, {self.m}, {self.n})"
            )
        if self.channel not in ('sym', 'antisym'):
            raise ValueError(
                f"channel must be 'sym' or 'antisym'; got {self.channel!r}"
            )

    def label(self) -> str:
        suffix = '+' if self.channel == 'sym' else '-'
        return f"P{suffix}_s{self.l}t{2 * self.m}u{self.n}"


def hylleraas_pstate_basis_total_degree(
    omega_max: int, channel: str = 'sym',
) -> List[HylleraasPStateBasisFn]:
    """All P-state basis functions in one channel with l + 2m + n <= omega_max."""
    if channel not in ('sym', 'antisym'):
        raise ValueError(f"channel must be 'sym' or 'antisym'; got {channel!r}")
    if omega_max < 0:
        raise ValueError(f"omega_max must be >= 0; got {omega_max}")

    basis: List[HylleraasPStateBasisFn] = []
    for l in range(omega_max + 1):
        for m in range((omega_max - l) // 2 + 1):
            for n in range(omega_max - l - 2 * m + 1):
                basis.append(HylleraasPStateBasisFn(l, m, n, channel=channel))
    return basis


# ---------------------------------------------------------------------------
# Angular reduction factor: <(z_1 + z_2)^2>_SO(3) = (s^2 + t^2 - u^2) / 3
#
# The pi^2 / 3 prefactor combines the standard S-state volume factor pi^2
# with the 1/3 from SO(3) averaging.
# ---------------------------------------------------------------------------

P_STATE_VOLUME_FACTOR = math.pi ** 2 / 3.0


# ---------------------------------------------------------------------------
# Symmetric-channel matrix elements
# (Both basis functions are sym; cosh(beta_p t) cosh(beta_q t) hyperbolic
# product gives the standard singlet cosh combination.)
# ---------------------------------------------------------------------------

def _eckart_cosh_combination_pstate(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    master_fn,  # I_HE_cosh_polynomial, J_HE_cosh_polynomial, or K_HE_cosh_polynomial
    L_extra: int = 0, M_extra: int = 0, N_extra: int = 0,
) -> float:
    """Singlet cosh combination of the master integral with optional
    (L, M, N) shifts:

        (1/2) [ master(L+L_extra, M+M_extra, N+N_extra; alpha, 2 beta)
                + master(L+L_extra, M+M_extra, N+N_extra; alpha, 0) ]

    Used to assemble P-state matrix elements where the basis function
    indices come from bf_p and bf_q and the geometric (s^2, t^2, u^2)
    factors of the angular reduction add (L_extra, M_extra, N_extra)
    shifts.

    Both bf_p and bf_q must be channel='sym' (this helper does not handle
    the antisymmetric channel).
    """
    L = bf_p.l + bf_q.l + L_extra
    M = bf_p.m + bf_q.m + M_extra
    N = bf_p.n + bf_q.n + N_extra
    if L < 0 or M < 0 or N < 0:
        raise ValueError(
            f"Negative master-integral index after shift: "
            f"({L}, {M}, {N}) with shifts ({L_extra}, {M_extra}, {N_extra})"
        )
    cf = master_fn(L, M, N)
    v_plus = cf.evaluate(alpha, 2.0 * beta)
    v_minus = cf.evaluate(alpha, 0.0)
    return 0.5 * (v_plus + v_minus)


def overlap_element_pstate_eckart(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Symmetric-channel P-state overlap <(z_1+z_2)*chi_p | (z_1+z_2)*chi_q>.

    After SO(3) averaging, <(z_1+z_2)^2> = (s^2 + t^2 - u^2)/3, so the
    overlap reduces to a linear combination of I-master integrals:

      <Phi_p | Phi_q> = (pi^2/3) [I(L+2, M, N) + I(L, M+1, N) - I(L, M, N+2)]
                       * singlet cosh combination

    where (L, M, N) = (l_p+l_q, m_p+m_q, n_p+n_q).
    """
    assert bf_p.channel == 'sym' and bf_q.channel == 'sym'

    # The volume factor (s^2-t^2) is built into the I master integral
    # (which has integrand u^{N+1} t^{2M} (s^2-t^2) cosh(Bt)). Adding
    # the angular-reduction factor (s^2 + t^2 - u^2):
    #   s^2 contribution: shift L by +2 -> I(L+2, M, N)
    #   t^2 contribution: shift M by +1 (t^2 * t^{2M} = t^{2(M+1)}) -> I(L, M+1, N)
    #   -u^2 contribution: shift N by +2 (u^2 * u^{N+1} = u^{N+3}, which is u^{(N+2)+1}) -> -I(L, M, N+2)
    I_s2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, I_HE_cosh_polynomial,
        L_extra=2, M_extra=0, N_extra=0,
    )
    I_t2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, I_HE_cosh_polynomial,
        L_extra=0, M_extra=1, N_extra=0,
    )
    I_u2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, I_HE_cosh_polynomial,
        L_extra=0, M_extra=0, N_extra=2,
    )
    return P_STATE_VOLUME_FACTOR * (I_s2 + I_t2 - I_u2)


def potential_vne_element_pstate_eckart(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float, Z: float,
) -> float:
    """Symmetric-channel P-state nuclear attraction.

    <Phi_p | -Z (1/r_1 + 1/r_2) | Phi_q>
      = -4 Z (pi^2/3) [J(L+3, M, N) + J(L+1, M+1, N) - J(L+1, M, N+2)]
        * singlet cosh combination

    The J master integral already has the L+1 shift baked in (since
    -4Z appears with s/(s^2-t^2) cancellation); we add an extra L+2 (so
    total +3), M+1, or N+2 from the angular reduction.
    """
    assert bf_p.channel == 'sym' and bf_q.channel == 'sym'

    # The J master is the V_ne master where the L+1 shift (extra s factor)
    # is the caller's responsibility. Here we add L_extra in addition to
    # the +1 shift that potential_vne_element_eckart applies for S-state.
    # Looking at potential_vne_element_eckart, it calls J_HE_cosh_polynomial(L+1, M, N).
    # For P-state symmetric, we add (L_extra=2, M_extra=0, N_extra=0) from s^2,
    # (L_extra=0, M_extra=1, N_extra=0) from t^2,
    # and (L_extra=0, M_extra=0, N_extra=2) from -u^2.
    J_s2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, J_HE_cosh_polynomial,
        L_extra=3, M_extra=0, N_extra=0,
    )
    J_t2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, J_HE_cosh_polynomial,
        L_extra=1, M_extra=1, N_extra=0,
    )
    J_u2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, J_HE_cosh_polynomial,
        L_extra=1, M_extra=0, N_extra=2,
    )
    return -4.0 * Z * P_STATE_VOLUME_FACTOR * (J_s2 + J_t2 - J_u2)


def potential_vee_element_pstate_eckart(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Symmetric-channel P-state electron-electron repulsion.

    <Phi_p | 1/r_12 | Phi_q>
      = (pi^2/3) [K(L+2, M, N) + K(L, M+1, N) - K(L, M, N+2)]
        * singlet cosh combination
    """
    assert bf_p.channel == 'sym' and bf_q.channel == 'sym'

    K_s2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, K_HE_cosh_polynomial,
        L_extra=2, M_extra=0, N_extra=0,
    )
    K_t2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, K_HE_cosh_polynomial,
        L_extra=0, M_extra=1, N_extra=0,
    )
    K_u2 = _eckart_cosh_combination_pstate(
        bf_p, bf_q, alpha, beta, K_HE_cosh_polynomial,
        L_extra=0, M_extra=0, N_extra=2,
    )
    return P_STATE_VOLUME_FACTOR * (K_s2 + K_t2 - K_u2)


# ---------------------------------------------------------------------------
# Antisymmetric-channel matrix elements
# ---------------------------------------------------------------------------
#
# Antisymmetric channel basis: phi^(-) = (z_1 - z_2) e^{-alpha s} sinh(beta t)
#                                          * s^l t^{2m} u^n
#
# Angular reduction:
#   <(z_1 - z_2)^2>_SO(3) = (1/3)(r_1^2 + r_2^2 - 2 r_1 r_2 cos(theta_12))
#                         = r_12^2 / 3
#                         = u^2 / 3
#
# So antisym × antisym matrix elements use master integrals with N
# shifted by +2 (from the u^2 factor), and the "triplet" cosh combination
# (sinh × sinh = (1/2)[cosh(B+) - cosh(B-)]).
#
# At beta = 0, sinh(0) = 0 so phi^(-) vanishes identically, and the
# triplet cosh combo gives zero — consistent with the basis collapsing
# at beta = 0.

def _eckart_triplet_combination_pstate(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    master_fn,
    L_extra: int = 0, M_extra: int = 0, N_extra: int = 0,
) -> float:
    """Triplet (sinh*sinh) cosh combination of master integrals with optional shifts.

        (1/2) [ master(L+L_extra, M+M_extra, N+N_extra; alpha, 2 beta)
                - master(L+L_extra, M+M_extra, N+N_extra; alpha, 0) ]
    """
    L = bf_p.l + bf_q.l + L_extra
    M = bf_p.m + bf_q.m + M_extra
    N = bf_p.n + bf_q.n + N_extra
    if L < 0 or M < 0 or N < 0:
        raise ValueError(
            f"Negative master-integral index after shift: ({L}, {M}, {N})"
        )
    cf = master_fn(L, M, N)
    v_plus = cf.evaluate(alpha, 2.0 * beta)
    v_minus = cf.evaluate(alpha, 0.0)
    return 0.5 * (v_plus - v_minus)


def overlap_element_pstate_eckart_antisym(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Antisymmetric-channel P-state overlap.

    <Phi^(-)_p | Phi^(-)_q> = (pi^2/3) * I(L, M, N+2) * triplet cosh combo

    where the u^2 angular factor shifts N by +2 in the master integral
    indexing (which has u^{N+1} already from the Jacobian; adding u^2
    gives u^{N+3} = u^{(N+2)+1}).
    """
    assert bf_p.channel == 'antisym' and bf_q.channel == 'antisym'
    val = _eckart_triplet_combination_pstate(
        bf_p, bf_q, alpha, beta, I_HE_cosh_polynomial,
        L_extra=0, M_extra=0, N_extra=2,
    )
    return P_STATE_VOLUME_FACTOR * val


def potential_vne_element_pstate_eckart_antisym(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float, Z: float,
) -> float:
    """Antisymmetric-channel P-state nuclear attraction.

    <Phi^(-)_p | -Z(1/r_1 + 1/r_2) | Phi^(-)_q>
      = -4 Z (pi^2/3) * J(L+1, M, N+2) * triplet cosh combo
    """
    assert bf_p.channel == 'antisym' and bf_q.channel == 'antisym'
    val = _eckart_triplet_combination_pstate(
        bf_p, bf_q, alpha, beta, J_HE_cosh_polynomial,
        L_extra=1, M_extra=0, N_extra=2,
    )
    return -4.0 * Z * P_STATE_VOLUME_FACTOR * val


def potential_vee_element_pstate_eckart_antisym(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Antisymmetric-channel P-state electron-electron repulsion.

    <Phi^(-)_p | 1/r_12 | Phi^(-)_q>
      = (pi^2/3) * K(L, M, N+2) * triplet cosh combo
    """
    assert bf_p.channel == 'antisym' and bf_q.channel == 'antisym'
    val = _eckart_triplet_combination_pstate(
        bf_p, bf_q, alpha, beta, K_HE_cosh_polynomial,
        L_extra=0, M_extra=0, N_extra=2,
    )
    return P_STATE_VOLUME_FACTOR * val


# ---------------------------------------------------------------------------
# Cross-sector matrix elements (sym × antisym and antisym × sym)
# ---------------------------------------------------------------------------
#
# For phi_p^(+) phi_q^(-) and phi_p^(-) phi_q^(+) cross products:
#   (z_1 + z_2)(z_1 - z_2) = z_1^2 - z_2^2
#   <z_1^2 - z_2^2>_SO(3) = (r_1^2 - r_2^2)/3 = s t / 3
#
# Hyperbolic content (shared beta):
#   cosh(beta_p t) sinh(beta_q t) = (1/2)[sinh(B_+ t) - sinh(B_- t)]
#                                 = (1/2) sinh(2 beta t)   (for shared beta)
#
# So cross-sector matrix elements have:
#   <Phi^(+)_p | O | Phi^(-)_q>  ~  s * t * sinh(2 beta t) * Q_p Q_q * O
#
# The t-power in the integrand becomes ODD (t^{2M+1}), and the hyperbolic
# kernel is sinh — so cross-sector matrix elements use master_S_gen.
#
# At beta = 0, sinh(0) = 0 and all cross-sector matrix elements vanish.

def _cross_sector_sinh_combination(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    L_arg: int, u_pow: int, M_arg: int,
) -> float:
    """Evaluate master_S_gen(L_arg, u_pow, M_arg; alpha, 2 beta).

    The cross-sector hyperbolic content cosh(beta_p t) sinh(beta_q t)
    equals (1/2) sinh(2 beta t) for shared beta, which combines with the
    1/2 from the angular average (s * t / 3 -> s * t / 6 prefactor)
    elsewhere. This helper just returns the bare sinh master.
    """
    if L_arg < 0 or u_pow < 0 or M_arg < 0:
        raise ValueError(
            f"Negative master-integral index: ({L_arg}, {u_pow}, {M_arg})"
        )
    cf = master_S_gen(L_arg, u_pow, M_arg)
    return cf.evaluate(alpha, 2.0 * beta)


# Cross-sector "PI/6" volume factor: pi^2 / 3 (P-state factor) * 1/2 (sinh
# combination at shared beta) = pi^2 / 6.
P_STATE_CROSS_VOLUME_FACTOR = math.pi ** 2 / 6.0


def overlap_element_pstate_eckart_cross(
    bf_p_plus: HylleraasPStateBasisFn,
    bf_q_minus: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Cross-sector overlap <Phi^(+)_p | Phi^(-)_q>.

    The integrand has factor s * t * sinh(2 beta t) from the angular
    average + hyperbolic combination, multiplied by the standard
    u^{N+1} t^{2M} (s^2 - t^2) volume content.

    Expanding (s^2 - t^2):
      s^2 piece: master_S_gen(L+3, N+1, M; alpha, 2 beta)
      -t^2 piece: -master_S_gen(L+1, N+1, M+1; alpha, 2 beta)
        (where t^2 * t^{2M+1} = t^{2(M+1)+1})

    The extra L+1 shift comes from the "s" factor of the s*t angular term.

    At beta = 0, master_S_gen(.., .., ..; 0, 0) = 0 (sinh master is odd in B
    and B=0), so cross-sector overlap is identically zero at beta = 0.
    """
    assert bf_p_plus.channel == 'sym' and bf_q_minus.channel == 'antisym'
    L = bf_p_plus.l + bf_q_minus.l
    M = bf_p_plus.m + bf_q_minus.m
    N = bf_p_plus.n + bf_q_minus.n

    S_s2 = _cross_sector_sinh_combination(
        bf_p_plus, bf_q_minus, alpha, beta,
        L_arg=L + 3, u_pow=N + 1, M_arg=M,
    )
    S_t2 = _cross_sector_sinh_combination(
        bf_p_plus, bf_q_minus, alpha, beta,
        L_arg=L + 1, u_pow=N + 1, M_arg=M + 1,
    )
    return P_STATE_CROSS_VOLUME_FACTOR * (S_s2 - S_t2)


def potential_vne_element_pstate_eckart_cross(
    bf_p_plus: HylleraasPStateBasisFn,
    bf_q_minus: HylleraasPStateBasisFn,
    alpha: float, beta: float, Z: float,
) -> float:
    """Cross-sector nuclear attraction <Phi^(+)_p | -Z(1/r_1+1/r_2) | Phi^(-)_q>.

    The 1/r_1 + 1/r_2 = 4 s / (s^2-t^2) cancels the (s^2-t^2) volume factor.
    The integrand is then -4Z * s^2 * t * sinh(2 beta t) * Q_p Q_q * u^{N+1}.

    Result: -4Z * (pi^2/6) * master_S_gen(L+2, N+1, M; alpha, 2 beta).

    (The extra L+2 = L+1 (from V_ne overall s) + L+1 (from s of s*t) = L+2.)
    """
    assert bf_p_plus.channel == 'sym' and bf_q_minus.channel == 'antisym'
    L = bf_p_plus.l + bf_q_minus.l
    M = bf_p_plus.m + bf_q_minus.m
    N = bf_p_plus.n + bf_q_minus.n

    val = _cross_sector_sinh_combination(
        bf_p_plus, bf_q_minus, alpha, beta,
        L_arg=L + 2, u_pow=N + 1, M_arg=M,
    )
    return -4.0 * Z * P_STATE_CROSS_VOLUME_FACTOR * val


def potential_vee_element_pstate_eckart_cross(
    bf_p_plus: HylleraasPStateBasisFn,
    bf_q_minus: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Cross-sector electron-electron repulsion <Phi^(+)_p | 1/r_12 | Phi^(-)_q>.

    1/r_12 = 1/u cancels one u from the volume's Jacobian. The integrand is
    then s * t * sinh(2 beta t) * Q_p Q_q * (s^2-t^2) (no extra u).

    Result: (pi^2/6) * [master_S_gen(L+3, N, M; 2 beta) - master_S_gen(L+1, N, M+1; 2 beta)]
    """
    assert bf_p_plus.channel == 'sym' and bf_q_minus.channel == 'antisym'
    L = bf_p_plus.l + bf_q_minus.l
    M = bf_p_plus.m + bf_q_minus.m
    N = bf_p_plus.n + bf_q_minus.n

    V_s2 = _cross_sector_sinh_combination(
        bf_p_plus, bf_q_minus, alpha, beta,
        L_arg=L + 3, u_pow=N, M_arg=M,
    )
    V_t2 = _cross_sector_sinh_combination(
        bf_p_plus, bf_q_minus, alpha, beta,
        L_arg=L + 1, u_pow=N, M_arg=M + 1,
    )
    return P_STATE_CROSS_VOLUME_FACTOR * (V_s2 - V_t2)


# ---------------------------------------------------------------------------
# P-state kinetic energy (sym × sym channel, singlet)
# ---------------------------------------------------------------------------
#
# Hartree-form decomposition <Phi_p | T | Phi_q> = T_1 + T_mid_q + T_mid_p + T_3
#
# where:
#   T_1     = (s^2 + t^2 - u^2)/3 multiplied INTO the S-state kinetic integrand
#             -- apply 3 index shifts (s^2, t^2, -u^2) to each S-state kinetic
#             contribution and multiply by 1/3
#   T_mid_q = (pi^2/3) * integral chi_p * u * [s(s^2-u^2)*chi_{q,s} +
#                                                t(u^2-t^2)*chi_{q,t}] dV
#             (the (s^2-t^2) of the SO(3)-averaged middle term cancels with
#             the (s^2-t^2) of the volume factor, leaving the clean integrand)
#   T_mid_p = same as T_mid_q with (p, q) swapped
#   T_3     = (pi^2) * integral chi_p chi_q u(s^2-t^2) dV  ==  S-state overlap
#
# Validated at (l=m=n=0, beta=0) against the analytical reference 2 pi^2/alpha^6.

def _apply_pstate_angular_shifts(base_contribs):
    """For each S-state kinetic contribution (coeff, L, u_pow, M, kind, B),
    produce three P-state contributions corresponding to multiplying by
    (s^2 + t^2 - u^2)/3:

      s^2:  shift L by +2
      t^2:  shift M by +1 (since t^2 * t^{2M} = t^{2(M+1)})
      -u^2: shift u_pow by +2, sign flip
    """
    out = []
    for (c, L, u_pow, M, kind, B) in base_contribs:
        c_over_3 = c / 3.0
        out.append((c_over_3, L + 2, u_pow, M, kind, B))     # s^2
        out.append((c_over_3, L, u_pow, M + 1, kind, B))     # t^2
        out.append((-c_over_3, L, u_pow + 2, M, kind, B))    # -u^2
    return out


def _ke_pstate_mid_contributions(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    spin: str = 'singlet',
):
    """Contributions to T_mid_q = (pi^2/3) integral chi_p [s(s^2-u^2) chi_{q,s}
                                                            + t(u^2-t^2) chi_{q,t}] u dV.

    Returns list of (coeff, L, u_pow, M, kind, B) tuples. The (pi^2/3)
    factor is NOT included; multiply the summed tuples by pi^2/3 outside.

    For singlet (cosh basis):
      chi = exp(-alpha s) cosh(beta t) Q   with Q = s^l t^{2m} u^n
      chi_s = chi * (-alpha + l/s) = exp(-alpha s) cosh(beta t) (l s^{l-1} - alpha s^l) t^{2m} u^n
      chi_t = exp(-alpha s) [beta sinh(beta t) Q + cosh(beta t) Q_t]
            (Q_t = 2 m t^{2m-1} s^l u^n)
      chi_u = chi * n/u

    The hyperbolic products for shared beta (singlet) are:
      F_p F_q     = cosh^2(beta t) = (1 + cosh(2 beta t))/2  --> FF_pieces
      F_p F'_q    = (beta / 2) sinh(2 beta t)                 --> FpF_pieces
                  (F' = derivative of F; here we use F=cosh, F'=beta sinh)
    """
    if spin != 'singlet':
        raise NotImplementedError(
            f"Triplet P-state kinetic not yet implemented; got spin={spin!r}"
        )

    l_p, m_p, n_p = bf_p.l, bf_p.m, bf_p.n
    l_q, m_q, n_q = bf_q.l, bf_q.m, bf_q.n
    L = l_p + l_q
    M = m_p + m_q
    N = n_p + n_q
    two_beta = 2.0 * beta

    # Hyperbolic pieces.
    FF_pieces = [(0.5, 'C', 0.0), (0.5, 'C', two_beta)]
    # F_p F'_q = (beta/2) sinh(2 beta t)
    FpF_pieces = [(0.5 * beta, 'S', two_beta)]

    out = []

    def add(coeff, L_, u_pow, M_, kind, B):
        if M_ < 0 or L_ < 0 or u_pow < 0:
            return
        if coeff == 0.0:
            return
        out.append((coeff, L_, u_pow, M_, kind, B))

    # =====================================================================
    # Piece A: u s(s^2 - u^2) * chi_p * chi_{q,s}
    #
    # chi_p chi_{q,s} = cosh^2 * Q_p * (-alpha Q_q + Q_{q,s})  * e^{-2 alpha s}
    #
    # Sub-piece A1: (-alpha Q_p Q_q) part. Q_p Q_q = s^L t^{2M} u^N.
    #   Times u s(s^2-u^2) = u (s^3 - s u^2):
    #     u^{N+1} t^{2M} [s^{L+3} - s^{L+1} u^2]
    #     = (u^{N+1} s^{L+3} - u^{N+3} s^{L+1}) t^{2M}
    #   Contribution = (-alpha) * sum_{FF} h_coeff * [C(L+3, N+1, M; h_B) - C(L+1, N+3, M; h_B)]
    #
    # Sub-piece A2: (Q_p Q_{q,s} = l_q s^{L-1} t^{2M} u^N) -- only when l_q >= 1.
    #   Times u s(s^2-u^2) gives:
    #     l_q * [u^{N+1} s^{L+2} - u^{N+3} s^L] t^{2M}
    #   Contribution = l_q * sum_{FF} h_coeff * [C(L+2, N+1, M; h_B) - C(L, N+3, M; h_B)]
    # =====================================================================
    for (h_c, h_kind, h_B) in FF_pieces:
        # A1
        add(-alpha * h_c, L + 3, N + 1, M, h_kind, h_B)
        add(+alpha * h_c, L + 1, N + 3, M, h_kind, h_B)
        # A2 (when l_q > 0)
        if l_q > 0:
            add(l_q * h_c, L + 2, N + 1, M, h_kind, h_B)
            add(-l_q * h_c, L, N + 3, M, h_kind, h_B)

    # =====================================================================
    # Piece B: u t(u^2 - t^2) * chi_p * chi_{q,t}
    #
    # chi_p chi_{q,t} = e^{-2 alpha s} * [F_p F'_q * Q_p Q_q + F_p F_q * Q_p Q_{q,t}]
    #
    # Sub-piece B1: F_p F'_q * Q_p Q_q
    #   F_p F'_q = (beta/2) sinh(2 beta t)
    #   t * (u^2 - t^2) * t^{2M} = (u^2 t - t^3) t^{2M} = u^2 t^{2M+1} - t^{2M+3}
    #   Times u^{N+1} = u^{N+3} t^{2M+1} - u^{N+1} t^{2M+3}
    #   Contribution = (beta/2) * [S(L, N+3, M; 2 beta) - S(L, N+1, M+1; 2 beta)]
    #   (S = sinh master with t^{2M+1}; t^{2M+3} = t^{2(M+1)+1} -> M_arg = M+1)
    #
    # Sub-piece B2: F_p F_q * Q_p Q_{q,t} -- only when m_q >= 1
    #   Q_{q,t} = 2 m_q t^{2m_q - 1} s^{l_q} u^{n_q}
    #   Q_p Q_{q,t} = 2 m_q s^L t^{2M-1} u^N
    #   t * Q_p Q_{q,t} = 2 m_q s^L t^{2M} u^N (t * t^{2M-1} = t^{2M})
    #   Combined with u(u^2-t^2): u^{N+1}(u^2-t^2) gives (u^{N+3} - u^{N+1} t^2) terms
    #     so integrand = 2 m_q s^L t^{2M} [u^{N+3} - u^{N+1} t^2]
    #                  = 2 m_q [s^L u^{N+3} t^{2M} - s^L u^{N+1} t^{2(M+1)}]
    #   Contribution = 2 m_q * sum_{FF} h_coeff * [C(L, N+3, M; h_B) - C(L, N+1, M+1; h_B)]
    # =====================================================================
    # B1 (always present at beta > 0)
    for (h_c, h_kind, h_B) in FpF_pieces:
        add(h_c, L, N + 3, M, h_kind, h_B)
        add(-h_c, L, N + 1, M + 1, h_kind, h_B)
    # B2 (only when m_q > 0)
    if m_q > 0:
        for (h_c, h_kind, h_B) in FF_pieces:
            add(2 * m_q * h_c, L, N + 3, M, h_kind, h_B)
            add(-2 * m_q * h_c, L, N + 1, M + 1, h_kind, h_B)

    return out


def _evaluate_pstate_contributions(contribs, alpha):
    """Sum coeff * master_X(L, u_pow, M).evaluate(alpha, B) over a contributions list."""
    total = 0.0
    for (c, L, u_pow, M, kind, B) in contribs:
        if kind == 'C':
            cf = master_C_gen(L, u_pow, M)
        elif kind == 'S':
            cf = master_S_gen(L, u_pow, M)
        else:
            raise ValueError(f"unknown kind {kind!r}")
        total += c * cf.evaluate(alpha, B)
    return total


def kinetic_element_pstate_eckart_sym_sym(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
) -> float:
    """Symmetric-channel P-state kinetic energy <Phi^(+)_p | T | Phi^(+)_q>.

    Algebraic implementation via Hartree decomposition T = T_1 + T_mid_q
    + T_mid_p + T_3:

      T_1     = (pi^2) * sum_{shifted contribs} ; angular factor
                (s^2+t^2-u^2)/3 multiplied into S-state kinetic
      T_mid_q = (pi^2/3) * sum over chi_p * (derivative terms of chi_q)
      T_mid_p = (pi^2/3) * sum over chi_q * (derivative terms of chi_p)
      T_3     = (pi^2) * S-state overlap of chi_p, chi_q   [no angular factor]

    At (l=m=n=0, beta=0), returns the analytical value 2 pi^2 / alpha^6.
    """
    assert bf_p.channel == 'sym' and bf_q.channel == 'sym'

    # S-state HylleraasBasisFn equivalents.
    s_bf_p = HylleraasBasisFn(bf_p.l, bf_p.m, bf_p.n)
    s_bf_q = HylleraasBasisFn(bf_q.l, bf_q.m, bf_q.n)

    # T_1: (s^2+t^2-u^2)/3 multiplied into S-state kinetic.
    s_contribs = _ke_term_contributions(s_bf_p, s_bf_q, alpha, beta, spin='singlet')
    pstate_T1_contribs = _apply_pstate_angular_shifts(s_contribs)
    T_1 = math.pi ** 2 * _evaluate_pstate_contributions(pstate_T1_contribs, alpha)

    # T_mid_q: chi_p * M(chi_q)
    mid_q_contribs = _ke_pstate_mid_contributions(bf_p, bf_q, alpha, beta, spin='singlet')
    T_mid_q = (math.pi ** 2 / 3.0) * _evaluate_pstate_contributions(mid_q_contribs, alpha)

    # T_mid_p: chi_q * M(chi_p) — same as T_mid_q with (p, q) swapped.
    mid_p_contribs = _ke_pstate_mid_contributions(bf_q, bf_p, alpha, beta, spin='singlet')
    T_mid_p = (math.pi ** 2 / 3.0) * _evaluate_pstate_contributions(mid_p_contribs, alpha)

    # T_3: S-state overlap of chi_p, chi_q.
    T_3 = overlap_element_eckart(s_bf_p, s_bf_q, alpha, beta, spin='singlet')

    return T_1 + T_mid_q + T_mid_p + T_3


# ---------------------------------------------------------------------------
# Universal P-state kinetic via 3D quadrature with analytical SO(3) reduction
# ---------------------------------------------------------------------------
#
# Handles ALL four channel combinations (sym/antisym × sym/antisym) by
# evaluating the SO(3)-averaged kinetic density at each (r_1, r_2, cos θ_12)
# quadrature point.
#
# Derivation (see notes 2026-05-18). For Φ_p^{(a)} = X^{(a)} χ_p^{(a)} with
# X^{(+)} = z_1+z_2, X^{(-)} = z_1-z_2:
#
#   ⟨Φ_p^{(a)}|T|Φ_q^{(b)}⟩ = (1/2) ∫ {T_1 + mid_p + mid_q + T_3}_{SO(3)} dV
#
# where T_1 = ⟨X_p X_q⟩·∑_i ∇_iχ_p·∇_iχ_q is the "no-X-gradient" piece;
# mid_p is the cross term from gradient on χ_p; mid_q likewise on χ_q;
# T_3 = ⟨X·X⟩-equivalent factor on χ_p χ_q (kills cross-sector by
# (−1+1) cancellation).
#
# Channel-resolved SO(3) averages:
#   ⟨(z_1+z_2)²⟩ = (2r_1²+2r_2²-r_12²)/3
#   ⟨(z_1-z_2)²⟩ = r_12²/3
#   ⟨(z_1+z_2)(z_1-z_2)⟩ = (r_1²-r_2²)/3
#
# Mid-piece SO(3) averages (4 types):
#   ⟨(z_1+z_2)\hat z·(∇_1+∇_2)χ⟩ =
#     (3r_1²+r_2²-r_12²)/(6r_1)·χ_{r_1} + (r_1²+3r_2²-r_12²)/(6r_2)·χ_{r_2}
#   ⟨(z_1+z_2)\hat z·(∇_1-∇_2)χ⟩ =
#     (3r_1²+r_2²-r_12²)/(6r_1)·χ_{r_1} - (r_1²+3r_2²-r_12²)/(6r_2)·χ_{r_2}
#     + 2(r_1²-r_2²)/(3 r_12)·χ_{r_12}
#   ⟨(z_1-z_2)\hat z·(∇_1+∇_2)χ⟩ =
#     (r_1²-r_2²+r_12²)/(6r_1)·χ_{r_1} + (r_1²-r_2²-r_12²)/(6r_2)·χ_{r_2}
#   ⟨(z_1-z_2)\hat z·(∇_1-∇_2)χ⟩ =
#     (r_1²-r_2²+r_12²)/(6r_1)·χ_{r_1} - (r_1²-r_2²-r_12²)/(6r_2)·χ_{r_2}
#     + (2 r_12/3)·χ_{r_12}
#
# The choice of gradient combination ∇_1±∇_2 for mid_p is determined by the
# OTHER basis function's channel (b for mid_p; a for mid_q).


def _kinetic_via_quadrature_pstate(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    n_r: int = 32, n_theta: int = 16,
) -> float:
    """Universal P-state kinetic energy matrix element via 3D quadrature.

    Handles all four channel combinations: sym×sym, antisym×antisym,
    sym×antisym (cross), antisym×sym (cross).

    At (l=m=n=0, beta=0) with channel='sym' both, returns the analytic
    reference 2 π² / α⁶ (validated 2026-05-18). For antisym channels at
    β=0, returns 0 identically (basis vanishes).
    """
    from scipy.special import roots_genlaguerre, roots_legendre

    nodes_r, weights_r = roots_genlaguerre(n_r, 0.0)
    nodes_c, weights_c = roots_legendre(n_theta)
    two_alpha = 2.0 * alpha

    chan_p = bf_p.channel
    chan_q = bf_q.channel

    # T_3 coefficient: 2 for same channel, 0 for cross
    t3_coeff = 2.0 if chan_p == chan_q else 0.0

    # Gradient combination for mid pieces — based on the OTHER basis fn's channel
    # mid_p uses combination based on chan_q; mid_q uses combination based on chan_p
    eps_p_for_mid_q = +1.0 if chan_p == 'sym' else -1.0  # ε_2 for mid_q
    eps_q_for_mid_p = +1.0 if chan_q == 'sym' else -1.0  # ε_2 for mid_p
    # X kind for each mid piece
    X_kind_p = '+' if chan_p == 'sym' else '-'  # X factor in mid_p
    X_kind_q = '+' if chan_q == 'sym' else '-'  # X factor in mid_q

    total = 0.0

    for ir1 in range(n_r):
        x1 = nodes_r[ir1]
        wr1 = weights_r[ir1]
        r1 = x1 / two_alpha

        for ir2 in range(n_r):
            x2 = nodes_r[ir2]
            wr2 = weights_r[ir2]
            r2 = x2 / two_alpha

            s = r1 + r2
            t = r1 - r2

            for ic in range(n_theta):
                c = nodes_c[ic]
                wc = weights_c[ic]

                r12_sq = max(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * c, 1e-30)
                r12 = math.sqrt(r12_sq)

                exp_factor = math.exp(-alpha * s)

                # Hyperbolic factors and their t-derivatives
                if chan_p == 'sym':
                    bt_p = math.cosh(beta * t)
                    bt_p_dt = beta * math.sinh(beta * t)
                else:
                    bt_p = math.sinh(beta * t)
                    bt_p_dt = beta * math.cosh(beta * t)
                if chan_q == 'sym':
                    bt_q = math.cosh(beta * t)
                    bt_q_dt = beta * math.sinh(beta * t)
                else:
                    bt_q = math.sinh(beta * t)
                    bt_q_dt = beta * math.cosh(beta * t)

                # Polynomial parts Q = s^l t^{2m} u^n
                Q_p = (s ** bf_p.l) \
                    * (t ** (2 * bf_p.m) if bf_p.m > 0 else 1.0) \
                    * (r12 ** bf_p.n if bf_p.n > 0 else 1.0)
                Q_q = (s ** bf_q.l) \
                    * (t ** (2 * bf_q.m) if bf_q.m > 0 else 1.0) \
                    * (r12 ** bf_q.n if bf_q.n > 0 else 1.0)

                phi_p = exp_factor * bt_p * Q_p
                phi_q = exp_factor * bt_q * Q_q

                # ∂χ/∂s = phi · (-α + l/s)
                if s > 0:
                    chi_p_s = phi_p * (-alpha + (bf_p.l / s if bf_p.l > 0 else 0.0))
                    chi_q_s = phi_q * (-alpha + (bf_q.l / s if bf_q.l > 0 else 0.0))
                else:
                    chi_p_s = -alpha * phi_p
                    chi_q_s = -alpha * phi_q

                # ∂χ/∂t — from bt_factor derivative and Q_t
                if bf_p.m > 0 and abs(t) > 1e-15:
                    chi_p_t_Q = phi_p * (2 * bf_p.m / t)
                else:
                    chi_p_t_Q = 0.0
                chi_p_t_bt = exp_factor * Q_p * bt_p_dt
                chi_p_t = chi_p_t_Q + chi_p_t_bt

                if bf_q.m > 0 and abs(t) > 1e-15:
                    chi_q_t_Q = phi_q * (2 * bf_q.m / t)
                else:
                    chi_q_t_Q = 0.0
                chi_q_t_bt = exp_factor * Q_q * bt_q_dt
                chi_q_t = chi_q_t_Q + chi_q_t_bt

                # ∂χ/∂u
                if bf_p.n > 0 and r12 > 1e-15:
                    chi_p_u = phi_p * (bf_p.n / r12)
                else:
                    chi_p_u = 0.0
                if bf_q.n > 0 and r12 > 1e-15:
                    chi_q_u = phi_q * (bf_q.n / r12)
                else:
                    chi_q_u = 0.0

                chi_p_r1 = chi_p_s + chi_p_t
                chi_p_r2 = chi_p_s - chi_p_t
                chi_p_r12 = chi_p_u
                chi_q_r1 = chi_q_s + chi_q_t
                chi_q_r2 = chi_q_s - chi_q_t
                chi_q_r12 = chi_q_u

                # cos_1, cos_2 for ∇_i·∇_i' inner products
                if r1 > 1e-15 and r12 > 1e-15:
                    cos_a = (r1 ** 2 - r2 ** 2 + r12_sq) / (2 * r1 * r12)
                else:
                    cos_a = 0.0
                if r2 > 1e-15 and r12 > 1e-15:
                    cos_b = (r2 ** 2 - r1 ** 2 + r12_sq) / (2 * r2 * r12)
                else:
                    cos_b = 0.0

                # Inner products ∇_i χ_p · ∇_i χ_q (SO(3) scalars)
                grad1_dot = (
                    chi_p_r1 * chi_q_r1
                    + chi_p_r12 * chi_q_r12
                    + cos_a * (chi_p_r1 * chi_q_r12 + chi_p_r12 * chi_q_r1)
                )
                grad2_dot = (
                    chi_p_r2 * chi_q_r2
                    + chi_p_r12 * chi_q_r12
                    + cos_b * (chi_p_r2 * chi_q_r12 + chi_p_r12 * chi_q_r2)
                )

                # T_1: ⟨X_p X_q⟩ · sum of inner products
                if chan_p == 'sym' and chan_q == 'sym':
                    X_pq_avg = (2 * r1 ** 2 + 2 * r2 ** 2 - r12_sq) / 3.0
                elif chan_p == 'antisym' and chan_q == 'antisym':
                    X_pq_avg = r12_sq / 3.0
                else:
                    X_pq_avg = (r1 ** 2 - r2 ** 2) / 3.0

                T_1 = X_pq_avg * (grad1_dot + grad2_dot)

                # Mid-piece SO(3) averages (4 types based on (X_kind, grad_kind))
                inv_r1 = (1.0 / r1) if r1 > 1e-15 else 0.0
                inv_r2 = (1.0 / r2) if r2 > 1e-15 else 0.0
                inv_r12 = (1.0 / r12) if r12 > 1e-15 else 0.0

                def mid_avg_coeffs(X_kind: str, grad_kind: str):
                    r"""Return (c_r1, c_r2, c_r12) for the SO(3)-averaged
                    ⟨X^{(X_kind)} \hat z · (∇_1 + grad_kind*∇_2) χ⟩."""
                    if X_kind == '+':
                        c_r1 = (3 * r1 ** 2 + r2 ** 2 - r12_sq) / 6 * inv_r1
                        if grad_kind == 'plus':
                            c_r2 = (r1 ** 2 + 3 * r2 ** 2 - r12_sq) / 6 * inv_r2
                            c_r12 = 0.0
                        else:  # 'minus'
                            c_r2 = -(r1 ** 2 + 3 * r2 ** 2 - r12_sq) / 6 * inv_r2
                            c_r12 = 2 * (r1 ** 2 - r2 ** 2) / 3 * inv_r12
                    else:  # X_kind == '-'
                        c_r1 = (r1 ** 2 - r2 ** 2 + r12_sq) / 6 * inv_r1
                        if grad_kind == 'plus':
                            c_r2 = (r1 ** 2 - r2 ** 2 - r12_sq) / 6 * inv_r2
                            c_r12 = 0.0
                        else:
                            c_r2 = -(r1 ** 2 - r2 ** 2 - r12_sq) / 6 * inv_r2
                            c_r12 = 2 * r12 / 3
                    return c_r1, c_r2, c_r12

                # Mid_p: X factor of p, gradient combination determined by q's channel
                grad_kind_q = 'plus' if chan_q == 'sym' else 'minus'
                cp_r1, cp_r2, cp_r12 = mid_avg_coeffs(X_kind_p, grad_kind_q)
                mid_p = (
                    cp_r1 * chi_p_r1 + cp_r2 * chi_p_r2 + cp_r12 * chi_p_r12
                ) * phi_q

                # Mid_q: X factor of q, gradient combination determined by p's channel
                grad_kind_p = 'plus' if chan_p == 'sym' else 'minus'
                cq_r1, cq_r2, cq_r12 = mid_avg_coeffs(X_kind_q, grad_kind_p)
                mid_q = (
                    cq_r1 * chi_q_r1 + cq_r2 * chi_q_r2 + cq_r12 * chi_q_r12
                ) * phi_p

                # T_3: scalar product of χ_p χ_q with the constant ∑ε_p^i ε_q^i
                T_3 = t3_coeff * phi_p * phi_q

                # Total kinetic density (Hartree form)
                kinetic_density = 0.5 * (T_1 + mid_p + mid_q + T_3)

                # Strip exp(-2αs)
                if exp_factor > 0:
                    k_stripped = kinetic_density / (exp_factor * exp_factor)
                else:
                    k_stripped = 0.0

                vol = 8.0 * math.pi ** 2 * r1 ** 2 * r2 ** 2
                jacobian = 1.0 / (two_alpha * two_alpha)

                total += k_stripped * wr1 * wr2 * wc * vol * jacobian

    return total


def kinetic_element_pstate(
    bf_p: HylleraasPStateBasisFn,
    bf_q: HylleraasPStateBasisFn,
    alpha: float, beta: float,
    n_r: int = 32, n_theta: int = 16,
) -> float:
    """Universal P-state kinetic dispatch.

    For sym×sym, routes to the fast algebraic implementation
    (kinetic_element_pstate_eckart_sym_sym). For other channel
    combinations, falls back to quadrature.
    """
    if bf_p.channel == 'sym' and bf_q.channel == 'sym':
        return kinetic_element_pstate_eckart_sym_sym(bf_p, bf_q, alpha, beta)
    return _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, n_r, n_theta)


# ---------------------------------------------------------------------------
# Sym-only 2^1P variational driver
# ---------------------------------------------------------------------------
#
# A single-channel (symmetric only) 2^1P trial wavefunction:
#   Psi_{2^1P} = sum_q c_q (z_1+z_2) e^{-alpha s} cosh(beta t) s^{l_q} t^{2m_q} u^{n_q}
#
# This is an approximation that drops the antisymmetric channel
# (z_1-z_2) sinh(beta t) Q. The variational principle gives an upper
# bound; the antisymmetric channel would lower it.
#
# Drake handbook: E(2^1P) = -2.123843 Ha; we expect to come in within a
# few mHa with omega = 3 or 4.


def assemble_p_state_matrices_sym(
    basis: List[HylleraasPStateBasisFn],
    alpha: float, Z: float, beta: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build (H, S) matrices for the sym-only P-state trial."""
    n = len(basis)
    S = np.zeros((n, n))
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            S[i, j] = overlap_element_pstate_eckart(basis[i], basis[j], alpha, beta)
            T_ij = kinetic_element_pstate_eckart_sym_sym(
                basis[i], basis[j], alpha, beta
            )
            Vne_ij = potential_vne_element_pstate_eckart(
                basis[i], basis[j], alpha, beta, Z=Z
            )
            Vee_ij = potential_vee_element_pstate_eckart(
                basis[i], basis[j], alpha, beta
            )
            H[i, j] = T_ij + Vne_ij + Vee_ij
            S[j, i] = S[i, j]
            H[j, i] = H[i, j]
    H = 0.5 * (H + H.T)
    S = 0.5 * (S + S.T)
    return H, S


def solve_2p1_state_sym(
    basis: List[HylleraasPStateBasisFn],
    alpha: float, Z: float, beta: float,
) -> Tuple[float, np.ndarray]:
    """Solve generalized eigenvalue problem for sym-only 2^1P.
    Returns (E(2^1P), eigenvector c)."""
    from scipy.linalg import eigh
    H, S = assemble_p_state_matrices_sym(basis, alpha, Z, beta)
    eigvals, eigvecs = eigh(H, b=S)
    return float(eigvals[0]), eigvecs[:, 0]


def optimize_2p1_sym(
    basis: List[HylleraasPStateBasisFn],
    Z: float,
    alpha_init: float = 1.35,
    beta_init: float = 0.35,
    method: str = 'Nelder-Mead',
    tol: float = 1e-7,
    max_iter: int = 200,
) -> Dict[str, Any]:
    """2D variational optimization over (alpha, beta) for sym-only 2^1P.

    Returns dict with E_2p, alpha_opt, beta_opt, coeffs, basis, opt_n_iter.
    """
    from scipy.optimize import minimize

    def obj(x):
        a, b = float(x[0]), float(x[1])
        if a <= 0 or abs(b) >= a:
            return 1e10
        try:
            E, _c = solve_2p1_state_sym(basis, a, Z, b)
            return E
        except Exception:
            return 1e10

    x0 = np.array([alpha_init, beta_init])
    result = minimize(
        obj, x0=x0, method=method,
        options={'xatol': tol, 'fatol': tol, 'maxiter': max_iter, 'adaptive': True},
    )
    alpha_opt, beta_opt = float(result.x[0]), float(result.x[1])
    E_2p, coeffs = solve_2p1_state_sym(basis, alpha_opt, Z, beta_opt)
    return {
        'E_2p': E_2p,
        'alpha_opt': alpha_opt,
        'beta_opt': beta_opt,
        'coeffs': coeffs,
        'basis': list(basis),
        'opt_n_iter': int(result.nit),
        'opt_n_eval': int(result.nfev),
        'opt_success': bool(result.success),
    }


# ---------------------------------------------------------------------------
# Cross-basis dipole matrix element <1^1S | z_1 + z_2 | 2^1P>
# ---------------------------------------------------------------------------
#
# The dipole operator z_1 + z_2 turns an S-state basis function
# phi^S = e^{-alpha_S s} cosh(beta_S t) Q_S into a (z_1+z_2)-multiplied
# function, which has the same angular structure as the P-state basis.
#
# Matrix element between an S-state phi^S_p and a P-state phi^{P+}_q
# (sym channel):
#   <phi^S_p | (z_1+z_2) | phi^{P+}_q>
#     = int (z_1+z_2)^2 chi^S_p · chi^P_q · d^6 r
#
# After SO(3) averaging:
#   = (pi^2/3) int (s^2 + t^2 - u^2) chi^S_p chi^P_q u(s^2-t^2) ds dt du
#
# This is the same form as the P-state overlap, but with TWO DIFFERENT
# (alpha, beta) pairs. We evaluate by combining:
#   alpha_eff = (alpha_S + alpha_P) / 2
#   B_+ = beta_S + beta_P
#   B_- = beta_S - beta_P
# and using the singlet cosh-cosh combo (1/2)[I_at(B_+) + I_at(B_-)]
# with the three (L+2, M+1, N+2) shifts as in P-state overlap.

def dipole_element_1s_2p_sym(
    bf_s: HylleraasBasisFn,
    bf_p: HylleraasPStateBasisFn,
    alpha_s: float, beta_s: float,
    alpha_p: float, beta_p: float,
) -> float:
    """Cross-basis dipole matrix element <phi^S | z_1+z_2 | phi^{P,sym}>.

    Parameters
    ----------
    bf_s : HylleraasBasisFn
        S-state basis function (l_s, m_s, n_s).
    bf_p : HylleraasPStateBasisFn
        P-state symmetric-channel basis function (l_p, m_p, n_p).
    alpha_s, beta_s : float
        S-state nonlinear parameters.
    alpha_p, beta_p : float
        P-state nonlinear parameters.

    Returns
    -------
    float : the dipole matrix element.
    """
    assert bf_p.channel == 'sym'
    L = bf_s.l + bf_p.l
    M = bf_s.m + bf_p.m
    N = bf_s.n + bf_p.n

    # alpha_eff so that 2*alpha_eff = alpha_s + alpha_p (matching the
    # e^{-(alpha_s + alpha_p) s} = e^{-2 alpha_eff s} factor expected by
    # the I master).
    alpha_eff = 0.5 * (alpha_s + alpha_p)
    B_plus = beta_s + beta_p
    B_minus = beta_s - beta_p

    def cosh_combo(L_extra, M_extra, N_extra):
        cf = I_HE_cosh_polynomial(L + L_extra, M + M_extra, N + N_extra)
        return 0.5 * (cf.evaluate(alpha_eff, B_plus)
                      + cf.evaluate(alpha_eff, B_minus))

    val = cosh_combo(2, 0, 0) + cosh_combo(0, 1, 0) - cosh_combo(0, 0, 2)
    return P_STATE_VOLUME_FACTOR * val


def compute_dipole_1s_to_2p_sym(
    s_state: Dict[str, Any],
    p_state: Dict[str, Any],
) -> float:
    """Compute <Psi_{1^1S} | z_1 + z_2 | Psi_{2^1P, sym}> by summing over
    basis-function contributions.

    Parameters
    ----------
    s_state : dict
        Must have 'basis' (list of HylleraasBasisFn), 'coeffs', 'alpha', 'beta'.
    p_state : dict
        Must have 'basis' (list of HylleraasPStateBasisFn), 'coeffs',
        'alpha_opt', 'beta_opt'.

    Returns
    -------
    float : the scalar dipole matrix element D = <Psi_S | (z_1+z_2) | Psi_P>.
    """
    basis_S = s_state['basis']
    c_S = s_state['coeffs']
    alpha_S = s_state['alpha']
    beta_S = s_state.get('beta', 0.0)

    basis_P = p_state['basis']
    c_P = p_state['coeffs']
    alpha_P = p_state['alpha_opt']
    beta_P = p_state['beta_opt']

    D = 0.0
    for i, bf_s in enumerate(basis_S):
        for j, bf_p in enumerate(basis_P):
            D += c_S[i] * c_P[j] * dipole_element_1s_2p_sym(
                bf_s, bf_p, alpha_S, beta_S, alpha_P, beta_P
            )
    return D


def dipole_element_1s_2p_antisym(
    bf_s: HylleraasBasisFn,
    bf_p: HylleraasPStateBasisFn,
    alpha_s: float, beta_s: float,
    alpha_p: float, beta_p: float,
) -> float:
    """Cross-basis dipole <phi^S | (z_1+z_2) | phi^{P,antisym}>.

    The dipole operator (z_1+z_2) times the antisym basis (z_1-z_2)·χ
    gives (z_1+z_2)(z_1-z_2) = z_1²-z_2², which SO(3)-averages to st/3.

    Integrand: (π²/3) · st · χ^S · χ^{P,antisym} · u(s²-t²) ds dt du
             = (π²/3) [s^{L+3} u^{N+1} t^{2M+1} − s^{L+1} u^{N+1} t^{2(M+1)+1}]·sinh-content

    Hyperbolic content: cosh(β_s t) sinh(β_p t) = (1/2)[sinh(B_+ t) − sinh(B_- t)],
    so we use the master_S generators at B_± with α_eff = (α_s+α_p)/2.

    At β_p = 0 (no antisym basis), this is identically zero.
    """
    assert bf_p.channel == 'antisym'
    L = bf_s.l + bf_p.l
    M = bf_s.m + bf_p.m
    N = bf_s.n + bf_p.n

    alpha_eff = 0.5 * (alpha_s + alpha_p)
    B_plus = beta_s + beta_p
    B_minus = beta_s - beta_p

    def sinh_combo(L_extra, M_extra, N_extra):
        L_arg = L + L_extra
        M_arg = M + M_extra
        N_arg = N + N_extra
        if L_arg < 0 or M_arg < 0 or N_arg < 0:
            return 0.0
        cf = master_S_gen(L_arg, N_arg, M_arg)
        return 0.5 * (cf.evaluate(alpha_eff, B_plus)
                      - cf.evaluate(alpha_eff, B_minus))

    val = sinh_combo(3, 0, 1) - sinh_combo(1, 1, 1)
    # P_STATE_VOLUME_FACTOR = π²/3 (P-state factor 1/3 × π² overlap factor)
    return P_STATE_VOLUME_FACTOR * val


# ---------------------------------------------------------------------------
# Full two-channel 2^1P solver
# ---------------------------------------------------------------------------
#
# The Schwartz 1961 trial wavefunction has BOTH sym and antisym channels:
#   Psi_{2^1P}^{M=0} = Σ_q c^sym_q (z_1+z_2) cosh(βt) Q_q
#                    + Σ_q c^anti_q (z_1-z_2) sinh(βt) Q_q
#
# All basis functions share the same (alpha, beta) nonlinear parameters.
# The H, S matrices block-decompose as
#       [ H^{++}  H^{+-} ]      [ S^{++}  S^{+-} ]
#       [ H^{-+}  H^{--} ]      [ S^{-+}  S^{--} ]
# with cross-block (sym × antisym) entries vanishing identically at β=0.


def assemble_p_state_matrices_full(
    basis_sym: List[HylleraasPStateBasisFn],
    basis_antisym: List[HylleraasPStateBasisFn],
    alpha: float, Z: float, beta: float,
    n_r: int = 32, n_theta: int = 16,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build full (H, S) matrices for the two-channel 2^1P trial.

    The basis ordering is [sym basis, antisym basis], so the returned
    matrices have block structure [[++,+-], [-+,--]].
    """
    n_sym = len(basis_sym)
    n_anti = len(basis_antisym)
    n = n_sym + n_anti
    H = np.zeros((n, n))
    S = np.zeros((n, n))

    # === Sym × Sym block ===
    for i in range(n_sym):
        for j in range(i, n_sym):
            bp, bq = basis_sym[i], basis_sym[j]
            S[i, j] = overlap_element_pstate_eckart(bp, bq, alpha, beta)
            T_ij = kinetic_element_pstate_eckart_sym_sym(bp, bq, alpha, beta)
            Vne_ij = potential_vne_element_pstate_eckart(bp, bq, alpha, beta, Z=Z)
            Vee_ij = potential_vee_element_pstate_eckart(bp, bq, alpha, beta)
            H[i, j] = T_ij + Vne_ij + Vee_ij
            S[j, i] = S[i, j]
            H[j, i] = H[i, j]

    # === Antisym × Antisym block ===
    for i in range(n_anti):
        ii = n_sym + i
        for j in range(i, n_anti):
            jj = n_sym + j
            bp, bq = basis_antisym[i], basis_antisym[j]
            S[ii, jj] = overlap_element_pstate_eckart_antisym(bp, bq, alpha, beta)
            T_ij = _kinetic_via_quadrature_pstate(
                bp, bq, alpha, beta, n_r=n_r, n_theta=n_theta
            )
            Vne_ij = potential_vne_element_pstate_eckart_antisym(
                bp, bq, alpha, beta, Z=Z
            )
            Vee_ij = potential_vee_element_pstate_eckart_antisym(
                bp, bq, alpha, beta
            )
            H[ii, jj] = T_ij + Vne_ij + Vee_ij
            S[jj, ii] = S[ii, jj]
            H[jj, ii] = H[ii, jj]

    # === Sym × Antisym cross blocks ===
    # At β=0, all cross-sector elements vanish identically (sinh→0).
    if beta != 0.0:
        for i in range(n_sym):
            for j in range(n_anti):
                jj = n_sym + j
                bp_plus = basis_sym[i]
                bq_minus = basis_antisym[j]
                S[i, jj] = overlap_element_pstate_eckart_cross(
                    bp_plus, bq_minus, alpha, beta
                )
                T_cross = _kinetic_via_quadrature_pstate(
                    bp_plus, bq_minus, alpha, beta, n_r=n_r, n_theta=n_theta
                )
                Vne_cross = potential_vne_element_pstate_eckart_cross(
                    bp_plus, bq_minus, alpha, beta, Z=Z
                )
                Vee_cross = potential_vee_element_pstate_eckart_cross(
                    bp_plus, bq_minus, alpha, beta
                )
                H[i, jj] = T_cross + Vne_cross + Vee_cross
                S[jj, i] = S[i, jj]
                H[jj, i] = H[i, jj]

    H = 0.5 * (H + H.T)
    S = 0.5 * (S + S.T)
    return H, S


def solve_2p1_state_full(
    basis_sym: List[HylleraasPStateBasisFn],
    basis_antisym: List[HylleraasPStateBasisFn],
    alpha: float, Z: float, beta: float,
    n_r: int = 32, n_theta: int = 16,
) -> Tuple[float, np.ndarray]:
    """Solve generalized eigenvalue problem for full two-channel 2^1P."""
    from scipy.linalg import eigh
    H, S = assemble_p_state_matrices_full(
        basis_sym, basis_antisym, alpha, Z, beta, n_r, n_theta
    )
    eigvals, eigvecs = eigh(H, b=S)
    return float(eigvals[0]), eigvecs[:, 0]


def optimize_2p1_full(
    basis_sym: List[HylleraasPStateBasisFn],
    basis_antisym: List[HylleraasPStateBasisFn],
    Z: float,
    alpha_init: float = 1.35,
    beta_init: float = 0.35,
    method: str = 'Nelder-Mead',
    tol: float = 1e-6,
    max_iter: int = 80,
    n_r: int = 24, n_theta: int = 12,
) -> Dict[str, Any]:
    """2D variational optimization over (alpha, beta) for full two-channel 2^1P.

    Returns dict with E_2p, alpha_opt, beta_opt, coeffs, basis info.
    """
    from scipy.optimize import minimize

    def obj(x):
        a, b = float(x[0]), float(x[1])
        if a <= 0 or abs(b) >= a:
            return 1e10
        try:
            E, _c = solve_2p1_state_full(
                basis_sym, basis_antisym, a, Z, b, n_r=n_r, n_theta=n_theta
            )
            return E
        except Exception:
            return 1e10

    x0 = np.array([alpha_init, beta_init])
    result = minimize(
        obj, x0=x0, method=method,
        options={'xatol': tol, 'fatol': tol, 'maxiter': max_iter, 'adaptive': True},
    )
    alpha_opt, beta_opt = float(result.x[0]), float(result.x[1])
    E_2p, coeffs = solve_2p1_state_full(
        basis_sym, basis_antisym, alpha_opt, Z, beta_opt, n_r=n_r, n_theta=n_theta
    )
    return {
        'E_2p': E_2p,
        'alpha_opt': alpha_opt,
        'beta_opt': beta_opt,
        'coeffs': coeffs,
        'basis_sym': list(basis_sym),
        'basis_antisym': list(basis_antisym),
        'n_sym': len(basis_sym),
        'n_antisym': len(basis_antisym),
        'opt_n_iter': int(result.nit),
        'opt_n_eval': int(result.nfev),
        'opt_success': bool(result.success),
    }


def compute_dipole_1s_to_2p_full(
    s_state: Dict[str, Any],
    p_state_full: Dict[str, Any],
) -> Dict[str, float]:
    """Compute <Psi_{1^1S} | z_1+z_2 | Psi_{2^1P}> for the FULL two-channel
    2^1P trial.

    Returns dict with D_total, D_sym (contribution from sym channel) and
    D_antisym (contribution from antisym channel).
    """
    basis_S = s_state['basis']
    c_S = s_state['coeffs']
    alpha_S = s_state['alpha']
    beta_S = s_state.get('beta', 0.0)

    basis_sym = p_state_full['basis_sym']
    basis_antisym = p_state_full['basis_antisym']
    c_full = p_state_full['coeffs']
    alpha_P = p_state_full['alpha_opt']
    beta_P = p_state_full['beta_opt']

    n_sym = len(basis_sym)
    n_antisym = len(basis_antisym)
    c_sym = c_full[:n_sym]
    c_antisym = c_full[n_sym:n_sym + n_antisym]

    # Sym contribution
    D_sym = 0.0
    for i, bf_s in enumerate(basis_S):
        for j, bf_p in enumerate(basis_sym):
            D_sym += c_S[i] * c_sym[j] * dipole_element_1s_2p_sym(
                bf_s, bf_p, alpha_S, beta_S, alpha_P, beta_P
            )

    # Antisym contribution
    D_antisym = 0.0
    for i, bf_s in enumerate(basis_S):
        for j, bf_p in enumerate(basis_antisym):
            D_antisym += c_S[i] * c_antisym[j] * dipole_element_1s_2p_antisym(
                bf_s, bf_p, alpha_S, beta_S, alpha_P, beta_P
            )

    return {
        'D_total': D_sym + D_antisym,
        'D_sym': D_sym,
        'D_antisym': D_antisym,
    }


def oscillator_strength_2p_to_1s_full(
    s_state: Dict[str, Any],
    p_state_full: Dict[str, Any],
) -> Dict[str, Any]:
    """Compute the He 2^1P → 1^1S oscillator strength with FULL two-channel
    2^1P trial.

    For an L=0 → L'=1 transition, Wigner-Eckart relates the single-component
    dipole to the reduced matrix element:
        |<L'||r||L>|^2 = (2L'+1) |<L', M=0 | r_z | L, 0>|^2 = 3 |D_z|^2
    so the standard absorption oscillator strength becomes
        f = (2/3) * Delta_E * |<L'||r||L>|^2 = 2 * Delta_E * |D_z|^2

    where D_z = <Psi_{1^1S} | z_1+z_2 | Psi_{2^1P}> summing over both sym and
    antisym channels. Validated against hydrogen 1S->2P (f = 0.4162).
    """
    dip = compute_dipole_1s_to_2p_full(s_state, p_state_full)
    D = dip['D_total']
    E_S = s_state['energy']
    E_P = p_state_full['E_2p']
    delta_E = E_P - E_S
    f = 2.0 * delta_E * D * D
    return {
        'f': f,
        'D': D,
        'D_sym': dip['D_sym'],
        'D_antisym': dip['D_antisym'],
        'D_squared': D * D,
        'delta_E': delta_E,
        'E_1s': E_S,
        'E_2p': E_P,
        'alpha_1s': s_state['alpha'],
        'beta_1s': s_state.get('beta', 0.0),
        'alpha_2p': p_state_full['alpha_opt'],
        'beta_2p': p_state_full['beta_opt'],
    }


def oscillator_strength_2p_to_1s_sym(
    s_state: Dict[str, Any],
    p_state: Dict[str, Any],
) -> Dict[str, Any]:
    """Compute the He 2^1P -> 1^1S oscillator strength (sym-only 2^1P trial).

    f = 2 * (E_P - E_S) * |D_z|^2 with D_z = <1^1S | z_1+z_2 | 2^1P, sym>

    The factor of 2 (rather than the bare 2/3 in some texts) absorbs the
    Wigner-Eckart sum over final M_L states for L=0 -> L'=1 transitions.
    Validated against hydrogen 1S -> 2P f = 0.4162.

    For absorption from ground state, E_P > E_S and f > 0.

    Returns dict with f, dipole D, transition energy Delta_E, and intermediate
    quantities.
    """
    D = compute_dipole_1s_to_2p_sym(s_state, p_state)
    E_S = s_state['energy']
    E_P = p_state['E_2p']
    delta_E = E_P - E_S
    f = 2.0 * delta_E * D * D
    return {
        'f': f,
        'D': D,
        'D_squared': D * D,
        'delta_E': delta_E,
        'E_1s': E_S,
        'E_2p': E_P,
        'alpha_1s': s_state['alpha'],
        'beta_1s': s_state.get('beta', 0.0),
        'alpha_2p': p_state['alpha_opt'],
        'beta_2p': p_state['beta_opt'],
    }
