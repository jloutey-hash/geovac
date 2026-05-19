"""
Hylleraas-Eckart master integral via algebraic recurrence
=========================================================

Implements the closed-form Hylleraas-Eckart master integral

    I_HE^cosh(L, M, N; alpha, B)
        = integral_0^infty ds exp(-2 alpha s) s^L
          * integral_0^s du u^{N+1}
          * integral_{-u}^u dt t^{2M} (s^2 - t^2) cosh(B t)

WITHOUT sympy in the hot loop. The three-stage integration

  Stage 1 (t-integration):  closed via integration-by-parts recurrence in M
  Stage 2 (u-integration):  closed via IBP recurrence on int_0^s u^a {sinh,cosh}(Bu) du
  Stage 3 (s-integration):  closed via int_0^inf s^k e^{-2as} {1,sinh,cosh}(Bs) ds

each have integer-rational-coefficient closed forms. Composing them gives
the same ClosedForm output as `hylleraas_eckart_closed_forms.I_HE_cosh_polynomial`
but in MICROSECONDS per cell instead of SECONDS-TO-MINUTES (sympy bottleneck
is `cancel`/`simplify`, not the math itself).

This module is the GeoVac algebraic-first realization of Track 1.5 of the
Hylleraas-Eckart sprint, per CLAUDE.md feedback_algebraic_first.md.

Stage 1: G_{2M}(u, B) recurrence
--------------------------------

Define G_k(u, B) = int_{-u}^u t^k cosh(B t) dt. For even k = 2M, two IBPs
in t give

    G_0(u, B)     = 2 sinh(Bu) / B
    G_{2M}(u, B)  = (2 u^{2M} sinh(Bu)) / B
                  + (4 M u^{2M-1} cosh(Bu)) / B^2
                  - (2 M (2M-1) / B^2) G_{2M-2}(u, B)            (M >= 1)

The t-integrated factor in the master integral is

    s^2 G_{2M}(u, B) - G_{2M+2}(u, B)

(from int_{-u}^u t^{2M} (s^2 - t^2) cosh(Bt) dt = s^2 G_{2M} - G_{2M+2}).

Multiplied by u^{N+1} (the Jacobian u-factor plus basis u^N) and ready for
stage 2.

Stage 2: u-integration via IBP
------------------------------

Let

    F^sinh_a(s, B) = int_0^s u^a sinh(Bu) du
    F^cosh_a(s, B) = int_0^s u^a cosh(Bu) du

Base:
    F^sinh_0 = (cosh(Bs) - 1) / B
    F^cosh_0 = sinh(Bs) / B

Recurrence (IBP):
    F^sinh_a = s^a cosh(Bs) / B - (a/B) F^cosh_{a-1}    (a >= 1)
    F^cosh_a = s^a sinh(Bs) / B - (a/B) F^sinh_{a-1}    (a >= 1)

These telescope to a polynomial in s, with sinh(Bs), cosh(Bs), and the
constant -1/B (from F^sinh_0) as kernels.

Stage 3: s-integration (textbook)
---------------------------------

    int_0^inf s^k e^{-2 alpha s} ds         = k! / (2 alpha)^{k+1}
    int_0^inf s^k e^{-2 alpha s} sinh(Bs) ds = (k!/2)[1/(2a-B)^{k+1} - 1/(2a+B)^{k+1}]
    int_0^inf s^k e^{-2 alpha s} cosh(Bs) ds = (k!/2)[1/(2a-B)^{k+1} + 1/(2a+B)^{k+1}]

Putting over common denominator (4 alpha^2 - B^2)^{k+1}:

    [1/(2a-B)^{k+1} +- 1/(2a+B)^{k+1}]
        = [(2a+B)^{k+1} +- (2a-B)^{k+1}] / (4 alpha^2 - B^2)^{k+1}

By symmetry the "minus" (sinh) version has only odd-B numerator powers
and the "plus" (cosh) version has only even-B numerator powers. Combined
with the 1/B^k factors from stages 1+2, every term contributes an
even-B-power numerator term over (4 alpha^2 - B^2)^k denominator, matching
the t-parity of the original integral.

Architecture
------------

* `UExpr` : list of terms `(coeff, b_exp, s_pow, u_pow, kind)` where
  kind in {'sinh_Bu', 'cosh_Bu'}, representing
  sum_i coeff_i * B^{b_exp_i} * s^{s_pow_i} * u^{u_pow_i} * kind_i(Bu).

* `SExpr` : list of terms `(coeff, b_exp, s_pow, kind)` where
  kind in {'one', 'sinh_Bs', 'cosh_Bs'}, representing
  sum_i coeff_i * B^{b_exp_i} * s^{s_pow_i} * kind_i(Bs).

* `_build_G(M_max)` : list of UExpr for G_0, G_2, ..., G_{2 M_max}.
* `_t_integrated(L, M, N)` : UExpr for u^{N+1} (s^2 G_{2M} - G_{2M+2}).
* `_integrate_u(uexpr)` : SExpr from integrating u from 0 to s.
* `_integrate_s(sexpr, L)` : ClosedForm from integrating s from 0 to inf.
* `I_HE_cosh_polynomial(L, M, N)` : main API, composes the three stages.

Each step is integer-rational polynomial arithmetic. The whole thing is
deterministic and < 1 ms per cell at moderate (L, M, N).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Tuple

# Import ClosedForm + cache machinery from the existing module.
from geovac.hylleraas_eckart_closed_forms import (
    ClosedForm,
    _CACHE_TABLE,
)


# ---------------------------------------------------------------------------
# Intermediate symbolic expressions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UTerm:
    """Term: coeff * B^b_exp * s^s_pow * u^u_pow * kind(B*u).

    kind in {'sinh_Bu', 'cosh_Bu'}.
    """
    coeff: Fraction
    b_exp: int
    s_pow: int
    u_pow: int
    kind: str  # 'sinh_Bu' or 'cosh_Bu'


@dataclass(frozen=True)
class STerm:
    """Term: coeff * B^b_exp * s^s_pow * kind(B*s).

    kind in {'one', 'sinh_Bs', 'cosh_Bs'}.
    'one' has no hyperbolic factor (constant boundary terms from u=0).
    """
    coeff: Fraction
    b_exp: int
    s_pow: int
    kind: str  # 'one', 'sinh_Bs', 'cosh_Bs'


class UExpr:
    """Bag of UTerms representing a polynomial in (u, s, B) with sinh/cosh of Bu."""

    __slots__ = ("terms",)

    def __init__(self, terms: List[UTerm] = None):
        self.terms = list(terms) if terms else []

    def __add__(self, other: "UExpr") -> "UExpr":
        return UExpr(self.terms + other.terms)

    def scaled(self, factor: Fraction, b_shift: int) -> "UExpr":
        return UExpr([
            UTerm(t.coeff * factor, t.b_exp + b_shift, t.s_pow, t.u_pow, t.kind)
            for t in self.terms
        ])

    def shift_s(self, ds: int) -> "UExpr":
        return UExpr([
            UTerm(t.coeff, t.b_exp, t.s_pow + ds, t.u_pow, t.kind)
            for t in self.terms
        ])

    def shift_u(self, du: int) -> "UExpr":
        return UExpr([
            UTerm(t.coeff, t.b_exp, t.s_pow, t.u_pow + du, t.kind)
            for t in self.terms
        ])


class SExpr:
    """Bag of STerms representing a polynomial in (s, B) with sinh/cosh of Bs."""

    __slots__ = ("terms",)

    def __init__(self, terms: List[STerm] = None):
        self.terms = list(terms) if terms else []

    def __add__(self, other: "SExpr") -> "SExpr":
        return SExpr(self.terms + other.terms)

    def scaled(self, factor: Fraction, b_shift: int) -> "SExpr":
        return SExpr([
            STerm(t.coeff * factor, t.b_exp + b_shift, t.s_pow, t.kind)
            for t in self.terms
        ])

    def shift_s(self, ds: int) -> "SExpr":
        return SExpr([
            STerm(t.coeff, t.b_exp, t.s_pow + ds, t.kind)
            for t in self.terms
        ])


# ---------------------------------------------------------------------------
# Stage 1: G_{2M}(u, B) recurrence
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _G_2M(M: int) -> UExpr:
    """G_{2M}(u, B) = int_{-u}^u t^{2M} cosh(Bt) dt.

    Stored as UExpr in (u, B), with s_pow=0 throughout.
    """
    if M == 0:
        # G_0 = 2 sinh(Bu) / B  =>  coeff=2, b_exp=-1, s_pow=0, u_pow=0, kind=sinh
        return UExpr([UTerm(Fraction(2), -1, 0, 0, 'sinh_Bu')])

    # G_{2M} = 2 u^{2M} sinh(Bu) / B
    #        - 4 M u^{2M-1} cosh(Bu) / B^2
    #        + 2 M (2M-1) / B^2 * G_{2(M-1)}
    # (Derived via two IBPs in t over the symmetric interval [-u, u];
    # the cosh boundary term picks up a minus sign because the inner
    # integrand t^{2M-1} sinh(Bt) is even in t, so both endpoints
    # contribute additively after the t -> -t substitution.)
    term_a = UExpr([UTerm(Fraction(2), -1, 0, 2 * M, 'sinh_Bu')])
    term_b = UExpr([UTerm(Fraction(-4 * M), -2, 0, 2 * M - 1, 'cosh_Bu')])
    prev_G = _G_2M(M - 1)
    term_c = prev_G.scaled(Fraction(2 * M * (2 * M - 1)), -2)
    return term_a + term_b + term_c


@lru_cache(maxsize=None)
def _H_2M_plus_1(M: int) -> UExpr:
    """H_{2M+1}(u, B) = int_{-u}^u t^{2M+1} sinh(Bt) dt.

    Closes via single-step IBP:
        H_{2M+1}(u, B) = 2 u^{2M+1} cosh(Bu) / B - (2M+1)/B * G_{2M}(u, B)

    Note this is a ONE-step recurrence (not a chain): each H value directly
    expressible via the existing G_{2M} from `_G_2M`. So no separate cache
    needed beyond the lru_cache on this function itself.

    Returns a UExpr with sinh_Bu and cosh_Bu terms (same data type as G).
    """
    if M < 0:
        raise ValueError(f"M must be >= 0, got {M}")
    # 2 u^{2M+1} cosh(Bu) / B
    head = UExpr([UTerm(Fraction(2), -1, 0, 2 * M + 1, 'cosh_Bu')])
    # -(2M+1)/B * G_{2M}(u, B)
    G_2M = _G_2M(M)
    tail = G_2M.scaled(Fraction(-(2 * M + 1)), -1)
    return head + tail


def _t_integrated(L: int, M: int, N: int) -> UExpr:
    """u^{N+1} * (s^2 G_{2M}(u,B) - G_{2M+2}(u,B)).

    Used for OVERLAP master integral I_HE^cosh: includes (s^2 - t^2)
    factor from the volume element and u^{N+1} factor (basis u^N times
    Jacobian u).
    """
    return _t_integrated_general(M=M, u_pow=N + 1, include_s2_minus_t2=True)


def _t_integrated_general(M: int, u_pow: int, include_s2_minus_t2: bool) -> UExpr:
    """Generalized t-integrated factor for Eckart-class master integrals.

    Parameters
    ----------
    M : int
        Half-power of t in the basis (basis has t^{2*M}, where M = m_p + m_q).
    u_pow : int
        Final u-power on the integrand. For overlap (I) and V_ne (J), this is
        N+1 (basis u^N times Jacobian u). For V_ee (K), this is N (the 1/r_12
        factor cancels the Jacobian u).
    include_s2_minus_t2 : bool
        True for overlap (I) and V_ee (K) where the (s^2-t^2) volume factor
        is present. False for V_ne (J) where 4s/(s^2-t^2) cancels the volume
        (s^2-t^2), leaving the integrand without it.
    """
    G_2M = _G_2M(M)
    if include_s2_minus_t2:
        G_2Mp2 = _G_2M(M + 1)
        # s^2 G_{2M} - G_{2M+2}
        integrand = G_2M.shift_s(2) + G_2Mp2.scaled(Fraction(-1), 0)
    else:
        integrand = G_2M
    if u_pow != 0:
        integrand = integrand.shift_u(u_pow)
    return integrand


# ---------------------------------------------------------------------------
# Stage 2: u-integration via IBP recurrence
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _F_sinh(a: int) -> SExpr:
    """F^sinh_a(s, B) = int_0^s u^a sinh(Bu) du as SExpr in (s, B).

    Base a=0:  (cosh(Bs) - 1) / B
              = (1/B) cosh(Bs) - (1/B) * one
    Recurrence a>=1: s^a cosh(Bs)/B - (a/B) F^cosh_{a-1}
    """
    if a == 0:
        return SExpr([
            STerm(Fraction(1), -1, 0, 'cosh_Bs'),
            STerm(Fraction(-1), -1, 0, 'one'),
        ])
    # s^a cosh(Bs)/B
    head = SExpr([STerm(Fraction(1), -1, a, 'cosh_Bs')])
    # -(a/B) F^cosh_{a-1}
    tail = _F_cosh(a - 1).scaled(Fraction(-a), -1)
    return head + tail


@lru_cache(maxsize=None)
def _F_cosh(a: int) -> SExpr:
    """F^cosh_a(s, B) = int_0^s u^a cosh(Bu) du as SExpr in (s, B).

    Base a=0:  sinh(Bs) / B
    Recurrence a>=1: s^a sinh(Bs)/B - (a/B) F^sinh_{a-1}
    """
    if a == 0:
        return SExpr([STerm(Fraction(1), -1, 0, 'sinh_Bs')])
    head = SExpr([STerm(Fraction(1), -1, a, 'sinh_Bs')])
    tail = _F_sinh(a - 1).scaled(Fraction(-a), -1)
    return head + tail


def _integrate_u(uexpr: UExpr) -> SExpr:
    """Integrate UExpr's u-dependence from 0 to s.

    For each term coeff * B^b * s^p * u^q * K(Bu), substitute
    int_0^s u^q K(Bu) du = F^K_q(s, B) and combine.
    """
    out_terms: List[STerm] = []
    for t in uexpr.terms:
        if t.kind == 'sinh_Bu':
            F = _F_sinh(t.u_pow)
        elif t.kind == 'cosh_Bu':
            F = _F_cosh(t.u_pow)
        else:
            raise ValueError(f"unknown UTerm kind: {t.kind}")
        # F is SExpr in (s, B). Multiply each F-term by coeff * B^b * s^p.
        for st in F.terms:
            out_terms.append(STerm(
                coeff=t.coeff * st.coeff,
                b_exp=t.b_exp + st.b_exp,
                s_pow=t.s_pow + st.s_pow,
                kind=st.kind,
            ))
    return SExpr(out_terms)


# ---------------------------------------------------------------------------
# Stage 3: s-integration -> ClosedForm
# ---------------------------------------------------------------------------

def _binomial_expand_sum(j_plus_one: int, sign: int
                         ) -> Dict[Tuple[int, int], Fraction]:
    """Compute (2a + B)^{j+1} + sign * (2a - B)^{j+1} as a polynomial in (a, B).

    Returns dict mapping (alpha_power, B_power) -> Fraction coefficient.

    sign = +1 -> only even B-powers survive (cosh-type numerator).
    sign = -1 -> only odd B-powers survive (sinh-type numerator).
    """
    out: Dict[Tuple[int, int], Fraction] = {}
    n = j_plus_one
    for k in range(n + 1):
        c_plus = math.comb(n, k)  # for (2a + B)^n: coefficient of (2a)^{n-k} B^k
        c_minus = math.comb(n, k) * ((-1) ** k)  # for (2a - B)^n
        c_total = c_plus + sign * c_minus
        if c_total == 0:
            continue
        # (2a)^{n-k} B^k = 2^{n-k} a^{n-k} B^k
        a_pow = n - k
        b_pow = k
        coeff = Fraction(c_total * (2 ** (n - k)))
        out[(a_pow, b_pow)] = out.get((a_pow, b_pow), Fraction(0)) + coeff
    return out


def _scale_poly2d(poly: Dict[Tuple[int, int], Fraction], factor: Fraction
                  ) -> Dict[Tuple[int, int], Fraction]:
    return {k: v * factor for k, v in poly.items()}


def _add_poly2d_inplace(target: Dict[Tuple[int, int], Fraction],
                         poly: Dict[Tuple[int, int], Fraction]) -> None:
    for k, v in poly.items():
        target[k] = target.get(k, Fraction(0)) + v


def _shift_b_poly2d(poly: Dict[Tuple[int, int], Fraction], db: int
                    ) -> Dict[Tuple[int, int], Fraction]:
    """Multiply polynomial by B^db (shift all B-exponents by db)."""
    out: Dict[Tuple[int, int], Fraction] = {}
    for (a, b), c in poly.items():
        new_b = b + db
        # Negative B-power should never survive in the FINAL result, but
        # may appear in intermediates as we accumulate over a common
        # denominator. We allow them here and verify at the end.
        out[(a, new_b)] = c
    return out


def _integrate_s_to_closed_form(sexpr: SExpr, L: int,
                                 b_parity: str = 'even') -> ClosedForm:
    """Integrate SExpr over s in [0, inf) with extra s^L factor.

    b_parity ('even' default, or 'odd'): the parity of B-powers expected
    in the result. For cosh-type masters (integrand even in t and B),
    b_parity='even' and the result has only even B-powers. For sinh-type
    masters (integrand even in t but odd in B), b_parity='odd' and the
    result has only odd B-powers, which we factor out as B^1 *
    (even-B-power polynomial) and return with b_factor=1.

    Each SExpr term contributes:
      - kind='one':       coeff * B^{b} * s^{p+L} -> coeff * B^b * (p+L)! / (2a)^{p+L+1}
      - kind='sinh_Bs':   ... (p+L)!/2 [1/(2a-B)^{p+L+1} - 1/(2a+B)^{p+L+1}]
      - kind='cosh_Bs':   ... (p+L)!/2 [1/(2a-B)^{p+L+1} + 1/(2a+B)^{p+L+1}]

    We use a common denominator (4 alpha^2 - B^2)^k_max for the sinh/cosh
    pieces, and for the 'one' pieces we multiply by (2 alpha)^{...} as
    appropriate. Final ClosedForm has numerator/denominator in (alpha, B^2).
    """
    # First, determine the max denominator power so we can put everything
    # over a common denominator (4 a^2 - B^2)^{k_max}.
    # For 'one' kind with s^{p+L}, denominator is (2a)^{p+L+1}; we'll need
    # to convert that to (4 a^2 - B^2)^{k_max} as well. Approach: separately
    # accumulate "rational over (4 a^2 - B^2)^k_max" and "rational over
    # (2a)^k" pieces, then combine at the end.

    # Compute k_max from sinh/cosh contributions.
    k_max_hyp = 0
    for t in sexpr.terms:
        if t.kind in ('sinh_Bs', 'cosh_Bs'):
            k = t.s_pow + L + 1  # exponent of (2a +- B) in the closed form
            k_max_hyp = max(k_max_hyp, k)

    # For 'one' terms, we have (2a)^{p+L+1} denominators. We'll put them
    # over (2a)^{k_max_one}, then convert to (4a^2 - B^2)^{k_max_one} by
    # multiplying num and den by (2a + B)^{k_max_one} -- wait, that's not
    # quite right. (2a)^k != (4a^2 - B^2)^{k/2}. We need a separate
    # accumulator for the 'one' contributions, and combine at the end.

    # Easier: compute each contribution as a (numerator dict, denominator power)
    # pair, then put everything over a common denominator at the end.

    # numerator polynomial in (alpha, B^2-index)... actually let me store as
    # (alpha_power, B_power) in numerator + denominator factors (which we
    # track as: alpha_power, B_power, denom_kind in {(2a-B), (2a+B), (2a)}).

    # Simplest: accumulate ALL terms over the common denominator
    # D = (4 a^2 - B^2)^{k_max_hyp} * (2 a)^{k_max_one}. After collecting
    # everything, polynomially divide to a canonical reduced form, or just
    # output as the unreduced rational.

    # Determine k_max_one.
    k_max_one = 0
    for t in sexpr.terms:
        if t.kind == 'one':
            k = t.s_pow + L + 1
            k_max_one = max(k_max_one, k)

    # Numerator polynomial in (alpha, B), with arbitrary B-powers
    # (we'll combine at the end). Stored as {(a_pow, b_pow): Fraction}.
    num: Dict[Tuple[int, int], Fraction] = {}

    # Common denominator: (4 a^2 - B^2)^{k_max_hyp} * (2 a)^{k_max_one}.

    for t in sexpr.terms:
        j_plus_one = t.s_pow + L + 1
        # The s^{p+L} factor in the integrand gives factorial (p+L)! = (j+1-1)! = j!
        # where j = s_pow + L. So (p+L)! = (j_plus_one - 1)!
        factorial = Fraction(math.factorial(j_plus_one - 1))

        if t.kind == 'one':
            # Contribution: coeff * B^{b_exp} * factorial / (2a)^{j_plus_one}.
            # Over common denom (4 a^2 - B^2)^{k_max_hyp} * (2a)^{k_max_one},
            # we multiply num by (4 a^2 - B^2)^{k_max_hyp} * (2a)^{k_max_one - j_plus_one}.
            multiplier_4a2mB2 = _pow_4a2_minus_B2(k_max_hyp)
            extra_2a_pow = k_max_one - j_plus_one
            factor = Fraction(2 ** extra_2a_pow)
            # (2a)^extra means a^extra * 2^extra. Multiply by a^{extra} on num.
            local_num: Dict[Tuple[int, int], Fraction] = {}
            for (a_pow, b_pow), c in multiplier_4a2mB2.items():
                key = (a_pow + extra_2a_pow, b_pow + t.b_exp)
                local_num[key] = local_num.get(key, Fraction(0)) + c * factor * t.coeff * factorial
            _add_poly2d_inplace(num, local_num)

        elif t.kind in ('sinh_Bs', 'cosh_Bs'):
            # Closed form:
            #   sinh: (p+L)!/2 [1/(2a-B)^{j+1} - 1/(2a+B)^{j+1}]
            #   cosh: (p+L)!/2 [1/(2a-B)^{j+1} + 1/(2a+B)^{j+1}]
            # Over common denom (4 a^2 - B^2)^{k_max_hyp} * (2a)^{k_max_one},
            # this contributes
            #   factorial/2 * [(2a+B)^{j+1} +- (2a-B)^{j+1}] * (4a^2-B^2)^{k_max_hyp - j_plus_one} * (2a)^{k_max_one}
            sign = -1 if t.kind == 'sinh_Bs' else +1
            binom_poly = _binomial_expand_sum(j_plus_one, sign)
            half_fact_coeff = factorial * Fraction(1, 2)

            # Pad the (4 a^2 - B^2)^{k_max_hyp - j_plus_one} factor.
            pad_4a2mB2 = _pow_4a2_minus_B2(k_max_hyp - j_plus_one)
            # And the (2a)^{k_max_one} factor.
            pad_2a_factor = Fraction(2 ** k_max_one)  # 2^{k_max_one} for the constant
            pad_2a_a_pow = k_max_one  # extra alpha power

            # Multiply binom_poly * pad_4a2mB2 * pad_2a -> contribution to num.
            # Then multiply by coeff * B^{b_exp} * half_fact_coeff * sign-adjust.
            # (the sign was already baked into binom_poly via the sign argument)
            intermediate: Dict[Tuple[int, int], Fraction] = {}
            for (a1, b1), c1 in binom_poly.items():
                for (a2, b2), c2 in pad_4a2mB2.items():
                    key = (a1 + a2 + pad_2a_a_pow, b1 + b2 + t.b_exp)
                    intermediate[key] = intermediate.get(key, Fraction(0)) + c1 * c2 * pad_2a_factor

            multiplier = t.coeff * half_fact_coeff
            for k, v in intermediate.items():
                num[k] = num.get(k, Fraction(0)) + v * multiplier

        else:
            raise ValueError(f"unknown STerm kind: {t.kind}")

    # Now num is a dict {(alpha_pow, B_pow): Fraction}. Verify all NON-ZERO
    # entries have non-negative B-powers of the expected parity.
    expected_parity = 0 if b_parity == 'even' else 1
    for (a_pow, b_pow), c in num.items():
        if c == 0:
            continue  # zero entries are bookkeeping artifacts, ignore
        if b_pow < 0:
            raise ValueError(
                f"Negative B-power {b_pow} in numerator: c={c}, a_pow={a_pow}. "
                "Indicates cancellation failure -- check the recurrence."
            )
        if b_pow % 2 != expected_parity:
            raise ValueError(
                f"B-power {b_pow} has wrong parity (expected {b_parity}) "
                f"in numerator: c={c}, a_pow={a_pow}."
            )

    # If b_parity is 'odd', factor out one explicit B from the numerator
    # so the remaining polynomial has even B-powers (storable in our nested
    # (alpha, B^2) representation).
    if b_parity == 'odd':
        num = {(a_pow, b_pow - 1): c for (a_pow, b_pow), c in num.items()}
        b_factor = 1
    else:
        b_factor = 0

    # Build denominator polynomial: (4 a^2 - B^2)^{k_max_hyp} * (2 a)^{k_max_one}.
    den_poly = _pow_4a2_minus_B2(k_max_hyp)
    # Multiply by (2a)^{k_max_one} = 2^{k_max_one} * a^{k_max_one}.
    den_poly = _multiply_by_2a_pow(den_poly, k_max_one)

    # Numerator and denominator are polynomials in (alpha, B^2). Convert
    # to the nested-tuple Fraction representation expected by ClosedForm.
    num_coeffs = _poly2d_to_nested_alpha_b2(num)
    den_coeffs = _poly2d_to_nested_alpha_b2(den_poly)

    return ClosedForm(
        L=L, M=0, N=0,  # M, N set by caller
        num_coeffs=num_coeffs,
        den_coeffs=den_coeffs,
        b_factor=b_factor,
    )


# Helpers for stage 3 polynomial accumulation
# -------------------------------------------

@lru_cache(maxsize=None)
def _pow_4a2_minus_B2(k: int) -> Dict[Tuple[int, int], Fraction]:
    """(4 a^2 - B^2)^k as polynomial dict {(a_pow, b_pow): coeff}."""
    if k == 0:
        return {(0, 0): Fraction(1)}
    if k == 1:
        return {(2, 0): Fraction(4), (0, 2): Fraction(-1)}
    # Iterative multiplication.
    out = {(0, 0): Fraction(1)}
    base = {(2, 0): Fraction(4), (0, 2): Fraction(-1)}
    for _ in range(k):
        new_out: Dict[Tuple[int, int], Fraction] = {}
        for (a1, b1), c1 in out.items():
            for (a2, b2), c2 in base.items():
                key = (a1 + a2, b1 + b2)
                new_out[key] = new_out.get(key, Fraction(0)) + c1 * c2
        out = new_out
    return out


def _multiply_by_2a_pow(poly: Dict[Tuple[int, int], Fraction], k: int
                        ) -> Dict[Tuple[int, int], Fraction]:
    """Multiply polynomial by (2a)^k = 2^k * a^k."""
    if k == 0:
        return dict(poly)
    factor = Fraction(2 ** k)
    return {(a_pow + k, b_pow): c * factor for (a_pow, b_pow), c in poly.items()}


def _poly2d_to_nested_alpha_b2(poly: Dict[Tuple[int, int], Fraction]
                                ) -> Tuple[Tuple[Fraction, ...], ...]:
    """Convert {(alpha_pow, B_pow_even): Fraction} to nested
    coeffs[i][j] = coefficient of alpha^i * (B^2)^j.

    Assumes B_pow is even on every key (enforced by upstream parity check).
    """
    if not poly:
        return ((Fraction(0),),)

    # Bucket by alpha-power and j := b_pow // 2.
    max_alpha = 0
    max_b2 = 0
    bucket: Dict[Tuple[int, int], Fraction] = {}
    for (a_pow, b_pow), c in poly.items():
        if c == 0:
            continue
        j = b_pow // 2
        bucket[(a_pow, j)] = c
        max_alpha = max(max_alpha, a_pow)
        max_b2 = max(max_b2, j)

    if not bucket:
        return ((Fraction(0),),)

    rows = []
    for i in range(max_alpha + 1):
        row = []
        for j in range(max_b2 + 1):
            row.append(bucket.get((i, j), Fraction(0)))
        # Trim trailing zeros.
        while len(row) > 1 and row[-1] == 0:
            row.pop()
        rows.append(tuple(row))
    while len(rows) > 1 and all(c == 0 for c in rows[-1]):
        rows.pop()
    return tuple(rows)


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def I_HE_cosh_polynomial_recurrence(L: int, M: int, N: int) -> ClosedForm:
    """Compute I_HE^cosh(L, M, N) ClosedForm via algebraic recurrence (no sympy).

    OVERLAP master integral:
        I_HE^cosh(L, M, N; a, B)
            = integral exp(-2 a s) s^L u^{N+1} t^{2M} (s^2 - t^2) cosh(B t) dt du ds

    Equivalent to the sympy-based I_HE_cosh_polynomial but ~1000x-100000x
    faster. Result evaluates identically to the sympy oracle (different
    polynomial reduction is possible).
    """
    if L < 0 or M < 0 or N < 0:
        raise ValueError(f"L, M, N must be >= 0, got ({L}, {M}, {N})")

    uexpr = _t_integrated_general(M=M, u_pow=N + 1, include_s2_minus_t2=True)
    sexpr = _integrate_u(uexpr)
    cf = _integrate_s_to_closed_form(sexpr, L)
    return ClosedForm(L=L, M=M, N=N,
                      num_coeffs=cf.num_coeffs,
                      den_coeffs=cf.den_coeffs)


@lru_cache(maxsize=None)
def J_HE_cosh_polynomial_recurrence(L: int, M: int, N: int) -> ClosedForm:
    """Compute J_HE^cosh(L, M, N) ClosedForm via algebraic recurrence.

    NUCLEAR ATTRACTION master integral (V_ne):
        J_HE^cosh(L, M, N; a, B)
            = integral exp(-2 a s) s^L u^{N+1} t^{2M} cosh(B t) dt du ds

    (No (s^2 - t^2) factor: 1/r_1 + 1/r_2 = 4s/(s^2-t^2) cancels the volume.)
    Use J_HE^cosh(L+1, M, N) in the V_ne matrix element (because the cancellation
    leaves an extra factor of s).

    At B=0 reduces to the single-alpha V_ne master:
        J(L+1, M, N; a, 0)
            = 2 (L + N + 2M + 4)! / [(2M+1)(N+2M+3)(2a)^{L+N+2M+5}]
    """
    if L < 0 or M < 0 or N < 0:
        raise ValueError(f"L, M, N must be >= 0, got ({L}, {M}, {N})")

    uexpr = _t_integrated_general(M=M, u_pow=N + 1, include_s2_minus_t2=False)
    sexpr = _integrate_u(uexpr)
    cf = _integrate_s_to_closed_form(sexpr, L)
    return ClosedForm(L=L, M=M, N=N,
                      num_coeffs=cf.num_coeffs,
                      den_coeffs=cf.den_coeffs)


@lru_cache(maxsize=None)
def master_C_gen(L: int, u_pow: int, M: int) -> ClosedForm:
    """Generalized cosh master integral (for kinetic-element assembly).

        master_C_gen(L, u_pow, M; alpha, B)
            = int_0^infty ds exp(-2 alpha s) s^L
              * int_0^s du u^{u_pow}
              * int_{-u}^u dt t^{2M} cosh(B t)

    No (s^2 - t^2) factor, no implicit u-Jacobian shift; you control L,
    u_pow, M directly. Use this as a building block for kinetic terms.

    Returns a ClosedForm that can be evaluated cheaply at any (alpha, B).
    """
    if L < 0 or u_pow < 0 or M < 0:
        raise ValueError(f"L, u_pow, M must be >= 0; got ({L}, {u_pow}, {M})")
    uexpr = _G_2M(M).shift_u(u_pow) if u_pow != 0 else _G_2M(M)
    sexpr = _integrate_u(uexpr)
    cf = _integrate_s_to_closed_form(sexpr, L)
    return ClosedForm(L=L, M=M, N=u_pow,
                      num_coeffs=cf.num_coeffs, den_coeffs=cf.den_coeffs)


@lru_cache(maxsize=None)
def master_S_gen(L: int, u_pow: int, M: int) -> ClosedForm:
    """Generalized sinh master integral (with odd t-power; for kinetic).

        master_S_gen(L, u_pow, M; alpha, B)
            = int_0^infty ds exp(-2 alpha s) s^L
              * int_0^s du u^{u_pow}
              * int_{-u}^u dt t^{2M+1} sinh(B t)

    The integrand is EVEN in t (odd t-power times odd sinh = even), so the
    [-u, u] integral is nonzero. Required for kinetic energy matrix elements
    in the Eckart basis (φ_t · φ_t and φ_t · φ_u cross-derivatives produce
    sinh(2β t) factors that combine with explicit t powers to give odd-t
    sinh integrands).
    """
    if L < 0 or u_pow < 0 or M < 0:
        raise ValueError(f"L, u_pow, M must be >= 0; got ({L}, {u_pow}, {M})")
    uexpr = _H_2M_plus_1(M).shift_u(u_pow) if u_pow != 0 else _H_2M_plus_1(M)
    sexpr = _integrate_u(uexpr)
    cf = _integrate_s_to_closed_form(sexpr, L, b_parity='odd')
    return ClosedForm(L=L, M=M, N=u_pow,
                      num_coeffs=cf.num_coeffs, den_coeffs=cf.den_coeffs,
                      b_factor=cf.b_factor)


@lru_cache(maxsize=None)
def K_HE_cosh_polynomial_recurrence(L: int, M: int, N: int) -> ClosedForm:
    """Compute K_HE^cosh(L, M, N) ClosedForm via algebraic recurrence.

    ELECTRON-ELECTRON REPULSION master integral (V_ee):
        K_HE^cosh(L, M, N; a, B)
            = integral exp(-2 a s) s^L u^N t^{2M} (s^2 - t^2) cosh(B t) dt du ds

    (u^N instead of u^{N+1}: 1/r_12 = 1/u cancels the Jacobian u, leaving
    one less power of u in the integrand.)

    At B=0 reduces to the single-alpha V_ee master:
        K(L, M, N; a, 0)
            = 2 [1/((2M+1)(N+2M+2)) - 1/((2M+3)(N+2M+4))]
              * (L + N + 2M + 4)! / (2a)^{L+N+2M+5}

    Note N can be 0 (most common case), giving u^0 = 1 in the integrand;
    the u-integration starts from u^0 K(Bu) which is the F^K_0 base case
    of the u-IBP recurrence.
    """
    if L < 0 or M < 0 or N < 0:
        raise ValueError(f"L, M, N must be >= 0, got ({L}, {M}, {N})")

    uexpr = _t_integrated_general(M=M, u_pow=N, include_s2_minus_t2=True)
    sexpr = _integrate_u(uexpr)
    cf = _integrate_s_to_closed_form(sexpr, L)
    return ClosedForm(L=L, M=M, N=N,
                      num_coeffs=cf.num_coeffs,
                      den_coeffs=cf.den_coeffs)


def precompute_table_recurrence(max_total_degree: int = 8,
                                 verbose: bool = False
                                 ) -> Dict[Tuple[int, int, int], ClosedForm]:
    """Drop-in replacement for hylleraas_eckart_closed_forms.precompute_table
    using the algebraic recurrence engine.
    """
    import time
    table: Dict[Tuple[int, int, int], ClosedForm] = {}
    n_total = 0
    for L in range(max_total_degree + 1):
        for M in range((max_total_degree - L) // 2 + 1):
            for N in range(max_total_degree - L - 2 * M + 1):
                n_total += 1
    if verbose:
        print(f"[precompute_table_recurrence] computing {n_total} cells "
              f"(max_total_degree={max_total_degree}) ...")

    t_start = time.time()
    count = 0
    for L in range(max_total_degree + 1):
        for M in range((max_total_degree - L) // 2 + 1):
            for N in range(max_total_degree - L - 2 * M + 1):
                key = (L, M, N)
                t0 = time.time()
                cf = I_HE_cosh_polynomial_recurrence(L, M, N)
                dt = time.time() - t0
                table[key] = cf
                count += 1
                if verbose and (count % 100 == 0 or n_total < 50):
                    print(f"  [{count}/{n_total}] (L={L}, M={M}, N={N}): "
                          f"{dt * 1000:.2f} ms")
    elapsed = time.time() - t_start
    if verbose:
        print(f"[precompute_table_recurrence] total {elapsed:.2f} s "
              f"for {n_total} cells.")
    return table
