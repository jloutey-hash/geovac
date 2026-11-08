"""
Sprint Q5'-Stage1-Sub-Sprint-2b — Formal Mellin transform of the
JLO degree-0 cochain on the truncated Camporesi-Higuchi spectral triple
at n_max = 2 (and n_max = 3).

Goal
----
Compute the formal Mellin transform of

    phi_0^odd(1; t) = Tr(e^{-t D^2})

against t^{s-1} explicitly in Q[Gamma(s)] at finite cutoff n_max in {2, 3}.
Identify the s -> 0 residue structure and the explicit pi-coefficient
appearance — which is where the M1 Hopf-base-measure mechanism injects
in the master Mellin engine (Paper 18 §III.7, Paper 32 §VIII, Paper 55 §3).

Method
------
At finite n_max, the Dirac D acts on a finite-dimensional Hilbert space
H of dim N_{n_max}. The eigenvalues of D^2 are bit-exact rationals (CH
spectrum is half-integer; Dirac perturbation D = Lambda + kappa A modifies
them slightly but the structural symmetry preserves rationality of D^2
moments via Tr(D^{2j})).

For each eigenvalue mu > 0 of D^2, the Mellin transform of e^{-t mu}
against t^{s-1} is

    integral_0^inf t^{s-1} e^{-t mu} dt = Gamma(s) * mu^{-s}     (Re s > 0)

Hence by linearity,

    M[phi_0^odd(1; t)](s) = sum_{eigenvalues mu_i of D^2} Gamma(s) * mu_i^{-s}
                          = Gamma(s) * zeta_{D^2}(s)

where zeta_{D^2}(s) = sum_i mu_i^{-s} is the FINITE-CUTOFF spectral zeta
function (literal sum, no analytic continuation needed because finite).

For the diagonal CH Lambda at n_max = 2:
  eigenvalues of Lambda^2: (3/2)^2 = 9/4 with multiplicity 4 (sectors
  (1,0), (1,1) at n=1, 4 states each chirality; sum: 4 states total at n=1
  by mult 2*n*(n+1) = 4, then... wait, dim=16, with n=1 mult 4 and n=2
  mult 12, so mu = 9/4 with mult 4 and mu = 25/4 with mult 12).

  Check: Tr(Lambda^0) = 4 + 12 = 16. OK.
  Tr(Lambda^2) = 4 * 9/4 + 12 * 25/4 = 9 + 75 = 84. OK.

The full Dirac D = Lambda + kappa A has perturbed eigenvalues; we use the
CH-1 trace data (Tr(D^{2j}) bit-exact) for the j-moment expansion, but
for the Mellin transform we need the actual eigenvalues of D^2.

Approach: build D^2 explicitly as a sympy Rational matrix, diagonalize
(eigenvalues over algebraic closure), then form the spectral zeta as a
sum over distinct eigenvalues with multiplicities.

If diagonalization is intractable in exact form for the full D, we work
with the diagonal Lambda case (which IS bit-exact and clean) and note
the kappa A perturbation as a structural extension.

The Mellin transform at integer s = 1, 2, 3, ...
-----------------------------------------------
At positive integer s, Gamma(s) = (s-1)! is rational, and zeta_{D^2}(s)
is a sum of inverse-rational eigenvalues raised to integer powers, so
M[phi_0](s) is bit-exact rational at finite n_max.

At s -> 0:
  Gamma(s) ~ 1/s - gamma_Euler + O(s)
  zeta_{D^2}(s) -> zeta_{D^2}(0) = (sum of mu_i^0) = N (the dim of H!)

So at s -> 0:
  M[phi_0](s) ~ Gamma(s) * zeta_{D^2}(0) = Gamma(s) * N

  This means the residue at s = 0 is exactly N = dim H (a bit-exact integer).
  And in the Laurent expansion around s = 0:
    M[phi_0](s) = N/s + N * (-gamma_E) + O(s)

The M1 mechanism (Hopf-base measure pi = Vol(S^2)/4) appears not at this
finite-cutoff stage but in the CONTINUUM LIMIT n_max -> infty, where
zeta_{D^2}(s) at fixed s acquires a pole structure (and pi-power values
at integer s, the M2 ring), and at s -> 0 the regularised dimension
relates to the heat-kernel zero-mode count, which on a continuum S^3
introduces volume factors of (2 pi^2) (Vol(S^3) = 2 pi^2).

So the s -> 0 RESIDUE structure has the form:
  Finite cutoff: residue = N (rational integer, no pi).
  Continuum limit (n_max -> infty): residue = (regularised zeta(0))
    which has a structural pi-coefficient through the Weyl law.

This sprint computes the finite-cutoff version exactly and explicitly,
and identifies HOW the continuum-limit pi appears (through Vol(S^3) and
the spectral-action zeroth coefficient).

Discipline
----------
- bit-exact sympy.Rational throughout for the finite-cutoff zeta values.
- symbolic sympy.gamma for Gamma(s).
- no PSLQ.
- transcendentals (pi, gamma_E) appear only at the s -> 0 residue stage
  and at the continuum-limit comparison.

Output
------
- debug/data/sprint_q5p_2b_phi0_mellin_data.json
- debug/sprint_q5p_2b_phi0_mellin_memo.md

Author: PM session, 2026-06-05.
"""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Symbol, gamma, pi, zeta, factorial, simplify
from sympy import N as numeric

from geovac.spectral_triple import FockSpectralTriple


# =====================================================================
# Eigenvalue extraction for Lambda^2 and D^2 at small n_max
# =====================================================================


def lambda_squared_spectrum(n_max: int) -> List[Tuple[Rational, int]]:
    """Eigenvalues of Lambda^2 with multiplicities at given n_max.

    The diagonal CH Dirac Lambda has eigenvalues +/- (n + 1/2) for
    n = 1, ..., n_max, each with multiplicity 2n(n+1) split equally
    between + and -. So Lambda^2 has eigenvalue (n + 1/2)^2 with
    multiplicity 2n(n+1).

    Returns
    -------
    List of (mu, mult) pairs, mu = (n+1/2)^2 in sympy.Rational.
    """
    spectrum = []
    for n in range(1, n_max + 1):
        mu = Rational(2 * n + 1, 2) ** 2
        mult = 2 * n * (n + 1)
        spectrum.append((mu, mult))
    return spectrum


def D_squared_spectrum_exact(n_max: int) -> List[Tuple[sp.Expr, int]]:
    """Exact eigenvalues of D^2 = (Lambda + kappa A)^2 at given n_max.

    Builds D from FockSpectralTriple and computes eigenvalues via
    sympy. The CH-1 structure shows Tr(D^{2j}) is rational with
    denominators powers of 2 (powers of 16 from kappa = -1/16).

    For full D at n_max = 2, sympy eigvals returns degree-4 Galois roots
    involving complex cube roots — algebraic but messy. We do NOT use this
    route. The cleaner route for the Mellin transform at positive integer s
    is to compute zeta_{D^2}(s) = Tr((D^2)^{-s}) directly as a rational
    matrix trace; see finite_cutoff_zeta_via_matrix_power().

    Returns
    -------
    List of (mu, mult) pairs.
    """
    st = FockSpectralTriple(n_max=n_max)
    D = st.dirac_operator
    D2 = D * D
    eigvals_dict = D2.eigenvals()
    spectrum = [(mu, mult) for mu, mult in eigvals_dict.items()]
    return spectrum


def finite_cutoff_zeta_via_matrix_power(n_max: int, s_val: int) -> sp.Rational:
    """Compute zeta_{D^2}(s) = Tr((D^2)^{-s}) for positive integer s
    by direct matrix inversion. Bit-exact rational for full D = Lambda
    + kappa A.

    This bypasses eigenvalue extraction and gives a clean rational result
    for any positive integer s. Internally:
      Tr((D^2)^{-s}) = Tr( (D^{-2})^s ) = Tr( (D^2)^{-1} . (D^2)^{-1} . ... )
    sympy matrix inversion of a rational matrix gives a rational matrix.
    """
    assert s_val >= 1, "matrix-power route requires positive integer s"
    st = FockSpectralTriple(n_max=n_max)
    D = st.dirac_operator
    D2 = D * D
    D2_inv = D2.inv()
    power = D2_inv ** s_val
    tr = sp.trace(power)
    return sp.simplify(tr)


def finite_cutoff_log_det_D2_via_charpoly(n_max: int) -> sp.Expr:
    """Compute sum_i m_i log(mu_i) = log det(D^2) bit-exactly via the
    characteristic polynomial.

    log det(D^2) = log(prod mu_i) = log(det D^2).
    sympy returns det(D^2) as an exact rational, so log det = log(rational).

    For full D at n_max = 2, this gives a single log of an exact rational
    (= log(numerator) - log(denominator)) without needing eigenvalues.
    """
    st = FockSpectralTriple(n_max=n_max)
    D = st.dirac_operator
    D2 = D * D
    det_D2 = D2.det()
    return sp.log(det_D2)


# =====================================================================
# Finite-cutoff spectral zeta and Mellin transform
# =====================================================================


def finite_cutoff_zeta(
    spectrum: List[Tuple[sp.Expr, int]],
    s_val,
) -> sp.Expr:
    """Spectral zeta at finite cutoff: zeta_{D^2}(s_val) = sum mu^{-s} * mult.

    Excludes zero eigenvalues (if any). Returns exact sympy expression.

    s_val: integer (for bit-exact rational at finite cutoff) OR sympy Symbol
    (for symbolic expansion).
    """
    if isinstance(s_val, int):
        s = Rational(s_val)
    else:
        s = s_val
    total = sp.Integer(0)
    for mu, mult in spectrum:
        if mu == 0:
            continue
        total += Integer(mult) * mu ** (-s)
    return simplify(total)


def mellin_transform_phi0(
    spectrum: List[Tuple[sp.Expr, int]],
    s_val,
) -> sp.Expr:
    """Formal Mellin transform M[Tr(e^{-t D^2})](s) at integer s_val.

    M[Tr(e^{-t D^2})](s) = Gamma(s) * zeta_{D^2}(s)

    at finite cutoff, this is bit-exact rational times symbolic Gamma(s).
    """
    s = Symbol("s") if not isinstance(s_val, int) else Rational(s_val)
    zeta_val = finite_cutoff_zeta(spectrum, s_val)
    if isinstance(s_val, int):
        gamma_val = factorial(s_val - 1) if s_val >= 1 else gamma(s)
        return simplify(gamma_val * zeta_val)
    return gamma(s) * zeta_val


def laurent_expansion_at_zero(
    spectrum: List[Tuple[sp.Expr, int]],
    n_terms: int = 4,
) -> Dict:
    """Laurent expansion of M[phi_0](s) at s = 0, up to O(s^{n_terms-1}).

    Around s = 0:
      Gamma(s) = 1/s - gamma_E + (gamma_E^2/2 + pi^2/12) s + ...
      zeta_{D^2}(s) at finite cutoff is analytic at s = 0:
        zeta_{D^2}(s) = N - s log(prod mu_i^{m_i}) + O(s^2)
                     = N - s sum_i m_i log(mu_i) + O(s^2)
      where N = sum_i m_i = dim H.

    So:
      M[phi_0](s) * s = (1 - gamma_E s + ...) * (N + s * (- sum m log mu) + ...)
                     = N + s (-gamma_E N - sum m log mu) + ...

    Hence:
      residue at s = 0 = N (the dim H)
      next coefficient = -gamma_E * N - sum_i m_i log(mu_i)

    The sum sum_i m_i log(mu_i) is the log-determinant of D^2 (the
    spectral-zeta derivative at zero, up to sign):
      log det D^2 = -zeta'_{D^2}(0)
    i.e., -zeta'_{D^2}(0) = sum_i m_i log(mu_i).

    Returns dict with: residue, next-order coefficient (in terms of gamma_E,
    log mu values), and the explicit pi-content (none at finite cutoff).
    """
    s = Symbol("s", positive=True)
    # Build symbolic zeta_{D^2}(s) and expand around s = 0
    zeta_sym = finite_cutoff_zeta(spectrum, s)
    # Compute zeta(0) and zeta'(0) exactly (finite sum)
    N = sum(mult for _, mult in spectrum)  # = dim H
    log_det_D2_neg = sum(
        Integer(mult) * sp.log(mu) for mu, mult in spectrum if mu != 0
    )
    # zeta'_{D^2}(0) = -log det D^2 = -log_det_D2_neg = -(sum m log mu)
    zeta_prime_at_zero = -log_det_D2_neg

    # Laurent: Gamma(s) zeta(s) = (1/s) * zeta(s) * s * Gamma(s)
    # s * Gamma(s) = Gamma(s + 1) (sympy convention check)
    # In fact s * Gamma(s) = Gamma(s+1), so
    # Gamma(s) = Gamma(s+1) / s, and Gamma(s+1) is regular at 0 with value 1.
    # Expanding Gamma(s+1) = 1 - gamma_E s + (gamma_E^2/2 + pi^2/12) s^2 - ...
    gamma_E = sp.EulerGamma

    # Series of Gamma(s+1) around s=0:
    g_series = [
        Integer(1),
        -gamma_E,
        gamma_E ** 2 / 2 + pi ** 2 / 12,
        -gamma_E ** 3 / 6 - gamma_E * pi ** 2 / 12 - sp.zeta(3) / 3,
    ]
    # zeta_{D^2}(s) series: zeta_{D^2}(0) = N, zeta'(0) = -log det,
    # higher: need to compute via series expansion of sum mu^{-s} = sum exp(-s log mu)
    # zeta_{D^2}(s) = sum_i m_i exp(-s log mu_i)
    #              = sum_i m_i (1 - s log mu_i + (s log mu_i)^2/2 - ...)
    #              = N + s * (- sum m log mu) + s^2/2 * (sum m (log mu)^2) - ...

    log_powers = []
    for k in range(n_terms + 2):
        coef = sum(
            Integer(mult) * (-sp.log(mu)) ** k / factorial(k)
            for mu, mult in spectrum
            if mu != 0
        )
        log_powers.append(simplify(coef))
    # log_powers[k] = coefficient of s^k in zeta_{D^2}(s) series at s=0

    # Now: Gamma(s) zeta(s) = (1/s) * (Gamma(s+1) * zeta(s))
    # Compute the product series (Gamma(s+1) * zeta(s)) to n_terms + 1 powers
    n_keep = n_terms + 1
    product_series = []
    for n in range(n_keep + 1):
        c = sum(
            (g_series[i] if i < len(g_series) else Integer(0))
            * log_powers[n - i]
            for i in range(n + 1)
        )
        product_series.append(simplify(c))

    # Laurent: M[phi_0](s) = (1/s) * sum_n product_series[n] * s^n
    # Coefficient of s^{n-1} is product_series[n].
    laurent_coeffs = product_series  # idx k -> coeff of s^{k-1}

    result = {
        "dim_H": N,
        "residue_at_s_zero": product_series[0],
        "next_order_coeff": product_series[1] if len(product_series) > 1 else None,
        "log_det_D_squared": -log_det_D2_neg if log_det_D2_neg != 0 else Integer(0),
        "zeta_prime_at_zero": zeta_prime_at_zero,
        "laurent_around_zero": product_series,
        "log_powers_of_zeta": log_powers,
    }
    return result


# =====================================================================
# Continuum-limit comparison: M1 mechanism appears in continuum spectral zeta
# =====================================================================


def continuum_spectral_zeta_at_integer(s_val: int) -> sp.Expr:
    """Continuum CH spectral zeta zeta_{D^2}(s) at positive integer s
    (Hurwitz-zeta reduction, same as Sprint Q5'-CH-2).

    For s = 0, the literal sum diverges; the regularised value at the
    continuum is obtained via analytic continuation and at s = 0 on a
    smooth 3-manifold equals (Weyl-law leading) zeta(0) = ... (specific
    closed forms).
    """
    s = Rational(s_val)
    if s == 0:
        # Direct evaluation of the regularised value zeta_{D^2}(0):
        # Use Hurwitz form: zeta_{D^2}(s) = 2 zeta(2s-2, 3/2) - (1/2) zeta(2s, 3/2)
        # At s = 0: 2 * zeta(-2, 3/2) - (1/2) * zeta(0, 3/2)
        # zeta(0, a) = 1/2 - a; zeta(-2, a) = -B_3(a)/3 where B_3 is the Bernoulli polynomial
        # B_3(a) = a^3 - (3/2) a^2 + (1/2) a
        # B_3(3/2) = 27/8 - 27/8 + 3/4 = 3/4
        # So zeta(-2, 3/2) = -(3/4)/3 = -1/4
        # 2 * (-1/4) - (1/2) * (1/2 - 3/2) = -1/2 + 1/2 = 0
        # i.e., zeta_{D^2}(0) = 0 in the continuum. Interesting: the regularised
        # dimension on continuum CH-S^3 is 0, NOT infinite (this is the
        # zeta-regularised statement of the Weyl law on a 3-manifold for
        # the Laplace-Beltrami; D^2 of CH-Dirac has same structure).
        two_s_minus_2 = -2
        two_s = 0
        # zeta(-2, 3/2):
        # B_3(3/2) = (3/2)^3 - 3*(3/2)^2/2 + (3/2)/2 = 27/8 - 27/8 + 3/4 = 3/4
        # zeta(-2, a) = -B_3(a)/3
        z_low = Rational(-3, 4) / Integer(3)  # = -1/4
        # zeta(0, 3/2) = 1/2 - 3/2 = -1
        z_high = Rational(1, 2) - Rational(3, 2)  # = -1
        return simplify(Integer(2) * z_low - Rational(1, 2) * z_high)
    # Use the same Hurwitz reduction as compute_ch_spectral_zeta_continuum.py
    two_s_minus_2 = 2 * s - 2
    two_s = 2 * s
    if two_s_minus_2 == 0:
        # s = 1: zeta(0, 3/2) = -1
        hurwitz_low = Rational(1, 2) - Rational(3, 2)
    elif two_s_minus_2 < 0:
        hurwitz_low = sp.zeta(two_s_minus_2, Rational(3, 2))
    else:
        # zeta(s, 1/2) = (2^s - 1) zeta(s); zeta(s, 3/2) = zeta(s, 1/2) - 2^s
        z_half = (2 ** two_s_minus_2 - 1) * sp.zeta(two_s_minus_2)
        hurwitz_low = z_half - 2 ** two_s_minus_2
    z_half_high = (2 ** two_s - 1) * sp.zeta(two_s)
    hurwitz_high = z_half_high - 2 ** two_s
    expr = 2 * hurwitz_low - Rational(1, 2) * hurwitz_high
    return simplify(expr)


def continuum_mellin_transform(s_val: int) -> sp.Expr:
    """Continuum M[Tr(e^{-t D^2})](s) = Gamma(s) * zeta_{D^2}^{continuum}(s)."""
    if s_val == 0:
        return None  # pole / s -> 0 expansion needed
    z = continuum_spectral_zeta_at_integer(s_val)
    g = factorial(s_val - 1) if s_val >= 1 else gamma(Rational(s_val))
    return simplify(g * z)


def hopf_base_measure_identification() -> Dict:
    """Identify where M1 pi = Vol(S^2)/4 injects in the Mellin transform.

    Sources (from Paper 18 §III.7, Paper 55 §3, Paper 25 §II):
      - Vol(S^3) = 2 pi^2
      - Vol(S^2) = 4 pi
      - Hopf-base measure factor: Vol(S^2)/4 = pi  [M1 fingerprint]
      - 4/pi = Vol(S^2)/pi^2 = 2 Vol(S^1)/Vol(SU(2))  [L2 propinquity rate]

    For the continuum spectral zeta zeta_{D^2}^{cont}(s):
      - Each value at integer s is in the M2 ring oplus_k pi^{2k} . Q
      - The Mellin transform Gamma(s) * zeta(s) at integer s is rational
        * (M2 ring element), since Gamma(integer) is integer.

    For the s -> 0 expansion (regularised dimension):
      - In the CONTINUUM: zeta_{D^2}(0) = 0 (computed above)
        and zeta'_{D^2}(0) = some closed form involving log 2, pi, gamma_E
        (the log-determinant on a continuum compact 3-manifold).
      - In the FINITE CUTOFF: zeta_{D^2}(0) = N = dim H (rational integer)
        and zeta'_{D^2}(0) = -sum m log mu (transcendentals from logs, but
        no pi).

    So at the FINITE-CUTOFF stage of the Mellin transform, the M1
    pi-coefficient is ABSENT. It appears only after the continuum limit
    AND at the s -> 0 regularised dimension stage, through the Vol(S^3)
    Weyl law leading.

    Returns dict with the structural identification.
    """
    Vol_S2 = 4 * pi
    Vol_S3 = 2 * pi ** 2
    Hopf_base_measure = Vol_S2 / 4  # = pi (M1 leading)
    # Weyl law leading on continuum 3-manifold:
    # Tr(e^{-t Delta}) ~ Vol/(4 pi t)^{3/2} as t -> 0+
    # for Laplace-Beltrami; for Dirac D^2 same leading. Multiplicity of
    # spinor bundle on S^3 = 2 (Camporesi-Higuchi rank).
    spinor_mult = Integer(2)
    weyl_leading_coeff = spinor_mult * Vol_S3 / (4 * pi) ** Rational(3, 2)
    return {
        "Vol_S3": str(Vol_S3),
        "Vol_S2": str(Vol_S2),
        "M1_Hopf_base_measure": str(Hopf_base_measure),
        "spinor_multiplicity": str(spinor_mult),
        "weyl_leading_coeff_continuum": str(weyl_leading_coeff),
        "weyl_leading_simplified": str(simplify(weyl_leading_coeff)),
        "interpretation": (
            "At finite cutoff, M[phi_0](s) at s in {1,2,3} is bit-exact "
            "rational * Gamma(integer); no pi appears. The s -> 0 residue is "
            "N = dim H (rational integer). M1 Hopf-base measure pi enters only "
            "at the continuum limit, through the Weyl-law leading coefficient "
            "Vol(S^3) = 2 pi^2 which sets the heat-kernel small-t leading "
            "as t^{-3/2} * pi^{-1/2}. The M1 mechanism is therefore the "
            "continuum-limit M[Tr(e^{-t D^2})] short-time asymptotic, "
            "NOT the finite-cutoff Laurent expansion at s = 0."
        ),
    }


# =====================================================================
# Sub-Sprint 1 cross-check: small-t Taylor expansion vs Mellin direct
# =====================================================================


def small_t_taylor_phi0(spectrum: List[Tuple[sp.Expr, int]], M_max: int) -> List[sp.Expr]:
    """Small-t Taylor expansion of Tr(e^{-t D^2}) = sum_i m_i e^{-t mu_i}.

    c_m = (-1)^m / m! * sum_i m_i mu_i^m

    Returns list [c_0, c_1, ..., c_{M_max}].
    """
    coeffs = []
    for m in range(M_max + 1):
        c = (-Integer(1)) ** m * sum(
            Integer(mult) * mu ** m for mu, mult in spectrum
        ) / factorial(m)
        coeffs.append(simplify(c))
    return coeffs


# =====================================================================
# Driver
# =====================================================================


def run_sprint() -> Dict:
    out: Dict = {
        "sprint": "Q5'-Stage1-Sub-Sprint-2b",
        "title": "Formal Mellin transform of phi_0^odd at n_max in {2, 3}",
        "discipline": "bit-exact sympy.Rational; symbolic Gamma(s)",
        "results": {},
    }

    # ----------------------- n_max = 2, Lambda only -----------------------
    print("=" * 72)
    print("Sub-Sprint 2b — n_max = 2, diagonal Lambda only (bit-exact)")
    print("=" * 72)
    spec_L2 = lambda_squared_spectrum(2)
    print(f"\nLambda^2 spectrum at n_max=2: {[(str(mu), m) for mu, m in spec_L2]}")
    dim_H_2 = sum(m for _, m in spec_L2)
    print(f"dim H = {dim_H_2}")

    # Verify small-t Taylor matches Sub-Sprint 1 phi_0^odd(1; t)
    taylor_L2 = small_t_taylor_phi0(spec_L2, M_max=4)
    print(f"\nSmall-t Taylor of Tr(e^{{-t Lambda^2}}) on Lambda (n_max=2):")
    for m, c in enumerate(taylor_L2):
        print(f"  c_{m} = {c}")
    expected_sub1 = ["16", "-84", "489/2", "-3967/8", "98203/128"]
    print("\nCross-check vs Sub-Sprint 1 phi_0^odd(1; t) on Lambda:")
    cross_check = []
    for m, c in enumerate(taylor_L2):
        match = str(c) == expected_sub1[m]
        cross_check.append({"m": m, "computed": str(c), "sub1": expected_sub1[m], "match": match})
        print(f"  m={m}: computed={c}, sub1={expected_sub1[m]}, match={match}")

    out["results"]["nmax2_Lambda"] = {
        "dim_H": dim_H_2,
        "spectrum": [(str(mu), m) for mu, m in spec_L2],
        "small_t_taylor": [str(c) for c in taylor_L2],
        "cross_check_sub_sprint_1": cross_check,
    }

    # Compute Mellin transform at integer s
    print("\nMellin transform M[phi_0^odd(1; t)](s) = Gamma(s) zeta_{Lambda^2}(s)")
    print("at finite cutoff n_max = 2:\n")
    mellin_at_int = {}
    for s_val in [1, 2, 3, 4, 5]:
        zeta_val = finite_cutoff_zeta(spec_L2, s_val)
        gamma_val = factorial(s_val - 1)
        mellin = simplify(gamma_val * zeta_val)
        print(f"  s = {s_val}:")
        print(f"    zeta_{{Lambda^2}}({s_val}) = {zeta_val}")
        print(f"    Gamma({s_val}) = {gamma_val}")
        print(f"    M[phi_0]({s_val}) = {mellin}")
        # Numeric for context
        zeta_n = float(numeric(zeta_val, 20))
        mellin_n = float(numeric(mellin, 20))
        # Compare to continuum
        try:
            cont_z = continuum_spectral_zeta_at_integer(s_val)
            cont_z_n = float(numeric(cont_z, 20))
            print(f"    continuum zeta = {cont_z}  (num: {cont_z_n:.6f}); ratio trunc/cont = {zeta_n/cont_z_n:.4f}")
        except Exception as e:
            cont_z = None
            cont_z_n = None
            print(f"    continuum zeta: error {e}")
        mellin_at_int[str(s_val)] = {
            "Gamma_s": str(gamma_val),
            "zeta_D2_truncated": str(zeta_val),
            "zeta_D2_truncated_numeric": zeta_n,
            "M_phi0_truncated": str(mellin),
            "M_phi0_truncated_numeric": mellin_n,
            "continuum_zeta_D2": str(cont_z) if cont_z is not None else None,
            "continuum_zeta_D2_numeric": cont_z_n,
        }

    out["results"]["nmax2_Lambda"]["mellin_at_integer_s"] = mellin_at_int

    # s -> 0 Laurent expansion
    print("\ns -> 0 Laurent expansion of M[phi_0^odd(1; t)](s) at n_max = 2, Lambda:")
    laurent_2 = laurent_expansion_at_zero(spec_L2, n_terms=4)
    print(f"  dim H = {laurent_2['dim_H']}")
    print(f"  Residue at s = 0 (coeff of 1/s): {laurent_2['residue_at_s_zero']}")
    print(f"  Next-order coeff (s^0): {laurent_2['next_order_coeff']}")
    print(f"  log det D^2 = sum m log mu = {laurent_2['log_det_D_squared']}")
    print(f"  Numeric log det = {float(numeric(laurent_2['log_det_D_squared'], 20))}")
    print(f"  zeta'(0) = -log det = {laurent_2['zeta_prime_at_zero']}")

    out["results"]["nmax2_Lambda"]["laurent_at_zero"] = {
        "dim_H": int(laurent_2["dim_H"]),
        "residue": str(laurent_2["residue_at_s_zero"]),
        "next_order": str(laurent_2["next_order_coeff"]),
        "log_det_D_squared": str(laurent_2["log_det_D_squared"]),
        "log_det_D_squared_numeric": float(numeric(laurent_2["log_det_D_squared"], 30)),
        "zeta_prime_at_zero": str(laurent_2["zeta_prime_at_zero"]),
        "laurent_coefficients_around_zero": [str(c) for c in laurent_2["laurent_around_zero"]],
    }

    # ----------------------- n_max = 2, full D = Lambda + kappa A -----------------------
    print("\n" + "=" * 72)
    print("Sub-Sprint 2b — n_max = 2, full D = Lambda + kappa A (bit-exact)")
    print("=" * 72)

    # The full D = Lambda + kappa A has eigenvalues that are algebraic
    # but ugly (degree-4 irreducible factors over Q, Galois roots with
    # complex cube roots, recalling Sprint 3 RH-F S_6 non-solvable Galois
    # finding for related polynomials). Direct closed-form eigenvalues are
    # not the right route.
    #
    # Cleanest: zeta_{D^2}(s) at positive integer s is a bit-exact rational
    # via Tr((D^2)^{-s}) (matrix inversion + power, all entries rational).
    # log det(D^2) is a single log of a rational (the determinant).
    try:
        dim_H_2_full = 16
        mellin_int_full = {}
        for s_val in [1, 2, 3]:
            t0 = time.time()
            zeta_val = finite_cutoff_zeta_via_matrix_power(2, s_val)
            gamma_val = factorial(s_val - 1)
            mellin = simplify(gamma_val * zeta_val)
            t_one = time.time() - t0
            print(f"\n  s = {s_val} (computed in {t_one:.2f}s):")
            print(f"    zeta_{{D^2}}({s_val}) = {zeta_val}")
            print(f"    Gamma({s_val}) = {gamma_val}")
            print(f"    M[phi_0]({s_val}) = {mellin}")
            zeta_n = float(numeric(zeta_val, 30))
            mellin_n = float(numeric(mellin, 30))
            mellin_int_full[str(s_val)] = {
                "Gamma_s": str(gamma_val),
                "zeta_D2_truncated": str(zeta_val),
                "zeta_D2_truncated_numeric": zeta_n,
                "M_phi0_truncated": str(mellin),
                "M_phi0_truncated_numeric": mellin_n,
                "matrix_power_seconds": t_one,
            }

        t0 = time.time()
        log_det_full = finite_cutoff_log_det_D2_via_charpoly(2)
        t_logdet = time.time() - t0
        log_det_full_n = float(numeric(log_det_full, 30))
        print(f"\n  log det(D^2) via det matrix (computed in {t_logdet:.2f}s):")
        print(f"    log det(D^2) = {log_det_full}")
        print(f"    log det(D^2) numeric: {log_det_full_n}")

        # Laurent at s = 0:
        # M[phi_0](s) = (1/s) * Gamma(s+1) * zeta_{D^2}(s)
        # Gamma(s+1) = 1 - gamma_E s + (gamma_E^2/2 + pi^2/12) s^2 - ...
        # zeta_{D^2}(s) = N - s log det(D^2) + (s^2/2)(sum m (log mu)^2) - ...
        # So:
        # M[phi_0](s) = (1/s) [ N + s * (-log det - N gamma_E) + O(s^2) ]
        #             = N/s + (-log det - N gamma_E) + O(s)
        gamma_E = sp.EulerGamma
        residue = Integer(dim_H_2_full)
        next_order = -log_det_full - Integer(dim_H_2_full) * gamma_E
        print(f"\n  s -> 0 Laurent (full D, n_max=2):")
        print(f"    Residue at s = 0 (coeff of 1/s): {residue}")
        print(f"    Next-order coeff (s^0): {next_order}")
        print(f"    Next-order numeric: {float(numeric(next_order, 30))}")

        out["results"]["nmax2_full_D"] = {
            "dim_H": dim_H_2_full,
            "mellin_at_integer_s": mellin_int_full,
            "log_det_D_squared": str(log_det_full),
            "log_det_D_squared_numeric": log_det_full_n,
            "log_det_compute_seconds": t_logdet,
            "laurent_at_zero": {
                "residue": str(residue),
                "next_order_coeff": str(next_order),
                "next_order_numeric": float(numeric(next_order, 30)),
                "note": "(1/s) . Gamma(s+1) . zeta_D2(s); pole = N = dim_H; next = -log det - N gamma_E",
            },
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"\n  Full D Mellin at n_max=2 FAILED: {e}")
        out["results"]["nmax2_full_D"] = {"error": str(e)}

    # ----------------------- n_max = 3, Lambda only -----------------------
    print("\n" + "=" * 72)
    print("Sub-Sprint 2b — n_max = 3, diagonal Lambda only (bit-exact)")
    print("=" * 72)
    spec_L3 = lambda_squared_spectrum(3)
    print(f"\nLambda^2 spectrum at n_max=3: {[(str(mu), m) for mu, m in spec_L3]}")
    dim_H_3 = sum(m for _, m in spec_L3)
    print(f"dim H = {dim_H_3}")

    taylor_L3 = small_t_taylor_phi0(spec_L3, M_max=4)
    print("\nSmall-t Taylor of Tr(e^{{-t Lambda^2}}) on Lambda (n_max=3):")
    for m, c in enumerate(taylor_L3):
        print(f"  c_{m} = {c}")

    mellin_int_L3 = {}
    print("\nMellin transform M[phi_0](s) at integer s (n_max=3, Lambda):")
    for s_val in [1, 2, 3, 4, 5]:
        zeta_val = finite_cutoff_zeta(spec_L3, s_val)
        gamma_val = factorial(s_val - 1)
        mellin = simplify(gamma_val * zeta_val)
        print(f"  s = {s_val}: zeta = {zeta_val}, M[phi_0] = {mellin}")
        cont_z = continuum_spectral_zeta_at_integer(s_val)
        ratio = float(numeric(zeta_val / cont_z, 20)) if cont_z != 0 else None
        mellin_int_L3[str(s_val)] = {
            "Gamma_s": str(gamma_val),
            "zeta_D2_truncated": str(zeta_val),
            "M_phi0_truncated": str(mellin),
            "continuum_zeta": str(cont_z),
            "ratio_trunc_to_cont": ratio,
        }

    laurent_3 = laurent_expansion_at_zero(spec_L3, n_terms=3)
    print("\ns -> 0 Laurent (n_max=3, Lambda):")
    print(f"  dim H = {laurent_3['dim_H']}")
    print(f"  Residue at s = 0: {laurent_3['residue_at_s_zero']}")
    print(f"  log det Lambda^2 (numeric): {float(numeric(laurent_3['log_det_D_squared'], 30))}")

    out["results"]["nmax3_Lambda"] = {
        "dim_H": dim_H_3,
        "spectrum": [(str(mu), m) for mu, m in spec_L3],
        "small_t_taylor": [str(c) for c in taylor_L3],
        "mellin_at_integer_s": mellin_int_L3,
        "laurent_at_zero": {
            "dim_H": int(laurent_3["dim_H"]),
            "residue": str(laurent_3["residue_at_s_zero"]),
            "next_order": str(laurent_3["next_order_coeff"]),
            "log_det_numeric": float(numeric(laurent_3["log_det_D_squared"], 30)),
        },
    }

    # ----------------------- Continuum-limit M1 identification -----------------------
    print("\n" + "=" * 72)
    print("Continuum-limit M1 Hopf-base-measure identification")
    print("=" * 72)
    m1_id = hopf_base_measure_identification()
    for k, v in m1_id.items():
        print(f"\n  {k}: {v}")
    out["m1_identification"] = m1_id

    # Continuum spectral-zeta value at s = 0 (regularised dimension)
    # zeta_{D^2}^{cont}(0) — computed above to be 0 (exact continuum).
    print("\n" + "=" * 72)
    print("Continuum regularised dimension zeta_{D^2}^{cont}(0)")
    print("=" * 72)
    cont_zero = continuum_spectral_zeta_at_integer(0)
    print(f"\n  zeta_{{D^2}}^{{cont}}(0) = {cont_zero}")
    # This is the regularised dimension in the continuum limit.
    # In the FINITE CUTOFF, zeta_{D^2}(0) = dim H = N (a rational integer).
    # In the CONTINUUM (n_max -> infty), zeta_{D^2}^{cont}(0) = 0.
    # The DIFFERENCE is the meaningful "anomaly":
    #   anomaly = lim_{n_max -> infty} [zeta_{n_max}(0) - zeta_{cont}(0)]
    #          = N - 0 = infty
    # which is the Weyl-law divergence of dim H.
    # The Hopf-base measure pi enters via the Weyl-law coefficient
    # 2 Vol(S^3)/(4 pi)^{3/2} setting the t^{-3/2} short-time leading.
    out["continuum_results"] = {
        "regularised_zeta_at_s_zero": str(cont_zero),
        "comment": (
            "Continuum regularised zeta_{D^2}(0) = 0 exactly (computed via "
            "Hurwitz Bernoulli polynomial reduction). The finite-cutoff value "
            "is dim H. The Weyl-law divergence dim H -> infinity matches the "
            "continuum t^{-3/2} short-time leading of Tr(e^{-t D^2})."
        ),
    }

    return out


def main() -> None:
    t_start = time.time()
    out = run_sprint()
    out["wall_seconds"] = time.time() - t_start

    out_path = Path("debug/data/sprint_q5p_2b_phi0_mellin_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\n\nOutput written: {out_path}")
    print(f"Total wall: {out['wall_seconds']:.2f} s")


if __name__ == "__main__":
    main()
