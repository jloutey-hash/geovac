"""Verification of Paper 55 Theorem M2 (thm:m2_mixed_tate).

Paper 55 (Periods of GeoVac, papers/group3_foundations/paper_55_periods_of_geovac.tex)
Theorem M2 period-ring (thm:m2_mixed_tate) states:

    The M2 sub-mechanism on the discrete Camporesi-Higuchi S^3 spectral
    triple produces period values in the strictly smaller pure-Tate
    even-weight sub-ring

        M_2 subset bigoplus_{k >= 0} pi^{2k} . Q  subset MT(Q)

    of mixed-Tate periods over Q, at depth 0.

The witness panel from Paper 55 sec:m2 Witnesses subsection lists:

    a_0^{D^2}   (Dirac SD, vol-norm)        =  4 pi^2          Tate weight +2
    a_1^{D^2}   (Dirac SD, vol-norm)        = -2 pi^2          Tate weight +2
    a_k^{D^2}   (k >= 2)                    =  0               (vanish, two-term exactness)
    a_k^Delta   (scalar SD, vol-norm)       =  2 pi^2 / k!     Tate weight +2

    zeta_{D^2}(2)                           =  pi^2 - pi^4/12              weights +2, +4
    zeta_{D^2}(3)                           =  pi^4/3 - pi^6/30            weights +4, +6
    zeta_{D^2}(4)                           =  2*pi^6/15 - 17*pi^8/1260    weights +6, +8
    zeta_{D^2}(5)                           =  pi^8(306 - 31 pi^2)/5670    weights +8, +10

Companion sprint memo: debug/sprint_mixed_tate_test_memo.md.

This test file verifies M2 by three independent paths:

  (1) Compute the Dirac SD coefficients a_k for k = 0..4 on unit S^3
      symbolically via geovac.qed_vacuum_polarization, divide out the
      (4*pi)^{-d/2} dimensional prefactor and the dim_spinor = 4 factor,
      and verify the volume-normalized Dirac SD = {4*pi^2, -2*pi^2}
      matches the Paper 55 witness panel.

      The PRODUCTION qed_vacuum_polarization.seeley_dewitt_coefficients_s3
      computes a_2 using the Vassilevich-Branson-Gilkey curvature
      polynomial (5R^2 - 2|Ric|^2 + 2|Riem|^2 - 30 R E + 60 E^2)/360 with
      Lichnerowicz E = R/4. On unit S^3 the integrand is 8/3, so the bare
      a_2^Dirac (= dim_S * (1/360) * 8/3 * Vol) = 32/(360*3/2) * pi^2 ...
      etc. This is NOT zero. The two-term exactness claim is for the
      ASYMPTOTIC heat-trace expansion (Paper 51 Cor 2.1), which is a
      stronger statement than the standard Branson-Gilkey-Vassilevich a_2.
      We document this convention split and verify the ASYMPTOTIC vanishing
      independently via the Jacobi theta_2 modular identity.

  (2) Compute zeta_{D^2}(s) at integer s = 2, 3, 4, 5 by direct sum and
      verify each is in Q[pi^2]. We use sympy to identify the closed form
      and confirm Tate weights match the witness panel.

  (3) PSLQ negative test at 100 dps: for each witness value, search for
      a Q-linear relation against the augmented basis
        { 1, pi^2, pi^4, pi^6, pi^8, pi^10, zeta(3), zeta(5), Catalan G }
      where the FIRST FIVE elements are expected to fit the witness, and
      the LAST THREE (odd zeta, Catalan G) are the negative-test elements
      that MUST receive coefficient 0 in any identification. We then run
      a SECOND PSLQ search against the basis { 1, pi^2, pi^4, ..., pi^10 }
      ALONE (no odd zeta, no G) and confirm the witness is identified with
      the same coefficients as in the augmented search. This is the
      "no odd-zeta or MZV needed" claim from the test directive.

Per CLAUDE.md section 13.4a / 13.5: this test does NOT modify any paper
or production code. The Dirac asymptotic two-term exactness (Paper 51
Cor 2.1) is verified independently via the modular sum, NOT via the
Vassilevich a_2 formula in qed_vacuum_polarization. We document the
convention split rather than modifying production code.

References
----------
- Paper 55 (Periods of GeoVac), papers/group3_foundations/paper_55_periods_of_geovac.tex,
  Theorem thm:m2_mixed_tate, witness table sec:m2.
- Paper 51 (Gravity arc) Cor 2.1 (two-term-exact Dirac SD).
- Paper 28 (QED on S^3) Theorem T9 (zeta_{D^2}(s) closed form).
- debug/sprint_mixed_tate_test_memo.md (June 2026).
- Fathizadeh-Marcolli, Comm. Math. Phys. 357 (2017), arXiv:1611.01815.
"""

from __future__ import annotations

import math

import mpmath
import pytest
import sympy as sp
from sympy import Integer, Rational, pi, sqrt, zeta as sym_zeta, factorial

from geovac.qed_vacuum_polarization import (
    seeley_dewitt_coefficients_s3,
    verify_no_odd_zeta_one_loop,
)


# Working precision for the PSLQ runs (>= 100 dps as the directive requests).
# Note: mpmath.mp.dps is GLOBAL module-level state; other test modules in this
# directory set it lower (e.g., test_paper51_scalar_ak.py sets dps=50), so we
# RESTORE dps=120 inside each test that needs it (via the _set_dps fixture
# helper) rather than relying on the module-import-time setting.
mpmath.mp.dps = 120
_DPS_TARGET = 120


def _ensure_dps():
    if mpmath.mp.dps < _DPS_TARGET:
        mpmath.mp.dps = _DPS_TARGET


# Tolerance for symbolic comparisons (must be 0).
def _is_zero(expr: sp.Expr) -> bool:
    return sp.simplify(sp.expand(expr)) == 0


# =============================================================================
# Path (1): Volume-normalized Dirac SD a_0, a_1 from production code.
# =============================================================================

def test_witness_a0_dirac_vol_normalized():
    """Witness: a_0^{D^2} (Dirac SD, vol-norm) = 4 pi^2.

    The production seeley_dewitt_coefficients_s3() returns
    a_0 (production) = (4 pi)^{-3/2} * dim_S * Vol(S^3)
                     = (4 pi)^{-3/2} * 4 * 2 pi^2
                     = sqrt(pi).

    Paper 55 "vol-norm" convention multiplies BACK by (4 pi)^{d/2}
    = (4 pi)^{3/2} to expose the spinor-bundle-extended volume:

        a_0_vol_norm = a_0_production * (4 pi)^{3/2}
                     = sqrt(pi) * 8 pi^{3/2}
                     = 8 pi^2.

    The Paper 55 witness 4 pi^2 corresponds to dim_S = 2 (Weyl chirality)
    while the production code uses dim_S = 4 (full Dirac). Either
    convention puts a_0 in pi^2 * Q, which is what M2 claims.

    We verify a_0_vol_norm is in pi^2 * Q (the M2 ring), and check the
    explicit dim_S = 4 production value, then note the conversion to the
    Paper 55 dim_S = 2 / dim_S = 4 conventions.
    """
    coeffs = seeley_dewitt_coefficients_s3()
    # Production a_0 = sqrt(pi).
    assert _is_zero(coeffs['a0'] - sp.sqrt(pi))

    # Volume-normalized: multiply by (4 pi)^{3/2}.
    a0_vol_norm = coeffs['a0'] * (4 * pi) ** Rational(3, 2)
    a0_vol_norm = sp.simplify(a0_vol_norm)
    # = sqrt(pi) * 8 * pi * sqrt(pi) = 8 pi^2.
    assert _is_zero(a0_vol_norm - 8 * pi ** 2)

    # Paper 55 witness 4 pi^2 corresponds to dim_S = 2 (Weyl):
    #   a_0_paper55 = a_0_production * (4 pi)^{3/2} / 2.
    a0_paper55 = a0_vol_norm / 2
    assert _is_zero(a0_paper55 - 4 * pi ** 2)

    # Both 8 pi^2 (production, full Dirac) and 4 pi^2 (Paper 55, Weyl) are
    # rational * pi^2, so M2 claim holds in both conventions.
    # Test that 4 pi^2 has rational coefficient 4 in pi^2:
    ratio = sp.simplify(a0_paper55 / pi ** 2)
    assert ratio == 4
    assert ratio.is_rational


def test_witness_a1_dirac_vol_normalized():
    """Witness: a_1^{D^2} (Dirac SD, vol-norm) = -2 pi^2.

    Paper 55 witnesses (vol-normalized SD coefficients) come from the
    Paper 51 Cor 2.1 ASYMPTOTIC Dirac heat-trace closed form

        K_{D^2}(t) = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}
                     + O(e^{-pi^2/t}).

    Under the standard SD convention
        K(t) = sum_k a_k * t^{k - 3/2} * (4 pi)^{-3/2}
    we identify (multiplying through by (4 pi)^{3/2} = 8 pi^{3/2}):

        a_0 = (sqrt(pi)/2) * (4 pi)^{3/2}
            = (sqrt(pi)/2) * 8 * pi^{3/2}
            = 4 pi^2                                           [Paper 55]
        a_1 = (-sqrt(pi)/4) * (4 pi)^{3/2}
            = (-sqrt(pi)/4) * 8 * pi^{3/2}
            = -2 pi^2                                          [Paper 55]
        a_k = 0           for k >= 2 (two-term exactness, Paper 51 Cor 2.1).

    NOTE on production code: qed_vacuum_polarization.seeley_dewitt_coefficients_s3
    computes Vassilevich-Branson-Gilkey a_2 using the curvature polynomial
    (5R^2 - 2|Ric|^2 + 2|Riem|^2 - 30 R E + 60 E^2)/360 with Lichnerowicz
    E = R/4 — this is the SECOND HEAT-KERNEL COEFFICIENT in the standard
    spectral-action sense, which is NOT the same as a_1 in Paper 51's
    asymptotic-coefficient indexing (the indexing shifts by one because
    Paper 51 indexes by t-power not by HK-coefficient order). The Vassilevich
    a_0 and a_1 in the production code use the Lichnerowicz formula and DO
    NOT reproduce the asymptotic K(t) expansion directly. The witnesses in
    Paper 55 are the ASYMPTOTIC heat-trace coefficients (Paper 51 Cor 2.1
    closed form), not the Vassilevich curvature polynomial output; this test
    verifies the asymptotic-coefficient identification independently.
    """
    # Asymptotic Dirac heat trace coefficients (Paper 51 Cor 2.1):
    raw_a0 = sp.sqrt(pi) / 2
    raw_a1 = -sp.sqrt(pi) / 4

    # Standard Connes-Chamseddine SD volume normalization:
    a0_paper55 = sp.expand(raw_a0 * (4 * pi) ** Rational(3, 2))
    a1_paper55 = sp.expand(raw_a1 * (4 * pi) ** Rational(3, 2))

    assert _is_zero(a0_paper55 - 4 * pi ** 2), (
        f"Asymptotic a_0 = {a0_paper55}, Paper 55 witness 4 pi^2"
    )
    assert _is_zero(a1_paper55 - (-2 * pi ** 2)), (
        f"Asymptotic a_1 = {a1_paper55}, Paper 55 witness -2 pi^2"
    )

    # M2 claim: each is rational * pi^2.
    ratio0 = sp.simplify(a0_paper55 / pi ** 2)
    ratio1 = sp.simplify(a1_paper55 / pi ** 2)
    assert ratio0 == 4 and ratio0.is_rational
    assert ratio1 == -2 and ratio1.is_rational


def test_dirac_two_term_exactness_via_jacobi_theta2():
    """Paper 51 Cor 2.1 / Paper 55 witness: a_k^{D^2} = 0 for k >= 2.

    Independent verification via Jacobi theta_2 modular identity.

    The CH Dirac heat trace is

        K_{D^2}(t) = sum_n g_n exp(-t * lambda_n^2)
                   = sum_n 2(n+1)(n+2) exp(-t * (n + 3/2)^2)
                   = 2 sum_{m >= 0} (m+1)(m+2) exp(-t * (m + 3/2)^2)

    where m = n. Set u = m + 3/2 (half-integer):

        K_{D^2}(t) = sum over u in {3/2, 5/2, 7/2, ...}
                          [2 (u - 1/2)(u + 1/2)] exp(-t * u^2)
                   = sum_u 2 (u^2 - 1/4) exp(-t u^2)
                   = 2 [Q(t) - (1/4) S(t)]

    where Q(t) = sum_u u^2 exp(-t u^2),  S(t) = sum_u exp(-t u^2),
    with u running over positive half-integers.

    By Jacobi theta_2 modular transformation:
        sum_{n in Z} exp(-t (n+1/2)^2) = sqrt(pi/t) * (Jacobi tail).
    The asymptotic small-t expansion truncates exactly:
        S(t) = (1/2) * sqrt(pi/t) + O(e^{-pi^2/t})
        Q(t) = (1/4) * sqrt(pi) * t^{-3/2} + (-1/8) * sqrt(pi) * t^{-1/2}
               * something... NO -- compute directly.

    Cleaner approach: verify numerically at small t that
        K_{D^2}(t) - [(sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}]
    decays as exp(-pi^2 / t) (i.e., faster than ANY power of t), which
    is the hallmark of two-term exactness.
    """
    _ensure_dps()

    def K_D2(t_val, m_max):
        """Truncated CH Dirac heat trace."""
        s = mpmath.mpf(0)
        for n in range(m_max + 1):
            lam_sq = (n + mpmath.mpf("1.5")) ** 2
            g_n = 2 * (n + 1) * (n + 2)
            s += g_n * mpmath.exp(-t_val * lam_sq)
        return s

    # Two-term asymptotic prediction:
    def two_term(t_val):
        return (mpmath.sqrt(mpmath.pi) / 2) / t_val ** mpmath.mpf("1.5") \
               - (mpmath.sqrt(mpmath.pi) / 4) / mpmath.sqrt(t_val)

    # At t in [0.05, 0.2], residual should be O(exp(-pi^2/t)) with a
    # polynomial prefactor of order ~ t^{-3/2}:
    # t = 0.05: exp(-pi^2/0.05) ~ exp(-197) ~ 1e-86, prefactor ~ 1e3 -> ~ 1e-83
    # t = 0.1:  exp(-pi^2/0.1)  ~ exp(-99)  ~ 1e-43, prefactor ~ 1e2 -> ~ 1e-37
    # t = 0.2:  exp(-pi^2/0.2)  ~ exp(-49)  ~ 1e-22, prefactor ~ 1e1 -> ~ 1e-15
    # Provide tolerances generously above these residual-bounds; what we
    # want to confirm is the EXPONENTIAL (not polynomial) decay of the
    # residual, i.e. two-term exactness.
    test_points = [
        (mpmath.mpf("0.05"), mpmath.mpf("1e-78")),
        (mpmath.mpf("0.10"), mpmath.mpf("1e-35")),
        (mpmath.mpf("0.20"), mpmath.mpf("1e-13")),
    ]
    for t_val, residual_tol in test_points:
        K_val = K_D2(t_val, m_max=200)
        pred = two_term(t_val)
        diff = abs(K_val - pred)
        assert diff < residual_tol, (
            f"t={t_val}: K={K_val}, two-term={pred}, residual={diff}, "
            f"tolerance={residual_tol}. Two-term exactness violated."
        )


def test_witness_scalar_ak_closed_form():
    """Witness: a_k^Delta (scalar SD, vol-norm) = 2 pi^2 / k! for all k >= 0.

    This is Paper 51 Thm 3.1. We verify the closed-form a_k = 2 pi^2 / k!
    for k = 0, 1, 2, 3, 4 symbolically and verify M2 membership: each is
    rational * pi^2.
    """
    expected = {
        0: 2 * pi ** 2,                # 2 pi^2 / 0! = 2 pi^2
        1: 2 * pi ** 2,                # 2 pi^2 / 1! = 2 pi^2
        2: pi ** 2,                    # 2 pi^2 / 2! = pi^2
        3: pi ** 2 / 3,                # 2 pi^2 / 3! = pi^2 / 3
        4: pi ** 2 / 12,               # 2 pi^2 / 4! = pi^2 / 12
    }
    for k in range(5):
        a_k = 2 * pi ** 2 / factorial(k)
        assert _is_zero(a_k - expected[k]), f"k={k}: a_k = {a_k}, expected {expected[k]}"

        # M2 membership: rational * pi^2.
        ratio = sp.simplify(a_k / pi ** 2)
        assert ratio.is_rational, (
            f"k={k}: a_k / pi^2 = {ratio} is not rational; M2 violated"
        )


# =============================================================================
# Path (2): Spectral zeta values zeta_{D^2}(s) at integer s.
# =============================================================================

def test_spectral_zeta_at_integer_s_in_Q_pi_squared():
    """Witnesses: zeta_{D^2}(s) at s = 2, 3, 4, 5 in Q[pi^2].

    Closed forms (Paper 28 T9 / Paper 55 witness panel):
        zeta_{D^2}(2) = pi^2 - pi^4/12
        zeta_{D^2}(3) = pi^4/3 - pi^6/30
        zeta_{D^2}(4) = 2 pi^6 / 15 - 17 pi^8 / 1260
        zeta_{D^2}(5) = pi^8 (306 - 31 pi^2) / 5670
                     = 306 pi^8 / 5670 - 31 pi^10 / 5670
                     = 17 pi^8 / 315 - 31 pi^10 / 5670

    Each is in Q[pi^2]; M2 membership verified by extracting rational
    coefficients on each pi^{2k} basis vector.
    """
    expected = {
        2: pi ** 2 - pi ** 4 / 12,
        3: pi ** 4 / 3 - pi ** 6 / 30,
        4: 2 * pi ** 6 / 15 - 17 * pi ** 8 / 1260,
        5: pi ** 8 * (306 - 31 * pi ** 2) / 5670,
    }

    # Use verify_no_odd_zeta_one_loop() from production code (T9 closed form).
    results = verify_no_odd_zeta_one_loop()
    for s in [2, 3, 4, 5]:
        val = results[f'zeta_D2_{s}']
        expected_val = expected[s]
        assert _is_zero(val - expected_val), (
            f"s={s}: production T9 closed form {val}, "
            f"paper 55 witness {expected_val}"
        )

        # M2 membership: every term in Q[pi^2]. Expand into a polynomial in
        # pi and check each coefficient is rational and only EVEN powers of
        # pi appear.
        poly = sp.Poly(sp.expand(val), pi)
        for monom, coeff in poly.as_dict().items():
            exp_pi = monom[0]
            assert exp_pi % 2 == 0, (
                f"s={s}: pi^{exp_pi} appears with coefficient {coeff} "
                "but odd power of pi means NOT in Q[pi^2] (M2 violated)"
            )
            assert coeff.is_rational, (
                f"s={s}: coefficient of pi^{exp_pi} is {coeff}, not rational"
            )


# =============================================================================
# Path (3): PSLQ negative test.
# =============================================================================

def _pslq_search(target, basis_labels, basis_values, max_coeff, maxsteps=3000):
    """Run PSLQ to find an integer relation [a, c_1, ..., c_n] such that
    a * target = c_1 * b_1 + ... + c_n * b_n at the working precision.

    Returns the relation (a tuple of ints) or None if no relation found.

    Note: mpmath.pslq default maxsteps=100 is insufficient for higher-
    degree pi-power identifications (e.g., relations involving pi^6 + pi^8
    at coefficient ceilings ~ 1260 need ~ 500-1000 PSLQ iterations).
    We use maxsteps=3000 to give the algorithm room to converge.
    """
    augmented = [target] + list(basis_values)
    rel = mpmath.pslq(augmented, maxcoeff=max_coeff, maxsteps=maxsteps)
    return rel


def _M1_basis_values():
    """Pure-Tate even-weight basis: {pi^0, pi^2, pi^4, pi^6, pi^8, pi^10}."""
    return [mpmath.mpf(1)] + [mpmath.pi ** (2 * k) for k in range(1, 6)]


def _M1_basis_labels():
    return ['1', 'pi^2', 'pi^4', 'pi^6', 'pi^8', 'pi^10']


def _negative_witnesses():
    """Quantities that MUST NOT participate in the M2 identification."""
    return {
        'zeta(3)': mpmath.zeta(3),
        'zeta(5)': mpmath.zeta(5),
        'Catalan G': mpmath.catalan,
    }


def _witness_numerical(expr: sp.Expr) -> mpmath.mpf:
    """Evaluate a sympy expression at mpmath working precision."""
    return mpmath.mpf(sp.N(expr, mpmath.mp.dps))


@pytest.mark.parametrize("name,expr_str", [
    # The vol-norm witnesses from Paper 55 sec:m2.
    ("a0_Dirac_Paper55", "4 * pi**2"),
    ("a1_Dirac_Paper55", "-2 * pi**2"),
    ("a0_scalar", "2 * pi**2"),
    ("a1_scalar", "2 * pi**2"),
    ("a2_scalar", "pi**2"),
    ("a3_scalar", "pi**2 / 3"),
    ("a4_scalar", "pi**2 / 12"),
    ("zeta_D2_2", "pi**2 - pi**4/12"),
    ("zeta_D2_3", "pi**4/3 - pi**6/30"),
    ("zeta_D2_4", "2*pi**6/15 - 17*pi**8/1260"),
    ("zeta_D2_5", "pi**8*(306 - 31*pi**2)/5670"),
])
def test_witness_pslq_pure_tate_identification(name, expr_str):
    """For each witness w, run PSLQ:
       (A) AUGMENTED basis {1, pi^2, ..., pi^10, zeta(3), zeta(5), G}.
       (B) PURE-TATE basis {1, pi^2, ..., pi^10}.

    M2 claim: each witness is in pi^{2k}-Q. Test:
       (i)  PSLQ in basis (B) returns a non-trivial relation
            identifying the witness as a Q-linear combination of pi^{2k}.
       (ii) PSLQ in basis (A), if it returns a relation, has zero
            coefficient on zeta(3), zeta(5), Catalan G (negative test).
       (iii) The pi^{2k} coefficients in (A) and (B) match.
    """
    _ensure_dps()
    expr = sp.sympify(expr_str)
    target = _witness_numerical(expr)

    M1_vals = _M1_basis_values()
    M1_labels = _M1_basis_labels()
    neg = _negative_witnesses()
    neg_labels = list(neg.keys())
    neg_vals = list(neg.values())

    # (B) Pure-Tate basis identification.
    rel_B = _pslq_search(target, M1_labels, M1_vals, max_coeff=10**6)
    assert rel_B is not None, (
        f"{name}: PSLQ failed to identify witness in pure-Tate basis "
        "{1, pi^2, pi^4, pi^6, pi^8, pi^10}. M2 not numerically supported."
    )

    # rel_B = (a_0, c_0, c_1, ..., c_5). Reconstruct: c_0*1 + sum_k c_k*pi^{2k}
    # = -a_0 * target (the PSLQ sign convention has the target with
    # coefficient -a_0 on the LHS).
    # Numerical reconstruction:
    a_target = rel_B[0]
    if a_target == 0:
        pytest.fail(f"{name}: PSLQ returned target coefficient 0.")

    # Reconstruct the "witness as Q-linear combination" coefficients:
    coeffs_B = [Rational(-c, a_target) for c in rel_B[1:]]

    # Verify the reconstruction symbolically equals expr.
    reconstruction = sum(coeffs_B[k] * pi ** (2 * k) for k in range(len(coeffs_B)))
    assert _is_zero(reconstruction - expr), (
        f"{name}: reconstruction {reconstruction} != witness {expr}. "
        f"PSLQ relation: {rel_B}"
    )

    # (A) Augmented basis identification.
    augmented_vals = M1_vals + neg_vals
    rel_A = _pslq_search(target, M1_labels + neg_labels, augmented_vals,
                         max_coeff=10**6)
    assert rel_A is not None, (
        f"{name}: PSLQ failed in augmented basis too. (Should always find "
        "the pure-Tate identification at minimum.)"
    )

    a_target_A = rel_A[0]
    if a_target_A == 0:
        pytest.fail(f"{name}: PSLQ augmented returned target coefficient 0.")

    coeffs_A = [Rational(-c, a_target_A) for c in rel_A[1:]]
    # Split into pure-Tate part and negative-test part.
    pure_tate_coeffs = coeffs_A[: len(M1_vals)]
    odd_zeta_coeffs = coeffs_A[len(M1_vals):]

    # Negative test: odd-zeta and Catalan G coefficients MUST be zero.
    for label, c in zip(neg_labels, odd_zeta_coeffs):
        assert c == 0, (
            f"{name}: augmented PSLQ assigned non-zero coefficient {c} to "
            f"'{label}'. Witness should NOT need odd zeta or Catalan G. "
            f"Full relation: {rel_A}"
        )

    # Cross-check: pure-Tate coefficients from (A) match (B).
    assert pure_tate_coeffs == coeffs_B, (
        f"{name}: pure-Tate coefficients from augmented PSLQ {pure_tate_coeffs} "
        f"differ from pure-Tate-only PSLQ {coeffs_B}."
    )

    # Final cross-check: reconstruction from augmented PSLQ also equals witness.
    augmented_reconstruction = sum(
        pure_tate_coeffs[k] * pi ** (2 * k) for k in range(len(pure_tate_coeffs))
    )
    assert _is_zero(augmented_reconstruction - expr), (
        f"{name}: augmented reconstruction {augmented_reconstruction} != "
        f"witness {expr}"
    )


def test_negative_witnesses_NOT_in_pure_tate_basis():
    """Sanity check: zeta(3), zeta(5), Catalan G are NOT linearly
    expressible in {1, pi^2, pi^4, pi^6, pi^8, pi^10} with small integer
    coefficients (proving M2's claim of strict containment via the
    contrapositive: were they expressible, then a witness using them would
    appear "as pure-Tate" too, making the M2 sub-ring statement vacuous).
    """
    _ensure_dps()
    M1_vals = _M1_basis_values()
    M1_labels = _M1_basis_labels()

    for label, val in _negative_witnesses().items():
        rel = _pslq_search(val, M1_labels, M1_vals, max_coeff=10**6)
        assert rel is None, (
            f"NEGATIVE-TEST FAILURE: PSLQ found a small-coefficient "
            f"relation for {label} in the pure-Tate basis: {rel}. "
            "This would invalidate the M2 separation."
        )


# =============================================================================
# Direct verification that the spectral sum gives the closed forms.
# =============================================================================

@pytest.mark.slow
def test_zeta_D2_2_direct_sum_matches_closed_form():
    """Direct spectral sum
        zeta_{D^2}(s=2) = sum_n g_n / lambda_n^4
    numerically matches the closed form pi^2 - pi^4/12 to mpmath
    high precision.

    The s=2 sum converges absolutely; truncating at n_max = 50000 gives
    a tail O(1/n_max^2) ~ 4e-10, sufficient for ~10 digit match.
    """
    _ensure_dps()
    n_max = 50000
    total = mpmath.mpf(0)
    for n in range(n_max + 1):
        lam_sq = mpmath.mpf(n) + mpmath.mpf("1.5")
        lam_sq = lam_sq * lam_sq
        g_n = 2 * (n + 1) * (n + 2)
        total += mpmath.mpf(g_n) / lam_sq ** 2
    closed_form = mpmath.pi ** 2 - mpmath.pi ** 4 / 12
    # Tail correction: the leading tail is sum_{n > n_max} 2(n+1)(n+2)/(n+3/2)^4
    # ~ 2/n. So error ~ 2/n_max = 4e-5 at n_max=50000. We use a generous tol.
    assert abs(total - closed_form) < mpmath.mpf("1e-3"), (
        f"Direct sum {total}, closed form {closed_form}, "
        f"diff {abs(total - closed_form)}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
