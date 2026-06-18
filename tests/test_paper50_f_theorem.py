"""Genuine F-theorem test for Paper 50 (CFT3 partition function), Thm 3.4
(``thm:scalar_bitexact`` + ``thm:dirac_bitexact``): the framework's S^3
spectral-zeta zeta'(0) reproduces the Klebanov-Pufu-Safdi (2011) continuum
free energies bit-exactly.

WHY THIS EXISTS.  Before v4.21.x the only "F-theorem" coverage was the LITERAL
cross-product ``F_D - 2 F_s`` in ``tests/test_paper55_m1_pure_tate.py`` (a Paper
55 test), which types the KPS F-values in as constants and never computes them
from the framework spectrum -- a false-positive backing flagged by the group1
Bite-A code review.  This test COMPUTES F from the GeoVac S^3 spectrum
(degeneracy + eigenvalue -> analytic continuation -> zeta'(0)), so a wrong
framework spectrum would now be caught.

FRAMEWORK SPECTRUM (the input under test):
  conformally coupled scalar on unit S^3:
      lambda_l = (l+1/2)(l+3/2) = (l+1)^2 - 1/4,  deg = (l+1)^2,  l >= 0
      [ -Delta eigenvalue l(l+2) plus conformal mass R/8 = 3/4 ]
  Dirac on unit S^3 (Camporesi-Higuchi):
      |lambda_n| = n + 3/2,  deg = 2(n+1)(n+2),  n >= 0

KPS continuum (Klebanov-Pufu-Safdi, arXiv:1105.4598):
  F_s = log2/8 - 3 zeta(3)/(16 pi^2)  =  -1/2 zeta'_scalar(0)   ~ 0.0638071
  F_D = log2/4 + 3 zeta(3)/( 8 pi^2)  =       zeta'_Dirac(0)    ~ 0.2189595
"""

import mpmath as mp
import pytest

from geovac.qed_two_loop import (
    conformal_scalar_zeta_s3,
    dirac_dirichlet_series_hurwitz,
    dirac_dirichlet_series_numerical,
    dirac_F_theorem,
    scalar_F_theorem,
)

mp.mp.dps = 60


# --- KPS reference closed forms ----------------------------------------------

def _F_s_kps() -> mp.mpf:
    return mp.log(2) / 8 - 3 * mp.zeta(3) / (16 * mp.pi ** 2)


def _F_D_kps() -> mp.mpf:
    return mp.log(2) / 4 + 3 * mp.zeta(3) / (8 * mp.pi ** 2)


# ---------------------------------------------------------------------------
# Spectrum ties: the analytic continuations reproduce the direct spectral sums
# where those converge (Re s large).  This pins the continuation to the
# framework's S^3 spectrum, not an abstract zeta.
# ---------------------------------------------------------------------------

def test_scalar_continuation_matches_direct_spectral_sum():
    """conformal_scalar_zeta_s3(s) == sum_l (l+1)^2 [(l+1)^2-1/4]^{-s} at s=4."""
    s = mp.mpf(4)
    direct = mp.nsum(
        lambda l: ((l + 1) ** 2) * (((l + 1) ** 2 - mp.mpf(1) / 4) ** (-s)),
        [0, mp.inf],
    )
    cont = conformal_scalar_zeta_s3(s)
    assert abs(direct - cont) < mp.mpf(10) ** (-25)


def test_dirac_continuation_matches_direct_spectral_sum():
    """The Hurwitz Dirac zeta (2 zeta(s-2,3/2) - 1/2 zeta(s,3/2)) reproduces the
    extrapolated direct sum 2(n+1)(n+2)(n+3/2)^{-s} at s=4 -- the spectrum tie
    that the s=0 derivative below inherits."""
    s = mp.mpf(4)
    direct = mp.nsum(
        lambda n: 2 * (n + 1) * (n + 2) * ((n + mp.mpf(3) / 2) ** (-s)),
        [0, mp.inf],
    )
    hurwitz = (2 * mp.zeta(s - 2, mp.mpf(3) / 2)
               - mp.mpf(1) / 2 * mp.zeta(s, mp.mpf(3) / 2))
    assert abs(direct - hurwitz) < mp.mpf(10) ** (-25)
    # cross-check the shipped fixed-N summation converges toward the same value
    res = dirac_dirichlet_series_numerical(4, N=2000)
    assert res["rel_error_vs_hurwitz"] < 1e-3


# ---------------------------------------------------------------------------
# F-theorem: framework zeta'(0) == KPS continuum, bit-exact
# ---------------------------------------------------------------------------

def test_scalar_F_theorem_bit_exact_kps():
    """F_s = -1/2 zeta'_scalar(0) computed from the framework spectrum equals
    the KPS value log2/8 - 3 zeta(3)/(16 pi^2) to the precision floor."""
    F_s = scalar_F_theorem()
    assert abs(F_s - _F_s_kps()) < mp.mpf(10) ** (-40)


def test_dirac_F_theorem_bit_exact_kps():
    """F_D = zeta'_Dirac(0) computed from the framework spectrum equals the KPS
    value log2/4 + 3 zeta(3)/(8 pi^2) to the precision floor."""
    F_D = dirac_F_theorem()
    assert abs(F_D - _F_D_kps()) < mp.mpf(10) ** (-40)


def test_scalar_F_theorem_independent_second_method():
    """Independent cross-check of zeta'_scalar(0): numerical differentiation of
    the binomial continuation at s=0 must agree with the term-by-term series in
    scalar_F_theorem (guards against an error in either derivation)."""
    zp0_diff = mp.diff(lambda s: conformal_scalar_zeta_s3(s), mp.mpf(0))
    F_s_diff = -zp0_diff / 2
    assert abs(F_s_diff - scalar_F_theorem()) < mp.mpf(10) ** (-25)


def test_dirac_F_theorem_independent_second_method():
    """Independent cross-check of zeta'_Dirac(0): numerical differentiation of
    the closed-form coefficients (2(2^{s-2}-1) zeta(s-2) - (2^s-1)/2 zeta(s),
    the coeff_a/coeff_b of dirac_dirichlet_series_hurwitz) at s=0 must agree
    with the Hurwitz-derivative form in dirac_F_theorem."""
    def zeta_dirac(s):
        return (2 * (mp.mpf(2) ** (s - 2) - 1) * mp.zeta(s - 2)
                - (mp.mpf(2) ** s - 1) / 2 * mp.zeta(s))

    zp0 = mp.diff(zeta_dirac, mp.mpf(0))
    assert abs(zp0 - dirac_F_theorem()) < mp.mpf(10) ** (-25)


def test_closed_form_coeffs_match_shipped_function():
    """The closed-form coefficients used for the s=0 continuation are the same
    coeff_a/coeff_b the shipped dirac_dirichlet_series_hurwitz reports at s=4,
    so the F_D derivation rides on the already-tested T9 representation."""
    row = dirac_dirichlet_series_hurwitz(4)
    assert row["coeff_zeta_R_s_minus_2"] == 2 * (2 ** (4 - 2) - 1)
    assert row["coeff_zeta_R_s"] == -(2 ** 4 - 1) / 2


def test_F_theorem_cross_product_consistency():
    """Sanity: the framework F-values reproduce the Paper 55 W5 cross-product
    F_D - 2 F_s = 3 zeta(3)/(4 pi^2) (the 1/pi^2 M1 witness) -- now from
    COMPUTED F-values rather than literals."""
    cross = dirac_F_theorem() - 2 * scalar_F_theorem()
    assert abs(cross - 3 * mp.zeta(3) / (4 * mp.pi ** 2)) < mp.mpf(10) ** (-40)
