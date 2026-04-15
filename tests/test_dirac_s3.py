"""Tests for geovac/dirac_s3.py (Track D1, Dirac-on-S³ Tier 1 sprint).

All tests use exact sympy arithmetic. No floating-point comparisons.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Rational

from geovac.dirac_s3 import (
    SpinorHarmonicLabel,
    count_spinor_labels,
    delta_inverse_identity,
    dirac_degeneracy,
    dirac_eigenvalue_abs,
    fock_to_ch,
    ch_to_fock,
    spinor_labels_at_n,
    verify_pi_free,
)


# ---------------------------------------------------------------------------
# Camporesi–Higuchi spectrum
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(7)))
def test_eigenvalue_abs_symbolic(n):
    """|λ_n| = n + 3/2 exactly (sympy Rational equality)."""
    got = dirac_eigenvalue_abs(n, convention="ch")
    expected = Rational(2 * n + 3, 2)
    assert got == expected
    assert isinstance(got, Rational)
    # And not a float
    assert got.q in (1, 2)  # denominator 1 or 2 (half-integer)


@pytest.mark.parametrize("n", list(range(7)))
def test_degeneracy_dirac_symbolic(n):
    """g_n^Dirac = 2 (n+1)(n+2) symbolically for n = 0..6."""
    got = dirac_degeneracy(n, sector="dirac", convention="ch")
    expected = 2 * (n + 1) * (n + 2)
    assert got == expected
    assert isinstance(got, int)


@pytest.mark.parametrize("n", list(range(7)))
def test_degeneracy_weyl_symbolic(n):
    """g_n^Weyl = (n+1)(n+2) symbolically for n = 0..6."""
    got = dirac_degeneracy(n, sector="weyl", convention="ch")
    expected = (n + 1) * (n + 2)
    assert got == expected
    assert isinstance(got, int)


# ---------------------------------------------------------------------------
# Phase 4H SM-D identity: g_3^Dirac = 40 = Δ^{-1}
# ---------------------------------------------------------------------------

def test_g3_dirac_equals_40_exact():
    """Phase 4H SM-D reproduction: g_3^Dirac = 40 as exact integer."""
    g3 = dirac_degeneracy(3, sector="dirac", convention="ch")
    assert g3 == 40
    # And as a sympy Rational equality (not approximate):
    assert sp.Integer(g3) == sp.Integer(40)


def test_delta_inverse_identity_exact():
    """delta_inverse_identity() returns (40, 1/40) exactly."""
    g3, delta = delta_inverse_identity()
    assert g3 == 40
    assert delta == Rational(1, 40)
    assert isinstance(delta, Rational)


# ---------------------------------------------------------------------------
# π-free rational certification (Paper 24 §III pattern)
# ---------------------------------------------------------------------------

def test_verify_pi_free_dirac_nmax_6():
    """verify_pi_free(6) returns True for the full Dirac sector."""
    assert verify_pi_free(6, sector="dirac") is True


def test_verify_pi_free_weyl_nmax_6():
    """verify_pi_free(6) returns True for the Weyl sector."""
    assert verify_pi_free(6, sector="weyl") is True


def test_verify_pi_free_small_n():
    """Certify n_max = 0, 1, 2 explicitly."""
    for n in (0, 1, 2):
        assert verify_pi_free(n, sector="dirac") is True
        assert verify_pi_free(n, sector="weyl") is True


# ---------------------------------------------------------------------------
# Label generator: counts match degeneracies
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(7)))
def test_label_count_matches_degeneracy_dirac(n):
    """spinor_labels_at_n yields exactly g_n^Dirac labels."""
    labels = spinor_labels_at_n(n, sector="dirac", convention="ch")
    assert len(labels) == dirac_degeneracy(n, sector="dirac", convention="ch")


@pytest.mark.parametrize("n", list(range(7)))
def test_label_count_matches_degeneracy_weyl(n):
    """spinor_labels_at_n yields exactly g_n^Weyl labels."""
    labels = spinor_labels_at_n(n, sector="weyl", convention="ch")
    assert len(labels) == dirac_degeneracy(n, sector="weyl", convention="ch")


def test_label_fields_types():
    """Labels have correct field types (integers, sympy Rational for σ)."""
    for lab in spinor_labels_at_n(3, sector="dirac", convention="ch"):
        assert isinstance(lab, SpinorHarmonicLabel)
        assert isinstance(lab.n_fock, int) and lab.n_fock >= 1
        assert isinstance(lab.l, int) and lab.l >= 0
        assert isinstance(lab.m, int)
        assert isinstance(lab.sigma, Rational)
        assert lab.sigma in (Rational(1, 2), Rational(-1, 2))
        assert lab.chirality in (+1, -1)


def test_labels_unique():
    """All labels at a given (n, sector) are distinct."""
    labels = spinor_labels_at_n(4, sector="dirac", convention="ch")
    assert len(set(labels)) == len(labels)


def test_weyl_chirality_is_plus_only():
    """Weyl sector labels have chirality = +1."""
    for lab in spinor_labels_at_n(3, sector="weyl", convention="ch"):
        assert lab.chirality == +1
        assert lab.sigma == Rational(1, 2)


def test_dirac_has_both_chiralities():
    """Full Dirac sector populates both chiralities equally."""
    labels = spinor_labels_at_n(3, sector="dirac", convention="ch")
    plus = [l for l in labels if l.chirality == +1]
    minus = [l for l in labels if l.chirality == -1]
    assert len(plus) == len(minus) == 20  # (3+1)(3+2) = 20


# ---------------------------------------------------------------------------
# Convention translation
# ---------------------------------------------------------------------------

def test_fock_ch_roundtrip():
    for n in range(1, 8):
        assert ch_to_fock(fock_to_ch(n)) == n
    for n in range(7):
        assert fock_to_ch(ch_to_fock(n)) == n


def test_fock_convention_gives_same_spectrum():
    """Using Fock convention (n_Fock = n_CH + 1) gives the same spectrum."""
    for n_ch in range(7):
        n_fock = n_ch + 1
        lam_ch = dirac_eigenvalue_abs(n_ch, convention="ch")
        lam_fock = dirac_eigenvalue_abs(n_fock, convention="fock")
        assert lam_ch == lam_fock
        g_ch = dirac_degeneracy(n_ch, sector="dirac", convention="ch")
        g_fock = dirac_degeneracy(n_fock, sector="dirac", convention="fock")
        assert g_ch == g_fock


def test_count_up_to_n_max():
    """Cumulative count up to n_max = 3 for full Dirac is 2·Σ(n+1)(n+2)."""
    # n=0: 2, n=1: 6, n=2: 12, n=3: 20 → Weyl cumulative = 40, Dirac = 80
    assert count_spinor_labels(3, sector="weyl", convention="ch") == (
        2 + 6 + 12 + 20)
    assert count_spinor_labels(3, sector="dirac", convention="ch") == 2 * (
        2 + 6 + 12 + 20)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def test_invalid_sector():
    with pytest.raises(ValueError):
        dirac_degeneracy(3, sector="majorana", convention="ch")  # type: ignore[arg-type]


def test_invalid_convention():
    with pytest.raises(ValueError):
        dirac_eigenvalue_abs(3, convention="weird")  # type: ignore[arg-type]


def test_negative_n_rejected():
    with pytest.raises(ValueError):
        dirac_eigenvalue_abs(-1, convention="ch")
    with pytest.raises(ValueError):
        dirac_degeneracy(-1, sector="dirac", convention="ch")
    with pytest.raises(ValueError):
        verify_pi_free(-1)


def test_fock_must_be_positive():
    with pytest.raises(ValueError):
        fock_to_ch(0)
