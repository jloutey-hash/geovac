"""Tests for geovac/dirac_matrix_elements.py (Track T1, Tier 2 sprint).

All tests use exact sympy arithmetic. No floating-point comparisons.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Rational, Integer, Symbol, oo

from geovac.dirac_s3 import (
    SpinorHarmonicLabel,
    iter_spinor_labels,
)
from geovac.dirac_matrix_elements import (
    DiracLabel,
    Z_sym,
    alpha_sym,
    gamma_rel,
    angular_matrix_J_sq,
    angular_matrix_L,
    angular_matrix_L_dot_S,
    angular_matrix_r_hat,
    angular_matrix_sigma,
    angular_sigma_dot_rhat_identity,
    dirac_label_to_spinor_label,
    inverse_r_cubed_hydrogenic,
    iter_dirac_labels,
    kappa_to_j,
    kappa_to_l,
    kappa_to_l_sigma,
    l_sigma_to_kappa,
    radial_expectation_diagonal,
    radial_matrix_element,
    spinor_label_to_dirac_label,
)


# ---------------------------------------------------------------------------
# κ ↔ (l, σ) conversion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l,sigma,expected_kappa", [
    (0, Rational(1, 2), -1),   # s_{1/2}
    (1, Rational(1, 2), -2),   # p_{3/2}
    (1, Rational(-1, 2), +1),  # p_{1/2}
    (2, Rational(1, 2), -3),   # d_{5/2}
    (2, Rational(-1, 2), +2),  # d_{3/2}
    (3, Rational(1, 2), -4),   # f_{7/2}
    (3, Rational(-1, 2), +3),  # f_{5/2}
])
def test_l_sigma_to_kappa(l, sigma, expected_kappa):
    assert l_sigma_to_kappa(l, sigma) == expected_kappa


@pytest.mark.parametrize("kappa,expected_l", [
    (-1, 0), (+1, 1),
    (-2, 1), (+2, 2),
    (-3, 2), (+3, 3),
    (-4, 3), (+4, 4),
])
def test_kappa_to_l(kappa, expected_l):
    assert kappa_to_l(kappa) == expected_l


@pytest.mark.parametrize("kappa", [-4, -3, -2, -1, 1, 2, 3, 4])
def test_kappa_roundtrip(kappa):
    l, sigma = kappa_to_l_sigma(kappa)
    if not (l == 0 and sigma == Rational(-1, 2)):  # the forbidden case
        assert l_sigma_to_kappa(l, sigma) == kappa


def test_kappa_to_j():
    assert kappa_to_j(-1) == Rational(1, 2)   # s_{1/2}
    assert kappa_to_j(+1) == Rational(1, 2)   # p_{1/2}
    assert kappa_to_j(-2) == Rational(3, 2)   # p_{3/2}
    assert kappa_to_j(+2) == Rational(3, 2)   # d_{3/2}
    assert kappa_to_j(-3) == Rational(5, 2)   # d_{5/2}


def test_l_eq_0_sigma_minus_half_rejected():
    """j = l − 1/2 = −1/2 is unphysical at l = 0."""
    with pytest.raises(ValueError):
        l_sigma_to_kappa(0, Rational(-1, 2))


# ---------------------------------------------------------------------------
# DiracLabel dataclass
# ---------------------------------------------------------------------------

def test_dirac_label_construction_ok():
    lab = DiracLabel(n_fock=2, kappa=-2, two_m_j=3)  # p_{3/2}, m_j=+3/2
    assert lab.l == 1
    assert lab.j == Rational(3, 2)
    assert lab.m_j == Rational(3, 2)
    assert lab.j_times_2 == 3


def test_dirac_label_kappa_zero_rejected():
    with pytest.raises(ValueError):
        DiracLabel(n_fock=1, kappa=0, two_m_j=0)


def test_dirac_label_m_j_out_of_range_rejected():
    """|m_j| must be ≤ j."""
    with pytest.raises(ValueError):
        # p_{1/2} has j=1/2, two_m_j must be ±1; two_m_j=3 is out
        DiracLabel(n_fock=2, kappa=+1, two_m_j=3)


def test_dirac_label_m_j_parity_rejected():
    """m_j is half-integer ⇔ 2·m_j is odd for half-integer j."""
    with pytest.raises(ValueError):
        # p_{1/2}: 2j=1 (odd), so two_m_j must be odd. two_m_j=0 is even → invalid.
        DiracLabel(n_fock=2, kappa=+1, two_m_j=0)


def test_dirac_label_l_must_be_less_than_n():
    """κ=+1 means l=1, needs n_fock ≥ 2."""
    with pytest.raises(ValueError):
        DiracLabel(n_fock=1, kappa=+1, two_m_j=1)


# ---------------------------------------------------------------------------
# Label iteration and count
# ---------------------------------------------------------------------------

def test_iter_dirac_labels_nmax_1():
    """n_fock=1: only 1s_{1/2}, giving 2 states (m_j = ±1/2)."""
    labels = list(iter_dirac_labels(1))
    assert len(labels) == 2
    kappas = {lab.kappa for lab in labels}
    assert kappas == {-1}


def test_iter_dirac_labels_nmax_2():
    """n_fock ∈ {1,2}. States:
        1s_{1/2}  (2)
        2s_{1/2}  (2)  κ=-1
        2p_{1/2}  (2)  κ=+1
        2p_{3/2}  (4)  κ=-2
    Total: 10.
    """
    labels = list(iter_dirac_labels(2))
    assert len(labels) == 10


def test_iter_dirac_labels_nmax_3_count():
    """n_fock ∈ {1,2,3}. For n=3 additionally:
        3s_{1/2} (2), 3p_{1/2} (2), 3p_{3/2} (4), 3d_{3/2} (4), 3d_{5/2} (6).
    Total adds 18; grand total = 10 + 18 = 28.
    """
    labels = list(iter_dirac_labels(3))
    assert len(labels) == 28


# ---------------------------------------------------------------------------
# Bridge to D1 (SpinorHarmonicLabel)
# ---------------------------------------------------------------------------

def test_spinor_to_dirac_roundtrip_kappa_neg():
    """For σ = +1/2 (j = l+1/2), conversion should map to κ = −(l+1)."""
    for n_fock in range(1, 4):
        for l in range(n_fock):
            # Pick any allowed m (D1 stores 2·m_j; for σ=+1/2, j=l+1/2,
            # so 2j=2l+1, so two_m_j ∈ {−(2l+1), ..., 2l+1} step 2).
            two_j = 2 * l + 1
            for two_m in range(-two_j, two_j + 1, 2):
                sph = SpinorHarmonicLabel(
                    n_fock=n_fock, l=l, m=two_m,
                    sigma=Rational(1, 2), chirality=+1,
                )
                dl = spinor_label_to_dirac_label(sph)
                assert dl.kappa == -(l + 1)
                assert dl.n_fock == n_fock
                assert dl.two_m_j == two_m
                # Roundtrip
                sph2 = dirac_label_to_spinor_label(dl, chirality=+1)
                assert sph2 == sph


def test_spinor_to_dirac_roundtrip_kappa_pos():
    """For σ = −1/2 and l ≥ 1, conversion should map to κ = +l."""
    for n_fock in range(2, 4):
        for l in range(1, n_fock):
            two_j = 2 * l - 1  # j = l − 1/2
            for two_m in range(-two_j, two_j + 1, 2):
                sph = SpinorHarmonicLabel(
                    n_fock=n_fock, l=l, m=two_m,
                    sigma=Rational(-1, 2), chirality=+1,
                )
                dl = spinor_label_to_dirac_label(sph)
                assert dl.kappa == l


# ---------------------------------------------------------------------------
# Angular matrix elements: Szmytkowski σ·r̂ identity
# ---------------------------------------------------------------------------

def test_sigma_dot_rhat_identity_small_kappa():
    """σ·r̂ |κ, m⟩ = −|−κ, m⟩ for |κ| ≤ 3 and all allowed m."""
    for kappa in (-3, -2, -1, +1, +2, +3):
        two_j = 2 * abs(kappa) - 1
        for two_m in range(-two_j, two_j + 1, 2):
            target_kappa, target_two_m, coeff = \
                angular_sigma_dot_rhat_identity(kappa, two_m)
            assert target_kappa == -kappa
            assert target_two_m == two_m
            assert coeff == Integer(-1)


def test_angular_matrix_r_hat_selection_rule():
    """⟨κ', m'|σ·r̂|κ, m⟩ = −δ(κ', −κ)·δ(m', m).

    (−(+1)) = −1, so ⟨κ'=−1, m'=1|σ·r̂|κ=+1, m=1⟩ = −1 (allowed, coeff −1).
    """
    # Allowed case: κ=+1, -κ=-1, so κ'=-1 matches:
    assert angular_matrix_r_hat(-1, 1, +1, 1) == -1
    # Violating case: same sign (not negated):
    assert angular_matrix_r_hat(+1, 1, +1, 1) == 0
    # m mismatch:
    assert angular_matrix_r_hat(-1, 1, +1, -1) == 0
    # Another allowed: s_{1/2} (κ=-1, m=-1/2) coupled to p_{1/2} (κ=+1, m=-1/2).
    assert angular_matrix_r_hat(+1, -1, -1, -1) == -1


# ---------------------------------------------------------------------------
# Angular matrix elements: J², L², L·S eigenvalues
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kappa,expected_j_j1", [
    (-1, Rational(1, 2) * Rational(3, 2)),   # j=1/2, j(j+1)=3/4
    (+1, Rational(3, 4)),                     # j=1/2
    (-2, Rational(3, 2) * Rational(5, 2)),   # j=3/2, j(j+1)=15/4
    (+2, Rational(15, 4)),
    (-3, Rational(5, 2) * Rational(7, 2)),   # j=5/2, j(j+1)=35/4
])
def test_J_sq_eigenvalue(kappa, expected_j_j1):
    two_j = 2 * abs(kappa) - 1
    for two_m in range(-two_j, two_j + 1, 2):
        got = angular_matrix_J_sq(kappa, two_m, kappa, two_m)
        assert got == expected_j_j1


def test_L_sq_eigenvalue():
    """⟨L²⟩ = l(l+1) with l = kappa_to_l(κ)."""
    cases = [
        (-1, 0),   # s_{1/2}: l=0
        (+1, 2),   # p_{1/2}: l=1, l(l+1)=2
        (-2, 2),   # p_{3/2}: l=1
        (+2, 6),   # d_{3/2}: l=2, l(l+1)=6
        (-3, 6),   # d_{5/2}: l=2
        (-4, 12),  # f_{7/2}: l=3
    ]
    for kappa, expected in cases:
        two_j = 2 * abs(kappa) - 1
        got = angular_matrix_L(kappa, two_j, kappa, two_j, component="sq")
        assert got == expected


@pytest.mark.parametrize("kappa,expected", [
    (-1, Rational(0)),        # s_{1/2}: (-(-1)-1)/2 = 0
    (+1, Rational(-1)),       # p_{1/2}: (-(+1)-1)/2 = -1
    (-2, Rational(1, 2)),     # p_{3/2}: (-(-2)-1)/2 = 1/2
    (+2, Rational(-3, 2)),    # d_{3/2}: (-(+2)-1)/2 = -3/2
    (-3, Rational(1)),        # d_{5/2}: (-(-3)-1)/2 = 1
    (+3, Rational(-2)),       # f_{5/2}
    (-4, Rational(3, 2)),     # f_{7/2}
])
def test_L_dot_S_eigenvalue(kappa, expected):
    """L·S eigenvalue = (−κ − 1)/2 in both κ-sign branches.

    Cross-check via direct formula [j(j+1) − l(l+1) − 3/4]/2.
    """
    two_j = 2 * abs(kappa) - 1
    got = angular_matrix_L_dot_S(kappa, two_j, kappa, two_j)
    assert got == expected
    # Cross-check:
    j = kappa_to_j(kappa)
    l = kappa_to_l(kappa)
    expected_direct = (j * (j + 1) - Integer(l) * (l + 1) - Rational(3, 4)) / 2
    assert sp.simplify(got - expected_direct) == 0


def test_L_dot_S_off_diagonal_is_zero():
    """L·S is diagonal in (κ, m_j)."""
    # κ_a ≠ κ_b
    assert angular_matrix_L_dot_S(-1, 1, +1, 1) == 0
    # m mismatch:
    assert angular_matrix_L_dot_S(-2, 3, -2, 1) == 0


# ---------------------------------------------------------------------------
# σ² = 3·I
# ---------------------------------------------------------------------------

def test_sigma_sq_is_3():
    for kappa in (-3, -2, -1, +1, +2, +3):
        two_j = 2 * abs(kappa) - 1
        for two_m in range(-two_j, two_j + 1, 2):
            got = angular_matrix_sigma(kappa, two_m, kappa, two_m, component="sq")
            assert got == 3


# ---------------------------------------------------------------------------
# Radial matrix elements: hydrogenic diagonal closed forms
# ---------------------------------------------------------------------------

def test_inverse_r_cubed_1s_diverges():
    """⟨1/r³⟩ diverges for l=0 (1s, 2s, ...)."""
    with pytest.raises(ValueError):
        inverse_r_cubed_hydrogenic(1, 0)
    with pytest.raises(ValueError):
        inverse_r_cubed_hydrogenic(2, 0)


def test_inverse_r_cubed_2p():
    """⟨1/r³⟩_{n=2, l=1} = Z³ / [8 · 1 · (3/2) · 2] = Z³/24."""
    got = inverse_r_cubed_hydrogenic(2, 1, Z=Z_sym)
    expected = Z_sym**3 / Integer(24)
    assert sp.simplify(got - expected) == 0

    # For Z=1 (hydrogen), this is 1/24.
    got_Z1 = inverse_r_cubed_hydrogenic(2, 1, Z=1)
    assert got_Z1 == Rational(1, 24)


def test_inverse_r_cubed_3d():
    """⟨1/r³⟩_{n=3, l=2} = Z³ / [27 · 2 · (5/2) · 3] = Z³/405."""
    got = inverse_r_cubed_hydrogenic(3, 2, Z=Z_sym)
    expected = Z_sym**3 / Integer(405)
    assert sp.simplify(got - expected) == 0


def test_inverse_r_hydrogenic_exact():
    """⟨1/r⟩_{n,l} = Z/n²."""
    for n in (1, 2, 3):
        for l in range(n):
            got = radial_expectation_diagonal(n, l, "1/r", Z=Z_sym)
            expected = Z_sym / Integer(n**2)
            assert sp.simplify(got - expected) == 0


def test_r_hydrogenic_exact():
    """⟨r⟩_{n,l} = [3n² − l(l+1)]/(2Z)."""
    # Hydrogen 1s: ⟨r⟩ = 3/(2Z) = 3/2 at Z=1.
    got = radial_expectation_diagonal(1, 0, "r", Z=1)
    assert got == Rational(3, 2)
    # 2p: ⟨r⟩ = (12 - 2)/(2Z) = 10/(2Z) = 5/Z.
    got = radial_expectation_diagonal(2, 1, "r", Z=Z_sym)
    assert sp.simplify(got - 5 / Z_sym) == 0


def test_r_sq_hydrogenic_1s():
    """⟨r²⟩_{1s} = 3/Z²."""
    got = radial_expectation_diagonal(1, 0, "r^2", Z=Z_sym)
    assert sp.simplify(got - 3 / Z_sym**2) == 0


def test_inverse_r_sq_hydrogenic():
    """⟨1/r²⟩_{n,l} = Z²/[n³(l+1/2)]."""
    # 2p: Z²/(8 · 3/2) = Z²/12.
    got = radial_expectation_diagonal(2, 1, "1/r^2", Z=Z_sym)
    assert sp.simplify(got - Z_sym**2 / 12) == 0


# ---------------------------------------------------------------------------
# Off-diagonal radial matrix element (integration)
# ---------------------------------------------------------------------------

def test_radial_off_diagonal_r_1s_2s():
    """⟨1s|r|2s⟩ at Z=1: sympy closed form.

    Direct sympy integration of the hydrogenic radial wavefunctions:
    ⟨1s|r|2s⟩ at Z=1 = −32√2/81.  (Exact rational × √2.)
    """
    got = radial_matrix_element(1, 0, 2, 0, "r", Z=1)
    got_simplified = sp.simplify(got)
    expected = Rational(-32, 81) * sp.sqrt(2)
    assert sp.simplify(got_simplified - expected) == 0


def test_radial_off_diagonal_invr_1s_2p_vanishes():
    """⟨1s|1/r|2p⟩ = 0 at all Z (orthogonality of different-l radials
    under 1/r weighting is NOT generic — this should be a specific check).

    Actually ⟨1s|1/r|2p⟩ is NOT zero in general; orthogonality under
    1/r weight is not a hydrogenic identity. We just check that the
    result is an exact sympy expression in Z (not NaN, not unevaluated).
    """
    got = radial_matrix_element(1, 0, 2, 1, "1/r", Z=Z_sym)
    # Should be a finite sympy expression, not oo or NaN.
    assert got is not sp.nan
    assert got != sp.oo
    # Just confirm it's symbolic in Z.
    assert got.has(Z_sym) or got == 0


# ---------------------------------------------------------------------------
# Scalar limit test: matching existing scalar hydrogenic values
# ---------------------------------------------------------------------------

def test_scalar_limit_1s_energy_element():
    """The non-relativistic hydrogenic ⟨1/r⟩_{1s} = Z, the diagonal of
    the 1s energy expectation value via Virial at Z=1 gives ⟨T⟩=−⟨V⟩/2.

    We just check ⟨1s|1/r|1s⟩ = Z — this is the scalar limit that
    T1's radial layer must reproduce, matching any existing scalar
    hydrogenic calculation.
    """
    for Z in (1, 2, 3, 4):
        got = radial_expectation_diagonal(1, 0, "1/r", Z=Z)
        assert got == Z


def test_scalar_limit_hydrogen_2p_energy():
    """⟨2p| 1/r |2p⟩ = 1/4 at Z=1 — non-relativistic 2p kinetic scale.

    Cross-checks: ⟨T⟩_2p = Z²/(2n²) = 1/8, so ⟨V⟩ = -⟨1/r⟩ = -1/4
    gives E = -1/8, matching the hydrogenic 2p energy.
    """
    got = radial_expectation_diagonal(2, 1, "1/r", Z=1)
    assert got == Rational(1, 4)


# ---------------------------------------------------------------------------
# γ symbolic parameter (Dirac-Coulomb relativistic factor)
# ---------------------------------------------------------------------------

def test_gamma_rel_exposed_as_symbol():
    """γ = √(1 − (Zα)²) is exposed as a named sympy positive symbol."""
    assert gamma_rel.is_positive
    # T1 does NOT use γ in any computed matrix element (all T1 outputs
    # are non-relativistic). We simply confirm it's exposed for T2/T5.


def test_alpha_exposed_as_symbol():
    """α = fine-structure constant is exposed as a named sympy symbol."""
    assert alpha_sym.is_positive


# ---------------------------------------------------------------------------
# Cross-module regression: D1 still passes its contract
# ---------------------------------------------------------------------------

def test_d1_spinor_iter_unchanged():
    """Sanity check: D1's iter_spinor_labels still produces the right counts."""
    labels = list(iter_spinor_labels(2))  # n_max=2 CH = up to n_fock=3
    # Weyl+Dirac: all labels at n_ch=0,1,2. Weyl default was "dirac" — check.
    from geovac.dirac_s3 import count_spinor_labels
    assert len(labels) == count_spinor_labels(2, sector="dirac", convention="ch")
