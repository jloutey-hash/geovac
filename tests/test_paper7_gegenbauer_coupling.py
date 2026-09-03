"""Paper 7 eq:fock_coupling (corrected 2026-09-03) -- the inter-shell coupling of
cos(chi) in the Gegenbauer eigenbasis of the unit-S^3 Laplacian.

    c^2(n, l) := |<n+1, l | cos chi | n, l>|^2 = (1/4) [1 - l(l+1) / (n(n+1))]

with radial functions  R_{nl}(chi) ~ sin^l(chi) C^{(l+1)}_{n-l-1}(cos chi)
on the measure sin^2(chi) dchi.  Closed form: Gegenbauer three-term
recurrence + norm ratio.  [SYMBOLIC + MEASURED]

WHY THIS FILE EXISTS (trunk FULL run #3, 2026-09-03, carryforward I.1.1).
Paper 7 (and Papers 0, 2, 18, 32 and the group3 synthesis) printed the
formula with prefactor 1/16 and read it as "(1/4)^2"; the Chebyshev-U
recurrence  cos(chi) U_{n-1} = U_n/2 + U_{n-2}/2  has amplitude 1/2, so
c^2(n, 0) = 1/4.  The 1/16 that the matching kappa = -1/16 coincides with
is a DIFFERENT quantity, the inverse Fock Jacobian 1/Omega^4(0) = c^2(n,0)/4.
The old backing test reached 1/16 through an unexplained `/2`; this file
replaces it with an independent quadrature route.  The paper's special
value 1/40 = Delta is the COMPOSITE  [1 - l(l+1)/(n(n+1))]_{(4,3)} / Omega^4(0)
= (2/5)(1/16); the actual coupling is c^2(4,3) = 1/10.

Tiers: closed form is a symbolic statement (INTERNAL THEOREM via the
recurrence); the quadrature panel is its numerical certificate.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_gegenbauer

PAIRS = [(1, 0), (2, 0), (3, 0), (4, 3), (5, 2), (7, 4), (6, 1)]


def c2_derived(n: int, l: int) -> Fraction:
    """The derived closed form (exact rational)."""
    return Fraction(1, 4) * (1 - Fraction(l * (l + 1), n * (n + 1)))


def c2_retired(n: int, l: int) -> Fraction:
    """The formula printed before 2026-09-03 (prefactor 1/16) -- kept only so
    the guard below can prove it is NOT the matrix element."""
    return Fraction(1, 16) * (1 - Fraction(l * (l + 1), n * (n + 1)))


def _radial(n: int, l: int):
    f = lambda c: eval_gegenbauer(n - l - 1, l + 1, np.cos(c)) * np.sin(c) ** l
    norm = np.sqrt(quad(lambda c: f(c) ** 2 * np.sin(c) ** 2, 0.0, np.pi)[0])
    return lambda c: f(c) / norm


def coupling_by_quadrature(n: int, l: int) -> float:
    a, b = _radial(n, l), _radial(n + 1, l)
    me = quad(lambda c: a(c) * b(c) * np.cos(c) * np.sin(c) ** 2, 0.0, np.pi)[0]
    return me ** 2


@pytest.mark.parametrize("n,l", PAIRS)
def test_closed_form_matches_quadrature(n, l):
    """Independent route: numerical Gegenbauer quadrature vs the closed form."""
    assert abs(coupling_by_quadrature(n, l) - float(c2_derived(n, l))) < 1e-9


@pytest.mark.parametrize("n,l", PAIRS)
def test_retired_one_sixteenth_formula_is_not_the_matrix_element(n, l):
    """Guard (fires on the retired formula): the 1/16-prefactor version is
    exactly four times too small at every (n, l)."""
    q = coupling_by_quadrature(n, l)
    assert abs(q - float(c2_retired(n, l))) > 1e-3
    assert abs(q / float(c2_retired(n, l)) - 4.0) < 1e-8


def test_s_wave_amplitude_is_one_half():
    """l = 0: R_{n0} ~ sin(n chi)/sin(chi); cos(chi) sin(n chi) =
    [sin((n+1)chi) + sin((n-1)chi)]/2 with equal norms, so the amplitude is
    exactly 1/2 and c^2(n, 0) = 1/4 for every n."""
    for n in range(1, 8):
        assert c2_derived(n, 0) == Fraction(1, 4)
        assert abs(coupling_by_quadrature(n, 0) - 0.25) < 1e-10


def test_inverse_fock_jacobian_is_a_different_quantity():
    """1/Omega^4(0) = 1/16 with Omega(0) = 2 is c^2(n,0)/4, not c^2(n,0):
    the kappa = -1/16 coincidence is with the Jacobian, not the coupling."""
    omega0 = Fraction(2, 1)
    inv_jacobian = 1 / omega0 ** 4
    assert inv_jacobian == Fraction(1, 16)
    assert c2_derived(3, 0) / 4 == inv_jacobian
    assert c2_derived(3, 0) != inv_jacobian


def test_delta_is_the_composite_not_the_coupling():
    """Delta = 1/40 = [1 - l(l+1)/(n(n+1))]_{(4,3)} / Omega^4(0) = (2/5)(1/16);
    the actual coupling at (4,3) is 1/10."""
    casimir_factor = 1 - Fraction(12, 20)
    assert casimir_factor == Fraction(2, 5)
    assert casimir_factor * Fraction(1, 16) == Fraction(1, 40)
    assert c2_derived(4, 3) == Fraction(1, 10)
    assert abs(coupling_by_quadrature(4, 3) - 0.1) < 1e-9
