"""
TRUNK QA — Claim 4: Delta = 1/40 = (2/5) / Omega^4(0)  (Paper 7 -> Paper 2).

CORRECTED 2026-09-03 (trunk FULL run #3).  Paper 7 had printed
c^2(n,l) = (1/16)[1 - l(l+1)/(n(n+1))] and c^2(4,3) = 1/40 = Delta.  The
derived inter-shell coupling is (1/4)[1 - l(l+1)/(n(n+1))] (Chebyshev
amplitude 1/2; tests/test_paper7_gegenbauer_coupling.py), so c^2(4,3) = 1/10.
The number that equals Delta is the COMPOSITE
    [1 - l(l+1)/(n(n+1))]_{(4,3)} / Omega^4(0) = (2/5) (1/16) = 1/40,
the Casimir factor times the inverse Fock Jacobian -- not a matrix element
of anything.  Delta itself is independently derived in production code from
g_3^Dirac = 2(n+1)(n+2)|_{n=3} = 40.  The coincidence is an OBSERVATION about
a composite; this file records exactly that status.
"""

from __future__ import annotations

import sympy as sp
from sympy import Rational

from geovac.dirac_s3 import delta_inverse_identity, dirac_degeneracy


def casimir_factor(n: int, l: int) -> sp.Rational:
    """1 - l(l+1)/(n(n+1)); at (4,3) this is 2/5."""
    return 1 - Rational(l * (l + 1), n * (n + 1))


def c2_formula(n: int, l: int) -> sp.Rational:
    """The DERIVED coupling |<n+1,l|cos chi|n,l>|^2 = (1/4) casimir_factor
    (tests/test_paper7_gegenbauer_coupling.py verifies it by quadrature)."""
    return Rational(1, 4) * casimir_factor(n, l)


INVERSE_FOCK_JACOBIAN = Rational(1, 16)      # 1/Omega^4(0), Omega(0) = 2


def delta_composite(n: int, l: int) -> sp.Rational:
    """The composite that equals Delta at (4,3): casimir_factor / Omega^4(0)."""
    return casimir_factor(n, l) * INVERSE_FOCK_JACOBIAN


def test_composite_4_3_is_one_fortieth_and_coupling_is_one_tenth():
    assert delta_composite(4, 3) == Rational(1, 40)
    assert c2_formula(4, 3) == Rational(1, 10)
    assert c2_formula(4, 3) != Rational(1, 40)      # the coupling is NOT Delta


def test_composite_4_3_matches_independent_delta():
    """The composite (2/5)/Omega^4(0) coincides with Delta = 1/40 from the
    independent g_3^Dirac degeneracy count -- an Observation, since the
    composite is not a matrix element (the coupling is 1/10)."""
    comp = delta_composite(4, 3)              # 1/40 from (2/5)(1/16)
    g3, delta = delta_inverse_identity()      # (40, 1/40) from g_3^Dirac count
    assert g3 == 40
    assert delta == Rational(1, 40)
    assert comp == delta


def test_g3_dirac_is_genuinely_independent():
    """g_3^Dirac = 2(n+1)(n+2)|_{n=3} = 40 contains NO 1/16 and NO c^2 formula
    — it is a pure spinor-shell degeneracy count. Confirms non-circularity.
    """
    g3 = dirac_degeneracy(3, sector="dirac", convention="ch")
    assert g3 == 2 * (3 + 1) * (3 + 2) == 40
    # And the composite route does not reference degeneracy: it is built from
    # the Casimir ratio and the Fock Jacobian. Different objects, same 1/40.
    assert delta_composite(4, 3) == Rational(1, g3)


def test_c2_formula_not_trivially_constant():
    """Guard against tautology: c^2 varies with (n,l); it is 1/4 at l=0 and
    1/10 at (4,3). A constant formula would make the composite vacuous."""
    vals = {(n, l): c2_formula(n, l) for n in range(1, 5) for l in range(n)}
    assert vals[(1, 0)] == Rational(1, 4)
    assert vals[(4, 3)] == Rational(1, 10)
    assert len(set(vals.values())) > 1        # genuinely non-constant
