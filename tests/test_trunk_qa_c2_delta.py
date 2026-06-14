"""
TRUNK QA — Claim 4: c^2(4,3) = 1/40 = Delta (Paper 7 -> Paper 2).

Paper 7 Eq fock_coupling: c^2(n,l) = (1/16)[1 - l(l+1)/(n(n+1))].
At the Paper 2 cutoff (highest l in the n=4 shell, l=3):
  c^2(4,3) = (1/16)(1 - 12/20) = (1/16)(2/5) = 1/40.
Paper 7 identifies this with Delta = 1/40, the boundary term in Paper 2's
K = pi(B + F - Delta), and with the Dirac degeneracy g_3^Dirac = 40.

Could-have-failed content:
  - c^2(4,3) computed from the (n,l) formula must equal 1/40 (exact rational),
  - and this must equal the INDEPENDENTLY-defined Delta from geovac.dirac_s3
    (which derives 1/40 from g_3^Dirac = 2(n+1)(n+2)|_{n=3} = 40, a spinor
    degeneracy count with no reference to the c^2 formula).
The match of two independent rational routes to 1/40 is the content; if the
c^2 formula gave anything else, the bridge would fail.
"""

from __future__ import annotations

import sympy as sp
from sympy import Rational

from geovac.dirac_s3 import delta_inverse_identity, dirac_degeneracy


def c2_formula(n: int, l: int) -> sp.Rational:
    return Rational(1, 16) * (1 - Rational(l * (l + 1), n * (n + 1)))


def test_c2_4_3_is_one_fortieth_exact():
    assert c2_formula(4, 3) == Rational(1, 40)
    assert isinstance(c2_formula(4, 3), Rational)


def test_c2_4_3_matches_independent_delta():
    """c^2(4,3) (Fock coupling route) == Delta (spinor degeneracy route)."""
    c2 = c2_formula(4, 3)                     # 1/40 from (1/16)(1 - 12/20)
    g3, delta = delta_inverse_identity()      # (40, 1/40) from g_3^Dirac count
    assert g3 == 40
    assert delta == Rational(1, 40)
    assert c2 == delta                        # two independent routes agree


def test_g3_dirac_is_genuinely_independent():
    """g_3^Dirac = 2(n+1)(n+2)|_{n=3} = 40 contains NO 1/16 and NO c^2 formula
    — it is a pure spinor-shell degeneracy count. Confirms non-circularity.
    """
    g3 = dirac_degeneracy(3, sector="dirac", convention="ch")
    assert g3 == 2 * (3 + 1) * (3 + 2) == 40
    # And the c^2 route does not reference degeneracy: it is built from the
    # Casimir ratio l(l+1)/(n(n+1)). Different mathematical objects, same 1/40.
    assert c2_formula(4, 3) == Rational(1, g3)


def test_c2_formula_not_trivially_constant():
    """Guard against tautology: c^2 varies with (n,l); it is 1/16 only at l=0
    and 1/40 at (4,3). A constant formula would make the bridge vacuous."""
    vals = {(n, l): c2_formula(n, l) for n in range(1, 5) for l in range(n)}
    assert vals[(1, 0)] == Rational(1, 16)
    assert vals[(4, 3)] == Rational(1, 40)
    assert len(set(vals.values())) > 1        # genuinely non-constant
