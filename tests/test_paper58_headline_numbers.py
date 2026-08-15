"""Pin the two OWED headline numbers Paper 58 leans on.

The 2026-08-15 full /qa certifying run found Paper 58's two most-quotable numbers
were backed only by exploratory ``debug/`` drivers (matrix-logged OWED), not by a
pytest assertion -- which is exactly how a single-digit corruption could reach the
paper text and survive the whole suite (the run's two number-seeds did precisely
that). These two ``@slow`` tests close that gap:

  (AA|BB) decided census  195/195 nonzero, 0 Gaunt-zero, 0 accidental
                          (Z_A=3, Z_B=1, n_max=2, R=3; cf. debug/step2_decided_census.py)
  H2 native FCI energy    -1.106556606091 Ha
                          (R=1.4, 1s-per-centre; cf. debug/step1_native_molecule.py)

The census leg is geovac-only (stable regardless of debug/ pruning). The energy leg
drives the native two-center pipeline through the N3b census helpers and skips
gracefully if those exploratory modules are absent.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

import pytest

sp = pytest.importorskip("sympy")

REPO = Path(__file__).resolve().parents[1]

from geovac.two_center_eri import (  # noqa: E402
    R_s, aabb_closed_form, multipole_decomposition,
)


# ---------------------------------------------------- (AA|BB) DECIDED census: 195
# geovac-only replication of debug/step2_decided_census.py on Paper 58's own
# config, pinning the DECIDED (not merely counted) headline: 195 permitted, all
# 195 proven genuinely nonzero by Lindemann separation, 0 missed symmetry zeros,
# 0 accidental zeros.

_ZA, _ZB, _RCEN = Fraction(3), Fraction(1), sp.Integer(3)
_ORBS = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]


def _side_permitted(a, b, Mv):
    l1, l2 = a[1], b[1]
    return any((l1 + l2 + L) % 2 == 0 and abs(Mv) <= L
               for L in range(abs(l1 - l2), l1 + l2 + 1))


def _gaunt_all_zero(a, b, c, d):
    tA = multipole_decomposition(_ZA, *a, _ZA, *b)
    tB = multipole_decomposition(_ZB, *c, _ZB, *d)
    for _LA, MA, gA, _r, _bb in tA:
        for _LB, MB, gB, _r2, _b2 in tB:
            if MA + MB == 0 and sp.simplify(gA * gB) != 0:
                return False
    return True


def _decide_zero(expr, Rval):
    """Lindemann decision: group by exponential rate, require every A_j(R) = 0."""
    groups: dict = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        rate, coeff = sp.Integer(0), sp.Integer(1)
        for f in sp.Mul.make_args(term):
            if isinstance(f, sp.exp):
                arg = sp.expand(f.args[0])
                d = -sp.diff(arg, R_s)
                rate += d
                coeff *= sp.exp(sp.expand(arg + d * R_s))
            else:
                coeff *= f
        key = sp.nsimplify(rate)
        groups[key] = groups.get(key, 0) + coeff
    return all(sp.simplify(c.subs(R_s, Rval)) == 0 for c in groups.values())


@pytest.mark.slow
def test_paper58_aabb_decided_census_is_195_of_195():
    permitted = []
    for p, q, r, s in product(range(len(_ORBS)), repeat=4):
        a, b, c, d = _ORBS[p], _ORBS[q], _ORBS[r], _ORBS[s]
        MA, MB = b[2] - a[2], d[2] - c[2]
        if MA + MB != 0:
            continue
        if _side_permitted(a, b, MA) and _side_permitted(c, d, MB):
            permitted.append((a, b, c, d))
    assert len(permitted) == 195, "permitted count drifted from Paper 58's 195"

    n_nonzero = n_gaunt = n_accidental = 0
    for a, b, c, d in permitted:
        if _gaunt_all_zero(a, b, c, d):
            n_gaunt += 1
            continue
        e = aabb_closed_form(_ZA, a, b, _ZB, c, d)
        if _decide_zero(e, _RCEN):
            n_accidental += 1
        else:
            n_nonzero += 1
    assert (n_nonzero, n_gaunt, n_accidental) == (195, 0, 0)


# ---------------------------------------------- native H2 end-to-end FCI energy
@pytest.mark.slow
def test_paper58_native_h2_fci_energy():
    """Paper 58 Sec. II compose check: H2 at R=1.4, 1s-per-centre, native S+h+g
    -> FCI + V_NN = -1.106556606091 Ha (the value the Gaussian-fit path converges
    onto). Skips if the exploratory native-engine drivers are absent."""
    sys.path.insert(0, str(REPO / "debug"))
    step1 = pytest.importorskip("step1_native_molecule")
    S, h = step1.build_S_h()
    g = step1.build_g(verbose=False)
    energy = step1.solve(S, h, g)
    assert energy == pytest.approx(-1.106556606091, abs=1e-9)
