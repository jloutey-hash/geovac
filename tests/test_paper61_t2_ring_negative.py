r"""Paper 61 sec:modular (eq:kw paragraph) and the conclusion: the guarded integer-relation
search on the 66-digit collinear observable T2 is DECISIVELY NEGATIVE in the corrected
weight-<=3 ring {pi, K(1/2)^{+-1}, G} at 64 working digits, two precisions.

Owed guard (docs/qa/group3.debt.md, open since 2026-09-07): the negative had no test.

What this test pins (self-contained; the value is the certified 66-digit entry of
benchmarks/certified_reference, quoted there, not recomputed):

  ring R = {1} U {pi^a vp^b G^c : 1 <= a + |b| + 2c <= 3, a,c >= 0}   (20 elements,
           vp = K(1/2) = Gamma(1/4)^2 / (4 sqrt pi), G = Catalan = beta(2)),
  exactly the pre-registered "T-2 ring" of the driver (leg J of debug/beta2_t2_pslq.py);
  targets T2 and the natural period W = T2 * pi / 8 (both legs the driver ran).

  (1) POSITIVE CONTROL: a height-<=10 element of the SAME 20-element ring, built at the same
      precision, is recovered by the same PSLQ call (so a "no relation" is not a search that
      cannot see this ring).
  (2) NEGATIVE: at dps 60 and 64, PSLQ finds no relation of height <= 100 for T2 or W.  The
      paper states height <= 10; the reviewer found it holds to 100, so the guard checks the
      stronger bound (still inside the false-positive budget 10^(dps/n) = 10^3.0 .. 10^3.2).
  (3) The height budget is real: 10^(dps/n) > 100 at both precisions, so (2) is not
      UNDERPOWERED by construction.

Wrong answer it excludes: "T2 (or W) is a low-height element of the weight-<=3 Catalan-
completed CM ring" -- planting T2 := a ring element makes (2) fail.
"""
from __future__ import annotations

import itertools

import mpmath as mp
import pytest

from benchmarks.certified_reference.entries_t2 import T2_VALUE_66

HEIGHT = 100
PRECISIONS = (60, 64)


def _ring(dps: int) -> dict:
    mp.mp.dps = dps + 40
    pi = +mp.pi
    vp = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
    G = +mp.catalan
    mp.mp.dps = dps
    ring = {"1": mp.mpf(1)}
    for a, b, c in itertools.product(range(0, 4), range(-3, 4), range(0, 2)):
        wt = a + abs(b) + 2 * c
        if 1 <= wt <= 3:
            ring[f"pi^{a}*vp^{b}*G^{c}"] = +(pi ** a * vp ** b * G ** c)
    return ring


def _pslq(target, ring: dict, dps: int):
    mp.mp.dps = dps
    vec = [+target] + list(ring.values())
    return mp.pslq(vec, tol=mp.mpf(10) ** -(dps - 3), maxcoeff=HEIGHT, maxsteps=10 ** 6)


def _targets(dps: int) -> dict:
    mp.mp.dps = dps + 20
    t2 = mp.mpf(T2_VALUE_66)
    return {"T2": t2, "W": t2 * mp.pi / 8}


def test_ring_has_twenty_elements_and_budget():
    for dps in PRECISIONS:
        ring = _ring(dps)
        assert len(ring) == 20
        # false-positive height of an n-element ring at this precision exceeds the search height
        assert mp.mpf(10) ** (mp.mpf(dps) / len(ring)) > HEIGHT


def test_ring_generators_pinned_independently():
    """The control is built from the ring itself, so pin each generator by an
    independent route: vp = K(1/2) via the AGM (mp.ellipk, m = k^2 = 1/2) and
    G = beta(2) via the Dirichlet series sum (-1)^n/(2n+1)^2 (mp.nsum)."""
    dps = max(PRECISIONS)
    ring = _ring(dps)
    mp.mp.dps = dps + 10
    vp_ref = mp.ellipk(mp.mpf(1) / 2)
    G_ref = mp.nsum(lambda n: (-1) ** n / (2 * n + 1) ** 2, [0, mp.inf])
    tol = mp.mpf(10) ** -(dps - 5)
    assert abs(ring["pi^0*vp^1*G^0"] - vp_ref) < tol
    assert abs(ring["pi^0*vp^0*G^1"] - G_ref) < tol
    assert abs(ring["pi^1*vp^0*G^0"] - 4 * mp.atan(1)) < tol


@pytest.mark.parametrize("dps", PRECISIONS)
def test_positive_control_recovered(dps):
    ring = _ring(dps)
    keys = list(ring)
    # a height-<=10 combination using a weight-3 G-requiring monomial and a 1/vp term
    coeffs = {"1": -3, "pi^1*vp^0*G^1": 5, "pi^0*vp^-1*G^1": -2, "pi^2*vp^1*G^0": 7}
    mp.mp.dps = dps
    ctrl = sum(c * ring[k] for k, c in coeffs.items()) / 9
    rel = _pslq(ctrl, ring, dps)
    assert rel is not None, "positive control not recovered: the search cannot see this ring"
    # the recovered relation must be the planted one (up to sign)
    expect = [9] + [-coeffs.get(k, 0) for k in keys]
    assert rel == expect or rel == [-x for x in expect], rel


@pytest.mark.slow
@pytest.mark.parametrize("dps", PRECISIONS)
@pytest.mark.parametrize("name", ["T2", "W"])
def test_t2_no_low_height_relation(dps, name):
    ring = _ring(dps)
    target = _targets(dps)[name]
    rel = _pslq(target, ring, dps)
    assert rel is None, (name, dps, rel)
