"""Validation of the two-center Coulomb--Sturmian integral engine
(``geovac/sturmian_integrals.py``), promoted from the Paper-60 sprint drivers.

Each check pins a numerical engine output against a known closed form, so a
regression in the multipole/grid machinery is caught.  Backs the Paper-60
``sec:manyelectron`` claim that the enabling two-center ERIs are computable and
validated against the framework's exact closed forms.
"""
import numpy as np
import pytest

from geovac.sturmian_integrals import GoscinskianIntegrals, validate


@pytest.mark.slow
def test_sturmian_integrals_validate_against_closed_forms():
    """Engine matches all five closed-form references to <= ~2e-4."""
    res = validate()
    tol = {
        "self_A": 5e-4,          # <1s|1/rA|1s> = a
        "cross_B": 5e-4,         # <1s_A|1/rB|1s_A> Coulomb-of-1s closed form
        "eri_1c": 5e-4,          # (1s1s|1s1s) = 5a/8
        "overlap_mixed": 1e-6,   # mixed-scale one-center overlap (2 sqrt(ab)/(a+b))^3
        "eri_aabb": 5e-4,        # (AA|BB) vs geovac.two_center_eri.aabb_value
    }
    for name, (num, exact) in res.items():
        assert abs(num - exact) < tol[name], f"{name}: numeric {num} vs exact {exact}"
    # the mixed-scale overlap is essentially analytic on the grid
    assert abs(res["overlap_mixed"][0] - res["overlap_mixed"][1]) < 1e-8


def test_sturmian_integrals_one_center_eri_is_five_eighths():
    """Fast spot check: one-center (1s1s|1s1s) = 5a/8 at a=1 (=5/8)."""
    g = GoscinskianIntegrals(R=1.5, Lmax=16, nr=2000, nth=120, rmax=55.0)
    a = 1.0
    v = g.eri((1, "A", a), (1, "A", a), (1, "A", a), (1, "A", a))
    assert abs(v - 5 * a / 8) < 3e-3, f"(1s1s|1s1s) = {v}, expected {5*a/8}"


def test_sturmian_integrals_overlap_diagonal_normalized():
    """Self-overlap of an L2-normalized orbital is 1 to grid precision."""
    g = GoscinskianIntegrals(R=1.5, Lmax=12, nr=2500, nth=120, rmax=55.0)
    for n, a in ((1, 1.0), (2, 0.7), (1, 1.3)):
        s = g.overlap((n, "A", a), (n, "A", a))
        assert abs(s - 1.0) < 2e-3, f"<{n}|{n}> = {s} (a={a})"
