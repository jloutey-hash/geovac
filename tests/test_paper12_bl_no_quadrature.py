"""Guards for Paper 12: the B_l table is quadrature-free, and the numerical
V_ee comparator column of tab:convergence.

1. Quadrature spy.  Paper 12 states the Neumann V_ee route evaluates the
   auxiliary integral B_l(p, alpha) = int_1^inf xi^p e^{-alpha xi} Q_l(xi)
   in closed form (closed-form seeds + recurrence), with no numerical
   integration.  geovac/neumann_vee.compute_Bl_table previously built the
   table by scipy.integrate.quad of the lqmn integrand.  This guard replaces
   scipy.integrate.quad (and mpmath.quad, the other integrator in reach) with
   a function that raises, empties the B_l cache, and builds the table.
   Wrong answer rejected: any quadrature call on the B_l path.  Fire-tested
   by planting a quad call inside compute_Bl_table.
   The closed-form values themselves are checked against an independent quad
   in tests/test_neumann_vee.py::test_Bl_against_quad.

2. Numerical comparator (tab:convergence, numerical column, N = 27 row,
   j = 2, l = 2, p = 0, alpha = 1.0, R = 1.4011 bohr, 30 x 20 (xi, eta)
   quadrature grid): E_num = -1.138836 Ha, D_e = 79.6% of 0.174475 Ha.
   The Neumann column (-1.160961, 92.2%) is pinned elsewhere; this pins the
   comparator the paper's "algebraic beats numerical" sentence rests on.
   Slow (the numerical 5D V_ee build takes several minutes).
"""
from __future__ import annotations

import mpmath
import numpy as np
import pytest
import scipy.integrate

import geovac.neumann_vee as nv


class _QuadCalled(AssertionError):
    pass


def _forbidden(*args, **kwargs):
    raise _QuadCalled("quadrature called on the closed-form B_l path")


def test_compute_Bl_table_calls_no_quadrature(monkeypatch):
    monkeypatch.setattr(scipy.integrate, "quad", _forbidden)
    monkeypatch.setattr(mpmath, "quad", _forbidden)
    monkeypatch.setattr(mpmath.mp, "quad", _forbidden)
    monkeypatch.setattr(nv, "_Bl_cache", {})
    # a few (l_max, p_max, alpha) shapes, including a high-l table
    for l_max, p_max, alpha in ((4, 3, 2.0), (12, 6, 1.3), (20, 4, 0.7)):
        Bl = nv.compute_Bl_table(l_max, p_max, alpha)
        assert Bl.shape == (l_max + 1, p_max + 1)
        assert np.all(np.isfinite(Bl))
        assert np.all(Bl[0] > 0.0)  # Q_0 > 0 on (1, inf)


def _integrator_frames(fn):
    """Run fn() under a profiler and return the integrator frames it entered.

    The attribute spy above cannot see a quadrature routine bound by name before
    the monkeypatch (e.g. a module-level `from scipy.integrate import quad as q`).
    A profiler sees every Python frame entered, however the name was bound, so it
    catches any call into scipy.integrate or mpmath's quadrature module.
    """
    import sys
    hits = []

    def prof(frame, event, arg):
        if event == "call":
            f = frame.f_code.co_filename.replace("\\", "/")
            if "/scipy/integrate/" in f or f.endswith("/mpmath/calculus/quadrature.py"):
                hits.append((f, frame.f_code.co_name))

    sys.setprofile(prof)
    try:
        fn()
    finally:
        sys.setprofile(None)
    return hits


def test_compute_Bl_table_enters_no_integrator_frame(monkeypatch):
    monkeypatch.setattr(nv, "_Bl_cache", {})

    def build():
        for l_max, p_max, alpha in ((4, 3, 2.0), (12, 6, 1.3)):
            nv.compute_Bl_table(l_max, p_max, alpha)

    hits = _integrator_frames(build)
    assert not hits, f"integrator frames entered on the B_l path: {sorted(set(hits))[:5]}"


def test_profiler_spy_is_live():
    """The profiler must see a quad bound by name before any patching."""
    from scipy.integrate import quad as bound_early
    hits = _integrator_frames(lambda: bound_early(lambda x: x, 0.0, 1.0))
    assert hits


def test_spy_is_live(monkeypatch):
    """The spy must actually intercept a quad call made the way the old code
    did it (a function-local `from scipy.integrate import quad`)."""
    monkeypatch.setattr(scipy.integrate, "quad", _forbidden)

    def old_style():
        from scipy.integrate import quad
        return quad(lambda x: x, 0.0, 1.0)

    with pytest.raises(_QuadCalled):
        old_style()


@pytest.mark.slow
def test_numerical_vee_comparator_N27():
    """tab:convergence numerical column, N = 27: 79.6% D_e (E = -1.138836)."""
    import geovac.hylleraas as h
    basis = h.generate_basis(j_max=2, l_max=2, p_max=0, alpha=1.0)
    assert len(basis) == 27
    grids = h.build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)
    r = h.solve_hylleraas(basis, 1.4011, grids, verbose=False,
                          vee_method="numerical")
    E = float(r["E_total"])
    pct = (-1.0 - E) / 0.174475 * 100.0
    # Recomputed 2026-11-04 on the 30 x 20 grid: -1.1388188 (79.56%); the
    # table prints -1.138836, 1.7e-5 Ha lower (same 79.6% after rounding).
    assert abs(E - (-1.138836)) < 3e-5, f"E_num = {E:.7f}"
    assert abs(pct - 79.6) < 0.1, f"D_e = {pct:.3f}%"
