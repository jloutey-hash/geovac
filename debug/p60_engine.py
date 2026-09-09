"""COMPATIBILITY SHIM -- the Paper-60 grid harness now lives in ``geovac/``.

The machinery formerly defined here (``set_grid`` and the two non-uniform-safe
cumulative integrators it installs, ``family``, ``norms``) was promoted verbatim
to :mod:`geovac.sturmian_variational` on 2026-09-08, because
``tests/test_paper60_scale_lock.py`` is the only independent second route
backing Paper 60's ``eq:scale_lock`` and ``debug/`` is prunable by the CLAUDE.md
SS9 clean-room rule (gate C22 check D).

This file is kept, per the SS14 redirect-before-archive rule, so the existing
``debug/p60_*.py`` drivers keep working unchanged.  It re-exports only; add
nothing here.
"""
from __future__ import annotations

import geovac.sturmian_secular as S  # noqa: F401  (drivers reach through E.S)
from geovac.sturmian_variational import (  # noqa: F401
    _ORIG_FWD,
    _ORIG_REV,
    _fwd_nonuniform,
    _rev_nonuniform,
    family,
    norms,
    set_grid,
)

__all__ = ["set_grid", "family", "norms", "S"]
