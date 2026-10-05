"""DEPRECATED source-level overrides: the geminal range now has a real API; only the kernel GRID
(NXI, NETA, NPHI) and the triangle cutoffs (MMAX, LMAX) still need a source override.

Geminal range.  Use ``geovac.lih_r12ci.configure_geminal(gamma)`` (= ``kernels.configure_geminal``):
every consumer reads ``kernels.GAM`` at call / construction time, and ``gVee``'s gamma tables are built
at import for the current gamma or rebuilt on change.  ``install_kernels(gam, None)`` is kept as a thin
wrapper that just calls ``configure_geminal(gam)``; it reproduces the old source-override results
bit-identically (tests/data/lih_golden items ``pi_marriage_exp1.5_*``).

Grid / triangle cutoffs (still import-time bindings).  ``NXI, NETA, NPHI`` size the grid arrays of
kernels.py that every consumer from-imports and builds tables from at import; ``MMAX, LMAX`` size
triangle.py's import-time Legendre tables and are bound by value into enumerators / pi_m16.
``install_kernels(gam, grid)`` with a grid and ``install_triangle`` therefore still read the module
SOURCE, substitute the constant line, and exec it as the module BEFORE any consumer imports it:
  * must be called before ANY import of the overridden module (asserted);
  * it REPLACES ``sys.modules[...]`` (and the package attribute) for the rest of the process; there is
    no uninstall -- use a fresh process per grid.

PROMOTED from debug/lih_exp_geminal_run.py in the LiH productionization step 2c (2026-10-04).
Paper 19 result it backs: the geminal-range table (tab:geminal_gamma) and the corpus-best
pi-marriage energy with the Slater geminal e^{-1.5 r}.
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
from types import ModuleType
from typing import Optional, Sequence

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))


def install_kernels(gam: float, grid: Optional[Sequence[int]]) -> ModuleType:
    """Set the geminal range to ``gam`` and (optionally) source-override the kernel grid.

    ``grid is None``: equivalent to ``configure_geminal(gam)`` (callable at any time) and returns the
    ordinary ``geovac.lih_r12ci.kernels`` module.  ``grid = (NXI, NETA, NPHI)``: exec kernels.py with
    the grid line substituted as ``geovac.lih_r12ci.kernels`` (must precede every kernels import),
    then ``configure_geminal(gam)`` on it.  Returns the kernels module in use."""
    if grid is None:
        import geovac.lih_r12ci.kernels as KG
        KG.configure_geminal(gam)
        return KG
    import geovac.lih_r12ci as pkg
    name = "geovac.lih_r12ci.kernels"
    assert name not in sys.modules, "kernels already imported -- grid override would not propagate"
    path = os.path.join(_PKG_DIR, "kernels.py")
    src = open(path, encoding="utf-8").read()
    nxi, neta, nphi = grid
    src, n2 = re.subn(r"^NXI, NETA, NPHI = 72, 44, 28\b",
                      f"NXI, NETA, NPHI = {nxi}, {neta}, {nphi}", src, flags=re.M)
    assert n2 == 1, "grid line not found exactly once"
    spec = importlib.util.spec_from_loader(name, loader=None, origin=path)
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = path
    mod.__package__ = "geovac.lih_r12ci"
    sys.modules[name] = mod
    exec(compile(src, path, "exec"), mod.__dict__)
    pkg.kernels = mod
    mod.configure_geminal(gam)
    return mod


# the geminal-range override under its descriptive name
install_geminal = install_kernels


def install_triangle(mmax: Optional[int], lmax: Optional[int]) -> ModuleType:
    """Source-override triangle.py's MMAX / LMAX (module constants that size import-time tables and
    are bound BY VALUE into the enumerators / pi_m16 at their import) -- must run before those import."""
    import geovac.lih_r12ci as pkg
    name = "geovac.lih_r12ci.triangle"
    assert name not in sys.modules
    path = os.path.join(_PKG_DIR, "triangle.py")
    src = open(path, encoding="utf-8").read()
    if mmax is not None:
        src, n = re.subn(r"^MMAX = 4", f"MMAX = {int(mmax)}", src, flags=re.M); assert n == 1
    if lmax is not None:
        src, n = re.subn(r"^LMAX = 34", f"LMAX = {int(lmax)}", src, flags=re.M); assert n == 1
    spec = importlib.util.spec_from_loader(name, loader=None, origin=path)
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = path; mod.__package__ = "geovac.lih_r12ci"
    sys.modules[name] = mod
    exec(compile(src, path, "exec"), mod.__dict__)
    pkg.triangle = mod
    return mod
