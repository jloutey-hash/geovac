"""Source-level overrides of the geminal range GAM / the kernel grid / the triangle cutoffs.

PROMOTED from debug/lih_exp_geminal_run.py (``install_kernels``, ``install_triangle``) in the LiH
productionization step 2c (2026-10-04); the mechanism is unchanged and bit-identical against
``tests/data/lih_golden/golden.json`` (items ``pi_marriage_exp1.5_*``).

Why a source override and not a parameter.  ``GAM`` (and ``NXI, NETA, NPHI``) are module constants
of ``geovac/lih_r12ci/kernels.py`` that every consumer binds BY VALUE at its own import
(``from .kernels import GAM, ...``), and ``gVee`` even builds its phi-kernels ``Kf``/``_KS`` at import.
A post-import monkeypatch therefore would not reach them.  ``install_kernels`` reads kernels.py's
SOURCE, substitutes the constant line(s), and execs the result as the module
``geovac.lih_r12ci.kernels`` BEFORE any consumer imports it; every later ``from .kernels import ...``
then binds the new value.  (Owed refactor: a real ``gamma`` parameter threaded through
kernels/hT/gT/gVee/triangle/enumerators/pi_m16 -- not done here because every import-time
binding would have to move, which is not a minimal, bit-identity-preserving change.)

Contract / side effects:
  * must be called before ANY import of ``geovac.lih_r12ci.kernels`` (asserted) -- i.e. before
    importing gVee, hT, gT, hVee, triangle, enumerators, pi_m16, gvee_sparse, ...;
  * it REPLACES ``sys.modules['geovac.lih_r12ci.kernels']`` (and the package attribute) for the
    rest of the process; there is no uninstall -- use a fresh process per gamma;
  * the package ``__init__`` (``assembly``) imports kernels lazily, so importing
    ``geovac.lih_r12ci`` itself does not defeat the assertion.

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
    """Exec kernels.py with ``GAM = gam`` (and optionally ``NXI, NETA, NPHI = grid``) as
    ``geovac.lih_r12ci.kernels``.  Returns the installed module."""
    import geovac.lih_r12ci as pkg
    name = "geovac.lih_r12ci.kernels"
    assert name not in sys.modules, "kernels already imported -- override would not propagate"
    path = os.path.join(_PKG_DIR, "kernels.py")
    src = open(path, encoding="utf-8").read()
    src, n1 = re.subn(r"^GAM = 0\.5\b", f"GAM = {float(gam)!r}", src, flags=re.M)
    assert n1 == 1, "GAM line not found exactly once"
    if grid is not None:
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
