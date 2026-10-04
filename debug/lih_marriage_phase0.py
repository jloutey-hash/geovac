"""Shim: moved to ``geovac.lih_r12ci.ref_tensor`` (LiH productionization step 2c, 2026-10-04).

Importing this module returns the promoted module object itself, so old drivers, private
names, module globals and monkeypatches keep working; running it as a script runs the
promoted module's ``__main__`` block (same argv).
"""
import os as _os
import sys as _sys

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _ROOT not in _sys.path:
    _sys.path.insert(0, _ROOT)

if __name__ == "__main__":
    import runpy as _runpy
    _runpy.run_module("geovac.lih_r12ci.ref_tensor", run_name="__main__", alter_sys=True)
else:
    import geovac.lih_r12ci.ref_tensor as _m
    _sys.modules[__name__] = _m
