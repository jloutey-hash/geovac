"""Shim: moved to ``geovac.lih_r12ci.pi_m16`` (LiH productionization step 2c, 2026-10-04).

Importing this module returns the promoted module object itself, so old drivers, private
names, module globals and monkeypatches keep working; running it as a script runs the
promoted module's ``__main__`` block (same argv).
"""
import os as _os
import sys as _sys

# The original pinned BLAS/OpenMP to one thread BEFORE numpy was imported; keep that for drivers
# that import this shim first (the promoted module's own setdefault runs after geovac imports numpy).
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_v, "1")

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _ROOT not in _sys.path:
    _sys.path.insert(0, _ROOT)

if __name__ == "__main__":
    import runpy as _runpy
    _runpy.run_module("geovac.lih_r12ci.pi_m16", run_name="__main__", alter_sys=True)
else:
    import geovac.lih_r12ci.pi_m16 as _m
    _sys.modules[__name__] = _m
