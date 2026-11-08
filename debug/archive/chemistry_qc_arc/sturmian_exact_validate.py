"""Validation drivers for debug/sturmian_exact_slater.py.

(a) equal-exponent case vs geovac/hypergeometric_slater.py (independent algorithm);
(b) mixed-exponent case vs geovac/sturmian_secular.py as the radial box is opened.
"""
from __future__ import annotations

import time

import numpy as np

import debug.sturmian_exact_slater as ex
import geovac.sturmian_secular as ss


def regrid(rmax: float, ngrid: int) -> None:
    """Rebuild the sturmian_secular module grid in place."""
    ss.R_MAX = rmax
    ss.N_GRID = ngrid
    ss.r = np.linspace(1e-7, rmax, ngrid)
    ss.dr = ss.r[1] - ss.r[0]
    ss.r2 = ss.r * ss.r
    ss.reset_caches()


def grid_M(tuples, rmax: float, ngrid: int, Z: float = 2.0):
    regrid(rmax, ngrid)
    cfgs = ss.build_configs(tuples)
    return ss.build_M(cfgs, Z), cfgs


def exact_M(tuples, Z: float = 2.0):
    cfgs = ex.build_exact_configs(tuples)
    return ex.build_exact_M(cfgs, Z), cfgs


def compare(tuples, boxes, label: str = "") -> None:
    Me, _ = exact_M(tuples)
    print("=== %s   K=%d ===" % (label, len(tuples)))
    print("  exact:  ||M||_1 = %.10f   E0 = %.10f"
          % (np.abs(Me).sum(), -np.sort(np.linalg.eigvalsh(Me))[-1] ** 2 / 2))
    for (rmax, ngrid) in boxes:
        t0 = time.time()
        Mg, _ = grid_M(tuples, rmax, ngrid)
        d = np.abs(Mg - Me)
        rel = np.abs(Mg - Me) / np.maximum(np.abs(Me), 1e-30)
        print("  grid rmax=%6.0f n=%6d: ||M||_1=%.10f  E0=%.10f  maxabs=%.3e "
              "maxrel=%.3e  1-norm rel err=%.3e  (%.1fs)"
              % (rmax, ngrid, np.abs(Mg).sum(),
                 -np.sort(np.linalg.eigvalsh(Mg))[-1] ** 2 / 2,
                 d.max(), rel.max(),
                 abs(np.abs(Mg).sum() - np.abs(Me).sum()) / np.abs(Me).sum(),
                 time.time() - t0))
