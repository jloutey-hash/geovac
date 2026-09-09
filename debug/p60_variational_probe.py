"""COMPATIBILITY SHIM + driver -- the variational assembly now lives in ``geovac/``.

DIAGNOSTIC: is Paper 60's energy floor a property of the SPAN or of the POSING?

Builds a genuine variational CI over the EXACT same Goscinskian configuration
span the isoenergetic secular equation uses, with the global scale lambda
optimized (the variational analogue of the isoenergetic p_kappa).

The assembly itself (``build``, ``var_energy``, and the ``radial_1r`` /
``u_terms`` primitives beneath them) was promoted verbatim to
:mod:`geovac.sturmian_variational` on 2026-09-08, because
``tests/test_paper60_scale_lock.py`` is the only independent second route
backing Paper 60's ``eq:scale_lock`` and ``debug/`` is prunable by the CLAUDE.md
SS9 clean-room rule (gate C22 check D).  This file re-exports it, per the SS14
redirect-before-archive rule, so the sibling ``debug/p60_*.py`` drivers keep
working unchanged -- and keeps its own CLI below.

Usage: python debug/p60_variational_probe.py NMAX LMAX [NPTS] [BOXC]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import debug.p60_engine as E                                    # noqa: E402,F401
import geovac.sturmian_secular as S                             # noqa: E402
from geovac.sturmian_variational import (                       # noqa: E402,F401
    build,
    radial_1r,
    radial_overlap_c,
    u_terms,
    var_energy,
    var_levels,
)

EXACT = -2.903724377
S_LIMIT = -2.879028767


if __name__ == "__main__":
    nmax = int(sys.argv[1]); lmax = int(sys.argv[2])
    npts = int(sys.argv[3]) if len(sys.argv) > 3 else 24000
    boxc = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    Z = 2.0
    E.set_grid(max(80.0, boxc * nmax * nmax), npts, "grade", 2.0)
    Smat, T, W, G, K, asym = build(nmax, lmax)
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda L: var_energy(Smat, T, W, G, Z, L),
                          bounds=(0.3, 40.0), method="bounded",
                          options=dict(xatol=1e-7))
    Evar, lam = res.fun, res.x
    e_iso, onenorm, _, _ = S.solve(E.family(nmax, lmax), Z=Z)
    ref = S_LIMIT if lmax == 0 else EXACT
    lab = "s-limit" if lmax == 0 else "exact"
    print("nmax=%2d lmax=%d K=%4d  cond(S)=%8.1f  bra/ket asym=%.2e" % (nmax, lmax, K, np.linalg.cond(Smat), asym))
    print("   isoenergetic (metric-free) E = %.7f   gap to %s = %8.3f mHa" % (e_iso, lab, (e_iso - ref) * 1000))
    print("   VARIATIONAL CI, same span    E = %.7f   gap to %s = %8.3f mHa   (lambda*=%.4f)" % (Evar, lab, (Evar - ref) * 1000, lam))
    print("   posing cost = %.3f mHa" % ((e_iso - Evar) * 1000))
