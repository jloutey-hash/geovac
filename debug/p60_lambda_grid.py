"""Integrity check: the lambda optimum by GRID scan, not Brent.

If the isoenergetic solution is variational (lowest root of H(lam)C = E S C at
lam = p_kappa), then min_lam E_var(lam) <= E_var(p_kappa) = E_iso ALWAYS.
A positive "posing cost" is therefore mandatory; the n_max=4 negative value in
the earlier sweep must be an optimizer failure.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S
from debug.p60_variational_probe import build, var_energy

EXACT, S_LIMIT = -2.903724377, -2.879028767
lmax = int(sys.argv[1]); ns = [int(x) for x in sys.argv[2:]]
Z = 2.0
for nmax in ns:
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    Smat, T, W, G, K, asym = build(nmax, lmax)
    e_iso, _, _, _ = S.solve(E.family(nmax, lmax), Z=Z)
    pk = np.sqrt(-2 * e_iso)
    grid = np.concatenate([np.linspace(0.5, 40, 800), [pk]])
    vals = np.array([var_energy(Smat, T, W, G, Z, L) for L in grid])
    i = vals.argmin()
    ref = S_LIMIT if lmax == 0 else EXACT
    print("nmax=%2d lmax=%d K=%4d | p_k=%.4f  E_iso=%.7f (gap %6.3f) | "
          "grid lam*=%6.3f E_var=%.7f (gap %6.3f) | posing cost %+7.3f mHa"
          % (nmax, lmax, K, pk, e_iso, (e_iso - ref) * 1000,
             grid[i], vals[i], (vals[i] - ref) * 1000, (e_iso - vals[i]) * 1000))
    assert vals[i] <= e_iso + 1e-9, "VARIATIONAL BOUND VIOLATED"
