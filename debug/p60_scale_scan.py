"""Is Paper 60's floor the FROZEN BASIS SCALE?

The isoenergetic posing pins the basis scale to Q_nu = p_kappa / R_nu, with
p_kappa = sqrt(-2E) the eigenvalue itself.  The variational optimum over the
same span wants a different scale.  This scans E_var(lambda) and reports it at
lambda = p_kappa against the isoenergetic answer.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S
from debug.p60_variational_probe import build, var_energy

EXACT = -2.903724377
S_LIMIT = -2.879028767

nmax = int(sys.argv[1]); lmax = int(sys.argv[2])
npts = int(sys.argv[3]) if len(sys.argv) > 3 else 24000
Z = 2.0
E.set_grid(max(80.0, 5.0 * nmax * nmax), npts, "grade", 2.0)
Smat, T, W, G, K, asym = build(nmax, lmax)
e_iso, _, _, _ = S.solve(E.family(nmax, lmax), Z=Z)
pk = np.sqrt(-2 * e_iso)
ref = S_LIMIT if lmax == 0 else EXACT
lab = "s-limit" if lmax == 0 else "exact"

from scipy.optimize import minimize_scalar
res = minimize_scalar(lambda L: var_energy(Smat, T, W, G, Z, L),
                      bounds=(0.3, 40.0), method="bounded", options=dict(xatol=1e-7))
print("nmax=%2d lmax=%d K=%4d   p_kappa = %.4f" % (nmax, lmax, K, pk))
print("  E_iso                    = %.7f   gap to %s = %7.3f mHa" % (e_iso, lab, (e_iso - ref) * 1000))
print("  E_var at lambda = p_kappa= %.7f   gap to %s = %7.3f mHa" % (var_energy(Smat, T, W, G, Z, pk), lab, (var_energy(Smat, T, W, G, Z, pk) - ref) * 1000))
print("  E_var at lambda* = %6.3f  = %.7f   gap to %s = %7.3f mHa" % (res.x, res.fun, lab, (res.fun - ref) * 1000))
print("  --> scale freedom is worth %.3f mHa; residual span deficit %.3f mHa"
      % ((var_energy(Smat, T, W, G, Z, pk) - res.fun) * 1000, (res.fun - ref) * 1000))
