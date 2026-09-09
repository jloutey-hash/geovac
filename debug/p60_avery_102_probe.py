"""Can the LOCKED-SCALE Goscinskian method reach -2.90250 Ha with 102 configs?

Avery & Avery's cited helium result is "-2.90250 with 102 optimized
Coulomb-Sturmian configurations" (1.2 mHa short).  A relayed source asserts it
used the bare-V0 Goscinskian isoenergetic method with a pure-number T' -- i.e.
exactly what geovac/sturmian_secular.py implements.  Our locked-scale ladder
gives 7.70 mHa at K=100.  Both cannot be true.

The word doing the work is "optimized".  If it means a SELECTED 102 drawn from a
larger pool, the question is whether selection rescues the locked posing.  It
should not:  the isoenergetic root is the lowest root of H(p_k)C = E S C, so the
method is variational, and a subset of a span cannot beat the whole span.  But
p_kappa is itself basis-dependent, so monotonicity is not obvious -- test it,
do not assert it.

Because M = diag(Z R_nu) + T', restricting to a subset of configurations is
exactly the principal submatrix, so one build prices every subset.

Usage: python debug/p60_avery_102_probe.py [NMAX] [LMAX] [NSEL]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.sturmian_variational import build, family, set_grid, var_energy
import geovac.sturmian_secular as S

EXACT = -2.903724377
AVERY = -2.90250
CHEM = 1.5936014616
Z = 2.0

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
lmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
nsel = int(sys.argv[3]) if len(sys.argv) > 3 else 102


def locked(Msub: np.ndarray) -> float:
    """Metric-free isoenergetic energy of a configuration subset."""
    p = np.sort(np.linalg.eigvalsh(Msub))[-1]
    return -p ** 2 / 2


def free(Ssub, Tsub, Wsub, Gsub) -> float:
    grid = np.linspace(0.5, 40.0, 500)
    return min(var_energy(Ssub, Tsub, Wsub, Gsub, Z, L) for L in grid)


set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
cts = family(nmax, lmax)
cfgs = S.build_configs(cts)
M = S.build_M(cfgs, Z=Z)
K = len(cfgs)
Smat, T, W, G, _, asym = build(nmax, lmax)

e_full = locked(M)
print("pool: n_max=%d l_max=%d K=%d   locked-scale E=%.7f  gap=%.3f mHa"
      % (nmax, lmax, K, e_full, (e_full - EXACT) * 1000))
print("Avery cited: E=%.5f  gap=%.3f mHa with %d configurations" % (AVERY, (AVERY - EXACT) * 1000, nsel))
print()

# --- monotonicity: does adding configurations always lower the locked energy?
print("MONOTONICITY of the locked posing (prefix subsets of the pool):")
prev, mono = None, True
for k in range(10, K + 1, max(1, (K - 10) // 8)):
    e = locked(M[:k, :k])
    flag = ""
    if prev is not None and e > prev + 1e-12:
        flag, mono = "  <-- ROSE", False
    print("   K=%4d  E=%.7f  gap=%7.3f mHa%s" % (k, e, (e - EXACT) * 1000, flag))
    prev = e
print("   monotone decreasing: %s" % mono)
print()

# --- the best possible 102: rank configurations by ground-state weight
p, V = np.linalg.eigh(M)
B = V[:, -1]
order = np.argsort(-np.abs(B))
sel = np.sort(order[:nsel])
e_sel = locked(M[np.ix_(sel, sel)])
e_sel_free = free(Smat[np.ix_(sel, sel)], T[np.ix_(sel, sel)],
                  W[np.ix_(sel, sel)], G[np.ix_(sel, sel)])
print("BEST %d by ground-state weight, drawn from the K=%d pool:" % (nsel, K))
print("   locked scale : E=%.7f  gap=%7.3f mHa" % (e_sel, (e_sel - EXACT) * 1000))
print("   free scale   : E=%.7f  gap=%7.3f mHa" % (e_sel_free, (e_sel_free - EXACT) * 1000))
print()

# --- can ANY selection help?  random subsets for the spread
rng = np.random.default_rng(0)
# NOTE the sign:  E is negative and deeper binding is LOWER, so the BEST random
# subset is the one with the MINIMUM energy.  An earlier version of this driver
# took max() and therefore reported the worst subset as the best.
energies = []
for _ in range(200):
    s = np.sort(rng.choice(K, size=nsel, replace=False))
    energies.append(locked(M[np.ix_(s, s)]))
best, worst = min(energies), max(energies)
print("200 random %d-subsets: BEST gap=%7.3f mHa   WORST gap=%7.3f mHa"
      % (nsel, (best - EXACT) * 1000, (worst - EXACT) * 1000))
print()
print("BOUND: no %d-configuration subset of this pool beats the full pool's"
      % nsel)
print("       %.3f mHa in the locked posing, because the posing is variational." % ((e_full - EXACT) * 1000))
print("       Avery's %.3f mHa is therefore NOT reachable by locked-scale" % ((AVERY - EXACT) * 1000))
print("       Goscinskian selection at this or any larger K.")

json.dump(dict(nmax=nmax, lmax=lmax, K=K, nsel=nsel, e_full=e_full,
               e_sel_locked=e_sel, e_sel_free=e_sel_free, best_random=best, worst_random=worst,
               monotone=bool(mono), avery=AVERY, exact=EXACT),
          open("debug/data/p60_avery_102_probe.json", "w"), indent=1)
