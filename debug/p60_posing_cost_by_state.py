"""Reference-free test of the split-shell mechanism.

Posing cost = E_iso(state k) - min_lambda E_var(state k) in the SAME span.
Needs no known limit.  If the ground state is starved because both electrons
are pinned to one exponent, its posing cost should dwarf the excited state's.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S
from debug.p60_variational_probe import build

lmax = int(sys.argv[1]); Z = 2.0


def var_levels(Smat, T, W, G, lam, nlev, tol=1e-10):
    H = lam ** 2 * T + lam * (-Z * W + G)
    w, V = np.linalg.eigh(Smat)
    keep = w > tol * w.max()
    X = V[:, keep] / np.sqrt(w[keep])
    return np.linalg.eigvalsh(X.T @ H @ X)[:nlev]


print("nmax   K  | state | E_iso        E_var(lam*)   posing cost (mHa)")
for nmax in [int(x) for x in sys.argv[2:]]:
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    Smat, T, W, G, K, asym = build(nmax, lmax)
    M = S.build_M(S.build_configs(E.family(nmax, lmax)), Z=Z)
    p = np.sort(np.linalg.eigvalsh(M))[::-1]
    grid = np.linspace(0.5, 40, 700)
    NL = int(os.environ.get("NLEV", "2"))
    lev = np.array([var_levels(Smat, T, W, G, L, NL) for L in grid])
    for k in range(NL):
        e_iso = -p[k] ** 2 / 2
        e_var = lev[:, k].min()
        print("%4d %4d  |  %s  | %.7f   %.7f   %+8.3f"
              % (nmax, K, "gnd " if k == 0 else "%d^1S"%(k+1), e_iso, e_var, (e_iso - e_var) * 1000))
