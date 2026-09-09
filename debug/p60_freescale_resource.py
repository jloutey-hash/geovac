"""What does FREEING THE SCALE cost in block-encoding terms?

Metric-free isoenergetic:  M C = p_kappa C,     ||M||_1 ~ K^0.82, no metric.
Free-scale variational:    H(lam) C = E S C,    metric S returns.

Measures, per K: E and gap for both; lambda*; cond(S); the entrywise 1-norm of
M, of H(lam*), and of the whitened Hhat = S^-1/2 H S^-1/2 whose eigenvalues are
the physical energies (the object a qubitized algorithm actually encodes).
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S
from debug.p60_variational_probe import build, var_energy
from scipy.optimize import minimize_scalar

EXACT, S_LIMIT = -2.903724377, -2.879028767
lmax = int(sys.argv[1]); ns = [int(x) for x in sys.argv[2:]]
Z, ref = 2.0, None
rows = []
for nmax in ns:
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    Smat, T, Wm, G, K, asym = build(nmax, lmax)
    Rn = np.array([c.Rnu for c in S.build_configs(E.family(nmax, lmax))])
    e_iso, m1, _, M = S.solve(E.family(nmax, lmax), Z=Z)
    res = minimize_scalar(lambda L: var_energy(Smat, T, Wm, G, Z, L),
                          bounds=(0.3, 40.0), method="bounded", options=dict(xatol=1e-7))
    lam = res.x
    H = lam ** 2 * T + lam * (-Z * Wm + G)
    w, V = np.linalg.eigh(Smat)
    Xi = V @ np.diag(w ** -0.5) @ V.T          # S^-1/2
    Hh = Xi @ H @ Xi
    ref = S_LIMIT if lmax == 0 else EXACT
    r = dict(K=K, nmax=nmax, e_iso=e_iso, e_var=res.fun, lam=lam,
             gap_iso=(e_iso - ref) * 1000, gap_var=(res.fun - ref) * 1000,
             condS=float(np.linalg.cond(Smat)),
             M1=float(np.abs(M).sum()), H1=float(np.abs(H).sum()),
             Hh1=float(np.abs(Hh).sum()))
    rows.append(r)
    print("K=%4d  gap_iso=%7.3f  gap_var=%7.3f mHa | lam*=%6.3f condS=%7.2f | "
          "||M||1=%9.2f ||Hhat||1=%9.2f" % (K, r["gap_iso"], r["gap_var"], lam, r["condS"], r["M1"], r["Hh1"]))
    sys.stdout.flush()
json.dump(rows, open("debug/data/p60_freescale_l%d.json" % lmax, "w"), indent=1)
K = np.array([r["K"] for r in rows], float)
for key in ("M1", "Hh1", "condS"):
    y = np.array([r[key] for r in rows], float)
    p = np.polyfit(np.log(K), np.log(y), 1)[0]
    print("  scaling  %-6s ~ K^%.3f" % (key, p))
