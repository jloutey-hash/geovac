"""Transfer test, done the way the flat basin demands.

The (k, gamma) scan showed gamma_opt is ill-determined -- the basin is flat (0.9-2.2 mHa
across a 24x range in gamma) and the sign of dE/dgamma flips with k.  Chasing gamma_opt is
therefore the wrong question.  The question F12 actually answers, and the one that matters
for "solve-and-tabulate", is:

    does ONE fixed gamma give near-optimal energy for BOTH a 2-electron and a 3-electron
    system at the same Z?

Reported as a PENALTY vs each system's own per-gamma optimum, at each system's best k.
"""
import io, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
os.chdir(ROOT); sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import r12ci_ne_engine as E

NG, LMAX, NX = 320, 10, 48
GAMMAS = [0.10, 0.20, 0.35, 0.50, 0.65, 0.85, 1.10, 1.45, 1.90, 2.50, 3.30]
CASES = [("Li+", 2, 0, 3.0, 4, 2.4), ("Li", 3, 1, 3.0, 4, 1.7)]

res = {}
for name, N, ms2, Z, ns, k in CASES:
    t0 = time.time()
    S0, H0, _ = E.build(ns, N, Z, k, 0.5, ms2, Ng=NG, with_geminal=False, Lmax=LMAX, nx=NX)
    e0, _ = E.solve(S0, H0)
    es = []
    for g in GAMMAS:
        S1, H1, d1 = E.build(ns, N, Z, k, g, ms2, Ng=NG, with_geminal=True, Lmax=LMAX, nx=NX)
        e, _ = E.solve(S1, H1, nd=len(d1))
        es.append(e)
    es = np.array(es)
    res[name] = dict(k=k, Z=Z, N=N, plain=e0, E=es.tolist(),
                     best=float(es.min()), gbest=GAMMAS[int(es.argmin())])
    print(f"{name}: k={k}  plain={e0:.6f}  best={es.min():.6f} at gamma={GAMMAS[int(es.argmin())]}"
          f"  ({time.time()-t0:.0f}s)")

print()
print("=" * 76)
print("PENALTY (mHa) of a single fixed gamma vs each system's own optimum")
print("=" * 76)
print(f"{'gamma':>7}{'Li+ (2e)':>12}{'Li (3e)':>12}{'max':>10}{'gam/Z':>9}")
best_c = None
for i, g in enumerate(GAMMAS):
    pa = (res["Li+"]["E"][i] - res["Li+"]["best"]) * 1000
    pb = (res["Li"]["E"][i] - res["Li"]["best"]) * 1000
    m = max(pa, pb)
    if best_c is None or m < best_c[0]:
        best_c = (m, g)
    print(f"{g:>7.2f}{pa:>12.3f}{pb:>12.3f}{m:>10.3f}{g/3.0:>9.3f}")
print()
print(f"best single gamma = {best_c[1]:.2f}  (gamma/Z = {best_c[1]/3.0:.3f}), "
      f"worst-case penalty {best_c[0]:.3f} mHa")
edge = []
for nm in res:
    if res[nm]["gbest"] in (GAMMAS[0], GAMMAS[-1]):
        edge.append(nm)
print("grid-edge optima:", edge if edge else "none -- both interior")
res["transfer"] = dict(best_gamma=best_c[1], worst_penalty_mHa=best_c[0], edge=edge)
with open("debug/data/r12ci_li_transfer_test.json", "w") as f:
    json.dump({"gammas": GAMMAS, **res}, f, indent=2)
print("\nwrote debug/data/r12ci_li_transfer_test.json")
