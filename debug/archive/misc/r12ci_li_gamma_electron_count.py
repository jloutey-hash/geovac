"""THE EXPERIMENT: does the variational geminal width survive ELECTRON COUNT?

The isoelectronic series (debug/r12ci_gamma_isoelectronic.py) showed gamma_opt ∝ Z to
within one grid step -- but all four systems were 1s^2, TWO electrons.  Electron count is
exactly where the cheap non-Hermitian route died (Be 33x worse than Li, v5.0.12).

So: run Li+ (2 electrons, Z=3) and Li (3 electrons, Z=3) in the SAME engine with the SAME
conventions, at the same nuclear charge, and ask whether gamma_opt moves.

  gamma_opt unchanged  =>  the tabulation unit is the PAIR, spectator-independent
                           -> a determined law survives into many-electron atoms
  gamma_opt shifts     =>  screening changes the pair scale -> tabulate against Z_eff,
                           and the isoelectronic law is 2-electron-only

Engine: debug/r12ci_ne_engine.py, validated by debug/r12ci_ne_gates.py (gamma=0 exact
reduction at 1e-13; N=2 vs an independent engine at 1e-15; references converge onto it).
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import r12ci_ne_engine as E  # noqa: E402

EXACT = {"Li+": -7.2799134126, "Li": -7.4780603236}
NG, LMAX, NX = 320, 10, 48
KS = [1.1, 1.4, 1.7, 2.0, 2.4, 2.85]
GAMMAS = [0.04, 0.08, 0.13, 0.20, 0.30, 0.45, 0.65, 0.95]

SYSTEMS = [("Li+", 2, 0, 3.0, 4), ("Li", 3, 1, 3.0, 4)]

out = {"grid": dict(Ng=NG, Lmax=LMAX, nx=NX), "ks": KS, "gammas": GAMMAS,
       "exact": EXACT, "systems": {}}

for name, N, ms2, Z, ns in SYSTEMS:
    t0 = time.time()
    print("=" * 84)
    print(f"{name}:  N={N} electrons, Z={Z}, ns={ns}, ms2={ms2}")
    print("=" * 84)
    rows = []
    # plain (no geminal) baseline at each k
    base = {}
    for k in KS:
        S0, H0, _ = E.build(ns, N, Z, k, 0.5, ms2, Ng=NG, with_geminal=False,
                            Lmax=LMAX, nx=NX)
        base[k], _ = E.solve(S0, H0)
    print("  plain (no geminal):  " +
          "  ".join(f"k={k}: {base[k]:.6f}" for k in KS))
    print()
    hdr = f"{'gamma':>7}" + "".join(f"{'k=' + str(k):>14}" for k in KS)
    print(hdr)
    best = None
    for gam in GAMMAS:
        line = f"{gam:>7.2f}"
        for k in KS:
            S1, H1, d1 = E.build(ns, N, Z, k, gam, ms2, Ng=NG, with_geminal=True,
                                 Lmax=LMAX, nx=NX)
            e, _ = E.solve(S1, H1, nd=len(d1))
            line += f"{e:>14.6f}"
            rows.append(dict(gamma=gam, k=k, E=e, gain_mHa=(base[k] - e) * 1000))
            if best is None or e < best[0]:
                best = (e, k, gam)
        print(line)
    E_best, k_best, g_best = best
    print()
    print(f"  optimum: k={k_best}  gamma={g_best}  E={E_best:.6f}  "
          f"gamma/Z={g_best / Z:.3f}  geminal gain={1000 * (base[k_best] - E_best):.2f} mHa")
    print(f"  ({time.time() - t0:.0f} s)")
    print()
    out["systems"][name] = dict(N=N, Z=Z, ns=ns, base=base, rows=rows,
                                E_best=E_best, k_best=k_best, gamma_best=g_best,
                                gamma_over_Z=g_best / Z)

a, b = out["systems"]["Li+"], out["systems"]["Li"]
print("=" * 84)
print("VERDICT")
print("=" * 84)
print(f"  Li+  (2e, Z=3):  gamma_opt = {a['gamma_best']:.2f}   (gamma/Z = {a['gamma_over_Z']:.3f})")
print(f"  Li   (3e, Z=3):  gamma_opt = {b['gamma_best']:.2f}   (gamma/Z = {b['gamma_over_Z']:.3f})")
ratio = b["gamma_best"] / a["gamma_best"]
print(f"  shift with electron count: {ratio:.2f}x")
print(f"  isoelectronic law predicted gamma/Z ~ 0.18-0.23 (Z>=2)")
# EDGE CHECK -- an optimum at a grid boundary is not a determination
for nm, d in out["systems"].items():
    edge = []
    if d["gamma_best"] in (GAMMAS[0], GAMMAS[-1]):
        edge.append("gamma")
    if d["k_best"] in (KS[0], KS[-1]):
        edge.append("k")
    d["at_grid_edge"] = edge
    flag = ("AT GRID EDGE in " + "+".join(edge) + " -- NOT a determination") if edge else "interior"
    print("  %4s: %s" % (nm, flag))
out["verdict"] = dict(gamma_Liplus=a["gamma_best"], gamma_Li=b["gamma_best"], ratio=ratio,
                      edge_Liplus=a.get("at_grid_edge"), edge_Li=b.get("at_grid_edge"))

os.makedirs("debug/data", exist_ok=True)
with open("debug/data/r12ci_li_gamma_electron_count.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_li_gamma_electron_count.json")
