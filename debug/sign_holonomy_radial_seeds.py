"""Addendum to sign_holonomy_probe.py: radial-only (Gaunt-preserving) rotations,
5 seeds per n_max, physical He singlet.  Is the native hydrogenic RADIAL basis
special for the stoquastic gap, or only the angular (Gaunt) structure?"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sign_holonomy_probe as P

for nmax in (2, 3, 4):
    h, V, orbs = P.he_tensors(nmax)
    ns = len(orbs)
    mv = np.array([o[2] for o in orbs])
    keep = lambda a, b: mv[a] + mv[b] == 0
    H0, _ = P.singlet_H(h, V, keep=keep)
    d0 = P.graph_stats(H0, n_rand=0)["delta"]
    ds = []
    for sd in range(5):
        rng = np.random.default_rng(2000 + sd)
        U = np.zeros((ns, ns))
        groups = {}
        for i, (n, l, m) in enumerate(orbs):
            groups.setdefault((l, m), []).append(i)
        ru = {}
        for (l, m), ids in groups.items():
            if l not in ru:
                ru[l] = P.random_orth(len(ids), rng)
            U[np.ix_(ids, ids)] = ru[l]
        hr, Vr = P.rotate(h, V, U)
        Hr, _ = P.singlet_H(hr, Vr, keep=keep)
        ds.append(P.graph_stats(Hr, n_rand=0)["delta"])
    ds = np.array(ds)
    print(f"n_max={nmax}: M_L=0 Delta native={d0:.4e}; radial-rotated seeds "
          f"{np.round(ds, 4).tolist()} ratio min/median/max = "
          f"{ds.min()/d0:.2f}/{np.median(ds)/d0:.2f}/{ds.max()/d0:.2f}", flush=True)
