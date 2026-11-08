"""Does xTC (the cheap native non-Hermitian TC) buy per-atom accuracy, and does it
SCALE with electron count?  The v5.0.2-5.0.8 arc answered this for He (2e) and Li (3e);
this extends the SAME tracked engine to Be (4e) -- the untested point -- to see whether
the cheap-TC win holds, degrades, or the gamma-fragility worsens as electrons pile up.

For each atom, at a matched (ns, k) s-only Coulomb-Sturmian basis:
  E_plain      = plain s-only Coulomb FCI (gamma-independent)
  E_TC(gamma)  = native non-Hermitian TC (D + K [+ xTC-L3]) ground, over a gamma sweep
We report the plain error, the BEST TC error + its gamma, and whether E_TC crosses exact
over the sweep (monotone-through = gamma-fragile: the good gamma is a non-variational
crossing, not an internal optimum -- the v5.0.7 STOP).

Non-relativistic exact energies (Ha): He -2.903724, Li -7.478060, Be -14.667356.
Everything is the CHEAP route (keeps sparsity, O(1) conditioning); the ACCURATE route
(R12-CI, 0.45-0.80 mHa on He) is a separate, conditioning-limited engine (He-only).
"""
from __future__ import annotations

import json
import numpy as np
from scipy.optimize import minimize_scalar

from geovac import transcorrelated_sturmian as TC

EXACT = {"He": -2.903724, "Li": -7.478060, "Be": -14.667356}
# (label, Z, n_elec, with_L3, k_seed)   He/Li k from the v5.0.x arc; Be optimized here.
ATOMS = [
    ("He", 2, 2, False, 1.7),
    ("Li", 3, 3, True, 1.5),
    ("Be", 4, 4, True, 1.35),
]
NS, NG, NX = 3, 500, 96
GAMMAS = [0.4, 0.6, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6]


def plain_E(Z, ne, k):
    sys = TC.build_atomic_system(NS, k, 1.0, Z=Z, n_elec=ne, Ng=NG, nx=NX,
                                 with_L3=(ne >= 3))
    return TC.plain_energy(sys)


def opt_k(Z, ne, k0):
    r = minimize_scalar(lambda k: plain_E(Z, ne, k), bounds=(max(0.8, k0 - 0.6), k0 + 0.6),
                        method="bounded", options={"xatol": 1e-3})
    return float(r.x), float(r.fun)


def tc_E(Z, ne, k, g, with_L3):
    sys = TC.build_atomic_system(NS, k, g, Z=Z, n_elec=ne, Ng=NG, nx=NX, with_L3=with_L3)
    op = sys.xtc_full() if sys.v2 is not None else sys.tc2()   # 2e has no genuine 3-body
    H = TC.build_fci_matrix(sys, op)
    E, im = TC.ground(H, hermitian=False)
    return float(E), float(im)


def main():
    out = []
    print(f"{'atom':>4} {'k*':>5} {'E_plain':>11} {'dplain':>8} | "
          f"{'bestTC':>8} {'@g':>4} {'cross?':>6} {'improve':>8}")
    print("-" * 64)
    for label, Z, ne, with_L3, k0 in ATOMS:
        k, E_plain = opt_k(Z, ne, k0)
        ex = EXACT[label]
        dplain = (E_plain - ex) * 1e3
        rows = []
        for g in GAMMAS:
            E, im = tc_E(Z, ne, k, g, with_L3)
            rows.append((g, E, (E - ex) * 1e3, im))
        dTC = np.array([r[2] for r in rows])
        i_best = int(np.argmin(np.abs(dTC)))
        g_best, d_best = rows[i_best][0], dTC[i_best]
        crosses = bool(dTC.min() < 0 < dTC.max())   # E_TC passes through exact -> fragile
        improve = abs(dplain) - abs(d_best)          # mHa of error removed at oracle gamma
        print(f"{label:>4} {k:>5.2f} {E_plain:>11.5f} {dplain:>+8.2f} | "
              f"{d_best:>+8.2f} {g_best:>4.2f} {str(crosses):>6} {improve:>+8.2f}")
        out.append(dict(atom=label, Z=Z, n_elec=ne, k=k, E_plain=E_plain, exact=ex,
                        dplain_mHa=dplain, best_TC_mHa=float(d_best), best_gamma=g_best,
                        crosses_exact=crosses, oracle_improve_mHa=float(improve),
                        rows=[dict(gamma=g, E_TC=E, dTC_mHa=d, imag=im) for g, E, d, im in rows]))
    with open("debug/data/xtc_accuracy_per_atom.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nlegend: dplain/bestTC in mHa vs non-rel exact; cross?=E_TC passes through exact")
    print("        over the gamma sweep (True => the accurate gamma is a non-variational")
    print("        crossing, not an internal optimum = the v5.0.7 gamma-fragility STOP).")
    print("wrote debug/data/xtc_accuracy_per_atom.json")


if __name__ == "__main__":
    main()
