"""End-to-end ALL-QUADRATURE pipeline, as an independent check on the energy.

The closed-form pipeline (`qfd_h2.py` / `qfd_lih.py`) is validated integral by
integral in those drivers.  This module does the complementary thing: it builds
S, h and g ENTIRELY from `qfd_quad` numerics, runs the identical Loewdin + FCI
machinery, and compares the resulting total energy with the closed-form one.  A
per-integral agreement can in principle hide an assembly/indexing error; an
energy agreement cannot.

Run:  python -u debug/qfd_quad_pipeline.py [h2|lih|both]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from mpmath import mp

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import qfd_assemble as AS  # noqa: E402
import qfd_core as Q  # noqa: E402
import qfd_quad as V  # noqa: E402

OUT = REPO / "debug" / "data"


def quad_S_h(orbs, ZA, ZB, R):
    n = len(orbs)
    Sm, hm = mp.matrix(n, n), mp.matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            s = V.one_electron_quad(orbs[i], orbs[j], "S", R)
            t = V.one_electron_quad(orbs[i], orbs[j], "T", R)
            va = V.one_electron_quad(orbs[i], orbs[j], "inv_ra", R)
            vb = V.one_electron_quad(orbs[i], orbs[j], "inv_rb", R)
            Sm[i, j] = Sm[j, i] = s
            hm[i, j] = hm[j, i] = t - ZA * va - ZB * vb
    return Sm, hm


def quad_g(orbs, R, tau_max):
    n = len(orbs)
    canon: dict = {}
    gm: dict = {}
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    key = AS._canon_key(p, q, r, s)
                    if key not in canon:
                        oa, ob, oc, od = (orbs[i] for i in key)
                        cls, _ = AS.classify((oa, ob, oc, od))
                        if cls == "one-center":
                            v = V.one_center_eri_quad(oa, ob, oc, od)
                        elif cls == "aabb":
                            v = V.aabb_quad(oa, ob, oc, od, R)
                        elif cls == "hybrid":
                            v = (V.hybrid_quad(oa, ob, oc, od, R)
                                 if oa[0] == ob[0]
                                 else V.hybrid_quad(oc, od, oa, ob, R))
                        else:
                            a, b, c, d = oa, ob, oc, od
                            if a[0] != "A":
                                a, b = b, a
                            if c[0] != "A":
                                c, d = d, c
                            v, _ = V.exchange_numeric_neumann(
                                a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                                (c[2], 0, 0), (d[2], 0, 0), R,
                                tau_max=tau_max(key) if callable(tau_max)
                                else tau_max)
                        canon[key] = v
                        print(f"    quad ({key[0]}{key[1]}|{key[2]}{key[3]}) "
                              f"{cls:11s} {mp.nstr(v, 18)}")
                    gm[(p, q, r, s)] = canon[key]
    return gm


def energy_from(Sm, hm, gm, n, n_elec, ZA, ZB, R):
    X = AS.lowdin(Sm, n)
    ht, gt = AS.transform(X, hm, gm, n)
    e = AS.fci_ground(ht, gt, n, n_elec)
    return e + mp.mpf(str(sp.N(sp.nsimplify(ZA) * sp.nsimplify(ZB)
                               / sp.nsimplify(R), 40)))


def run(which):
    res = {}
    if which in ("h2", "both"):
        print("H2 -- all-quadrature pipeline (dps 20)")
        R = sp.Rational(7, 5)
        orbs = [("A", Fraction(1), 1), ("B", Fraction(1), 1)]
        t = time.time()
        mp.dps = 20
        Sm, hm = quad_S_h(orbs, 1, 1, R)
        gm = quad_g(orbs, R, tau_max=2)      # H2 tau series terminates at 2
        e_q = energy_from(Sm, hm, gm, 2, 2, 1, 1, R)
        mp.dps = 60
        S, h = AS.build_S_h(orbs, 1, 1, R)
        gmap, _c = AS.build_g(orbs, R, tau_max=6)
        e_c, _ee, _v = AS.total_energy(S, h, gmap, orbs, 2, 1, 1, R, 40)
        d = abs(mp.mpf(e_q) - mp.mpf(e_c))
        print(f"  E(all-quadrature) = {mp.nstr(mp.mpf(e_q), 20)}")
        print(f"  E(closed form)    = {mp.nstr(mp.mpf(e_c), 20)}")
        print(f"  |difference|      = {mp.nstr(d, 3)}   "
              f"({time.time()-t:.0f}s)")
        res["h2"] = {"E_quadrature": mp.nstr(mp.mpf(e_q), 20),
                     "E_closed_form": mp.nstr(mp.mpf(e_c), 20),
                     "difference": mp.nstr(d, 3)}
    if which in ("lih", "both"):
        print("\nLiH -- all-quadrature pipeline (dps 20, exchange tau <= 10)")
        R = sp.Rational(603, 200)
        orbs = [("A", Fraction(3), 1), ("A", Fraction(3), 2),
                ("B", Fraction(1), 1)]
        t = time.time()
        mp.dps = 20
        Sm, hm = quad_S_h(orbs, 3, 1, R)
        gm = quad_g(orbs, R, tau_max=10)
        e_q = energy_from(Sm, hm, gm, 3, 4, 3, 1, R)
        mp.dps = 60
        S, h = AS.build_S_h(orbs, 3, 1, R)
        gmap, _c = AS.build_g(orbs, R, tau_max=10, exchange_hp_dps=25)
        e_c, _ee, _v = AS.total_energy(S, h, gmap, orbs, 4, 3, 1, R, 25)
        d = abs(mp.mpf(e_q) - mp.mpf(e_c))
        print(f"  E(all-quadrature) = {mp.nstr(mp.mpf(e_q), 20)}")
        print(f"  E(closed form, same tau_max = 10) = "
              f"{mp.nstr(mp.mpf(e_c), 20)}")
        print(f"  |difference|      = {mp.nstr(d, 3)}   "
              f"({time.time()-t:.0f}s)")
        res["lih"] = {"E_quadrature": mp.nstr(mp.mpf(e_q), 20),
                      "E_closed_form_tau10": mp.nstr(mp.mpf(e_c), 20),
                      "difference": mp.nstr(d, 3)}
    (OUT / "qfd_quad_pipeline.json").write_text(json.dumps(res, indent=2),
                                                encoding="utf-8")
    print(f"\nwrote {OUT / 'qfd_quad_pipeline.json'}")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "both")
