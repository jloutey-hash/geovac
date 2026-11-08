"""Probe B -- independent validation of the momentum-reconstructed ERIs and rho~.

Three independent ground truths:
  (1) one-centre (1s1s|1s1s) = 5a/8 exactly;
  (2) the corpus's EXACT closed forms geovac.two_center_eri.{aabb_value,
      hybrid_closed_form, exchange_value} for the (AA|BB) / (AA|AB) / (AB|AB)
      two-centre classes (Paper 58 engine);
  (3) the direct prolate-spheroidal FT for rho~ itself at large k (the regime
      where the Feynman t-quadrature phase is hardest).
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from debug.probeB_thc_density import (feynman_grid, rho_tilde,  # noqa: E402
                                      rho_tilde_direct, system)
from debug.probeB_thc_lambda import build_grid, node_count  # noqa: E402


def eri_raw(orbs, bands, tgrid):
    n = len(orbs)
    om_all, C_all, S_all = [], [], []
    for (kb, wkb, mub, wmub) in bands:
        raw = np.zeros((n, n, kb.size, mub.size), dtype=complex)
        for i in range(n):
            for j in range(i, n):
                v = rho_tilde(orbs[i], orbs[j], kb, mub, tgrid)
                raw[i, j] = v
                raw[j, i] = v
        rot = np.transpose(raw, (2, 3, 0, 1)).reshape(-1, n, n)
        om = ((2.0 / math.pi) * wkb[:, None] * wmub[None, :]).reshape(-1)
        om_all.append(om)
        C_all.append(rot.real.copy())
        S_all.append(rot.imag.copy())
    om = np.concatenate(om_all)
    C = np.concatenate(C_all, 0)
    S = np.concatenate(S_all, 0)
    return (np.einsum("n,npq,nrs->pqrs", om, C, C, optimize=True)
            + np.einsum("n,npq,nrs->pqrs", om, S, S, optimize=True))


if __name__ == "__main__":
    tg = feynman_grid(y_max=26.0, panel=0.30, n_g=12)
    print(f"Feynman t-grid: {tg[0].size} nodes\n")

    # ---------- (3) rho~ at large k, two-centre pairs (hardest phase regime)
    print("--- rho~ vs direct prolate FT at LARGE k (two-centre pairs)")
    for sysname, pair in (("H2", (0, 1)), ("LiH", (0, 2))):
        orbs, _ = system(sysname)
        p, q = orbs[pair[0]], orbs[pair[1]]
        for kv, mv in ((8.0, 0.6), (18.0, 1.0), (25.0, 0.35)):
            fast = rho_tilde(p, q, np.array([kv]), np.array([mv]), tg)[0, 0]
            ref = rho_tilde_direct(p, q, kv, mv, n_g=48, n_pan=48, n_eta=420)
            print(f"    {sysname:6} {pair} k={kv:>5.0f} mu={mv:+.2f}  "
                  f"|rho~|={abs(ref):.3e}  abs err {abs(fast-ref):.2e}  "
                  f"rel err {abs(fast-ref)/max(abs(ref),1e-300):.2e}")

    # ---------- (1)+(2) ERI tensor vs exact closed forms
    print("\n--- reconstructed ERIs vs EXACT references")
    from geovac.two_center_eri import aabb_value, exchange_value, hybrid_closed_form
    import sympy as sp

    CASES = {
        "H2": dict(R=1.4, ZA=Fraction(1), ZB=Fraction(1),
                   oA=[(1, 0, 0)], oB=[(1, 0, 0)], iA=[0], iB=[1]),
        "LiH": dict(R=3.015, ZA=Fraction(3), ZB=Fraction(1),
                    oA=[(1, 0, 0), (2, 0, 0)], oB=[(1, 0, 0)], iA=[0, 1], iB=[2]),
    }
    for sysname, spec in CASES.items():
        orbs, _ = system(sysname)
        dP = max(abs(o1.z - o2.z) for o1 in orbs for o2 in orbs)
        for (Kmax, npan, ng, mpr) in ((40., 44, 10, 4), (80., 80, 10, 8), (200., 170, 12, 12)):
            t0 = time.time()
            bands = build_grid(Kmax, npan, ng, dP, pad=mpr)
            e = eri_raw(orbs, bands, tg)
            print(f"\n  {sysname}  Kmax={Kmax:g}  M={node_count(bands)}  "
                  f"({time.time()-t0:.1f}s)")
            R = spec["R"]
            ZA, ZB = spec["ZA"], spec["ZB"]
            # (AA|BB)
            for ia, oa in zip(spec["iA"], spec["oA"]):
                for ib, ob in zip(spec["iB"], spec["oB"]):
                    ex = aabb_value(ZA, oa, oa, ZB, ob, ob, R, prec=25)
                    got = e[ia, ia, ib, ib]
                    print(f"     (AA|BB) [{ia}{ia}|{ib}{ib}]  exact {ex:.12f}  "
                          f"got {got:.12f}  d={abs(got-ex):.2e}")
            # (AA|AB) hybrid  -- a,b,c on A, d on B
            oa = spec["oA"][0]
            ob = spec["oB"][0]
            ia, ib = spec["iA"][0], spec["iB"][0]
            ex = float(sp.re(sp.N(hybrid_closed_form(ZA, oa, oa, oa, ZB, ob, R), 25)))
            got = e[ia, ia, ia, ib]
            print(f"     (AA|AB) [{ia}{ia}|{ia}{ib}]  exact {ex:.12f}  "
                  f"got {got:.12f}  d={abs(got-ex):.2e}")
            # (AB|AB) exchange
            ex = exchange_value(ZA, oa, ob, ZB, oa, ob, R, tau_max=14)
            got = e[ia, ib, ia, ib]
            print(f"     (AB|AB) [{ia}{ib}|{ia}{ib}]  exact {ex:.12f}  "
                  f"got {got:.12f}  d={abs(got-ex):.2e}")
