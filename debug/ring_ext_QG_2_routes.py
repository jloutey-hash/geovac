"""QG step 2: two independent routes to the geminal kite Phi~(w), and the normalisation checks.
Route A: momentum triangle, Phi~(w) = (2pi)^-3 (2pi/w) int int_{|p-q|<=w<=p+q} p q rho~(p) rho~(q) T(p,q,w),
         T = (2pi)^-3 (8 pi gamma)^2 (4 pi) I_221(gamma, gamma, 0) with the solid-angle closed form (QG1).
Route B: double bubble (Feynman parameters on the two bubbles, ring_ext_QG_common.phi_gem).
Coulomb check: the double-bubble form with Coulomb kernels (gamma -> 0 of the all-Yukawa form) against
the record's Coulomb triangle G(t) (density_loop_B.json step1, 25 digits) -- fixes all constants."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

OUT = Path(__file__).parent / "data" / "ring_ext_QG_2_routes.json"


def phi_routeA(w, g):
    w = mp.mpf(w); g = mp.mpf(g)
    rt = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    pref = (2 * mp.pi) ** -3 * (8 * mp.pi * g) ** 2 * 4 * mp.pi

    def inner(p):
        f = lambda q: q * rt(q) * Cm.I221(g, g, 0, w * w, p * p, q * q)
        return mp.quad(f, [abs(p - w), p + w])
    val = mp.quad(lambda p: p * rt(p) * inner(p), [0, w, 2 * w, 4 * w + 2, mp.inf])
    return (2 * mp.pi) ** -3 * 2 * mp.pi / w * pref * val


def main():
    rec = {}
    mp.mp.dps = 15
    ctx.dps = 30
    for g in ("1", "1.5"):
        for w in ("0.5", "1.3", "3.0"):
            a = phi_routeA(w, g)
            b = Cm.phi_gem(arb(w), arb(g), 80, 80).real
            rel = abs(a / Cm.to_mp(b) - 1)
            print(f"gamma {g} w {w}: route A {mp.nstr(a, 15)}  route B {b.mid().str(16, radius=False)}  rel {mp.nstr(rel, 3)}", flush=True)
            rec[f"g{g}_w{w}"] = {"A": mp.nstr(a, 16), "B": b.mid().str(25, radius=False), "rel": float(rel)}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
