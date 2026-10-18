"""QG secondary placement (t on density i): route A (momentum triangle, closed-form loop) against the
double-bubble form phi1_gem_z, at three momenta and gamma = 1, 3/2.
Route A: Phi~_1(p) = (2pi)^-3 (2pi/p) int int_{|q-w|<=p<=q+w} q w rho~(q) rho~(w) T(p,q,w)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

OUT = Path(__file__).parent / "data" / "ring_ext_QG_2b_routes_i.json"


def phi1_routeA(p, g):
    p = mp.mpf(p); g = mp.mpf(g)
    rt = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    pref = (2 * mp.pi) ** -3 * (8 * mp.pi * g) ** 2 * 4 * mp.pi

    def inner(q):
        f = lambda w: w * rt(w) * Cm.I221(g, g, 0, w * w, p * p, q * q)
        return mp.quad(f, [abs(p - q), p + q])
    val = mp.quad(lambda q: q * rt(q) * inner(q), [0, p, 2 * p, 4 * p + 2, mp.inf])
    return (2 * mp.pi) ** -3 * 2 * mp.pi / p * pref * val


def main():
    rec = {}
    mp.mp.dps = 15
    ctx.dps = 30
    for g in ("1", "1.5"):
        for p in ("0.5", "1.3", "3.0"):
            a = phi1_routeA(p, g)
            b = Cm.phi1_gem_z(acb(arb(p) ** 2), arb(g), 120, 120).real
            rel = abs(a / Cm.to_mp(b) - 1)
            print(f"gamma {g} p {p}: route A {mp.nstr(a, 15)}  bubble {b.mid().str(16, radius=False)}  rel {mp.nstr(rel, 3)}", flush=True)
            rec[f"g{g}_p{p}"] = {"A": mp.nstr(a, 16), "B": b.mid().str(25, radius=False), "rel": float(rel)}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
