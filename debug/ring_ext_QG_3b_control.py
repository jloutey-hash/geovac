"""QG control (registered): Yukawa kernels exp(-gamma r)/r on all three edges; gamma -> 0 must
reproduce the Coulomb triangle G(t) of PT2 at two t.  Same double-bubble pipeline as Y:
  Y_Yuk(t; gamma) = (1/(2 pi^2)) int w^2 t^3 (w^2+t^2)^(-3/2) Phi~_Y(w) dw.
gamma = 0 member computed exactly (v = (1 - cos phi)/2 removes the 1/N endpoint behaviour), and a
sequence gamma = 0.1, 0.03, 0.01 to show the approach.  Reference G: density_loop_B.json (25 digits)."""
from __future__ import annotations
import json, sys
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

DATA = Path(__file__).parent / "data"
G = {"0.7": "0.05376960760532172894474858", "1.6": "0.1199894528338774458832898"}


def wnodes(digits):
    mp.mp.dps = digits + 10
    h = mp.mpf(1) / 16
    out = []
    for k in range(-110, 111):
        tau = k * h
        w = mp.exp(mp.pi / 2 * mp.sinh(tau))
        dw = h * mp.pi / 2 * mp.cosh(tau) * w
        if w ** 2 * dw > mp.mpf(10) ** (-digits - 6) and dw / w ** 4 > mp.mpf(10) ** (-digits - 6):
            out.append((w, dw))
    return out


def job(args):
    kind, w, g, digits = args
    ctx.dps = digits + 10
    if kind == "coul":
        n = 60 + int(4 * float(w))
        a = Cm.phi_coul_bubble(arb(w), min(n, 400), min(n, 400))
        b = Cm.phi_coul_bubble(arb(w), min(2 * n, 800), min(2 * n, 800))
        return b.mid().str(digits + 5, radius=False), float(abs(a - b) / abs(b))
    h = arb(1) / 32
    a = Cm.phi_yuk3_ts(arb(w), g, h, digits + 6).real
    b = Cm.phi_yuk3_ts(arb(w), g, h / 2, digits + 6).real
    return b.mid().str(digits + 5, radius=False), float(abs(a - b) / abs(b))


def Y_of(nodes, phis, t):
    t = mp.mpf(t)
    s = sum(w * w * mp.mpf(p) * dw * t ** 3 * (w * w + t * t) ** mp.mpf(-1.5) for (w, dw), p in zip(nodes, phis))
    return s / (2 * mp.pi ** 2)


def main():
    rec = {}
    with Pool(15) as pool:
        for kind, g, digits in [("coul", "0", 20), ("yuk", "0.1", 10), ("yuk", "0.03", 10), ("yuk", "0.01", 10)]:
            nodes = wnodes(digits)
            res = pool.map(job, [(kind, mp.nstr(w, digits + 8), g, digits) for w, _ in nodes], chunksize=1)
            mp.mp.dps = digits + 8
            phis = [r[0] for r in res]
            err = max(r[1] for r in res)
            ys = {t: Y_of(nodes, phis, t) for t in G}
            rec[f"{kind}_g{g}"] = {t: {"Y": mp.nstr(ys[t], digits), "G": G[t][:digits + 2],
                                       "Y/G-1": mp.nstr(ys[t] / mp.mpf(G[t]) - 1, 3)} for t in G}
            rec[f"{kind}_g{g}"]["max_node_selfcheck"] = err
            print(kind, g, rec[f"{kind}_g{g}"], flush=True)
    (DATA / "ring_ext_QG_3b_control.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
