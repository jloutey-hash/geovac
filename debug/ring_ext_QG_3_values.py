"""QG step 3: the geminal triangle Y(t) and its theta_t-derivatives at the fingerprint fibre, by the
double-bubble form on an exp-sinh w-grid (Phi~(w) cached), and the gamma -> 0 control.

Y(t) = (1/(2 pi^2)) int_0^inf w^2 Phi~(w) k(t,w) dw,  k = t^3 (w^2+t^2)^(-3/2) = (1 + x)^(-3/2), x = w^2/t^2.
theta_t^j k = (-2 x d/dx)^j (1+x)^(-3/2).

usage: python ring_ext_QG_3_values.py grid <kind> <gamma> <digits> [h_tau]   (kind = gem | yuk3)
       python ring_ext_QG_3_values.py control
"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
import sympy as sp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

DATA = Path(__file__).parent / "data"


def phi_adaptive(args):
    kind, w, g, digits = args
    ctx.dps = digits + 12
    w = arb(w)
    f = Cm.phi_gem_ts if kind == "gem" else Cm.phi_yuk3_ts
    h = arb(1) / 16
    prev = f(w, g, h, digits + 10).real
    while True:
        h = h / 2
        cur = f(w, g, h, digits + 10).real
        if abs(cur - prev) < abs(cur) * arb(10) ** (-digits) or h < arb(1) / 1024:
            return cur.mid().str(digits + 8, radius=False), str(h)
        prev = cur


def wgrid(htau, digits):
    """exp-sinh nodes on (0, inf): w = exp(pi/2 sinh tau), weight dw."""
    mp.mp.dps = digits + 10
    out = []
    k = -int(6 / htau)
    while k <= int(6 / htau):
        tau = k * htau
        w = mp.exp(mp.pi / 2 * mp.sinh(tau))
        dw = htau * mp.pi / 2 * mp.cosh(tau) * w
        out.append((w, dw))
        k += 1
    return out


def theta_kernels(J):
    x = sp.Symbol("x")
    out, cur = [], (1 + x) ** sp.Rational(-3, 2)
    for j in range(J + 1):
        out.append(sp.lambdify(x, sp.simplify(cur), "mpmath"))
        cur = sp.simplify(-2 * x * sp.diff(cur, x))
    return out


def build_grid(kind, g, digits, htau):
    """Phi~ on the w-grid, cached in data/ring_ext_QG_3_grid_<kind>_g<g>_d<digits>.json."""
    fn = DATA / f"ring_ext_QG_3_grid_{kind}_g{g}_d{digits}_h{htau}.json"
    if fn.exists():
        return json.loads(fn.read_text())
    grid = wgrid(mp.mpf(htau), digits)
    # integrand cut: w^2 Phi~ ~ w^3 at small w, ~ w^-6 at large w (k ~ w^-3 is not included)
    keep = [(w, dw) for w, dw in grid if w ** 3 * dw > mp.mpf(10) ** (-digits - 4) and dw / w ** 5 > mp.mpf(10) ** (-digits - 4)]
    tasks = [(kind, mp.nstr(w, digits + 8), g, digits) for w, dw in keep]
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(phi_adaptive, tasks, chunksize=1)
    rec = {"kind": kind, "gamma": g, "digits": digits, "htau": htau,
           "w": [mp.nstr(w, digits + 8) for w, _ in keep], "dw": [mp.nstr(dw, digits + 8) for _, dw in keep],
           "phi": [r[0] for r in res], "h_used": [r[1] for r in res], "seconds": time.time() - t0}
    fn.write_text(json.dumps(rec, indent=0))
    return rec


def Y_theta(rec, t, J):
    mp.mp.dps = rec["digits"] + 8
    t = mp.mpf(t)
    ks = theta_kernels(J)
    acc = [mp.mpf(0)] * (J + 1)
    for w, dw, ph in zip(rec["w"], rec["dw"], rec["phi"]):
        w, dw, ph = mp.mpf(w), mp.mpf(dw), mp.mpf(ph)
        x = w * w / (t * t)
        base = w * w * ph * dw
        for j in range(J + 1):
            acc[j] += base * ks[j](x)
    return [a / (2 * mp.pi ** 2) for a in acc]


def main():
    if sys.argv[1] == "grid":
        kind, g, digits = sys.argv[2], sys.argv[3], int(sys.argv[4])
        htau = sys.argv[5] if len(sys.argv) > 5 else "0.0625"
        rec = build_grid(kind, g, digits, htau)
        print(kind, g, "nodes", len(rec["w"]), "seconds", round(rec["seconds"], 1))
        for t in ("0.7", "1.6"):
            print("t", t, "Y =", mp.nstr(Y_theta(rec, t, 0)[0], digits))
    elif sys.argv[1] == "control":
        G = {"0.7": mp.mpf("0.05376960760532172894474858"), "1.6": mp.mpf("0.1199894528338774458832898")}
        out = {}
        for g in ("1e-3", "1e-4", "1e-5"):
            rec = build_grid("yuk3", g, 16, "0.0625")
            out[g] = {t: Y_theta(rec, t, 0)[0] for t in G}
        res = {}
        for t in G:
            y3, y4, y5 = (out[g][t] for g in ("1e-3", "1e-4", "1e-5"))
            rich = (10 * y5 - y4) / 9
            res[t] = {"Y_1e-3": mp.nstr(y3, 14), "Y_1e-4": mp.nstr(y4, 14), "Y_1e-5": mp.nstr(y5, 14),
                      "Richardson": mp.nstr(rich, 14), "G": mp.nstr(G[t], 14),
                      "rel_Y1e-5": float(abs(y5 / G[t] - 1)), "rel_rich": float(abs(rich / G[t] - 1)),
                      "slope_ratio": mp.nstr((y3 - y4) / (y4 - y5), 6)}
            print(t, res[t], flush=True)
        (DATA / "ring_ext_QG_3_control.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
