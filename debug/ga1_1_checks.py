"""GA1 stage 1: checks of the new loop and kite before any series (record: relative 1e-8 or better at
three random points).
 (a) TD2 algebra control: I_221 through TD2 equals ring_ext_QG_common.I221 (bi-dual) at three triples.
 (b) lin loop closed form (analytic 4th mixed derivative) against direct 3-D quadrature at three random
     triples, gamma = 3/2; also the exp (2,2,1) loop as the QG1 control.
 (c) phi_z(shape='exp') equals ring_ext_QG_common.phi_gem_z at three z (one complex), 50 digits.
 (d) lin kite: route A (momentum triangle, closed-form loop, 2-D quadrature) against route B (double
     bubble, Feynman parameters) at three random w, gamma = 3/2.
 (e) phi_ts('lin') against phi_z('lin') at one real w.
Gate: every relative difference < 1e-8 (the 50-digit identities far below)."""
from __future__ import annotations
import json, random, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm
import ga1_common as GA

OUT = Path(__file__).parent / "data" / "ga1_1_checks.json"


def direct(x1, x2, g, shape, Nang=160, Nr=240):
    """int d^3 l g(l-x1) g(l-x2) / |l|^2 (x3 = 0), float64 Gauss rules, independent of the closed form.
    shape exp: g(q) = (q^2+g^2)^-2 ; lin: g(q) = (3 g^2 - q^2)/(q^2+g^2)^3."""
    x1, x2 = np.array(x1), np.array(x2)
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, Nang, endpoint=False); wph = 2 * np.pi / Nang
    st = np.sqrt(1 - uu**2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    d1, d2 = n @ x1, n @ x2
    R = max(np.linalg.norm(x1), np.linalg.norm(x2)) + 2
    sr, wr = np.polynomial.legendre.leggauss(Nr)

    def gfun(q2):
        if shape == "exp":
            return (q2 + g * g) ** -2
        return (3 * g * g - q2) / (q2 + g * g) ** 3

    def f_at(ri):
        q1 = ri**2 - 2 * ri * d1 + x1 @ x1
        q2 = ri**2 - 2 * ri * d2 + x2 @ x2
        return gfun(q1) * gfun(q2)          # times ri^2 / ri^2 (the 1/|l|^2 cancels the Jacobian)
    tot = 0.0
    for lo, hi in [(0, R / 3), (R / 3, 2 * R / 3), (2 * R / 3, R)]:
        r = (hi - lo) / 2 * sr + (hi + lo) / 2; w = (hi - lo) / 2 * wr
        for ri, wi in zip(r, w):
            tot += wi * np.sum(wu[:, None] * f_at(ri)) * wph
    s = (sr + 1) / 2; ws = wr / 2
    for si, wi in zip(s, ws):
        ri = R / si
        tot += wi * R / si**2 * np.sum(wu[:, None] * f_at(ri)) * wph
    return tot


def main():
    mp.mp.dps = 30
    rec = {}
    g = mp.mpf(3) / 2
    rng = random.Random(20261022)
    worst = 0.0
    # (a) TD2 vs bi-dual on I_221
    rows = []
    for k in range(3):
        x1 = [rng.uniform(-2, 2) for _ in range(3)]; x2 = [rng.uniform(-2, 2) for _ in range(3)]
        s12 = sum((mp.mpf(a) - b) ** 2 for a, b in zip(x1, x2)); s13 = sum(mp.mpf(a) ** 2 for a in x1); s23 = sum(mp.mpf(b) ** 2 for b in x2)
        a = GA.loop_exp_TD2(g, g, 0, s12, s13, s23); b = Cm.I221(g, g, 0, s12, s13, s23)
        rel = float(abs(a / b - 1)); rows.append({"x1": x1, "x2": x2, "TD2": mp.nstr(a, 25), "bidual": mp.nstr(b, 25), "rel": rel})
        print(f"(a) I_221 TD2 vs bi-dual, triple {k}: rel {rel:.2e}", flush=True)
        worst = max(worst, rel)
    rec["a_I221_TD2_vs_bidual"] = rows
    # (b) loops vs direct 3-D
    for shape in ("exp", "lin"):
        rows = []
        for k in range(3):
            x1 = [rng.uniform(-2, 2) for _ in range(3)]; x2 = [rng.uniform(-2, 2) for _ in range(3)]
            s12 = sum((mp.mpf(a) - b) ** 2 for a, b in zip(x1, x2)); s13 = sum(mp.mpf(a) ** 2 for a in x1); s23 = sum(mp.mpf(b) ** 2 for b in x2)
            closed = Cm.I221(g, g, 0, s12, s13, s23) if shape == "exp" else GA.loop_lin(g, g, 0, s12, s13, s23)
            dd = direct(x1, x2, float(g), shape)
            rel = abs(dd / float(closed) - 1)
            rows.append({"x1": x1, "x2": x2, "direct": dd, "closed": float(closed), "rel": rel})
            print(f"(b) {shape} loop triple {k}: direct {dd:.15e} closed {float(closed):.15e} rel {rel:.2e}", flush=True)
            worst = max(worst, rel)
        rec[f"b_loop_{shape}_vs_direct3d"] = rows
    # (c) phi_z('exp') vs Cm.phi_gem_z
    ctx.dps = 50
    rows = []
    for z in (acb("0.7"), acb("-1.3"), acb("0.9", "1.4")):
        a = GA.phi_z(z, arb(3) / 2, 120, 120, "exp"); b = Cm.phi_gem_z(z, arb(3) / 2, 120, 120)
        rel = float(abs(a - b) / abs(b)); rows.append({"z": str(z), "rel": rel})
        print(f"(c) phi exp new vs QG at z = {z}: rel {rel:.2e}", flush=True)
        worst = max(worst, rel)
    rec["c_phi_exp_new_vs_QG"] = rows
    # (d) lin kite: route A vs route B
    rows = []
    t0 = time.time()
    mp.mp.dps = 15
    ctx.dps = 30
    for w in ("0.5", "1.3", "3.0"):
        a = GA.phi_routeA(w, "1.5", "lin")
        b = GA.phi_z(acb(w) * acb(w), arb(3) / 2, 80, 80, "lin").real
        rel = float(abs(a / Cm.to_mp(b) - 1))
        rows.append({"w": w, "routeA": mp.nstr(a, 16), "routeB": b.mid().str(25, radius=False), "rel": rel})
        print(f"(d) lin kite w = {w}: route A {mp.nstr(a, 15)}  route B {b.mid().str(16, radius=False)}  rel {rel:.2e}  ({time.time()-t0:.0f}s)", flush=True)
        worst = max(worst, rel)
    rec["d_lin_kite_routeA_vs_routeB"] = rows
    # exp kite route A vs B with the same code (control of route A's normalisation)
    a = GA.phi_routeA("1.3", "1.5", "exp"); b = Cm.phi_gem(arb("1.3"), arb(3) / 2, 80, 80).real
    rel = float(abs(a / Cm.to_mp(b) - 1)); rec["d_exp_kite_control"] = {"w": "1.3", "rel": rel}
    print(f"(d) exp kite control w = 1.3: rel {rel:.2e}", flush=True)
    worst = max(worst, rel)
    # (e) tanh-sinh form
    ctx.dps = 40
    a = GA.phi_ts(arb("1.7"), arb(3) / 2, arb(1) / 32, 30, "lin").real
    b = GA.phi_z(acb("1.7") * acb("1.7"), arb(3) / 2, 120, 120, "lin").real
    rel = float(abs(a - b) / abs(b)); rec["e_phi_ts_vs_phi_z_lin"] = {"w": "1.7", "rel": rel}
    print(f"(e) lin phi_ts vs phi_z at w = 1.7: rel {rel:.2e}", flush=True)
    worst = max(worst, rel)
    rec["worst"] = worst
    rec["gate_1e-8"] = "PASS" if worst < 1e-8 else "FAIL"
    print("WORST rel:", worst, "GATE", rec["gate_1e-8"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
