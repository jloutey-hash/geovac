"""PF1-PF4 stage 1 (docs/predictions/2026-10-23_price_formula_slots.md): checks of the off-slice kite before any series.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions,
pf_common).
 (a) the closed-form loop I_221(gamma_1, gamma_2, 0) against direct 3-D quadrature at five random triples (gate 1e-8);
 (b) route A (momentum triangle, closed-form loop) against route B (kite, Feynman parameters, k2 = (1-u-v)^2 z) at
     three w (gate 1e-14 at mp.dps 20 / 40 arb digits); the slice-form assignment k2 = (u-v)^2 z at the same w (must
     disagree: the discriminating control);
 (c) slice members: phi_z at (1,1;2,2) and (1,1;3,3) against ring_ext_QG_common.phi_gem_z at three z (one complex),
     50 digits;
 (d) phi_ts (tanh-sinh, real w) against phi_z at one real w;
 (e) phi(0) is not in pi Q (PSLQ on [phi(0), pi, 1] at 60 digits with a height bound) -- the reason the exact mod-p
     route of the slice does not apply off the slice.
Gate: (a) < 1e-8, (b) < 1e-14 for the derived assignment, (c), (d) < 1e-40."""
from __future__ import annotations
import json, random, sys, time
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm
import pf_common as P

OUT = Path(__file__).parent / "data" / "pf_1_checks.json"


def main():
    rec = {"object": "a = 1, b = 2, gamma_1 = 3/2, gamma_2 = 3, t on the apex"}
    rng = random.Random(20261023)
    worst_a = 0.0
    mp.mp.dps = 30
    g1, g2 = P.mpq(P.G1), P.mpq(P.G2)
    rows = []
    for k in range(5):
        x1 = [rng.uniform(-2, 2) for _ in range(3)]
        x2 = [rng.uniform(-2, 2) for _ in range(3)]
        s12 = sum((mp.mpf(p) - q) ** 2 for p, q in zip(x1, x2))
        s13 = sum(mp.mpf(p) ** 2 for p in x1)
        s23 = sum(mp.mpf(q) ** 2 for q in x2)
        closed = Cm.I221(g1, g2, 0, s12, s13, s23)
        dd = P.direct3d(x1, x2, float(g1), float(g2))
        rel = abs(dd / float(closed) - 1)
        rows.append({"x1": x1, "x2": x2, "direct": dd, "closed": float(closed), "rel": rel})
        print(f"(a) loop (3/2, 3, 0) triple {k}: direct {dd:.15e} closed {float(closed):.15e} rel {rel:.2e}", flush=True)
        worst_a = max(worst_a, rel)
    rec["a_loop_vs_direct3d"] = {"rows": rows, "worst": worst_a, "gate_1e-8": "PASS" if worst_a < 1e-8 else "FAIL"}
    # (b) route A vs route B
    mp.mp.dps = 20
    ctx.dps = 40
    rows = []
    worst_b, best_wrong = 0.0, 1.0
    t0 = time.time()
    for w in ("0.5", "0.9", "2.1"):
        a = P.phi_routeA(w)
        zz = acb(w) * acb(w)
        b = P.phi_z(zz, 200, 200).real
        bw = P.phi_z(zz, 200, 200, assign="wrong").real
        rel = float(abs(Cm.to_mp(b) / a - 1))
        relw = float(abs(Cm.to_mp(bw) / a - 1))
        rows.append({"w": w, "routeA": mp.nstr(a, 18), "routeB_correct": b.mid().str(18, radius=False),
                     "routeB_wrong_assignment": bw.mid().str(18, radius=False), "rel_correct": rel, "rel_wrong": relw})
        print(f"(b) w = {w}: route A {mp.nstr(a, 17)}  route B {b.mid().str(17, radius=False)}  rel {rel:.2e};  "
              f"slice-form assignment rel {relw:.2e}  ({time.time()-t0:.0f}s)", flush=True)
        worst_b = max(worst_b, rel)
        best_wrong = min(best_wrong, relw)
    rec["b_routeA_vs_routeB"] = {"rows": rows, "worst_correct": worst_b, "gate_1e-14": "PASS" if worst_b < 1e-14 else "FAIL",
                                 "wrong_assignment_smallest_rel": best_wrong,
                                 "wrong_assignment_disagrees": best_wrong > 1e-3}
    # (c) slice members against the stored QG kite
    ctx.dps = 50
    rows = []
    worst_c = 0.0
    for g in ("2", "3"):
        for z in (acb("0.7"), acb("-1.3"), acb("0.9", "1.4")):
            a = P.phi_z(z, 160, 160, a="1", b="1", g1=g, g2=g)
            b = Cm.phi_gem_z(z, arb(g), 160, 160)
            rel = float(abs(a - b) / abs(b))
            rows.append({"gamma": g, "z": str(z), "rel": rel})
            print(f"(c) slice (1,1;{g},{g}) at z = {z}: rel {rel:.2e}", flush=True)
            worst_c = max(worst_c, rel)
    rec["c_slice_vs_QG_kite"] = {"rows": rows, "worst": worst_c}
    # (d) tanh-sinh form
    ctx.dps = 45
    a = P.phi_ts(arb("1.7"), arb(1) / 32, 35).real
    b = P.phi_z(acb("1.7") * acb("1.7"), 200, 200).real
    rel_d = float(abs(a - b) / abs(b))
    rec["d_phi_ts_vs_phi_z"] = {"w": "1.7", "rel": rel_d}
    print(f"(d) phi_ts vs phi_z at w = 1.7: rel {rel_d:.2e}", flush=True)
    # (e) phi(0) against pi Q
    ctx.dps = 80
    p0 = P.phi_z(acb(0), 300, 200).real
    p0s = Cm.phi_gem_z(acb(0), arb(2), 300, 200).real
    mp.mp.dps = 70
    v0 = Cm.to_mp(p0)
    v0s = Cm.to_mp(p0s)
    rel_pslq = mp.pslq([v0, mp.pi], maxcoeff=10**20, maxsteps=10**6)
    rel_pslq_s = mp.pslq([v0s, mp.pi], maxcoeff=10**20, maxsteps=10**6)
    rec["e_phi0"] = {"phi0": mp.nstr(v0, 60), "phi0_over_pi": mp.nstr(v0 / mp.pi, 60),
                     "pslq_phi0_pi_height_1e20": rel_pslq,
                     "slice_control_gamma2_phi0_over_pi": mp.nstr(v0s / mp.pi, 40), "slice_control_pslq": rel_pslq_s}
    print("(e) phi(0) =", mp.nstr(v0, 30), "; phi(0)/pi =", mp.nstr(v0 / mp.pi, 30), "; PSLQ [phi0, pi]:", rel_pslq,
          "; slice control gamma = 2: phi0/pi =", mp.nstr(v0s / mp.pi, 25), "PSLQ:", rel_pslq_s, flush=True)
    rec["gate"] = "PASS" if (worst_a < 1e-8 and worst_b < 1e-14 and worst_c < 1e-40 and rel_d < 1e-30 and best_wrong > 1e-3) else "FAIL"
    print("GATE:", rec["gate"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
