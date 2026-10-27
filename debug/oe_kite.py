"""OE1-OE3 (docs/predictions/2026-10-26_operator_evaluator.md, committed 0a98d7e3 before computing): the off-slice
kite phi(z) = Phi~(sqrt z) at the 247 nodes of the 60-digit Y grid, from its exact order-11 IBP operator
(data/ib_pf_op_Q.json) and its first six Taylor coefficients at z = 0, then Y(t) by the grid's own w-sum.

Input modes: 'ref' takes phi_0..phi_5 from the 280-digit series (control); 'fresh' recomputes them with
pf_4b_series_direct.py at <input dps> digits on 80 x 80 nodes on ONE process, timed; 'pert' is 'ref' with one operator
coefficient times 201/200 (negative control).

Path. The series at 0 serves |z| <= 1/16 (radius 1/4, set by -1/4). Then a walker: from centre c towards a target T,
step to c + (r/2)(T - c)/|T - c|, r = distance from c to the nearest root of the leading coefficient (genuine or
apparent: rounding excites solutions of the recurrence that blow up at apparent points too). Waypoints along the
positive axis with a detour above each real apparent root rho (rho(1-0.35) -> rho(1+0.35i) -> rho(1+0.35)). A node
is evaluated from the main-path centre where |node - c| / r is smallest if that ratio is <= 1/2; otherwise by a
branch walk from that centre. Each disk is evaluated at |v| <= 1/2 of its radius, so the series converge at ratio
1/2 at worst.
usage: python oe_kite.py <ref|fresh|pert> <work digits> [input dps]"""
from __future__ import annotations
import json, math, subprocess, sys, time
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx, fmpz_poly
sys.path.insert(0, str(Path(__file__).parent))
import oe_engine as O

DATA = Path(__file__).parent / "data"
TS = ("0.8", "1.3", "1.7", "2.5")


def s(x, d=4):
    if isinstance(x, acb):
        x = x.real
    return x.mid().str(d, radius=False)


def rel(x, y):
    q = abs(x - y) / abs(y)
    try:
        return float(mp.mpf(s(q, 6)))
    except Exception:
        return float("inf")


def cplx(x):
    return complex(float(mp.mpf(x.real.mid().str(25, radius=False))), float(mp.mpf(x.imag.mid().str(25, radius=False))))


class Walker:
    def __init__(self, op, roots, D):
        self.op, self.roots, self.D = op, roots, D
        self.order = op.order
        self.nterms = int(D * 3.33) + 30
        self.steps = 0

    def dist(self, cf):
        return min(abs(cf - r) for r in self.roots)

    def series(self, c, der):
        """Local series at c (scale h = r/2) from the derivatives der[k] = f^(k)(c)."""
        r = self.dist(cplx(c))
        h = acb(arb(r / 2))
        loc = [mid_mul(der[k], h ** k) / math.factorial(k) for k in range(self.order)]
        R = self.op.local(c, h, max(4 * self.D, 1200))
        a = O.step_series(R, self.order, loc, self.nterms)
        self.steps += 1
        return a, h

    def walk(self, c, der, T):
        """From (c, der) towards T; returns the list of (c, der, a, h) disks visited, the last containing T."""
        out = []
        while True:
            a, h = self.series(c, der)
            out.append((c, der, a, h))
            cf, Tf = cplx(c), cplx(T)
            hf = abs(cplx(h))
            if abs(Tf - cf) <= hf:          # T within |v| <= 1, i.e. half the radius
                return out
            u = (Tf - cf) / abs(Tf - cf)
            v = acb(arb(repr(u.real)), arb(repr(u.imag)))
            ders = O.eval_series(a, v, self.order - 1)
            c = c + h * v
            der = [ders[k] / h ** k for k in range(self.order)]


def mid_mul(x, y):
    return O.mid(x * y)


def main():
    mode, D = sys.argv[1], int(sys.argv[2])
    din = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    ctx.dps = D
    opd = json.loads((DATA / "ib_pf_op_Q.json").read_text())
    op = O.Operator(opd["polys"])
    if mode == "pert":
        op = op.perturbed(5, 40, 201, 200)
    ref = [arb(v) for v in json.loads((DATA / "pf_series_direct_explore_N300_d280_460x460_sub.json").read_text())["phi_n"]]
    log = {"mode": mode, "work_digits": D}
    T0 = time.time()
    if mode == "fresh":
        t = time.time()
        r = subprocess.run([sys.executable, str(Path(__file__).parent / "pf_4b_series_direct.py"), "run", "6", str(din),
                            "80", "80", "oe", "1", "sub"], capture_output=True, text=True)
        log["input_seconds"] = time.time() - t
        log["input_stdout"] = r.stdout[-300:] + r.stderr[-300:]
        fn = DATA / f"pf_series_direct_oe_N6_d{din}_80x80_sub.json"
        free = [arb(v) for v in json.loads(fn.read_text())["phi_n"][:6]]
        log["input_vs_ref"] = [rel(free[n], ref[n]) for n in range(6)]
    else:
        free = ref[:6]
    t = time.time()
    lead = fmpz_poly([int(v) for v in opd["polys"][-1]])
    roots = []
    for f, e in lead.factor()[1]:
        roots += [cplx(rr[0]) for rr in f.complex_roots()]
    real_app = sorted(r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0)
    log["real_apparent_roots"] = real_app

    grid = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    ws = [arb(w) for w in grid["w"]]
    zs = [w * w for w in ws]
    zf = [float(mp.mpf(s(z, 20))) for z in zs]
    zmax = max(zf)
    order = op.order

    a0, resid = op.series_at_zero({n: free[n] for n in range(6)}, int(D * 3.33) + 40)
    log["series0_free_residuals"] = {n: s(abs(v), 3) for n, v in resid.items()}
    log["recurrence_vs_ref_6_20"] = max(rel(a0[n], ref[n]) for n in range(6, 21)) if mode != "pert" else None
    a0c = [acb(x) for x in a0]
    phi = [None] * len(zs)
    for i, z in enumerate(zs):
        if zf[i] <= 1 / 16:
            phi[i] = O.eval_series(a0c, acb(z))[0].real
    # main path
    W = Walker(op, roots, D)
    wps = []
    for rho in real_app:
        wps += [rho * 0.65, rho * complex(1, 0.35), rho * 1.35]
    wps.append(zmax * 1.01)
    c = acb(arb(1) / 16)
    der = O.eval_series(a0c, c, order - 1)
    disks = []
    for T in wps:
        Ta = acb(arb(repr(T.real)), arb(repr(T.imag))) if isinstance(T, complex) else acb(arb(repr(T)))
        seg = W.walk(c, der, Ta)
        disks += seg
        cl, dl, al, hl = seg[-1]
        v = (Ta - cl) / hl
        ders = O.eval_series(al, v, order - 1)
        c, der = Ta, [ders[k] / hl ** k for k in range(order)]
    log["main_steps"] = W.steps
    # nodes
    nbranch = 0
    for i, z in enumerate(zs):
        if phi[i] is not None:
            continue
        best = min(disks, key=lambda dk: abs(zf[i] - cplx(dk[0])) / abs(cplx(dk[3])))
        cb, db, ab, hb = best
        if abs(zf[i] - cplx(cb)) <= abs(cplx(hb)):
            phi[i] = O.eval_series(ab, (acb(z) - cb) / hb)[0].real
        else:
            seg = W.walk(cb, db, acb(z))
            cl, dl, al, hl = seg[-1]
            phi[i] = O.eval_series(al, (acb(z) - cl) / hl)[0].real
            nbranch += 1
    log["branch_nodes"] = nbranch
    log["steps_total"] = W.steps
    log["seconds_kite_total"] = time.time() - T0
    log["seconds_operator_part"] = time.time() - t
    gp = [arb(v) for v in grid["phi"]]
    errs = [rel(phi[i], gp[i]) for i in range(len(zs))]
    byd = {}
    for i, e in enumerate(errs):
        dec = int(math.floor(math.log10(float(mp.mpf(s(ws[i], 10))))))
        byd[dec] = max(byd.get(dec, 0.0), e)
    log["node_err_by_decade_of_w"] = byd
    log["node_err_max"] = max(errs)
    import ring_ext_QG_3_values as V
    recO = dict(grid)
    recO["phi"] = [s(p, 75) for p in phi]
    log["Y_rel_err"] = {}
    for tt in TS:
        yg = V.Y_theta(grid, tt, 0)[0]
        yo = V.Y_theta(recO, tt, 0)[0]
        log["Y_rel_err"][tt] = mp.nstr(abs(yo / yg - 1), 3)
    t = time.time()
    for tt in TS:
        V.Y_theta(recO, tt, 0)
    log["seconds_Y_sums_4t"] = time.time() - t
    print(json.dumps(log, indent=1), flush=True)
    (DATA / f"oe_kite_{mode}_d{D}.json").write_text(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
