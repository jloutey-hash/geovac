"""KS route O (record docs/predictions/2026-10-28_kite_speed_test.md): the kite's exact order-11 operator
(data/ib_pf_op_Q.json) walked by oe_engine + oe_kite.Walker (both unchanged, imported) to the nodes of the shared
outer grid only, then stage 2 in float64 (ks_common).
Segments timed by process_time inside this one process: load (operator JSON, polynomials, leading-coefficient roots),
input (phi_0..phi_5 by pf_4b_series_direct.u_chunk in-process, one process, 80 x 80 'sub' nodes, at D + 10 digits),
walk (series at 0, main path, nodes), stage 2.
Path rule as oe_kite.main, with one change forced by the shorter grid: the detour around a real apparent root rho is
used only if rho * 0.65 < zmax (waypoints past the last node are not walked).
usage: python ks_route_o.py <T12|T9> <D> <ok|pert>"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json, sys, time
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx, fmpz_poly
import oe_engine as O
from oe_kite import Walker, cplx, s
import pf_4b_series_direct as PF
import pf_common as P
import ks_common as C

DATA = Path(__file__).parent / "data"


def run(target, D, pert=False):
    nd = C.outer_nodes(target)
    tm = {}
    t0 = time.process_time()
    tm["before_run"] = t0
    ctx.dps = D
    opd = json.loads((DATA / "ib_pf_op_Q.json").read_text())
    op = O.Operator(opd["polys"])
    if pert:
        op = op.perturbed(5, 40, 201, 200)
    lead = fmpz_poly([int(v) for v in opd["polys"][-1]])
    roots = []
    for f, e in lead.factor()[1]:
        roots += [cplx(rr[0]) for rr in f.complex_roots()]
    real_app = sorted(r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0)
    t1 = time.process_time(); tm["load"] = t1 - t0
    # input: six Taylor coefficients, one process, in-process
    din = D + 10
    cs, _rad = PF.u_chunk((6, din, 80, 80, list(range(80)), "sub"))
    ctx.dps = din
    a, b, g1, g2 = P.params()
    pref = 2 * g1 * g2 * a * a * a * b * b * b
    free = [arb(c) * pref for c in cs[:6]]
    ctx.dps = D
    free = [arb(x.mid().str(D + 5, radius=False)) for x in free]
    t2 = time.process_time(); tm["input"] = t2 - t1
    ws = [arb(w) for w, _ in nd]
    zs = [w * w for w in ws]
    zf = [float(mp.mpf(s(z, 20))) for z in zs]
    zmax = max(zf)
    order = op.order
    a0, resid = op.series_at_zero({n: free[n] for n in range(6)}, int(D * 3.33) + 40)
    a0c = [acb(x) for x in a0]
    phi = [None] * len(zs)
    for i, z in enumerate(zs):
        if zf[i] <= 1 / 16:
            phi[i] = O.eval_series(a0c, acb(z))[0].real
    W = Walker(op, roots, D)
    wps = []
    for rho in real_app:
        if rho * 0.65 < zmax:
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
    main_steps = W.steps
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
    phif = [float(mp.mpf(s(p, 20))) for p in phi]
    t3 = time.process_time(); tm["walk"] = t3 - t2
    Y = C.stage2([float(w) for w, _ in nd], [float(d) for _, d in nd], phif)
    t4 = time.process_time(); tm["stage2"] = t4 - t3
    info = {"main_steps": main_steps, "steps_total": W.steps, "branch_nodes": nbranch, "waypoints": len(wps)}
    return Y, tm, info, phif


if __name__ == "__main__":
    tg, D, mode = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    Y, tm, info, phif = run(tg, D, mode == "pert")
    print(json.dumps({"route": "O", "target": tg, "D": D, "mode": mode, "Y": list(map(float, Y)),
                      "rel": C.rel_errs(Y), "cpu": tm, "info": info}))
