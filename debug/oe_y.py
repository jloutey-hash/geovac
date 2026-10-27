"""OE4 (docs/predictions/2026-10-26_operator_evaluator.md, committed 0a98d7e3 before computing): Y's own minimal
operator M13 (tests/data/paper59_y_2026_10_26.json; sum_k c_k(x) theta_t^k, x = t^2) as an evaluator of Y in t, from
theta_t^j Y (j = 0..12) at t0 = 1.7 on the 60-digit grid, against one grid w-sum per t.

M13 is rewritten in d/dx form, P_i(x) = sum_k c_k(x) 2^k S2(k, i) x^i (theta_t = 2 theta_x,
theta_x^k = sum_i S2(k, i) x^i D^i). The walker of oe_kite.py steps at half the distance to the nearest root of the
leading coefficient; the path rises to height H = 3/4 above the axis, runs across, and comes down, avoiding the cut
points and the apparent points on the axis.
usage: python oe_y.py <work digits> [pert]"""
from __future__ import annotations
import json, math, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
from flint import arb, acb, ctx, fmpz_poly
sys.path.insert(0, str(Path(__file__).parent))
import oe_engine as O
from oe_kite import Walker, cplx, s

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).parent.parent
T0S = "1.7"
TARGETS = ("0.8", "1.3", "2.5", "2.6", "3.8")
H = 0.75


def dform(M):
    """theta_t-form integer coefficient lists -> d/dx-form integer coefficient lists."""
    order = len(M) - 1
    P = []
    for i in range(order + 1):
        acc = {}
        for k in range(i, order + 1):
            w = (2 ** k) * int(sp.functions.combinatorial.numbers.stirling(k, i, kind=2))
            for j, cf in enumerate(M[k]):
                acc[j + i] = acc.get(j + i, 0) + w * int(cf)
        deg = max(acc) if acc else 0
        P.append([acc.get(j, 0) for j in range(deg + 1)])
    return P


def main():
    D = int(sys.argv[1])
    pert = len(sys.argv) > 2 and sys.argv[2] == "pert"
    ctx.dps = D
    mp.mp.dps = 70
    M = json.loads((ROOT / "tests" / "data" / "paper59_y_2026_10_26.json").read_text())["M13"]
    P = dform(M)
    op = O.Operator(P)
    if pert:
        op = op.perturbed(7, 30, 201, 200)
    order = op.order
    lead = fmpz_poly(P[-1])
    roots = []
    for f, e in lead.factor()[1]:
        roots += [cplx(rr[0]) for rr in f.complex_roots()]
    log = {"work_digits": D, "pert": pert, "order": order,
           "real_roots": sorted(round(r.real, 6) for r in roots if abs(r.imag) < 1e-12)}
    import ring_ext_QG_8_onY as E
    grid = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    tq = time.time()
    th = E.thetaY(grid, T0S, order - 1)
    log["seconds_grid_initial_data"] = time.time() - tq
    x0 = arb(T0S) ** 2
    thx = [arb(mp.nstr(th[j], 66)) / 2 ** j for j in range(order)]
    der = []
    for k in range(order):
        sacc = arb(0)
        for i in range(k + 1):
            s1 = int(sp.functions.combinatorial.numbers.stirling(k, i, kind=1, signed=True))
            sacc += s1 * thx[i]
        der.append(acb(sacc / x0 ** k))
    W = Walker(op, roots, D)
    out = {}
    for tt in TARGETS:
        xt = arb(tt) ** 2
        xtf = float(mp.mpf(tt)) ** 2
        x0f = float(mp.mpf(T0S)) ** 2
        st = time.time()
        steps0 = W.steps
        c, d = acb(x0), list(der)
        for T in (complex(x0f, H), complex(xtf, H), complex(xtf, 0)):
            Ta = acb(arb(repr(T.real)), arb(repr(T.imag))) if T.imag else acb(xt)
            seg = W.walk(c, d, Ta)
            cl, dl, al, hl = seg[-1]
            v = (Ta - cl) / hl
            ders = O.eval_series(al, v, order - 1)
            c, d = Ta, [ders[k] / hl ** k for k in range(order)]
        yo = d[0]
        sec = time.time() - st
        tg = time.time()
        yg = E.thetaY(grid, tt, 0)[0]
        secg = time.time() - tg
        err = abs(yo.real - arb(mp.nstr(yg, 66))) / abs(arb(mp.nstr(yg, 66)))
        out[tt] = {"rel_err": s(err, 3), "imag": s(abs(yo.imag), 3), "seconds_operator": round(sec, 3),
                   "steps": W.steps - steps0, "seconds_grid_sum": round(secg, 4)}
        print(tt, out[tt], flush=True)
    log["targets"] = out
    tag = "_pert" if pert else ""
    (DATA / f"oe_y_d{D}{tag}.json").write_text(json.dumps(log, indent=1))
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
