"""NP record, Stage C: z9 as a function of all four rates (a^2 = 1, X = (b^2, g1^2, g2^2) random residues).
  IB_PRIME=p0 python np_stageC.py line            -> C1: random line u + t v, 36 fit + 6 held out; degrees (dn, dd)
  IB_PRIME=p0 python np_stageC.py dense dn dd     -> C2: dense ansatz from U + 12 points, 12 held out, lift to Q
  IB_PRIME=p1 python np_stageC.py check           -> 6 random points at p1 vs the lifted formula"""
from __future__ import annotations

import itertools
import json
import random
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

from flint import nmod, nmod_mat

from np_engine import P, run, readout, ratrec, fr
from ib_kite_de import fit_rational

DATA = Path(__file__).parent / "data"


def measure(points, rng):
    """points: list of (b2, g12, g22) residues.  Returns dict i -> z9 (4-sample readout only) and failures."""
    args, idx = [], []
    for i, X in enumerate(points):
        for _ in range(4):
            args.append((rng.randrange(10**6, P - 1), [1] + [x % P for x in X] + [0], False, "top"))
            idx.append(i)
    res, wall = run(args)
    groups = {}
    for r, i in zip(res, idx):
        groups.setdefault(i, []).append(r)
    vals, fails = {}, []
    for i in range(len(points)):
        g = groups[i]
        if any(len(r["masters"]) != 18 for r in g):
            fails.append((i, "masters"))
            continue
        r0, ok, note = readout(g)
        if ok:
            vals[i] = r0
        else:
            fails.append((i, note))
    return vals, fails, wall, [r["cpu"] for r in res]


def monos(d):
    return [m for m in itertools.product(range(d + 1), repeat=3) if sum(m) <= d]


def ev(mon, X):
    v = 1
    for e, x in zip(mon, X):
        v = v * pow(x, e, P) % P
    return v


def _main():
    mode = sys.argv[1]
    rng = random.Random(f"np_stageC|{mode}|{P}")
    if mode == "line":
        out = DATA / "np_stageC_line.json"
        assert not out.exists()
        u = [rng.randrange(10**6, P - 1) for _ in range(3)]
        v = [rng.randrange(10**6, P - 1) for _ in range(3)]
        ts = [rng.randrange(10**6, P - 1) for _ in range(42)]
        pts = [[(u[k] + t * v[k]) % P for k in range(3)] for t in ts]
        vals, fails, wall, cpu = measure(pts, rng)
        print(f"wall {wall:.0f}s, readout failures {fails}")
        assert not fails
        xs, ys = ts[:36], [vals[i] for i in range(36)]
        f = fit_rational(xs, ys, len(xs) - 4)
        rec = {"P": P, "wall": wall, "cpu": cpu, "fails": fails}
        if f is None:
            print("NO reconstruction on the line within total degree 32")
            rec["deg"] = None
        else:
            N, D = f
            held = [(N(t) - nmod(vals[36 + i], P) * D(t)) == 0 for i, t in enumerate(ts[36:])]
            dn, dd = N.degree(), D.degree()
            U = comb(dn + 3, 3) + comb(dd + 3, 3) - 1
            print(f"line degrees num {dn} den {dd}; held-out {held}; dense unknowns U = {U}")
            rec.update({"deg": [dn, dd], "held": held, "U": U})
        out.write_text(json.dumps(rec))
    elif mode == "dense":
        out = DATA / "np_stageC_dense.json"
        assert not out.exists()
        dn, dd = int(sys.argv[2]), int(sys.argv[3])
        MN, MD = monos(dn), monos(dd)
        U = len(MN) + len(MD) - 1
        assert U <= 250, "cap"
        nfit, nhold = U + 12, 12
        pts = [[rng.randrange(10**6, P - 1) for _ in range(3)] for _ in range(nfit + nhold)]
        vals, fails, wall, cpu = measure(pts, rng)
        print(f"{len(pts)} points, wall {wall:.0f}s, readout failures {len(fails)} {fails[:5]}", flush=True)
        fit_idx = [i for i in range(nfit) if i in vals]
        rows = []
        for i in fit_idx:
            X, z = pts[i], vals[i]
            rows.append([ev(m, X) for m in MN] + [(-z * ev(m, X)) % P for m in MD])
        M = nmod_mat(rows, P)
        rk = M.rank()
        ncol = len(MN) + len(MD)
        print("rank", rk, "of", ncol, "columns; rows", len(rows), "(need rank = cols - 1)")
        rec = {"P": P, "wall": wall, "cpu": cpu, "fails": fails, "dn": dn, "dd": dd, "U": U, "rank": rk}
        (DATA / "np_stageC_dense_samples.json").write_text(json.dumps({"P": P, "pts": pts, "vals": vals}))
        if rk == ncol - 1:
            R, rr = M.rref()
            piv = []
            for i in range(rr):
                piv.append(next(j for j in range(ncol) if int(R[i, j]) != 0))
            free = [j for j in range(ncol) if j not in piv]
            assert len(free) == 1
            full = [0] * ncol
            full[free[0]] = 1
            for i, j in enumerate(piv):
                full[j] = (-int(R[i, free[0]])) % P
            # normalisation rule (fixed in advance of the lift): first nonzero D coefficient in MD order = 1
            jn = next(j for j in range(len(MN), ncol) if full[j])
            inv = pow(full[jn], P - 2, P)
            full = [x * inv % P for x in full]
            Nco, Dco = full[:len(MN)], full[len(MN):]

            def zf(X):
                n = sum(c * ev(m, X) for c, m in zip(Nco, MN)) % P
                d = sum(c * ev(m, X) for c, m in zip(Dco, MD)) % P
                return n * pow(d, P - 2, P) % P if d else None
            held = [zf(pts[i]) == vals.get(i) for i in range(nfit, nfit + nhold)]
            print("held-out (12):", held)
            Nq = [ratrec(c) for c in Nco]
            Dq = [ratrec(c) for c in Dco]
            print("lift complete:", None not in Nq + Dq)
            rec.update({"held": held, "N": {str(m): str(c) for m, c in zip(MN, Nq)},
                        "D": {str(m): str(c) for m, c in zip(MD, Dq)}})
            if None not in Nq + Dq:
                X2 = (F(9, 4), F(49, 9), F(64, 9))
                n2 = sum(c * X2[0] ** m[0] * X2[1] ** m[1] * X2[2] ** m[2] for c, m in zip(Nq, MN))
                d2 = sum(c * X2[0] ** m[0] * X2[1] ** m[1] * X2[2] ** m[2] for c, m in zip(Dq, MD))
                print("formula at P2:", n2 / d2 if d2 else "pole", "(measured -111/1120)")
                rec["P2"] = str(n2 / d2) if d2 else "pole"
        out.write_text(json.dumps(rec))
    elif mode == "check":
        out = DATA / "np_stageC_check.json"
        assert not out.exists()
        dz = json.loads((DATA / "np_stageC_dense.json").read_text())
        MN = [tuple(json.loads(k.replace("(", "[").replace(")", "]"))) for k in dz["N"]]
        MD = [tuple(json.loads(k.replace("(", "[").replace(")", "]"))) for k in dz["D"]]
        Nq = [F(v) for v in dz["N"].values()]
        Dq = [F(v) for v in dz["D"].values()]
        pts = [[rng.randrange(10**6, P - 1) for _ in range(3)] for _ in range(6)]
        vals, fails, wall, cpu = measure(pts, rng)

        def zf(X):
            n = sum(fr(c) * ev(m, X) for c, m in zip(Nq, MN)) % P
            d = sum(fr(c) * ev(m, X) for c, m in zip(Dq, MD)) % P
            return n * pow(d, P - 2, P) % P
        chk = [zf(pts[i]) == vals.get(i) for i in range(6)]
        print("p1 check (6):", chk, "failures", fails)
        out.write_text(json.dumps({"P": P, "wall": wall, "check": chk, "fails": fails}))


if __name__ == "__main__":
    _main()
