"""NP record, Stage B: the massive-Coulomb limit at P2.  IB_PRIME = p0.
At z in {z9, -121/9, -25/9, two random}: 64 random eps = m5^2 (58 fit, 6 held out); every entry of (c, A) rebuilt
as a rational function of eps; order of vanishing at eps = 0 of each denominator, and its rational eps-roots.
Code control first: the engine (massless) vs the saved P2 p0 top row at 2 random z."""
from __future__ import annotations

import json
import random
import statistics
from fractions import Fraction as F
from pathlib import Path

from flint import nmod_poly, nmod

from np_engine import P, run, ratrec, TOP, fr
from ib_kite_de import fit_rational

DATA = Path(__file__).parent / "data"


def _main():
    out = DATA / "np_stageB.json"
    assert not out.exists()
    rng = random.Random(f"np_stageB|{P}")
    P2 = [fr(F(x) ** 2) for x in (1, F(3, 2), F(7, 3), F(8, 3))]

    # code control
    sv = json.loads((DATA / "ib_pko_p2_p0_op.json").read_text())
    assert sv["P"] == P
    zc = [rng.randrange(10**6, P - 1) for _ in range(2)]
    resc, _ = run([(z, P2 + [0], False, "top") for z in zc], nproc=2)
    ms = sv["masters"]
    i17 = ms.index(TOP)
    okc = True
    for r in resc:
        for j, lab in enumerate(ms):
            nu, de = sv["A"][i17][j]
            N, D = nmod_poly(nu, P), nmod_poly(de, P)
            val = int(N(r["z"]) * (nmod(1, P) / D(r["z"]))) if not N.is_zero() else 0
            if val != r["A"][TOP].get(lab, 0) % P:
                okc = False
    print("code control (engine vs saved P2 top row, 2 z):", okc, flush=True)

    zpts = {"z9": fr(F(-111, 1120)), "-121/9": fr(F(-121, 9)), "-25/9": fr(F(-25, 9)),
            "rand1": rng.randrange(10**6, P - 1), "rand2": rng.randrange(10**6, P - 1)}
    NE = 64
    eps = {k: [rng.randrange(10**6, P - 1) for _ in range(NE)] for k in zpts}
    args = [(zpts[k], P2 + [e], True, "all") for k in zpts for e in eps[k]]
    # pre-check: first 7 at z9 for timing + master count
    pre, wpre = run(args[:7])
    print("pre-check at z9: masters", [len(r["masters"]) for r in pre], "median cpu",
          f"{statistics.median(r['cpu'] for r in pre):.1f}s", flush=True)
    if statistics.median(r["cpu"] for r in pre) > 120:
        raise SystemExit("STOP: massive per-sample cost above cap")
    res, wall = run(args[7:])
    res = pre + res
    wall += wpre
    print(f"{len(res)} massive samples, wall {wall:.0f}s", flush=True)

    result = {"P": P, "code_control": okc, "wall": wall, "cpu": [r["cpu"] for r in res], "zpts": zpts}
    per = {}
    pos = 0
    for k in zpts:
        grp = res[pos:pos + NE]
        pos += NE
        mcount = sorted({len(r["masters"]) for r in grp})
        labs = grp[0]["masters"]
        if any(r["masters"] != labs for r in grp):
            per[k] = {"masters": mcount, "error": "basis differs"}
            print(k, "basis differs", mcount)
            continue
        fitp, testp = grp[:58], grp[58:]
        xs = [r["m2"][4] for r in fitp]
        ents = {}
        fails = 0
        for m in labs:
            for key, get in [(f"c|{m}", lambda r, m=m: r["c"][m])] + \
                    [(f"A|{m}|{n}", (lambda r, m=m, n=n: r["A"][m].get(n, 0))) for n in labs]:
                ys = [get(r) % P for r in fitp]
                if all(y == 0 for y in ys):
                    continue
                f = fit_rational(xs, ys, len(xs) - 3)
                if f is None:
                    fails += 1
                    continue
                N, D = f
                if not all((N(r["m2"][4]) - nmod(get(r), P) * D(r["m2"][4])) == 0 for r in testp):
                    fails += 1
                    continue
                x = nmod_poly([0, 1], P)
                v = 0
                DD = D
                while int(DD(0)) == 0:
                    DD, v = DD // x, v + 1
                vn = 0
                NN = N
                while not NN.is_zero() and int(NN(0)) == 0:
                    NN, vn = NN // x, vn + 1
                roots = []
                for fct, e in DD.factor()[1]:
                    if fct.degree() == 1:
                        co = [int(c) for c in fct.coeffs()]
                        roots.append(str(ratrec((-co[0]) * pow(co[1], P - 2, P) % P)))
                ents[key] = {"v": v, "vnum": vn, "dn": N.degree(), "dd": D.degree(), "roots": roots}
        maxdeg = max((e["dn"] + e["dd"] for e in ents.values()), default=0)
        per[k] = {"masters": mcount, "fails": fails, "entries": ents, "maxdeg": maxdeg}
        print(k, "masters", mcount, "entries", len(ents), "failed fits", fails, "max eps-degree", maxdeg, flush=True)
    result["per"] = per
    # comparison
    if all("entries" in per[k] for k in zpts):
        gen = {}
        for key in set(per["rand1"]["entries"]) | set(per["rand2"]["entries"]):
            gen[key] = max(per[r]["entries"].get(key, {"v": 0})["v"] for r in ("rand1", "rand2"))
        vsame = all(per["rand1"]["entries"].get(k, {"v": 0})["v"] == per["rand2"]["entries"].get(k, {"v": 0})["v"]
                    for k in gen)
        print("random z agree on every v:", vsame)
        comp = {}
        for k in ("z9", "-121/9", "-25/9"):
            up = [(key, e["v"], gen.get(key, 0)) for key, e in per[k]["entries"].items() if e["v"] > gen.get(key, 0)]
            rts = sorted({r for e in per[k]["entries"].values() for r in e["roots"]})
            comp[k] = {"raised": up, "roots": rts}
            print(k, "entries with v above generic:", len(up), up[:6])
            print("   eps-roots seen:", rts[:40])
        rts_gen = sorted({r for z in ("rand1", "rand2") for e in per[z]["entries"].values() for r in e["roots"]})
        print("eps-roots at random z:", rts_gen[:40])
        result["comparison"] = comp
        result["random_agree"] = vsame
    out.write_text(json.dumps(result))



if __name__ == "__main__":
    _main()
