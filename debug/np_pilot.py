"""NP record: pilot + readout controls.  Run with IB_PRIME = p0.
 - readout identity at P2 (4 samples; must give -111/1120)
 - scaling control at 2*P2 (4 samples; must give -111/280)
 - pilot: one random s on F1 = (1, 2, 3/2, g2^2 = s), 40 samples, full system rebuilt in z (35 fit, 5 held out);
   masters, common denominator factors, residual factor, carrier forms.  Writes debug/data/np_pilot.json."""
from __future__ import annotations

import json
import random
import statistics
from fractions import Fraction as F
from pathlib import Path

from flint import nmod_poly, nmod

from np_engine import P, run, readout, ratrec, TOP, CARRIERS, fr
from ib_kite_de import fit_rational

DATA = Path(__file__).parent / "data"


def _main():
    out = DATA / "np_pilot.json"
    assert not out.exists()
    rng = random.Random(f"np_pilot|{P}")


    def sq(*xs):
        return [fr(F(x) ** 2) for x in xs]


    P2 = sq(1, F(3, 2), F(7, 3), F(8, 3)) + [0]
    P2x2 = sq(2, 3, F(14, 3), F(16, 3)) + [0]
    s_pil = rng.randrange(10**6, P - 1)
    F1 = sq(1, 2, F(3, 2)) + [s_pil, 0]

    args = [(rng.randrange(10**6, P - 1), P2, False, "top") for _ in range(4)]
    args += [(rng.randrange(10**6, P - 1), P2x2, False, "top") for _ in range(4)]
    args += [(rng.randrange(10**6, P - 1), F1, False, "all") for _ in range(40)]
    res, wall = run(args)
    cpus = [r["cpu"] for r in res]
    print(f"{len(res)} samples, wall {wall:.0f}s, median cpu {statistics.median(cpus):.1f}s, max {max(cpus):.1f}s")
    r0, ok, note = readout(res[0:4])
    print("readout identity P2:", ok, note, ratrec(r0) if ok else None)
    r1, ok1, note1 = readout(res[4:8])
    print("scaling control 2*P2:", ok1, note1, ratrec(r1) if ok1 else None)

    pil = res[8:]
    ms = pil[0]["masters"]
    assert all(r["masters"] == ms for r in pil), "master basis differs"
    print("pilot masters", len(ms), "same labels as P2:", ms == json.loads((DATA / "ib_pko_p2_p0_op.json").read_text())["masters"])
    fitp, testp = pil[:35], pil[35:]
    xs = [r["z"] for r in fitp]


    def rebuild(get):
        ys = [get(r) % P for r in fitp]
        if all(y == 0 for y in ys):
            return None
        f = fit_rational(xs, ys, len(xs) - 3)
        assert f is not None, "fit failed"
        N, D = f
        for r in testp:
            assert (N(r["z"]) - nmod(get(r), P) * D(r["z"])) == 0, "held-out fails"
        return N, D


    ents = {}
    for m in ms:
        f = rebuild(lambda r, m=m: r["c"][m])
        if f:
            ents[("c", m)] = f
        for n in ms:
            f = rebuild(lambda r, m=m, n=n: r["A"][m].get(n, 0))
            if f:
                ents[("A", m, n)] = f
    print("entries", len(ents), "max total degree", max(N.degree() + D.degree() for N, D in ents.values()))
    a2, b2, g12 = fr(1), fr(4), fr(F(9, 4))
    g22 = s_pil
    X = {"a": a2, "b": b2, "g1": g12, "g2": g22}
    known = {"z": nmod_poly([0, 1], P), "sys(a2+g22)": nmod_poly([(a2 + g22) % P, 1], P),
             "sys(b2+g12)": nmod_poly([(b2 + g12) % P, 1], P)}
    for u, v in [("a", "b"), ("g1", "g2"), ("a", "g2"), ("b", "g1")]:
        s0 = (X[u] + X[v]) % P
        known[f"cut({u},{v})"] = nmod_poly([s0, 1], P) ** 2 - nmod_poly([4 * X[u] * X[v] % P], P)
    for u, v in [("a", "g1"), ("b", "g2")]:
        s0 = (X[u] + X[v]) % P
        known[f"adj({u},{v})"] = nmod_poly([s0, 1], P) ** 2 - nmod_poly([4 * X[u] * X[v] % P], P)


    def strip(Dp):
        rep = {}
        rest = Dp
        for k, q in known.items():
            e = 0
            while True:
                qq, rr = divmod(rest, q)
                if rr.is_zero() and q.degree() > 0:
                    rest, e = qq, e + 1
                else:
                    break
            if e:
                rep[k] = e
        # quadratic cut factors may split: also strip linear factors of them
        return rep, rest


    lcmA, lcmc = nmod_poly([1], P), nmod_poly([1], P)
    for k, (N, D) in ents.items():
        if k[0] == "A":
            lcmA = lcmA * D // lcmA.gcd(D)
        else:
            lcmc = lcmc * D // lcmc.gcd(D)
    info = {}
    for name, L in [("A", lcmA), ("c", lcmc)]:
        rep, rest = strip(L)
        fac = rest.factor()[1]
        info[name] = {"known": rep, "residual": [(f.degree(), e, [int(c) for c in f.coeffs()]) for f, e in fac]}
        print(name, "lcm deg", L.degree(), "known", rep, "residual factors (deg, mult)", [(f.degree(), e) for f, e in fac])
    car = {}
    for lab in CARRIERS:
        N, D = ents[("A", TOP, lab)]
        car[lab] = (N.degree(), D.degree(), [int(c) for c in D.coeffs()])
        print("carrier", lab, "num deg", N.degree(), "den deg", D.degree())
    rres = [(f, e) for f, e in lcmA.factor()[1]]
    rp, okp, notep = readout(pil[:4])
    print("4-sample readout on pilot:", okp, notep, rp)
    # which entries carry the residual root
    z9 = None
    resid = info["A"]["residual"]
    if len(resid) == 1 and resid[0][0] == 1:
        co = resid[0][2]
        z9 = (-co[0]) * pow(co[1], P - 2, P) % P
        print("residual root mod P", z9, "readout agrees:", z9 == rp)
        carriers_all = [k for k, (N, D) in ents.items() if int(D(z9)) == 0]
        print("entries carrying it:", len(carriers_all))
        for k in carriers_all:
            N, D = ents[k]
            if k[0] == "c":
                e = 0
                q = nmod_poly([(-z9) % P, 1], P)
                DD = D
                while (DD % q).is_zero():
                    DD, e = DD // q, e + 1
                print("  c entry", k[1], "power", e)
    out.write_text(json.dumps({"P": P, "wall": wall, "cpus": cpus, "s_pilot": s_pil, "readout_P2": r0,
                               "readout_P2_ratrec": str(ratrec(r0)) if ok else None, "scaling": str(ratrec(r1)) if ok1 else None,
                               "masters": ms, "lcm": info, "carriers": car, "pilot_readout": rp, "z9_resid": z9,
                               "samples": pil}))



if __name__ == "__main__":
    _main()
