"""NP record, Stage A: z9 along F1 = (a, b, g1) = (1, 2, 3/2), s = g2^2.
  IB_PRIME=p0 python np_stageA.py p0     -> 36 fit + 6 held-out random s, rationals g2 = 5/2, 7/2, 4; fit, lift
  IB_PRIME=p1 python np_stageA.py p1     -> 4 random s + the rationals at p1, checked against the lifted formula
Readout: 4 samples per point; on failure 8 more and a general fit of the carriers (fallback)."""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

from flint import nmod_poly, nmod

from np_engine import P, run, readout, ratrec, TOP, CARRIERS, fr
from ib_kite_de import fit_rational

DATA = Path(__file__).parent / "data"


def _main():
    mode = sys.argv[1]
    out = DATA / f"np_stageA_{mode}.json"
    assert not out.exists()
    rng = random.Random(f"np_stageA|{mode}|{P}")
    BASE3 = [fr(1), fr(4), fr(F(9, 4))]
    RATS = [F(5, 2), F(7, 2), F(4)]


    def pts_for(svals, nper=4):
        args, idx = [], []
        for i, s in enumerate(svals):
            for _ in range(nper):
                args.append((rng.randrange(10**6, P - 1), BASE3 + [s % P, 0], False, "top"))
                idx.append(i)
        return args, idx


    def known_polys(s):
        X = {"a": BASE3[0], "b": BASE3[1], "g1": BASE3[2], "g2": s % P}
        ks = [nmod_poly([0, 1], P), nmod_poly([(X["a"] + X["g2"]) % P, 1], P), nmod_poly([(X["b"] + X["g1"]) % P, 1], P)]
        for u, v in [("a", "b"), ("g1", "g2"), ("a", "g2"), ("b", "g1")]:
            ks.append(nmod_poly([(X[u] + X[v]) % P, 1], P) ** 2 - nmod_poly([4 * X[u] * X[v] % P], P))
        return ks


    def fallback(samps, s):
        """general fit of each carrier from all samples (last 2 held out); strip known factors; residual roots."""
        xs = [r["z"] for r in samps[:-2]]
        L = nmod_poly([1], P)
        for lab in CARRIERS:
            ys = [r["A"][TOP].get(lab, 0) % P for r in samps[:-2]]
            f = fit_rational(xs, ys, len(xs) - 3)
            if f is None:
                return None, "fit failed"
            N, D = f
            for r in samps[-2:]:
                if (N(r["z"]) - nmod(r["A"][TOP].get(lab, 0), P) * D(r["z"])) != 0:
                    return None, "held-out fails"
            L = L * D // L.gcd(D)
        for q in known_polys(s):
            while True:
                qq, rr = divmod(L, q)
                if rr.is_zero():
                    L = qq
                else:
                    break
        fac = L.factor()[1]
        if len(fac) == 1 and fac[0][0].degree() == 1:
            co = [int(c) for c in fac[0][0].coeffs()]
            return (-co[0]) * pow(co[1], P - 2, P) % P, "fallback ok"
        return None, f"residual factors {[(f.degree(), e) for f, e in fac]}"


    def measure(svals):
        args, idx = pts_for(svals)
        res, wall = run(args)
        groups = {}
        for r, i in zip(res, idx):
            groups.setdefault(i, []).append(r)
        vals, notes, nfall = {}, {}, 0
        bad = []
        for i in range(len(svals)):
            g = groups[i]
            assert all(len(r["masters"]) == 18 for r in g), f"master count at point {i}: {[len(r['masters']) for r in g]}"
            r0, ok, note = readout(g)
            if ok:
                vals[i] = r0
            else:
                bad.append(i)
                notes[i] = note
        if bad:
            nfall = len(bad)
            args2, idx2 = pts_for([svals[i] for i in bad], nper=8)
            res2, w2 = run(args2)
            wall += w2
            for j, i in enumerate(bad):
                g = groups[i] + [r for r, k in zip(res2, idx2) if k == j]
                v, note = fallback(g, svals[i])
                notes[i] += " | " + note
                if v is not None:
                    vals[i] = v
        cpu = [r["cpu"] for r in res]
        return vals, notes, nfall, wall, cpu


    def eval_formula(Nc, Dc, s, p):
        num = sum(F(c) * F(s) ** k for k, c in enumerate(Nc))
        den = sum(F(c) * F(s) ** k for k, c in enumerate(Dc))
        return num, den


    if mode == "p0":
        sfit = [rng.randrange(10**6, P - 1) for _ in range(36)]
        shold = [rng.randrange(10**6, P - 1) for _ in range(6)]
        srat = [fr(g * g) for g in RATS]
        allv = sfit + shold + srat
        vals, notes, nfall, wall, cpu = measure(allv)
        print(f"samples wall {wall:.0f}s, fallbacks {nfall}, notes {notes}", flush=True)
        assert all(i in vals for i in range(len(allv))), "readout failed at some point"
        xs, ys = sfit, [vals[i] for i in range(36)]
        f = fit_rational(xs, ys, len(xs) - 4)  # total degree <= 32 (nunk <= 34 = 36 - 2)
        rec = {"P": P, "wall": wall, "cpu": cpu, "fallbacks": nfall, "notes": {str(k): v for k, v in notes.items()},
               "s": allv, "z9": [vals[i] for i in range(len(allv))]}
        if f is None:
            print("NO reconstruction within total degree 32")
            rec["fit"] = None
        else:
            N, D = f
            lc = D.coeffs()[-1]
            il = nmod(1, P) / lc
            N, D = N * il, D * il
            print("fit: num deg", N.degree(), "den deg", D.degree(), flush=True)
            held = [(N(s) - nmod(vals[36 + i], P) * D(s)) == 0 for i, s in enumerate(shold)]
            print("held-out (6) pass:", held)
            Nq = [ratrec(int(c)) for c in N.coeffs()]
            Dq = [ratrec(int(c)) for c in D.coeffs()]
            print("lift:", Nq, Dq)
            rat = []
            for g, i in zip(RATS, range(42, 45)):
                meas = ratrec(vals[i])
                pred = None
                if None not in Nq + Dq:
                    n_, d_ = eval_formula(Nq, Dq, g * g, P)
                    pred = n_ / d_ if d_ != 0 else "pole"
                rat.append((str(g), str(meas), str(pred)))
                print(f"g2 = {g}: measured z9 = {meas}, formula = {pred}")
            rec.update({"fit": [[int(c) for c in N.coeffs()], [int(c) for c in D.coeffs()]], "held": held,
                        "lift": [[str(x) for x in Nq], [str(x) for x in Dq]], "rationals": rat})
            if None not in Nq + Dq:
                n9, d9 = eval_formula(Nq, Dq, 9, P)
                print("formula at s = 9 (g2 = 3): num", n9, "den", d9)
                rec["at9"] = [str(n9), str(d9)]
            # negative control: permuted labels
            perm = ys[:]
            random.Random(99).shuffle(perm)
            fp = fit_rational(xs, perm, len(xs) - 4)
            if fp is None:
                pc = "no reconstruction (control fails as required)"
            else:
                Np, Dp = fp
                hp = [(Np(s) - nmod(vals[36 + i], P) * Dp(s)) == 0 for i, s in enumerate(shold)]
                pc = f"reconstructed deg {Np.degree()},{Dp.degree()}; held-out {hp}"
            print("permuted-label control:", pc)
            rec["perm_control"] = pc
        out.write_text(json.dumps(rec))
    elif mode == "p1":
        a0 = json.loads((DATA / "np_stageA_p0.json").read_text())
        Nq = [F(x) for x in a0["lift"][0]]
        Dq = [F(x) for x in a0["lift"][1]]
        srand = [rng.randrange(10**6, P - 1) for _ in range(4)]
        srat = [fr(g * g) for g in RATS]
        vals, notes, nfall, wall, cpu = measure(srand + srat)
        print(f"samples wall {wall:.0f}s, fallbacks {nfall}, notes {notes}")

        def ev(s):
            n = sum(fr(c) * pow(s, k, P) for k, c in enumerate(Nq)) % P
            d = sum(fr(c) * pow(s, k, P) for k, c in enumerate(Dq)) % P
            return n * pow(d, P - 2, P) % P

        chk = [ev(s) == vals[i] for i, s in enumerate(srand)]
        print("p1 random s (4) agree:", chk)
        rat = [(str(g), str(ratrec(vals[4 + i]))) for i, g in enumerate(RATS)]
        print("p1 rationals:", rat)
        out.write_text(json.dumps({"P": P, "wall": wall, "cpu": cpu, "check": chk, "rationals": rat,
                                   "s": srand + srat, "z9": [vals.get(i) for i in range(7)]}))



if __name__ == "__main__":
    _main()
