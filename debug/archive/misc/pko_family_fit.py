"""PK record, stage 2b: reconstruct the family system (c, A) as rational functions of (z, s), s = g2^2, and use it.

  python pko_family_fit.py fit <grid json> <out json> <nzfit> <nzhold> <nshold>
      one prime: per s, each entry as N(z)/D(z) with D monic (ib_kite_de.fit_rational, held-out z); then every
      z-coefficient as a rational function of s with monic denominator (held-out s).  Writes the residues.
  python pko_family_fit.py lift <out json> <fit json> <fit json> ... [--check <grid json>]
      CRT + rational reconstruction of every s-coefficient over the primes; optional check against raw samples of
      a further prime (each sample's entries from the lifted family at its (z, s)).
  python pko_family_fit.py special <family json> <s rational> <op json out> [--compare <ib op json> ...]
      specialize at s, rebuild the system over Q(z), compare with saved systems mod their primes, build the
      operator by the cyclic vector over Q(z), time it.
  python pko_family_fit.py points <family json> <points json>
      compare the family with direct samples at their (z, s), mod the samples' prime.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import lcm
from pathlib import Path

sys.set_int_max_str_digits(10**6)


def fit_one_prime(grid, nzfit, nzhold, nshold):
    import os
    os.environ["IB_PRIME"] = str(grid["P"])
    from ib_kite_de import fit_rational
    from flint import nmod_poly, nmod
    P = grid["P"]
    by_s = defaultdict(list)
    for r in grid["samples"]:
        by_s[r[1]].append(r)
    ms = grid["samples"][0][2]
    keys = [("c", m, None) for m in ms] + [("A", m, n) for m in ms for n in ms]

    def get(r, k):
        return r[3][k[1]] if k[0] == "c" else r[4][k[1]].get(k[2], 0)

    zfits = {}   # s -> key -> (num coeffs, den coeffs) monic den
    for s, rs in by_s.items():
        assert all(r[2] == ms for r in rs), "master basis differs"
        fitp, testp = rs[:nzfit], rs[nzfit:nzfit + nzhold]
        xs = [r[0] for r in fitp]
        out = {}
        for k in keys:
            ys = [get(r, k) for r in fitp]
            if all(y == 0 for y in ys) and all(get(r, k) == 0 for r in testp):
                out[k] = ([], [1])
                continue
            f = fit_rational(xs, ys, len(xs) - 3)
            assert f is not None, f"z-fit failed at s={s} entry {k}"
            N, D = f
            assert all((N(r[0]) - nmod(get(r, k), P) * D(r[0])) == 0 for r in testp), "held-out z fails"
            lc = D.coeffs()[-1]
            il = nmod(1, P) / lc
            N, D = N * il, D * il
            out[k] = ([int(v) for v in N.coeffs()], [int(v) for v in D.coeffs()])
        zfits[s] = out
    # shape per entry: the most common (deg N, deg D); drop s values where an entry has a lower shape
    shape = {}
    for k in keys:
        cnt = Counter((len(zfits[s][k][0]), len(zfits[s][k][1])) for s in zfits)
        shape[k] = max(cnt, key=lambda t: (t[0] + t[1], cnt[t]))
    good = [s for s in zfits if all((len(zfits[s][k][0]), len(zfits[s][k][1])) == shape[k] for k in keys)]
    print(f"prime {P}: {len(zfits)} s values, {len(good)} with the generic z-shape", flush=True)
    sfit, shold = good[:len(good) - nshold], good[len(good) - nshold:]
    res = {}
    maxdeg = 0
    for k in keys:
        dn, dd = shape[k]
        if dn == 0:
            continue
        ent = {"num": [], "den": []}
        for part, ln in (("num", dn), ("den", dd)):
            for i in range(ln):
                ys = [zfits[s][k][0 if part == "num" else 1][i] for s in sfit]
                if all(y == ys[0] for y in ys) and all(zfits[s][k][0 if part == "num" else 1][i] == ys[0] for s in shold):
                    ent[part].append(([ys[0]], [1]))
                    continue
                f = fit_rational(sfit, ys, len(sfit) - 3)
                assert f is not None, f"s-fit failed, entry {k} {part} {i}"
                N, D = f
                for s in shold:
                    assert (N(s) - nmod(zfits[s][k][0 if part == "num" else 1][i], P) * D(s)) == 0, "held-out s fails"
                lc = D.coeffs()[-1]
                il = nmod(1, P) / lc
                N, D = N * il, D * il
                maxdeg = max(maxdeg, N.degree() + D.degree())
                ent[part].append(([int(v) for v in N.coeffs()], [int(v) for v in D.coeffs()]))
        res["|".join(str(x) for x in k)] = ent
    print(f"prime {P}: max s-degree (num + den) of a z-coefficient: {maxdeg}", flush=True)
    return {"P": P, "masters": ms, "entries": res, "max_s_degree": maxdeg, "n_s_used": len(sfit),
            "n_s_held": len(shold)}


def crt2(res, mods):
    x, M = 0, 1
    for r, m in zip(res, mods):
        t = ((r - x) * pow(M, -1, m)) % m
        x, M = x + M * t, M * m
    return x, M


def ratrec(a, M):
    a %= M
    r0, r1, s0, s1 = M, a, 0, 1
    bound = int((M // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return F(r1, s1)


def lift(fits):
    mods = [f["P"] for f in fits]
    out = {}
    hb = 0
    for k in fits[0]["entries"]:
        ent = {}
        for part in ("num", "den"):
            lst = []
            for i in range(len(fits[0]["entries"][k][part])):
                pair = []
                for j in (0, 1):
                    ln = max(len(f["entries"][k][part][i][j]) for f in fits)
                    co = []
                    for t in range(ln):
                        res = [(f["entries"][k][part][i][j][t] if t < len(f["entries"][k][part][i][j]) else 0) for f in fits]
                        x, M = crt2(res, mods)
                        q = ratrec(x, M)
                        if q is None:
                            return None, None
                        hb = max(hb, abs(q.numerator).bit_length(), q.denominator.bit_length())
                        co.append(q)
                    pair.append(co)
                lst.append(pair)
            ent[part] = lst
        out[k] = ent
    return out, hb


def ev_poly(co, x):
    acc = 0
    for c in reversed(co):
        acc = acc * x + c
    return acc


def entry_at(ent, s):
    """the z-rational function (num coeffs, den coeffs) in Q at rational s"""
    num = [ev_poly(n, s) / ev_poly(d, s) for n, d in ent["num"]]
    den = [ev_poly(n, s) / ev_poly(d, s) for n, d in ent["den"]]
    return num, den


def to_str(fam):
    return {k: {p: [[[str(c) for c in a], [str(c) for c in b]] for a, b in v[p]] for p in v} for k, v in fam.items()}


def from_str(d):
    return {k: {p: [[[F(c) for c in a], [F(c) for c in b]] for a, b in v[p]] for p in v} for k, v in d.items()}


def main():
    mode = sys.argv[1]
    if mode == "fit":
        grid = json.loads(Path(sys.argv[2]).read_text())
        out = Path(sys.argv[3])
        assert not out.exists()
        t0 = time.time()
        r = fit_one_prime(grid, int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
        r["seconds"] = time.time() - t0
        out.write_text(json.dumps(r))
    elif mode == "lift":
        out = Path(sys.argv[2])
        assert not out.exists()
        args = sys.argv[3:]
        check = None
        if "--check" in args:
            i = args.index("--check")
            check = args[i + 1]
            args = args[:i]
        fits = [json.loads(Path(a).read_text()) for a in args]
        fam, hb = lift(fits)
        if fam is None:
            print("lift failed: more primes needed")
            return
        print(f"family lifted from {len(fits)} primes, largest height {hb} bits", flush=True)
        if check:
            g = json.loads(Path(check).read_text())
            p = g["P"]
            bad = 0
            tot = 0
            for r in g["samples"]:
                z, s = r[0], r[1]
                for k, ent in fam.items():
                    typ, m, n = k.split("|")
                    val = r[3][m] if typ == "c" else r[4][m].get(n, 0)

                    def red(co):
                        acc = 0
                        for c in reversed(co):
                            acc = (acc * s + c.numerator * pow(c.denominator, -1, p)) % p
                        return acc
                    numz = [red(a) * pow(red(b), -1, p) % p for a, b in ent["num"]]
                    denz = [red(a) * pow(red(b), -1, p) % p for a, b in ent["den"]]
                    N = sum(c * pow(z, i, p) for i, c in enumerate(numz)) % p
                    D = sum(c * pow(z, i, p) for i, c in enumerate(denz)) % p
                    tot += 1
                    bad += (N - val * D) % p != 0
            print(f"check prime {p}: {tot - bad} of {tot} entry values agree", flush=True)
        out.write_text(json.dumps({"primes": [f["P"] for f in fits], "height_bits": hb, "masters": fits[0]["masters"],
                                   "entries": to_str(fam)}))
    elif mode == "special":
        from flint import fmpq_poly, fmpq, nmod_poly
        from ib_lift_system import RQ, cyclic
        fam = json.loads(Path(sys.argv[2]).read_text())
        ms = fam["masters"]
        ents = from_str(fam["entries"])
        s = F(sys.argv[3])
        out = Path(sys.argv[4])
        assert not out.exists()
        rest = sys.argv[5:]
        compare = rest[rest.index("--compare") + 1:] if "--compare" in rest else []
        t0 = time.time()
        zero = RQ(fmpq_poly([0]))

        def rq(k):
            if k not in ents:
                return zero
            num, den = entry_at(ents[k], s)
            return RQ(fmpq_poly([fmpq(q.numerator, q.denominator) for q in num]),
                      fmpq_poly([fmpq(q.numerator, q.denominator) for q in den]))
        c = [rq(f"c|{m}|None") for m in ms]
        A = [[rq(f"A|{m}|{n}") for n in ms] for m in ms]
        t1 = time.time()
        for path in compare:
            o = json.loads(Path(path).read_text())
            p = o["P"]
            assert o["masters"] == ms

            def red(poly):
                return [(int(q.p) * pow(int(q.q), -1, p)) % p for q in poly.coeffs()]
            ok = 0
            tot = 0
            for i, m in enumerate(ms):
                tot += 1
                ok += red(c[i].n) == [v % p for v in o["c"][i][0]] and red(c[i].d) == [v % p for v in o["c"][i][1]]
                for j in range(len(ms)):
                    tot += 1
                    ok += red(A[i][j].n) == [v % p for v in o["A"][i][j][0]] and red(A[i][j].d) == [v % p for v in o["A"][i][j][1]]
            print(f"specialization at s = {s} vs saved system mod {p}: {ok} of {tot} entries identical", flush=True)
        t2 = time.time()
        a = cyclic(c, A)
        den = fmpq_poly([1])
        for x in a:
            den = den * x.d // den.gcd(x.d)
        polys = [x.n * (den // x.d) for x in a]
        Lc = 1
        for pp in polys:
            for q in pp.coeffs():
                Lc = lcm(Lc, int(q.q))
        ipolys = [[int((q * Lc).p) for q in pp.coeffs()] for pp in polys]
        t3 = time.time()
        print(f"s = {s}: specialization {t1-t0:.1f} s, cyclic vector {t3-t2:.1f} s; operator order {len(a)-1}, "
              f"degrees {[len(pp)-1 for pp in ipolys]}", flush=True)
        out.write_text(json.dumps({"form": "sum_k polys[k](z) (d/dz)^k phi = 0, integer coefficients",
                                   "s": str(s), "polys": [[str(v) for v in pp] for pp in ipolys],
                                   "seconds_specialize": t1 - t0, "seconds_cyclic": t3 - t2}))
    elif mode == "points":
        fam = json.loads(Path(sys.argv[2]).read_text())
        ents = from_str(fam["entries"])
        g = json.loads(Path(sys.argv[3]).read_text())
        p = g["P"]
        bad = tot = 0
        for r in g["samples"]:
            z, si = r[0], r[1]
            for k, ent in ents.items():
                typ, m, n = k.split("|")
                val = r[3][m] if typ == "c" else r[4][m].get(n, 0)

                def red(co):
                    acc = 0
                    for cc in reversed(co):
                        acc = (acc * si + cc.numerator * pow(cc.denominator, -1, p)) % p
                    return acc
                numz = [red(a) * pow(red(b), -1, p) % p for a, b in ent["num"]]
                denz = [red(a) * pow(red(b), -1, p) % p for a, b in ent["den"]]
                N = sum(cc * pow(z, i, p) for i, cc in enumerate(numz)) % p
                D = sum(cc * pow(z, i, p) for i, cc in enumerate(denz)) % p
                tot += 1
                bad += (N - val * D) % p != 0
        print(f"direct samples mod {p}: {tot - bad} of {tot} entry values agree with the family")


if __name__ == "__main__":
    main()
