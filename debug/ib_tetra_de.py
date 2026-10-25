"""IB3: Y's differential equation in x = t^2 from IBP on the tetrahedron, mod P; operator for Y/t^3.
usage: python ib_tetra_de.py <tag> <a> <b> <g1> <g2> [R S nsamp workers]"""
from __future__ import annotations
import json, random, sys, time
from fractions import Fraction as F
from multiprocessing import Pool
from pathlib import Path
from flint import nmod_poly, nmod
from ib_reduce import P, Reducer, shifts_box, fr
from ib_tetra import tetra_family, TARGET, dx_terms
from ib_kite_de import fit_rational, RF, cyclic_operator

DATA = Path(__file__).parent / "data"


def sample(args):
    xi, masses, R_, S_, sym = args
    a, b, g1, g2 = masses
    fam = tetra_family(a, b, g1, g2)
    fam.setup({}, 3, masses2=[fr(a * a), fr(b * b), xi, fr(g1 * g1), fr(g2 * g2), 0])
    seeds = [n for n in shifts_box(6, R_, S_, fam.integer) if not fam.is_zero(n)]
    Rd = Reducer(fam)
    Rd.seed(seeds, sym=(lambda n: (n[1], n[0], n[2], n[4], n[3], n[5])) if sym else None)
    c = Rd.reduce({TARGET: 1})
    masters, rows, todo = set(c), {}, list(c)
    while todo:
        m = todo.pop()
        r = Rd.reduce(dx_terms(fam, m))
        rows[m] = r
        for n in r:
            if n not in masters:
                masters.add(n); todo.append(n)
    ms = sorted(masters, key=Rd.key)
    return xi, [str(m) for m in ms], {str(m): c.get(m, 0) for m in ms}, \
        {str(m): {str(n): v for n, v in rows[m].items()} for m in ms}


def main():
    tag = sys.argv[1]
    sys.argv, flags = [a for a in sys.argv if not a.startswith("--")], [a for a in sys.argv if a.startswith("--")]
    sys.argv += flags
    masses = [F(v) for v in sys.argv[2:6]]
    R_ = int(sys.argv[6]) if len(sys.argv) > 6 else 5
    S_ = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    nsamp = int(sys.argv[8]) if len(sys.argv) > 8 else 80
    workers = int(sys.argv[9]) if len(sys.argv) > 9 else 8
    sym = masses[0] == masses[1] and masses[2] == masses[3]
    rng = random.Random(4242)
    xs_all = [rng.randrange(10**6, P - 1) for _ in range(nsamp)]
    t0 = time.time()
    saved = DATA / f"ib_{tag}_samples.json"
    if "--from-samples" in sys.argv:
        res = json.loads(saved.read_text())
        print(f"loaded {len(res)} saved samples", flush=True)
    else:
        with Pool(workers, maxtasksperchild=1) as pool:
            res = pool.map(sample, [(x, masses, R_, S_, sym) for x in xs_all], chunksize=1)
        print(f"{nsamp} samples in {time.time()-t0:.0f}s", flush=True)
        saved.write_text(json.dumps(res))
    bases = {}
    for r in res:
        bases[tuple(r[1])] = bases.get(tuple(r[1]), 0) + 1
    best = max(bases, key=bases.get)
    if len(bases) > 1:
        print("WARNING: basis differs at", nsamp - bases[best], "points")
    res = [r for r in res if tuple(r[1]) == best]
    ms = list(best)
    print("masters", len(ms), ms, flush=True)
    fitp, testp = res[:-10], res[-10:]
    xs = [r[0] for r in fitp]

    def rebuild(get):
        ys = [get(r) for r in fitp]
        if not any(ys):
            return RF(nmod_poly([0], P)), 0
        f = fit_rational(xs, ys, len(xs) - 3)
        if f is None:
            raise RuntimeError("rational fit failed: too few samples")
        N, D = f
        for r in testp:
            assert (N(r[0]) - nmod(get(r), P) * D(r[0])) == 0, "held-out point fails"
        return RF(N, D), N.degree() + D.degree()

    c, A, maxdeg = [], [], 0
    for m in ms:
        f, d = rebuild(lambda r, m=m: r[2][m]); c.append(f); maxdeg = max(maxdeg, d)
    for m in ms:
        row = []
        for n in ms:
            f, d = rebuild(lambda r, m=m, n=n: r[3][m].get(n, 0)); row.append(f); maxdeg = max(maxdeg, d)
        A.append(row)
    print("rebuilt; max total degree", maxdeg, "held-out points pass", flush=True)
    a = cyclic_operator(c, A)
    den = nmod_poly([1], P)
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [[int(v) for v in (x.n * (den // x.d)).coeffs()] for x in a]
    out = {"tag": tag, "masses": [str(m) for m in masses], "P": P, "masters": ms, "order": len(a) - 1,
           "form": "sum_k polys[k](x) (d/dx)^k (Y/t^3) = 0, x = t^2, mod P", "polys": polys,
           "c": [[[int(v) for v in x.n.coeffs()], [int(v) for v in x.d.coeffs()]] for x in c],
           "A": [[[[int(v) for v in x.n.coeffs()], [int(v) for v in x.d.coeffs()]] for x in row] for row in A],
           "seconds": time.time() - t0}
    (DATA / f"ib_{tag}_op.json").write_text(json.dumps(out))
    print("order", len(a) - 1, "degrees", [len(p) - 1 for p in polys], f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
