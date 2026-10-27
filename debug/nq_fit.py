"""Phase-2 analysis (record 2026-10-27_kite_ninth_point_origin): rebuild the six traces T_v(z) at fixed eps from
nq_scan jsonl samples (rational reconstruction by the extended Euclidean algorithm, 6 held-out points), form the
drop polynomial D_eps = lcm of the denominators, and test listed values.
usage: python nq_fit.py <tag_p0> <generic_count> [--p1 <tag_p1>] [--controls z1,z2,...] [--test z1,...]
       [--quad c2,c1,c0]  (test divisibility of D_eps by c2 z^2 + c1 z + c0, coefficients rational)
Control verdicts are printed BEFORE any other root is printed."""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

from flint import nmod_poly

DATA = Path(__file__).parent / "data"
HOLD = 6


def load(tag):
    return [json.loads(l) for l in open(DATA / f"nq_scan_{tag}.jsonl")]


def interp(xs, ys, p):
    # Newton divided differences
    n = len(xs)
    c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % p, -1, p) % p
    P = nmod_poly([c[-1]], p)
    for i in range(n - 2, -1, -1):
        P = P * nmod_poly([(-xs[i]) % p, 1], p) + nmod_poly([c[i]], p)
    return P


def ratrec_poly(xs, ys, hx, hy, p):
    n = len(xs)
    M = nmod_poly([1], p)
    for x in xs:
        M = M * nmod_poly([(-x) % p, 1], p)
    P = interp(xs, ys, p)
    r0, r1, t0, t1 = M, P, nmod_poly([0], p), nmod_poly([1], p)
    cands = [(r1, t1)]
    while not r1.is_zero():
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, t0 - q * t1
        if not r1.is_zero() or True:
            cands.append((r1, t1))
    cands.sort(key=lambda c: max(c[0].degree(), 0) + c[1].degree())
    for N, D in cands:
        if D.is_zero() or N.degree() + D.degree() > n - 3:
            continue
        g = N.gcd(D) if not N.is_zero() else D
        N, D = N // g, D // g
        lc = int(D.leading_coefficient())
        inv = pow(lc, -1, p)
        N, D = N * inv, D * inv
        if any(int(D(x)) == 0 for x in hx):
            continue
        if all((int(N(x)) - y * int(D(x))) % p == 0 for x, y in zip(hx, hy)):
            return N, D
    return None


def rr(a, p):
    """rational reconstruction of a residue, |num|, den < sqrt(p/2)."""
    import math
    bound = math.isqrt(p // 2)
    r0, r1, s0, s1 = p, a % p, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1)


def res(q, p):
    q = Fr(q)
    return q.numerator * pow(q.denominator, -1, p) % p


def fit_tag(tag, generic):
    rows = load(tag)
    p = rows[0]["p"]
    bad = [r for r in rows if r["count"] != generic]
    rows = [r for r in rows if r["count"] == generic]
    print(f"{tag}: {len(rows)} samples with count {generic}, {len(bad)} with other counts {[r['count'] for r in bad]}")
    xs = [r["z"] for r in rows]
    fits = []
    for v in range(6):
        ys = [r["traces"][v] for r in rows]
        f = ratrec_poly(xs[:-HOLD], ys[:-HOLD], xs[-HOLD:], ys[-HOLD:], p)
        if f is None:
            print(f"  trace {v}: NO FIT with {len(xs) - HOLD} points")
            return p, None, None
        print(f"  trace {v}: deg N {f[0].degree()}, deg D {f[1].degree()}, held-out {HOLD}/{HOLD} pass")
        fits.append(f)
    D = nmod_poly([1], p)
    for _, d in fits:
        g = D.gcd(d)
        D = D * (d // g)
    return p, fits, D


if __name__ == "__main__":
    args = sys.argv[1:]
    tag, generic = args[0], int(args[1])
    opts = {args[i]: args[i + 1] for i in range(2, len(args) - 1, 2)}
    p, fits, D = fit_tag(tag, generic)
    if D is None:
        sys.exit("no fit")
    print(f"D_eps: degree {D.degree()} (mod p0)")
    if "--controls" in opts:
        print("CONTROLS (read first):")
        allok = True
        for zs in opts["--controls"].split(","):
            ok = int(D(res(Fr(zs), p))) == 0
            allok &= ok
            print(f"  control {zs}: {'ROOT' if ok else 'NOT a root'}")
        print("CONTROL VERDICT:", "PASS" if allok else "FAIL")
        if not allok:
            sys.exit("controls failed: nothing else read")
    if "--test" in opts:
        for zs in opts["--test"].split(","):
            print(f"  test {zs}: {'ROOT' if int(D(res(Fr(zs), p))) == 0 else 'not a root'}")
    if "--quad" in opts:
        cs = [Fr(c) for c in opts["--quad"].split(",")]
        Qp = nmod_poly([res(cs[2], p), res(cs[1], p), res(cs[0], p)], p)
        print("  quadratic divides D mod p0:", (D % Qp).is_zero())
    lc, facs = D.factor()
    print("factorization of D mod p0 (degree: multiplicity, root as rational if linear):")
    for f, e in sorted(facs, key=lambda fe: fe[0].degree()):
        if f.degree() == 1:
            root = (-int(f.coeffs()[0]) * pow(int(f.coeffs()[1]), -1, p)) % p
            print(f"  linear, power {e}, root {rr(root, p)}")
        else:
            print(f"  degree {f.degree()}, power {e}")
    json.dump({"p": p, "D": [int(c) for c in D.coeffs()],
               "fits": [[[int(c) for c in N.coeffs()], [int(c) for c in Dd.coeffs()]] for N, Dd in fits]},
              open(DATA / f"nq_fit_{tag}.json", "w"))
