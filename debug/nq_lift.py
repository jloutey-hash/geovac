"""Lift the p0 trace fits (nq_fit_<tag>.json) to Q by rational reconstruction (monic denominators), confirm on the
p1 samples (nq_scan_<tag_p1>.jsonl), then test, over Q, divisibility of D_eps by the model quadratic and by the
perturbed (negative-control) quadratic.  Record 2026-10-27_kite_ninth_point_origin, NM1.
usage: python nq_lift.py <tag_p0> <tag_p1> <eps>"""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp

from nq_fit import rr, res

DATA = Path(__file__).parent / "data"
z = sp.Symbol("z")


def quad(e):
    e = sp.Rational(e.numerator, e.denominator)
    return sp.Poly(1296 * e * z**2 + (20484 * e - 1296 * e**2 - 28000) * z + 2700 * e - 2775, z)


if __name__ == "__main__":
    t0, t1, eps = sys.argv[1], sys.argv[2], Fr(sys.argv[3])
    fit = json.load(open(DATA / f"nq_fit_{t0}.json"))
    p0 = fit["p"]
    lifted = []
    for N, D in fit["fits"]:
        Nq = [rr(c, p0) for c in N]
        Dq = [rr(c, p0) for c in D]
        if any(c is None for c in Nq + Dq):
            sys.exit("rational reconstruction failed for some coefficient: lift not available")
        lifted.append((Nq, Dq))
    rows = [json.loads(l) for l in open(DATA / f"nq_scan_{t1}.jsonl")]
    p1 = rows[0]["p"]
    ok = True
    for r in rows:
        for v, (Nq, Dq) in enumerate(lifted):
            n = sum(res(c, p1) * pow(r["z"], k, p1) for k, c in enumerate(Nq)) % p1
            d = sum(res(c, p1) * pow(r["z"], k, p1) for k, c in enumerate(Dq)) % p1
            if (n - r["traces"][v] * d) % p1:
                ok = False
    print(f"lift from p0 confirmed on {len(rows)} p1 samples x 6 traces:", ok)
    if not ok:
        sys.exit("lift not confirmed")
    Ds = [sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * z**k for k, c in enumerate(Dq)), z) for _, Dq in lifted]
    Dl = Ds[0]
    for d in Ds[1:]:
        Dl = sp.lcm(Dl, d)
    print("D_eps over Q, degree", Dl.degree())
    print("factors over Q:", [(sp.factor(f.as_expr()) if f.degree() > 1 else f.as_expr(), m) for f, m in sp.factor_list(Dl)[1]])
    q = quad(eps)
    qn = quad(eps * Fr(1001, 1000))
    print("model quadratic", q.as_expr(), "divides D_eps over Q:", Dl.rem(q).is_zero)
    print("negative control (eps x 1001/1000) divides:", Dl.rem(qn).is_zero)
    print("real roots of D_eps:", [sp.N(r, 15) for r in sp.Poly(Dl, z).real_roots()])
