"""Hypothesis after TG (records docs/predictions/2026-10-25_tetrahedron_route_y.md, disposition; PF addendum 4): F_8, the
order-12 factor of Y's operator, contains the rank-4 system INDUCED from Bailey's curve E+ over the conic
w^2 = (x - 1)(x - 9), not the tensor-induced K3 (R_V).  At x = 1 and 9 the induced system swaps E+ and E-, so it has
two eigenvalues -1 there, as F_8 does; R_V has one.

Source: f = g * F(lambda+(x)), F a solution of Legendre's equation lambda(1-lambda)F'' + (1-2 lambda)F' - F/4 = 0,
lambda+ = (x - 3 + w)/(2x) (Bailey's modulus for scales 1, 2, t; debug/geo_y_route_memo.md, Computation 2), g a twist
with g^2 in Q(x)(w): 1, sqrt(lambda+), sqrt(1 - lambda+), sqrt(lambda+ (1 - lambda+)).  f satisfies an order-two
equation over K = Q(x)(w); with its conjugate it gives the order-four operator over Q(x) computed here exactly
(d/dx-form, then theta form).  Checks: the operator annihilates f and its conjugate numerically; its local exponents.
Test: the fixed map search of pf_12_hom.py into F_8 and F_8* (both ends), twelve primes if anything is found.
usage: python pf_13_induced.py <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import pf_lib as A
import pf_12_hom as P

DATA = Path(__file__).parent / "data"
x = sp.Symbol("x")
W = (x - 1) * (x - 9)


class K:
    """a + b w, w^2 = W, a and b rational functions of x"""

    def __init__(self, a, b=0):
        self.a, self.b = sp.cancel(sp.sympify(a)), sp.cancel(sp.sympify(b))

    def __add__(self, o):
        o = o if isinstance(o, K) else K(o)
        return K(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, K) else K(o)))

    def __rsub__(self, o):
        return K(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, K) else K(o)
        return K(self.a * o.a + self.b * o.b * W, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inv(self):
        d = sp.cancel(self.a ** 2 - self.b ** 2 * W)
        return K(self.a / d, -self.b / d)

    def __truediv__(self, o):
        o = o if isinstance(o, K) else K(o)
        return self * o.inv()

    def d(self):
        return K(sp.diff(self.a, x), sp.diff(self.b, x) + self.b * (x - 5) / W)


def induced_op(twist):
    lam = K((x - 3) / (2 * x), 1 / (2 * x))
    l1, l2 = lam.d(), lam.d().d()
    one = K(1)
    p = l2 / l1 - (one - 2 * lam) * l1 / (lam * (one - lam))
    q = l1 * l1 / (4 * lam * (one - lam))
    if twist == "1":
        ell = K(0)
    elif twist == "sqrt(lam)":
        ell = l1 / (2 * lam)
    elif twist == "sqrt(1-lam)":
        ell = -l1 / (2 * (one - lam))
    elif twist == "sqrt(lam(1-lam))":
        ell = l1 / (2 * lam) - l1 / (2 * (one - lam))
    p2 = p + 2 * ell
    q2 = q - p * ell - ell * ell + ell.d()
    # d^k f = alpha_k f + beta_k f'
    al, be = [K(1)], [K(0)]
    for k in range(4):
        a, b = al[-1], be[-1]
        al.append(a.d() + b * q2)
        be.append(a + b.d() + b * p2)
    M = sp.Matrix([[al[k].a for k in range(5)], [al[k].b for k in range(5)],
                   [be[k].a for k in range(5)], [be[k].b for k in range(5)]])
    ns = M.nullspace(simplify=True)
    assert len(ns) == 1, f"nullspace {len(ns)}"
    v = [sp.cancel(c) for c in ns[0]]
    den = sp.lcm([sp.fraction(sp.together(c))[1] for c in v])
    v = [sp.expand(sp.cancel(c * den)) for c in v]
    g = sp.gcd_list(v)
    return [sp.expand(sp.cancel(c / g)) for c in v], (p2, q2, ell)


def numeric_check(b, twist, x0="12"):
    mp.mp.dps = 40
    X0 = mp.mpf(x0)
    out = []
    for sgn in (1, -1):
        def f(xx, sgn=sgn):
            w = mp.sqrt((xx - 1) * (xx - 9))
            lam = (xx - 3 + sgn * w) / (2 * xx)
            gv = {"1": 1, "sqrt(lam)": mp.sqrt(lam), "sqrt(1-lam)": mp.sqrt(1 - lam),
                  "sqrt(lam(1-lam))": mp.sqrt(lam * (1 - lam))}[twist]
            return gv * mp.hyp2f1(0.5, 0.5, 1, lam)
        ds = [mp.diff(f, X0, k) for k in range(len(b))]
        terms = [mp.mpf(sp.N(bj.subs(x, sp.Rational(x0)), 50)) * ds[j] for j, bj in enumerate(b)]
        out.append(abs(mp.fsum(terms)) / max(abs(v) for v in terms))
    return out


def main():
    tag = sys.argv[1]
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    out = {}
    t0 = time.time()
    for twist in ("1", "sqrt(lam)", "sqrt(1-lam)", "sqrt(lam(1-lam))"):
        b, _ = induced_op(twist)
        chk = numeric_check(b, twist)
        S = P.theta_from_dform(b)
        rep = C.local_report(S, A.XPTS, verbose=False)
        loc = {"t0": rep["t0"], "inf": rep["inf_mu"],
               **{k: v["exps"] for k, v in rep["points"].items() if v["exps"] != [str(i) for i in range(4)]},
               "lead": [(f, e) for f, e in rep["lead"]]}
        print(f"twist {twist}: order {len(S)-1}, x-degree {G.op_degree(S)}; numeric |S f|/max term (both conjugates) "
              f"{[mp.nstr(c, 3) for c in chk]}; local {loc} ({time.time()-t0:.0f}s)", flush=True)
        r = P.hom(S, F8, 8, 6, 300)
        l = P.hom(G.adjoint(S), G.adjoint(F8), 8, 6, 300)
        print(f"   Hom(S, F_8) {({k: v['nullity_two_primes'] for k, v in r.items() if k.startswith('nu')})}; "
              f"Hom(S*, F_8*) {({k: v['nullity_two_primes'] for k, v in l.items() if k.startswith('nu')})} ({time.time()-t0:.0f}s)", flush=True)
        out[twist] = {"dform": [str(c) for c in b], "theta_form": G.op_to_strings(S), "numeric_check": [mp.nstr(c, 3) for c in chk],
                      "local": loc, "Hom(S,F_8)": r, "Hom(S*,F_8*)": l}
        (DATA / f"pf_13_induced_{tag}.json").write_text(json.dumps(out, indent=1, default=str))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
