"""QG step 7b: exact rational gauge map from a reference family into an operator's solution space.
Reference: f3 (L3 = 4 th^3 - t^2 (th+1)^3, the K3 rung) or f2 (Legendre th^2 - t^2 (th+1)^2).
Side R: A with L (A f) = 0 (reference = sub-module); side L: the same on the adjoint (quotient).
Exact rational linear algebra on the reference t-series, then the operator identity L o A = 0 modulo
the reference operator, in exact arithmetic.
usage: python ring_ext_QG_7b_map.py <op-json> <key> <f3|f2> <R|L> <den> <Dg> <a> <out-tag>
  op-json holds a list of theta-polynomials (sum_i t^(2i) Q_i(u)) under <key>."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
import ring_ext_QG_7_k3map as K
from ring_ext_QG_ops import u, tt as t

DATA = Path(__file__).parent / "data"


def ref_series(name, N):
    out, term = {}, Fraction(1)
    for n in range(N // 2):
        out[2 * n] = term
        if name == "f3":
            term = term * Fraction(2 * n + 1, 2) ** 3 / Fraction(n + 1) ** 3 / 4
        else:
            term = term * Fraction(2 * n + 1, 2) ** 2 / Fraction(n + 1) ** 2
    return out


def main():
    path, key, ref, side, dens, Dg, a, tag = sys.argv[1:9]
    Dg, a = int(Dg), int(a)
    Q = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(path).read_text())[key]]
    den = sp.sympify(dens, locals={"t": t})
    if side == "L":
        Q = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
    rr = 3 if ref == "f3" else 2
    NT = max(200, 2 * (3 * (Dg + a + 1)) + 4 * len(Q) + 80)
    w = K.inv_poly_series(sp.expand(den), NT)
    f = ref_series(ref, NT)
    top = NT - 2 * len(Q) - Dg - 6
    cols, labels = [], []
    t0 = time.time()
    for j in range(rr):
        gj = K.mul(w, K.theta_s(f, j), NT)
        for dd in range(-a, Dg + 1):
            cols.append(K.apply_Q(Q, {e + dd: v for e, v in gj.items()}))
            labels.append((j, dd))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[sp.Rational(col.get(e, 0).numerator, col.get(e, 0).denominator) if col.get(e, 0) else sp.Integer(0)
             for col in cols] for e in exps]
    M = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    ns = M.nullspace().to_Matrix()
    print("exact nullity", ns.shape[0], f"({time.time()-t0:.1f}s, {len(rows)}x{len(cols)})", flush=True)
    if ns.shape[0] == 0:
        return
    A = [sp.Integer(0)] * rr
    for v, (j, dd) in zip(list(ns.row(0)), labels):
        A[j] += v * t**dd
    A = [sp.cancel(sp.together(x / den)) for x in A]
    Lref = O.Q_to_coeffs([4 * u**3, -(u + 1) ** 3] if ref == "f3" else [u**2, -(u + 1) ** 2])
    S = [sp.cancel(-c / Lref[rr]) for c in Lref[:rr]]

    def th_vec(p):
        out = []
        for i in range(rr):
            v = O.th_fn(p[i]) + p[rr - 1] * S[i]
            if i > 0:
                v += p[i - 1]
            out.append(sp.cancel(v))
        return out
    Lc = O.Q_to_coeffs(Q)
    acc, cur = [0] * rr, A
    for m in Lc:
        acc = [sp.cancel(x + m * y) for x, y in zip(acc, cur)]
        cur = th_vec(cur)
    ok = all(sp.simplify(x) == 0 for x in acc)
    print("operator identity L o A = 0 mod reference:", ok, flush=True)
    (DATA / f"ring_ext_QG_7b_map_{tag}.json").write_text(json.dumps(
        {"op": path, "key": key, "ref": ref, "side": side, "den": dens, "Dg": Dg, "a": a,
         "nullity": ns.shape[0], "A": [str(x) for x in A], "identity_exact": ok}, indent=1))


if __name__ == "__main__":
    main()
