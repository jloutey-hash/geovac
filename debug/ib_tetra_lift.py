"""PF (record IB follow-up): lift the off-slice tetrahedron system (c, A) to Q from primes, check it on a held-out
prime, build the operator for S = Y/t^3 in x = t^2 exactly over Q(x), check it against each prime's operator, and
write it in the Q-form of the ga3/pf scripts: sum_i t^(2i) Q_i(u) on Y, u = theta_t, via Q_i(u) = R_i((u - 3)/2)
where x^order L = sum_i x^i R_i(theta_x) on S.
usage: python ib_tetra_lift.py <out tag> <op json> <op json> ... (last one held out)"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from math import lcm
from pathlib import Path
import sympy as sp
from flint import fmpq_poly
from ib_lift_system import lift_rf, check_rf, cyclic

sys.set_int_max_str_digits(10**7)
DATA = Path(__file__).parent / "data"
FLAGS = [a for a in sys.argv if a.startswith("--")]


def d_to_Q(ip):
    """integer d/dx form on S = Y/t^3 -> Q-form on Y in u = theta_t; returns (list of sympy exprs, x shift)."""
    order = len(ip) - 1
    # Q-form: x^order sum_k p_k(x) D^k = sum_k p_k(x) x^(order-k) ff(theta_x, k)
    th, u = sp.symbols("th u")
    R = {}
    for k, p in enumerate(ip):
        ff = sp.Integer(1)
        for j in range(k):
            ff *= (th - j)
        ffp = sp.Poly(sp.expand(ff), th)
        for i, v in enumerate(p):
            if v:
                R[i + order - k] = R.get(i + order - k, 0) + v * ffp
    lo = min(R)
    Q = []
    for i in range(lo, max(R) + 1):
        Ri = R.get(i, sp.Poly(0, th))
        Q.append(sp.expand(Ri.as_expr().subs(th, (u - 3) / 2) * 2 ** order))
    return Q, lo


def main():
    tag = sys.argv[1]
    ops = [json.loads(Path(a).read_text()) for a in sys.argv[2:] if not a.startswith("--")]
    assert all(o["masters"] == ops[0]["masters"] for o in ops), "master sets differ between primes"
    use, hold = ops[:-1], ops[-1]
    mods = [o["P"] for o in use]
    nm = len(use[0]["c"])
    t0 = time.time()
    A = [[lift_rf([o["A"][i][j] for o in use], mods) for j in range(nm)] for i in range(nm)]
    if any(x is None for row in A for x in row):
        print("system reconstruction failed: more primes needed"); return
    ok = all(check_rf(A[i][j], hold["A"][i][j], hold["P"]) for i in range(nm) for j in range(nm))
    c = [lift_rf([o["c"][i] for o in use], mods) for i in range(nm)]
    if all(x is not None for x in c):
        ok = ok and all(check_rf(c[i], hold["c"][i], hold["P"]) for i in range(nm))
        print("c lifted from the same primes and checked on the held-out one", flush=True)
    elif "--c-all" in FLAGS:
        allm = [o["P"] for o in ops]
        c = [lift_rf([o["c"][i] for o in ops], allm) for i in range(nm)]
        if any(x is None for x in c):
            print("c reconstruction failed from all primes"); return
        hc = max(max(abs(int(q.p)).bit_length(), int(q.q).bit_length()) for x in c
                 for poly in (x.n, x.d) for q in poly.coeffs())
        bound = sum(m.bit_length() for m in allm) // 2
        print(f"c lifted from all {len(ops)} primes, NOT held out (largest height {hc} bits, ratrec bound ~{bound})",
              flush=True)
    else:
        print("c reconstruction failed: more primes needed (or --c-all)"); return
    hb = max(max(abs(int(q.p)).bit_length(), int(q.q).bit_length()) for row in A + [c] for x in row
             for poly in (x.n, x.d) for q in poly.coeffs())
    print(f"system lifted from {len(use)} primes (largest height {hb} bits); held-out prime agrees: {ok}", flush=True)
    if not ok:
        return
    a = cyclic(c, A)
    order = len(a) - 1
    den = fmpq_poly([1])
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [x.n * (den // x.d) for x in a]
    L = 1
    for p in polys:
        for q in p.coeffs():
            L = lcm(L, int(q.q))
    ip = [[int((q * L).p) for q in p.coeffs()] for p in polys]
    from math import gcd
    g = 0
    for p in ip:
        for v in p:
            g = gcd(g, v)
    ip = [[v // g for v in p] for p in ip]
    print("operator over Q: order", order, "degrees", [len(p) - 1 for p in ip], "largest coefficient",
          max(abs(v).bit_length() for p in ip for v in p), "bits", f"{time.time()-t0:.0f}s", flush=True)
    agree = 0
    for o in ops:
        p = o["P"]
        il = pow(ip[-1][-1] % p, -1, p)
        agree += [[v * il % p for v in poly] for poly in ip] == [[v % p for v in poly] for poly in o["polys"]]
    print("agrees with the operator at", agree, "of", len(ops), "primes", flush=True)
    Q, lo = d_to_Q(ip)
    out = {"tag": tag, "form_D": "sum_k polys[k](x) (d/dx)^k (Y/t^3) = 0, x = t^2, integer coefficients",
           "order": order, "polys": [[str(v) for v in p] for p in ip],
           "form_Q": "sum_i t^(2i) Q[i](u) Y = 0, u = theta_t (x-power shifted by %d)" % lo,
           "Q": [str(q) for q in Q]}
    (DATA / f"ib_{tag}_op_Q.json").write_text(json.dumps(out))
    print("wrote", f"ib_{tag}_op_Q.json", "Q-form x-degree", len(Q) - 1)


if __name__ == "__main__":
    main()
