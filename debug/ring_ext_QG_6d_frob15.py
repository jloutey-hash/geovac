"""QG step 6d (gamma = 3/2): factor lattice of L9 by Frobenius solutions at t = 0 (mod p): for each
non-resonant exponent, the smallest order k <= 8 of a polynomial-coefficient operator (t-degree <= D)
killing that solution.  A solution with k < 9 spans a right factor of order k."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D

DATA = Path(__file__).parent / "data"
u = sp.Symbol("u")
P = C.P


def min_order(ser, maxorder, Dg, NT):
    for r in range(1, maxorder + 1):
        fs = [D.theta(ser, j) for j in range(r + 1)]
        cols = [D.shift(f, d) for f in fs for d in range(Dg + 1)]
        lo = min(ser)
        exps = [e for e in sorted(set().union(*cols)) if e < lo + NT - 2 * Dg - 4]
        if len(exps) < len(cols) + 30:
            return None
        rows = [[col.get(e, 0) for col in cols] for e in exps[: len(cols) + 60]]
        if len(cols) - C.rank_mod(rows) > 0:
            return r
    return f">{maxorder}"


def main():
    rec = json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())
    Q = [sp.sympify(q, locals={"u": u}) for q in rec["Q"]]
    if len(sys.argv) > 1 and sys.argv[1] == "adj":
        Q = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
    NT = 1400
    out = {}
    roots = sorted(set(int(r) for r in sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    for lam in roots:
        try:
            y = D.frob_t(Q, lam, NT=NT)
        except ZeroDivisionError as e:
            print("exponent", lam, "resonant", e, flush=True)
            continue
        for Dg in (20, 40, 60):
            k = min_order(y, 8, Dg, NT)
            print(f"exponent {lam}: minimal order with t-degree <= {Dg}: {k}", flush=True)
            out[f"{lam}|{Dg}"] = k
    (DATA / ("ring_ext_QG_6d_frob_g1.5_" + (sys.argv[1] if len(sys.argv) > 1 else "") + ".json")).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
