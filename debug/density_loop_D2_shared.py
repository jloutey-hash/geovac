"""Loop test PT2, step D2: shared composition factors, tested by relations
  sum_{j<ry} q_j(t) theta^j y = sum_{j<rf} p_j(t) theta^j f,   (q, p) both nonzero,
between a cyclic solution y of one operator and a reference f (series at t = 0, mod p).
A relation with both sides nonzero means the two D-modules share a composition factor.
Runs: chain control (y = its t^7 solution, order 7) vs f3 (K3 rung) and vs f2 (elliptic rung);
triangle (y = its t^3 solution, order 3) vs f3 and f2; F_4 climb control (its t^0 regular
solution at t = 0 if non-resonant, else the first non-resonant) vs f3 and f2.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402

u = sp.Symbol("u")
OUT = Path(__file__).parent / "data" / "density_loop_D2.json"


def shared(y, ry, f, rf, Dg):
    ys = [D.theta(y, j) for j in range(ry)]
    fs = [D.theta(f, j) for j in range(rf)]
    cols = [D.shift(s, d) for s in ys for d in range(Dg + 1)]
    ncq = len(cols)
    cols += [{e: (-c) % D.P for e, c in D.shift(s, d).items()} for s in fs for d in range(Dg + 1)]
    lo = min(min(y), min(f))
    hi = lo + D.NT - 2 * Dg - 8
    rows = [[col.get(e, 0) for col in cols] for e in range(lo, hi)]
    tot = len(cols) - C.rank_mod(rows)
    only_q = ncq - C.rank_mod([r[:ncq] for r in rows])
    only_p = (len(cols) - ncq) - C.rank_mod([r[ncq:] for r in rows])
    return tot, only_q, only_p, len(cols), hi - lo


def main() -> None:
    B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
    QG = [sp.sympify(q) for q in B["QG"]]
    QC = [sp.sympify(q) for q in B["QC"]]
    LC = json.loads((Path(__file__).parent / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4 = [sp.sympify(q) for q in LC["4"]["Q"]]
    from fractions import Fraction
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1])
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], scale=Fraction(1, 4))
    rec = {}
    cases = [("chain", D.frob_t(QC, 7), 7), ("triangle", D.frob_t(QG, 3), 3)]
    roots4 = sorted(set(int(r) for r in sp.roots(sp.Poly(Q4[0], u).as_expr(), u)))
    print("F_4 operator t^0 coefficient:", sp.factor(Q4[0]), flush=True)
    for lam in sorted(roots4, reverse=True):
        try:
            cases.append((f"F4(lam={lam})", D.frob_t(Q4, lam), 4))
            break
        except ZeroDivisionError:
            continue
    cases.append(("control f3 vs itself", f3, 3))
    for name, y, ry in cases:
        for ref, f, rf in [("f3", f3, 3), ("f2", f2, 2)]:
            for Dg in (6, 10):
                tot, oq, op, nc, nr = shared(y, ry, f, rf, Dg)
                mixed = tot - oq - op
                print(f"{name:22s} vs {ref}: deg {Dg}: relations {tot} (y-only {oq}, f-only {op}) -> shared {mixed}   [{nc} unknowns, {nr} rows]", flush=True)
                rec[f"{name}|{ref}|{Dg}"] = [tot, oq, op]
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
