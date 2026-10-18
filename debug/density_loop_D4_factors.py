"""Loop test PT2, step D4: is a reference family a right factor (sub) or a left factor (quotient)
of an operator?  L = sum_i t^{2i} Q_i(theta).
  (R) exists A = sum_{j<r} a_j(t) theta^j, a_j polynomials of degree <= D, with L (A f) = 0, A f != 0;
      then A maps Sol(ref) into Sol(L): the reference module is a submodule (right factor).
  (L) the same with the adjoint L* = sum_i t^{2i} Q_i(-theta - 2i): the reference is a quotient.
References f2 (Legendre, two-density elliptic) and f3 (F_3, K3) are self-adjoint up to sign, so the
same series serve for (L).  Series at t = 0 mod p.  Controls: L = the reference operator itself
must give (R) and (L); F_4's operator must give neither.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402

u = sp.Symbol("u")
P = C.P
NT = 520
OUT = Path(__file__).parent / "data" / "density_loop_D4.json"


def apply_L(Q, ser):
    Qp = [sp.Poly(q, u) for q in Q]
    out = {}
    for i, qp in enumerate(Qp):
        cs = [(m[0], D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)))) for m, c in zip(qp.monoms(), qp.coeffs())]
        for e, c in ser.items():
            val = sum(cc * pow(e % P, m, P) for m, cc in cs) % P
            if val:
                out[e + 2 * i] = (out.get(e + 2 * i, 0) + val * c) % P
    return {e: c for e, c in out.items() if c}


def adjoint(Q):
    return [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]


def factor_test(Q, f, r, Dg):
    cols = [apply_L(Q, D.shift(D.theta(f, j), d)) for j in range(r) for d in range(Dg + 1)]
    top = min(f) + NT - 2 * len(Q) - Dg - 6
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows), len(cols), len(rows)


def main() -> None:
    B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
    QG = [sp.sympify(q) for q in B["QG"]]
    QC = [sp.sympify(q) for q in B["QC"]]
    QC6 = [sp.cancel(q / (u - 2)) for q in QC]          # right factor theta - 2 removed
    LCd = json.loads((Path(__file__).parent / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4 = [sp.sympify(q) for q in LCd["4"]["Q"]]
    Q3 = [4 * u**3, -(u + 1) ** 3]
    Q2 = [u**2, -(u + 1) ** 2]
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    assert not apply_L(Q3, f3) or min(apply_L(Q3, f3)) > NT - 10, "f3 is not a solution of F_3's operator"
    assert not apply_L(Q2, f2) or min(apply_L(Q2, f2)) > NT - 10, "f2 is not a solution of Legendre"
    rec = {}
    ops = [("control F_3 op", Q3), ("control Legendre", Q2), ("control F_4 op", Q4),
           ("triangle L_G", QG), ("chain L_C", QC), ("chain L_C/(th-2)", QC6)]
    for name, Q in ops:
        for ref, f, r in [("f3", f3, 3), ("f2", f2, 2)]:
            for side, QQ in [("R", Q), ("L", adjoint(Q))]:
                for Dg in (8, 16):
                    n, nc, nr = factor_test(QQ, f, r, Dg)
                    rec[f"{name}|{ref}|{side}|{Dg}"] = n
                    print(f"{name:20s} {ref} ({side}) deg {Dg:2d}: {n}  [{nc} unknowns, {nr} rows]", flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
