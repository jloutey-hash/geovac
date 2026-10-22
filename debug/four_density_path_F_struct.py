"""FP step F: structure of W12 (step D): leading coefficient, exponents at t = 0, simple right factors
(theta - 2, theta - 3, E2 o (theta - 3) = the one-leaf moment's operator)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_F.json"


def load(key="W12"):
    d = json.loads((HERE / "data" / "four_density_path_D.json").read_text())
    return [sp.sympify(q).subs(sp.Symbol("u"), u) for q in d[key]["Q"]]


def main():
    rec = {}
    Q = load()
    r = max(sp.Poly(q, u).degree() for q in Q)
    lead = sp.factor(sum(sp.Poly(q, u).coeff_monomial(u**r) * t**(2*i) for i, q in enumerate(Q)))
    print("order", r, "t-degree", 2*(len(Q)-1), "leading coefficient:", lead, flush=True)
    rec["lead"] = str(lead)
    print("Q_0:", sp.factor(Q[0]), "\nQ_top:", sp.factor(Q[-1]), flush=True)
    for a in [1, 2, 3]:
        rec[f"t^{a} solution"] = all(sp.Poly(q, u).eval(a) == 0 for q in Q)
    print({k: v for k, v in rec.items() if 'solution' in k}, flush=True)
    W = O.from_Q(Q)
    for nm, R in [("theta-2", O.from_Q([u - 2])), ("theta-3", O.from_Q([u - 3])),
                  ("E2o(theta-3)", O.from_Q([(u - 3) ** 3, -u * (u - 3) * (u - 2)]))]:
        _, rem = O.rdiv(W, R)
        rec[f"right factor {nm}"] = all(c == 0 for c in rem)
        print(nm, "right factor of W12:", rec[f"right factor {nm}"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
