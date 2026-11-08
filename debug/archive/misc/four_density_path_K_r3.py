"""FP step K: the exponent-0 log-free solution y0 of W11 has an order-3 minimal annihilator (step J, mod p,
t-degree 14).  Exact R3 from the exact series (nullspace over QQ), exact right division W11 = M8 o R3."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_G_peel1 as G  # noqa: E402
import four_density_star_H_order4 as H  # noqa: E402
import four_density_star_I_w9 as I9  # noqa: E402
import four_density_star_ops as O  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_K.json"


def main():
    rec = {}
    t0 = time.time()
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_path_I.json").read_text())["quotient_Q"]]
    y0 = G.frob_exact(Q, 0, 2 * 70)
    print(f"y0: {len(y0)} terms; first {[str(y0[e]) for e in sorted(y0)[:5]]} ({time.time()-t0:.0f}s)", flush=True)
    sols = H.guess(y0, 3, 7, 60)
    print("nullity order 3, t-degree 14:", len(sols), flush=True)
    R3Q = I9.primitive(sols[0])
    while R3Q[-1] == 0:
        R3Q.pop()
    print("R3:", [sp.factor(q) for q in R3Q], flush=True)
    W11 = O.from_Q(Q)
    R3 = O.from_Q(R3Q)
    M8, rem = O.rdiv(W11, R3)
    rec["W11_div_R3_zero"] = all(c == 0 for c in rem)
    print("W11 = M8 o R3 exact:", rec["W11_div_R3_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    rec["R3_Q"] = [str(q) for q in R3Q]
    lead = sp.factor(sum(sp.Poly(q, u).coeff_monomial(u**3) * t**(2*i) for i, q in enumerate(R3Q)))
    rec["R3_lead"] = str(lead)
    print("R3 leading coefficient:", lead, flush=True)
    M8Q, D = O.to_Q(M8)
    rec["M8_Q"] = [str(q) for q in I9.primitive(M8Q)]
    print("M8 order", len(M8) - 1, "t-terms", len(M8Q), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
