"""Loop test PT2, step D12: does the chain's K3 piece split off its elementary pieces?
W = R3 o R1 o (theta - 2) annihilates C (D11).  A split (K3 as a sub, not only a quotient)
would mean a rational gauge B with W(B f) = 0 for f in Sol(L3); registered outcomes (A) and (B)
both need one.  Exact search over denominators t^a (t^2-1)^b (t^2-4)^c (4t^2+11)^e,
a <= 8, b, c <= 3, e <= 2, numerators of degree <= 16.  Positive control: the reversed
composition (elementary) o R3 has the split by construction and must be found."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D9_exact_k3 as X
import density_loop_D10_k3_exact_check as D10
import density_loop_D11_chain_minimal as D11

t, u = D10.t, D10.u


def scan(Q):
    X.M5 = [sp.Poly(q, u) for q in Q]
    X.NT = 240
    hits = []
    for e in range(3):
        for b in range(4):
            for c in range(4):
                ns, _, _ = X.search(b, c, 8, 16, e)
                if ns.shape[0]:
                    hits.append((e, b, c, ns.shape[0]))
    return hits


def main():
    QW, _ = D11.to_Q(D11.W)
    hits = scan(QW)
    print("W (chain):", hits or "no split in the searched bounds", flush=True)
    E2 = D10.compose([-(3 + t**2), 1 - t**2], [sp.Integer(-2), sp.Integer(1)])
    Qs, _ = D11.to_Q(D10.compose(E2, D11.R3))
    X.M5 = [sp.Poly(q, u) for q in Qs]; X.NT = 240
    ctrl = X.search(0, 2, 8, 16)[0].shape[0]
    print("positive control (elementary o R3), (t^2-4)^2:", ctrl, flush=True)
    (Path(__file__).parent / "data" / "density_loop_D12.json").write_text(
        json.dumps({"chain_hits": hits, "control_nullity": ctrl}, indent=1))


if __name__ == "__main__":
    main()
