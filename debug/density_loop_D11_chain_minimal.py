"""Loop test PT2, step D11: is the chain control's top piece N2 part of the chain's D-module?
D10 gives L_C = N2 o R3 o R1 o (theta - 2) exactly (R1 = (theta-3) - t^2 (theta+1)), with N2 of
order 2, irregular at infinity.  Evaluate W C for the order-5 right part W = R3 o R1 o (theta-2)
on the chain integral (step B2's evaluator).  W C = 0 means N2 is spurious: the chain is
annihilated by K3 (R3, equivalent to F_3's operator) above two elementary pieces."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_B2_inhom as B2
import density_loop_D10_k3_exact_check as X

t, u = X.t, X.u


def to_Q(L):
    """function-left theta-coefficients (even in t, rational) -> polynomial t^{2i} Q_i(u) form,
    after clearing denominators on the left."""
    den = sp.lcm([sp.denom(sp.together(c)) for c in L])
    cs = [sp.expand(sp.cancel(c * den)) for c in L]
    z = sp.Symbol("z")
    Q = {}
    for k, c in enumerate(cs):
        pc = sp.Poly(c, t)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            assert e % 2 == 0
            Q[e // 2] = Q.get(e // 2, 0) + v * u**k
    return [sp.expand(Q.get(i, 0)) for i in range(max(Q) + 1)], den


R3 = [sp.sympify(x) for x in json.loads((Path(__file__).parent / "data" / "density_loop_D10.json").read_text())["R3"]]
R1 = [-(3 + t**2), 1 - t**2]
W = X.compose(X.compose(R3, R1), [sp.Integer(-2), sp.Integer(1)])
QW, den = to_Q(W)
# control: the same composition with R3's theta^2 coefficient 198 -> 199 must not annihilate C
R3p = [R3[0], R3[1], sp.cancel(R3[2] * sp.Rational(199, 198)), R3[3]]
QWp, _ = to_Q(X.compose(X.compose(R3p, R1), [sp.Integer(-2), sp.Integer(1)]))
print("W order", len(W) - 1, " cleared denominator", sp.factor(den), flush=True)

def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    mp.mp.dps = dps
    rec = {"dps": dps, "QW": [str(q) for q in QW], "WC": {}, "C": {}, "control": {}}
    for ts in ["0.5", "1.3", "2.5", "3.5"]:
        wc = B2.LC(QW, ts)
        c0 = B2.LC([sp.Integer(1)], ts)
        cp = B2.LC(QWp, ts)
        rec["WC"][ts], rec["C"][ts], rec["control"][ts] = mp.nstr(wc, 12), mp.nstr(c0, 12), mp.nstr(cp, 12)
        print(f"t {ts}:  C = {mp.nstr(c0, 15)}   (den W) C = {mp.nstr(wc, 6)}   control = {mp.nstr(cp, 6)}", flush=True)
    (Path(__file__).parent / "data" / "density_loop_D11.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
