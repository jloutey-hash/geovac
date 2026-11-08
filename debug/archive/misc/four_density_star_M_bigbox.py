"""FD1-FD3 step M: maximal-box mod-p sub tests (step J machinery) of reference families in the
quotients M15 = L/R4 and M14 = M15/(theta - r) (step L).  One maximal denominator
(t^2-1)^b (t^2-4)^c (t^2-9)^e and Laurent numerators t^[-a..Dg], per parity: a maximal box contains
every smaller one.  Controls: F_4's operator composed on the left/right of known pieces must show
F_4 at the same box (see __main__ control block)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_J_gauge as J  # noqa: E402
import four_density_star_ops as O  # noqa: E402

u = O.u
OUT = HERE / "data" / "four_density_star_M.json"


def bigbox(Q, f, r, b, c, e, a, Dg, NT, extra=1):
    out = []
    for par in (0, 1):
        out.append(J.test(Q, f, r, Dg, b, c, e, a=a, NT=NT, par=par, extra=extra))
    return out


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "M14"
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    Ld = json.loads((HERE / "data" / "four_density_star_L.json").read_text())
    if which == "M14":
        QQ = [sp.sympify(q) for q in Ld["M14_Q"]]
    elif which in ("M14sun", "M14x"):
        QQ = [sp.sympify(q) for q in Ld["M14_Q"]]
    elif which == "control":
        # planted: Q4 conjugated by (t^2 - 9)^2 (solutions f/(t^2-9)^2, inside the box), then composed
        # on the right with the order-1 piece theta - r of step L and divided back (exact), and
        # left-multiplied by (5 t^2 + 4)^3 (spurious coefficient factor, as in M14's cleared form)
        g = (O.W - 9) ** 2
        Q4op = O.from_Q(J.Q4)
        tw = O.compose(O.compose([1 / g], Q4op), [g])
        tw = [c * (O.K(5) * O.W + 4) ** 3 for c in tw]
        QQ, _ = O.to_Q(tw)
    else:
        raise SystemExit(which)
    NT = 900
    if which == "M14x":
        NT = 1100
    R = J.refs(NT)
    box = dict(b=4, c=3, e=4, a=10, Dg=50)
    if which == "M14x":
        tt = sp.Symbol("t")
        box = dict(b=4, c=3, e=4, a=10, Dg=74, extra=sp.expand((5 * tt**2 + 4) ** 3 * (5 * tt**4 + 2 * tt**2 + 2) ** 3))
        NT = 1100
    t0 = time.time()
    names = {"M14": ["F4", "F3", "Leg(t^2)", "Leg(t^2/4)", "Leg(t^2/9)"], "M14sun": ["Sun"],
             "control": ["F4", "F3"], "M14x": ["Sun", "F4", "F3"]}[which]
    ops = [(which, QQ)]
    for opname, Qop in ops:
        print(f"{opname}: order {max(sp.Poly(q, u).degree() for q in Qop)}, {len(Qop)} t-terms", flush=True)
        for nm in names:
            f, r_, _ = R[nm]
            res = bigbox(Qop, f, r_, NT=NT, **box)
            rec[f"{opname}|{nm}|R|box{box}"] = res
            print(f"  {opname} sub {nm}: (nullity, unknowns, rows) even/odd = {res}  ({time.time()-t0:.0f}s)", flush=True)
            OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
