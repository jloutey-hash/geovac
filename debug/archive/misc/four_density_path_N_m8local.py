"""FP step N: log-free Frobenius solutions of M8 = W11/R3 at t = 0 and the orders of their minimal annihilators
(mod p).  Usage: python ... NT Imax"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
OUT = HERE / "data" / "four_density_path_N.json"


def run(Q, name, NT, Imax, kmax):
    qc = T.qcoeffs(Q)
    q0 = sp.sympify(Q[0]).subs(sp.Symbol('u'), T.u)
    roots = sorted(set(r for r in sp.roots(sp.Poly(q0, T.u).as_expr(), T.u)))
    print(name, "t^0 coefficient:", sp.factor(q0), flush=True)
    rec = {}
    for lam in roots:
        t0 = time.time()
        ser, free = T.frob(qc, int(lam), NT)
        if ser is None:
            print(f"  lam={lam}: log forced (resonance at {free})", flush=True)
            rec[str(lam)] = "log"
            continue
        k, I = T.min_annihilator_order([(ser, int(lam) + 2 * NT)], kmax, Imax)
        print(f"  lam={lam}: free {free[:6]}; minimal annihilator order {k}, t-degree {None if not isinstance(I, int) else 2*I} ({time.time()-t0:.0f}s)", flush=True)
        rec[str(lam)] = {"free": free, "order": k, "tdeg": None if not isinstance(I, int) else 2 * I}
    return rec


def main():
    NT, Imax = int(sys.argv[1]), int(sys.argv[2])
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    rec = run([sp.sympify(q) for q in K["M8_Q"]], "M8", NT, Imax, 7)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
