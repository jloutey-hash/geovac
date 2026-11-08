"""FP step J: log-free Frobenius solutions of W11 = W12/(theta-2) at t = 0 (free slots 0) and the orders of
their minimal annihilators (mod p guessing, t-degree <= 2*Imax).  A solution with minimal annihilator of
order k < 11 gives a right factor of order k."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
OUT = HERE / "data" / "four_density_path_J.json"


def main():
    Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_path_I.json").read_text())["quotient_Q"]]
    qc = T.qcoeffs(Q)
    roots = sorted(set(r for r in sp.roots(sp.Poly(sp.sympify(Q[0]).subs(sp.Symbol('u'), T.u), T.u).as_expr(), T.u)))
    print("exponents at 0:", sp.factor(sp.sympify(Q[0])), flush=True)
    rec = {}
    NT = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    Imax = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    for lam in roots:
        t0 = time.time()
        ser, free = T.frob(qc, int(lam), NT)
        if ser is None:
            print(f"lam={lam}: no log-free solution (inconsistent resonance at {free})", flush=True)
            rec[str(lam)] = "log"
            continue
        top = int(lam) + 2 * NT
        k, I = T.min_annihilator_order([(ser, top)], 11, Imax)
        print(f"lam={lam}: free slots {free[:6]}; minimal annihilator order {k}, t-degree {None if I is None else 2*I}"
              f" ({time.time()-t0:.0f}s)", flush=True)
        rec[str(lam)] = {"free": free, "order": k, "tdeg": None if I is None else 2 * I}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
