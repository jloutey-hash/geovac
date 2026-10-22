"""FP step G: exact right division L = M62 o W12 over Q(t^2) (L: step A, W12: step D)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
u = O.u
OUT = HERE / "data" / "four_density_path_G.json"


def main():
    t0 = time.time()
    QL = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in json.loads((HERE / "data" / "four_density_path_A.json").read_text())["L"]["Q"]]
    QW = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in json.loads((HERE / "data" / "four_density_path_D.json").read_text())["W12"]["Q"]]
    L, W = O.from_Q(QL), O.from_Q(QW)
    print(f"loaded ({time.time()-t0:.0f}s)", flush=True)
    M, rem = O.rdiv(L, W)
    ok = all(c == 0 for c in rem)
    print("L = M o W12 exact:", ok, f"({time.time()-t0:.0f}s)", flush=True)
    rec = {"L_div_W12_zero": ok, "quotient_order": len(M) - 1}
    if not ok:
        rec["rem_nonzero"] = [i for i, c in enumerate(rem) if c != 0]
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
