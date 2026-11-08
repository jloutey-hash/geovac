"""Does Y need each piece of L15?  Right factors of L15 that drop one piece, applied to Y (pf_9_onY.py).
R_S  = quotient by both sunrise sectors              (drops R_V; keeps z1, z2 and R_7; order 11)
R_D  = quotient by the sunrise sectors and {g1,C}     (drops R_V and the z1 piece; order 9 expected)
R_E  = quotient by the sunrise sectors and {g2,C}     (drops R_V and the z2 piece; order 9 expected)
R_C  = quotient by the sunrise sectors, {g1,C}, {g2,C} (= R_7, order 7)
Writes theta_t forms on Y to pf_23_needs.json.
usage: python pf_23_needs.py"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
from ib_lift_system import lift_rf, check_rf, cyclic
import correlated_triangle_g32_lib as G
import pf_12_hom as PH
from pf_16_blocks_Q import to_int_dform
from pf_21_pieces import sector, FILES, SUN

DATA = Path(__file__).parent / "data"
IDEALS = {"R_S": SUN, "R_D": SUN + [(3, 5)], "R_E": SUN + [(4, 5)], "R_C": SUN + [(3, 5), (4, 5)]}


def main():
    t0 = time.time()
    ops = [json.loads((DATA / f).read_text()) for f in FILES]
    ms = ops[0]["masters"]; nm = len(ms)
    A = [[lift_rf([o["A"][i][j] for o in ops[:2]], [o["P"] for o in ops[:2]]) for j in range(nm)] for i in range(nm)]
    assert all(check_rf(A[i][j], ops[2]["A"][i][j], ops[2]["P"]) for i in range(nm) for j in range(nm))
    c = [lift_rf([o["c"][i] for o in ops], [o["P"] for o in ops]) for i in range(nm)]
    sec = [sector(m) for m in ms]
    x = sp.Symbol("x")
    out = {}
    for name, ideal in IDEALS.items():
        idx = [i for i in range(nm) if sec[i] not in ideal]
        a = cyclic([c[i] for i in idx], [[A[i][j] for j in idx] for i in idx])
        ip = to_int_dform(a)
        S = PH.theta_from_dform([sum(v * x ** k for k, v in enumerate(p)) for p in ip])
        out[name] = G.op_to_strings(G.content_free(G.shift_theta(S, -3)))
        print(f"{name}: order {len(a)-1} ({time.time()-t0:.0f}s)", flush=True)
    (DATA / "pf_23_needs.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
