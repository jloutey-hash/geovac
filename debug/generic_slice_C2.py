"""GS step C2: step C's sub/quotient test of R_V in M14 with a wider box (the b = 1 slice's exact K3 map
had poles of order 5 at t = 0 and t = +-2, beyond step C's box).  Box: t^-8 (t^2-1)^-6 (t^2-4)^-6
(t^2-9)^-6, numerator degree 60.  Controls rerun at the same box."""
import json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_C as GC  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402

Bp, a, Dg, NT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
GC.NT = NT
F4.NT = NT
HERE = Path(__file__).parent
Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
QM = [sp.sympify(q, locals={"u": GC.u}) for q in Cd["M14_Q"]]
fV = D.frob_t(GC.QV, 2, NT=NT)
LCd = json.loads((HERE / "data" / "linking_ceiling_B_ode.json").read_text())
Q4 = [sp.sympify(q, locals={"u": GC.u}) for q in LCd["4"]["Q"]]
rec = {"box": [Bp, a, Dg, NT]}
for name, Qo in [("control R_V", GC.QV), ("control F_4", Q4), ("M14", QM)]:
    for side, QQ in [("R", Qo), ("L", F4.adjoint(Qo))]:
        t1 = time.time()
        n, nc, nr = GC.test(QQ, fV, 4, Bp, a, Dg)
        rec[f"{name}|{side}"] = n
        print(f"{name:12s} ({side}): nullity {n}  [{nc} cols, {nr} rows, {time.time()-t1:.0f}s]", flush=True)
        (HERE / "data" / f"generic_slice_C2_{Bp}_{a}_{Dg}.json").write_text(json.dumps(rec, indent=1))
