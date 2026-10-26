"""The six sector blocks of the off-slice tetrahedron system, exactly over Q(x).

Lift the system A from two primes, check it on the third (as ib_tetra_lift.py), then for each sector take the
diagonal block A_ss (the sector's own D-module, its maximal cut) and the minimal operator of a vector with small
integer entries (cyclic if its order equals the block size).  Write each as an integer d/dx form in x = t^2 and in
the theta_t form of pf_12_hom.py, with its local exponents at 0, the eight cuts and infinity.
usage: python pf_16_blocks_Q.py"""
from __future__ import annotations
import json, sys, time
from math import lcm, gcd
from pathlib import Path
from flint import fmpq_poly
sys.path.insert(0, str(Path(__file__).parent))
from ib_lift_system import lift_rf, check_rf, cyclic, RQ
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import pf_lib as PL
import pf_12_hom as PH

sys.set_int_max_str_digits(10**7)
DATA = Path(__file__).parent / "data"
FILES = ["ib_tetra_pf_op.json", "ib_tetra_pf_p1_op.json", "ib_tetra_pf_p2_op.json"]
SECTORS = [(3,), (4,), (3, 5), (4, 5), (3, 4), (3, 4, 5)]


def sector(m):
    n = eval(m)
    return tuple(i for i in (3, 4, 5) if n[i] >= 0)


def to_int_dform(a):
    den = fmpq_poly([1])
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [x.n * (den // x.d) for x in a]
    L = 1
    for p in polys:
        for q in p.coeffs():
            L = lcm(L, int(q.q))
    ip = [[int((q * L).p) for q in p.coeffs()] for p in polys]
    g = 0
    for p in ip:
        for v in p:
            g = gcd(g, v)
    return [[v // g for v in p] for p in ip]


def main():
    t0 = time.time()
    ops = [json.loads((DATA / f).read_text()) for f in FILES]
    assert all(o["masters"] == ops[0]["masters"] for o in ops)
    use, hold = ops[:2], ops[2]
    mods = [o["P"] for o in use]
    nm = len(ops[0]["masters"])
    A = [[lift_rf([o["A"][i][j] for o in use], mods) for j in range(nm)] for i in range(nm)]
    assert all(x is not None for row in A for x in row), "lift failed"
    ok = all(check_rf(A[i][j], hold["A"][i][j], hold["P"]) for i in range(nm) for j in range(nm))
    print(f"system lifted from 2 primes; held-out prime agrees: {ok} ({time.time()-t0:.0f}s)", flush=True)
    assert ok
    sec = [sector(m) for m in ops[0]["masters"]]
    xs = __import__("sympy").Symbol("x")
    out = {}
    for s in SECTORS:
        idx = [i for i in range(nm) if sec[i] == s]
        AA = [[A[i][j] for j in idx] for i in idx]
        best = None
        for trial in range(4):
            v = [RQ(fmpq_poly([1 + (k * (trial + 2)) % 7])) for k in range(len(idx))]
            a = cyclic(v, AA)
            if best is None or len(a) > len(best):
                best = a
            if len(best) - 1 == len(idx):
                break
        ip = to_int_dform(best)
        b = [sum(c * xs ** k for k, c in enumerate(p)) for p in ip]
        S = PH.theta_from_dform(b)
        rep = C.local_report(S, PL.XPTS, verbose=False)
        loc = {"t0": rep["t0"], "inf": rep["inf_mu"],
               **{k: v["exps"] for k, v in rep["points"].items() if v["exps"] != [str(i) for i in range(len(S) - 1)]},
               "logs": {k: {cl["class"]: cl["dim"] - cl["log_free_solutions"] for cl in v["mono"]
                            if cl["dim"] > cl["log_free_solutions"]} for k, v in rep["points"].items()},
               "logs_t0": {cl["class"]: cl["dim"] - cl["log_free_solutions"] for cl in rep["mono_t0"]},
               "lead": rep["lead"]}
        print(f"sector {s}: block {len(idx)}, cyclic order {len(best)-1}, x-degree {G.op_degree(S)}; local {loc} "
              f"({time.time()-t0:.0f}s)", flush=True)
        out[str(s)] = {"masters": [ops[0]["masters"][i] for i in idx], "order": len(best) - 1,
                       "dform": [[str(v) for v in p] for p in ip], "theta_form": G.op_to_strings(S),
                       "local": loc, "mono": {k: v["mono"] for k, v in rep["points"].items()},
                       "mono_t0": rep["mono_t0"], "mono_inf": rep["mono_inf"]}
        (DATA / "pf_16_blocks_Q.json").write_text(json.dumps(out, indent=1, default=str))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
