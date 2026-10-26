"""Exact order-15 operator for Y: the IBP system over Q, quotiented by the sunrise relation of pf_18_quotient.py.

Lift A from two primes (checked on the third) and c from all three (as ib_tetra_lift.py --c-all), identify the {g1}
masters with lam times the {g2} masters, and build Y's minimal operator in the quotient by the cyclic vector.
lam = 1/2 is the tadpole ratio g1/g2; lam = 1/3 is the control (also order 15, but it should not annihilate Y).
Writes theta_t forms on Y (pf_9_onY.py keys "L15_lam_1/2", "L15_lam_1/3") and the integer d/dx forms on S = Y/t^3.
usage: python pf_19_quotient_Q.py"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
from flint import fmpq_poly, fmpq
sys.path.insert(0, str(Path(__file__).parent))
from ib_lift_system import lift_rf, check_rf, cyclic, RQ
import correlated_triangle_g32_lib as G
import pf_12_hom as PH
from pf_16_blocks_Q import to_int_dform

sys.set_int_max_str_digits(10**7)
DATA = Path(__file__).parent / "data"
FILES = ["ib_tetra_pf_op.json", "ib_tetra_pf_p1_op.json", "ib_tetra_pf_p2_op.json"]


def main():
    t0 = time.time()
    ops = [json.loads((DATA / f).read_text()) for f in FILES]
    ms = ops[0]["masters"]
    nm = len(ms)
    mods2 = [o["P"] for o in ops[:2]]
    A = [[lift_rf([o["A"][i][j] for o in ops[:2]], mods2) for j in range(nm)] for i in range(nm)]
    assert all(check_rf(A[i][j], ops[2]["A"][i][j], ops[2]["P"]) for i in range(nm) for j in range(nm))
    c = [lift_rf([o["c"][i] for o in ops], [o["P"] for o in ops]) for i in range(nm)]
    assert all(x is not None for x in c)
    print(f"lifted ({time.time()-t0:.0f}s)", flush=True)

    def key(m):
        return eval(m)[:3]
    g1 = {key(m): i for i, m in enumerate(ms) if eval(m)[3] >= 0 and eval(m)[4] < 0 and eval(m)[5] < 0}
    g2 = {key(m): i for i, m in enumerate(ms) if eval(m)[3] < 0 and eval(m)[4] >= 0 and eval(m)[5] < 0}
    pair = {g1[k]: g2[k] for k in g1}
    keep = [i for i in range(nm) if i not in pair]
    pos = {i: k for k, i in enumerate(keep)}
    x = sp.Symbol("x")
    out = {}
    for lam_q in (Fraction(1, 2), Fraction(1, 3)):
        lam = RQ(fmpq_poly([fmpq(lam_q.numerator, lam_q.denominator)]))
        cn = [c[i] for i in keep]
        for i, j in pair.items():
            cn[pos[j]] = cn[pos[j]] + lam * c[i]
        An = [[A[i][j] for j in keep] for i in keep]
        for r, i in enumerate(keep):
            for a, b in pair.items():
                if not A[i][a].is_zero():
                    An[r][pos[b]] = An[r][pos[b]] + lam * A[i][a]
        a = cyclic(cn, An)
        ip = to_int_dform(a)
        S = PH.theta_from_dform([sum(v * x ** k for k, v in enumerate(p)) for p in ip])
        LY = G.shift_theta(S, -3)
        nmk = f"L15_lam_{lam_q}"
        out[nmk] = G.op_to_strings(LY)
        out[nmk + "_dform_S"] = [[str(v) for v in p] for p in ip]
        print(f"lam = {lam_q}: order {len(a)-1}, d/dx degrees {[len(p)-1 for p in ip]}, largest coefficient "
              f"{max(abs(v).bit_length() for p in ip for v in p)} bits ({time.time()-t0:.0f}s)", flush=True)
        (DATA / "pf_19_quotient_Q.json").write_text(json.dumps(out))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
