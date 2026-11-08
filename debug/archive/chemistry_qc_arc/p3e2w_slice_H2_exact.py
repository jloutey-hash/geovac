"""P3e-3 step H2 (faster replacement of step H's back end; same mathematics):
take the exact gauge vector of step H (integerized, debug/data/p3e2w_slice_H_rawgauge.json),
A = sum_j a_j theta^j, a_j = sum_d w_{j,d} t^d / (t^2-4)^5;
(2) check M5 o A = 0 modulo L3 = 4 th^3 - t^2 (th+1)^3 exactly (sp.cancel arithmetic);
(3) R3 = monic annihilator of A(Sol L3) (3x4 rational nullspace); N2 = M5 / R3, zero remainder.
"""
from __future__ import annotations
import json, sys, time
sys.set_int_max_str_digits(0)
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D10_k3_exact_check as D10  # noqa: E402

u, t = sp.symbols("u t")
HERE = Path(__file__).parent


def main():
    raw = json.loads((HERE / "data" / "p3e2w_slice_H_rawgauge.json").read_text())
    den = sp.sympify(raw["den"])
    a = [sp.Integer(0)] * 3
    for (j, d), w in zip(raw["labels"], raw["w"]):
        a[j] += sp.Integer(w) * t**d
    A = [sp.cancel(x / den) for x in a]
    for j, x in enumerate(A):
        print(f"a{j} =", sp.factor(x), flush=True)
    M5 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_F.json").read_text())["M5"]]
    M5L = D10.from_Q(M5)
    t0 = time.time()
    rem = D10.apply_mod_L3(M5L, A)
    rem = [sp.cancel(x) for x in rem]
    print("M5 o A mod L3 =", rem, f"({time.time()-t0:.0f}s)", flush=True)
    vecs = [A]
    for _ in range(3):
        vecs.append(D10.th_vec(vecs[-1]))
    Mx = sp.Matrix([[vecs[k][i] for k in range(4)] for i in range(3)])
    nsp = Mx.nullspace(simplify=sp.cancel)
    assert len(nsp) == 1
    R3 = [sp.factor(sp.cancel(x / nsp[0][3])) for x in nsp[0]]
    print("R3 (monic):", R3, f"({time.time()-t0:.0f}s)", flush=True)
    N2, rem3 = D10.rdiv(M5L, R3)
    rem3 = [sp.cancel(x) for x in rem3]
    print("N2:", N2, " remainder:", rem3, f"({time.time()-t0:.0f}s)", flush=True)
    lead = sp.factor(sp.denom(sp.together(R3[0])))
    rec = {"A": [str(x) for x in A], "M5A_mod_L3_zero": all(x == 0 for x in rem),
           "R3": [str(x) for x in R3], "N2": [str(x) for x in N2], "rem_zero": all(x == 0 for x in rem3),
           "R3_den": str(lead)}
    (HERE / "data" / "p3e2w_slice_H.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
