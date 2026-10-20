"""P3e-3 step H: make step G's mod-p finding exact (density_loop_D9/D10 machinery, imported unmodified).
(1) exact rational gauge A = a0 + a1 th + a2 th^2, denominators t^-10 .. (t^2-4)^5, numerator degree
    <= 24, with M5 (A f) = 0 for f = 3F2(1/2,1/2,1/2;1,1;t^2/4);
(2) exact operator identity M5 o A = 0 modulo L3 = 4 th^3 - t^2 (th+1)^3 (no series);
(3) R3 = annihilator of A(Sol L3); N2 = M5 / R3 exactly (zero remainder).
"""
from __future__ import annotations
import json, sys, time
sys.set_int_max_str_digits(0)
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D9_exact_k3 as X  # noqa: E402
import density_loop_D10_k3_exact_check as D10  # noqa: E402

u, t = sp.symbols("u t")
HERE = Path(__file__).parent


def main():
    M5 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_F.json").read_text())["M5"]]
    X.M5 = [sp.Poly(q, u) for q in M5]
    X.NT = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    t0 = time.time()
    ns, labels, den = X.search(0, 5, 10, 24)
    print(f"exact nullity {ns.shape[0]}  ({time.time()-t0:.0f}s)", flush=True)
    assert ns.shape[0] >= 1
    A = X.gauge_from(list(ns.row(0)), labels, den)
    for j, aj in enumerate(A):
        nu, de = sp.fraction(sp.together(aj))
        print(f"A[{j}]: num deg {sp.degree(nu, t)}, den {sp.factor(de)}, max coeff digits {max(len(str(abs(c))) for c in sp.Poly(nu, t).coeffs())}", flush=True)
    rem = X.apply_op_mod_L3(X.M5_as_theta_coeffs(), A)
    print("M5 o A mod L3 =", rem, flush=True)
    M5L = D10.from_Q(M5)
    vecs = [[sp.cancel(a) for a in A]]
    for _ in range(3):
        vecs.append(D10.th_vec(vecs[-1]))
    Mx = sp.Matrix([[vecs[k][i] for k in range(4)] for i in range(3)])
    nsp = Mx.nullspace()
    assert len(nsp) == 1
    R3 = [sp.factor(sp.cancel(x / nsp[0][3])) for x in nsp[0]]
    print("R3 (monic):", R3, flush=True)
    N2, rem3 = D10.rdiv(M5L, R3)
    rem3 = [sp.simplify(x) for x in rem3]
    print("N2:", N2, " remainder:", rem3, flush=True)
    rec = {"nullity": ns.shape[0], "A": [str(a) for a in A], "M5A_mod_L3": [str(x) for x in rem],
           "R3": [str(x) for x in R3], "N2": [str(x) for x in N2], "rem": [str(x) for x in rem3]}
    (HERE / "data" / "p3e2w_slice_H.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
