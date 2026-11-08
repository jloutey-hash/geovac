"""GS4 (record docs/predictions/2026-10-21_generic_slice.md): no split.  Does Sol(R_V) map into Sol(W7)
itself, W7 = R o E2 o (theta+2)?  A hit would make the K3 piece a sub of W7, i.e. G_b a rational
combination of V-type functions plus elementary ones.  Box from the leading-term bound of step F's local
exponents (pole order at t0 <= max exp R_V - min exp W7; numerator degree <= max rho_V - min rho_W at oo),
enlarged by a margin; mod-p series test of generic_slice_C at that box.
Positive control at the same box: X = E3 o R_V (split by construction).  Negative control: W7 with R_V
replaced by F_4's operator as reference."""
from __future__ import annotations
import json, math, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_C as GC  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
from generic_slice_B import compose, coeff_list, t, u  # noqa: E402

HERE = Path(__file__).parent


def exps(F, name, t0):
    return [sp.Rational(k) for k, m in F[f"{name}|{t0}"].items() for _ in range(m)]


def main(margin=4):
    F = json.loads((HERE / "data" / "generic_slice_F.json").read_text())
    rec = {}
    bounds = {}
    for t0 in ["0", "1", "2", "3"]:
        eV = exps(F, "R_V", t0)
        eW = exps(F, "W7", t0)
        bounds[t0] = int(math.floor(max(eV) - min(eW)))
    bounds["oo"] = int(math.floor(max(exps(F, "R_V", "oo")) - min(exps(F, "W7", "oo"))))
    print("leading-term bounds (pole order at t0; numerator excess at oo):", bounds, flush=True)
    rec["bounds"] = bounds
    Bp = max(bounds["1"], bounds["2"], bounds["3"], 0) + margin
    a = max(bounds["0"], 0) + margin
    Dd = json.loads((HERE / "data" / "generic_slice_D2.json").read_text())
    Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
    R = [sp.sympify(c, locals={"t": t}) for c in Dd["R"]]
    E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
    E3 = compose(E2, [sp.Integer(2), sp.Integer(1)])
    W7 = compose(R, E3)
    QW = GC.to_Q(W7)
    # numerator degree: denominator degree 6 Bp, plus excess at infinity, plus margin
    Dg = 6 * Bp + max(bounds["oo"], 0) + margin
    NT = 2 * (4 * (Dg + a + 1)) + 300
    GC.NT = NT
    F4.NT = NT
    rec["box"] = {"Bp": Bp, "a": a, "Dg": Dg, "NT": NT}
    print("box", rec["box"], flush=True)
    fV = D.frob_t(GC.QV, 2, NT=NT)
    X = compose(E3, coeff_list(GC.QV))
    QX = GC.to_Q(X)
    LCd = json.loads((HERE / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4 = [sp.sympify(q, locals={"u": u}) for q in LCd["4"]["Q"]]
    for name, Qo, f in [("control E3 o R_V (split by construction)", QX, fV), ("W7", QW, fV)]:
        t1 = time.time()
        n, nc, nr = GC.test(Qo, f, 4, Bp, a, Dg)
        rec[name] = n
        print(f"{name}: nullity {n}  [{nc} cols, {nr} rows, {time.time()-t1:.0f}s]", flush=True)
        (HERE / "data" / "generic_slice_G.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
