"""The pieces of Y's order-15 operator, cut out by sector quotients.

For a downward-closed set I of sectors containing both sunrise sectors, Y's image in M/M_I has a minimal operator L_I,
a right factor of L15 (pf_19) and of L19.  Quotients by the top-only right factor isolate the middle pieces:
  I_top = all but top                  -> R_top (top sector alone)
  I_A   = sunrise, {g1,C}, {g1,g2}     -> R_A = E_A' o R_top  ({g2,C} piece, logs at 4 and 16 expected)
  I_B   = sunrise, {g2,C}, {g1,g2}     -> R_B = E_B' o R_top  ({g1,C} piece, logs at 1/4 and 49/4 expected)
  I_C   = sunrise, {g1,C}, {g2,C}      -> R_C = G_C o R_top   ({g1,g2} piece)
  I_S   = sunrise                      -> R_S (order 11 = L15 without its R_V left factor)
Each quotient by exact right division; local data; exact divisibility checks against L15, L19 and chain_R7;
Hom tests of the order-two pieces against the four Legendre pullbacks of PF1 (both directions).
usage: python pf_21_pieces.py"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
from ib_lift_system import lift_rf, check_rf, cyclic
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import pf_lib as PL
import pf_12_hom as PH
from pf_16_blocks_Q import to_int_dform

sys.set_int_max_str_digits(10**7)
DATA = Path(__file__).parent / "data"
FILES = ["ib_tetra_pf_op.json", "ib_tetra_pf_p1_op.json", "ib_tetra_pf_p2_op.json"]
SUN = [(3,), (4,)]
IDEALS = {"R_top": SUN + [(3, 5), (4, 5), (3, 4)], "R_A": SUN + [(3, 5), (3, 4)], "R_B": SUN + [(4, 5), (3, 4)],
          "R_C": SUN + [(3, 5), (4, 5)], "R_S": SUN}


def sector(m):
    n = eval(m)
    return tuple(i for i in (3, 4, 5) if n[i] >= 0)


def local(S):
    rep = C.local_report(S, PL.XPTS, verbose=False)
    return {"order": len(S) - 1, "t0": rep["t0"], "inf": rep["inf_mu"],
            **{k: v["exps"] for k, v in rep["points"].items() if v["exps"] != [str(i) for i in range(len(S) - 1)]},
            "logs": {k: {cl["class"]: cl["dim"] - cl["log_free_solutions"] for cl in v["mono"]
                         if cl["dim"] > cl["log_free_solutions"]} for k, v in rep["points"].items()},
            "logs_t0": {cl["class"]: cl["dim"] - cl["log_free_solutions"] for cl in rep["mono_t0"]}}


def nul(r):
    return {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")}


def main():
    t0 = time.time()
    ops = [json.loads((DATA / f).read_text()) for f in FILES]
    ms = ops[0]["masters"]; nm = len(ms)
    A = [[lift_rf([o["A"][i][j] for o in ops[:2]], [o["P"] for o in ops[:2]]) for j in range(nm)] for i in range(nm)]
    assert all(check_rf(A[i][j], ops[2]["A"][i][j], ops[2]["P"]) for i in range(nm) for j in range(nm))
    c = [lift_rf([o["c"][i] for o in ops], [o["P"] for o in ops]) for i in range(nm)]
    sec = [sector(m) for m in ms]
    x = sp.Symbol("x")
    R, out = {}, {}
    for name, ideal in IDEALS.items():
        idx = [i for i in range(nm) if sec[i] not in ideal]
        a = cyclic([c[i] for i in idx], [[A[i][j] for j in idx] for i in idx])
        ip = to_int_dform(a)
        S = PH.theta_from_dform([sum(v * x ** k for k, v in enumerate(p)) for p in ip])
        R[name] = G.content_free(G.shift_theta(S, -3))
        out[name] = {"order": len(a) - 1, "op": G.op_to_strings(R[name])}
        print(f"{name}: order {len(a)-1} ({time.time()-t0:.0f}s)", flush=True)
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    L19, R7 = G.parse_op(pc["chain_R8"]), G.parse_op(pc["chain_R7"])
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    div = {}
    for nm_, Rr in list(R.items()) + [("L15", L15)]:
        for tn, T in (("L15", L15), ("L19", L19)):
            if nm_ != tn:
                div[f"{tn} / {nm_}"] = C.divides(T, Rr)[0]
        div[f"{nm_} / chain_R7"] = C.divides(Rr, R7)[0]
    out["divisions_exact"] = div
    print("exact right divisions:", div, flush=True)
    Q = {}
    for nm_ in ("R_A", "R_B", "R_C"):
        ok, q = C.divides(R[nm_], R["R_top"])
        assert ok, nm_
        Q[nm_ + " / R_top"] = q
    ok, q = C.divides(L15, R["R_S"]); assert ok
    Q["L15 / R_S"] = q
    for k, q in Q.items():
        out[k] = {"op": G.op_to_strings(q), "local": local(q)}
        print(f"{k}: {out[k]['local']} ({time.time()-t0:.0f}s)", flush=True)
        (DATA / "pf_21_pieces.json").write_text(json.dumps(out, indent=1, default=str))
    RV = __import__("pf_10_rvmap").RV
    out["Hom(R_V, L15/R_S)"] = PH.hom(RV, Q["L15 / R_S"], 8, 6, 200)
    out["Hom(L15/R_S, R_V)"] = PH.hom(Q["L15 / R_S"], RV, 8, 6, 200)
    print("R_V vs the left order-4 quotient of L15:", nul(out["Hom(R_V, L15/R_S)"]), nul(out["Hom(L15/R_S, R_V)"]), flush=True)
    for k in ("R_A / R_top", "R_B / R_top"):
        if len(Q[k]) - 1 != 2:
            continue
        for z in ("z1", "z2", "d1", "d2"):
            Sz = PH.legendre_source(z)
            r = PH.hom(Sz, Q[k], 8, 6, 200)
            out[f"Hom({z}, {k})"] = r
            print(f"Hom({z}, {k}): {nul(r)} ({time.time()-t0:.0f}s)", flush=True)
    (DATA / "pf_21_pieces.json").write_text(json.dumps(out, indent=1, default=str))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
