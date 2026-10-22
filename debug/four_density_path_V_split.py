"""FP step V: exact splitting W12 = Y4 o R8 (R8 from step U), R8 = X4 o R3 o (theta - 2), and the equivalence
Y4 ~ A4 (= Mellin operator of P_h4) through C (step U: R8 P = C P_h4): Y4 o C = 0 mod A4 and C injective on Sol(A4)
(the annihilator of C(Sol A4) equals Y4).  Then A4 ~ F_4 through G = [t^-2 (theta-3)^2]^3 o t^3 (P_{z h} = t^-2 (theta-3)^2 P_h):
A4 o G = 0 mod Q4, G injective."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_V.json"


def ld(path, key):
    return [sp.sympify(q).subs(sp.Symbol("u"), u) for q in json.loads((HERE / "data" / path).read_text())[key]]


def image_annihilator(Aop_delta, X):
    """annihilator (monic, delta form) of A(Sol X), A in delta form (even); None if A not injective."""
    r = len(X) - 1
    S = [-c / X[r] for c in X[:r]]

    def dvec(vv):
        out = [O.d(vv[0]) + vv[-1] * S[0]]
        for k in range(1, r):
            out.append(vv[k - 1] + O.d(vv[k]) + vv[-1] * S[k])
        return out
    # reduce A mod X first
    _, Ared = O.rdiv(Aop_delta, X) if len(Aop_delta) > r else (None, list(Aop_delta))
    v0 = list(Ared) + [O.K(0)] * (r - len(Ared))
    vecs = [v0]
    for _ in range(r):
        vecs.append(dvec(vecs[-1]))
    Mx = [[vecs[k][i] for k in range(r + 1)] for i in range(r)]
    piv, rk = [], 0
    for col in range(r + 1):
        pr = next((i for i in range(rk, r) if Mx[i][col] != 0), None)
        if pr is None:
            continue
        Mx[rk], Mx[pr] = Mx[pr], Mx[rk]
        iv = 1 / Mx[rk][col]
        Mx[rk] = [x * iv for x in Mx[rk]]
        for i in range(r):
            if i != rk and Mx[i][col] != 0:
                fct = Mx[i][col]
                Mx[i] = [x - fct * y for x, y in zip(Mx[i], Mx[rk])]
        piv.append(col)
        rk += 1
    free = [c for c in range(r + 1) if c not in piv]
    if free != [r]:
        return None
    row = [O.K(0)] * (r + 1)
    row[r] = O.K(1)
    for i, pc in enumerate(piv):
        row[pc] = -Mx[i][r]
    return row


def main():
    rec = {}
    t0 = time.time()
    W12 = O.from_Q(ld("four_density_path_D.json", "W12") if False else [sp.sympify(q).subs(sp.Symbol("u"), u) for q in json.loads((HERE / "data" / "four_density_path_D.json").read_text())["W12"]["Q"]])
    U = json.loads((HERE / "data" / "four_density_path_U.json").read_text())
    RQ = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in U["R_Q"]]
    while RQ and RQ[-1] == 0:
        RQ.pop()
    CQ = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in U["C_Q"]]
    while CQ and CQ[-1] == 0:
        CQ.pop()
    R8, Cop = O.from_Q(RQ), O.from_Q(CQ)
    print("R8 order", len(R8) - 1, " C order", len(Cop) - 1, flush=True)
    Y4, rem = O.rdiv(W12, R8)
    rec["W12_eq_Y4_o_R8"] = all(c == 0 for c in rem)
    print("W12 = Y4 o R8:", rec["W12_eq_Y4_o_R8"], f"({time.time()-t0:.0f}s)", flush=True)
    S = json.loads((HERE / "data" / "four_density_path_S.json").read_text())
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    R3 = O.from_Q([sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]])
    R4low = O.compose(R3, O.from_Q([u - 2]))
    X4, rem = O.rdiv(R8, R4low)
    rec["R8_eq_X4_o_R3_o_th2"] = all(c == 0 for c in rem)
    print("R8 = X4 o R3 o (theta-2):", rec["R8_eq_X4_o_R3_o_th2"], f"({time.time()-t0:.0f}s)", flush=True)
    X4Q, Dx = O.to_Q(X4)
    rec["X4_Q"] = [str(q) for q in X4Q]
    Y4Q, Dy = O.to_Q(Y4)
    rec["Y4_Q"] = [str(q) for q in Y4Q]
    rec["R8_Q"] = [str(q) for q in RQ]
    rec["C_Q"] = [str(q) for q in CQ]
    OUT.write_text(json.dumps(rec, indent=1))
    O2 = json.loads((HERE / "data" / "four_density_path_O2.json").read_text())
    A4Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in O2["A4_Q"]]
    A4 = O.from_Q(A4Q)
    comp = O.compose(Y4, Cop)
    _, rem = O.rdiv(comp, A4)
    rec["Y4_o_C_mod_A4_zero"] = all(c == 0 for c in rem)
    img = image_annihilator(Cop, A4)
    Y4m = [c / Y4[-1] for c in Y4]
    rec["C_injective_image_ann_eq_Y4"] = img is not None and all(a == b for a, b in zip(img, Y4m))
    print("Y4 o C = 0 mod A4:", rec["Y4_o_C_mod_A4_zero"], " C injective, image annihilator = Y4:", rec["C_injective_image_ann_eq_Y4"], f"({time.time()-t0:.0f}s)", flush=True)
    # A4 ~ Q4 via G = [t^-2 (theta-3)^2]^3 o t^3 = t^-6 (theta-7)^2 (theta-5)^2 (theta-3)^2 o t^3 = t^-3 (theta-4)^2 (theta-2)^2 theta^2
    Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]
    # t^-3 h(theta): conjugation: A4 o t^-3 = t^-3 o A4(theta - 3); work with A4sh = A4(u - 3) and Gs = (theta-4)^2 (theta-2)^2 theta^2
    A4sh = O.from_Q([sp.expand(q.subs(u, u - 3)) for q in A4Q])
    Gs = O.from_Q([sp.expand((u - 4) ** 2 * (u - 2) ** 2 * u**2)])
    comp = O.compose(A4sh, Gs)
    _, rem = O.rdiv(comp, O.from_Q(Q4))
    rec["A4_o_G_mod_Q4_zero"] = all(c == 0 for c in rem)
    img = image_annihilator(Gs, O.from_Q(Q4))
    A4m = [c / A4sh[-1] for c in A4sh]
    rec["G_injective_image_ann_eq_A4"] = img is not None and all(a == b for a, b in zip(img, A4m))
    print("A4 o G = 0 mod Q4:", rec["A4_o_G_mod_Q4_zero"], " G injective, image annihilator = A4:", rec["G_injective_image_ann_eq_A4"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
