"""Identify the pieces of L15 (pf_21_pieces.py) with the predicted ones, in a box large enough for them.

Box: x^-a with a 16, cuts^-b with b 12 (the quotients have exponents down to -21/2 at the cuts), the source's extra
points to power b, and for the adjoint (left-end) searches also the target's extra points to power 2 (pf_20_left_fix:
adjoint solutions have poles at the target's apparent points); numerator degree M.
Tests: R_V into the left order-4 quotient L15/R_S; the four Legendre pullbacks of PF1 into, and (adjoint) out of,
the order-3 quotients R_A/R_top ({g2,C}: logs at 4, 16) and R_B/R_top ({g1,C}: logs at 1/4, 49/4) and the middle
order-4 piece X = R_S / R_7.  Controls at this box: planted F o S(theta + 1) and S o F for z1 with F = R_A/R_top.
usage: python pf_22_pieces_hom.py <tag> [M]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import pf_lib  # noqa: F401
import pf_12_hom as P
import pf_20_left_fix as LF

DATA = Path(__file__).parent / "data"
AMAX, BMAX = 16, 12


def nul(r):
    return {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")}


def right(S, F, M):
    P.den_coeffs = LF._orig_den
    return P.hom(S, F, AMAX, BMAX, M)


def left(S, F, M):
    P.den_coeffs = LF.den_with_target
    LF.TEXTRA.update(poly=P.source_extra(F), k=2)
    r = P.hom(G.adjoint(S), G.adjoint(F), AMAX, BMAX, M)
    P.den_coeffs = LF._orig_den
    return r


def main():
    tag = sys.argv[1]
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    t0 = time.time()
    pc = json.loads((DATA / "pf_21_pieces.json").read_text())
    QA = G.parse_op(pc["R_A / R_top"]["op"]); QB = G.parse_op(pc["R_B / R_top"]["op"])
    LS = G.parse_op(pc["L15 / R_S"]["op"])
    RS = G.parse_op(pc["R_S"]["op"])
    R7 = G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["chain_R7"])
    ok, X = C.divides(RS, R7)
    assert ok
    RV = __import__("pf_10_rvmap").RV
    out = {"box": {"amax": AMAX, "bmax": BMAX, "M": M}, "X_order": len(X) - 1}

    def rec(name, r):
        out[name] = r
        print(f"{name}: {nul(r)} ({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_22_pieces_hom_{tag}.json").write_text(json.dumps(out, indent=1, default=str))

    out["X"] = G.op_to_strings(X)
    z1 = P.legendre_source("z1")
    rec("planted right z1 into QA o z1(theta+1)", right(z1, G.compose(QA, G.shift_theta(z1, 1)), M))
    rec("planted left z1 o QA", left(z1, G.compose(z1, QA), M))
    rec("Hom(R_V, L15/R_S)", right(RV, LS, M))
    for z in ("z1", "z2", "d1", "d2"):
        S = P.legendre_source(z)
        for nm, T in (("R_A/R_top", QA), ("R_B/R_top", QB), ("X", X)):
            rec(f"right Hom({z}, {nm})", right(S, T, M))
            rec(f"left Hom({z}*, {nm}*)", left(S, T, M))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
