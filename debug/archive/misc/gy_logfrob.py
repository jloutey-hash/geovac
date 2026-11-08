"""GY: the Jordan structure of the local monodromy at x = 0 in an exponent class e0 mod Z, from the dimensions
d_J = number of independent local solutions of the form sum_{j<=J} log(x)^j x^(e0) P_j(x) (mod p), J = 0, 1, 2, ...:
d_J = dim ker (M - lambda)^(J+1).  The block sizes follow (number of blocks of size >= J+1 is d_J - d_{J-1}).
A composition factor equivalent to F_3's operator twisted by t needs a block of size 3 (its MUM point at x = 0).
theta = theta_t = 2 x d/dx acts by theta(x^s log^j x) = 2 s x^s log^j x + 2 j x^s log^(j-1) x.
Controls: F_3's operator: d = (1, 2, 3) (one block of size 3); the Slater G4 at gamma = 1, class 1/2: a block of size 3
(it contains the MUM K3 piece R3) plus one of size 1: d = (2, 3, 4).
usage: python gy_logfrob.py controls | python gy_logfrob.py <op json> <key|Q|theta_form> <e0> <N> <tag>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import gy_logfree as LF

DATA = Path(__file__).parent / "data"


def dims(L, e0, N, Jmax=3, p=H.P0):
    Lp = G.op_mod(L, p)
    D = G.op_degree(L)
    e0p = G.frac_mod(e0, p)
    out = []
    for J in range(Jmax + 1):
        # basis index (n, j) -> column; target index (m, j) rows for m <= N
        cols = []
        for n in range(N + 1):
            for j in range(J + 1):
                # apply L to x^(e0+n) log^j: first theta^k on the log-polynomial, then multiply by c_k(x)
                vec = {}
                for k, c in enumerate(Lp):
                    ck = [int(v) for v in c.coeffs()]
                    if not any(ck):
                        continue
                    # theta^k (x^s log^j) = sum_i a_i x^s log^(j-i): iterate
                    cur = {j: 1}
                    s = (e0p + n) % p
                    for _ in range(k):
                        nxt = {}
                        for jj, a in cur.items():
                            nxt[jj] = (nxt.get(jj, 0) + 2 * s % p * a) % p
                            if jj >= 1:
                                nxt[jj - 1] = (nxt.get(jj - 1, 0) + 2 * jj % p * a) % p
                        cur = nxt
                    for i, cv in enumerate(ck):
                        if cv:
                            for jj, a in cur.items():
                                key = (n + i, jj)
                                vec[key] = (vec.get(key, 0) + cv * a) % p
                cols.append(vec)
        rows = []
        for m in range(N + 1):
            for jj in range(J + 1):
                rows.append([col.get((m, jj), 0) for col in cols])
        out.append(len(G.nullspace_mod(rows, len(cols), p)))
    return out


def blocks(d):
    """Jordan block sizes from d_J = dim ker (M-lambda)^(J+1)"""
    dd = [0] + d
    ge = [dd[J + 1] - dd[J] for J in range(len(d))]      # number of blocks of size >= J+1
    sizes = []
    for J in range(len(d)):
        cnt = ge[J] - (ge[J + 1] if J + 1 < len(ge) else 0)
        sizes += [J + 1] * cnt
    return sorted(sizes, reverse=True)


if __name__ == "__main__":
    if sys.argv[1] == "controls":
        d3 = dims(G.L3_op(), Fraction(0), 40)
        d4 = dims(G.load_G4(), Fraction(1, 2), 60)
        print("F_3 at x=0, class 0: d_J =", d3, "blocks", blocks(d3), "(expect [3])")
        print("Slater G4 (gamma=1) at x=0, class 1/2: d_J =", d4, "blocks", blocks(d4), "(expect [3, 1])")
        (DATA / "gy_logfrob_controls.json").write_text(json.dumps({"F3": {"d": d3, "blocks": blocks(d3)}, "G4": {"d": d4, "blocks": blocks(d4)}}))
    else:
        L = LF.load(sys.argv[1], sys.argv[2])
        e0, N, tag = Fraction(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        d = dims(L, e0, N)
        print(f"{tag}: order {len(L)-1}, class x^{e0} at x = 0: d_J = {d}, Jordan blocks {blocks(d)}")
        (DATA / f"gy_logfrob_{tag}.json").write_text(json.dumps({"order": len(L) - 1, "class": str(e0), "d_J": d, "blocks": blocks(d)}))
