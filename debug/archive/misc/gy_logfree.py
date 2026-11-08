"""GY: the number of independent LOG-FREE local solutions of an operator at x = 0 in a given exponent class e0 mod Z
(mod p): the nullity of L acting on x^(e0) * (polynomials of degree <= N), truncated consistently.  This is the number
of Jordan blocks of the local monodromy at x = 0 with eigenvalue exp(2 pi i e0); with n solutions in the class, a Jordan
block of size >= 3 (maximal unipotency, as at F_3's x = 0) needs nullity <= n - 2.  Positive control: F_3's operator
(class 0, 3 solutions) must give nullity 1; the Slater G4 at gamma = 1 (class 1/2 in x: t-exponents 1,3,3,5) must give
nullity <= 2 (it contains the K3 piece R3, MUM at t = 0).
usage: python gy_logfree.py <op json> <key|Q|theta_form> <e0 (Fraction)> <N> [tag]"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

DATA = Path(__file__).parent / "data"


def logfree_nullity(L, e0, N, p=H.P0):
    """columns: L applied to x^(e0 + n), n = 0..N; L = sum_k c_k(x) theta^k with theta x^s = 2 s x^s (theta_t)."""
    Lp = G.op_mod(L, p)
    D = G.op_degree(L)
    e0p = G.frac_mod(e0, p)
    cols = []
    for n in range(N + 1):
        s = (e0p + n) % p
        col = [0] * (N + D + 1)
        for k, c in enumerate(Lp):
            ck = [int(v) for v in c.coeffs()]
            w = pow(2 * s % p, k, p)
            for i, cv in enumerate(ck):
                if cv:
                    col[n + i] = (col[n + i] + cv * w) % p
        cols.append(col)
    # rows m = 0..N (the coefficient of x^(e0+m) only involves n <= m, complete for m <= N)
    rows = [[col[m] for col in cols] for m in range(N + 1)]
    return len(G.nullspace_mod(rows, N + 1, p))


def load(path, key):
    d = json.loads(Path(path).read_text())
    if key == "Q":
        return G.from_Q([sp.sympify(q, locals={"u": G.u_sym}) for q in d["Q"]])
    if key == "theta_form":
        return G.parse_op(d["theta_form"])
    return G.parse_op(d[key])


if __name__ == "__main__":
    if sys.argv[1] == "controls":
        f3 = G.L3_op()
        n3 = logfree_nullity(f3, Fraction(0), 60)
        g4 = G.load_G4()
        n4 = logfree_nullity(g4, Fraction(1, 2), 80)
        print("F_3 at x = 0, class 0 (3 solutions): log-free nullity", n3, "(expect 1: MUM)")
        print("Slater G4 (gamma = 1) at x = 0, class 1/2 (4 solutions): log-free nullity", n4, "(expect <= 2: contains the MUM K3 piece)")
        (DATA / "gy_logfree_controls.json").write_text(json.dumps({"F3_class0": n3, "G4_class_half": n4}))
    else:
        L = load(sys.argv[1], sys.argv[2])
        e0 = Fraction(sys.argv[3]); N = int(sys.argv[4])
        tag = sys.argv[5] if len(sys.argv) > 5 else sys.argv[2]
        nul = logfree_nullity(L, e0, N)
        ind = G.exps_of(G.indicial_t0(L))[0]
        ncls = sum(1 for e in ind if (Fraction(str(e)) / 2 - e0).denominator == 1)
        print(f"{tag}: order {len(L)-1}; t-exponents at 0: {[str(e) for e in ind]}; class x^{e0}: {ncls} solutions; log-free nullity {nul}; "
              f"largest possible Jordan block <= {ncls - nul + 1}")
        (DATA / f"gy_logfree_{tag}.json").write_text(json.dumps({"order": len(L) - 1, "t_exponents": [str(e) for e in ind], "class": str(e0),
                                                                 "solutions_in_class": ncls, "logfree_nullity": nul, "max_block_bound": ncls - nul + 1}))
