"""Jordan block sizes of local monodromy, exactly (mod two primes), for the order-4 piece of Y at t = 0 and infinity.

Method.  L = sum_k c_k(x) theta_t^k, x = t^2 (the correlated_triangle_g32_lib convention); write L = sum_j x^j Q_j(theta_t).
For each class of t-exponents s mod 2 (= x-exponent class mod 1, i.e. one eigenvalue of the x-monodromy T), take the
smallest root s0 of the class and solve L y = 0 for y = sum_{n<=N} sum_{k<=K} c[n][k] t^(s0+2n) (log t)^k over F_p.
The kernel dimension d_K is dim ker (T - lambda)^(K+1) on that class (log x = 2 log t; monodromy shifts log by 2 pi i).
Blocks of size > K number d_K - d_(K-1).  Infinity: correlated_triangle_g32_lib3.at_infinity (w = 1/x).

Controls (known answers): theta_x^2 (theta_x - 1)^2 -> 2+2 at 0; theta_x^3 (theta_x - 1) -> 3+1 at 0;
AESZ #1 (quintic) -> one block of 4 at 0, four blocks of size 1 at infinity.
Targets: R_V (pf_10_rvmap.RV) and the order-4 left quotient L15 / R_S (pf_21_pieces.json): equivalent by an exact map,
so they must agree (two routes).
usage: python pf_27_jordan.py"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc

DATA = Path(__file__).parent / "data"
u = G.u_sym
PRIMES = [1000000007, 998244353]


def qforms(L):
    """Q_j(u) as lists of Fractions (coefficient of u^k), j = 0..deg"""
    D = G.op_degree(L)
    out = []
    for j in range(D + 1):
        row = []
        for k in range(len(L)):
            cs = L[k].coeffs()
            row.append(Fraction(int(cs[j])) if j < len(cs) else Fraction(0))
        out.append(row)
    return out


def roots_t0(Q0):
    poly = sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * u ** k for k, c in enumerate(Q0)), u)
    rts = sp.roots(poly)
    if sum(rts.values()) != poly.degree():
        raise ValueError("non-explicit indicial roots")
    out = []
    for r, m in rts.items():
        if not r.is_rational:
            raise ValueError(f"irrational exponent {r}")
        out += [Fraction(int(r.p), int(r.q))] * m
    return sorted(out)


def classes(rts):
    cl = {}
    for s in rts:
        key = s - 2 * (s.numerator // (2 * s.denominator))      # s mod 2 in [0, 2)
        cl.setdefault(key, []).append(s)
    return cl


def matpoly(coeffs, a, K, p):
    """Q(theta) on the log-coefficient vector at t-power a: theta v = a v + shift, (theta v)[k] = a v[k] + (k+1) v[k+1]"""
    n = K + 1
    A = [[0] * n for _ in range(n)]
    for k in range(n):
        A[k][k] = a % p
        if k + 1 < n:
            A[k][k + 1] = (k + 1) % p
    R = [[0] * n for _ in range(n)]
    for c in reversed(coeffs):                                 # Horner: R = R A + c I
        R = [[sum(R[i][m] * A[m][j] for m in range(n)) % p for j in range(n)] for i in range(n)]
        cc = c.numerator * pow(c.denominator, -1, p) % p
        for i in range(n):
            R[i][i] = (R[i][i] + cc) % p
    return R


def kernel_dim(Q, s0, N, K, p):
    s0p = s0.numerator * pow(s0.denominator, -1, p) % p
    n = K + 1
    nunk = (N + 1) * n
    rows = []
    for ne in range(N + 1):
        block = [[0] * nunk for _ in range(n)]
        for j in range(min(ne, len(Q) - 1) + 1):
            if all(c == 0 for c in Q[j]):
                continue
            a = (s0p + 2 * (ne - j)) % p
            Mj = matpoly(Q[j], a, K, p)
            col0 = (ne - j) * n
            for i in range(n):
                for m in range(n):
                    block[i][col0 + m] = (block[i][col0 + m] + Mj[i][m]) % p
        rows += block
    M = flint.nmod_mat(len(rows), nunk, [v for r in rows for v in r], p)
    return nunk - M.rank()


def jordan(L, where, name):
    Lw = L if where == "0" else Lc.at_infinity(L)
    Q = qforms(Lw)
    rts = roots_t0(Q[0])
    r = len(L) - 1
    res = {}
    for key, lst in sorted(classes(rts).items()):
        s0 = min(lst)
        N = int((max(lst) - s0) / 2) + 25
        per_p = []
        for p in PRIMES:
            d = [kernel_dim(Q, s0, N, K, p) for K in range(r)]
            per_p.append(d)
        d = per_p[0]
        blocks = []
        for K in range(r):
            more = d[K] - (d[K - 1] if K else 0)
            nxt = (d[K + 1] - d[K]) if K + 1 < r else 0
            blocks += [K + 1] * (more - nxt)
        res[str(key)] = {"t_exponents": [str(s) for s in lst], "d_K (log degree <= K)": d,
                         "primes_agree": per_p[0] == per_p[1], "blocks": sorted(blocks, reverse=True)}
    print(f"{name} at {where}: " + "; ".join(f"class {k}: exps {v['t_exponents']} d={v['d_K (log degree <= K)']} "
                                             f"blocks {v['blocks']} agree={v['primes_agree']}" for k, v in res.items()),
          flush=True)
    return res


def theta_x_op(expr_in_thx):
    """an operator given as a polynomial in theta_x = u/2 times x^0 (pure indicial), plus optional x-terms"""
    return G.from_Q([sp.expand(sp.Poly(e.subs(sp.Symbol("T"), u / 2) * 2 ** 4, u).as_expr()) for e in expr_in_thx])


def main():
    T = sp.Symbol("T")
    out = {"controls": {}, "targets": {}}
    ctl = {
        "theta_x^2 (theta_x - 1)^2 [2+2 at 0]": [T ** 2 * (T - 1) ** 2],
        "theta_x^3 (theta_x - 1) [3+1 at 0]": [T ** 3 * (T - 1)],
        "AESZ #1 [4 at 0; 1+1+1+1 at inf]": [T ** 4, -3125 * (T + sp.Rational(1, 5)) * (T + sp.Rational(2, 5))
                                              * (T + sp.Rational(3, 5)) * (T + sp.Rational(4, 5))],
    }
    for name, Qs in ctl.items():
        L = theta_x_op(Qs)
        out["controls"][name] = {"0": jordan(L, "0", name)}
        if len(Qs) > 1:
            out["controls"][name]["inf"] = jordan(L, "inf", name)
    RV = __import__("pf_10_rvmap").RV
    pieces = json.loads((DATA / "pf_21_pieces.json").read_text())
    Q4 = G.parse_op(pieces["L15 / R_S"]["op"])
    for name, L in (("R_V", RV), ("L15 / R_S", Q4)):
        out["targets"][name] = {w: jordan(L, w, name) for w in ("0", "inf")}
    (DATA / "pf_27_jordan.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
