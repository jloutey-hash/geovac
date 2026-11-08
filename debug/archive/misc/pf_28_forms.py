"""Irreducibility of R_V and the type of its invariant form, by exterior and symmetric squares (mod two primes).

V = Sol(L) at a regular point x0, basis y_1..y_4 (Taylor series mod p).  Wedge map: y_i ^ y_j -> w_ij =
y_i theta y_j - theta y_i y_j (theta = theta_t = 2 x d/dx).  Sym map: y_i y_j -> y_i y_j.
For each map we report
  dep   = dim of its kernel (identically vanishing combinations): a monodromy-invariant subspace;
  hyp   = for each exponent-admissible hyperexponential shape  s = x^F0 prod_a (x - a)^E_a P(x), deg P <= N,
          the number of independent P with s in the image (P != 0).
Order-two right factor => its theta-Wronskian is a NONZERO hyperexponential element of the wedge image
(a decomposable vector never maps to 0), so wedge hyp = 0 over every admissible shape excludes order-two factors.
Order-one right factors: hyp of V itself; order-three right factors: order-one right factors of the adjoint.
Negatives are rigorous: mod-p rank <= rank over Q, and truncation only drops equations, so a solution over Q
always survives.  Positives are checked at two primes and two base points.
Invariant form (irreducible V): a one-dimensional kernel or hyperexponential line of the sym map = a symmetric form
on V* preserved up to a character; of the wedge map = an alternating one.

Exponent bounds (x-exponents; t-exponents halved at 0 and infinity): order one, class minima; sym, sums
e_i + e_j (i <= j); wedge, e_i + e_j (i < j) at 0 and infinity and e_i + e_j - 1 at finite a.

Controls: AESZ #1 (symplectic, CY: wedge dep 1, sym dep 0); A (x) B with A = Legendre theta_x^2 - x (theta_x+1/2)^2,
B = A(x/9) (orthogonal tensor: sym dep 1, wedge dep 0, irreducible); A (+) B (planted order-two factor: wedge hyp > 0).
usage: python pf_28_forms.py"""
from __future__ import annotations
import itertools, json, sys
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
NT = 100
NTV = NT - 10          # rows used in ranks: theta and products are exact below NT - 1; keep a margin


def fm(q, p):
    q = Fraction(q)
    return q.numerator % p * pow(q.denominator % p, -1, p) % p


def mul(a, b, p):
    return [int(v) for v in (flint.nmod_poly(a, p) * flint.nmod_poly(b, p)).coeffs()][:NT] + [0] * NT


def trunc(a):
    return (list(a) + [0] * NT)[:NT]


def theta(s, x0, p):
    """theta_t = 2 (x0 + h) d/dh on a power series in h"""
    d = [(n + 1) * s[n + 1] % p for n in range(NT - 1)] + [0]
    x0p = fm(x0, p)
    return [(2 * (x0p * d[n] + (d[n - 1] if n else 0))) % p for n in range(NT)]


def basis(L, x0, p):
    return [trunc([0] * j + G.frob_mod(L, Fraction(x0), j, NT, p)[: NT - j]) for j in range(len(L) - 1)]


def binom_series(c, E, p):
    """(1 + c h)^E mod p"""
    out, coef = [], 1
    Ep, cp = fm(E, p), fm(c, p)
    for n in range(NT):
        out.append(coef * pow(cp, n, p) % p)
        coef = coef * ((Ep - n) % p) % p * pow(n + 1, -1, p) % p
    return out


def rank(vecs, p):
    if not vecs:
        return 0
    M = flint.nmod_mat(len(vecs), NTV, [v % p for vec in vecs for v in vec[:NTV]], p)
    return M.rank()


def hyp_count(funcs, shape, x0, p):
    """number of independent P (deg <= N) with prefactor * P in span(funcs)"""
    F0, E, N = shape
    pre = binom_series(Fraction(1) / Fraction(x0), F0, p)
    for a, Ea in E.items():
        pre = trunc(mul(pre, binom_series(Fraction(1) / (Fraction(x0) - a), Ea, p), p))
    cols = [trunc([0] * k + pre[: NT - k]) for k in range(N + 1)]
    r_all = rank(cols + funcs, p)
    r_f = rank(funcs, p)
    # dim{P : pre*P in span} = (N+1) - (rank(all) - rank(funcs))
    return (N + 1) - (r_all - r_f)


def class_mins(vals):
    out = {}
    for e in vals:
        key = e - (e.numerator // e.denominator)
        if key not in out or e < out[key]:
            out[key] = e
    return sorted(out.values())


def class_maxs(vals):
    out = {}
    for e in vals:
        key = e - (e.numerator // e.denominator)
        if key not in out or e > out[key]:
            out[key] = e
    return out


def shapes(ex, kind):
    """ex = {'0': [...], 'inf': [...growth...], a: [...]} x-exponents; returns admissible (F0, {a: E_a}, N)"""
    def comb(vals, k):
        if k == "one":
            return list(vals)
        if k == "sym":
            return [vals[i] + vals[j] for i in range(len(vals)) for j in range(i, len(vals))]
        return [vals[i] + vals[j] for i in range(len(vals)) for j in range(i + 1, len(vals))]
    f0s = class_mins(comb(ex["0"], kind))
    mus = class_maxs(comb(ex["inf"], kind))
    fin = [a for a in ex if a not in ("0", "inf")]
    per = []
    for a in fin:
        vals = comb(ex[a], kind)
        if kind == "wedge":
            vals = [v - 1 for v in vals]
        per.append(class_mins(vals))
    out = []
    for f0 in f0s:
        for combo in itertools.product(*per):
            s = f0 + sum(combo, Fraction(0))
            for key, mu in mus.items():
                N = mu - s
                if N >= 0 and N.denominator == 1:
                    out.append((f0, dict(zip(fin, combo)), int(N)))
    return out


def analyse(name, sols_fn, ex, x0s=(Fraction(2), Fraction(5))):
    res = {}
    for x0 in x0s:
        for p in PRIMES:
            Y = sols_fn(x0, p)
            n = len(Y)
            W = [trunc([(a - b) % p for a, b in zip(mul(Y[i], theta(Y[j], x0, p), p),
                                                     mul(theta(Y[i], x0, p), Y[j], p))])
                 for i in range(n) for j in range(i + 1, n)]
            S = [trunc(mul(Y[i], Y[j], p)) for i in range(n) for j in range(i, n)]
            r = {"wedge dep": len(W) - rank(W, p), "sym dep": len(S) - rank(S, p)}
            for kind, funcs in (("one", Y), ("wedge", W), ("sym", S)):
                hs = [(sh, hyp_count(funcs, sh, x0, p)) for sh in shapes(ex, kind)]
                r[f"{kind} shapes"] = len(hs)
                r[f"{kind} hyp"] = [[str(sh[0]), {str(k): str(v) for k, v in sh[1].items()}, sh[2], c]
                                    for sh, c in hs if c]
            res[f"x0={x0},p={p}"] = r
    vals = list(res.values())
    same = all(v == vals[0] for v in vals)
    print(f"{name}: {vals[0]}  (same at both primes and base points: {same})", flush=True)
    return {"by_run": res, "agree": same}


def xexps(L):
    """x-exponents of an operator in the G convention at 0, infinity and its finite singular points"""
    lc = sp.Poly([int(c) for c in reversed(L[-1].coeffs())], sp.Symbol("x"))
    pts = [Fraction(int(r.p), int(r.q)) for r in sp.roots(lc).keys() if r != 0]
    ex = {"0": [Fraction(e) / 2 for e in Lc.exps_t0(L)[0]], "inf": [Fraction(m) / 2 for m in Lc.exps_inf_growth(L)[0]]}
    for a in pts:
        ex[a] = [Fraction(e) for e in Lc.exps_finite_t(L, a)[0]]
    return ex


def main():
    out = {}
    A = G.from_Q([u ** 2, -(u + 1) ** 2])
    B = G.from_Q([9 * u ** 2, -(u + 1) ** 2])
    Q5 = G.from_Q([u ** 4, sp.expand(-5 * (5 * u + 2) * (5 * u + 4) * (5 * u + 6) * (5 * u + 8))])
    h = Fraction(-1, 2)
    exA = {"0": [0, 0], "inf": [h, h], Fraction(1): [0, 0], Fraction(9): [0, 1]}
    exB = {"0": [0, 0], "inf": [h, h], Fraction(1): [0, 1], Fraction(9): [0, 0]}
    F = lambda v: [Fraction(x) for x in v]
    ex_sum = {k: F(exA[k] + exB[k]) for k in exA}
    ex_ten = {k: sorted(F([a + b for a in exA[k] for b in exB[k]])) for k in exA}
    out["control AESZ #1"] = analyse("control AESZ #1", lambda x0, p: basis(Q5, x0, p), xexps(Q5))
    out["control A (x) B"] = analyse("control A (x) B", lambda x0, p: [trunc(mul(f, g, p)) for f in basis(A, x0, p)
                                                                      for g in basis(B, x0, p)], ex_ten)
    out["control A (+) B"] = analyse("control A (+) B", lambda x0, p: basis(A, x0, p) + basis(B, x0, p), ex_sum)
    RV = __import__("pf_10_rvmap").RV
    print("R_V exponents:", {str(k): [str(e) for e in v] for k, v in xexps(RV).items()}, flush=True)
    out["R_V"] = analyse("R_V", lambda x0, p: basis(RV, x0, p), xexps(RV))
    RVa = G.adjoint(RV)
    print("R_V* exponents:", {str(k): [str(e) for e in v] for k, v in xexps(RVa).items()}, flush=True)
    out["R_V adjoint"] = analyse("R_V adjoint", lambda x0, p: basis(RVa, x0, p), xexps(RVa))
    (DATA / "pf_28_forms.json").write_text(json.dumps(out, indent=1, default=str))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
